# -*- coding: utf-8 -*-
"""
한반도 LMM (Local Magnetic Model) 구축 스크립트
=================================================

오석훈 교수 설명자료(lmm_korea_explainer_kr.pdf)의 4-층 결합 구조를
현재 보유 자료로 구현한다.

    B_LMM(r,t) = B_IGRF + B_Regional + B_Crustal + B_External

| Layer    | 자료원                                  | 보유 |
|----------|-----------------------------------------|------|
| Core     | IGRF-14 (ppigrf / IGRF14.shc)           |  O   |
| Regional | 지자기측량 성과정리(22_25).xlsx 절대측정  |  O   |
| Crustal  | KIGAM 자력이상도 mag_1982-2018_1.5min    |  O   |
| External | CYG 1분 자료 (INTERMAGNET)               |  X   |

External 층이 없으므로 본 모델은 **정온시(Kp<=2) 기준 baseline 모델**이다.
CYG dmin 자료를 확보하면 add_external_layer() 를 통해 확장한다.

출력:
    docs/data/lmm_model.json   — 웹 계산기용 계수 묶음
    docs/output/YYYYMMDD_HHMMSS_LMM_validation.csv — 검증표
"""

import json
import math
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd
from ppigrf import igrf

BASE = Path(__file__).parent
DATA = BASE / "data"
DOCS_DATA = BASE / "docs" / "data"
DOCS_OUT = BASE / "docs" / "output"

SURVEY_XLSX = DATA / "지자기측량 성과정리(22_25).xlsx"
KIGAM_DAT = DATA / "mag_1982-2018_1.5min_ed.dat"

EPOCH = 2027.0  # 목표 Epoch

# 2019 야장 성과를 Regional 층에 포함할지 (lmm_verify2019.py 검증 통과분)
INCLUDE_2019 = True

# Regional 층 다항식 차수. 측점 26개 기준 2차(6계수)까지가 안전.
REGIONAL_DEGREE = 2

# 지각층을 편각·복각에도 기여시킬지 (스칼라 dF -> 벡터 b 복원, 아래 crustal_vector).
# False 면 종전처럼 지각층은 F 에만 더해진다.
CRUSTAL_VECTOR = True

# 정규화 기준점 (한반도 중심)
LAT0, LON0 = 36.0, 127.5

KM_PER_DEG = 111.32


# ---------------------------------------------------------------- 자료 적재
def load_survey_points() -> pd.DataFrame:
    """지표 절대측정(F·D·I) 성과를 정리해 반환."""
    df = pd.read_excel(SURVEY_XLSX, sheet_name="Sheet1")
    df = df.rename(
        columns={
            "실": "lat",
            "실.1": "lon",
            "실.2": "D_deg",
            "실.3": "I_deg",
            "표고": "elev_m",
            "총자력": "F_nT",
            "관측연도": "year",
            "도엽명": "name",
        }
    )
    df["name"] = df["name"].ffill()
    df = df[["name", "year", "lat", "lon", "elev_m", "D_deg", "I_deg", "F_nT"]]
    df = df.dropna(subset=["lat", "lon", "D_deg", "I_deg", "F_nT"])

    # 양산/언양처럼 동일 지점이 두 이름으로 중복 기재된 행 제거.
    # 좌표가 소수점 끝자리만 다르므로 약 10 m 격자로 반올림해 비교한다.
    df = df.assign(
        _k_lat=df["lat"].round(4),
        _k_lon=df["lon"].round(4),
    ).drop_duplicates(subset=["_k_lat", "_k_lon", "year"], keep="first")
    df = df.drop(columns=["_k_lat", "_k_lon"])

    # 대표 관측일시 — 야장에서 복원한 실제 일시를 우선한다.
    # 성과표에는 연도만 있어 종전에는 7월 1일로 대체했는데, 그러면 ① 영년변화
    # 환산이 최대 반년 어긋나고 ② 관측시각을 몰라 External 보정을 걸 수 없었다.
    # 2026-08 지리원 야장을 확보해 (측점, 연도) 로 맞춰 실제 일시를 붙인다.
    df["date"] = df["year"].apply(lambda y: dt.datetime(int(y), 7, 1))
    df["date_src"] = "연도대체(7/1)"
    df = attach_fieldbook_datetime(df)
    df = apply_external_correction(df)
    df = apply_quality_filter(df)
    return df.reset_index(drop=True)


# ── Regional 성과 신뢰도 필터 (audit_regional_reliability.py, 2026-08-19) ──────
#
# 재방문 잔여(관측 ΔD − IGRF 영년변화)로 방문별 신뢰도를 판정했다. 같은 표석을
# 다시 재면 0 이어야 하므로, 크게 벗어난 방문은 그 회차의 방위 기준이 틀렸다는 뜻이다.
#
#   불일치   연속 구간의 잔여 부호가 뒤집혀 **그 방문이 원인으로 특정**된 것
#   판별불가 2회 방문뿐이라 둘 중 어느 쪽이 틀렸는지 못 가리는 것
#
# 내부 산술(재계산차 0.16″)과 정밀도(세션 산포 1.41′)는 통과했으므로 측정·계산의
# 문제가 아니다. 방문마다 마크 참방위각이 바뀌는 것과 맞물린 **기준의 문제**다.
REGIONAL_QUALITY = {
    ("부안", 2023): "불일치",     # 앞구간 +62.2′ / 뒷구간 −70.6′ 로 부호 반전
    ("포천", 2023): "불일치",     # −32.2′ / +22.3′ 로 부호 반전
    ("가야", 2022): "판별불가",   # 2회 방문 · 잔여 −35.5′
    ("가야", 2024): "판별불가",
    ("거제", 2022): "판별불가",   # 야장 2019~2024 잔여 −19.7′ (주의)
    ("거제", 2024): "판별불가",
}
#   "none"         : 전건 사용 (종전 동작)
#   "strict"       : 원인이 특정된 「불일치」만 제외
#   "conservative" : 「판별불가」까지 제외
#
# 채택 = "strict". 감사가 **모델 잔차를 보지 않고** 야장 원시 판독값과 IGRF
# 영년변화만으로 지목한 2개 방문을 빼면 LOO 편각이 0.7691 -> 0.5675° 로 26%
# 개선된다(F 60.26 -> 58.41 nT). 순환 논증이 아니라는 점이 중요하다.
# ⚠ 대가 — 복각이 0.1787 -> 0.2600° 로 악화되고, 부안이 성과표에서 2023 한 건뿐이라
#   측점이 17 -> 16 으로 줄어 서해안 공간 커버리지를 잃는다.
REGIONAL_FILTER = "strict"


def apply_quality_filter(df: pd.DataFrame) -> pd.DataFrame:
    if REGIONAL_FILTER == "none":
        df["quality"] = "미적용"
        return df
    drop = {"strict": {"불일치"},
            "conservative": {"불일치", "판별불가"}}[REGIONAL_FILTER]
    key = list(zip(df["name"], df["year"].astype(int)))
    tag = [REGIONAL_QUALITY.get(k, "정상") for k in key]
    df = df.assign(quality=tag)
    keep = ~df["quality"].isin(drop)
    n = int((~keep).sum())
    if n:
        rm = df.loc[~keep, ["name", "year", "quality"]]
        print(f"  [신뢰도] {REGIONAL_FILTER} - {n}행 제외: "
              + ", ".join(f"{r['name']}{int(r['year'])}({r['quality']})"
                          for _, r in rm.iterrows()))
    return df[keep]


# ── ④ External 층 — 성과표 값에 외부장 보정 적용 ──────────────────────────────
#
# 성과표가 **원시 세션평균**임이 확정됐으므로(audit_survey_table.py, 2026-08-19)
# 뺄 외부장이 그대로 남아 있다. lmm_external.py 가 세션별 편차를 산출해 두었고,
# 성과값이 세션 평균이므로 보정도 **세션 평균**으로 걸어야 짝이 맞는다.
#
#   "none"          : 미적용 (종전 동작)
#   "subtract"      : 전 세션 평균 보정량을 뺀다
#   "subtract_quiet": Kp<=2 세션만으로 평균한 보정량을 뺀다
#   "drop_storm"    : 보정은 하지 않고, 교란시(Kp>2) 관측 성과를 표본에서 뺀다
#   "subtract_D"    : **편각만** 보정한다 (총자력은 손대지 않는다)
#
# ⚠ 균일 V 근사의 한계가 남아 있다 — CYG 한 곳의 변동을 전국에 동일 적용한다.
#
# 채택 = "subtract_D" (2026-08-19). 신뢰도 감사로 편각 잔차가 35→20.7′ 로 줄자
# **성분별로 결과가 갈렸다** — 같은 보정을 걸어도 편각은 개선되고 총자력은 악화된다.
#   D : 보정량 5.4′ / 잔차 20.7′ → 비중이 커서 이득 (LOO 0.5675 -> 0.5403°)
#   F : 보정량 38 nT / 잔차 약 190 nT → 지각장 불일치가 지배해 잡음만 더해짐 (+8.3)
# 그래서 **편각에만** 적용한다. 감사 전 표본에서는 D 도 악화됐었다(+0.079) —
# 방위 기준 오류가 외부장 신호를 덮고 있었기 때문이다.
# 2026-08-19 갱신 — 관측소 4소 공간보간(multi) + 정온세션 한정 + 편각 전용.
#   생산 설정 LOO D:  미적용 0.5899 / subtract_D·cyg 0.3653 /
#                     subtract_D·multi 0.4028 / **subtract_quiet_D·multi 0.3232**
# 다중 관측소는 «정온 세션에 한정할 때만» 이득이 난다. 교란시에는 4소 평면적합이
# 크게 발산해(|dF| 최대 567 nT) 잡음을 키운다 — 공간투영 근사가 폭풍 중에 먼저
# 깨진다는 뜻이다. ⚠ 8개 설정 중 LOO 최소를 고른 것이므로 0.365→0.323 차이는
# 선택편의를 포함한다. 방향(정온 한정·공간보간)이 물리적으로 타당한 것이 근거이고,
# 차이 자체를 성과로 읽지 말 것.
EXTERNAL_MODE = "subtract_quiet_D"

# 보정량 자료원 — "cyg"(청양 단독, 균일 V 근사) | "multi"(관측소 4소 공간보간).
# 자문 권고는 대체관측소 활용이므로 자료가 있으면 multi 를 쓴다. 파일이 없으면
# 자동으로 cyg 로 내려간다(조용히 무보정이 되지 않도록).
EXTERNAL_SOURCE = "multi"
EXTERNAL_CSV_CYG = BASE / "docs" / "data" / "external_corrections.csv"
EXTERNAL_CSV_MULTI = BASE / "docs" / "data" / "external_corrections_multi.csv"
EXTERNAL_CSV = (EXTERNAL_CSV_MULTI
                if EXTERNAL_SOURCE == "multi" and EXTERNAL_CSV_MULTI.exists()
                else EXTERNAL_CSV_CYG)


def apply_external_correction(df: pd.DataFrame) -> pd.DataFrame:
    """(측점, 연도) 세션 평균 외부장 편차를 성과값에서 뺀다."""
    if EXTERNAL_MODE == "none" or not EXTERNAL_CSV.exists():
        df["ext_src"] = "미적용"
        return df
    ec = pd.read_csv(EXTERNAL_CSV, encoding="utf-8-sig")
    ec = ec[ec["상태"] == "정상"].copy()
    ec["연도"] = pd.to_datetime(ec["날짜"]).dt.year
    if EXTERNAL_MODE in ("subtract_quiet", "subtract_quiet_D"):
        # 교란시 세션의 보정량은 신뢰하지 않는다(균일·평면 근사가 폭풍 중에 가장
        # 크게 깨진다). 값을 만지는 대신 그 세션의 보정만 쓰지 않는 방식.
        ec = ec[ec["Kp"] <= FB2019_MAX_KP]
    if ec.empty:
        df["ext_src"] = "미적용(표본없음)"
        return df

    g = (ec.groupby(["측점", "연도"])
           .agg(dF=("dF", "mean"), dD=("dD_arcmin", "mean"),
                kp_max=("Kp", "max"), n=("dF", "size")).reset_index())
    out = df.merge(g, left_on=["name", "year"], right_on=["측점", "연도"],
                   how="left")
    hit = out["dF"].notna()

    if EXTERNAL_MODE == "drop_storm":
        # 보정하지 않고 교란시 성과를 표본에서 뺀다(값을 만지지 않는 방식)
        storm = hit & (out["kp_max"] > FB2019_MAX_KP)
        out["ext_src"] = np.where(storm, "교란배제", "유지")
        print(f"  [External] drop_storm - 교란시 성과 {int(storm.sum())}행 배제")
        out = out[~storm]
    else:
        # ⚠ 성분별로 결과가 갈린다. 감사 후 표본에서 편각은 개선되지만
        #   총자력은 악화된다 — F 잔차는 외부장(38 nT)이 아니라 **지각장 불일치**
        #   (약 190 nT)가 지배하므로, 균일한 시간 보정을 빼면 잡음만 더해진다.
        #   편각 잔차는 20.7′ 이고 보정량이 5.4′ 라 비중이 커서 이득이 난다.
        if EXTERNAL_MODE not in ("subtract_D", "subtract_quiet_D"):
            out.loc[hit, "F_nT"] = out.loc[hit, "F_nT"] - out.loc[hit, "dF"]
        out.loc[hit, "D_deg"] = out.loc[hit, "D_deg"] - out.loc[hit, "dD"] / 60.0
        out["ext_src"] = np.where(hit, EXTERNAL_MODE, "미적용")
        print(f"  [External] {EXTERNAL_MODE} - {int(hit.sum())}/{len(out)}행 보정 "
              f"(|dF| 평균 {out.loc[hit, 'dF'].abs().mean():.1f} nT · "
              f"|dD| 평균 {out.loc[hit, 'dD'].abs().mean():.2f}')")
    return out.drop(columns=[c for c in ("측점", "연도", "dF", "dD", "kp_max", "n")
                             if c in out])


# 야장에서 복원한 세션표 (lmm_fieldbook.py 산출). 없으면 종전 7/1 대체로 동작한다.
FIELDBOOK_SESSIONS = BASE / "docs" / "data" / "fieldbook_sessions.csv"


def attach_fieldbook_datetime(df: pd.DataFrame) -> pd.DataFrame:
    """(측점, 연도) 로 야장 세션을 찾아 대표 관측일시를 붙인다.

    한 해에 여러 세션이 있으므로 **세션 시작시각의 중앙값**을 대표로 쓴다.
    성과 자체가 세션 평균이라 특정 세션 시각을 고르는 것은 근거가 없고,
    중앙값이 외부장 보정의 기준시각으로 가장 무난하다.
    """
    if not FIELDBOOK_SESSIONS.exists():
        print("  [야장] 세션표 없음 - 연도 대체일 유지"
              " (python lmm_fieldbook.py 로 생성)")
        return df
    fb = pd.read_csv(FIELDBOOK_SESSIONS, encoding="utf-8-sig")
    fb = fb.dropna(subset=["측점", "날짜"])
    fb["날짜"] = pd.to_datetime(fb["날짜"])
    fb["연도"] = fb["날짜"].dt.year
    fb = fb[fb["시작"].notna()]
    if fb.empty:
        return df

    def _dt(r):
        return pd.Timestamp(f"{r['날짜'].date()} {str(r['시작'])[:8]}")

    fb["일시"] = fb.apply(_dt, axis=1)
    rep_dt = (fb.groupby(["측점", "연도"])["일시"]
                .median().rename("fb_dt").reset_index())

    out = df.merge(rep_dt, left_on=["name", "year"],
                   right_on=["측점", "연도"], how="left")
    hit = out["fb_dt"].notna()
    out.loc[hit, "date"] = out.loc[hit, "fb_dt"]
    out.loc[hit, "date_src"] = "야장복원"
    out = out.drop(columns=[c for c in ("측점", "연도", "fb_dt") if c in out])
    print(f"  [야장] 관측일시 복원 {hit.sum()}/{len(out)}행"
          f" - 나머지는 연도 대체일(7/1) 유지")
    return out


# 2019 야장 투입 정책 (lmm_verify2019.py 검증 결과 반영)
#   · 성산 — 관측시각 미기입 + 원시 판독값 배열이 달라 역산 검증 불가 -> 제외
#   · 겹치는 측점은 성과표 좌표를 기준으로 삼는다 (미원은 353 m 불일치 확인됨)
#   · 신규 측점(장흥·남지) 표고는 개략값. IGRF 표고기울기 -26 nT/km 이므로
#     수십 m 오차는 1~2 nT 수준으로 영향이 미미하다.
FB2019_EXCLUDE = ("성산",)
FB2019_NEW_SITE_COORD = {
    "장흥": (34 + 37 / 60 + 38 / 3600, 126 + 58 / 60 + 40 / 3600, 100.0),
    "남지": (35 + 26 / 60 + 10 / 3600, 128 + 28 / 60 + 25 / 3600, 50.0),
}

# ────────────────────────────────────────────────────────────────────────
# ④ External 층 정책
#
# ⚠️ 잠정 전제 (2026-07, 미확인) — 성과에 일변화 보정이 적용되지 **않았다**고 본다.
#
# 「지구물리측량 작업규정」 제21조는 상시 기준점 대비 일변화 보정을 법정 절차로
# 정하고 있으나, 실제 적용 여부가 확인되지 않았다. 다음 정황으로 미적용으로 본다:
#   · 동일 측점 재방문 F 산포가 최대 138 nT — 보정되었다면 나올 수 없는 크기
#   · 성과표에 관측시각이 없음 (보정하려면 시각이 필요하므로 산출 근거가 남지 않음)
# 확인되면 이 전제와 아래 플래그를 수정할 것.
ASSUME_NO_DIURNAL_CORRECTION = True

# CYG 외부장 값을 직접 차감(V-subtraction)할지. **적용하지 않는다.** 두 가지 이유:
#   ① 「V 가 전역 균일」이라는 1차 근사가 CYG 로부터 76~229 km 떨어진 측점에서
#      성립하지 않는다. 실측 LOO D 가 0.751° -> 0.858° 로 악화된다.
#      (설명자료가 NOC 공간투영을 요구하는 이유)
#   ② 야장의 F 는 **일자별 단일값**이다(거제·장흥은 측점당 1개). 세션 시각 기준
#      보정량을 일자 대표값에 빼는 것은 원리적으로 성립하지 않는다.
#      제19조①6 은 세션마다 총자기장 측정을 요구하나 야장 양식이 이를 담지 않는다.
FB2019_APPLY_EXTERNAL = False

# 자기폭풍 중 관측된 세션은 배제한다(Kp). 설명자료가 Kp<=2 정온시간대 관측을
# 권고하며, 미원 2019-05-14 두 세션이 Kp 5.7~6.3 자기폭풍 중 관측되었다.
# 이는 값을 조작하지 않고 오염된 표본만 제외하므로 위 ①②의 문제가 없다.
FB2019_MAX_KP = 2.0


def external_corrections_2019(raw: pd.DataFrame) -> pd.DataFrame:
    """
    2019 세션별 CYG 외부장 보정량과 Kp 를 계산해 raw 에 덧붙인다.

    반환 컬럼: cor_D(deg), cor_F(nT), Kp — 측정값에서 **빼야 할** 양.
    """
    import numpy as np

    from lmm_cyg import fetch_kp, fetch_range
    from lmm_external import derive_dif, quiet_baseline, to_utc

    r = raw[raw["시작"].notna()].copy()
    if r.empty:
        return raw.assign(cor_D=np.nan, cor_F=np.nan, Kp=np.nan)

    days = sorted(set(r["날짜"]))
    parts = [fetch_range(d - dt.timedelta(days=10), d + dt.timedelta(days=10),
                         quiet=True) for d in days]
    cyg = pd.concat(parts, ignore_index=True)
    cyg["time"] = pd.to_datetime(cyg["time"], utc=True)
    cyg = cyg.sort_values("time").drop_duplicates("time")
    kp = fetch_kp()
    cygd = derive_dif(cyg).set_index("time")

    from lmm_cyg import kp_for

    base_cache, rows = {}, []
    for _, s in raw.iterrows():
        if pd.isna(s["시작"]):
            rows.append((np.nan, np.nan, np.nan))
            continue
        day = s["날짜"]
        if day not in base_cache:
            base_cache[day] = quiet_baseline(cyg, kp, day)
        b = base_cache[day]
        t0 = to_utc(day, s["시작"])
        seg = cygd[(cygd.index >= t0 - pd.Timedelta(minutes=5))
                   & (cygd.index <= t0 + pd.Timedelta(minutes=20))].dropna(subset=["X"])
        if b is None or seg.empty:
            rows.append((np.nan, np.nan, np.nan))
            continue
        vX, vY, vZ = seg.X.mean() - b["X"], seg.Y.mean() - b["Y"], seg.Z.mean() - b["Z"]
        bH = np.hypot(b["X"], b["Y"])
        bF = np.hypot(bH, b["Z"])
        mH = np.hypot(b["X"] + vX, b["Y"] + vY)
        mF = np.hypot(mH, b["Z"] + vZ)
        dD = np.degrees(np.arctan2(b["Y"] + vY, b["X"] + vX)
                        - np.arctan2(b["Y"], b["X"]))
        rows.append((dD, mF - bF, float(kp_for([t0], kp).iloc[0])))

    raw = raw.copy()
    raw[["cor_D", "cor_F", "Kp"]] = pd.DataFrame(rows, index=raw.index)
    return raw


def load_fieldbook_2019(exclude=FB2019_EXCLUDE, apply_external=None,
                        max_kp=None) -> pd.DataFrame:
    """
    2019 야장 성과를 성과표와 동일한 스키마로 반환.

    세션 단위 성과를 그대로 쓰면 같은 날 여러 세션이 과대 가중되므로
    측점·일자별로 평균해 1건으로 축약한다.

    apply_external : CYG 외부장 보정 적용 여부 (None 이면 모듈 기본값)
    max_kp         : 이 값을 넘는 Kp 세션을 배제 (None 이면 배제 안 함)
    """
    from lmm_verify2019 import collect_raw

    apply_external = (FB2019_APPLY_EXTERNAL if apply_external is None
                      else apply_external)
    max_kp = FB2019_MAX_KP if max_kp is None else max_kp

    raw = collect_raw()
    # ⚠ collect_raw() 는 2026-08 부터 지리원 야장 전체(2019~2025)를 훑는다.
    #   그 확장은 **신뢰도 감사용**이고, 여기서는 이름 그대로 2019 분만 쓴다.
    #   걸러내지 않으면 성과표와 같은 (측점,연도)가 중복 투입된다.
    raw = raw[pd.to_datetime(raw["날짜"]).dt.year == 2019]
    raw = raw[~raw["측점"].isin(exclude)]
    raw = raw.dropna(subset=["D_야장", "I_야장", "F_야장"])

    if apply_external or max_kp is not None:
        raw = external_corrections_2019(raw)
        if max_kp is not None:
            keep = raw["Kp"].isna() | (raw["Kp"] <= max_kp)
            raw = raw[keep]
        if apply_external:
            raw = raw.assign(
                D_야장=raw["D_야장"] - raw["cor_D"].fillna(0.0),
                F_야장=raw["F_야장"] - raw["cor_F"].fillna(0.0),
            )

    # 겹치는 측점은 성과표 좌표, 신규 측점은 관측현황 좌표를 쓴다
    ref = load_survey_points().groupby("name").first()
    rows = []
    for (site, day), g in raw.groupby(["측점", "날짜"]):
        if site in ref.index:
            lat, lon, elev = (ref.loc[site, "lat"], ref.loc[site, "lon"],
                              ref.loc[site, "elev_m"])
        elif site in FB2019_NEW_SITE_COORD:
            lat, lon, elev = FB2019_NEW_SITE_COORD[site]
        else:
            continue
        rows.append({
            "name": site,
            "year": day.year,
            "lat": lat, "lon": lon, "elev_m": elev,
            "D_deg": g["D_야장"].mean(),
            "I_deg": g["I_야장"].mean(),
            "F_nT": g["F_야장"].mean(),
            "date": dt.datetime.combine(day, dt.time(12, 0)),
        })
    return pd.DataFrame(rows)


def load_all_points(include_2019=True, **fb_kw) -> pd.DataFrame:
    """성과표(2022~2025) + 선택적으로 2019 야장 성과."""
    pts = load_survey_points()
    if not include_2019:
        return pts
    fb = load_fieldbook_2019(**fb_kw)
    if fb.empty:
        return pts
    return pd.concat([pts, fb], ignore_index=True)


def repeatability(pts, residuals):
    """
    동일 측점 재방문 성과의 산포를 집계한다.

    현장 배리오미터·CYG 보정이 없으면 방문 시각의 외부장 잡음이 그대로
    성과에 남는다(설명자료 p.9-10). 재방문 산포는 그 잡음의 실측 하한이다.
    """
    t = pts.copy()
    t[["dD", "dI", "dF"]] = residuals[["dD", "dI", "dF"]].values

    rows = []
    for name, g in t.groupby("name"):
        if len(g) < 2:
            continue
        rows.append(
            {
                "측점": name,
                "방문수": len(g),
                "연도": "/".join(str(int(y)) for y in g["year"]),
                "F_산포_nT": float(g["dF"].max() - g["dF"].min()),
                "D_산포_arcmin": float((g["dD"].max() - g["dD"].min()) * 60),
                "I_산포_arcmin": float((g["dI"].max() - g["dI"].min()) * 60),
            }
        )
    return pd.DataFrame(rows).sort_values("F_산포_nT", ascending=False)


# 재방문 산포 허용 한계.
# 2~3년간 영년변화(secular variation)로 설명 가능한 D 변화는 최대 15' 수준이므로,
# 그보다 큰 산포는 방위 기준 오류 등 조사 블런더로 본다.
#
# 참고 — 「지구물리측량 작업규정」(국토지리정보원고시 제2021-2985호) 제20조는
# 편각·복각의 정수차 한계를 30' 로 정한다. 그러나 이는 측정 품질관리 기준이며
# 모델 KPI(0.1°=6')보다 5배 느슨하다. 실제로 30' 를 적용하면 포천(22.3')이
# 포함되어 LOO D-RMS 가 0.751° -> 0.807° 로 악화된다. 따라서 법정 한계보다
# 엄격한 20' 를 쓴다. (법정 기준 충족이 곧 모델 투입 적합을 뜻하지 않는다)
MAX_D_SPREAD_ARCMIN = 20.0
MAX_I_SPREAD_ARCMIN = 15.0


def quality_flags(pts, residuals):
    """재방문 산포로 측점별 성분 신뢰도를 판정한다."""
    rep = repeatability(pts, residuals).set_index("측점")
    names = sorted(pts["name"].unique())

    flags = {}
    for n in names:
        if n in rep.index:
            flags[n] = {
                "D_ok": rep.loc[n, "D_산포_arcmin"] <= MAX_D_SPREAD_ARCMIN,
                "I_ok": rep.loc[n, "I_산포_arcmin"] <= MAX_I_SPREAD_ARCMIN,
                "repeat": True,
            }
        else:
            # 단일 방문 측점은 검증 불가 -> 일단 채택하고 MAD 클리핑에 맡긴다
            flags[n] = {"D_ok": True, "I_ok": True, "repeat": False}
    return flags


def aggregate_sites(pts, residuals):
    """
    측점별로 잔차를 평균해 대표값 1개로 축약한다.

    재방문 평균은 외부장 잡음을 sqrt(n)만큼 낮추는 근사 보정이며,
    CYG 자료 확보 전까지의 차선책이다.
    """
    t = pts.copy()
    t[["dD", "dI", "dF"]] = residuals[["dD", "dI", "dF"]].values

    agg = (
        t.groupby("name")
        .agg(
            lat=("lat", "mean"),
            lon=("lon", "mean"),
            elev_m=("elev_m", "mean"),
            n_visit=("year", "size"),
            dD=("dD", "mean"),
            dI=("dI", "mean"),
            dF=("dF", "mean"),
        )
        .reset_index()
    )

    fl = quality_flags(pts, residuals)
    agg["D_ok"] = agg["name"].map(lambda n: fl[n]["D_ok"])
    agg["I_ok"] = agg["name"].map(lambda n: fl[n]["I_ok"])
    agg["repeat"] = agg["name"].map(lambda n: fl[n]["repeat"])
    return agg


def load_kigam_grid():
    """
    KIGAM 자력이상도(잔차, IGRF 제거됨)를 **균일 간격** 규칙격자로 반환.

    주의: 원자료의 위도축에는 제주해협 구간(33.55°~34.1°)에 0.55° 공백이 있어
    np.unique() 결과가 균일 간격이 아니다. 이를 그대로 축으로 쓰면
    (값-원점)/간격 방식의 색인이 공백 위쪽 전 구간에서 21행씩 어긋난다.
    따라서 결측 위도행을 NaN 으로 채운 완전 균일축을 구성한다.
    """
    a = np.loadtxt(KIGAM_DAT, skiprows=9)
    lon, lat, anom = a[:, 0], a[:, 1], a[:, 2]

    step = 0.025
    lons = np.round(np.arange(lon.min(), lon.max() + step / 2, step), 4)
    lats = np.round(np.arange(lat.min(), lat.max() + step / 2, step), 4)

    grid = np.full((lats.size, lons.size), np.nan)
    ix = np.rint((lon - lons[0]) / step).astype(int)
    iy = np.rint((lat - lats[0]) / step).astype(int)
    grid[iy, ix] = anom

    return lons, lats, grid


# ------------------------------------------------------------- Layer 1 IGRF
def igrf_dif(lat, lon, elev_m, date):
    """
    IGRF-14 를 평가해 (D, I, F, X, Y, Z) 반환. 각도는 도(deg), 자력은 nT.

    date 는 단일 datetime 또는 측점별 datetime 리스트를 받는다.
    ppigrf 는 (날짜 x 측점) 외적을 반환하므로 측점별 날짜는 개별 평가한다.
    """
    lat = np.atleast_1d(np.asarray(lat, float))
    lon = np.atleast_1d(np.asarray(lon, float))
    h_km = np.atleast_1d(np.asarray(elev_m, float)) / 1000.0

    if isinstance(date, (list, tuple, np.ndarray, pd.Series)):
        parts = [igrf(lon[i:i + 1], lat[i:i + 1], h_km[i:i + 1], d)
                 for i, d in enumerate(date)]
        Be = np.concatenate([np.ravel(p[0]) for p in parts])
        Bn = np.concatenate([np.ravel(p[1]) for p in parts])
        Bu = np.concatenate([np.ravel(p[2]) for p in parts])
    else:
        Be, Bn, Bu = igrf(lon, lat, h_km, date)

    X = np.ravel(Bn)          # 북
    Y = np.ravel(Be)          # 동
    Z = -np.ravel(Bu)         # 하향(+)

    H = np.hypot(X, Y)
    F = np.hypot(H, Z)
    D = np.degrees(np.arctan2(Y, X))
    I = np.degrees(np.arctan2(Z, H))
    return D, I, F, X, Y, Z


# ---------------------------------------------------------- Layer 3 Crustal
class CrustalGrid:
    """KIGAM 자력이상 격자의 쌍선형 보간자."""

    def __init__(self, lons, lats, grid):
        self.lons, self.lats, self.grid = lons, lats, grid
        self.dlon = float(np.median(np.diff(lons)))
        self.dlat = float(np.median(np.diff(lats)))

    def __call__(self, lat, lon):
        lat = np.atleast_1d(np.asarray(lat, float))
        lon = np.atleast_1d(np.asarray(lon, float))

        fx = (lon - self.lons[0]) / self.dlon
        fy = (lat - self.lats[0]) / self.dlat
        i0 = np.clip(np.floor(fx).astype(int), 0, self.lons.size - 2)
        j0 = np.clip(np.floor(fy).astype(int), 0, self.lats.size - 2)
        tx = np.clip(fx - i0, 0, 1)
        ty = np.clip(fy - j0, 0, 1)

        g = self.grid
        c00, c10 = g[j0, i0], g[j0, i0 + 1]
        c01, c11 = g[j0 + 1, i0], g[j0 + 1, i0 + 1]

        # NaN(자료 공백) 이웃은 유효 이웃 평균으로 대체
        stack = np.stack([c00, c10, c01, c11])
        w = np.stack(
            [
                (1 - tx) * (1 - ty),
                tx * (1 - ty),
                (1 - tx) * ty,
                tx * ty,
            ]
        )
        valid = ~np.isnan(stack)
        w = np.where(valid, w, 0.0)
        wsum = w.sum(axis=0)
        vals = np.where(valid, np.nan_to_num(stack), 0.0)
        out = np.where(wsum > 0, (w * vals).sum(axis=0) / np.where(wsum > 0, wsum, 1), np.nan)
        return out



def crustal_ref_dir(lats, ref_lat=None, ref_lon=127.5, epoch=None):
    """벡터 역산에 쓰는 기준 주자기장 방향 (l, m, n) — 격자 중앙 위도에서 평가."""
    if ref_lat is None:
        ref_lat = float(np.mean(lats))
    if epoch is None:
        epoch = dt.datetime(int(EPOCH), 1, 1)
    Dm, Im = igrf_dif(np.array([ref_lat]), np.array([ref_lon]),
                      np.array([0.0]), epoch)[:2]
    Dm, Im = float(np.ravel(Dm)[0]), float(np.ravel(Im)[0])
    return {"D": Dm, "I": Im,
            "l": math.cos(math.radians(Im)) * math.cos(math.radians(Dm)),
            "m": math.cos(math.radians(Im)) * math.sin(math.radians(Dm)),
            "n": math.sin(math.radians(Im)),
            "ref_lat": ref_lat, "ref_lon": ref_lon}


def crustal_vector(lons, lats, grid, ref_lat=None, ref_lon=127.5, epoch=None):
    """
    스칼라 자기이상 dF 격자에서 이상벡터 3성분 (북·동·하) 을 복원한다.

    항공자력은 총자력 이상만 준다. 이는 이상벡터를 주자기장 방향
    (l, m, n) 으로 투영한 값이므로 dF = l*b_north + m*b_east + n*b_down 이다.
    포텐셜장은 조화함수라 파수영역에서 이 관계가 닫히고, 역산하면 3성분을
    모두 되찾을 수 있다 — 지각층이 편각·복각에 기여하는 통로가 이것이다.

    반환: (b_north, b_east, b_down)  격자와 같은 shape, 단위 nT.

    주의: 파수 0 성분(격자 평균)은 복원되지 않아 상수 오프셋이 남는다.
          Regional 상수항이 흡수하므로 무해하다.
    """
    from scipy.signal.windows import tukey

    ref = crustal_ref_dir(lats, ref_lat, ref_lon, epoch)
    ref_lat = ref["ref_lat"]
    l, m, n = ref["l"], ref["m"], ref["n"]

    dlat = float(np.median(np.diff(lats)))
    dlon = float(np.median(np.diff(lons)))
    dx_km = dlat * KM_PER_DEG
    dy_km = dlon * KM_PER_DEG * math.cos(math.radians(ref_lat))

    g0 = np.nan_to_num(grid, nan=0.0)
    ny, nx = g0.shape
    py, px = ny // 2, nx // 2
    gp = np.pad(g0, ((py, py), (px, px)), mode="reflect")
    gp = gp * np.outer(tukey(gp.shape[0], alpha=2 * py / gp.shape[0]),
                       tukey(gp.shape[1], alpha=2 * px / gp.shape[1]))

    kx = 2 * np.pi * np.fft.fftfreq(gp.shape[0], d=dx_km)[:, None]
    ky = 2 * np.pi * np.fft.fftfreq(gp.shape[1], d=dy_km)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    Ft = np.fft.fft2(gp)
    den = 1j * (l * kx + m * ky) + n * k
    den[0, 0] = 1.0
    inv = np.where(k > 0, 1.0 / den, 0.0)

    cy, cx = slice(py, py + ny), slice(px, px + nx)
    bx = np.real(np.fft.ifft2(1j * kx * Ft * inv))[cy, cx]
    by = np.real(np.fft.ifft2(1j * ky * Ft * inv))[cy, cx]
    bz = np.real(np.fft.ifft2(k * Ft * inv))[cy, cx]
    return bx, by, bz


class CrustalVector:
    """복원된 지각 이상벡터의 쌍선형 보간자. 반환 (b_north, b_east, b_down)."""

    def __init__(self, lons, lats, bx, by, bz):
        self.lons, self.lats = lons, lats
        self.b = (bx, by, bz)
        self.dlon = float(np.median(np.diff(lons)))
        self.dlat = float(np.median(np.diff(lats)))

    def __call__(self, lat, lon):
        lat = np.atleast_1d(np.asarray(lat, float))
        lon = np.atleast_1d(np.asarray(lon, float))
        fx = (lon - self.lons[0]) / self.dlon
        fy = (lat - self.lats[0]) / self.dlat
        i0 = np.clip(np.floor(fx).astype(int), 0, self.lons.size - 2)
        j0 = np.clip(np.floor(fy).astype(int), 0, self.lats.size - 2)
        tx = np.clip(fx - i0, 0, 1)
        ty = np.clip(fy - j0, 0, 1)
        out = []
        for a in self.b:
            out.append((1 - tx) * (1 - ty) * a[j0, i0]
                       + tx * (1 - ty) * a[j0, i0 + 1]
                       + (1 - tx) * ty * a[j0 + 1, i0]
                       + tx * ty * a[j0 + 1, i0 + 1])
        return tuple(out)


def crustal_di(lat, lon, elev_m, cvec, epoch=None):
    """
    지각 이상벡터가 만드는 편각·복각 기여량(도)을 계산한다.

    IGRF 벡터에 b 를 더한 뒤 각을 다시 재는 방식이라 소각 근사를 쓰지 않는다.
    """
    if epoch is None:
        epoch = dt.datetime(int(EPOCH), 1, 1)
    lat = np.atleast_1d(np.asarray(lat, float))
    lon = np.atleast_1d(np.asarray(lon, float))
    elev_m = np.atleast_1d(np.asarray(elev_m, float))
    D_i, I_i, _, X, Y, Z = igrf_dif(lat, lon, elev_m, epoch)
    bx, by, bz = cvec(lat, lon)
    bx, by, bz = (np.nan_to_num(v, nan=0.0) for v in (bx, by, bz))
    Xc, Yc, Zc = X + bx, Y + by, Z + bz
    dD = np.degrees(np.arctan2(Yc, Xc)) - D_i
    dI = np.degrees(np.arctan2(Zc, np.hypot(Xc, Yc))) - I_i
    return dD, dI


def attach_crustal_di(sites, cvec):
    """측점표에 지각층의 D·I 기여량 열(crD, crI)을 붙인다. cvec 이 None 이면 0."""
    sites = sites.copy()
    if cvec is None:
        sites["crD"] = 0.0
        sites["crI"] = 0.0
    else:
        dD, dI = crustal_di(sites["lat"].values, sites["lon"].values,
                            sites["elev_m"].values, cvec)
        sites["crD"] = dD
        sites["crI"] = dI
    return sites


# --------------------------------------------------------- Layer 2 Regional
def poly_terms(lat, lon, degree=REGIONAL_DEGREE):
    """정규화 좌표 다항식 기저. 반환 shape (n, n_terms)."""
    u = (np.asarray(lat, float) - LAT0)
    v = (np.asarray(lon, float) - LON0)
    cols = []
    for total in range(degree + 1):
        for p in range(total + 1):
            cols.append((u ** (total - p)) * (v ** p))
    return np.column_stack(cols)


def poly_labels(degree=REGIONAL_DEGREE):
    out = []
    for total in range(degree + 1):
        for p in range(total + 1):
            out.append(f"u^{total - p} v^{p}")
    return out


def igrf_residuals(pts):
    """관측 성과에서 IGRF-14 를 뺀 잔차(방문 단위)."""
    D_i, I_i, F_i, *_ = igrf_dif(
        pts["lat"].values, pts["lon"].values, pts["elev_m"].values, pts["date"].tolist()
    )
    return pd.DataFrame(
        {
            "dD": pts["D_deg"].values - D_i,
            "dI": pts["I_deg"].values - I_i,
            "dF": pts["F_nT"].values - F_i,
        }
    )


def robust_lstsq(A, y, n_sigma=2.5, n_iter=10):
    """
    MAD 기반 반복 재가중 최소제곱.

    Lim2011 이 32점 중 27 inlier 만 사용한 것과 같은 취지로,
    조사 블런더·국지 자성체 영향을 받은 측점을 반복적으로 배제한다.

    표준편차 대신 MAD(중앙값 절대편차)로 산포를 추정하는 이유:
    블런더가 표준편차 자체를 부풀려 시그마 클리핑이 아무것도
    걸러내지 못하는 마스킹 효과를 피하기 위함이다.
    반환: (계수, inlier 불리언 마스크)
    """
    y = np.asarray(y, float)
    finite = np.isfinite(y)
    mask = finite.copy()

    for _ in range(n_iter):
        c, *_ = np.linalg.lstsq(A[mask], y[mask], rcond=None)
        resid = y - A @ c
        # 1.4826 = 정규분포에서 MAD -> 표준편차 환산계수
        scale = 1.4826 * np.median(np.abs(resid[mask] - np.median(resid[mask])))
        if scale <= 0:
            break
        new_mask = finite & (np.abs(resid - np.median(resid[mask])) <= n_sigma * scale)
        # 자유도 확보를 위해 최소 측점 수는 유지
        if new_mask.sum() < A.shape[1] + 3:
            break
        if (new_mask == mask).all():
            break
        mask = new_mask

    c, *_ = np.linalg.lstsq(A[mask], y[mask], rcond=None)
    return c, mask


def fit_regional(sites, crustal, degree=REGIONAL_DEGREE):
    """
    측점별 평균 잔차에 다항식을 적합(시그마 클리핑 적용).

    F 는 Crustal 층이 단파장을 먼저 설명하므로 그 잔차에 적합한다.
    D·I 도 CRUSTAL_VECTOR 를 쓰면(sites 에 crD·crI 열이 있으면) 같은 방식으로
    지각 기여를 먼저 빼고 적합한다. 열이 없으면 종전처럼 전체 잔차에 적합한다.
    """
    crust = np.nan_to_num(
        crustal(sites["lat"].values, sites["lon"].values), nan=0.0
    )
    A = poly_terms(sites["lat"].values, sites["lon"].values, degree)

    # 지각 벡터 기여를 먼저 뺀다 (열이 없으면 0 — 종전 동작과 동일)
    crD = sites["crD"].values if "crD" in sites else 0.0
    crI = sites["crI"].values if "crI" in sites else 0.0

    # 재방문 산포로 불합격한 측점은 해당 성분에서 NaN 처리해 적합에서 제외
    y_D = np.where(sites["D_ok"].values, sites["dD"].values - crD, np.nan)
    y_I = np.where(sites["I_ok"].values, sites["dI"].values - crI, np.nan)

    # inlier 판정은 IGRF 잔차만으로 수행한다.
    # Crustal 보정 후 잔차로 판정하면 "Crustal 이 잘 맞는 측점"만 남게 되어
    # Crustal 층의 기여도가 부풀려지는 선택편향이 생긴다.
    inliers = {}
    _, inliers["D"] = robust_lstsq(A, y_D)
    _, inliers["I"] = robust_lstsq(A, y_I)
    _, inliers["F"] = robust_lstsq(A, sites["dF"].values)

    # 계수는 각 성분의 inlier 집합 위에서 산출
    coef = {}
    coef["D"] = np.linalg.lstsq(A[inliers["D"]], y_D[inliers["D"]], rcond=None)[0]
    coef["I"] = np.linalg.lstsq(A[inliers["I"]], y_I[inliers["I"]], rcond=None)[0]
    mF = inliers["F"]
    coef["F"] = np.linalg.lstsq(
        A[mF], (sites["dF"].values - crust)[mF], rcond=None
    )[0]

    return coef, crust, inliers


# ------------------------------------------------------------------ 검증
def validate(sites, coef, crust, inliers, degree=REGIONAL_DEGREE):
    """IGRF 단독 / +Crustal / +Regional 단계별 RMS 를 inlier 기준으로 산출."""
    A = poly_terms(sites["lat"].values, sites["lon"].values, degree)
    rows = []

    crD = sites["crD"].values if "crD" in sites else np.zeros(len(sites))
    crI = sites["crI"].values if "crI" in sites else np.zeros(len(sites))
    has_cv = bool(np.any(crD != 0) or np.any(crI != 0))

    for comp, resid, cr, reg in [
        ("D_deg", sites["dD"].values, crD, A @ coef["D"]),
        ("I_deg", sites["dI"].values, crI, A @ coef["I"]),
    ]:
        m = inliers[comp[0]]
        rows.append((comp, "IGRF", rms(resid[m]), int(m.sum())))
        if has_cv:
            rows.append((comp, "+Crustal", rms((resid - cr)[m]), int(m.sum())))
            rows.append((comp, "+Crustal+Regional",
                         rms((resid - cr - reg)[m]), int(m.sum())))
        else:
            rows.append((comp, "+Regional", rms((resid - reg)[m]), int(m.sum())))

    m = inliers["F"]
    dF = sites["dF"].values
    rows.append(("F_nT", "IGRF", rms(dF[m]), int(m.sum())))
    rows.append(("F_nT", "+Crustal", rms((dF - crust)[m]), int(m.sum())))
    rows.append(
        ("F_nT", "+Crustal+Regional", rms((dF - crust - A @ coef["F"])[m]), int(m.sum()))
    )

    return pd.DataFrame(rows, columns=["성분", "단계", "RMS", "inlier수"])


def rms(x):
    x = np.asarray(x, float)
    return float(np.sqrt(np.nanmean(x ** 2)))


def loo_cv(sites, crustal, degree=REGIONAL_DEGREE):
    """
    Leave-one-out 교차검증.

    적합에 쓰이지 않은 측점에서의 예측오차이므로 과적합에 속지 않는다.
    inlier 로 판정된 측점만 평가 대상으로 삼는다.
    """
    _, _, inl = fit_regional(sites, crustal, degree)
    errs = {"D": [], "I": [], "F": []}

    for i in range(len(sites)):
        train = sites.drop(sites.index[i])
        test = sites.iloc[[i]]
        c, _, _ = fit_regional(train, crustal, degree)

        A = poly_terms(test["lat"].values, test["lon"].values, degree)
        cr = np.nan_to_num(crustal(test["lat"].values, test["lon"].values), nan=0.0)

        crD = float(test["crD"].values[0]) if "crD" in test else 0.0
        crI = float(test["crI"].values[0]) if "crI" in test else 0.0
        if inl["D"][i]:
            errs["D"].append(test["dD"].values[0] - crD - (A @ c["D"])[0])
        if inl["I"][i]:
            errs["I"].append(test["dI"].values[0] - crI - (A @ c["I"])[0])
        if inl["F"][i]:
            errs["F"].append(test["dF"].values[0] - cr[0] - (A @ c["F"])[0])

    return {k: rms(v) for k, v in errs.items()}


def crustal_diagnostics(sites, crust, inliers):
    """
    Crustal 층이 지표 점 잔차를 실제로 얼마나 설명하는지 정량화.

    설명자료의 사례연구(mag_1.txt, 측선간격 250 m)에서는 점 잔차와
    항공자력 평균이 5 nT 이내로 일치했다. 우리가 보유한 자료는
    전국 1.5분(약 2.8 km) 격자 컴필레이션이므로 그보다 훨씬 성긴다.
    """
    m = inliers["F"] & np.isfinite(crust)
    dF = sites["dF"].values

    r = float(np.corrcoef(dF[m], crust[m])[0, 1]) if m.sum() > 2 else np.nan

    # 잔차 자체를 0 으로 줄이는 것이 목표이므로 평균 주위 분산이 아니라
    # 0 주위 제곱평균(RMS)의 감소율로 설명력을 정의한다.
    rms_before = rms(dF[m])
    rms_after = rms(dF[m] - crust[m])
    reduction = 100 * (1 - rms_after / rms_before) if rms_before > 0 else np.nan

    return {
        "n": int(m.sum()),
        "corr": r,
        "rms_reduction_pct": reduction,
        "rms_before_nT": rms_before,
        "rms_after_nT": rms_after,
    }


def select_degree(sites, crustal, candidates=(1, 2, 3)):
    """LOO-CV 오차가 가장 작은 다항식 차수를 고른다."""
    scores = {}
    for d in candidates:
        n_terms = (d + 1) * (d + 2) // 2
        if len(sites) < n_terms + 3:
            continue
        cv = loo_cv(sites, crustal, d)
        # D 와 F 를 각 KPI 로 정규화해 합산
        scores[d] = cv["D"] / 0.1 + cv["F"] / 50.0
        print(f"   degree {d} ({n_terms}계수): "
              f"LOO D={cv['D']:.4f}deg  I={cv['I']:.4f}deg  F={cv['F']:.1f}nT")
    best = min(scores, key=scores.get)
    return best


# ------------------------------------------------------------------ 내보내기
def export_model(coef, lons, lats, grid, sites, val, cv, degree, rep, cd,
                 epoch_label="", years=(), vec=None):
    DOCS_DATA.mkdir(parents=True, exist_ok=True)

    # 격자를 정수 nT 로 양자화하고 NaN 은 null 로 (JSON 용량 절감)
    g = np.where(np.isnan(grid), None, np.round(grid).astype(object))
    flat = [None if v is None else int(v) for v in g.ravel()]

    # 지각 이상벡터 — 북·동 두 성분만 싣고 하향은 클라이언트에서 유도한다.
    #   b_down = (dF + dc - l*b_north - m*b_east) / n
    # 스칼라 격자와 벡터가 정의상 정합하고 용량도 한 성분 아낀다.
    vec_out = None
    if vec is not None:
        bx, by, bz = vec
        ref = crustal_ref_dir(lats)
        proj = ref["l"] * bx + ref["m"] * by + ref["n"] * bz
        dc = float(np.median(proj - np.nan_to_num(grid, nan=0.0)))
        vec_out = {
            "unit": "nT",
            "note": ("anomaly vector recovered from scalar dF by Fourier "
                     "potential-field inversion; b_down is derived on the client "
                     "as (dF + dc - l*b_north - m*b_east)/n"),
            "l": ref["l"], "m": ref["m"], "n": ref["n"],
            "ref_D": ref["D"], "ref_I": ref["I"], "dc_nT": dc,
            "b_north": [int(v) for v in np.round(bx).ravel()],
            "b_east": [int(v) for v in np.round(by).ravel()],
        }

    model = {
        "name": "Korea LMM v1 (quiet-time baseline)",
        "epoch": EPOCH,
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "layers": {
            "core": "IGRF-14 (degree 13, ppigrf/IGRF14.shc)",
            "regional": (f"polynomial degree {degree} on {len(sites)} absolute sites "
                         f"({epoch_label})"),
            "crustal": ("KIGAM aeromagnetic anomaly grid 1.5 arcmin (1982-2018)"
                        + (" + anomaly vector (b_north, b_east, b_down) recovered by "
                           "Fourier potential-field inversion, applied to D and I"
                           if vec_out is not None else
                           " (scalar dF only, applied to F)")),
            "external": (
                f"APPLIED TO DECLINATION ONLY (mode={EXTERNAL_MODE}, "
                f"source={EXTERNAL_SOURCE}). Session-time disturbance is estimated "
                "from four observatories (Cheongyang CYG, Jeju, Gangneung, Icheon) "
                "by first-order spatial interpolation of the quiet-night baseline "
                "deviation, and subtracted from the survey declination. Storm "
                "sessions (Kp>2) are left uncorrected because the spatial "
                "approximation breaks down first under disturbance. Total field F "
                "is NOT corrected: its residual is dominated by crustal-field "
                "mismatch (~190 nT) rather than external field (~38 nT), so a "
                "spatially smooth time correction only adds noise. The model "
                "itself remains a quiet-time baseline - no external term is "
                "evaluated at prediction time."),
        },
        "observation_years": sorted(int(y) for y in years),
        "epoch_label": epoch_label,
        "sites": sites[["name", "lat", "lon", "n_visit"]].to_dict(orient="records"),
        "repeatability": rep.to_dict(orient="records"),
        "regional": {
            "lat0": LAT0,
            "lon0": LON0,
            "degree": degree,
            "terms": poly_labels(degree),
            "D": list(map(float, coef["D"])),
            "I": list(map(float, coef["I"])),
            "F": list(map(float, coef["F"])),
        },
        "crustal": {
            "lon0": float(lons[0]),
            "lat0": float(lats[0]),
            "dlon": float(np.median(np.diff(lons))),
            "dlat": float(np.median(np.diff(lats))),
            "nlon": int(lons.size),
            "nlat": int(lats.size),
            "unit": "nT",
            "values": flat,
            "vector": vec_out,
        },
        "validation": val.to_dict(orient="records"),
        "loo_cv": cv,
        "crustal_diagnostics": cd,
    }

    out = DOCS_DATA / "lmm_model.json"
    out.write_text(json.dumps(model, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"[저장] {out}  ({out.stat().st_size/1024:.0f} KB)")
    return out


def main():
    print("=" * 60)
    print("한반도 LMM v1 구축")
    print("=" * 60)

    pts = load_all_points(include_2019=INCLUDE_2019)
    src = "성과표 2022~2025 + 2019 야장" if INCLUDE_2019 else "성과표 2022~2025"
    print(f"[Layer 2] {src}")
    print(f"          {len(pts)}회 관측 / {pts['name'].nunique()}개 측점 "
          f"({int(pts['year'].min())}~{int(pts['year'].max())})")

    lons, lats, grid = load_kigam_grid()
    crustal = CrustalGrid(lons, lats, grid)
    print(f"[Layer 3] KIGAM 자력이상 격자 {lons.size}x{lats.size} "
          f"(유효 {np.isfinite(grid).sum()}점, 간격 "
          f"{np.median(np.diff(lons))*111:.1f} km, "
          f"커버리지 {100*np.isfinite(grid).sum()/grid.size:.0f}%)")
    vec = None
    if CRUSTAL_VECTOR:
        vec = crustal_vector(lons, lats, grid)
        _ref = crustal_ref_dir(lats)
        print(f"[Layer 3] 지각 이상벡터 복원 (주자기장 D={_ref['D']:.2f}deg "
              f"I={_ref['I']:.2f}deg) -> 편각·복각에도 기여")
    print("[Layer 4] CYG 1분 자료 없음 -> 외부장 보정 미적용 (정온시 baseline 모델)")

    res = igrf_residuals(pts)

    rep = repeatability(pts, res)
    print("\n--- 동일 측점 재방문 산포 (외부장 잡음의 실측 하한) ---")
    print(rep.round(1).to_string(index=False))
    print(f"\n  F 산포 중앙값 {rep['F_산포_nT'].median():.0f} nT "
          f"-> KPI 50 nT 대비 {'초과' if rep['F_산포_nT'].median() > 50 else '이내'}")

    sites = aggregate_sites(pts, res)
    print(f"\n[집계] 측점 {len(sites)}개 (재방문 평균으로 외부장 잡음 부분 완화)")

    cvec = CrustalVector(lons, lats, *vec) if vec is not None else None
    sites = attach_crustal_di(sites, cvec)
    if cvec is not None:
        print(f"          지각 벡터 기여 D {rms(sites['crD'])*60:.2f}분 · "
              f"I {rms(sites['crI'])*60:.2f}분")

    print("\n--- Regional 다항식 차수 선정 (LOO-CV) ---")
    degree = select_degree(sites, crustal)
    print(f"   -> 채택: degree {degree}")

    coef, crust, inl = fit_regional(sites, crustal, degree)
    print(f"[Layer 2] inlier D={inl['D'].sum()}/{len(sites)}  "
          f"I={inl['I'].sum()}/{len(sites)}  F={inl['F'].sum()}/{len(sites)}")

    val = validate(sites, coef, crust, inl, degree)
    print("\n--- 단계별 RMS (inlier 기준) ---")
    print(val.round(4).to_string(index=False))

    cd = crustal_diagnostics(sites, crust, inl)
    print("\n--- Crustal 층 설명력 진단 ---")
    print(f"  상관계수 r = {cd['corr']:+.3f}  (n={cd['n']})")
    print(f"  RMS 감소율 = {cd['rms_reduction_pct']:.1f} %")
    print(f"  RMS {cd['rms_before_nT']:.1f} -> {cd['rms_after_nT']:.1f} nT")
    if cd["rms_reduction_pct"] < 40:
        print("  [경고] 전국 1.5분 격자는 지표 점 잔차를 거의 설명하지 못한다.")
        print("         원측선 자료는 존재하지 않으므로 현 격자가 해상도 상한이다.")

    cv = loo_cv(sites, crustal, degree)
    print("\n--- Leave-one-out 교차검증 (실제 예측성능) ---")
    print(f"D: {cv['D']:.4f} deg   I: {cv['I']:.4f} deg   F: {cv['F']:.1f} nT")

    print("\n--- KPI 판정 (D < 0.1 deg, F < 50 nT) ---")
    d_fin = float(val.query("성분=='D_deg'")["RMS"].values[-1])
    f_fin = float(val.query("성분=='F_nT'")["RMS"].values[-1])
    print(f"적합 D-RMS {d_fin:.4f} deg -> {'통과' if d_fin < 0.1 else '위반'}")
    print(f"적합 F-RMS {f_fin:.1f} nT   -> {'통과' if f_fin < 50 else '위반'}")
    print(f"LOO  D-RMS {cv['D']:.4f} deg -> {'통과' if cv['D'] < 0.1 else '위반'}")
    print(f"LOO  F-RMS {cv['F']:.1f} nT   -> {'통과' if cv['F'] < 50 else '위반'}")

    # 관측연도는 하드코딩하지 않는다 — 2019 야장 포함 여부에 따라 달라진다
    years = sorted({int(y) for y in pts["year"]})
    groups, run = [], [years[0], years[0]]
    for y in years[1:]:
        if y == run[1] + 1:
            run[1] = y
        else:
            groups.append(tuple(run))
            run = [y, y]
    groups.append(tuple(run))
    epoch_label = ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in groups)

    export_model(coef, lons, lats, grid, sites, val, cv, degree, rep, cd,
                 epoch_label=epoch_label, years=years, vec=vec)

    DOCS_OUT.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    val.to_csv(DOCS_OUT / f"{stamp}_LMM_validation.csv", index=False, encoding="utf-8-sig")
    rep.to_csv(DOCS_OUT / f"{stamp}_LMM_repeatability.csv", index=False, encoding="utf-8-sig")
    print(f"[저장] {DOCS_OUT / f'{stamp}_LMM_validation.csv'}")
    print(f"[저장] {DOCS_OUT / f'{stamp}_LMM_repeatability.csv'}")


if __name__ == "__main__":
    main()
