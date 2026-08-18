# -*- coding: utf-8 -*-
"""
성과표 원본 감사 — 성과표 값이 원시 세션평균인가, 기준시점 환산값인가
======================================================================

13절 작업순서의 1단계. 성과표(지자기측량 성과정리 22_25)의 D·I·F 가

  (a) 야장 세션의 **원시 평균** 인지,
  (b) 기준시점으로 **영년변화 환산**된 값인지,
  (c) 일변화(외부장) **보정이 적용된** 값인지

를 가린다. 이 판정에 따라 External 층 적용 방식이 달라진다 — 이미 환산된 값에
실제 관측일시로 영년변화를 다시 걸면 이중보정이 된다.

판정 방법
---------
① 야장 대조 : 야장에 세션별 편각·복각·자력과 「전자력(평균)」이 그대로 있다.
              성과표와 직접 비교하면 (a)/(b) 를 한 번에 가른다.
② 영년변화  : 같은 측점을 2년 간격으로 재측정한 쌍에서 성과표 차이와 IGRF 영년변화
              차이를 비교한다. 원시값이면 둘이 같고, 환산값이면 성과표 차가 0 에 가깝다.

⚠ 야장 시각은 **KST** 다. IGRF 는 시각 의존성이 거의 없어 이 판정에는 영향이 없지만,
   External 단계에서 CYG(UTC)·Kp(UT) 와 결합할 때는 반드시 -9h 변환할 것.

    python audit_survey_table.py
"""
import datetime as dt
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from lmm_build import load_survey_points
from lmm_fieldbook import ngii_files, SKIP_SHEETS, parse_file

ROOT = Path(__file__).parent
OUT = ROOT / "docs" / "output"

KST = dt.timezone(dt.timedelta(hours=9))


# ── 야장에서 세션별 성과(D·I·F)와 평균 추출 ─────────────────────────────────
def _norm(v):
    return re.sub(r"\s+", "", str(v)) if v is not None else ""


def dms_to_deg(s):
    """'-8°38′45.5″' → 십진도. 실패 시 None."""
    if s is None:
        return None
    t = str(s).strip()
    if not t or t in ("nan", "-"):
        return None
    neg = t.lstrip().startswith("-")
    nums = re.findall(r"\d+\.?\d*", t)
    if not nums:
        return None
    v = float(nums[0])
    if len(nums) >= 2:
        v += float(nums[1]) / 60
    if len(nums) >= 3:
        v += float(nums[2]) / 3600
    return -v if neg else v


def label_row(df, key, start=0):
    """A열에서 라벨 행을 찾는다(공백 무시). 없으면 None."""
    for r in range(start, df.shape[0]):
        if _norm(df.iat[r, 0]) == key:
            return r
    return None


def parse_results(path):
    """야장 1파일 → [{측점, 날짜, 회차, D, I, F}] + 전자력(평균).

    관측회차는 가로 3열 간격(c1, c4, c7, …)이다. 세로로 블록이 반복될 수 있어
    '편각' 라벨을 만날 때마다 블록으로 처리한다.
    """
    out, means, bad = [], [], []
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return out, means, bad
    for sh in xl.sheet_names:
        if sh in SKIP_SHEETS:
            continue
        try:
            df = xl.parse(sh, header=None)
        except Exception:
            continue
        if df.shape[0] < 20 or df.shape[1] < 8:
            continue
        site = str(df.iat[2, 1]).strip() if df.shape[0] > 2 else ""
        if not site or site == "nan":
            continue

        r = 0
        while True:
            rD = label_row(df, "편각", r)
            if rD is None:
                break
            rI = label_row(df, "복각", rD)
            rF = label_row(df, "자력", max(0, rD - 6))
            rDate = label_row(df, "일자", max(0, rD - 12))
            # 세션 열: c1 부터 3칸 간격
            for c in range(1, min(df.shape[1], 20), 3):
                D = dms_to_deg(df.iat[rD, c]) if rD is not None else None
                I = dms_to_deg(df.iat[rI, c]) if rI is not None else None
                F = df.iat[rF, c] if rF is not None else None
                F = float(F) if isinstance(F, (int, float, np.floating)) else None
                d = df.iat[rDate, c] if rDate is not None else None
                d = pd.Timestamp(d).normalize() if isinstance(
                    d, (dt.datetime, pd.Timestamp)) else None
                # ⚠ 야장에 각도 셀이 날짜로 저장된 사례가 있다
                #   (미원 2024 5회차 복각 = '1900-01-04'). 한반도 범위로 배제한다.
                if D is not None and not (-20 <= D <= 20):
                    bad.append((site, "편각", str(df.iat[rD, c])[:24]))
                    D = None
                if I is not None and not (30 <= I <= 70):
                    bad.append((site, "복각", str(df.iat[rI, c])[:24]))
                    I = None
                if F is not None and not (30000 <= F <= 60000):
                    F = None
                if D is None and I is None and F is None:
                    continue
                out.append(dict(측점=site, 날짜=d, D=D, I=I, F=F,
                                파일=Path(path).name, 시트=sh))
            # 전자력(평균)
            rM = label_row(df, "자력/분력", rD)
            if rM is not None:
                for c in range(1, min(df.shape[1], 20)):
                    if "전자력" in str(df.iat[rM, c]):
                        v = df.iat[rM, c + 1] if c + 1 < df.shape[1] else None
                        if isinstance(v, (int, float, np.floating)):
                            means.append(dict(측점=site, F_평균=float(v),
                                              파일=Path(path).name))
                        break
            r = rD + 1
    return out, means, bad


def file_year(path):
    """야장 1파일의 관측 연도. 검증된 lmm_fieldbook 파서의 날짜를 재사용한다.

    성과 블록의 '일자' 라벨은 파일마다 위치·표기가 흔들려 직접 찾으면 놓친다.
    파일 하나에 여러 해가 섞이는 경우는 없으므로 최빈 연도를 쓴다.
    """
    try:
        recs = parse_file(Path(path))
    except Exception:
        return None
    ys = [pd.Timestamp(r["날짜"]).year for r in recs
          if r.get("날짜") is not None and pd.notna(r.get("날짜"))]
    if not ys:
        return None
    return int(pd.Series(ys).mode().iloc[0])


def collect_results():
    rows, means, bads = [], [], []
    for f in ngii_files():
        a, b, c = parse_results(f)
        bads += [(Path(f).name, *x) for x in c]
        y = file_year(f) if a else None
        for r in a:
            r["연도"] = y
        rows += a
        means += b
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.dropna(subset=["측점", "연도"])
        df["연도"] = df["연도"].astype(int)
        # 빈 세션열(F=0, D/I 없음)은 양식의 미사용 칸이다
        df = df[~((df["F"].fillna(0) == 0) & df["D"].isna() & df["I"].isna())]
    return df, pd.DataFrame(means), bads


# ── ② 영년변화 검정 ─────────────────────────────────────────────────────────
def igrf_sv(lat, lon, elev_m, y0, y1):
    """두 시점 사이 IGRF 예측 변화량 (ΔD도, ΔI도, ΔF nT)."""
    import ppigrf

    def elems(year):
        d = dt.datetime(int(year), 7, 1)
        Be, Bn, Bu = ppigrf.igrf(lon, lat, elev_m / 1000.0, d)
        X, Y, Z = float(Bn), float(Be), -float(Bu)
        H = np.hypot(X, Y)
        return (np.degrees(np.arctan2(Y, X)), np.degrees(np.arctan2(Z, H)),
                np.sqrt(X * X + Y * Y + Z * Z))
    a, b = elems(y0), elems(y1)
    return b[0] - a[0], b[1] - a[1], b[2] - a[2]


def main():
    print("=" * 78)
    print("성과표 원본 감사 — 원시 세션평균인가, 기준시점 환산값인가")
    print("=" * 78)

    sp = load_survey_points()
    sp["year"] = sp["year"].astype(int)
    fb, fmean, bads = collect_results()
    print(f"\n야장 세션성과 {len(fb)}건 · 전자력평균 {len(fmean)}건")
    if bads:
        print(f"\n⚠ 야장 셀 서식 오류 {len(bads)}건 (각도가 날짜 등으로 저장됨)")
        for fn, site, kind, raw in bads:
            print(f"   {site:6} {kind}  '{raw}'   ({fn})")

    # ── ① 야장 세션평균 ↔ 성과표 ──
    g = (fb.dropna(subset=["연도"])
           .groupby(["측점", "연도"])
           .agg(D_야장=("D", "mean"), I_야장=("I", "mean"), F_야장=("F", "mean"),
                n=("D", "size")).reset_index())
    g["연도"] = g["연도"].astype(int)
    m = sp.merge(g, left_on=["name", "year"], right_on=["측점", "연도"], how="inner")
    print(f"\n■ ① 야장 세션평균 대조 — 매칭 {len(m)}행")
    print(f"{'측점':8} {'연도':5} {'ΔD(′)':>9} {'ΔI(′)':>9} {'ΔF(nT)':>9} {'n':>3}")
    print("-" * 50)
    for _, r in m.sort_values(["name", "year"]).iterrows():
        dD = (r["D_deg"] - r["D_야장"]) * 60 if pd.notna(r["D_야장"]) else np.nan
        dI = (r["I_deg"] - r["I_야장"]) * 60 if pd.notna(r["I_야장"]) else np.nan
        dF = r["F_nT"] - r["F_야장"] if pd.notna(r["F_야장"]) else np.nan
        print(f"{r['name'][:7]:8} {r['year']:5} {dD:9.2f} {dI:9.2f} {dF:9.2f} "
              f"{int(r['n']):3}")
    for lab, col, unit, scale in (("ΔD", "D", "′", 60), ("ΔI", "I", "′", 60),
                                  ("ΔF", "F", " nT", 1)):
        obs = m[f"{col}_deg" if col != "F" else "F_nT"]
        fbv = m[f"{col}_야장"]
        d = (obs - fbv) * scale
        d = d.dropna()
        if len(d):
            print(f"  {lab} 중앙 {d.median():+.2f}{unit} · RMS {np.sqrt((d**2).mean()):.2f}"
                  f"{unit} · 최대 {d.abs().max():.2f}{unit}")

    # ── ② 영년변화 검정 ──
    print("\n■ ② 반복측정 영년변화 검정 (성과표 차 vs IGRF 예측 차)")
    print(f"{'측점':8} {'구간':11} {'성과ΔD(′)':>10} {'IGRFΔD(′)':>10} "
          f"{'성과ΔF':>9} {'IGRFΔF':>9}  판정")
    print("-" * 74)
    verdicts = []
    for name, grp in sp.groupby("name"):
        if len(grp) < 2:
            continue
        grp = grp.sort_values("year")
        a, b = grp.iloc[0], grp.iloc[-1]
        if a["year"] == b["year"]:
            continue
        sD, sI, sF = igrf_sv(a["lat"], a["lon"], a["elev_m"], a["year"], b["year"])
        oD = (b["D_deg"] - a["D_deg"]) * 60
        oF = b["F_nT"] - a["F_nT"]
        # 원시값이면 관측차 ≈ IGRF차, 환산값이면 관측차 ≈ 0
        raw_err = abs(oD - sD * 60)
        red_err = abs(oD)
        v = "원시값" if raw_err < red_err else "환산값"
        verdicts.append(v)
        print(f"{name[:7]:8} {a['year']}~{b['year']:<6} {oD:10.2f} {sD*60:10.2f} "
              f"{oF:9.1f} {sF:9.1f}  {v}")
    if verdicts:
        raw = verdicts.count("원시값")
        print(f"\n  판정: 원시값 {raw} · 환산값 {len(verdicts)-raw} / {len(verdicts)}쌍")

    OUT.mkdir(parents=True, exist_ok=True)
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    p = OUT / f"{ts}_성과표_원본감사.csv"
    m.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")
    return m


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
