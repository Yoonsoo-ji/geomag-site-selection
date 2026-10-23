# -*- coding: utf-8 -*-
"""
기존 측정점 절대관측 통합 시계열 — 다섯 출처를 한 표로
========================================================

D·I·F 측정값을 담은 자료가 저장소에 **다섯 갈래**로 흩어져 있다. 측점마다
시기별 변동을 보려면 그때마다 다섯 파일을 다시 열어야 했고, 같은 관측이
두 출처에 겹쳐 실리는 일도 있었다. 이 스크립트가 **하나로 모은다.**

| 출처 | 기간 | 단위 | 담긴 것 |
|---|---|---|---|
| `register` 지자기 점조서 관측이력 | 1984~2025 | 방문 | D·I·F·장비 |
| `legacy1019` '10~'19 관측현황 | 2010~2019 | 방문 | D·I·F |
| `survey2225` 22~25 성과정리 | 2022~2025 | 방문 | D·I·F·좌표·표고 |
| `kriss` 2012~2016 절대관측 원자료 | 2012~2016 | **세션** | D·I·F·시각·참방위각 |
| `fieldbook` 2019~2025 야장 | 2019~2025 | **세션** | D·I(원시 역산)·F·시각 |

산출은 둘이다.

  `docs/data/obs_archive_sessions.csv`  세션·방문 한 줄씩 — 원자료에 가장 가깝다
  `docs/data/obs_archive_visits.csv`    (측점, 일자) 로 묶은 방문 대표값

두 표 모두 **IGRF-14 잔차**(rD·rI·rF)를 함께 싣는다. 영년변화가 빠지므로
시기를 가로질러 견줄 수 있는 유일한 양이다.

규약
----
⚠️ **값을 고치지 않는다.** 부호 교정처럼 되돌릴 수 있는 변환만 걸고, 원값은
   `D_원본` 열에 남긴다. 범위 밖·블런더는 **표시만** 한다(현장 탈락 금지 규약).

⚠️ **측점 이름은 현재 이름으로 통일한다.** 옛 이름(임계·경주·양산·곡성·인제)은
   `SITE_CANON` 으로 바꾸되 원문 이름을 `측점_원문` 에 남긴다. 이름을 합치지
   않으면 같은 계보가 두 측점으로 갈려 시계열이 끊긴다.

⚠️ **좌표는 그 관측의 좌표를 쓴다.** 현재 관측망 좌표로 전부 덮으면 이설 전
   관측의 잔차가 틀어진다(삼척 2012 는 30 km 떨어진 폐지점에서 쟀다).
   관측 좌표가 없는 출처는 관측망 좌표로 채우고 `좌표출처` 에 적는다.

    python build_obs_archive.py
"""
from __future__ import annotations

import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"
DOCS = ROOT / "docs" / "data"
OUT_S = DOCS / "obs_archive_sessions.csv"
OUT_V = DOCS / "obs_archive_visits.csv"

import existing_network as EN                        # noqa: E402

# 옛 이름 → 현재 이름. 같은 계보가 두 측점으로 갈리지 않게 한다.
SITE_CANON = {"임계": "삼척", "경주": "영천", "양산": "언양",
              "곡성": "순창", "인제": "설악"}

# 한반도 범위 — 넘으면 「범위밖」으로 표시만 한다(버리지 않는다)
RANGE = {"D": (-20.0, 20.0), "I": (30.0, 70.0), "F": (30000.0, 60000.0)}


def parse_day(txt):
    """
    출처마다 다른 일자 표기를 (날짜, 정밀도) 로 돌려준다.

    ⚠️ **`pd.to_datetime(…, errors="coerce")` 를 그대로 쓰면 안 된다.** 판본에
       따라 첫 비어 있지 않은 값에서 형식을 하나 추정하고 나머지를 전부 NaT 로
       버린다 — 실제로 587행 중 482행이 날아갔다. 점조서가 「1993.10.19」·
       「1995.0927」·「1996.00.00」·「2010연구성과」를 섞어 쓰기 때문이다.

    정밀도는 **IGRF 영년변화를 어느 해상도로 걸 수 있는지**를 뜻한다.
    「연」이면 7월 1일로 대체하므로 최대 반년 어긋나고, 편각으로는 약 ±2.5′ 다
    (잔차 20~30′ 에 견주면 작지만 적어 두는 편이 낫다).
    """
    t = str(txt or "").strip()
    m = re.search(r"(19|20)\d{2}", t)
    if not m:
        return None, ""
    year = int(m.group(0))
    rest = t[m.end():]
    nums = [int(x) for x in re.findall(r"\d+", rest)]
    # 「1995.0927」처럼 월일이 붙은 경우
    if len(nums) == 1 and 100 <= nums[0] <= 1231:
        nums = [nums[0] // 100, nums[0] % 100]
    mo = nums[0] if nums and 1 <= nums[0] <= 12 else 0
    dy = nums[1] if len(nums) > 1 and 1 <= nums[1] <= 31 else 0
    if mo and dy:
        try:
            return dt.date(year, mo, dy), "일"
        except ValueError:
            return dt.date(year, mo, 1), "월"
    if mo:
        return dt.date(year, mo, 1), "월"
    return dt.date(year, 7, 1), "연"


def canon(name):
    n = str(name).strip()
    return SITE_CANON.get(n, n)


def flag(v, comp):
    if v is None or not np.isfinite(v):
        return "결측"
    lo, hi = RANGE[comp]
    return "정상" if lo <= v <= hi else "범위밖"


# ── 출처 ① 지자기 점조서 관측이력 ────────────────────────────────────
def src_register():
    """측점마다 1984~2025 관측이력. 가장 긴 시계열이다."""
    try:
        from parse_point_register import (REGISTER_XLSX, load_register,
                                          current_points)
    except Exception as e:
        print(f"  ! register 건너뜀: {e}")
        return pd.DataFrame()
    if not REGISTER_XLSX.exists():
        print(f"  ! register 원본이 없다 (28 MB · git 제외): {REGISTER_XLSX}")
        return pd.DataFrame()
    try:
        # ⚠️ `current_points` 는 **시트 목록을 인자로 받는다** — 계보마다
        #    현재점 하나를 매설일 기준으로 고른다. 폐지점 시트까지 넣으면
        #    같은 계보가 두 측점으로 갈린다.
        pts = current_points(load_register())
    except Exception as e:
        print(f"  ! register 원본을 못 읽음: {e}")
        return pd.DataFrame()

    # ⚠️ 측점 이름 열은 `지점명` 이 아니라 **`도엽명`** 이다. `지점명` 을 읽으면
    #    전부 None 이 되어 180 행이 한 측점으로 뭉친다(실제로 그랬다).
    # ⚠️ 좌표는 `EN.COORD_OVERRIDE` 를 거친 값을 쓴다 — 점조서 제천·남양 시트가
    #    현장과 어긋나므로 원값을 그대로 실으면 잔차가 통째로 틀어진다.
    ov = EN.load_register()
    OVC = {str(r["지점명"]).strip(): (float(r["위도"]), float(r["경도"]))
           for _, r in ov.iterrows()}

    rows = []
    for p in pts:
        nm = canon(p.get("도엽명"))
        la, lo = p.get("위도"), p.get("경도")
        csrc = "점조서 현재점"
        if nm in OVC:
            la, lo = OVC[nm]
            csrc = "점조서 현재점(교정 반영)"
        bury = _bury_year(p.get("매설일"))
        for r in p.get("관측기록") or []:
            day = str(r.get("일자", ""))[:10]
            # 현재 표석을 묻기 «전» 관측은 폐지된 다른 표석에서 한 것이다.
            # 좌표를 현재점으로 붙이면 그 방문의 잔차가 틀어진다 — 표시해 둔다.
            yr = None
            try:
                yr = int(day[:4])
            except Exception:
                pass
            pre = (bury is not None and yr is not None and yr < bury)
            rows.append(dict(
                측점=nm, 측점_원문=str(p.get("도엽명", "")).strip(),
                일자=day, 출처="register", 단위="방문", 회차=np.nan,
                위도=la, 경도=lo,
                좌표출처=(f"{csrc} · ⚠ {bury} 매설 이전 관측" if pre else csrc),
                좌표신뢰=("이설 이전 — 좌표 불확실" if pre else "현재점"),
                D=r.get("편각"), I=r.get("복각"), F=r.get("전자력"),
                장비=r.get("장비", ""), 점의부호=p.get("점의부호", "")))
    return pd.DataFrame(rows)


def _bury_year(v):
    """매설일에서 연도만. 「2022.10.12」·「2022-10-12」 등 꼴이 섞여 있다."""
    import re
    m = re.search(r"(19|20)[0-9]{2}", str(v or ""))
    return int(m.group(0)) if m else None


# ── 출처 ② '10~'19 관측현황 ──────────────────────────────────────────
def src_legacy1019():
    f = DATA / "'10~'19년 지자기점 관측현황(최종).xls"
    if not f.exists():
        print("  ! legacy1019 파일 없음")
        return pd.DataFrame()
    try:
        from geomag_site_selection import load_existing_sites
        d = load_existing_sites()
    except Exception as e:
        print(f"  ! legacy1019 건너뜀: {e}")
        return pd.DataFrame()
    col = {c.lower(): c for c in d.columns}

    def pick(*names):
        for n in names:
            if n in d.columns:
                return n
            if n.lower() in col:
                return col[n.lower()]
        return None

    cn, cy = pick("도엽명", "name"), pick("관측연도", "year")
    cd, ci, cf = pick("편각", "D_deg"), pick("복각", "I_deg"), pick("총자력", "F_nT")
    cla, clo = pick("위도", "lat"), pick("경도", "lon")
    if not cn:
        return pd.DataFrame()
    rows = []
    for _, r in d.iterrows():
        y = r.get(cy) if cy else None
        rows.append(dict(
            측점=canon(r[cn]), 측점_원문=str(r[cn]).strip(),
            일자=(f"{int(y)}-07-01" if pd.notna(y) else ""),
            출처="legacy1019", 단위="방문", 회차=np.nan,
            위도=r.get(cla), 경도=r.get(clo), 좌표출처="관측현황/점조서 병합",
            D=r.get(cd), I=r.get(ci), F=r.get(cf), 장비="", 점의부호=""))
    return pd.DataFrame(rows)


# ── 출처 ③ 22~25 성과정리 ────────────────────────────────────────────
def src_survey2225():
    f = DATA / "지자기측량 성과정리(22_25).xlsx"
    if not f.exists():
        print("  ! survey2225 파일 없음")
        return pd.DataFrame()
    d = pd.read_excel(f, sheet_name="Sheet1")
    d = d.rename(columns={"실": "lat", "실.1": "lon", "실.2": "D_deg",
                          "실.3": "I_deg", "표고": "elev_m",
                          "총자력": "F_nT", "관측연도": "year",
                          "도엽명": "name"})
    d["name"] = d["name"].ffill()
    d = d.dropna(subset=["lat", "lon"])
    rows = []
    for _, r in d.iterrows():
        y = r.get("year")
        rows.append(dict(
            측점=canon(r["name"]), 측점_원문=str(r["name"]).strip(),
            일자=(f"{int(y)}-07-01" if pd.notna(y) else ""),
            출처="survey2225", 단위="방문", 회차=np.nan,
            위도=r["lat"], 경도=r["lon"], 좌표출처="성과표 실측",
            D=r.get("D_deg"), I=r.get("I_deg"), F=r.get("F_nT"),
            장비="", 점의부호=""))
    return pd.DataFrame(rows)


# ── 출처 ④ 2012~2016 절대관측 원자료 ─────────────────────────────────
def src_kriss():
    f = DOCS / "kriss_legacy_sessions.csv"
    if not f.exists():
        print("  ! kriss 세션표 없음 — parse_kriss_legacy.py 를 먼저")
        return pd.DataFrame()
    d = pd.read_csv(f, encoding="utf-8-sig")
    d = d[d["잔재"].isna() | (d["잔재"].astype(str).str.strip() == "")]
    rows = []
    for _, r in d.iterrows():
        rows.append(dict(
            측점=canon(r["측점"]), 측점_원문=str(r["측점"]).strip(),
            일자=str(r.get("날짜", ""))[:10], 출처="kriss", 단위="세션",
            회차=r.get("회차"),
            위도=r.get("위도"), 경도=r.get("경도"), 좌표출처="관측 야장 머리글",
            D=r.get("편각"), I=r.get("복각"), F=r.get("총자력"),
            장비="", 점의부호=str(r.get("코드", "")),
            참방위각=r.get("FM_True"), 편각시작=r.get("편각시작")))
    return pd.DataFrame(rows)


# ── 출처 ⑤ 2019~2025 야장 (원시 판독값 역산) ─────────────────────────
def src_fieldbook():
    try:
        import lmm_verify2019 as V
        raw = V.collect_raw()
    except Exception as e:
        print(f"  ! fieldbook 건너뜀: {e}")
        return pd.DataFrame()
    if raw is None or not len(raw):
        return pd.DataFrame()
    # 역산 열 이름이 판본마다 달라 후보를 훑는다
    # ⚠️ `collect_raw()` 는 야장 기재값을 `D_야장`·`I_야장`·`F_야장` 으로
    #    낸다(ED/WD/EU/WU 원시 판독값도 함께 싣는다). 열 이름을 짐작하지 말고
    #    실제 이름부터 찾을 것 — 처음에 D_deg 로 찾아 197행이 전부 결측이 됐다.
    dcol = next((c for c in ("D_야장", "D_deg", "편각", "D")
                 if c in raw.columns), None)
    icol = next((c for c in ("I_야장", "I_deg", "복각", "I")
                 if c in raw.columns), None)
    fcol = next((c for c in ("F_야장", "F_nT", "총자력", "F", "전자력")
                 if c in raw.columns), None)
    if not dcol:
        print(f"  ! fieldbook 에서 편각 열을 못 찾음: {list(raw.columns)[:12]}")
    rows = []
    for _, r in raw.iterrows():
        rows.append(dict(
            측점=canon(r["측점"]), 측점_원문=str(r["측점"]).strip(),
            일자=str(r.get("날짜", ""))[:10], 출처="fieldbook", 단위="세션",
            회차=r.get("회차", np.nan),
            위도=r.get("위도"), 경도=r.get("경도"), 좌표출처="야장 머리글",
            D=(r.get(dcol) if dcol else None),
            I=(r.get(icol) if icol else None),
            F=(r.get(fcol) if fcol else None),
            장비="", 점의부호="",
            참방위각=r.get("FM_True", np.nan)))
    return pd.DataFrame(rows)


# ── IGRF 잔차 ────────────────────────────────────────────────────────
def attach_igrf(d):
    """
    관측 일시·좌표에서 IGRF-14 를 계산해 잔차를 붙인다.

    ⚠️ 영년변화가 빠지므로 **시기를 가로질러 견줄 수 있는 유일한 양**이다.
       같은 측점을 다시 재면 잔차가 같아야 한다 — 달라지면 그 방문에 문제가
       있다는 뜻이다(13절 재방문 감사와 같은 논리).
    """
    import ppigrf
    cache, out = {}, {"rD": [], "rI": [], "rF": []}
    for _, r in d.iterrows():
        la, lo, day = r["위도"], r["경도"], r["_date"]
        if (not np.isfinite(la) or not np.isfinite(lo) or pd.isna(day)):
            out["rD"].append(np.nan)
            out["rI"].append(np.nan)
            out["rF"].append(np.nan)
            continue
        key = (round(la, 4), round(lo, 4), day.year, day.month)
        if key not in cache:
            e, n, u = ppigrf.igrf(lo, la, 0.0,
                                  dt.datetime(day.year, day.month, 15))
            e, n, u = (float(np.ravel(x)[0]) for x in (e, n, u))
            h = np.hypot(e, n)
            cache[key] = (np.degrees(np.arctan2(e, n)),
                          np.degrees(np.arctan2(-u, h)),
                          float(np.sqrt(e * e + n * n + u * u)))
        d0, i0, f0 = cache[key]
        out["rD"].append((r["D"] - d0) * 60 if np.isfinite(r["D"]) else np.nan)
        out["rI"].append((r["I"] - i0) * 60 if np.isfinite(r["I"]) else np.nan)
        out["rF"].append((r["F"] - f0) if np.isfinite(r["F"]) else np.nan)
    for k, v in out.items():
        d[k + ("분" if k != "rF" else "_nT")] = np.round(v, 2)
    return d


def main():
    print("■ 출처별 수집")
    frames = []
    for name, fn in (("register", src_register),
                     ("legacy1019", src_legacy1019),
                     ("survey2225", src_survey2225),
                     ("kriss", src_kriss),
                     ("fieldbook", src_fieldbook)):
        f = fn()
        print(f"  {name:<12}{len(f):>5} 행"
              + (f" · 측점 {f['측점'].nunique()}" if len(f) else ""))
        if len(f):
            frames.append(f)
    if not frames:
        print("수집된 자료가 없다")
        return 1

    d = pd.concat(frames, ignore_index=True)
    for c in ("D", "I", "F", "위도", "경도", "회차", "참방위각"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        else:
            d[c] = np.nan

    # ── 편각 부호 교정 — 한반도는 전부 서편(음)이다 ──────────────────
    # 원값은 남긴다. 설악 2012 는 파일에 +7.9989 로 적혀 있어 그대로 쓰면
    # IGRF 잔차가 977′ 로 나온다(음수화하면 20.5′).
    d["D_원본"] = d["D"]
    d["D부호교정"] = d["D"] > 0
    d.loc[d["D"].notna(), "D"] = -d.loc[d["D"].notna(), "D"].abs()

    for comp in ("D", "I", "F"):
        d[f"{comp}_판정"] = [flag(v, comp) for v in d[comp]]

    # ── 좌표가 없는 행은 관측망 좌표로 채운다 ───────────────────────
    reg = EN.load_register()
    R = {str(r["지점명"]).strip(): (float(r["위도"]), float(r["경도"]))
         for _, r in reg.iterrows()}
    miss = d["위도"].isna() | d["경도"].isna()
    for i in d.index[miss]:
        c = R.get(d.at[i, "측점"])
        if c:
            d.at[i, "위도"], d.at[i, "경도"] = c
            d.at[i, "좌표출처"] = "관측망 좌표로 보충"
    print(f"■ 좌표 보충 {int(miss.sum())} 행 "
          f"(남은 결측 {int((d['위도'].isna()).sum())})")

    if "좌표신뢰" not in d.columns:
        d["좌표신뢰"] = np.nan
    d["좌표신뢰"] = d["좌표신뢰"].fillna(
        d["좌표출처"].map(lambda t: "관측 당시 좌표"
                          if "야장" in str(t) or "실측" in str(t)
                          else "현재점 또는 병합값"))

    parsed = [parse_day(v) for v in d["일자"]]
    d["_date"] = pd.to_datetime([x[0] for x in parsed], errors="coerce")
    d["일자정밀도"] = [x[1] for x in parsed]
    d["연도"] = d["_date"].dt.year
    d["일자"] = [(x[0].isoformat() if x[0] else "") for x in parsed]
    d["일자_원문"] = [str(v).strip() for v in d["일자_원문"]]         if "일자_원문" in d.columns else [str(v).strip() for v in d["일자"]]
    nd = int(d["_date"].isna().sum())
    print(f"■ 일자 해석 — 일 {sum(1 for x in parsed if x[1] == '일')}"
          f" · 월 {sum(1 for x in parsed if x[1] == '월')}"
          f" · 연 {sum(1 for x in parsed if x[1] == '연')}"
          f" · 못 읽음 {nd}")
    d = d.sort_values(["측점", "_date", "출처", "회차"], na_position="last")
    print("■ IGRF-14 잔차 계산 중…")
    d = attach_igrf(d)

    cols = ["측점", "측점_원문", "연도", "일자", "출처", "단위", "회차",
            "위도", "경도", "좌표출처", "D", "I", "F",
            "D_판정", "I_판정", "F_판정", "D_원본", "D부호교정",
            "rD분", "rI분", "rF_nT", "일자정밀도", "좌표신뢰",
            "참방위각", "장비", "점의부호"]
    S = d[[c for c in cols if c in d.columns]].reset_index(drop=True)

    # ── 방문 대표값 — 세션은 중앙값으로 묶는다 ──────────────────────
    vis = []
    for (s, day, src), g in S.groupby(["측점", "일자", "출처"], dropna=False):
        rec = dict(측점=s, 일자=day, 출처=src,
                   연도=(int(g["연도"].iloc[0])
                         if pd.notna(g["연도"].iloc[0]) else None),
                   세션=len(g), 위도=g["위도"].median(), 경도=g["경도"].median(),
                   좌표출처=g["좌표출처"].iloc[0])
        for comp, col in (("D", "D"), ("I", "I"), ("F", "F")):
            ok = g[g[f"{comp}_판정"] == "정상"]
            rec[col] = round(float(ok[col].median()), 6) if len(ok) else None
            rec[f"{comp}_정상세션"] = len(ok)
            rec[f"{comp}_산포"] = (
                round(float((ok[col].max() - ok[col].min())
                            * (60 if comp != "F" else 1)), 2)
                if len(ok) > 1 else None)
        for k, col in (("rD분", "rD분"), ("rI분", "rI분"), ("rF_nT", "rF_nT")):
            v = g[col].dropna()
            rec[k] = round(float(v.median()), 2) if len(v) else None
        rec["참방위각"] = (round(float(g["참방위각"].dropna().median()), 4)
                        if g["참방위각"].notna().any() else None)
        vis.append(rec)
    V = pd.DataFrame(vis).sort_values(["측점", "일자", "출처"])

    DOCS.mkdir(parents=True, exist_ok=True)
    S.to_csv(OUT_S, index=False, encoding="utf-8-sig")
    V.to_csv(OUT_V, index=False, encoding="utf-8-sig")

    print(f"\n■ 통합 결과")
    print(f"  세션·방문 행 {len(S)} · 방문 {len(V)} · 측점 {S['측점'].nunique()}")
    print(f"  연도 범위 {int(S['연도'].min())} ~ {int(S['연도'].max())}")
    print("\n  출처 × 단위")
    print(S.groupby(["출처", "단위"]).agg(
        행=("측점", "size"), 측점=("측점", "nunique"),
        D정상=("D_판정", lambda x: (x == "정상").sum()),
        F정상=("F_판정", lambda x: (x == "정상").sum())).to_string())
    print(f"\n  편각 부호 교정 {int(S['D부호교정'].sum())} 행")
    print(f"\n[저장] {OUT_S}")
    print(f"[저장] {OUT_V}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
