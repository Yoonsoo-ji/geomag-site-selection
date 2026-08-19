# -*- coding: utf-8 -*-
"""
2019년 야장 신뢰도 검증
=======================

2019 성과를 Regional 층에 투입하기 전에, 야장 자체로 검증 가능한 항목과
외부 자료와 대조 가능한 항목을 모두 점검한다.

각도 단위
---------
MAG-01H 야장의 세오돌라이트 판독값은 **곤(gon, 1회전 = 400)** 단위다.
따라서 대척은 180°가 아니라 200 gon 이고, 야장 수식의 180/360 은 실제로
200/400 이다. 도(°) 환산은 x0.9.

    편각 D[gon] = 참방위각[gon] - 마크평균 + (DEC1 - 200) - 100
    복각 I[gon] = [(200-SU) + (400-ND) + (SD-200) + NU] / 4

검증 항목
---------
[내부] A. 마크 폐합      Mark1-1 - Mark1-0 이 200 gon 인가
       B. 편각 재계산     DEC1 = (ED+WD+EU+WU)/4 및 최종 편각
       C. 복각 재계산     위 수식으로 최종 복각
       D. 분력 정합       F^2 = H^2 + Z^2
       E. 자오선 대칭성   ED/WU, WD/EU 쌍의 대척 관계
[외부] F. IGRF-14 잔차
       G. 세션 반복성     CYG 외부장 보정 후 세션 간 수렴 여부
       H. 성과표 대조     겹치는 측점의 2022~2025 성과와 영년변화 정합

실행:
    python lmm_verify2019.py
"""

import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lmm_build import igrf_dif, load_survey_points
from lmm_fieldbook import FIELDBOOK_DIR, PATTERNS, SKIP_SHEETS, _is_date, _is_time, _norm

GON2DEG = 0.9
FULL_GON = 400.0
HALF_GON = 200.0

# 판정 임계값
TOL_MARK_GON = 0.02      # 마크 폐합 허용 (0.02 gon = 65")
TOL_RECALC_DEG = 0.002   # 재계산 일치 허용 (7.2")
TOL_F_NT = 2.0           # 분력 정합 허용


def dms_to_deg(s):
    """-8°53′21.78″ 형태를 도(십진)로."""
    if not isinstance(s, str):
        return np.nan
    m = re.findall(r"-?\d+(?:\.\d+)?", s)
    if len(m) < 3:
        return np.nan
    d, mi, se = float(m[0]), float(m[1]), float(m[2])
    sign = -1 if s.strip().startswith("-") else 1
    return sign * (abs(d) + mi / 60 + se / 3600)


def _num(v):
    return float(v) if isinstance(v, (int, float)) and pd.notna(v) else np.nan


def _val_at(df, row, col):
    """라벨행에 값이 있으면 그것, 시각이면 다음 행의 값."""
    if row is None or row + 1 >= df.shape[0]:
        return np.nan
    v = df.iat[row, col]
    if isinstance(v, (int, float)) and pd.notna(v):
        return float(v)
    return _num(df.iat[row + 1, col])


def parse_sessions(df: pd.DataFrame, sheet: str, fname: str):
    nrow, ncol = df.shape
    labels = {r: _norm(df.iat[r, 0]) for r in range(nrow)}

    site = None
    for r in range(min(6, nrow)):
        row = [_norm(df.iat[r, c]) for c in range(min(ncol, 12))]
        for key in ("지자기측점명", "도엽(측점)"):
            if key in row and r + 1 < nrow:
                site = df.iat[r + 1, row.index(key)]
        if site is not None:
            break

    block_rows = [r for r in range(nrow)
                  if labels[r] in ("일자", "측정일자", "관측일자")]
    if not block_rows:
        return []

    def find(lo, hi, *tokens):
        for r in range(lo, hi):
            for t in tokens:
                if labels[r].startswith(t):
                    return r
        return None

    out = []
    for bi, br in enumerate(block_rows):
        end = block_rows[bi + 1] if bi + 1 < len(block_rows) else nrow

        col_date = {c: _is_date(df.iat[br, c]) for c in range(1, ncol)}
        col_date = {c: d for c, d in col_date.items() if d is not None}
        if not col_date:
            continue

        r_gps = find(br, end, "GPS")
        r_m0 = find(br, end, "Mark1-0")
        r_m1 = find(br, end, "Mark1-1")
        r_mn = find(br, end, "1/n")
        r_ed = find(br, end, "ED")
        r_wd = find(br, end, "WD")
        r_eu = find(br, end, "EU")
        r_wu = find(br, end, "WU")
        r_su = find(br, end, "SU")
        r_nd = find(br, end, "ND")
        r_sd = find(br, end, "SD")
        r_nu = find(br, end, "NU")
        r_dec1 = find(br, end, "DEC1")
        r_F = find(br, end, "자력")
        r_dec = find(br, end, "편각")
        r_inc = find(br, end, "복각")
        r_H = find(br, end, "수평분력")
        r_Z = find(br, end, "수직분력")

        for c, d in col_date.items():
            rec = {
                "파일": fname, "시트": sheet, "측점": site, "날짜": d.date(),
                "Mark1_0": _num(df.iat[r_m0, c]) if r_m0 else np.nan,
                "Mark1_1": _num(df.iat[r_m1, c]) if r_m1 else np.nan,
                "mark_mean": _num(df.iat[r_mn, c]) if r_mn else np.nan,
                "DEC1": _num(df.iat[r_dec1, c]) if r_dec1 else np.nan,
                "F_야장": _num(df.iat[r_F, c]) if r_F else np.nan,
                "H_야장": _num(df.iat[r_H, c]) if r_H else np.nan,
                "Z_야장": _num(df.iat[r_Z, c]) if r_Z else np.nan,
                "D_야장": dms_to_deg(df.iat[r_dec, c]) if r_dec else np.nan,
                "I_야장": dms_to_deg(df.iat[r_inc, c]) if r_inc else np.nan,
            }
            # 참방위각 (도/분/초가 c, c+1, c+2 에 배치)
            if r_gps is not None and c + 2 < ncol:
                dd, mm, ss = (_num(df.iat[r_gps, c]), _num(df.iat[r_gps, c + 1]),
                              _num(df.iat[r_gps, c + 2]))
                rec["방위각_deg"] = (abs(dd) + mm / 60 + ss / 3600) if np.isfinite(dd) else np.nan
            for key, rr in [("ED", r_ed), ("WD", r_wd), ("EU", r_eu), ("WU", r_wu),
                            ("SU", r_su), ("ND", r_nd), ("SD", r_sd), ("NU", r_nu)]:
                rec[key] = _val_at(df, rr, c)

            # 시각(있으면)
            ts = [df.iat[r, c] for r in range(br, end)
                  if c < ncol and _is_time(df.iat[r, c])]
            rec["시작"] = min(ts) if ts else None
            out.append(rec)
    return out


def collect_raw() -> pd.DataFrame:
    import glob
    rows = []
    # docs/data 의 2019 야장 + 지리원 야장 트리(2020~2025) 전체
    from lmm_fieldbook import ngii_files
    files = sorted({f for p in PATTERNS for f in glob.glob(str(FIELDBOOK_DIR / p))}
                   | set(ngii_files()))
    for f in files:
        try:
            xl = pd.ExcelFile(f)
        except Exception:
            continue
        for s in xl.sheet_names:
            if s in SKIP_SHEETS:
                continue
            try:
                df = xl.parse(s, header=None)
            except Exception:
                continue
            if df.shape[0] < 10:
                continue
            rows += parse_sessions(df, s, Path(f).name)

    d = pd.DataFrame(rows).dropna(subset=["측점", "날짜"])
    d = d[d[["ED", "WD", "EU", "WU"]].notna().all(axis=1)]
    d = d.drop_duplicates(subset=["측점", "날짜", "ED", "WD", "SU"], keep="first")
    return d.sort_values(["측점", "날짜", "시작"]).reset_index(drop=True)


# ------------------------------------------------------------- 내부 검증
def internal_checks(d: pd.DataFrame) -> pd.DataFrame:
    r = d.copy()

    # A. 마크 폐합
    r["A_마크폐합_gon"] = (r.Mark1_0 - r.Mark1_1).abs() - HALF_GON

    # B. 편각 재계산
    dec1 = (r.ED + r.WD + r.EU + r.WU) / 4
    r["DEC1_재계산"] = dec1
    mark_mean = r[["Mark1_0", "Mark1_1"]].mean(axis=1)
    D_gon = r.방위각_deg / GON2DEG - mark_mean + (dec1 - HALF_GON) - 100
    # 마크를 정/반 어느 쪽으로 시준했는지에 따라 결과가 180° 갈린다
    # (DEC1 / DEC2, Mark1-0 / Mark1-1 이 각각 200 gon 차이).
    # 편각은 물리적으로 |D| < 90° 이므로 그 분기를 택한다.
    dd = (D_gon * GON2DEG + 180) % 360 - 180
    r["D_재계산"] = np.where(np.abs(dd) > 90, dd - np.sign(dd) * 180, dd)
    r["B_편각차_deg"] = r.D_재계산 - r.D_야장

    # C. 복각 재계산
    I_gon = ((HALF_GON - r.SU) + (FULL_GON - r.ND)
             + (r.SD - HALF_GON) + r.NU) / 4
    r["I_재계산"] = I_gon * GON2DEG
    r["C_복각차_deg"] = r.I_재계산 - r.I_야장

    # D. 분력 정합
    r["D_분력차_nT"] = np.hypot(r.H_야장, r.Z_야장) - r.F_야장

    # E. 자오선 대칭성: ED/WU, WD/EU 는 대척(200 gon) 관계
    r["E_ED_WU_gon"] = (r.WU - r.ED).abs() - HALF_GON
    r["E_WD_EU_gon"] = (r.EU - r.WD).abs()
    return r


# ------------------------------------------------------------- 외부 검증
# '10~'19년 지자기점 관측현황(최종).xls 에서 확인한 2019 측점 좌표
SITE_COORD = {
    "미원": (36 + 43/60 + 19/3600, 127 + 35/60 + 7/3600, 97.24),
    "이원": (36 + 4/60 + 15/3600, 127 + 41/60 + 59/3600, 210.7),
    "남지": (35 + 26/60 + 10/3600, 128 + 28/60 + 25/3600, 50.0),
    "거제": (34 + 52/60 + 13/3600, 128 + 35/60 + 34/3600, 91.0),
    "장흥": (34 + 37/60 + 38/3600, 126 + 58/60 + 40/3600, 100.0),
    "성산": (33 + 31/60 + 28/3600, 126 + 53/60 + 39/3600, 30.0),
}


def external_checks(r: pd.DataFrame):
    """F. IGRF 잔차 · G. CYG 보정 후 세션 반복성 · H. 성과표 대조"""
    from lmm_cyg import fetch_kp, fetch_range
    from lmm_external import KST, derive_dif, quiet_baseline, to_utc

    r = r[r["시작"].notna()].copy()
    r["lat"] = r.측점.map(lambda s: SITE_COORD.get(s, (np.nan,)*3)[0])
    r["lon"] = r.측점.map(lambda s: SITE_COORD.get(s, (np.nan,)*3)[1])
    r["elev"] = r.측점.map(lambda s: SITE_COORD.get(s, (np.nan,)*3)[2])
    r = r.dropna(subset=["lat"])

    # --- F. IGRF-14 잔차 ---
    dates = [dt.datetime.combine(d, dt.time(12, 0)) for d in r.날짜]
    Di, Ii, Fi, *_ = igrf_dif(r.lat.values, r.lon.values, r.elev.values, dates)
    r["F_잔차_D"] = r.D_야장.values - Di
    r["F_잔차_I"] = r.I_야장.values - Ii
    r["F_잔차_F"] = r.F_야장.values - Fi

    # --- G. CYG 외부장 보정 ---
    days = sorted(set(r.날짜))
    lo, hi = min(days) - dt.timedelta(days=10), max(days) + dt.timedelta(days=10)
    cyg_parts = []
    for d in days:
        cyg_parts.append(fetch_range(d - dt.timedelta(days=10),
                                     d + dt.timedelta(days=10), quiet=True))
    cyg = pd.concat(cyg_parts, ignore_index=True)
    cyg["time"] = pd.to_datetime(cyg["time"], utc=True)
    cyg = cyg.sort_values("time").drop_duplicates("time")
    kp = fetch_kp()
    cygd = derive_dif(cyg).set_index("time")

    base_cache, corr = {}, []
    for _, s in r.iterrows():
        day = s["날짜"]
        if day not in base_cache:
            base_cache[day] = quiet_baseline(cyg, kp, day)
        b = base_cache[day]
        t0 = to_utc(day, s["시작"])
        seg = cygd[(cygd.index >= t0 - pd.Timedelta(minutes=5))
                   & (cygd.index <= t0 + pd.Timedelta(minutes=20))].dropna(subset=["X"])
        if b is None or seg.empty:
            corr.append((np.nan, np.nan, np.nan))
            continue
        vX, vY, vZ = seg.X.mean()-b["X"], seg.Y.mean()-b["Y"], seg.Z.mean()-b["Z"]
        bH = np.hypot(b["X"], b["Y"]); bF = np.hypot(bH, b["Z"])
        mH = np.hypot(b["X"]+vX, b["Y"]+vY); mF = np.hypot(mH, b["Z"]+vZ)
        dD = np.degrees(np.arctan2(b["Y"]+vY, b["X"]+vX) - np.arctan2(b["Y"], b["X"]))
        dI = np.degrees(np.arctan2(b["Z"]+vZ, mH) - np.arctan2(b["Z"], bH))
        corr.append((dD, dI, mF - bF))

    r[["cor_D", "cor_I", "cor_F"]] = pd.DataFrame(corr, index=r.index)
    r["D_보정"] = r.D_야장 - r.cor_D
    r["I_보정"] = r.I_야장 - r.cor_I
    return r


def main():
    print("=" * 78)
    print("2019년 야장 신뢰도 검증")
    print("=" * 78)

    d = collect_raw()
    print(f"\n원시 판독값이 온전한 세션 {len(d)}건 / 측점 {d.측점.nunique()}개")

    r = internal_checks(d)

    print("\n" + "-" * 78)
    print("[내부검증] 야장 자체 정합성")
    print("-" * 78)

    checks = [
        ("A 마크폐합", "A_마크폐합_gon", TOL_MARK_GON, "gon"),
        ("B 편각재계산", "B_편각차_deg", TOL_RECALC_DEG, "°"),
        ("C 복각재계산", "C_복각차_deg", TOL_RECALC_DEG, "°"),
        ("D 분력정합", "D_분력차_nT", TOL_F_NT, "nT"),
    ]
    for name, col, tol, unit in checks:
        v = r[col].abs().dropna()
        if v.empty:
            print(f"  {name:12s} 자료 없음")
            continue
        bad = (v > tol).sum()
        print(f"  {name:12s} 최대 {v.max():9.4f} {unit:3s}  "
              f"중앙 {v.median():8.4f}  허용 {tol} → "
              f"{'통과' if bad == 0 else f'{bad}건 초과'}")

    # 개별 불합격 목록
    print("\n  ── 항목별 이상 세션 ──")
    any_bad = False
    for name, col, tol, unit in checks:
        bad = r[r[col].abs() > tol]
        for _, b in bad.iterrows():
            any_bad = True
            t = b["시작"].strftime("%H:%M") if b["시작"] else "—"
            print(f"   [{name}] {b['측점']} {b['날짜']} {t}  "
                  f"편차 {b[col]:+.4f} {unit}")
    if not any_bad:
        print("   없음 — 모든 세션이 야장 내부 정합성 통과")

    # ---------------------------------------------------------- 외부 검증
    print("\n" + "-" * 78)
    print("[외부검증] 독립 자료와의 대조")
    print("-" * 78)
    e = external_checks(r)

    print("\n  ── F. IGRF-14 잔차 ──")
    for site, g in e.groupby("측점"):
        print(f"   {site:4s} ΔD {g.F_잔차_D.mean():+7.3f}°  "
              f"ΔI {g.F_잔차_I.mean():+7.3f}°  ΔF {g.F_잔차_F.mean():+8.1f} nT")

    print("\n  ── G. CYG 외부장 보정 후 세션 간 산포 (같은 날·같은 측점) ──")
    print("     보정으로 산포가 줄면 야장 시각·측정값이 함께 신뢰할 만하다는 뜻")
    from lmm_cyg import CYG_LAT, CYG_LON, fetch_kp, kp_for
    from lmm_external import haversine, to_utc
    kp_tab = fetch_kp()

    rows = []
    for (site, day), g in e.groupby(["측점", "날짜"]):
        if len(g) < 2 or g.D_보정.isna().any():
            continue
        km = haversine(g.lat.iloc[0], g.lon.iloc[0], CYG_LAT, CYG_LON)
        kps = [float(kp_for([to_utc(day, t)], kp_tab).iloc[0]) for t in g.시작]
        rows.append({
            "측점": site, "날짜": day, "세션": len(g),
            "CYG거리_km": km, "Kp최대": max(kps),
            "D원_arcmin": (g.D_야장.max()-g.D_야장.min())*60,
            "D보정_arcmin": (g.D_보정.max()-g.D_보정.min())*60,
            "I원_arcmin": (g.I_야장.max()-g.I_야장.min())*60,
            "I보정_arcmin": (g.I_보정.max()-g.I_보정.min())*60,
        })
    sp = pd.DataFrame(rows)
    if len(sp):
        print(sp.round(2).to_string(index=False))
        for label, sub in [("정온 Kp<=2", sp[sp.Kp최대 <= 2]),
                           ("교란 Kp>2", sp[sp.Kp최대 > 2])]:
            if sub.empty:
                continue
            dg = sub.D원_arcmin.mean() - sub.D보정_arcmin.mean()
            ig = sub.I원_arcmin.mean() - sub.I보정_arcmin.mean()
            print(f"\n     [{label}] {len(sub)}건")
            print(f"       D {sub.D원_arcmin.mean():.2f}' -> {sub.D보정_arcmin.mean():.2f}'"
                  f"  ({'개선' if dg > 0 else '악화'} {abs(dg):.2f}')")
            print(f"       I {sub.I원_arcmin.mean():.2f}' -> {sub.I보정_arcmin.mean():.2f}'"
                  f"  ({'개선' if ig > 0 else '악화'} {abs(ig):.2f}')")

        print(f"\n     보정 전 원자료 세션산포: D 최대 {sp.D원_arcmin.max():.2f}'"
              f" = {sp.D원_arcmin.max()/60:.4f}° (KPI 0.1°)")

    print("\n  ── H. 2022~2025 성과표와 대조 (IGRF 잔차 기준, 영년변화 제거됨) ──")
    pts = load_survey_points()
    Dp, Ip, Fp, *_ = igrf_dif(pts.lat.values, pts.lon.values,
                              pts.elev_m.values, pts.date.tolist())
    pts = pts.assign(rD=pts.D_deg.values-Dp, rI=pts.I_deg.values-Ip,
                     rF=pts.F_nT.values-Fp)
    for site, g19 in e.groupby("측점"):
        g22 = pts[pts.name == site]
        if g22.empty:
            print(f"   {site:4s} 성과표에 없음 (신규 측점)")
            continue
        print(f"   {site:4s} ΔD 2019 {g19.F_잔차_D.mean():+7.3f}° vs "
              f"{'/'.join(str(int(y)) for y in g22.year)} "
              f"{g22.rD.mean():+7.3f}°   차이 {g19.F_잔차_D.mean()-g22.rD.mean():+.3f}°")
    return e


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
