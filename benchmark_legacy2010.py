# -*- coding: utf-8 -*-
"""
2010 국가성과 대조 — 독립 표본 검증 (설계안 v2 검토 D)
=========================================================

`data/'10~'19년 지자기점 관측현황(최종).xls` 에는 30개 지자기점의 **공표 성과**
(편각·복각·전자력·수평/수직분력)가 들어 있다. 이 값은 2010 기준시점으로 환산된
것이고, 우리 LMM 은 2019~2025 관측 16측점으로 적합했다. 즉 **겹치지 않는 표본**이다.

두 가지를 본다.

  ① IGRF 단독 대비 LMM 이 그 30점에서 실제로 더 맞는가
  ② 남는 잔차의 크기·공간 구조가 우리 16점 잔차와 같은가

⚠ 한계 — 2010 원 관측자료(야장)와 그 시기 청양 자료가 우리에게 없다. INTERMAGNET
   CYG 는 2013-10-30 부터다. 따라서 2010 보고서의 **환산 절차 자체를 재현하는 것은
   불가능**하고, 여기서는 공표 성과를 독립 표본으로 쓰는 데 그친다.

    python benchmark_legacy2010.py
"""
import datetime as dt
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
import lmm_build as LB      # noqa: E402

XLS = ROOT / "data" / "'10~'19년 지자기점 관측현황(최종).xls"
EPOCH_2010 = dt.datetime(2010, 1, 1)


def dms(txt):
    """「7° 26′ 52.0000」 → 7.4478. 값이 없으면 NaN."""
    if not isinstance(txt, str):
        return np.nan
    n = re.findall(r"[-+]?\d+(?:\.\d+)?", txt.replace("\xa0", " "))
    if len(n) < 2:
        return np.nan
    v = abs(float(n[0])) + float(n[1]) / 60 + (float(n[2]) / 3600 if len(n) > 2 else 0)
    return -v if txt.strip().startswith("-") else v


def load_published():
    d = pd.ExcelFile(XLS).parse("point_excel", header=None, dtype=object)
    rows = []
    for r in range(3, d.shape[0]):
        name = str(d.iat[r, 2]).strip()
        if not name or name == "nan":
            continue
        lon, lat = dms(str(d.iat[r, 6])), dms(str(d.iat[r, 7]))
        if not (120 < lon < 132 and 32 < lat < 40):
            continue
        D = dms(str(d.iat[r, 9]))
        I = dms(str(d.iat[r, 10]))
        F = d.iat[r, 11]
        try:
            F = float(F)
        except Exception:
            F = np.nan
        h = str(d.iat[r, 8]).replace("\xa0", "").replace("m", "").strip()
        try:
            elev = float(h)
        except Exception:
            elev = 0.0
        rows.append(dict(측점=name, lat=lat, lon=lon, elev=elev,
                         # 한반도 편각은 서편각이므로 공표 도분초는 음수로 읽는다
                         D_pub=-D if D == D else np.nan,
                         I_pub=I, F_pub=F))
    return pd.DataFrame(rows)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 84)
    print("2010 공표 성과 대조 — LMM 독립 표본 검증")
    print("=" * 84)

    pub = load_published()
    ok = pub["D_pub"].notna() | pub["I_pub"].notna() | pub["F_pub"].notna()
    pub = pub[ok].reset_index(drop=True)
    print(f"공표 성과 {len(pub)}점 "
          f"(편각 {int(pub['D_pub'].notna().sum())} · "
          f"복각 {int(pub['I_pub'].notna().sum())} · "
          f"전자력 {int(pub['F_pub'].notna().sum())})")

    M = __import__("json").loads(
        (ROOT / "docs" / "data" / "lmm_model.json").read_text(encoding="utf-8"))
    R, C = M["regional"], M["crustal"]
    alpha = float(C.get("alpha", 1.0))
    gv = np.array([np.nan if x is None else x for x in C["values"]],
                  float).reshape(C["nlat"], C["nlon"])

    def crustal_at(lat, lon):
        fi = (lon - C["lon0"]) / C["dlon"]
        fj = (lat - C["lat0"]) / C["dlat"]
        i = int(np.clip(round(fi), 0, C["nlon"] - 1))
        j = int(np.clip(round(fj), 0, C["nlat"] - 1))
        v = gv[j, i]
        return 0.0 if np.isnan(v) else float(v)

    def reg(k, lat, lon):
        c = R[k]
        u, v = lat - R["lat0"], lon - R["lon0"]
        basis = [1.0, u, v]
        return float(np.dot(c, basis[:len(c)]))

    D_i, I_i, F_i, *_ = LB.igrf_dif(pub["lat"].values, pub["lon"].values,
                                    pub["elev"].values, EPOCH_2010)
    pub["D_igrf"], pub["I_igrf"], pub["F_igrf"] = D_i, I_i, F_i
    pub["D_lmm"] = [d + reg("D", la, lo)
                    for d, la, lo in zip(D_i, pub["lat"], pub["lon"])]
    pub["I_lmm"] = [i + reg("I", la, lo)
                    for i, la, lo in zip(I_i, pub["lat"], pub["lon"])]
    pub["F_lmm"] = [f + reg("F", la, lo) + alpha * crustal_at(la, lo)
                    for f, la, lo in zip(F_i, pub["lat"], pub["lon"])]

    for k, unit, sc in (("D", "′", 60), ("I", "′", 60), ("F", " nT", 1)):
        pub[f"r{k}_igrf"] = (pub[f"{k}_pub"] - pub[f"{k}_igrf"]) * sc
        pub[f"r{k}_lmm"] = (pub[f"{k}_pub"] - pub[f"{k}_lmm"]) * sc

    print()
    print(f"{'성분':>6} {'표본':>4} {'IGRF 잔차 RMS':>14} {'LMM 잔차 RMS':>14} {'개선':>8}")
    for k, unit in (("D", "′"), ("I", "′"), ("F", " nT")):
        a = pub[f"r{k}_igrf"].dropna()
        b = pub[f"r{k}_lmm"].dropna()
        if a.empty:
            continue
        ra, rb = float(np.sqrt((a ** 2).mean())), float(np.sqrt((b ** 2).mean()))
        print(f"{k:>6} {len(a):>4} {ra:>12.2f}{unit} {rb:>12.2f}{unit} "
              f"{100*(1-rb/ra):>6.1f} %")

    print("\n측점별 (편각, 분)")
    show = pub[["측점", "lat", "lon", "D_pub", "rD_igrf", "rD_lmm"]].dropna(
        subset=["rD_igrf"]).copy()
    show["D_pub"] = show["D_pub"].round(3)
    print(show.round(2).to_string(index=False))

    print("\n⚠ 해석 주의")
    print("  · 공표 성과는 2010 기준시점 환산값이고 우리 Regional 은 2019~2025")
    print("    관측에서 얻은 상수항이다. 두 시기 사이의 계통 차이가 남아 있으면")
    print("    그대로 잔차에 들어간다.")
    print("  · 2010 원 관측·당시 청양 자료가 없어 환산 절차 재현은 불가능하다.")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_2010성과_대조.csv"
    pub.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")
    return pub


if __name__ == "__main__":
    main()
