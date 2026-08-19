# -*- coding: utf-8 -*-
"""완전 공정 비교 — 평가 측점 집합을 «기준 설정»에서 한 번만 정하고 고정한다.

기준 설정: 지각 벡터 미적용 · External 미적용 · 0차. 이 집합은 어느 설정에도
유리하지 않으며, 모든 설정을 같은 측점에서 견준다.
"""
import warnings; warnings.filterwarnings("ignore")
import importlib, numpy as np, pandas as pd
import lmm_build as LB

DEG = 0
lons, lats, grid = LB.load_kigam_grid()
crustal = LB.CrustalGrid(lons, lats, grid)
cvec = LB.CrustalVector(lons, lats, *LB.crustal_vector(lons, lats, grid))


def build(mode, source, vec):
    importlib.reload(LB)
    LB.EXTERNAL_MODE = mode
    LB.EXTERNAL_CSV = (LB.EXTERNAL_CSV_MULTI if source == "multi"
                       else LB.EXTERNAL_CSV_CYG)
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    return LB.attach_crustal_di(sites, cvec if vec else None)


base = build("none", "cyg", False)
_, _, INL = LB.fit_regional(base, crustal, DEG)
names = list(base["name"])
print("고정 평가 집합:", {k: int(v.sum()) for k, v in INL.items()})
print("  D 평가 측점:", ", ".join(n for n, m in zip(names, INL["D"]) if m))


def loo(sites):
    errs = {"D": [], "I": [], "F": []}
    for i in range(len(sites)):
        tr = sites.drop(sites.index[i]); te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, crustal, DEG)
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEG)
        cr = np.nan_to_num(crustal(te["lat"].values, te["lon"].values), nan=0.0)
        if INL["D"][i]:
            errs["D"].append(te["dD"].values[0]-te["crD"].values[0]-(A@c["D"])[0])
        if INL["I"][i]:
            errs["I"].append(te["dI"].values[0]-te["crI"].values[0]-(A@c["I"])[0])
        if INL["F"][i]:
            errs["F"].append(te["dF"].values[0]-cr[0]-(A@c["F"])[0])
    return {k: LB.rms(v) for k, v in errs.items()}


rows = []
for vec in (False, True):
    for mode, src in [("none", "cyg"), ("subtract_D", "cyg"),
                      ("subtract_D", "multi"), ("subtract_quiet_D", "cyg"),
                      ("subtract_quiet_D", "multi")]:
        r = loo(build(mode, src, vec))
        rows.append(dict(지각벡터="적용" if vec else "F만", External=mode,
                         자료원=src, **{k: round(v, 4) for k, v in r.items()}))
df = pd.DataFrame(rows)
print("\n" + df.to_string(index=False))
print("\nD 최소:", df.loc[df.D.idxmin()].to_dict())
print("F 최소:", df.loc[df.F.idxmin()].to_dict())
