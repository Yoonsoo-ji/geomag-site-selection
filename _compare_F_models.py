# -*- coding: utf-8 -*-
"""F0 / F1 / Fα 비교 — 평가 측점을 한 번만 정해 고정한다(순환 방지)."""
import warnings; warnings.filterwarnings("ignore")
import importlib, numpy as np, pandas as pd, lmm_build as LB

DEG = 0
lons, lats, grid = LB.load_kigam_grid()
crustal = LB.CrustalGrid(lons, lats, grid)
cvec = LB.CrustalVector(lons, lats, *LB.crustal_vector(lons, lats, grid))

def build(mode, timeterm=False, vec=True):
    importlib.reload(LB)
    LB.CRUSTAL_SCALE_MODE = mode
    LB.REGIONAL_TIME_TERM = timeterm
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    return LB.attach_crustal_di(sites, cvec if vec else None)

base = build("one")
_, _, INL = LB.fit_regional(base, crustal, DEG)

def loo(sites):
    errs = {"D": [], "I": [], "F": []}
    for i in range(len(sites)):
        tr = sites.drop(sites.index[i]); te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, crustal, DEG)
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEG)
        cr = np.nan_to_num(crustal(te["lat"].values, te["lon"].values), nan=0.0)
        AD = LB.design_DI(te, A); ct = LB.crustal_term(te, cr, c)
        crD = float(te["crD"].values[0]); crI = float(te["crI"].values[0])
        if INL["D"][i]: errs["D"].append(te["dD"].values[0]-crD-(AD@c["D"])[0])
        if INL["I"][i]: errs["I"].append(te["dI"].values[0]-crI-(AD@c["I"])[0])
        if INL["F"][i]: errs["F"].append(te["dF"].values[0]-ct[0]-(AD@c["F"])[0])
    return {k: LB.rms(v) for k, v in errs.items()}

rows = []
for lab, mode, tt, vec in [("지각 벡터 적용 (현행)", "one", False, True),
                           ("지각 벡터 미적용", "one", False, False)]:
    s = build(mode, tt, vec)
    c, _, _ = LB.fit_regional(s, crustal, DEG)
    r = loo(s)
    rows.append(dict(모델=lab, alpha=round(c.get("alpha", float("nan")), 3),
                     LOO_D=round(r["D"], 4), LOO_I=round(r["I"], 4),
                     LOO_F=round(r["F"], 1)))
print("\n평가 측점 고정 · degree 0")
print(pd.DataFrame(rows).to_string(index=False))
