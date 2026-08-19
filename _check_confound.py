# -*- coding: utf-8 -*-
"""β≫1 이 시기(캠페인) 교락 때문인지 확인."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, glob
f = sorted(glob.glob("docs/output/*일변화보정_적용여부.csv"))[-1]
w = pd.read_csv(f, encoding="utf-8-sig")
w["yr_w"] = w["year"] - w.groupby("station")["year"].transform("mean")
def c(a, b):
    m = np.isfinite(w[a]) & np.isfinite(w[b])
    return float(np.corrcoef(w.loc[m, a], w.loc[m, b])[0, 1])
print("측점 내 차분 상관 (n=%d)" % len(w))
for x in ("dD_ext_w", "dI_ext_w", "dF_ext_w"):
    print(f"  {x:10} vs 연도차 : r = {c(x,'yr_w'):+.3f}")
for y in ("rD_w", "rI_w", "rF_w"):
    print(f"  {y:10} vs 연도차 : r = {c(y,'yr_w'):+.3f}")
print("\n연도차를 뺀 뒤(부분회귀) 편각 β")
for xc, yc in (("dD_ext_w","rD_w"), ("dI_ext_w","rI_w"), ("dF_ext_w","rF_w")):
    m = np.isfinite(w[xc]) & np.isfinite(w[yc])
    X = np.column_stack([w.loc[m, xc], w.loc[m, "yr_w"]])
    y = w.loc[m, yc].to_numpy(float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(f"  {yc}: β(외부장) = {b[0]:+.2f} · β(연도) = {b[1]:+.2f}/년")
print("\n세션 편차 자체의 크기")
for x, u in (("dD_ext_w","′"), ("dI_ext_w","′"), ("dF_ext_w"," nT")):
    print(f"  {x:10} 측점 내 산포 {w[x].std():.2f}{u}")
