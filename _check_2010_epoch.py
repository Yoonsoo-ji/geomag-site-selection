# -*- coding: utf-8 -*-
"""공표 성과의 기준시점 추정 — IGRF 와 일치하는 연도를 역산한다."""
import warnings; warnings.filterwarnings("ignore")
import datetime as dt, numpy as np, pandas as pd
from benchmark_legacy2010 import load_published
import lmm_build as LB

pub = load_published()
pub = pub.dropna(subset=["D_pub"]).reset_index(drop=True)
print(f"편각 공표값 {len(pub)}점")
rows=[]
for _, r in pub.iterrows():
    best, bestd = None, 1e9
    for y in range(2008, 2027):
        D,_,_,*_ = LB.igrf_dif(np.array([r.lat]), np.array([r.lon]),
                               np.array([r.elev]), dt.datetime(y,7,1))
        d = abs(float(D[0]) - r.D_pub)*60
        if d < bestd: best, bestd = y, d
    rows.append(dict(측점=r["측점"], D공표=round(r.D_pub,4),
                     최적연도=best, 잔차분=round(bestd,1)))
d=pd.DataFrame(rows)
print(d.to_string(index=False))
print(f"\n최적연도 중앙 {d.최적연도.median():.0f} · 범위 {d.최적연도.min()}~{d.최적연도.max()}")
print(f"그 연도에서의 잔차 중앙 {d.잔차분.median():.1f}′")
