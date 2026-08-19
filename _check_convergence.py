# -*- coding: utf-8 -*-
"""가설: 편각 잔차의 동서 기울기 = 자오선 수차(grid convergence) 미적용/부호오류.

   γ = (λ − λ0)·sin φ   (평면직각좌표 격자북 ↔ 진북 차이)
   원점 후보: 통일원점 127.5° / 중부 127° / 권역별(서 125·중 127·동 129)
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, lmm_build as LB
from scipy.stats import pearsonr

pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
res = LB.igrf_residuals(pts)
sites = LB.aggregate_sites(pts, res)
lons, lats, grid = LB.load_kigam_grid()
cvec = LB.CrustalVector(lons, lats, *LB.crustal_vector(lons, lats, grid))
sites = LB.attach_crustal_di(sites, cvec)
y = (sites["dD"].values - sites["crD"].values) * 60      # arcmin

def gamma(lat, lon, lon0):
    return (lon - lon0) * np.sin(np.radians(lat)) * 60   # arcmin

la, lo = sites["lat"].values, sites["lon"].values
zone = np.where(lo < 126, 125.0, np.where(lo < 128, 127.0, 129.0))
cands = {"통일원점 127.5°": gamma(la, lo, 127.5),
         "중부원점 127°": gamma(la, lo, 127.0),
         "권역별 원점(125/127/129)": gamma(la, lo, zone),
         "경도만(수차 아님)": (lo - 127.5) * 60}
print(f"{'가설':>24} {'상관 r':>8} {'p':>8} {'회귀기울기':>10} {'잔차RMS 감소':>12}")
rms0 = np.sqrt((y ** 2).mean())
for k, g in cands.items():
    r, p = pearsonr(g, y)
    a = float((g @ y) / (g @ g))
    rms1 = np.sqrt(((y - a * g) ** 2).mean())
    print(f"{k:>24} {r:+8.3f} {p:8.4f} {a:10.3f} {rms0:6.1f}→{rms1:5.1f}′")

g = cands["권역별 원점(125/127/129)"]
d = pd.DataFrame({"측점": sites["name"], "위도": la.round(2), "경도": lo.round(2),
                  "원점": zone, "수차γ′": np.round(g, 1),
                  "잔차dD′": np.round(y, 1),
                  "γ차감후′": np.round(y - g, 1)})
print("\n" + d.sort_values("경도").to_string(index=False))
print(f"\n잔차 RMS {rms0:.1f}′  →  권역별 γ 1:1 차감 후 {np.sqrt(((y-g)**2).mean()):.1f}′")
