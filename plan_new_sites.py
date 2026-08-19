# -*- coding: utf-8 -*-
"""
신규 선점 배치 계획 — 어디가 비어 있고, 어느 후보가 그 구멍을 메우나
======================================================================

Regional 층의 약점은 측점 수(16)만이 아니라 **배치**다. 신뢰도 감사로 부안이 빠지면서
서해안이 통째로 비었다. 이 스크립트는

  ① 현재 측점망의 공간 공백을 격자로 계산하고(최근접 측점까지 거리),
  ② 현장조사 103개 후보 중 그 공백을 가장 크게 줄이는 순서로 순위를 매긴다.

순위는 **탐욕적 최대최소 거리 감소**로 매긴다 — 후보를 하나 넣을 때마다 국토 전역의
「최근접 측점 거리」 최댓값이 얼마나 줄어드는지로 평가한다. 등급(A/B/C/D)과 도엽
대표 여부를 함께 실어 자기환경 조건과 배치 효과를 같이 보게 한다.

    python plan_new_sites.py
"""
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs" / "output"

GRID_STEP = 0.1          # 평가 격자 (약 11 km)
TOP_N = 20               # 출력할 상위 후보 수


def geo_km(lat1, lon1, lat2, lon2):
    """WGS84 국소 평면 근사 거리(km). 한반도 규모에서 충분하다."""
    la = np.radians((lat1 + lat2) / 2)
    return np.hypot((lon2 - lon1) * 111.320 * np.cos(la),
                    (lat2 - lat1) * 110.574)


def korea_mask(lats, lons):
    """국토 경계 안쪽 격자만 남긴다(해상 제외)."""
    p = DATA / "korea_boundary.geojson"
    if not p.exists():
        return np.ones((lats.size, lons.size), bool)
    gj = json.load(open(p, encoding="utf-8"))
    from matplotlib.path import Path as MPath
    polys = []
    for f in gj["features"]:
        g = f["geometry"]
        for poly in ([g["coordinates"]] if g["type"] == "Polygon"
                     else g["coordinates"]):
            polys.append(MPath(np.asarray(poly[0])))
    LO, LA = np.meshgrid(lons, lats)
    pts = np.column_stack([LO.ravel(), LA.ravel()])
    m = np.zeros(len(pts), bool)
    for pp in polys:
        m |= pp.contains_points(pts)
    return m.reshape(LA.shape)


def main():
    # ── 현재 측점망 ──
    model = json.load(open(DATA / "lmm_model.json", encoding="utf-8"))
    sites = pd.DataFrame(model["sites"])
    print(f"현재 Regional 측점 {len(sites)}개")

    # ── 평가 격자 ──
    lats = np.arange(33.1, 38.7, GRID_STEP)
    lons = np.arange(125.0, 129.7, GRID_STEP)
    mask = korea_mask(lats, lons)
    LO, LA = np.meshgrid(lons, lats)
    gl, gn = LO[mask], LA[mask]
    print(f"평가 격자 {mask.sum():,}점 (국토 내부, {GRID_STEP}° 간격)")

    def nearest(site_lat, site_lon):
        d = np.full(gl.size, np.inf)
        for la, lo in zip(site_lat, site_lon):
            d = np.minimum(d, geo_km(gn, gl, la, lo))
        return d

    base = nearest(sites["lat"].values, sites["lon"].values)
    print(f"\n■ 현재 공백 — 최근접 측점까지 거리")
    print(f"   중앙 {np.median(base):.0f} km · 90백분위 {np.percentile(base,90):.0f} km "
          f"· 최대 {base.max():.0f} km")
    worst = np.argmax(base)
    print(f"   최대 공백 지점 {gn[worst]:.2f}°N {gl[worst]:.2f}°E "
          f"({base[worst]:.0f} km)")
    for th in (50, 75, 100):
        print(f"   {th} km 초과 지역 {100*(base>th).mean():5.1f}%")

    # ── 후보지 ──
    p = DATA / "survey_sites.geojson"
    if not p.exists():
        print("\n[중단] survey_sites.geojson 없음 — export_globe_sites.py 실행 필요")
        return
    gj = json.load(open(p, encoding="utf-8"))
    cand = pd.DataFrame([{**f["properties"],
                          "lon": f["geometry"]["coordinates"][0],
                          "lat": f["geometry"]["coordinates"][1]}
                         for f in gj["features"]])
    print(f"\n현장조사 후보 {len(cand)}개 "
          f"(등급 {cand['grade'].value_counts().to_dict()})")

    # ── 탐욕적 선택 ──
    cur = base.copy()
    chosen, rows = [], []
    pool = cand.copy().reset_index(drop=True)
    for step in range(TOP_N):
        best, bi = None, None
        for i, r in pool.iterrows():
            if i in chosen:
                continue
            d = np.minimum(cur, geo_km(gn, gl, r["lat"], r["lon"]))
            score = cur.max() - d.max()          # 최대공백 감소량
            mean_gain = cur.mean() - d.mean()
            if best is None or (score, mean_gain) > best:
                best, bi = (score, mean_gain), i
        r = pool.loc[bi]
        cur = np.minimum(cur, geo_km(gn, gl, r["lat"], r["lon"]))
        chosen.append(bi)
        rows.append({"순위": step + 1, "관리번호": r.get("mid", ""),
                     "후보지명": r.get("name", ""), "등급": r.get("grade", ""),
                     "위도": round(r["lat"], 4), "경도": round(r["lon"], 4),
                     "최대공백_감소_km": round(best[0], 1),
                     "평균거리_감소_km": round(best[1], 2),
                     "적용후_최대공백_km": round(cur.max(), 1)})

    res = pd.DataFrame(rows)
    print(f"\n■ 공백을 가장 크게 줄이는 후보 {TOP_N}개 (탐욕적 순서)")
    print(res.to_string(index=False))
    print(f"\n   {TOP_N}개 모두 반영 시 최대공백 {base.max():.0f} → {cur.max():.0f} km "
          f"· 평균 {base.mean():.0f} → {cur.mean():.0f} km")

    # 등급 A/B 만으로 제한했을 때
    ab = cand[cand["grade"].isin(["A", "B"])]
    print(f"\n   ※ 위 목록 중 등급 A·B 는 "
          f"{int(res['등급'].isin(['A','B']).sum())}/{TOP_N}개 "
          f"(전체 후보 중 A·B 는 {len(ab)}/{len(cand)})")

    OUT.mkdir(parents=True, exist_ok=True)
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    f = OUT / f"{ts}_신규선점_배치계획.csv"
    res.to_csv(f, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {f}")
    return res


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
