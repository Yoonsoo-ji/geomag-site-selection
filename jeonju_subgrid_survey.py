#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
전주(완주) P1 후보지 1km 서브격자 현장 답사 세부 후보 선정
site_id=11  (35.592°N, 127.158°E)  10km 셀 → 1km 서브격자

흐름:
  1. 1km 서브격자 생성 (±5km 반경)
  2. Open-Elevation API 경사도 산정
  3. 제외구역 WKB 캐시 로드 → 불적합점 제거
  4. 도로 접근성 필터 (≤ 2km)
  5. 점수 산정 (s2~s6, 100점 정규화)
  6. Folium 지도 + CSV 저장
"""
import json, warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point, LineString, box
from shapely.prepared import prep
from pathlib import Path
import requests
import folium
from folium.plugins import MarkerCluster
from datetime import datetime

warnings.filterwarnings("ignore")

# ─── 경로 ─────────────────────────────────────────────────────
MAIN_DIR  = Path("C:/LG_gram_backup_users/LX/2026_geomag")
DATA_DIR  = MAIN_DIR / "data"
CACHE_DIR = DATA_DIR / "zone_cache"
OUT_DIR   = Path("docs/output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

UTM_CRS   = "EPSG:5179"
WGS84_CRS = "EPSG:4326"
TS        = datetime.now().strftime("%Y%m%d_%H%M%S")

# ─── 분석 파라미터 ────────────────────────────────────────────
CENTER_LAT, CENTER_LON = 35.59196, 127.15795   # site_id=11
RADIUS_M   = 5_000    # 탐색 반경
GRID_STEP  = 1_000    # 서브격자 간격
ROAD_MAX_M = 2_000    # 차량 접근 최대 거리

print("=" * 62)
print("  전주(완주) P1 후보지 — 1km 서브격자 현장 답사 세부 선정")
print(f"  중심: {CENTER_LAT}°N, {CENTER_LON}°E  /  반경 {RADIUS_M//1000}km  간격 {GRID_STEP//1000}km")
print("=" * 62)


# ═══════════════════════════════════════════════════════════════
# 1. 서브격자 생성
# ═══════════════════════════════════════════════════════════════
print("\n[1] 서브격자 생성...")
center_gdf = gpd.GeoDataFrame(geometry=[Point(CENTER_LON, CENTER_LAT)],
                               crs=WGS84_CRS).to_crs(UTM_CRS)
cx, cy = center_gdf.geometry.iloc[0].x, center_gdf.geometry.iloc[0].y

xs = np.arange(cx - RADIUS_M, cx + RADIUS_M + 1, GRID_STEP)
ys = np.arange(cy - RADIUS_M, cy + RADIUS_M + 1, GRID_STEP)
center_pt = Point(cx, cy)
pts_utm   = [Point(x, y) for x in xs for y in ys
             if Point(x, y).distance(center_pt) <= RADIUS_M]
n_all     = len(pts_utm)

gdf_all = gpd.GeoDataFrame({"gid": range(n_all)}, geometry=pts_utm, crs=UTM_CRS)
gdf_wgs = gdf_all.to_crs(WGS84_CRS)
gdf_all["lat"] = gdf_wgs.geometry.y.round(6)
gdf_all["lon"] = gdf_wgs.geometry.x.round(6)
print(f"  격자점: {n_all}개")


# ═══════════════════════════════════════════════════════════════
# 2. 지형 경사도 (Open-Elevation API)
# ═══════════════════════════════════════════════════════════════
print("\n[2] 지형 경사도 산정 (Open-Elevation API)...")
ELEV_CACHE = OUT_DIR / "jeonju_elev_cache.json"
ELEV_OFFSET = 500   # 4방향 오프셋(m) — 중앙 차분 경사도

def fetch_elev_batch(latlons):
    url = "https://api.open-elevation.com/api/v1/lookup"
    results = [None] * len(latlons)
    for start in range(0, len(latlons), 1000):
        batch  = latlons[start:start + 1000]
        payload = {"locations": [{"latitude": la, "longitude": lo}
                                  for la, lo in batch]}
        try:
            r = requests.post(url, json=payload, timeout=90)
            if r.status_code == 200:
                for j, item in enumerate(r.json()["results"]):
                    results[start + j] = item.get("elevation")
        except Exception as e:
            print(f"  ⚠ API 오류: {e}")
    return results

DIRS = {"C": (0, 0), "N": (0, ELEV_OFFSET),
        "S": (0, -ELEV_OFFSET), "E": (ELEV_OFFSET, 0), "W": (-ELEV_OFFSET, 0)}

if ELEV_CACHE.exists():
    print(f"  캐시 로드: {ELEV_CACHE.name}")
    with open(ELEV_CACHE) as f:
        elev_dict = json.load(f)
    elev_dict = {int(k): v for k, v in elev_dict.items()}
else:
    # 모든 격자점 × 5방향 UTM → WGS84 배치 변환
    all_pts_utm, gi_list, dir_list = [], [], []
    for gi in range(n_all):
        xc = gdf_all.geometry.iloc[gi].x
        yc = gdf_all.geometry.iloc[gi].y
        for dname, (dx, dy) in DIRS.items():
            all_pts_utm.append(Point(xc + dx, yc + dy))
            gi_list.append(gi)
            dir_list.append(dname)

    pts_conv = gpd.GeoDataFrame(geometry=all_pts_utm, crs=UTM_CRS).to_crs(WGS84_CRS)
    latlons_batch = [(row.geometry.y, row.geometry.x)
                     for _, row in pts_conv.iterrows()]

    print(f"  API 요청: {len(latlons_batch)}개 지점...")
    elevs = fetch_elev_batch(latlons_batch)

    elev_dict = {}
    for idx, (gi, dname) in enumerate(zip(gi_list, dir_list)):
        if gi not in elev_dict:
            elev_dict[gi] = {}
        elev_dict[gi][dname] = elevs[idx]

    with open(ELEV_CACHE, "w") as f:
        json.dump(elev_dict, f)
    print(f"  캐시 저장: {ELEV_CACHE.name}")

slopes = np.full(n_all, np.nan)
for gi in range(n_all):
    t = elev_dict.get(gi, {})
    ec = t.get("C")
    grads = []
    for dname in ["N", "S", "E", "W"]:
        ev = t.get(dname)
        if ec is not None and ev is not None:
            grads.append(np.degrees(np.arctan(abs(ev - ec) / ELEV_OFFSET)))
    if grads:
        slopes[gi] = np.mean(grads)

gdf_all["slope_deg"] = np.round(slopes, 2)
print(f"  경사도 취득: {(~np.isnan(slopes)).sum()}/{n_all}개, "
      f"평균 {np.nanmean(slopes):.1f}°")


# ═══════════════════════════════════════════════════════════════
# 3. 제외구역 WKB 캐시 로드 → 필터링
# ═══════════════════════════════════════════════════════════════
print("\n[3] 제외구역 필터링 (WKB 캐시 로드)...")
LOCAL_BOX = box(cx - RADIUS_M * 1.2, cy - RADIUS_M * 1.2,
                cx + RADIUS_M * 1.2, cy + RADIUS_M * 1.2)

ZONE_NAMES = ["power", "railway", "urban_dense", "urban_resid",
              "pipeline", "comm", "wind", "quarry",
              "water", "geology", "fault", "anomaly"]

def load_zone(name):
    p = CACHE_DIR / f"{name}.wkb"
    if not p.exists():
        return None
    g = shapely.from_wkb(p.read_bytes())
    try:
        local = g.intersection(LOCAL_BOX)
        return None if local.is_empty else local
    except Exception:
        return g   # fallback: 전체 사용

zones = {}
for nm in ZONE_NAMES:
    zones[nm] = load_zone(nm)
    status = "없음" if zones[nm] is None else "로드"
    print(f"  [{nm:12s}] {status}")

excl_mask  = np.zeros(n_all, dtype=bool)
excl_reason = [""] * n_all

for nm in ZONE_NAMES:
    z = zones[nm]
    if z is None or z.is_empty:
        continue
    pz = prep(z)
    cnt = 0
    for i, pt in enumerate(pts_utm):
        if not excl_mask[i] and pz.contains(pt):
            excl_mask[i] = True
            excl_reason[i] = nm
            cnt += 1
    if cnt > 0:
        print(f"  → [{nm}] {cnt}개 제외")

gdf_all["excl"]   = excl_mask
gdf_all["excl_r"] = excl_reason
candidates = gdf_all[~excl_mask].copy().reset_index(drop=True)
print(f"\n  필터 후 후보점: {len(candidates)}개 / {n_all}개")


# ═══════════════════════════════════════════════════════════════
# 4. 도로 접근성 (access_roads.json)
# ═══════════════════════════════════════════════════════════════
print("\n[4] 도로 접근성 계산...")
road_json = DATA_DIR / "access_roads.json"
road_dist = np.full(n_all, np.nan)

if road_json.exists():
    with open(road_json, encoding="utf-8") as f:
        rd = json.load(f)
    nodes_rd = {e["id"]: (e["lon"], e["lat"]) for e in rd.get("elements", [])
                if e["type"] == "node"}
    road_lines = []
    for e in rd.get("elements", []):
        if e["type"] == "way":
            coords = [nodes_rd[n] for n in e.get("nodes", []) if n in nodes_rd]
            if len(coords) >= 2:
                road_lines.append(LineString(coords))

    if road_lines:
        road_gdf = gpd.GeoDataFrame(geometry=road_lines,
                                     crs=WGS84_CRS).to_crs(UTM_CRS)
        # 로컬 영역 내 도로만 clip
        local_roads_gdf = road_gdf[road_gdf.geometry.intersects(LOCAL_BOX)]
        print(f"  로컬 도로 세그먼트: {len(local_roads_gdf)}개")

        if len(local_roads_gdf) > 0:
            road_union = local_roads_gdf.geometry.unary_union
            p_road = prep(road_union)
            for i, pt in enumerate(pts_utm):
                road_dist[i] = pt.distance(road_union)
        else:
            print("  ⚠ 로컬 영역 내 도로 없음 — 접근성 필터 건너뜀")
else:
    print("  ⚠ access_roads.json 없음")

gdf_all["road_dist_m"] = np.round(road_dist, 0)

# 후보점에 도로 거리 병합
c_idx = candidates.index
candidates["road_dist_m"] = gdf_all.loc[c_idx, "road_dist_m"].values

# 도로 접근성 필터 적용
if not np.all(np.isnan(road_dist)):
    before = len(candidates)
    candidates = candidates[
        candidates["road_dist_m"].isna() |
        (candidates["road_dist_m"] <= ROAD_MAX_M)
    ].copy().reset_index(drop=True)
    print(f"  도로 ≤ {ROAD_MAX_M}m 필터: {before}개 → {len(candidates)}개")
else:
    print("  도로 거리 미산정 — 접근성 필터 생략")


# ═══════════════════════════════════════════════════════════════
# 5. 점수 산정 (s2 지형 / s3 전력철도 / s4 인구 / s5 KIGAM / s6 암상)
# ═══════════════════════════════════════════════════════════════
print("\n[5] 점수 산정...")
nc = len(candidates)
pts_c = candidates.geometry.tolist()

# ── ② 지형 경사도 (15점) ─────────────────────────────────
sl = candidates["slope_deg"].values.astype(float)
sl_filled = np.where(np.isnan(sl), np.nanmedian(sl) if not np.all(np.isnan(sl)) else 5.0, sl)
s2 = np.clip((1.0 - sl_filled / 30.0) * 15.0, 0.0, 15.0)
candidates["s2_지형"] = np.round(s2, 1)
print(f"  ② 지형: 평균 {s2.mean():.1f}/15점")

# ── ③ 전력/철도 이격도 (15점) ────────────────────────────
def zone_dist(pts, z, label=""):
    if z is None or z.is_empty:
        return np.full(len(pts), 50_000.0)
    zs = z.simplify(300)
    d  = np.array([pt.distance(zs) for pt in pts])
    return d

d_power   = zone_dist(pts_c, zones["power"],   "전력")
d_railway = zone_dist(pts_c, zones["railway"],  "철도")
lp = np.log1p(d_power)
lr = np.log1p(d_railway)
s3 = (lp / (lp.max() or 1) * 0.5 + lr / (lr.max() or 1) * 0.5) * 15
candidates["s3_전력철도"] = np.round(s3, 1)
candidates["d_power_km"]   = np.round(d_power   / 1000, 1)
candidates["d_railway_km"] = np.round(d_railway / 1000, 1)
print(f"  ③ 전력/철도: 평균 {s3.mean():.1f}/15점")

# ── ④ 인구 밀집 이격도 (15점) ────────────────────────────
d_urban = zone_dist(pts_c, zones["urban_resid"], "도시")
lu = np.log1p(d_urban)
s4 = (lu / (lu.max() or 1)) * 15
candidates["s4_인구이격"] = np.round(s4, 1)
candidates["d_urban_km"]  = np.round(d_urban / 1000, 1)
print(f"  ④ 인구이격: 평균 {s4.mean():.1f}/15점")

# ── ⑤ 자기이상·지질경계 모델 기여도 (10점) ──────────────
print("  ⑤ KIGAM 자기이상 모델 기여도...")
kigam_path = DATA_DIR / "mag_1982-2018_1.5min_ed.dat"
s5       = np.full(nc, 6.0)
p90p10_arr = np.full(nc, np.nan)
if kigam_path.exists():
    kdf = pd.read_csv(kigam_path, sep=r"\s+", skiprows=9, header=None,
                      names=["lon", "lat", "anomaly_nT"],
                      na_values=["99999", "-99999"], on_bad_lines="skip")
    kdf = kdf.apply(pd.to_numeric, errors="coerce").dropna()
    klat = kdf.lat.values; klon = kdf.lon.values; kanom = kdf.anomaly_nT.values
    cands_wgs_c = candidates.to_crs(WGS84_CRS)
    for i in range(nc):
        clat = cands_wgs_c.geometry.iloc[i].y
        clon = cands_wgs_c.geometry.iloc[i].x
        mask = (np.abs(klat - clat) <= 0.05) & (np.abs(klon - clon) <= 0.05)
        pts_k = kanom[mask]
        if len(pts_k) >= 3:
            v = np.percentile(pts_k, 90) - np.percentile(pts_k, 10)
            p90p10_arr[i] = round(v, 1)
            if   v <  30:  s5[i] = 5.0
            elif v <= 150: s5[i] = 8.0
            elif v <= 400: s5[i] = 10.0
            elif v <= 800: s5[i] = 7.0
            else:          s5[i] = 0.0
candidates["s5_모델기여도"] = np.round(s5, 1)
candidates["mag_p90p10_nT"] = np.round(p90p10_arr, 1)
print(f"  ⑤ 모델기여도: 평균 {s5.mean():.1f}/10점")

# ── ⑥ 암상 적합성 (5점) ─────────────────────────────────
s6 = np.full(nc, 2.5)
gz   = zones.get("geology")
fz   = zones.get("fault")
have_rock  = gz   is not None and not gz.is_empty
have_fault = fz   is not None and not fz.is_empty

if have_rock:
    d_rock  = zone_dist(pts_c, gz.simplify(300) if have_rock else None, "암상")
    lr_rock = np.log1p(d_rock)
    candidates["d_geo_rock_km"] = np.round(d_rock / 1000, 1)
else:
    lr_rock = np.ones(nc)
    candidates["d_geo_rock_km"] = np.nan

if have_fault:
    d_fault  = zone_dist(pts_c, fz.simplify(300) if have_fault else None, "단층")
    lr_fault = np.log1p(d_fault)
    candidates["d_geo_fault_km"] = np.round(d_fault / 1000, 1)
else:
    lr_fault = np.ones(nc)
    candidates["d_geo_fault_km"] = np.nan

if have_rock and have_fault:
    s6 = (lr_rock / (lr_rock.max() or 1) * 0.5 +
          lr_fault / (lr_fault.max() or 1) * 0.5) * 5
elif have_rock:
    s6 = (lr_rock / (lr_rock.max() or 1)) * 5
elif have_fault:
    s6 = (lr_fault / (lr_fault.max() or 1)) * 5

candidates["s6_암상"] = np.round(s6, 1)
print(f"  ⑥ 암상: 평균 {s6.mean():.1f}/5점")

# ── 종합 점수 (60점 만점 → 100점 정규화) ─────────────────
raw = s2 + s3 + s4 + s5 + s6    # 최대 60점
candidates["score_raw"] = np.round(raw, 1)
candidates["score"]     = np.round(raw / 60 * 100, 1)
candidates = candidates.sort_values("score", ascending=False).reset_index(drop=True)
candidates["rank"] = candidates.index + 1

print(f"\n  최종 후보: {len(candidates)}개")
print(f"  점수 범위: {candidates['score'].min():.1f} ~ {candidates['score'].max():.1f}점")
print(f"\n  ── 상위 10개 ──────────────────────────────────────")
top10 = candidates.head(10)[["rank","lat","lon","score",
                              "s2_지형","s3_전력철도","s4_인구이격",
                              "s5_모델기여도","s6_암상",
                              "slope_deg","d_power_km","d_railway_km",
                              "d_urban_km","road_dist_m"]]
print(top10.to_string(index=False))


# ═══════════════════════════════════════════════════════════════
# 6. Folium 지도 생성
# ═══════════════════════════════════════════════════════════════
print("\n[6] Folium 지도 생성...")

# 점수 → 색상 (red-yellow-green)
def score_color(s, mn, mx):
    t = (s - mn) / (mx - mn + 1e-9)
    r = int((1 - t) * 220)
    g = int(t * 200 + 40)
    return f"#{r:02x}{g:02x}32"

smin, smax = candidates["score"].min(), candidates["score"].max()

m = folium.Map(location=[CENTER_LAT, CENTER_LON], zoom_start=12,
               tiles="CartoDB Positron")

# 5km 반경 원
folium.Circle(
    location=[CENTER_LAT, CENTER_LON],
    radius=RADIUS_M,
    color="#4472C4", fill=True, fill_opacity=0.05,
    weight=2, dash_array="8",
    tooltip="탐색 반경 5km"
).add_to(m)

# P1 원점 마커
folium.Marker(
    location=[CENTER_LAT, CENTER_LON],
    icon=folium.Icon(color="blue", icon="star", prefix="fa"),
    tooltip="P1 원점 (site_id=11, 도상 점수 61.1점)",
    popup=folium.Popup(
        f"<b>P1 원점 (site_id=11)</b><br>"
        f"도상 종합점수: 61.1점<br>"
        f"위도: {CENTER_LAT}°N<br>"
        f"경도: {CENTER_LON}°E",
        max_width=220
    )
).add_to(m)

# 제외된 격자점 (회색 소형 원)
excluded = gdf_all[gdf_all["excl"]].copy()
excl_wgs = excluded.to_crs(WGS84_CRS)
for _, row in excl_wgs.iterrows():
    folium.CircleMarker(
        location=[row.geometry.y, row.geometry.x],
        radius=3, color="#BBBBBB", fill=True, fill_opacity=0.4,
        weight=0.5,
        tooltip=f"제외: {row['excl_r']}"
    ).add_to(m)

# 도로 접근성 불만족 (도달 불가)
if "road_dist_m" in gdf_all.columns:
    road_far = gdf_all[
        (~gdf_all["excl"]) &
        (gdf_all["road_dist_m"].notna()) &
        (gdf_all["road_dist_m"] > ROAD_MAX_M)
    ].to_crs(WGS84_CRS)
    for _, row in road_far.iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=4, color="#FF8800", fill=True, fill_opacity=0.4,
            weight=1,
            tooltip=f"도로 {row['road_dist_m']:.0f}m (>{ROAD_MAX_M}m 접근 불리)"
        ).add_to(m)

# 후보점 (점수별 색상)
cands_wgs = candidates.to_crs(WGS84_CRS)
top5_rank = set(range(1, 6))

for idx, row in cands_wgs.iterrows():
    lat, lon = row.geometry.y, row.geometry.x
    sc  = row["score"]
    rnk = row["rank"]
    col = score_color(sc, smin, smax)

    # 슬로프, 거리 정보 포맷
    sl_str   = f"{row['slope_deg']:.1f}°" if not np.isnan(row['slope_deg']) else "-"
    rd_str   = f"{row['road_dist_m']:.0f}m" if not np.isnan(row.get('road_dist_m', np.nan)) else "-"
    p90_str  = f"{row['mag_p90p10_nT']:.0f}nT" if not np.isnan(row.get('mag_p90p10_nT', np.nan)) else "-"
    rock_str = f"{row.get('d_geo_rock_km', np.nan):.1f}km" if not np.isnan(row.get('d_geo_rock_km', np.nan)) else "-"
    flt_str  = f"{row.get('d_geo_fault_km', np.nan):.1f}km" if not np.isnan(row.get('d_geo_fault_km', np.nan)) else "-"

    popup_html = f"""
<div style='font-family:맑은고딕,sans-serif;font-size:12px;min-width:240px'>
<b style='font-size:14px'>#{rnk} 후보점</b>&nbsp;
<span style='background:{col};color:#fff;padding:2px 8px;border-radius:10px;font-weight:bold'>{sc:.1f}점</span>
<hr style='margin:4px 0'>
<table style='border-collapse:collapse;width:100%'>
<tr><td style='color:#555'>위도</td><td><b>{lat:.5f}°N</b></td>
    <td style='color:#555'>경도</td><td><b>{lon:.5f}°E</b></td></tr>
<tr><td colspan=4><hr style='margin:3px 0'></td></tr>
<tr style='background:#f5f5f5'><td colspan=4 style='font-weight:bold;padding:2px'>점수 내역</td></tr>
<tr><td>② 지형경사</td><td>{row['s2_지형']:.1f}/15</td>
    <td>경사도</td><td>{sl_str}</td></tr>
<tr><td>③ 전력·철도</td><td>{row['s3_전력철도']:.1f}/15</td>
    <td>전력/철도</td><td>{row['d_power_km']:.1f}/{row['d_railway_km']:.1f}km</td></tr>
<tr><td>④ 인구이격</td><td>{row['s4_인구이격']:.1f}/15</td>
    <td>도시</td><td>{row['d_urban_km']:.1f}km</td></tr>
<tr><td>⑤ 모델기여</td><td>{row['s5_모델기여도']:.1f}/10</td>
    <td>P90-P10</td><td>{p90_str}</td></tr>
<tr><td>⑥ 암상</td><td>{row['s6_암상']:.1f}/5</td>
    <td>암상/단층</td><td>{rock_str}/{flt_str}</td></tr>
<tr style='background:#EEF4FF'><td colspan=2><b>합계</b></td>
    <td colspan=2><b>{row['score_raw']:.1f}/60점 → {sc:.1f}점</b></td></tr>
</table>
<hr style='margin:4px 0'>
<span style='color:#1a5276'>🚗 도로접근: {rd_str}</span>
</div>"""

    if rnk in top5_rank:
        folium.Marker(
            location=[lat, lon],
            icon=folium.DivIcon(
                html=f"""<div style='
                    background:{col};color:#fff;font-weight:bold;
                    border-radius:50%;width:32px;height:32px;
                    display:flex;align-items:center;justify-content:center;
                    font-size:14px;border:2px solid #333;
                    box-shadow:0 2px 6px rgba(0,0,0,0.4)'>#{rnk}</div>""",
                icon_size=(32, 32), icon_anchor=(16, 16)
            ),
            popup=folium.Popup(popup_html, max_width=300),
            tooltip=f"#{rnk} 종합 {sc:.1f}점"
        ).add_to(m)
    else:
        folium.CircleMarker(
            location=[lat, lon],
            radius=7, color=col, fill=True, fill_color=col,
            fill_opacity=0.8, weight=1.5,
            popup=folium.Popup(popup_html, max_width=300),
            tooltip=f"#{rnk} {sc:.1f}점"
        ).add_to(m)

# 범례
legend_html = """
<div style='position:fixed;bottom:28px;left:20px;z-index:9999;
     background:white;padding:12px 16px;border-radius:8px;
     border:1px solid #ccc;font-size:12px;font-family:맑은고딕,sans-serif;
     box-shadow:0 2px 8px rgba(0,0,0,0.2)'>
<b>범례</b><br>
<span style='color:#1F5C96'>●</span> 탐색 반경 5km<br>
<span style='color:#BBBBBB'>●</span> 제외 격자점 (간섭·지질)<br>
<span style='color:#FF8800'>●</span> 도로접근 불리 (>2km)<br>
<span style='color:#3a3'>●</span> 후보점 (점수 高)<br>
<span style='color:#d33'>●</span> 후보점 (점수 低)<br>
<b style='color:#1a5276'>#1~#5</b> 최우선 답사 후보<br>
<span style='color:#4472C4'>★</span> P1 원점 (도상 선정)
</div>"""
m.get_root().html.add_child(folium.Element(legend_html))

# 제목
title_html = f"""
<div style='position:fixed;top:12px;left:50%;transform:translateX(-50%);
     z-index:9999;background:rgba(31,56,100,0.9);color:white;
     padding:8px 20px;border-radius:6px;font-size:14px;font-weight:bold;
     font-family:맑은고딕,sans-serif;white-space:nowrap'>
전주(완주) P1 후보지 — 1km 서브격자 현장 답사 후보 |
후보 {len(candidates)}개 | 상위 5개 번호 표시
</div>"""
m.get_root().html.add_child(folium.Element(title_html))

# 레이어 컨트롤
folium.LayerControl().add_to(m)

# ─── 저장 ─────────────────────────────────────────────────────
html_path = OUT_DIR / f"{TS}_jeonju_subgrid_survey.html"
csv_path  = OUT_DIR / f"{TS}_jeonju_subgrid_candidates.csv"

m.save(str(html_path))
print(f"  지도 저장: {html_path.name}")

# CSV
save_cols = ["rank", "lat", "lon", "score", "score_raw",
             "s2_지형", "s3_전력철도", "s4_인구이격",
             "s5_모델기여도", "s6_암상", "mag_p90p10_nT",
             "slope_deg", "d_power_km", "d_railway_km", "d_urban_km",
             "d_geo_rock_km", "d_geo_fault_km", "road_dist_m"]
save_cols = [c for c in save_cols if c in candidates.columns]
candidates[save_cols].to_csv(csv_path, index=False, encoding="utf-8-sig")
print(f"  CSV  저장: {csv_path.name}")

print("\n" + "=" * 62)
print("  완료!")
print(f"  전체 격자: {n_all}개 → 제외 후: {n_all - excl_mask.sum()}개 "
      f"→ 도로접근 필터: {len(candidates)}개")
print("  ── 최종 상위 5개 답사 후보 ──────────────────────────────")
for _, r in candidates.head(5).iterrows():
    sl_v = f"{r['slope_deg']:.1f}°" if not np.isnan(r['slope_deg']) else "-"
    print(f"  #{int(r['rank']):2d} | {r['lat']:.5f}°N {r['lon']:.5f}°E | "
          f"점수 {r['score']:.1f}점 | 경사 {sl_v} | "
          f"전력 {r['d_power_km']:.1f}km | 철도 {r['d_railway_km']:.1f}km")
print("=" * 62)
