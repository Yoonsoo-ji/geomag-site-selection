#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
도엽 위치 정합성 진단 스크립트
================================
SHP 파일의 실제 WGS84 좌표를 여러 EPSG 코드로 시험 변환하여
NGII 실제 도엽 위치와 비교합니다.

사용법:
    python _check_topo_alignment.py

기준 확인 지점 (알려진 위치):
  - 예진 (36806 / NJ52-11-10): 약 36.65°N, 128.45°E (예천군 인근)
  - 안동 (36807 / NJ52-11-11): 약 36.57°N, 128.72°E
"""

import geopandas as gpd
from pathlib import Path

SHP_PATH = Path("data/국가기본도_도엽인덱스50K/TN_MAPINDX_50K.shp")

if not SHP_PATH.exists():
    print(f"[ERROR] SHP 파일 없음: {SHP_PATH}")
    raise SystemExit(1)

print("=" * 68)
print("  NGII 도엽인덱스50K — CRS 진단")
print("=" * 68)

# ── 1. 파일 로드 및 기본 정보 출력 ──────────────────────────────
raw = gpd.read_file(str(SHP_PATH), encoding="cp949")
print(f"\n[파일 정보]")
print(f"  행 수     : {len(raw)}")
print(f"  컬럼      : {raw.columns.tolist()}")
print(f"  GPD CRS   : {raw.crs}")
print(f"  EPSG 추출 : {raw.crs.to_epsg() if raw.crs else 'None'}")
print(f"  CRS name  : {raw.crs.name if raw.crs else 'None'}")

# SHP 원시 좌표 범위 (투영 좌표계 단위)
xmin, ymin, xmax, ymax = raw.total_bounds
print(f"\n[원시 좌표 범위] (투영 단위)")
print(f"  X: {xmin:.0f} ~ {xmax:.0f}")
print(f"  Y: {ymin:.0f} ~ {ymax:.0f}")

# 좌표 범위로 CRS 추정:
# EPSG:5179  — X ≈ 800,000~1,200,000 / Y ≈ 1,600,000~2,200,000
# EPSG:5186  — X ≈ 0~400,000         / Y ≈ 0~1,000,000
# EPSG:4326  — X ≈ 124~130           / Y ≈ 33~38 (도 단위)
print("\n  [CRS 추정]")
if 700_000 < xmin < 1_300_000 and 1_500_000 < ymin < 2_300_000:
    print("  → 좌표 범위가 EPSG:5179 (Korea Unified, FE=1,000,000 / FN=2,000,000) 에 해당")
elif 0 < xmin < 500_000 and 0 < ymin < 1_200_000:
    print("  → 좌표 범위가 EPSG:5186 (Korea Central Belt, FE=200,000 / FN=600,000) 에 해당")
elif 120 < xmin < 135 and 30 < ymin < 42:
    print("  → 좌표가 이미 WGS84 도 단위 (EPSG:4326)")
else:
    print(f"  → 알 수 없는 좌표 범위 — 수동 확인 필요")

# ── 2. 대상 도엽 찾기 ───────────────────────────────────────────
TARGET_MAPIDCD = ["36806", "36807"]  # 예진, 안동

id_col = None
for col in ["MAPIDCD_NO", "MAPID_NO", "MAPID_NM"]:
    if col in raw.columns:
        id_col = col
        break

print(f"\n[대상 도엽 진단] — 기준 컬럼: {id_col}")

# ── 3. 여러 EPSG로 변환하여 centroid 출력 ─────────────────────────
test_epsg = [5179, 5186, 5187, 5185, 5174, 5173]

print("\n  EPSG 별 변환 centroid (도엽 36806 / 36807):")
print(f"  {'EPSG':<8}  {'코드':>7}  {'위도(°N)':>10}  {'경도(°E)':>10}  {'평가'}")
print("  " + "-" * 60)

def _try_convert(df, epsg_in):
    try:
        r = df.set_crs(f"EPSG:{epsg_in}", allow_override=True).to_crs("EPSG:4326")
        return r
    except Exception as e:
        return None

# 기준 위치 (알려진 지점)
known = {
    "36806": (36.65, 128.45),  # 예진 (예천 인근 추정)
    "36807": (36.57, 128.72),  # 안동
}

for epsg in test_epsg:
    converted = _try_convert(raw, epsg)
    if converted is None:
        continue
    if id_col is None:
        continue
    for tgt in TARGET_MAPIDCD:
        mask = converted[id_col].astype(str) == tgt
        sub = converted[mask]
        if len(sub) == 0:
            continue
        geom = sub.iloc[0].geometry
        cx = geom.centroid.x
        cy = geom.centroid.y
        ref_lat, ref_lon = known.get(tgt, (None, None))
        if ref_lat and ref_lon:
            dlat = abs(cy - ref_lat)
            dlon = abs(cx - ref_lon)
            dist_km = ((dlat * 111)**2 + (dlon * 89)**2)**0.5
            ok = "✅ 양호" if dist_km < 5 else ("⚠ 주의" if dist_km < 20 else "❌ 오차 큼")
            note = f"참조({ref_lat:.2f}°N, {ref_lon:.2f}°E)와 {dist_km:.1f}km 차이 {ok}"
        else:
            note = ""
        print(f"  EPSG:{epsg:<4}  {tgt:>7}  {cy:>10.4f}  {cx:>10.4f}  {note}")

# ── 4. 네이티브 CRS (override 없이) ──────────────────────────────
print("\n  [네이티브 PRJ CRS 변환] (set_crs 없이 바로 to_crs):")
try:
    converted_native = raw.to_crs("EPSG:4326")
    for tgt in TARGET_MAPIDCD:
        if id_col is None:
            break
        mask = converted_native[id_col].astype(str) == tgt
        sub = converted_native[mask]
        if len(sub) == 0:
            continue
        geom = sub.iloc[0].geometry
        cx = geom.centroid.x
        cy = geom.centroid.y
        ref_lat, ref_lon = known.get(tgt, (None, None))
        if ref_lat and ref_lon:
            dlat = abs(cy - ref_lat)
            dlon = abs(cx - ref_lon)
            dist_km = ((dlat * 111)**2 + (dlon * 89)**2)**0.5
            ok = "✅ 양호" if dist_km < 5 else ("⚠ 주의" if dist_km < 20 else "❌ 오차 큼")
            print(f"  네이티브  {tgt:>7}  {cy:>10.4f}  {cx:>10.4f}  {dist_km:.1f}km 차이 {ok}")
        else:
            print(f"  네이티브  {tgt:>7}  {cy:>10.4f}  {cx:>10.4f}")
except Exception as e:
    print(f"  네이티브 변환 실패: {e}")

# ── 5. 36806 도엽 bbox 상세 출력 ─────────────────────────────────
print("\n[도엽 36806 — bbox 상세 (EPSG:5179 기준)]")
c5179 = _try_convert(raw, 5179)
if c5179 is not None and id_col is not None:
    mask = c5179[id_col].astype(str) == "36806"
    sub = c5179[mask]
    if len(sub) > 0:
        from shapely.geometry import box as sbox
        g = sub.iloc[0].geometry
        mn_lon, mn_lat, mx_lon, mx_lat = g.bounds
        print(f"  WGS84 bbox:  {mn_lat:.4f}°N ~ {mx_lat:.4f}°N, "
              f"{mn_lon:.4f}°E ~ {mx_lon:.4f}°E")
        print(f"  가로(경도): {(mx_lon - mn_lon) * 60:.2f}분 = {(mx_lon - mn_lon) * 89.5:.1f} km")
        print(f"  세로(위도): {(mx_lat - mn_lat) * 60:.2f}분 = {(mx_lat - mn_lat) * 111:.1f} km")
        print(f"  centroid:    {g.centroid.y:.4f}°N, {g.centroid.x:.4f}°E")

print("\n" + "=" * 68)
print("  결론: 위 결과에서 참조 좌표(예진 36.65°N 128.45°E)와 가장 가까운")
print("  EPSG 코드가 실제 SHP CRS입니다.")
print("  그 EPSG 코드를 create_topo_sheet_grid()의 set_crs()에 사용하세요.")
print("=" * 68)
