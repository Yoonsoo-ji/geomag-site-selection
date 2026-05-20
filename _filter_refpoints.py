#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KML 내 기준점 중 제외구역에 포함된 것을 제거
- 제외구역: power(1km), railway(5km), urban_dense(500m), urban_resid(300m),
            pipeline(500m), comm(500m), wind(500m), quarry(1km), anomaly(WKB)
- 입력 : 20260520_170218_P1_63sites_subgrid_기준점_측정점.kml
- 출력 : {TS}_P1_63sites_filtered.kml
"""
import sys, json, math, warnings, re
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime

MAIN_DIR  = Path("C:/LG_gram_backup_users/LX/2026_geomag")
DATA_DIR  = MAIN_DIR / "data"
CACHE_DIR = DATA_DIR / "zone_cache"
OUT_DIR   = MAIN_DIR / "docs" / "output"

KML_IN  = OUT_DIR / "20260520_170218_P1_63sites_subgrid_기준점_측정점.kml"
TS      = datetime.now().strftime("%Y%m%d_%H%M%S")
KML_OUT = OUT_DIR / f"{TS}_P1_63sites_filtered.kml"

WGS84 = "EPSG:4326"
UTM   = "EPSG:5179"

# 버퍼 거리 (m)
BUFFERS = {
    "power":       1_000,
    "railway":     5_000,
    "urban_dense":   500,
    "urban_resid":   300,
    "pipeline":      500,
    "comm":          500,
    "wind":          500,
    "quarry":      1_000,
    "anomaly":         0,  # WKB는 이미 버퍼 포함
}

# JSON 파일 매핑
JSON_FILES = {
    "railway":     DATA_DIR / "railways.json",
    "urban_dense": DATA_DIR / "urban_dense.json",
    "urban_resid": DATA_DIR / "urban_residential_v2.json",
    "pipeline":    DATA_DIR / "pipelines.json",
    "comm":        DATA_DIR / "comm_towers.json",
    "wind":        DATA_DIR / "wind_turbines.json",
    "quarry":      DATA_DIR / "quarries_mines.json",
}

NS  = "http://www.opengis.net/kml/2.2"
GX  = "http://www.google.com/kml/ext/2.2"
ATOM= "http://www.w3.org/2005/Atom"

print("=" * 64)
print("  기준점 제외구역 필터링")
print("=" * 64)


# ════════════════════════════════════════════════════════════════
# 1. KML 파싱 — 기준점 Placemark 추출
# ════════════════════════════════════════════════════════════════
print("\n[1] KML 파싱 중...")
ET.register_namespace("", NS)
ET.register_namespace("gx", GX)
ET.register_namespace("atom", ATOM)

def ktag(n): return f"{{{NS}}}{n}"

tree = ET.parse(KML_IN)
root = tree.getroot()

# 기준점 폴더 찾기 (📐 근접 기준점) — iter()로 전체 트리 탐색
ref_folders = []  # folder_el 목록
for folder_el in root.iter(ktag("Folder")):
    nm_el = folder_el.find(ktag("name"))
    nm = (nm_el.text or "") if nm_el is not None else ""
    if "근접 기준점" in nm:
        ref_folders.append(folder_el)
print(f"  기준점 폴더 수: {len(ref_folders)}개")

# Placemark 정보 수집 (좌표 + 元素 참조)
pm_records = []  # {folder_el, pm_el, lon, lat, idx}
for folder_el in ref_folders:
    for pm_el in folder_el.findall(ktag("Placemark")):
        coord_el = pm_el.find(f".//{ktag('coordinates')}")
        if coord_el is None or not coord_el.text: continue
        parts = coord_el.text.strip().split(",")
        if len(parts) < 2: continue
        try:
            lon, lat = float(parts[0]), float(parts[1])
        except:
            continue
        pm_records.append({
            "folder_el": folder_el,
            "pm_el":     pm_el,
            "lon": lon,
            "lat": lat,
        })

total_pm = len(pm_records)
print(f"  기준점 Placemark 총수: {total_pm}개")

# UTM GeoDataFrame 생성
pts_gdf = gpd.GeoDataFrame(
    {"idx": range(total_pm)},
    geometry=[Point(r["lon"], r["lat"]) for r in pm_records],
    crs=WGS84
).to_crs(UTM)

excluded = set()  # 제외할 pm_records 인덱스


# ════════════════════════════════════════════════════════════════
# 2. 제외구역 로드 및 필터
# ════════════════════════════════════════════════════════════════
print("\n[2] 제외구역 로드 및 교차 검사...")

def check_zone_wkb(name, wkb_path):
    """WKB 캐시 기반 교차 검사"""
    if not wkb_path.exists():
        print(f"  {name}: WKB 없음 → 건너뜀")
        return set()
    zone_geom = shapely.from_wkb(wkb_path.read_bytes())
    zone_utm = gpd.GeoDataFrame(geometry=[zone_geom], crs=UTM)
    try:
        joined = gpd.sjoin(pts_gdf[["idx","geometry"]], zone_utm, how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  {name:15s}: WKB  → 제외 {len(excl)}개")
        return excl
    except Exception as e:
        print(f"  {name}: WKB 검사 오류: {e}")
        return set()


def json_to_gdf(json_path):
    """Overpass JSON → GeoDataFrame (WGS84)"""
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    elements = data.get("elements", [])
    nodes = {e["id"]: (e["lon"], e["lat"])
             for e in elements if e["type"] == "node" and "lat" in e}
    geoms = []
    for e in elements:
        if e["type"] == "node" and "lat" in e:
            geoms.append(Point(e["lon"], e["lat"]))
        elif e["type"] == "way" and "nodes" in e:
            coords = [nodes[n] for n in e["nodes"] if n in nodes]
            if len(coords) < 2: continue
            if coords[0] == coords[-1] and len(coords) >= 4:
                try:   geoms.append(Polygon(coords))
                except: geoms.append(LineString(coords))
            else:
                geoms.append(LineString(coords))
    if not geoms:
        return gpd.GeoDataFrame(geometry=[], crs=WGS84)
    return gpd.GeoDataFrame(geometry=geoms, crs=WGS84)


def check_zone_json(name, json_path, buffer_m):
    """JSON 파일 기반 교차 검사
    기준점 점들을 buffer_m 만큼 팽창시켜 zone 피처와 교차 여부 확인"""
    if not json_path.exists():
        print(f"  {name}: JSON 없음 → 건너뜀")
        return set()
    try:
        zone_gdf = json_to_gdf(json_path)
        if len(zone_gdf) == 0:
            print(f"  {name:15s}: 피처 없음 → 건너뜀")
            return set()
        zone_utm = zone_gdf.to_crs(UTM)
        # 기준점을 buffer_m 팽창 → zone 피처와 교차 여부
        pts_buf = pts_gdf[["idx","geometry"]].copy()
        pts_buf["geometry"] = pts_buf.geometry.buffer(buffer_m)
        joined = gpd.sjoin(pts_buf, zone_utm[["geometry"]], how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  {name:15s}: {len(zone_gdf):7,}개 피처 → 제외 {len(excl)}개")
        return excl
    except Exception as e:
        print(f"  {name}: 오류 → {e}")
        return set()


# ── WKB 캐시 존 ─────────────────────────────────────────────
excluded |= check_zone_wkb("power",   CACHE_DIR / "power.wkb")
excluded |= check_zone_wkb("anomaly", CACHE_DIR / "anomaly.wkb")

# ── JSON 존 ─────────────────────────────────────────────────
for zone_name, json_path in JSON_FILES.items():
    buf = BUFFERS[zone_name]
    excluded |= check_zone_json(zone_name, json_path, buf)

print(f"\n  ▶ 총 제외 대상: {len(excluded)} / {total_pm}개 기준점")
print(f"  ▶ 잔존 기준점 : {total_pm - len(excluded)}개")


# ════════════════════════════════════════════════════════════════
# 3. 제외 Placemark 제거 + 폴더명 업데이트
# ════════════════════════════════════════════════════════════════
print("\n[3] KML 수정 중...")

# 폴더별 pm 리스트 구성
from collections import defaultdict
folder_pms = defaultdict(list)  # id(folder_el) → [pm_record_idx, ...]
for i, rec in enumerate(pm_records):
    folder_pms[id(rec["folder_el"])].append(i)
print(f"  폴더 수: {len(folder_pms)}개 (각 폴더의 기준점 합: {sum(len(v) for v in folder_pms.values())}개)")

removed_total = 0
for folder_el in ref_folders:
    fid = id(folder_el)
    pm_idxs = folder_pms[fid]

    # 제거 대상
    to_remove = [i for i in pm_idxs if i in excluded]
    for i in to_remove:
        pm_el = pm_records[i]["pm_el"]
        try:
            folder_el.remove(pm_el)
            removed_total += 1
        except ValueError:
            pass

    # 폴더명 업데이트
    remaining = len(pm_idxs) - len(to_remove)
    nm_el = folder_el.find(ktag("name"))
    if nm_el is not None and nm_el.text:
        old_nm = nm_el.text
        # "📐 근접 기준점 (13개)" → 숫자만 교체
        new_nm = re.sub(r"\(\d+개\)", f"({remaining}개)", old_nm)
        nm_el.text = new_nm

print(f"  제거 완료: {removed_total}개")


# ════════════════════════════════════════════════════════════════
# 4. KML 저장
# ════════════════════════════════════════════════════════════════
print("\n[4] KML 저장...")

# ET.tostring → 선언 없이 출력되므로 직접 붙이기
kml_str = ET.tostring(root, encoding="unicode")
# 선언 추가
output = '<?xml version="1.0" encoding="UTF-8"?>\n' + kml_str

KML_OUT.write_text(output, encoding="utf-8")
sz = KML_OUT.stat().st_size / 1024
print(f"  ✅ {KML_OUT.name}")
print(f"     크기: {sz:.1f} KB")

# ════════════════════════════════════════════════════════════════
# 5. 요약
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 64)
print("  ✅ 완료!")
print()
print("  ── 기준점 필터링 결과 ───────────────────────────────")
print(f"  원본 기준점 : {total_pm}개")
print(f"  제외 (제외구역 내): {len(excluded)}개")
print(f"  최종 잔존 : {total_pm - len(excluded)}개")
print()
print("  ── 제외구역 기준 ───────────────────────────────────")
for nm, buf in BUFFERS.items():
    src = "WKB캐시" if nm in ("power","anomaly") else "JSON"
    print(f"    {nm:15s}  {buf:5,}m  [{src}]")
print("=" * 64)
