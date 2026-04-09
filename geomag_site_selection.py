#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
대한민국 지구자기장 모델 구축을 위한 측정 입지 선정 시스템
Korea Geomagnetic Field Model — Measurement Site Selection Tool

━━━ 법적 근거 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  국립기상과학원 지자기관측 업무 규정
  제13조(도상선점)
    ② 지자기점은 가능한 전국에 균일하게 분포되게 하고, 시가지·철도·
       송전탑 등으로부터 다음과 같이 충분한 거리를 확보하여야 한다.
       ┌──────────────────┬─────────────┐
       │  구  분           │ 거리(km)    │
       ├──────────────────┼─────────────┤
       │ 직류철도          │ 5.0 이상    │
       │ 교류 및 일반철도  │ 2.0 이상    │
       │ 고압철탑          │ 1.0 이상    │
       │ 송전탑            │ 0.5 이상    │
       └──────────────────┴─────────────┘
  제14조(선점실시)
    ① 인공적인 지자기 잡음이 발생하지 않고 자연자장을 측정할 수
       있는 조건을 갖추고 있는 지역을 선점한다.
    ② 장래에도 영구적으로 유지될 수 있는지를 고려한다.
    ③ 지하의 자기발생원이 있을 수 있으므로 선점지역 주변에서
       전자력을 측정하여 자장분포를 조사한다.
    ④ 태양 또는 북극성 관측에 의한 진북관측과 방위표 설치 등에
       지장이 없어야 한다.

━━━ 적용 버퍼 (OSM 데이터 기반 단순화) ━━━━━━━━━━━━━━━━━━━
  [1] 전력 인프라 (고압철탑·송전탑 통합)  반경  1.0 km
        ▸ OSM에서 고압·일반 철탑 구분 어려움 → 보수적으로 1.0 km 적용
  [2] 철도 (직류·교류·일반 통합)          반경  5.0 km
        ▸ 직류/교류 구분 어려움 → 직류 기준(최대값) 5.0 km 적용
  [3] 도시·주거 지역                       전체 제외
  [4] 파이프라인                           반경  0.5 km
  [5] 통신탑·기지국                        반경  0.5 km
  [6] 풍력발전기                           반경  0.5 km
  [7] 채석장·광산                          반경  1.0 km  (제14조③ 지하 자기발생원)

━━━ 지자기 이상 기준 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [8] 부지 내 자기장 변화폭 (KIGAM 1.5분 격자, 반경 0.05°≈5.5km, P90-P10)
        ▸ ≤ 100 nT : 우수 후보
        ▸ 100~200 nT: 현장 검토 필요
        ▸ > 200 nT : 제외 (hard cut)
        ▸ ※ 광역 자력이상도 예비선정이며, 최종 확정은 현장 정밀 자력측량 필요

━━━ 모델 구축 우선순위 가중치 ━━━━━━━━━━━━━━━━━━━━━━━━━━━
  - 후보점 공간 희소성 (전국 균일 분포, 제13조① 기준)
  - 지형 대표성, 전자기 이격도, 자기이상 균일도

사용법:
  python geomag_site_selection.py

출력:
  output/geomag_site_selection.html   - 대화형 지도 (OSM 기반)
  output/candidate_sites.csv          - 후보 지점 목록 (위경도·우선순위)
  data/                               - Overpass API 캐시 파일
"""

import os
import sys
import json
import time
import warnings
import requests
import numpy as np
import pandas as pd
import geopandas as gpd
import folium
from folium import plugins
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, shape
from shapely.ops import unary_union
from pathlib import Path

warnings.filterwarnings("ignore")

# ============================================================
# 설정 (Configuration)
# ============================================================

OUTPUT_DIR = Path("output")
DATA_DIR   = Path("data")
OUTPUT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

# 대한민국 경계 bbox (minlon, minlat, maxlon, maxlat)
KOREA_BBOX = (124.5, 33.0, 129.6, 38.9)

# ── 인공 간섭 제외 버퍼 거리 (제13조·제14조 기준) ───────────
#   제13조 원문:  직류철도 5 km / 교류·일반철도 2 km
#                고압철탑 1 km / 송전탑 0.5 km
#   OSM 단순화:  철도 종류 구분 불가 → 직류 기준 5 km (보수적)
#                철탑 종류 구분 불가 → 고압철탑 기준 1 km (보수적)
POWER_BUFFER_M      = 1_000   # [1] 고압철탑·송전탑 통합  1.0 km  (제13조)
RAILWAY_BUFFER_M    = 5_000   # [2] 철도 통합             5.0 km  (제13조 직류철도 기준)
# [3] 도시·주거       → 폴리곤 직접 제외             (제14조①)
PIPELINE_BUFFER_M   =   500   # [4] 파이프라인             0.5 km
# [5] 군사시설        → 폴리곤 직접 제외             (제14조①)
COMM_TOWER_BUFFER_M =   500   # [6] 통신탑·기지국          0.5 km
WIND_BUFFER_M       =   500   # [7] 풍력발전기              0.5 km
QUARRY_BUFFER_M     = 1_000   # [8] 채석장·광산            1.0 km  (제14조③)
# 도시·주거 지역 버퍼 (완화된 기준)
URBAN_DENSE_BUFFER_M    =   500   # 상업/공업/건설 지역 (핵심 도심)
URBAN_RESIDENT_BUFFER_M =   300   # 주거 지역 (농촌 친화적 완화)

# 후보 격자 간격 (m)
GRID_SPACING_M = 10_000   # 10 km

# 좌표 참조계: Korea 2000 / Korea Unified Coordinate System
# EPSG:5179 — 한반도 전역 분석 기준 (UTM 52N 대비 한국 전역 등거리 보장)
UTM_CRS   = "EPSG:5179"
WGS84_CRS = "EPSG:4326"

# Overpass API
OVERPASS_URL = "https://overpass-api.de/api/interpreter"

# 지자기 이상 기준 (KIGAM 1.5분 격자, 반경 0.05°≈5.5km 내 P90-P10)
# ■ 평가 스케일: 반경 0.05°(≈5.5 km) — KIGAM 1.5' 최소 실용 반경
# ■ 지표: P90-P10 (이상치 강건 범위, robust range)
#   ≤ 100 nT : 우수 (녹색)
#   100~200 nT: 현장 검토 필요 (주황)
#   > 200 nT : 제외 (red, hard cut)
# ※ KIGAM 1.5분(≈2.8 km) 데이터는 1 km 스케일을 직접 해상하지 못함.
#   광역 자력이상도 기반 예비선정이며, 최종 확정은 현장 정밀 자력측량 필요.
ANOMALY_EXCLUDE_THRESHOLD_NT = 200.0  # nT  hard cut (제외)
ANOMALY_CAUTION_THRESHOLD_NT = 100.0  # nT  현장검토 권고
ANOMALY_SITE_RADIUS_DEG      =  0.05  # °   ≈ 5.5 km
# 하위 호환성 유지
ANOMALY_VARIATION_THRESHOLD  = ANOMALY_EXCLUDE_THRESHOLD_NT

# ── 1:50,000 지형도 도엽 셰이프파일 (NGII 국가기본도) ────────────
# 출처: 국토지리정보원 국가기본도_도엽인덱스50K (TN_MAPINDX_50K.shp)
# CRS: Korea 2000 Korea Unified Coordinate System (EPSG:5186 유사)
# 주요 컬럼: MAPID_NM(도엽명), MAPID_NO(도엽번호, 예: NI52-7-01)
TOPO50K_SHP = DATA_DIR / "국가기본도_도엽인덱스50K" / "TN_MAPINDX_50K.shp"
NGII_SHEETS_CSV = DATA_DIR / "ngii_50k_sheets.csv"  # 백업용 CSV (SHP 없을 때)
KIGAM_MAG_DAT = DATA_DIR / "mag_1982-2018_1.5min_ed.dat"


# ============================================================
# Overpass API 유틸리티
# ============================================================

def _overpass_request(query: str, timeout: int = 300) -> dict | None:
    """단일 Overpass API POST 요청"""
    resp = requests.post(
        OVERPASS_URL,
        data={"data": query},
        timeout=timeout,
        headers={"Accept-Charset": "utf-8"},
    )
    resp.raise_for_status()
    return resp.json()


def query_overpass(query: str, cache_file: Path, timeout: int = 300) -> dict | None:
    """Overpass API 쿼리 (디스크 캐싱 지원, 3회 재시도)"""
    if cache_file.exists():
        print(f"    캐시 로드: {cache_file.name}")
        with open(cache_file, "r", encoding="utf-8") as f:
            return json.load(f)

    print(f"    Overpass API 쿼리 전송 중 ({cache_file.stem})...")
    for attempt in range(3):
        try:
            data = _overpass_request(query, timeout)
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False)
            print(f"    완료 — {len(data.get('elements', []))}개 요소")
            return data
        except Exception as exc:
            print(f"    시도 {attempt + 1}/3 실패: {exc}")
            if attempt < 2:
                wait = 60 * (attempt + 1)
                print(f"    {wait}초 대기 후 재시도...")
                time.sleep(wait)
    print("    ✗ 데이터 취득 실패")
    return None


def elements_to_gdf(
    data: dict | None,
    keep_tags: list | None = None,
) -> gpd.GeoDataFrame:
    """Overpass JSON 응답 → GeoDataFrame (WGS84)

    keep_tags: 보존할 OSM 태그 키 목록.
               None(기본)이면 geometry + osm_id 만 저장.
               태그를 모두 컬럼으로 펼치면 수백만 행에서 메모리 폭발이 발생하므로
               필요한 태그만 명시적으로 지정한다.
    """
    empty = gpd.GeoDataFrame(geometry=[], crs=WGS84_CRS)
    if not data or "elements" not in data:
        return empty

    elements = data["elements"]
    # node 좌표 인덱스 구축
    nodes = {
        e["id"]: (e["lon"], e["lat"])
        for e in elements
        if e["type"] == "node" and "lat" in e
    }

    geoms:    list = []
    osm_ids:  list = []
    tag_cols: dict = {k: [] for k in (keep_tags or [])}

    for elem in elements:
        tags = elem.get("tags", {})
        geom = None

        if elem["type"] == "node" and "lat" in elem:
            geom = Point(elem["lon"], elem["lat"])

        elif elem["type"] == "way" and "nodes" in elem:
            coords = [nodes[n] for n in elem["nodes"] if n in nodes]
            if len(coords) >= 2:
                if coords[0] == coords[-1] and len(coords) >= 4:
                    try:
                        geom = Polygon(coords)
                    except Exception:
                        geom = LineString(coords)
                else:
                    geom = LineString(coords)

        if geom is not None:
            geoms.append(geom)
            osm_ids.append(elem["id"])
            for k in tag_cols:
                tag_cols[k].append(tags.get(k))

    if not geoms:
        return empty

    gdf = gpd.GeoDataFrame(
        {"osm_id": osm_ids, **tag_cols},
        geometry=geoms,
        crs=WGS84_CRS,
    )
    return gdf


# ============================================================
# 데이터 취득
# ============================================================

def _get_naturalearth_kor_land():
    """
    naturalearth **50m** 해상도 한국(KOR) 육지 폴리곤을 GeoDataFrame으로 반환.
    110m 해상도는 제주도·서해 도서를 누락하므로 50m 사용.
    우선순위:
      1. data/naturalearth_kor_50m.geojson (로컬 캐시)
      2. GitHub nvkelso/natural-earth-vector 50m 다운로드 후 캐시
      3. 실패 시 None 반환
    """
    ne_cache = DATA_DIR / "naturalearth_kor_50m.geojson"
    # 110m 구버전 캐시 삭제 (해상도 부족)
    old_cache = DATA_DIR / "naturalearth_kor.geojson"
    if old_cache.exists():
        old_cache.unlink()
    if ne_cache.exists():
        try:
            gdf = gpd.read_file(str(ne_cache))
            if "iso_a3" not in gdf.columns:
                gdf["iso_a3"] = "KOR"
            return gdf
        except Exception:
            pass
    try:
        import io
        url = (
            "https://raw.githubusercontent.com/nvkelso/natural-earth-vector"
            "/master/geojson/ne_50m_admin_0_countries.geojson"
        )
        print("    naturalearth 50m 다운로드 중 (정밀 해안선)...")
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        world = gpd.read_file(io.BytesIO(resp.content))
        # 컬럼명이 버전마다 다를 수 있음
        iso_col = next(
            (c for c in world.columns if c.upper() in ("ISO_A3", "ADM0_ISO", "ISO_A3_EH")),
            None,
        )
        if iso_col:
            kor = world[world[iso_col] == "KOR"].copy()
        else:
            name_col = next((c for c in world.columns if "NAME" in c.upper()), None)
            kor = world[world[name_col].str.contains("Korea", na=False)].copy() if name_col else gpd.GeoDataFrame()
        if len(kor) > 0:
            kor[["geometry"]].to_file(str(ne_cache), driver="GeoJSON")
            kor["iso_a3"] = "KOR"
            bounds = unary_union(kor.geometry).bounds
            print(f"    naturalearth 50m KOR 저장: {ne_cache.name}")
            print(f"    범위: lon {bounds[0]:.2f}~{bounds[2]:.2f}, lat {bounds[1]:.2f}~{bounds[3]:.2f}")
            return kor[["iso_a3", "geometry"]]
    except Exception as exc:
        print(f"    ⚠ naturalearth 다운로드 실패 ({exc})")
    return None


def _clip_to_land(poly):
    """
    shapely polygon을 naturalearth 육지 폴리곤과 교차(intersection)하여
    해양·연안 바다 영역을 제거한다.

    OSM admin_level=2 경계는 영해(territorial sea)를 포함하는 경우가 있어
    격자점이 바다 위에 생성되는 문제를 이 함수로 해결한다.
    naturalearth 폴리곤을 0.05도(~5 km) 외측 버퍼하여 소규모 연안 섬이
    잘리지 않도록 보정한다.
    """
    try:
        kor_land = _get_naturalearth_kor_land()
        if kor_land is None or len(kor_land) == 0:
            return poly
        land_geom = unary_union(kor_land.geometry).buffer(0.05)  # ~5 km 외측 여유
        clipped = poly.intersection(land_geom)
        if clipped is None or clipped.is_empty:
            print("    ⚠ 육지 클리핑 결과 비어있음 — 원본 사용")
            return poly
        print(f"    육지 클리핑: {poly.area:.3f}° → {clipped.area:.3f}° (해양 제거)")
        return clipped
    except Exception as exc:
        print(f"    ⚠ 육지 클리핑 실패 ({exc}) — 원본 사용")
        return poly


def _build_korea_polygon_from_overpass(data: dict):
    """
    Overpass `out geom` 응답(relation)에서 outer way 선분을 모아
    shapely polygonize로 대한민국 육지 폴리곤을 재구성한다.
    섬을 포함한 MultiPolygon을 반환.
    """
    from shapely.ops import polygonize, linemerge

    if not data or not data.get("elements"):
        return None

    outer_lines = []
    for elem in data["elements"]:
        if elem["type"] != "relation":
            continue
        for member in elem.get("members", []):
            if member.get("role", "outer") in ("outer", "") and "geometry" in member:
                coords = [(p["lon"], p["lat"]) for p in member["geometry"]]
                if len(coords) >= 2:
                    outer_lines.append(LineString(coords))

    if not outer_lines:
        return None

    merged = linemerge(outer_lines)
    polys  = list(polygonize(merged))
    if not polys:
        return None

    return unary_union(polys)   # MultiPolygon (본토 + 섬 포함)


def get_korea_boundary() -> gpd.GeoDataFrame:
    """
    대한민국 육지 경계 (GeoDataFrame, WGS84).

    우선순위:
      1. data/korea_boundary.geojson 캐시 (정밀 경계)
      2. Overpass API — admin_level=2 relation의 outer way 재구성
      3. naturalearth_lowres fallback (저해상도, 해양 오염 가능)
    """
    print("\n[1/6] 대한민국 경계 로드...")

    # ── 1. 캐시 우선 로드 ────────────────────────────────────
    geojson_cache = DATA_DIR / "korea_boundary.geojson"
    if geojson_cache.exists():
        gdf = gpd.read_file(str(geojson_cache))
        print(f"    정밀 경계 캐시 로드: {len(gdf)} 피처")
        return gdf.to_crs(WGS84_CRS)

    # ── 2. Overpass에서 정밀 경계 취득 ───────────────────────
    print("    Overpass API에서 정밀 경계 취득 중 (최초 1회)...")
    json_cache = DATA_DIR / "korea_boundary_raw.json"
    query = """
    [out:json][timeout:120];
    relation["ISO3166-1"="KR"]["admin_level"="2"];
    out geom;
    """
    data = query_overpass(query, json_cache, timeout=120)
    korea_poly = _build_korea_polygon_from_overpass(data)

    if korea_poly is not None and not korea_poly.is_empty:
        # ── naturalearth 육지 폴리곤과 교차 → 해양 제거 ────────
        korea_poly = _clip_to_land(korea_poly)
        gdf = gpd.GeoDataFrame(geometry=[korea_poly], crs=WGS84_CRS)
        gdf.to_file(str(geojson_cache), driver="GeoJSON")
        print(f"    정밀 경계 취득 성공 (육지 클리핑 완료) → {geojson_cache.name} 저장")
        return gdf

    # ── 3. naturalearth fallback ──────────────────────────────
    print("    Overpass 실패, naturalearth fallback 사용")
    try:
        korea = _get_naturalearth_kor_land()
        if korea is not None and len(korea) > 0:
            return korea.to_crs(WGS84_CRS)
    except Exception:
        pass

    # ── 4. 최후 수단: 근사 폴리곤 ────────────────────────────
    approx = Polygon([
        (124.6, 34.0), (125.5, 33.2), (126.9, 33.1),
        (129.5, 34.6), (129.4, 35.6), (129.5, 37.0),
        (129.2, 37.7), (128.6, 38.3), (127.4, 38.6),
        (126.5, 38.3), (125.0, 38.0), (124.4, 37.0),
        (124.6, 36.0), (125.3, 34.7), (124.6, 34.0),
    ])
    print("    ⚠ 근사 경계 사용")
    return gpd.GeoDataFrame(geometry=[approx], crs=WGS84_CRS)


def get_power_infrastructure() -> gpd.GeoDataFrame:
    """
    송전탑 (power=tower/pole) 및 송전선 (power=line/cable) 취득
    OSM 대한민국 전역 쿼리
    """
    print("\n[2/6] 송전 인프라 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      node["power"="tower"]({s},{w},{n},{e});
      node["power"="pole"]({s},{w},{n},{e});
      way["power"="line"]({s},{w},{n},{e});
      way["power"="cable"]({s},{w},{n},{e});
      way["power"="minor_line"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "power_infra.json"),
        keep_tags=["power"],
    )
    print(f"    송전 인프라: {len(gdf)}개 요소")
    return gdf


def get_railways() -> gpd.GeoDataFrame:
    """
    철도 (rail, subway, light_rail, narrow_gauge, tram) 취득
    """
    print("\n[3/6] 철도 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      way["railway"~"^(rail|subway|light_rail|narrow_gauge|tram)$"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "railways.json"),
        keep_tags=["railway"],
    )
    print(f"    철도: {len(gdf)}개 요소")
    return gdf


def get_urban_dense() -> gpd.GeoDataFrame:
    """
    [3-A] 핵심 도심·산업 지역 취득 (완전 제외 대상)
    상업/공업/유통/건설 지역 — 강한 전자기 잡음원
    """
    print("\n[4a] 핵심 도심·산업 지역 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      way["landuse"~"^(commercial|industrial|retail|construction|garages)$"]({s},{w},{n},{e});
      relation["landuse"~"^(commercial|industrial|retail)$"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "urban_dense.json"),
        keep_tags=["landuse"],
    )
    print(f"    핵심도심/산업: {len(gdf)}개 요소")
    return gdf


def get_urban_residential() -> gpd.GeoDataFrame:
    """
    [3-B] 주거 지역 취득 (완화 버퍼 적용 대상)
    주거/마을 지역 — 농촌 지역도 포함되므로 작은 버퍼만 적용
    """
    print("\n[4b] 주거·취락 지역 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      way["landuse"="residential"]({s},{w},{n},{e});
      relation["landuse"="residential"]({s},{w},{n},{e});
      way["place"~"^(city|town|village|suburb|quarter)$"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "urban_residential_v2.json"),
        keep_tags=["landuse", "place"],
    )
    print(f"    주거/취락: {len(gdf)}개 요소")
    return gdf


def get_pipelines() -> gpd.GeoDataFrame:
    """
    [4] 파이프라인 취득 (man_made=pipeline)
    금속 파이프는 잔류 자화로 수십~수백 nT 이상을 일으킬 수 있음.
    """
    print("\n[추가-1] 파이프라인 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:180];
    (
      way["man_made"="pipeline"]({s},{w},{n},{e});
      node["man_made"="pipeline"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "pipelines.json"),
        keep_tags=["man_made", "substance"],
    )
    print(f"    파이프라인: {len(gdf)}개 요소")
    return gdf


def get_military_areas() -> gpd.GeoDataFrame:
    """
    [5] 군사시설 취득 (landuse=military, military=*)
    레이더·통신·무기 시스템에 의한 강한 전자기 간섭.
    """
    print("\n[추가-2] 군사시설 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    # way 노드 좌표 재구성(>;) 없이 out geom 으로 직접 수신 → 타임아웃 방지
    query = f"""
    [out:json][timeout:120];
    (
      way["landuse"="military"]({s},{w},{n},{e});
      relation["landuse"="military"]({s},{w},{n},{e});
    );
    out geom qt;
    """
    raw = query_overpass(query, DATA_DIR / "military.json", timeout=120)

    # out geom 포맷: way에 geometry 배열이 직접 포함
    geoms, ids = [], []
    if raw and raw.get("elements"):
        for elem in raw["elements"]:
            if elem["type"] == "way" and "geometry" in elem:
                coords = [(p["lon"], p["lat"]) for p in elem["geometry"]]
                if len(coords) >= 4 and coords[0] == coords[-1]:
                    try:
                        geoms.append(Polygon(coords))
                        ids.append(elem["id"])
                    except Exception:
                        pass
    if geoms:
        gdf = gpd.GeoDataFrame({"osm_id": ids, "landuse": "military"},
                               geometry=geoms, crs=WGS84_CRS)
    else:
        gdf = gpd.GeoDataFrame(geometry=[], crs=WGS84_CRS)
    print(f"    군사시설: {len(gdf)}개 요소")
    return gdf


def get_comm_towers() -> gpd.GeoDataFrame:
    """
    [6] 통신탑·기지국 취득 (man_made=mast/tower, communication=*)
    전자기파 방사로 수 nT ~ 수십 nT 간섭 가능.
    """
    print("\n[추가-3] 통신탑·기지국 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:180];
    (
      node["man_made"="mast"]({s},{w},{n},{e});
      node["man_made"="tower"]["tower:type"~"communication|radio|television"]({s},{w},{n},{e});
      way["man_made"="mast"]({s},{w},{n},{e});
      node["communication"~"mobile_phone|television|radio|microwave"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "comm_towers.json"),
        keep_tags=["man_made", "tower:type", "communication"],
    )
    print(f"    통신탑/기지국: {len(gdf)}개 요소")
    return gdf


def get_wind_turbines() -> gpd.GeoDataFrame:
    """
    [7] 풍력발전기 취득 (generator:source=wind)
    회전하는 대형 금속 날개·발전기 자장 (수 nT ~ 수십 nT).
    """
    print("\n[추가-4] 풍력발전기 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:180];
    (
      node["generator:source"="wind"]({s},{w},{n},{e});
      way["generator:source"="wind"]({s},{w},{n},{e});
      node["power"="generator"]["generator:source"="wind"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "wind_turbines.json"),
        keep_tags=["power", "generator:source"],
    )
    print(f"    풍력발전기: {len(gdf)}개 요소")
    return gdf


def get_quarries_mines() -> gpd.GeoDataFrame:
    """
    [8] 채석장·광산 취득 (landuse=quarry, man_made=mineshaft/adit)
    지표·지하 광물(철광석, 자철광 등)에 의한 국소 자기 이상.
    """
    print("\n[추가-5] 채석장·광산 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:180];
    (
      way["landuse"="quarry"]({s},{w},{n},{e});
      relation["landuse"="quarry"]({s},{w},{n},{e});
      node["man_made"~"mineshaft|adit"]({s},{w},{n},{e});
      way["man_made"~"mineshaft|adit"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "quarries_mines.json"),
        keep_tags=["landuse", "man_made"],
    )
    print(f"    채석장/광산: {len(gdf)}개 요소")
    return gdf


def get_water_bodies() -> gpd.GeoDataFrame:
    """
    [수계] 대형 수계·수면 취득 (자력계 설치 불가 지역)
    호수, 저수지, 강 수면 등 — 관측소 부지로 부적합
    """
    print("\n[추가-6] 수계·수면 데이터 취득...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:180];
    (
      way["natural"="water"]({s},{w},{n},{e});
      relation["natural"="water"]({s},{w},{n},{e});
      way["water"~"^(lake|reservoir|pond|oxbow)$"]({s},{w},{n},{e});
      way["waterway"~"^(riverbank|dock)$"]({s},{w},{n},{e});
      relation["waterway"="riverbank"]({s},{w},{n},{e});
      way["landuse"="reservoir"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "water_bodies.json"),
        keep_tags=["natural", "water", "waterway", "landuse"],
    )
    print(f"    수계/수면: {len(gdf)}개 요소")
    return gdf


def get_public_land() -> gpd.GeoDataFrame:
    """
    [⑦] 국공유지·공공 관리 토지 취득 (부지 지속성 산정용)
    국립공원, 자연보호구역, 국유림, 정부 소유 토지 등
    측정 부지로서 장기 안정적 유지가 가능한 공공토지 우선
    """
    print("\n국공유지 데이터 취득 (부지 지속성용)...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      way["boundary"~"^(national_park|protected_area)$"]({s},{w},{n},{e});
      relation["boundary"~"^(national_park|protected_area)$"]({s},{w},{n},{e});
      way["leisure"~"^(nature_reserve|park)$"]({s},{w},{n},{e});
      way["landuse"~"^(forest|recreation_ground|grass|meadow)$"]
         ["access"~"^(public|yes)$"]({s},{w},{n},{e});
      way["ownership"~"^(government|public|national)$"]({s},{w},{n},{e});
      way["operator:type"~"^(government|public)$"]({s},{w},{n},{e});
      way["landuse"="forest"]["operator"~"산림청|국유"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "public_land.json"),
        keep_tags=["boundary", "leisure", "landuse", "ownership", "operator"],
    )
    print(f"    국공유지: {len(gdf)}개 요소")
    return gdf


def get_access_roads() -> gpd.GeoDataFrame:
    """
    [⑧] 일반 도로 취득 (관리 접근성 산정용, 고속도로 제외)
    국도·지방도·군도 등 — 반경 1 km 이내 접근 가능 여부 확인
    고속국도(motorway·trunk)는 접근이 어려우므로 제외
    """
    print("\n일반 도로 데이터 취득 (접근성용)...")
    s, w, n, e = KOREA_BBOX[1], KOREA_BBOX[0], KOREA_BBOX[3], KOREA_BBOX[2]
    query = f"""
    [out:json][timeout:300];
    (
      way["highway"~"^(primary|secondary|tertiary|unclassified|residential|service)$"]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;
    """
    gdf = elements_to_gdf(
        query_overpass(query, DATA_DIR / "access_roads.json"),
        keep_tags=["highway"],
    )
    print(f"    일반 도로: {len(gdf)}개 요소")
    return gdf


# ============================================================
# 자기이상도 (Magnetic Anomaly) — EMAG2 선택적 처리
# ============================================================

def load_emag2_korea(emag2_path: str | Path | None = None) -> np.ndarray | None:
    """
    EMAG2 데이터 로드 및 한반도 영역 추출.

    emag2_path: EMAG2 ASCII grid 파일 경로
                (https://www.ngdc.noaa.gov/geomag/emag2.shtml 에서 다운로드)
                컬럼: lon, lat, anomaly_nT
                예: EMAG2_V3_20170530.csv

    반환: (N,3) ndarray — [lon, lat, anomaly_nT]
          파일이 없으면 None 반환
    """
    if emag2_path is None:
        # 기본 위치 탐색
        candidates = [
            DATA_DIR / "EMAG2_V3_20170530.csv",
            DATA_DIR / "emag2.csv",
            DATA_DIR / "emag2.xyz",
        ]
        for p in candidates:
            if p.exists():
                emag2_path = p
                break

    if emag2_path is None or not Path(emag2_path).exists():
        return None

    print(f"    EMAG2 로드: {emag2_path}")

    # ── EMAG2 파일 포맷 자동 감지 ────────────────────────────
    # EMAG2 V3 (EMAG2_V3_20170530.csv) — 8컬럼 CSV:
    #   col0: row_idx  col1: col_idx  col2: lon  col3: lat
    #   col4: ?        col5: anomaly_nT  col6: data_type  col7: uncertainty
    # EMAG2 V2 / 단순 xyz — 3컬럼 (공백 구분):
    #   col0: lon  col1: lat  col2: anomaly_nT
    # ─────────────────────────────────────────────────────────
    with open(emag2_path, "r") as _f:
        # 헤더·주석 건너뛰고 첫 데이터 라인 확인
        first_line = ""
        for line in _f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                first_line = stripped
                break

    ncols = len(first_line.split(",")) if "," in first_line else len(first_line.split())

    if ncols >= 6:
        # V3: comma-separated 8컬럼
        print(f"    EMAG2 포맷: V3 ({ncols}컬럼 CSV) — col2=lon, col3=lat, col5=anomaly")
        df = pd.read_csv(
            emag2_path,
            comment="#",
            header=None,
            sep=",",
            usecols=[2, 3, 5],
            names=["lon", "lat", "anomaly"],
            na_values=["99999", "-99999", "999999", "-999999"],
            on_bad_lines="skip",
        )
    else:
        # V2 / XYZ: 공백 구분 3컬럼
        print(f"    EMAG2 포맷: V2/XYZ ({ncols}컬럼) — col0=lon, col1=lat, col2=anomaly")
        df = pd.read_csv(
            emag2_path,
            comment="#",
            header=None,
            sep=r"\s+",
            usecols=[0, 1, 2],
            names=["lon", "lat", "anomaly"],
            na_values=["99999", "-99999"],
            on_bad_lines="skip",
        )

    # 숫자 변환 + 결측 제거
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    df = df.astype(np.float32)

    # 한반도 영역 필터
    mask = (
        (df["lon"] >= KOREA_BBOX[0]) & (df["lon"] <= KOREA_BBOX[2]) &
        (df["lat"] >= KOREA_BBOX[1]) & (df["lat"] <= KOREA_BBOX[3])
    )
    sub = df[mask].reset_index(drop=True)
    sub.columns = ["lon", "lat", "anomaly_nT"]
    print(f"    한반도 영역 EMAG2 포인트: {len(sub)}개")
    return sub


def load_kigam_anomaly(path: Path | None = None) -> "pd.DataFrame | None":
    """
    한국지질자원연구원(KIGAM) 수치 자력 이상도 로드.
    파일: mag_1982-2018_1.5min_ed.dat
    포맷: ASCII, 9줄 헤더 건너뜀, 공백 구분 3컬럼
          Col0: Longitude(deg), Col1: Latitude(deg), Col2: MagAnomaly(nT)
    해상도: 1.5분 간격, WGS84 기준, 한반도 전역 18,593개 포인트
    출처: 박영수 외, 2019, 한국의 자력 이상도, 지구물리와 물리탐사, 22(1), p.29-36
    """
    if path is None:
        path = KIGAM_MAG_DAT
    if not Path(path).exists():
        return None

    print(f"    KIGAM 자력이상도 로드: {path.name}")
    try:
        df = pd.read_csv(
            path, sep=r"\s+", skiprows=9, header=None,
            names=["lon", "lat", "anomaly_nT"],
            na_values=["99999", "-99999"],
            on_bad_lines="skip",
        )
        df = df.apply(pd.to_numeric, errors="coerce").dropna()
        # 한반도 영역 필터
        mask = (
            (df["lon"] >= KOREA_BBOX[0]) & (df["lon"] <= KOREA_BBOX[2]) &
            (df["lat"] >= KOREA_BBOX[1]) & (df["lat"] <= KOREA_BBOX[3])
        )
        sub = df[mask].reset_index(drop=True)
        print(f"    KIGAM 포인트: {len(sub)}개 (1.5분 격자, {sub['anomaly_nT'].min():.0f}~{sub['anomaly_nT'].max():.0f} nT)")
        return sub
    except Exception as exc:
        print(f"    ⚠ KIGAM 로드 실패: {exc}")
        return None


def compute_anomaly_variation_zones(
    mag_data,
    exclude_threshold_nT: float = ANOMALY_EXCLUDE_THRESHOLD_NT,
    caution_threshold_nT: float = ANOMALY_CAUTION_THRESHOLD_NT,
    site_radius_deg:      float = ANOMALY_SITE_RADIUS_DEG,
    threshold_nT:         float | None = None,  # 하위 호환성
) -> gpd.GeoDataFrame | None:
    """
    부지 내 자기장 변화폭 기반 제외 구역 생성.

    ■ 평가 기준 (제14조③) — 단일 공간 스케일 통일:
      반경 0.05°(≈5.5 km) 내 P90-P10 (robust range, 이상치 강건)
      ≤ 100 nT   : 우수
      100~200 nT : 현장 검토 필요
      > 200 nT   : 제외 (hard cut) → 이 함수에서 제외 구역 생성

    ■ 데이터 한계 고지:
      KIGAM 1.5분(≈2.8 km) 데이터는 1 km 스케일을 직접 해상하지 못함.
      광역 자력이상도 기반 예비선정이며, 최종 확정은 현장 정밀 자력측량 필요.

    Parameters
    ----------
    mag_data             : KIGAM/EMAG2 DataFrame (lon, lat, anomaly_nT 컬럼)
    exclude_threshold_nT : 제외 임계값 (기본 200 nT, hard cut)
    caution_threshold_nT : 현장검토 임계값 (기본 100 nT, 점수 감점용)
    site_radius_deg      : 탐색 반경 (°, 기본 0.05° ≈ 5.5 km)

    Returns
    -------
    GeoDataFrame (WGS84) with column 'variability_nT' — 제외 구역(>200 nT)만 포함
    """
    # 하위 호환성: 구 threshold_nT 인수 → exclude_threshold_nT
    if threshold_nT is not None:
        exclude_threshold_nT = threshold_nT

    if mag_data is None or len(mag_data) < 4:
        return None

    from shapely.geometry import box as _box

    print(f"    자기이상 변화폭(P90-P10) 계산 중 (반경 {site_radius_deg}° ≈ "
          f"{site_radius_deg * 111:.0f} km, 제외 기준 {exclude_threshold_nT} nT)...")

    if hasattr(mag_data, "columns"):
        lons = mag_data["lon"].values.astype(float)
        lats = mag_data["lat"].values.astype(float)
        vals = mag_data["anomaly_nT"].values.astype(float)
    else:
        lons = mag_data[:, 0].astype(float)
        lats = mag_data[:, 1].astype(float)
        vals = mag_data[:, 2].astype(float)

    res  = site_radius_deg
    glon = np.arange(KOREA_BBOX[0], KOREA_BBOX[2] + res, res)
    glat = np.arange(KOREA_BBOX[1], KOREA_BBOX[3] + res, res)

    exclude_cells  = []
    variabilities  = []
    n_valid        = 0
    n_caution      = 0

    for lat_c in glat:
        for lon_c in glon:
            mask = (
                (np.abs(lats - lat_c) <= site_radius_deg) &
                (np.abs(lons - lon_c) <= site_radius_deg)
            )
            pts = vals[mask]
            if len(pts) < 3:
                continue
            n_valid += 1
            # P90-P10: robust range (이상치에 강건)
            variability = float(np.percentile(pts, 90) - np.percentile(pts, 10))
            if variability > caution_threshold_nT:
                n_caution += 1
            if variability > exclude_threshold_nT:
                exclude_cells.append(_box(
                    lon_c - res / 2, lat_c - res / 2,
                    lon_c + res / 2, lat_c + res / 2,
                ))
                variabilities.append(round(variability, 1))

    n_caution_only = n_caution - len(exclude_cells)
    pct_excl = len(exclude_cells) / n_valid * 100 if n_valid else 0
    print(f"    자기이상 평가: {n_valid}개 격자 평가")
    print(f"      ≤ {caution_threshold_nT} nT (우수): {n_valid - n_caution}개")
    print(f"      {caution_threshold_nT}~{exclude_threshold_nT} nT (현장검토): {n_caution_only}개")
    print(f"      > {exclude_threshold_nT} nT (제외): {len(exclude_cells)}개 ({pct_excl:.1f}%)")

    if not exclude_cells:
        print(f"    ℹ 제외 구역 없음 (P90-P10 기준 {exclude_threshold_nT} nT 초과 없음)")
        return None

    gdf = gpd.GeoDataFrame(
        {"variability_nT": variabilities},
        geometry=exclude_cells,
        crs=WGS84_CRS,
    )
    print(f"    ✅ 자기이상 제외 구역 생성: {len(exclude_cells)}개 셀 "
          f"(P90-P10 > {exclude_threshold_nT} nT)")
    print(f"    ⚠  광역 자력이상도(KIGAM 1.5분) 예비선정 결과 — 최종 확정은 현장 정밀 자력측량 필요")
    return gdf


# ============================================================
# 1:50,000 지형도 도엽 격자
# ============================================================

def create_topo_sheet_grid(korea_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    1:50,000 지형도 도엽 격자 생성.
    우선순위:
      1) data/국가기본도_도엽인덱스50K/TN_MAPINDX_50K.shp (NGII 셰이프파일, 248개)
      2) data/ngii_50k_sheets.csv (CSV 백업)
    반환: WGS84 GeoDataFrame, 컬럼: sheet_name(도엽명), sheet_code(도엽번호)
    """
    from shapely.geometry import box

    korea_wgs  = korea_gdf.to_crs(WGS84_CRS)
    korea_geom = unary_union(korea_wgs.geometry)

    # ── ① 셰이프파일 우선 로드 ───────────────────────────────
    if TOPO50K_SHP.exists():
        try:
            raw = gpd.read_file(str(TOPO50K_SHP), encoding="cp949")

            # ── CRS 처리 (네이티브 PRJ WKT 우선) ─────────────────────────
            # 문제 분석:
            #   pyproj의 EPSG:5179 정의가 버전에 따라 TOWGS84 파라미터를 포함하면
            #   WGS84 변환 시 수십~수백m 오프셋이 발생하여 도엽 위치가 어긋날 수 있음.
            # 해결책:
            #   SHP PRJ 파일의 네이티브 WKT("Korea_2000_Korea_Unified_CS")를 직접 사용.
            #   pyproj가 CRS를 읽어오지 못할 때만 EPSG:5179 폴백.
            #
            # ⚠ 만약 변환 후에도 위치가 맞지 않으면 _check_topo_alignment.py 실행하여
            #   실제 SHP 좌표 범위와 올바른 EPSG 코드를 확인할 것.
            native_crs = raw.crs
            if native_crs is None:
                print("    ⚠ PRJ CRS 없음 → EPSG:5179 폴백")
                raw = raw.set_crs("EPSG:5179", allow_override=True)
            else:
                print(f"    CRS: 네이티브 PRJ 사용 ({native_crs.name})")
                # EPSG:5179 override 하지 않음 — 네이티브 WKT로 정확 변환
            topo_raw = raw.to_crs(WGS84_CRS)

            # ★ TM 투영 곡선 보정: WGS84 변환 후 bbox 직사각형으로 교체
            # 지리 좌표계에서 도엽은 경위선 직교 직사각형이어야 함.
            # 변환 잔류 곡선을 제거해 NGII 도엽 경계와 정렬.
            topo_raw["geometry"] = topo_raw.geometry.apply(
                lambda g: box(*g.bounds)
            )
            topo = topo_raw

            # 컬럼 정리
            name_col   = "MAPID_NM"   if "MAPID_NM"   in topo.columns else topo.columns[1]
            code_col   = "MAPID_NO"   if "MAPID_NO"   in topo.columns else topo.columns[2]
            mapidcd_col = "MAPIDCD_NO" if "MAPIDCD_NO" in topo.columns else None

            # 남한 범위 필터 (WGS84 변환 후 위도 33~38.65°)
            centroids_y = topo.geometry.centroid.y
            topo = topo[centroids_y.between(33.0, 38.65)].copy()

            # 도엽명 정리:
            # MAPID_NM이 순수 숫자(코드)인 경우 MAPID_NO 코드로 대체
            def _clean_name(row):
                nm = str(row[name_col]).strip()
                no = str(row[code_col]).strip()
                if not nm or nm.isdigit() or nm in ("nan", "None"):
                    return no
                return nm

            topo["sheet_name"] = topo.apply(_clean_name, axis=1)
            topo["sheet_code"] = topo[code_col].astype(str)
            # MAPIDCD_NO: NGII 5자리 도엽번호 (예: 36807 = 안동)
            if mapidcd_col:
                topo["sheet_mapidcd"] = topo[mapidcd_col].astype(str)
            else:
                topo["sheet_mapidcd"] = ""

            # 한국 경계와 교차하는 도엽만
            topo = topo[topo.geometry.intersects(korea_geom)].reset_index(drop=True)
            keep_cols = ["sheet_name", "sheet_code", "sheet_mapidcd", "geometry"]
            result = topo[keep_cols].copy()
            sample_names = result["sheet_name"].head(5).tolist()
            print(f"  1:50,000 도엽 격자: {len(result)}개 셀 생성 (SHP, bbox 직사각형) "
                  f"— 도엽명 예: {sample_names}")
            return result
        except Exception as exc:
            print(f"  ⚠ 셰이프파일 로드 실패: {exc} — CSV 백업으로 전환")

    # ── ② CSV 백업 ────────────────────────────────────────────
    if not NGII_SHEETS_CSV.exists():
        raise FileNotFoundError(
            f"도엽 데이터 없음. 다음 중 하나를 준비하세요:\n"
            f"  {TOPO50K_SHP}\n  {NGII_SHEETS_CSV}"
        )

    df = pd.read_csv(NGII_SHEETS_CSV, encoding="utf-8-sig", dtype=str)
    for col in ["min_lon", "min_lat", "max_lon", "max_lat"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["min_lon", "min_lat", "max_lon", "max_lat"])
    df = df[
        (df["min_lat"] >= 33.0) & (df["max_lat"] <= 39.5) &
        (df["min_lon"] >= 124.0) & (df["max_lon"] <= 131.0)
    ].drop_duplicates(subset=["code"]).reset_index(drop=True)

    cells, names, codes = [], [], []
    for _, row in df.iterrows():
        cell = box(row["min_lon"], row["min_lat"], row["max_lon"], row["max_lat"])
        if korea_geom.intersects(cell):
            cells.append(cell)
            names.append(str(row["name"]))
            codes.append(str(row["code"]))

    result = gpd.GeoDataFrame(
        {"sheet_name": names, "sheet_code": codes, "geometry": cells},
        crs=WGS84_CRS,
    )
    print(f"  1:50,000 도엽 격자: {len(result)}개 셀 생성 (CSV 백업)")
    return result


# ============================================================
# 공간 분석
# ============================================================

def create_candidate_grid(korea_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """대한민국 영역 내 정규 격자점 생성 (UTM, GRID_SPACING_M 간격)"""
    print("\n[5/6] 후보 격자점 생성...")

    korea_utm = korea_gdf.to_crs(UTM_CRS)
    bounds = korea_utm.total_bounds  # (minx, miny, maxx, maxy)

    xs = np.arange(bounds[0], bounds[2], GRID_SPACING_M)
    ys = np.arange(bounds[1], bounds[3], GRID_SPACING_M)
    pts = [Point(x, y) for x in xs for y in ys]

    grid = gpd.GeoDataFrame(geometry=pts, crs=UTM_CRS)

    # 남한 경계 내부만
    korea_geom = unary_union(korea_utm.geometry)
    mask = grid.geometry.within(korea_geom)
    grid = grid[mask].reset_index(drop=True)
    print(f"    격자 범위: {len(xs)} × {len(ys)} = {len(xs)*len(ys)} 셀 → 경계내 {len(grid)}개")

    # ── 육지 마스크 (해양 격자점 이중 제거) ─────────────────────
    # naturalearth 육지 폴리곤으로 최종 확인 — 경계 재구성 오류 보완
    try:
        kor_land = _get_naturalearth_kor_land()
        if kor_land is not None and len(kor_land) > 0:
            land_utm = kor_land.to_crs(UTM_CRS)
            land_geom = unary_union(land_utm.geometry).buffer(2_000)  # 2 km 외측 여유 (도서 포함, 해양 오염 최소화)
            sea_mask = grid.geometry.within(land_geom)
            removed = (~sea_mask).sum()
            if removed > 0:
                grid = grid[sea_mask].reset_index(drop=True)
                print(f"    육지 마스크 적용: 해양 {removed}개 제거 → 잔여 {len(grid)}개")
    except Exception as exc:
        print(f"    ⚠ 육지 마스크 실패 ({exc}) — 건너뜀")

    return grid


def _make_valid_geom(g):
    """비유효 geometry를 shapely.make_valid 또는 buffer(0)으로 복구"""
    try:
        import shapely
        if hasattr(shapely, "make_valid"):
            v = shapely.make_valid(g)
        else:
            v = g.buffer(0)
        return v if v is not None and not v.is_empty else None
    except Exception:
        return None


def _chunked_union(geoms: list, chunk_size: int = 2_000):
    """
    대량의 shapely geometry 리스트를 청크 단위로 나누어 계층적으로 union.
    각 청크 처리 전에 비유효 geometry를 make_valid로 복구하여
    TopologyException을 방지한다.
    """
    # 1단계: 유효화
    print(f"      geometry 유효화 중 ({len(geoms)}개)...")
    valid = []
    for g in geoms:
        if g is None or g.is_empty:
            continue
        if not g.is_valid:
            g = _make_valid_geom(g)
        if g is not None and not g.is_empty:
            valid.append(g)
    print(f"      유효 geometry: {len(valid)}개 (제외 {len(geoms)-len(valid)}개)")
    geoms = valid

    # 2단계: 계층적 청크 유니온
    while len(geoms) > 1:
        next_level = []
        for i in range(0, len(geoms), chunk_size):
            chunk = geoms[i : i + chunk_size]
            try:
                merged = unary_union(chunk)
            except Exception:
                # 청크 내 문제 geometry 개별 처리 후 재시도
                fixed = [_make_valid_geom(g) or g for g in chunk]
                try:
                    merged = unary_union([g for g in fixed if g and not g.is_empty])
                except Exception as exc:
                    print(f"      ⚠ 청크 유니온 실패, 건너뜀: {exc}")
                    merged = None
            if merged is not None and not merged.is_empty:
                next_level.append(merged)
        geoms = next_level
        if len(geoms) > 1:
            print(f"      유니온 중간 결과: {len(geoms)}개 → 다음 단계...")
    return geoms[0] if geoms else None


def _build_zone(gdf: gpd.GeoDataFrame, buffer_m: float | None, label: str):
    """
    GeoDataFrame → UTM 유니온 폴리곤.
    buffer_m=None 이면 폴리곤 직접 사용 (군사시설·도시처럼 제외 구역이 이미 면적인 경우).
    buffer_m > 0 이면 선·점 geometry를 버퍼링 후 유니온.
    """
    if len(gdf) == 0:
        print(f"    {label}: 데이터 없음 (건너뜀)")
        return None

    utm = gdf.to_crs(UTM_CRS)

    if buffer_m is not None and buffer_m > 0:
        geoms = list(utm.geometry.buffer(buffer_m))
    else:
        # 폴리곤 직접 사용
        polys  = utm[utm.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].geometry.tolist()
        pts    = utm[utm.geometry.geom_type == "Point"].geometry.tolist()
        geoms  = polys + [p.buffer(1_500) for p in pts]  # 포인트는 1.5km 버퍼

    if not geoms:
        return None

    if len(geoms) > 5_000:
        print(f"    {label}: {len(geoms)}개 유니온 중...")
        result = _chunked_union(geoms, chunk_size=2_000)
    else:
        try:
            result = unary_union(geoms)
        except Exception:
            fixed = [_make_valid_geom(g) or g for g in geoms]
            result = unary_union([g for g in fixed if g and not g.is_empty])

    if result is None or result.is_empty:
        return None

    km2 = result.area / 1e6
    buf_str = f"  (버퍼 {buffer_m/1000:.1f} km)" if buffer_m else ""
    print(f"    {label}: {km2:.0f} km²{buf_str}")
    return result


def build_exclusion_zones(
    power_gdf:        gpd.GeoDataFrame,
    railway_gdf:      gpd.GeoDataFrame,
    urban_dense_gdf:  gpd.GeoDataFrame,
    urban_resid_gdf:  gpd.GeoDataFrame,
    pipeline_gdf:     gpd.GeoDataFrame,
    comm_gdf:         gpd.GeoDataFrame,
    wind_gdf:         gpd.GeoDataFrame,
    quarry_gdf:       gpd.GeoDataFrame,
    anomaly_gdf:      gpd.GeoDataFrame | None,
    water_gdf:        gpd.GeoDataFrame | None = None,
) -> dict:
    """
    각 제외 조건별 단일 유니온 폴리곤 반환 (UTM CRS, EPSG:5179)
    키: power, railway, urban_dense, urban_resid, pipeline, comm, wind, quarry,
        water, anomaly
    """
    print("\n제외 구역 구축 중...")
    zones = {}

    zones["power"]    = _build_zone(power_gdf,    POWER_BUFFER_M,
                                     "[1] 송전 인프라")
    zones["railway"]  = _build_zone(
        railway_gdf[railway_gdf.geometry.geom_type.isin(
            ["LineString", "MultiLineString"])
        ] if len(railway_gdf) > 0 else railway_gdf,
        RAILWAY_BUFFER_M, "[2] 철도")
    zones["urban_dense"] = _build_zone(urban_dense_gdf, URBAN_DENSE_BUFFER_M,
                                        "[3a] 핵심도심·산업")
    zones["urban_resid"] = _build_zone(urban_resid_gdf, URBAN_RESIDENT_BUFFER_M,
                                        "[3b] 주거·취락")
    zones["pipeline"] = _build_zone(pipeline_gdf, PIPELINE_BUFFER_M,
                                     "[4] 파이프라인")
    zones["comm"]     = _build_zone(comm_gdf,     COMM_TOWER_BUFFER_M,
                                     "[6] 통신탑/기지국")
    zones["wind"]     = _build_zone(wind_gdf,     WIND_BUFFER_M,
                                     "[7] 풍력발전기")
    zones["quarry"]   = _build_zone(quarry_gdf,   QUARRY_BUFFER_M,
                                     "[8] 채석장/광산")

    # 수계 제외 (폴리곤 직접 사용 — 수면 위에 관측소 설치 불가)
    if water_gdf is not None and len(water_gdf) > 0:
        zones["water"] = _build_zone(water_gdf, None, "[수계] 호수·저수지·강수면")
    else:
        zones["water"] = None
        print("    [수계] 수계 데이터 없음 (건너뜀)")

    if anomaly_gdf is not None and len(anomaly_gdf) > 0:
        an_utm = anomaly_gdf.to_crs(UTM_CRS)
        zones["anomaly"] = unary_union(an_utm.geometry)
        print(f"    [9] 자기이상도: {zones['anomaly'].area/1e6:.0f} km²  "
              f"(P90-P10 >{ANOMALY_EXCLUDE_THRESHOLD_NT} nT 제외)")
    else:
        zones["anomaly"] = None
        print("    [9] 자기이상도: KIGAM/EMAG2 파일 배치 시 활성화")

    return zones


def filter_candidates(
    grid:  gpd.GeoDataFrame,
    zones: dict,
) -> gpd.GeoDataFrame:
    """제외 구역 적용 후 최종 후보점 반환 (UTM)"""
    print("\n후보점 필터링...")
    candidates = grid.copy()

    labels = {
        "power":       f"[1] 송전 인프라 ({POWER_BUFFER_M//1000} km)",
        "railway":     f"[2] 철도 ({RAILWAY_BUFFER_M//1000} km)",
        "urban_dense": f"[3a] 핵심도심·산업 ({URBAN_DENSE_BUFFER_M}m)",
        "urban_resid": f"[3b] 주거·취락 ({URBAN_RESIDENT_BUFFER_M}m)",
        "pipeline":    f"[4] 파이프라인 ({PIPELINE_BUFFER_M//1000} km)",
        "comm":        f"[6] 통신탑 ({COMM_TOWER_BUFFER_M//1000} km)",
        "wind":        f"[7] 풍력발전기 ({WIND_BUFFER_M//1000} km)",
        "quarry":      f"[8] 채석장/광산 ({QUARRY_BUFFER_M//1000} km)",
        "water":       "[수계] 호수·저수지·강수면",
        "anomaly":     f"[9] 자기이상도 P90-P10 >{ANOMALY_EXCLUDE_THRESHOLD_NT} nT",
    }
    for key, geom in zones.items():
        if geom is None or geom.is_empty:
            continue
        before = len(candidates)
        mask = ~candidates.geometry.intersects(geom)
        candidates = candidates[mask].reset_index(drop=True)
        removed = before - len(candidates)
        if removed > 0:
            print(f"  {labels.get(key, key)}: -{removed}개 → 잔여 {len(candidates)}개")

    print(f"\n  ✅ 최종 후보점: {len(candidates)}개 / 전체 격자 {len(grid)}개")
    return candidates


def _slope_from_five(center, north, south, east, west,
                     lat_c: float, radius_deg: float) -> float:
    """
    중심(center) + 4방향(N/S/E/W) 표고 5점으로 경사도(°) 계산.
    중앙차분(central difference) 기반.

    dz/dy = (north - south) / (2 * d_lat_m)
    dz/dx = (east  - west)  / (2 * d_lon_m)
    slope  = arctan(sqrt(dz/dy² + dz/dx²))  [°]
    """
    lat_m = 111_000.0                                     # 위도 1° ≈ 111 km
    lon_m = 111_000.0 * np.cos(np.radians(lat_c))        # 경도 1° — 위도 보정
    d_lat_m = radius_deg * lat_m
    d_lon_m = radius_deg * lon_m
    if d_lat_m == 0 or d_lon_m == 0:
        return np.nan
    dz_dy = (north - south) / (2.0 * d_lat_m)
    dz_dx = (east  - west)  / (2.0 * d_lon_m)
    return float(np.degrees(np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))))


def _fetch_dem_via_srtm(lats, lons, radius_deg) -> np.ndarray | None:
    """
    로컬 SRTM GeoTIFF에서 경사도(°) 계산.
    data/srtm/ 디렉토리에 한반도 SRTM3 타일 배치 시 자동 사용.
    예: data/srtm/N33E124.hgt, N34E125.hgt, ...  (HGT 또는 GeoTIFF 가능)

    반환: shape (n,) 경사도 배열(°), 또는 rasterio 미설치/파일 없음 시 None
    """
    try:
        import rasterio
        from rasterio.merge import merge as rio_merge
        from rasterio.transform import rowcol
    except ImportError:
        return None

    srtm_dir = DATA_DIR / "srtm"
    if not srtm_dir.exists():
        return None

    tif_files = list(srtm_dir.glob("*.tif")) + list(srtm_dir.glob("*.TIF")) \
              + list(srtm_dir.glob("*.hgt")) + list(srtm_dir.glob("*.HGT"))
    if not tif_files:
        return None

    print(f"    SRTM 로컬 파일 사용 ({len(tif_files)}개 타일)...")
    try:
        src_list = [rasterio.open(str(f)) for f in tif_files]
        mosaic, transform = rio_merge(src_list)
        elev_arr = mosaic[0].astype(float)
        elev_arr[elev_arr < -9000] = np.nan  # SRTM nodata

        n = len(lats)
        r = radius_deg
        slopes = np.full(n, np.nan)
        for i in range(n):
            lat_c, lon_c = lats[i], lons[i]
            sample_pts = [
                (lat_c,      lon_c),       # center
                (lat_c + r,  lon_c),       # north
                (lat_c - r,  lon_c),       # south
                (lat_c,      lon_c + r),   # east
                (lat_c,      lon_c - r),   # west
            ]
            elevs = []
            for slat, slon in sample_pts:
                try:
                    row, col = rowcol(transform, slon, slat)
                    if 0 <= row < elev_arr.shape[0] and 0 <= col < elev_arr.shape[1]:
                        v = elev_arr[row, col]
                        elevs.append(float(v) if not np.isnan(v) else np.nan)
                    else:
                        elevs.append(np.nan)
                except Exception:
                    elevs.append(np.nan)
            if len(elevs) == 5 and not any(np.isnan(e) for e in elevs):
                slopes[i] = _slope_from_five(*elevs, lat_c=lat_c, radius_deg=r)

        for src in src_list:
            src.close()
        valid = (~np.isnan(slopes)).sum()
        print(f"    SRTM 경사도 계산 완료: {valid}/{n}개 유효")
        return slopes
    except Exception as exc:
        print(f"    ⚠ SRTM 처리 실패: {exc}")
        return None


def fetch_dem_slopes(
    candidates_wgs: gpd.GeoDataFrame,
    radius_deg: float = 0.045,   # ~5 km at Korea latitudes
    cache_file: Path = DATA_DIR / "dem_slopes.json",
) -> np.ndarray:
    """
    후보 지점의 지형 경사도(°) 취득.

    우선순위:
      1. 로컬 SRTM GeoTIFF (data/srtm/*.hgt 또는 *.tif)
         → rasterio 설치 및 data/srtm/ 타일 배치 시 사용
      2. Open-Elevation API (무료, 느림, 불안정)
         → SRTM 로컬 파일 없을 때 폴백

    각 후보점에 대해 중심 + 4방향(N/S/E/W, ~5 km) 총 5개 지점의 표고를 조회.
    중앙차분으로 경사도(°) 계산:
        dz/dx = (E - W) / (2 × d_lon_m)
        dz/dy = (N - S) / (2 × d_lat_m)
        slope = arctan(√(dz/dx² + dz/dy²))  [°]

    반환: shape (n,) numpy array — 후보점별 경사도 (°)
    """
    import json

    pts_wgs = candidates_wgs
    n = len(pts_wgs)
    lats = pts_wgs.geometry.y.values
    lons = pts_wgs.geometry.x.values

    # ── 캐시 확인 (재현성 강화: 좌표 해시 + 격자 간격 포함) ───────────
    import hashlib
    coord_bytes = (lats.tobytes() + lons.tobytes()
                   + np.array([GRID_SPACING_M]).tobytes())
    coord_hash = hashlib.md5(coord_bytes).hexdigest()[:12]
    cache_key = {
        "n": n,
        "radius_deg": radius_deg,
        "grid_spacing_m": GRID_SPACING_M,
        "coord_hash": coord_hash,
    }
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                cached = json.load(f)
            if (cached.get("n") == n
                    and cached.get("radius_deg") == radius_deg
                    and cached.get("coord_hash") == coord_hash):
                print(f"    DEM 경사도 캐시 로드: {cache_file.name} ({n}개 지점)")
                return np.array(cached["slopes"])
            else:
                print(f"    DEM 캐시 무효 (좌표/설정 변경) — 재취득")
        except Exception:
            pass

    # ── 1단계: 로컬 SRTM ─────────────────────────────────────
    slopes = _fetch_dem_via_srtm(lats, lons, radius_deg)

    # ── 2단계: Open-Elevation API fallback ───────────────────
    if slopes is None:
        print(f"    Open-Elevation API 조회 중 ({n}개 후보점 × 5방향)...")
        print(f"    ※ 로컬 SRTM 사용 시 data/srtm/ 에 한반도 HGT/GeoTIFF 타일 배치 권장")
        r = radius_deg
        all_locations = []
        for i in range(n):
            all_locations += [
                {"latitude": round(lats[i],      5), "longitude": round(lons[i],      5)},  # center
                {"latitude": round(lats[i] + r,  5), "longitude": round(lons[i],      5)},  # north
                {"latitude": round(lats[i] - r,  5), "longitude": round(lons[i],      5)},  # south
                {"latitude": round(lats[i],       5), "longitude": round(lons[i] + r, 5)},  # east
                {"latitude": round(lats[i],       5), "longitude": round(lons[i] - r, 5)},  # west
            ]
        chunk_size = 1000
        elevations = []
        url = "https://api.open-elevation.com/api/v1/lookup"
        total_locs = len(all_locations)

        for start in range(0, total_locs, chunk_size):
            chunk = all_locations[start : start + chunk_size]
            for attempt in range(3):
                try:
                    resp = requests.post(
                        url,
                        json={"locations": chunk},
                        timeout=60,
                        headers={"Content-Type": "application/json"},
                    )
                    resp.raise_for_status()
                    results = resp.json().get("results", [])
                    elevations += [r_["elevation"] for r_ in results]
                    pct = min(start + chunk_size, total_locs)
                    print(f"      {pct}/{total_locs} 지점 완료", end="\r")
                    break
                except Exception as exc:
                    if attempt < 2:
                        print(f"\n      재시도 {attempt+1}/3 ({exc})")
                        time.sleep(10)
                    else:
                        print(f"\n      ✗ 청크 실패 — 0m 대체")
                        elevations += [0] * len(chunk)

        print(f"\n    DEM 표고 취득 완료: {len(elevations)}개")
        slopes_list = []
        elev_arr = np.array(elevations, dtype=float)
        for i in range(n):
            s = i * 5
            five = elev_arr[s : s + 5]
            if len(five) == 5 and not np.any(np.isnan(five)):
                # center=five[0], north=five[1], south=five[2], east=five[3], west=five[4]
                deg = _slope_from_five(five[0], five[1], five[2], five[3], five[4],
                                       lat_c=lats[i], radius_deg=r)
                slopes_list.append(deg)
            else:
                slopes_list.append(np.nan)
        slopes = np.array(slopes_list)

    # ── 캐시 저장 ─────────────────────────────────────────────
    try:
        with open(cache_file, "w") as f:
            json.dump({
                **cache_key,
                "slopes": [float(v) for v in slopes],
            }, f)
        print(f"    DEM 경사도 캐시 저장: {cache_file.name}")
    except Exception as exc:
        print(f"    ⚠ DEM 캐시 저장 실패: {exc}")

    return slopes


def compute_priority(
    candidates:   gpd.GeoDataFrame,
    zones:        dict,
    emag2_data=None,
) -> gpd.GeoDataFrame:
    """
    입지 점수 산정 (0~100점, UTM CRS 입력).

    ┌─────────────┬──────────────────┬────┬──────────────────────────────┐
    │ 평가 항목    │ 세부 지표         │배점│ 가용 여부                     │
    ├─────────────┼──────────────────┼────┼──────────────────────────────┤
    │ 공간 대표성  │ 격자 데이터 희소성 │ 25 │ ✅ 후보점 분포 분석            │
    │             │ 지형 경사도        │ 15 │ ⚠  Open-Elevation API        │
    │ 환경 정온도  │ 전력/철도 이격도   │ 15 │ ✅ OSM 데이터                 │
    │             │ 인구 밀집 이격도   │ 15 │ ✅ OSM 데이터                 │
    │ 지질 안정성  │ 자기 이상 균일도   │ 10 │ ⚠  KIGAM/EMAG2 배치 시       │
    │             │ 암상 적합성       │  5 │ ❌ 지질도 미확보               │
    │ 운영 인프라  │ 부지 지속성       │ 10 │ ※ 최종 선정 후 육안 확인       │
    │             │ 관리 접근성       │  5 │ ※ 최종 선정 후 지도 확인       │
    └─────────────┴──────────────────┴────┴──────────────────────────────┘
    ⑦ 부지 지속성 / ⑧ 관리 접근성은 최종 선정 후보지에서 지도·현장 육안 확인.
    가용 항목 점수를 100점 만점으로 정규화하여 등급 분류.
    """
    if len(candidates) == 0:
        return candidates

    result = candidates.copy()
    pts = candidates.geometry
    n = len(pts)
    print("  입지 점수 산정 중...")

    # ── ① 격자 데이터 희소성 (25점) ──────────────────────────
    # 주변 후보점 밀도의 역수 — 고립된(이웃이 먼) 지점일수록 높은 점수
    try:
        from scipy.spatial import KDTree
        xy = np.column_stack([pts.x, pts.y])
        tree = KDTree(xy)
        k = min(6, n)
        dists_k, _ = tree.query(xy, k=k)
        avg_nn = dists_k[:, 1:].mean(axis=1)   # self 제외
        mx = avg_nn.max() if avg_nn.max() > 0 else 1
        s1 = (avg_nn / mx) * 25
    except ImportError:
        s1 = np.full(n, 12.5)
        print("    ⚠ scipy 미설치 — 희소성 중간값 적용")
    result["s1_희소성"] = np.round(s1, 1)
    print(f"    ① 희소성: 평균 {s1.mean():.1f} / 25점")

    # ── ② 지형 경사도 (15점) — Open-Elevation API / SRTM ────
    # 후보점 중심 + 4방향(~5km) 표고 5점으로 중앙차분 경사도(°) 산정
    # 경사도 낮음(평탄) → 설치·측정 유리 → 높은 점수
    # 경사도 높음(급경사) → 불리 → 낮은 점수
    # 정규화: s2 = clip(1 - slope / SLOPE_MAX, 0, 1) × 15
    #   SLOPE_MAX = 30° (한국 산지 기준 실질적 상한)
    SLOPE_MAX_DEG = 30.0
    s2 = np.full(n, np.nan)
    dem_available = False
    try:
        cands_wgs = candidates.to_crs(WGS84_CRS)
        dem_slopes = fetch_dem_slopes(cands_wgs)
        valid_dem = ~np.isnan(dem_slopes)
        if valid_dem.sum() > 0:
            s2[valid_dem] = np.clip(
                (1.0 - dem_slopes[valid_dem] / SLOPE_MAX_DEG) * 15.0,
                0.0, 15.0,
            )
            s2[~valid_dem] = 7.5    # 취득 실패 지점 중간값 대체
            dem_available = True
            result["dem_slope_deg"] = np.round(dem_slopes, 2)
            print(f"    ② 지형 경사도: 평균 {np.nanmean(s2):.1f} / 15점  "
                  f"(경사도 평균 {np.nanmean(dem_slopes):.1f}°)")
    except Exception as exc:
        print(f"    ⚠ DEM 취득 실패 ({exc}) — ② 지형 미산정")
    result["s2_지형"] = np.round(s2, 1)

    # ── ③ 전력/철도 이격도 (15점) ──────────────────────────
    # log(dist) 기반 점수화: 제외 구역 경계로부터의 이격 거리
    def _zone_distances(zone_geom, label):
        """제외 구역 경계까지의 거리(m) 배열 반환"""
        if zone_geom is None or zone_geom.is_empty:
            return np.full(n, 50_000.0)
        # 복잡 geometry 간소화 (500m 허용오차) — 성능 최적화
        simp = zone_geom.simplify(500)
        print(f"    {label} 이격 거리 계산 ({n}개)...")
        return np.array([pt.distance(simp) for pt in pts])

    d_power   = _zone_distances(zones.get("power"),   "전력")
    d_railway = _zone_distances(zones.get("railway"),  "철도")

    lp = np.log1p(d_power)
    lr = np.log1p(d_railway)
    mx_p = lp.max() if lp.max() > 0 else 1
    mx_r = lr.max() if lr.max() > 0 else 1
    s3 = (lp / mx_p * 0.5 + lr / mx_r * 0.5) * 15
    result["s3_전력철도"] = np.round(s3, 1)
    result["d_power_km"]  = np.round(d_power  / 1000, 1)
    result["d_railway_km"] = np.round(d_railway / 1000, 1)
    print(f"    ③ 전력/철도: 평균 {s3.mean():.1f} / 15점")

    # ── ④ 인구 밀집 이격도 (15점) ──────────────────────────
    d_urban = _zone_distances(zones.get("urban_resid"), "도시·주거")
    lu = np.log1p(d_urban)
    mx_u = lu.max() if lu.max() > 0 else 1
    s4 = (lu / mx_u) * 15
    result["s4_인구이격"] = np.round(s4, 1)
    result["d_urban_km"]  = np.round(d_urban / 1000, 1)
    print(f"    ④ 인구이격: 평균 {s4.mean():.1f} / 15점")

    # ── ⑤ 자기 이상 균일도 (10점) ─ 고정 임계값 기반 점수화 ──────
    # 기준: KIGAM 1.5분 격자, 반경 0.05°(≈5.5km) 내 P90-P10
    # 고정 임계값 점수화 (상대정규화 X):
    #   ≤ 50 nT  → 10점 (최우수)
    #   ≤ 100 nT → 8점  (우수)
    #   ≤ 150 nT → 5점  (보통)
    #   ≤ 200 nT → 2점  (현장검토)
    #   > 200 nT → 0점  (제외 — filter_candidates에서 이미 제거됨)
    s5 = np.full(n, np.nan)
    emag_available = False
    if emag2_data is not None:
        try:
            e_lons = emag2_data["lon"].values
            e_lats = emag2_data["lat"].values
            e_anom = emag2_data["anomaly_nT"].values
            cw = candidates.to_crs(WGS84_CRS)
            c_lats = cw.geometry.y.values
            c_lons = cw.geometry.x.values
            radius = ANOMALY_SITE_RADIUS_DEG  # 0.05° — 동일 스케일 통일
            p90p10 = np.full(n, np.nan)
            for i, (clat, clon) in enumerate(zip(c_lats, c_lons)):
                mask = (np.abs(e_lats - clat) <= radius) & (np.abs(e_lons - clon) <= radius)
                pts = e_anom[mask]
                if len(pts) >= 3:
                    p90p10[i] = np.percentile(pts, 90) - np.percentile(pts, 10)
            valid = ~np.isnan(p90p10)
            if valid.sum() > 0:
                # 고정 임계값 기반 점수 (상대정규화 대신 절대 기준)
                def _score_variability(v):
                    if v <= 50:   return 10.0
                    if v <= 100:  return 8.0
                    if v <= 150:  return 5.0
                    if v <= 200:  return 2.0
                    return 0.0
                s5[valid] = np.array([_score_variability(v) for v in p90p10[valid]])
                s5[~valid] = 5.0   # 데이터 없는 지점 중간값
                emag_available = True
                result["mag_p90p10_nT"] = np.round(p90p10, 1)
                print(f"    ⑤ 자기균일: 평균 {np.nanmean(s5):.1f} / 10점  "
                      f"(P90-P10 반경 {radius}°, 고정임계값 기반)")
                print(f"    ⚠  KIGAM 1.5분 광역 자력이상도 예비선정 — 최종 확정은 현장 정밀 자력측량 필요")
        except Exception as exc:
            print(f"    ⚠ 자기균일 계산 실패: {exc}")
    result["s5_자기균일"] = np.round(s5, 1)

    # ── ⑥ 암상 적합성 (5점) — 지질도 미확보 ─────────────────
    result["s6_암상"] = np.nan
    # ⑦ 부지 지속성 (10점): 최종 선정 후보지에서 국공유지 여부 육안 확인
    # ⑧ 관리 접근성  (5점): 최종 선정 후보지에서 도로망 지도 육안 확인

    # ── 종합 점수 (가용 항목 합산 → 100점 정규화) ────────────
    total = s1 + s3 + s4                       # 기본 55점 만점
    available_max = 25 + 15 + 15               # = 55
    if dem_available:
        total = total + np.nan_to_num(s2, nan=0)
        available_max += 15                    # +15 = 70
    if emag_available:
        total = total + np.nan_to_num(s5, nan=0)
        available_max += 10                    # +10 (최대 80점)

    result["score"] = np.round(total / available_max * 100, 1)
    result["score_max"] = available_max

    # 3등급 분류
    scores = result["score"].values
    p33, p66 = np.percentile(scores, [33, 66])
    result["priority"] = np.where(scores >= p66, 1,
                          np.where(scores >= p33, 2, 3))

    for p in [1, 2, 3]:
        cnt = (result["priority"] == p).sum()
        label = ["최우선", "우선", "일반"][p - 1]
        print(f"  우선순위 {p}등급({label}): {cnt}개")

    unavail_list = ["⑥암상(5)", "⑦부지(10,육안확인)", "⑧접근성(5,육안확인)"]
    if not dem_available:
        unavail_list.insert(0, "②지형(15)")
    print(f"  가용 항목 합계: {available_max}/100점 → 100점 정규화")
    print(f"  미산정/후확인 항목: {', '.join(unavail_list)}")

    return result


# ============================================================
# 지도 시각화 (Folium)
# ============================================================

def _geom_to_wgs84(geom_utm, crs=UTM_CRS) -> gpd.GeoDataFrame:
    """단일 shapely geometry → WGS84 GeoDataFrame"""
    return gpd.GeoDataFrame(geometry=[geom_utm], crs=crs).to_crs(WGS84_CRS)


def load_existing_sites() -> "pd.DataFrame | None":
    """
    지자기측량 성과정리(22_25).xlsx 에서 기존 15개 측정점 정보를 로드.
    각 측점별 최신 관측연도 행을 반환한다.
    """
    xlsx_path = DATA_DIR / "지자기측량 성과정리(22_25).xlsx"
    if not xlsx_path.exists():
        print(f"  ⚠ 기존 측정점 파일 없음: {xlsx_path}")
        return None
    try:
        df = pd.read_excel(xlsx_path, sheet_name="Sheet1", header=0, engine="openpyxl")
        # 열 인덱스(0-based) 기준 필요 컬럼 추출
        # A=0연번, B=1도엽명, C=2주소, D=3최초설치, E=4관측연도
        # I=8위도실수, M=12경도실수, N=13표고
        # R=17편각실수, V=21복각실수, W=22총자력
        df_work = pd.DataFrame({
            "연번":     df.iloc[:, 0],
            "도엽명":   df.iloc[:, 1],
            "주소":     df.iloc[:, 2],
            "최초설치": df.iloc[:, 3],
            "관측연도": df.iloc[:, 4],
            "위도":     df.iloc[:, 8],
            "경도":     df.iloc[:, 12],
            "표고":     df.iloc[:, 13],
            "편각":     df.iloc[:, 17],
            "복각":     df.iloc[:, 21],
            "총자력":   df.iloc[:, 22],
        })
        # 도엽명·주소·최초설치는 첫 행에만 존재 → forward fill
        df_work["도엽명"]   = df_work["도엽명"].ffill()
        df_work["주소"]     = df_work["주소"].ffill()
        df_work["최초설치"] = df_work["최초설치"].ffill()
        # 위도·경도·관측연도가 유효한 행만 유지
        df_work["위도"]     = pd.to_numeric(df_work["위도"],     errors="coerce")
        df_work["경도"]     = pd.to_numeric(df_work["경도"],     errors="coerce")
        df_work["관측연도"] = pd.to_numeric(df_work["관측연도"], errors="coerce")
        df_work = df_work.dropna(subset=["위도", "경도", "관측연도"])
        # 각 도엽명별 최신 관측연도 행 선택
        idx_max    = df_work.groupby("도엽명")["관측연도"].idxmax()
        df_latest  = df_work.loc[idx_max].reset_index(drop=True)
        print(f"  ✅ 기존 측정점 {len(df_latest)}개 로드 완료 (최신 관측연도 기준)")
        return df_latest
    except Exception as exc:
        print(f"  ⚠ 기존 측정점 로드 실패: {exc}")
        return None


def save_map_data(
    zones:          dict,
    grid:           gpd.GeoDataFrame,
    final_cands:    gpd.GeoDataFrame,
    existing_sites: "pd.DataFrame | None",
    korea_gdf:      "gpd.GeoDataFrame | None",
    data_dir:       Path,
) -> None:
    """
    지도 레이어 데이터를 output/data/ 에 GeoJSON/JSON 파일로 저장.
    HTML 의 JS 가 fetch() 로 로드하여 지도에 표시한다.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"  데이터 파일 저장 중: {data_dir.resolve()}")

    # ── 제외 구역 (UTM → WGS84 변환 후 저장) ─────────────────
    for key in ("power", "railway", "urban_dense", "urban_resid",
                "pipeline", "comm", "wind", "quarry", "anomaly"):
        geom = zones.get(key)
        if geom is None or geom.is_empty:
            continue
        gdf = gpd.GeoDataFrame(geometry=[geom], crs=UTM_CRS).to_crs(WGS84_CRS)
        out = data_dir / f"zone_{key}.geojson"
        gdf.to_file(str(out), driver="GeoJSON")

    # ── 격자점 샘플 ─────────────────────────────────────────
    grid_wgs = grid.to_crs(WGS84_CRS)
    sample   = grid_wgs.sample(min(len(grid_wgs), 600), random_state=42)
    gpd.GeoDataFrame(geometry=sample.geometry, crs=WGS84_CRS).to_file(
        str(data_dir / "grid_sample.geojson"), driver="GeoJSON"
    )

    # ── 후보 지점 (우선순위별, 속성 포함) ────────────────────
    final_wgs    = final_cands.to_crs(WGS84_CRS).copy().reset_index(drop=True)
    has_priority = "priority" in final_wgs.columns
    final_wgs["idx"] = final_wgs.index + 1
    final_wgs["lat"] = final_wgs.geometry.y
    final_wgs["lon"] = final_wgs.geometry.x
    rename_map = {
        "s1_희소성": "s1", "s2_지형": "s2", "s3_전력철도": "s3",
        "s4_인구이격": "s4", "s5_자기균일": "s5",
        "d_power_km": "dp", "d_railway_km": "dr", "d_urban_km": "du",
        "dem_slope_deg": "dem", "d_public_km": "dpub", "d_road_km": "drd",
    }
    prop_cols = ["idx", "lat", "lon"]
    for c in ["priority", "score"] + list(rename_map.keys()):
        if c in final_wgs.columns:
            prop_cols.append(c)
    for p in (1, 2, 3):
        if has_priority:
            sub = final_wgs[final_wgs["priority"] == p][prop_cols].copy()
        else:
            sub = final_wgs[prop_cols].copy() if p == 3 else pd.DataFrame()
        if len(sub) == 0:
            continue
        sub = sub.rename(columns={k: v for k, v in rename_map.items() if k in sub.columns})
        geom_col = final_wgs.geometry.iloc[sub.index]
        out_gdf  = gpd.GeoDataFrame(sub, geometry=geom_col.values, crs=WGS84_CRS)
        out_gdf.to_file(str(data_dir / f"candidates_p{p}.geojson"), driver="GeoJSON")

    # ── 기존 측정점 ─────────────────────────────────────────
    if existing_sites is not None and len(existing_sites) > 0:
        features = []
        for _, r in existing_sites.iterrows():
            props = {
                "name":      str(r["도엽명"]),
                "address":   str(r["주소"])     if pd.notna(r["주소"])     else None,
                "inst_year": int(r["최초설치"]) if pd.notna(r["최초설치"]) else None,
                "obs_year":  int(r["관측연도"]),
                "lat":       float(r["위도"]),
                "lon":       float(r["경도"]),
                "elev":      float(r["표고"])   if pd.notna(r["표고"])   else None,
                "decl":      float(r["편각"])   if pd.notna(r["편각"])   else None,
                "incl":      float(r["복각"])   if pd.notna(r["복각"])   else None,
                "total":     float(r["총자력"]) if pd.notna(r["총자력"]) else None,
            }
            features.append({
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [props["lon"], props["lat"]]},
                "properties": props,
            })
        with open(str(data_dir / "existing_sites.geojson"), "w", encoding="utf-8") as fh:
            json.dump({"type": "FeatureCollection", "features": features}, fh,
                      ensure_ascii=False, indent=None)

    # ── 1:50,000 도엽 격자 ──────────────────────────────────
    if korea_gdf is not None:
        try:
            print("  1:50,000 도엽 GeoJSON 저장 중...")
            topo_gdf = create_topo_sheet_grid(korea_gdf)
            topo_gdf = topo_gdf.copy()
            topo_gdf["cx"] = topo_gdf.geometry.centroid.x
            topo_gdf["cy"] = topo_gdf.geometry.centroid.y
            for col in ["sheet_mapidcd", "sheet_name", "sheet_code"]:
                if col not in topo_gdf.columns:
                    topo_gdf[col] = ""
            topo_gdf["sheet_mapidcd"] = (
                topo_gdf["sheet_mapidcd"].astype(str).str.strip()
                .replace({"nan": "", "None": ""})
            )
            keep = [c for c in ["sheet_name", "sheet_code", "sheet_mapidcd", "cx", "cy"]
                    if c in topo_gdf.columns]
            topo_gdf[keep + ["geometry"]].to_file(
                str(data_dir / "topo_sheets.geojson"), driver="GeoJSON"
            )
            print(f"  ✅ 도엽 GeoJSON 저장: {len(topo_gdf)}개")
        except Exception as exc:
            print(f"  ⚠ 도엽 GeoJSON 저장 실패: {exc}")

    print(f"  ✅ 데이터 저장 완료 ({data_dir.resolve()})")


def create_folium_map(
    zones:          dict,
    grid:           gpd.GeoDataFrame,
    final_cands:    gpd.GeoDataFrame,   # UTM, priority 컬럼 포함
    korea_gdf:      gpd.GeoDataFrame | None = None,
    existing_sites: "pd.DataFrame | None" = None,
    data_subdir:    str = "data",
) -> folium.Map:
    """대화형 Folium 지도 생성 (OSM 기반)"""
    print("\n지도 생성 중...")

    m = folium.Map(
        location=[36.5, 127.5],
        zoom_start=7,
        tiles=None,
        prefer_canvas=True,
    )

    # ── 기본 타일 레이어 ──────────────────────────────────────
    folium.TileLayer(
        tiles="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
        attr='© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name="OpenStreetMap",
        overlay=False,
        control=True,
        max_zoom=19,
    ).add_to(m)

    folium.TileLayer(
        tiles=(
            "https://server.arcgisonline.com/ArcGIS/rest/services"
            "/World_Imagery/MapServer/tile/{z}/{y}/{x}"
        ),
        attr="Esri World Imagery",
        name="위성 이미지 (Esri)",
        overlay=False,
        control=True,
    ).add_to(m)

    # ── 제외 구역 레이어 (외부 GeoJSON fetch) ───────────────────
    _excl = {
        "power":       ("#FF3300", "#CC0000", "⚡ [1] 고압철탑·송전탑 제외 (1.0 km)",    True),
        "railway":     ("#FF7700", "#CC4400", "🚆 [2] 철도 제외 (5.0 km)",               True),
        "urban_dense": ("#AA00FF", "#7700CC", "🏙️ [3a] 핵심도심·산업 제외 (500m)",       True),
        "urban_resid": ("#CC66FF", "#9933CC", "🏘️ [3b] 주거·취락 제외 (300m)",           True),
        "pipeline":    ("#FF00AA", "#CC0088", "🔩 [4] 파이프라인 제외 (0.5 km)",         True),
        "comm":        ("#0088FF", "#0055CC", "📡 [5] 통신탑·기지국 제외 (0.5 km)",      True),
        "wind":        ("#00CCAA", "#009977", "💨 [6] 풍력발전기 제외 (0.5 km)",         True),
        "quarry":      ("#AA6600", "#774400", "⛏️ [7] 채석장·광산 제외 (1.0 km)",        True),
        "anomaly":     ("#0044FF", "#0022CC", "🧲 [8] 자기이상도 고변화 제외 (KIGAM)",    True),
    }
    for key, (fc, ec, lname, show) in _excl.items():
        geom = zones.get(key)
        if geom is None or geom.is_empty:
            continue
        layer = folium.FeatureGroup(name=lname, show=show)
        layer.add_to(m)
        lv  = layer.get_name()
        tip = lname.replace("'", "\\'")
        js  = (
            "<script>(function(){"
            "fetch('%s/zone_%s.geojson')"
            ".then(function(r){if(!r.ok)throw new Error(r.status);return r.json();})"
            ".then(function(d){"
            "L.geoJSON(d,{"
            "style:function(){return{fillColor:'%s',color:'%s',weight:1,fillOpacity:0.32}},"
            "onEachFeature:function(f,l){l.bindTooltip('%s');}"
            "}).addTo(%s);"
            "}).catch(function(e){console.warn('zone_%s:',e.message);});"
            "})();</script>"
        ) % (data_subdir, key, fc, ec, tip, lv, key)
        m.get_root().html.add_child(folium.Element(js))

    # ── 전체 격자점 (참고용, 외부 GeoJSON fetch) ─────────────
    grid_layer = folium.FeatureGroup(name="격자점 전체 (참고)", show=False)
    grid_layer.add_to(m)
    gv  = grid_layer.get_name()
    js  = (
        "<script>(function(){"
        "fetch('%s/grid_sample.geojson')"
        ".then(function(r){return r.json();})"
        ".then(function(d){"
        "L.geoJSON(d,{pointToLayer:function(f,ll){"
        "return L.circleMarker(ll,{radius:2,color:'#888888',"
        "fillColor:'#AAAAAA',fillOpacity:0.3,weight:0});"
        "}}).addTo(%s);"
        "}).catch(function(e){console.warn('grid:',e.message);});"
        "})();</script>"
    ) % (data_subdir, gv)
    m.get_root().html.add_child(folium.Element(js))

    # ── 최종 후보 지점 (우선순위별, 외부 GeoJSON fetch) ────────
    priority_cfg = {
        1: ("#FF0000", "#CC0000", "🔴 우선순위 1등급 (최우선, 모델 공백 지역)"),
        2: ("#FF8800", "#CC5500", "🟠 우선순위 2등급 (우선)"),
        3: ("#00BB00", "#008800", "🟢 우선순위 3등급 (일반)"),
    }
    final_wgs    = final_cands.to_crs(WGS84_CRS)
    has_priority = "priority" in final_wgs.columns

    for p, (fc, ec, pname) in priority_cfg.items():
        if has_priority:
            has_data = int((final_wgs["priority"] == p).sum()) > 0
        else:
            has_data = (p == 3) and len(final_wgs) > 0
        if not has_data:
            continue
        plabel = "최우선" if p == 1 else ("우선" if p == 2 else "일반")
        radius = 7 if p == 1 else (6 if p == 2 else 5)
        p_layer = folium.FeatureGroup(name=pname, show=True)
        p_layer.add_to(m)
        pv  = p_layer.get_name()
        js  = (
            "<script>(function(){"
            "fetch('%s/candidates_p%d.geojson')"
            ".then(function(r){return r.json();})"
            ".then(function(d){"
            "L.geoJSON(d,{"
            "pointToLayer:function(f,ll){"
            "return L.circleMarker(ll,{radius:%d,color:'%s',"
            "fillColor:'%s',fillOpacity:0.85,weight:1.5});"
            "},"
            "onEachFeature:function(f,l){"
            "var p=f.properties;"
            "var vn=function(x){return(x===null||x===undefined||isNaN(+x))?'-':(+x).toFixed(1);};"
            "var lat=(p.lat!=null)?p.lat:f.geometry.coordinates[1];"
            "var lon=(p.lon!=null)?p.lon:f.geometry.coordinates[0];"
            "var html='<div style=\"font-family:sans-serif;font-size:12.5px;min-width:250px;\">'"
            "+'<b style=\"color:%s;\">측정 후보지 #'+p.idx+'</b>&nbsp;'"
            "+'<span style=\"background:%s;color:white;padding:1px 5px;border-radius:3px;font-size:11px;\">%s</span><br>'"
            "+'<hr style=\"margin:4px 0;\">'"
            "+'<b>위도:</b> '+lat.toFixed(5)+'° N &nbsp; <b>경도:</b> '+lon.toFixed(5)+'° E<br>'"
            "+'<hr style=\"margin:4px 0;border-color:#ddd;\">'"
            "+'<b>입지 점수: '+vn(p.score)+' / 100</b><br>'"
            "+'<span style=\"font-size:11.5px;color:#333;\">'"
            "+'&nbsp;① 희소성: '+vn(p.s1)+' / 25<br>'"
            "+'&nbsp;② 지형적 대표성: '+vn(p.s2)+' / 15<br>'"
            "+'&nbsp;③ 전력·철도 이격: '+vn(p.s3)+' / 15 '"
            "+'<span style=\"color:#888;\">(전력 '+vn(p.dp)+'km, 철도 '+vn(p.dr)+'km)</span><br>'"
            "+'&nbsp;④ 인구 이격: '+vn(p.s4)+' / 15 '"
            "+'<span style=\"color:#888;\">(도시 '+vn(p.du)+'km)</span><br>'"
            "+'&nbsp;⑤ 자기균일: '+vn(p.s5)+' / 10<br>'"
            "+'&nbsp;<span style=\"color:#999;\">⑦ 부지 지속성: 현장검토 / 10</span><br>'"
            "+'&nbsp;<span style=\"color:#999;\">⑧ 관리 접근성: 현장검토 / 5</span><br>'"
            "+'&nbsp;<span style=\"color:#999;\">⑥암상: 미산정</span>'"
            "+'</span></div>';"
            "l.bindPopup(html,{maxWidth:270});"
            "l.bindTooltip('후보지 #'+p.idx+' P%d ('+lat.toFixed(4)+'°N, '+lon.toFixed(4)+'°E)');"
            "}"
            "}).addTo(%s);"
            "}).catch(function(e){console.warn('candidates_p%d:',e.message);});"
            "})();</script>"
        ) % (data_subdir, p, radius, ec, fc, ec, fc, plabel, p, pv, p)
        m.get_root().html.add_child(folium.Element(js))

    # ── 기존 측정점 레이어 (외부 GeoJSON fetch) ──────────────
    exist_layer = folium.FeatureGroup(name="⭐ 기존 측정점 (22-25년)", show=True)
    exist_layer.add_to(m)
    ev  = exist_layer.get_name()
    js  = (
        "<script>(function(){"
        "fetch('%s/existing_sites.geojson')"
        ".then(function(r){return r.json();})"
        ".then(function(d){"
        "L.geoJSON(d,{"
        "pointToLayer:function(f,ll){"
        "var icon=L.divIcon({"
        "html:'<div style=\"font-size:22px;line-height:1;"
              "text-shadow:1px 1px 2px rgba(0,0,0,0.4);\">⭐</div>',"
        "className:'',iconAnchor:[11,11]});"
        "return L.marker(ll,{icon:icon});"
        "},"
        "onEachFeature:function(f,l){"
        "var p=f.properties;"
        "var vf=function(x,fmt){return(x===null||x===undefined)?'-':fmt(x);};"
        "var html='<div style=\"font-family:sans-serif;font-size:12.5px;min-width:260px;\">'"
        "+'<b style=\"color:#8B4513;\">⭐ 기존 측정점: '+p.name+'</b><br>'"
        "+'<hr style=\"margin:4px 0;\">'"
        "+'<b>위도:</b> '+p.lat.toFixed(6)+'° N &nbsp; <b>경도:</b> '+p.lon.toFixed(6)+'° E<br>'"
        "+'<b>주소:</b> <span style=\"font-size:11px;\">'+(p.address||'-')+'</span><br>'"
        "+'<hr style=\"margin:4px 0;border-color:#ddd;\">'"
        "+'<b>최초설치:</b> '+(p.inst_year?p.inst_year+'년':'-')+' &nbsp; <b>최신관측:</b> '+p.obs_year+'년<br>'"
        "+'<b>표고:</b> '+vf(p.elev,function(x){return x.toFixed(1)+' m'})+'<br>'"
        "+'<hr style=\"margin:4px 0;border-color:#ddd;\">'"
        "+'<b>측정값 ('+p.obs_year+'년)</b><br>'"
        "+'<span style=\"font-size:11.5px;color:#333;\">'"
        "+'&nbsp;편각: '+vf(p.decl,function(x){return x.toFixed(4)+'°'})+'<br>'"
        "+'&nbsp;복각: '+vf(p.incl,function(x){return x.toFixed(4)+'°'})+'<br>'"
        "+'&nbsp;총자력: '+vf(p.total,function(x){"
        "return x.toLocaleString(undefined,{minimumFractionDigits:1,maximumFractionDigits:1})+' nT';"
        "})"
        "+'</span></div>';"
        "l.bindPopup(html,{maxWidth:290});"
        "l.bindTooltip('⭐ '+p.name+' ('+p.obs_year+'년 관측)');"
        "}"
        "}).addTo(%s);"
        "}).catch(function(e){console.warn('existing_sites:',e.message);});"
        "})();</script>"
    ) % (data_subdir, ev)
    m.get_root().html.add_child(folium.Element(js))

    # ── 범례 ─────────────────────────────────────────────────
    n_cands = len(final_wgs)
    n_p1 = int((final_wgs["priority"] == 1).sum()) if has_priority else 0
    n_p2 = int((final_wgs["priority"] == 2).sum()) if has_priority else 0
    n_p3 = int((final_wgs["priority"] == 3).sum()) if has_priority else n_cands

    def _swatch(color, label):
        return (f'<span style="background:{color};display:inline-block;'
                f'width:12px;height:12px;opacity:0.8;vertical-align:middle;'
                f'margin-right:4px;"></span>{label}<br>')

    def _dot(color, label):
        return (f'<span style="background:{color};display:inline-block;'
                f'width:12px;height:12px;border-radius:50%;opacity:0.85;'
                f'vertical-align:middle;margin-right:4px;"></span>{label}<br>')

    legend_html = f"""
    <div style="
        position:fixed; bottom:30px; left:30px; width:280px;
        background:rgba(255,255,255,0.96); border:2px solid #555;
        z-index:9999; padding:12px 14px; border-radius:8px;
        font-family:sans-serif; font-size:11.5px; line-height:1.65;
        box-shadow:2px 2px 8px rgba(0,0,0,0.25);">
      <b style="font-size:13px;">🗺️ 대한민국 지구자기장 모델<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;측정 입지 선정</b>
      <hr style="margin:6px 0;border-color:#ccc;">
      <b>▸ 제외 구역</b><br>
      {_swatch('#FF3300','[1] 고압철탑·송전탑 1.0 km')}
      {_swatch('#FF7700','[2] 철도 5.0 km')}
      {_swatch('#AA00FF','[3a] 핵심도심·산업 500m')}
      {_swatch('#CC66FF','[3b] 주거·취락 300m')}
      {_swatch('#FF00AA','[4] 파이프라인 0.5 km')}
      {_swatch('#0088FF','[5] 통신탑·기지국 0.5 km')}
      {_swatch('#00CCAA','[6] 풍력발전기 0.5 km')}
      {_swatch('#AA6600','[7] 채석장·광산 1.0 km')}
      {_swatch('#0044FF','[8] 자기이상도 P90-P10>200nT *')}
      <hr style="margin:6px 0;border-color:#ccc;">
      <b>▸ 측정 후보지 (총 {n_cands}개)</b><br>
      {_dot('#FF0000',f'1등급 최우선 (데이터 공백): {n_p1}개')}
      {_dot('#FF8800',f'2등급 우선: {n_p2}개')}
      {_dot('#00BB00',f'3등급 일반: {n_p3}개')}
      <hr style="margin:6px 0;border-color:#ccc;">
      <b>▸ 기존 측정점</b><br>
      <span style="display:inline-block;vertical-align:middle;margin-right:4px;">⭐</span>기존 측정점 15개 (최신 관측연도)<br>
      <hr style="margin:6px 0;border-color:#ccc;">
      <small style="color:#555;">
        격자 간격: {GRID_SPACING_M//1000} km | 좌표계: WGS84/EPSG:5179<br>
        데이터: OpenStreetMap (Overpass API)<br>
        * KIGAM 배치 시 활성화<br>
        ⚠ 광역 자력이상도 예비선정.<br>
        &nbsp;&nbsp;최종 확정은 현장 정밀<br>
        &nbsp;&nbsp;자력측량으로 검증 필요.
      </small>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    # ── 타이틀 ────────────────────────────────────────────────
    m.get_root().html.add_child(folium.Element("""
    <div style="
        position:fixed; top:10px; left:50%; transform:translateX(-50%);
        z-index:9999; background:rgba(255,255,255,0.93);
        padding:8px 22px; border-radius:6px; border:1px solid #999;
        font-family:sans-serif; font-size:15px; font-weight:bold;
        box-shadow:1px 1px 5px rgba(0,0,0,0.2);">
      🧭 대한민국 지구자기장 총자력 측정 입지 선정 지도
    </div>
    """))

    # ── 입지 점수 산정 기준 패널 (지도 우측 하단 고정) ────────────
    scoring_html = """
    <div style="
        position:fixed; bottom:30px; right:20px; width:290px;
        background:rgba(255,255,255,0.96); border:2px solid #446;
        z-index:9999; padding:11px 13px; border-radius:8px;
        font-family:sans-serif; font-size:11px; line-height:1.65;
        box-shadow:2px 2px 8px rgba(0,0,0,0.22);">
      <b style="font-size:12.5px;">📊 입지 점수 산정 기준</b>
      <hr style="margin:5px 0;border-color:#ccc;">
      <table style="width:100%;border-collapse:collapse;font-size:10.5px;">
        <tr style="border-bottom:1px solid #ddd;">
          <th style="text-align:left;padding:2px;">평가 항목</th>
          <th style="text-align:left;padding:2px;">세부 지표</th>
          <th style="text-align:center;padding:2px;">배점</th>
          <th style="text-align:center;padding:2px;">산정</th>
        </tr>
        <tr><td rowspan="2">공간 대표성</td>
            <td>① 격자 데이터 희소성</td>
            <td style="text-align:center;">25</td>
            <td style="text-align:center;">✅</td></tr>
        <tr><td>② 지형 경사도</td>
            <td style="text-align:center;">15</td>
            <td style="text-align:center;">✅</td></tr>
        <tr style="border-top:1px solid #eee;">
            <td rowspan="2">환경 정온도</td>
            <td>③ 전력/철도 이격도</td>
            <td style="text-align:center;">15</td>
            <td style="text-align:center;">✅</td></tr>
        <tr><td>④ 인구 밀집 이격도</td>
            <td style="text-align:center;">15</td>
            <td style="text-align:center;">✅</td></tr>
        <tr style="border-top:1px solid #eee;">
            <td rowspan="2">지질 안정성</td>
            <td>⑤ 자기 이상 균일도</td>
            <td style="text-align:center;">10</td>
            <td style="text-align:center;">⚠</td></tr>
        <tr><td style="color:#999;">⑥ 암상 적합성</td>
            <td style="text-align:center;color:#999;">5</td>
            <td style="text-align:center;color:#999;">-</td></tr>
        <tr style="border-top:1px solid #eee;">
            <td rowspan="2">운영 인프라</td>
            <td style="color:#888;">⑦ 부지 지속성</td>
            <td style="text-align:center;color:#888;">10</td>
            <td style="text-align:center;color:#888;">※</td></tr>
        <tr><td style="color:#888;">⑧ 관리 접근성</td>
            <td style="text-align:center;color:#888;">5</td>
            <td style="text-align:center;color:#888;">※</td></tr>
      </table>
      <hr style="margin:5px 0;border-color:#ccc;">
      <span style="color:#444;font-size:10.5px;">
      <b>산정 방식:</b><br>
      ① 후보점 간 평균 이격 거리 역산 (KDTree K=5)<br>
      ② 중심+4방향(5km) 표고 → 중앙차분 경사도(°) (Open-Elevation)<br>
      ③ log(전력·철도 이격 거리) 정규화<br>
      ④ log(도시·주거 이격 거리) 정규화<br>
      ⑤ KIGAM 반경 20km 내 자기이상 표준편차 역산<br>
      ⑦ 부지 지속성: 최종 선정 후 현장·지도 육안 확인<br>
      ⑧ 관리 접근성: 최종 선정 후 도로망 지도 확인<br>
      <br>
      가용 항목 합산 → 100점 정규화<br>
      </span>
      <hr style="margin:5px 0;border-color:#ccc;">
      <span style="color:#555;font-size:10.5px;">
      🔴 상위 34% → 1등급 최우선<br>
      🟠 34~67% → 2등급 우선<br>
      🟢 하위 33% → 3등급 일반<br>
      <span style="color:#999;">✅ 산정 가능 &nbsp; ⚠ KIGAM 필요 &nbsp; - 미확보 &nbsp; ※ 현장확인</span>
      </span>
    </div>
    """
    m.get_root().html.add_child(folium.Element(scoring_html))

    # ── 1:50,000 도엽 격자 레이어 (외부 GeoJSON fetch) ──────────
    topo_file = OUTPUT_DIR / data_subdir / "topo_sheets.geojson"
    if topo_file.exists():
        topo_layer = folium.FeatureGroup(
            name="📐 1:50,000 지형도 도엽 (15'×15')", show=False
        )
        topo_layer.add_to(m)
        tv  = topo_layer.get_name()
        # fetch GeoJSON, draw each polygon + zoom-responsive label
        topo_js = (
            "<script>(function(){"
            "fetch('%s/topo_sheets.geojson')"
            ".then(function(r){return r.json();})"
            ".then(function(fc){"
            "fc.features.forEach(function(f){"
            "var props=f.properties;"
            "var nm=props.sheet_name||'';"
            "var code=props.sheet_code||'';"
            "var mcd=props.sheet_mapidcd||'';"
            "var tipCode=mcd?mcd+'  /  '+code:code;"
            "var tip='<b>'+nm+'</b><br><span style=\"color:#555;font-size:11px;\">NGII No.&nbsp;'+tipCode+'</span>';"
            "var coords=f.geometry.coordinates[0].map(function(c){return[c[1],c[0]];});"
            "L.polygon(coords,{color:'#1A4A8A',weight:1.2,fill:true,"
            "fillColor:'#4A90D9',fillOpacity:0.06}).bindTooltip(tip,{sticky:true}).addTo(%s);"
            "var cx=props.cx||0,cy=props.cy||0;"
            "var sub=mcd?'<br><span style=\"font-weight:normal;font-size:0.85em;color:#336;\">'+mcd+'</span>':'';"
            "var lbl='<div class=\"topo-label\" style=\"display:none;pointer-events:none;"
            "text-align:center;font-family:Malgun Gothic,sans-serif;color:#1A3A6A;"
            "white-space:nowrap;font-weight:bold;line-height:1.3;"
            "text-shadow:1px 1px 0 white,-1px -1px 0 white,1px -1px 0 white,-1px 1px 0 white;\">'"
            "+nm+sub+'</div>';"
            "L.marker([cy,cx],{icon:L.divIcon({html:lbl,className:'',iconAnchor:[45,15],"
            "iconSize:[90,30]})}).addTo(%s);"
            "});"
            "});"
            "})();</script>"
        ) % (data_subdir, tv, tv)
        m.get_root().html.add_child(folium.Element(topo_js))

        # 줌 레벨별 도엽명 라벨 크기·표시 동적 조절
        zoom_js = """
        <script>
        (function() {
            function applyTopoLabels(zoom) {
                var labels = document.querySelectorAll('.topo-label');
                var display, size, weight;
                if (zoom <= 7) {
                    display = 'none';
                } else if (zoom === 8) {
                    display = 'block'; size = '7px'; weight = 'normal';
                } else if (zoom === 9) {
                    display = 'block'; size = '9px'; weight = 'bold';
                } else if (zoom === 10) {
                    display = 'block'; size = '12px'; weight = 'bold';
                } else if (zoom === 11) {
                    display = 'block'; size = '14px'; weight = 'bold';
                } else {
                    display = 'block'; size = '17px'; weight = 'bold';
                }
                labels.forEach(function(el) {
                    el.style.display = display;
                    if (display !== 'none') {
                        el.style.fontSize   = size;
                        el.style.fontWeight = weight;
                    }
                });
            }
            document.addEventListener('DOMContentLoaded', function() {
                var mapVarName = Object.keys(window).find(
                    function(k) { return /^map_[a-f0-9]+$/.test(k); }
                );
                if (!mapVarName) return;
                var leafletMap = window[mapVarName];
                setTimeout(function() { applyTopoLabels(leafletMap.getZoom()); }, 300);
                leafletMap.on('zoomend', function() {
                    applyTopoLabels(leafletMap.getZoom());
                });
                leafletMap.on('overlayadd overlayremove', function() {
                    setTimeout(function() { applyTopoLabels(leafletMap.getZoom()); }, 100);
                });
            });
        })();
        </script>
        """
        m.get_root().html.add_child(folium.Element(zoom_js))

    # ── 주소 검색 창 — 좌상단 고정 ──────────────────────────────
    # Nominatim OSM geocoder — 한국 도로명·지번·지명 검색
    # 3단계 폴백: ① countrycodes=kr ② +대한민국 ③ 상위 행정구역만
    geocoder_html = """
    <div id="geocoder-box" style="
        position:fixed; top:130px; left:10px; z-index:9999;
        background:rgba(255,255,255,0.97); border:2px solid #336699;
        border-radius:8px; padding:9px 11px;
        font-family:'Malgun Gothic',sans-serif; font-size:12px;
        box-shadow:2px 2px 8px rgba(0,0,0,0.25); width:278px;">
      <b style="font-size:12px;color:#224;">&#128269; 주소 / 지명 / 좌표 검색</b>
      <div style="display:flex;gap:4px;margin-top:6px;">
        <input id="gc-input" type="text"
          placeholder="예: 태백시 / 37.5665,126.9780"
          style="flex:1;padding:5px 7px;border:1px solid #aac;
                 border-radius:5px;font-size:11.5px;outline:none;
                 font-family:'Malgun Gothic',sans-serif;">
        <button id="gc-btn"
          style="padding:5px 10px;background:#336699;color:white;
                 border:none;border-radius:5px;cursor:pointer;
                 font-size:11.5px;white-space:nowrap;flex-shrink:0;">검색</button>
      </div>
      <div style="margin-top:3px;font-size:9.5px;color:#aaa;">
        위경도 직접 입력 가능: <i>37.5665, 126.9780</i> (위도, 경도)
      </div>
      <div id="gc-result" style="margin-top:5px;font-size:11px;color:#444;
           max-height:200px;overflow-y:auto;border-top:1px solid #eee;
           padding-top:3px;display:none;"></div>
    </div>

    <script>
    (function(){
      var gcMarker = null;
      var NOM = 'https://nominatim.openstreetmap.org/search';

      /* ── Folium 지도 객체 탐색 ─────────────────────────────── */
      function getMap() {
        var keys = Object.keys(window);
        for (var i = 0; i < keys.length; i++) {
          if (/^map_[a-f0-9]{8,}$/.test(keys[i]) && window[keys[i]] &&
              typeof window[keys[i]].setView === 'function') {
            return window[keys[i]];
          }
        }
        return null;
      }

      /* ── 결과 렌더링 ─────────────────────────────────────── */
      /* note: 상단에 표시할 안내 문자열 (optional) */
      function gcRender(items, query, note) {
        var div = document.getElementById('gc-result');
        if (!div) return;
        div.style.display = 'block';
        if (!items || items.length === 0) {
          div.innerHTML =
            '<div style="color:#c00;padding:4px 0;">'
            + '&#10060; <b>"' + query + '"</b> 결과 없음</div>'
            + '<div style="color:#888;font-size:10px;margin-top:3px;">'
            + '&#128161; 시&#183;군&#183;구 단위로 검색해 보세요<br>'
            + '예) <i>안동시</i>, <i>경상북도 안동</i>, <i>설악산</i></div>';
          return;
        }
        var typeMap = {
          administrative:'행정', road:'도로', amenity:'시설',
          place:'장소', natural:'자연', tourism:'관광',
          building:'건물', highway:'도로', residential:'주거',
          peak:'산', water:'수계', forest:'산림'
        };
        var html = note
          ? '<div style="color:#886;font-size:10px;padding:2px 0 4px;">'
            + '&#128161; ' + note + '</div>'
          : '';
        items.slice(0, 8).forEach(function(item, i) {
          var parts  = item.display_name.split(',');
          var label  = parts.slice(0, Math.min(3, parts.length)).join(', ').trim();
          var badge  = typeMap[item.type] || typeMap[item.class] || '';
          var lat    = parseFloat(item.lat);
          var lon    = parseFloat(item.lon);
          var safe   = item.display_name.replace(/\\\\/g,'\\\\\\\\').replace(/'/g,"\\\\'");
          html += '<div style="padding:5px 2px;border-bottom:1px solid #f0f0f0;cursor:pointer;" '
               +  'onmouseover="this.style.background=\'#eef3ff\'"'
               +  ' onmouseout="this.style.background=\'\'"'
               +  ' onclick="gcJump(' + lat + ',' + lon + ',\'' + safe + '\')">'
               +  '<span style="color:#336699;font-weight:bold;">' + (i+1) + '.</span> ' + label;
          if (badge) html += ' <span style="background:#ddeeff;color:#336;font-size:9px;'
                           + 'padding:1px 3px;border-radius:3px;">' + badge + '</span>';
          html += '<br><span style="color:#bbb;font-size:10px;">'
               +  lat.toFixed(5) + '&deg;N &nbsp;' + lon.toFixed(5) + '&deg;E</span></div>';
        });
        div.innerHTML = html;
      }

      /* ── fetch 래퍼 (에러 시 null 반환) ─────────────────────── */
      function nomFetch(url, cb) {
        fetch(url)
          .then(function(r) {
            if (!r.ok) { console.warn('[GC] HTTP', r.status, url); cb(null); return; }
            return r.json();
          })
          .then(function(d) { if (d !== undefined) cb(d || null); })
          .catch(function(e) { console.warn('[GC] fetch error:', e.message); cb(null); });
      }

      /* ── 좌표 파싱 헬퍼 ──────────────────────────────────────── */
      /* "33.524, 126.894" / "33.524480,126.894162" 등 다양한 구분자 지원 */
      function parseCoords(q) {
        /* 숫자 두 개를 쉼표·공백·세미콜론 중 하나로 구분 */
        var re = /^(-?[0-9]+[.]?[0-9]*)[,; ]+(-?[0-9]+[.]?[0-9]*)$/;
        var m = q.trim().match(re);
        if (!m) return null;
        var a = parseFloat(m[1]), b = parseFloat(m[2]);
        if (isNaN(a) || isNaN(b)) return null;
        /* 한반도 범위: 위도 33-38°, 경도 124-132° */
        if (a >= 33 && a <= 39 && b >= 124 && b <= 132) return {lat:a, lon:b};
        /* 입력이 경도·위도 순서로 뒤바뀐 경우 */
        if (b >= 33 && b <= 39 && a >= 124 && a <= 132) return {lat:b, lon:a};
        /* 한반도 범위 밖이지만 위경도 형식으로 판단되면 그대로 사용 */
        if (a >= -90 && a <= 90 && b >= -180 && b <= 180) return {lat:a, lon:b};
        return null;
      }

      /* ── 검색 메인 (좌표 우선 → 3단계 Nominatim 폴백) ───────── */
      window.gcSearch = function() {
        var inp = document.getElementById('gc-input');
        var div = document.getElementById('gc-result');
        if (!inp || !div) return;
        var q = inp.value.trim();
        if (!q) return;

        /* ① 좌표 직접 입력 감지 ───────────────────────────────── */
        var coords = parseCoords(q);
        if (coords) {
          div.style.display = 'block';
          div.innerHTML = '<div style="color:#226;padding:3px 0;">'
            + '&#128205; 좌표 이동: ' + coords.lat.toFixed(5) + '°N, '
            + coords.lon.toFixed(5) + '°E</div>';
          gcJump(coords.lat, coords.lon,
            coords.lat.toFixed(5) + ', ' + coords.lon.toFixed(5));
          return;
        }

        div.style.display = 'block';
        div.innerHTML = '<span style="color:#888;"><i>&#128269; 검색 중...</i></span>';
        console.log('[GC] query:', q);

        var base = NOM + '?format=json&limit=8&addressdetails=1&accept-language=ko,en';

        /* 시도1: countrycodes=kr, 원문 그대로 */
        var u1 = base + '&countrycodes=kr&q=' + encodeURIComponent(q);
        nomFetch(u1, function(d1) {
          console.log('[GC] try1 count:', d1 ? d1.length : 'null');
          if (d1 && d1.length > 0) { gcRender(d1, q); return; }

          /* 시도2: countrycodes 없이, "대한민국" 덧붙임 */
          var u2 = base + '&q=' + encodeURIComponent(q + ' 대한민국');
          nomFetch(u2, function(d2) {
            console.log('[GC] try2 count:', d2 ? d2.length : 'null');
            if (d2 && d2.length > 0) { gcRender(d2, q); return; }

            /* 시도3: 상위 행정구역만 (공백 기준 앞 2토큰) */
            var tokens = q.trim().split(' ').filter(function(t){return t.length>0;});
            var short  = tokens.slice(0, Math.min(2, tokens.length)).join(' ');
            if (short === q) { gcRender([], q); return; }   // 이미 단어 1개
            console.log('[GC] try3 simplified:', short);
            var u3 = base + '&countrycodes=kr&q=' + encodeURIComponent(short);
            nomFetch(u3, function(d3) {
              console.log('[GC] try3 count:', d3 ? d3.length : 'null');
              if (d3 && d3.length > 0) {
                gcRender(d3, short, '<i>"' + short + '"</i> 기준 결과 (주소 단순화)');
              } else {
                gcRender([], q);
              }
            });
          });
        });
      };

      /* ── 지도 이동 + 마커 ───────────────────────────────────── */
      window.gcJump = function(lat, lon, name) {
        var map = getMap();
        if (!map) { console.warn('[GC] map not found'); return; }
        map.setView([lat, lon], 13);
        if (gcMarker) { try { gcMarker.remove(); } catch(e){} gcMarker = null; }
        var label = name.split(',')[0].trim();
        gcMarker = L.marker([lat, lon], {
          icon: L.divIcon({
            html: '<div style="background:#336699;color:#fff;padding:3px 9px;'
                + 'border-radius:5px;font-size:12px;white-space:nowrap;'
                + 'font-family:Malgun Gothic,sans-serif;'
                + 'box-shadow:1px 2px 6px rgba(0,0,0,0.4);">&#128205; ' + label + '</div>',
            className: '', iconAnchor:[12, 32]
          })
        }).addTo(map);
        gcMarker.bindPopup(
          '<b>' + name.split(',').slice(0,3).join(', ') + '</b><br>'
          + '<span style="color:#666;font-size:11px;">'
          + lat.toFixed(5) + '&deg;N, ' + lon.toFixed(5) + '&deg;E</span>'
        ).openPopup();
        var div = document.getElementById('gc-result');
        if (div) div.innerHTML =
          '<div style="color:#226;padding:4px 0;">&#10004; <b>' + label + '</b> 이동 완료</div>';
      };

      /* ── Enter 키 + 버튼 바인딩 ────────────────────────────── */
      function bindUI() {
        var inp = document.getElementById('gc-input');
        var btn = document.getElementById('gc-btn');
        if (inp) inp.addEventListener('keydown', function(e) {
          if (e.key === 'Enter') { e.preventDefault(); window.gcSearch(); }
        });
        if (btn) btn.addEventListener('click', function() { window.gcSearch(); });
      }
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', bindUI);
      } else {
        bindUI();
      }
    })();
    </script>
    """
    m.get_root().html.add_child(folium.Element(geocoder_html))

    folium.LayerControl(collapsed=False).add_to(m)
    return m


# ============================================================
# 예상 소요 시간 추정
# ============================================================

def _fmt_time(secs: float) -> str:
    """초 → 사람이 읽기 쉬운 시간 문자열"""
    if secs < 60:
        return f"약 {int(secs)}초"
    elif secs < 3600:
        m, s = divmod(int(secs), 60)
        return f"약 {m}분 {s}초" if s else f"약 {m}분"
    else:
        h, rem = divmod(int(secs), 3600)
        return f"약 {h}시간 {rem // 60}분"


def estimate_runtime() -> None:
    """
    실행 전 단계별 예상 소요 시간 출력.

    추정 기준:
      - OSM 캐시 있음  : 항목당 ~1초 (파일 읽기)
      - OSM 캐시 없음  : 항목당 ~60초 (Overpass API, 네트워크 상황에 따라 변동)
      - DEM 캐시 있음  : 후보점당 ~0.05초
      - DEM 캐시 없음  : 후보점당 ~1.5초 (Open-Elevation API 5점 배치)
      - 점수 계산       : 후보점당 ~0.05초 (거리 행렬)
    """
    osm_items = {
        "전력(power_infra.json)":             DATA_DIR / "power_infra.json",
        "철도(railways.json)":               DATA_DIR / "railways.json",
        "도시핵심(urban_dense.json)":         DATA_DIR / "urban_dense.json",
        "도시주거(urban_residential_v2.json)": DATA_DIR / "urban_residential_v2.json",
        "파이프라인(pipelines.json)":          DATA_DIR / "pipelines.json",
        "군사시설(military.json)":             DATA_DIR / "military.json",
        "통신탑(comm_towers.json)":           DATA_DIR / "comm_towers.json",
        "풍력(wind_turbines.json)":           DATA_DIR / "wind_turbines.json",
        "채석장(quarries_mines.json)":        DATA_DIR / "quarries_mines.json",
        "수계(water_bodies.json)":            DATA_DIR / "water_bodies.json",
    }

    cached   = [k for k, p in osm_items.items() if p.exists()]
    uncached = [k for k, p in osm_items.items() if not p.exists()]

    osm_time = len(cached) * 1 + len(uncached) * 60  # seconds

    # 한반도 육상 면적 약 100,000 km², 격자 간격 기준 후보점 추산 (통과율 ~30%)
    grid_km      = GRID_SPACING_M / 1000
    raw_pts      = 100_000 / (grid_km ** 2)
    est_cands    = int(raw_pts * 0.30)

    dem_cache    = DATA_DIR / "dem_slopes.json"
    dem_per_pt   = 0.05 if dem_cache.exists() else 1.5
    dem_time     = est_cands * dem_per_pt

    score_time   = est_cands * 0.05
    total_time   = osm_time + dem_time + score_time + 30  # 기타 여유 30초

    print("\n━━━ 예상 소요 시간 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  격자 간격 : {grid_km:.0f} km  →  예상 후보점 : 약 {est_cands}개")
    print(f"  OSM 데이터: {len(cached)}개 캐시됨, {len(uncached)}개 신규 취득 필요")
    if uncached:
        print(f"    미캐시 항목 ({len(uncached)}개, 각 ~1분):")
        for item in uncached:
            print(f"      · {item}")
    print(f"  ① OSM 취득  : {_fmt_time(osm_time)}")
    dem_label = "캐시" if dem_cache.exists() else "API"
    print(f"  ② DEM 고도  : {_fmt_time(dem_time)}  ({dem_label})")
    print(f"  ③ 점수 계산 : {_fmt_time(score_time)}")
    print(f"  ──────────────────────────────────────────────────────────")
    print(f"  총 예상 시간: {_fmt_time(total_time)}")
    if len(uncached) == 0 and dem_cache.exists():
        print(f"  (모든 캐시 존재 — 빠른 실행 예상)")
    elif len(uncached) > 0:
        print(f"  ※ 네트워크 상태에 따라 OSM 취득 시간 변동 가능")
    print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")


# ============================================================
# 메인
# ============================================================

def main():
    t_total_start = time.time()

    print("=" * 66)
    print("  대한민국 지구자기장 모델 구축 — 측정 입지 선정 시스템")
    print("  Korea Geomagnetic Field Model — Site Selection Tool")
    print("=" * 66)

    estimate_runtime()

    # ── 1. 대한민국 경계 ─────────────────────────────────────
    korea = get_korea_boundary()

    # ── 2. OSM 데이터 취득 (캐싱) ────────────────────────────
    t_step = time.time()
    print("\n▶ 인공 간섭 요소 데이터 취득 (OSM Overpass API)")
    power_gdf       = get_power_infrastructure()  # [1] 송전탑/선
    railway_gdf     = get_railways()              # [2] 철도
    urban_dense_gdf = get_urban_dense()           # [3a] 핵심도심·산업
    urban_resid_gdf = get_urban_residential()     # [3b] 주거·취락
    pipeline_gdf    = get_pipelines()             # [4] 파이프라인
    comm_gdf        = get_comm_towers()           # [5] 통신탑·기지국
    wind_gdf        = get_wind_turbines()         # [6] 풍력발전기
    quarry_gdf      = get_quarries_mines()        # [7] 채석장·광산
    water_gdf       = get_water_bodies()          # [수계] 호수·저수지·강수면
    # ⑦ 국공유지 / ⑧ 일반도로 → 최종 선정 후 현장·지도 육안 확인 (계산 제외)
    print(f"  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 3. 자기이상도 — KIGAM 우선, EMAG2 폴백 ───────────────
    t_step = time.time()
    print("\n▶ 자기이상도 데이터 확인 [9]")
    mag_data    = load_kigam_anomaly()       # KIGAM 1.5분 격자 우선
    if mag_data is None:
        print("  ℹ  KIGAM 데이터 없음 — EMAG2 시도")
        mag_data = load_emag2_korea()        # 폴백: EMAG2 V3
    anomaly_gdf = None
    if mag_data is not None:
        anomaly_gdf = compute_anomaly_variation_zones(mag_data)
        print(f"  ✅ 자기이상도 변화폭 제외 구역 생성 완료 (기준: {ANOMALY_VARIATION_THRESHOLD} nT)")
    else:
        print(f"  ℹ  {DATA_DIR}/mag_1982-2018_1.5min_ed.dat 필요")
    print(f"  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 4. 후보 격자 생성 ─────────────────────────────────────
    t_step = time.time()
    print("\n▶ 후보 격자 생성")
    grid = create_candidate_grid(korea)
    print(f"  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 5. 제외 구역 구축 + 필터링 ───────────────────────────
    t_step = time.time()
    print("\n▶ 제외 구역 구축 및 후보점 필터링")
    zones = build_exclusion_zones(
        power_gdf, railway_gdf,
        urban_dense_gdf, urban_resid_gdf,
        pipeline_gdf, comm_gdf,
        wind_gdf, quarry_gdf, anomaly_gdf,
        water_gdf=water_gdf,
    )
    final_candidates = filter_candidates(grid, zones)
    print(f"  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 6. 입지 점수 산출 ──────────────────────────────────────
    t_step = time.time()
    print("\n▶ 입지 점수 산출 (다기준 평가)")
    final_candidates = compute_priority(final_candidates, zones, mag_data)
    print(f"  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 7. 지도 생성 ─────────────────────────────────────────
    t_step = time.time()
    print("\n▶ 기존 측정점 로드")
    existing_sites = load_existing_sites()

    print("\n▶ 지도 데이터 파일 저장 (output/data/)")
    data_dir = OUTPUT_DIR / "data"
    save_map_data(zones, grid, final_candidates, existing_sites, korea, data_dir)

    m = create_folium_map(
        zones, grid, final_candidates,
        korea_gdf=korea,
        existing_sites=existing_sites,
        data_subdir="data",
    )

    # ── 8. HTML 저장 ─────────────────────────────────────────
    html_path = OUTPUT_DIR / "geomag_site_selection.html"
    m.save(str(html_path))
    print(f"\n✅ 지도 저장: {html_path.resolve()}  [소요 {_fmt_time(time.time() - t_step)}]")

    # ── 9. CSV 저장 ──────────────────────────────────────────
    if len(final_candidates) > 0:
        result = final_candidates.to_crs(WGS84_CRS).copy()
        result["lat"]     = result.geometry.y.round(5)
        result["lon"]     = result.geometry.x.round(5)
        result["site_id"] = range(1, len(result) + 1)

        cols = ["site_id", "lat", "lon"]
        for c in ["priority", "score",
                  "s1_희소성", "s2_지형", "s3_전력철도", "s4_인구이격",
                  "s5_자기균일",
                  "dem_slope_deg", "d_power_km", "d_railway_km", "d_urban_km"]:
            if c in result.columns:
                cols.append(c)

        csv_path = OUTPUT_DIR / "candidate_sites.csv"
        result[cols].to_csv(str(csv_path), index=False, encoding="utf-8-sig")
        print(f"✅ 후보지 목록: {csv_path.resolve()}")

        # 우선순위별 요약
        print("\n━━━ 최종 결과 요약 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        print(f"  총 후보지: {len(final_candidates)}개 (격자 {GRID_SPACING_M//1000} km 기준)")
        if "priority" in result.columns:
            for p, label in [(1,"최우선"),(2,"우선"),(3,"일반")]:
                n = int((result["priority"]==p).sum())
                print(f"  우선순위 {p}등급 ({label}): {n}개")
        avail = int(result.get("score_max", pd.Series([55])).iloc[0])
        print(f"  입지 점수: 가용 {avail}/100점 → 100점 정규화")
        print(f"  제외 항목: {sum(1 for v in zones.values() if v is not None and not v.is_empty)}개 조건 적용")
        print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

    elapsed = time.time() - t_total_start
    print(f"\n✅ 완료!  총 소요 시간: {_fmt_time(elapsed)}")

    # ── 로컬 서버 실행 안내 ──────────────────────────────────
    print("\n" + "=" * 66)
    print("  ⚠  HTML 파일은 직접 열기(file://) 대신 로컬 HTTP 서버로 열어야")
    print("     fetch() 가 정상 동작합니다.")
    print()
    print("  방법 1 — Python 내장 서버 (권장):")
    print(f"    cd \"{OUTPUT_DIR.resolve()}\"")
    print(f"    python -m http.server 8000")
    print(f"    → 브라우저: http://localhost:8000/geomag_site_selection.html")
    print()
    print("  방법 2 — 스크립트 폴더에서 실행:")
    print(f"    python -m http.server 8000 --directory \"{OUTPUT_DIR.resolve()}\"")
    print("=" * 66)

    return final_candidates


if __name__ == "__main__":
    main()
