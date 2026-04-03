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

━━━ 지자기 이상 제외 기준 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [8] 자기이상도 gradient > 5 nT/km  (EMAG2 사용 시, 선택 / 제14조③)

━━━ 모델 구축 우선순위 가중치 ━━━━━━━━━━━━━━━━━━━━━━━━━━━
  - 기존 관측소와의 거리 (멀수록 데이터 공백 → 높은 우선순위)
  - 남북 위도 분포 균형 (제13조① 전국 균일 분포)

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

# 좌표 참조계: UTM Zone 52N (한반도 표준)
UTM_CRS   = "EPSG:32652"
WGS84_CRS = "EPSG:4326"

# Overpass API
OVERPASS_URL = "https://overpass-api.de/api/interpreter"

# 자기이상도 gradient 임계값 (nT/km) — EMAG2 사용 시 적용
ANOMALY_GRADIENT_THRESHOLD = 5.0   # nT/km  (조정 가능)

# ── 한국 지자기 관측소 (기준점) ───────────────────────────────
# 출처: KMA (기상청), INTERMAGNET
KOREA_OBSERVATORIES = [
    {
        "code": "SLE", "name": "서울",  "name_en": "Seoul",
        "lat": 37.160, "lon": 127.074, "org": "KMA",
        "established": 1957, "status": "운영중",
    },
    {
        "code": "ICH", "name": "인천",  "name_en": "Incheon",
        "lat": 37.271, "lon": 126.622, "org": "KMA",
        "established": 2006, "status": "운영중",
    },
    {
        "code": "USN", "name": "울산",  "name_en": "Ulsan",
        "lat": 35.527, "lon": 129.244, "org": "KMA",
        "established": 2014, "status": "운영중",
    },
]

# 모델 구축 시 신규 관측소 최소 이격 거리 (기존 관측소로부터)
OBS_MIN_SPACING_KM = 100.0

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


def compute_anomaly_gradient_zones(
    emag2_data: np.ndarray,
    threshold_nT_per_km: float = ANOMALY_GRADIENT_THRESHOLD,
) -> gpd.GeoDataFrame | None:
    """
    EMAG2 anomaly 수평 gradient 계산 →
    임계값 초과 지역을 제외 구역 GeoDataFrame으로 반환.

    Returns GeoDataFrame (WGS84) or None
    """
    if emag2_data is None or len(emag2_data) < 4:
        return None

    try:
        from scipy.interpolate import griddata
        from scipy.ndimage import sobel
    except ImportError:
        print("    scipy 없음 — 자기이상도 gradient 계산 생략")
        return None

    print("    자기이상도 gradient 계산 중...")

    # DataFrame 또는 numpy array 모두 지원
    if hasattr(emag2_data, "columns"):
        lons = emag2_data["lon"].values.astype(float)
        lats = emag2_data["lat"].values.astype(float)
        vals = emag2_data["anomaly_nT"].values.astype(float)
    else:
        lons = emag2_data[:, 0].astype(float)
        lats = emag2_data[:, 1].astype(float)
        vals = emag2_data[:, 2].astype(float)

    # 정규 격자 보간 (0.05도 간격)
    res = 0.05
    glon = np.arange(KOREA_BBOX[0], KOREA_BBOX[2] + res, res)
    glat = np.arange(KOREA_BBOX[1], KOREA_BBOX[3] + res, res)
    GLON, GLAT = np.meshgrid(glon, glat)

    grid_vals = griddata(
        (lons, lats), vals,
        (GLON, GLAT),
        method="linear",
        fill_value=np.nan,
    )

    # Sobel 필터로 gradient 추정 (nT / grid_cell)
    # 0.05도 ≈ 5.5 km → nT/grid_cell ÷ 5.5 = nT/km
    deg_km = 111.0 * res   # ~5.55 km per 0.05 deg
    gx = sobel(np.nan_to_num(grid_vals), axis=1) / (8 * deg_km)
    gy = sobel(np.nan_to_num(grid_vals), axis=0) / (8 * deg_km)
    grad = np.sqrt(gx**2 + gy**2)  # nT/km

    # 임계값 초과 → 제외 구역 폴리곤 생성
    from shapely.geometry import box
    high_grad_cells = []
    for i in range(len(glat)):
        for j in range(len(glon)):
            if grad[i, j] >= threshold_nT_per_km:
                lon0, lon1 = glon[j] - res / 2, glon[j] + res / 2
                lat0, lat1 = glat[i] - res / 2, glat[i] + res / 2
                high_grad_cells.append(box(lon0, lat0, lon1, lat1))

    if not high_grad_cells:
        return None

    gdf = gpd.GeoDataFrame(
        geometry=high_grad_cells, crs=WGS84_CRS
    )
    print(f"    고이상도 격자셀: {len(high_grad_cells)}개 (임계값 {threshold_nT_per_km} nT/km)")
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

            # CRS 처리:
            # 파일에 내장된 Korea 2000 UCS CRS는 EPSG 코드가 없지만 정확함.
            # → crs가 None일 때만 5179 강제 지정 (None이 아니면 내장 CRS 그대로 사용)
            if raw.crs is None:
                raw = raw.set_crs("EPSG:5179", allow_override=True)
            # WGS84로 변환 (내장 CRS 기반 정확 변환)
            topo = raw.to_crs(WGS84_CRS)

            # 컬럼 정리
            name_col = "MAPID_NM" if "MAPID_NM" in topo.columns else topo.columns[1]
            code_col = "MAPID_NO" if "MAPID_NO" in topo.columns else topo.columns[2]

            # 남한 범위 필터 (WGS84 변환 후 위도 33~38.65°)
            centroids_y = topo.geometry.centroid.y
            topo = topo[centroids_y.between(33.0, 38.65)].copy()

            # 도엽명 정리:
            # MAPID_NM이 순수 숫자(코드)인 경우 MAPID_NO 코드로 대체
            def _clean_name(row):
                nm = str(row[name_col]).strip()
                no = str(row[code_col]).strip()
                # 숫자만 있거나 비어있으면 도엽번호로 대체
                if not nm or nm.isdigit() or nm in ("nan", "None"):
                    return no
                return nm

            topo["sheet_name"] = topo.apply(_clean_name, axis=1)
            topo["sheet_code"] = topo[code_col].astype(str)

            # 한국 경계와 교차하는 도엽만
            topo = topo[topo.geometry.intersects(korea_geom)].reset_index(drop=True)
            result = topo[["sheet_name", "sheet_code", "geometry"]].copy()
            # 이름 샘플 출력 (디버깅용)
            sample_names = result["sheet_name"].head(5).tolist()
            print(f"  1:50,000 도엽 격자: {len(result)}개 셀 생성 (SHP) — 도엽명 예: {sample_names}")
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
) -> dict:
    """
    각 제외 조건별 단일 유니온 폴리곤 반환 (UTM CRS)
    키: power, railway, urban_dense, urban_resid, pipeline, comm, wind, quarry, anomaly
    도시 분류:
      - urban_dense: 상업/공업/건설 → URBAN_DENSE_BUFFER_M(500m) 버퍼
      - urban_resid: 주거/취락     → URBAN_RESIDENT_BUFFER_M(300m) 버퍼
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

    if anomaly_gdf is not None and len(anomaly_gdf) > 0:
        an_utm = anomaly_gdf.to_crs(UTM_CRS)
        zones["anomaly"] = unary_union(an_utm.geometry)
        print(f"    [9] 자기이상도: {zones['anomaly'].area/1e6:.0f} km²  (임계 {ANOMALY_GRADIENT_THRESHOLD} nT/km)")
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
        "anomaly":     "[9] 자기이상도 고변화",
    }
    for key, geom in zones.items():
        if geom is None or geom.is_empty:
            continue
        before = len(candidates)
        mask = ~candidates.geometry.intersects(geom)
        candidates = candidates[mask].reset_index(drop=True)
        removed = before - len(candidates)
        if removed > 0:
            print(f"  {labels[key]}: -{removed}개 → 잔여 {len(candidates)}개")

    print(f"\n  ✅ 최종 후보점: {len(candidates)}개 / 전체 격자 {len(grid)}개")
    return candidates


def fetch_dem_elevations(
    candidates_wgs: gpd.GeoDataFrame,
    radius_deg: float = 0.045,   # ~5 km at Korea latitudes
    cache_file: Path = DATA_DIR / "dem_elevations.json",
) -> np.ndarray:
    """
    Open-Elevation API를 통해 후보 지점의 지형적 대표성 산정용 표고 데이터 취득.

    각 후보점에 대해 중심 + 4방향(N/S/E/W, ~5 km) 총 5개 지점의 표고를 배치 조회.
    표고 표준편차(std) → 지형 기복 지표 (높을수록 지형 다양성 높음 → 대표성 높음).

    반환: shape (n,) numpy array — 각 후보점의 주변 5점 표고 표준편차 (m)
    """
    import json

    pts_wgs = candidates_wgs
    n = len(pts_wgs)

    # ── 캐시 확인 ─────────────────────────────────────────────
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                cached = json.load(f)
            if cached.get("n") == n and cached.get("radius_deg") == radius_deg:
                print(f"    DEM 캐시 로드: {cache_file.name} ({n}개 지점)")
                return np.array(cached["stds"])
        except Exception:
            pass

    print(f"    Open-Elevation API 조회 중 ({n}개 후보점 × 5방향)...")
    lats = pts_wgs.geometry.y.values
    lons = pts_wgs.geometry.x.values
    r = radius_deg

    # 각 후보점마다 5개 쿼리 포인트 생성: [center, N, S, E, W]
    all_locations = []
    for i in range(n):
        all_locations += [
            {"latitude": round(lats[i],          5), "longitude": round(lons[i],      5)},
            {"latitude": round(lats[i] + r,      5), "longitude": round(lons[i],      5)},
            {"latitude": round(lats[i] - r,      5), "longitude": round(lons[i],      5)},
            {"latitude": round(lats[i],           5), "longitude": round(lons[i] + r, 5)},
            {"latitude": round(lats[i],           5), "longitude": round(lons[i] - r, 5)},
        ]

    # ── 배치 분할 요청 (최대 1000포인트/요청) ────────────────────
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

    print(f"\n    DEM 취득 완료: {len(elevations)}개 표고값")

    # ── 후보점별 5포인트 표고 표준편차 계산 ──────────────────────
    stds = []
    elev_arr = np.array(elevations, dtype=float)
    for i in range(n):
        s = i * 5
        five = elev_arr[s : s + 5]
        if len(five) == 5 and not np.any(np.isnan(five)):
            stds.append(float(np.std(five)))
        else:
            stds.append(float(np.nanstd(five)) if len(five) > 0 else np.nan)
    stds = np.array(stds)

    # ── 캐시 저장 ─────────────────────────────────────────────
    try:
        with open(cache_file, "w") as f:
            json.dump({
                "n": n,
                "radius_deg": radius_deg,
                "stds": [float(v) for v in stds],
            }, f)
        print(f"    DEM 캐시 저장: {cache_file.name}")
    except Exception as exc:
        print(f"    ⚠ DEM 캐시 저장 실패: {exc}")

    return stds


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
    │             │ 지형적 대표성      │ 15 │ ⚠  Open-Elevation API        │
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

    # ── ② 지형적 대표성 (15점) — Open-Elevation API ──────────
    # 후보점 중심 + 4방향(~5km) 표고 표준편차 → 지형 기복 지표
    # std 높음 = 고저 변화 큼 = 지형 대표성 높음 (산지 등)
    # std 낮음 = 평탄지 — 측정은 유리하지만 지형 대표성 낮음
    s2 = np.full(n, np.nan)
    dem_available = False
    try:
        cands_wgs = candidates.to_crs(WGS84_CRS)
        dem_stds = fetch_dem_elevations(cands_wgs)
        valid_dem = ~np.isnan(dem_stds)
        if valid_dem.sum() > 0:
            mx_std = np.nanmax(dem_stds)
            if mx_std > 0:
                s2[valid_dem] = (dem_stds[valid_dem] / mx_std) * 15
            else:
                s2[valid_dem] = 7.5  # 모두 평탄지인 경우 중간값
            s2[~valid_dem] = 7.5    # 취득 실패 지점 중간값 대체
            dem_available = True
            result["dem_std_m"] = np.round(dem_stds, 1)
            print(f"    ② 지형적대표성: 평균 {np.nanmean(s2):.1f} / 15점  "
                  f"(표고 std 평균 {np.nanmean(dem_stds):.0f} m)")
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

    # ── ⑤ 자기 이상 균일도 (10점) ──────────────────────────
    # EMAG2 데이터 반경 0.2°(~20 km) 내 자기이상 표준편차 역산
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
            radius = 0.2
            stds = np.full(n, np.nan)
            for i, (clat, clon) in enumerate(zip(c_lats, c_lons)):
                mask = (np.abs(e_lats - clat) < radius) & (np.abs(e_lons - clon) < radius)
                if mask.sum() > 2:
                    stds[i] = np.std(e_anom[mask])
            valid = ~np.isnan(stds)
            if valid.sum() > 0:
                mx_s = np.nanmax(stds)
                if mx_s > 0:
                    s5[valid] = (1 - stds[valid] / mx_s) * 10
                else:
                    s5[valid] = 10
                s5[~valid] = 5
                emag_available = True
                print(f"    ⑤ 자기균일: 평균 {np.nanmean(s5):.1f} / 10점")
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


def create_folium_map(
    zones:        dict,
    grid:         gpd.GeoDataFrame,
    final_cands:  gpd.GeoDataFrame,   # UTM, priority 컬럼 포함
    korea_gdf:    gpd.GeoDataFrame | None = None,
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

    # ── 제외 구역 레이어 ─────────────────────────────────────
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
    for key, (fc, ec, name, show) in _excl.items():
        geom = zones.get(key)
        if geom is None or geom.is_empty:
            continue
        layer = folium.FeatureGroup(name=name, show=show)
        folium.GeoJson(
            _geom_to_wgs84(geom).__geo_interface__,
            style_function=lambda x, f=fc, e=ec: {
                "fillColor": f, "color": e, "weight": 1, "fillOpacity": 0.32,
            },
            tooltip=name,
        ).add_to(layer)
        layer.add_to(m)

    # ── 전체 격자점 (참고용) ──────────────────────────────────
    grid_layer = folium.FeatureGroup(name="격자점 전체 (참고)", show=False)
    grid_wgs   = grid.to_crs(WGS84_CRS)
    sample     = grid_wgs.sample(min(len(grid_wgs), 600), random_state=42)
    for _, row in sample.iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=2, color="#888888", fill=True,
            fill_color="#AAAAAA", fill_opacity=0.3, weight=0,
        ).add_to(grid_layer)
    grid_layer.add_to(m)

    # ── 최종 후보 지점 (우선순위별 3개 레이어) ──────────────────
    priority_cfg = {
        1: ("#FF0000", "#CC0000", "🔴 우선순위 1등급 (최우선, 모델 공백 지역)"),
        2: ("#FF8800", "#CC5500", "🟠 우선순위 2등급 (우선)"),
        3: ("#00BB00", "#008800", "🟢 우선순위 3등급 (일반)"),
    }

    final_wgs = final_cands.to_crs(WGS84_CRS)
    has_priority = "priority" in final_wgs.columns

    for p, (fc, ec, pname) in priority_cfg.items():
        if has_priority:
            subset = final_wgs[final_wgs["priority"] == p]
        else:
            subset = final_wgs if p == 3 else gpd.GeoDataFrame()

        if len(subset) == 0:
            continue

        p_layer = folium.FeatureGroup(name=pname, show=True)
        for idx, row in subset.iterrows():
            lat, lon = row.geometry.y, row.geometry.x
            score_str = f"{row.get('score', 0):.0f}" if has_priority else "-"
            v1 = f"{row.get('s1_희소성',   0):.1f}" if has_priority else "-"
            v2 = row.get('s2_지형', np.nan) if has_priority else np.nan
            v2s = f"{v2:.1f}" if not (isinstance(v2, float) and np.isnan(v2)) else "-"
            v3 = f"{row.get('s3_전력철도', 0):.1f}" if has_priority else "-"
            v4 = f"{row.get('s4_인구이격', 0):.1f}" if has_priority else "-"
            v5 = row.get('s5_자기균일', np.nan) if has_priority else np.nan
            v5s = f"{v5:.1f}" if not (isinstance(v5, float) and np.isnan(v5)) else "-"
            v7 = row.get('s7_부지', np.nan) if has_priority else np.nan
            v7s = f"{v7:.1f}" if not (isinstance(v7, float) and np.isnan(v7)) else "-"
            v8 = row.get('s8_접근성', np.nan) if has_priority else np.nan
            v8s = f"{v8:.1f}" if not (isinstance(v8, float) and np.isnan(v8)) else "-"
            dp = f"{row.get('d_power_km',  0):.1f}" if has_priority else "-"
            dr = f"{row.get('d_railway_km',0):.1f}" if has_priority else "-"
            du = f"{row.get('d_urban_km',  0):.1f}" if has_priority else "-"
            dpub = row.get('d_public_km', np.nan) if has_priority else np.nan
            dpub_str = (f" <span style='color:#888;'>({dpub:.1f}km)</span>"
                        if not (isinstance(dpub, float) and np.isnan(dpub)) else "")
            drd = row.get('d_road_km', np.nan) if has_priority else np.nan
            drd_str = (f" <span style='color:#888;'>({drd:.1f}km)</span>"
                       if not (isinstance(drd, float) and np.isnan(drd)) else "")
            dem_std = row.get('dem_std_m', np.nan) if has_priority else np.nan
            dem_str = (f" <span style='color:#888;'>({dem_std:.0f}m std)</span>"
                       if not (isinstance(dem_std, float) and np.isnan(dem_std)) else "")
            popup_html = (
                f'<div style="font-family:sans-serif;font-size:12.5px;min-width:250px;">'
                f'<b style="color:{ec};">측정 후보지 #{idx+1}</b>'
                f'&nbsp;<span style="background:{fc};color:white;padding:1px 5px;'
                f'border-radius:3px;font-size:11px;">'
                f'{"최우선" if p==1 else "우선" if p==2 else "일반"}</span><br>'
                f'<hr style="margin:4px 0;">'
                f'<b>위도:</b> {lat:.5f}° N &nbsp; '
                f'<b>경도:</b> {lon:.5f}° E<br>'
                f'<hr style="margin:4px 0;border-color:#ddd;">'
                f'<b>입지 점수: {score_str} / 100</b><br>'
                f'<span style="font-size:11.5px;color:#333;">'
                f'&nbsp;① 희소성: {v1} / 25<br>'
                f'&nbsp;② 지형적 대표성: {v2s} / 15{dem_str}<br>'
                f'&nbsp;③ 전력·철도 이격: {v3} / 15'
                f' <span style="color:#888;">(전력 {dp}km, 철도 {dr}km)</span><br>'
                f'&nbsp;④ 인구 이격: {v4} / 15'
                f' <span style="color:#888;">(도시 {du}km)</span><br>'
                f'&nbsp;⑤ 자기균일: {v5s} / 10<br>'
                f'&nbsp;⑦ 부지 지속성: {v7s} / 10{dpub_str}<br>'
                f'&nbsp;⑧ 관리 접근성: {v8s} / 5{drd_str}<br>'
                f'&nbsp;<span style="color:#999;">⑥암상: 미산정</span>'
                f'</span></div>'
            )
            folium.CircleMarker(
                location=[lat, lon],
                radius=7 if p == 1 else (6 if p == 2 else 5),
                color=ec,
                fill=True,
                fill_color=fc,
                fill_opacity=0.85,
                weight=1.5,
                popup=folium.Popup(popup_html, max_width=270),
                tooltip=f"후보지 #{idx+1}  P{p}  ({lat:.4f}°N, {lon:.4f}°E)",
            ).add_to(p_layer)
        p_layer.add_to(m)

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
      {_swatch('#0044FF','[8] 자기이상도 고변화 *')}
      <hr style="margin:6px 0;border-color:#ccc;">
      <b>▸ 측정 후보지 (총 {n_cands}개)</b><br>
      {_dot('#FF0000',f'1등급 최우선 (데이터 공백): {n_p1}개')}
      {_dot('#FF8800',f'2등급 우선: {n_p2}개')}
      {_dot('#00BB00',f'3등급 일반: {n_p3}개')}
      <hr style="margin:6px 0;border-color:#ccc;">
      <small style="color:#555;">
        격자 간격: {GRID_SPACING_M//1000} km | 좌표계: WGS84/UTM 52N<br>
        데이터: OpenStreetMap (Overpass API)<br>
        * EMAG2 배치 시 활성화
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
        <tr><td>② 지형적 대표성</td>
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
            <td>⑦ 부지 지속성</td>
            <td style="text-align:center;">10</td>
            <td style="text-align:center;">✅</td></tr>
        <tr><td>⑧ 관리 접근성</td>
            <td style="text-align:center;">5</td>
            <td style="text-align:center;">✅</td></tr>
      </table>
      <hr style="margin:5px 0;border-color:#ccc;">
      <span style="color:#444;font-size:10.5px;">
      <b>산정 방식:</b><br>
      ① 후보점 간 평균 이격 거리 역산 (KDTree K=5)<br>
      ② 중심+4방향(5km) 표고 표준편차 (Open-Elevation)<br>
      ③ log(전력·철도 이격 거리) 정규화<br>
      ④ log(도시·주거 이격 거리) 정규화<br>
      ⑤ KIGAM 반경 20km 내 자기이상 표준편차 역산<br>
      ⑦ 국공유지 거리: 내부=10, 1km=7, 2km=4, 외부=1<br>
      ⑧ 일반도로 거리: ≤500m=5, ≤1km=3, >1km=0<br>
      <br>
      가용 항목 합산 → 100점 정규화<br>
      </span>
      <hr style="margin:5px 0;border-color:#ccc;">
      <span style="color:#555;font-size:10.5px;">
      🔴 상위 34% → 1등급 최우선<br>
      🟠 34~67% → 2등급 우선<br>
      🟢 하위 33% → 3등급 일반<br>
      <span style="color:#999;">✅ 산정 가능 &nbsp; ⚠ KIGAM/EMAG2 필요 &nbsp; - 데이터 미확보</span>
      </span>
    </div>
    """
    m.get_root().html.add_child(folium.Element(scoring_html))

    # ── 1:50,000 도엽 격자 레이어 ────────────────────────────
    if korea_gdf is not None:
        try:
            print("  1:50,000 도엽 격자 생성 중...")
            topo_gdf = create_topo_sheet_grid(korea_gdf)
            topo_layer = folium.FeatureGroup(
                name="📐 1:50,000 지형도 도엽 (15'×15')", show=False
            )
            for _, row in topo_gdf.iterrows():
                geom = row.geometry
                coords = list(geom.exterior.coords)
                # 폴리곤 경계선
                folium.Polygon(
                    locations=[(lat, lon) for lon, lat in coords],
                    color="#1A4A8A",
                    weight=1.2,
                    fill=True,
                    fill_color="#4A90D9",
                    fill_opacity=0.06,
                    tooltip=(
                        f"<b>{row['sheet_name']}</b><br>"
                        f"<span style='color:#666;font-size:11px;'>{row['sheet_code']}</span>"
                    ),
                ).add_to(topo_layer)
                # 셀 중심에 도엽명 라벨 (줌 반응형 — CSS 클래스로 제어)
                cx = geom.centroid.x
                cy = geom.centroid.y
                name_escaped = row["sheet_name"].replace("'", "\\'").replace('"', "&quot;")
                folium.Marker(
                    location=[cy, cx],
                    icon=folium.DivIcon(
                        html=(
                            f'<div class="topo-label"'
                            f' data-name="{name_escaped}"'
                            f' style="display:none;pointer-events:none;'
                            f'font-family:\'Malgun Gothic\',sans-serif;'
                            f'color:#1A3A6A;white-space:nowrap;font-weight:bold;'
                            f'text-shadow:1px 1px 0 white,-1px -1px 0 white,'
                            f'1px -1px 0 white,-1px 1px 0 white,'
                            f'0 0 4px white,0 0 4px white;">'
                            f'{row["sheet_name"]}</div>'
                        ),
                        icon_size=(80, 18),
                        icon_anchor=(40, 9),
                    ),
                ).add_to(topo_layer)
            topo_layer.add_to(m)

            # 줌 레벨별 도엽명 라벨 크기·표시 동적 조절
            zoom_js = """
            <script>
            (function() {
                // 줌 레벨별 라벨 스타일 설정
                // zoom  7 이하: 숨김
                // zoom  8~9  : 7px — 전국 뷰
                // zoom 10    : 10px
                // zoom 11    : 13px
                // zoom 12+   : 16px — 시/군 뷰
                function applyTopoLabels(zoom) {
                    var labels = document.querySelectorAll('.topo-label');
                    var display, size, weight, shadow;
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
                    // 초기 적용
                    setTimeout(function() {
                        applyTopoLabels(leafletMap.getZoom());
                    }, 300);
                    // 줌 변경 시 적용
                    leafletMap.on('zoomend', function() {
                        applyTopoLabels(leafletMap.getZoom());
                    });
                    // 레이어 토글 시에도 적용
                    leafletMap.on('overlayadd overlayremove', function() {
                        setTimeout(function() {
                            applyTopoLabels(leafletMap.getZoom());
                        }, 100);
                    });
                });
            })();
            </script>
            """
            m.get_root().html.add_child(folium.Element(zoom_js))

        except Exception as exc:
            print(f"  ⚠ 도엽 격자 생성 실패: {exc}")

    # ── 주소 검색 창 — 좌상단 고정 ──────────────────────────────
    # Nominatim OSM geocoder + V-World 좌표 검색 통합
    # 한국 도로명·지번·지명 모두 지원
    geocoder_html = """
    <div id="geocoder-box" style="
        position:fixed; top:55px; left:10px; z-index:9999;
        background:rgba(255,255,255,0.97); border:2px solid #336699;
        border-radius:8px; padding:9px 11px;
        font-family:'Malgun Gothic',sans-serif; font-size:12px;
        box-shadow:2px 2px 8px rgba(0,0,0,0.25); width:272px;">
      <b style="font-size:12px;color:#224;">🔍 주소 / 지명 검색</b>
      <span style="font-size:10px;color:#888;margin-left:3px;">도로명·지번·지명</span>
      <div style="display:flex;gap:4px;margin-top:6px;">
        <input id="gc-input" type="text"
          placeholder="예: 태백시, 홍성군, 설악산..."
          style="flex:1;padding:5px 7px;border:1px solid #aac;
                 border-radius:5px;font-size:11.5px;outline:none;
                 font-family:'Malgun Gothic',sans-serif;"
          onkeydown="if(event.key==='Enter'){event.preventDefault();gcSearch();}">
        <button onclick="gcSearch()"
          style="padding:5px 10px;background:#336699;color:white;
                 border:none;border-radius:5px;cursor:pointer;
                 font-size:11.5px;white-space:nowrap;flex-shrink:0;">검색</button>
      </div>
      <div style="margin-top:3px;font-size:9.5px;color:#aaa;">
        예) 세종대로 110 &nbsp;|&nbsp; 북한산 &nbsp;|&nbsp; 충청남도 홍성군
      </div>
      <div id="gc-result" style="margin-top:5px;font-size:11px;color:#444;
           max-height:180px;overflow-y:auto;border-top:1px solid #eee;
           padding-top:3px;display:none;"></div>
    </div>

    <script>
    (function(){
      var gcMarker = null;

      // 지도 객체 가져오기
      function getMap() {
        var k = Object.keys(window).find(function(x){return /^map_[a-f0-9]+$/.test(x);});
        return k ? window[k] : null;
      }

      // 결과 렌더링
      function gcRender(items) {
        var div = document.getElementById('gc-result');
        div.style.display = 'block';
        if (!items || items.length === 0) {
          div.innerHTML = '<span style="color:#c00;">검색 결과 없음 — 지명이나 행정구역으로 검색해 보세요</span>';
          return;
        }
        var typeMap = {
          'administrative':'행정','road':'도로','amenity':'시설',
          'place':'장소','natural':'자연','tourism':'관광',
          'building':'건물','highway':'도로','residential':'주거',
          'peak':'산','water':'수계','forest':'산림'
        };
        var html = '';
        items.slice(0,8).forEach(function(item, i){
          var parts = item.display_name.split(',');
          var label = parts.slice(0, Math.min(3, parts.length)).join(', ').trim();
          var badge = typeMap[item.type] || typeMap[item.class] || '';
          var lat = parseFloat(item.lat), lon = parseFloat(item.lon);
          // escape quotes in name for onclick
          var safeName = item.display_name.replace(/\\/g,'\\\\').replace(/'/g,"\\'");
          html += '<div style="padding:4px 2px;border-bottom:1px solid #f0f0f0;cursor:pointer;" '
               + 'onmouseover="this.style.background=\'#eef3ff\'" onmouseout="this.style.background=\'\'" '
               + 'onclick="gcJump(' + lat + ',' + lon + ',\'' + safeName + '\')">'
               + '<span style="color:#336699;font-weight:bold;">' + (i+1) + '.</span> ' + label;
          if(badge) html += ' <span style="background:#ddeeff;color:#336;font-size:9px;'
                          + 'padding:1px 3px;border-radius:3px;">' + badge + '</span>';
          html += '<br><span style="color:#bbb;font-size:10px;">'
               + lat.toFixed(5) + '°N &nbsp;' + lon.toFixed(5) + '°E</span></div>';
        });
        div.innerHTML = html;
      }

      // 검색 실행 (Nominatim API — 한국 한정, 뷰박스 설정)
      window.gcSearch = function() {
        var q = document.getElementById('gc-input').value.trim();
        if (!q) return;
        var div = document.getElementById('gc-result');
        div.style.display = 'block';
        div.innerHTML = '<span style="color:#888;"><i>검색 중...</i></span>';

        // viewbox: 한국 전체 (서→동, 남→북)
        var vb = '124.5,33.0,130.0,38.9';
        var url = 'https://nominatim.openstreetmap.org/search?format=json&limit=8'
                + '&countrycodes=kr&viewbox=' + vb + '&bounded=0'
                + '&accept-language=ko,en'
                + '&q=' + encodeURIComponent(q);

        fetch(url, {method:'GET', mode:'cors'})
          .then(function(r){
            if(!r.ok) throw new Error('HTTP ' + r.status);
            return r.json();
          })
          .then(function(data){
            if(!data || data.length === 0){
              // 결과 없으면 '대한민국' 추가 재시도
              var url2 = 'https://nominatim.openstreetmap.org/search?format=json&limit=8'
                       + '&accept-language=ko,en'
                       + '&q=' + encodeURIComponent(q + ' 대한민국');
              return fetch(url2, {method:'GET',mode:'cors'}).then(function(r2){return r2.json();});
            }
            return data;
          })
          .then(function(data){ gcRender(data); })
          .catch(function(err){
            div.innerHTML = '<span style="color:#c00;">API 오류: ' + err.message
              + '<br><small>인터넷 연결을 확인하세요</small></span>';
          });
      };

      // 지도 이동 + 마커
      window.gcJump = function(lat, lon, name) {
        var map = getMap();
        if (!map) return;
        map.setView([lat, lon], 13);
        if (gcMarker) { gcMarker.remove(); gcMarker = null; }
        var label = name.split(',')[0].trim();
        gcMarker = L.marker([lat, lon], {
          icon: L.divIcon({
            html: '<div style="background:#336699;color:#fff;padding:3px 8px;'
                + 'border-radius:5px;font-size:12px;white-space:nowrap;'
                + 'font-family:Malgun Gothic,sans-serif;'
                + 'box-shadow:1px 2px 6px rgba(0,0,0,0.4);">📍 ' + label + '</div>',
            className: '', iconAnchor:[10,30]
          })
        }).addTo(map);
        gcMarker.bindPopup(
          '<b>' + name.split(',').slice(0,3).join(', ') + '</b><br>'
          + '<span style="color:#666;font-size:11px;">'
          + lat.toFixed(5) + '°N, ' + lon.toFixed(5) + '°E</span>'
        ).openPopup();
        document.getElementById('gc-result').innerHTML =
          '<div style="color:#226;padding:3px 0;">✔ <b>' + label + '</b> 이동 완료</div>';
      };

      // Enter 키 이벤트 (DOMContentLoaded 후 바인딩)
      document.addEventListener('DOMContentLoaded', function(){
        var inp = document.getElementById('gc-input');
        if(inp) inp.addEventListener('keydown', function(e){
          if(e.key === 'Enter'){ e.preventDefault(); gcSearch(); }
        });
      });
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
        "전력(power_infrastructure.json)":   DATA_DIR / "power_infrastructure.json",
        "철도(railways.json)":               DATA_DIR / "railways.json",
        "도시핵심(urban_dense.json)":         DATA_DIR / "urban_dense.json",
        "도시주거(urban_residential_v2.json)": DATA_DIR / "urban_residential_v2.json",
        "파이프라인(pipelines.json)":          DATA_DIR / "pipelines.json",
        "통신탑(comm_towers.json)":           DATA_DIR / "comm_towers.json",
        "풍력(wind_turbines.json)":           DATA_DIR / "wind_turbines.json",
        "채석장(quarries.json)":              DATA_DIR / "quarries.json",
    }

    cached   = [k for k, p in osm_items.items() if p.exists()]
    uncached = [k for k, p in osm_items.items() if not p.exists()]

    osm_time = len(cached) * 1 + len(uncached) * 60  # seconds

    # 한반도 육상 면적 약 100,000 km², 격자 간격 기준 후보점 추산 (통과율 ~30%)
    grid_km      = GRID_SPACING_M / 1000
    raw_pts      = 100_000 / (grid_km ** 2)
    est_cands    = int(raw_pts * 0.30)

    dem_cache    = DATA_DIR / "dem_elevations.json"
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
    # ⑦ 국공유지 / ⑧ 일반도로 → 최종 선정 후 육안 확인 (계산 제외)
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
        anomaly_gdf = compute_anomaly_gradient_zones(mag_data)
        print(f"  ✅ 자기이상도 gradient 제외 구역 생성 완료")
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
    m = create_folium_map(zones, grid, final_candidates, korea_gdf=korea)

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
                  "dem_std_m", "d_power_km", "d_railway_km", "d_urban_km"]:
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
    return final_candidates


if __name__ == "__main__":
    main()
