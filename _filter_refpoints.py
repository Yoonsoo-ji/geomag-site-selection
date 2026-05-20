#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KML 내 기준점 중 제외구역에 포함된 것을 제거 (지질·암상 필터 포함)
- 제외구역: power(1km), railway(5km), urban_dense(500m), urban_resid(300m),
            pipeline(500m), comm(500m), wind(500m), quarry(1km), anomaly(WKB),
            자성암종 폴리곤(직접 제외), 단층선(500m)
- 입력 : 20260520_170218_P1_63sites_subgrid_기준점_측정점.kml  (원본, 733개)
- 출력 : {TS}_P1_63sites_filtered.kml
"""
import sys, json, warnings, re
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

import geopandas as gpd
import shapely
from shapely.geometry import Point, LineString, Polygon
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime
from collections import defaultdict

MAIN_DIR  = Path("C:/LG_gram_backup_users/LX/2026_geomag")
DATA_DIR  = MAIN_DIR / "data"
CACHE_DIR = DATA_DIR / "zone_cache"
GEO_DIR   = DATA_DIR / "수치지질도_25만축척_전국"
OUT_DIR   = MAIN_DIR / "docs" / "output"

KML_IN  = OUT_DIR / "20260520_170218_P1_63sites_subgrid_기준점_측정점.kml"
TS      = datetime.now().strftime("%Y%m%d_%H%M%S")
KML_OUT = OUT_DIR / f"{TS}_P1_63sites_filtered.kml"

WGS84 = "EPSG:4326"
UTM   = "EPSG:5179"

# ── 버퍼 거리 (m) ────────────────────────────────────────────
BUFFERS = {
    "power":       1_000,   # 고압철탑·송전탑
    "railway":     5_000,   # 철도
    "urban_dense":   500,   # 핵심도심·산업
    "urban_resid":   300,   # 주거·취락
    "pipeline":      500,   # 파이프라인
    "comm":          500,   # 통신탑·기지국
    "wind":          500,   # 풍력발전기
    "quarry":      1_000,   # 채석장·광산
    "anomaly":         0,   # 자기이상 WKB (이미 버퍼 포함)
    "geology":         0,   # 자성암종 폴리곤 — 내부 직접 제외
    "fault":         500,   # 단층선
}

# ── JSON 파일 매핑 ────────────────────────────────────────────
JSON_FILES = {
    "railway":     DATA_DIR / "railways.json",
    "urban_dense": DATA_DIR / "urban_dense.json",
    "urban_resid": DATA_DIR / "urban_residential_v2.json",
    "pipeline":    DATA_DIR / "pipelines.json",
    "comm":        DATA_DIR / "comm_towers.json",
    "wind":        DATA_DIR / "wind_turbines.json",
    "quarry":      DATA_DIR / "quarries_mines.json",
}

# ── 자성 암종 코드 (geomag_site_selection.py 동일 기준) ──────
MAGNETIC_ROCK_CODES = {
    # 현무암 (Quaternary basalt)
    "Qb","Qb(I)","Qb(I)(T)","Qb(II)","Qb(II)(S)","Qb(III)","Qb(III)(S)",
    "Qtb","Qtb(I)","Qtb(I)(T)","Qtb(II)","Qtb(II)(S)","Qtb(III)","Qtb(III)(S)",
    "Qtb(IV)","Qtb(IV)(S)","Qtb(V)","Qtb(V)(S)","Qtb(VI)","Qtb(VI)(S)",
    "Qtb(VII)","Qtb(VII)(S)","Qtb(VIII)","Qtb(VIII)(S)","Qtb(VIII)(T)",
    "Qtb(VIII)((T)","Qtt",
    # 신생대 화산암 (Tertiary)
    "Tb","Tbv","Tav","Teo","Te1","Te2","Te3",
    # 쥐라기 현무암
    "Jb",
    # 반려암 (Gabbro)
    "Jgb","Jga","Kgb",
    # 휘록암·반려암질 (Diabase/Doleritic)
    "Jda","Jdi","Kdi","Kdd",
    # 안산암류 (Andesite)
    "Kan","Kav",
    # 백악기 화산암류 (Cretaceous volcanic)
    "Kv","Kkv","Kjv","Kgsv","Khsv","Kpsv","Kssv","Ktsv","Kchv","Kcv",
    # 미분류 염기성암·화산암
    "bg","v",
    # 각섬암 (Amphibolite)
    "am",
}

NS   = "http://www.opengis.net/kml/2.2"
GX   = "http://www.google.com/kml/ext/2.2"
ATOM = "http://www.w3.org/2005/Atom"

print("=" * 66)
print("  기준점 제외구역 필터링 (지질·암상 포함)")
print("=" * 66)


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

ref_folders = []
for folder_el in root.iter(ktag("Folder")):
    nm_el = folder_el.find(ktag("name"))
    nm = (nm_el.text or "") if nm_el is not None else ""
    if "근접 기준점" in nm:
        ref_folders.append(folder_el)
print(f"  기준점 폴더 수: {len(ref_folders)}개")

pm_records = []
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
        pm_records.append({"folder_el": folder_el, "pm_el": pm_el,
                            "lon": lon, "lat": lat})

total_pm = len(pm_records)
print(f"  기준점 Placemark 총수: {total_pm}개")

pts_gdf = gpd.GeoDataFrame(
    {"idx": range(total_pm)},
    geometry=[Point(r["lon"], r["lat"]) for r in pm_records],
    crs=WGS84
).to_crs(UTM)

excluded = set()
excl_reason = {}   # idx → 제외 이유


# ════════════════════════════════════════════════════════════════
# 2. 제외구역 로드 및 교차 검사
# ════════════════════════════════════════════════════════════════
print("\n[2] 제외구역 로드 및 교차 검사...")

def check_and_record(name, new_excl):
    """제외 집합 누적 + 이유 기록"""
    fresh = new_excl - excluded   # 이번에 새로 추가되는 것만
    for i in fresh:
        excl_reason[i] = name
    excluded.update(new_excl)


def check_zone_wkb(name, wkb_path):
    if not wkb_path.exists():
        print(f"  {name:15s}: WKB 없음 → 건너뜀"); return set()
    zone_geom = shapely.from_wkb(wkb_path.read_bytes())
    zone_utm  = gpd.GeoDataFrame(geometry=[zone_geom], crs=UTM)
    try:
        joined = gpd.sjoin(pts_gdf[["idx","geometry"]], zone_utm,
                           how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  {name:15s}: WKB          → 제외 {len(excl):3d}개")
        return excl
    except Exception as e:
        print(f"  {name}: WKB 오류: {e}"); return set()


def json_to_gdf(json_path):
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
    return gpd.GeoDataFrame(geometry=geoms or [], crs=WGS84)


def check_zone_json(name, json_path, buffer_m):
    if not json_path.exists():
        print(f"  {name:15s}: JSON 없음 → 건너뜀"); return set()
    try:
        zone_gdf = json_to_gdf(json_path)
        if len(zone_gdf) == 0:
            print(f"  {name:15s}: 피처 없음 → 건너뜀"); return set()
        zone_utm = zone_gdf.to_crs(UTM)
        pts_buf  = pts_gdf[["idx","geometry"]].copy()
        pts_buf["geometry"] = pts_buf.geometry.buffer(buffer_m)
        joined = gpd.sjoin(pts_buf, zone_utm[["geometry"]],
                           how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  {name:15s}: {len(zone_gdf):8,}개 피처 → 제외 {len(excl):3d}개")
        return excl
    except Exception as e:
        print(f"  {name}: 오류: {e}"); return set()


def check_zone_geology(shp_path):
    """자성 암종 폴리곤 내부 기준점 제외 (lithoidx 기준)"""
    if not shp_path.exists():
        print(f"  geology        : SHP 없음 → 건너뜀"); return set()
    try:
        for enc in ("utf-8", "cp949", "euc-kr"):
            try:
                litho = gpd.read_file(str(shp_path), encoding=enc)
                break
            except: continue
        else:
            print("  geology        : 인코딩 오류"); return set()

        litho = litho[litho.geometry.notna()].copy()
        if litho.crs is None: litho = litho.set_crs(WGS84)
        elif str(litho.crs) != WGS84: litho = litho.to_crs(WGS84)

        # 자성 암종만 필터
        if "lithoidx" not in litho.columns:
            # 컬럼명 확인
            print(f"  geology        : 컬럼 = {litho.columns.tolist()[:8]}")
            return set()
        mag = litho[litho["lithoidx"].isin(MAGNETIC_ROCK_CODES)].copy()
        print(f"  geology        : 전체 {len(litho):,}개 → 자성암종 {len(mag):,}개")
        if len(mag) == 0:
            print(f"  geology        : 자성암종 없음 → 건너뜀"); return set()

        mag_utm = mag.to_crs(UTM)
        # 폴리곤 내부 포함 여부 (buffer=0 → 경계 포함)
        joined = gpd.sjoin(pts_gdf[["idx","geometry"]], mag_utm[["geometry"]],
                           how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  geology        : 자성암종 내부 → 제외 {len(excl):3d}개")
        return excl
    except Exception as e:
        print(f"  geology        : 오류: {e}"); return set()


def check_zone_fault(shp_path, buffer_m=500):
    """단층선 buffer_m 이내 기준점 제외"""
    if not shp_path.exists():
        print(f"  fault          : SHP 없음 → 건너뜀"); return set()
    try:
        for enc in ("cp949", "euc-kr", "utf-8"):
            try:
                fault = gpd.read_file(str(shp_path), encoding=enc)
                break
            except: continue
        else:
            print("  fault          : 인코딩 오류"); return set()

        fault = fault[fault.geometry.notna()].copy()
        if fault.crs is None: fault = fault.set_crs(WGS84)
        elif str(fault.crs) != WGS84: fault = fault.to_crs(WGS84)

        fault_utm = fault.to_crs(UTM)
        # 기준점을 buffer_m 팽창 → 단층선과 교차 여부
        pts_buf = pts_gdf[["idx","geometry"]].copy()
        pts_buf["geometry"] = pts_buf.geometry.buffer(buffer_m)
        joined = gpd.sjoin(pts_buf, fault_utm[["geometry"]],
                           how="left", predicate="intersects")
        excl = set(joined[joined.index_right.notna()]["idx"].values)
        print(f"  fault          : {len(fault):,}개 단층선 ({buffer_m}m 버퍼) → 제외 {len(excl):3d}개")
        return excl
    except Exception as e:
        print(f"  fault          : 오류: {e}"); return set()


# ── WKB 캐시 ───────────────────────────────────────────────
check_and_record("power",   check_zone_wkb("power",   CACHE_DIR / "power.wkb"))
check_and_record("anomaly", check_zone_wkb("anomaly", CACHE_DIR / "anomaly.wkb"))

# ── JSON 존 ────────────────────────────────────────────────
for zone_name, json_path in JSON_FILES.items():
    check_and_record(zone_name, check_zone_json(zone_name, json_path, BUFFERS[zone_name]))

# ── 지질·암상 (SHP) ────────────────────────────────────────
check_and_record("geology", check_zone_geology(GEO_DIR / "Geology_250K_Litho.shp"))
check_and_record("fault",   check_zone_fault(GEO_DIR / "Geology_250K_Fault.shp", BUFFERS["fault"]))

print(f"\n  ▶ 총 제외 대상: {len(excluded)} / {total_pm}개 기준점")
print(f"  ▶ 잔존 기준점 : {total_pm - len(excluded)}개")


# ════════════════════════════════════════════════════════════════
# 3. 제외 Placemark 제거 + 폴더명 업데이트
# ════════════════════════════════════════════════════════════════
print("\n[3] KML 수정 중...")

folder_pms = defaultdict(list)
for i, rec in enumerate(pm_records):
    folder_pms[id(rec["folder_el"])].append(i)

removed_total = 0
for folder_el in ref_folders:
    fid    = id(folder_el)
    pm_idxs = folder_pms[fid]
    to_remove = [i for i in pm_idxs if i in excluded]
    for i in to_remove:
        try:
            folder_el.remove(pm_records[i]["pm_el"])
            removed_total += 1
        except ValueError:
            pass
    remaining = len(pm_idxs) - len(to_remove)
    nm_el = folder_el.find(ktag("name"))
    if nm_el is not None and nm_el.text:
        nm_el.text = re.sub(r"\(\d+개\)", f"({remaining}개)", nm_el.text)

print(f"  제거 완료: {removed_total}개")


# ════════════════════════════════════════════════════════════════
# 4. KML 저장
# ════════════════════════════════════════════════════════════════
print("\n[4] KML 저장...")
kml_str = ET.tostring(root, encoding="unicode")
output  = '<?xml version="1.0" encoding="UTF-8"?>\n' + kml_str
KML_OUT.write_text(output, encoding="utf-8")
sz = KML_OUT.stat().st_size / 1024
print(f"  ✅ {KML_OUT.name}")
print(f"     크기: {sz:.1f} KB")


# ════════════════════════════════════════════════════════════════
# 5. 요약
# ════════════════════════════════════════════════════════════════
# 제외 이유별 집계
from collections import Counter
reason_cnt = Counter(excl_reason.values())

print("\n" + "=" * 66)
print("  ✅ 완료!")
print()
print("  ── 기준점 필터링 결과 ─────────────────────────────────")
print(f"  원본 기준점     : {total_pm}개")
print(f"  제외 합계       : {len(excluded)}개")
print(f"  최종 잔존       : {total_pm - len(excluded)}개")
print()
print("  ── 제외구역별 건수 (중복 제거 후 첫 번째 이유 기준) ───")
zone_labels = {
    "power":       f"고압철탑·송전탑     ({BUFFERS['power']//1000}km)",
    "anomaly":      "자기이상 >800nT     (WKB)",
    "railway":     f"철도                ({BUFFERS['railway']//1000}km)",
    "urban_dense": f"핵심도심·산업       ({BUFFERS['urban_dense']}m)",
    "urban_resid": f"주거·취락           ({BUFFERS['urban_resid']}m)",
    "pipeline":    f"파이프라인          ({BUFFERS['pipeline']}m)",
    "comm":        f"통신탑·기지국       ({BUFFERS['comm']}m)",
    "wind":        f"풍력발전기          ({BUFFERS['wind']}m)",
    "quarry":      f"채석장·광산         ({BUFFERS['quarry']//1000}km)",
    "geology":      "자성암종 폴리곤     (직접 제외)",
    "fault":       f"단층선              ({BUFFERS['fault']}m)",
}
for k, label in zone_labels.items():
    cnt = reason_cnt.get(k, 0)
    if cnt > 0:
        print(f"    {label}: {cnt}개")
print("=" * 66)
