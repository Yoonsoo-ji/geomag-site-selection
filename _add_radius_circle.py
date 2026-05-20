#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P1 63개 사이트 중심에 5km 반경 원(Circle Polygon) 추가
입력: 20260520_172945_P1_63sites_filtered.kml
출력: {TS}_P1_63sites_with_radius.kml
"""
import sys, re, warnings
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime

MAIN_DIR = Path("C:/LG_gram_backup_users/LX/2026_geomag")
OUT_DIR  = MAIN_DIR / "docs" / "output"
KML_IN   = OUT_DIR / "20260520_174034_P1_63sites_filtered.kml"
TS       = datetime.now().strftime("%Y%m%d_%H%M%S")
KML_OUT  = OUT_DIR / f"{TS}_P1_63sites_with_radius.kml"

RADIUS_M = 5_000   # 5 km
N_PTS    = 90      # 원 꼭짓점 수 (90 = 4도 간격, 충분히 부드러움)
WGS84    = "EPSG:4326"
UTM      = "EPSG:5179"

NS   = "http://www.opengis.net/kml/2.2"
GX   = "http://www.google.com/kml/ext/2.2"
ATOM = "http://www.w3.org/2005/Atom"
def ktag(n): return f"{{{NS}}}{n}"

print("=" * 60)
print("  P1 사이트 5km 반경 원 추가")
print("=" * 60)

# ── 1. KML 파싱 ─────────────────────────────────────────────
print(f"\n[1] KML 로드...")
ET.register_namespace("", NS)
ET.register_namespace("gx", GX)
ET.register_namespace("atom", ATOM)

tree = ET.parse(KML_IN)
root = tree.getroot()

# ── 2. 반경 원 스타일 정의 ─────────────────────────────────
# KML color: AABBGGRR
# 원 외곽선: 빨간(ff0000ff), 채움: 10% 투명 빨간(1a0000ff)
STYLE_XML = """
  <!-- ── 5km 반경 원 스타일 ── -->
  <Style id="circle_5km_style">
    <LineStyle>
      <color>cc0000ff</color>
      <width>2</width>
    </LineStyle>
    <PolyStyle>
      <color>150000ff</color>
      <fill>1</fill>
      <outline>1</outline>
    </PolyStyle>
  </Style>
"""

# Document 첫 번째 Style 앞에 삽입 (string 방식)
# → 나중에 ET.tostring()에 포함시키기 위해 스타일 요소를 직접 트리에 추가
doc_el = root.find(ktag("Document"))
if doc_el is None:
    # <kml> 바로 아래에 Document 없으면 root 자체
    doc_el = root

# 스타일을 Document 첫 번째 자식으로 삽입
style_el = ET.fromstring(
    '<Style xmlns="http://www.opengis.net/kml/2.2" id="circle_5km_style">'
    '<LineStyle><color>cc0000ff</color><width>2</width></LineStyle>'
    '<PolyStyle><color>150000ff</color><fill>1</fill><outline>1</outline></PolyStyle>'
    '</Style>'
)
doc_el.insert(0, style_el)
print("  스타일 정의 추가 완료")

# ── 3. 원 좌표 생성 함수 ────────────────────────────────────
def make_circle_coords(clat, clon, radius_m=5000, n=90):
    """
    UTM에서 원 생성 후 WGS84로 역변환 → [(lon, lat), ...] 반환
    n+1개 (마지막은 첫점 반복으로 폴리곤 닫기)
    """
    center_gdf = gpd.GeoDataFrame(
        geometry=[Point(clon, clat)], crs=WGS84
    ).to_crs(UTM)
    cx, cy = center_gdf.geometry.iloc[0].x, center_gdf.geometry.iloc[0].y

    # UTM에서 원 → 다각형 (buffer)
    circle_utm = Point(cx, cy).buffer(radius_m, resolution=n // 4)

    # WGS84 역변환
    circle_gdf = gpd.GeoDataFrame(
        geometry=[circle_utm], crs=UTM
    ).to_crs(WGS84)
    coords = list(circle_gdf.geometry.iloc[0].exterior.coords)
    return coords  # [(lon, lat), ...]


def make_circle_pm(sid, clat, clon, radius_m=5000):
    """원 Polygon Placemark ET 요소 생성"""
    coords = make_circle_coords(clat, clon, radius_m)
    coords_txt = "\n".join(
        f"          {lon:.7f},{lat:.7f},0" for lon, lat in coords
    )

    pm = ET.Element(ktag("Placemark"))
    nm = ET.SubElement(pm, ktag("name"))
    nm.text = f"◯ #{sid} 반경 {radius_m//1000}km"
    desc = ET.SubElement(pm, ktag("description"))
    desc.text = (f"<![CDATA[후보지 <b>#{sid}</b> 서브격자 분석 반경<br>"
                 f"반경: <b>{radius_m//1000} km</b><br>"
                 f"중심: {clat:.5f}°N, {clon:.5f}°E]]>")
    vis = ET.SubElement(pm, ktag("visibility"))
    vis.text = "1"
    su = ET.SubElement(pm, ktag("styleUrl"))
    su.text = "#circle_5km_style"
    poly = ET.SubElement(pm, ktag("Polygon"))
    alt = ET.SubElement(poly, ktag("altitudeMode"))
    alt.text = "clampToGround"
    outer = ET.SubElement(poly, ktag("outerBoundaryIs"))
    ring = ET.SubElement(outer, ktag("LinearRing"))
    coords_el = ET.SubElement(ring, ktag("coordinates"))
    coords_el.text = "\n" + coords_txt + "\n        "
    return pm


# ── 4. 각 사이트 폴더에 원 삽입 ─────────────────────────────
print("\n[2] 사이트 폴더 탐색 및 원 생성...")

# 사이트 폴더 패턴: name이 "#숫자 (" 로 시작
site_re = re.compile(r'^#(\d+)\s+\((\d+\.\d+)°N[,\s]+(\d+\.\d+)°E\)')

added = 0
site_info = []   # (sid, clat, clon) 확인용

for folder_el in doc_el.iter(ktag("Folder")):
    nm_el = folder_el.find(ktag("name"))
    if nm_el is None or not nm_el.text:
        continue
    m = site_re.match(nm_el.text.strip())
    if not m:
        continue

    sid  = int(m.group(1))
    clat = float(m.group(2))
    clon = float(m.group(3))

    # 원 Placemark 생성
    circle_pm = make_circle_pm(sid, clat, clon, RADIUS_M)

    # 폴더의 첫 번째 자식으로 삽입 (중심 마커 앞)
    folder_el.insert(0, circle_pm)
    added += 1
    site_info.append((sid, clat, clon))
    print(f"  #{sid:3d}  {clat:.4f}°N {clon:.4f}°E  → 원 추가")

print(f"\n  ✅ 총 {added}개 사이트 원 추가 완료")

# ── 5. KML 저장 ─────────────────────────────────────────────
print(f"\n[3] KML 저장...")

# CDATA 처리: ET.tostring은 CDATA를 지원하지 않으므로
# description을 포함한 내용은 그냥 텍스트로 써도 구글어스에서 렌더링됨
kml_str = ET.tostring(root, encoding="unicode")

# 선언 추가
output = '<?xml version="1.0" encoding="UTF-8"?>\n' + kml_str
KML_OUT.write_text(output, encoding="utf-8")
sz = KML_OUT.stat().st_size / 1024

print(f"  ✅ {KML_OUT.name}")
print(f"     크기: {sz:.1f} KB")

print("\n" + "=" * 60)
print("  ✅ 완료!")
print()
print(f"  원 추가: {added}개 사이트 × 반경 {RADIUS_M//1000}km")
print(f"  스타일: 빨간 외곽선(2px) + 10% 투명 채움")
print(f"  출력: {KML_OUT.name}")
print("=" * 60)
