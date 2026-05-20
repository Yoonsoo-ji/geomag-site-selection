#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
site16_18.kml + 기준점.kml(전라북도) 합치기
전라북도 행정 bbox 안에 있는 기준점만 추출하여 하나의 KML로 생성
"""
import sys, re
sys.stdout.reconfigure(encoding="utf-8")

from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime

OUT_DIR = Path("C:/LG_gram_backup_users/LX/2026_geomag/docs/output")
KML_SITE = OUT_DIR / "20260519_113151_subgrid_site16_18.kml"
KML_REF  = OUT_DIR / "기준점.kml"
TS       = datetime.now().strftime("%Y%m%d_%H%M%S")
KML_OUT  = OUT_DIR / f"{TS}_site16_18_with_기준점_전북.kml"

# ── 전라북도 bbox (군산~무주 전체 포함) ──────────────────────
JB_SOUTH =  35.38
JB_NORTH =  36.16
JB_WEST  = 126.28
JB_EAST  = 127.90

print("=" * 60)
print("  site16_18.kml + 기준점.kml(전북) 합치기")
print(f"  전라북도 bbox: {JB_SOUTH}~{JB_NORTH}°N, {JB_WEST}~{JB_EAST}°E")
print("=" * 60)


# ════════════════════════════════════════════════════════════
# 1. 기준점.kml 파싱 — 전라북도 범위 Placemark 추출
#    구조: Folder/Folder(통합기준점|삼각점)/Placemark
# ════════════════════════════════════════════════════════════
print("\n[1] 기준점.kml 파싱 중 (18,854개 전국)...")

# ElementTree로 파싱 (namespace 처리)
ET.register_namespace("", "http://www.opengis.net/kml/2.2")
ET.register_namespace("gx", "http://www.google.com/kml/ext/2.2")
ET.register_namespace("atom", "http://www.w3.org/2005/Atom")

NS = "http://www.opengis.net/kml/2.2"

tree_ref = ET.parse(KML_REF)
root_ref = tree_ref.getroot()

def tag(name):
    return f"{{{NS}}}{name}"

jb_placemarks = []   # (folder_name, name, lon, lat, raw_element)

for folder_top in root_ref.iter(tag("Folder")):
    folder_nm_el = folder_top.find(tag("name"))
    folder_nm = folder_nm_el.text.strip() if folder_nm_el is not None else ""

    for subfolder in folder_top.findall(tag("Folder")):
        sf_nm_el = subfolder.find(tag("name"))
        sf_nm = sf_nm_el.text.strip() if sf_nm_el is not None else ""
        category = sf_nm  # "통합기준점" | "삼각점"

        for pm in subfolder.findall(tag("Placemark")):
            pt_el = pm.find(tag("Point"))
            if pt_el is None:
                continue
            coord_el = pt_el.find(tag("coordinates"))
            if coord_el is None or not coord_el.text:
                continue
            parts = coord_el.text.strip().split(",")
            if len(parts) < 2:
                continue
            try:
                lon, lat = float(parts[0]), float(parts[1])
            except ValueError:
                continue

            if not (JB_SOUTH <= lat <= JB_NORTH and JB_WEST <= lon <= JB_EAST):
                continue

            # 기준점 이름 추출 (ExtendedData > SchemaData > SimpleData)
            ctrl_nm = ""
            ed = pm.find(tag("ExtendedData"))
            if ed is not None:
                for sd in ed.iter(tag("SimpleData")):
                    if sd.get("name") == "CTRLPNT_NM":
                        ctrl_nm = (sd.text or "").strip()
                        break
            # name 태그가 있으면 우선 사용
            nm_el = pm.find(tag("name"))
            if nm_el is not None and nm_el.text:
                ctrl_nm = nm_el.text.strip() or ctrl_nm

            jb_placemarks.append({
                "category": category,
                "name":     ctrl_nm,
                "lon":      lon,
                "lat":      lat,
            })

cnt_tri  = sum(1 for p in jb_placemarks if p["category"] == "삼각점")
cnt_uni  = sum(1 for p in jb_placemarks if p["category"] == "통합기준점")
print(f"  전라북도 기준점: {len(jb_placemarks)}개")
print(f"    삼각점: {cnt_tri}개 / 통합기준점: {cnt_uni}개")


# ════════════════════════════════════════════════════════════
# 2. site16_18.kml 원본 읽기
# ════════════════════════════════════════════════════════════
print("\n[2] site16_18.kml 로드...")
site_kml = KML_SITE.read_text(encoding="utf-8")
print(f"  로드 완료 ({len(site_kml):,} bytes)")


# ════════════════════════════════════════════════════════════
# 3. 합본 KML 작성
# ════════════════════════════════════════════════════════════
print("\n[3] 합본 KML 생성...")

# 기준점 스타일 정의
style_block = """
  <!-- ── 기준점 스타일 ── -->
  <Style id="삼각점_style">
    <IconStyle>
      <color>ff00ff00</color>
      <scale>0.9</scale>
      <Icon><href>http://maps.google.com/mapfiles/kml/shapes/triangle.png</href></Icon>
    </IconStyle>
    <LabelStyle><scale>0.7</scale><color>ff00cc00</color></LabelStyle>
  </Style>
  <Style id="통합기준점_style">
    <IconStyle>
      <color>ffff8800</color>
      <scale>0.85</scale>
      <Icon><href>http://maps.google.com/mapfiles/kml/shapes/donut.png</href></Icon>
    </IconStyle>
    <LabelStyle><scale>0.7</scale><color>ffff8800</color></LabelStyle>
  </Style>
"""

# 기준점 Folder 블록 생성
def make_ref_folder(placemarks):
    lines = []
    lines.append("  <Folder>")
    lines.append(f"    <name>📐 전라북도 기준점 ({len(placemarks)}개)</name>")
    lines.append(f"    <description>삼각점 {cnt_tri}개 / 통합기준점 {cnt_uni}개</description>")
    lines.append("    <visibility>1</visibility>")

    # 삼각점 서브폴더
    tri_pts = [p for p in placemarks if p["category"] == "삼각점"]
    lines.append("    <Folder>")
    lines.append(f"      <name>삼각점 ({len(tri_pts)}개)</name>")
    lines.append("      <visibility>1</visibility>")
    for p in tri_pts:
        nm  = p["name"] or "삼각점"
        lon = p["lon"]; lat = p["lat"]
        naver_url = f"https://map.naver.com/v5/search/{lat},{lon}"
        lines.append("      <Placemark>")
        lines.append(f"        <name>△ {nm}</name>")
        lines.append(f"        <description><![CDATA["
                     f"<b>삼각점</b><br>"
                     f"명칭: {nm}<br>"
                     f"위도: {lat:.6f}°N<br>"
                     f"경도: {lon:.6f}°E<br><br>"
                     f'<a href="{naver_url}">📍 네이버 지도</a>'
                     f"]]></description>")
        lines.append("        <styleUrl>#삼각점_style</styleUrl>")
        lines.append("        <Point>")
        lines.append(f"          <coordinates>{lon},{lat},0</coordinates>")
        lines.append("        </Point>")
        lines.append("      </Placemark>")
    lines.append("    </Folder>")

    # 통합기준점 서브폴더
    uni_pts = [p for p in placemarks if p["category"] == "통합기준점"]
    lines.append("    <Folder>")
    lines.append(f"      <name>통합기준점 ({len(uni_pts)}개)</name>")
    lines.append("      <visibility>1</visibility>")
    for p in uni_pts:
        nm  = p["name"] or "통합기준점"
        lon = p["lon"]; lat = p["lat"]
        naver_url = f"https://map.naver.com/v5/search/{lat},{lon}"
        lines.append("      <Placemark>")
        lines.append(f"        <name>◎ {nm}</name>")
        lines.append(f"        <description><![CDATA["
                     f"<b>통합기준점</b><br>"
                     f"명칭: {nm}<br>"
                     f"위도: {lat:.6f}°N<br>"
                     f"경도: {lon:.6f}°E<br><br>"
                     f'<a href="{naver_url}">📍 네이버 지도</a>'
                     f"]]></description>")
        lines.append("        <styleUrl>#통합기준점_style</styleUrl>")
        lines.append("        <Point>")
        lines.append(f"          <coordinates>{lon},{lat},0</coordinates>")
        lines.append("        </Point>")
        lines.append("      </Placemark>")
    lines.append("    </Folder>")

    lines.append("  </Folder>")
    return "\n".join(lines)

ref_folder_block = make_ref_folder(jb_placemarks)

# site16_18.kml의 </Document> 직전에 스타일 + 기준점 폴더 삽입
merged = site_kml.replace(
    "  <!-- 스타일 정의 -->",
    style_block + "  <!-- ── 후보지 스타일 ── -->"
).replace(
    "</Document>",
    ref_folder_block + "\n\n</Document>"
)

# Document name 수정
merged = merged.replace(
    "<name>지구자기장 측정 후보지 서브격자 분석</name>",
    "<name>지구자기장 후보지 + 전북 기준점</name>"
).replace(
    "<description>후보지 #16·#18 1km 서브격자 현장 답사 세부 선정</description>",
    f"<description>후보지 #16·#18 서브격자 + 전라북도 기준점 {len(jb_placemarks)}개 (삼각점·통합기준점)</description>"
)

KML_OUT.write_text(merged, encoding="utf-8")
print(f"  ✅ 저장: {KML_OUT.name}")
print(f"     크기: {KML_OUT.stat().st_size/1024:.1f} KB")


# ════════════════════════════════════════════════════════════
# 4. 요약
# ════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("  ✅ 완료!")
print(f"  파일: {KML_OUT.name}")
print(f"  저장 위치: {OUT_DIR}")
print()
print("  ── 포함 레이어 ─────────────────────────────")
print(f"  📍 후보지 #16 서브격자   7개 (빨강별·주황)")
print(f"  📍 후보지 #18 서브격자  10개 (빨강별·주황)")
print(f"  △  삼각점 (전라북도)   {cnt_tri}개 (초록 삼각형)")
print(f"  ◎  통합기준점 (전라북도) {cnt_uni}개 (주황 도넛)")
print()
print("  KML 열기:")
print("  · Google Earth: 파일 → KML 열기")
print("  · 네이버 지도 앱: 내 장소 → KML 가져오기")
print("=" * 60)
