#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20260520_163338_P1_63sites_subgrid_기준점.kml 에
기존 측정점 (20260415_geomag_site_results.xlsx > '기존 측정점' 시트, 33개) 추가
"""
import sys, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")

from pathlib import Path
from datetime import datetime

OUT_DIR  = Path("C:/LG_gram_backup_users/LX/2026_geomag/docs/output")
KML_BASE = OUT_DIR / "20260520_163338_P1_63sites_subgrid_기준점.kml"
XLSX_SRC = OUT_DIR / "20260415_geomag_site_results.xlsx"
TS       = datetime.now().strftime("%Y%m%d_%H%M%S")
KML_OUT  = OUT_DIR / f"{TS}_P1_63sites_subgrid_기준점_측정점.kml"

print("=" * 62)
print("  기존 측정점 추가 → P1 63사이트 KML")
print("=" * 62)

# ── 1. 기존 측정점 로드 ────────────────────────────────────────
print(f"\n[1] {XLSX_SRC.name} '기존 측정점' 시트 로드...")
df = pd.read_excel(XLSX_SRC, sheet_name="기존 측정점")
print(f"  총 {len(df)}개 측정점")

pts = []
for _, row in df.iterrows():
    total_raw = row["총 자력 (nT)"]
    try:
        total_val = float(total_raw)
        total_str = f"{total_val:.1f} nT"
    except (ValueError, TypeError):
        total_val = None
        total_str = "미측정"

    pts.append({
        "sheet_name": str(row["도엽 이름"]).strip(),
        "name":       str(row["측정점 이름"]).strip(),
        "lat":        float(row["위도 (°N)"]),
        "lon":        float(row["경도 (°E)"]),
        "total":      total_val,
        "total_str":  total_str,
        "year":       int(row["관측 연도"]),
    })

recent = [p for p in pts if p["year"] >= 2022]
old    = [p for p in pts if p["year"] <  2022]
print(f"  최근 관측(2022~): {len(recent)}개 / 구 관측(~2021): {len(old)}개")

# ── 2. 스타일 ──────────────────────────────────────────────────
# KML color = AABBGGRR  (A=alpha, B=blue, G=green, R=red)
# 노란 별  ff00ffff  (2022~ 최근)
# 연회색 별 ffbbbbaa  (~2021 구)
STYLE_BLOCK = """
  <!-- ── 기존 측정점 스타일 ── -->
  <Style id="meas_recent_style">
    <IconStyle>
      <color>ff00ffff</color>
      <scale>1.1</scale>
      <Icon><href>http://maps.google.com/mapfiles/kml/shapes/star.png</href></Icon>
    </IconStyle>
    <LabelStyle><scale>0.85</scale><color>ff00cccc</color></LabelStyle>
  </Style>
  <Style id="meas_old_style">
    <IconStyle>
      <color>ffbbbbaa</color>
      <scale>0.9</scale>
      <Icon><href>http://maps.google.com/mapfiles/kml/shapes/star.png</href></Icon>
    </IconStyle>
    <LabelStyle><scale>0.75</scale><color>ffaaaaaa</color></LabelStyle>
  </Style>
"""

# ── 3. Placemark 생성 함수 ─────────────────────────────────────
def pm_block(p, style_id, indent="      "):
    nm        = p["name"]
    sheet_nm  = p["sheet_name"]
    lat       = p["lat"]
    lon       = p["lon"]
    yr        = p["year"]
    total_str = p["total_str"]
    naver_url = f"https://map.naver.com/v5/search/{lat},{lon}"
    kakao_url = f"https://map.kakao.com/link/map/{nm},{lat},{lon}"

    desc = (
        f"<b>📡 기존 지자기 측정점: {nm}</b><br>"
        f"도엽명: {sheet_nm}<br>"
        f"위도: {lat:.6f}°N | 경도: {lon:.6f}°E<br>"
        f"총 자력: <b>{total_str}</b><br>"
        f"관측 연도: <b>{yr}년</b><br><br>"
        f'<a href="{naver_url}">📍 네이버 지도</a>'
        f" &nbsp; "
        f'<a href="{kakao_url}">🗺 카카오맵</a>'
    )

    lines = [
        f"{indent}<Placemark>",
        f"{indent}  <name>★ {nm} ({yr})</name>",
        f"{indent}  <description><![CDATA[{desc}]]></description>",
        f"{indent}  <styleUrl>#{style_id}</styleUrl>",
        f"{indent}  <Point>",
        f"{indent}    <coordinates>{lon},{lat},0</coordinates>",
        f"{indent}  </Point>",
        f"{indent}</Placemark>",
    ]
    return "\n".join(lines)


# ── 4. Folder 블록 생성 ────────────────────────────────────────
def make_folder(pts_all):
    lines = []
    lines.append("  <Folder>")
    lines.append(f"    <name>★ 기존 지자기 측정점 ({len(pts_all)}개)</name>")
    lines.append(
        f"    <description><![CDATA["
        f"국토지리정보원 지자기 측정점 (관측 연도 2012~2025)<br>"
        f"<b>노란별(★)</b>: 2022년 이후 최근 관측 ({len(recent)}개)<br>"
        f"<b>회색별(★)</b>: 2021년 이전 구 관측 ({len(old)}개)"
        f"]]></description>"
    )
    lines.append("    <visibility>1</visibility>")

    # 최근 관측 서브폴더
    lines.append("    <Folder>")
    lines.append(f"      <name>최근 관측 (2022~, {len(recent)}개)</name>")
    lines.append("      <visibility>1</visibility>")
    for p in sorted(recent, key=lambda x: x["name"]):
        lines.append(pm_block(p, "meas_recent_style"))
    lines.append("    </Folder>")

    # 구 관측 서브폴더
    lines.append("    <Folder>")
    lines.append(f"      <name>구 관측 (~2021, {len(old)}개)</name>")
    lines.append("      <visibility>1</visibility>")
    for p in sorted(old, key=lambda x: x["name"]):
        lines.append(pm_block(p, "meas_old_style"))
    lines.append("    </Folder>")

    lines.append("  </Folder>")
    return "\n".join(lines)


# ── 5. KML 합체 ───────────────────────────────────────────────
print("\n[2] KML 로드 및 합체...")
base_kml = KML_BASE.read_text(encoding="utf-8")
print(f"  원본: {len(base_kml):,} bytes")

folder_block = make_folder(pts)

# 스타일 삽입 위치: 기존 "<!-- 스타일 정의 -->" 앞
# (없으면 첫 번째 <Style> 앞에 삽입)
if "<!-- 스타일 정의 -->" in base_kml:
    merged = base_kml.replace(
        "<!-- 스타일 정의 -->",
        "<!-- 스타일 정의 -->" + STYLE_BLOCK
    )
elif "<!-- ── 기준점 스타일" in base_kml:
    merged = base_kml.replace(
        "<!-- ── 기준점 스타일",
        STYLE_BLOCK.strip() + "\n  <!-- ── 기준점 스타일"
    )
else:
    # 첫 번째 <Style id= 앞에 삽입
    merged = base_kml.replace(
        "  <Style id=",
        STYLE_BLOCK + "  <Style id=",
        1  # 첫 번째 하나만
    )

# 기존 측정점 폴더: </Document> 직전
merged = merged.replace(
    "</Document>",
    folder_block + "\n\n</Document>"
)

# Document name 업데이트
old_name = "<name>1등급 P1 63개 후보지 — 서브격자 분석 + 기준점</name>"
new_name = "<name>1등급 P1 63개 후보지 — 서브격자 + 기준점 + 기존 측정점</name>"
merged = merged.replace(old_name, new_name)

old_desc = "<description>P1 우선순위 63개 사이트: 1km 서브격자 후보 · 삼각점·통합기준점 매칭</description>"
new_desc = (
    f"<description>P1 우선순위 63개 사이트: 서브격자 후보 · 기준점 매칭 · "
    f"기존 측정점 {len(pts)}개 (2012~2025)</description>"
)
merged = merged.replace(old_desc, new_desc)

KML_OUT.write_text(merged, encoding="utf-8")
sz = KML_OUT.stat().st_size / 1024

print(f"\n[3] 저장 완료")
print(f"  ✅ {KML_OUT.name}")
print(f"     크기: {sz:.1f} KB")

# ── 6. 요약 ───────────────────────────────────────────────────
print("\n" + "=" * 62)
print("  ✅ 완료!")
print()
print("  ── 포함 레이어 ────────────────────────────────────")
print(f"  📍 P1 후보지 63개  (1km 서브격자 후보 314개)")
print(f"  △/◎ 근접 기준점  733개 (삼각점·통합기준점)")
print(f"  ★  기존 측정점  {len(pts)}개 (2012~2025년 관측)")
print(f"      ├ 최근(2022~): {len(recent)}개  → 노란별")
print(f"      └ 구(~2021):  {len(old)}개  → 회색별")
print()
print("  ── 기존 측정점 목록 ─────────────────────────────────")
for p in sorted(pts, key=lambda x: x["name"]):
    print(f"    {p['name']:6s}  {p['lat']:.5f}N {p['lon']:.5f}E  "
          f"{p['year']}년  {p['total_str']}")
print("=" * 62)
