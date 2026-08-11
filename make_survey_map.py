# -*- coding: utf-8 -*-
"""
현장조사 선점 검토 결과 → 지도(Folium) 마커 페이지
===================================================

aggregate_survey_xlsx 의 파서·검토 로직을 재사용해 103개 후보지를 등급별
(A 선점가능 / B 조건부 / C 부적합·재검토 / 미완료) 마커로 지도에 표출한다.

  · 타일: OSM + Esri 위성 (기존 geomag 지도와 동일)
  · 레이어: 등급별 FeatureGroup 토글, 범례, 팝업(검토 결론·조사자 의견)
  · 출력: docs/survey_review.html  (GitHub Pages 에서 /survey_review.html 로 표출)
    index.html 은 건드리지 않는 독립 페이지.

사용:
    python make_survey_map.py
    python make_survey_map.py --dir <조사서폴더>
"""
import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import folium

from aggregate_survey_xlsx import parse_workbook, review, key_disturb, survey_files

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs" / "survey_review.html"
PHOTO_DIR = ROOT / "docs" / "survey_photos"

# 기존 지자기 측정점(24) — 현황표 1차 비고 초록(15) ∪ 2차 비고 파랑(9). (사용자 확정)
EXISTING_USE = [
    # 초록(15)
    "춘양", "청송", "서산", "청양", "이원", "남지", "거제", "장흥",
    "성산", "여주", "화천", "봉평", "삼척", "제천", "와도",
    # 파랑(9)
    "상주", "미원", "부안", "강화", "함양", "가야", "양산", "순천", "포천",
]
EXISTING_ALIAS = {}   # 확정 이름이 existing_sites.geojson 과 직접 일치

GRADE = {   # 등급 → (색, 라벨)
    "A": ("#2E8B57", "A · 선점 가능 (자기구배 조사)"),
    "B": ("#E0A400", "B · 조건부 (현장 재확인)"),
    "C": ("#CC3333", "C · 부적합 (대체 후보지 검토)"),
    "미완료": ("#888888", "미완료 · 재조사 필요"),
}


def fnum(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            if s is not None else "")


def az_text(d):
    """방위표지 1·2 방위각·거리 요약."""
    det = d.get("방위표지상세", {})
    out = []
    for i, k in enumerate(("표지1", "표지2"), 1):
        c = det.get(k, {})
        az, dist = c.get("방위각", ""), c.get("거리", "")
        if az and az != "-":
            seg = f"표지{i} 방위각 {esc(az)}"
            if dist and dist != "-":
                seg += f" · {esc(dist)}m"
            out.append(seg)
    return " / ".join(out)


def photo_html(d):
    ph = d.get("사진", {})
    if not ph:
        return ""
    cells = ""
    for slot in ("중심", "동", "서", "남", "북"):
        f = ph.get(slot)
        if not f:
            continue
        url = "survey_photos/" + f
        cells += (
            f"<a href='{url}' target='_blank' title='{slot}측 (클릭 확대)' "
            f"style='text-decoration:none;display:inline-block;text-align:center;margin:2px'>"
            f"<img src='{url}' loading='lazy' style='width:66px;height:66px;"
            f"object-fit:cover;border:1px solid #ccc;border-radius:3px'/>"
            f"<div style='font-size:10px;color:#666'>{slot}</div></a>")
    if not cells:
        return ""
    return ("<div style='margin-top:7px'>"
            "<div style='color:#666;font-size:11px;margin-bottom:2px'>주변 사진 (클릭 확대)</div>"
            f"<div style='display:flex;flex-wrap:wrap'>{cells}</div></div>")


def _fmt_az(s):
    nums = re.findall(r"\d+\.?\d*", str(s))
    if len(nums) >= 3:
        return f"{int(float(nums[0]))}°{int(float(nums[1]))}′{float(nums[2]):.0f}″"
    return str(s) if s else "-"


def svg_sketch(d):
    """기준점+방위표지1·2 를 실제 방위·거리로 배치한 약도(SVG). 좌표 없으면 안내문."""
    import math
    base, m1, m2 = d.get("기준점ll"), d.get("표지1ll"), d.get("표지2ll")
    if not (base and m1 and m2):
        return ("<div style='margin-top:7px;color:#999;font-size:11px'>"
                "방위표지 좌표 미확정 — 약도 없음</div>")
    det = d.get("방위표지상세", {})
    W, H = 296, 236
    cx, cy = W / 2, H / 2 + 6
    lat0 = base[1]

    def en(ll):
        de = (ll[0] - base[0]) * 111320 * math.cos(math.radians(lat0))
        dn = (ll[1] - base[1]) * 111320
        return de, dn
    e1, e2 = en(m1), en(m2)
    maxr = max(math.hypot(*e1), math.hypot(*e2), 1.0)
    scale = (min(W, H) / 2 - 46) / maxr

    def xy(e):
        return cx + e[0] * scale, cy - e[1] * scale
    x1, y1 = xy(e1)
    x2, y2 = xy(e2)

    def mark(x, y, color, glyph):
        return (f"<circle cx='{x:.1f}' cy='{y:.1f}' r='7' fill='{color}' "
                f"stroke='#fff' stroke-width='1.5'/>"
                f"<text x='{x:.1f}' y='{y+3.5:.1f}' font-size='9' fill='#fff' "
                f"text-anchor='middle' font-weight='bold'>{glyph}</text>")

    def label(x, y, txt, color, dy=-12):
        return (f"<text x='{x:.1f}' y='{y+dy:.1f}' font-size='10' fill='{color}' "
                f"text-anchor='middle' font-weight='bold' "
                f"style='paint-order:stroke;stroke:#fff;stroke-width:3px'>{txt}</text>")

    def leg(tag, xm, ym):   # 방위각·거리 중점 라벨
        c = det.get(tag, {})
        az = _fmt_az(c.get("방위각", ""))
        dist = c.get("거리", "")
        t = f"{az}" + (f" · {dist}m" if dist and dist != "-" else "")
        mx, my = (cx + xm) / 2, (cy + ym) / 2
        return (f"<text x='{mx:.1f}' y='{my:.1f}' font-size='8.5' fill='#222' "
                f"text-anchor='middle' style='paint-order:stroke;stroke:#fff;"
                f"stroke-width:3px'>{t}</text>")

    parts = [f"<svg viewBox='0 0 {W} {H}' width='100%' style='max-width:{W}px;"
             "background:#f7f8fa;border:1px solid #ddd;border-radius:5px'>"]
    # 북 화살표
    parts.append(f"<line x1='{W-20}' y1='30' x2='{W-20}' y2='14' stroke='#444' "
                 "stroke-width='1.6' marker-end='url(#arr)'/>"
                 "<defs><marker id='arr' markerWidth='7' markerHeight='7' refX='3' refY='3' "
                 "orient='auto'><path d='M0,0 L6,3 L0,6 Z' fill='#444'/></marker></defs>"
                 f"<text x='{W-20}' y='42' font-size='10' fill='#444' text-anchor='middle' "
                 "font-weight='bold'>N</text>")
    # 라인
    parts.append(f"<line x1='{cx}' y1='{cy}' x2='{x1:.1f}' y2='{y1:.1f}' stroke='#1D4ED8' stroke-width='2'/>")
    parts.append(f"<line x1='{cx}' y1='{cy}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='#0E8A6B' stroke-width='2'/>")
    # 중점 방위각·거리
    parts.append(leg("표지1", x1, y1))
    parts.append(leg("표지2", x2, y2))
    # 마커 + 이름
    parts.append(mark(x1, y1, "#1D4ED8", "1"))
    parts.append(mark(x2, y2, "#0E8A6B", "2"))
    parts.append(mark(cx, cy, "#E8531F", "★"))
    parts.append(label(x1, y1, "방위표지 1", "#1D4ED8"))
    parts.append(label(x2, y2, "방위표지 2", "#0E8A6B"))
    parts.append(label(cx, cy, "기준점", "#E8531F"))
    parts.append("</svg>")
    return ("<div style='margin-top:8px'>"
            "<div style='color:#666;font-size:11px;margin-bottom:3px'>방위표지 약도 (진북 기준)</div>"
            + "".join(parts) + "</div>")


def popup_html(d, grade, concl, note):
    color = GRADE[grade][0]
    dist = key_disturb(d) or "없음"
    bang = d["방위표지"] or "-"
    az = az_text(d)
    rows = [
        ("종합 판정", esc(d["종합판정"])),
        ("핵심 교란요인", esc(dist)),
        ("방위표지", esc(bang) + (f" <span style='color:#888'>({az})</span>" if az else "")),
        ("검토 결론", f"<b>{esc(concl)}</b>"),
        ("조사자 의견", esc(note)),
        ("조사", f"{esc(d['조사일'])} · {esc(d['조사자'])}"),
        ("좌표", f"{esc(d['위도'])}, {esc(d['경도'])} (표고 {esc(d['표고'])} m)"),
        ("연계 기존점", esc(d["연계기존점"])),
    ]
    body = "".join(
        f"<tr><td style='color:#666;padding:2px 8px 2px 0;white-space:nowrap;"
        f"vertical-align:top'>{k}</td><td style='padding:2px 0'>{v}</td></tr>"
        for k, v in rows)
    return (
        f"<div style='font-family:\"맑은 고딕\",sans-serif;font-size:12.5px;"
        f"width:330px;line-height:1.45'>"
        f"<div style='background:{color};color:#fff;padding:6px 10px;margin:-2px -2px 6px;"
        f"border-radius:4px 4px 0 0;font-weight:bold'>"
        f"[{grade}] {esc(d['관리번호'])} · {esc(d['후보지명'])}</div>"
        f"<table style='border-collapse:collapse'>{body}</table>"
        f"{svg_sketch(d)}"
        f"{photo_html(d)}"
        f"<div style='margin-top:6px;color:#999;font-size:11px'>{esc(d['관할본부'])}</div>"
        f"</div>")


def legend_html(counts, total):
    items = "".join(
        f"<div style='display:flex;align-items:center;margin:3px 0'>"
        f"<span style='width:13px;height:13px;border-radius:50%;background:{c};"
        f"display:inline-block;margin-right:7px;border:1px solid #fff;"
        f"box-shadow:0 0 0 1px #999'></span>"
        f"<span style='font-size:12px'>{lab} "
        f"<b>{counts.get(g,0)}</b></span></div>"
        for g, (c, lab) in GRADE.items())
    return (
        "<div style='position:fixed;bottom:22px;left:22px;z-index:9999;"
        "background:rgba(255,255,255,.95);padding:11px 14px;border-radius:8px;"
        "box-shadow:0 1px 6px rgba(0,0,0,.3);font-family:\"맑은 고딕\",sans-serif'>"
        "<div style='font-weight:bold;font-size:13px;margin-bottom:6px'>"
        f"선점 검토 등급 <span style='color:#888;font-weight:normal'>(총 {total}점)</span></div>"
        f"{items}</div>")


def add_topo_layer(m):
    """도엽(1:50,000) 경계 폴리곤 토글 — index.html 과 동일 스타일, 기본 꺼짐."""
    p = DATA / "topo_sheets.geojson"
    if not p.exists():
        return 0
    fc = json.load(open(p, encoding="utf-8"))
    fg = folium.FeatureGroup(name="🗺️ 도엽 경계 (1:50,000)", show=False)
    for f in fc["features"]:
        pr = f["properties"]
        nm, code = pr.get("sheet_name", ""), pr.get("sheet_code", "")
        mcd = pr.get("sheet_mapidcd", "")
        tip = (f"<b>{nm}</b><br><span style='color:#555;font-size:11px;'>"
               f"NGII No.&nbsp;{(mcd + '  /  ' + code) if mcd else code}</span>")
        coords = [[c[1], c[0]] for c in f["geometry"]["coordinates"][0]]
        folium.Polygon(coords, color="#1A4A8A", weight=1.2, fill=True,
                       fill_color="#4A90D9", fill_opacity=0.06,
                       tooltip=folium.Tooltip(tip, sticky=True)).add_to(fg)
    fg.add_to(m)
    return len(fc["features"])


def add_existing_layer(m):
    """기존 지자기 측정점(EXISTING_USE) 토글 — ⭐ 마커, 기본 꺼짐."""
    p = DATA / "existing_sites.geojson"
    if not p.exists():
        return 0
    fc = json.load(open(p, encoding="utf-8"))
    want = {EXISTING_ALIAS.get(n, n) for n in EXISTING_USE}
    by_name = {f["properties"].get("name"): f for f in fc["features"]}
    fg = folium.FeatureGroup(name=f"⭐ 기존 측정점 ({len(EXISTING_USE)})", show=False)
    n = 0
    for name in EXISTING_USE:
        gname = EXISTING_ALIAS.get(name, name)
        f = by_name.get(gname)
        if not f:
            print(f"    ! 기존측정점 미매칭: {name}", file=sys.stderr)
            continue
        pr = f["properties"]
        lat, lon = pr["lat"], pr["lon"]
        vf = lambda x, u="": f"{x:.4f}{u}" if isinstance(x, (int, float)) else "-"
        pop = (
            f"<div style='font-family:\"맑은 고딕\",sans-serif;font-size:12.5px;min-width:250px'>"
            f"<b style='color:#8B4513'>⭐ 기존 측정점: {esc(name)}</b><hr style='margin:4px 0'>"
            f"<b>위도:</b> {lat:.6f}° N &nbsp; <b>경도:</b> {lon:.6f}° E<br>"
            f"<b>주소:</b> <span style='font-size:11px'>{esc(pr.get('address','-'))}</span><br>"
            f"<hr style='margin:4px 0;border-color:#ddd'>"
            f"<b>최초설치:</b> {pr.get('inst_year') or '-'}년 &nbsp; "
            f"<b>최신관측:</b> {pr.get('obs_year','-')}년<br>"
            f"<b>편각:</b> {vf(pr.get('decl'),'°')} &nbsp; <b>복각:</b> {vf(pr.get('incl'),'°')} &nbsp; "
            f"<b>총자력:</b> {vf(pr.get('total'),' nT')}</div>")
        folium.Marker(
            [lat, lon],
            icon=folium.DivIcon(html="<div style='font-size:22px;line-height:1;"
                                "text-shadow:1px 1px 2px rgba(0,0,0,.4)'>⭐</div>",
                                icon_anchor=(11, 11)),
            tooltip=f"기존 측정점: {name}",
            popup=folium.Popup(pop, max_width=320)).add_to(fg)
        n += 1
    fg.add_to(m)
    return n


def build(records):
    m = folium.Map(location=[36.3, 127.8], zoom_start=7, tiles=None, prefer_canvas=True)
    folium.TileLayer(
        tiles="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
        attr='© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name="OpenStreetMap", overlay=False, control=True, max_zoom=19).add_to(m)
    folium.TileLayer(
        tiles=("https://server.arcgisonline.com/ArcGIS/rest/services"
               "/World_Imagery/MapServer/tile/{z}/{y}/{x}"),
        attr="Esri World Imagery", name="위성 이미지 (Esri)",
        overlay=False, control=True).add_to(m)

    groups = {g: folium.FeatureGroup(name=f"{lab}", show=True)
              for g, (c, lab) in GRADE.items()}
    counts = {}
    for d in records:
        lat, lon = fnum(d["위도"]), fnum(d["경도"])
        if lat is None or lon is None:
            continue
        grade, concl, note = review(d)
        counts[grade] = counts.get(grade, 0) + 1
        color = GRADE[grade][0]
        folium.CircleMarker(
            location=[lat, lon], radius=7, color="#333", weight=1.2,
            fill=True, fill_color=color, fill_opacity=0.92,
            tooltip=f"[{grade}] {d['관리번호']} {d['후보지명']}",
            popup=folium.Popup(popup_html(d, grade, concl, note), max_width=360),
        ).add_to(groups[grade])
    for g in GRADE:
        groups[g].add_to(m)

    add_topo_layer(m)
    n_exist = add_existing_layer(m)

    folium.LayerControl(collapsed=False).add_to(m)
    total = sum(counts.values())
    m.get_root().html.add_child(folium.Element(legend_html(counts, total)))

    title = (
        "<div style='position:fixed;top:12px;left:50%;transform:translateX(-50%);"
        "z-index:9999;background:rgba(31,56,100,.95);color:#fff;padding:8px 20px;"
        "border-radius:8px;box-shadow:0 1px 6px rgba(0,0,0,.3);"
        "font-family:\"맑은 고딕\",sans-serif;font-size:15px;font-weight:bold'>"
        "지자기 도상선점 — 현장조사 선점 검토"
        f"<span style='font-size:11px;font-weight:normal;opacity:.8;margin-left:10px'>"
        f"{datetime.now():%Y-%m-%d} 기준 · {total}개 후보지</span></div>")
    m.get_root().html.add_child(folium.Element(title))
    return m, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=None)
    a = ap.parse_args()
    files = survey_files(a.dir)
    PHOTO_DIR.mkdir(parents=True, exist_ok=True)
    recs = []
    for f in files:
        cards = parse_workbook(f, photo_dir=PHOTO_DIR)
        nph = sum(len(c.get("사진", {})) for c in cards)
        print(f"  {f.name:28} 카드 {len(cards):3}건  사진 {nph:3}장")
        recs += cards
    recs.sort(key=lambda d: (d["관할본부"] or "", d["관리번호"] or ""))
    m, counts = build(recs)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(OUT))
    print(f"저장: {OUT}")
    print("등급 분포:", counts, " 합계:", sum(counts.values()))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
