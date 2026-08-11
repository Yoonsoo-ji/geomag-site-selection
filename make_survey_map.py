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
import re
import sys
from datetime import datetime
from pathlib import Path

import folium

from aggregate_survey_xlsx import (DEF_DIR, parse_workbook, review, key_disturb)

ROOT = Path(__file__).parent
OUT = ROOT / "docs" / "survey_review.html"

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


def popup_html(d, grade, concl, note):
    color = GRADE[grade][0]
    dist = key_disturb(d) or "없음"
    bang = d["방위표지"] or "-"
    rows = [
        ("종합 판정", esc(d["종합판정"])),
        ("핵심 교란요인", esc(dist)),
        ("방위표지", esc(bang)),
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
        f"width:320px;line-height:1.45'>"
        f"<div style='background:{color};color:#fff;padding:6px 10px;margin:-2px -2px 6px;"
        f"border-radius:4px 4px 0 0;font-weight:bold'>"
        f"[{grade}] {esc(d['관리번호'])} · {esc(d['후보지명'])}</div>"
        f"<table style='border-collapse:collapse'>{body}</table>"
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
    ap.add_argument("--dir", default=str(DEF_DIR))
    a = ap.parse_args()
    files = sorted(Path(a.dir).glob("*.xlsx"),
                   key=lambda p: int(re.match(r"(\d+)", p.name).group(1))
                   if re.match(r"(\d+)", p.name) else 999)
    recs = []
    for f in files:
        recs += parse_workbook(f)
    recs.sort(key=lambda d: (d["관할본부"] or "", d["관리번호"] or ""))
    m, counts = build(recs)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(OUT))
    print(f"저장: {OUT}")
    print("등급 분포:", counts, " 합계:", sum(counts.values()))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
