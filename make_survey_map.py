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

from aggregate_survey_xlsx import (parse_workbook, review, key_disturb, survey_files,
                                   sheet_priority, mark_max_dist)

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs" / "survey_review.html"
PHOTO_DIR = ROOT / "docs" / "survey_photos"

# 기존 지자기 측정점(24) — 현황표 1차 비고 초록(15) ∪ 2차 비고 파랑(9). (사용자 확정)
EXISTING_USE = [
    # 기존 지자기 측정점 30점 (2026-08-19 확정)
    "남양", "춘양", "청송", "상주", "미원", "서산", "청양", "이원", "부안",
    "와도", "순창", "강화", "함양", "가야", "영천", "양산", "남지", "거제",
    "남해", "순천", "장흥", "조도", "포천", "성산", "여주", "화천", "설악",
    "봉평", "삼척", "제천",
]

# 확정 이름이 existing_sites.geojson 과 직접 일치하므로 별칭이 필요 없다.
# (경주=영천 · 언양=양산 · 임계=삼척 도 확정 이름 쪽이 geojson 에 함께 실려 있다.)
EXISTING_ALIAS = {}

# 기존 30점 중 **선점 대상 15점**.
# 나머지 15점은 저수지(제방)에 설치된 점으로 판정되어 선점에서 뺀다 —
# 지반 거동·수위 변동으로 표석 지속성이 확보되지 않는다.
EXISTING_TARGET = [
    "상주", "미원", "부안", "강화", "양산", "남지", "거제", "조도",
    "포천", "성산", "여주", "화천", "봉평", "삼척", "제천",
]

GRADE = {   # 등급 → (색, 라벨) — 4단계
    "A": ("#2E8B57", "A · 선점 가능 (자기구배 조사)"),
    "B": ("#2E86C1", "B · 조건부 선점 가능 (방위표지 거리 확보 필요)"),
    "C": ("#E0A400", "C · 현장 확인 필요 (자기교란 재확인)"),
    "D": ("#CC3333", "D · 부적합 (대체 후보지 검토)"),
    "미완료": ("#888888", "미완료 · 재조사 필요"),
}

# 표출 레이어 — 자침편각 표기는 도엽당 1점이면 충분하므로 B 를 대표/예비로 분리.
LAYERS = [
    ("A",      "#2E8B57", "A · 선점 가능 (자기구배 조사)"),
    ("B대표",  "#2E86C1", "B · 조건부 선점 가능 (도엽 대표)"),
    ("B예비",  "#AFCBE3", "└ B 예비 (도엽 중복·A 확보 도엽)"),
    ("C",      "#E0A400", "C · 현장 확인 필요 (자기교란 재확인)"),
    ("D",      "#CC3333", "D · 부적합 (대체 후보지 검토)"),
    ("미완료", "#888888", "미완료 · 재조사 필요"),
]
LAYER_COL = {k: c for k, c, _ in LAYERS}


def layer_of(grade, prio):
    """등급 + 도엽 구분 → 표출 레이어 키."""
    if grade == "B" and prio:
        return "B대표" if prio[0] == "대표" else "B예비"
    return grade


def fnum(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            if s is not None else "")


def az_text(d):
    """방위표지 1·2 방위각·거리 요약.

    값은 `fix_azimuth()` 가 좌표 기준으로 정정한 것이다. 정정된 항목은 ✔ 를 달고
    카드 기재값을 회색으로 병기해 현장 재확인이 가능하게 둔다.
    """
    det = d.get("방위표지상세", {})
    out = []
    for i, k in enumerate(("표지1", "표지2"), 1):
        c = det.get(k, {})
        az, dist = c.get("방위각", ""), c.get("거리", "")
        if not (az and az != "-"):
            continue
        seg = f"표지{i} 방위각 {esc(az)}"
        if c.get("방위각정정"):
            seg += (f"<span style='color:#B8860B'> ✔</span>"
                    f"<span style='color:#999'>(카드 {esc(c.get('방위각기재') or '-')})</span>")
        if dist and dist != "-":
            seg += f" · {esc(dist)}m"
            if c.get("거리정정"):
                seg += (f"<span style='color:#B8860B'> ✔</span>"
                        f"<span style='color:#999'>(카드 {esc(c.get('거리기재') or '-')}m)</span>")
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


def bangwi_entry(d):
    """후보지 클릭 시 지도에 찍을 방위표지 데이터. 좌표 없으면 None."""
    base, m1, m2 = d.get("기준점ll"), d.get("표지1ll"), d.get("표지2ll")
    if not (base and m1 and m2):
        return None
    det = d.get("방위표지상세", {})

    def leg(tag):
        c = det.get(tag, {})
        return _fmt_az(c.get("방위각", "")), (c.get("거리", "") or "-")
    az1, d1 = leg("표지1")
    az2, d2 = leg("표지2")
    return {
        "name": d["후보지명"], "mid": d["관리번호"],
        "base": [round(base[1], 7), round(base[0], 7)],   # [lat, lon]
        "m1": [round(m1[1], 7), round(m1[0], 7)], "az1": az1, "d1": d1,
        "m2": [round(m2[1], 7), round(m2[0], 7)], "az2": az2, "d2": d2,
    }


def bangwi_script(data):
    """후보지 마커 클릭 → 기준점·방위표지1·2 를 지도에 찍고 방위각·거리 표시."""
    payload = json.dumps(data, ensure_ascii=False)
    return """
<script>
(function(){
  var DATA = __DATA__;
  function pin(color, ch){
    return L.divIcon({className:'', iconAnchor:[13,30], html:
      "<div style='position:relative'><div style='width:22px;height:22px;border-radius:50% 50% 50% 0;"
      +"background:"+color+";border:2px solid #fff;transform:rotate(-45deg);box-shadow:0 1px 3px rgba(0,0,0,.5)'></div>"
      +"<div style='position:absolute;top:3px;left:0;width:22px;text-align:center;color:#fff;"
      +"font:bold 12px sans-serif'>"+ch+"</div></div>"});
  }
  function init(){
    var mk = Object.keys(window).find(function(k){return k.indexOf('map_')===0 && window[k] instanceof L.Map;});
    if(!mk){ return setTimeout(init,200); }
    var map = window[mk];
    var grp = L.layerGroup().addTo(map);
    function plot(b){
      grp.clearLayers();
      var base=b.base, m1=b.m1, m2=b.m2;
      L.polyline([base,m1],{color:'#1D4ED8',weight:2.5}).addTo(grp);
      L.polyline([base,m2],{color:'#0E8A6B',weight:2.5}).addTo(grp);
      L.circleMarker(base,{radius:7,color:'#fff',weight:2,fillColor:'#E8531F',fillOpacity:1})
        .bindTooltip('기준점 · '+b.name,{direction:'top'}).addTo(grp);
      L.marker(m1,{icon:pin('#1D4ED8','1')})
        .bindTooltip('<b>방위표지 1</b><br>방위각 '+b.az1+'<br>거리 '+b.d1+' m',
          {permanent:true,direction:'right',offset:[8,-14],className:'bw-tip'}).addTo(grp);
      L.marker(m2,{icon:pin('#0E8A6B','2')})
        .bindTooltip('<b>방위표지 2</b><br>방위각 '+b.az2+'<br>거리 '+b.d2+' m',
          {permanent:true,direction:'right',offset:[8,-14],className:'bw-tip'}).addTo(grp);
      map.fitBounds(L.latLngBounds([base,m1,m2]).pad(0.8),{maxZoom:18});
    }
    map.eachLayer(function(l){
      if(l instanceof L.CircleMarker && l.getLatLng){
        var ll=l.getLatLng();
        var key=ll.lat.toFixed(5)+','+ll.lng.toFixed(5);
        if(DATA[key]){ l.on('click', function(){ plot(DATA[key]); }); }
      }
    });
  }
  init();
})();
</script>
<style>.bw-tip{font-size:11px;line-height:1.3;border-color:#888}</style>
""".replace("__DATA__", payload)


def popup_html(d, grade, concl, note, rank=None):
    color = GRADE[grade][0]
    dist = key_disturb(d) or "없음"
    bang = d["방위표지"] or "-"
    az = az_text(d)
    rows = [
        ("종합 판정", esc(d["종합판정"])),
        ("핵심 교란요인", esc(dist)),
        ("방위표지", esc(bang) + (f" <span style='color:#888'>({az})</span>" if az else "")),
        ("검토 결론", f"<b>{esc(concl)}</b>"),
    ]
    # 도엽 대표/예비 (자침편각 표기는 도엽당 1점 → 방위표지 최장거리 1순위가 대표)
    if rank and grade == "B":
        kind, rk, n, why = rank
        md = mark_max_dist(d) or 0
        if kind == "대표":
            txt = (f"<b style='color:#2E8B57'>★ {esc(d['도엽명'])} 도엽 대표</b>"
                   + (f" <span style='color:#888'>({rk}/{n}위)</span>" if n > 1 else ""))
        else:
            txt = (f"<b style='color:#888'>예비 — {esc(why)}</b>"
                   + (f" <span style='color:#888'>({rk}/{n}위)</span>" if n > 1 else ""))
        rows.append(("도엽 구분",
                     txt + f" <span style='color:#888'>· 방위표지 최장 {md:.0f}m</span>"))
    if d.get("방위정정"):
        rows.append(("방위 정정",
                     f"<span style='color:#B8860B'>✔ 좌표 기준으로 정정</span>"
                     f"<div style='color:#777;font-size:11px;margin-top:2px'>"
                     f"{esc(d['방위정정'])}</div>"))
    gap = d.get("후보지이격")
    if gap is not None and gap >= 50:
        # 후보지는 도상에서 고른 점이고 성과 기준은 현장 기준점이다.
        # 이격이 크다는 것은 도상 후보지로 **진입이 어려워 옮겼다**는 뜻(접근성 지표).
        rows.append(("도상 후보지 이격",
                     f"{gap:,.0f} m <span style='color:#888'>— 도상 후보지 진입이 어려워 "
                     f"현장에서 이동{' (1km 초과 · 사실상 별개 지점)' if gap >= 1000 else ''}"
                     f"</span>"))
    rows += [
        ("조사자 의견", esc(note)),
        ("조사", f"{esc(d['조사일'])} · {esc(d['조사자'])}"),
        ("도상 후보지 좌표", f"{esc(d['위도'])}, {esc(d['경도'])} "
         f"(표고 {esc(d['표고'])} m)"),
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
        f"{photo_html(d)}"
        f"<div style='margin-top:6px;color:#999;font-size:11px'>{esc(d['관할본부'])}</div>"
        f"</div>")


def legend_html(counts, total, n_tgt=0, n_oth=0):
    items = "".join(
        f"<div style='display:flex;align-items:center;margin:3px 0"
        f"{';padding-left:12px' if k == 'B예비' else ''}'>"
        f"<span style='width:13px;height:13px;border-radius:50%;background:{c};"
        f"display:inline-block;margin-right:7px;border:1px solid #fff;"
        f"box-shadow:0 0 0 1px #999'></span>"
        f"<span style='font-size:12px{';color:#666' if k == 'B예비' else ''}'>{lab} "
        f"<b>{counts.get(k,0)}</b></span></div>"
        for k, c, lab in LAYERS)
    return (
        "<div style='position:fixed;bottom:22px;left:22px;z-index:9999;"
        "background:rgba(255,255,255,.95);padding:11px 14px;border-radius:8px;"
        "box-shadow:0 1px 6px rgba(0,0,0,.3);font-family:\"맑은 고딕\",sans-serif'>"
        "<div style='font-weight:bold;font-size:13px;margin-bottom:6px'>"
        f"선점 검토 등급 <span style='color:#888;font-weight:normal'>(총 {total}점)</span></div>"
        f"{items}"
        "<div style='border-top:1px solid #ddd;margin:7px 0 5px'></div>"
        "<div style='display:flex;align-items:center;margin:3px 0'>"
        "<span style='width:15px;height:15px;border-radius:50%;background:#D7A400;"
        "border:2px solid #fff;box-shadow:0 0 0 1px #8B6A00;color:#fff;"
        "font:bold 10px/13px sans-serif;text-align:center;display:inline-block;"
        "margin-right:5px'>★</span>"
        f"<span style='font-size:12px'>기존점 — 선점 대상 <b>{n_tgt}</b></span></div>"
        "<div style='display:flex;align-items:center;margin:3px 0;padding-left:12px'>"
        "<span style='display:inline-block;width:15px;margin-right:5px;"
        "font-size:13px;text-align:center'>⭐</span>"
        f"<span style='font-size:12px;color:#666'>기존점 — 기타 <b>{n_oth}</b> "
        "<span style='color:#999'>(저수지·제방)</span></span></div>"
        "<div style='border-top:1px solid #ddd;margin:7px 0 5px'></div>"
        "<div style='font-size:11px;color:#777;max-width:210px;line-height:1.45'>"
        "<b style='color:#B8860B'>✔</b> 방위각·거리는 기준점·방위표지 좌표로 "
        "재계산한 값. 카드 기재값은 팝업에 병기.</div>"
        "</div>")


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


def _existing_popup(name, pr, target):
    """기존 측정점 팝업."""
    vf = lambda x, u="": f"{x:.4f}{u}" if isinstance(x, (int, float)) else "-"
    head = ("🎯 선점 대상 기존점" if target else "⭐ 기존 측정점")
    col = "#B8860B" if target else "#8B4513"
    tail = ("<div style='margin-top:5px;padding:4px 6px;background:#FFF6DC;"
            "border-left:3px solid #D7A400;font-size:11px;color:#6B5200'>"
            f"기존 {len(EXISTING_USE)}점 중 <b>선점 대상 {len(EXISTING_TARGET)}점</b>"
            " — 표석 지속성이 확보되어 재설치 없이 선점 가능.</div>"
            ) if target else (
            "<div style='margin-top:5px;padding:4px 6px;background:#F4F4F4;"
            "border-left:3px solid #AAA;font-size:11px;color:#555'>"
            "저수지(제방) 설치점 — 지반 거동·수위 변동으로 표석 지속성이 "
            "확보되지 않아 선점 대상에서 제외.</div>")
    return (
        f"<div style='font-family:\"맑은 고딕\",sans-serif;font-size:12.5px;min-width:250px'>"
        f"<b style='color:{col}'>{head}: {esc(name)}</b><hr style='margin:4px 0'>"
        f"<b>위도:</b> {pr['lat']:.6f}° N &nbsp; <b>경도:</b> {pr['lon']:.6f}° E<br>"
        f"<b>주소:</b> <span style='font-size:11px'>{esc(pr.get('address','-'))}</span><br>"
        f"<hr style='margin:4px 0;border-color:#ddd'>"
        f"<b>최초설치:</b> {pr.get('inst_year') or '-'}년 &nbsp; "
        f"<b>최신관측:</b> {pr.get('obs_year','-')}년<br>"
        f"<b>편각:</b> {vf(pr.get('decl'),'°')} &nbsp; <b>복각:</b> {vf(pr.get('incl'),'°')} &nbsp; "
        f"<b>총자력:</b> {vf(pr.get('total'),' nT')}"
        f"{tail}</div>")


def add_existing_layer(m):
    """기존 지자기 측정점(EXISTING_USE 30) 토글 — 두 레이어로 분리.

      · 🎯 선점 대상 (15) — 표석 지속성이 확보되는 점. 강조 마커, 기본 켬.
      · ⭐ 기타 기존점 (15) — 저수지(제방) 설치점. 종전 ⭐ 마커, 기본 끔.

    반환 (선점대상 수, 기타 수).
    """
    p = DATA / "existing_sites.geojson"
    if not p.exists():
        return 0, 0
    fc = json.load(open(p, encoding="utf-8"))
    by_name = {f["properties"].get("name"): f for f in fc["features"]}
    tgt = {EXISTING_ALIAS.get(n, n) for n in EXISTING_TARGET}
    rest = [n for n in EXISTING_USE if EXISTING_ALIAS.get(n, n) not in tgt]

    for n in EXISTING_TARGET:                     # 목록 자체의 오타·누락 방어
        if EXISTING_ALIAS.get(n, n) not in by_name:
            print(f"    ! 선점 대상 미매칭: {n}", file=sys.stderr)

    fg_t = folium.FeatureGroup(
        name=f"🎯 기존 측정점 — 선점 대상 ({len(EXISTING_TARGET)})", show=True)
    fg_o = folium.FeatureGroup(
        name=f"⭐ 기존 측정점 — 기타 ({len(rest)})", show=False)

    def put(name, fg, target):
        gname = EXISTING_ALIAS.get(name, name)
        f = by_name.get(gname)
        if not f:
            print(f"    ! 기존측정점 미매칭: {name}", file=sys.stderr)
            return 0
        pr = f["properties"]
        if target:
            # 선점 대상: 금색 원판 + ★ — ⭐ 이모지보다 눈에 띄고 등급 마커와도 구분된다.
            icon = folium.DivIcon(html=(
                "<div style='width:26px;height:26px;border-radius:50%;"
                "background:#D7A400;border:2.5px solid #fff;color:#fff;"
                "font:bold 15px/22px sans-serif;text-align:center;"
                "box-shadow:0 0 0 1.5px #8B6A00,0 1px 4px rgba(0,0,0,.5)'>★</div>"),
                icon_anchor=(13, 13))
            tip = f"🎯 선점 대상 기존점: {name}"
        else:
            icon = folium.DivIcon(html=("<div style='font-size:22px;line-height:1;"
                                        "text-shadow:1px 1px 2px rgba(0,0,0,.4)'>⭐</div>"),
                                  icon_anchor=(11, 11))
            tip = f"기존 측정점: {name}"
        folium.Marker([pr["lat"], pr["lon"]], icon=icon, tooltip=tip,
                      popup=folium.Popup(_existing_popup(name, pr, target),
                                         max_width=320)).add_to(fg)
        return 1

    nt = sum(put(n, fg_t, True) for n in EXISTING_TARGET)
    no = sum(put(n, fg_o, False) for n in rest)
    fg_t.add_to(m)
    fg_o.add_to(m)
    return nt, no



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

    groups = {k: folium.FeatureGroup(name=lab, show=True) for k, c, lab in LAYERS}
    counts = {}
    bangwi = {}
    prio = sheet_priority(records)   # B 도엽 대표/예비
    for d in records:
        lat, lon = fnum(d["위도"]), fnum(d["경도"])
        if lat is None or lon is None:
            continue
        grade, concl, note = review(d)
        rk = prio.get(d["관리번호"])
        key = layer_of(grade, rk)
        counts[key] = counts.get(key, 0) + 1
        color = LAYER_COL[key]
        tip = f"[{grade}] {d['관리번호']} {d['후보지명']}"
        if rk and grade == "B":
            tip += (f" · ★{d['도엽명']}도엽 대표" if rk[0] == "대표"
                    else f" · {d['도엽명']}도엽 예비({rk[3]})")
        tip += " (클릭: 방위표지 표시)"
        is_top = bool(rk and grade == "B" and rk[0] == "대표")
        folium.CircleMarker(
            location=[lat, lon], radius=8 if is_top else 7,
            color="#fff" if is_top else "#333", weight=2.5 if is_top else 1.2,
            fill=True, fill_color=color, fill_opacity=0.92,
            tooltip=tip,
            popup=folium.Popup(popup_html(d, grade, concl, note, rk), max_width=360),
        ).add_to(groups[key])
        be = bangwi_entry(d)
        if be:
            bangwi[f"{lat:.5f},{lon:.5f}"] = be
    for k, _c, _lab in LAYERS:
        groups[k].add_to(m)

    add_topo_layer(m)
    n_tgt, n_oth = add_existing_layer(m)

    folium.LayerControl(collapsed=False).add_to(m)
    total = sum(counts.values())
    m.get_root().html.add_child(folium.Element(legend_html(counts, total, n_tgt, n_oth)))
    m.get_root().html.add_child(folium.Element(bangwi_script(bangwi)))

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
