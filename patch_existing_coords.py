# -*- coding: utf-8 -*-
"""
기존 측정점 좌표 교정을 **생성기를 못 돌리는 산출물**에 반영한다.

`docs/survey_review.html` 은 Folium 이 좌표를 **HTML 안에 박아** 낸다(fetch 아님).
그래서 `existing_sites.geojson` 을 고쳐도 자동으로 따라오지 않는다. 그런데
그 생성기(`make_survey_map.py`)는 현장조사 회신본 폴더
(`20260811_사전 현장 조사서 취합`)를 요구하는데 저장소에 없다 — 재생성이 불가하다.

그래서 **이미 만들어진 HTML 의 좌표만 이름으로 찾아 교정**한다. 마커 좌표와
팝업에 표시된 위도·경도 문자열 둘 다 고친다(둘이 어긋나면 더 나쁘다).

`docs/data/existing_pts.geojson`(지구본용)은 정적 파일이라 그냥 다시 만든다 —
`export_globe_sites.py` 는 조사서까지 읽지만 기존점 부분만 떼어 쓰면 된다.

    python patch_existing_coords.py [--check]

⚠️ 원본 회신본이 확보되면 `make_survey_map.py` 를 정상 재실행하는 것이 옳다.
   이 스크립트는 그때까지의 임시 경로이며, **멱등**이다(이미 맞으면 건드리지 않음).
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
HTML = ROOT / "docs" / "survey_review.html"

import existing_network as EN


def load_truth():
    """이름 → (lat, lon). `existing_sites.geojson` 이 기준(대조표 교정 반영본)."""
    gj = json.load(open(DATA / "existing_sites.geojson", encoding="utf-8"))
    return {f["properties"]["name"]: (float(f["properties"]["lat"]),
                                      float(f["properties"]["lon"]))
            for f in gj["features"]}


# ======================================================================
# survey_review.html
# ======================================================================

def patch_html(truth, check_only=False):
    if not HTML.exists():
        print(f"[건너뜀] {HTML.name} 없음")
        return 0
    s = HTML.read_text(encoding="utf-8")

    # marker → popup → html 연결을 되짚어 이름을 찾는다
    m2p = dict(re.findall(r"marker_(\w+)\.bindPopup\(popup_(\w+)\)", s))
    p2h = dict(re.findall(r"popup_(\w+)\.setContent\(html_(\w+)\)", s))
    html_txt = {mid: body for mid, body in
                re.findall(r"var html_(\w+) = \$\(`(.*?)`\)", s, re.S)}

    n_fix = 0
    report = []
    for mm in list(re.finditer(
            r"var marker_(\w+) = L\.marker\(\s*\n?\s*\[([-\d.]+), ([-\d.]+)\]", s)):
        mid, la_s, lo_s = mm.group(1), mm.group(2), mm.group(3)
        body = html_txt.get(p2h.get(m2p.get(mid, ""), ""), "")
        # 두 갈래 팝업을 모두 받는다 — 「⭐ 기존 측정점: 이름」(기타 15)과
        # 「🎯 선점 대상 기존점: 이름」(선점 대상 15). 하나만 보면 절반을 놓친다.
        nm = re.search(r"(?:기존 측정점|선점 대상 기존점): ([^<]+)</b>", body)
        if not nm:
            continue                                   # 선점 후보 마커 등
        name = nm.group(1).strip()
        if name not in truth:
            report.append((name, "정본에 없음", None))
            continue
        la, lo = truth[name]
        if abs(float(la_s) - la) < 1e-7 and abs(float(lo_s) - lo) < 1e-7:
            continue                                   # 이미 맞음 (멱등)
        report.append((name, f"{float(la_s):.6f},{float(lo_s):.6f}",
                       f"{la:.6f},{lo:.6f}"))
        n_fix += 1
        if not check_only:
            s = s.replace(mm.group(0),
                          f"var marker_{mid} = L.marker(\n                [{la}, {lo}]", 1)
            # 팝업에 찍힌 좌표 문자열도 함께
            old_disp = (f"<b>위도:</b> {float(la_s):.6f}° N &nbsp; "
                        f"<b>경도:</b> {float(lo_s):.6f}° E")
            new_disp = (f"<b>위도:</b> {la:.6f}° N &nbsp; "
                        f"<b>경도:</b> {lo:.6f}° E")
            if old_disp in s:
                s = s.replace(old_disp, new_disp, 1)

    print(f"■ survey_review.html — 교정 대상 {n_fix}점")
    for nm, a, b in report:
        print(f"    {nm:<6} {a}  →  {b}")
    if n_fix and not check_only:
        HTML.write_text(s, encoding="utf-8")
        print(f"    [저장] {HTML.name}")
    return n_fix


# ======================================================================
# existing_pts.geojson (지구본)
# ======================================================================

def rebuild_globe_pts(check_only=False):
    src = json.load(open(DATA / "existing_sites.geojson", encoding="utf-8"))
    by = {f["properties"]["name"]: f for f in src["features"]}
    feats = []
    missing = []
    for name in EN.EXISTING_NETWORK:
        g = EN.EXISTING_ALIAS.get(name, name)
        f = by.get(g)
        if f is None:
            missing.append(name)
            continue
        p = f["properties"]
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point",
                         "coordinates": [round(p["lon"], 6), round(p["lat"], 6)]},
            "properties": {k: p.get(k) for k in
                           ("name", "address", "inst_year", "obs_year",
                            "elev", "decl", "incl", "total", "net")},
        })
    out = DATA / "existing_pts.geojson"
    old = json.load(open(out, encoding="utf-8")) if out.exists() else None
    same = old is not None and old.get("features") == feats
    print(f"■ existing_pts.geojson — {len(feats)}점"
          + (f" · 누락 {missing}" if missing else "")
          + ("  (변경 없음)" if same else "  (갱신 필요)"))
    if not same and not check_only:
        out.write_text(json.dumps({"type": "FeatureCollection",
                                   "features": feats}, ensure_ascii=False),
                       encoding="utf-8")
        print(f"    [저장] {out.name}")
    return 0 if same else 1



# ======================================================================
# 「점조서 좌표(현장 확인용)」 레이어
# ======================================================================
#
# ⚠️ **좌표를 덮었다고 옛 자리를 지우지 않는다.** 어느 쪽에 표석이 있는지는
# 현장에 나가 봐야 알 수 있으므로 두 자리를 **함께** 찍어 둔다. Folium 이
# 좌표를 HTML 에 박아 내므로 여기서도 마커·팝업·레이어 컨트롤 항목을 직접
# 끼워 넣는다. **멱등**이다 — 표식(`PREV_MARK`)이 이미 있으면 손대지 않는다.
#
# ⚠️ `make_survey_map.py` 를 다시 돌리면 이 레이어는 사라진다. 그때는 이
#    스크립트를 다시 실행하면 된다.

PREV_LAYER_LABEL = "▫️ 점조서 좌표 (현장 확인용)"
PREV_MARK = "prevcoord"            # 멱등 판정용 표식


def _km(la0, lo0, la1, lo1):
    import math
    m = math.radians((la0 + la1) / 2)
    return math.hypot((lo1 - lo0) * 111.320 * math.cos(m),
                      (la1 - la0) * 110.574)


def _tpl(t):
    """JS 템플릿 리터럴에 넣을 문자열 이스케이프."""
    return t.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")


def patch_prev_markers(check_only=False):
    """`COORD_OVERRIDE` 측점마다 **점조서 원좌표** 마커를 새 레이어로 얹는다."""
    if not HTML.exists():
        print("[건너뜀] survey_review.html 없음")
        return 0
    s = HTML.read_text(encoding="utf-8")
    if PREV_MARK in s:
        print("■ 점조서 좌표 레이어 — 이미 있음")
        return 0
    if not EN.COORD_OVERRIDE:
        print("■ 점조서 좌표 레이어 — 대상 없음")
        return 0

    mv = re.search(r"var (map_\w+) = L\.map\(", s)
    ctrl = re.search(r"\n\s*var layer_control_\w+_layers = \{", s)
    if not (mv and ctrl):
        print("■ 점조서 좌표 레이어 — HTML 구조를 못 찾음(건너뜀)")
        return 1

    raw = EN.load_register(override=False).set_index("지점명")
    js = ["", f"        // ---- {PREV_MARK}: 점조서 원좌표(현장 확인용) ----",
          f"        var fg_{PREV_MARK} = L.featureGroup({{}})"
          f".addTo({mv.group(1)});"]
    n = 0
    for nm, o in EN.COORD_OVERRIDE.items():
        if nm not in raw.index:
            continue
        la, lo = float(raw.at[nm, "위도"]), float(raw.at[nm, "경도"])
        d = _km(la, lo, o["lat"], o["lon"])
        body = (
            "<div style='font-family:&quot;맑은 고딕&quot;,sans-serif;"
            "font-size:12.5px;min-width:250px'>"
            f"<b style='color:#555'>▫️ {nm} — "
            "점조서 좌표 (현장 확인용)</b><hr style='margin:4px 0'>"
            f"<b>위도:</b> {la:.6f}° N &nbsp; "
            f"<b>경도:</b> {lo:.6f}° E<br>"
            f"<b>현 표기 위치까지:</b> {d:.1f} km<br>"
            "<div style='margin-top:6px;padding:6px;background:#F4F4F4;"
            "border-left:3px solid #AAA;font-size:11px;color:#555'>"
            f"{o.get('prev_note', '')}</div></div>")
        icon = ("<div style='width:26px;height:26px;border-radius:50%;"
                "background:#fff;border:2.5px solid #666;color:#444;"
                "font:bold 15px/22px sans-serif;text-align:center;"
                "box-shadow:0 1px 4px rgba(0,0,0,.4)'>?</div>")
        v = f"mk_{PREV_MARK}_{n}"
        js += [
            f"        var {v} = L.marker([{la}, {lo}], {{}})"
            f".addTo(fg_{PREV_MARK});",
            f"        {v}.setIcon(L.divIcon("
            + json.dumps({"html": icon, "iconAnchor": [13, 13],
                          "className": "empty"}, ensure_ascii=False) + "));",
            f"        {v}.bindTooltip(`<div>▫️ {nm} — "
            "점조서 좌표 (현장 확인)</div>`, "
            + json.dumps({"sticky": True}) + ");",
            f"        {v}.bindPopup(L.popup("
            + json.dumps({"maxWidth": 340})
            + f").setContent(`{_tpl(body)}`));",
        ]
        n += 1
    if not n:
        print("■ 점조서 좌표 레이어 — 대상 없음")
        return 0

    # 그룹은 레이어 컨트롤 «생성 전»에 있어야 한다
    s = s[:ctrl.start()] + "\n" + "\n".join(js) + "\n" + s[ctrl.start():]
    # 기존 측정점 항목 «뒤»에 붙인다 — overlays 블록의 끝자리
    ov = re.search(r"overlays :  \{\n", s)
    end = s.index("\n                },", ov.end()) + 1
    s = (s[:end]
         + f'                    "{PREV_LAYER_LABEL}" : fg_{PREV_MARK},\n'
         + s[end:])

    print(f"■ 점조서 좌표 레이어 — {n}점 추가")
    if not check_only:
        HTML.write_text(s, encoding="utf-8")
        print(f"    [저장] {HTML.name}")
    return n


def main():
    check = "--check" in sys.argv
    truth = load_truth()
    print(f"기준 {len(truth)}점 (docs/data/existing_sites.geojson)\n")
    a = patch_html(truth, check)
    b = rebuild_globe_pts(check)
    c = patch_prev_markers(check)
    if check:
        print(f"\n[점검] 교정 필요: survey_review {a}점 · globe {b}"
              f" · 점조서 좌표 레이어 {c}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
