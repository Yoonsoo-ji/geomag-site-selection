# -*- coding: utf-8 -*-
"""
3D 지구본(geomag_globe.html) 3단계용 정적 데이터 생성.
  · docs/data/survey_sites.geojson  — 선점 후보 103(등급 A/B/C·판정·좌표)
  · docs/data/existing_pts.geojson  — 기존 측정점(EXISTING_USE)

선점 데이터는 조사서(aggregate_survey_xlsx)에서, 기존점은 existing_sites.geojson 에서.
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"

from aggregate_survey_xlsx import survey_files, parse_workbook, review, key_disturb
from make_survey_map import EXISTING_USE, EXISTING_ALIAS


def fnum(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


# ── 선점 후보 103 ──
feats = []
for f in survey_files():
    for d in parse_workbook(f):
        lat, lon = fnum(d["위도"]), fnum(d["경도"])
        if lat is None or lon is None:
            continue
        grade, concl, note = review(d)
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [round(lon, 6), round(lat, 6)]},
            "properties": {
                "mid": d["관리번호"], "name": d["후보지명"], "hq": d["관할본부"],
                "grade": grade, "verdict": d["종합판정"], "bang": d["방위표지"] or "-",
                "disturb": key_disturb(d) or "없음", "concl": concl,
                "surveyor": d["조사자"], "date": d["조사일"], "elev": d["표고"],
            },
        })
feats.sort(key=lambda x: x["properties"]["mid"])
out1 = {"type": "FeatureCollection", "features": feats}
(DATA / "survey_sites.geojson").write_text(
    json.dumps(out1, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
gc = {}
for x in feats:
    gc[x["properties"]["grade"]] = gc.get(x["properties"]["grade"], 0) + 1
print(f"survey_sites.geojson: {len(feats)}건  등급 {gc}")

# ── 기존 측정점 24 ──
ex = json.load(open(DATA / "existing_sites.geojson", encoding="utf-8"))
by_name = {f["properties"].get("name"): f for f in ex["features"]}
efeats = []
for name in EXISTING_USE:
    gname = EXISTING_ALIAS.get(name, name)
    f = by_name.get(gname)
    if not f:
        print("  ! 기존점 미매칭:", name)
        continue
    p = f["properties"]
    efeats.append({
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [round(p["lon"], 6), round(p["lat"], 6)]},
        "properties": {
            "name": name, "address": p.get("address", "-"),
            "inst_year": p.get("inst_year"), "obs_year": p.get("obs_year"),
            "decl": p.get("decl"), "incl": p.get("incl"), "total": p.get("total"),
        },
    })
out2 = {"type": "FeatureCollection", "features": efeats}
(DATA / "existing_pts.geojson").write_text(
    json.dumps(out2, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
print(f"existing_pts.geojson: {len(efeats)}건")
