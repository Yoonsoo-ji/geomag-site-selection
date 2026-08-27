#!/usr/bin/env python3
"""전문가·발주자용 LMM 시네마틱 웹 브리핑을 단일 HTML로 생성한다.

수치와 지점은 배포 JSON/GeoJSON을 읽어 삽입한다. 결과물은 외부 라이브러리나
네트워크 fetch 없이 ``docs/lmm_cinematic.html`` 하나로 동작한다.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
DATA = DOCS / "data"
OUT = DOCS / "lmm_cinematic.html"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def flatten_coords(value):
    """GeoJSON 좌표 중 [lon, lat] 쌍만 순회한다."""
    if (
        isinstance(value, list)
        and len(value) >= 2
        and isinstance(value[0], (int, float))
        and isinstance(value[1], (int, float))
    ):
        yield float(value[0]), float(value[1])
        return
    if isinstance(value, list):
        for item in value:
            yield from flatten_coords(item)


def compact_boundary(doc):
    """화면용 해안선을 단순화하되 shapely가 없으면 원자료를 그대로 쓴다."""
    geom = doc["features"][0]["geometry"]
    try:
        from shapely.geometry import mapping, shape

        geom = mapping(shape(geom).simplify(0.012, preserve_topology=True))
    except Exception:
        pass
    return geom


def point_rows(doc, fields):
    rows = []
    for feature in doc.get("features", []):
        geom = feature.get("geometry") or {}
        if geom.get("type") != "Point":
            continue
        lon, lat = geom["coordinates"][:2]
        props = feature.get("properties", {})
        row = {"lon": round(float(lon), 6), "lat": round(float(lat), 6)}
        for key in fields:
            row[key] = props.get(key)
        rows.append(row)
    return rows


def density_rows():
    rows = []
    files = {
        "critical": DATA / "density_gap_critical.geojson",
        "short": DATA / "density_gap_short.geojson",
        "mild": DATA / "density_gap_mild.geojson",
        "ok": DATA / "density_gap_ok.geojson",
    }
    counts = {}
    missing_sigma = 0
    for tier, path in files.items():
        doc = read_json(path)
        counts[tier] = len(doc.get("features", []))
        for feature in doc.get("features", []):
            pts = list(flatten_coords(feature["geometry"]["coordinates"]))
            if not pts:
                continue
            xs, ys = zip(*pts)
            props = feature.get("properties", {})
            observed = bool(props.get("sigma_obs", 0))
            missing_sigma += int(not observed)
            rows.append(
                {
                    "lon": round((min(xs) + max(xs)) / 2, 5),
                    "lat": round((min(ys) + max(ys)) / 2, 5),
                    "tier": tier,
                    "observed": observed,
                    "sigma": props.get("sigma_nT"),
                    "required": props.get("spacing_req"),
                    "relative": props.get("deficit_rel"),
                }
            )
    return rows, counts, missing_sigma


def validation_value(rows, component, stage):
    for row in rows:
        if row.get("성분") == component and row.get("단계") == stage:
            return float(row["RMS"])
    raise KeyError((component, stage))


def build_payload():
    model = read_json(DATA / "lmm_model.json")
    diagnosis = read_json(DATA / "lmm_diagnosis.json")
    existing_doc = read_json(DATA / "existing_pts.geojson")
    survey_doc = read_json(DATA / "survey_sites.geojson")
    boundary_doc = read_json(DATA / "korea_boundary.geojson")

    existing = point_rows(existing_doc, ["name", "address", "obs_year"])
    survey = point_rows(survey_doc, ["mid", "name", "grade", "verdict"])
    density, density_counts, sigma_missing = density_rows()

    candidate_counts = {}
    for priority in ("p1", "p2", "p3"):
        candidate_counts[priority] = len(
            read_json(DATA / f"candidates_{priority}.geojson").get("features", [])
        )

    grade_counts = Counter(str(row.get("grade") or "미분류") for row in survey)
    crust_values = model["crustal"]["values"]
    crust_valid = sum(v is not None for v in crust_values)
    validation = model["validation"]
    dataset = diagnosis["dataset"]

    return {
        "generated": model["generated"],
        "epoch": model["epoch"],
        "model_name": model["name"],
        "boundary": compact_boundary(boundary_doc),
        "existing": existing,
        "model_sites": [
            {
                "name": s["name"],
                "lat": round(float(s["lat"]), 6),
                "lon": round(float(s["lon"]), 6),
                "visits": int(s["n_visit"]),
            }
            for s in model["sites"]
        ],
        "survey": survey,
        "density": density,
        "stats": {
            "network_sites": len(existing),
            "model_sites": len(model["sites"]),
            "observations": int(dataset["n_obs"]),
            "survey_sites": len(survey),
            "candidate_sites": sum(candidate_counts.values()),
            "candidate_counts": candidate_counts,
            "survey_grades": dict(grade_counts),
            "density_counts": density_counts,
            "land_cells": sum(density_counts.values()),
            "sigma_missing": sigma_missing,
            "target_scenario": len(existing) + 50,
            "crust_coverage_pct": round(100 * crust_valid / len(crust_values), 1),
            "crust_corr": model["crustal_diagnostics"]["corr"],
            "crust_rms_reduction_pct": model["crustal_diagnostics"][
                "rms_reduction_pct"
            ],
        },
        "model": {
            "observation_years": model["observation_years"],
            "epoch_label": model["epoch_label"],
            "regional_degree": int(model["regional"]["degree"]),
            "regional": {
                "D": float(model["regional"]["D"][0]),
                "I": float(model["regional"]["I"][0]),
                "F": float(model["regional"]["F"][0]),
            },
            "loo": {key: float(value) for key, value in model["loo_cv"].items()},
            "stage": {
                "D_igrf": validation_value(validation, "D_deg", "IGRF"),
                "D_final": validation_value(validation, "D_deg", "+Regional"),
                "I_igrf": validation_value(validation, "I_deg", "IGRF"),
                "I_final": validation_value(validation, "I_deg", "+Regional"),
                "F_igrf": validation_value(validation, "F_nT", "IGRF"),
                "F_crust": validation_value(validation, "F_nT", "+Crustal"),
                "F_final": validation_value(
                    validation, "F_nT", "+Crustal+Regional"
                ),
            },
            "inliers": {
                "D": int(dataset["inlier_D"]),
                "I": int(dataset["inlier_I"]),
                "F": int(dataset["inlier_F"]),
            },
            "external_mode": "D·I 관측값 전처리 보정 / 예측 시 외부장 미평가",
        },
        "sources": [
            {
                "group": "Core",
                "name": "IGRF-14 Gauss 계수",
                "detail": "degree 13 · ppigrf/IGRF14.shc",
            },
            {
                "group": "Regional",
                "name": "지표 절대측정 D·I·F",
                "detail": "2019 야장 + 2022–2025 성과 · 품질감사 후 16측점",
            },
            {
                "group": "Crustal",
                "name": "KIGAM 항공자력이상 격자",
                "detail": "1982–2018 · 1.5분(약 2.8 km) · F에 스칼라 적용",
            },
            {
                "group": "External",
                "name": "CYG·제주·강릉·이천 1분 자료",
                "detail": "현재 subtract_DI: 상태=정상 세션의 D·I 환산 · QC=QUIET 전용 필터 아님 · 예측층 미적용",
            },
            {
                "group": "Network",
                "name": "1등 지자기점 30점 좌표 대조표",
                "detail": "출처 확인 필요 · 원장/정본으로 간주하지 않음",
            },
            {
                "group": "Survey",
                "name": "현장조사 일괄취합 103건",
                "detail": "20260818_211750 · A 6 / B 35 / C 49 / D 13",
            },
        ],
    }


HTML = r'''<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<meta name="color-scheme" content="dark">
<title>한반도 LMM — 관측망에서 국가 자기장 기준으로</title>
<style>
:root{
  --ink:#eaf1f6;--muted:#9eacb9;--dim:#687887;--bg:#071018;--panel:#0d1822;
  --line:#20303d;--cyan:#51d2d7;--blue:#427df4;--orange:#ff7048;--gold:#f0bc54;
  --good:#6cc49a;--danger:#ff775c;--pad:clamp(22px,5.6vw,92px);--max:1480px;
}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-snap-type:y mandatory;background:var(--bg)}
body{margin:0;color:var(--ink);background:var(--bg);font-family:"Pretendard","Noto Sans KR","Malgun Gothic",system-ui,sans-serif;overflow-x:hidden}
body::before{content:"";position:fixed;inset:0;pointer-events:none;z-index:40;opacity:.17;background-image:radial-gradient(rgba(255,255,255,.22) .45px,transparent .6px);background-size:5px 5px;mix-blend-mode:soft-light}
a{color:inherit}.mono{font-family:"IBM Plex Mono","Cascadia Mono",Consolas,monospace;letter-spacing:.02em}
.topbar{position:fixed;left:0;right:0;top:0;z-index:60;height:58px;display:flex;align-items:center;gap:18px;padding:0 clamp(18px,3vw,48px);background:rgba(7,16,24,.84);backdrop-filter:blur(16px);border-bottom:1px solid rgba(255,255,255,.07)}
.brand{font-size:12px;font-weight:800;letter-spacing:.15em;white-space:nowrap}.brand b{color:var(--cyan)}
.progress{height:2px;flex:1;background:#172632;position:relative;overflow:hidden}.progress span{position:absolute;inset:0 auto 0 0;width:0;background:var(--orange);transition:width .35s ease}
.counter{font-size:12px;color:var(--muted);min-width:58px;text-align:right}
.scene-nav{position:fixed;right:18px;top:50%;transform:translateY(-50%);z-index:55;display:grid;gap:10px}
.scene-nav button{width:8px;height:8px;border:0;border-radius:50%;padding:0;background:#506170;opacity:.46;cursor:pointer;transition:.25s}
.scene-nav button.active{height:28px;border-radius:10px;background:var(--orange);opacity:1}
.controls{position:fixed;right:clamp(18px,3vw,48px);bottom:22px;z-index:55;display:flex;gap:8px}
.controls button,.chip{border:1px solid rgba(255,255,255,.16);background:rgba(11,23,33,.82);color:var(--ink);backdrop-filter:blur(12px);padding:10px 14px;font:700 12px inherit;cursor:pointer}
.controls button:hover,.chip:hover{border-color:var(--cyan)}
.scene{position:relative;min-height:100svh;scroll-snap-align:start;padding:clamp(92px,11vh,138px) var(--pad) clamp(74px,8vh,104px);display:flex;align-items:center;isolation:isolate;overflow:hidden}
.inner{width:min(100%,var(--max));margin:auto;position:relative;z-index:2}
.scene::after{content:attr(data-no);position:absolute;right:4vw;bottom:-.09em;font:900 clamp(130px,22vw,360px)/1 "Arial Narrow",sans-serif;color:rgba(255,255,255,.018);z-index:-1}
.eyebrow{display:flex;align-items:center;gap:12px;color:var(--cyan);font:700 12px/1.3 inherit;letter-spacing:.16em;text-transform:uppercase;margin-bottom:19px}
.eyebrow::before{content:"";width:38px;height:2px;background:currentColor}
h1,h2,h3,p{margin-top:0}h1{font-size:clamp(47px,7vw,112px);line-height:.93;letter-spacing:-.065em;max-width:1120px;margin-bottom:28px}
h2{font-size:clamp(36px,5vw,72px);line-height:1.04;letter-spacing:-.045em;margin-bottom:26px;max-width:1080px}
h3{font-size:clamp(18px,1.7vw,28px);line-height:1.25;letter-spacing:-.025em}
.lead{font-size:clamp(17px,1.65vw,25px);line-height:1.62;color:#c5d1da;max-width:930px}.lead strong{color:#fff}
.tag{display:inline-flex;align-items:center;gap:8px;padding:8px 10px;border:1px solid rgba(255,255,255,.13);color:var(--muted);font-size:12px}.tag::before{content:"";width:6px;height:6px;background:var(--cyan)}
.warning{border-left:3px solid var(--orange);background:rgba(255,112,72,.08);padding:15px 18px;color:#f2c7bb;line-height:1.55;font-size:14px}
.reveal{opacity:0;transform:translateY(28px);transition:opacity .72s ease,transform .72s cubic-bezier(.2,.8,.2,1);transition-delay:var(--delay,0ms)}
.scene.active .reveal{opacity:1;transform:none}
.accent{color:var(--orange)}.cyan{color:var(--cyan)}.muted{color:var(--muted)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:clamp(22px,4vw,64px);align-items:center}.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}.grid4{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}
.hero{background:radial-gradient(circle at 72% 42%,rgba(47,125,244,.13),transparent 36%),linear-gradient(135deg,#071018 0%,#09141e 60%,#050b11 100%)}
#fieldCanvas{position:absolute;inset:0;width:100%;height:100%;z-index:0;opacity:.72}
.hero-copy{position:relative;z-index:2}.hero .kicker{display:flex;flex-wrap:wrap;gap:9px;margin-bottom:36px}
.hero h1 span{display:block}.hero h1 .ghost{color:transparent;-webkit-text-stroke:1px rgba(234,241,246,.48)}
.hero-foot{display:flex;gap:28px;align-items:center;margin-top:44px;color:var(--muted);font-size:13px}.scroll-cue{display:flex;align-items:center;gap:10px}.scroll-cue i{display:block;width:24px;height:38px;border:1px solid #577080;border-radius:16px;position:relative}.scroll-cue i::after{content:"";position:absolute;top:7px;left:10px;width:2px;height:8px;background:var(--orange);animation:wheel 1.8s infinite}@keyframes wheel{0%{transform:translateY(0);opacity:0}30%{opacity:1}100%{transform:translateY(14px);opacity:0}}
.flow-track{display:grid;grid-template-columns:repeat(7,1fr);gap:0;margin-top:54px;position:relative}.flow-track::before{content:"";position:absolute;top:28px;left:6%;right:6%;height:1px;background:linear-gradient(90deg,var(--cyan),var(--blue),var(--orange))}
.flow-node{position:relative;padding:0 9px;text-align:center}.flow-node .dot{width:56px;height:56px;margin:0 auto 16px;display:grid;place-items:center;background:var(--bg);border:1px solid #3b596d;color:var(--cyan);font:bold 14px monospace;position:relative;z-index:1}.flow-node:nth-child(n+5) .dot{border-color:#8c4937;color:var(--orange)}.flow-node b{display:block;font-size:14px;margin-bottom:8px}.flow-node span{display:block;color:var(--muted);font-size:12px;line-height:1.45}
.loop-note{margin:34px auto 0;max-width:980px;padding:17px 22px;border:1px solid #294052;text-align:center;color:#c2d0da;font-size:14px}.loop-note b{color:var(--orange)}
.map-shell{position:relative;min-height:600px;border:1px solid rgba(255,255,255,.1);background:#09141d;overflow:hidden}.map-shell canvas{display:block;width:100%;height:600px}.map-toolbar{position:absolute;left:14px;top:14px;display:flex;flex-wrap:wrap;gap:7px;max-width:75%}.chip{padding:8px 10px}.chip[aria-pressed="true"]{border-color:var(--cyan);color:var(--cyan);background:rgba(81,210,215,.09)}
.map-legend{position:absolute;right:14px;bottom:14px;padding:13px 15px;background:rgba(7,16,24,.86);border:1px solid rgba(255,255,255,.11);font-size:11px;color:var(--muted);line-height:1.8}.legend-row{display:flex;gap:8px;align-items:center}.swatch{width:9px;height:9px;display:inline-block}.map-note{position:absolute;left:14px;bottom:14px;max-width:50%;font-size:11px;line-height:1.5;color:#9fb0bd;background:rgba(7,16,24,.82);padding:9px 11px}
.stat-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:11px}.stat{padding:18px;border-top:2px solid var(--line);background:rgba(255,255,255,.025)}.stat.orange{border-color:var(--orange)}.stat.cyanline{border-color:var(--cyan)}.stat .num{font:700 clamp(31px,4vw,58px)/1 "Arial Narrow",sans-serif;letter-spacing:-.03em}.stat .unit{font-size:14px;color:var(--muted);margin-left:4px}.stat small{display:block;margin-top:8px;color:var(--muted);line-height:1.45}
.axis-pair{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:28px}.axis-card{padding:22px;border:1px solid var(--line);background:rgba(255,255,255,.02)}.axis-card .symbol{font:700 36px/1 monospace;margin-bottom:12px}.axis-card.point .symbol{color:var(--cyan)}.axis-card.region .symbol{color:var(--orange)}.axis-card p{color:var(--muted);line-height:1.6;font-size:14px}.axis-card b{color:#fff}
.layer-stack{display:grid;gap:10px;perspective:900px}.layer{display:grid;grid-template-columns:58px 1fr auto;align-items:center;gap:14px;padding:19px 20px;border:1px solid var(--line);background:rgba(13,24,34,.9);transform:rotateX(2deg) translateX(var(--shift));box-shadow:0 14px 30px rgba(0,0,0,.16)}.layer:nth-child(1){--shift:0px}.layer:nth-child(2){--shift:12px}.layer:nth-child(3){--shift:24px}.layer:nth-child(4){--shift:36px}.layer .n{font:700 22px monospace;color:var(--cyan)}.layer:nth-child(4) .n{color:var(--orange)}.layer b{display:block;margin-bottom:4px}.layer small{color:var(--muted);line-height:1.45}.layer .status{font-size:11px;color:var(--muted);text-align:right}
.formula{font:500 clamp(22px,2.3vw,38px)/1.5 "Cambria Math","Times New Roman",serif;padding:25px 28px;border-left:3px solid var(--orange);background:rgba(255,255,255,.026);margin:28px 0}.formula .sub{font-size:.62em;vertical-align:sub}.formula .sup{font-size:.62em;vertical-align:super}.equation-note{font-size:13px;color:var(--muted);line-height:1.65}
.algo{display:grid;grid-template-columns:repeat(6,1fr);gap:9px;margin-top:34px}.algo-step{min-height:185px;padding:17px 15px;border:1px solid var(--line);position:relative;background:#0a151f}.algo-step::after{content:"→";position:absolute;right:-15px;top:44%;z-index:3;color:#557181}.algo-step:last-child::after{display:none}.algo-step .k{font:700 11px monospace;color:var(--cyan);margin-bottom:30px}.algo-step b{display:block;font-size:15px;margin-bottom:10px}.algo-step p{font-size:12px;line-height:1.55;color:var(--muted)}
.ledger{margin-top:28px;border-top:1px solid var(--line)}.ledger-row{display:grid;grid-template-columns:110px 1.2fr 2fr;gap:16px;padding:12px 6px;border-bottom:1px solid rgba(255,255,255,.07);font-size:13px;align-items:center}.ledger-row .group{color:var(--cyan);font:700 11px monospace}.ledger-row span:last-child{color:var(--muted)}
.perf-grid{display:grid;grid-template-columns:1fr 1fr;gap:28px}.perf-panel{padding:24px;border:1px solid var(--line);background:rgba(255,255,255,.02)}.perf-panel h3{display:flex;justify-content:space-between;align-items:baseline}.perf-panel h3 small{font-size:11px;color:var(--muted);font-weight:500}.bar-row{display:grid;grid-template-columns:118px 1fr 68px;gap:10px;align-items:center;margin:14px 0;font-size:12px}.bar-track{height:12px;background:#15232e;position:relative}.bar-fill{height:100%;background:var(--blue);width:var(--w)}.bar-fill.orange{background:var(--orange)}.bar-fill.cyanbar{background:var(--cyan)}.bar-value{text-align:right;font:700 12px monospace}.gate{margin-top:26px;display:grid;grid-template-columns:auto 1fr;gap:17px;padding:20px;border:1px solid rgba(255,112,72,.35);background:rgba(255,112,72,.065)}.gate .seal{width:68px;height:68px;display:grid;place-items:center;border:2px solid var(--orange);color:var(--orange);font:900 14px monospace;transform:rotate(-5deg)}.gate p{margin:0;color:#dfbdb4;line-height:1.58;font-size:14px}
.feedback{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;align-items:center;margin-top:45px}.feedback-card{min-height:170px;padding:20px;border:1px solid var(--line);background:#0a151f;position:relative}.feedback-card:not(:last-child)::after{content:"→";position:absolute;right:-18px;top:46%;z-index:3;color:var(--orange);font-size:22px}.feedback-card .k{color:var(--cyan);font:700 11px monospace;margin-bottom:22px}.feedback-card b{display:block;margin-bottom:9px}.feedback-card p{font-size:12px;line-height:1.55;color:var(--muted)}
.case{padding:25px 24px;border-top:3px solid var(--blue);background:rgba(255,255,255,.025);min-height:340px;display:flex;flex-direction:column}.case.jp{border-color:var(--orange)}.case.uk{border-color:var(--cyan)}.case.us{border-color:var(--gold)}.case .country{font:700 12px monospace;color:var(--muted);letter-spacing:.12em;margin-bottom:25px}.case p{color:#b8c6d0;font-size:14px;line-height:1.65}.case .verdict{margin-top:auto;padding-top:18px;border-top:1px solid var(--line);font-weight:700;font-size:13px}.case a{margin-top:12px;font-size:11px;color:var(--cyan);text-underline-offset:3px}
.usecase{padding:23px;border:1px solid var(--line);min-height:220px;background:rgba(255,255,255,.018)}.usecase .icon{font:700 13px monospace;color:var(--cyan);margin-bottom:36px}.usecase h3{font-size:19px}.usecase p{font-size:13px;color:var(--muted);line-height:1.55}.usecase .when{display:inline-block;margin-top:15px;color:var(--orange);font-size:11px;font-weight:800}
.roadmap{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:35px}.phase{padding:24px;border:1px solid var(--line);min-height:280px;background:rgba(255,255,255,.018)}.phase .year{font:700 12px monospace;color:var(--cyan)}.phase h3{margin:22px 0 17px}.phase ul{padding-left:18px;margin:0;color:var(--muted);font-size:13px;line-height:1.75}.phase:last-child{border-color:#684235}.phase:last-child .year{color:var(--orange)}
.decision-strip{display:grid;grid-template-columns:repeat(3,1fr);gap:1px;background:var(--line);margin-top:20px;border:1px solid var(--line)}.decision{background:#0a151f;padding:16px 18px;font-size:13px;line-height:1.55}.decision b{color:var(--orange)}
.closing{background:radial-gradient(circle at 50% 60%,rgba(81,210,215,.08),transparent 35%),#071018;text-align:center}.closing h2{font-size:clamp(42px,6.8vw,100px);margin-inline:auto}.closing .lead{margin-inline:auto}.link-row{display:flex;justify-content:center;flex-wrap:wrap;gap:10px;margin-top:35px}.primary-link{display:inline-flex;padding:13px 18px;border:1px solid var(--cyan);color:var(--cyan);text-decoration:none;font-size:13px;font-weight:800}.primary-link.orange-link{border-color:var(--orange);color:var(--orange)}
.source-grid{display:grid;grid-template-columns:1.1fr .9fr;gap:40px}.sources-list{display:grid;gap:10px}.source-item{padding:15px 17px;border:1px solid var(--line);font-size:13px;line-height:1.55}.source-item b{display:block;margin-bottom:4px}.source-item span{color:var(--muted)}.source-item a{color:var(--cyan);word-break:break-all;font-size:11px}.footnote{font-size:11px;line-height:1.65;color:var(--dim)}
@media(max-width:1050px){.grid4{grid-template-columns:repeat(2,1fr)}.flow-track{grid-template-columns:repeat(4,1fr);gap:22px}.flow-track::before{display:none}.algo{grid-template-columns:repeat(3,1fr)}.algo-step:nth-child(3)::after,.algo-step:last-child::after{display:none}.feedback{grid-template-columns:repeat(2,1fr)}.feedback-card::after{display:none!important}}
@media(max-width:780px){html{scroll-snap-type:y proximity}.topbar{height:52px}.scene{align-items:flex-start;padding-top:82px;overflow:visible}.scene-nav{display:none}.grid2,.grid3,.perf-grid,.source-grid,.roadmap{grid-template-columns:1fr}.grid4{grid-template-columns:1fr}.flow-track{grid-template-columns:repeat(2,1fr)}.algo{grid-template-columns:1fr}.algo-step::after{display:none}.feedback{grid-template-columns:1fr}.map-shell{min-height:480px}.map-shell canvas{height:480px}.map-note{max-width:72%}.controls{bottom:12px}.controls button{padding:8px 11px}.ledger-row{grid-template-columns:72px 1fr}.ledger-row span:last-child{grid-column:2}.decision-strip{grid-template-columns:1fr}.hero-foot{align-items:flex-start;flex-direction:column}.layer{grid-template-columns:45px 1fr}.layer .status{display:none}.bar-row{grid-template-columns:90px 1fr 58px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}.reveal{opacity:1;transform:none;transition:none}.scroll-cue i::after{animation:none}}
</style>
</head>
<body>
<header class="topbar"><div class="brand"><b>KOREA</b> LMM / DECISION BRIEF</div><div class="progress"><span id="progress"></span></div><div class="counter mono" id="counter">01 / 13</div></header>
<nav class="scene-nav" id="sceneNav" aria-label="장면 이동"></nav>
<div class="controls"><button id="prev" aria-label="이전 장면">↑ 이전</button><button id="next" aria-label="다음 장면">다음 ↓</button></div>

<main>
<section class="scene hero" data-no="01" data-title="표지"><canvas id="fieldCanvas" aria-hidden="true"></canvas><div class="inner hero-copy">
  <div class="kicker reveal"><span class="tag">전문가·발주자 브리핑</span><span class="tag">단일 파일 · 오프라인</span><span class="tag" id="modelStamp">MODEL —</span></div>
  <h1 class="reveal" style="--delay:80ms"><span>30개의 점을</span><span class="ghost">하나의 국가 자기장 기준으로.</span></h1>
  <p class="lead reveal" style="--delay:160ms">선점의 목표는 점을 늘리는 일이 아닙니다. <strong>전국 어디서나 설명 가능한 자기장</strong>을 만들고, 설명하지 못하는 곳을 다음 관측으로 되돌려 보내는 일입니다.</p>
  <div class="hero-foot reveal" style="--delay:240ms"><div class="scroll-cue"><i></i><span>스크롤 · 화살표 키로 진행</span></div><span class="mono">LMM = LOCAL MAGNETIC MODEL</span></div>
</div></section>

<section class="scene" data-no="02" data-title="전체 흐름"><div class="inner">
  <div class="eyebrow reveal">THE WHOLE LOOP</div><h2 class="reveal">점에서 모델로, 모델에서 다시 점으로.</h2>
  <p class="lead reveal">이 사업은 후보지를 한 번 고르고 끝나는 선점이 아니라, <strong>관측망과 모델을 함께 개선하는 순환 절차</strong>입니다.</p>
  <div class="flow-track reveal" style="--delay:120ms">
    <div class="flow-node"><div class="dot">01</div><b>자료 감사</b><span>좌표·이설·야장·시각·장비 이력</span></div>
    <div class="flow-node"><div class="dot">02</div><b>기존망 진단</b><span>공백·담당면적·재방문 품질</span></div>
    <div class="flow-node"><div class="dot">03</div><b>도상 후보</b><span>제외구역·환경·저구배 입지</span></div>
    <div class="flow-node"><div class="dot">04</div><b>현장 검증</b><span>103개소 접근·교란·방위표지</span></div>
    <div class="flow-node"><div class="dot">05</div><b>절대 관측</b><span>D·I·F·시간 동시 관측</span></div>
    <div class="flow-node"><div class="dot">06</div><b>LMM 결합</b><span>Core + Regional + Crustal</span></div>
    <div class="flow-node"><div class="dot">07</div><b>다음 보강</b><span>잔차·불확실성 큰 권역 재설계</span></div>
  </div>
  <div class="loop-note reveal" style="--delay:210ms">핵심 산출물은 후보점 목록이 아니라 <b>검증 가능한 관측–모델–보강의 폐루프</b>입니다.</div>
</div></section>

<section class="scene" data-no="03" data-title="왜 추가 선점인가"><div class="inner grid2">
  <div>
    <div class="eyebrow reveal">WHY MORE SITES</div><h2 class="reveal">30점이 있어도, 모델에 바로 쓸 수 있는 점은 16점입니다.</h2>
    <p class="lead reveal">위치는 존재하지만 모든 성과가 같은 시기·품질·표석 이력을 갖지는 않습니다. <strong>관측망의 점 수</strong>와 <strong>현재 모델의 유효 입력 수</strong>를 구분해야 합니다.</p>
    <div class="stat-grid reveal" style="--delay:120ms">
      <div class="stat cyanline"><span class="num" data-bind="network_sites">30</span><span class="unit">점</span><small>1등 지자기점 공간망<br>좌표 대조표 출처는 추가 확인 필요</small></div>
      <div class="stat orange"><span class="num" data-bind="model_sites">16</span><span class="unit">점</span><small>현재 LMM 품질감사 통과 입력<br>2019·2022–2025</small></div>
      <div class="stat"><span class="num" data-bind="survey_sites">103</span><span class="unit">개소</span><small>현장조사 취합본<br>A 6 · B 35 · C 49 · D 13</small></div>
      <div class="stat"><span class="num" data-bind="candidate_sites">178</span><span class="unit">후보</span><small>도상선점 후보<br>P1 62 · P2 57 · P3 59</small></div>
    </div>
    <div class="warning reveal" style="--delay:190ms;margin-top:18px">목표 80점 = 기존 30 + 신규 50은 <strong>설계 시나리오</strong>입니다. 승인된 최종 점수나 해외 표준이 아닙니다.</div>
  </div>
  <div class="map-shell reveal" style="--delay:100ms"><canvas id="networkMap" aria-label="기존망·모델 입력·현장조사 지점 지도"></canvas><div class="map-toolbar">
    <button class="chip" data-map="network" data-layer="existing" aria-pressed="true">기존망 30</button><button class="chip" data-map="network" data-layer="model" aria-pressed="true">LMM 16</button><button class="chip" data-map="network" data-layer="survey" aria-pressed="false">현장 103</button><button class="chip" data-map="network" data-layer="density" aria-pressed="false">상대 부족도</button>
  </div><div class="map-legend"><div class="legend-row"><i class="swatch" style="background:#eaf1f6"></i>기존망</div><div class="legend-row"><i class="swatch" style="background:#ff7048"></i>LMM 입력</div><div class="legend-row"><i class="swatch" style="background:#51d2d7"></i>현장조사</div></div><div class="map-note">버튼으로 레이어를 켜고 끌 수 있습니다.</div></div>
</div></section>

<section class="scene" data-no="04" data-title="적응형 선점"><div class="inner grid2">
  <div class="map-shell reveal"><canvas id="densityMap" aria-label="항공자력 복잡도 대비 상대 부족도 지도"></canvas><div class="map-legend"><div class="legend-row"><i class="swatch" style="background:#ff7048"></i>보강 1순위</div><div class="legend-row"><i class="swatch" style="background:#f0bc54"></i>평균보다 부족</div><div class="legend-row"><i class="swatch" style="background:#427df4"></i>평균보다 양호</div><div class="legend-row"><i class="swatch" style="background:#244257"></i>충분</div></div><div class="map-note"><span data-bind="land_cells">1,091</span> 국토셀 · 보강 1순위 <span data-bind="critical_cells">76</span> · σ 대체 <span data-bind="sigma_missing">57</span></div></div>
  <div>
    <div class="eyebrow reveal">TWO SCALES, ONE RULE</div><h2 class="reveal">복잡한 권역에 더 두고, 각 점은 조용한 곳에 둡니다.</h2>
    <p class="lead reveal">항공자력은 한 점의 적합성과 한 권역의 필요한 밀도를 서로 다른 척도로 판단합니다.</p>
    <div class="axis-pair reveal" style="--delay:120ms">
      <article class="axis-card point"><div class="symbol">|∇ΔT|</div><h3>점 축 · 국소 구배</h3><p><span class="mono">s5 = f(|∇ΔT|)</span>. 표석은 주변 자기장이 급변하지 않는 <b>저구배 지점</b>을 우대합니다. 현장 정밀 자력측량이 최종 확인합니다.</p></article>
      <article class="axis-card region"><div class="symbol">σ<sub>25 km</sub></div><h3>권역 축 · 공간 복잡도</h3><p><span class="mono">n<sub>h</sub> ∝ A<sub>h</sub>·σ<sub>h</sub>, L<sub>h</sub> ∝ 1/√σ<sub>h</sub></span>. 복잡한 권역은 같은 오차를 내려면 <b>더 촘촘한 관측</b>이 필요하다고 가정합니다.</p></article>
    </div>
    <div class="warning reveal" style="--delay:210ms;margin-top:16px">Neyman 배분과 s5 저구배 곡선은 이 연구의 <strong>잠정 설계안</strong>입니다. σ 결측 셀과 소표본 검증 한계가 있어 최종 규칙으로 확정하지 않았습니다.</div>
  </div>
</div></section>

<section class="scene" data-no="05" data-title="LMM의 역할"><div class="inner grid2">
  <div>
    <div class="eyebrow reveal">WHY LMM</div><h2 class="reveal">관측점은 값 하나를 주고, LMM은 전국의 값을 설명합니다.</h2>
    <p class="lead reveal"><strong>LMM(Local Magnetic Model)</strong>은 전 지구 기준장에 한반도 관측·지각 이상을 결합해 위치와 시점별 편각 D, 복각 I, 총자력 F를 추정하는 지역 모델입니다.</p>
    <div class="formula reveal" style="--delay:120ms">B<span class="sub">obs, quiet</span> = B<span class="sub">obs</span> − δB<span class="sub">external</span><br> B̂<span class="sub">LMM</span>(r,t) = B<span class="sub">IGRF-14</span> + R<span class="sub">Korea</span> + C<span class="sub">KIGAM</span></div>
    <p class="equation-note reveal" style="--delay:180ms">현재판은 외부장을 미래 시점 예측에 더하지 않습니다. 네 관측소 자료로 <strong>관측 세션의 D·I 시간변화만 환산</strong>한 뒤, Core + Regional + Crustal의 baseline을 만듭니다. 현재 <span class="mono">subtract_DI</span>는 QC=QUIET만 고르는 엄격 필터가 아닙니다.</p>
  </div>
  <div class="layer-stack reveal" style="--delay:80ms">
    <div class="layer"><span class="n">01</span><div><b>Core · 전 지구 핵 기원장</b><small>IGRF-14, 구면조화 degree 13</small></div><span class="status">기준장</span></div>
    <div class="layer"><span class="n">02</span><div><b>Regional · 한반도 장파장 잔차</b><small>16측점 · LOO 선택 결과 현재 0차(상수)</small></div><span class="status">D·I·F</span></div>
    <div class="layer"><span class="n">03</span><div><b>Crustal · 지각 단파장 이상</b><small>KIGAM 1.5분 격자 · 현재 F에만 스칼라 적용</small></div><span class="status">F only</span></div>
    <div class="layer"><span class="n">04</span><div><b>External · 시간 변화 환산</b><small>CYG·제주·강릉·이천 1분 자료 · D·I 전처리</small></div><span class="status">전처리</span></div>
  </div>
</div></section>

<section class="scene" data-no="06" data-title="데이터와 알고리듬"><div class="inner">
  <div class="eyebrow reveal">DATA → ALGORITHM → MODEL</div><h2 class="reveal">수식보다 중요한 것은, 어떤 오차를 어느 단계에서 제거했는가입니다.</h2>
  <div class="algo reveal" style="--delay:100ms">
    <div class="algo-step"><div class="k">01 / AUDIT</div><b>성과 감사</b><p>좌표·이설·야장·관측시각·재방문 산포로 입력 신뢰도를 분리합니다.</p></div>
    <div class="algo-step"><div class="k">02 / REDUCE</div><b>외부장 환산</b><p>네 관측소의 시간 편차를 공간보간해 세션 D·I에서 뺍니다.</p></div>
    <div class="algo-step"><div class="k">03 / RESIDUAL</div><b>IGRF 잔차</b><p>ΔD, ΔI, ΔF = 관측 − IGRF-14를 방문별로 계산합니다.</p></div>
    <div class="algo-step"><div class="k">04 / ROBUST</div><b>강건 선별</b><p>MAD 2.5σ 반복 선별로 블런더·국지 자성체 영향을 제한합니다.</p></div>
    <div class="algo-step"><div class="k">05 / FIT</div><b>층별 적합</b><p>F에는 KIGAM 이상을 1:1 결합하고 남은 장파장 잔차를 적합합니다.</p></div>
    <div class="algo-step"><div class="k">06 / LOO</div><b>독립 예측</b><p>한 점씩 빼고 예측해 차수와 성능을 평가합니다. 현재 최적은 0차입니다.</p></div>
  </div>
  <div class="ledger reveal" id="ledger" style="--delay:190ms"></div>
</div></section>

<section class="scene" data-no="07" data-title="현재 성능"><div class="inner">
  <div class="eyebrow reveal">QUALITY GATE</div><h2 class="reveal">지각층은 F를 개선했지만, 현재 LMM은 아직 정식 편각 산출 단계가 아닙니다.</h2>
  <div class="perf-grid reveal" style="--delay:100ms"><div class="perf-panel"><h3>편각 D <small>RMS · degree</small></h3><div id="perfD"></div></div><div class="perf-panel"><h3>총자력 F <small>RMS · nT</small></h3><div id="perfF"></div></div></div>
  <div class="gate reveal" style="--delay:180ms"><div class="seal">HOLD</div><p><strong>LOO D <span data-bind="loo_D">0.499</span>°</strong>는 공학 KPI 0.1°의 약 <span data-bind="loo_D_ratio">5.0</span>배, <strong>LOO F <span data-bind="loo_F">62.8</span> nT</strong>는 KPI 50 nT의 약 <span data-bind="loo_F_ratio">1.26</span>배입니다. 현재판은 진단·설계 참고용이며, 품질 게이트 통과 전 국가기본도 자편각 등 정식 성과에 쓰지 않습니다. KPI와 법정 관측 정수차 30′는 서로 다른 기준입니다.</p></div>
  <p class="footnote reveal" style="--delay:240ms;margin-top:15px">표본내 단계 RMS: D IGRF → Regional, F IGRF → Crustal → Regional. LOO는 적합에 쓰지 않은 점을 예측한 값이므로 최종 판정은 LOO를 우선합니다. <span class="mono">lmm_model.json</span>의 External 설명은 2026-08-27에 실제 <span class="mono">subtract_DI</span> 구현(D·I 전처리 보정 · F 미보정 · strict quiet 필터 없음 · 예측 시 미평가)에 맞춰 정정했습니다.</p>
</div></section>

<section class="scene" data-no="08" data-title="LMM과 선점의 되먹임"><div class="inner">
  <div class="eyebrow reveal">THE FEEDBACK LOOP</div><h2 class="reveal">LMM은 첫 후보를 고르는 도구가 아니라, 다음 관측의 가치를 계산하는 도구입니다.</h2>
  <p class="lead reveal">초기 후보는 공간 공백·간섭원·저구배로 고릅니다. 관측 뒤에는 LMM 잔차와 불확실성을 더해 다음 보강 권역을 좁힙니다. 이 구분이 순환논리를 막습니다.</p>
  <div class="feedback reveal" style="--delay:120ms">
    <div class="feedback-card"><div class="k">A / BASE</div><b>기본 공간균형</b><p>도엽·담당면적·접근 가능성으로 전국 공백을 먼저 줄입니다.</p></div>
    <div class="feedback-card"><div class="k">B / REGION</div><b>복잡권역 보강</b><p>σ(25 km)가 큰 권역은 더 촘촘한 설계를 검토합니다.</p></div>
    <div class="feedback-card"><div class="k">C / POINT</div><b>저구배 점 선택</b><p>각 권역 안에서는 |∇ΔT|가 낮고 인공 교란이 적은 곳을 고릅니다.</p></div>
    <div class="feedback-card"><div class="k">D / MODEL</div><b>LMM 재적합</b><p>새 D·I·F를 넣고 독립 검증으로 실제 개선을 확인합니다.</p></div>
    <div class="feedback-card"><div class="k">E / NEXT</div><b>잔차 기반 다음 점</b><p>모델이 설명하지 못한 권역만 다시 보강합니다.</p></div>
  </div>
  <div class="loop-note reveal" style="--delay:220ms">추가점의 가치는 <b>점수</b>가 아니라 “그 점을 넣었을 때 독립 예측오차가 얼마나 줄었는가”로 최종 검증합니다.</div>
</div></section>

<section class="scene" data-no="09" data-title="해외 근거"><div class="inner">
  <div class="eyebrow reveal">FOREIGN EVIDENCE, PRECISELY</div><h2 class="reveal">해외는 답을 복사해 주지 않았지만, 설계 원리를 검증해 줍니다.</h2>
  <div class="grid3 reveal" style="--delay:100ms">
    <article class="case jp"><div class="country">JAPAN / GSI</div><h3>구배 대비 밀도가 낮으면 모델 오차가 커집니다.</h3><p>평균 약 20 km 지상망은 더 짧은 파장의 이상을 표현하지 못해 항공자력으로 보완했습니다. 야쓰가타케·홋카이도 동부의 큰 오차는 큰 자기구배에 비해 1·2등 자기점 밀도가 낮은 데서 왔다고 분석했습니다.</p><div class="verdict">우리 설계의 가장 강한 실증 근거<br><span class="muted">단, 자동 밀도배분 선례는 아님</span></div><a href="https://www.gsi.go.jp/common/000071404.pdf" target="_blank" rel="noreferrer">GSI 2010.0 모델 정확도 평가 ↗</a><a href="https://www.gsi.go.jp/common/000024739.pdf" target="_blank" rel="noreferrer">GSI 항공자력·20 km 지상망 설명 ↗</a></article>
    <article class="case uk"><div class="country">UNITED KINGDOM / BGS</div><h3>시간 변화와 지각 이상은 역할을 나눕니다.</h3><p>41개 반복점과 3개 관측소로 UK Regional Geomagnetic Model을 갱신하며, 반복점은 자기 간섭이 적은 곳을 택합니다. 반복점은 주로 core field의 시간변화, 항공자력은 crustal field를 담당합니다.</p><div class="verdict">현재 LMM의 Regional + Crustal 분리와 정합</div><a href="https://www.geomag.bgs.ac.uk/operations/uksurvey.html" target="_blank" rel="noreferrer">BGS UK Magnetic Survey ↗</a><a href="https://geomag.bgs.ac.uk/education/earthmag.html" target="_blank" rel="noreferrer">BGS 자료 역할 설명 ↗</a></article>
    <article class="case us"><div class="country">US / UK JOINT MODEL</div><h3>관측 공백은 다른 자료로 보완하되, 역할을 혼동하지 않습니다.</h3><p>1995 개정 모델은 관측소·반복점을 주 자료로 쓰고, 이 자료가 없는 영역에서 항공자력·위성 자료로 영년변화 정보를 보완했습니다.</p><div class="verdict">공백 보완의 개념적 선례<br><span class="muted">신규 반복점 밀도 알고리듬은 아님</span></div><a href="https://pubs.usgs.gov/publication/70019514" target="_blank" rel="noreferrer">USGS Publications Warehouse ↗</a></article>
  </div>
</div></section>

<section class="scene" data-no="10" data-title="최종 활용"><div class="inner">
  <div class="eyebrow reveal">FROM MODEL TO PUBLIC VALUE</div><h2 class="reveal">최종 성과는 지도 한 장이 아니라, 위치·시점별 자기장 기준 서비스입니다.</h2>
  <div class="grid4 reveal" style="--delay:100ms">
    <article class="usecase"><div class="icon">01 / MAP</div><h3>국가기본도 자편각 갱신</h3><p>도폭과 기준년별 편각·변화율을 일관된 모델에서 산출하고, 현장 측량의 진북–자북 관계를 제공합니다.</p><span class="when">KPI 통과 후 공식 활용</span></article>
    <article class="usecase"><div class="icon">02 / FIELD</div><h3>측량·관측 품질관리</h3><p>관측값과 예상값의 차이로 방위표지 오류, 이설, 국지 교란, 시각 보정 누락을 조기에 탐지합니다.</p><span class="when">현재도 진단 활용 가능</span></article>
    <article class="usecase"><div class="icon">03 / NAV</div><h3>항법·헤딩 보정</h3><p>차량·드론·해양·자기센서의 지역 편차를 보정하는 고해상도 기준으로 확장할 수 있습니다.</p><span class="when">독립 검증·배포 체계 필요</span></article>
    <article class="usecase"><div class="icon">04 / GEO</div><h3>지질·인프라 기준장</h3><p>항공·지상 자력탐사, 방향성 시추, 지하 구조 해석에서 core·crustal·external 성분을 분리하는 기준이 됩니다.</p><span class="when">목적별 해상도 추가 필요</span></article>
  </div>
  <p class="footnote reveal" style="--delay:190ms;margin-top:22px">BGS는 지역모델을 지도 magnetic north 갱신에 사용하고, WMM/BGGM 계열은 항법·헤딩·방향성 시추 등에 활용합니다. 국내 적용은 정확도·책임주체·갱신주기를 별도 확정해야 합니다.</p>
</div></section>

<section class="scene" data-no="11" data-title="의사결정 로드맵"><div class="inner">
  <div class="eyebrow reveal">DECISION ROADMAP</div><h2 class="reveal">다음 예산은 점 수가 아니라, 품질 게이트를 통과시키는 데 배분해야 합니다.</h2>
  <div class="roadmap reveal" style="--delay:100ms">
    <article class="phase"><div class="year">GATE 01 / INPUT</div><h3>자료 계보를 닫습니다.</h3><ul><li>30점 좌표 대조표의 원본파일명·작성일·문서번호 확인</li><li>표석 이설·개명은 지점코드로 추적</li><li>모든 관측에 시각·장비·방위표지 메타데이터 결합</li></ul></article>
    <article class="phase"><div class="year">GATE 02 / CAMPAIGN</div><h3>가치가 큰 점부터 관측합니다.</h3><ul><li>기본 공간 공백 + 항공자력 복잡권역을 단계 배치</li><li>각 점은 저구배·저인공교란 현장 검증</li><li>신규점 수는 순차 LOO 개선량으로 재결정</li></ul></article>
    <article class="phase"><div class="year">GATE 03 / RELEASE</div><h3>독립 검증 뒤 공개합니다.</h3><ul><li>외부장 D·I·F 전 성분 환산 체계 확정</li><li>항공자력 결측·벡터 지각장 개선</li><li>D&lt;0.1°·F&lt;50 nT 독립 검증 통과 후 공식 서비스</li></ul></article>
  </div>
  <div class="decision-strip reveal" style="--delay:190ms"><div class="decision"><b>결정 1</b><br>80점을 확정할 것인가가 아니라, 적응형 보강 원칙을 승인할 것인가</div><div class="decision"><b>결정 2</b><br>관측시각·상시관측소 연계를 필수 성과로 둘 것인가</div><div class="decision"><b>결정 3</b><br>고복잡도 권역의 정밀 항공/지상 자력자료를 추가 확보할 것인가</div></div>
</div></section>

<section class="scene closing" data-no="12" data-title="결론"><div class="inner">
  <div class="eyebrow reveal" style="justify-content:center">THE ONE SENTENCE</div><h2 class="reveal">좋은 선점은 조용한 점을 찾고,<br><span class="accent">좋은 모델은 설명 못한 곳을 다시 찾습니다.</span></h2>
  <p class="lead reveal">관측망·항공자력·상시관측·IGRF를 하나의 검증 순환으로 묶을 때, 신규점은 비용이 아니라 국가 자기장 기준의 해상도가 됩니다.</p>
  <div class="link-row reveal" style="--delay:140ms"><a class="primary-link" href="index.html">입지 선정 지도 열기</a><a class="primary-link orange-link" href="lmm.html">LMM 계산기 열기</a></div>
  <p class="footnote reveal" style="--delay:210ms;margin:28px auto 0;max-width:760px">현재 모델은 시범 구축판이며 품질 게이트 미통과입니다. 이 브리핑은 의사결정과 연구설계를 위한 자료입니다.</p>
</div></section>

<section class="scene" data-no="13" data-title="근거와 한계"><div class="inner source-grid">
  <div><div class="eyebrow reveal">PROVENANCE</div><h2 class="reveal">수치의 출처와 한계를 같은 화면에 둡니다.</h2><div class="sources-list reveal" id="sourceLedger" style="--delay:100ms"></div></div>
  <div class="reveal" style="--delay:140ms"><h3>공식 해외 근거</h3><div class="sources-list">
    <div class="source-item"><b>GSI · 지자기 시공간모델 정확도 평가</b><span>큰 자기구배 대비 낮은 측점밀도와 모델오차의 관계</span><br><a href="https://www.gsi.go.jp/common/000071404.pdf" target="_blank" rel="noreferrer">gsi.go.jp/common/000071404.pdf</a></div>
    <div class="source-item"><b>GSI · 항공자력과 지상망 공간해상도</b><span>약 20 km 지상망이 표현하지 못하는 단파장 이상 보완</span><br><a href="https://www.gsi.go.jp/common/000024739.pdf" target="_blank" rel="noreferrer">gsi.go.jp/common/000024739.pdf</a></div>
    <div class="source-item"><b>BGS · UK Magnetic Survey / Earth magnetic field</b><span>반복점·관측소와 항공자력의 역할 분리, 지도·모델 활용</span><br><a href="https://www.geomag.bgs.ac.uk/operations/uksurvey.html" target="_blank" rel="noreferrer">geomag.bgs.ac.uk/operations/uksurvey.html</a></div>
    <div class="source-item"><b>USGS · 1995 joint US/UK model</b><span>관측소·반복점 공백에서 항공·위성자료 보완</span><br><a href="https://pubs.usgs.gov/publication/70019514" target="_blank" rel="noreferrer">pubs.usgs.gov/publication/70019514</a></div>
  </div><p class="footnote" style="margin-top:18px">해외 자료는 설계 원리의 근거입니다. 항공자력 구배로 반복점 밀도를 자동 배분한 직접 선례로 인용하지 않습니다.</p></div>
</div></section>
</main>

<script>
const DATA=__PAYLOAD__;
const scenes=[...document.querySelectorAll('.scene')];
const fmt=(v,d=0)=>Number(v).toLocaleString('ko-KR',{minimumFractionDigits:d,maximumFractionDigits:d});

function bindData(){
  const S=DATA.stats,M=DATA.model;
  const values={...S,critical_cells:S.density_counts.critical,loo_D:fmt(M.loo.D,3),loo_F:fmt(M.loo.F,1),loo_D_ratio:fmt(M.loo.D/.1,1),loo_F_ratio:fmt(M.loo.F/50,2)};
  document.querySelectorAll('[data-bind]').forEach(el=>{const k=el.dataset.bind;if(k in values){const v=values[k];el.textContent=typeof v==='number'?fmt(v,0):v;}});
  document.querySelector('#modelStamp').textContent=`MODEL ${DATA.generated.slice(0,10)} · EPOCH ${DATA.epoch.toFixed(1)}`;
  const ledger=document.querySelector('#ledger');
  DATA.sources.forEach(s=>{ledger.insertAdjacentHTML('beforeend',`<div class="ledger-row"><span class="group">${s.group}</span><b>${s.name}</b><span>${s.detail}</span></div>`)});
  const sourceLedger=document.querySelector('#sourceLedger');
  DATA.sources.forEach(s=>{sourceLedger.insertAdjacentHTML('beforeend',`<div class="source-item"><b>${s.group} · ${s.name}</b><span>${s.detail}</span></div>`)});
}

function perfRows(target,rows,max){
  const box=document.querySelector(target);box.innerHTML='';
  rows.forEach(r=>{const w=Math.min(100,100*r.value/max);box.insertAdjacentHTML('beforeend',`<div class="bar-row"><span>${r.label}</span><div class="bar-track"><div class="bar-fill ${r.cls||''}" style="--w:${w}%"></div></div><span class="bar-value">${r.text}</span></div>`)});
}
function buildPerformance(){const m=DATA.model;
  perfRows('#perfD',[{label:'IGRF 단독',value:m.stage.D_igrf,text:fmt(m.stage.D_igrf,3)+'°'},{label:'+ Regional',value:m.stage.D_final,text:fmt(m.stage.D_final,3)+'°',cls:'cyanbar'},{label:'LOO 예측',value:m.loo.D,text:fmt(m.loo.D,3)+'°',cls:'orange'},{label:'공학 KPI',value:.1,text:'0.100°'}],.6);
  perfRows('#perfF',[{label:'IGRF 단독',value:m.stage.F_igrf,text:fmt(m.stage.F_igrf,1)},{label:'+ Crustal',value:m.stage.F_crust,text:fmt(m.stage.F_crust,1),cls:'cyanbar'},{label:'+ Regional',value:m.stage.F_final,text:fmt(m.stage.F_final,1),cls:'cyanbar'},{label:'LOO 예측',value:m.loo.F,text:fmt(m.loo.F,1),cls:'orange'},{label:'공학 KPI',value:50,text:'50.0'}],100);
}

function drawGeometry(ctx,geom,project){ctx.beginPath();
  const ring=r=>{r.forEach((p,i)=>{const [x,y]=project(p[0],p[1]);if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)});ctx.closePath()};
  if(geom.type==='Polygon')geom.coordinates.forEach(ring);else if(geom.type==='MultiPolygon')geom.coordinates.forEach(p=>p.forEach(ring));
}
const colors={critical:'#ff7048',short:'#f0bc54',mild:'#427df4',ok:'#244257'};
function makeMap(canvas,mode){const ctx=canvas.getContext('2d');const state={existing:true,model:true,survey:false,density:mode==='density'};
  function render(){const rect=canvas.getBoundingClientRect(),dpr=Math.min(2,devicePixelRatio||1);canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);const W=rect.width,H=rect.height;ctx.clearRect(0,0,W,H);ctx.fillStyle='#09141d';ctx.fillRect(0,0,W,H);
    const b={w:124.25,e:130.15,s:32.9,n:38.9},pad=32,scale=Math.min((W-pad*2)/(b.e-b.w),(H-pad*2)/(b.n-b.s)),ox=(W-(b.e-b.w)*scale)/2,oy=(H-(b.n-b.s)*scale)/2;const project=(lon,lat)=>[ox+(lon-b.w)*scale,oy+(b.n-lat)*scale];
    if(state.density){DATA.density.forEach(c=>{const [x,y]=project(c.lon,c.lat);ctx.globalAlpha=c.observed?.62:.22;ctx.fillStyle=colors[c.tier];const s=Math.max(3.2,scale*.095);ctx.fillRect(x-s/2,y-s/2,s,s);if(!c.observed){ctx.strokeStyle='#b8c3cb';ctx.lineWidth=.7;ctx.beginPath();ctx.moveTo(x-s/2,y-s/2);ctx.lineTo(x+s/2,y+s/2);ctx.stroke()}});ctx.globalAlpha=1}
    drawGeometry(ctx,DATA.boundary,project);ctx.fillStyle='rgba(13,27,38,.46)';ctx.fill('evenodd');ctx.strokeStyle='#5d798b';ctx.lineWidth=1.1;ctx.stroke();
    if(state.survey)drawPoints(DATA.survey,'#51d2d7',2.3,.52);
    if(state.existing)drawPoints(DATA.existing,'#eaf1f6',3.2,.9,true);
    if(state.model)drawPoints(DATA.model_sites,'#ff7048',4.4,1);
    function drawPoints(rows,color,r,alpha,ring=false){ctx.globalAlpha=alpha;rows.forEach(p=>{const [x,y]=project(p.lon,p.lat);ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ring?(ctx.strokeStyle=color,ctx.lineWidth=1.4,ctx.stroke()):(ctx.fillStyle=color,ctx.fill())});ctx.globalAlpha=1}
  }
  new ResizeObserver(render).observe(canvas);render();return{state,render};
}
const maps={};
function initMaps(){maps.network=makeMap(document.querySelector('#networkMap'),'network');maps.density=makeMap(document.querySelector('#densityMap'),'density');document.querySelectorAll('[data-map]').forEach(btn=>btn.addEventListener('click',()=>{const m=maps[btn.dataset.map],k=btn.dataset.layer;m.state[k]=!m.state[k];btn.setAttribute('aria-pressed',String(m.state[k]));m.render()}));}

function initNavigation(){const nav=document.querySelector('#sceneNav');scenes.forEach((s,i)=>{const b=document.createElement('button');b.title=`${String(i+1).padStart(2,'0')} ${s.dataset.title}`;b.setAttribute('aria-label',b.title);b.addEventListener('click',()=>s.scrollIntoView({behavior:'smooth'}));nav.appendChild(b)});const dots=[...nav.children];
  const setActive=i=>{scenes.forEach((s,j)=>s.classList.toggle('active',i===j));dots.forEach((d,j)=>d.classList.toggle('active',i===j));document.querySelector('#counter').textContent=`${String(i+1).padStart(2,'0')} / ${String(scenes.length).padStart(2,'0')}`;document.querySelector('#progress').style.width=`${100*(i+1)/scenes.length}%`;window.activeScene=i};
  const io=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting&&e.intersectionRatio>.42)setActive(scenes.indexOf(e.target))}),{threshold:[.42,.65]});scenes.forEach(s=>io.observe(s));setActive(0);
  const go=d=>scenes[Math.max(0,Math.min(scenes.length-1,(window.activeScene||0)+d))].scrollIntoView({behavior:'smooth'});document.querySelector('#prev').onclick=()=>go(-1);document.querySelector('#next').onclick=()=>go(1);
  addEventListener('keydown',e=>{if(/INPUT|TEXTAREA|SELECT/.test(e.target.tagName))return;if(['ArrowDown','PageDown',' '].includes(e.key)){e.preventDefault();go(1)}else if(['ArrowUp','PageUp'].includes(e.key)){e.preventDefault();go(-1)}else if(e.key==='Home'){scenes[0].scrollIntoView()}else if(e.key==='End'){scenes.at(-1).scrollIntoView()}});
}

function initField(){const c=document.querySelector('#fieldCanvas'),ctx=c.getContext('2d');let W=0,H=0,t=0;const reduced=matchMedia('(prefers-reduced-motion:reduce)').matches;function resize(){const r=c.getBoundingClientRect(),d=Math.min(2,devicePixelRatio||1);W=r.width;H=r.height;c.width=W*d;c.height=H*d;ctx.setTransform(d,0,0,d,0,0)}new ResizeObserver(resize).observe(c);resize();
  function frame(){ctx.clearRect(0,0,W,H);const cx=W*.71,cy=H*.47;for(let i=-9;i<=9;i++){const off=i*H*.043;ctx.beginPath();for(let x=0;x<=W;x+=14){const dx=x-cx;const y=cy+off+Math.sin(dx/W*5+i*.23)*22+(dx*dx)/(W*W)*off*.65;x===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.strokeStyle=i%3===0?'rgba(81,210,215,.16)':'rgba(82,116,140,.11)';ctx.lineWidth=i===0?1.2:.7;ctx.stroke()}
    for(let i=0;i<24;i++){const q=(t*.00005+i/24)%1,x=q*W,dx=x-cx,y=cy+(i-12)*H*.032+Math.sin(dx/W*5+i*.2)*22+(dx*dx)/(W*W)*(i-12)*H*.021;ctx.fillStyle=i%5===0?'#ff7048':'rgba(81,210,215,.66)';ctx.beginPath();ctx.arc(x,y,i%5===0?2.5:1.35,0,Math.PI*2);ctx.fill()}ctx.beginPath();ctx.arc(cx,cy,5,0,Math.PI*2);ctx.fillStyle='#ff7048';ctx.fill();t+=16;if(!reduced)requestAnimationFrame(frame)}frame();}

bindData();buildPerformance();initNavigation();initMaps();initField();
</script>
</body></html>'''


def main():
    payload = build_payload()
    text = HTML.replace(
        "__PAYLOAD__", json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding="utf-8")
    print(f"[saved] {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")
    print(
        "[data] "
        f"network {payload['stats']['network_sites']} · "
        f"model {payload['stats']['model_sites']} · "
        f"survey {payload['stats']['survey_sites']} · "
        f"land cells {payload['stats']['land_cells']}"
    )


if __name__ == "__main__":
    main()
