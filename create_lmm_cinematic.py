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
THREE_RUNTIME = ROOT / "vendor" / "three.r149.min.js"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_three_runtime() -> str:
    """공식 Three.js 런타임을 단일 HTML 안에 인라인한다."""
    if not THREE_RUNTIME.exists():
        raise FileNotFoundError(
            f"Three.js runtime not found: {THREE_RUNTIME}\n"
            "Download the official r149 build before generating the cinematic."
        )
    runtime = THREE_RUNTIME.read_text(encoding="utf-8")
    required = ("SPDX-License-Identifier: MIT", 'const e="149"')
    missing = [token for token in required if token not in runtime]
    if missing:
        raise RuntimeError(
            f"Unexpected Three.js runtime ({THREE_RUNTIME}); missing {missing}"
        )
    if "</script" in runtime.lower():
        raise RuntimeError("Three.js runtime contains a closing script tag")
    return runtime


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
        "east_sea": [
            {"name": "울릉도", "lat": 37.4845, "lon": 130.9057},
            {"name": "독도", "lat": 37.2429, "lon": 131.8664},
        ],
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
                "detail": "1982–2018 · 1.5분(약 2.8 km) · F에 적용",
            },
            {
                "group": "External",
                "name": "CYG·제주·강릉·이천 1분 자료",
                "detail": "subtract_DI: 정상 세션 D·I 환산 · QUIET 전용 필터 아님 · 예측층 미적용",
            },
            {
                "group": "Network",
                "name": "1등 지자기점 30점 좌표 대조표",
                "detail": "출처 확인 필요 · 원장/정본으로 간주하지 않음",
            },
            {
                "group": "Survey",
                "name": "현장조사 일괄취합 103건",
                "detail": "20260818_211750 · A/B/C/D=6/35/49/13",
            },
        ],
    }


HTML = r'''<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<meta name="color-scheme" content="dark">
<link rel="icon" href="data:,">
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
h3{font-size:clamp(20px,1.9vw,30px);line-height:1.25;letter-spacing:-.025em}
.lead{font-size:clamp(19px,1.9vw,28px);line-height:1.62;color:#c5d1da;max-width:930px}.lead strong{color:#fff}
.tag{display:inline-flex;align-items:center;gap:8px;padding:8px 10px;border:1px solid rgba(255,255,255,.13);color:var(--muted);font-size:12px}.tag::before{content:"";width:6px;height:6px;background:var(--cyan)}
.warning{border-left:3px solid var(--orange);background:rgba(255,112,72,.08);padding:16px 20px;color:#f2c7bb;line-height:1.6;font-size:16px}
.reveal{opacity:0;transform:translateY(28px);transition:opacity .72s ease,transform .72s cubic-bezier(.2,.8,.2,1);transition-delay:var(--delay,0ms)}
.scene.active .reveal{opacity:1;transform:none}
.accent{color:var(--orange)}.cyan{color:var(--cyan)}.muted{color:var(--muted)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:clamp(22px,4vw,64px);align-items:center}.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}.grid4{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}
.hero{background:radial-gradient(circle at 72% 42%,rgba(47,125,244,.13),transparent 36%),linear-gradient(135deg,#071018 0%,#09141e 60%,#050b11 100%)}
#heroWebGL,#fieldCanvas{position:absolute;inset:0;width:100%;height:100%;z-index:0}.webgl-surface{pointer-events:none}.webgl-surface canvas{display:block;width:100%;height:100%}#fieldCanvas{opacity:.72;transition:opacity .4s ease}body.webgl-ready:not(.webgl-paused) #fieldCanvas{opacity:0}
.hero-copy{position:relative;z-index:2}.hero .kicker{display:flex;flex-wrap:wrap;gap:9px;margin-bottom:36px}
.hero h1{max-width:none;font-size:clamp(30px,6vw,96px)}.hero h1 span{display:block;white-space:nowrap}@media(max-width:620px){.hero h1 span{white-space:normal}.hero h1{font-size:clamp(34px,8.4vw,54px)}}.hero h1 .ghost{color:transparent;-webkit-text-stroke:1px rgba(234,241,246,.48)}
.hero-foot{display:flex;gap:28px;align-items:center;margin-top:44px;color:var(--muted);font-size:13px}.scroll-cue{display:flex;align-items:center;gap:10px}.scroll-cue i{display:block;width:24px;height:38px;border:1px solid #577080;border-radius:16px;position:relative}.scroll-cue i::after{content:"";position:absolute;top:7px;left:10px;width:2px;height:8px;background:var(--orange);animation:wheel 1.8s infinite}@keyframes wheel{0%{transform:translateY(0);opacity:0}30%{opacity:1}100%{transform:translateY(14px);opacity:0}}
.flow-track{display:grid;grid-template-columns:repeat(7,1fr);gap:0;margin-top:54px;position:relative}.flow-track::before{content:"";position:absolute;top:28px;left:6%;right:6%;height:1px;background:linear-gradient(90deg,var(--cyan),var(--blue),var(--orange))}
.flow-node{position:relative;padding:0 9px;text-align:center}.flow-node .dot{width:56px;height:56px;margin:0 auto 16px;display:grid;place-items:center;background:var(--bg);border:1px solid #3b596d;color:var(--cyan);font:bold 14px monospace;position:relative;z-index:1}.flow-node:nth-child(n+5) .dot{border-color:#8c4937;color:var(--orange)}.flow-node b{display:block;font-size:14px;margin-bottom:8px}.flow-node span{display:block;color:var(--muted);font-size:12px;line-height:1.45}
.loop-note{margin:34px auto 0;max-width:980px;padding:17px 22px;border:1px solid #294052;text-align:center;color:#c2d0da;font-size:16px}.loop-note b{color:var(--orange)}
.map-shell{position:relative;min-height:600px;border:1px solid rgba(255,255,255,.1);background:#09141d;overflow:hidden}.map-shell canvas{display:block;width:100%;height:600px;transition:opacity .35s ease}.map-toolbar{position:absolute;z-index:8;left:14px;top:14px;display:flex;flex-wrap:wrap;gap:7px;max-width:82%}.chip{padding:8px 10px}.chip[aria-pressed="true"]{border-color:var(--cyan);color:var(--cyan);background:rgba(81,210,215,.09)}
.map-legend{position:absolute;z-index:8;right:14px;bottom:14px;padding:13px 15px;background:rgba(7,16,24,.86);border:1px solid rgba(255,255,255,.11);font-size:11px;color:var(--muted);line-height:1.8}.legend-row{display:flex;gap:8px;align-items:center}.swatch{width:9px;height:9px;display:inline-block}.map-note{position:absolute;z-index:8;left:14px;bottom:14px;max-width:54%;font-size:13px;line-height:1.5;color:#9fb0bd;background:rgba(7,16,24,.82);padding:9px 11px}
.stat-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:11px}.stat{padding:18px;border-top:2px solid var(--line);background:rgba(255,255,255,.025)}.stat.orange{border-color:var(--orange)}.stat.cyanline{border-color:var(--cyan)}.stat .num{font:700 clamp(31px,4vw,58px)/1 "Arial Narrow",sans-serif;letter-spacing:-.03em}.stat .unit{font-size:15px;color:var(--muted);margin-left:4px}.stat small{display:block;margin-top:8px;color:var(--muted);line-height:1.5;font-size:14.5px}
.axis-pair{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:28px}.axis-card{padding:22px;border:1px solid var(--line);background:rgba(255,255,255,.02)}.axis-card .symbol{font:700 36px/1 monospace;margin-bottom:12px}.axis-card.point .symbol{color:var(--cyan)}.axis-card.region .symbol{color:var(--orange)}.axis-card p{color:var(--muted);line-height:1.68;font-size:16px}.axis-card b{color:#fff}
.layer-stack{display:grid;gap:10px;perspective:900px}.layer{display:grid;grid-template-columns:58px 1fr auto;align-items:center;gap:14px;padding:19px 20px;border:1px solid var(--line);background:rgba(13,24,34,.9);transform:rotateX(2deg) translateX(var(--shift));box-shadow:0 14px 30px rgba(0,0,0,.16)}.layer:nth-child(1){--shift:0px}.layer:nth-child(2){--shift:12px}.layer:nth-child(3){--shift:24px}.layer:nth-child(4){--shift:36px}.layer .n{font:700 22px monospace;color:var(--cyan)}.layer:nth-child(4) .n{color:var(--orange)}.layer b{display:block;margin-bottom:4px}.layer small{color:var(--muted);line-height:1.45}.layer .status{font-size:11px;color:var(--muted);text-align:right}
.formula{font:500 clamp(22px,2.3vw,38px)/1.5 "Cambria Math","Times New Roman",serif;padding:25px 28px;border-left:3px solid var(--orange);background:rgba(255,255,255,.026);margin:28px 0}.formula .sub{font-size:.62em;vertical-align:sub}.formula .sup{font-size:.62em;vertical-align:super}.equation-note{font-size:15px;color:var(--muted);line-height:1.65}
.algo{display:grid;grid-template-columns:repeat(6,1fr);gap:9px;margin-top:34px}.algo-step{min-height:185px;padding:17px 15px;border:1px solid var(--line);position:relative;background:#0a151f}.algo-step::after{content:"→";position:absolute;right:-15px;top:44%;z-index:3;color:#557181}.algo-step:last-child::after{display:none}.algo-step .k{font:700 11px monospace;color:var(--cyan);margin-bottom:30px}.algo-step b{display:block;font-size:17px;margin-bottom:10px}.algo-step p{font-size:14px;line-height:1.55;color:var(--muted)}
.ledger{margin-top:28px;border-top:1px solid var(--line)}.ledger-row{display:grid;grid-template-columns:110px 1.2fr 2fr;gap:16px;padding:12px 6px;border-bottom:1px solid rgba(255,255,255,.07);font-size:13px;align-items:center}.ledger-row .group{color:var(--cyan);font:700 11px monospace}.ledger-row span:last-child{color:var(--muted)}
.perf-grid{display:grid;grid-template-columns:1fr 1fr;gap:28px}.perf-panel{padding:24px;border:1px solid var(--line);background:rgba(255,255,255,.02)}.perf-panel h3{display:flex;justify-content:space-between;align-items:baseline}.perf-panel h3 small{font-size:11px;color:var(--muted);font-weight:500}.bar-row{display:grid;grid-template-columns:118px 1fr 68px;gap:10px;align-items:center;margin:14px 0;font-size:12px}.bar-track{height:12px;background:#15232e;position:relative}.bar-fill{height:100%;background:var(--blue);width:var(--w)}.bar-fill.orange{background:var(--orange)}.bar-fill.cyanbar{background:var(--cyan)}.bar-value{text-align:right;font:700 12px monospace}.gate{margin-top:26px;display:grid;grid-template-columns:auto 1fr;gap:17px;padding:20px;border:1px solid rgba(255,112,72,.35);background:rgba(255,112,72,.065)}.gate .seal{width:68px;height:68px;display:grid;place-items:center;border:2px solid var(--orange);color:var(--orange);font:900 14px monospace;transform:rotate(-5deg)}.gate p{margin:0;color:#dfbdb4;line-height:1.58;font-size:14px}
.feedback{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;align-items:center;margin-top:45px}.feedback-card{min-height:170px;padding:20px;border:1px solid var(--line);background:#0a151f;position:relative}.feedback-card:not(:last-child)::after{content:"→";position:absolute;right:-18px;top:46%;z-index:3;color:var(--orange);font-size:22px}.feedback-card .k{color:var(--cyan);font:700 11px monospace;margin-bottom:22px}.feedback-card b{display:block;margin-bottom:9px}.feedback-card p{font-size:14.5px;line-height:1.55;color:var(--muted)}
.case{padding:25px 24px;border-top:3px solid var(--blue);background:rgba(255,255,255,.025);min-height:340px;display:flex;flex-direction:column}.case.jp{border-color:var(--orange)}.case.uk{border-color:var(--cyan)}.case.us{border-color:var(--gold)}.case .country{font:700 12px monospace;color:var(--muted);letter-spacing:.12em;margin-bottom:25px}.case p{color:#b8c6d0;font-size:16px;line-height:1.65}.case .verdict{margin-top:auto;padding-top:18px;border-top:1px solid var(--line);font-weight:700;font-size:15px}.case a{margin-top:12px;font-size:13px;color:var(--cyan);text-underline-offset:3px}
.usecase{padding:23px;border:1px solid var(--line);min-height:220px;background:rgba(255,255,255,.018)}.usecase .icon{font:700 13px monospace;color:var(--cyan);margin-bottom:36px}.usecase h3{font-size:21px}.usecase p{font-size:15px;color:var(--muted);line-height:1.55}.usecase .when{display:inline-block;margin-top:15px;color:var(--orange);font-size:13px;font-weight:800}
.roadmap{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:35px}.phase{padding:24px;border:1px solid var(--line);min-height:280px;background:rgba(255,255,255,.018)}.phase .year{font:700 12px monospace;color:var(--cyan)}.phase h3{margin:22px 0 17px}.phase ul{padding-left:18px;margin:0;color:var(--muted);font-size:15px;line-height:1.75}.phase:last-child{border-color:#684235}.phase:last-child .year{color:var(--orange)}
.decision-strip{display:grid;grid-template-columns:repeat(3,1fr);gap:1px;background:var(--line);margin-top:20px;border:1px solid var(--line)}.decision{background:#0a151f;padding:17px 20px;font-size:15px;line-height:1.55}.decision b{color:var(--orange)}
.closing{background:radial-gradient(circle at 50% 60%,rgba(81,210,215,.08),transparent 35%),#071018;text-align:center}.closing h2{font-size:clamp(42px,6.8vw,100px);margin-inline:auto}.closing .lead{margin-inline:auto}.link-row{display:flex;justify-content:center;flex-wrap:wrap;gap:10px;margin-top:35px}.primary-link{display:inline-flex;padding:13px 18px;border:1px solid var(--cyan);color:var(--cyan);text-decoration:none;font-size:15px;font-weight:800}.primary-link.orange-link{border-color:var(--orange);color:var(--orange)}
.source-grid{display:grid;grid-template-columns:1.1fr .9fr;gap:40px}.sources-list{display:grid;gap:10px}.source-item{padding:15px 17px;border:1px solid var(--line);font-size:13px;line-height:1.55}.source-item b{display:block;margin-bottom:4px}.source-item span{color:var(--muted)}.source-item a{color:var(--cyan);word-break:break-all;font-size:11px}.footnote{font-size:13.5px;line-height:1.65;color:var(--dim)}
@media(max-width:1050px){.grid4{grid-template-columns:repeat(2,1fr)}.flow-track{grid-template-columns:repeat(4,1fr);gap:22px}.flow-track::before{display:none}.algo{grid-template-columns:repeat(3,1fr)}.algo-step:nth-child(3)::after,.algo-step:last-child::after{display:none}.feedback{grid-template-columns:repeat(2,1fr)}.feedback-card::after{display:none!important}}
@media(max-width:780px){html{scroll-snap-type:y proximity}.topbar{height:52px}.scene{align-items:flex-start;padding-top:82px;overflow:visible}.scene-nav{display:none}.grid2,.grid3,.perf-grid,.source-grid,.roadmap{grid-template-columns:1fr}.grid4{grid-template-columns:1fr}.flow-track{grid-template-columns:repeat(2,1fr)}.algo{grid-template-columns:1fr}.algo-step::after{display:none}.feedback{grid-template-columns:1fr}.map-shell{min-height:480px}.map-shell canvas{height:480px}.map-note{max-width:72%}.controls{bottom:12px}.controls button{padding:8px 11px}.ledger-row{grid-template-columns:72px 1fr}.ledger-row span:last-child{grid-column:2}.decision-strip{grid-template-columns:1fr}.hero-foot{align-items:flex-start;flex-direction:column}.layer{grid-template-columns:45px 1fr}.layer .status{display:none}.bar-row{grid-template-columns:90px 1fr 58px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}.reveal{opacity:1;transform:none;transition:none}.scroll-cue i::after{animation:none}}

/* BEGIN AURORA SPECTRAL */
/* Aurora Spectral — 한글 조판·720p 밀도·지도 색 토큰 */
:root{
  --ink:#f8f5ff;--muted:#b5acd0;--dim:#786f98;--bg:#050515;--panel:#111026;
  --line:rgba(189,171,255,.2);--cyan:#5af7f2;--blue:#7867ff;--orange:#ff5c9c;--gold:#ffd36a;
  --map-bg:#09081d;--map-land:rgba(27,26,53,.46);--map-stroke:#82e8e4;
  --inset-bg:#0d0b25;--inset-line:#655a9a;--inset-text:#f8f5ff;
}
html,body{max-width:100%;overflow-x:hidden}
body{word-break:keep-all;overflow-wrap:normal;line-break:strict;hyphens:none;background:radial-gradient(ellipse at 15% 20%,rgba(90,247,242,.12),transparent 31%),radial-gradient(ellipse at 82% 18%,rgba(120,103,255,.2),transparent 34%),radial-gradient(ellipse at 60% 84%,rgba(255,92,156,.12),transparent 32%),var(--bg)}
body::before{opacity:.28;background-image:radial-gradient(rgba(255,255,255,.45) .55px,transparent .8px);background-size:7px 7px;mix-blend-mode:screen}
h1,h2,h3,.loop-note,.warning,.decision{text-wrap:balance;word-break:keep-all;overflow-wrap:normal;line-break:strict}
p,li,.lead,.axis-card p,.case p,.usecase p,.feedback-card p,.phase li,.flow-node span,.stat small,.layer small{text-wrap:pretty;word-break:keep-all;overflow-wrap:normal;line-break:strict}
.no-break{white-space:nowrap}
h1{font-size:clamp(44px,6.4vw,96px);line-height:.96;letter-spacing:-.055em}
h2{font-size:clamp(34px,4.35vw,62px);line-height:1.08;letter-spacing:-.038em}
h3{font-size:clamp(18px,1.55vw,24px);line-height:1.32}
.lead{font-size:clamp(17px,1.55vw,23px);line-height:1.58}
.grid2 h2{font-size:clamp(34px,3.65vw,52px);line-height:1.1}.grid2 .lead{font-size:clamp(17px,1.42vw,21px)}
.axis-pair{margin-top:20px}.axis-card{padding:18px}.axis-card .symbol{font-size:30px}.axis-card h3{margin-bottom:10px}.axis-card p{font-size:14px;line-height:1.58;margin-bottom:0}
.loop-note{font-size:14px;line-height:1.55}.warning{font-size:13.5px;line-height:1.55}.map-note{font-size:11.5px}.stat small{font-size:13px}.algo-step b{font-size:15px}.algo-step p{font-size:12.5px}.feedback-card p{font-size:13px}.case p{font-size:14px}
.adaptive-scene .grid2{grid-template-columns:minmax(430px,.9fr) minmax(540px,1.1fr)}
.adaptive-scene h2{font-size:clamp(32px,3.3vw,47px);margin-bottom:13px}.adaptive-scene .lead{margin-bottom:0}
.adaptive-scene .axis-pair{gap:12px;margin-top:16px}.adaptive-scene .axis-card{padding:16px}.adaptive-scene .axis-card p{font-size:13.2px;line-height:1.52}
.adaptive-scene .loop-note{margin-top:9px!important;padding:11px 15px}.adaptive-scene .warning{margin-top:8px!important;padding:11px 15px}
.topbar{background:rgba(7,6,24,.62);border-color:rgba(184,164,255,.18);backdrop-filter:blur(22px) saturate(145%)}
.progress{background:rgba(255,255,255,.08)}.progress span{background:linear-gradient(90deg,var(--cyan),var(--blue),var(--orange));box-shadow:0 0 18px var(--blue)}
.scene::before{content:"";position:absolute;width:52vw;height:52vw;right:-20vw;top:-25vw;border-radius:50%;background:conic-gradient(from 80deg,transparent,rgba(90,247,242,.08),rgba(120,103,255,.13),rgba(255,92,156,.08),transparent);filter:blur(32px);animation:auroraSpin 22s linear infinite;pointer-events:none}
@keyframes auroraSpin{to{transform:rotate(360deg)}}
.scene::after{color:rgba(185,166,255,.035);text-shadow:0 0 50px rgba(120,103,255,.2)}
.eyebrow{color:var(--cyan)}h2 .accent,.accent{color:var(--orange)}
.hero h1 .ghost{background:linear-gradient(90deg,var(--cyan),#9e8cff,var(--orange));-webkit-background-clip:text;color:transparent;-webkit-text-stroke:0}
.tag,.axis-card,.stat,.layer,.algo-step,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item,.map-legend,.map-note{background:linear-gradient(145deg,rgba(31,28,65,.7),rgba(11,10,34,.58));border-color:rgba(203,188,255,.2);backdrop-filter:blur(18px) saturate(135%);box-shadow:inset 0 1px rgba(255,255,255,.06),0 22px 55px rgba(0,0,0,.18)}
.axis-card,.stat,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item,.warning,.loop-note,.gate{border-radius:18px}
.map-shell{border-radius:26px;border-color:rgba(90,247,242,.25);box-shadow:0 35px 100px rgba(0,0,0,.48),0 0 80px rgba(120,103,255,.13)}.map-shell canvas{border-radius:26px}
.chip,.controls button{border-radius:999px;background:rgba(20,17,52,.72)}
.webgl-toggle{margin-left:auto;flex:0 0 auto;border:1px solid rgba(90,247,242,.34);border-radius:999px;background:rgba(12,12,38,.74);color:var(--cyan);padding:7px 11px;font:800 10px/1 Consolas,monospace;letter-spacing:.08em;cursor:pointer}.webgl-toggle[hidden]{display:none}.webgl-toggle[aria-pressed="false"]{color:var(--muted);border-color:rgba(255,255,255,.16)}
#heroWebGL{opacity:.94;filter:saturate(1.15)}#heroWebGL::after{content:"";position:absolute;inset:0;background:linear-gradient(90deg,rgba(5,5,21,.94) 0%,rgba(5,5,21,.52) 42%,rgba(5,5,21,.06) 72%);pointer-events:none}.webgl-paused .webgl-surface{visibility:hidden}.webgl-paused #fieldCanvas{opacity:.72}
.webgl-badge{position:absolute;right:18px;top:18px;z-index:4;padding:8px 11px;border:1px solid rgba(90,247,242,.25);border-radius:999px;background:rgba(7,6,24,.64);color:var(--cyan);font:700 9px/1 Consolas,monospace;letter-spacing:.1em;pointer-events:none}
.density-webgl{position:absolute;inset:0;z-index:2;opacity:0;visibility:hidden;transition:opacity .4s ease,visibility .4s;cursor:grab;touch-action:none;user-select:none}.density-webgl.is-dragging{cursor:grabbing}.density-webgl canvas{width:100%!important;height:100%!important}.map-shell.webgl-view .density-webgl{opacity:1;visibility:visible}.map-shell.webgl-view>#densityMap{opacity:0}.map-shell:not(.webgl-view) .density-webgl{pointer-events:none}.map-shell:not(.webgl-view) .density-spin{display:none}.map-view-label{display:none}.map-shell.webgl-view .map-view-label.webgl-copy,.map-shell:not(.webgl-view) .map-view-label.map-copy{display:inline}.webgl-unavailable [data-density-view="3d"],.webgl-unavailable .density-spin{display:none}
.layer-composite{display:grid;gap:9px}.layer-webgl{height:238px;position:relative;overflow:hidden;border:1px solid rgba(90,247,242,.2);border-radius:22px;background:radial-gradient(circle at 52% 38%,rgba(120,103,255,.16),rgba(7,6,24,.58) 64%);box-shadow:0 28px 80px rgba(0,0,0,.28)}.layer-webgl canvas{display:block;width:100%;height:100%}.layer-webgl .webgl-badge{top:12px;right:12px}.layer-focus{position:absolute;z-index:4;left:12px;top:11px;max-width:62%;padding:8px 10px;border-left:2px solid var(--cyan);background:rgba(7,6,24,.68);pointer-events:none}.layer-focus b,.layer-focus span{display:block}.layer-focus b{color:#f8f5ff;font-size:11px;margin-bottom:3px}.layer-focus span{color:#9fb7ca;font-size:9px;line-height:1.35}.layer-hint{position:absolute;z-index:4;right:12px;bottom:9px;color:#89a2b4;font:700 8.5px/1.2 Consolas,monospace;letter-spacing:.04em;pointer-events:none}.layer-composite .layer-stack{gap:6px}.layer-composite .layer{appearance:none;width:100%;color:inherit;text-align:left;font:inherit;cursor:pointer;padding:10px 14px;grid-template-columns:42px 1fr auto;border-radius:12px;transform:translateX(var(--shift));box-shadow:inset 0 1px rgba(255,255,255,.05);transition:opacity .22s ease,border-color .22s ease,background .22s ease}.layer-composite .layer[data-field-layer="core"]{--layer-accent:#5af7f2}.layer-composite .layer[data-field-layer="regional"]{--layer-accent:#9e8cff}.layer-composite .layer[data-field-layer="crustal"]{--layer-accent:#ff5c9c}.layer-composite .layer[data-field-layer="external"]{--layer-accent:#ffd36a}.layer-composite.has-layer-focus .layer[aria-pressed="false"]{opacity:.42}.layer-composite .layer[aria-pressed="true"]{opacity:1;border-color:var(--layer-accent);background:linear-gradient(100deg,color-mix(in srgb,var(--layer-accent) 11%,transparent),rgba(11,10,34,.72))}.layer-composite .layer:focus-visible{outline:2px solid var(--layer-accent);outline-offset:2px}.layer-composite .layer:nth-child(2){--shift:7px}.layer-composite .layer:nth-child(3){--shift:14px}.layer-composite .layer:nth-child(4){--shift:21px}.layer-composite .layer .n{font-size:16px}.layer-composite .layer b{font-size:14px;margin-bottom:2px}.layer-composite .layer small{font-size:11.5px}.layer-composite .layer .status{font-size:10px}
.layer-composite .layer .n,.layer-composite .layer .status{color:var(--layer-accent)}
.axis-card.point{box-shadow:inset 0 1px rgba(255,255,255,.08),0 0 48px rgba(90,247,242,.06)}.axis-card.region{box-shadow:inset 0 1px rgba(255,255,255,.08),0 0 48px rgba(255,92,156,.08)}
.formula{border:1px solid rgba(185,166,255,.24);border-left:3px solid var(--cyan);border-radius:18px;background:linear-gradient(100deg,rgba(90,247,242,.08),rgba(120,103,255,.08),rgba(255,92,156,.06))}
.bar-track{border-radius:8px;overflow:hidden}.bar-fill{background:linear-gradient(90deg,var(--blue),var(--cyan))}.bar-fill.orange{background:linear-gradient(90deg,var(--orange),var(--gold))}
.warning{background:rgba(255,92,156,.09);color:#f2bfd7}.loop-note{background:rgba(120,103,255,.08);border-color:rgba(185,166,255,.28)}
.closing{background:radial-gradient(circle at 50% 55%,rgba(90,247,242,.11),transparent 30%),radial-gradient(circle at 62% 45%,rgba(255,92,156,.1),transparent 28%),var(--bg)}
@media(min-width:781px) and (max-height:820px){
  .scene{padding-top:64px;padding-bottom:12px}.eyebrow{margin-bottom:10px}h2{margin-bottom:17px}.hero .kicker{margin-bottom:24px}.hero-foot{margin-top:30px}
  .map-shell{min-height:510px}.map-shell canvas{height:510px}.flow-track{margin-top:32px}.algo{margin-top:22px}.algo-step{min-height:156px}.feedback{margin-top:28px}.feedback-card{min-height:148px}.roadmap{margin-top:22px}
  .phase{min-height:238px}.case{min-height:300px}.usecase{min-height:190px}.formula{margin:18px 0;padding:18px 22px}.ledger{margin-top:18px}.perf-panel{padding:18px}.bar-row{margin:10px 0}.gate{margin-top:16px;padding:15px}
  .adaptive-scene .map-shell,.adaptive-scene .map-shell canvas{height:500px;min-height:500px}.layer-webgl{height:220px}.layer-composite .layer{padding-block:8px}.closing h2{font-size:clamp(44px,5.2vw,70px)}.source-grid{gap:28px}.sources-list{gap:6px}.source-item{padding:10px 13px;font-size:12px;line-height:1.42}.source-item a{font-size:10px}.footnote{font-size:12px;line-height:1.5}
}
@media(max-width:1050px){.adaptive-scene .grid2{grid-template-columns:1fr}.adaptive-scene .map-shell{order:2}}
@media(max-width:780px){h1{font-size:clamp(38px,10vw,58px)}h2,.grid2 h2{font-size:clamp(31px,8.4vw,46px)}.lead,.grid2 .lead{font-size:17px}.scene{padding-inline:22px}.adaptive-scene .map-shell{order:2}.axis-pair{grid-template-columns:1fr}.layer-webgl{height:220px}.map-toolbar{max-width:calc(100% - 28px)}.webgl-toggle{padding:6px 8px}.brand{max-width:42vw;overflow:hidden;text-overflow:ellipsis}.map-shell.webgl-view{touch-action:pan-y}}
@media(max-width:520px){.layer-webgl .webgl-badge{display:none}.layer-focus{max-width:78%}.layer-hint{font-size:8px}}
@media(max-width:520px){.flow-track,.stat-grid{grid-template-columns:1fr}.flow-node{text-align:left;display:grid;grid-template-columns:56px 1fr;column-gap:14px;align-items:center}.flow-node .dot{grid-row:1/3;margin:0}.flow-node b{align-self:end;margin:0 0 3px}.flow-node span{align-self:start}.stat small{font-size:13.5px}}
@media(prefers-reduced-motion:reduce){.scene::before{animation:none}.density-webgl,.map-shell canvas{transition:none}}
/* END AURORA SPECTRAL */
</style>
</head>
<body>
<header class="topbar"><div class="brand"><b>KOREA</b> LMM / AURORA SPECTRAL</div><div class="progress"><span id="progress"></span></div><button class="webgl-toggle" id="webglToggle" type="button" aria-pressed="true" hidden>WEBGL ON</button><div class="counter mono" id="counter">01 / 13</div></header>
<nav class="scene-nav" id="sceneNav" aria-label="장면 이동"></nav>
<div class="controls"><button id="prev" aria-label="이전 장면">↑ 이전</button><button id="next" aria-label="다음 장면">다음 ↓</button></div>

<main>
<section class="scene hero" data-no="01" data-title="표지"><div id="heroWebGL" class="webgl-surface" role="img" aria-label="한반도 관측망과 자기력선을 입체적으로 표현한 WebGL 장면"></div><canvas id="fieldCanvas" aria-hidden="true"></canvas><div class="webgl-badge" aria-hidden="true">THREE.JS · LIVE FIELD</div><div class="inner hero-copy">
  <div class="kicker reveal"><span class="tag">전문가·발주자 브리핑</span><span class="tag">단일 파일 · 오프라인</span><span class="tag" id="modelStamp">MODEL —</span></div>
  <h1 class="reveal" style="--delay:80ms"><span>30개의 점을</span><span class="ghost">하나의 국가 자기장 기준으로.</span></h1>
  <p class="lead reveal" style="--delay:160ms">선점의 목표는 점을 늘리는 일이 아닙니다. <strong><span class="no-break">전국 어디서나</span> 설명 가능한 </strong><span class="no-break"><strong>자기장</strong>을</span> <span class="no-break">만들고,</span> 설명하지 못하는 곳을 다음 관측으로 되돌려 보내는 일입니다.</p>
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
    <p class="lead reveal">표석은 30곳에 있지만, 관측 시기·품질·표석 이력이 모두 같지는 않습니다. 품질감사를 통과해 <strong>지금 모델에 실제로 들어간 점은 16곳</strong>입니다. 103개소는 그 다음 표석을 어디에 놓을지 정하기 위해 <strong>미리 현장을 확인한 곳</strong>입니다.</p>
    <div class="stat-grid reveal" style="--delay:120ms">
      <div class="stat cyanline"><span class="num" data-bind="network_sites">30</span><span class="unit">점</span><small>1등 지자기점 공간망<br>전국에 이미 설치된 표석</small></div>
      <div class="stat orange"><span class="num" data-bind="model_sites">16</span><span class="unit">점</span><small>현재 LMM 품질감사 통과 입력<br>2019·2022–2025</small></div>
      <div class="stat"><span class="num" data-bind="survey_sites">103</span><span class="unit">개소</span><small><b>선점 전 사전 현장조사</b><br>접근·자기교란·방위표지<br>등급 A/B/C/D<br>6 / 35 / 49 / 13</small></div>
      <div class="stat"><span class="num" data-bind="candidate_sites">178</span><span class="unit">후보</span><small>도상선점 후보<br>P1 62 · P2 57 · P3 59</small></div>
    </div>
  </div>
  <div class="map-shell reveal" style="--delay:100ms"><canvas id="networkMap" aria-label="기존망·모델 입력·현장조사 지점 지도"></canvas><div class="map-toolbar">
    <button class="chip" data-map="network" data-layer="existing" aria-pressed="true">기존망 30</button><button class="chip" data-map="network" data-layer="model" aria-pressed="true">LMM 16</button><button class="chip" data-map="network" data-layer="survey" aria-pressed="false">현장 103</button><button class="chip" data-map="network" data-layer="density" aria-pressed="false">상대 부족도</button>
  </div><div class="map-legend"><div class="legend-row"><i class="swatch" style="background:#eaf1f6"></i>기존망</div><div class="legend-row"><i class="swatch" style="background:#ff7048"></i>LMM 입력</div><div class="legend-row"><i class="swatch" style="background:#51d2d7"></i>현장조사</div></div><div class="map-note">버튼으로 레이어를 켜고 끌 수 있습니다.</div></div>
</div></section>

<section class="scene adaptive-scene" data-no="04" data-title="적응형 선점"><div class="inner grid2">
  <div class="map-shell reveal" id="densityShell"><canvas id="densityMap" aria-label="자기이상 공간복잡도 대비 상대 부족도 2차원 지도"></canvas><div id="densityWebGL" class="density-webgl" role="img" aria-label="상대 밀도 부족도를 높이로 표현한 한반도 3차원 WebGL 지도. 마우스나 터치로 회전할 수 있습니다."></div><div class="map-toolbar" aria-label="밀도 지도 보기 방식"><button class="chip" type="button" data-density-view="2d" aria-pressed="true">2D 지도</button><button class="chip" type="button" data-density-view="3d" aria-pressed="false">3D WebGL</button><button class="chip density-spin" id="densitySpin" type="button" aria-pressed="true">자동 회전 ON</button></div><div class="map-legend"><div class="legend-row"><i class="swatch" style="background:#ff7048"></i>보강 1순위</div><div class="legend-row"><i class="swatch" style="background:#f0bc54"></i>평균보다 부족</div><div class="legend-row"><i class="swatch" style="background:#427df4"></i>평균보다 양호</div><div class="legend-row"><i class="swatch" style="background:#244257"></i>충분</div></div><div class="map-note"><span class="map-view-label map-copy">색 = 상대 부족도</span><span class="map-view-label webgl-copy">색·기둥 높이 = 상대 부족도 · 드래그로 회전</span><br><span data-bind="land_cells">1,091</span> 국토셀 · 보강 1순위 <span data-bind="critical_cells">76</span> · σ 대체 <span data-bind="sigma_missing">57</span><br>울릉도·독도 위치 포함</div></div>
  <div>
    <div class="eyebrow reveal">TWO SCALES, ONE RULE</div><h2 class="reveal">권역 밀도는 공간복잡도로,<br>후보점은 국소 구배로 설계합니다.</h2>
    <p class="lead reveal">항공자력 자료를 <strong>점 대표성</strong>과 <strong>권역 밀도</strong>라는 서로 다른 질문에 사용합니다.</p>
    <div class="axis-pair reveal" style="--delay:120ms">
      <article class="axis-card point"><div class="symbol">|∇ΔT|</div><h3>① 후보점 대표성 · <b>저구배 우대</b></h3><p>KIGAM 1.5분(약 2.8 km) 격자를 중앙차분한 수평 공간구배(nT/km)를 후보 178점의 s5에 실제 반영했습니다. <b>1 km 해상도 분석은 아니며</b>, 최종 확정은 현장 정밀 자력측량으로 합니다.</p></article>
      <article class="axis-card region"><div class="symbol">σ<sub>25 km</sub></div><h3>② 권역 밀도 · <b>고복잡도 보강</b></h3><p>반경 25 km 자기이상에서 1차 지역추세를 제거한 표준편차 σ를 <b>자기이상 공간복잡도</b>로 사용합니다. σ가 클수록 권장 간격을 줄입니다.</p></article>
    </div>
    <div class="loop-note reveal" style="--delay:190ms"><b>고복잡도 권역을 보강하되, 그 안에서는 저구배 후보점을 고릅니다.</b></div><div class="warning reveal" style="--delay:230ms">Neyman 배분과 s5 곡선은 결측·소표본 한계가 있는 <strong>잠정 설계안</strong>입니다.</div>
  </div>
</div></section>

<section class="scene" data-no="05" data-title="LMM의 역할"><div class="inner grid2">
  <div>
    <div class="eyebrow reveal">WHY LMM</div><h2 class="reveal">관측점은 값 하나를 주고, LMM은 전국의 값을 설명합니다.</h2>
    <p class="lead reveal"><strong>LMM(Local Magnetic Model)</strong>은 전 지구 기준장에 한반도 관측·지각 이상을 결합해 위치와 시점별 편각 D, 복각 I, 총자력 F를 추정하는 지역 모델입니다.</p>
    <div class="formula reveal" style="--delay:120ms">B<span class="sub">obs, quiet</span> = B<span class="sub">obs</span> − δB<span class="sub">external</span><br> B̂<span class="sub">LMM</span>(r,t) = B<span class="sub">IGRF-14</span> + R<span class="sub">Korea</span> + C<span class="sub">KIGAM</span></div>
    <p class="equation-note reveal" style="--delay:180ms">현재판은 외부장을 미래 시점 예측에 더하지 않습니다. 네 관측소 자료로 <strong>관측 세션의 D·I 시간변화만 환산</strong>한 뒤, Core + Regional + Crustal의 baseline을 만듭니다. 현재 <span class="mono">subtract_DI</span>는 QC=QUIET만 고르는 엄격 필터가 아닙니다.</p>
  </div>
  <div class="layer-composite reveal" id="layerComposite" style="--delay:80ms">
    <div id="layersWebGL" class="layer-webgl webgl-surface" role="img" aria-label="전 지구 핵장, 한반도 지역보정, 국지 지각 자기이상, 관측시각 외부장 전처리를 서로 다른 물리 형태로 구분한 WebGL 모식도"><span class="webgl-badge" aria-hidden="true">FIELD SOURCES</span><div class="layer-focus" id="layerFocus" aria-live="polite"><b>결합 보기</b><span>공간 3층 + 관측시각 전처리 1층</span></div><div class="layer-hint" aria-hidden="true">아래 층 카드를 눌러 분리 보기</div></div>
    <div class="layer-stack">
      <button class="layer" type="button" data-field-layer="core" aria-pressed="false"><span class="n">01</span><span><b>Core · 전 지구 핵 기원장</b><small>IGRF-14, 구면조화 degree 13</small></span><span class="status">기준장</span></button>
      <button class="layer" type="button" data-field-layer="regional" aria-pressed="false"><span class="n">02</span><span><b>Regional · 한반도 장파장 잔차</b><small>16측점 · LOO 선택 결과 현재 0차(상수)</small></span><span class="status">D·I·F</span></button>
      <button class="layer" type="button" data-field-layer="crustal" aria-pressed="false"><span class="n">03</span><span><b>Crustal · 지각 단파장 이상</b><small>KIGAM 1.5분 격자 · F 전용</small></span><span class="status">F only</span></button>
      <button class="layer" type="button" data-field-layer="external" aria-pressed="false"><span class="n">04</span><span><b>External · 시간 변화 환산</b><small>4개 관측소 1분 자료 · D·I 전처리</small></span><span class="status">전처리</span></button>
    </div>
  </div>
</div></section>

<section class="scene" data-no="06" data-title="데이터와 알고리듬"><div class="inner">
  <div class="eyebrow reveal">DATA → ALGORITHM → MODEL</div><h2 class="reveal">수식보다 중요한 것은, 어떤 오차를 어느 단계에서 제거했는가입니다.</h2>
  <div class="algo reveal" style="--delay:100ms">
    <div class="algo-step"><div class="k">01 / AUDIT</div><b>성과 감사</b><p>좌표·이설·야장·관측시각·재방문 산포로 입력 신뢰도를 분리합니다.</p></div>
    <div class="algo-step"><div class="k">02 / REDUCE</div><b>외부장 환산</b><p>네 관측소의 시간 편차를 공간보간해 세션 D·I에서 뺍니다.</p></div>
    <div class="algo-step"><div class="k">03 / RESIDUAL</div><b>IGRF 잔차</b><p>실제 관측이 이뤄진 <b>측점·관측시기마다</b> ΔD, ΔI, ΔF = 관측값 − IGRF-14를 계산합니다.</p></div>
    <div class="algo-step"><div class="k">04 / ROBUST</div><b>강건 선별</b><p>MAD 2.5σ 반복 선별로 블런더·국지 자성체 영향을 제한합니다.</p></div>
    <div class="algo-step"><div class="k">05 / FIT</div><b>층별 적합</b><p>F에는 KIGAM 이상을 1:1 결합하고 남은 장파장 잔차를 적합합니다.</p></div>
    <div class="algo-step"><div class="k">06 / LOO</div><b>독립 예측</b><p>한 점씩 빼고 예측해 차수와 성능을 평가합니다. 현재 최적은 0차입니다.</p></div>
  </div>
  <div class="ledger reveal" id="ledger" style="--delay:190ms"></div>
</div></section>

<section class="scene" data-no="07" data-title="현재 성능"><div class="inner">
  <div class="eyebrow reveal">QUALITY GATE</div><h2 class="reveal">지각층이 F를 크게 줄였습니다. 다만 <b>정확도 목표를 정하기엔 자료가 더 필요합니다.</b></h2>
  <div class="perf-grid reveal" style="--delay:100ms"><div class="perf-panel"><h3>편각 D <small>RMS · degree</small></h3><div id="perfD"></div></div><div class="perf-panel"><h3>총자력 F <small>RMS · nT</small></h3><div id="perfF"></div></div></div>
  <div class="gate reveal" style="--delay:180ms"><div class="seal">HOLD</div><p>현재 <strong>LOO D <span data-bind="loo_D">0.499</span>°</strong> · <strong>LOO F <span data-bind="loo_F">62.8</span> nT</strong> 입니다. <b>이 모델의 정확도 목표는 아직 확정되지 않았습니다.</b> 자주 인용되는 D&lt;0.1°·F&lt;50 nT는 설명자료의 <b>공학적 목표치</b>일 뿐 발주기관이 정한 값이 아니며, 참고로 두면 D는 약 <span data-bind="loo_D_ratio">5.0</span>배 · F는 약 <span data-bind="loo_F_ratio">1.26</span>배 남았습니다. <b>목표를 정하려면 먼저 «어느 수준까지 도달할 수 있는가»를 알아야 하고, 그러려면 모델 개선과 관측자료 확충이 더 필요합니다.</b> 그 전까지 현재판은 진단·설계 참고용이며 국가기본도 자편각 등 정식 성과에 쓰지 않습니다. 참고로 법정 관측 정수차 30′와는 성격이 다른 기준입니다.</p></div>
  <p class="footnote reveal" style="--delay:240ms;margin-top:15px">표본내 단계 RMS: D IGRF → Regional, F IGRF → Crustal → Regional. LOO는 적합에 쓰지 않은 점을 예측한 값이므로 최종 판정은 LOO를 우선합니다. <span class="mono">lmm_model.json</span>의 External 설명은 2026-08-27에 실제 <span class="mono">subtract_DI</span> 구현(D·I 전처리 보정 · F 미보정 · strict quiet 필터 없음 · 예측 시 미평가)에 맞춰 정정했습니다.</p>
</div></section>

<section class="scene" data-no="08" data-title="LMM과 선점의 되먹임"><div class="inner">
  <div class="eyebrow reveal">THE FEEDBACK LOOP</div><h2 class="reveal">첫 후보는 모델 없이 고르고, <b>다음 후보는 모델이 알려줍니다.</b></h2>
  <p class="lead reveal">처음에는 모델을 쓰지 않습니다 — 전국 공백·인공 간섭원·국소 자기이상 구배가 작은 자리만 보고 고릅니다. 관측을 마친 뒤에야 LMM이 <strong>“아직 설명하지 못하는 지역”</strong>을 알려주고, 거기가 다음 관측 대상이 됩니다. 이 순서를 지켜야 <strong>모델이 고른 점으로 그 모델을 검증하는 순환논리</strong>를 피할 수 있습니다.</p>
  <div class="feedback reveal" style="--delay:120ms">
    <div class="feedback-card"><div class="k">1단계</div><b>전국 공백 메우기</b><p>모델을 쓰지 않습니다. 도엽·담당면적·접근성만 보고 비어 있는 곳부터 채웁니다.</p></div>
    <div class="feedback-card"><div class="k">2단계</div><b>어느 지역을 촘촘히</b><p>자기이상 공간복잡도가 큰 권역을 항공자력으로 골라 더 촘촘한 배치를 검토합니다.</p></div>
    <div class="feedback-card"><div class="k">3단계</div><b>그 안에서 어느 자리</b><p>고른 권역 안에서 국소 자기이상 구배가 작고 인공 교란이 없는 자리를 찾습니다.</p></div>
    <div class="feedback-card"><div class="k">4단계</div><b>관측하고 다시 맞추기</b><p>새로 잰 D·I·F를 모델에 넣고, 독립 검증으로 실제로 좋아졌는지 확인합니다.</p></div>
    <div class="feedback-card"><div class="k">5단계</div><b>남은 곳만 다시</b><p>그래도 모델이 설명하지 못한 지역만 골라 1단계로 돌아갑니다.</p></div>
  </div>
  <div class="loop-note reveal" style="--delay:220ms">새 점이 좋은 점이었는지는 <b>점수표로 정하지 않습니다.</b> 그 점을 <b>넣기 전과 넣은 뒤의 예측오차 차이</b>로 확인합니다.</div>
</div></section>

<section class="scene" data-no="09" data-title="해외 근거"><div class="inner">
  <div class="eyebrow reveal">FOREIGN EVIDENCE, PRECISELY</div><h2 class="reveal">그대로 가져올 해외 규칙은 없습니다. 다만 <b>설계의 전제는 이미 확인됐습니다.</b></h2>
  <div class="grid3 reveal" style="--delay:100ms">
    <article class="case jp"><div class="country">JAPAN / GSI</div><h3>구배 대비 밀도가 낮으면 모델 오차가 커집니다.</h3><p>평균 약 20 km 지상망은 더 짧은 파장의 이상을 표현하지 못해 항공자력으로 보완했습니다. 야쓰가타케·홋카이도 동부의 큰 오차는 큰 자기구배에 비해 1·2등 자기점 밀도가 낮은 데서 왔다고 분석했습니다.</p><div class="verdict">우리 설계의 가장 강한 실증 근거<br><span class="muted">단, 자동 밀도배분 선례는 아님</span></div><a href="https://www.gsi.go.jp/common/000071404.pdf" target="_blank" rel="noreferrer">GSI 2010.0 모델 정확도 평가 ↗</a><a href="https://www.gsi.go.jp/common/000024739.pdf" target="_blank" rel="noreferrer">GSI 항공자력·20 km 지상망 설명 ↗</a></article>
    <article class="case uk"><div class="country">UNITED KINGDOM / BGS</div><h3>시간 변화와 지각 이상은 역할을 나눕니다.</h3><p>41개 반복점과 3개 관측소로 UK Regional Geomagnetic Model을 갱신하며, 반복점은 자기 간섭이 적은 곳을 택합니다. 반복점은 주로 core field의 시간변화, 항공자력은 crustal field를 담당합니다.</p><div class="verdict">현재 LMM의 Regional + Crustal 분리와 정합</div><a href="https://www.geomag.bgs.ac.uk/operations/uksurvey.html" target="_blank" rel="noreferrer">BGS UK Magnetic Survey ↗</a><a href="https://geomag.bgs.ac.uk/education/earthmag.html" target="_blank" rel="noreferrer">BGS 자료 역할 설명 ↗</a></article>
    <article class="case us"><div class="country">US / UK JOINT MODEL</div><h3>관측 공백은 다른 자료로 보완하되, 역할을 혼동하지 않습니다.</h3><p>1995 개정 모델은 관측소·반복점을 주 자료로 쓰고, 이 자료가 없는 영역에서 항공자력·위성 자료로 영년변화 정보를 보완했습니다.</p><div class="verdict">공백 보완의 개념적 선례<br><span class="muted">신규 반복점 밀도 알고리듬은 아님</span></div><a href="https://pubs.usgs.gov/publication/70019514" target="_blank" rel="noreferrer">USGS Publications Warehouse ↗</a></article>
  </div>
</div></section>

<section class="scene" data-no="10" data-title="최종 활용"><div class="inner">
  <div class="eyebrow reveal">FROM MODEL TO PUBLIC VALUE</div><h2 class="reveal">최종 성과는 <b>자기도(지도)</b>와 <b>위치·시점별 계산 서비스</b>, 두 형태로 나갑니다.</h2>
  <p class="lead reveal">일본 국토지리원은 같은 성격의 결과를 <strong>전국 자기도(磁気図) 웹 지도</strong>로 공개하고 있습니다. 지도는 최종 산출물의 한 형태이고, 같은 모델이 도폭 자편각·항법 보정·탐사 기준장으로도 이어집니다.</p>
  <div class="grid4 reveal" style="--delay:100ms">
    <article class="usecase"><div class="icon">01 / MAP</div><h3>자기도·국가기본도 자편각</h3><p>전국 편각·복각·총자력 분포도를 내고, 도폭별 자편각과 연변화율을 같은 모델에서 산출합니다. 일본 GSI 자기도와 같은 형태의 공개 서비스가 가능합니다.</p><span class="when">정확도 목표 확정·검증 후</span></article>
    <article class="usecase"><div class="icon">02 / FIELD</div><h3>측량·관측 품질관리</h3><p>관측값과 예상값의 차이로 방위표지 오류, 이설, 국지 교란, 시각 보정 누락을 조기에 탐지합니다.</p><span class="when">현재도 진단 활용 가능</span></article>
    <article class="usecase"><div class="icon">03 / NAV</div><h3>항법·헤딩 보정</h3><p>차량·드론·해양·자기센서의 지역 편차를 보정하는 고해상도 기준으로 확장할 수 있습니다.</p><span class="when">독립 검증·배포 체계 필요</span></article>
    <article class="usecase"><div class="icon">04 / GEO</div><h3>지질·인프라 기준장</h3><p>항공·지상 자력탐사, 방향성 시추, 지하 구조 해석에서 core·crustal·external 성분을 분리하는 기준이 됩니다.</p><span class="when">목적별 해상도 추가 필요</span></article>
  </div>
  <p class="footnote reveal" style="--delay:190ms;margin-top:22px">일본 GSI는 자기도를 웹 지도로 공개하고, 영국 BGS는 지역모델을 지도 magnetic north 갱신에 사용합니다. WMM/BGGM 계열은 항법·헤딩·방향성 시추 등에 쓰입니다. 국내 적용은 정확도 목표·책임주체·갱신주기를 별도로 확정해야 합니다.</p>
</div></section>

<section class="scene" data-no="11" data-title="의사결정 로드맵"><div class="inner">
  <div class="eyebrow reveal">RESEARCH ROADMAP</div><h2 class="reveal">추가 연구는 <b>점을 늘리는 일</b>이 아니라, <b>품질 게이트를 통과시키는 일</b>입니다.</h2>
  <div class="roadmap reveal" style="--delay:100ms">
    <article class="phase"><div class="year">GATE 01 / INPUT</div><h3>자료 계보를 닫습니다.</h3><ul><li>30점 좌표 대조표의 원본파일명·작성일·문서번호 확인</li><li>표석 이설·개명은 지점코드로 추적</li><li>모든 관측에 시각·장비·방위표지 메타데이터 결합</li></ul></article>
    <article class="phase"><div class="year">GATE 02 / CAMPAIGN</div><h3>가치가 큰 점부터 관측합니다.</h3><ul><li>기본 공간 공백 + 항공자력 복잡권역을 단계 배치</li><li>각 점은 저구배·저인공교란 현장 검증</li><li>신규점 수는 순차 LOO 개선량으로 재결정</li></ul></article>
    <article class="phase"><div class="year">GATE 03 / RELEASE</div><h3>독립 검증 뒤 공개합니다.</h3><ul><li>외부장 D·I·F 전 성분 환산 체계 확정</li><li>항공자력 결측·벡터 지각장 개선</li><li>달성 가능한 정확도 목표를 먼저 확정</li><li>독립 검증 통과 후 공식 서비스</li></ul></article>
  </div>
  <div class="decision-strip reveal" style="--delay:190ms"><div class="decision"><b>연구과제 1</b><br>신규점 수와 적응형 보강 원칙의 <span class="no-break">확정 순서</span></div><div class="decision"><b>연구과제 2</b><br>관측시각·상시관측소 연계의 필수 성과 <span class="no-break">지정 여부</span></div><div class="decision"><b>연구과제 3</b><br>고복잡도 권역 정밀 항공·지상 자력자료의 추가 <span class="no-break">확보 여부</span></div></div>
</div></section>

<section class="scene closing" data-no="12" data-title="결론"><div class="inner">
  <div class="eyebrow reveal" style="justify-content:center">THE ONE SENTENCE</div><h2 class="reveal">좋은 선점은 저구배 후보점을 찾고,<br><span class="accent">좋은 모델은 설명 못한 곳을 다시 찾습니다.</span></h2>
  <p class="lead reveal">관측망·항공자력·상시관측·IGRF를 하나로 묶어 검증하면, 새로 놓는 점 하나하나가 <strong>국가 자기장 기준을 얼마나 정밀하게 만드는지</strong> 숫자로 확인할 수 있습니다. 점을 늘리는 것 자체가 목표가 아니라, <strong>설명되지 않는 지역을 없애는 것</strong>이 목표입니다.</p>
  <div class="link-row reveal" style="--delay:140ms"><a class="primary-link" href="index.html">입지 선정 지도 열기</a><a class="primary-link orange-link" href="lmm.html">LMM 계산기 열기</a></div>
  <p class="footnote reveal" style="--delay:210ms;margin:28px auto 0;max-width:760px">현재 모델은 시범 구축판이고 <b>정확도 목표는 아직 확정되지 않았습니다.</b> 이 브리핑은 연구설계와 의사결정을 위한 자료입니다.</p>
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
__THREE_JS__
</script>
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
  perfRows('#perfD',[{label:'IGRF 단독',value:m.stage.D_igrf,text:fmt(m.stage.D_igrf,3)+'°'},{label:'+ Regional',value:m.stage.D_final,text:fmt(m.stage.D_final,3)+'°',cls:'cyanbar'},{label:'LOO 예측',value:m.loo.D,text:fmt(m.loo.D,3)+'°',cls:'orange'},{label:'참고 목표치',value:.1,text:'0.100° (미확정)'}],.6);
  perfRows('#perfF',[{label:'IGRF 단독',value:m.stage.F_igrf,text:fmt(m.stage.F_igrf,1)},{label:'+ Crustal',value:m.stage.F_crust,text:fmt(m.stage.F_crust,1),cls:'cyanbar'},{label:'+ Regional',value:m.stage.F_final,text:fmt(m.stage.F_final,1),cls:'cyanbar'},{label:'LOO 예측',value:m.loo.F,text:fmt(m.loo.F,1),cls:'orange'},{label:'참고 목표치',value:50,text:'50.0 (미확정)'}],100);
}

function drawGeometry(ctx,geom,project){ctx.beginPath();
  const ring=r=>{r.forEach((p,i)=>{const [x,y]=project(p[0],p[1]);if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)});ctx.closePath()};
  if(geom.type==='Polygon')geom.coordinates.forEach(ring);else if(geom.type==='MultiPolygon')geom.coordinates.forEach(p=>p.forEach(ring));
}
const css=(name,fallback)=>getComputedStyle(document.documentElement).getPropertyValue(name).trim()||fallback;
function drawEastSeaInset(ctx,W,H){
  const iw=Math.min(184,Math.max(142,W*.34)),ih=112,x=W-iw-14,y=14;
  const b={w:130.62,e:132.08,s:37.05,n:37.68},pad=10;
  const scale=Math.min((iw-pad*2)/(b.e-b.w),(ih-pad*2)/(b.n-b.s));
  const ox=x+(iw-(b.e-b.w)*scale)/2,oy=y+(ih-(b.n-b.s)*scale)/2;
  const project=(lon,lat)=>[ox+(lon-b.w)*scale,oy+(b.n-lat)*scale];
  ctx.save();ctx.fillStyle=css('--inset-bg','#071018');ctx.fillRect(x,y,iw,ih);
  ctx.strokeStyle=css('--inset-line','#355267');ctx.lineWidth=1;ctx.strokeRect(x+.5,y+.5,iw-1,ih-1);
  ctx.beginPath();ctx.rect(x+1,y+1,iw-2,ih-2);ctx.clip();
  drawGeometry(ctx,DATA.boundary,project);ctx.fillStyle=css('--map-land','#0d1b26');ctx.fill('evenodd');
  ctx.strokeStyle=css('--map-stroke','#66869a');ctx.lineWidth=.8;ctx.stroke();
  DATA.east_sea.forEach(site=>{const [px,py]=project(site.lon,site.lat);ctx.beginPath();ctx.arc(px,py,site.name==='독도'?3.2:2.6,0,Math.PI*2);
    ctx.fillStyle=site.name==='독도'?css('--orange','#ff7048'):css('--cyan','#51d2d7');ctx.fill();
    if(site.name==='독도'){ctx.beginPath();ctx.arc(px,py,7,0,Math.PI*2);ctx.strokeStyle=css('--orange','#ff7048');ctx.lineWidth=.8;ctx.stroke()}
    ctx.font='700 10px "Noto Sans KR","Malgun Gothic",sans-serif';ctx.fillStyle=css('--inset-text','#dce9f1');
    ctx.textAlign=site.name==='독도'?'right':'left';ctx.textBaseline='middle';ctx.fillText(site.name,px+(site.name==='독도'?-7:7),py-1)});
  ctx.restore();ctx.font='700 8px Consolas,monospace';ctx.fillStyle=css('--inset-text','#dce9f1');ctx.textAlign='left';ctx.textBaseline='top';ctx.fillText('EAST SEA ISLANDS',x+8,y+7);
}
const colors={critical:'#ff7048',short:'#f0bc54',mild:'#427df4',ok:'#244257'};
function makeMap(canvas,mode){const ctx=canvas.getContext('2d');const state={existing:true,model:true,survey:false,density:mode==='density'};
  function render(){const rect=canvas.getBoundingClientRect(),dpr=Math.min(2,devicePixelRatio||1);canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);const W=rect.width,H=rect.height;ctx.clearRect(0,0,W,H);ctx.fillStyle=css('--map-bg','#09141d');ctx.fillRect(0,0,W,H);
    const b={w:124.25,e:130.15,s:32.9,n:38.9},pad=32,scale=Math.min((W-pad*2)/(b.e-b.w),(H-pad*2)/(b.n-b.s)),ox=(W-(b.e-b.w)*scale)/2,oy=(H-(b.n-b.s)*scale)/2;const project=(lon,lat)=>[ox+(lon-b.w)*scale,oy+(b.n-lat)*scale];
    if(state.density){DATA.density.forEach(c=>{const [x,y]=project(c.lon,c.lat);ctx.globalAlpha=c.observed?.62:.22;ctx.fillStyle=colors[c.tier];const s=Math.max(3.2,scale*.095);ctx.fillRect(x-s/2,y-s/2,s,s);if(!c.observed){ctx.strokeStyle='#b8c3cb';ctx.lineWidth=.7;ctx.beginPath();ctx.moveTo(x-s/2,y-s/2);ctx.lineTo(x+s/2,y+s/2);ctx.stroke()}});ctx.globalAlpha=1}
    drawGeometry(ctx,DATA.boundary,project);ctx.fillStyle=css('--map-land','rgba(13,27,38,.46)');ctx.fill('evenodd');ctx.strokeStyle=css('--map-stroke','#5d798b');ctx.lineWidth=1.1;ctx.stroke();
    if(state.survey)drawPoints(DATA.survey,'#51d2d7',2.3,.52);
    if(state.existing)drawPoints(DATA.existing,'#eaf1f6',3.2,.9,true);
    if(state.model)drawPoints(DATA.model_sites,'#ff7048',4.4,1);
    drawEastSeaInset(ctx,W,H);
    function drawPoints(rows,color,r,alpha,ring=false){ctx.globalAlpha=alpha;rows.forEach(p=>{const [x,y]=project(p.lon,p.lat);ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ring?(ctx.strokeStyle=color,ctx.lineWidth=1.4,ctx.stroke()):(ctx.fillStyle=color,ctx.fill())});ctx.globalAlpha=1}
  }
  new ResizeObserver(render).observe(canvas);render();return{state,render};
}
const maps={};
function initMaps(){maps.network=makeMap(document.querySelector('#networkMap'),'network');maps.density=makeMap(document.querySelector('#densityMap'),'density');document.querySelectorAll('[data-map]').forEach(btn=>btn.addEventListener('click',()=>{const m=maps[btn.dataset.map],k=btn.dataset.layer;m.state[k]=!m.state[k];btn.setAttribute('aria-pressed',String(m.state[k]));m.render()}));}

const WEBGL={available:false,enabled:true,reduced:matchMedia('(prefers-reduced-motion:reduce)').matches,engines:[],raf:0,densityView:'2d'};
function canUseWebGL(){try{const c=document.createElement('canvas');return !!(window.WebGLRenderingContext&&(c.getContext('webgl2')||c.getContext('webgl')||c.getContext('experimental-webgl')))}catch(_){return false}}
function geometryRings(geom){if(!geom)return[];if(geom.type==='Polygon')return geom.coordinates;if(geom.type==='MultiPolygon')return geom.coordinates.flat();return[]}
function geoVector(lon,lat,y=0,scale=.9,centerLon=127.45,centerLat=35.65){return new THREE.Vector3((lon-centerLon)*scale,y,-(lat-centerLat)*scale*1.06)}
function addBoundary3D(parent,{y=.02,scale=.9,color=0x5af7f2,opacity=.72,centerLon=127.45,centerLat=35.65}={}){const group=new THREE.Group();geometryRings(DATA.boundary).forEach(ring=>{const pts=ring.map(p=>geoVector(p[0],p[1],y,scale,centerLon,centerLat));if(pts.length<2)return;const g=new THREE.BufferGeometry().setFromPoints(pts);const line=new THREE.Line(g,new THREE.LineBasicMaterial({color,transparent:true,opacity}));group.add(line)});parent.add(group);return group}
function pointCloud(rows,{color=0xff5c9c,size=.08,y=.12,scale=.9,centerLon=127.45,centerLat=35.65,opacity=1}={}){const pts=rows.map(p=>geoVector(p.lon,p.lat,y,scale,centerLon,centerLat));const geometry=new THREE.BufferGeometry().setFromPoints(pts);return new THREE.Points(geometry,new THREE.PointsMaterial({color,size,transparent:opacity<1,opacity,sizeAttenuation:true}))}
function textSprite(text,color='#f8f5ff'){const c=document.createElement('canvas');c.width=320;c.height=80;const ctx=c.getContext('2d');ctx.clearRect(0,0,c.width,c.height);ctx.font='700 30px "Noto Sans KR","Malgun Gothic",sans-serif';ctx.textAlign='center';ctx.textBaseline='middle';ctx.shadowColor='rgba(5,5,21,.95)';ctx.shadowBlur=10;ctx.lineWidth=7;ctx.strokeStyle='rgba(5,5,21,.9)';ctx.strokeText(text,160,40);ctx.fillStyle=color;ctx.fillText(text,160,40);const texture=new THREE.CanvasTexture(c);texture.minFilter=THREE.LinearFilter;const sprite=new THREE.Sprite(new THREE.SpriteMaterial({map:texture,transparent:true,depthTest:false}));sprite.scale.set(1.25,.31,1);return sprite}
function makeRenderer(host,scene,camera,sceneIndex,updateScene){const renderer=new THREE.WebGLRenderer({antialias:true,alpha:true,powerPreference:'high-performance'});renderer.setPixelRatio(Math.min(1.75,devicePixelRatio||1));renderer.setClearColor(0x000000,0);renderer.outputEncoding=THREE.sRGBEncoding;renderer.domElement.setAttribute('aria-hidden','true');host.appendChild(renderer.domElement);const draw=(now,force)=>{updateScene(now,force);renderer.render(scene,camera)},engine={host,scene,camera,renderer,sceneIndex,renderScene:draw,resize(){const w=Math.max(1,host.clientWidth),h=Math.max(1,host.clientHeight);renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();draw(performance.now(),true)}};new ResizeObserver(()=>engine.resize()).observe(host);renderer.domElement.addEventListener('webglcontextlost',e=>{e.preventDefault();document.body.classList.add('webgl-paused')});engine.resize();WEBGL.engines.push(engine);return engine}
function buildHeroWebGL(){const host=document.querySelector('#heroWebGL'),scene=new THREE.Scene();scene.fog=new THREE.FogExp2(0x050515,.047);const camera=new THREE.PerspectiveCamera(38,1,.1,80);camera.position.set(7.2,5.4,9.3);const root=new THREE.Group();root.position.x=2.15;root.rotation.y=-.1;scene.add(root);addBoundary3D(root,{y:.04,scale:.9,color:0x5af7f2,opacity:.82});root.add(pointCloud(DATA.existing,{color:0xf8f5ff,size:.065,y:.11,opacity:.72}));root.add(pointCloud(DATA.model_sites,{color:0xff5c9c,size:.12,y:.16}));
  const field=new THREE.Group();for(let i=-8;i<=8;i++){const curvePts=[];for(let k=0;k<=42;k++){const x=-5.8+k*(11.6/42),z=i*.29+Math.sin(x*.72+i*.31)*.12,y=.48+Math.cos(x*.38+i*.17)*.18+Math.abs(i)*.018;curvePts.push(new THREE.Vector3(x,y,z))}const curve=new THREE.CatmullRomCurve3(curvePts);const g=new THREE.BufferGeometry().setFromPoints(curve.getPoints(120));const hot=i%4===0;field.add(new THREE.Line(g,new THREE.LineBasicMaterial({color:hot?0xff5c9c:0x5af7f2,transparent:true,opacity:hot?.3:.14})))}root.add(field);
  const particleCount=720,pos=new Float32Array(particleCount*3),col=new Float32Array(particleCount*3),cyan=new THREE.Color(0x5af7f2),violet=new THREE.Color(0x7867ff),mix=new THREE.Color();for(let i=0;i<particleCount;i++){pos[i*3]=(Math.random()-.5)*13;pos[i*3+1]=Math.random()*5-1;pos[i*3+2]=(Math.random()-.5)*9;mix.copy(cyan).lerp(violet,Math.random());col.set([mix.r,mix.g,mix.b],i*3)}const pg=new THREE.BufferGeometry();pg.setAttribute('position',new THREE.BufferAttribute(pos,3));pg.setAttribute('color',new THREE.BufferAttribute(col,3));const particles=new THREE.Points(pg,new THREE.PointsMaterial({size:.025,vertexColors:true,transparent:true,opacity:.66,depthWrite:false,blending:THREE.AdditiveBlending}));scene.add(particles);
  DATA.east_sea.forEach(site=>{const p=geoVector(site.lon,site.lat,.17,.9);const marker=new THREE.Mesh(new THREE.SphereGeometry(site.name==='독도'?.075:.06,14,10),new THREE.MeshBasicMaterial({color:site.name==='독도'?0xff5c9c:0x5af7f2}));marker.position.copy(p);root.add(marker)});let px=0,py=0;addEventListener('pointermove',e=>{px=e.clientX/innerWidth-.5;py=e.clientY/innerHeight-.5},{passive:true});return makeRenderer(host,scene,camera,0,(now)=>{const t=now*.0001;root.rotation.y=-.1+px*.07+Math.sin(t)*.018;root.rotation.x=-.12+py*.025;field.position.y=Math.sin(t*1.4)*.025;particles.rotation.y=t*.08;camera.lookAt(1.25,.15,0);})}
function buildDensityWebGL(){const host=document.querySelector('#densityWebGL'),scene=new THREE.Scene();scene.fog=new THREE.FogExp2(0x050515,.032);const camera=new THREE.PerspectiveCamera(39,1,.1,80);const root=new THREE.Group();root.rotation.y=-.08;scene.add(root);const geo=new THREE.BoxGeometry(.071,1,.071),matrix=new THREE.Matrix4(),position=new THREE.Vector3(),quaternion=new THREE.Quaternion(),scaleV=new THREE.Vector3();['critical','short','mild','ok'].forEach(tier=>{[true,false].forEach(observed=>{const cells=DATA.density.filter(cell=>cell.tier===tier&&cell.observed===observed);if(!cells.length)return;const mat=new THREE.MeshBasicMaterial({color:colors[tier]||'#7867ff',transparent:true,opacity:observed?.92:.3}),mesh=new THREE.InstancedMesh(geo,mat,cells.length);cells.forEach((cell,i)=>{const p=geoVector(cell.lon,cell.lat,0,.82,127.55,35.65),r=Number.isFinite(cell.relative)?cell.relative:1,h=.07+Math.max(0,Math.min(1.28,(r-.35)*.82));position.set(p.x,h/2,p.z);scaleV.set(1,h,1);matrix.compose(position,quaternion,scaleV);mesh.setMatrixAt(i,matrix)});mesh.instanceMatrix.needsUpdate=true;root.add(mesh)})});addBoundary3D(root,{y:.045,scale:.82,color:0x82e8e4,opacity:1,centerLon:127.55,centerLat:35.65});root.add(pointCloud(DATA.model_sites,{color:0xff8db8,size:.105,y:1.02,scale:.82,centerLon:127.55,centerLat:35.65}));const grid=new THREE.GridHelper(7.7,22,0x594b8d,0x211d46);grid.position.y=-.015;grid.material.transparent=true;grid.material.opacity=.28;root.add(grid);
  DATA.east_sea.forEach(site=>{const p=geoVector(site.lon,site.lat,.1,.82,127.55,35.65),marker=new THREE.Mesh(new THREE.SphereGeometry(site.name==='독도'?.08:.065,14,10),new THREE.MeshBasicMaterial({color:site.name==='독도'?0xff5c9c:0x5af7f2}));marker.position.copy(p);root.add(marker);const label=textSprite(site.name,site.name==='독도'?'#ff8db8':'#8ffffb');label.position.set(p.x,p.y+.34,p.z);root.add(label)});
  const target=new THREE.Vector3(.15,.18,0),radius=13.45,spinButton=document.querySelector('#densitySpin');let yaw=.62,pitch=.59,spinning=!WEBGL.reduced,dragging=false,lastX=0,lastY=0,lastNow=performance.now();
  const syncSpinButton=()=>{if(!spinButton)return;spinButton.setAttribute('aria-pressed',String(spinning));spinButton.textContent=spinning?'자동 회전 ON':'자동 회전 OFF';spinButton.setAttribute('aria-label',spinning?'자동 회전 끄기':'자동 회전 켜기')};
  const setSpinning=value=>{spinning=!!value;syncSpinButton()};syncSpinButton();if(spinButton)spinButton.addEventListener('click',()=>setSpinning(!spinning));
  host.addEventListener('pointerdown',e=>{if(e.pointerType==='mouse'&&e.button!==0)return;dragging=true;lastX=e.clientX;lastY=e.clientY;host.classList.add('is-dragging');if(host.setPointerCapture)host.setPointerCapture(e.pointerId)});
  host.addEventListener('pointermove',e=>{if(!dragging)return;const dx=e.clientX-lastX,dy=e.clientY-lastY;yaw-=dx*.008;pitch=Math.max(.28,Math.min(1.02,pitch+dy*.006));lastX=e.clientX;lastY=e.clientY;e.preventDefault()},{passive:false});
  const finishDrag=e=>{if(!dragging)return;dragging=false;host.classList.remove('is-dragging');if(host.hasPointerCapture&&host.hasPointerCapture(e.pointerId))host.releasePointerCapture(e.pointerId)};host.addEventListener('pointerup',finishDrag);host.addEventListener('pointercancel',finishDrag);
  let statusStamp=0;camera.position.set(6.6,7.7,9.1);const engine=makeRenderer(host,scene,camera,3,(now)=>{const dt=Math.min(40,Math.max(0,now-lastNow));lastNow=now;if(spinning&&!dragging)yaw+=dt*.000085;const horizontal=radius*Math.cos(pitch);camera.position.set(target.x+Math.sin(yaw)*horizontal,target.y+Math.sin(pitch)*radius,target.z+Math.cos(yaw)*horizontal);camera.lookAt(target);if(now-statusStamp>250||dragging){host.dataset.orbitYaw=yaw.toFixed(3);host.dataset.orbitPitch=pitch.toFixed(3);statusStamp=now}});engine.setSpinning=setSpinning;engine.getOrbitState=()=>({yaw,pitch,spinning,dragging});return engine}
function waveGeometry(width,depth,nx,nz,amplitude,fx,fz,phase){const p=[],idx=[];for(let z=0;z<=nz;z++){for(let x=0;x<=nx;x++){const px=(x/nx-.5)*width,pz=(z/nz-.5)*depth,py=amplitude*(Math.sin(px*fx+phase)+.45*Math.cos(pz*fz-phase));p.push(px,py,pz)}}for(let z=0;z<nz;z++){for(let x=0;x<nx;x++){const a=z*(nx+1)+x,b=a+1,c=a+nx+1,d=c+1;idx.push(a,c,b,b,c,d)}}const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(p,3));g.setIndex(idx);g.computeVertexNormals();return g}
function buildLayersWebGL(){
  const host=document.querySelector('#layersWebGL'),composite=document.querySelector('#layerComposite'),focusEl=document.querySelector('#layerFocus'),buttons=[...document.querySelectorAll('[data-field-layer]')],scene=new THREE.Scene(),camera=new THREE.PerspectiveCamera(37,1,.1,60);scene.fog=new THREE.FogExp2(0x050515,.045);camera.position.set(4.65,2.9,5.85);const root=new THREE.Group();root.rotation.y=-.28;scene.add(root);
  const mapFrame=new THREE.Group();root.add(mapFrame);addBoundary3D(mapFrame,{y:.08,scale:.73,color:0xeaf1f6,opacity:.8,centerLon:127.55,centerLat:35.65});const mapGrid=new THREE.GridHelper(5.5,14,0x5a6c86,0x211d46);mapGrid.scale.z=.67;mapGrid.position.y=.01;mapGrid.material.transparent=true;mapGrid.material.opacity=.14;mapFrame.add(mapGrid);
  const core=new THREE.Group(),regional=new THREE.Group(),crustal=new THREE.Group(),external=new THREE.Group(),groups={core,regional,crustal,external};Object.values(groups).forEach(group=>root.add(group));

  const earthMat=new THREE.MeshBasicMaterial({color:0x227b8b,transparent:true,opacity:.16,wireframe:true,depthWrite:false}),earth=new THREE.Mesh(new THREE.SphereGeometry(2.65,40,18,0,Math.PI*2,0,Math.PI*.5),earthMat);earth.position.y=-2.48;core.add(earth);const innerCore=new THREE.Mesh(new THREE.SphereGeometry(.72,24,14),new THREE.MeshBasicMaterial({color:0xff7048,transparent:true,opacity:.34,depthWrite:false}));innerCore.position.y=-2.48;core.add(innerCore);const coreLineMats=[];[-1.15,-.78,-.39,0,.39,.78,1.15].forEach((z,i)=>{const lift=1.02-Math.abs(z)*.12,curve=new THREE.CatmullRomCurve3([new THREE.Vector3(-2.65,-.05,z),new THREE.Vector3(-1.35,lift*.74,z*.78),new THREE.Vector3(0,lift,z*.58),new THREE.Vector3(1.35,lift*.74,z*.78),new THREE.Vector3(2.65,-.05,z)]),geometry=new THREE.BufferGeometry().setFromPoints(curve.getPoints(80)),material=new THREE.LineBasicMaterial({color:0x5af7f2,transparent:true,opacity:i===3?.42:.22,depthWrite:false});material.userData.baseOpacity=material.opacity;coreLineMats.push(material);core.add(new THREE.Line(geometry,material))});

  const plateGeometry=new THREE.PlaneGeometry(5.25,3.35),plateMaterial=new THREE.MeshBasicMaterial({color:0x7867ff,transparent:true,opacity:.13,side:THREE.DoubleSide,depthWrite:false}),regionalPlate=new THREE.Mesh(plateGeometry,plateMaterial);regionalPlate.rotation.x=-Math.PI/2;regionalPlate.position.y=.14;regional.add(regionalPlate);const plateEdge=new THREE.LineSegments(new THREE.EdgesGeometry(plateGeometry),new THREE.LineBasicMaterial({color:0x9e8cff,transparent:true,opacity:.55}));plateEdge.rotation.x=-Math.PI/2;plateEdge.position.y=.145;regional.add(plateEdge);regional.add(pointCloud(DATA.model_sites,{color:0xb8a9ff,size:.085,y:.23,scale:.73,centerLon:127.55,centerLat:35.65,opacity:.92}));

  const anomalyHeads=[];[{lon:126.25,lat:34.75,h:.58,c:0xff5c9c},{lon:127.25,lat:36.05,h:-.43,c:0x427df4},{lon:129.05,lat:35.35,h:.82,c:0xff5c9c},{lon:128.2,lat:37.05,h:.46,c:0xff8db8},{lon:126.75,lat:37.35,h:-.34,c:0x427df4}].forEach(a=>{const p=geoVector(a.lon,a.lat,.18,.73,127.55,35.65),height=Math.abs(a.h),material=new THREE.MeshBasicMaterial({color:a.c,transparent:true,opacity:.8}),stem=new THREE.Mesh(new THREE.CylinderGeometry(.045,.1,height,14),material);stem.position.set(p.x,.18+a.h/2,p.z);crustal.add(stem);const head=new THREE.Mesh(new THREE.SphereGeometry(.12,16,10),material);head.position.set(p.x,.18+a.h,p.z);head.userData.baseY=head.position.y;anomalyHeads.push(head);crustal.add(head);const ring=new THREE.Mesh(new THREE.TorusGeometry(.22,.014,6,40),new THREE.MeshBasicMaterial({color:a.c,transparent:true,opacity:.45}));ring.rotation.x=Math.PI/2;ring.position.set(p.x,.19,p.z);crustal.add(ring)});

  const externalPulse=new THREE.Group(),externalMats=[];external.add(externalPulse);[0,1].forEach(i=>{const material=new THREE.MeshBasicMaterial({color:0xffd36a,transparent:true,opacity:i?.16:.28,depthWrite:false}),ring=new THREE.Mesh(new THREE.TorusGeometry(2.15-i*.42,.025,8,96),material);ring.rotation.x=Math.PI/2;ring.position.y=1.5+i*.19;material.userData.baseOpacity=material.opacity;externalMats.push(material);externalPulse.add(ring)});[[-1.65,-.65],[0,.25],[1.65,-.55]].forEach(([x,z],i)=>{const arrow=new THREE.ArrowHelper(new THREE.Vector3(0,-1,0),new THREE.Vector3(x,1.72,z),1.16,0xffd36a,.18,.11);arrow.traverse(obj=>{if(obj.material){obj.material.transparent=true;obj.material.opacity=i===1?.65:.38;obj.material.userData.baseOpacity=obj.material.opacity;externalMats.push(obj.material)}});external.add(arrow)});

  const copy={all:['결합 보기','Core + Regional + Crustal · External은 관측시각 전처리'],core:['Core · 전 지구 기준장','지구 핵 기원 장파장 · IGRF-14'],regional:['Regional · 한반도 보정','16측점으로 정한 현재 0차 상수 보정'],crustal:['Crustal · 국지 자기이상','KIGAM 항공자력 · 총자력 F에 적용'],external:['External · 관측시각 변화','현재 예측층이 아닌 관측 D·I 전처리']};let active='all',engine=null;
  const applyFocus=key=>{active=key in groups?key:'all';Object.entries(groups).forEach(([name,group])=>{group.visible=active==='all'||active===name});buttons.forEach(btn=>btn.setAttribute('aria-pressed',String(btn.dataset.fieldLayer===active)));composite.classList.toggle('has-layer-focus',active!=='all');composite.dataset.fieldFocus=active;host.dataset.fieldFocus=active;root.scale.setScalar(({all:1.1,core:1.1,regional:1.4,crustal:1.5,external:1.32})[active]);const [title,note]=copy[active];focusEl.innerHTML=`<b>${title}</b><span>${note}</span>`;if(engine)engine.renderScene(performance.now(),true)};buttons.forEach(btn=>btn.addEventListener('click',()=>applyFocus(active===btn.dataset.fieldLayer?'all':btn.dataset.fieldLayer)));
  engine=makeRenderer(host,scene,camera,4,(now)=>{const t=now*.001;root.rotation.y=-.28+Math.sin(t*.16)*.025;coreLineMats.forEach((material,i)=>{material.opacity=material.userData.baseOpacity*(.82+.18*Math.sin(t*1.1+i*.7))});innerCore.scale.setScalar(1+Math.sin(t*.8)*.035);anomalyHeads.forEach((head,i)=>{head.position.y=head.userData.baseY+Math.sin(t*1.35+i)*.025});externalPulse.scale.setScalar(1+Math.sin(t*1.2)*.045);externalMats.forEach((material,i)=>{material.opacity=material.userData.baseOpacity*(.72+.28*Math.sin(t*1.45+i*.6))});camera.lookAt(0,({all:.02,core:-.5,regional:.12,crustal:.2,external:.68})[active],0)});engine.setLayerFocus=applyFocus;engine.getLayerFocus=()=>active;applyFocus('all');return engine
}
function setDensityView(view){const shell=document.querySelector('#densityShell');if(!shell)return;if(view==='3d'&&(!WEBGL.available||!WEBGL.enabled))view='2d';WEBGL.densityView=view;shell.classList.toggle('webgl-view',view==='3d');document.querySelectorAll('[data-density-view]').forEach(btn=>btn.setAttribute('aria-pressed',String(btn.dataset.densityView===view)));if(maps.density)maps.density.render();renderCinematicWebGL()}
function renderCinematicWebGL(){if(!WEBGL.available||!WEBGL.enabled)return;const active=window.activeScene||0;WEBGL.engines.forEach(engine=>{if(engine.sceneIndex===3&&WEBGL.densityView!=='3d')return;if(Math.abs(engine.sceneIndex-active)<=1)engine.renderScene(performance.now(),true)})}
window.renderCinematicWebGL=renderCinematicWebGL;
function startWebGLLoop(){cancelAnimationFrame(WEBGL.raf);if(!WEBGL.available||!WEBGL.enabled){return}const tick=now=>{if(!WEBGL.enabled)return;const active=window.activeScene||0;WEBGL.engines.forEach(engine=>{if(engine.sceneIndex===3&&WEBGL.densityView!=='3d')return;if(Math.abs(engine.sceneIndex-active)<=1)engine.renderScene(now,false)});if(!WEBGL.reduced)WEBGL.raf=requestAnimationFrame(tick)};tick(performance.now())}
function initWebGL(){const toggle=document.querySelector('#webglToggle');if(typeof THREE==='undefined'||!canUseWebGL()){document.body.classList.add('webgl-unavailable');return}try{buildHeroWebGL();buildDensityWebGL();buildLayersWebGL();WEBGL.available=WEBGL.engines.length>0;document.body.classList.add('webgl-ready');toggle.hidden=false;document.querySelectorAll('[data-density-view]').forEach(btn=>btn.addEventListener('click',()=>setDensityView(btn.dataset.densityView)));toggle.addEventListener('click',()=>{WEBGL.enabled=!WEBGL.enabled;document.body.classList.toggle('webgl-paused',!WEBGL.enabled);toggle.setAttribute('aria-pressed',String(WEBGL.enabled));toggle.textContent=WEBGL.enabled?'WEBGL ON':'WEBGL OFF';if(!WEBGL.enabled)setDensityView('2d');else if(innerWidth>780)setDensityView('3d');startWebGLLoop()});setDensityView(!WEBGL.reduced&&innerWidth>780?'3d':'2d');startWebGLLoop()}catch(error){console.warn('WebGL enhancement disabled:',error);WEBGL.available=false;document.body.classList.add('webgl-unavailable','webgl-paused');toggle.hidden=true;setDensityView('2d')}}

function initNavigation(){const nav=document.querySelector('#sceneNav');scenes.forEach((s,i)=>{const b=document.createElement('button');b.title=`${String(i+1).padStart(2,'0')} ${s.dataset.title}`;b.setAttribute('aria-label',b.title);b.addEventListener('click',()=>s.scrollIntoView({behavior:'smooth'}));nav.appendChild(b)});const dots=[...nav.children];
  const setActive=i=>{scenes.forEach((s,j)=>s.classList.toggle('active',i===j));dots.forEach((d,j)=>d.classList.toggle('active',i===j));document.querySelector('#counter').textContent=`${String(i+1).padStart(2,'0')} / ${String(scenes.length).padStart(2,'0')}`;document.querySelector('#progress').style.width=`${100*(i+1)/scenes.length}%`;window.activeScene=i;if(window.renderCinematicWebGL)window.renderCinematicWebGL()};
  const io=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting&&e.intersectionRatio>.42)setActive(scenes.indexOf(e.target))}),{threshold:[.42,.65]});scenes.forEach(s=>io.observe(s));setActive(0);
  const go=d=>scenes[Math.max(0,Math.min(scenes.length-1,(window.activeScene||0)+d))].scrollIntoView({behavior:'smooth'});document.querySelector('#prev').onclick=()=>go(-1);document.querySelector('#next').onclick=()=>go(1);
  addEventListener('keydown',e=>{if(/INPUT|TEXTAREA|SELECT/.test(e.target.tagName))return;if(['ArrowDown','PageDown',' '].includes(e.key)){e.preventDefault();go(1)}else if(['ArrowUp','PageUp'].includes(e.key)){e.preventDefault();go(-1)}else if(e.key==='Home'){scenes[0].scrollIntoView()}else if(e.key==='End'){scenes.at(-1).scrollIntoView()}});
}

function initField(){const c=document.querySelector('#fieldCanvas'),ctx=c.getContext('2d');let W=0,H=0,t=0;const reduced=matchMedia('(prefers-reduced-motion:reduce)').matches;function resize(){const r=c.getBoundingClientRect(),d=Math.min(2,devicePixelRatio||1);W=r.width;H=r.height;c.width=W*d;c.height=H*d;ctx.setTransform(d,0,0,d,0,0)}new ResizeObserver(resize).observe(c);resize();
  function frame(){ctx.clearRect(0,0,W,H);const cx=W*.71,cy=H*.47;for(let i=-9;i<=9;i++){const off=i*H*.043;ctx.beginPath();for(let x=0;x<=W;x+=14){const dx=x-cx;const y=cy+off+Math.sin(dx/W*5+i*.23)*22+(dx*dx)/(W*W)*off*.65;x===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.strokeStyle=i%3===0?'rgba(81,210,215,.16)':'rgba(82,116,140,.11)';ctx.lineWidth=i===0?1.2:.7;ctx.stroke()}
    for(let i=0;i<24;i++){const q=(t*.00005+i/24)%1,x=q*W,dx=x-cx,y=cy+(i-12)*H*.032+Math.sin(dx/W*5+i*.2)*22+(dx*dx)/(W*W)*(i-12)*H*.021;ctx.fillStyle=i%5===0?'#ff7048':'rgba(81,210,215,.66)';ctx.beginPath();ctx.arc(x,y,i%5===0?2.5:1.35,0,Math.PI*2);ctx.fill()}ctx.beginPath();ctx.arc(cx,cy,5,0,Math.PI*2);ctx.fillStyle='#ff7048';ctx.fill();t+=16;if(!reduced)requestAnimationFrame(frame)}frame();}

bindData();buildPerformance();initNavigation();initMaps();initField();initWebGL();
</script>
</body></html>'''


def render_template(template: str, payload: dict) -> str:
    """데이터와 Three.js를 주입하고 미치환 토큰이 없는지 확인한다."""
    text = template.replace("__THREE_JS__", read_three_runtime())
    text = text.replace(
        "__PAYLOAD__", json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    )
    leftovers = [token for token in ("__THREE_JS__", "__PAYLOAD__") if token in text]
    if leftovers:
        raise RuntimeError(f"Unreplaced cinematic placeholders: {leftovers}")
    return text


def main():
    payload = build_payload()
    text = render_template(HTML, payload)
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
