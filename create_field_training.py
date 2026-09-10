# -*- coding: utf-8 -*-
"""
현장 자료취득 요원 교육 — 웹 시네마틱 (`docs/field_training.html`)
====================================================================

    python create_field_training.py

현장 요원에게 **지자기 측량이 무엇이고 왜 당신의 기록 하나가 중요한가**를
전달하는 스크롤형 교육 자료. 22장면 · 약 35분(질의 별도).

## 이 자료가 다른 발표자료와 다른 점

`create_lmm_cinematic.py` 는 전문가·발주자용이라 「모델이 이만큼 나왔다」를 말한다.
이 자료는 **현장 요원용**이라 「당신이 이렇게 기록하지 않으면 이런 일이 벌어진다」를
말한다. 그래서 가상 사례를 쓰지 않고 **이 프로젝트가 실제로 겪은 일**만 싣는다:

| 장면 | 실화 |
|---|---|
| 시각을 안 적으면 | 야장 68개 전수조사에서 **F 측정시각 0건** → 외부장 보정이 원리적으로 불가능해졌다 |
| 좌표가 어긋나면 | 미원 353 m · 남양은 세 자료가 전부 다르고 **현재 좌표 미상** · 서산 시트에 남양 좌표가 통째로 복사 |
| 한 방문이 틀리면 | 부안·포천 2023 두 방문을 빼니 **전국 LOO 편각 −26 %** (0.7691 → 0.5675°) |
| 나쁜 자리를 고르면 | 이원 여의저수지 동쪽 **−250.72 nT/m** — 참고 기준 3 nT/m 의 83배 |

## ⚠️ 올해 요원이 하는 것은 «절대측정이 아니다»

발주자 확정(2026-09-07) — 올해 50점에서 하는 것은 **자기교란(구배) 측정**이고
DI-flux 절대측정은 선점이 끝난 뒤의 별도 단계다. 그래서 편각 4방향(E Up/W Up…)·
복각 4방향·Null 측정은 **「앞으로 이 자리에서 이런 관측을 합니다」라는 미래 맥락**
으로만 짧게 보이고, **올해 실제 작업**(P0 → 4방향 × 11점 → P0 · 수직 20~200 cm ·
사진 9칸 · 카드 7구역)을 3부 본론에 놓는다.

⚠️ 이 원칙을 어기면 요원이 「내가 DI-flux 를 배워야 하는구나」로 이해한다.

## 수치는 전부 파일에서 읽는다

`trial_survey_points`(50점) · `trial_survey_spec`(측선·참고값·발주자 결정) ·
`lmm_model.json`(LOO) · `existing_pts.geojson`(관측망 30) · `korea_boundary.geojson`.
**하드코딩 금지** — 재적합·명단 변경이 자동으로 따라온다.

⚠️ 숫자 네 갈래를 섞지 말 것 — **관측망 30** / 2022~25 측량 15 / LMM 투입 16 /
**올해 선점 검토 50**. 교육 대상이 헷갈리는 첫 번째 지점이다.

## 오프라인 단일 파일

Three.js r149(MIT)를 인라인하므로 **외부 script·CDN·fetch 가 없다**. 현장 교육장의
인터넷이 끊겨도 파일 하나로 뜬다. 지도는 타일 없이 해안선 GeoJSON 을 캔버스로 그린다.
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs" / "field_training.html"
THREE_RUNTIME = ROOT / "vendor" / "three.r149.min.js"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_three_runtime() -> str:
    """공식 Three.js 런타임을 단일 HTML 안에 인라인한다."""
    runtime = THREE_RUNTIME.read_text(encoding="utf-8")
    for token in ("SPDX-License-Identifier: MIT", 'const e="149"'):
        if token not in runtime:
            raise RuntimeError(f"예상과 다른 Three.js 런타임: {token} 없음")
    if "</script" in runtime.lower():
        raise RuntimeError("Three.js 런타임에 닫는 script 태그가 있다")
    return runtime


def flatten_rings(value, out):
    if not value:
        return
    if isinstance(value[0], (int, float)):
        return
    if isinstance(value[0][0], (int, float)):
        out.append([[round(float(x), 3), round(float(y), 3)] for x, y in value])
        return
    for part in value:
        flatten_rings(part, out)


def compact_boundary(doc):
    rings = []
    for feat in doc.get("features", []):
        flatten_rings((feat.get("geometry") or {}).get("coordinates"), rings)
    rings.sort(key=len, reverse=True)
    return rings[:6]


# ══════════════════════════════════════════════════════════════
def build_payload():
    import trial_survey_points as TP
    import trial_survey_spec as SP

    pts = TP.load_points()
    for i, p in enumerate(pts, 1):
        p.setdefault("_sid", f"S{i:02d}")
    kinds = Counter(p["구분"] for p in pts)

    model = read_json(DATA / "lmm_model.json")
    loo = {k: float(v) for k, v in model["loo_cv"].items()}
    boundary = compact_boundary(read_json(DATA / "korea_boundary.geojson"))

    net = read_json(DATA / "existing_pts.geojson")
    network = [
        {
            "name": f["properties"].get("name"),
            "lat": round(float(f["geometry"]["coordinates"][1]), 5),
            "lon": round(float(f["geometry"]["coordinates"][0]), 5),
        }
        for f in net.get("features", [])
    ]

    sites = [
        {
            "id": p["_sid"],
            "k": p["구분"],
            "name": p["지점명"],
            "lat": round(p["위도"], 5),
            "lon": round(p["경도"], 5),
            "grad": None if p["예측구배"] is None else round(p["예측구배"], 1),
        }
        for p in pts
    ]

    pilot = [
        {"name": nm, "dirs": v, "note": note}
        for nm, v, note in SP.PILOT if v
    ]
    worst = min(
        (abs(x), k, nm) for nm, v, _ in SP.PILOT if v for k, x in v.items()
    )
    worst_max = max(
        (abs(x), k, nm) for nm, v, _ in SP.PILOT if v for k, x in v.items()
    )

    nH = len(SP.H_DIRECTIONS) * (len(SP.H_OFFSETS_M) + 2)
    nV = len(SP.V_HEIGHTS_CM) + 2

    return {
        "generated": f"{dt.date.today():%Y년 %m월 %d일}",
        "boundary": boundary,
        "network": network,
        "sites": sites,
        "n_sites": len(pts),
        "n_new": kinds["A"] + kinds["B"],
        "n_exist": kinds["기존점"],
        "n_grade_a": kinds["A"],
        "n_grade_b": kinds["B"],
        "n_network": len(network),
        "target_network": len(network) + 50,
        "loo_d": round(loo["D"], 4),
        "loo_d_min": round(loo["D"] * 60, 1),
        "loo_i": round(loo["I"], 4),
        "loo_f": round(loo["F"], 1),
        "offsets": SP.H_OFFSETS_M,
        "dirs": SP.H_DIRECTIONS,
        "heights": SP.V_HEIGHTS_CM,
        "n_read_h": nH,
        "n_read_v": nV,
        "n_read": nH + nV,
        "iaga_range": SP.IAGA_RANGE_NT,
        "iaga_radius": SP.IAGA_RADIUS_M,
        "euro_grad": SP.EURO_GRAD_NT_PER_M,
        "pilot": pilot,
        "worst_val": round(worst_max[0], 2),
        "worst_dir": worst_max[1],
        "worst_site": worst_max[2],
        "worst_ratio": round(worst_max[0] / SP.EURO_GRAD_NT_PER_M),
        "decisions": [
            {"k": k, "v": v.replace("**", "")} for k, v, _ in SP.DECISIONS
        ],
        "gap_sites": [p["지점명"] for p in pts if p["예측구배"] is None],
    }


# ══════════════════════════════════════════════════════════════
TEMPLATE = r"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>지자기 측량 현장교육 — 당신의 기록 하나에서 시작합니다</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
:root{
 --bg:#050a12; --bg2:#08121e; --ink:#eaf2f8; --muted:#8ba3b8;
 --cyan:#3fd8d0; --blue:#4a7fe8; --violet:#9b7fe8; --orange:#ff7048;
 --red:#e8503f; --line:rgba(140,180,210,.18); --good:#3fd8a0;
}
html{scroll-behavior:smooth}
body{background:var(--bg);color:var(--ink);
 font-family:"Pretendard","Noto Sans KR","맑은 고딕",system-ui,sans-serif;
 line-height:1.7;word-break:keep-all;overflow-wrap:normal;line-break:strict;
 -webkit-font-smoothing:antialiased}
.no-break{white-space:nowrap}
/* ⚠️ 배경 자기력선은 «분위기»지 내용이 아니다. 진하면 본문을 덮어
   현장 교육장 빔프로젝터에서 글자가 안 읽힌다 — 실제로 그랬다. */
#bg{position:fixed;inset:0;z-index:0;opacity:.34}
#bg canvas{display:block;width:100%;height:100%}
#bg::after{content:"";position:absolute;inset:0;
 background:radial-gradient(ellipse at 30% 45%,rgba(5,10,18,.92) 0%,
 rgba(5,10,18,.62) 45%,rgba(5,10,18,.25) 100%)}
main{position:relative;z-index:1}
section{min-height:100vh;display:flex;align-items:center;
 padding:clamp(48px,7vh,110px) clamp(20px,6vw,96px);position:relative}
.wrap{width:100%;max-width:1180px;margin:0 auto}
.tag{font:700 12px/1 ui-monospace,Consolas,monospace;letter-spacing:.18em;
 color:var(--cyan);margin-bottom:20px;display:block}
.tag.o{color:var(--orange)}
h1{font-size:clamp(34px,6.2vw,82px);line-height:1.16;font-weight:800;
 letter-spacing:-.02em;text-wrap:balance}
h2{font-size:clamp(26px,3.9vw,50px);line-height:1.24;font-weight:800;
 letter-spacing:-.015em;margin-bottom:22px;text-wrap:balance}
h3{font-size:clamp(18px,2vw,25px);font-weight:700;margin-bottom:12px}
.lead{font-size:clamp(16px,1.55vw,21px);color:var(--muted);
 max-width:64ch;margin-top:18px;text-wrap:pretty}
p+p{margin-top:14px}
.hl{color:var(--cyan)}.hl-o{color:var(--orange)}.hl-r{color:var(--red)}
.hl-g{color:var(--good)}
b{color:#fff;font-weight:700}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:clamp(22px,4vw,58px);
 align-items:center}
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}
.grid4{display:grid;grid-template-columns:repeat(4,1fr);gap:13px}
/* ⚠️ 그리드 자식은 기본이 min-width:auto 라 «내용 폭»만큼 벌어진다.
   그래서 안에 넣은 overflow-x 상자가 스크롤하지 않고 페이지를 민다 —
   좁은 화면에서 표가 실제로 그랬다. 0 으로 눌러 줘야 상자가 일한다. */
.grid2>*,.grid3>*,.grid4>*{min-width:0}
@media(max-width:860px){.grid2,.grid3,.grid4{grid-template-columns:1fr}}
.card{padding:22px 20px;border:1px solid var(--line);
 background:rgba(255,255,255,.025);border-radius:2px}
.card.warn{border-color:rgba(232,80,63,.45);background:rgba(232,80,63,.07)}
.card.good{border-color:rgba(63,216,160,.4);background:rgba(63,216,160,.06)}
.card h3{margin-bottom:8px}
.card p{color:var(--muted);font-size:15px;line-height:1.62}
.num{font:800 clamp(30px,4.4vw,64px)/1 "Arial Narrow","Pretendard",sans-serif;
 letter-spacing:-.03em;display:block}
.num.c{color:var(--cyan)}.num.o{color:var(--orange)}.num.r{color:var(--red)}
.unit{font-size:.42em;color:var(--muted);margin-left:5px;font-weight:600}
.stat small{display:block;margin-top:9px;color:var(--muted);font-size:14px;
 line-height:1.5}
.quote{border-left:3px solid var(--cyan);padding:6px 0 6px 24px;
 font-size:clamp(19px,2.4vw,32px);line-height:1.42;font-weight:700;
 margin:30px 0;text-wrap:balance}
.quote.o{border-color:var(--orange)}
.quote.r{border-color:var(--red)}
.tblwrap{overflow-x:auto;-webkit-overflow-scrolling:touch;margin-top:20px}
table{width:100%;border-collapse:collapse;font-size:15px;min-width:420px}
th,td{padding:11px 13px;text-align:left;border-bottom:1px solid var(--line);
 vertical-align:top}
th{font-size:12px;letter-spacing:.08em;color:var(--muted);font-weight:700;
 text-transform:uppercase}
td b{color:#fff}
.bad td{color:#f0b5ad}
canvas.fig{width:100%;height:auto;display:block;border:1px solid var(--line);
 background:#04080e;border-radius:2px}
.steps{counter-reset:s;margin-top:26px;display:grid;gap:11px}
.step{display:grid;grid-template-columns:44px 1fr;gap:16px;align-items:start;
 padding:15px 17px;border:1px solid var(--line);background:rgba(255,255,255,.02)}
.step::before{counter-increment:s;content:counter(s);
 font:800 20px/1.5 ui-monospace,monospace;color:var(--cyan);text-align:center}
.step b{display:block;margin-bottom:3px}
.step small{color:var(--muted);font-size:14px;line-height:1.55}
.rules{display:grid;gap:12px;margin-top:34px;counter-reset:r}
.rule{display:grid;grid-template-columns:52px 1fr;gap:18px;align-items:center;
 padding:19px 22px;border:1px solid var(--line);
 background:linear-gradient(90deg,rgba(63,216,208,.07),transparent 70%)}
.rule::before{counter-increment:r;content:counter(r);
 font:800 27px/1 ui-monospace,monospace;color:var(--cyan);text-align:center}
.rule b{font-size:clamp(17px,1.9vw,23px)}
.rule small{display:block;color:var(--muted);font-size:14.5px;margin-top:4px}
.chip{display:inline-block;padding:5px 12px;border:1px solid var(--line);
 font:700 12px/1.5 ui-monospace,monospace;color:var(--muted);margin:0 7px 7px 0;
 border-radius:2px}
.chip.on{border-color:var(--cyan);color:var(--cyan)}
.chip.off{border-color:rgba(232,80,63,.5);color:#f0b5ad}
.prog{position:fixed;left:0;top:0;height:2px;background:var(--cyan);z-index:50;
 width:0;transition:width .12s linear}
.nav{position:fixed;right:16px;top:50%;transform:translateY(-50%);z-index:40;
 display:flex;flex-direction:column;gap:7px}
.nav i{width:7px;height:7px;border:1px solid var(--muted);border-radius:50%;
 cursor:pointer;transition:.2s;display:block}
.nav i.on{background:var(--cyan);border-color:var(--cyan);transform:scale(1.35)}
@media(max-width:860px){.nav{display:none}}
.reveal{opacity:0;transform:translateY(22px);
 transition:opacity .7s ease,transform .7s cubic-bezier(.2,.7,.3,1)}
.reveal.in{opacity:1;transform:none}
.foot{padding:64px clamp(20px,6vw,96px) 90px;color:var(--muted);font-size:13px;
 border-top:1px solid var(--line);text-align:center;line-height:1.9}
@media (prefers-reduced-motion:reduce){
 .reveal{opacity:1;transform:none;transition:none}
 html{scroll-behavior:auto}
}
</style></head><body>
<div class="prog" id="prog"></div>
<div class="nav" id="nav"></div>
<div id="bg"></div>
<main>

<!-- ══ 프롤로그 ══ -->
<section data-t="프롤로그"><div class="wrap reveal">
 <span class="tag">현장 자료취득 요원 교육 · {{GENERATED}}</span>
 <h1>우리는 보이지 않는<br>국가 기준을 측정합니다</h1>
 <p class="lead">지구자기장은 눈에 보이지 않지만 지도의 북쪽, 나침반의 방향,
 항법 장비의 기준이 됩니다. 그 값은 <b>현장에서 사람이 직접 재야만</b>
 알 수 있습니다.</p>
 <p class="lead">오늘 이야기는 하나입니다 —
 <b class="hl">전국 지자기 모델의 품질은 현장에서 기록한 단 한 번의 관측에서
 시작됩니다.</b></p>
</div></section>

<section data-t="편각"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">01</span>
  <h2>지도의 북쪽과<br>나침반의 북쪽은 다릅니다</h2>
  <p class="lead">우리나라에서 나침반 바늘은 진북보다 <b>서쪽으로 약 8도</b>
  기울어 있습니다. 이 각도를 <b class="hl">편각(D)</b> 이라 부릅니다.</p>
  <p class="lead">8도는 작아 보이지만 1 km 를 가면 140 m 가 어긋납니다.
  그래서 국가기본도에 이 값을 인쇄하고, 정확히 관리하려면 전국에서
  실제로 재야 합니다.</p>
 </div>
 <div class="reveal"><canvas class="fig" id="figCompass" width="620" height="620"
  aria-label="진북과 자북의 차이"></canvas></div>
</div></section>

<section data-t="올해 50점"><div class="wrap">
 <div class="reveal">
  <span class="tag">02</span>
  <h2>올해 여러분이 갈 곳은<br><span class="hl">{{N_SITES}}곳</span>입니다</h2>
 </div>
 <div class="grid2" style="margin-top:34px">
  <div class="reveal"><canvas class="fig" id="figSites" width="640" height="800"
   aria-label="올해 조사 대상 50지점 분포"></canvas></div>
  <div class="reveal">
   <div class="grid2" style="gap:13px">
    <div class="card stat"><span class="num c">{{N_NEW}}</span>
     <small>신규 후보지 — 도상선점 103곳에서 현장조사를 거쳐 좁힌 자리</small></div>
    <div class="card stat"><span class="num o">{{N_EXIST}}</span>
     <small>기존 지자기점 — 관측망 {{N_NETWORK}}점 중 표석이 계속 유지될
     수 있는 자리</small></div>
   </div>
   <p class="lead" style="margin-top:24px">이 {{N_SITES}}곳에서 여러분이 할 일은
   <b>지자기 측량이 아닙니다.</b> 그 자리가 자기적으로 조용한지를 재는
   <b class="hl">자기교란 조사</b> 입니다. 왜 그것부터 하는지는 3부에서
   말씀드리겠습니다.</p>
  </div>
 </div>
</div></section>

<!-- ══ 1부 왜 ══ -->
<section data-t="1부 · 네 겹"><div class="wrap">
 <div class="reveal">
  <span class="tag">1부 — 왜 이 일을 하는가 · 03</span>
  <h2>우리가 재는 값은<br>한 가지 자기장이 아닙니다</h2>
  <p class="lead">현장에서 자력계에 뜨는 숫자 하나는 네 가지가 겹쳐진
  결과입니다. 이걸 갈라내는 것이 이 연구의 핵심입니다.</p>
 </div>
 <div class="grid4 reveal" style="margin-top:38px">
  <div class="card"><h3 class="hl">① 주 자기장</h3>
   <p>지구 외핵의 유체 운동에서 생깁니다. 전체의 <b>99 %</b> 를 차지하고
   해마다 조금씩 움직입니다.</p></div>
  <div class="card"><h3 class="hl">② 지역장</h3>
   <p>전지구 모델이 못 담는 우리나라만의 차이입니다.
   <b>이 값을 만드는 것이 이 연구</b> 입니다.</p></div>
  <div class="card"><h3 class="hl">③ 지각 자기장</h3>
   <p>땅속 암석이 만듭니다. 몇 m 만 옮겨도 값이 달라지는 것이
   대부분 이것 때문입니다.</p></div>
  <div class="card"><h3 class="hl-o">④ 외부 자기장</h3>
   <p>태양 활동으로 <b>하루 사이에도 변합니다.</b> 그래서 «몇 시 몇 분»에
   쟀는지가 반드시 필요합니다.</p></div>
 </div>
 <div class="quote reveal">④ 때문에 <span class="hl-o">시각을 안 적으면 그 값은
 쓸 수 없습니다.</span> 뒤에서 실제로 그런 일이 있었습니다.</div>
</div></section>

<section data-t="변한다"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">04</span>
  <h2>한 번 재고 끝낼 수 없습니다</h2>
  <p class="lead">지구자기장은 해마다 변합니다. 편각은 우리나라에서
  <b>연간 약 3분씩</b> 서쪽으로 움직입니다. 10년이면 0.5도입니다.</p>
  <p class="lead">그래서 과거 성과만으로는 오늘의 자기장을 설명할 수 없고,
  <b class="hl">주기적으로 다시 재야</b> 합니다. 여러분이 고른 자리에서
  <b>앞으로 5년마다</b> 관측이 이어집니다.</p>
 </div>
 <div class="reveal"><canvas class="fig" id="figDrift" width="640" height="440"
  aria-label="편각의 연간 변화"></canvas></div>
</div></section>

<section data-t="관측망"><div class="wrap">
 <div class="reveal">
  <span class="tag">05</span>
  <h2>지금 관측망에는<br>빈 곳이 있습니다</h2>
  <p class="lead">전국 1등 지자기점은 <b>{{N_NETWORK}}점</b> 입니다.
  그런데 최근 측량은 그중 15점에만 집중됐고, 강원·충청 내륙 축이
  빠지면서 <b class="hl-o">가장 먼 빈 곳이 129 km</b> 까지 벌어졌습니다.</p>
 </div>
 <div class="grid3 reveal" style="margin-top:32px">
  <div class="card stat"><span class="num">{{N_NETWORK}}</span>
   <small>1등 지자기점 관측망 — 국가가 관리하는 반복관측점 전체</small></div>
  <div class="card stat"><span class="num o">129<span class="unit">km</span></span>
   <small>가장 가까운 관측점까지 이만큼 떨어진 곳이 있습니다(강원 북동)</small></div>
  <div class="card stat"><span class="num c">{{TARGET_NETWORK}}</span>
   <small>목표 관측망 — 현 {{N_NETWORK}}점에 신규 50점을 더한 설계</small></div>
 </div>
 <div class="quote reveal o">관측점이 부족하면 그 지역의 값은
 <b>측정이 아니라 추정</b>이 됩니다.</div>
</div></section>

<section data-t="현재 성능"><div class="wrap">
 <div class="reveal">
  <span class="tag">06</span>
  <h2>지금 모델은<br>얼마나 맞습니까</h2>
  <p class="lead">현재 우리 모델의 편각 예측 오차는
  <b class="hl-o">{{LOO_D_MIN}}분</b> 입니다. 목표는 6분인데,
  아직 다섯 배 떨어져 있습니다.</p>
 </div>
 <div class="grid3 reveal" style="margin-top:30px">
  <div class="card stat"><span class="num o">{{LOO_D_MIN}}<span class="unit">분</span></span>
   <small>편각 예측 오차 — 모델에 넣지 않은 지점에서 재어 본 값</small></div>
  <div class="card stat"><span class="num">{{LOO_F}}<span class="unit">nT</span></span>
   <small>총자력 예측 오차</small></div>
  <div class="card stat"><span class="num c">16</span>
   <small>모델에 실제로 들어간 관측점 수 — 이것이 부족의 핵심입니다</small></div>
 </div>
 <p class="lead reveal" style="margin-top:26px">원인을 찾아보니 장비도 계산도
 문제가 없었습니다. 같은 날 두 번 재면 <b class="hl-g">1.4분</b> 안에 들어옵니다.
 문제는 <b class="hl-r">몇 년 뒤 다시 가서 재면 34분이 어긋난다</b>는 것이었습니다.
 그 이야기를 4부에서 하겠습니다.</p>
</div></section>

<!-- ══ 2부 무엇을 ══ -->
<section data-t="2부 · 벡터"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">2부 — 무엇을 재는가 · 07</span>
  <h2>자기장은 방향과 세기를<br>함께 가진 화살표입니다</h2>
  <p class="lead">한 지점의 자기장은 3차원 화살표 하나로 나타납니다.
  이 화살표를 완전히 알려면 <b>세 가지</b>만 재면 됩니다.</p>
  <div style="margin-top:22px">
   <span class="chip on">D 편각</span><span class="chip on">I 복각</span>
   <span class="chip on">F 총자력</span>
   <span class="chip">→</span>
   <span class="chip">X 북</span><span class="chip">Y 동</span>
   <span class="chip">Z 연직</span><span class="chip">H 수평</span>
  </div>
  <p class="lead">셋을 재면 나머지 넷은 계산으로 나옵니다.
  <b class="hl">그래서 이 셋을 정확히 재는 것이 지자기 측량입니다.</b></p>
 </div>
 <div class="reveal"><canvas class="fig" id="figVector" width="640" height="620"
  aria-label="자기장 벡터의 성분 분해"></canvas></div>
</div></section>

<section data-t="D·I·F"><div class="wrap">
 <div class="reveal"><span class="tag">08</span>
 <h2>세 가지를 하나씩 보면</h2></div>
 <div class="grid3 reveal" style="margin-top:30px">
  <div class="card"><h3 class="hl">편각 D</h3>
   <p>위에서 내려다봤을 때 <b>진북과 자북 사이의 각도</b> 입니다.
   우리나라는 서쪽으로 약 8도. 재기 가장 어렵고, 오차도 여기서 제일 많이
   납니다.</p></div>
  <div class="card"><h3 class="hl">복각 I</h3>
   <p>옆에서 봤을 때 자기장이 <b>수평면 아래로 기울어진 각도</b> 입니다.
   우리나라는 약 53도. 적도에서 0도, 극에서 90도입니다.</p></div>
  <div class="card good"><h3 class="hl-g">총자력 F</h3>
   <p>화살표 <b>전체의 세기</b> 입니다. 약 5만 nT.
   <b class="hl-g">올해 여러분이 재는 것이 바로 이 F 하나입니다.</b></p></div>
 </div>
 <div class="quote reveal">D 와 I 는 비자성 경위의와 플럭스게이트로 재는
 <b>절대측정</b>입니다. 그건 <span class="hl-o">선점이 끝난 뒤</span>의
 일입니다.</div>
</div></section>

<!-- ══ 3부 올해 할 일 ══ -->
<section data-t="3부 · 올해 할 일"><div class="wrap">
 <div class="reveal">
  <span class="tag o">3부 — 올해 여러분이 하는 일 · 09</span>
  <h2>측량이 아니라<br><span class="hl-o">자리를 고르는 일</span> 입니다</h2>
  <p class="lead">지자기 측량은 그 자리가 <b>조용한 자리일 때만</b> 뜻이 있습니다.
  땅속 암석이나 주변 철구조물 때문에 몇 m 만 움직여도 값이 뛰는 자리라면,
  아무리 정밀하게 재도 그 값은 그 지역을 대표하지 못합니다.</p>
  <p class="lead">그래서 <b class="hl">먼저 자리가 조용한지를 재고</b>,
  그다음에 그 자리를 국가 기준점으로 삼습니다. 올해가 그 첫 단계입니다.</p>
 </div>
 <div class="grid2 reveal" style="margin-top:34px">
  <div class="card good"><h3>올해 하는 것</h3>
   <p><b>총자력 F 만</b> 재서 그 자리 주변이 얼마나 균질한지 봅니다.
   오버하우저 자력계 한 대와 GNSS 면 됩니다.</p></div>
  <div class="card"><h3 style="color:var(--muted)">올해 하지 않는 것</h3>
   <p>편각·복각 절대측정(DI-flux). 4방향 관측이나 Null 측정은
   <b>이번 범위가 아닙니다.</b> 선점이 확정된 뒤 별도로 진행합니다.</p></div>
 </div>
</div></section>

<section data-t="측선 설계"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag o">10</span>
  <h2>중심점에서 네 방향,<br>그리고 위로</h2>
  <p class="lead"><b>수평</b> — 중심점(P0)에서 재고, 동쪽으로
  {{OFFSETS_TXT}} m 를 차례로 잰 뒤 <b class="hl">다시 중심점으로 돌아와
  한 번 더</b> 잽니다. 이걸 동·서·남·북 네 방향 반복합니다.</p>
  <p class="lead"><b>수직</b> — 중심점에서 센서 높이만 바꿔가며
  지상 {{H_MIN}}~{{H_MAX}} cm 를 {{H_STEP}} cm 간격으로 잽니다.</p>
  <p class="lead">한 지점에서 모두 <b class="hl">{{N_READ}}번</b>
  재게 됩니다(수평 {{N_READ_H}} + 수직 {{N_READ_V}}).</p>
 </div>
 <div class="reveal"><canvas class="fig" id="figProfile" width="640" height="640"
  aria-label="수평 4방향 측선과 수직 측선 배치"></canvas></div>
</div></section>

<section data-t="P0 재측정"><div class="wrap">
 <div class="reveal">
  <span class="tag o">11</span>
  <h2>왜 중심점으로<br>자꾸 돌아옵니까</h2>
  <p class="lead">자력계가 <b>한 대뿐</b>이기 때문입니다. 한 대로 여러 자리를 재면
  <b class="hl-o">위치가 바뀌는 동안 시간도 함께 흘러갑니다.</b> 그 사이 자기장은
  저절로 변합니다.</p>
  <p class="lead">중심점을 앞뒤로 두 번 재두면, 그 두 값을 직선으로 이어
  «각 측점을 잰 순간의 중심점 값»을 추정할 수 있습니다. 그만큼을 빼면
  <b class="hl">순수하게 자리 때문에 생긴 차이</b>만 남습니다.</p>
 </div>
 <div class="grid2 reveal" style="margin-top:30px">
  <div><canvas class="fig" id="figCorr" width="640" height="380"
   aria-label="시간변화 보정 예시"></canvas></div>
  <div>
   <div class="tblwrap"><table>
    <tr><th>순서</th><th>시각</th><th>값</th></tr>
    <tr><td>중심점 P0</td><td>10:00</td><td><b>100 nT</b></td></tr>
    <tr><td>동쪽 2 m</td><td>10:02</td><td><b>106 nT</b></td></tr>
    <tr><td>중심점 P0</td><td>10:04</td><td><b>105 nT</b></td></tr>
   </table></div>
   <div class="card good" style="margin-top:20px">
    <p style="color:var(--ink)">10:02 의 중심점 추정값 = <b>102.5</b> nT<br>
    차이 ΔF = 106 − 102.5 = <b>3.5</b> nT<br>
    구배 = 3.5 ÷ 2 m = <b class="hl-g">1.75 nT/m</b></p>
   </div>
   <div class="card warn" style="margin-top:13px">
    <p>보정하지 않으면 (106−100)÷2 = <b class="hl-r">3.0 nT/m</b> —
    <b>71 % 과대평가</b> 됩니다. 자리 탓이 아니라 시간이 흐른 몫까지
    자리 탓으로 돌리게 됩니다.</p>
   </div>
  </div>
 </div>
</div></section>

<section data-t="실측 사례"><div class="wrap">
 <div class="reveal">
  <span class="tag o">12</span>
  <h2>실제로 재어 봤더니</h2>
  <p class="lead">기존 지자기점 두 곳에서 시범으로 재어 본 결과입니다.
  중심점에서 1 m 떨어진 네 방향 값입니다.</p>
 </div>
 <div class="reveal" style="margin-top:26px">{{PILOT_TABLE}}</div>
 <div class="quote reveal r">한 곳은 동쪽으로 1 m 만 가도
 <span class="hl-r">{{WORST_VAL}} nT/m</span> 씩 값이 뜁니다 —
 국제 권고 기준 {{EURO_GRAD}} nT/m 의 <b>{{WORST_RATIO}}배</b> 입니다.</div>
 <p class="lead reveal">이런 자리에서는 정밀하게 재도 소용이 없습니다.
 표석을 몇 cm 만 잘못 찾아도 값이 달라지니까요.
 <b class="hl">그 자리가 이런 자리인지 아닌지를 가려내는 것</b>이
 올해 여러분의 일입니다.</p>
</div></section>

<section data-t="기준이 없다"><div class="wrap">
 <div class="reveal">
  <span class="tag o">13</span>
  <h2>그런데 아직<br><span class="hl-o">합격 기준이 없습니다</span></h2>
  <p class="lead">국제 권고에는 「반경 {{IAGA_RADIUS}} m 안에서
  {{IAGA_RANGE}} nT 이내」「구배 {{EURO_GRAD}} nT/m 미만」이라는 값이 있습니다.
  그런데 <b class="hl-r">국내 기준은 아직 정해져 있지 않습니다.</b></p>
  <p class="lead">앞 장의 시범관측을 보십시오. 두 곳 모두 네 방향 전부가
  {{EURO_GRAD}} nT/m 를 넘었습니다. 이 기준을 그대로 쓰면 상당수가 탈락합니다.
  그 기준이 우리 땅에 맞는지 아직 아무도 모릅니다.</p>
 </div>
 <div class="quote reveal">그래서 이번 조사는 기준을 <b>지키는</b> 일이 아니라
 <span class="hl">기준을 만들 자료를 모으는</span> 일입니다.</div>
 <div class="card warn reveal" style="margin-top:8px">
  <p style="font-size:17px;color:var(--ink)">값이 크게 나와도
  <b>현장에서 그 자리를 접지 마십시오.</b> 나온 대로 적어 주셔야 합니다.
  큰 값만 골라 빼거나 다시 재면, 나중에 기준을 정할 때 자료가 한쪽으로
  치우쳐 <b class="hl-r">엉뚱한 기준</b>이 만들어집니다.</p>
 </div>
</div></section>

<section data-t="자기청정"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag o">14</span>
  <h2>측정하는 사람도<br>자기장을 만듭니다</h2>
  <p class="lead">센서 가까이에 쇠붙이가 있으면 그 영향이 값에 섞여 들어옵니다.
  문제는 나중에 보면 <b class="hl-r">그것이 사람 때문인지 그 자리의 성질인지
  가려낼 방법이 없다</b>는 것입니다.</p>
  <p class="lead">측정 전에 몸에서 빼 주십시오. 사정상 못 뺀 것이 있으면
  무엇을 왜 못 뺐는지 적어 주시면 됩니다.</p>
 </div>
 <div class="reveal"><canvas class="fig" id="figClean" width="620" height="560"
  aria-label="측정자가 지닌 자성 물품"></canvas></div>
</div></section>

<section data-t="카드"><div class="wrap">
 <div class="reveal">
  <span class="tag o">15</span>
  <h2>현장에서는 이 카드<br>한 장을 채웁니다</h2>
  <p class="lead">점마다 카드가 한 장씩 있습니다. 한 장을 다 채우면
  그 지점의 기록이 끝납니다. 순서대로 되어 있으니 위에서부터 채워 가시면 됩니다.</p>
 </div>
 <div class="steps reveal">
  <div class="step"><div><b>측정 전 확인</b>
   <small>몸에 지닌 자성 물품, 차량 위치, 전자기기 전원</small></div></div>
  <div class="step"><div><b>도착 · 관측 조건 · 원시자료 연결</b>
   <small>일자·시각·관측자·기상·Kp, 그리고 <b>장비 원시파일명과 레코드 범위</b></small></div></div>
  <div class="step"><div><b>기준점 좌표</b>
   <small>GNSS 실측 위도·경도를 십진도로</small></div></div>
  <div class="step"><div><b>수평 자기구배</b>
   <small>네 방향 × {{N_OFF}}점, 방향마다 P0 앞뒤로 한 번씩</small></div></div>
  <div class="step"><div><b>수직 자기구배</b>
   <small>중심점에서 높이별, 실제 센서 «가운데» 높이를 적습니다</small></div></div>
  <div class="step"><div><b>방위표지 시준</b>
   <small>정·반 시준과 표지 좌표·취득방법</small></div></div>
  <div class="step"><div><b>현장 사진 9칸</b>
   <small>중심점·측정 모습·표지 2장·네 방향 전경·교란 요소</small></div></div>
  <div class="step"><div><b>중심점 최종 결정 · 소견</b>
   <small>더 조용한 자리로 옮겼다면 사유와 함께</small></div></div>
 </div>
</div></section>

<!-- ══ 4부 왜 당신이 중요한가 ══ -->
<section data-t="4부 · 실화"><div class="wrap">
 <div class="reveal">
  <span class="tag" style="color:var(--red)">4부 — 실제로 있었던 일 · 16</span>
  <h2>시각을 안 적어서<br><span class="hl-r">보정을 통째로 못 한 일</span></h2>
  <p class="lead">지난 자료를 정리하면서 야장 <b>68권을 전부</b> 열어 봤습니다.
  총자력을 몇 시 몇 분에 쟀는지 적힌 것이
  <b class="hl-r">한 건도 없었습니다.</b></p>
 </div>
 <div class="grid3 reveal" style="margin-top:30px">
  <div class="card stat"><span class="num">68</span>
   <small>확인한 야장 권수</small></div>
  <div class="card stat warn"><span class="num r">0</span>
   <small>총자력 측정시각이 적힌 건수</small></div>
  <div class="card stat"><span class="num o">불가</small></span>
   <small>외부장 보정 — 시각이 없으니 어느 시점 값을 빼야 할지 알 수 없습니다</small></div>
 </div>
 <div class="quote reveal r">각도를 아무리 정확히 재도,
 <b>몇 시에 쟀는지 모르면 그 값은 쓸 수 없습니다.</b></div>
 <p class="lead reveal">규정에도 「편각·복각을 잴 때 총자기장과 시간을 함께
 측정한다」고 되어 있었습니다. 다만 <b>야장 양식에 그 칸이 없었습니다.</b>
 그래서 이번 카드에는 <b class="hl">F 를 잰 모든 줄에 시각 칸</b>을 넣었습니다.</p>
</div></section>

<section data-t="좌표"><div class="wrap">
 <div class="reveal">
  <span class="tag" style="color:var(--red)">17</span>
  <h2>좌표가 어긋나서<br><span class="hl-r">점 하나를 잃은 일</span></h2>
 </div>
 <div class="reveal tblwrap" style="margin-top:24px"><table>
  <tr><th>점</th><th>무슨 일이 있었나</th><th>결과</th></tr>
  <tr><td><b>미원</b></td>
   <td>두 문서의 좌표가 <b>353 m</b> 어긋납니다</td>
   <td>같은 표석인지 확인될 때까지 과거 성과와 대조 보류</td></tr>
  <tr class="bad"><td><b>남양</b></td>
   <td>세 자료의 좌표가 <b>전부 다릅니다.</b> 어느 것도 기재된 소재지와
   맞지 않습니다</td>
   <td><b>현재 좌표 미상</b> — 공간 계산에서 제외</td></tr>
  <tr><td><b>서산</b></td>
   <td>점조서에 <b>남양의 좌표가 통째로 복사</b>돼 있었습니다(52 km 차이)</td>
   <td>다른 자료 다섯 갈래로 대조해 참값 확정</td></tr>
 </table></div>
 <div class="quote reveal">반복측량은 <b>비슷한 곳</b>을 재는 것이 아니라
 <span class="hl">가능한 한 같은 지점</span>을 다시 재는 것입니다.</div>
 <p class="lead reveal">남양은 표석이 어디 있는지 아직 모릅니다.
 40년치 관측 기록이 있는데도 <b class="hl-r">지금 그 자리를 찾아갈 수가
 없습니다.</b> 좌표와 소재지를 정확히 남기는 일이 이만큼 중요합니다.</p>
</div></section>

<section data-t="두 번의 방문"><div class="wrap">
 <div class="reveal">
  <span class="tag" style="color:var(--red)">18</span>
  <h2>두 번의 방문이<br>전국 모델을 <span class="hl-r">26 % 망가뜨렸습니다</span></h2>
  <p class="lead">같은 표석을 몇 년 뒤 다시 가서 재면, 자연스러운 변화만큼만
  달라져야 합니다. 그런데 <b class="hl-r">34분</b> 씩 어긋나는 구간이
  있었습니다.</p>
  <p class="lead">원인을 찾아보니 <b>방문할 때마다 시준한 방위표지가
  달라져 있었습니다.</b> 표지를 바꾼 것 자체는 잘못이 아닙니다.
  바뀐 표지의 참방위각이 틀리면 그 방문의 편각이 통째로 그만큼 틀어집니다.</p>
 </div>
 <div class="reveal" style="margin-top:30px"><canvas class="fig" id="figAudit"
  width="1100" height="300" aria-label="두 방문 제거 전후의 모델 오차"></canvas></div>
 <div class="grid2 reveal" style="margin-top:22px">
  <div class="card warn"><h3>빼기 전</h3>
   <p>전국 편각 예측 오차 <b class="hl-r">0.7691°</b><br>
   편각 잔차 <b class="hl-r">35.0분</b></p></div>
  <div class="card good"><h3>부안·포천 2023 두 방문을 빼면</h3>
   <p>전국 편각 예측 오차 <b class="hl-g">0.5675°</b> (−26 %)<br>
   편각 잔차 <b class="hl-g">20.7분</b> (−41 %)</p></div>
 </div>
 <div class="quote reveal r">모델은 <b>잘못된 현장값도 정답이라고 믿습니다.</b>
 두 번의 방문이 전국의 값을 끌고 갔습니다.</div>
</div></section>

<section data-t="당신의 자리"><div class="wrap">
 <div class="reveal">
  <span class="tag">19</span>
  <h2>여러분이 고른 자리가<br><span class="hl">앞으로 수십 년</span>의 기준입니다</h2>
 </div>
 <div class="grid2 reveal" style="margin-top:30px">
  <div class="card"><h3 class="hl">① 그 자리에서 계속 잽니다</h3>
   <p>선점이 확정되면 표석을 세우고 <b>5년마다</b> 반복 관측을 합니다.
   지금 조용한 자리를 못 고르면 그 오차가 수십 년 따라갑니다.</p></div>
  <div class="card"><h3 class="hl">② 여러분의 자료가 기준이 됩니다</h3>
   <p>국내 판정 기준이 아직 없습니다. <b>이번에 {{N_SITES}}곳에서 모은 값으로
   그 기준을 만듭니다.</b> 여러분의 측정이 규정이 됩니다.</p></div>
 </div>
 <p class="lead reveal" style="margin-top:26px">그래서 「기록을 빠뜨리지 마세요」가
 형식적인 당부가 아닙니다. 시각 하나, 좌표 하나, 차수 하나가 빠지면
 그 줄은 <b class="hl-r">계산에서 통째로 빠집니다.</b> 나중에 사무실에서
 되살릴 방법이 없습니다.</p>
</div></section>

<section data-t="활용"><div class="wrap">
 <div class="reveal">
  <span class="tag">5부 — 이 자료가 어디에 쓰이는가 · 20</span>
  <h2>여러분의 기록은<br>여기에 쓰입니다</h2>
 </div>
 <div class="grid3 reveal" style="margin-top:30px">
  <div class="card good"><h3 class="hl-g">지금 바로</h3>
   <p>국가기본도의 자침편차 표기 · 측량과 방위 기준 ·
   국가 지자기 성과 고시 · 상시관측소 기선값 관리</p></div>
  <div class="card"><h3 class="hl">연구에</h3>
   <p>대한민국 지역 자기장 모델 · 장기 영년변화 분석 ·
   전지구 모델(IGRF)과 국내 실측값 대조 · 지각 자기이상 해석</p></div>
  <div class="card"><h3 style="color:var(--muted)">앞으로 넓혀 갈 수 있는 곳</h3>
   <p>항법 보조 · 우주환경 영향 분석 · 지하자원 탐사 기준 ·
   전력망·철도의 지자기 유도 영향 — <b>가능성이지 확정된 계획은
   아닙니다</b></p></div>
 </div>
</div></section>

<!-- ══ 에필로그 ══ -->
<section data-t="에필로그"><div class="wrap">
 <div class="reveal">
  <span class="tag">에필로그 · 21</span>
  <h2>한 점에서 전국으로</h2>
  <p class="lead">여러분이 현장에서 적은 각도 하나, 시각 하나, 좌표 하나가
  데이터베이스로, 전국 지도로, 대한민국 지자기 기준으로 이어집니다.</p>
 </div>
 <div class="quote reveal" style="font-size:clamp(22px,3.2vw,42px);margin:44px 0">
  대한민국 지자기 기준은<br><span class="hl">여러분의 기록 하나에서
  시작됩니다.</span></div>
 <div class="rules reveal">
  <div class="rule"><div><b>같은 지점에서 잽니다</b>
   <small>좌표와 소재지를 정확히 남겨야 다음 사람이 그 자리를 찾아갑니다</small></div></div>
  <div class="rule"><div><b>자기적 영향을 없앱니다</b>
   <small>몸에 지닌 쇠붙이, 차량, 전자기기 — 못 뺀 것은 적어 둡니다</small></div></div>
  <div class="rule"><div><b>시각을 시·분·초까지 적습니다</b>
   <small>F 를 잰 모든 줄에. 이것 하나가 빠지면 그 줄은 계산에서 빠집니다</small></div></div>
  <div class="rule"><div><b>이상한 값은 원인을 확인합니다</b>
   <small>값이 커서가 아니라 «절차가 깨져서» 다시 잽니다. 시작한 방향은
   끝까지 재 주세요</small></div></div>
  <div class="rule"><div><b>재현할 수 있게 남깁니다</b>
   <small>차수, 장비 원시파일명, 사진 파일명 — 나중에 되짚을 수 있어야
   자료입니다</small></div></div>
 </div>
</div></section>

</main>
<div class="foot">
 지자기 측량 현장교육 · {{GENERATED}}<br>
 수치는 <code>lmm_model.json</code> · <code>trial_survey_spec</code> ·
 <code>existing_pts.geojson</code> 에서 읽습니다 — 자료가 갱신되면 이 자료도
 따라갑니다.<br>
 Three.js r149 (MIT) 인라인 · 외부 요청 없음 · 오프라인 실행
</div>

<script>{{THREE}}</script>
<script>
const P = {{PAYLOAD}};

/* ── 스크롤 연출 ─────────────────────────────────── */
const secs = [...document.querySelectorAll("section")];
const nav = document.getElementById("nav");
secs.forEach((s, i) => {
  const b = document.createElement("i");
  b.title = s.dataset.t || ("장면 " + (i + 1));
  b.onclick = () => s.scrollIntoView({behavior: "smooth"});
  nav.appendChild(b);
});
const dots = [...nav.children];
const io = new IntersectionObserver(es => {
  es.forEach(e => { if (e.isIntersecting) e.target.classList.add("in"); });
}, {threshold: .15});
document.querySelectorAll(".reveal").forEach(el => io.observe(el));
const prog = document.getElementById("prog");
function onScroll() {
  const h = document.body.scrollHeight - innerHeight;
  prog.style.width = (h > 0 ? (scrollY / h) * 100 : 0) + "%";
  let cur = 0;
  secs.forEach((s, i) => { if (s.getBoundingClientRect().top < innerHeight * .5) cur = i; });
  dots.forEach((d, i) => d.classList.toggle("on", i === cur));
}
addEventListener("scroll", onScroll, {passive: true});
onScroll();

/* 발표 모드 — 방향키·스페이스로 장면 이동 */
let idx = 0;
addEventListener("keydown", e => {
  if (["ArrowDown", "PageDown", " ", "ArrowRight"].includes(e.key)) {
    e.preventDefault(); idx = Math.min(idx + 1, secs.length - 1);
    secs[idx].scrollIntoView({behavior: "smooth"});
  } else if (["ArrowUp", "PageUp", "ArrowLeft"].includes(e.key)) {
    e.preventDefault(); idx = Math.max(idx - 1, 0);
    secs[idx].scrollIntoView({behavior: "smooth"});
  }
});

/* ── 지도 도우미 ─────────────────────────────────── */
function mapProj(w, h, pad) {
  let x0 = 999, x1 = -999, y0 = 999, y1 = -999;
  P.boundary.forEach(r => r.forEach(([x, y]) => {
    if (x < x0) x0 = x; if (x > x1) x1 = x;
    if (y < y0) y0 = y; if (y > y1) y1 = y;
  }));
  const sx = (w - pad * 2) / (x1 - x0), sy = (h - pad * 2) / (y1 - y0);
  const s = Math.min(sx, sy);
  const ox = pad + ((w - pad * 2) - (x1 - x0) * s) / 2;
  const oy = pad + ((h - pad * 2) - (y1 - y0) * s) / 2;
  return (lon, lat) => [ox + (lon - x0) * s, h - oy - (lat - y0) * s];
}
function drawCoast(ctx, pj, color, lw) {
  ctx.strokeStyle = color; ctx.lineWidth = lw || 1; ctx.lineJoin = "round";
  P.boundary.forEach(r => {
    ctx.beginPath();
    r.forEach(([x, y], i) => {
      const [a, b] = pj(x, y);
      if (i === 0) ctx.moveTo(a, b); else ctx.lineTo(a, b);
    });
    ctx.stroke();
  });
}
function hidpi(cv) {
  const r = Math.min(devicePixelRatio || 1, 2);
  const w = cv.width, h = cv.height;
  cv.width = w * r; cv.height = h * r;
  cv.style.aspectRatio = w + " / " + h;
  const ctx = cv.getContext("2d");
  ctx.scale(r, r);
  return {ctx, w, h};
}

/* ── 나침반 ──────────────────────────────────────── */
(function () {
  const cv = document.getElementById("figCompass"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const cx = w / 2, cy = h / 2, R = Math.min(w, h) * .36;
  let ang = 0, t0 = null;
  function draw(ts) {
    if (t0 === null) t0 = ts;
    const k = Math.min((ts - t0) / 2200, 1);
    ang = -8.2 * (k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2);
    ctx.clearRect(0, 0, w, h);
    ctx.strokeStyle = "rgba(140,180,210,.22)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
    ctx.beginPath(); ctx.arc(cx, cy, R * .74, 0, 7); ctx.stroke();
    /* 진북 */
    ctx.strokeStyle = "#8ba3b8"; ctx.lineWidth = 1.6;
    ctx.setLineDash([5, 5]);
    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx, cy - R * 1.06); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "700 15px sans-serif";
    ctx.textAlign = "center"; ctx.fillText("진북 (지도의 북)", cx, cy - R * 1.14);
    /* 자북 바늘 */
    const rad = (ang - 90) * Math.PI / 180;
    ctx.save(); ctx.translate(cx, cy); ctx.rotate(rad + Math.PI / 2);
    ctx.fillStyle = "#3fd8d0";
    ctx.beginPath(); ctx.moveTo(0, -R * .96);
    ctx.lineTo(11, 0); ctx.lineTo(0, 16); ctx.lineTo(-11, 0); ctx.closePath();
    ctx.fill();
    ctx.fillStyle = "rgba(255,255,255,.22)";
    ctx.beginPath(); ctx.moveTo(0, R * .74);
    ctx.lineTo(8, 0); ctx.lineTo(0, -12); ctx.lineTo(-8, 0); ctx.closePath();
    ctx.fill();
    ctx.restore();
    /* 각 표시 */
    ctx.strokeStyle = "rgba(63,216,208,.55)"; ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(cx, cy, R * .42, -Math.PI / 2 + rad + Math.PI / 2 - Math.PI / 2, -Math.PI / 2);
    ctx.stroke();
    ctx.fillStyle = "#3fd8d0"; ctx.font = "800 30px 'Arial Narrow',sans-serif";
    ctx.textAlign = "left";
    ctx.fillText(ang.toFixed(1) + "°", cx - R * .62, cy - R * .28);
    ctx.font = "600 13px sans-serif"; ctx.fillStyle = "#8ba3b8";
    ctx.fillText("서편각", cx - R * .62, cy - R * .28 + 20);
    ctx.fillStyle = "#3fd8d0"; ctx.font = "700 15px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("자북 (나침반)", cx + Math.sin(rad + Math.PI / 2) * R * 1.18,
                 cy - Math.cos(rad + Math.PI / 2) * R * 1.18);
    ctx.beginPath(); ctx.arc(cx, cy, 5, 0, 7); ctx.fillStyle = "#eaf2f8"; ctx.fill();
    if (k < 1) requestAnimationFrame(draw);
  }
  new IntersectionObserver((es, ob) => {
    es.forEach(e => { if (e.isIntersecting) { requestAnimationFrame(draw); ob.disconnect(); } });
  }, {threshold: .3}).observe(cv);
})();

/* ── 50점 분포 ───────────────────────────────────── */
(function () {
  const cv = document.getElementById("figSites"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const pj = mapProj(w, h, 26);
  ctx.clearRect(0, 0, w, h);
  drawCoast(ctx, pj, "rgba(140,180,210,.35)", 1.1);
  const col = {"A": "#3fd8a0", "B": "#4a7fe8", "기존점": "#ff7048"};
  let n = 0;
  const tick = () => {
    for (let i = 0; i < 2 && n < P.sites.length; i++, n++) {
      const s = P.sites[n], [x, y] = pj(s.lon, s.lat);
      ctx.beginPath(); ctx.arc(x, y, 5.2, 0, 7);
      ctx.fillStyle = col[s.k] || "#fff"; ctx.fill();
      ctx.strokeStyle = "rgba(5,10,18,.8)"; ctx.lineWidth = 1.4; ctx.stroke();
    }
    if (n < P.sites.length) requestAnimationFrame(tick);
    else {
      ctx.font = "700 13px sans-serif"; ctx.textAlign = "left";
      const L = [["A 등급 " + P.n_grade_a, col.A],
                 ["B 등급 " + P.n_grade_b, col.B],
                 ["기존점 " + P.n_exist, col["기존점"]]];
      L.forEach((it, i) => {
        const y = h - 74 + i * 22;
        ctx.beginPath(); ctx.arc(30, y - 4, 5.2, 0, 7);
        ctx.fillStyle = it[1]; ctx.fill();
        ctx.fillStyle = "#8ba3b8"; ctx.fillText(it[0], 44, y);
      });
    }
  };
  new IntersectionObserver((es, ob) => {
    es.forEach(e => { if (e.isIntersecting) { tick(); ob.disconnect(); } });
  }, {threshold: .25}).observe(cv);
})();

/* ── 편각 영년변화 ───────────────────────────────── */
(function () {
  const cv = document.getElementById("figDrift"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const L = 58, R = w - 24, T = 26, B = h - 44;
  ctx.clearRect(0, 0, w, h);
  ctx.strokeStyle = "rgba(140,180,210,.22)"; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, B); ctx.lineTo(R, B); ctx.stroke();
  const y0 = 1990, y1 = 2030, d0 = -6.0, d1 = -9.6;
  const px = y => L + (y - y0) / (y1 - y0) * (R - L);
  const py = d => B - (d - d0) / (d1 - d0) * (B - T);
  ctx.font = "600 12px sans-serif"; ctx.fillStyle = "#8ba3b8";
  ctx.textAlign = "center";
  for (let y = 1990; y <= 2030; y += 10) {
    ctx.fillText(y, px(y), B + 20);
    ctx.strokeStyle = "rgba(140,180,210,.1)";
    ctx.beginPath(); ctx.moveTo(px(y), T); ctx.lineTo(px(y), B); ctx.stroke();
  }
  ctx.textAlign = "right";
  for (let d = -6; d >= -9.5; d -= 1) ctx.fillText(d.toFixed(0) + "°", L - 9, py(d) + 4);
  ctx.strokeStyle = "#3fd8d0"; ctx.lineWidth = 2.6;
  ctx.beginPath();
  for (let y = y0; y <= y1; y += .5) {
    const d = -6.0 - (y - y0) * 0.055;
    const X = px(y), Y = py(d);
    if (y === y0) ctx.moveTo(X, Y); else ctx.lineTo(X, Y);
  }
  ctx.stroke();
  const yn = 2026, dn = -6.0 - (yn - y0) * 0.055;
  ctx.beginPath(); ctx.arc(px(yn), py(dn), 6, 0, 7);
  ctx.fillStyle = "#ff7048"; ctx.fill();
  ctx.fillStyle = "#ff7048"; ctx.font = "700 13px sans-serif";
  ctx.textAlign = "left"; ctx.fillText("현재", px(yn) + 11, py(dn) - 6);
  ctx.fillStyle = "#8ba3b8"; ctx.font = "600 12.5px sans-serif";
  ctx.fillText("편각은 해마다 약 3분씩 서쪽으로 움직입니다", L + 6, T + 16);
})();

/* ── 벡터 분해 (WebGL) ───────────────────────────── */
(function () {
  const cv = document.getElementById("figVector"); if (!cv || !window.THREE) return;
  const W = cv.width, H = cv.height;
  const rn = new THREE.WebGLRenderer({canvas: cv, antialias: true, alpha: true});
  rn.setPixelRatio(Math.min(devicePixelRatio || 1, 2));
  rn.setSize(W, H, false);
  const sc = new THREE.Scene();
  const cam = new THREE.PerspectiveCamera(42, W / H, .1, 100);
  cam.position.set(3.1, 2.3, 3.6); cam.lookAt(0, .1, 0);
  sc.add(new THREE.AmbientLight(0xffffff, .85));
  const dl = new THREE.DirectionalLight(0xffffff, .55);
  dl.position.set(3, 5, 2); sc.add(dl);
  /* 지면 */
  const g = new THREE.GridHelper(4, 8, 0x2c4759, 0x1b2f3d);
  g.position.y = 0; sc.add(g);
  function arrow(dir, len, col, lw) {
    const a = new THREE.ArrowHelper(dir.clone().normalize(),
      new THREE.Vector3(0, 0, 0), len, col, len * .17, len * .09);
    a.line.material.linewidth = lw || 2;
    sc.add(a); return a;
  }
  const D = 8.2 * Math.PI / 180, I = 53 * Math.PI / 180;
  const Fv = new THREE.Vector3(
    Math.cos(I) * Math.cos(-D), -Math.sin(I), Math.cos(I) * Math.sin(-D));
  const items = [
    {v: Fv, l: 2.2, c: 0x3fd8d0},                                  /* F */
    {v: new THREE.Vector3(Fv.x, 0, Fv.z), l: 1.35, c: 0x4a7fe8},   /* H */
    {v: new THREE.Vector3(0, Fv.y, 0), l: 1.75, c: 0x9b7fe8},      /* Z */
    {v: new THREE.Vector3(1, 0, 0), l: 1.3, c: 0x55707f},          /* 북 */
  ];
  const arrows = items.map(it => { const a = arrow(it.v, .01, it.c); a.userData = it; return a; });
  let t = 0, run = false;
  function loop() {
    requestAnimationFrame(loop);
    if (!run) return;
    t = Math.min(t + .012, 1);
    arrows.forEach((a, i) => {
      const k = Math.max(0, Math.min((t - i * .16) / .5, 1));
      const e = k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
      a.setLength(Math.max(a.userData.l * e, .001),
                  a.userData.l * e * .17, a.userData.l * e * .09);
    });
    sc.rotation.y = Math.sin(t * Math.PI) * .18 + (t >= 1 ? performance.now() * .00012 : 0);
    rn.render(sc, cam);
  }
  loop();
  new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) run = true; }),
    {threshold: .25}).observe(cv);
})();

/* ── 측선 배치 ───────────────────────────────────── */
(function () {
  const cv = document.getElementById("figProfile"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const cx = w / 2, cy = h * .46, S = Math.min(w, h) * .036;
  ctx.clearRect(0, 0, w, h);
  ctx.strokeStyle = "rgba(140,180,210,.12)"; ctx.lineWidth = 1;
  [2, 5, 10].forEach(r => {
    ctx.beginPath(); ctx.arc(cx, cy, r * S, 0, 7); ctx.stroke();
    ctx.fillStyle = "#55707f"; ctx.font = "600 11px sans-serif";
    ctx.textAlign = "left"; ctx.fillText(r + " m", cx + r * S + 4, cy - 4);
  });
  const DIR = [[1, 0, "동"], [-1, 0, "서"], [0, 1, "남"], [0, -1, "북"]];
  let step = 0;
  const seq = [];
  DIR.forEach(([dx, dy, nm]) => {
    seq.push({p0: true, dx, dy, nm});
    P.offsets.forEach(o => seq.push({dx, dy, o, nm}));
    seq.push({p0: true, dx, dy, nm});
  });
  function tick() {
    const it = seq[step];
    if (!it) {
      ctx.fillStyle = "#3fd8d0"; ctx.font = "700 14px sans-serif";
      ctx.textAlign = "center";
      ctx.fillText("수평 " + P.n_read_h + "회", cx, h - 44);
      ctx.fillStyle = "#9b7fe8";
      ctx.fillText("+ 수직 " + P.n_read_v + "회 = 한 지점 " + P.n_read + "회", cx, h - 24);
      return;
    }
    if (it.p0) {
      ctx.beginPath(); ctx.arc(cx, cy, 8, 0, 7);
      ctx.fillStyle = "#ff7048"; ctx.fill();
    } else {
      const x = cx + it.dx * it.o * S, y = cy + it.dy * it.o * S;
      ctx.strokeStyle = "rgba(63,216,208,.3)"; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(x, y); ctx.stroke();
      ctx.beginPath(); ctx.arc(x, y, 3.6, 0, 7);
      ctx.fillStyle = "#3fd8d0"; ctx.fill();
      if (it.o === P.offsets[P.offsets.length - 1]) {
        ctx.fillStyle = "#8ba3b8"; ctx.font = "700 13px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(it.nm, x + it.dx * 22, y + it.dy * 22 + 4);
      }
    }
    step++;
    setTimeout(() => requestAnimationFrame(tick), 26);
  }
  ctx.fillStyle = "#ff7048"; ctx.font = "700 13px sans-serif";
  ctx.textAlign = "center"; ctx.fillText("P0 중심점", cx, cy - 20);
  new IntersectionObserver((es, ob) => {
    es.forEach(e => { if (e.isIntersecting) { tick(); ob.disconnect(); } });
  }, {threshold: .25}).observe(cv);
})();

/* ── 시간변화 보정 ───────────────────────────────── */
(function () {
  const cv = document.getElementById("figCorr"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const L = 52, R = w - 20, T = 26, B = h - 40;
  ctx.clearRect(0, 0, w, h);
  ctx.strokeStyle = "rgba(140,180,210,.22)";
  ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, B); ctx.lineTo(R, B); ctx.stroke();
  const px = m => L + m / 4 * (R - L);
  const py = v => B - (v - 96) / 12 * (B - T);
  ctx.font = "600 12px sans-serif"; ctx.fillStyle = "#8ba3b8";
  ctx.textAlign = "center";
  ["10:00", "10:01", "10:02", "10:03", "10:04"].forEach((s, i) =>
    ctx.fillText(s, px(i), B + 19));
  /* P0 선형보간 */
  ctx.strokeStyle = "#ff7048"; ctx.lineWidth = 2; ctx.setLineDash([6, 5]);
  ctx.beginPath(); ctx.moveTo(px(0), py(100)); ctx.lineTo(px(4), py(105)); ctx.stroke();
  ctx.setLineDash([]);
  [[0, 100], [4, 105]].forEach(([m, v]) => {
    ctx.beginPath(); ctx.arc(px(m), py(v), 6, 0, 7);
    ctx.fillStyle = "#ff7048"; ctx.fill();
  });
  /* 추정값 */
  ctx.beginPath(); ctx.arc(px(2), py(102.5), 5, 0, 7);
  ctx.fillStyle = "#ffd7c8"; ctx.fill();
  /* 측정값 */
  ctx.beginPath(); ctx.arc(px(2), py(106), 7, 0, 7);
  ctx.fillStyle = "#3fd8d0"; ctx.fill();
  /* ΔF */
  ctx.strokeStyle = "#3fd8a0"; ctx.lineWidth = 2.4;
  ctx.beginPath(); ctx.moveTo(px(2), py(106)); ctx.lineTo(px(2), py(102.5)); ctx.stroke();
  ctx.fillStyle = "#3fd8a0"; ctx.font = "800 15px sans-serif";
  ctx.textAlign = "left"; ctx.fillText("ΔF = 3.5 nT", px(2) + 12, py(104.2));
  ctx.fillStyle = "#ff7048"; ctx.font = "700 12px sans-serif";
  ctx.fillText("P0 앞뒤를 이은 선", px(0) + 8, py(100) - 12);
  ctx.fillStyle = "#3fd8d0";
  ctx.fillText("동쪽 2 m 측정값", px(2) + 12, py(106) - 10);
})();

/* ── 자기청정 ────────────────────────────────────── */
(function () {
  const cv = document.getElementById("figClean"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  ctx.clearRect(0, 0, w, h);
  const cx = w * .38, cy = h * .5;
  /* 사람 실루엣 */
  ctx.strokeStyle = "rgba(140,180,210,.5)"; ctx.lineWidth = 2.2;
  ctx.beginPath(); ctx.arc(cx, cy - 118, 26, 0, 7); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx, cy - 92); ctx.lineTo(cx, cy + 24); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx - 44, cy - 46); ctx.lineTo(cx + 44, cy - 46); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx, cy + 24); ctx.lineTo(cx - 30, cy + 128); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx, cy + 24); ctx.lineTo(cx + 30, cy + 128); ctx.stroke();
  const ITEM = [
    [cx + 46, cy - 44, "손목시계"], [cx - 46, cy - 30, "휴대전화"],
    [cx + 16, cy - 62, "금속 펜"], [cx, cy + 22, "혁대 버클"],
    [cx - 30, cy + 126, "안전화 철심"], [cx + 6, cy - 122, "안경 부품"],
    [cx - 20, cy - 20, "열쇠"],
  ];
  let n = 0;
  function tick() {
    if (n >= ITEM.length) return;
    const [x, y, nm] = ITEM[n];
    const gr = ctx.createRadialGradient(x, y, 0, x, y, 34);
    gr.addColorStop(0, "rgba(232,80,63,.55)");
    gr.addColorStop(1, "rgba(232,80,63,0)");
    ctx.fillStyle = gr;
    ctx.beginPath(); ctx.arc(x, y, 34, 0, 7); ctx.fill();
    ctx.beginPath(); ctx.arc(x, y, 4.5, 0, 7);
    ctx.fillStyle = "#e8503f"; ctx.fill();
    ctx.strokeStyle = "rgba(232,80,63,.5)"; ctx.lineWidth = 1;
    const tx = w * .74;
    ctx.beginPath(); ctx.moveTo(x + 5, y); ctx.lineTo(tx - 8, y); ctx.stroke();
    ctx.fillStyle = "#f0b5ad"; ctx.font = "600 13.5px sans-serif";
    ctx.textAlign = "left"; ctx.fillText(nm, tx, y + 4);
    n++;
    setTimeout(() => requestAnimationFrame(tick), 340);
  }
  new IntersectionObserver((es, ob) => {
    es.forEach(e => { if (e.isIntersecting) { tick(); ob.disconnect(); } });
  }, {threshold: .3}).observe(cv);
})();

/* ── 감사 전후 ───────────────────────────────────── */
(function () {
  const cv = document.getElementById("figAudit"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  ctx.clearRect(0, 0, w, h);
  const rows = [
    ["전국 편각 예측 오차", 0.7691, 0.5675, "°", 1.0],
    ["편각 잔차", 35.0, 20.73, "분", 40],
  ];
  const L = 200, R = w - 130;
  rows.forEach((r, i) => {
    const y = 82 + i * 108;
    ctx.fillStyle = "#8ba3b8"; ctx.font = "700 14px sans-serif";
    ctx.textAlign = "right"; ctx.fillText(r[0], L - 18, y + 5);
    const bw = v => (v / r[4]) * (R - L);
    ctx.fillStyle = "rgba(232,80,63,.75)";
    ctx.fillRect(L, y - 24, bw(r[1]), 20);
    ctx.fillStyle = "rgba(63,216,160,.8)";
    ctx.fillRect(L, y + 6, bw(r[2]), 20);
    ctx.textAlign = "left"; ctx.font = "800 14px ui-monospace,monospace";
    ctx.fillStyle = "#f0b5ad";
    ctx.fillText(r[1] + " " + r[3], L + bw(r[1]) + 10, y - 9);
    ctx.fillStyle = "#8ce8c4";
    ctx.fillText(r[2] + " " + r[3], L + bw(r[2]) + 10, y + 21);
  });
  ctx.fillStyle = "#8ba3b8"; ctx.font = "600 12.5px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText("■ 두 방문 포함", L, h - 22);
  ctx.fillStyle = "#8ce8c4";
  ctx.fillText("■ 두 방문 제외", L + 140, h - 22);
})();

/* ── 배경 자기력선 (WebGL) ───────────────────────── */
(function () {
  const host = document.getElementById("bg");
  if (!window.THREE || matchMedia("(prefers-reduced-motion:reduce)").matches) return;
  if (innerWidth < 860) return;
  const rn = new THREE.WebGLRenderer({antialias: true, alpha: true});
  rn.setPixelRatio(Math.min(devicePixelRatio || 1, 1.6));
  rn.setSize(innerWidth, innerHeight);
  host.appendChild(rn.domElement);
  const sc = new THREE.Scene();
  const cam = new THREE.PerspectiveCamera(50, innerWidth / innerHeight, .1, 200);
  cam.position.set(0, 0, 46);
  const grp = new THREE.Group(); sc.add(grp);
  const mat = new THREE.LineBasicMaterial({
    color: 0x1e5570, transparent: true, opacity: .30});
  for (let i = 0; i < 18; i++) {
    const pts = [], a = (i / 18) * Math.PI * 2, r0 = 5 + (i % 4) * 2.2;
    for (let t = -1; t <= 1; t += .025) {
      const r = r0 * (1 - t * t * .72);
      pts.push(new THREE.Vector3(
        Math.cos(a) * r, t * 17, Math.sin(a) * r));
    }
    grp.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat));
  }
  let raf;
  function loop() {
    raf = requestAnimationFrame(loop);
    grp.rotation.y += .0011;
    grp.rotation.x = Math.sin(performance.now() * .00012) * .1;
    const s = 1 - Math.min(scrollY / (document.body.scrollHeight || 1), .35);
    grp.scale.setScalar(s);
    rn.render(sc, cam);
  }
  loop();
  addEventListener("resize", () => {
    cam.aspect = innerWidth / innerHeight; cam.updateProjectionMatrix();
    rn.setSize(innerWidth, innerHeight);
  });
})();
</script>
</body></html>
"""


def pilot_table(payload):
    rows = []
    for p in payload["pilot"]:
        cells = "".join(
            f'<td><b class="{"hl-r" if abs(v) >= 50 else ""}">{v:+.2f}</b></td>'
            for v in (p["dirs"].get(d) for d in ("동", "서", "남", "북"))
        )
        rows.append(f"<tr><td><b>{p['name']}</b></td>{cells}</tr>")
    return (
        '<div class="tblwrap"><table>'
        '<tr><th>지점</th><th>동</th><th>서</th><th>남</th><th>북</th></tr>'
        + "".join(rows)
        + '</table></div>'
        '<p style="color:var(--muted);font-size:13.5px;margin-top:10px">'
        '단위 nT/m · 중심점에서 1 m 떨어진 지점의 값</p>'
    )


def render(payload):
    html = TEMPLATE
    tokens = {
        "GENERATED": payload["generated"],
        "N_SITES": payload["n_sites"],
        "N_NEW": payload["n_new"],
        "N_EXIST": payload["n_exist"],
        "N_NETWORK": payload["n_network"],
        "TARGET_NETWORK": payload["target_network"],
        "LOO_D_MIN": payload["loo_d_min"],
        "LOO_F": payload["loo_f"],
        "OFFSETS_TXT": ", ".join(str(x) for x in payload["offsets"]),
        "H_MIN": payload["heights"][0],
        "H_MAX": payload["heights"][-1],
        "H_STEP": payload["heights"][1] - payload["heights"][0],
        "N_READ": payload["n_read"],
        "N_READ_H": payload["n_read_h"],
        "N_READ_V": payload["n_read_v"],
        "N_OFF": len(payload["offsets"]),
        "IAGA_RADIUS": f"{payload['iaga_radius']:.0f}",
        "IAGA_RANGE": f"{payload['iaga_range']:.0f}",
        "EURO_GRAD": f"{payload['euro_grad']:.0f}",
        "WORST_VAL": payload["worst_val"],
        "WORST_RATIO": payload["worst_ratio"],
        "PILOT_TABLE": pilot_table(payload),
    }
    for k, v in tokens.items():
        html = html.replace("{{" + k + "}}", str(v))
    html = html.replace("{{PAYLOAD}}", json.dumps(payload, ensure_ascii=False))
    html = html.replace("{{THREE}}", read_three_runtime())
    return html


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    sys.path.insert(0, str(ROOT))
    payload = build_payload()
    html = render(payload)
    OUT.write_text(html, encoding="utf-8")
    size = OUT.stat().st_size / 1e6
    print(f"대상 {payload['n_sites']}점 "
          f"(신규 {payload['n_new']} · 기존 {payload['n_exist']}) · "
          f"관측망 {payload['n_network']}점")
    print(f"측선 {payload['n_read']}회/점 · 최악 시범값 {payload['worst_val']} nT/m")
    print(f"[저장] {OUT}  ({size:.1f} MB)")
    return OUT


if __name__ == "__main__":
    main()
