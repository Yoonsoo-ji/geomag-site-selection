#!/usr/bin/env python3
"""LMM 시네마틱 브리핑의 A/B/C 디자인 비교본을 생성한다.

현재 배포본 ``docs/lmm_cinematic.html`` 은 건드리지 않는다. 세 비교본은 같은
데이터와 같은 교정 문구를 사용하며 디자인만 달라야 한다.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import create_lmm_cinematic as base


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "docs" / "design_previews"


COMMON_CSS = r"""
/* ── A/B/C 공통: 한글 조판·720p 압축·지도 색 토큰 ─────────────── */
:root{
  --map-bg:#09141d;--map-land:#0d1b26;--map-stroke:#66869a;
  --inset-bg:#071018;--inset-line:#355267;--inset-text:#dce9f1;
}
html,body{max-width:100%;overflow-x:hidden}
body{word-break:keep-all;overflow-wrap:normal;line-break:strict;hyphens:none}
h1,h2,h3,.loop-note,.warning,.decision{text-wrap:balance;word-break:keep-all;overflow-wrap:normal;line-break:strict}
p,li,.lead,.axis-card p,.case p,.usecase p,.feedback-card p,.phase li,
.flow-node span,.stat small,.layer small{text-wrap:pretty;word-break:keep-all;overflow-wrap:normal;line-break:strict}
.no-break{white-space:nowrap}
h1{font-size:clamp(44px,6.4vw,96px);line-height:.96;letter-spacing:-.055em}
h2{font-size:clamp(34px,4.35vw,62px);line-height:1.08;letter-spacing:-.038em}
h3{font-size:clamp(18px,1.55vw,24px);line-height:1.32}
.lead{font-size:clamp(17px,1.55vw,23px);line-height:1.58}
.grid2 h2{font-size:clamp(34px,3.65vw,52px);line-height:1.1}
.grid2 .lead{font-size:clamp(17px,1.42vw,21px)}
.axis-pair{margin-top:20px}.axis-card{padding:18px}.axis-card .symbol{font-size:30px}
.axis-card h3{margin-bottom:10px}.axis-card p{font-size:14px;line-height:1.58;margin-bottom:0}
.loop-note{font-size:14px;line-height:1.55}.warning{font-size:13.5px;line-height:1.55}
.map-note{font-size:11.5px}.stat small{font-size:13px}.algo-step b{font-size:15px}
.algo-step p{font-size:12.5px}.feedback-card p{font-size:13px}.case p{font-size:14px}
.theme-mark{font:800 10px/1.2 "IBM Plex Mono",Consolas,monospace;letter-spacing:.16em;
  color:var(--orange);white-space:nowrap}
.adaptive-scene .grid2{grid-template-columns:minmax(430px,.9fr) minmax(540px,1.1fr)}
.adaptive-scene h2{font-size:clamp(32px,3.3vw,47px);margin-bottom:13px}.adaptive-scene .lead{margin-bottom:0}
.adaptive-scene .axis-pair{gap:12px;margin-top:16px}.adaptive-scene .axis-card{padding:16px}
.adaptive-scene .axis-card p{font-size:13.2px;line-height:1.52}
.adaptive-scene .loop-note{margin-top:9px!important;padding:11px 15px}.adaptive-scene .warning{margin-top:8px!important;padding:11px 15px}
@media(min-width:781px) and (max-height:820px){
  .scene{padding-top:64px;padding-bottom:12px}.eyebrow{margin-bottom:10px}
  h2{margin-bottom:17px}.hero .kicker{margin-bottom:24px}.hero-foot{margin-top:30px}
  .map-shell{min-height:510px}.map-shell canvas{height:510px}
  .flow-track{margin-top:32px}.algo{margin-top:22px}.algo-step{min-height:156px}
  .feedback{margin-top:28px}.feedback-card{min-height:148px}.roadmap{margin-top:22px}
  .phase{min-height:238px}.case{min-height:300px}.usecase{min-height:190px}
  .formula{margin:18px 0;padding:18px 22px}.ledger{margin-top:18px}
  .perf-panel{padding:18px}.bar-row{margin:10px 0}.gate{margin-top:16px;padding:15px}
  .adaptive-scene .map-shell,.adaptive-scene .map-shell canvas{height:500px;min-height:500px}
  .closing h2{font-size:clamp(44px,5.2vw,70px)}
  .source-grid{gap:28px}.sources-list{gap:6px}.source-item{padding:10px 13px;font-size:12px;line-height:1.42}
  .source-item a{font-size:10px}.footnote{font-size:12px;line-height:1.5}
}
@media(max-width:1050px){.adaptive-scene .grid2{grid-template-columns:1fr}.adaptive-scene .map-shell{order:2}}
@media(max-width:780px){
  h1{font-size:clamp(38px,10vw,58px)}h2,.grid2 h2{font-size:clamp(31px,8.4vw,46px)}
  .lead,.grid2 .lead{font-size:17px}.scene{padding-inline:22px}
  .adaptive-scene .map-shell{order:2}.axis-pair{grid-template-columns:1fr}
  .theme-mark{display:none}
}
@media(max-width:520px){.flow-track,.stat-grid{grid-template-columns:1fr}.flow-node{text-align:left;display:grid;grid-template-columns:56px 1fr;column-gap:14px;align-items:center}.flow-node .dot{grid-row:1/3;margin:0}.flow-node b{align-self:end;margin:0 0 3px}.flow-node span{align-self:start}.stat small{font-size:13.5px}}
"""


THEMES = {
    "A": {
        "name": "MAGNETIC OBSERVATORY",
        "label": "A · Magnetic Observatory",
        "css": r"""
:root{--ink:#f1fbff;--muted:#9db4c2;--dim:#637d8e;--bg:#020812;--panel:#071521;
  --line:#17394a;--cyan:#43f2e8;--blue:#4e7dff;--orange:#ff6840;--gold:#ffd166;
  --map-bg:#020b14;--map-land:rgba(10,32,43,.46);--map-stroke:#69bac0;
  --inset-bg:#04111b;--inset-line:#2d6573;--inset-text:#e8fbff}
body{background:radial-gradient(circle at 78% 18%,rgba(54,115,255,.15),transparent 28%),
  radial-gradient(circle at 18% 84%,rgba(67,242,232,.08),transparent 30%),var(--bg)}
body::before{opacity:.24;background-image:linear-gradient(rgba(67,242,232,.025) 1px,transparent 1px),
  linear-gradient(90deg,rgba(67,242,232,.025) 1px,transparent 1px);background-size:34px 34px;mix-blend-mode:screen}
.topbar{background:rgba(2,8,18,.84);border-color:rgba(67,242,232,.16)}
.progress{background:#0c2431}.progress span{box-shadow:0 0 16px var(--orange)}
.scene::before{content:"";position:absolute;inset:0;pointer-events:none;background:
  linear-gradient(112deg,transparent 0 72%,rgba(67,242,232,.025) 72% 72.15%,transparent 72.15%)}
.scene::after{color:rgba(67,242,232,.028);-webkit-text-stroke:1px rgba(67,242,232,.06)}
.eyebrow{color:var(--cyan);text-shadow:0 0 18px rgba(67,242,232,.32)}
.tag,.axis-card,.stat,.layer,.algo-step,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item{
  background:linear-gradient(145deg,rgba(12,31,43,.9),rgba(4,14,24,.88));border-color:rgba(85,185,198,.23)}
.axis-card,.perf-panel,.map-shell{clip-path:polygon(0 0,calc(100% - 17px) 0,100% 17px,100% 100%,0 100%)}
.axis-card{box-shadow:inset 0 2px 0 rgba(67,242,232,.16)}
.axis-card.region{box-shadow:inset 0 2px 0 rgba(255,104,64,.35)}
.map-shell{border-color:rgba(67,242,232,.3);box-shadow:0 24px 80px rgba(0,0,0,.48),0 0 0 1px rgba(67,242,232,.04)}
.chip,.controls button{background:rgba(2,13,23,.86);border-color:rgba(67,242,232,.24)}
.hero h1 .ghost{-webkit-text-stroke-color:rgba(67,242,232,.45)}
.formula{border-color:var(--cyan);background:linear-gradient(90deg,rgba(67,242,232,.08),transparent)}
.gate{background:linear-gradient(90deg,rgba(255,104,64,.1),rgba(255,104,64,.02))}
""",
    },
    "B": {
        "name": "CARTOGRAPHIC EDITORIAL",
        "label": "B · Cartographic Editorial",
        "css": r"""
:root{--ink:#17242b;--muted:#59676c;--dim:#7a827f;--bg:#f2eee4;--panel:#fbf8f0;
  --line:#c9c1b2;--cyan:#0d7477;--blue:#2457a2;--orange:#c85234;--gold:#a87923;
  --good:#34775f;--danger:#b8422e;--map-bg:#dfe9e7;--map-land:#f7f2e7;
  --map-stroke:#245463;--map-land:rgba(247,242,231,.44);--inset-bg:#efe9dc;--inset-line:#79919a;--inset-text:#17242b}
html{background:var(--bg)}body{background-color:var(--bg);background-image:
  radial-gradient(rgba(23,36,43,.05) .55px,transparent .75px);background-size:5px 5px}
body::before{opacity:.55;background:linear-gradient(90deg,transparent 0 7.5%,rgba(36,87,162,.055) 7.5% 7.65%,transparent 7.65%);
  mix-blend-mode:multiply}
.topbar{background:rgba(242,238,228,.92);border-color:#bdb4a4;color:var(--ink);backdrop-filter:blur(9px)}
.progress{background:#d6cfc1}.counter,.brand{color:#586468}.scene-nav button{background:#7f8d90}.controls button,.chip{
  background:rgba(247,243,233,.94);color:var(--ink);border-color:#a9a091;backdrop-filter:none}
.scene::after{color:rgba(36,87,162,.045);font-family:"Bodoni 72","Times New Roman",serif}
.scene:nth-child(even){background:rgba(255,255,255,.2)}
.hero{background:linear-gradient(100deg,rgba(242,238,228,.98),rgba(229,235,230,.92))}
.hero h1,.closing h2,h2{font-family:"Noto Serif KR","Nanum Myeongjo","Batang",serif;font-weight:800;letter-spacing:-.045em}
.hero h1 .ghost{color:transparent;-webkit-text-stroke:1px rgba(36,87,162,.58)}
.eyebrow{color:var(--blue)}.eyebrow::before{height:1px;width:58px}.tag{background:#f8f4eb;border-color:#c8c0b1}
.lead{color:#45565c}.lead strong{color:#142329}
.axis-card b{color:var(--ink)}
.axis-card,.stat,.layer,.algo-step,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item{
  background:rgba(250,247,239,.86);border-color:#c9c1b2;box-shadow:4px 4px 0 rgba(35,64,74,.065)}
.axis-card{border-top:4px solid var(--blue)}.axis-card.region{border-top-color:var(--orange)}
.map-shell{background:var(--map-bg);border:1px solid #91a2a2;box-shadow:12px 14px 0 rgba(45,74,80,.09)}
.map-legend,.map-note{background:rgba(248,244,234,.9);border-color:#aba394;color:#495b60}
.loop-note{background:rgba(255,255,255,.22);border-color:#b8aea0;color:#34484f}.loop-note b{color:var(--orange)}
.warning{background:rgba(200,82,52,.07);color:#7c3728}.formula{background:rgba(36,87,162,.05);border-color:var(--blue)}
.layer{transform:none;box-shadow:4px 4px 0 rgba(35,64,74,.065)}.bar-track{background:#d9d3c7}
.gate{background:rgba(200,82,52,.06);border-color:rgba(200,82,52,.42)}.gate p{color:#713c31}
.case p,.usecase p,.phase ul,.feedback-card p,.source-item span,.footnote{color:#59676c}
.decision-strip{background:#bdb5a8;border-color:#bdb5a8}.decision{background:#f7f3ea}
.closing{background:radial-gradient(circle at 50% 60%,rgba(13,116,119,.1),transparent 38%),var(--bg)}
.primary-link{background:#f7f3ea}.source-item a{color:var(--blue)}
""",
    },
    "C": {
        "name": "AURORA SPECTRAL",
        "label": "C · Aurora Spectral",
        "css": r"""
:root{--ink:#f8f5ff;--muted:#b5acd0;--dim:#786f98;--bg:#050515;--panel:#111026;
  --line:rgba(189,171,255,.2);--cyan:#5af7f2;--blue:#7867ff;--orange:#ff5c9c;--gold:#ffd36a;
  --map-bg:#09081d;--map-land:rgba(27,26,53,.46);--map-stroke:#82e8e4;
  --inset-bg:#0d0b25;--inset-line:#655a9a;--inset-text:#f8f5ff}
body{background:radial-gradient(ellipse at 15% 20%,rgba(90,247,242,.12),transparent 31%),
  radial-gradient(ellipse at 82% 18%,rgba(120,103,255,.2),transparent 34%),
  radial-gradient(ellipse at 60% 84%,rgba(255,92,156,.12),transparent 32%),var(--bg)}
body::before{opacity:.28;background-image:radial-gradient(rgba(255,255,255,.45) .55px,transparent .8px);
  background-size:7px 7px;mix-blend-mode:screen}
.topbar{background:rgba(7,6,24,.62);border-color:rgba(184,164,255,.18);backdrop-filter:blur(22px) saturate(145%)}
.progress{background:rgba(255,255,255,.08)}.progress span{background:linear-gradient(90deg,var(--cyan),var(--blue),var(--orange));box-shadow:0 0 18px var(--blue)}
.scene::before{content:"";position:absolute;width:52vw;height:52vw;right:-20vw;top:-25vw;border-radius:50%;
  background:conic-gradient(from 80deg,transparent,rgba(90,247,242,.08),rgba(120,103,255,.13),rgba(255,92,156,.08),transparent);
  filter:blur(32px);animation:auroraSpin 22s linear infinite;pointer-events:none}
@keyframes auroraSpin{to{transform:rotate(360deg)}}
.scene::after{color:rgba(185,166,255,.035);text-shadow:0 0 50px rgba(120,103,255,.2)}
.eyebrow{color:var(--cyan)}h2 .accent,.accent{color:var(--orange)}
.hero h1 .ghost{background:linear-gradient(90deg,var(--cyan),#9e8cff,var(--orange));-webkit-background-clip:text;
  color:transparent;-webkit-text-stroke:0}
.tag,.axis-card,.stat,.layer,.algo-step,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item,
.map-legend,.map-note{background:linear-gradient(145deg,rgba(31,28,65,.7),rgba(11,10,34,.58));
  border-color:rgba(203,188,255,.2);backdrop-filter:blur(18px) saturate(135%);box-shadow:inset 0 1px rgba(255,255,255,.06),0 22px 55px rgba(0,0,0,.18)}
.axis-card,.stat,.perf-panel,.feedback-card,.case,.usecase,.phase,.source-item,.warning,.loop-note,.gate{border-radius:18px}
.map-shell{border-radius:26px;border-color:rgba(90,247,242,.25);box-shadow:0 35px 100px rgba(0,0,0,.48),0 0 80px rgba(120,103,255,.13)}
.map-shell canvas{border-radius:26px}.chip,.controls button{border-radius:999px;background:rgba(20,17,52,.72)}
.axis-card.point{box-shadow:inset 0 1px rgba(255,255,255,.08),0 0 48px rgba(90,247,242,.06)}
.axis-card.region{box-shadow:inset 0 1px rgba(255,255,255,.08),0 0 48px rgba(255,92,156,.08)}
.formula{border:1px solid rgba(185,166,255,.24);border-left:3px solid var(--cyan);border-radius:18px;
  background:linear-gradient(100deg,rgba(90,247,242,.08),rgba(120,103,255,.08),rgba(255,92,156,.06))}
.bar-track{border-radius:8px;overflow:hidden}.bar-fill{background:linear-gradient(90deg,var(--blue),var(--cyan))}
.bar-fill.orange{background:linear-gradient(90deg,var(--orange),var(--gold))}.warning{background:rgba(255,92,156,.09);color:#f2bfd7}
.loop-note{background:rgba(120,103,255,.08);border-color:rgba(185,166,255,.28)}
.closing{background:radial-gradient(circle at 50% 55%,rgba(90,247,242,.11),transparent 30%),
  radial-gradient(circle at 62% 45%,rgba(255,92,156,.1),transparent 28%),var(--bg)}
@media(prefers-reduced-motion:reduce){.scene::before{animation:none}}
""",
    },
}


SCENE_04 = r'''<section class="scene adaptive-scene" data-no="04" data-title="적응형 선점"><div class="inner grid2">
  <div class="map-shell reveal"><canvas id="densityMap" aria-label="자기이상 공간복잡도 대비 상대 부족도 지도"></canvas><div class="map-legend"><div class="legend-row"><i class="swatch" style="background:#ff7048"></i>보강 1순위</div><div class="legend-row"><i class="swatch" style="background:#f0bc54"></i>평균보다 부족</div><div class="legend-row"><i class="swatch" style="background:#427df4"></i>평균보다 양호</div><div class="legend-row"><i class="swatch" style="background:#244257"></i>충분</div></div><div class="map-note"><span data-bind="land_cells">1,091</span> 국토셀 · 보강 1순위 <span data-bind="critical_cells">76</span> · σ 대체 <span data-bind="sigma_missing">57</span><br>울릉도·독도 위치 인셋 포함</div></div>
  <div>
    <div class="eyebrow reveal">TWO SCALES, ONE RULE</div><h2 class="reveal">권역 밀도는 공간복잡도로,<br>후보점은 국소 구배로 설계합니다.</h2>
    <p class="lead reveal">항공자력 자료를 <strong>점 대표성</strong>과 <strong>권역 밀도</strong>라는 서로 다른 질문에 사용합니다.</p>
    <div class="axis-pair reveal" style="--delay:120ms">
      <article class="axis-card point"><div class="symbol">|∇ΔT|</div><h3>① 후보점 대표성 · <b>저구배 우대</b></h3><p>KIGAM 1.5분(약 2.8 km) 격자를 중앙차분한 수평 공간구배(nT/km)를 후보 178점의 s5에 실제 반영했습니다. <b>1 km 해상도 분석은 아니며</b>, 최종 확정은 현장 정밀 자력측량으로 합니다.</p></article>
      <article class="axis-card region"><div class="symbol">σ<sub>25 km</sub></div><h3>② 권역 밀도 · <b>고복잡도 보강</b></h3><p>반경 25 km 자기이상에서 1차 지역추세를 제거한 표준편차 σ를 <b>자기이상 공간복잡도</b>로 사용합니다. σ가 클수록 권장 간격을 줄입니다.</p></article>
    </div>
    <div class="loop-note reveal" style="--delay:190ms"><b>고복잡도 권역을 보강하되, 그 안에서는 저구배 후보점을 고릅니다.</b></div><div class="warning reveal" style="--delay:230ms">Neyman 배분과 s5 곡선은 결측·소표본 한계가 있는 <strong>잠정 설계안</strong>입니다.</div>
  </div>
</div></section>'''


MAP_HELPERS = r'''
const css=(name,fallback)=>getComputedStyle(document.documentElement).getPropertyValue(name).trim()||fallback;
function drawEastSeaInset(ctx,W,H){
  const iw=Math.min(184,Math.max(142,W*.34)),ih=112,x=W-iw-14,y=14;
  const b={w:130.62,e:132.08,s:37.05,n:37.68},pad=10;
  const scale=Math.min((iw-pad*2)/(b.e-b.w),(ih-pad*2)/(b.n-b.s));
  const ox=x+(iw-(b.e-b.w)*scale)/2,oy=y+(ih-(b.n-b.s)*scale)/2;
  const p=(lon,lat)=>[ox+(lon-b.w)*scale,oy+(b.n-lat)*scale];
  ctx.save();ctx.fillStyle=css('--inset-bg','#071018');ctx.fillRect(x,y,iw,ih);
  ctx.strokeStyle=css('--inset-line','#355267');ctx.lineWidth=1;ctx.strokeRect(x+.5,y+.5,iw-1,ih-1);
  ctx.beginPath();ctx.rect(x+1,y+1,iw-2,ih-2);ctx.clip();
  drawGeometry(ctx,DATA.boundary,p);ctx.fillStyle=css('--map-land','#0d1b26');ctx.fill('evenodd');
  ctx.strokeStyle=css('--map-stroke','#66869a');ctx.lineWidth=.8;ctx.stroke();
  DATA.east_sea.forEach(site=>{const [px,py]=p(site.lon,site.lat);ctx.beginPath();ctx.arc(px,py,site.name==='독도'?3.2:2.6,0,Math.PI*2);
    ctx.fillStyle=site.name==='독도'?css('--orange','#ff7048'):css('--cyan','#51d2d7');ctx.fill();
    if(site.name==='독도'){ctx.beginPath();ctx.arc(px,py,7,0,Math.PI*2);ctx.strokeStyle=css('--orange','#ff7048');ctx.lineWidth=.8;ctx.stroke()}
    ctx.font='700 10px "Noto Sans KR","Malgun Gothic",sans-serif';ctx.fillStyle=css('--inset-text','#dce9f1');
    ctx.textAlign=site.name==='독도'?'right':'left';ctx.textBaseline='middle';ctx.fillText(site.name,px+(site.name==='독도'?-7:7),py-1)});
  ctx.restore();ctx.font='700 8px Consolas,monospace';ctx.fillStyle=css('--inset-text','#dce9f1');ctx.textAlign='left';ctx.textBaseline='top';ctx.fillText('EAST SEA ISLANDS',x+8,y+7);
}
'''


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {text.count(old)}")
    return text.replace(old, new, 1)


def corrected_base() -> str:
    text = base.HTML
    text, n = re.subn(
        r'\n/\* BEGIN AURORA SPECTRAL \*/.*?/\* END AURORA SPECTRAL \*/\n',
        "\n",
        text,
        count=1,
        flags=re.S,
    )
    if n != 1:
        raise RuntimeError(f"production theme block removal failed: {n}")
    text = text.replace('href="index.html"', 'href="../index.html"')
    text = text.replace('href="lmm.html"', 'href="../lmm.html"')
    return text


def render_variant(theme_key: str, payload: dict) -> Path:
    theme = THEMES[theme_key]
    text = corrected_base()
    css = COMMON_CSS + "\n" + theme["css"]
    text = replace_once(text, "</style>", css + "\n</style>", f"theme {theme_key} CSS")
    text = replace_once(text, "<body>", f'<body data-theme="{theme_key}">', f"theme {theme_key} body")
    text = replace_once(
        text,
        '<div class="brand"><b>KOREA</b> LMM / AURORA SPECTRAL</div>',
        f'<div class="brand"><b>KOREA</b> LMM / {theme["name"]}</div><div class="theme-mark">{theme["label"]}</div>',
        f"theme {theme_key} brand",
    )
    text = text.replace(
        "<title>한반도 LMM — 관측망에서 국가 자기장 기준으로</title>",
        f"<title>{theme['label']} — 한반도 LMM</title>",
    )
    if theme_key == "B":
        text = text.replace('<meta name="color-scheme" content="dark">', '<meta name="color-scheme" content="light">')
    text = text.replace(
        "__PAYLOAD__", json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    )
    out = OUT_DIR / f"lmm_cinematic_{theme_key}.html"
    out.write_text(text, encoding="utf-8")
    return out


def comparison_page() -> str:
    cards = []
    for key, theme in THEMES.items():
        cards.append(f'''<article class="card"><div class="label">{theme["label"]}</div>
          <img src="{key}_scene04.png" alt="{theme["label"]} 4번 장면 미리보기">
          <div class="actions"><a href="lmm_cinematic_{key}.html">전체 13장 보기</a><a href="{key}_cover.png">표지 보기</a></div></article>''')
    return '''<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>LMM 디자인 A/B/C 비교</title><style>
      *{box-sizing:border-box}body{margin:0;background:#10141a;color:#eef3f7;font-family:"Noto Sans KR","Malgun Gothic",sans-serif;padding:40px}
      header{max-width:1500px;margin:0 auto 28px}h1{margin:0 0 8px;font-size:clamp(30px,4vw,58px);letter-spacing:-.05em}p{color:#9dacb7;margin:0}
      main{max-width:1500px;margin:auto;display:grid;grid-template-columns:repeat(3,1fr);gap:18px}.card{border:1px solid #2d3944;background:#171d24;padding:12px}
      .label{font:800 13px/1.3 Consolas,monospace;letter-spacing:.08em;margin:4px 4px 12px}.card img{display:block;width:100%;aspect-ratio:16/9;object-fit:cover;background:#050910;border:1px solid #2c3944}
      .actions{display:flex;gap:8px;margin-top:12px}.actions a{flex:1;text-align:center;color:#eaf2f7;text-decoration:none;border:1px solid #42515e;padding:10px;font-size:13px}.actions a:hover{border-color:#62dce1;color:#62dce1}
      @media(max-width:1000px){main{grid-template-columns:1fr}body{padding:22px}}
    </style></head><body><header><h1>같은 내용, 세 가지 시각 언어</h1><p>4번 장면을 먼저 비교한 뒤 각 안의 전체 13장 흐름을 열어보세요.</p></header><main>''' + "".join(cards) + "</main></body></html>"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = base.build_payload()
    payload["east_sea"] = [
        {"name": "울릉도", "lat": 37.4845, "lon": 130.9057},
        {"name": "독도", "lat": 37.2429, "lon": 131.8664},
    ]
    for key in THEMES:
        out = render_variant(key, payload)
        print(f"[saved] {out}")
    compare = OUT_DIR / "index.html"
    compare.write_text(comparison_page(), encoding="utf-8")
    print(f"[saved] {compare}")


if __name__ == "__main__":
    main()
