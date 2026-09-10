# -*- coding: utf-8 -*-
"""
현장 자료취득 «직원» 교육 — 웹 시네마틱
(`docs/지자기_현장교육_CINEMATIC_EDITION.html`)
====================================================================

    python create_field_training.py

지자기 측량이 무엇이고 왜 당신의 기록 하나가 중요한가를 전달하는 스크롤형
교육 자료. **22장면 · 약 35분**(질의 별도). 실시일 `TRAINING_DATE`.

## 이 자료가 다른 발표자료와 다른 점

`create_lmm_cinematic.py` 는 전문가·발주자용이라 「모델이 이만큼 나왔다」를 말한다.
이 자료는 **현장 직원용**이라 「당신이 이렇게 기록하지 않으면 이런 일이 벌어진다」를
말한다. 그래서 가상 사례를 쓰지 않고 **이 프로젝트가 실제로 겪은 일**만 싣는다:

| 장면 | 실화 |
|---|---|
| 시각을 안 적으면 | 야장 68**건** 전수조사 — 총자력 절 36건, 그중 F 측정시각 **0건** → 외부장 보정이 막혔다 |
| 좌표가 어긋나면 | 미원 353 m · 남양은 세 자료가 전부 다르고 **현재 좌표 미상** · 서산 시트에 남양 좌표가 통째로 복사 |
| 같은 자리를 다시 가면 | 같은 날 산포는 1.4분인데 **재방문 잔여는 34분** — 방문마다 방위표지가 달라져 있었다 |
| 나쁜 자리를 고르면 | 이원 여의저수지 동쪽 **250.72 nT/m** — 참고 기준 3 nT/m 의 83배 |

## ⚠️ 「모델」을 말하지 않는다 (2026-09-14 발주자 지적)

1차본은 3부 앞에 「지금 모델은 얼마나 맞습니까」 장면을 두고 LOO 편각 29.9분 ·
총자력 62.8 nT 를 실었다. **쓸 수 있는 지역 자기장 모델은 아직 없다** —
있는 것처럼 말하면 직원이 「이미 다 되어 있는데 왜 또 재나」로 듣는다.
그 장면을 통째로 빼고 `lmm_model.json` 의존도 끊었다. 4부의 재방문 장면도
「모델이 26 % 망가졌다」에서 **「같은 자리를 다시 갔더니 34분이 어긋났다」**
로 바꿨다 — 모델 성능이 아니라 관측 자체로 말한다.

## ⚠️ 올해 하는 것은 «절대측정이 아니다»

발주자 확정(2026-09-07) — 올해 50점에서 하는 것은 **자기교란(구배) 측정**이고
DI-flux 절대측정은 선점이 끝난 뒤의 별도 단계다. 편각 4방향·복각 4방향·Null
측정은 **「앞으로 이 자리에서 이런 관측을 합니다」라는 미래 맥락**으로만 짧게
보이고, **올해 실제 작업**(P0 → 4방향 × 11점 → P0 · 수직 20~200 cm · 사진 9칸 ·
카드 7구역)이 3부 본론이다.

⚠️ 이 원칙을 어기면 직원이 「내가 DI-flux 를 배워야 하는구나」로 이해한다.
같은 이유로 **4부의 좌표·재방문 장면에는 「올해가 아니라 내년」 띠**를 붙였다 —
그 둘은 절대측정 캠페인의 이야기다.

## 수치는 전부 파일에서 읽는다

`trial_survey_points`(50점) · `trial_survey_spec`(측선·참고값·발주자 결정) ·
`existing_pts.geojson`(관측망 30) · `korea_boundary.geojson`.
**하드코딩 금지** — 명단이 바뀌면 자동으로 따라온다.

⚠️ 숫자 네 갈래를 섞지 말 것 — **관측망 30** / 2022~25 측량 15 / LMM 투입 16 /
**올해 선점 검토 50**. 교육 대상이 헷갈리는 첫 번째 지점이다.

## 발표 연출 층 · 한 화면에 담기

HUD(장면 이름·전체화면·H/F 단축키)·모서리 프레임·장면 전환 섬광·포인터 조명은
**발표용 연출**이고 본문을 건드리지 않는다. ⚠️ 이 층은 손으로 고친 HTML 에만
있었다 — 생성기를 다시 돌리면 둘로 갈라지므로 2026-09-14 에 **여기로 접어
넣었다.** 연출을 손보려면 이 파일을 고치고 다시 돌린다.

⚠️ **발표 자료는 「장면 하나 = 화면 하나」여야 한다.** 1366×768 에서 재어 보니
23장 중 **14장이 화면을 넘었고 최악이 1.98배**였다. 셋으로 잡는다:

1. `canvas.fig,svg.figsvg{max-height:50vh(짧은 화면 44vh);width:auto}` —
   도판이 화면을 다 먹지 못하게 한다. 이것만으로 14 → 7 장으로 준다.
2. `@media(max-height:900px)` 에서 여백·글자·표를 한 단계 줄이고,
   `.steps`·`.rules` 를 넓은 화면에서 2열로 편다.
3. 그래도 넘치면 **그 장면만** `--fit` 배율로 축소하고 `height` 를 화면
   높이로 못박는다. ⚠️ `transform` 은 레이아웃 상자를 바꾸지 않으므로
   height 를 함께 고정하지 않으면 빈 공간이 남는다.

⚠️ **fit 예약에 `requestAnimationFrame` 을 쓰지 말 것** — 보이지 않는 탭에서
멈춘다. 배경 탭으로 열어 두었다가 넘어오면 장면이 안 맞은 채로 뜬다.
`setTimeout` + `visibilitychange` + 장면 진입 IntersectionObserver 로 다시 잰다.

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
OUT = ROOT / "docs" / "지자기_현장교육_CINEMATIC_EDITION.html"

# ⚠️ 교육 «실시일»이지 파일 생성일이 아니다. today() 로 두면 다시 돌릴 때마다
#    표지 날짜가 바뀐다 — 인쇄본·배포본과 어긋난다.
TRAINING_DATE = "2026년 9월 14일"
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
DRIFT_SITE = (36.5, 127.5, "국토 중앙 36.5°N 127.5°E")
DRIFT_SPAN = (1990, 2026)
DRIFT_CHECK = [(33.5, 126.5), (36.5, 127.5), (37.8, 128.9)]   # 남·중·북동


def declination(lat, lon, year):
    """IGRF-14 편각(도). 어림값을 쓰지 않고 그때그때 계산한다."""
    import datetime as _dt
    import numpy as _np
    import ppigrf
    e, n, _ = ppigrf.igrf(lon, lat, 0, _dt.datetime(year, 1, 1))
    return float(_np.degrees(_np.arctan2(_np.ravel(e)[0], _np.ravel(n)[0])))


def declination_drift():
    import math
    lat, lon, name = DRIFT_SITE
    y0, y1 = DRIFT_SPAN
    curve = [[y, round(declination(lat, lon, y), 4)]
             for y in range(y0, y1 + 1, 2)]
    if curve[-1][0] != y1:
        curve.append([y1, round(declination(lat, lon, y1), 4)])
    d0, d1 = curve[0][1], curve[-1][1]
    rates = []
    for la, lo in DRIFT_CHECK:
        rates.append(abs(declination(la, lo, y1) - declination(la, lo, y0))
                     * 60 / (y1 - y0))
    return {
        "site": name, "y0": y0, "y1": y1, "curve": curve,
        "d0": d0, "d1": d1,
        "delta_min": round(abs(d1 - d0) * 60, 1),
        "offset_m": round(math.tan(math.radians(abs(d1 - d0))) * 1000, 1),
        "rate_lo": round(min(rates), 2), "rate_hi": round(max(rates), 2),
    }


def declination_now(year=2026):
    """올해 편각 — 지역마다 다르므로 기준점 값과 «폭»을 함께 낸다."""
    lat, lon, name = DRIFT_SITE
    d = declination(lat, lon, year)
    vals = [declination(la, lo, year) for la, lo in DRIFT_CHECK]
    a = abs(d)
    return {
        "dec_now": round(d, 3),
        "dec_dms": f"{int(a)}\u00b0{round((a - int(a)) * 60):02d}\u2032W",
        "dec_lo": round(min(abs(v) for v in vals), 1),
        "dec_hi": round(max(abs(v) for v in vals), 1),
        "dec_site": name, "dec_year": year,
    }


def build_payload():
    import trial_survey_points as TP
    import trial_survey_spec as SP

    pts = TP.load_points()
    for i, p in enumerate(pts, 1):
        p.setdefault("_sid", f"S{i:02d}")
    kinds = Counter(p["구분"] for p in pts)

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
        "generated": TRAINING_DATE,
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
        "drift": declination_drift(),
        **declination_now(),
    }


# ══════════════════════════════════════════════════════════════
# 도판 삽화 — 실장비 형태를 딴 인라인 SVG
#
# ⚠️ 캔버스 선그림이 「너무 단순하다」는 지적을 받았다(2026-09-14).
#    사진을 쓰면 외부 파일이 붙어 «오프라인 단일 파일» 원칙이 깨지므로,
#    실제 장비(GEM GSM-19T · MinGeo 010B)의 형태를 그대로 딴 음영 삽화를
#    SVG 로 그려 인라인한다. 벡터라 빔프로젝터에서도 깨지지 않는다.
#
# ⚠️ GNSS 는 «모델이 확정되지 않았다» — 저장소에도
#    `trial_survey_spec.EQUIPMENT` 의 모델 칸이 비어 있다. 임의로 특정
#    제조사를 그리지 말고 일반형 로버로 두고 「모델 확인 중」이라 적는다.
# ══════════════════════════════════════════════════════════════
SVG_GEAR = r'''<svg class="figsvg" viewBox="0 0 1180 530" role="img"
 aria-label="올해 쓰는 GSM-19T 오버하우저 자력계와 GNSS 수신기, 내년 쓰는 MinGeo 010B DI-flux 자기경위의">
<defs>
 <linearGradient id="gPole" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#232b34"/><stop offset=".32" stop-color="#66737f"/>
  <stop offset=".62" stop-color="#38424c"/><stop offset="1" stop-color="#171d23"/></linearGradient>
 <linearGradient id="gSens" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#fdfaf1"/><stop offset=".42" stop-color="#e7e0ce"/>
  <stop offset="1" stop-color="#ada596"/></linearGradient>
 <linearGradient id="gBrass" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#f4dc96"/><stop offset="1" stop-color="#9c7c28"/></linearGradient>
 <linearGradient id="gBox" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#525d69"/><stop offset="1" stop-color="#1e242c"/></linearGradient>
 <linearGradient id="gYel" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#b98c22"/><stop offset=".28" stop-color="#f8d155"/>
  <stop offset=".62" stop-color="#e6b833"/><stop offset="1" stop-color="#9c7318"/></linearGradient>
 <linearGradient id="gWht" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#7f8b98"/><stop offset=".34" stop-color="#eef2f6"/>
  <stop offset="1" stop-color="#6f7a87"/></linearGradient>
 <linearGradient id="gOrg" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#f08a52"/><stop offset="1" stop-color="#b0501f"/></linearGradient>
 <radialGradient id="gShd"><stop offset="0" stop-color="#000" stop-opacity=".6"/>
  <stop offset="1" stop-color="#000" stop-opacity="0"/></radialGradient>
</defs>

<!-- ═══ ① GEM GSM-19T 오버하우저 자력계 ═══ -->
<g>
 <ellipse cx="205" cy="432" rx="92" ry="12" fill="url(#gShd)"/>
 <rect x="198" y="118" width="14" height="292" rx="4" fill="url(#gPole)"/>
 <path d="M198 410h14l-7 26z" fill="#39434f"/>
 <rect x="192" y="110" width="26" height="14" rx="3" fill="#222a33"/>
 <rect x="152" y="62" width="106" height="48" rx="24" fill="url(#gSens)"
  stroke="#8d8579" stroke-width="1.6"/>
 <rect x="166" y="72" width="78" height="11" rx="5.5" fill="#fff" opacity=".6"/>
 <rect x="140" y="74" width="15" height="24" rx="4" fill="url(#gBrass)"/>
 <circle cx="138" cy="86" r="6" fill="url(#gBrass)"/>
 <path d="M247 100 C 300 112 330 176 324 248" fill="none" stroke="#333b45"
  stroke-width="6" stroke-linecap="round"/>
 <rect x="272" y="248" width="106" height="98" rx="9" fill="url(#gBox)"
  stroke="#6d7986" stroke-width="1.6"/>
 <rect x="284" y="260" width="82" height="48" rx="3" fill="#9db5a6"/>
 <rect x="288" y="266" width="52" height="4" rx="2" fill="#5d7a68"/>
 <rect x="288" y="276" width="66" height="4" rx="2" fill="#5d7a68"/>
 <rect x="288" y="286" width="40" height="4" rx="2" fill="#5d7a68"/>
 <g fill="#414c58">
  <rect x="286" y="316" width="22" height="16" rx="3"/>
  <rect x="313" y="316" width="22" height="16" rx="3"/>
  <rect x="340" y="316" width="22" height="16" rx="3"/></g>
 <path d="M272 268 l-22 -10 M378 268 l20 -10" stroke="#3a434e" stroke-width="5"
  stroke-linecap="round" fill="none"/>
</g>

<!-- ═══ ② GNSS 수신기 ═══ -->
<g>
 <ellipse cx="592" cy="432" rx="78" ry="11" fill="url(#gShd)"/>
 <g fill="none" stroke="#3fd8a0" stroke-width="3.4" stroke-linecap="round" opacity=".8">
  <path d="M548 104a62 62 0 0 1 88 0"/>
  <path d="M532 82a86 86 0 0 1 120 0"/></g>
 <rect x="586" y="176" width="14" height="234" rx="4" fill="url(#gPole)"/>
 <path d="M586 410h14l-7 24z" fill="#39434f"/>
 <path d="M531 166a61 32 0 0 1 122 0z" fill="url(#gWht)" stroke="#5e6975" stroke-width="1.6"/>
 <ellipse cx="592" cy="166" rx="61" ry="16" fill="#9aa6b3" stroke="#5e6975" stroke-width="1.6"/>
 <ellipse cx="592" cy="142" rx="30" ry="9" fill="#dfe6ed" opacity=".65"/>
 <rect x="580" y="166" width="24" height="14" rx="3" fill="#2b3440"/>
 <rect x="628" y="238" width="92" height="88" rx="8" fill="url(#gBox)"
  stroke="#6d7986" stroke-width="1.6"/>
 <rect x="638" y="248" width="72" height="52" rx="3" fill="#8fa8bd"/>
 <rect x="642" y="254" width="46" height="4" rx="2" fill="#4f6a83"/>
 <rect x="642" y="264" width="58" height="4" rx="2" fill="#4f6a83"/>
 <g fill="#414c58">
  <rect x="640" y="306" width="20" height="12" rx="2"/>
  <rect x="666" y="306" width="20" height="12" rx="2"/>
  <rect x="692" y="306" width="20" height="12" rx="2"/></g>
 <path d="M600 274 h28" stroke="#4a545f" stroke-width="8" stroke-linecap="round"/>
</g>

<!-- ═══ ③ MinGeo 010B DI-flux 자기경위의 ═══ -->
<g>
 <ellipse cx="980" cy="440" rx="112" ry="13" fill="url(#gShd)"/>
 <g stroke="#8a6534" stroke-width="8" stroke-linecap="round" fill="none">
  <path d="M980 296 L 896 432"/><path d="M980 296 L 1064 432"/><path d="M980 296 L 994 426"/></g>
 <path d="M940 272 h80 l-11 26 h-58z" fill="url(#gYel)" stroke="#7d5f12" stroke-width="1.4"/>
 <rect x="928" y="244" width="104" height="30" rx="6" fill="url(#gYel)"
  stroke="#7d5f12" stroke-width="1.4"/>
 <circle cx="928" cy="259" r="11" fill="#e6b833" stroke="#7d5f12" stroke-width="1.4"/>
 <circle cx="1032" cy="259" r="11" fill="#e6b833" stroke="#7d5f12" stroke-width="1.4"/>
 <rect x="936" y="176" width="27" height="70" rx="6" fill="url(#gYel)"
  stroke="#7d5f12" stroke-width="1.4"/>
 <rect x="998" y="176" width="27" height="70" rx="6" fill="url(#gYel)"
  stroke="#7d5f12" stroke-width="1.4"/>
 <circle cx="949" cy="205" r="18" fill="#f7e8b4" stroke="#7d5f12" stroke-width="2"/>
 <circle cx="949" cy="205" r="6" fill="#c9a233"/>
 <rect x="898" y="186" width="166" height="30" rx="15" fill="url(#gYel)"
  stroke="#7d5f12" stroke-width="1.4"/>
 <rect x="884" y="190" width="18" height="22" rx="5" fill="#2b3138"/>
 <circle cx="1068" cy="201" r="14" fill="#2b3138"/>
 <circle cx="1068" cy="201" r="7" fill="#4c5a68"/>
 <rect x="938" y="148" width="86" height="32" rx="6" fill="url(#gOrg)"
  stroke="#8a3a12" stroke-width="1.6"/>
 <rect x="948" y="156" width="66" height="7" rx="3.5" fill="#ffd0b4" opacity=".7"/>
 <path d="M932 160 H 862" stroke="#ff9a6a" stroke-width="1.2" fill="none"/>
 <text x="856" y="164" fill="#ff9a6a" font-size="13" font-weight="700"
  text-anchor="end">플럭스게이트 센서</text>
 <path d="M924 252 H 872" stroke="#e6b833" stroke-width="1.2" fill="none"/>
 <text x="866" y="256" fill="#e6b833" font-size="13" font-weight="700"
  text-anchor="end">비자성 경위의</text>
</g>

<!-- ═══ 배지 · 이름 ═══ -->
<g font-family="inherit" text-anchor="middle">
 <rect x="151" y="14" width="108" height="28" rx="14" fill="rgba(63,216,160,.10)"
  stroke="#3fd8a0" stroke-width="1.4"/>
 <text x="205" y="33" fill="#3fd8a0" font-size="14" font-weight="700">올해 사용</text>
 <rect x="538" y="14" width="108" height="28" rx="14" fill="rgba(63,216,160,.10)"
  stroke="#3fd8a0" stroke-width="1.4"/>
 <text x="592" y="33" fill="#3fd8a0" font-size="14" font-weight="700">올해 사용</text>
 <rect x="926" y="14" width="108" height="28" rx="14" fill="rgba(155,127,232,.12)"
  stroke="#9b7fe8" stroke-width="1.4"/>
 <text x="980" y="33" fill="#9b7fe8" font-size="14" font-weight="700">내년 사용</text>

 <text x="205" y="474" fill="#eaf2f8" font-size="18" font-weight="700">오버하우저 자력계</text>
 <text x="205" y="496" fill="#8ba3b8" font-size="14" font-weight="600">GEM GSM-19T · 총자력만</text>
 <text x="592" y="474" fill="#eaf2f8" font-size="18" font-weight="700">GNSS 수신기</text>
 <text x="592" y="496" fill="#8ba3b8" font-size="14" font-weight="600">기관 보유 장비 · 모델 확인 중</text>
 <text x="980" y="474" fill="#eaf2f8" font-size="18" font-weight="700">DI-flux 자기경위의</text>
 <text x="980" y="496" fill="#8ba3b8" font-size="14" font-weight="600">MinGeo 010B · 편각·복각</text>
 <text x="590" y="522" fill="#8ba3b8" font-size="14.5" font-weight="600">올해는 세기(F)만 잽니다 — 방향(D·I)을 재는 장비는 내년에 들어옵니다</text>
</g>
</svg>'''

SVG_CLEAN = r'''<svg class="figsvg" viewBox="0 0 700 640" role="img"
 aria-label="현장 작업자와 몸에서 빼 두어야 할 자성 물품 일곱 가지">
<defs>
 <linearGradient id="cSkin" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#d9a173"/><stop offset=".45" stop-color="#f0c69c"/>
  <stop offset="1" stop-color="#c4915f"/></linearGradient>
 <linearGradient id="cShirt" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#1c2f47"/><stop offset=".4" stop-color="#31506f"/>
  <stop offset="1" stop-color="#182739"/></linearGradient>
 <linearGradient id="cVest" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#a9c72c"/><stop offset=".4" stop-color="#dcee5c"/>
  <stop offset="1" stop-color="#8fae1c"/></linearGradient>
 <linearGradient id="cPant" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#232e3f"/><stop offset=".4" stop-color="#3b4a61"/>
  <stop offset="1" stop-color="#1c2534"/></linearGradient>
 <linearGradient id="cHat" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#d29d16"/><stop offset=".35" stop-color="#ffd956"/>
  <stop offset="1" stop-color="#c08c10"/></linearGradient>
 <linearGradient id="cPole" x1="0" y1="0" x2="1" y2="0">
  <stop offset="0" stop-color="#232b34"/><stop offset=".35" stop-color="#66737f"/>
  <stop offset="1" stop-color="#1a2028"/></linearGradient>
 <radialGradient id="cShd"><stop offset="0" stop-color="#000" stop-opacity=".6"/>
  <stop offset="1" stop-color="#000" stop-opacity="0"/></radialGradient>
</defs>

<ellipse cx="212" cy="560" rx="96" ry="14" fill="url(#cShd)"/>

<!-- 자력계 봉 (오른손) -->
<rect x="296" y="196" width="11" height="344" rx="4" fill="url(#cPole)"/>
<path d="M296 540h11l-5.5 20z" fill="#39434f"/>
<rect x="272" y="166" width="60" height="28" rx="14" fill="#e7e0ce" stroke="#8d8579" stroke-width="1.4"/>
<rect x="281" y="172" width="42" height="7" rx="3.5" fill="#fff" opacity=".6"/>

<!-- 다리 · 안전화 -->
<path d="M170 350 h84 l10 164 h-40 l-12-106 -12 106 h-40z" fill="url(#cPant)"/>
<path d="M156 506 h46 v30 a9 9 0 0 1 -9 9 h-52 a9 9 0 0 1 -9 -9 v-9 q0 -9 12 -11z"
 fill="#18202b"/>
<path d="M268 506 h-46 v30 a9 9 0 0 0 9 9 h52 a9 9 0 0 0 9 -9 v-9 q0 -9 -12 -11z"
 fill="#222b39"/>
<rect x="132" y="544" width="70" height="8" rx="4" fill="#333e50"/>
<rect x="222" y="544" width="70" height="8" rx="4" fill="#333e50"/>

<!-- 팔 -->
<g stroke="url(#cShirt)" stroke-width="27" stroke-linecap="round" fill="none">
 <path d="M178 222 C 152 254 144 302 150 350"/>
 <path d="M246 222 C 274 250 288 288 294 322"/></g>
<circle cx="150" cy="352" r="15" fill="#6d7785" stroke="#4d5765" stroke-width="1.5"/>
<circle cx="301" cy="330" r="17" fill="#6d7785" stroke="#4d5765" stroke-width="1.5"/>
<rect x="290" y="322" width="24" height="16" rx="7" fill="#5f6975"/>

<!-- 몸통 -->
<path d="M174 212 c-4 -28 16 -46 38 -50 c22 4 42 22 38 50 l6 140 h-88z" fill="url(#cShirt)"/>
<!-- 안전조끼 -->
<path d="M176 218 h34 l4 134 h-32z" fill="url(#cVest)" opacity=".95"/>
<path d="M248 218 h-34 l-4 134 h32z" fill="url(#cVest)" opacity=".95"/>
<g fill="#e8eef3" opacity=".85">
 <rect x="177" y="256" width="33" height="11"/><rect x="214" y="256" width="33" height="11"/>
 <rect x="178" y="296" width="32" height="11"/><rect x="214" y="296" width="32" height="11"/></g>
<rect x="205" y="216" width="14" height="138" fill="#22374f" opacity=".55"/>

<!-- 목 · 머리 -->
<path d="M200 152 h24 v26 h-24z" fill="#c4915f"/>
<ellipse cx="212" cy="126" rx="31" ry="35" fill="url(#cSkin)"/>
<path d="M181 128 a31 35 0 0 0 62 0z" fill="#000" opacity=".07"/>
<g fill="#3a2c22">
 <ellipse cx="201" cy="132" rx="3.2" ry="4.4"/><ellipse cx="223" cy="132" rx="3.2" ry="4.4"/></g>
<path d="M194 122 h13 M217 122 h13" stroke="#3a2c22" stroke-width="2.6"
 stroke-linecap="round" fill="none"/>
<path d="M204 146 q8 5 16 0" stroke="#a5714a" stroke-width="2.2" fill="none" stroke-linecap="round"/>
<!-- 안전모 -->
<path d="M178 116 a34 34 0 0 1 68 0 z" fill="url(#cHat)"/>
<rect x="168" y="110" width="88" height="12" rx="6" fill="#d9a516"/>
<path d="M209 84 h6 v30 h-6z" fill="#b98a0c" opacity=".55"/>
<path d="M186 100 a28 28 0 0 1 24 -16" stroke="#fff" stroke-width="4"
 stroke-linecap="round" fill="none" opacity=".45"/>

<!-- ═══ 자성 물품 — 몸에서 «뺄 수 있는» 것만 ═══ -->
<defs>
 <radialGradient id="hotR"><stop offset="0" stop-color="#e8503f" stop-opacity=".62"/>
  <stop offset="1" stop-color="#e8503f" stop-opacity="0"/></radialGradient>
</defs>
<g font-family="inherit">
 <g stroke="#e8503f" stroke-width="1" opacity=".45" fill="none">
  <path d="M229 134 L 420 135"/><path d="M243 246 L 420 182"/>
  <path d="M191 258 L 420 229"/><path d="M305 334 L 420 276"/>
  <path d="M217 348 L 420 323"/><path d="M249 378 L 420 370"/>
  <path d="M183 392 L 420 417"/></g>
 <g fill="url(#hotR)">
  <circle cx="224" cy="134" r="30"/><circle cx="238" cy="246" r="30"/>
  <circle cx="186" cy="258" r="30"/><circle cx="300" cy="334" r="30"/>
  <circle cx="212" cy="348" r="30"/><circle cx="244" cy="378" r="30"/>
  <circle cx="178" cy="392" r="30"/></g>
 <g fill="#e8503f">
  <circle cx="224" cy="134" r="4.5"/><circle cx="238" cy="246" r="4.5"/>
  <circle cx="186" cy="258" r="4.5"/><circle cx="300" cy="334" r="4.5"/>
  <circle cx="212" cy="348" r="4.5"/><circle cx="244" cy="378" r="4.5"/>
  <circle cx="178" cy="392" r="4.5"/></g>
 <g fill="#f0b5ad" font-size="15.5" font-weight="600">
  <text x="428" y="140">안경테 · 나사</text>
  <text x="428" y="187">가슴 주머니 볼펜</text>
  <text x="428" y="234">무전기</text>
  <text x="428" y="281">손목시계</text>
  <text x="428" y="328">혁대 버클</text>
  <text x="428" y="375">열쇠꾸러미</text>
  <text x="428" y="422">휴대전화</text></g>
 <text x="428" y="470" fill="#8ba3b8" font-size="13.5" font-weight="600">손에서 내려놓을 수 있는 것들입니다</text>
 <text x="428" y="492" fill="#55707f" font-size="13" font-weight="600">센서에 가까울수록 영향이 큽니다</text>
</g>
</svg>'''

SVG_USE = r'''<svg class="figsvg" viewBox="0 0 1180 470" role="img"
 aria-label="현장 기록에서 국가 지자기 성과, 5만분의 1 지형도 자침편차 표기, 나침반과 항법으로 이어지는 흐름">
<defs>
 <linearGradient id="uPaper" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#f2f5f8"/><stop offset="1" stop-color="#c9d2db"/></linearGradient>
 <linearGradient id="uMap" x1="0" y1="0" x2="0" y2="1">
  <stop offset="0" stop-color="#dfe8dd"/><stop offset="1" stop-color="#b7c6b8"/></linearGradient>
 <linearGradient id="uBezel" x1="0" y1="0" x2="1" y2="1">
  <stop offset="0" stop-color="#586574"/><stop offset=".5" stop-color="#22282f"/>
  <stop offset="1" stop-color="#4a5561"/></linearGradient>
</defs>

<!-- ① 현장 기록 -->
<g>
 <rect x="40" y="70" width="248" height="184" rx="3" fill="rgba(255,255,255,.03)"
  stroke="#3fd8a0" stroke-width="1.6"/>
 <rect x="108" y="96" width="112" height="140" rx="3" fill="url(#uPaper)"/>
 <rect x="100" y="86" width="128" height="20" rx="4" fill="#5b6b7a"/>
 <rect x="150" y="80" width="28" height="14" rx="4" fill="#8996a3"/>
 <g fill="#94a3b1">
  <rect x="120" y="120" width="88" height="4" rx="2"/><rect x="120" y="136" width="88" height="4" rx="2"/>
  <rect x="120" y="152" width="70" height="4" rx="2"/><rect x="120" y="168" width="88" height="4" rx="2"/>
  <rect x="120" y="184" width="58" height="4" rx="2"/><rect x="120" y="200" width="80" height="4" rx="2"/></g>
 <g fill="#3fd8a0"><rect x="120" y="119" width="34" height="6" rx="3"/>
  <rect x="120" y="167" width="42" height="6" rx="3"/></g>
 <path d="M232 232 L 262 190 l10 7 l-30 42z" fill="#e0a34a"/>
 <path d="M232 232 l10 -14 l6 4z" fill="#39434f"/>
</g>

<!-- ② 국가 지자기 성과 -->
<g>
 <rect x="324" y="70" width="248" height="184" rx="3" fill="rgba(255,255,255,.03)"
  stroke="#3fd8d0" stroke-width="1.6"/>
 <rect x="386" y="98" width="112" height="134" rx="3" fill="#b9c4cf"/>
 <rect x="380" y="90" width="112" height="134" rx="3" fill="url(#uPaper)"/>
 <g fill="#94a3b1">
  <rect x="392" y="110" width="88" height="4" rx="2"/><rect x="392" y="124" width="88" height="4" rx="2"/>
  <rect x="392" y="138" width="66" height="4" rx="2"/><rect x="392" y="152" width="88" height="4" rx="2"/>
  <rect x="392" y="166" width="74" height="4" rx="2"/></g>
 <circle cx="466" cy="198" r="22" fill="none" stroke="#c0392b" stroke-width="3.4"/>
 <text x="466" y="204" fill="#c0392b" font-size="14" font-weight="700"
  text-anchor="middle" font-family="inherit">고시</text>
</g>

<!-- ③ 1:50,000 지형도 -->
<g>
 <rect x="608" y="70" width="248" height="184" rx="3" fill="rgba(255,255,255,.03)"
  stroke="#4a7fe8" stroke-width="1.6"/>
 <rect x="628" y="84" width="208" height="84" fill="url(#uMap)" stroke="#4a7fe8" stroke-width="1.4"/>
 <g fill="none" stroke="#8fa88f" stroke-width="1.3">
  <path d="M652 168 q26 -44 54 -22 q30 24 56 -10 q22 -28 50 -6"/>
  <path d="M652 152 q28 -36 56 -16 q28 20 52 -10 q20 -24 48 -4"/>
  <path d="M660 134 q26 -26 50 -12 q26 16 46 -8"/></g>
 <g stroke="#4a7fe8" stroke-width=".7" opacity=".45">
  <path d="M680 84v84M732 84v84M784 84v84"/><path d="M628 112h208M628 140h208"/></g>
 <g stroke="#7c8f7c" stroke-width="1.6" fill="none"><path d="M628 126 h96 l30 22 h82"/></g>
 <text x="832" y="98" fill="#2f4f7a" font-size="12" font-weight="700"
  text-anchor="end" font-family="inherit">1:50,000</text>
 <!-- 난외 자침편차 도식 -->
 <g>
  <path d="M676 248 V 204" stroke="#c9d6e2" stroke-width="2.6" fill="none"/>
  <path d="M676 196 l-6 11 h12z" fill="#c9d6e2"/>
  <text x="676" y="191" fill="#c9d6e2" font-size="11.5" text-anchor="middle">★ 진북</text>
  <path d="M676 248 L 654 206" stroke="#4a7fe8" stroke-width="3" fill="none"/>
  <path d="M650 198 l-1 12 l12 -5z" fill="#4a7fe8"/>
  <text x="641" y="192" fill="#4a7fe8" font-size="11.5" font-weight="700"
   text-anchor="middle">자북</text>
  <circle cx="676" cy="248" r="3.5" fill="#e8eef4"/>
  <text x="700" y="226" fill="#4a7fe8" font-size="16" font-weight="700">{{DEC_DMS}}</text>
  <text x="700" y="244" fill="#55707f" font-size="11" font-weight="600">각도는 과장해 그림</text></g>
</g>

<!-- ④ 나침반 -->
<g>
 <rect x="892" y="70" width="248" height="184" rx="3" fill="rgba(255,255,255,.03)"
  stroke="#ff7048" stroke-width="1.6"/>
 <circle cx="1016" cy="152" r="62" fill="url(#uBezel)"/>
 <circle cx="1016" cy="152" r="52" fill="#0d141c" stroke="#7d8894" stroke-width="1.4"/>
 <g stroke="#7d8894" stroke-width="1.6">
  <path d="M1016 104v9M1016 191v9M968 152h9M1055 152h9"/></g>
 <g stroke="#4a5561" stroke-width="1.2">
  <path d="M1050 118l6-6M982 118l-6-6M1050 186l6 6M982 186l-6 6"/></g>
 <text x="1016" y="120" fill="#e8eef4" font-size="12" font-weight="700"
  text-anchor="middle" font-family="inherit">N</text>
 <g transform="rotate(-8.2 1016 152)">
  <path d="M1016 108 l7 44 l-7 8 l-7 -8z" fill="#e8503f"/>
  <path d="M1016 196 l7 -44 l-7 -8 l-7 8z" fill="#dfe6ed"/></g>
 <circle cx="1016" cy="152" r="4.5" fill="#c9a233"/>
 <text x="1016" y="238" fill="#8ba3b8" font-size="13" font-weight="600"
  text-anchor="middle" font-family="inherit">측량 · 항법 · 등산</text>
</g>

<!-- 화살표 -->
<g fill="#8ba3b8" opacity=".7">
 <path d="M294 162 h20 l-6 -7 l14 8 l-14 8 l6 -7 h-20z"/>
 <path d="M578 162 h20 l-6 -7 l14 8 l-14 8 l6 -7 h-20z"/>
 <path d="M862 162 h20 l-6 -7 l14 8 l-14 8 l6 -7 h-20z"/></g>

<!-- 이름 -->
<g text-anchor="middle" font-family="inherit">
 <text x="164" y="292" fill="#eaf2f8" font-size="17.5" font-weight="700">현장 기록</text>
 <text x="164" y="313" fill="#8ba3b8" font-size="13.5" font-weight="600">올해 자리 고르기 → 내년 D·I·F 측량</text>
 <text x="448" y="292" fill="#eaf2f8" font-size="17.5" font-weight="700">국가 지자기 성과</text>
 <text x="448" y="313" fill="#8ba3b8" font-size="13.5" font-weight="600">검증하고 고시합니다</text>
 <text x="732" y="292" fill="#eaf2f8" font-size="17.5" font-weight="700">1:50,000 지형도</text>
 <text x="732" y="313" fill="#8ba3b8" font-size="13.5" font-weight="600">가장자리의 자침편차 표기</text>
 <text x="1016" y="292" fill="#eaf2f8" font-size="17.5" font-weight="700">나침반 · 측량 · 항법</text>
 <text x="1016" y="313" fill="#8ba3b8" font-size="13.5" font-weight="600">그 각도를 보고 방향을 잡습니다</text>
 <text x="590" y="386" fill="#8ba3b8" font-size="14.5" font-weight="600">⚠ 올해 재는 총자력만으로는 편각이 나오지 않습니다 — 선점을 확정하고 «따로» 편각·복각을 재야 이 흐름이 이어집니다.</text>
 <text x="590" y="410" fill="#55707f" font-size="13.5" font-weight="600">그 출발점이 올해 여러분이 고르는 «자리»입니다.</text>
</g>
</svg>'''


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
canvas.fig,svg.figsvg{width:100%;height:auto;display:block;
 border:1px solid var(--line);background:#04080e;border-radius:2px}
svg.figsvg text{font-family:"Pretendard","Noto Sans KR","맑은 고딕",system-ui,sans-serif}
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
/* ── CINEMATIC EDITION : presentation layer ─────────────── */
body{--mx:72%;--my:24%;background:
 radial-gradient(circle at var(--mx) var(--my),rgba(49,133,170,.13),transparent 24%),
 linear-gradient(180deg,#040911 0%,#07111c 52%,#03070c 100%)}
body::before{content:"";position:fixed;inset:0;z-index:0;pointer-events:none;opacity:.16;
 background-image:linear-gradient(rgba(90,165,190,.11) 1px,transparent 1px),linear-gradient(90deg,rgba(90,165,190,.08) 1px,transparent 1px);
 background-size:72px 72px;mask-image:linear-gradient(to bottom,black,transparent 72%)}
section{isolation:isolate;overflow:hidden}
section::after{content:attr(data-t);position:absolute;right:3vw;bottom:-.14em;z-index:-1;
 font:900 clamp(70px,12vw,180px)/1 "Arial Narrow","Pretendard",sans-serif;
 letter-spacing:-.06em;color:rgba(130,196,216,.027);white-space:nowrap;pointer-events:none}
.wrap{position:relative}
.wrap::before{content:"";position:absolute;left:-32px;top:3px;width:1px;height:44px;
 background:linear-gradient(var(--cyan),transparent);box-shadow:0 0 16px var(--cyan)}
.tag{display:flex;align-items:center;gap:12px;text-shadow:0 0 16px rgba(63,216,208,.45)}
.tag::before{content:"";width:28px;height:1px;background:currentColor;box-shadow:0 0 10px currentColor}
h1,h2{text-shadow:0 12px 40px rgba(0,0,0,.34)}
.hl,.hl-g{filter:drop-shadow(0 0 15px rgba(63,216,208,.2))}
.hl-o,.hl-r{filter:drop-shadow(0 0 15px rgba(255,112,72,.18))}
.card{position:relative;overflow:hidden;backdrop-filter:blur(16px);-webkit-backdrop-filter:blur(16px);
 box-shadow:0 18px 55px rgba(0,0,0,.16);transition:transform .3s ease,border-color .3s ease,background .3s ease}
.card::after{content:"";position:absolute;left:0;top:0;width:38%;height:1px;
 background:linear-gradient(90deg,var(--cyan),transparent);opacity:.55}
.card:hover{transform:translateY(-4px);border-color:rgba(63,216,208,.42);background:rgba(25,65,82,.17)}
.card.warn::after{background:linear-gradient(90deg,var(--orange),transparent)}
.quote{background:linear-gradient(90deg,rgba(63,216,208,.055),transparent 70%);
 padding-top:18px;padding-bottom:18px;box-shadow:-16px 0 38px rgba(63,216,208,.035)}
canvas.fig{box-shadow:0 28px 70px rgba(0,0,0,.27),inset 0 0 55px rgba(23,91,117,.05)}
.tblwrap{border:1px solid rgba(140,180,210,.13);background:rgba(4,11,18,.48);padding:5px 16px 12px}
th{color:#71bbb9}.step,.rule{backdrop-filter:blur(12px);transition:.28s ease}
.step:hover,.rule:hover{border-color:rgba(63,216,208,.42);transform:translateX(5px)}
.num{filter:drop-shadow(0 0 24px rgba(63,216,208,.16))}
main>section:first-child{background:radial-gradient(circle at 80% 42%,rgba(24,114,145,.2),transparent 34%)}
main>section:first-child::before{content:"MAGNETIC FIELD / KOREA";position:absolute;right:6vw;top:16vh;
 color:rgba(133,206,224,.18);font:700 11px/1 ui-monospace,monospace;letter-spacing:.35em;
 writing-mode:vertical-rl}
main>section:first-child h1{font-size:clamp(44px,7vw,98px);max-width:940px;letter-spacing:-.035em}
main>section:first-child .lead{max-width:62ch}
main>section:nth-of-type(4),main>section:nth-of-type(6),main>section:nth-of-type(9),main>section:nth-of-type(16),main>section:nth-of-type(20){
 background:linear-gradient(110deg,rgba(10,35,51,.82),rgba(4,10,17,.45));border-top:1px solid rgba(63,216,208,.12)}
.cinema-hud{position:fixed;z-index:60;left:22px;right:22px;top:18px;height:50px;display:flex;
 align-items:center;pointer-events:none;color:#b8cad5;font-size:12px;letter-spacing:.08em}
.hud-brand{display:flex;align-items:center;gap:11px;font-weight:800}.hud-sigil{width:25px;height:25px;border:1px solid var(--cyan);
 border-radius:50%;position:relative;box-shadow:0 0 22px rgba(63,216,208,.28)}
.hud-sigil::before,.hud-sigil::after{content:"";position:absolute;background:var(--cyan);left:50%;top:50%;transform:translate(-50%,-50%)}
.hud-sigil::before{width:13px;height:1px}.hud-sigil::after{height:13px;width:1px}
.hud-current{margin-left:auto;margin-right:18px;color:var(--muted);font-family:ui-monospace,monospace}
.hud-actions{display:flex;gap:7px;pointer-events:auto}.hud-actions button{height:34px;border:1px solid var(--line);
 background:rgba(5,16,26,.72);backdrop-filter:blur(14px);color:#dceaf1;padding:0 13px;cursor:pointer;
 font:700 11px/1 inherit;letter-spacing:.06em;transition:.2s}.hud-actions button:hover{border-color:var(--cyan);color:var(--cyan)}
.hud-actions .primary{border-color:rgba(63,216,208,.45);background:rgba(63,216,208,.08)}
.cinema-corners{position:fixed;inset:11px;z-index:55;pointer-events:none;border:1px solid rgba(126,183,207,.07)}
.cinema-corners::before,.cinema-corners::after{content:"";position:absolute;width:52px;height:52px;border-color:rgba(63,216,208,.35)}
.cinema-corners::before{left:-1px;top:-1px;border-left:1px solid;border-top:1px solid}.cinema-corners::after{right:-1px;bottom:-1px;border-right:1px solid;border-bottom:1px solid}
.hud-hidden .cinema-hud,.hud-hidden .cinema-corners,.hud-hidden .nav{display:none}
.section-flash{position:fixed;inset:0;z-index:90;pointer-events:none;background:var(--cyan);opacity:0}
.section-flash.go{animation:sceneFlash .5s ease-out}
@keyframes sceneFlash{0%{opacity:.055}100%{opacity:0}}
@media(max-width:860px){.cinema-hud{left:13px;right:13px;top:9px}.hud-current{display:none}.hud-actions button:not(.primary){display:none}
 .cinema-corners{display:none}.wrap::before{display:none}section::after{font-size:20vw}.card:hover{transform:none}}
/* ── 한 화면에 담기 ─────────────────────────────────────
   ⚠️ 발표 자료라 «장면 하나 = 화면 하나»여야 한다. 1366×768 에서 재어 보니
   23장 중 14장이 화면을 넘었고 최악은 1.98배였다. 아래 셋으로 잡는다:
   ① 도판이 화면을 다 먹지 않게 높이를 묶고 ② 낮은 화면에서 여백·글자를 줄이고
   ③ 그래도 넘치면 스크립트가 그 장면만 축소한다(--fit). */
canvas.fig,svg.figsvg{max-height:50vh;width:auto;max-width:100%;
 margin-left:auto;margin-right:auto}
.wrap{transform:scale(var(--fit,1));transform-origin:center center}
.wrap.reveal{transform:translateY(22px) scale(var(--fit,1))}
.wrap.reveal.in{transform:scale(var(--fit,1))}
@media(min-width:1100px){.steps,.rules{grid-template-columns:1fr 1fr}}
@media(max-height:900px){
 section{padding-top:clamp(24px,3.4vh,56px);padding-bottom:clamp(24px,3.4vh,56px)}
 h1{font-size:clamp(34px,5vw,68px)}
 h2{margin-bottom:13px;font-size:clamp(24px,3.2vw,40px)}
 .lead{margin-top:11px;font-size:clamp(14.5px,1.3vw,17.5px)}
 .card{padding:14px 15px}
 .card p{font-size:13.5px;line-height:1.5}
 .card h3{font-size:clamp(16px,1.6vw,20px);margin-bottom:6px}
 .quote{margin:16px 0;padding-top:9px;padding-bottom:9px;
  font-size:clamp(17px,2vw,26px)}
 .steps{margin-top:14px;gap:7px}
 .step{padding:9px 12px;grid-template-columns:34px 1fr;gap:11px}
 .step small{font-size:13px;line-height:1.45}
 .rules{margin-top:14px;gap:8px}
 .rule{padding:11px 15px;grid-template-columns:40px 1fr;gap:13px}
 .rule small{font-size:13.5px}
 .num{font-size:clamp(26px,3.4vw,48px)}
 .stat small{margin-top:6px;font-size:13px}
 canvas.fig,svg.figsvg{max-height:44vh}
 .tblwrap{margin-top:12px;padding:3px 12px 8px}
 table{font-size:14px}th,td{padding:8px 10px}
 p+p{margin-top:8px}
}
@media(max-width:860px){.wrap{transform:none}}
</style></head><body>
<div class="cinema-hud" aria-label="발표 제어">
 <div class="hud-brand"><span class="hud-sigil"></span><span>GEOMAGNETIC FIELD</span></div>
 <div class="hud-current" id="hudCurrent">프롤로그</div>
 <div class="hud-actions">
  <button type="button" id="hudToggle" title="화면 요소 숨기기 (H)">HUD</button>
  <button type="button" id="fullScreen" class="primary" title="전체화면 전환 (F)">전체화면&nbsp; ↗</button>
 </div>
</div>
<div class="cinema-corners" aria-hidden="true"></div>
<div class="section-flash" id="sectionFlash" aria-hidden="true"></div>
<div class="prog" id="prog"></div>
<div class="nav" id="nav"></div>
<div id="bg"></div>
<main>

<!-- ══ 프롤로그 ══ -->
<section data-t="프롤로그"><div class="wrap reveal">
 <span class="tag">현장 자료취득 직원 교육 · {{GENERATED}}</span>
 <h1>우리는 보이지 않는<br>국가 기준을 측정합니다</h1>
 <p class="lead">지구자기장은 눈에 보이지 않지만 지도의 북쪽, 나침반의 방향,
 항법 장비의 기준이 됩니다. 그 값이 지금 이 자리에서 얼마인지는
 <b>현장에서 실제로 재어야</b> 확인됩니다.</p>
 <p class="lead">오늘 이야기는 하나입니다 —
 <b class="hl">대한민국 지자기 기준은 현장에서 기록한 단 한 번의 관측에서
 시작됩니다.</b></p>
</div></section>

<section data-t="편각"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">01</span>
  <h2>지도의 북쪽과<br>나침반의 북쪽은 다릅니다</h2>
  <p class="lead">우리나라에서 나침반 바늘은 진북보다
  <b>서쪽으로 약 {{DEC_LO}}~{{DEC_HI}}도</b> 기울어 있습니다 — 하나의 값이
  아니라 <b>지역마다 다릅니다.</b> 이 각도를 <b class="hl">편각(D)</b> 이라
  부릅니다.</p>
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
   <p class="lead" style="margin-top:24px"><b>이 {{N_SITES}}곳이 전부 새 관측점이
   되는 것은 아닙니다</b> — {{N_EXIST}}곳은 이미 관측점이고, 나머지 {{N_NEW}}곳은
   «후보지»라 조사 결과를 보고 나중에 가릅니다.</p>
   <p class="lead">이 {{N_SITES}}곳에서 여러분이 할 일은
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
 <div class="reveal" style="margin-top:30px"><canvas class="fig" id="figLayers"
  width="1180" height="520"
  aria-label="지구 외핵·지각 암석·상공 태양활동이 만드는 자기장과 항공자력·지상 반복관측의 관계"></canvas></div>
 <div class="grid4 reveal" style="margin-top:26px">
  <div class="card"><h3 class="hl">① 주 자기장</h3>
   <p>지구 외핵의 유체 운동에서 생깁니다. 전체의
   <b>대부분(95 % 이상)</b> 을 차지하고 해마다 조금씩 움직입니다.</p></div>
  <div class="card"><h3 class="hl">② 지역장</h3>
   <p>전지구 표준값이 우리나라에서 실제와 얼마나 다른지입니다.
   <b>따로 있는 자기장이 아니라 «표준값과 실제의 차이»</b>이고, 그 차이를
   실제로 재어 두는 것이 이 연구입니다.</p></div>
  <div class="card"><h3 class="hl">③ 지각 자기장</h3>
   <p>땅속 암석이 만듭니다. 몇 m 만 옮겨도 값이 달라지는 것이
   대부분 이것 때문입니다.</p></div>
  <div class="card"><h3 class="hl-o">④ 외부 자기장</h3>
   <p>태양 활동으로 <b>하루 사이에도 변합니다.</b> 그래서 «몇 시 몇 분»에
   쟀는지가 반드시 필요합니다.</p></div>
 </div>
 <div class="quote reveal">④ 때문에 <span class="hl-o">시각을 안 적으면
 그 값을 되살릴 수 없습니다.</span></div>
</div></section>

<section data-t="변한다"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">04</span>
  <h2>한 번 재고 끝낼 수 없습니다</h2>
  <p class="lead">지구자기장은 해마다 변합니다. 우리나라의 편각도
  <b>매년 수분 수준으로 변하지만, 변화량과 방향은 지역과 시기에 따라
  달라집니다.</b> 그래서 「한 해에 몇 분」이라고 하나로 정해 말할 수
  없습니다.</p>
  <p class="lead">그래서 과거 성과만으로는 오늘의 자기장을 설명할 수 없고,
  <b class="hl">주기적으로 다시 재야</b> 합니다. 여러분이 고른 자리에서
  <b>앞으로 5년마다</b> 관측이 이어집니다.</p>
 </div>
 <div class="reveal"><canvas class="fig" id="figDrift" width="640" height="540"
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
   <small>장래 설계 시나리오 — 현 {{N_NETWORK}}점에 «신규» 50점을 더하면 이렇게
   된다는 안입니다. <b>올해 도는 {{N_SITES}}곳과 같은 50 이 아닙니다</b></small></div>
 </div>
 <div class="quote reveal o">관측점이 부족하면 그 지역의 값은
 <b>측정이 아니라 추정</b>이 됩니다.</div>
</div></section>

<!-- ══ 2부 무엇을 ══ -->
<section data-t="2부 · 벡터"><div class="wrap grid2">
 <div class="reveal">
  <span class="tag">2부 — 무엇을 재는가 · 06</span>
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
 <div class="reveal"><canvas class="fig" id="figVector" width="760" height="560"
  aria-label="자기장 벡터의 성분 분해"></canvas></div>
</div></section>

<section data-t="D·I·F"><div class="wrap">
 <div class="reveal"><span class="tag">07</span>
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

<section data-t="장비"><div class="wrap">
 <div class="reveal">
  <span class="tag">08</span>
  <h2>그럼 무엇으로 잽니까</h2>
  <p class="lead">방금 「비자성 경위의」 같은 말이 나왔습니다. 장비 이름이 낯설 수
  있으니 한 번 짚고 넘어가겠습니다. <b class="hl-g">올해 쓰는 것은 앞의 둘</b>
  이고, 셋째는 내년 이야기입니다.</p>
 </div>
 <div class="reveal" style="margin-top:26px">{{SVG_GEAR}}</div>
 <div class="grid3 reveal" style="margin-top:24px">
  <div class="card good"><h3 class="hl-g">① 오버하우저 자력계</h3>
   <p><b>올해 씁니다.</b> 센서 속 액체의 수소 원자핵이 도는 «진동수»를 셉니다.
   그 진동수가 자기장 세기에 비례하는 것이 물리 상수라 <b>눈금을 맞추지 않아도
   절대값이 나옵니다</b>(동작·시각 점검은 그대로 합니다).
   <b>세기(F)만 재고 방향은 못 잽니다.</b></p></div>
  <div class="card good"><h3 class="hl-g">② GNSS 수신기</h3>
   <p><b>올해 씁니다.</b> 중심점의 위도·경도를 <b>십진도 여섯 자리</b>로
   잽니다(방위표지 좌표는 <b>통보된 대상 지점만</b>). 다시 찾아오는
   <b>가장 확실한 단서</b>라 이 값이 틀리면 그 점을 잃습니다.</p></div>
  <div class="card"><h3 style="color:var(--muted)">③ DI-flux 자기경위의</h3>
   <p><b>내년에 씁니다.</b> 편각·복각을 재는 장비이고,
   <b class="hl">비자성 경위의</b> 위에 <b class="hl">플럭스게이트</b> 센서를
   얹은 것입니다. 아래에서 이 둘을 풀어 설명합니다.</p></div>
 </div>
 <div class="grid2 reveal" style="margin-top:16px">
  <div class="card"><h3 class="hl">비자성 경위의가 무엇입니까</h3>
   <p>각도를 재는 측량용 망원경이 «경위의»입니다. 그것을
   <b>자성 재료를 쓰지 않고</b> 만들어(놋쇠·알루미늄·세라믹 등) 자기 영향이
   없음을 검증해 둔 것입니다. 보통 측량기는 강철이 들어 있어 장비 자체가
   자기장을 흐트러뜨립니다.</p></div>
  <div class="card"><h3 class="hl">플럭스게이트가 무엇입니까</h3>
   <p><b>1축 플럭스게이트는 센서축 방향의 자기장 성분을 측정합니다.</b>
   DI-flux 측정에서는 그 성분이 <b>0 에 가까워지는 자세</b>를 찾아 경위의의
   각도를 읽음으로써 편각과 복각을 결정합니다.</p></div>
 </div>
</div></section>

<!-- ══ 3부 올해 할 일 ══ -->
<section data-t="3부 · 올해 할 일"><div class="wrap">
 <div class="reveal">
  <span class="tag o">3부 — 올해 여러분이 하는 일 · 09</span>
  <h2>측량이 아니라<br><span class="hl-o">자리를 고르는 일</span> 입니다</h2>
  <p class="lead">지자기 측량은 그 자리가 <b>조용한 자리일 때만</b> 뜻이 있습니다.
  땅속 암석이나 주변 철구조물 때문에 몇 m 만 움직여도 값이 뛰는 자리라면,
  아무리 정밀하게 재도 그 값은 그 지역을 대표하지 못합니다.</p>
  <p class="lead">그래서 <b class="hl">먼저 자리가 조용한지를 재 둡니다.</b>
  선점을 확정하고, 표지를 세우고, 별도로 편각·복각·총자력을 재는 일은
  <b>그 뒤에 오는 단계</b>입니다. 올해는 그 첫 단계입니다.</p>
 </div>
 <div class="grid2 reveal" style="margin-top:26px">
  <div class="card warn"><h3 class="hl-o">★ 이번 조사는 «법정 선점»이 아닙니다</h3>
   <p style="font-size:15.5px">법으로 정해진 선점 절차를 밟는 것이 아니고,
   <b>국내 합격 기준도 아직 없습니다.</b> 그래서 <b class="hl-o">현장에서
   합격·탈락을 판정하지 않습니다.</b> 값이 크게 나와도 그 자리를 접지 마시고
   나온 대로 적어 주십시오.</p></div>
  <div class="card"><h3 class="hl">확정된 작업 조건 넷</h3>
   <p style="font-size:15.5px">① <b>야간 관측은 하지 않습니다.</b>
   ② 원본은 <b>종이 카드</b>이고 기기 파일은 그 뒤에 붙입니다.
   ③ 판독은 <b>한 자리에 1회</b>입니다.
   ④ Kp·기상은 <b>적기만</b> 하고 «오늘은 안 되겠다»의 기준으로 쓰지 않습니다 —
   어느 값에서 걸러야 할지를 정하려고 이 조사를 합니다.</p></div>
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
  지상 {{H_MIN}}~{{H_MAX}} cm 를 {{H_STEP}} cm 간격으로 잽니다. 높이는
  <b>{{N_HEIGHT}}개</b> 인데, 여기에도 <b class="hl">시작할 때와 끝날 때
  「기준높이」에서 한 번씩</b> 더 재기 때문에 모두
  <b>{{N_READ_V}}줄</b>이 됩니다. 기준높이는 재기 편한 높이 하나(보통 100 cm)를
  골라 카드에 적어 두고 처음·끝에 같은 높이로 잽니다.</p>
  <p class="lead">한 지점에서 모두 <b class="hl">{{N_READ}}번</b>
  재게 됩니다(수평 {{N_READ_H}} + 수직 {{N_READ_V}}).</p>
  <div class="card warn" style="margin-top:20px">
   <p style="font-size:16px;color:var(--ink)"><b class="hl-o">★ 현장 상황에 따라
   달라질 수 있습니다.</b> 바위·물·급경사·사유지 때문에 10 m 를 다 못 채우거나
   방향을 틀어야 할 수 있습니다. 그럴 때 임의로 건너뛰지 마시고
   <b>실제로 잰 방위와 거리, 그리고 그렇게 한 이유</b>를 카드에 적어 주십시오 —
   그 셋이 있어야 사무실에서 그 줄을 살려 쓸 수 있습니다.
   못 잰 칸은 빈칸으로 두지 말고 <b>「측정불가」</b> 로 표시합니다 —
   빈칸은 «안 쟀다»인지 «못 쟀다»인지 나중에 구분되지 않습니다.</p>
  </div>
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
  <b class="hl">자리 때문에 생긴 차이에 훨씬 가까워집니다.</b></p>
  <p class="lead">다만 이것은 두 시각 사이가 <b>직선으로 변했다고 «본»</b> 것이라
  전부를 걷어내지는 못합니다. 그래서 시각을 정확히 적는 일이 더 중요합니다 —
  사무실에서 다르게 계산해 볼 여지를 남겨 두는 것이 시각입니다.</p>
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
 <div class="quote reveal r">한 곳은 동쪽으로 1 m 만 가도 구배의 «크기»가
 <span class="hl-r">{{WORST_VAL}} nT/m</span> 입니다 —
 참고 수준 {{EURO_GRAD}} nT/m 의 <b>{{WORST_RATIO}}배</b> 입니다.</div>
 <p class="lead reveal" style="margin-top:-8px;font-size:15px">표의 부호는
 <b>방향</b>입니다 — 「−」는 중심점보다 값이 «작아졌다»는 뜻이고, 판정에 쓰는 것은
 부호를 뗀 <b>크기(|ΔF| ÷ 거리)</b> 입니다. 그러니 「−250」을 「3보다 작다」로
 읽지 마십시오.</p>
 <p class="lead reveal">이런 자리에서는 정밀하게 재도 소용이 없습니다.
 표석을 몇 cm 만 잘못 찾아도 값이 달라지니까요.
 <b class="hl">그 자리가 이런 자리인지 아닌지를 가려내는 것</b>이
 올해 여러분의 일입니다.</p>
</div></section>

<section data-t="기준이 없다"><div class="wrap">
 <div class="reveal">
  <span class="tag o">13</span>
  <h2>그런데 아직<br><span class="hl-o">합격 기준이 없습니다</span></h2>
  <p class="lead">해외에는 「반경 {{IAGA_RADIUS}} m 안에서 {{IAGA_RANGE}} nT
  이내」「구배 {{EURO_GRAD}} nT/m 미만」이라는 값이 있습니다. 다만
  <b class="hl-r">{{EURO_GRAD}} nT/m 는 특정 해외 사례에서 제시된 참고
  수준이며, 국내 반복관측점의 법정 또는 확정 판정기준이 아닙니다.</b></p>
  <p class="lead">앞 장의 시범관측을 보십시오. 재어 본 <b>두 곳</b> 모두 네 방향
  전부가 {{EURO_GRAD}} nT/m 를 넘었습니다. 두 곳만으로 전국이 그렇다고 말할 수는
  없지만, 이 참고값을 그대로 합격선으로 쓰면 상당수가 걸린다는 것은 분명합니다.
  그 값이 우리 땅에 맞는지는 <b>이번 조사를 다 마치고 나서</b> 판단합니다.</p>
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
  <p class="lead">측정 전에 <b>뺄 수 있는 것</b>부터 몸에서 내려놓아 주십시오 —
  시계·휴대전화·열쇠·펜·무전기처럼 손에서 놓아도 되는 것들입니다.
  사정상 못 뺀 것이 있으면 무엇을 왜 못 뺐는지 적어 주시면 됩니다.</p>
 </div>
 <div class="reveal">{{SVG_CLEAN}}</div>
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
   <small>일자·관측자·기상·Kp·<b>장비 원시파일명과 레코드 범위</b>.
   시각은 <b>모든 F 줄에 시·분·초까지</b></small></div></div>
  <div class="step"><div><b>기준점 좌표</b>
   <small>GNSS 실측 위도·경도를 <b>십진도 여섯 자리</b>로.
   사전좌표와 많이 다르면 그 사유도 함께</small></div></div>
  <div class="step"><div><b>수평 자기구배</b>
   <small>네 방향 × {{N_OFF}}점, 방향마다 P0 앞뒤로 한 번씩</small></div></div>
  <div class="step"><div><b>수직 자기구배</b>
   <small>높이 {{N_HEIGHT}}개 + <b>시작·끝 기준높이 2회 = {{N_READ_V}}줄</b>.
   높이는 목표값이 아니라 실제 센서 «가운데» 높이를 적습니다</small></div></div>
  <div class="step" style="border-color:rgba(155,127,232,.42);
   background:rgba(155,127,232,.06)">
   <div><b>방위표지 시준 <span style="color:#b9a6f2;font-size:15px;
    font-weight:600">— 대상 지점만</span></b>
   <small>정·반 시준과 표지 좌표·취득방법. <b>전 지점에서 하지 않습니다</b> —
   대상 여부를 착수 전에 알려 드리고, 아니면 「해당없음」으로 둡니다</small></div></div>
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
  <p class="lead">지난 자료를 정리하면서 야장 <b>68건을 전부</b> 열어 봤습니다.
  그중 총자력 측정 결과가 적힌 것은 <b>36건</b> 이었는데, 그 36건 어디에도
  <b class="hl-r">몇 시 몇 분에 쟀는지는 없었습니다.</b></p>
 </div>
 <div class="grid3 reveal" style="margin-top:30px">
  <div class="card stat"><span class="num">68</span>
   <small>전수조사한 야장 건수</small></div>
  <div class="card stat"><span class="num">36</span>
   <small>그중 총자력 측정 결과가 적힌 건수</small></div>
  <div class="card stat warn"><span class="num r">0</span>
   <small>그 36건 가운데 총자력을 «몇 시에» 쟀는지 적힌 건수</small></div>
 </div>
 <p class="lead reveal" style="margin-top:24px">편각·복각은 사정이 달랐습니다.
 세션 시각이 <b class="hl-g">212세션 가운데 211세션</b>에 적혀 있어서, 나중에
 시각을 되살려 그 시점의 외부장 보정량을 낼 수 있었습니다.
 <b class="hl-r">총자력만 그러지 못했습니다.</b></p>
 <p class="lead reveal">그리고 지금까지 정리된 관측 가운데 실제로 계산에 들어간
 것은 <b>16개 지점 · 30개 관측행</b> 입니다. 야장에 적힌 것이 전부 쓰이지는
 않습니다 — <b class="hl">쓸 수 있게 적혀 있어야</b> 쓰입니다.</p>
 <div class="quote reveal r">아무리 정확히 재도,
 <b>몇 시에 쟀는지 모르면 그 값은 외부장 보정에 쓸 수 없습니다.</b></div>
 <p class="lead reveal">규정에도 「편각·복각을 잴 때 총자기장과 시간을 함께
 측정한다」고 되어 있었습니다. 다만 <b>야장 양식에 그 칸이 없었습니다.</b>
 그래서 이번 카드에는 <b class="hl">F 를 잰 모든 줄에 시각 칸</b>을 넣었습니다.</p>
</div></section>

<section data-t="좌표"><div class="wrap">
 <div class="reveal">
  <span class="tag" style="color:var(--red)">17</span>
  <h2>좌표가 어긋나서<br><span class="hl-r">점 하나를 잃은 일</span></h2>
 </div>
 <div class="card reveal" style="border-color:rgba(155,127,232,.5);
  background:rgba(155,127,232,.09);margin:0 0 26px">
  <p style="font-size:17px;color:var(--ink)">
  <b style="color:#b9a6f2">★ 과거 성과와 대조하는 일은 내년 몫입니다.</b>
  다만 <b>중심점 좌표를 정확히 남기는 것은 올해 여러분의 책임</b>입니다 —
  아래 사고들이 바로 그 한 줄에서 시작됐습니다.</p>
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
  <h2>같은 자리를 다시 갔더니<br><span class="hl-r">34분이 어긋났습니다</span></h2>
 </div>
 <div class="card reveal" style="border-color:rgba(155,127,232,.5);
  background:rgba(155,127,232,.09);margin:0 0 26px">
  <p style="font-size:16.5px;color:var(--ink)">
  <b style="color:#b9a6f2">★ 올해 하는 일은 아닙니다.</b>
  <b>내년에 실제 측량이 시작되면</b> 확인해야 하는 내용인데, 올해 고른 자리와
  남긴 좌표가 그때 그대로 쓰이기 때문에 지금 알아 두셔야 합니다.</p>
 </div>
 <div class="reveal">
  <p class="lead">같은 표석을 몇 년 뒤 다시 재면 자연스러운 변화만큼만
  달라져야 합니다. 그런데 <b class="hl-r">재방문 자료의 방위표지 보정 후
  잔여 RMS 는 33.7분</b> 이었습니다. 16구간 가운데 9구간에서 어긋남이
  6분을 넘었습니다.</p>
  <p class="lead">장비나 계산 탓이 아니었습니다. 야장 산술은
  <b class="hl-g">0.16초</b> 까지 맞았고 같은 날 두 번 재면
  <b class="hl-g">1.4분</b> 안에 들어옵니다 — <b>정밀한데 부정확했습니다.</b>
  원인은 <b class="hl-r">방문할 때마다 시준한 방위표지가 달라져 있었다</b>는
  것이었습니다. 바뀐 표지의 참방위각이 틀리면 그 방문의 편각이 통째로
  그만큼 틀어집니다.</p>
 </div>
 <div class="reveal" style="margin-top:30px"><canvas class="fig" id="figAudit"
  width="1100" height="252" aria-label="같은 날 재현성과 재방문 재현성 비교"></canvas></div>
 <div class="quote reveal r">다음 사람이 <b>무엇을 보고 쟀는지</b> 알 수 있게
 남겨 주십시오.</div>
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
  <div class="card"><h3 class="hl">② 기준을 정할 때 이 자료를 씁니다</h3>
   <p>국내 판정 기준이 아직 없습니다. <b>이번에 {{N_SITES}}곳에서 모은 값</b>이
   그 기준을 정할 때 쓰는 자료가 됩니다. 한쪽으로 치우친 자료를 모으면
   기준도 그만큼 치우칩니다.</p></div>
 </div>
 <p class="lead reveal" style="margin-top:26px">그래서 「기록을 빠뜨리지 마세요」가
 형식적인 당부가 아닙니다. 시각 하나, 좌표 하나, 차수 하나가 빠지면
 그 줄은 <b class="hl-r">보정 계산에서 빠집니다.</b> 나중에 사무실에서
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
  <div class="card"><h3 class="hl">지금 만들어 가는 것</h3>
   <p>전국의 편각을 설명할 수 있는 기준 자료 · 장기 영년변화 분석 ·
   전지구 표준값과 국내 실측 대조 — <b>아직 만드는 중이고,
   그 재료가 여러분의 관측입니다</b></p></div>
  <div class="card"><h3 style="color:var(--muted)">앞으로 넓혀 갈 수 있는 곳</h3>
   <p>항법 보조 · 우주환경 영향 분석 · 지하자원 탐사 기준 ·
   전력망·철도의 지자기 유도 영향 — <b>가능성이지 확정된 계획은
   아닙니다</b></p></div>
 </div>
 <div class="reveal" style="margin-top:34px">{{SVG_USE}}</div>
 <p class="lead reveal" style="margin-top:22px"><b>여러분의 교란조사는 향후
 편각·복각·총자기장을 안정적으로 측정할 관측점을 결정합니다.</b>
 그 관측점에서 축적된 성과가 <b class="hl">지역 지자기 모델과 국가기본도
 자편각의 기초자료</b>가 됩니다.</p>
 <p class="lead reveal">가장 눈에 보이는 쓰임은 지형도입니다 — 1:50,000
 도엽 가장자리의 자침편차 표기가 그것입니다.</p>
</div></section>

<!-- ══ 에필로그 ══ -->
<section data-t="에필로그"><div class="wrap">
 <div class="reveal">
  <span class="tag">에필로그 · 21</span>
  <h2>현장에서 지킬 다섯 가지</h2>
  <p class="lead">여러분이 현장에서 적은 값 하나, 시각 하나, 좌표 하나가
  데이터베이스로, 전국 지도로, 대한민국 지자기 기준으로 이어집니다.</p>
 </div>
 <div class="quote reveal" style="font-size:clamp(22px,3.2vw,42px);margin:44px 0">
  대한민국 지자기 기준은<br><span class="hl">여러분의 기록 하나에서
  시작됩니다.</span></div>
 <div class="rules reveal">
  <div class="rule"><div><b>같은 지점에서 잽니다</b>
   <small>좌표와 소재지를 정확히 남겨야 다음 사람이 그 자리를 찾아갑니다</small></div></div>
  <div class="rule"><div><b>뺄 수 있는 쇠붙이는 내려놓습니다</b>
   <small>시계·휴대전화·열쇠·펜, 차량과 전자기기 — 못 뺀 것은 적어 둡니다</small></div></div>
  <div class="rule"><div><b>시각을 시·분·초까지 적습니다</b>
   <small>F 를 잰 모든 줄에. 이것 하나가 빠지면 그 줄은 계산에서 빠집니다</small></div></div>
  <div class="rule"><div><b>다시 재는 것은 «절차가 깨졌을 때»입니다</b>
   <small>값이 커서가 아닙니다. 시작한 방향은 값이 어떻게 나오든 끝까지 재고,
   다시 잴 때는 <b>차수를 올려 새 줄</b>에 적습니다 — 앞서 잰 줄은 지우지
   않습니다</small></div></div>
  <div class="rule"><div><b>재현할 수 있게 남깁니다</b>
   <small>차수, 장비 원시파일명, 사진 파일명 — 나중에 되짚을 수 있어야
   자료입니다</small></div></div>
 </div>
</div></section>

<!-- ══ 마무리 ══ -->
<section data-t="마무리"><div class="wrap">
 <div class="reveal">
  <span class="tag">마무리 · 22</span>
  <h2>다 기억하지 않으셔도 됩니다</h2>
  <p class="lead">오늘 말씀드린 것을 전부 외우실 필요는 없습니다.
  현장에서 필요한 것은 <b>카드에 순서대로 다 적혀 있습니다.</b>
  위에서부터 채워 가시면 됩니다.</p>
  <p class="lead">다만 카드가 대신해 줄 수 없는 것이 하나 있습니다 —
  <b class="hl">적을까 말까 망설여지는 순간</b> 입니다.
  그때는 적어 주십시오. 현장에서 10초 걸리는 일이,
  사무실에서는 <b class="hl-r">되살릴 방법이 없는 일</b>이 됩니다.</p>
 </div>
 <div class="grid3 reveal" style="margin-top:34px">
  <div class="card"><h3 class="hl">망설여지면</h3>
   <p>적습니다. 필요 없는 기록은 나중에 빼면 되지만,
   없는 기록은 만들어 낼 수 없습니다.</p></div>
  <div class="card"><h3 class="hl">이상하면</h3>
   <p>그대로 적습니다. 값이 크다고 고치거나 빼지 않습니다 —
   그 값도 기준을 정하는 자료입니다.</p></div>
  <div class="card"><h3 class="hl">막히면</h3>
   <p>혼자 판단하지 마시고 알려 주십시오.
   못 잰 이유가 적혀 있으면 그것도 자료가 됩니다.</p></div>
 </div>
 <div class="quote reveal" style="font-size:clamp(21px,3vw,38px);margin:44px 0 30px">
  오늘 여러분이 남기는 한 줄이<br>
  <span class="hl">앞으로 수십 년을 다시 잴 자리</span>를 정합니다.</div>
 <p class="lead reveal" style="font-size:clamp(17px,1.7vw,23px);color:var(--ink)">
 안전하게 다녀오십시오. 고맙습니다.</p>
</div></section>

</main>
<div class="foot">
 지자기 측량 현장교육 · {{GENERATED}}<br>
 수치는 <code>trial_survey_points</code> · <code>trial_survey_spec</code> ·
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
function nearest() {
  let best = 0, min = Infinity;
  secs.forEach((s, i) => {
    const d = Math.abs(s.getBoundingClientRect().top);
    if (d < min) { min = d; best = i; }
  });
  return best;
}
addEventListener("keydown", e => {
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  if (["ArrowDown", "PageDown", " ", "ArrowRight"].includes(e.key)) {
    e.preventDefault();
    secs[Math.min(nearest() + 1, secs.length - 1)].scrollIntoView({behavior: "smooth"});
  } else if (["ArrowUp", "PageUp", "ArrowLeft"].includes(e.key)) {
    e.preventDefault();
    secs[Math.max(nearest() - 1, 0)].scrollIntoView({behavior: "smooth"});
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
function rrect(ctx, x, y, w, h, r, fill, line) {
  ctx.beginPath();
  ctx.moveTo(x + r, y); ctx.lineTo(x + w - r, y);
  ctx.quadraticCurveTo(x + w, y, x + w, y + r); ctx.lineTo(x + w, y + h - r);
  ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h); ctx.lineTo(x + r, y + h);
  ctx.quadraticCurveTo(x, y + h, x, y + h - r); ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y); ctx.closePath();
  if (fill) { ctx.fillStyle = fill; ctx.fill(); }
  if (line) { ctx.strokeStyle = line; ctx.stroke(); }
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

/* ── 자기장의 출처 단면도 ────────────────────────── */
/* 「지구 핵 · 지각 · 항공자력 · 반복측정」을 한 그림에 둔다. 카드 넷은
   무엇이 있는지만 말하고, 이 그림이 «어디서 와서 어떻게 재는지»를 말한다. */
(function () {
  const cv = document.getElementById("figLayers"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const GY = 322, FLY = 152;
  ctx.clearRect(0, 0, w, h);
  ctx.lineJoin = "round"; ctx.lineCap = "round";
  /* 맨틀 */
  ctx.fillStyle = "rgba(140,180,210,.05)";
  ctx.fillRect(0, GY + 58, w, h - GY - 58);
  /* 외핵 */
  const cg = ctx.createLinearGradient(0, 402, 0, h);
  cg.addColorStop(0, "rgba(63,216,208,.07)");
  cg.addColorStop(1, "rgba(63,216,208,.26)");
  ctx.fillStyle = cg;
  ctx.beginPath(); ctx.moveTo(0, h); ctx.lineTo(0, 470);
  ctx.quadraticCurveTo(w / 2, 402, w, 470); ctx.lineTo(w, h); ctx.closePath(); ctx.fill();
  ctx.strokeStyle = "rgba(63,216,208,.45)"; ctx.lineWidth = 1.6;
  ctx.beginPath(); ctx.moveTo(0, 470); ctx.quadraticCurveTo(w / 2, 402, w, 470); ctx.stroke();
  [620, 820, 1020].forEach(x => {                 /* 쇳물 순환 */
    const r = 30, a = Math.PI * 1.7;
    ctx.strokeStyle = "rgba(63,216,208,.55)"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(x, h - 42, r, Math.PI * .2, a); ctx.stroke();
    ctx.save(); ctx.translate(x + Math.cos(a) * r, h - 42 + Math.sin(a) * r);
    ctx.rotate(a + Math.PI / 2); ctx.fillStyle = "rgba(63,216,208,.75)";
    ctx.beginPath(); ctx.moveTo(0, -9); ctx.lineTo(6, 5); ctx.lineTo(-6, 5);
    ctx.closePath(); ctx.fill(); ctx.restore();
  });
  /* 지각 */
  ctx.fillStyle = "rgba(155,127,232,.13)"; ctx.fillRect(0, GY, w, 58);
  ctx.strokeStyle = "rgba(155,127,232,.55)"; ctx.lineWidth = 1.6;
  [[600, 28, 62, 17], [800, 34, 44, 13], [1010, 24, 54, 15]].forEach(([x, d, rx, ry]) => {
    ctx.beginPath(); ctx.ellipse(x, GY + d, rx, ry, .2, 0, 7); ctx.stroke();
  });
  /* 지표 */
  ctx.strokeStyle = "#8ba3b8"; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(0, GY); ctx.lineTo(w, GY); ctx.stroke();
  /* 항공자력 측선 */
  ctx.strokeStyle = "rgba(74,127,232,.6)"; ctx.lineWidth = 1.6;
  ctx.setLineDash([9, 7]);
  ctx.beginPath(); ctx.moveTo(60, FLY); ctx.lineTo(w - 60, FLY); ctx.stroke();
  ctx.setLineDash([]);
  ctx.strokeStyle = "rgba(74,127,232,.22)"; ctx.lineWidth = 1;
  for (let x = 90; x < w - 60; x += 46) {
    ctx.beginPath(); ctx.moveTo(x, FLY + 6); ctx.lineTo(x, GY - 4); ctx.stroke();
  }
  (function plane(x, y) {                          /* 항공기 */
    ctx.fillStyle = "#4a7fe8";
    ctx.beginPath(); ctx.moveTo(x + 28, y); ctx.lineTo(x - 16, y - 9);
    ctx.lineTo(x - 7, y); ctx.lineTo(x - 16, y + 9); ctx.closePath(); ctx.fill();
    ctx.beginPath(); ctx.moveTo(x + 6, y - 2); ctx.lineTo(x - 10, y - 22);
    ctx.lineTo(x - 1, y - 22); ctx.lineTo(x + 15, y - 2); ctx.closePath(); ctx.fill();
  })(700, FLY);
  /* 지상 반복관측점 */
  [600, 800, 1010].forEach(x => {
    ctx.strokeStyle = "#3fd8a0"; ctx.lineWidth = 2;
    ctx.fillStyle = "rgba(63,216,160,.25)";
    ctx.beginPath(); ctx.moveTo(x - 12, GY); ctx.lineTo(x + 12, GY);
    ctx.lineTo(x + 9, GY - 17); ctx.lineTo(x - 9, GY - 17); ctx.closePath();
    ctx.fill(); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x, GY - 17); ctx.lineTo(x, GY - 62); ctx.stroke();
    ctx.beginPath(); ctx.arc(x, GY - 71, 9, 0, 7);
    ctx.fillStyle = "rgba(63,216,160,.4)"; ctx.fill(); ctx.stroke();
  });
  /* 태양 */
  const SX = w - 62, SY = 60;
  ctx.fillStyle = "rgba(255,112,72,.28)";
  ctx.beginPath(); ctx.arc(SX, SY, 23, 0, 7); ctx.fill();
  ctx.strokeStyle = "#ff7048"; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(SX, SY, 23, 0, 7); ctx.stroke();
  for (let i = 0; i < 10; i++) {
    const a = i / 10 * Math.PI * 2;
    ctx.beginPath();
    ctx.moveTo(SX + Math.cos(a) * 29, SY + Math.sin(a) * 29);
    ctx.lineTo(SX + Math.cos(a) * 37, SY + Math.sin(a) * 37); ctx.stroke();
  }
  ctx.strokeStyle = "rgba(255,112,72,.42)"; ctx.lineWidth = 1.6;
  [0, 26, 52].forEach(off => {                     /* 흔들리는 외부장 */
    ctx.beginPath();
    for (let t = 0; t <= 1; t += .02) {
      const x = SX - 40 - t * 420, y = SY + 44 + off + t * 190 + Math.sin(t * 22) * 9;
      if (t === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  });
  /* 이름표 — 왼쪽 한 줄로 모아 겹치지 않게 한다 */
  /* ⚠️ 발생원(자기장이 «오는 곳»)과 관측방법(그것을 «재는 법»)을 같은
     모양으로 찍으면 항공자력·반복관측이 자기장의 원인처럼 보인다.
     채운 사각형 = 발생원 · 빈 사각형 = 재는 방법 으로 갈라 둔다. */
  function lab(y, col, t1, t2, how) {
    if (how) {
      ctx.strokeStyle = col; ctx.lineWidth = 2;
      ctx.strokeRect(73, y - 9, 8, 8);
    } else {
      ctx.fillStyle = col; ctx.fillRect(72, y - 10, 10, 10);
    }
    ctx.fillStyle = how ? "#c9d8e4" : "#eaf2f8";
    ctx.font = "700 14.5px sans-serif"; ctx.textAlign = "left";
    ctx.fillText((how ? "[재는 방법] " : "") + t1, 90, y);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "600 13px sans-serif";
    ctx.fillText(t2, 90, y + 19);
  }
  lab(46, "#ff7048", "태양 활동 — 상공에서 오는 자기장",
      "하루 사이에도 값을 흔듭니다. 그래서 «몇 시에 쟀는지»가 필요합니다");
  lab(112, "#4a7fe8", "항공자력 측량 — 하늘에서 넓게 훑습니다",
      "지각 자기장의 «모양»을 봅니다. 이미 확보되어 있는 자료입니다", true);
  lab(250, "#3fd8a0", "지상 반복관측 — 한 점에서 정확하게 잽니다",
      "5년마다 다시 갑니다. 여러분이 하는 일이 바로 이것입니다", true);
  lab(356, "#9b7fe8", "지각 암석 — 땅속 암석이 만드는 국지 자기장",
      "몇 m 만 옮겨도 값이 달라지는 원인입니다");
  lab(486, "#3fd8d0", "지구 외핵 — 쇳물의 흐름",
      "자기장의 대부분을 만들고, 해마다 조금씩 움직입니다");
})();

/* ── 나침반 ──────────────────────────────────────── */
(function () {
  const cv = document.getElementById("figCompass"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const cx = w / 2, cy = h / 2, R = Math.min(w, h) * .36;
  const DEC = P.dec_now;                       /* IGRF-14 계산값 */
  let ang = 0, t0 = null;
  function draw(ts) {
    if (t0 === null) t0 = ts;
    const k = Math.min((ts - t0) / 2200, 1);
    ang = DEC * (k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2);
    ctx.clearRect(0, 0, w, h);
    ctx.strokeStyle = "rgba(140,180,210,.22)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
    ctx.beginPath(); ctx.arc(cx, cy, R * .74, 0, 7); ctx.stroke();
    /* 진북 */
    ctx.strokeStyle = "#8ba3b8"; ctx.lineWidth = 1.6;
    ctx.setLineDash([5, 5]);
    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx, cy - R * 1.10); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "700 15px sans-serif";
    ctx.textAlign = "center"; ctx.fillText("진북 (지도의 북)", cx, cy - R * 1.21);
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
    ctx.fillText("서편각 · " + P.dec_site, cx - R * .62, cy - R * .28 + 20);
    /* ⚠️ 자북 이름을 바늘 «머리 위»에 두면 편각이 8도뿐이라 「진북」 위에
       겹쳐 찍힌다. 지시선을 왼쪽으로 빼서 따로 앉힌다. */
    const tipX = cx + Math.sin(rad + Math.PI / 2) * R * .98;
    const tipY = cy - Math.cos(rad + Math.PI / 2) * R * .98;
    const lbX = cx - R * .56, lbY = cy - R * .99;
    ctx.strokeStyle = "rgba(63,216,208,.5)"; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(tipX - 6, tipY + 3); ctx.lineTo(lbX + 8, lbY - 4);
    ctx.stroke();
    ctx.fillStyle = "#3fd8d0"; ctx.font = "700 15px sans-serif";
    ctx.textAlign = "right";
    ctx.fillText("자북 (나침반 바늘)", lbX, lbY);
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
/* ⚠️ 꺾은선 그래프로는 «움직인다»가 안 느껴진다 — 1차본이 그랬다.
   바늘을 실제로 돌리되, 36년치 차이가 1.8° 뿐이라 맨눈으로는 안 보이므로
   돋보기로 8배 키워 세 해의 방향을 나란히 보인다. */
(function () {
  const cv = document.getElementById("figDrift"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const cx = w * .41, cy = h * .58, R = Math.min(w, h) * .30;
  const DR = P.drift, CV = DR.curve;             /* IGRF-14 계산값 */
  const D90 = DR.d0, D26 = DR.d1;
  function dAt(k) {                              /* 곡선 위를 따라간다 */
    const t = DR.y0 + (DR.y1 - DR.y0) * k;
    for (let i = 1; i < CV.length; i++) {
      if (t <= CV[i][0]) {
        const f = (t - CV[i - 1][0]) / (CV[i][0] - CV[i - 1][0]);
        return CV[i - 1][1] + (CV[i][1] - CV[i - 1][1]) * f;
      }
    }
    return CV[CV.length - 1][1];
  }
  const MX = w * .78, MY = h * .28, MR = Math.min(w, h) * .17;
  function needle(deg, col, alpha, dash) {
    const rad = deg * Math.PI / 180;
    ctx.save(); ctx.translate(cx, cy); ctx.rotate(rad); ctx.globalAlpha = alpha;
    if (dash) {
      ctx.setLineDash([6, 5]); ctx.strokeStyle = col; ctx.lineWidth = 2.2;
      ctx.beginPath(); ctx.moveTo(0, 0); ctx.lineTo(0, -R * .95); ctx.stroke();
      ctx.setLineDash([]);
    } else {
      ctx.fillStyle = col;
      ctx.beginPath(); ctx.moveTo(0, -R * .95);
      ctx.lineTo(9, 0); ctx.lineTo(0, 13); ctx.lineTo(-9, 0); ctx.closePath(); ctx.fill();
    }
    ctx.restore(); ctx.globalAlpha = 1;
  }
  function ray(deg, col, dash, nm, lr) {
    const bx = MX, by = MY + MR * .86, a = (deg - 90) * Math.PI / 180;
    ctx.strokeStyle = col; ctx.lineWidth = dash ? 2 : 3;
    if (dash) ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(bx, by);
    ctx.lineTo(bx + Math.cos(a) * MR * 1.8, by + Math.sin(a) * MR * 1.8); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = col; ctx.font = "700 12.5px ui-monospace,Consolas,monospace";
    ctx.textAlign = "center";
    ctx.fillText(nm, bx + Math.cos(a) * MR * lr, by + Math.sin(a) * MR * lr - 6);
  }
  function frame(k) {
    ctx.clearRect(0, 0, w, h);
    ctx.strokeStyle = "rgba(140,180,210,.20)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
    ctx.beginPath(); ctx.arc(cx, cy, R * .72, 0, 7); ctx.stroke();
    ctx.strokeStyle = "#8ba3b8"; ctx.setLineDash([5, 5]); ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx, cy - R * 1.08); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "700 13px sans-serif";
    ctx.textAlign = "center"; ctx.fillText("진북", cx, cy - R * 1.17);
    needle(D90, "#6b8296", .8, false);                /* 시작 해 자취 */
    const cur = dAt(k);
    needle(cur, "#3fd8d0", 1, false);                 /* 지금 */
    ctx.beginPath(); ctx.arc(cx, cy, 5, 0, 7); ctx.fillStyle = "#eaf2f8"; ctx.fill();
    ctx.fillStyle = "#3fd8d0";
    ctx.font = "800 32px 'Arial Narrow',Consolas,sans-serif"; ctx.textAlign = "center";
    ctx.fillText(String(Math.round(DR.y0 + (DR.y1 - DR.y0) * k)), cx, cy + R * .5);
    ctx.font = "600 12.5px sans-serif"; ctx.fillStyle = "#8ba3b8";
    ctx.fillText(cur.toFixed(1) + "°", cx, cy + R * .5 + 19);
    /* 돋보기 — 실제 각도차가 작아 8배로 키운다 */
    const tx = cx + Math.sin(cur * Math.PI / 180) * R * .95;
    const ty = cy - Math.cos(cur * Math.PI / 180) * R * .95;
    ctx.strokeStyle = "rgba(140,180,210,.26)"; ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(tx + 7, ty - 5); ctx.lineTo(MX - MR * .95, MY - MR * .55);
    ctx.moveTo(tx + 7, ty + 9); ctx.lineTo(MX - MR * .95, MY + MR * .66); ctx.stroke();
    ctx.save();
    ctx.beginPath(); ctx.arc(MX, MY, MR, 0, 7); ctx.clip();
    ctx.fillStyle = "rgba(9,17,28,.9)";
    ctx.fillRect(MX - MR, MY - MR, MR * 2, MR * 2);
    ray((D90 - D26) * 8, "#6b8296", false, String(DR.y0), 1.64);
    ray(0, "#3fd8d0", false, String(DR.y1), 1.3);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "600 10.5px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("각도차만 키운 것", MX, MY + MR * .58);
    ctx.fillText("숫자는 실제값", MX, MY + MR * .78);
    ctx.restore();
    ctx.strokeStyle = "rgba(140,180,210,.45)"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(MX, MY, MR, 0, 7); ctx.stroke();
    ctx.fillStyle = "#55707f"; ctx.font = "700 11px ui-monospace,Consolas,monospace";
    ctx.textAlign = "center";
    ctx.fillText("각도차만 8배로 벌려 그림", MX, MY + MR + 17);
    ctx.fillStyle = "#eaf2f8"; ctx.font = "700 14.5px sans-serif";
    ctx.fillText(DR.site + " · " + (DR.y1 - DR.y0) + "년 동안 "
                 + DR.delta_min + "분 (1 km 앞에서 " + DR.offset_m + " m)",
                 w / 2, h - 32);
    ctx.fillStyle = "#8ba3b8"; ctx.font = "600 12.5px sans-serif";
    ctx.fillText("IGRF-14 계산값 · 변화 속도는 지역과 시기에 따라 다릅니다 ("
                 + DR.rate_lo + "~" + DR.rate_hi + "분/년)", w / 2, h - 12);
  }
  let t0 = null;
  function run(ts) {
    if (t0 === null) t0 = ts;
    const k = Math.min((ts - t0) / 2600, 1);
    frame(k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2);
    if (k < 1) requestAnimationFrame(run);
  }
  frame(0);
  new IntersectionObserver((es, ob) => {
    es.forEach(e => { if (e.isIntersecting) { requestAnimationFrame(run); ob.disconnect(); } });
  }, {threshold: .3}).observe(cv);
})();

/* ── 벡터 분해 (2D · 화살표마다 이름) ────────────── */
/* ⚠️ 1차본은 WebGL 로 돌아가는 화살표 셋이었는데 «무엇이 무엇인지»가
   없었다. 회전보다 이름이 중요하다 — 정지 그림에 이름표를 단다. */
(function () {
  const cv = document.getElementById("figVector"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  const O = [w * .33, h * .36], L = Math.min(w, h) * .45;
  const Dg = P.dec_now, Ig = 53;
  const D = Dg * Math.PI / 180, I = Ig * Math.PI / 180;
  /* 등각 투영 — 북은 오른쪽 위, 동은 오른쪽 아래, 연직은 곧게 아래 */
  const pr = (n, e, z) => [O[0] + (n * .87 + e * .87) * L,
                           O[1] + (-n * .5 + e * .5 + z) * L];
  const Ft = pr(Math.cos(I) * Math.cos(D), Math.cos(I) * Math.sin(D), Math.sin(I));
  const Ht = pr(Math.cos(I) * Math.cos(D), Math.cos(I) * Math.sin(D), 0);
  const Zt = pr(0, 0, Math.sin(I));
  const Nt = pr(.85, 0, 0), Et = pr(0, .85, 0);
  ctx.clearRect(0, 0, w, h);
  /* 지면 격자 */
  ctx.strokeStyle = "rgba(140,180,210,.16)"; ctx.lineWidth = 1;
  for (let k = -5; k <= 5; k++) {
    let a = pr(k / 10, -.5, 0), b = pr(k / 10, .5, 0);
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
    a = pr(-.5, k / 10, 0); b = pr(.5, k / 10, 0);
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
  }
  function arrow(p, col, lw) {
    ctx.strokeStyle = col; ctx.lineWidth = lw; ctx.lineCap = "round";
    ctx.beginPath(); ctx.moveTo(O[0], O[1]); ctx.lineTo(p[0], p[1]); ctx.stroke();
    const a = Math.atan2(p[1] - O[1], p[0] - O[0]), t = 16;
    ctx.fillStyle = col; ctx.beginPath(); ctx.moveTo(p[0], p[1]);
    ctx.lineTo(p[0] - Math.cos(a - .36) * t, p[1] - Math.sin(a - .36) * t);
    ctx.lineTo(p[0] - Math.cos(a + .36) * t, p[1] - Math.sin(a + .36) * t);
    ctx.closePath(); ctx.fill();
  }
  arrow(Nt, "#55707f", 1.6); arrow(Et, "#55707f", 1.6);
  /* 각도 호 */
  function arc(fn, col) {
    ctx.strokeStyle = col; ctx.lineWidth = 2;
    ctx.beginPath();
    for (let t = 0; t <= 1.0001; t += .04) {
      const p = fn(t);
      if (t === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    }
    ctx.stroke();
  }
  arc(t => pr(Math.cos(D * t) * .55, Math.sin(D * t) * .55, 0), "#8ba3b8");
  arc(t => pr(Math.cos(I * t) * Math.cos(D) * .7, Math.cos(I * t) * Math.sin(D) * .7,
              Math.sin(I * t) * .7), "#8ba3b8");
  /* 보조선 */
  ctx.strokeStyle = "rgba(140,180,210,.4)"; ctx.setLineDash([5, 5]); ctx.lineWidth = 1.2;
  ctx.beginPath(); ctx.moveTo(Ht[0], Ht[1]); ctx.lineTo(Ft[0], Ft[1]);
  ctx.moveTo(Zt[0], Zt[1]); ctx.lineTo(Ft[0], Ft[1]); ctx.stroke();
  ctx.setLineDash([]);
  arrow(Zt, "#9b7fe8", 3.2); arrow(Ht, "#4a7fe8", 3.2); arrow(Ft, "#3fd8d0", 4);
  ctx.beginPath(); ctx.arc(O[0], O[1], 5, 0, 7); ctx.fillStyle = "#eaf2f8"; ctx.fill();
  /* 이름표 — 오른쪽 한 줄로 모으고 지시선을 뺀다 */
  const LX = w * .70;
  const dm = pr(Math.cos(D * .5) * .55, Math.sin(D * .5) * .55, 0);
  const im = pr(Math.cos(I * .5) * Math.cos(D) * .7, Math.cos(I * .5) * Math.sin(D) * .7,
                Math.sin(I * .5) * .7);
  [[Nt, 80, "#8ba3b8", "북 (진북)"],
   [Ht, 126, "#4a7fe8", "H  수평성분"],
   [dm, 172, "#8ba3b8", "D  편각 — 북과 H 사이"],
   [im, 218, "#8ba3b8", "I  복각 — H와 F 사이"],
   [Et, 300, "#8ba3b8", "동"],
   [Ft, 346, "#3fd8d0", "F  총자력 — 화살표 전체"],
   [Zt, 400, "#9b7fe8", "Z  연직성분 (아래로)"]].forEach(([p, y, col, t]) => {
    ctx.strokeStyle = "rgba(140,180,210,.35)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(p[0] + 6, p[1]); ctx.lineTo(LX - 10, y - 4); ctx.stroke();
    ctx.fillStyle = col; ctx.font = "700 14px sans-serif"; ctx.textAlign = "left";
    ctx.fillText(t, LX, y);
  });
  ctx.fillStyle = "#55707f"; ctx.font = "600 12.5px sans-serif"; ctx.textAlign = "left";
  ctx.fillText("계산값 — " + P.dec_site + " · D "
               + Math.abs(P.dec_now).toFixed(1) + "° 서편 · I 약 53° 아래", 24, h - 18);
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
  ctx.fillStyle = "#55707f"; ctx.font = "600 11.5px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText("설명용으로 줄인 값입니다 — 실제 판독은 약 5만 nT 입니다", L, h - 8);
})();

/* ── 재현성 비교 ─────────────────────────────────── */
(function () {
  const cv = document.getElementById("figAudit"); if (!cv) return;
  const {ctx, w, h} = hidpi(cv);
  ctx.clearRect(0, 0, w, h);
  const L = 340, R = w - 150, MAX = 40;
  const bw = v => (v / MAX) * (R - L);
  [["같은 날 두 번 재면", 1.4, "#3fd8a0", "이만큼 잘 맞습니다 (산포 중앙값)"],
   ["몇 년 뒤 다시 가서 재면", 33.7, "#e8503f", "이만큼 어긋났습니다 (잔여 RMS)"]
  ].forEach((r, i) => {
    const y = 74 + i * 80;
    ctx.fillStyle = "#eaf2f8"; ctx.font = "700 16px sans-serif"; ctx.textAlign = "right";
    ctx.fillText(r[0], L - 22, y + 8);
    ctx.fillStyle = "rgba(255,255,255,.05)"; ctx.fillRect(L, y - 14, R - L, 40);
    ctx.fillStyle = r[2]; ctx.fillRect(L, y - 14, Math.max(bw(r[1]), 3), 40);
    ctx.textAlign = "left"; ctx.fillStyle = r[2];
    ctx.font = "800 22px 'Arial Narrow',Consolas,sans-serif";
    ctx.fillText(r[1].toFixed(1) + "분", L + bw(r[1]) + 14, y + 6);
    ctx.font = "600 13px sans-serif"; ctx.fillStyle = "#8ba3b8";
    ctx.fillText(r[3], L + bw(r[1]) + 14, y + 26);
  });
  ctx.fillStyle = "#8ba3b8"; ctx.font = "600 13.5px sans-serif"; ctx.textAlign = "left";
  ctx.fillText("달라진 것은 «무엇을 보고 방위를 잡았는가» 하나였습니다.", 42, h - 34);
  ctx.fillStyle = "#55707f"; ctx.font = "600 12px sans-serif";
  ctx.fillText("두 값은 재는 방식이 다릅니다(산포 중앙값 ↔ 잔여 RMS) — 자릿수를 견주는 그림입니다.",
               42, h - 14);
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
/* ── Cinematic presentation controls ───────────────────── */
(function(){
 const title=document.getElementById('hudCurrent');
 const flash=document.getElementById('sectionFlash');
 const sections=[...document.querySelectorAll('main>section')];
 let current='';
 const io=new IntersectionObserver(entries=>entries.forEach(entry=>{
   if(!entry.isIntersecting)return;
   const next=entry.target.dataset.t||'현장교육';
   if(next!==current){
     current=next;title.textContent=next;
     flash.classList.remove('go');void flash.offsetWidth;flash.classList.add('go');
   }
 }),{threshold:.55});
 sections.forEach(s=>io.observe(s));
 document.getElementById('fullScreen').addEventListener('click',async()=>{
   if(!document.fullscreenElement)await document.documentElement.requestFullscreen?.();
   else await document.exitFullscreen?.();
 });
 document.getElementById('hudToggle').addEventListener('click',()=>document.body.classList.toggle('hud-hidden'));
 addEventListener('keydown',async e=>{
   /* ⚠️ Ctrl+F(찾기)·Cmd+H 까지 삼키면 안 된다 */
   if(e.ctrlKey||e.metaKey||e.altKey||!e.key)return;
   if(e.key.toLowerCase()==='h')document.body.classList.toggle('hud-hidden');
   if(e.key.toLowerCase()==='f'){
     if(!document.fullscreenElement)await document.documentElement.requestFullscreen?.();
     else await document.exitFullscreen?.();
   }
 });
 /* ⚠️ 포인터마다 CSS 변수를 쓰면 배경 그라디언트가 매번 다시 칠해진다.
    빔프로젝터를 물린 노트북에서 눈에 띄게 끊긴다 — 프레임당 한 번으로 묶는다. */
 let px=0,py=0,queued=false;
 addEventListener('pointermove',e=>{
   px=e.clientX/innerWidth*100;py=e.clientY/innerHeight*100;
   if(queued)return;queued=true;
   requestAnimationFrame(()=>{queued=false;
     document.body.style.setProperty('--mx',px.toFixed(2)+'%');
     document.body.style.setProperty('--my',py.toFixed(2)+'%');});
 },{passive:true});
})();

/* ── 장면을 한 화면에 맞춘다 ───────────────────────────── */
/* ⚠️ transform 은 «레이아웃 상자»를 바꾸지 않는다. 축소만 하면 section 은
   원래 높이를 그대로 차지해 빈 공간이 남는다 — 넘치는 장면만 height 를
   화면 높이로 못박고 overflow:hidden 이 축소된 결과를 담게 한다. */
(function () {
  const secs = [...document.querySelectorAll("main>section")];
  let raf = 0;
  /* 장면이 화면에 들어올 때(=캔버스가 그려질 때) 그 장면을 다시 잰다 */
  const seen = new IntersectionObserver(es => {
    if (es.some(e => e.isIntersecting)) schedule();
  }, {threshold: .2});
  function fit() {
    const vh = innerHeight, wide = innerWidth > 860;
    secs.forEach(s => {
      const w = s.firstElementChild; if (!w) return;
      s.style.height = ""; w.style.removeProperty("--fit");
      if (!wide) return;
      const cs = getComputedStyle(s);
      const pad = parseFloat(cs.paddingTop) + parseFloat(cs.paddingBottom);
      const need = w.getBoundingClientRect().height + pad;
      if (need > vh - 6) {
        w.style.setProperty("--fit", Math.max(.72, (vh - 6) / need).toFixed(4));
        s.style.height = vh + "px";
      }
    });
  }
  /* ⚠️ requestAnimationFrame 으로 예약하면 «보이지 않는 탭»에서 멈춘다.
     배경 탭으로 열어 두었다가 넘어오면 장면이 안 맞은 채로 뜬다 —
     setTimeout 은 숨어 있어도 돈다. */
  function schedule() { clearTimeout(raf); raf = setTimeout(fit, 40); }
  schedule();
  addEventListener("resize", schedule);
  addEventListener("load", () => setTimeout(schedule, 60));
  addEventListener("visibilitychange", schedule);
  /* 도판이 그려지고 글꼴이 앉은 뒤 다시 — 캔버스는 화면에 들어와야 크기가 잡힌다 */
  [400, 900, 1800].forEach(t => setTimeout(schedule, t));
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(schedule);
  secs.forEach(x => seen.observe(x));
  window.__fit = fit;
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
        '단위 nT/m · 중심점에서 1 m 떨어진 지점의 «변화율». 부호는 방향이라 '
        '「−」는 중심점보다 값이 작아졌다는 뜻이고, 판정에 쓰는 것은 '
        '부호를 뗀 크기입니다.</p>'
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
        "OFFSETS_TXT": ", ".join(str(x) for x in payload["offsets"]),
        "H_MIN": payload["heights"][0],
        "H_MAX": payload["heights"][-1],
        "H_STEP": payload["heights"][1] - payload["heights"][0],
        "N_READ": payload["n_read"],
        "N_READ_H": payload["n_read_h"],
        "N_READ_V": payload["n_read_v"],
        "N_HEIGHT": len(payload["heights"]),
        "N_OFF": len(payload["offsets"]),
        "DEC_DMS": payload["dec_dms"],
        "DEC_LO": f"{payload['dec_lo']:.0f}",
        "DEC_HI": f"{payload['dec_hi']:.0f}",
        "IAGA_RADIUS": f"{payload['iaga_radius']:.0f}",
        "IAGA_RANGE": f"{payload['iaga_range']:.0f}",
        "EURO_GRAD": f"{payload['euro_grad']:.0f}",
        "WORST_VAL": payload["worst_val"],
        "WORST_RATIO": payload["worst_ratio"],
        "PILOT_TABLE": pilot_table(payload),
    }
    for k, v in (("SVG_GEAR", SVG_GEAR), ("SVG_CLEAN", SVG_CLEAN),
                 ("SVG_USE", SVG_USE)):
        html = html.replace("{{" + k + "}}", v)
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
