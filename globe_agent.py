# -*- coding: utf-8 -*-
"""
지자기 3D 지구본 LUI 에이전트 — 로컬 LLM(LM Studio) 질의 응답.
========================================================================

2026_ex_simul 의 LUI 패턴(서버가 데이터·수치를 계산해 컨텍스트로 주입 → LLM 은
서술만) 을 따른다. 소형 로컬 모델의 tool-calling 불안정성을 피하려는 설계.

  · ① 자연어 지자기 계산기 — ppigrf 로 임의 지점·연도 F/D/I/H 정확 계산(+영년변화)
  · ③ 선점 어드바이저 — survey_sites/existing24 검색·집계

반환: {"answer": str, "actions": [{type,...}], "tier": "local|fallback|data-only"}
actions: 프론트가 지구본을 조종 — focus(lat,lon)·highlight(mids)·open_detail

LM Studio 미실행이면 서버 계산 결과(수치)만 담아 답한다(여전히 유용).
"""
from __future__ import annotations

import datetime as dt
import json
import math
import os
import re
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"

# .env(루트) 자동 로드 — python-dotenv 없이 동작하는 간이 파서(KEY=VALUE)
def _load_env():
    p = ROOT / ".env"
    if not p.exists():
        return
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#") or "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        k, v = k.strip(), v.strip().strip('"').strip("'")
        if k and k not in os.environ:      # 이미 설정된 환경변수 우선
            os.environ[k] = v


_load_env()

LM_STUDIO_URL = os.getenv("LM_STUDIO_URL", "http://localhost:1234/v1")
LOCAL_MODEL = os.getenv("LOCAL_LLM_MODEL", "").strip()   # 빈 값 = 로드된 모델
TIMEOUT_S = int(os.getenv("LOCAL_LLM_TIMEOUT", "120"))
NOW_YEAR = 2026.6

# ── LLM 백엔드 선택 ────────────────────────────────────────────────────────────
#   GLOBE_LLM_BACKEND = auto(기본) | lmstudio | openai | gemini
#   auto: LM Studio 먼저 → 실패하면 Gemini(키 있을 때) → OpenAI(키 있을 때) → 서버계산.
BACKEND = os.getenv("GLOBE_LLM_BACKEND", "auto").lower().strip()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "").strip()
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o-mini").strip()
OPENAI_BASE_URL = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1").strip()
# Gemini — OpenAI 호환 엔드포인트(https://ai.google.dev/gemini-api/docs/openai)라
# 별도 SDK 없이 openai 클라이언트를 base_url 만 바꿔 그대로 재사용한다.
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY", "").strip()
GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.0-flash").strip()
GEMINI_BASE_URL = os.getenv("GEMINI_BASE_URL",
                             "https://generativelanguage.googleapis.com/v1beta/openai/").strip()

# ── 게이지티어(주요 지점 + 우리 측점 좌표는 데이터에서 로드) ──────────────────
GAZETTEER = {
    "서울": (37.5665, 126.9780), "부산": (35.1796, 129.0756), "인천": (37.4563, 126.7052),
    "대구": (35.8714, 128.6014), "대전": (36.3504, 127.3845), "광주": (35.1595, 126.8526),
    "울산": (35.5384, 129.3114), "제주": (33.4996, 126.5312), "독도": (37.2429, 131.8664),
    "울릉도": (37.4845, 130.9057), "백령도": (37.9660, 124.6100), "마라도": (33.1216, 126.2683),
    "강릉": (37.7519, 128.8761), "춘천": (37.8813, 127.7300), "청양": (36.4595, 126.8020),
    "포항": (36.0190, 129.3435), "여수": (34.7604, 127.6622), "속초": (38.2070, 128.5918),
    "남대서양": (-25.0, -45.0), "남대서양이상": (-25.0, -45.0),
}


def _load(name):
    p = DATA / name
    if p.exists():
        return json.load(open(p, encoding="utf-8")).get("features", [])
    return []


SITES = _load("survey_sites.geojson")     # 선점 후보 103
EXIST = _load("existing24.geojson")       # 기존점 24
# 측점 이름 → 좌표 게이지티어 보강
for f in SITES + EXIST:
    nm = f["properties"].get("name")
    c = f["geometry"]["coordinates"]
    if nm:
        GAZETTEER.setdefault(nm, (c[1], c[0]))


# ── IGRF (ppigrf) ─────────────────────────────────────────────────────────────
def igrf_elements(lat, lon, year, elev_m=0):
    import ppigrf
    date = dt.datetime(int(year), 1, 1) + dt.timedelta(days=(year - int(year)) * 365.25)
    Be, Bn, Bu = ppigrf.igrf(lon, lat, elev_m / 1000.0, date)
    Be, Bn, Bu = float(Be), float(Bn), float(Bu)
    X, Y, Z = Bn, Be, -Bu
    H = math.hypot(X, Y)
    F = math.hypot(H, Z)
    D = math.degrees(math.atan2(Y, X))
    I = math.degrees(math.atan2(Z, H))
    return {"F": F, "D": D, "I": I, "H": H, "X": X, "Y": Y, "Z": Z}


def dms(d):
    s = "-" if d < 0 else ""
    a = abs(d)
    deg = int(a); mn = int((a - deg) * 60); sc = ((a - deg) * 60 - mn) * 60
    return f"{s}{deg}°{mn}'{sc:.1f}\""


# ── 의도 파싱 ─────────────────────────────────────────────────────────────────
HQ_KEYS = ["강원", "경남", "경북", "대구경북", "충북", "충남", "대전세종충남", "전북",
           "전남", "광주전남", "제주", "부산", "울산", "부산울산", "서울", "경기",
           "인천", "서울경기북부", "인천경기남부"]


def find_places(msg):
    """메시지에서 게이지티어 지명 탐지(긴 이름 우선)."""
    hits = []
    for name in sorted(GAZETTEER, key=len, reverse=True):
        if name in msg:
            hits.append(name)
    # 부분 중복 제거(짧은 게 긴 이름 안에 들어간 경우)
    out = []
    for h in hits:
        if not any(h != o and h in o for o in hits):
            out.append(h)
    return out[:4]


def find_years(msg):
    ys = [int(y) for y in re.findall(r"20\d{2}", msg)]
    return [y for y in ys if 2015 <= y <= 2030]


def query_sites(msg):
    """등급·본부·방위표지·판정 키워드로 선점 후보 검색."""
    grades = []
    for g in ("A", "B", "C", "D"):
        if re.search(rf"(등급\s*{g}|{g}\s*등급|\b{g}급)", msg):
            grades.append(g)
    if "적합" in msg and "부적합" not in msg and "조건부" not in msg and not grades:
        grades.append("A")
    if "부적합" in msg and not grades:
        grades.append("D")
    hq = next((h for h in HQ_KEYS if h in msg), None)
    want_bang = "방위표지" in msg
    res = []
    for f in SITES:
        p = f["properties"]
        if grades and p["grade"] not in grades:
            continue
        if hq and hq not in (p.get("hq") or ""):
            continue
        res.append(p | {"lat": f["geometry"]["coordinates"][1],
                        "lon": f["geometry"]["coordinates"][0]})
    return res, {"grades": grades, "hq": hq, "want_bang": want_bang}


# ── 컨텍스트 빌드 ──────────────────────────────────────────────────────────────
def build_context(msg):
    """질의 → (데이터 컨텍스트 문자열, actions). 계산은 여기서 정확히 수행."""
    ctx, actions = [], []
    years = find_years(msg)
    yr = float(years[0]) if years else NOW_YEAR

    # 선점 신호는 명시적 키워드로만(지명·본부명 city 와 혼동 방지)
    grades_pre = [g for g in ("A", "B", "C", "D")
                  if re.search(rf"(등급\s*{g}|{g}\s*등급|\b{g}\s*급)", msg)]
    site_kw = bool(grades_pre) or any(k in msg for k in
                                      ("후보지", "후보", "선점", "부적합", "방위표지", "대체", "검토 등급")) \
        or "본부" in msg

    # ① 지자기 계산(지명 탐지) — 선점 질의가 아닐 때
    places = find_places(msg)
    if places and not site_kw:
        for nm in places[:3]:
            lat, lon = GAZETTEER[nm]
            e = igrf_elements(lat, lon, yr)
            ctx.append(f"[{nm} {yr:.0f} IGRF-14]  위도 {lat:.4f}·경도 {lon:.4f}  "
                       f"편각 D {e['D']:.3f}° ({dms(e['D'])}) · 복각 I {e['I']:.3f}° · "
                       f"총자력 F {e['F']:.0f} nT · 수평 H {e['H']:.0f} nT")
            actions.append({"type": "focus", "lat": lat, "lon": lon, "label": nm})
            # 영년변화(두 연도)
            if len(years) >= 2:
                e2 = igrf_elements(lat, lon, float(years[1]))
                dD = (e2["D"] - e["D"]) * 60
                ctx.append(f"[{nm} 영년변화 {years[0]}→{years[1]}]  편각 ΔD {dD:+.1f}′ "
                           f"({dD/60:+.3f}°) · 총자력 ΔF {e2['F']-e['F']:+.0f} nT · "
                           f"복각 ΔI {(e2['I']-e['I'])*60:+.1f}′")

    # ③ 선점 어드바이저
    if site_kw:
        res, flt = query_sites(msg)
        gc = {"A": 0, "B": 0, "C": 0, "D": 0}
        for r in res:
            gc[r["grade"]] = gc.get(r["grade"], 0) + 1
        head = (f"[선점 검색] 조건 {flt} → {len(res)}건 "
                f"(A{gc['A']}·B{gc['B']}·C{gc['C']}·D{gc['D']})")
        lines = [head]
        for r in res[:12]:
            lines.append(f"  {r['mid']} {r['name']}({r.get('hq','')}) 등급 {r['grade']}·"
                         f"{r['verdict']}·방위표지 {r.get('bang','-')}·교란 {r.get('disturb','-')}")
        if len(res) > 12:
            lines.append(f"  …외 {len(res)-12}건")
        ctx.append("\n".join(lines))
        if res:
            actions.append({"type": "highlight", "mids": [r["mid"] for r in res]})
            actions.append({"type": "focus_korea"})

    return "\n".join(ctx), actions


# ── LLM 호출(LM Studio) / 폴백 ────────────────────────────────────────────────
SYSTEM = """당신은 '지자기 3D 지구본'의 분석 보조 AI다. 한국 국가기준점 지자기측량 프로젝트를 돕는다.
아래 [데이터]에 서버가 계산한 정확한 수치(IGRF-14 편각 D·복각 I·총자력 F, 선점 후보 등급 등)가 주어진다.
규칙:
- [데이터]의 수치는 한 글자도 바꾸지 말고 그대로 인용한다. 데이터에 없으면 추측하지 않는다.
- 한국어로 간결·정확하게. 표보다 문장 위주. 핵심 결론을 먼저.
- 편각(D)=진북 대비 자침 편차, 복각(I)=수평면 대비 자기장 기울기, F=총자력. 영년변화는 시간에 따른 변화.
- 선점 등급(4단계): A 선점가능(자기 청정+방위표지≥100m) / B 조건부 선점가능(자기 청정이나 방위표지 거리<100m — 정밀 방위각·표지 이설 시 선점가능) / C 현장 확인 필요(자기교란 재확인) / D 부적합(방위표지 불가 포함, 대체 검토).
- 데이터가 비어 있으면(일반 개념 질문) 프로젝트 지식으로 간결히 설명한다."""


def _lmstudio(messages):
    from openai import OpenAI
    client = OpenAI(base_url=LM_STUDIO_URL, api_key="lm-studio",
                    timeout=TIMEOUT_S, max_retries=0)   # 미실행 시 즉시 폴백
    r = client.chat.completions.create(model=LOCAL_MODEL or "local-model",
                                       messages=messages, temperature=0.3, max_tokens=1024)
    return r.choices[0].message.content or ""


def _openai(messages):
    """진짜 OpenAI API(또는 호환 base_url). OPENAI_API_KEY 필요. 키는 서버(.env)에만 둔다."""
    if not OPENAI_API_KEY:
        raise RuntimeError("OPENAI_API_KEY 미설정")
    from openai import OpenAI
    client = OpenAI(api_key=OPENAI_API_KEY, base_url=OPENAI_BASE_URL, timeout=60, max_retries=1)
    r = client.chat.completions.create(model=OPENAI_MODEL, messages=messages,
                                       temperature=0.3, max_tokens=1024)
    return r.choices[0].message.content or ""


def _gemini(messages):
    """Gemini(OpenAI 호환 엔드포인트). GEMINI_API_KEY 필요. 키는 서버(.env)에만 둔다."""
    if not GEMINI_API_KEY:
        raise RuntimeError("GEMINI_API_KEY 미설정")
    from openai import OpenAI
    client = OpenAI(api_key=GEMINI_API_KEY, base_url=GEMINI_BASE_URL, timeout=60, max_retries=1)
    r = client.chat.completions.create(model=GEMINI_MODEL, messages=messages,
                                       temperature=0.3, max_tokens=1024)
    return r.choices[0].message.content or ""


def _tiers():
    """시도 순서 [(tier명, 호출함수)] — 백엔드 설정에 따라."""
    lm = ("local", _lmstudio)
    ge = ("gemini", _gemini)
    oa = ("openai", _openai)
    if BACKEND == "gemini":
        return [ge] if GEMINI_API_KEY else []
    if BACKEND == "openai":
        return [oa] if OPENAI_API_KEY else []
    if BACKEND == "lmstudio":
        return [lm]
    # auto: 로컬 먼저, 그다음 Gemini, 그다음 OpenAI(키 있는 것만)
    return [lm] + ([ge] if GEMINI_API_KEY else []) + ([oa] if OPENAI_API_KEY else [])


def answer(message):
    message = (message or "").strip()
    if not message:
        return {"answer": "질문을 입력해 주세요.", "actions": [], "tier": "data-only"}
    try:
        ctx, actions = build_context(message)
    except Exception as e:            # 계산 오류도 치명적이지 않게
        ctx, actions = f"(계산 오류: {e})", []

    user = (f"[데이터]\n{ctx}\n\n[질문]\n{message}\n\n위 데이터를 근거로 답해줘."
            if ctx else f"[질문]\n{message}")
    messages = [{"role": "system", "content": SYSTEM}, {"role": "user", "content": user}]
    last = ""
    for tier, fn in _tiers():
        try:
            ans = fn(messages)
            if ans.strip():
                return {"answer": ans, "actions": actions, "tier": tier}
        except Exception as e:            # noqa: BLE001
            last = str(e)
    # 모든 LLM 실패/미설정 → 서버 계산 결과만이라도 반환
    note = ("⚠ LLM 미응답 — 서버 계산 결과만 표시합니다. LM Studio(localhost:1234) 실행 또는 "
            ".env 에 GEMINI_API_KEY/OPENAI_API_KEY 설정 시 자연어 설명이 붙습니다.\n\n")
    body = ctx if ctx else f"(질의 '{message}' 에 대한 계산 데이터가 없습니다. 지명·등급 등을 포함해 주세요.)"
    return {"answer": note + body, "actions": actions, "tier": "data-only", "err": last}


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    for q in ["2027년 독도 편각과 복각 알려줘",
              "서울 편각이 2025에서 2030 얼마나 변해?",
              "강원 A등급 후보지 알려줘",
              "부적합 후보지는?"]:
        r = answer(q)
        print("Q:", q)
        print("  tier:", r["tier"], " actions:", [a["type"] for a in r["actions"]])
        print("  answer:", r["answer"][:200].replace("\n", " "))
        print()
