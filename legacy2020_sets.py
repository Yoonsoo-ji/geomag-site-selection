# -*- coding: utf-8 -*-
"""
2020 연구의 15점 ↔ 2022~25 측량 15점 ↔ 현 LMM 적합 측점 — 집합·배치 대조
========================================================================

발표자료·보고서가 같은 수치를 쓰도록 **단일 출처**로 모아 둔다. 어느 것도
하드코딩하지 않고 1차 자료에서 읽는다.

  ① 2020 분석 15점   2020 지구물리측량 연구보고서 p.117 이 이름으로 나열한
                     「2012~2019 편복각 측량회수 1회」 15점의 여집합.
                     ⚠ 보고서의 「1회」는 원시 횟수가 아니라 **유효 횟수**다 —
                       미원·강화는 2013 관측이 있으나 신뢰도 미확보, 남양·양산은
                       2016 성과 불량으로 각각 1회로 취급됐다.
  ② 2022~25 측량 15점 data/지자기측량 성과정리(22_25).xlsx 의 연번 1~15
  ③ 현 LMM 측점      docs/data/lmm_model.json 의 sites

배치 지표는 `plan_new_sites` 의 격자·거리 함수를 그대로 쓴다(두 산출물이
갈라지지 않게).
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"
DOCS = ROOT / "docs" / "data"

# 2020 보고서 p.117 — 측량회수 1회로 분류돼 분석에서 빠진 15점
ONCE_2020 = ("상주", "부안", "순창", "함양", "가야", "경주", "남해", "순천",
             "조도", "포천", "설악", "미원", "강화", "남양", "양산")

# 개명 — 같은 표석의 다른 이름
ALIAS = {"임계": "삼척", "경주": "영천"}

# 저수지·제방 설치 판별에 쓰는 소재지 키워드
_POND = re.compile(r"저수지|제방|뚝|못|지뚝")


def canon(name):
    return ALIAS.get(name, name)


def pts30():
    """1등 지자기점 30점 — {점명: {lat, lon, addr}}."""
    d = pd.read_excel(DATA / "'10~'19년 지자기점 관측현황(최종).xls", header=None)

    def dms(s):
        n = [float(x) for x in re.findall(r"[\d.]+", str(s))]
        return n[0] + n[1] / 60 + n[2] / 3600 if len(n) >= 3 else np.nan

    out = {}
    for i in range(3, 33):
        r = d.iloc[i]
        if pd.isna(r[2]):
            continue
        out[r[2]] = dict(lat=dms(r[7]), lon=dms(r[6]), addr=str(r[4]))
    return out


def survey_sites_2225():
    """2022~25 성과표의 측점 — 연번이 붙은 행만(언양 등 별칭 행 제외)."""
    d = pd.read_excel(DATA / "지자기측량 성과정리(22_25).xlsx", header=None)
    return [d.iloc[i, 1] for i in range(1, len(d))
            if not pd.isna(d.iloc[i, 0]) and not pd.isna(d.iloc[i, 1])]


def model_sites():
    """현 LMM 적합 측점 — {점명: (lat, lon)}."""
    m = json.load(open(DOCS / "lmm_model.json", encoding="utf-8"))
    return {s["name"]: (s["lat"], s["lon"]) for s in m["sites"]}


def sets():
    """세 집합과 그 교차 관계를 한 번에 낸다."""
    P = pts30()
    cur = model_sites()
    once = {canon(n) for n in ONCE_2020}
    s2020 = {canon(n) for n in P if n not in ONCE_2020}
    s2225 = {canon(n) for n in survey_sites_2225()}

    def pond(n):
        src = next((k for k in P if canon(k) == n), None)
        return bool(src and _POND.search(P[src]["addr"]))

    only20 = s2020 - s2225
    only22 = s2225 - s2020
    return dict(
        pts30=P, current=cur,
        s2020=s2020, s2225=s2225, once=once,
        common=s2020 & s2225, only2020=only20, only2225=only22,
        # 2022~25 신규 측점 가운데 2020 이 「자료 부족」으로 배제했던 점
        new_from_once=only22 & once,
        # 2020 집합에서 빠진 측점 가운데 저수지·제방 설치점
        dropped_pond={n for n in only20 if pond(n)},
    )


def coords(names, prefer_model=False):
    """점명 목록 → [(lat, lon)]. prefer_model 이면 현 성과 좌표를 우선한다."""
    P, cur = pts30(), model_sites()
    out = []
    for n in names:
        if prefer_model and n in cur:
            out.append(cur[n])
            continue
        src = next((k for k in P if canon(k) == n), None)
        if src:
            out.append((P[src]["lat"], P[src]["lon"]))
        elif n in cur:
            out.append(cur[n])
    return out


def coverage(pts):
    """국토 0.1° 격자에서의 최근접 측점 거리 통계 (plan_new_sites 와 동일 방법)."""
    import plan_new_sites as PN

    lats = np.arange(33.1, 38.7, 0.1)
    lons = np.arange(125.0, 129.7, 0.1)
    mask = PN.korea_mask(lats, lons)
    LO, LA = np.meshgrid(lons, lats)
    gl, gn = LO[mask], LA[mask]
    d = np.full(gl.size, np.inf)
    for la, lo in pts:
        d = np.minimum(d, PN.geo_km(gn, gl, la, lo))
    i = int(d.argmax())
    return dict(n=len(pts), median=float(np.median(d)),
                p90=float(np.percentile(d, 90)), max=float(d.max()),
                at_lat=float(gn[i]), at_lon=float(gl[i]),
                over50=float(100 * (d > 50).mean()))
