# -*- coding: utf-8 -*-
"""
2020 연구의 15점 ↔ 2022~25 측량 15점 ↔ 현 LMM 16점 — 대조 보고
==============================================================

명단·좌표·배치는 `legacy2020_sets` 에서 읽는다(발표자료와 같은 출처).
여기서는 그 위에 **2020 보고서 자체 수치의 재현 점검**을 얹는다.

    python _compare_2020_vs_2225.py
"""
import datetime as dt
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

import legacy2020_sets as LS                                   # noqa: E402

# ── 2020 보고서 본문의 연변화 서술값 (p.121~127) ────────────────────────────
#    (IGRF 연변화 ′/년, 관측 연평균 편각변화 ′/년, 원문 표기)
SV_2020 = {
    "화천": (5 + 0 / 60, 8 + 57.9 / 60, "56′48.9″/6년"),
    "제천": (8 + 29 / 60, 2 + 50.5 / 60, "34′32.8″/3년"),
    "춘양": (4 + 53 / 60, 7 + 16.3 / 60, "21′49″/3년"),
    "서산": (5 + 12 / 60, 2 + 38.7 / 60, "13′2.7″/5년"),
    "청양": (5 + 9 / 60, 6 + 45 / 60, "33′45.2″/5년"),
    "이원": (5 + 4 / 60, 0 + 44.4 / 60, "3′37″/5년"),
    "와도": (5 + 15 / 60, 1 + 19.9 / 60, "11′14.9″/6년"),
    "남지": (4 + 59 / 60, 4 + 42.3 / 60, "18′49.2″/4년"),
    "거제": (4 + 59 / 60, 3 + 46.8 / 60, "18′53.9″/5년"),
    "장흥": (5 + 12 / 60, 9 + 40.5 / 60, "48′22.6″/5년"),
    "성산": (5 + 14 / 60, 4 + 32.7 / 60, "22′43.3″/5년"),
}

# ── 2020 보고서 <표4-1> KRISS 절대측정 5점 ──────────────────────────────────
#    (측정일, D관측, I관측, F관측, D보고IGRF, I보고IGRF, F보고IGRF)
KRISS = {
    "가야": ("2020-09-15", (-8, 27, 4), (46, 45, 52), 49045.27,
            (-7, 25, 18), (46, 51, 42), 49843.06),
    "남해": ("2020-07-02", (-7, 35, 54), (45, 54, 58), 49054.63,
            (-7, 9, 45), (45, 48, 49), 49372.67),
    "부안": ("2020-11-03", (-7, 24, 7), (46, 57, 4), 50530.20,
            (-7, 15, 51), (46, 53, 58), 50107.07),
    "설악": ("2020-10-22", (-9, 8, 0), (54, 29, 36), 50513.50,
            (-8, 54, 44), (54, 46, 23), 51023.00),
    "순창": ("2020-06-17", (-7, 47, 38), (46, 9, 52), 50009.93,
            (-7, 13, 16), (46, 29, 0), 49836.87),
}


def dms2deg(t):
    s = -1 if t[0] < 0 else 1
    return s * (abs(t[0]) + t[1] / 60 + t[2] / 3600)


def rms(v):
    return math.sqrt(sum(x * x for x in v) / len(v))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    warnings.filterwarnings("ignore")
    S = LS.sets()
    cur = S["current"]

    print("═" * 78)
    print("① 세 집합의 구성")
    print("═" * 78)
    print(f"2020 분석 15점  ({len(S['s2020']):>2}) {sorted(S['s2020'])}")
    print(f"’22~25 측량     ({len(S['s2225']):>2}) {sorted(S['s2225'])}")
    print(f"현 LMM 적합     ({len(cur):>2}) {sorted(cur)}")
    print()
    print(f"공통          ({len(S['common'])}) {sorted(S['common'])}")
    print(f"2020 에만     ({len(S['only2020'])}) {sorted(S['only2020'])}")
    print(f"’22~25 에만   ({len(S['only2225'])}) {sorted(S['only2225'])}")
    print(f"\n  → ’22~25 신규 {len(S['only2225'])}점 중 «2020 이 측량 1회로 배제한 점» : "
          f"{len(S['new_from_once'])}점")
    print(f"  → 2020 에서 빠진 {len(S['only2020'])}점 중 저수지·제방 설치 : "
          f"{len(S['dropped_pond'])}점  {sorted(S['dropped_pond'])}")

    print("\n■ 공간 커버리지 (국토 0.1° 격자 · plan_new_sites 와 동일 방법)")
    print(f"{'집합':<20}{'n':>4}{'중앙':>8}{'90%ile':>8}{'최대공백':>10}"
          f"{'최대공백 위치':>20}{'50km초과':>10}")
    for tag, names, pm in (("2020 분석 15점", S["s2020"], False),
                           ("’22~25 측량 15점", S["s2225"], True),
                           ("현 LMM", cur, True),
                           ("1등 지자기점 30점", set(map(LS.canon, S["pts30"])), False)):
        c = LS.coverage(LS.coords(names, prefer_model=pm))
        print(f"{tag:<20}{c['n']:>4}{c['median']:>7.0f}km{c['p90']:>7.0f}"
              f"{c['max']:>9.0f}km   {c['at_lat']:.2f}N {c['at_lon']:.2f}E"
              f"{c['over50']:>9.1f}%")

    print("\n" + "═" * 78)
    print("② 2020 보고서 연변화 서술 — 관측 vs IGRF (p.121~127)")
    print("═" * 78)
    print(f"{'측점':<5}{'IGRF ′/년':>10}{'관측 ′/년':>10}{'차이':>9}   "
          f"{'원문':<14} 산술재현")
    diff = []
    for k, (g, o, src) in SV_2020.items():
        n = [float(x) for x in re.findall(r"[\d.]+", src)]
        calc = (n[0] + n[1] / 60 if len(n) == 3 else n[0]) / n[-1]
        flag = "" if abs(calc - o) < 0.1 else f"  ← 본문 {o:.2f} ≠ {calc:.2f}"
        print(f"{k:<5}{g:>10.2f}{o:>10.2f}{o - g:>+9.2f}   {src:<14}{flag}")
        diff.append(o - g)
    print(f"\n  관측−IGRF 연변화 차이 RMS {rms(diff):.2f} ′/년 (n={len(diff)}) "
          f"→ 5년 재방문 환산 {rms(diff) * 5:.1f}′")

    print("\n" + "═" * 78)
    print("③ 2020 <표4-1> KRISS 절대측정 5점 — 보고서 IGRF 열의 재현 점검")
    print("═" * 78)
    from lmm_build import igrf_dif
    P = S["pts30"]
    print(f"{'측점':<5}{'D관측':>9}{'D보고IGRF':>11}{'D재계산':>10}"
          f"{'I관측':>9}{'I보고IGRF':>11}{'I재계산':>10}")
    r_rep, r_new = [], []
    for k, (dd, Do, Io, _Fo, Dg, Ig, _Fg) in KRISS.items():
        la, lo = P[k]["lat"], P[k]["lon"]
        y, m, day = (int(x) for x in dd.split("-"))
        D, I, _ = (float(v[0]) for v in
                   igrf_dif(np.array([la]), np.array([lo]), np.zeros(1),
                            dt.datetime(y, m, day))[:3])
        print(f"{k:<5}{dms2deg(Do):>9.3f}{dms2deg(Dg):>11.3f}{D:>10.3f}"
              f"{dms2deg(Io):>9.3f}{dms2deg(Ig):>11.3f}{I:>10.3f}")
        r_rep.append((dms2deg(Do) - dms2deg(Dg)) * 60)
        r_new.append((dms2deg(Do) - D) * 60)
    print(f"\n  관측−보고서IGRF    편각 RMS {rms(r_rep):.1f}′")
    print(f"  관측−IGRF-14재계산 편각 RMS {rms(r_new):.1f}′")
    print("  ⚠ 복각은 설악을 뺀 4점에서 보고서 IGRF 열이 재계산과 약 5° 어긋난다 "
          "— 표 원본 확인 필요")

    print("\n" + "═" * 78)
    print("④ 현 LMM 의 대응 수치 (lmm_model.json)")
    print("═" * 78)
    model = json.load(open(ROOT / "docs" / "data" / "lmm_model.json",
                           encoding="utf-8"))
    V = {(r["성분"], r["단계"]): r["RMS"] for r in model["validation"]}
    reg = model["regional"]
    print(f"  IGRF 단독 편각 잔차 RMS  {V[('D_deg', 'IGRF')] * 60:.1f}′")
    print(f"  +Regional 후            {V[('D_deg', '+Regional')] * 60:.1f}′")
    print(f"  Regional 상수항          D {reg['D'][0] * 60:+.2f}′ · "
          f"I {reg['I'][0] * 60:+.2f}′ · F {reg['F'][0]:+.2f} nT "
          f"(차수 {reg['degree']})")
    print(f"  Station-LOSO            D {model['loo_cv']['D']:.4f}° "
          f"({model['loo_cv']['D'] * 60:.1f}′) · F {model['loo_cv']['F']:.1f} nT")


if __name__ == "__main__":
    main()
