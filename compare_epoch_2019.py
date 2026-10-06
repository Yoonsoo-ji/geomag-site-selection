# -*- coding: utf-8 -*-
"""
2019 자료 포함 여부 — 에폭 혼용 검정 (station-LOSO)
=====================================================

자문 권고 「에폭 혼용 금지」의 검정. 우리 Regional 층은 **시간 항이 없는**
공간 다항식이므로(REGIONAL_TIME_TERM=False), 2019~2025 를 한 상수항에 섞으면
IGRF 가 못 잡는 영년변화 잔차가 평균되어 들어간다. 그 크기를 실측한다.

⚠ 공정한 비교의 조건 — 평가 측점을 **양쪽 공통 집합**으로 고정한다.
  2019 를 빼면 남지·장흥이 사라지므로, 그 두 점을 평가에 넣으면 비교가 성립하지 않는다.
  또 거제·이원은 2019 행이 있어 «관측 대표값 자체»가 바뀌므로 따로 분리해 본다.

    python compare_epoch_2019.py
"""
import importlib
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DEGREE = 0

import lmm_build as LB      # noqa: E402

_lons, _lats, _grid = LB.load_kigam_grid()
CRUSTAL = LB.CrustalGrid(_lons, _lats, _grid)


def build(inc2019):
    importlib.reload(LB)
    pts = LB.load_all_points(include_2019=inc2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    return LB.attach_crustal_di(sites, None)


def loso(sites, keep, evalmask):
    """station-LOSO. keep = 평가 대상 측점명 집합, evalmask = 성분별 {측점명: bool}"""
    errs = {"D": [], "I": [], "F": []}
    for i in range(len(sites)):
        nm = sites["name"].values[i]
        if nm not in keep:
            continue
        tr = sites.drop(sites.index[i])
        te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, CRUSTAL, DEGREE)
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEGREE)
        cr = np.nan_to_num(CRUSTAL(te["lat"].values, te["lon"].values), nan=0.0)
        AD = LB.design_DI(te, A)
        ct = LB.crustal_term(te, cr, c)
        crD = float(te["crD"].values[0]) if "crD" in te else 0.0
        crI = float(te["crI"].values[0]) if "crI" in te else 0.0
        if evalmask["D"].get(nm):
            errs["D"].append(te["dD"].values[0] - crD - (AD @ c["D"])[0])
        if evalmask["I"].get(nm):
            errs["I"].append(te["dI"].values[0] - crI - (AD @ c["I"])[0])
        if evalmask["F"].get(nm):
            errs["F"].append(te["dF"].values[0] - ct[0] - (AD @ c["F"])[0])
    return {k: LB.rms(v) for k, v in errs.items()}, \
           {k: len(v) for k, v in errs.items()}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 74)
    print("2019 포함 여부 — degree 0 · Grade A · α=1 · 벡터 OFF · 평가집합 공통 고정")
    print("=" * 74)

    A = build(True)
    B = build(False)
    nA, nB = set(A["name"]), set(B["name"])
    print(f"2019 포함: {len(nA)}측점   2019 제외: {len(nB)}측점")
    print(f"2019 에만 있는 측점: {sorted(nA - nB)}")

    # 평가집합 — 생산 기준(2019 포함) 설정에서 1회 확정
    _, _, INL = LB.fit_regional(A, CRUSTAL, DEGREE)
    names = list(A["name"])
    mask = {k: {names[i]: bool(INL[k][i]) for i in range(len(names))}
            for k in ("D", "I", "F")}

    common = nA & nB
    # 2019 행이 섞여 관측 대표값 자체가 바뀌는 측점
    ptsA = LB.load_all_points(include_2019=True)
    mixed = {s for s, g in ptsA.groupby("station")
             if (g["year"] <= 2019).any() and (g["year"] > 2019).any()}
    pure = common - mixed
    print(f"공통 측점 {len(common)}개   그중 2019 행이 섞인 측점 {sorted(mixed)}")

    for label, keep in (("① 공통 전체", common), ("② 2019 무관 측점만", pure)):
        print(f"\n{label}  (평가 {len(keep)}측점)")
        print(f"{'설정':>12} {'D(도)':>9} {'D(분)':>8} {'I(도)':>9} "
              f"{'F(nT)':>9}   n(D/I/F)")
        for tag, S in (("2019 포함", A), ("2019 제외", B)):
            r, n = loso(S, keep, mask)
            print(f"{tag:>12} {r['D']:>9.4f} {r['D']*60:>8.2f} {r['I']:>9.4f} "
                  f"{r['F']:>9.2f}   {n['D']}/{n['I']}/{n['F']}")
        rA, _ = loso(A, keep, mask)
        rB, _ = loso(B, keep, mask)
        print(f"{'차이(제외-포함)':>12} {rB['D']-rA['D']:>+9.4f} "
              f"{(rB['D']-rA['D'])*60:>+8.2f} {rB['I']-rA['I']:>+9.4f} "
              f"{rB['F']-rA['F']:>+9.2f}")


if __name__ == "__main__":
    main()
