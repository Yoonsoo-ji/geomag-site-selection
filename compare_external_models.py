# -*- coding: utf-8 -*-
"""
External 모델 E0/E1/E2 end-to-end 비교 (station-LOSO)
=======================================================

설계안 v2 9항. 보정량 자료원과 공간보간 방식을 바꿔 가며 **처음부터 다시 적합**하고,
평가 측점을 한 번만 정해 고정한 채 LOSO 로 견준다.

    E0  청양 단독            docs/data/external_E0_cyg.csv
    E1  최근접 관측소        docs/data/external_E1_nearest.csv
    E2  4소 1차 평면 보간    docs/data/external_corrections_multi.csv
    E2b 4소 역거리가중       docs/data/external_E2b_idw.csv

적용 방식은 세 가지를 함께 본다.
    subtract_D        편각만
    subtract_DI       편각·복각을 각자 관측구간 보정량으로
    subtract_quiet_DI QC 가 QUIET 인 세션만

⚠ 비교 조건을 모두 같게 맞춘다 — degree 0(현 생산 차수) · Grade A 필터 ·
  지각 결합 α=1 · 지각 벡터 OFF · **평가 측점은 기준 설정에서 한 번만 정해 고정**.

    python compare_external_models.py
"""
import importlib
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "docs" / "data"
DEGREE = 0

SOURCES = {
    "E0 청양단독": DATA / "external_E0_cyg.csv",
    "E1 최근접": DATA / "external_E1_nearest.csv",
    "E2 평면보간": DATA / "external_corrections_multi.csv",
    "E2b 역거리": DATA / "external_E2b_idw.csv",
}
MODES = ["subtract_D", "subtract_DI", "subtract_quiet_DI"]

import lmm_build as LB      # noqa: E402

_lons, _lats, _grid = LB.load_kigam_grid()
CRUSTAL = LB.CrustalGrid(_lons, _lats, _grid)


def build(mode, csv):
    importlib.reload(LB)
    LB.EXTERNAL_MODE = mode
    if csv is not None:
        LB.EXTERNAL_CSV = Path(csv)
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    return LB.attach_crustal_di(sites, None)     # 지각 벡터 OFF (생산 설정)


def loso(sites, inl):
    """station 단위 leave-one-out. 평가 집합 inl 은 밖에서 고정해 넘긴다."""
    errs = {"D": [], "I": [], "F": []}
    for i in range(len(sites)):
        tr = sites.drop(sites.index[i])
        te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, CRUSTAL, DEGREE)
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEGREE)
        cr = np.nan_to_num(CRUSTAL(te["lat"].values, te["lon"].values), nan=0.0)
        AD = LB.design_DI(te, A)
        ct = LB.crustal_term(te, cr, c)
        crD = float(te["crD"].values[0]) if "crD" in te else 0.0
        crI = float(te["crI"].values[0]) if "crI" in te else 0.0
        if inl["D"][i]:
            errs["D"].append(te["dD"].values[0] - crD - (AD @ c["D"])[0])
        if inl["I"][i]:
            errs["I"].append(te["dI"].values[0] - crI - (AD @ c["I"])[0])
        if inl["F"][i]:
            errs["F"].append(te["dF"].values[0] - ct[0] - (AD @ c["F"])[0])
    return {k: LB.rms(v) for k, v in errs.items()}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 74)
    print("External 모델 비교 — degree 0 · Grade A · 지각 α=1 · 벡터 OFF")
    print("=" * 74)

    base = build("none", None)
    _, _, INL = LB.fit_regional(base, CRUSTAL, DEGREE)
    names = list(base["name"])
    print(f"평가 집합(기준 설정에서 1회 확정): "
          f"D {int(INL['D'].sum())} · I {int(INL['I'].sum())} · "
          f"F {int(INL['F'].sum())} / {len(names)}측점")

    rows = [dict(자료원="—", 방식="none", **{k: round(v, 4)
                                           for k, v in loso(base, INL).items()})]
    for label, csv in SOURCES.items():
        if not Path(csv).exists():
            print(f"  ! 없음: {csv}")
            continue
        for m in MODES:
            try:
                r = loso(build(m, csv), INL)
            except Exception as e:              # noqa: BLE001
                print(f"  ! {label}/{m} 실패: {type(e).__name__} {e}")
                continue
            rows.append(dict(자료원=label, 방식=m,
                             **{k: round(v, 4) for k, v in r.items()}))

    df = pd.DataFrame(rows)
    df["D"] = df["D"].astype(float)
    print()
    print(df.to_string(index=False))

    b = df.iloc[0]
    best = df.loc[df["D"].idxmin()]
    print(f"\n기준선(none) D {b['D']:.4f}°  ·  최소 D {best['자료원']}/"
          f"{best['방식']} {best['D']:.4f}°")
    if best["D"] >= b["D"] - 1e-6:
        print("→ 어떤 External 모델도 기준선을 넘지 못한다. 미적용을 유지한다.")
    else:
        print(f"→ {best['자료원']}/{best['방식']} 채택 검토 "
              f"(개선 {b['D']-best['D']:+.4f}°)")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    import datetime as dt
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_External_모델비교.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"[저장] {p}")
    return df


if __name__ == "__main__":
    main()
