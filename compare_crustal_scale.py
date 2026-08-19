# -*- coding: utf-8 -*-
"""
지각 결합계수 α 회귀검증 — F0 / F1 / Fα (설계안 v2 검토 4항)
===============================================================

    F0  F = IGRF + R_F                    (지각 미사용)
    F1  F = IGRF + R_F + A_KIGAM          (1:1)
    Fα  F = IGRF + R_F + α·A_KIGAM        (α 를 자료로 추정)

⚠ **정보 누출 점검이 이 스크립트의 요점이다.** α 를 전체 자료에서 한 번 구해
   LOSO 에 넣으면 그 자체가 누출이다. 여기서는 fold 마다 훈련 측점만으로
   α 를 다시 추정하고, fold 별 α 의 산포도 함께 출력한다.

평가 측점은 `docs/data/eval_set.json` 에 **동결된 집합**을 쓴다.

    python compare_crustal_scale.py
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

_lo, _la, _g = LB.load_kigam_grid()
CRUSTAL = LB.CrustalGrid(_lo, _la, _g)


def build(scale_mode, timeterm=False):
    importlib.reload(LB)
    LB.CRUSTAL_SCALE_MODE = scale_mode
    LB.REGIONAL_TIME_TERM = timeterm
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    return LB.attach_crustal_di(sites, None)      # 생산 설정: 지각 벡터 OFF


def loso(sites):
    """fold 마다 훈련 측점만으로 재적합 — α 도 fold 안에서 다시 추정된다."""
    inl = LB.eval_mask(sites, LB.fit_regional(sites, CRUSTAL, DEGREE)[2])
    errs = {"D": [], "I": [], "F": []}
    alphas = []
    for i in range(len(sites)):
        tr = sites.drop(sites.index[i])
        te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, CRUSTAL, DEGREE)
        alphas.append(c.get("alpha", np.nan))
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEGREE)
        cr = np.nan_to_num(CRUSTAL(te["lat"].values, te["lon"].values), nan=0.0)
        AD = LB.design_DI(te, A)
        ct = LB.crustal_term(te, cr, c)
        if inl["D"][i]:
            errs["D"].append(te["dD"].values[0] - (AD @ c["D"])[0])
        if inl["I"][i]:
            errs["I"].append(te["dI"].values[0] - (AD @ c["I"])[0])
        if inl["F"][i]:
            errs["F"].append(te["dF"].values[0] - ct[0] - (AD @ c["F"])[0])
    return ({k: LB.rms(v) for k, v in errs.items()}, np.array(alphas, float))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 70)
    print("지각 결합계수 α 회귀검증 — 평가 집합 동결 · fold 내 α 재추정")
    print("=" * 70)

    rows = []
    for lab, mode, tt in [("F0  지각 미사용", "none", False),
                          ("F1  1:1 결합", "one", False),
                          ("Fα  결합계수 추정", "estimate", False),
                          ("F1 + 시간항 R1T", "one", True),
                          ("Fα + 시간항 R1T", "estimate", True)]:
        sites = build(mode, tt)
        full = LB.fit_regional(sites, CRUSTAL, DEGREE)[0]
        r, al = loso(sites)
        a_fold = al[np.isfinite(al)]
        rows.append(dict(
            모델=lab,
            α_전체=round(full.get("alpha", np.nan), 3),
            α_fold평균=round(float(a_fold.mean()), 3) if len(a_fold) else np.nan,
            α_fold범위=(f"{a_fold.min():.3f}~{a_fold.max():.3f}"
                       if len(a_fold) else "—"),
            LOO_D=round(r["D"], 4), LOO_I=round(r["I"], 4),
            LOO_F=round(r["F"], 1)))

    df = pd.DataFrame(rows)
    print()
    print(df.to_string(index=False))

    f1 = df[df["모델"].str.startswith("F1  ")].iloc[0]
    fa = df[df["모델"].str.startswith("Fα  ")].iloc[0]
    print()
    if fa["LOO_F"] < f1["LOO_F"]:
        print(f"→ α 추정이 낫다 ({fa['LOO_F']} < {f1['LOO_F']} nT). α 채택 검토.")
    else:
        print(f"→ α=1 이 낫다 ({f1['LOO_F']} < {fa['LOO_F']} nT). "
              "α 를 자유롭게 두면 과적합이다.")
    print("  fold 별 α 산포가 좁을수록 추정이 안정적이다. 위 표의 α_fold범위 참조.")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    import datetime as dt
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_지각결합계수_비교.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"[저장] {p}")
    return df


if __name__ == "__main__":
    main()
