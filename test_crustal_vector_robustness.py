# -*- coding: utf-8 -*-
"""
지각 벡터(D·I) 강건성 시험 — 설계안 v2 검토 5항
==================================================

복각이 0.2430° → 0.09° 수준으로 줄어드는 신호는 우연으로 보기 어렵다. 그러나
현재 역산은 KIGAM 격자의 결측 39 % 를 **0 nT 로 채운** 뒤 FFT 를 돌린다. 결측과
「이상 0」은 다른 뜻이므로, 그 개선이 채움 방식의 인공물일 수 있다.

네 조건에서 부호와 크기가 유지되는지 본다.

    ① 채움 방식      zero / mean / harmonic(라플라스 내삽)
    ② 결측 경계 근접  경계에서 N 격자 이내 측점을 빼도 개선이 남는가
    ③ 창·패딩        FFT 앞의 tukey 창·반사 패딩은 이미 적용 — 창 폭을 바꿔 본다
    ④ 블록 홀드아웃   남/중/북, 동/서 블록을 통째로 빼고 예측

네 조건에서 모두 ~0.09° 수준이면 「조건부 생산후보」로 올릴 근거가 된다.
하나에서만 좋아지면 넣지 않는다.

    python test_crustal_vector_robustness.py
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DEGREE = 0

import lmm_build as LB      # noqa: E402


# ── 결측 채움 ────────────────────────────────────────────────────────────
def fill_zero(g):
    return np.nan_to_num(g, nan=0.0)


def fill_mean(g):
    out = g.copy()
    out[np.isnan(out)] = np.nanmean(g)
    return out


def fill_harmonic(g, iters=400):
    """
    라플라스 내삽 — 결측 화소를 이웃 평균으로 반복 완화한다.

    결측 영역에서 ∇²φ = 0 을 만족시키므로 경계에서 값이 부드럽게 이어진다.
    0 채움이 만드는 인공 계단이 사라진다.
    """
    out = np.where(np.isnan(g), np.nanmean(g), g).astype(float)
    mask = np.isnan(g)
    if not mask.any():
        return out
    for _ in range(iters):
        nb = np.zeros_like(out)
        cnt = np.zeros_like(out)
        for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            nb += np.roll(out, sh, axis=ax)
            cnt += 1
        new = nb / cnt
        out[mask] = new[mask]
    return out


FILLS = {"zero": fill_zero, "mean": fill_mean, "harmonic": fill_harmonic}


def gap_distance(g):
    """각 화소에서 가장 가까운 결측까지의 격자 거리(대략, 반복 팽창)."""
    gap = np.isnan(g)
    if not gap.any():
        return np.full(g.shape, 99)
    dist = np.where(gap, 0, 99).astype(float)
    for d in range(1, 8):
        cur = dist == 99
        nb = np.zeros_like(gap)
        for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            nb |= np.roll(dist <= d - 1, sh, axis=ax)
        dist[cur & nb] = d
    return dist


# ── 평가 ────────────────────────────────────────────────────────────────
def prep(fill=None, alpha_win=None):
    lo, la, g = LB.load_kigam_grid()
    crustal = LB.CrustalGrid(lo, la, g)
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    if fill is None:
        return sites, crustal, LB.attach_crustal_di(sites, None), (lo, la, g)
    gf = FILLS[fill](g)
    bx, by, bz = LB.crustal_vector(lo, la, gf)
    cvec = LB.CrustalVector(lo, la, bx, by, bz)
    return sites, crustal, LB.attach_crustal_di(sites, cvec), (lo, la, g)


def loso(sites, crustal, keep=None):
    inl = LB.eval_mask(sites, LB.fit_regional(sites, crustal, DEGREE)[2])
    errs = {"D": [], "I": []}
    for i in range(len(sites)):
        if keep is not None and not keep[i]:
            continue
        tr = sites.drop(sites.index[i])
        te = sites.iloc[[i]]
        c, _, _ = LB.fit_regional(tr, crustal, DEGREE)
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEGREE)
        AD = LB.design_DI(te, A)
        crD = float(te["crD"].values[0]); crI = float(te["crI"].values[0])
        if inl["D"][i]:
            errs["D"].append(te["dD"].values[0] - crD - (AD @ c["D"])[0])
        if inl["I"][i]:
            errs["I"].append(te["dI"].values[0] - crI - (AD @ c["I"])[0])
    return {k: LB.rms(v) for k, v in errs.items()}, \
           {k: len(v) for k, v in errs.items()}


def block_holdout(sites, crustal, blocks):
    """블록을 통째로 빼고 나머지로 적합해 그 블록을 예측."""
    inl = LB.eval_mask(sites, LB.fit_regional(sites, crustal, DEGREE)[2])
    errs = {"D": [], "I": []}
    for name, sel in blocks.items():
        tr = sites[~sel]
        if len(tr) < 4:
            continue
        c, _, _ = LB.fit_regional(tr, crustal, DEGREE)
        te = sites[sel]
        idx = np.where(sel)[0]
        A = LB.poly_terms(te["lat"].values, te["lon"].values, DEGREE)
        AD = LB.design_DI(te, A)
        for k, i in enumerate(idx):
            if inl["D"][i]:
                errs["D"].append(te["dD"].values[k] - te["crD"].values[k]
                                 - (AD @ c["D"])[k])
            if inl["I"][i]:
                errs["I"].append(te["dI"].values[k] - te["crI"].values[k]
                                 - (AD @ c["I"])[k])
    return {k: LB.rms(v) for k, v in errs.items()}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 76)
    print("지각 벡터 강건성 시험 — 평가 집합 동결 · degree 0 · Grade A")
    print("=" * 76)

    rows = []
    sites0, crustal, s_off, (lo, la, g) = prep(None)
    r, n = loso(s_off, crustal)
    rows.append(dict(조건="벡터 OFF (기준)", LOO_D=round(r["D"], 4),
                     LOO_I=round(r["I"], 4), n_D=n["D"], n_I=n["I"]))

    cache = {}
    for fill in ("zero", "mean", "harmonic"):
        _, _, s_on, _ = prep(fill)
        cache[fill] = s_on
        r, n = loso(s_on, crustal)
        rows.append(dict(조건=f"벡터 ON · {fill} 채움", LOO_D=round(r["D"], 4),
                         LOO_I=round(r["I"], 4), n_D=n["D"], n_I=n["I"]))

    # ② 결측 경계 근접 측점 배제
    dist = gap_distance(g)
    dlon = float(np.median(np.diff(lo))); dlat = float(np.median(np.diff(la)))
    def near_gap(sites, n_cell):
        i = np.clip(np.rint((sites["lon"].values - lo[0]) / dlon).astype(int),
                    0, len(lo) - 1)
        j = np.clip(np.rint((sites["lat"].values - la[0]) / dlat).astype(int),
                    0, len(la) - 1)
        return dist[j, i] <= n_cell
    for n_cell in (1, 2):
        s_on = cache["harmonic"]
        keep = ~near_gap(s_on, n_cell)
        r, n = loso(s_on, crustal, keep=keep)
        rows.append(dict(조건=f"harmonic · 결측 {n_cell}칸 이내 배제",
                         LOO_D=round(r["D"], 4), LOO_I=round(r["I"], 4),
                         n_D=n["D"], n_I=n["I"]))

    df = pd.DataFrame(rows)
    print()
    print(df.to_string(index=False))

    # ④ 블록 홀드아웃
    print("\n블록 홀드아웃 (블록을 통째로 빼고 예측)")
    for lab, s_ in (("벡터 OFF", s_off), ("벡터 ON(harmonic)", cache["harmonic"])):
        la_ = s_["lat"].values; lo_ = s_["lon"].values
        blocks = {
            "남부(<35.5N)": la_ < 35.5,
            "중부(35.5~36.9N)": (la_ >= 35.5) & (la_ < 36.9),
            "북부(>=36.9N)": la_ >= 36.9,
        }
        rb = block_holdout(s_, crustal, blocks)
        blocks_ew = {"서부(<127.6E)": lo_ < 127.6, "동부(>=127.6E)": lo_ >= 127.6}
        re = block_holdout(s_, crustal, blocks_ew)
        print(f"  {lab:18} 남중북 D {rb['D']:.4f}° I {rb['I']:.4f}°  |  "
              f"동서 D {re['D']:.4f}° I {re['I']:.4f}°")

    print("\n판정 기준: 세 채움 방식과 경계 배제, 블록 홀드아웃에서 모두 복각이")
    print("           기준(OFF)보다 뚜렷이 낮으면 조건부 생산후보로 올린다.")
    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    import datetime as dt
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_지각벡터_강건성.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"[저장] {p}")


if __name__ == "__main__":
    main()
