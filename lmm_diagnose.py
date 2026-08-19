# -*- coding: utf-8 -*-
"""
LMM 검증 진단 — 편각 잔차의 원인을 좁히기 위한 일괄 분석.

    python lmm_diagnose.py        -> docs/data/lmm_diagnosis.json

`create_lmm_validation_report.py` 가 이 JSON 을 읽어 보고서를 쓴다.
보고서에 수치를 하드코딩하지 않기 위한 중간 산출물이다.

분석 구성
    1. 계산기 구현 검증        한 지점에서 ppigrf 와 웹 계산기·공식 계산기 대조
    2. 층별 기여 분해          같은 지점에서 Core / Regional / Crustal 분리
    3. 지각 벡터 복원          스칼라 ΔF -> 벡터 b (포텐셜장 FFT 역산)
    4. 가설 A~D 검정           복원 벡터·측정품질·자기거칠기·공간구조
    5. 자료 오류 검정          좌표 민감도·기준시점·자료원
    6. 시기간 대조             2019(검증완료) vs 2022~25(미검증) 겹치는 측점
    7. 2019 단독 구축 가능성    측점 수·공간 분포·자유도
    8. CYG 보유 현황           캐시 실측

⚠️ 표본이 17측점(편각 inlier 14)뿐이므로 상관계수의 통계적 유의성은 대부분
   확보되지 않는다. p 값을 반드시 함께 읽을 것. 이 진단의 결론은 상관보다
   규모 논거와 반변량도에 의존한다.
"""
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal.windows import tukey
from scipy.stats import mannwhitneyu, pearsonr

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import lmm_build as LB   # noqa: E402

OUT_JSON = HERE / "docs" / "data" / "lmm_diagnosis.json"
CYG_CACHE = HERE / "data" / "cyg"

DEGREE = 1                      # 배포 모형(lmm_model.json)과 동일
EPOCH = datetime(2027, 1, 1)
KM_PER_DEG = 111.32

# 검증용 단일 지점 — 사용자가 공식 IGRF 계산기와 웹 계산기로 함께 조회한 지점
PT = dict(lat=36.696, lon=127.392, elev_m=0.0, date="2026-07-28")

# 공식 IGRF 계산기(외부 웹 서비스) 출력 — 화면 표시 자릿수 그대로 기록
OFFICIAL = dict(D=-8.784, I=53.484, X=29892.0, Y=-4619.0,
                H=30247.0, Z=40852.0, F=50831.0)
# docs/lmm.html 을 브라우저에서 직접 실행해 얻은 값
WEBCALC = dict(D=-(8 + 45 / 60 + 17.0 / 3600), I=53 + 32 / 60 + 40.3 / 3600,
               F=50679.0, H=30113.4, X=29762.5, Y=-4583.4, Z=40762.1,
               igrf_D=-(8 + 47 / 60 + 2.6 / 3600),
               igrf_I=53 + 29 / 60 + 2.3 / 3600, igrf_F=50830.7,
               igrf_H=30246.7, igrf_X=29891.9, igrf_Y=-4619.0, igrf_Z=40852.1,
               layer_core=50830.7, layer_regional=-3.89, layer_crustal=-147.76)


def jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return [jsonable(v) for v in o.tolist()]
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def corr(pred, resid):
    """상관·차감 전후 RMS·최적배율. 표본이 작아 p 값을 반드시 함께 본다."""
    p, r = np.asarray(pred, float), np.asarray(resid, float)
    g = np.isfinite(p) & np.isfinite(r)
    p, r = p[g], r[g]
    if len(p) < 4:
        return dict(n=int(len(p)))
    cc, pv = pearsonr(p, r)
    a = float((p @ r) / (p @ p)) if (p @ p) > 0 else 0.0
    return dict(n=int(len(p)), r=float(cc), p=float(pv),
                rms=float(np.sqrt((r ** 2).mean())),
                rms_sub=float(np.sqrt(((r - p) ** 2).mean())),
                rms_scaled=float(np.sqrt(((r - a * p) ** 2).mean())), a=float(a))


# ══════════════════════════════════════════════════════════════════════
def main():
    D = {"generated": datetime.now().isoformat(timespec="seconds"),
         "degree": DEGREE, "epoch_used": EPOCH.isoformat()}

    # ── 자료 적재 ─────────────────────────────────────────────────────
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    sites = LB.aggregate_sites(pts, res)
    lons, lats, grid = LB.load_kigam_grid()
    crustal = LB.CrustalGrid(lons, lats, grid)

    # 배포 모형과 같은 설정으로 진단한다 — 지각 벡터가 켜져 있으면 D·I 잔차도
    # 그 기여를 뺀 뒤의 값이어야 진단이 실제 모형을 가리킨다.
    bxg, byg, bzg = LB.crustal_vector(lons, lats, grid)
    cvec = LB.CrustalVector(lons, lats, bxg, byg, bzg) if LB.CRUSTAL_VECTOR else None
    sites = LB.attach_crustal_di(sites, cvec)
    coef, crust, inliers = LB.fit_regional(sites, crustal, degree=DEGREE)

    A = LB.poly_terms(sites["lat"].values, sites["lon"].values, DEGREE)
    # rD·rI 는 「지각 벡터까지 뺀」 잔차 (crD·crI 는 벡터 미적용 시 0)
    sites["rD"] = (sites["dD"].values - sites["crD"].values - A @ coef["D"]) * 60
    sites["rI"] = (sites["dI"].values - sites["crI"].values - A @ coef["I"]) * 60
    # 지각 벡터 적용 전 잔차 — 가설 A 의 효과를 재현하기 위해 함께 남긴다
    sites["rD_novec"] = (sites["dD"].values - A @ coef["D"]) * 60
    sites["rI_novec"] = (sites["dI"].values - A @ coef["I"]) * 60
    sites["rF"] = sites["dF"].values - A @ coef["F"] - crust
    sites["crustal_nT"] = crust
    mD, mI = inliers["D"], inliers["I"]

    D["crustal_vector_applied"] = bool(LB.CRUSTAL_VECTOR)

    D["dataset"] = dict(
        n_sites=int(len(sites)), n_obs=int(len(pts)),
        year_min=int(pts["year"].min()), year_max=int(pts["year"].max()),
        D_ok=int(sites["D_ok"].sum()), I_ok=int(sites["I_ok"].sum()),
        inlier_D=int(mD.sum()), inlier_I=int(mI.sum()),
        inlier_F=int(inliers["F"].sum()))

    # ── 1. 단일 지점 구현 검증 ────────────────────────────────────────
    import ppigrf
    when = datetime.strptime(PT["date"], "%Y-%m-%d")
    Be, Bn, Bu = ppigrf.igrf(PT["lon"], PT["lat"], PT["elev_m"] / 1000.0, when)
    X0, Y0, Z0 = float(Bn.squeeze()), float(Be.squeeze()), -float(Bu.squeeze())
    H0 = math.hypot(X0, Y0)
    F0 = math.hypot(H0, Z0)
    Dg = math.degrees(math.atan2(Y0, X0))
    Ig = math.degrees(math.atan2(Z0, H0))
    ppi = dict(D=Dg, I=Ig, X=X0, Y=Y0, H=H0, Z=Z0, F=F0)

    # Regional / Crustal 을 배포 모형 계수로 재현.
    # ⚠️ 배포 격자(lmm_model.json)는 정수로 반올림되어 있어 원 KIGAM 격자와
    #    약 0.2 nT 차이가 난다. 계산기 검증이므로 배포 격자를 써야 일치한다.
    MJ = json.loads((HERE / "docs" / "data" / "lmm_model.json")
                    .read_text(encoding="utf-8"))
    R = MJ["regional"]
    u, v = PT["lat"] - R["lat0"], PT["lon"] - R["lon0"]
    reg = {k: float(np.dot(R[k], [1.0, u, v])) for k in ("D", "I", "F")}

    C = MJ["crustal"]
    gv = np.array([np.nan if x is None else x for x in C["values"]],
                  dtype=float).reshape(C["nlat"], C["nlon"])
    fi = (PT["lon"] - C["lon0"]) / C["dlon"]
    fj = (PT["lat"] - C["lat0"]) / C["dlat"]
    i0 = int(np.clip(np.floor(fi), 0, C["nlon"] - 2))
    j0 = int(np.clip(np.floor(fj), 0, C["nlat"] - 2))
    ti, tj = fi - i0, fj - j0
    vv = np.array([gv[j0, i0], gv[j0, i0 + 1], gv[j0 + 1, i0], gv[j0 + 1, i0 + 1]])
    ww = np.array([(1 - ti) * (1 - tj), ti * (1 - tj), (1 - ti) * tj, ti * tj])
    ok_ = np.isfinite(vv)
    cr_pt = float((ww[ok_] * vv[ok_]).sum() / ww[ok_].sum())
    cr_raw = float(crustal(np.array([PT["lat"]]), np.array([PT["lon"]]))[0])

    D_l = Dg + reg["D"]
    I_l = Ig + reg["I"]
    F_l = F0 + reg["F"] + cr_pt
    H_l = F_l * math.cos(math.radians(I_l))
    lmm_pt = dict(D=D_l, I=I_l, F=F_l, H=H_l,
                  X=H_l * math.cos(math.radians(D_l)),
                  Y=H_l * math.sin(math.radians(D_l)),
                  Z=F_l * math.sin(math.radians(I_l)))

    D["point_check"] = dict(
        point=PT, official=OFFICIAL, webcalc=WEBCALC, ppigrf=ppi, lmm=lmm_pt,
        layers=dict(core=F0, regional_F=reg["F"], crustal_F=cr_pt,
                    crustal_F_raw_grid=cr_raw,
                    grid_quantization_nT=abs(cr_pt - cr_raw),
                    regional_D_deg=reg["D"], regional_I_deg=reg["I"]),
        diff_lmm_igrf={k: lmm_pt[k] - ppi[k] for k in ppi},
        # 공식 계산기는 강도를 정수로 표시하므로 반올림 후 비교해야 한다
        vs_official={k: float(ppi[k] - OFFICIAL[k]) for k in OFFICIAL},
        vs_webcalc=dict(D=float(ppi["D"] - WEBCALC["igrf_D"]),
                        I=float(ppi["I"] - WEBCALC["igrf_I"]),
                        F=float(ppi["F"] - WEBCALC["igrf_F"]),
                        lmm_D=float(lmm_pt["D"] - WEBCALC["D"]),
                        lmm_F=float(lmm_pt["F"] - WEBCALC["F"])))

    # ── 2~3. 스칼라 -> 벡터 복원 (포텐셜장 FFT 역산) ──────────────────
    lat0 = float(np.nanmean(sites["lat"]))
    Dm, Im, _, *_ = LB.igrf_dif(np.array([lat0]), np.array([127.5]),
                                np.array([0.0]), [EPOCH])
    Dm, Im = float(Dm[0]), float(Im[0])
    l = math.cos(math.radians(Im)) * math.cos(math.radians(Dm))
    m = math.cos(math.radians(Im)) * math.sin(math.radians(Dm))
    nn = math.sin(math.radians(Im))

    dlat = float(np.median(np.diff(lats)))
    dlon = float(np.median(np.diff(lons)))
    dx_km = dlat * KM_PER_DEG
    dy_km = dlon * KM_PER_DEG * math.cos(math.radians(lat0))

    g0 = np.nan_to_num(grid, nan=0.0)
    ny, nx = g0.shape
    # 역산은 lmm_build 구현을 그대로 쓴다 (사본을 두면 조용히 갈라진다)
    bx, by, bz = bxg, byg, bzg
    okg = np.isfinite(grid)
    err = (l * bx + m * by + nn * bz)[okg] - grid[okg]

    def samp(arr, la_, lo_):
        fi = (lo_ - lons[0]) / dlon
        fj = (la_ - lats[0]) / dlat
        i0 = np.clip(np.floor(fi).astype(int), 0, arr.shape[1] - 2)
        j0 = np.clip(np.floor(fj).astype(int), 0, arr.shape[0] - 2)
        ti, tj = fi - i0, fj - j0
        return ((1 - ti) * (1 - tj) * arr[j0, i0] + ti * (1 - tj) * arr[j0, i0 + 1]
                + (1 - ti) * tj * arr[j0 + 1, i0] + ti * tj * arr[j0 + 1, i0 + 1])

    la, lo = sites["lat"].values, sites["lon"].values
    bxs, bys, bzs = samp(bx, la, lo), samp(by, la, lo), samp(bz, la, lo)
    Di, Ii, Fi, *_ = LB.igrf_dif(la, lo, sites["elev_m"].values,
                                 [EPOCH] * len(sites))
    Xs = Fi * np.cos(np.radians(Ii)) * np.cos(np.radians(Di))
    Ys = Fi * np.cos(np.radians(Ii)) * np.sin(np.radians(Di))
    Zs = Fi * np.sin(np.radians(Ii))
    predD = (np.degrees(np.arctan2(Ys + bys, Xs + bxs)) - Di) * 60
    predI = (np.degrees(np.arctan2(Zs + bzs, np.hypot(Xs + bxs, Ys + bys)))
             - Ii) * 60
    sites["predD"] = predD
    sites["predI"] = predI

    D["vector_recovery"] = dict(
        main_field=dict(D=Dm, I=Im, l=l, m=nn * 0 + m, n=nn),
        grid=dict(nlat=int(ny), nlon=int(nx), dx_km=dx_km, dy_km=dy_km,
                  n_gap=int(np.isnan(grid).sum()),
                  gap_pct=float(100 * np.isnan(grid).sum() / grid.size)),
        recon_rms_nT=float(np.sqrt((err ** 2).mean())),
        recon_max_nT=float(np.abs(err).max()),
        grid_std_nT=float(grid[okg].std()),
        b_east_rms=float(np.sqrt((bys ** 2).mean())),
        b_north_rms=float(np.sqrt((bxs ** 2).mean())),
        scalar_rms=float(np.sqrt((crust ** 2).mean())),
        predD_rms_arcmin=float(np.sqrt((predD[mD] ** 2).mean())),
        predI_rms_arcmin=float(np.sqrt((predI[mD] ** 2).mean())))

    # ── 4. 가설 검정 ──────────────────────────────────────────────────
    gy, gx = np.gradient(g0, dx_km, dy_km)

    def rough(half_km=5.0):
        hj = max(1, int(round(half_km / dx_km)))
        hi = max(1, int(round(half_km / dy_km)))
        o = []
        for a_, b_ in zip(la, lo):
            j = int(round((a_ - lats[0]) / dlat))
            i = int(round((b_ - lons[0]) / dlon))
            w = grid[max(0, j - hj):j + hj + 1, max(0, i - hi):i + hi + 1]
            w = w[np.isfinite(w)]
            o.append(w.std() if w.size >= 4 else np.nan)
        return np.array(o)

    rep = LB.repeatability(pts, res)
    rep = rep.rename(columns={rep.columns[0]: "name"})
    mg = sites.merge(rep[["name", "D_산포_arcmin", "F_산포_nT"]],
                     on="name", how="left")
    sub = mg[np.isfinite(mg["D_산포_arcmin"])]

    D["hypotheses"] = dict(
        # 가설 A 는 「지각 벡터를 빼기 전」 잔차에 대해 검정한다.
        # 모형이 이미 벡터를 빼고 있으면 rD 에는 남아 있지 않기 때문이다.
        A_vector_D=corr(predD[mD], sites["rD_novec"].values[mD]),
        A_vector_I=corr(predI[mI], sites["rI_novec"].values[mI]),
        # 적용 후에도 남는 상관 (0 에 가까워야 정상 — 이중계상 점검)
        A_vector_D_after=corr(predD[mD], sites["rD"].values[mD]),
        B_precision=corr(sub["D_산포_arcmin"].values,
                         np.abs(sub["rD"].values)),
        B_nvisit=corr(mg["n_visit"].values, np.abs(mg["rD"].values)),
        C_roughness=corr(rough()[mD], np.abs(sites["rD"].values[mD])),
        ctrl_scalar=corr(crust[mD], sites["rD"].values[mD]),
        ctrl_grad_east=corr(samp(gx, la, lo)[mD], sites["rD"].values[mD]),
        ctrl_grad_north=corr(samp(gy, la, lo)[mD], sites["rD"].values[mD]))

    one = mg.loc[mg["n_visit"] == 1, "rD"].abs()
    many = mg.loc[mg["n_visit"] > 1, "rD"].abs()
    D["precision_split"] = dict(
        single_median=float(one.median()), single_n=int(len(one)),
        multi_median=float(many.median()), multi_n=int(len(many)))

    # D/I 비대칭 — 동일 측점 집합에서 비교해야 공정하다
    common = mD & np.isfinite(sites["rI"].values)
    rD_rms = float(np.sqrt((sites["rD"].values[common] ** 2).mean()))
    rI_rms = float(np.sqrt((sites["rI"].values[common] ** 2).mean()))
    # 자기 이상이 원인이면 ΔD/ΔI 는 대략 F/H 배 (편각은 H 로, 복각은 F 로 나뉜다)
    expected = float(Fi.mean() / np.hypot(Xs, Ys).mean())
    D["asymmetry"] = dict(
        n=int(common.sum()), rD_rms=rD_rms, rI_rms=rI_rms,
        ratio=rD_rms / rI_rms,
        expected_from_magnetic=expected,
        pred_ratio=float(D["vector_recovery"]["predD_rms_arcmin"]
                         / D["vector_recovery"]["predI_rms_arcmin"]))

    # 잔차를 만들려면 필요한 국소 수평 이상
    need = np.abs(np.tan(np.radians(sites["rD"].values / 60))) * np.hypot(Xs, Ys)
    D["required_local_anomaly"] = dict(
        median_nT=float(np.median(need[mD])), max_nT=float(need[mD].max()),
        recoverable_nT=D["vector_recovery"]["b_east_rms"],
        factor=float(np.median(need[mD]) / D["vector_recovery"]["b_east_rms"]))

    # ── 반변량도 ──────────────────────────────────────────────────────
    sl, so_, sr = la[mD], lo[mD], sites["rD"].values[mD]
    dist, semi = [], []
    for i in range(len(sl)):
        for j in range(i + 1, len(sl)):
            dyk = (sl[i] - sl[j]) * KM_PER_DEG
            dxk = (so_[i] - so_[j]) * KM_PER_DEG * math.cos(math.radians(lat0))
            dist.append(math.hypot(dxk, dyk))
            semi.append(0.5 * (sr[i] - sr[j]) ** 2)
    dist, semi = np.array(dist), np.array(semi)
    edges = [0, 60, 110, 160, 220, 300, 500]
    bins = []
    for a_, b_ in zip(edges[:-1], edges[1:]):
        s_ = (dist >= a_) & (dist < b_)
        if s_.sum() >= 3:
            bins.append(dict(lo=a_, hi=b_, n=int(s_.sum()),
                             semivariance=float(semi[s_].mean())))
    cc, pv = pearsonr(dist, semi)
    D["variogram"] = dict(n_pairs=int(len(dist)), dist_min=float(dist.min()),
                          dist_max=float(dist.max()), bins=bins,
                          sill=float(np.var(sr)), rms=float(np.std(sr)),
                          trend_r=float(cc), trend_p=float(pv))

    # ── 5. 자료 오류 검정 ─────────────────────────────────────────────
    la0, lo0 = 36.5, 127.5
    d0, i0_, _, *_ = LB.igrf_dif(np.array([la0]), np.array([lo0]),
                                 np.array([0.0]), [datetime(2024, 1, 1)])
    sens = {}
    for tag, dla, dlo in (("north_1km", 1 / KM_PER_DEG, 0.0),
                          ("east_1km", 0.0,
                           1 / (KM_PER_DEG * math.cos(math.radians(la0))))):
        d1, i1, _, *_ = LB.igrf_dif(np.array([la0 + dla]), np.array([lo0 + dlo]),
                                    np.array([0.0]), [datetime(2024, 1, 1)])
        sens[tag] = dict(dD_arcmin=float(d1[0] - d0[0]) * 60,
                         dI_arcmin=float(i1[0] - i0_[0]) * 60)
    worst = max(abs(s["dD_arcmin"]) for s in sens.values())
    med_r = float(np.median(np.abs(sites["rD"].values[mD])))
    D["coord_sensitivity"] = dict(
        per_km=sens, median_abs_rD=med_r, km_needed=med_r / worst,
        mion_353m_arcmin=0.353 * worst)

    yrmin = pts.groupby("name")["year"].min()
    m2 = sites.merge(yrmin.rename("y0").reset_index(), on="name")
    D["epoch_trend"] = dict(
        D=dict(**corr(m2["y0"].values[mD].astype(float), sites["rD"].values[mD]),
               slope=float(np.polyfit(m2["y0"].values[mD].astype(float),
                                      sites["rD"].values[mD], 1)[0])),
        I=dict(**corr(m2["y0"].values[mD].astype(float), sites["rI"].values[mD]),
               slope=float(np.polyfit(m2["y0"].values[mD].astype(float),
                                      sites["rI"].values[mD], 1)[0])))

    src19 = m2["y0"].values < 2022
    srcs = {}
    for col in ("rD", "rI", "rF"):
        a_ = sites[col].values[mD & src19]
        b_ = sites[col].values[mD & ~src19]
        a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
        if len(a_) >= 2 and len(b_) >= 2:
            u_, pv_ = mannwhitneyu(a_, b_)
            srcs[col] = dict(n_2019=int(len(a_)), med_2019=float(np.median(a_)),
                             rms_2019=float(np.sqrt((a_ ** 2).mean())),
                             n_survey=int(len(b_)), med_survey=float(np.median(b_)),
                             rms_survey=float(np.sqrt((b_ ** 2).mean())),
                             p=float(pv_))
    D["source_split"] = srcs

    # ── 6. 시기간 대조 ────────────────────────────────────────────────
    pts2 = pts.assign(dD=res["dD"].values * 60, dI=res["dI"].values * 60,
                      dF=res["dF"].values)
    pts2["src"] = np.where(pts2["year"] < 2022, "fb2019", "survey")
    both = sorted(set(pts2.loc[pts2["src"] == "fb2019", "name"])
                  & set(pts2.loc[pts2["src"] == "survey", "name"]))
    cross = []
    for nm in both:
        g = pts2[pts2["name"] == nm]
        row = dict(name=nm)
        for c in ("dD", "dI", "dF"):
            a_ = g.loc[g["src"] == "fb2019", c].mean()
            b_ = g.loc[g["src"] == "survey", c].mean()
            row[c] = dict(fb2019=float(a_), survey=float(b_),
                          diff=float(b_ - a_))
        cross.append(row)
    D["cross_epoch"] = dict(sites=cross,
                            fb2019_only=sorted(
                                set(pts2.loc[pts2["src"] == "fb2019", "name"])
                                - set(pts2.loc[pts2["src"] == "survey", "name"])))

    # ── 7. 2019 단독 구축 가능성 ──────────────────────────────────────
    p19 = pts2[pts2["src"] == "fb2019"]
    s19 = (p19.groupby("name").agg(lat=("lat", "mean"), lon=("lon", "mean"),
                                   n=("year", "size")).reset_index())
    KOREA_N = 38.6
    D["only2019"] = dict(
        n_sites=int(len(s19)), n_sessions=int(len(p19)),
        lat_min=float(s19["lat"].min()), lat_max=float(s19["lat"].max()),
        lon_min=float(s19["lon"].min()), lon_max=float(s19["lon"].max()),
        gap_north_deg=float(KOREA_N - s19["lat"].max()),
        gap_north_km=float((KOREA_N - s19["lat"].max()) * KM_PER_DEG),
        sites=[dict(name=r["name"], lat=float(r["lat"]), lon=float(r["lon"]),
                    n=int(r["n"])) for _, r in s19.sort_values("lat").iterrows()],
        dof=[dict(degree=d, n_coef=(d + 1) * (d + 2) // 2,
                  dof=int(len(s19) - (d + 1) * (d + 2) // 2),
                  robust_min=(d + 1) * (d + 2) // 2 + 3,
                  robust_ok=bool(len(s19) >= (d + 1) * (d + 2) // 2 + 3))
             for d in (1, 2)])

    # ── 8. CYG 보유 현황 ──────────────────────────────────────────────
    cyg = []
    if CYG_CACHE.exists():
        for sub_ in sorted(p for p in CYG_CACHE.iterdir() if p.is_dir()):
            files = sorted(sub_.rglob("cyg_*.csv"))
            ds = sorted(re.search(r"(\d{8})", f.name).group(1) for f in files
                        if re.search(r"(\d{8})", f.name))
            if ds:
                cyg.append(dict(group=sub_.name, n_days=len(ds),
                                first=ds[0], last=ds[-1]))
    D["cyg_cache"] = dict(exists=CYG_CACHE.exists(), groups=cyg,
                          extra=[f.name for f in CYG_CACHE.glob("*.txt")]
                          if CYG_CACHE.exists() else [])

    # 성과표 날짜 조작 여부 (7월 1일 대표값)
    survey_dates = sorted(set(pts2.loc[pts2["src"] == "survey", "date"]
                              .dt.strftime("%m-%d")))
    D["survey_date_stub"] = dict(unique_monthday=survey_dates,
                                 is_synthetic=survey_dates == ["07-01"])

    # ── 측점 표 ───────────────────────────────────────────────────────
    D["site_table"] = [
        dict(name=r["name"], lat=float(r["lat"]), lon=float(r["lon"]),
             n_visit=int(r["n_visit"]), crustal_nT=float(r["crustal_nT"]),
             rD=float(r["rD"]), rI=float(r["rI"]), rF=float(r["rF"]),
             predD=float(r["predD"]), D_ok=bool(r["D_ok"]),
             inlier_D=bool(mD[i]))
        for i, (_, r) in enumerate(sites.iterrows())]

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(jsonable(D), ensure_ascii=False, indent=1),
                        encoding="utf-8")
    print("saved:", OUT_JSON)
    print(f"  측점 {D['dataset']['n_sites']} / 편각 inlier {D['dataset']['inlier_D']}")
    print(f"  편각 잔차 RMS {D['asymmetry']['rD_rms']:.2f}′  "
          f"복원 벡터 예측 {D['vector_recovery']['predD_rms_arcmin']:.2f}′")
    print(f"  가설 A  r={D['hypotheses']['A_vector_D']['r']:+.3f} "
          f"p={D['hypotheses']['A_vector_D']['p']:.3f}")
    return D


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
