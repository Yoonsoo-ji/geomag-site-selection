# -*- coding: utf-8 -*-
"""
④ External 층 — **다중 관측소 공간보간** 판본
================================================

자문 권고(청양 단독 차감 지양, 이천·제주·강릉 등 대체관측소 활용)를 모형에
실제로 반영한다. `lmm_external.py` 는 청양(CYG) 한 곳의 편차를 전국에 그대로
적용했다 — 76~229 km 떨어진 측점에서 「V 가 전역 균일」은 성립하지 않는다.

여기서는 관측소 4곳(청양 CYG · 제주 JJ · 강릉 GN · 이천 ICH)의 세션 구간
편차를 각각 구한 뒤 **공간적으로 보간해 측점 위치의 편차를 추정**한다.

    V̂(r_site, t) = 공간보간{ V_k(t) : 관측소 k }

보간 방식 (`--mode`)
    plane   1차 평면 적합 (기본) — 관측소 3곳 이상일 때. NOC 공간투영의 최소형
    idw     역거리가중 (p=2)
    nearest 최근접 관측소 (거리 효과만 보는 대조군)

⚠ 관측소 좌표계 정합 — 강릉·제주의 Y 는 **계기(변화계) 기준**이라 Y≈0 이다.
   시간 변화만 쓰면 무방하다는 것이 종전 설명이었으나, 편차 벡터를 측점의
   지리 좌표계로 옮기려면 회전을 맞춰야 한다. 관측소별로
   α = D_IGRF(관측소) − median D_관측 을 구해 편차 벡터를 회전시킨다.

⚠ 이 산출물은 **관측 성과 정리 단계**의 보정량이다. 계산기(예측 시점)에는
   들어가지 않는다 — LMM 은 정온시 기준값 모형이다.

실행:
    python lmm_external_multi.py [--mode plane|idw|nearest]
"""

import argparse
import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lmm_cyg import CYG_LAT, CYG_LON, fetch_kp, fetch_range, kp_for
from lmm_external import (BASELINE_WINDOW_DAYS, KST, NIGHT_HOURS,
                          QUIET_KP_MAX, haversine, to_utc)
from lmm_fieldbook import collect

BASE = Path(__file__).parent
KASA_DIR = BASE / "data" / "kasa"
OUT_CSV = BASE / "docs" / "data" / "external_corrections_multi.csv"

# 관측소 좌표 — KASA 3소는 공개 위치 기준의 근사값이다(수 km 오차는 공간보간에
# 영향이 없다). 청양은 INTERMAGNET 등록 좌표.
STATIONS = {
    "CYG": {"name": "청양", "lat": CYG_LAT, "lon": CYG_LON, "src": "cyg"},
    "JJ":  {"name": "제주", "lat": 33.43, "lon": 126.30, "src": "kasa"},
    "GN":  {"name": "강릉", "lat": 37.77, "lon": 128.87, "src": "kasa"},
    "ICH": {"name": "이천", "lat": 37.14, "lon": 127.53, "src": "kasa"},
}


# ------------------------------------------------------------------ 자료 적재
def load_kasa(key: str) -> pd.DataFrame:
    """KASA 연도 캐시를 모두 읽어 UTC 시각의 X·Y·Z 로 반환."""
    files = sorted(KASA_DIR.glob(f"kasa_{key}_*.csv"))
    if not files:
        return pd.DataFrame(columns=["time", "X", "Y", "Z"])
    parts = []
    for f in files:
        d = pd.read_csv(f)
        if "time_kst" not in d.columns:
            continue
        t = pd.to_datetime(d["time_kst"], errors="coerce")
        d = d.assign(time=t.dt.tz_localize(KST).dt.tz_convert("UTC"))
        parts.append(d[["time", "X_nT", "Y_nT", "Z_nT"]]
                     .rename(columns={"X_nT": "X", "Y_nT": "Y", "Z_nT": "Z"}))
    if not parts:
        return pd.DataFrame(columns=["time", "X", "Y", "Z"])
    out = pd.concat(parts, ignore_index=True)
    # 결측 표식 및 상수 placeholder 제거
    out = out.replace(-99999.0, np.nan)
    return out.dropna(subset=["time"]).sort_values("time").drop_duplicates("time")


def align_to_geographic(df: pd.DataFrame, lat: float, lon: float) -> pd.DataFrame:
    """
    계기 좌표계로 기록된 X·Y 를 지리 좌표계로 회전한다.

    관측 편각의 중앙값과 IGRF 편각의 차를 회전각으로 본다. 회전은 상수이므로
    시간 변화 자체는 바뀌지 않지만, **편차 벡터를 측점 좌표계로 옮길 때**
    방향이 맞아야 편각 보정의 부호와 크기가 성립한다.
    """
    import lmm_build as LB
    d = df.dropna(subset=["X", "Y"])
    if d.empty:
        return df.assign(alpha_deg=0.0)
    D_obs = float(np.degrees(np.arctan2(d.Y.median(), d.X.median())))
    yr = int(pd.to_datetime(d.time.iloc[len(d) // 2]).year)
    D_ig = float(LB.igrf_dif(np.array([lat]), np.array([lon]),
                             np.array([0.0]), dt.datetime(yr, 7, 1))[0][0])
    a = math.radians(D_ig - D_obs)
    out = df.copy()
    ca, sa = math.cos(a), math.sin(a)
    out["X"] = df.X * ca - df.Y * sa
    out["Y"] = df.X * sa + df.Y * ca
    out["alpha_deg"] = math.degrees(a)
    return out


def quiet_baseline(obs: pd.DataFrame, kp: pd.DataFrame, center: dt.date):
    """정온야간(Kp<=2 · 지방시 23~03시) 평균. lmm_external 과 같은 규약."""
    lo = pd.Timestamp(center - dt.timedelta(days=BASELINE_WINDOW_DAYS), tz="UTC")
    hi = pd.Timestamp(center + dt.timedelta(days=BASELINE_WINDOW_DAYS + 1), tz="UTC")
    w = obs[(obs.time >= lo) & (obs.time < hi)].dropna(subset=["X", "Y", "Z"])
    if w.empty:
        return None
    h = w.time.dt.tz_convert(KST).dt.hour
    night = (h >= NIGHT_HOURS[0]) | (h < NIGHT_HOURS[1])
    quiet = kp_for(w.time, kp).to_numpy() <= QUIET_KP_MAX
    sel = w[night.to_numpy() & quiet]
    if len(sel) < 60:
        sel = w[night.to_numpy()]
    if len(sel) < 60:
        return None
    return {"X": sel.X.mean(), "Y": sel.Y.mean(), "Z": sel.Z.mean(), "n": len(sel)}


# ------------------------------------------------------------------ 공간보간
def interpolate(vals, slat, slon, lat, lon, mode="plane"):
    """
    관측소별 편차 vals[(lat, lon, v)] 를 측점 (lat, lon) 으로 보간.

    plane 은 1차 평면(3계수)이라 관측소 3곳 이상이 필요하다. 미달이면 idw 로
    자동 강등한다 — 자료가 없는 날 조용히 0 을 쓰지 않기 위해서다.
    """
    a = np.array([[p[0], p[1]] for p in vals], float)
    v = np.array([p[2] for p in vals], float)
    if mode == "nearest" or len(v) == 1:
        d = haversine(a[:, 0], a[:, 1], lat, lon)
        return float(v[np.argmin(d)]), "nearest"
    if mode == "plane" and len(v) >= 3:
        A = np.column_stack([np.ones(len(v)), a[:, 0] - lat, a[:, 1] - lon])
        c, *_ = np.linalg.lstsq(A, v, rcond=None)
        return float(c[0]), "plane"
    d = haversine(a[:, 0], a[:, 1], lat, lon)
    w = 1.0 / np.maximum(d, 1.0) ** 2
    return float((w * v).sum() / w.sum()), "idw"


def main(mode="plane"):
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 78)
    print(f"④ External 층 — 다중 관측소 공간보간 ({mode})")
    print("=" * 78)

    import lmm_build as LB

    fb = collect()
    fb = fb[fb["시작"].notna()].copy()
    days = sorted(set(fb["날짜"]))
    print(f"야장 유효 세션 {len(fb)}건 / {len(days)}일 "
          f"({min(days)} ~ {max(days)})")

    # 측점 좌표 — 성과표·야장 공통 명칭 기준
    pts = LB.load_all_points(include_2019=True)
    site_xy = (pts.groupby("name")[["lat", "lon"]].mean().to_dict("index"))

    kp = fetch_kp()

    # ── 관측소 자료 적재 ────────────────────────────────────────────────
    obs = {}
    for key, st in STATIONS.items():
        if st["src"] == "kasa":
            d = load_kasa(key)
        else:
            ranges = []
            for day in days:
                ranges.append((day - dt.timedelta(days=BASELINE_WINDOW_DAYS),
                               day + dt.timedelta(days=BASELINE_WINDOW_DAYS)))
            ranges.sort()
            merged = [list(ranges[0])]
            for s, e in ranges[1:]:
                if s <= merged[-1][1] + dt.timedelta(days=1):
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            parts = [fetch_range(s, e, quiet=True) for s, e in merged]
            d = pd.concat(parts, ignore_index=True)
            d["time"] = pd.to_datetime(d["time"], utc=True)
            d = d.sort_values("time").drop_duplicates("time")
        if d.empty:
            print(f"  · {st['name']}({key}): 자료 없음 — 제외")
            continue
        d = align_to_geographic(d, st["lat"], st["lon"])
        obs[key] = d.set_index("time")
        print(f"  · {st['name']}({key}): {len(d):,}행 "
              f"(유효 {int(d.X.notna().sum()):,}, 좌표계 회전 "
              f"{float(d['alpha_deg'].iloc[0]):+.2f}°)")

    if not obs:
        print("관측소 자료가 없다 — 중단")
        return None

    # ── 세션별 보정량 ───────────────────────────────────────────────────
    base_cache = {}
    rows = []
    for _, r in fb.iterrows():
        day, site = r["날짜"], r["측점"]
        t0, t1 = to_utc(day, r["시작"]), to_utc(day, r["종료"])
        xy = site_xy.get(site)
        if xy is None:
            rows.append({**r, "Kp": np.nan, "dF": np.nan, "dD_arcmin": np.nan,
                         "n_station": 0, "보간": "-", "상태": "측점 좌표 없음"})
            continue

        dev = []          # (lat, lon, vX, vY, vZ)
        for key, d in obs.items():
            ck = (key, day)
            if ck not in base_cache:
                base_cache[ck] = quiet_baseline(d.reset_index(), kp, day)
            b = base_cache[ck]
            seg = d[(d.index >= t0) & (d.index <= t1)].dropna(subset=["X"])
            if b is None or seg.empty:
                continue
            dev.append((STATIONS[key]["lat"], STATIONS[key]["lon"],
                        seg.X.mean() - b["X"], seg.Y.mean() - b["Y"],
                        seg.Z.mean() - b["Z"]))

        if not dev:
            rows.append({**r, "Kp": np.nan, "dF": np.nan, "dD_arcmin": np.nan,
                         "n_station": 0, "보간": "-", "상태": "관측소 자료 없음"})
            continue

        vX, how = interpolate([(a, b, c) for a, b, c, _, _ in dev],
                              None, None, xy["lat"], xy["lon"], mode)
        vY, _ = interpolate([(a, b, d_) for a, b, _, d_, _ in dev],
                            None, None, xy["lat"], xy["lon"], mode)
        vZ, _ = interpolate([(a, b, e) for a, b, _, _, e in dev],
                            None, None, xy["lat"], xy["lon"], mode)

        # 편차를 **측점의 기준장**(IGRF)에 얹어 F·D 변화로 환산한다.
        # 관측소 기준장을 쓰던 종전 방식보다 측점 위도차를 바르게 반영한다.
        Di, Ii, Fi, X0, Y0, Z0 = LB.igrf_dif(
            np.array([xy["lat"]]), np.array([xy["lon"]]), np.array([0.0]),
            dt.datetime(day.year, day.month, day.day))
        X0, Y0, Z0 = float(X0[0]), float(Y0[0]), float(Z0[0])
        bH = math.hypot(X0, Y0)
        bF = math.hypot(bH, Z0)
        mH = math.hypot(X0 + vX, Y0 + vY)
        mF = math.hypot(mH, Z0 + vZ)
        dF = mF - bF
        dD = math.degrees(math.atan2(Y0 + vY, X0 + vX)
                          - math.atan2(Y0, X0)) * 60

        rows.append({**r, "Kp": float(kp_for([t0], kp).iloc[0]),
                     "dF": dF, "dD_arcmin": dD, "n_station": len(dev),
                     "보간": how, "상태": "정상"})

    res = pd.DataFrame(rows)
    ok = res[res["상태"] == "정상"]
    print(f"\n보정 산출 {len(ok)}건 / 실패 {len(res)-len(ok)}건")
    if len(ok):
        print(f"  관측소 수 분포: "
              + ", ".join(f"{k}소 {v}건" for k, v in
                          sorted(ok['n_station'].value_counts().items())))
        print(f"  |F| 평균 {ok.dF.abs().mean():.1f} nT  최대 {ok.dF.abs().max():.1f}")
        print(f"  |D| 평균 {ok.dD_arcmin.abs().mean():.2f}′ 최대 "
              f"{ok.dD_arcmin.abs().max():.2f}′")

    cols = [c for c in ("측점", "날짜", "시작", "종료", "Kp", "dF", "dD_arcmin",
                        "n_station", "보간", "상태", "파일") if c in res.columns]
    out = res[cols].copy()
    out["날짜"] = out["날짜"].astype(str)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"[저장] {OUT_CSV}  ({len(out)}행, 정상 {len(ok)}행)")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="plane",
                    choices=["plane", "idw", "nearest"])
    a = ap.parse_args()
    main(a.mode)
