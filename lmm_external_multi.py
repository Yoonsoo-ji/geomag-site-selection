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

**성분별 시간창** — 야장에는 세션 전체(시작~종료) 외에 「편각구간」·「복각구간」이
따로 있다. 한 세션에서 두 구간의 중심시각은 중앙 12분·최대 53분 떨어져 있으므로,
D 는 편각 관측 구간, I 는 복각 관측 구간의 편차로 각각 보정한다. F 는 야장 기록이
일자 단위라 세션 전체 구간을 쓴다.

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


# ── QC ────────────────────────────────────────────────────────────────────
# Kp 는 전지구 지수이고 CYG·KASA 는 한반도에서 직접 잰 값이다. 그래서 Kp 를
# 「제외 기준」으로 쓰지 않고 플래그로만 남기고, 실제 채택 여부는 **관측소가
# 그 시간창에서 실제로 흔들렸는지**(변동폭)로 판단한다.
# 기준선 통계 — 세션 대표값과 같은 규약(중앙값)을 기본으로 한다.
BASELINE_STAT = "median"

KP_QUIET, KP_CAUTION = 2.0, 3.0
# 성분별 문턱 — 편각 KPI 6′ 를 기준으로 잡았다. 관측 창 안에서 관측소 편각이
# 이미 3′ 이상 흔들렸다면 그 세션의 「한 값 대표」가 성립하지 않는다.
RANGE_QUIET = {"D": 3.0, "I": 3.0, "F": 15.0}      # ′, ′, nT
RANGE_CAUTION = {"D": 8.0, "I": 8.0, "F": 40.0}


def kp_flag(kp):
    if kp != kp:
        return "미상"
    if kp <= KP_QUIET:
        return "QUIET"
    return "CAUTION" if kp <= KP_CAUTION else "DISTURBED"


def session_qc(kp, rng: dict):
    """
    세션 채택 등급 — 관측소가 그 창에서 실제로 흔들렸는지로 판정한다.

    rng 는 성분별 변동폭 {"D": 분, "I": 분, "F": nT}. 하나라도 CAUTION 문턱을
    넘으면 그 등급으로 내려간다. Kp 는 전지구 지수이므로 보조로만 쓴다.
    """
    if not rng or all(v != v for v in rng.values()):
        return kp_flag(kp)
    worst = "QUIET"
    for k, v in rng.items():
        if v != v:
            continue
        if v > RANGE_CAUTION.get(k, 1e9):
            worst = "DISTURBED"
        elif v > RANGE_QUIET.get(k, 1e9) and worst != "DISTURBED":
            worst = "CAUTION"
    # 실측이 조용해도 Kp 가 폭풍이면 한 단계만 낮춘다(폐기하지는 않는다)
    if worst == "QUIET" and kp_flag(kp) == "DISTURBED":
        worst = "CAUTION"
    return worst


def parse_window(txt):
    """야장의 「HH:MM~HH:MM」 구간 문자열 → (시작, 종료) time. 없으면 None."""
    if not isinstance(txt, str) or "~" not in txt:
        return None
    try:
        a, b = [x.strip() for x in txt.split("~")]
        f = lambda x: dt.time(*[int(v) for v in x.split(":")[:2]])
        return f(a), f(b)
    except Exception:      # noqa: BLE001
        return None


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
    f = (lambda v: float(v.median())) if BASELINE_STAT == "median" \
        else (lambda v: float(v.mean()))
    return {"X": f(sel.X), "Y": f(sel.Y), "Z": f(sel.Z), "n": len(sel)}


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


def main(mode="plane", only=None, out_csv=None):
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
    use = STATIONS if not only else {k: v for k, v in STATIONS.items()
                                     if k in only}
    for key, st in use.items():
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
        blank = {"Kp": np.nan, "Kp_flag": "미상",
                 "관측소_D변동_분": np.nan, "관측소_I변동_분": np.nan,
                 "관측소_F변동_nT": np.nan, "관측소_dB_nT": np.nan,
                 "QC": "미상", "D_utc": "", "I_utc": "",
                 "dF": np.nan, "dD_arcmin": np.nan, "dI_arcmin": np.nan,
                 "n_station": 0, "보간": "-", "창": "-"}
        if xy is None:
            rows.append({**r, **blank, "상태": "측점 좌표 없음"})
            continue

        # 성분별 관측 구간 — 야장에 「편각구간」·「복각구간」이 따로 있다.
        # 없으면 세션 전체로 물러선다(2019 일부 양식).
        wD = parse_window(r.get("편각구간"))
        wI = parse_window(r.get("복각구간"))
        win = {
            "D": (to_utc(day, wD[0]), to_utc(day, wD[1])) if wD else (t0, t1),
            "I": (to_utc(day, wI[0]), to_utc(day, wI[1])) if wI else (t0, t1),
            "F": (t0, t1),
        }
        used = ("성분별" if (wD and wI) else "세션전체")

        def deviations(a0, a1):
            """
            주어진 시간창에서 관측소별 (위도, 경도, vX, vY, vZ).

            대표값은 **중앙값**이다 — 1분 자료의 스파이크 한두 개에 평균이 끌려간다.
            창 안의 산포(range)는 QC 로 따로 모은다.
            """
            out, spread = [], []
            for key, dfo in obs.items():
                ck = (key, day)
                if ck not in base_cache:
                    base_cache[ck] = quiet_baseline(dfo.reset_index(), kp, day)
                b = base_cache[ck]
                seg = dfo[(dfo.index >= a0) & (dfo.index <= a1)].dropna(subset=["X"])
                if b is None or seg.empty:
                    continue
                out.append((STATIONS[key]["lat"], STATIONS[key]["lon"],
                            seg.X.median() - b["X"], seg.Y.median() - b["Y"],
                            seg.Z.median() - b["Z"]))
                # 관측소 자신의 D·I·F 시계열 변동폭 — X 만 보면 Y·Z 의 흔들림을
                # 놓친다(편각은 Y, 복각은 Z 에 민감하다).
                H = np.hypot(seg.X, seg.Y)
                Dd = np.degrees(np.arctan2(seg.Y, seg.X)) * 60      # 분
                Ii = np.degrees(np.arctan2(seg.Z, H)) * 60          # 분
                Ff = np.hypot(H, seg.Z)
                rng_ = lambda v: float(np.nanmax(v) - np.nanmin(v)) if len(v) else np.nan
                spread.append(dict(
                    key=key, n=int(len(seg)),
                    rD=rng_(Dd), rI=rng_(Ii), rF=rng_(Ff),
                    dB=float(np.hypot(np.hypot(rng_(seg.X), rng_(seg.Y)),
                                      rng_(seg.Z)))))
            return out, spread

        Di, Ii, Fi, X0, Y0, Z0 = LB.igrf_dif(
            np.array([xy["lat"]]), np.array([xy["lon"]]), np.array([0.0]),
            dt.datetime(day.year, day.month, day.day))
        X0, Y0, Z0 = float(X0[0]), float(Y0[0]), float(Z0[0])
        bH = math.hypot(X0, Y0)
        bF = math.hypot(bH, Z0)

        res, how, nst = {}, "-", 0
        qc = {}
        for comp in ("D", "I", "F"):
            dev, spread = deviations(*win[comp])
            if spread:
                qc[comp] = dict(
                    n_min=min(x["n"] for x in spread),
                    rD=max(x["rD"] for x in spread),
                    rI=max(x["rI"] for x in spread),
                    rF=max(x["rF"] for x in spread),
                    dB=max(x["dB"] for x in spread))
            if not dev:
                res[comp] = np.nan
                continue
            nst = max(nst, len(dev))
            vX, how = interpolate([(a, b, c) for a, b, c, _, _ in dev],
                                  None, None, xy["lat"], xy["lon"], mode)
            vY, _ = interpolate([(a, b, d_) for a, b, _, d_, _ in dev],
                                None, None, xy["lat"], xy["lon"], mode)
            vZ, _ = interpolate([(a, b, e) for a, b, _, _, e in dev],
                                None, None, xy["lat"], xy["lon"], mode)
            mH = math.hypot(X0 + vX, Y0 + vY)
            mF = math.hypot(mH, Z0 + vZ)
            if comp == "D":
                res["D"] = math.degrees(math.atan2(Y0 + vY, X0 + vX)
                                        - math.atan2(Y0, X0)) * 60
            elif comp == "I":
                res["I"] = math.degrees(math.atan2(Z0 + vZ, mH)
                                        - math.atan2(Z0, bH)) * 60
            else:
                res["F"] = mF - bF

        if nst == 0:
            rows.append({**r, **blank, "상태": "관측소 자료 없음"})
            continue

        # Kp 는 3시간 구간값이라 세션이 경계를 넘으면 시작시각만으로는 놓친다.
        # 성분 시간창 전체(양 끝 + 30분 간격)에서 최대값을 QC 에 쓴다.
        wlo = min(w[0] for w in win.values())
        whi = max(w[1] for w in win.values())
        probes = pd.date_range(wlo, whi, freq="30min").tolist()
        if probes[-1] != whi:
            probes.append(whi)
        kpv = float(np.nanmax(kp_for(probes, kp).to_numpy()))
        rD = qc.get("D", {}).get("rD", np.nan)
        rI = qc.get("I", {}).get("rI", np.nan)
        rF = qc.get("F", {}).get("rF", np.nan)
        dB = max((v.get("dB", np.nan) for v in qc.values()), default=np.nan)
        rows.append({**r,
                     "Kp": kpv,
                     "Kp_flag": kp_flag(kpv),
                     "관측소_D변동_분": round(rD, 2) if rD == rD else np.nan,
                     "관측소_I변동_분": round(rI, 2) if rI == rI else np.nan,
                     "관측소_F변동_nT": round(rF, 1) if rF == rF else np.nan,
                     "관측소_dB_nT": round(dB, 1) if dB == dB else np.nan,
                     "QC": session_qc(kpv, {"D": rD, "I": rI, "F": rF}),
                     "D_utc": f"{win['D'][0]:%Y-%m-%d %H:%M}~{win['D'][1]:%H:%M}",
                     "I_utc": f"{win['I'][0]:%Y-%m-%d %H:%M}~{win['I'][1]:%H:%M}",
                     "dF": res.get("F", np.nan),
                     "dD_arcmin": res.get("D", np.nan),
                     "dI_arcmin": res.get("I", np.nan),
                     "n_station": nst, "보간": how, "창": used,
                     "상태": "정상"})

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
        print(f"  |I| 평균 {ok.dI_arcmin.abs().mean():.2f}′ 최대 "
              f"{ok.dI_arcmin.abs().max():.2f}′")
        print("  시간창: " + ", ".join(f"{k} {v}건"
                                    for k, v in ok['창'].value_counts().items()))
        print("  QC: " + ", ".join(f"{k} {v}건"
                                   for k, v in ok['QC'].value_counts().items()))
        print("  Kp: " + ", ".join(f"{k} {v}건"
                                   for k, v in ok['Kp_flag'].value_counts().items()))
        # 성분별 창을 따로 쓴 효과 — 같은 세션에서 D 와 I 보정량이 얼마나 갈리나
        both = ok[ok['창'] == '성분별']
        if len(both):
            g = (both.dD_arcmin - both.dI_arcmin).abs()
            print(f"  같은 세션 D-I 보정량 차: 중앙 {g.median():.2f}′ "
                  f"최대 {g.max():.2f}′")

    cols = [c for c in ("측점", "날짜", "시작", "종료", "편각구간", "복각구간",
                        "D_utc", "I_utc", "Kp", "Kp_flag",
                        "관측소_D변동_분", "관측소_I변동_분", "관측소_F변동_nT",
                        "관측소_dB_nT", "QC", "dF", "dD_arcmin", "dI_arcmin",
                        "n_station", "보간", "창", "상태", "파일")
            if c in res.columns]
    out = res[cols].copy()
    out["날짜"] = out["날짜"].astype(str)
    dst = Path(out_csv) if out_csv else OUT_CSV
    dst.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dst, index=False, encoding="utf-8-sig")
    print(f"[저장] {dst}  ({len(out)}행, 정상 {len(ok)}행)")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="plane",
                    choices=["plane", "idw", "nearest"])
    ap.add_argument("--baseline", default="median", choices=["median", "mean"],
                    help="정온야간 기준선 통계(민감도 비교용)")
    ap.add_argument("--stations", default="", help="쉼표로 구분(예: CYG). 비우면 전부")
    ap.add_argument("--out", default="", help="산출 CSV 경로")
    a = ap.parse_args()
    BASELINE_STAT = a.baseline
    main(a.mode, only=[x.strip() for x in a.stations.split(",") if x.strip()] or None,
         out_csv=a.out or None)
