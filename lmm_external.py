# -*- coding: utf-8 -*-
"""
④ External 층 — 야장 일시 × CYG × Kp 결합
==========================================

야장에서 복원한 관측 일시에 대해 청양(CYG) 관측소의 외부장 변동을 구하고,
지표 측정 성과에서 그 변동을 제거한다(설명자료 p.10-11 "가상 배리오미터").

    B_baseline(r) ≈ B_측정(r, t_방문) − V̂(r, t_방문)

기준선(baseline) 정의
---------------------
V(t) 는 "정온 상태 대비 편차"이므로 기준선을 먼저 정해야 한다.
표준 방식대로 **정온야간 평균**을 쓴다:

    · 관측월 내에서 Kp ≤ 2 인 3시간 구간
    · 지방시 23:00~03:00 (야간, Sq 전류 최소)
    · 해당 조건의 CYG X·Y·Z 평균 = 기준선

한계
----
1차 근사로 V 가 한반도 전역에서 동일하다고 본다("가상 배리오미터" 최단순형).
실제로는 CYG 로부터 거리에 따라 감쇠·위상차가 있으며, 설명자료의 NOC
모드 분해는 이를 공간 투영으로 처리한다. 측점-CYG 거리를 함께 출력하므로
근사 타당성을 판단할 수 있다.

실행:
    python lmm_external.py
"""

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lmm_cyg import CYG_LAT, CYG_LON, best_dataset, fetch_kp, fetch_range, kp_for
from lmm_fieldbook import collect

KST = dt.timezone(dt.timedelta(hours=9))

# 정온야간 기준선 조건
QUIET_KP_MAX = 2.0
NIGHT_HOURS = (23, 3)   # 지방시 23시 ~ 익일 3시
BASELINE_WINDOW_DAYS = 10   # 관측일 전후 며칠까지 기준선 산정에 쓸지


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def to_utc(day: dt.date, t: dt.time) -> pd.Timestamp:
    """야장 시각은 한국 표준시(KST)로 본다."""
    return pd.Timestamp(dt.datetime.combine(day, t, tzinfo=KST)).tz_convert("UTC")


def derive_dif(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["H"] = np.hypot(d.X, d.Y)
    d["D_deg"] = np.degrees(np.arctan2(d.Y, d.X))
    d["I_deg"] = np.degrees(np.arctan2(d.Z, d.H))
    return d


def quiet_baseline(cyg: pd.DataFrame, kp: pd.DataFrame, center: dt.date):
    """관측일 전후 구간의 정온야간 평균을 기준선으로 산출."""
    lo = pd.Timestamp(center - dt.timedelta(days=BASELINE_WINDOW_DAYS), tz="UTC")
    hi = pd.Timestamp(center + dt.timedelta(days=BASELINE_WINDOW_DAYS + 1), tz="UTC")
    w = cyg[(cyg.time >= lo) & (cyg.time < hi)].dropna(subset=["X", "Y", "Z"])
    if w.empty:
        return None

    local = w.time.dt.tz_convert(KST)
    h = local.dt.hour
    night = (h >= NIGHT_HOURS[0]) | (h < NIGHT_HOURS[1])
    quiet = kp_for(w.time, kp).to_numpy() <= QUIET_KP_MAX

    sel = w[night.to_numpy() & quiet]
    if len(sel) < 60:      # 최소 1시간 분량은 있어야 신뢰
        sel = w[night.to_numpy()]
    if len(sel) < 60:
        return None

    return {
        "X": sel.X.mean(), "Y": sel.Y.mean(), "Z": sel.Z.mean(),
        "n": len(sel),
    }


def main():
    print("=" * 78)
    print("④ External 층 — 야장 일시 기반 외부장 보정")
    print("=" * 78)

    fb = collect()
    fb = fb[fb["시작"].notna()].copy()
    print(f"\n야장에서 복원한 유효 세션 {len(fb)}건")

    # 필요한 날짜 범위 (기준선 산정 여유 포함)
    days = sorted(set(fb["날짜"]))
    ranges = []
    for d in days:
        ranges.append((d - dt.timedelta(days=BASELINE_WINDOW_DAYS),
                       d + dt.timedelta(days=BASELINE_WINDOW_DAYS)))
    # 겹치는 구간 병합
    ranges.sort()
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + dt.timedelta(days=1):
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    total = sum((e - s).days + 1 for s, e in merged)
    grades = sorted({best_dataset(d).split('/')[1] for d in days})
    print(f"CYG 수집 대상 {total}일 / {len(merged)}구간  (등급: {', '.join(grades)})")

    parts = []
    for s, e in merged:
        print(f"  · {s} ~ {e}")
        parts.append(fetch_range(s, e, quiet=True))
    cyg = pd.concat(parts, ignore_index=True)
    cyg["time"] = pd.to_datetime(cyg["time"], utc=True)
    cyg = cyg.sort_values("time").drop_duplicates("time")
    print(f"CYG 1분 자료 {len(cyg):,}행 (유효 {cyg.X.notna().sum():,})")

    kp = fetch_kp()
    print(f"Kp 지수 {len(kp):,}행 ({kp.time.min():%Y-%m-%d} ~ {kp.time.max():%Y-%m-%d})")

    cygd = derive_dif(cyg).set_index("time")

    rows = []
    base_cache = {}
    for _, r in fb.iterrows():
        day = r["날짜"]
        if day not in base_cache:
            base_cache[day] = quiet_baseline(cyg, kp, day)
        base = base_cache[day]

        t0 = to_utc(day, r["시작"])
        t1 = to_utc(day, r["종료"])
        seg = cygd[(cygd.index >= t0) & (cygd.index <= t1)].dropna(subset=["X"])

        if base is None or seg.empty:
            rows.append({**r, "Kp": np.nan, "dF": np.nan, "dD_arcmin": np.nan,
                         "상태": "CYG 자료 없음"})
            continue

        # 세션 구간 평균 외부장 편차
        vX = seg.X.mean() - base["X"]
        vY = seg.Y.mean() - base["Y"]
        vZ = seg.Z.mean() - base["Z"]

        # F·D 로 환산 (기준선 상태의 벡터에 편차를 더해 차이를 구함)
        bH = np.hypot(base["X"], base["Y"])
        bF = np.hypot(bH, base["Z"])
        mH = np.hypot(base["X"] + vX, base["Y"] + vY)
        mF = np.hypot(mH, base["Z"] + vZ)
        dF = mF - bF
        dD = np.degrees(np.arctan2(base["Y"] + vY, base["X"] + vX)
                        - np.arctan2(base["Y"], base["X"])) * 60

        kpv = float(kp_for([t0], kp).iloc[0])
        rows.append({**r, "Kp": kpv, "dF": dF, "dD_arcmin": dD, "상태": "정상"})

    res = pd.DataFrame(rows)

    ok = res[res["상태"] == "정상"].copy()
    print(f"\n보정 산출 {len(ok)}건 / 실패 {len(res)-len(ok)}건")

    if len(ok):
        show = ok.copy()
        show["보정F"] = show["dF"].map(lambda v: f"{-v:+.1f}")
        show["보정D"] = show["dD_arcmin"].map(lambda v: f"{-v:+.2f}")
        show["Kp"] = show["Kp"].map(lambda v: f"{v:.1f}")
        print("\n--- 세션별 외부장 보정량 (측정값에서 빼야 할 양의 부호 반전) ---")
        print(show[["측점", "날짜", "편각구간", "Kp", "보정F", "보정D"]]
              .rename(columns={"보정F": "F보정_nT", "보정D": "D보정_arcmin"})
              .to_string(index=False))

        print("\n--- 보정량 규모 ---")
        print(f"  |F| 평균 {ok.dF.abs().mean():.1f} nT   최대 {ok.dF.abs().max():.1f} nT")
        print(f"  |D| 평균 {ok.dD_arcmin.abs().mean():.2f}'  최대 "
              f"{ok.dD_arcmin.abs().max():.2f}'  "
              f"({ok.dD_arcmin.abs().max()/60:.3f}°)")

        # 같은 날 같은 측점인데도 세션 시각에 따라 보정량이 달라지는 폭.
        # 이것이 곧 "측정 시각을 모를 때 성과에 남는 오차"의 실측치다.
        # (야장의 F 는 세션별이 아니라 일자별 단일값이므로 원자료 산포 비교는
        #  성립하지 않는다 — 보정량 자체의 변동폭으로 평가한다)
        print("\n--- 동일 측점·동일일 세션 간 보정량 변동 ---")
        print("    = 관측 시각을 모르면 성과에 그대로 남는 오차")
        spreads = []
        for (site, day), g in ok.groupby(["측점", "날짜"]):
            if len(g) < 2:
                continue
            fs = g["dF"].to_numpy().ptp()
            ds = g["dD_arcmin"].to_numpy().ptp()
            spreads.append((fs, ds))
            print(f"  {site} {day} ({len(g)}세션): "
                  f"F {fs:5.1f} nT   D {ds:4.2f}'")
        if spreads:
            fmax = max(s[0] for s in spreads)
            dmax = max(s[1] for s in spreads)
            print(f"  -> 최대 F {fmax:.1f} nT (KPI 50), "
                  f"D {dmax:.2f}' = {dmax/60:.3f}° (KPI 0.1°)")

        print(f"\n--- Kp 판정 (정온 기준 Kp<=2) ---")
        print(f"  범위 {ok.Kp.min():.1f} ~ {ok.Kp.max():.1f}   "
              f"정온 {int((ok.Kp<=2).sum())}/{len(ok)}세션")
        storm = ok[ok.Kp > 2]
        if len(storm):
            print(f"  [경고] 교란시 관측 {len(storm)}건 — 성과 신뢰도 재검토 필요:")
            for _, s in storm.iterrows():
                print(f"     · {s['측점']} {s['날짜']} {s['편각구간']}  "
                      f"Kp={s['Kp']:.1f}  F보정 {-s['dF']:+.1f} nT")

    save_corrections(res)
    return res


# 세션별 외부장 보정량 — lmm_build 가 (측점, 연도) 평균으로 소비한다.
CORRECTION_CSV = Path(__file__).parent / "docs" / "data" / "external_corrections.csv"


def save_corrections(res: pd.DataFrame) -> None:
    """세션별 보정량을 저장한다.

    부호 규약 — `dF`·`dD_arcmin` 은 **관측시각의 외부장 편차**(기준선 대비)다.
    측정값에서 이 값을 **빼야** 정온 기준값이 된다. 소비 측에서 부호를 뒤집지
    말고 그대로 빼도록 컬럼명을 `dF`·`dD_arcmin` 으로 유지한다.
    """
    if res is None or res.empty:
        return
    cols = [c for c in ("측점", "날짜", "시작", "종료", "Kp", "dF", "dD_arcmin",
                        "상태", "파일") if c in res.columns]
    out = res[cols].copy()
    out["날짜"] = out["날짜"].astype(str)
    CORRECTION_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(CORRECTION_CSV, index=False, encoding="utf-8-sig")
    ok = (out["상태"] == "정상").sum() if "상태" in out else len(out)
    print(f"\n[저장] {CORRECTION_CSV}  ({len(out)}행, 정상 {ok}행)")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
