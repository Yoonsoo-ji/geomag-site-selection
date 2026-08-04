# -*- coding: utf-8 -*-
"""
청양(CYG) 관측소 자료 수집 — INTERMAGNET HAPI
==============================================

LMM ④ External 층 입력자료를 INTERMAGNET HAPI 서버에서 내려받아 캐시한다.

    서버: https://imag-data.bgs.ac.uk/GIN_V1/hapi
    자료셋 ID 구성: <IAGA>/<등급>/<주기>/<성분>

등급별 CYG 보유 기간 (2026-07 확인)
-----------------------------------
    definitive   2014-01-01 ~ 2017-12-31   <- 측정기간(2022~2025) 미포함
    quasi-def    2015-04-25 ~ 2020-12-31   <- 측정기간 미포함
    adjusted     2013-10-30 ~ 현재
    best-avail   2013-10-30 ~ 현재         <- 일자별 최상위 등급 자동 선택

따라서 2022~2025 성과 보정에는 definitive 를 쓸 수 없고 best-avail 을 쓴다.
등급이 섞이므로 보고서에는 사용 등급을 반드시 명시해야 한다.

사용 예:
    python lmm_cyg.py --start 2024-06-01 --end 2024-06-30
    python lmm_cyg.py --start 2022-01-01 --end 2025-12-31   # 전체(느림)

    from lmm_cyg import fetch_range
    df = fetch_range("2024-06-10", "2024-06-12")
"""

import argparse
import csv
import datetime as dt
import io
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).parent
CACHE = BASE / "data" / "cyg"
HAPI = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"

DATASET = "cyg/best-avail/PT1M/xyzf"
FILL = 99999.0
EXPECTED_PER_DAY = 1440  # 1분 자료

# 관측소 제원 (INTERMAGNET 등록값)
CYG_LAT, CYG_LON = 36.37, 126.85


def best_dataset(day: dt.date) -> str:
    """
    해당 일자에 사용할 최상위 등급 자료셋.

    definitive 는 2017 년까지, quasi-def 는 2020 년까지만 제공되므로
    일자에 따라 등급이 달라진다. 보고서에 등급을 명시할 수 있도록
    선택 결과를 그대로 자료셋 ID 로 돌려준다.
    """
    if day <= dt.date(2017, 12, 31):
        return "cyg/definitive/PT1M/xyzf"
    if dt.date(2015, 4, 25) <= day <= dt.date(2020, 12, 31):
        return "cyg/quasi-def/PT1M/xyzf"
    return DATASET  # best-avail


def _day_path(day: dt.date, ds: str) -> Path:
    grade = ds.split("/")[1]
    return CACHE / grade / f"{day:%Y}" / f"cyg_{day:%Y%m%d}.csv"


def fetch_day(day: dt.date, force=False, dataset=None) -> pd.DataFrame:
    """하루치 1분 자료를 받아 캐시하고 DataFrame 으로 반환."""
    ds = dataset or best_dataset(day)
    p = _day_path(day, ds)
    if p.exists() and not force:
        return pd.read_csv(p, parse_dates=["time"])

    stop = day + dt.timedelta(days=1)
    q = urllib.parse.urlencode({
        "dataset": ds,
        "parameters": "Field_Vector,Field_Magnitude",
        "start": f"{day}T00:00:00Z",
        "stop": f"{stop}T00:00:00Z",
    })
    # HAPI 는 400(자료 없음/등급 경계)·IncompleteRead(응답 잘림)·타임아웃 등 다양한
    # 오류를 낸다. 어떤 오류든 해당 일만 건너뛰고(빈 DataFrame, 캐시 안 함) 계속한다.
    try:
        with urllib.request.urlopen(f"{HAPI}/data?{q}", timeout=180) as r:
            text = r.read().decode()
    except Exception as e:                          # noqa: BLE001 (네트워크 불안정 대응)
        print(f"    ! {day} 건너뜀: {type(e).__name__} {e}", file=sys.stderr)
        return pd.DataFrame(columns=["time", "X", "Y", "Z", "F"])

    rows = []
    for rec in csv.reader(l for l in text.strip().split("\n") if l):
        if len(rec) < 5:
            continue
        rows.append((rec[0], *map(float, rec[1:5])))

    df = pd.DataFrame(rows, columns=["time", "X", "Y", "Z", "F"])
    if len(df):
        df["time"] = pd.to_datetime(df["time"], format="ISO8601", utc=True)
        # 결측 표식을 NaN 으로 (99999 를 실수값으로 오인하면 치명적)
        df[["X", "Y", "Z", "F"]] = df[["X", "Y", "Z", "F"]].where(
            df[["X", "Y", "Z", "F"]] < FILL)

    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)
    return df


def fetch_range(start, end, quiet=False, dataset=None) -> pd.DataFrame:
    """[start, end] 구간(양끝 포함)의 1분 자료를 이어붙여 반환."""
    s = pd.Timestamp(start).date()
    e = pd.Timestamp(end).date()
    out = []
    n = (e - s).days + 1
    for k in range(n):
        d = s + dt.timedelta(days=k)
        out.append(fetch_day(d, dataset=dataset))
        if not quiet and (k % 10 == 0 or k == n - 1):
            print(f"  [{k+1:4d}/{n}] {d}", flush=True)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# ------------------------------------------------------------------ Kp 지수
KP_URL = "https://kp.gfz.de/fileadmin/files_for_gfz_cms/Kp_ap_since_1932.txt"
KP_CACHE = CACHE / "Kp_ap_since_1932.txt"


def fetch_kp(force=False) -> pd.DataFrame:
    """
    GFZ Kp/ap 지수(1932~현재, 3시간 간격)를 받아 DataFrame 으로 반환.

    형식: YYYY MM DD hh.h hh._m days days_m Kp ap D
    D=1 이면 확정치(definitive), 0 이면 잠정치.
    """
    KP_CACHE.parent.mkdir(parents=True, exist_ok=True)
    if not KP_CACHE.exists() or force:
        with urllib.request.urlopen(KP_URL, timeout=180) as r:
            KP_CACHE.write_bytes(r.read())

    rows = []
    for ln in KP_CACHE.read_text(errors="replace").split("\n"):
        if not ln or ln.startswith("#"):
            continue
        p = ln.split()
        if len(p) < 9:
            continue
        rows.append((int(p[0]), int(p[1]), int(p[2]), float(p[3]),
                     float(p[7]), int(p[8]), int(p[9]) if len(p) > 9 else 0))

    df = pd.DataFrame(rows, columns=["y", "m", "d", "h", "Kp", "ap", "definitive"])
    df["time"] = (pd.to_datetime(df[["y", "m", "d"]].rename(
        columns={"y": "year", "m": "month", "d": "day"}), utc=True)
        + pd.to_timedelta(df["h"], unit="h"))
    return df[["time", "Kp", "ap", "definitive"]]


def kp_for(times, kp: pd.DataFrame = None) -> pd.Series:
    """임의 시각들에 해당하는 3시간 구간의 Kp 값을 돌려준다."""
    if kp is None:
        kp = fetch_kp()
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    idx = kp["time"].searchsorted(t, side="right") - 1
    idx = idx.clip(0, len(kp) - 1)
    return pd.Series(kp["Kp"].to_numpy()[idx], index=range(len(t)))


def derive(df: pd.DataFrame) -> pd.DataFrame:
    """X·Y·Z 로부터 편각 D·복각 I·수평분력 H 를 유도."""
    d = df.copy()
    d["H"] = np.hypot(d.X, d.Y)
    d["D_deg"] = np.degrees(np.arctan2(d.Y, d.X))
    d["I_deg"] = np.degrees(np.arctan2(d.Z, d.H))
    return d


def daily_report(df: pd.DataFrame) -> pd.DataFrame:
    """
    일별 완전성과 외부장 변동폭.

    변동폭(peak-to-peak)은 "그날 측정했다면 시각에 따라 성과가 얼마나
    달라졌을지"의 상한이며, 곧 외부장 보정으로 제거 가능한 오차의 규모다.
    """
    d = derive(df)
    d["date"] = d.time.dt.date
    g = d.groupby("date")

    rep = pd.DataFrame({
        "수신": g.size(),
        "유효": g["F"].count(),
        "F변동_nT": g["F"].apply(lambda s: s.max() - s.min()),
        "D변동_arcmin": g["D_deg"].apply(lambda s: (s.max() - s.min()) * 60),
        "I변동_arcmin": g["I_deg"].apply(lambda s: (s.max() - s.min()) * 60),
    }).reset_index()
    rep["완전성_%"] = (100 * rep["유효"] / EXPECTED_PER_DAY).round(1)
    return rep


def main():
    ap = argparse.ArgumentParser(description="CYG 1분 자료 수집 (INTERMAGNET HAPI)")
    ap.add_argument("--start", required=True, help="시작일 YYYY-MM-DD")
    ap.add_argument("--end", required=True, help="종료일 YYYY-MM-DD (포함)")
    ap.add_argument("--report", action="store_true", help="일별 요약만 출력")
    a = ap.parse_args()

    print(f"CYG 자료 수집: {a.start} ~ {a.end}  (자료셋 {DATASET})")
    df = fetch_range(a.start, a.end)
    if df.empty:
        print("수신 자료 없음")
        return

    rep = daily_report(df)
    ok = rep[rep["유효"] > 0]
    bad = rep[rep["완전성_%"] < 90]

    print(f"\n총 {len(rep)}일 / 유효자료 있는 날 {len(ok)}일")
    if len(bad):
        print(f"\n[주의] 완전성 90% 미만 {len(bad)}일 — 해당 일자 측정성과는 보정 불가:")
        print(bad[["date", "유효", "완전성_%"]].to_string(index=False))

    if len(ok):
        print("\n--- 외부장 일변동 (보정으로 제거 가능한 오차 규모) ---")
        print(f"  F : 중앙값 {ok['F변동_nT'].median():.1f} nT   "
              f"최대 {ok['F변동_nT'].max():.1f} nT   (KPI 50 nT)")
        print(f"  D : 중앙값 {ok['D변동_arcmin'].median()/60:.3f}°  "
              f"최대 {ok['D변동_arcmin'].max()/60:.3f}°   (KPI 0.1°)")

    if a.report:
        print("\n--- 일별 상세 ---")
        print(rep.to_string(index=False))

    print(f"\n캐시 위치: {CACHE}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
