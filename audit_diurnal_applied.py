# -*- coding: utf-8 -*-
"""
성과표에 일변화 보정이 적용됐는가 — 상시관측 자료로 직접 판정
================================================================

작업규정 제21조는 「동일시간 상시 기준점과 비교」해 일변화·영년변화를 보정하도록
정한다. 성과표가 이미 보정된 값이면 우리가 다시 빼면 **이중보정**이다.

판정 원리 — 관측 시각의 외부장 편차를 그대로 «되풀이»하는지 본다.

    보정 안 됐다면   잔차 = (측점 계통오차) + (그 시각의 외부장 편차) + 잡음
                    → 편차를 설명변수로 회귀하면 **기울기 β ≈ 1**
    보정 됐다면      외부장 몫이 이미 빠져 있으므로 **β ≈ 0**

⚠ 그냥 회귀하면 검정력이 없다. 잔차는 측점마다 고정된 계통오차(약 30′)가 지배하고
   외부장 편차는 6′ 수준이라 묻힌다. 그래서 **측점 안에서 평균을 뺀다**(within-station
   차분). 측점 고정 오차가 소거되고 시각에 따른 변동만 남아 검정력이 올라간다.

    python audit_diurnal_applied.py
"""
import datetime as dt
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
import lmm_build as LB      # noqa: E402

EXT_CSV = ROOT / "docs" / "data" / "external_corrections_multi.csv"


def load_ext():
    d = pd.read_csv(EXT_CSV, encoding="utf-8-sig")
    d = d[d["상태"] == "정상"].copy()
    d["연도"] = pd.to_datetime(d["날짜"]).dt.year
    g = (d.groupby(["측점", "연도"])
           .agg(dD_ext=("dD_arcmin", "mean"), dI_ext=("dI_arcmin", "mean"),
                dF_ext=("dF", "mean"), n=("dF", "size"),
                Kp=("Kp", "max")).reset_index())
    return g


def within_station(df, cols):
    """측점 평균을 빼 측점 고정 오차를 소거한다."""
    out = df.copy()
    for c in cols:
        out[c + "_w"] = out[c] - out.groupby("station")[c].transform("mean")
    return out


def fit(x, y):
    """β 와 표준오차·상관. 표본이 적으므로 CI 를 함께 본다."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4 or np.allclose(x, 0):
        return dict(n=int(len(x)))
    b = float((x @ y) / (x @ x))
    resid = y - b * x
    se = float(np.sqrt((resid @ resid) / max(len(x) - 1, 1) / (x @ x)))
    r = float(np.corrcoef(x, y)[0, 1])
    return dict(n=int(len(x)), beta=b, se=se, r=r,
                lo=b - 1.96 * se, hi=b + 1.96 * se,
                sd_x=float(np.std(x)), sd_y=float(np.std(y)))


def session_test():
    """
    [3] 같은 날 여러 세션 — 검정력이 가장 높은 판정.

    한 측점·한 날의 세션들은 **측점 계통오차와 그날의 기준값을 정확히 공유**한다.
    따라서 세션 간 차이는 오로지 그 시각의 외부장(과 측정 잡음)이다.

        ⚠ 야장 D 는 **정의상 원시 판독값**이다(우리 코드가 ED/WD/EU/WU 에서 재계산).
      따라서 이 β 는 「보정 여부」가 아니라 **우리 외부장 추정이 실제 세션 간 변동을
      얼마나 맞히는가**를 잰다. β≈1 이면 잘 맞히는 것이고, β≈0 이면 추정이 실제와
      무관하다는 뜻이다.
    """
    import lmm_verify2019 as V

    print("\n[3] 같은 날 여러 세션 — 측점·일자 고정 (검정력 최상)")
    raw = V.collect_raw()
    if raw is None or raw.empty:
        print("  야장 원시 판독값을 읽지 못했다.")
        return
    raw = raw.dropna(subset=["D_야장", "시작"]).copy()
    raw["날짜"] = pd.to_datetime(raw["날짜"]).dt.date.astype(str)
    raw["시작"] = raw["시작"].astype(str).str.slice(0, 8)

    ext = pd.read_csv(EXT_CSV, encoding="utf-8-sig")
    ext = ext[ext["상태"] == "정상"].copy()
    ext["날짜"] = ext["날짜"].astype(str)
    ext["시작"] = ext["시작"].astype(str).str.slice(0, 8)

    j = raw.merge(ext[["측점", "날짜", "시작", "dD_arcmin", "dI_arcmin", "dF"]],
                  on=["측점", "날짜", "시작"], how="inner")
    j["D_am"] = j["D_야장"] * 60
    j["I_am"] = j["I_야장"] * 60
    j["key"] = j["측점"] + "|" + j["날짜"]
    g = j.groupby("key").filter(lambda x: len(x) >= 2)
    if g.empty:
        print("  같은 날 2세션 이상인 표본이 없다.")
        return
    print(f"  표본 {len(g)}세션 / {g['key'].nunique()}개 (측점·일자) 묶음")

    for lab, xc, yc, unit in (("편각 D", "dD_arcmin", "D_am", "′"),
                              ("복각 I", "dI_arcmin", "I_am", "′")):
        x = (g[xc] - g.groupby("key")[xc].transform("mean")).to_numpy(float)
        y = (g[yc] - g.groupby("key")[yc].transform("mean")).to_numpy(float)
        f = fit(x, y)
        if "beta" not in f:
            continue
        b, lo, hi = f["beta"], f["lo"], f["hi"]
        # 야장 값은 정의상 원시 판독값이다(우리 코드가 ED/WD/EU/WU 로 재계산).
        # 따라서 β≈1 이 나와야 «우리 외부장 추정이 실제 변동을 맞힌다»는 뜻이고,
        # β≈0 은 보정이 됐다는 뜻이 아니라 **추정이 실제와 무관하다**는 뜻이다.
        if lo <= 1 <= hi and not (lo <= 0 <= hi):
            v = "예측이 실제 변동을 잘 맞힌다 (β=1 포함 · 0 배제)"
        elif lo <= 0 <= hi and not (lo <= 1 <= hi):
            v = "⚠ 예측이 실제 세션 간 변동을 설명하지 못한다 (β=0 포함 · 1 배제)"
        else:
            v = "판별 불가"
        print(f"  {lab}: n={f['n']:3d}  β={b:+.2f} [{lo:+.2f}, {hi:+.2f}]  "
              f"r={f['r']:+.2f}  예측 편차 산포 {f['sd_x']:.2f}{unit} · "
              f"실제 세션간 산포 {f['sd_y']:.2f}{unit}")
        print(f"          → {v}")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 80)
    print("성과표 일변화 보정 적용 여부 — 상시관측 편차와의 회귀로 판정")
    print("=" * 80)

    # External 을 걸지 않은 «원시» 잔차가 필요하다
    LB.EXTERNAL_MODE = "none"
    pts = LB.load_all_points(include_2019=LB.INCLUDE_2019)
    res = LB.igrf_residuals(pts)
    d = pts[["name", "station", "year", "lat", "lon"]].copy()
    d[["rD", "rI", "rF"]] = res[["dD", "dI", "dF"]].values
    d["rD"] *= 60      # 분
    d["rI"] *= 60
    d["year"] = d["year"].astype(int)

    ext = load_ext()
    m = d.merge(ext, left_on=["station", "year"], right_on=["측점", "연도"],
                how="inner")
    print(f"매칭 {len(m)}행 / 측점 {m['station'].nunique()}개 "
          f"({m['year'].min()}~{m['year'].max()})")
    if m.empty:
        return

    rep = m.groupby("station").filter(lambda g: len(g) >= 2)
    print(f"재방문 2회 이상 측점만: {len(rep)}행 / {rep['station'].nunique()}측점")

    print("\n[1] 전체 표본 회귀 (측점 고정 오차가 섞여 검정력이 낮다)")
    print(f"{'성분':>6} {'n':>4} {'β':>8} {'95% CI':>18} {'r':>7}")
    for lab, xc, yc in (("편각 D", "dD_ext", "rD"), ("복각 I", "dI_ext", "rI"),
                        ("총자력 F", "dF_ext", "rF")):
        f = fit(m[xc].to_numpy(float), m[yc].to_numpy(float))
        if "beta" not in f:
            continue
        print(f"{lab:>6} {f['n']:>4} {f['beta']:>8.2f} "
              f"[{f['lo']:>7.2f}, {f['hi']:>6.2f}] {f['r']:>7.2f}")

    print("\n[2] 측점 내 차분 회귀 — 측점 고정 오차를 소거한 판정")
    w = within_station(rep, ["dD_ext", "dI_ext", "dF_ext", "rD", "rI", "rF"])
    print(f"{'성분':>6} {'n':>4} {'β':>8} {'95% CI':>18} {'r':>7} "
          f"{'편차 산포':>9} {'잔차 산포':>9}")
    verdict = {}
    for lab, xc, yc, unit in (("편각 D", "dD_ext_w", "rD_w", "′"),
                              ("복각 I", "dI_ext_w", "rI_w", "′"),
                              ("총자력 F", "dF_ext_w", "rF_w", " nT")):
        f = fit(w[xc].to_numpy(float), w[yc].to_numpy(float))
        if "beta" not in f:
            print(f"{lab:>6} {f.get('n', 0):>4}  표본 부족")
            continue
        print(f"{lab:>6} {f['n']:>4} {f['beta']:>8.2f} "
              f"[{f['lo']:>7.2f}, {f['hi']:>6.2f}] {f['r']:>7.2f} "
              f"{f['sd_x']:>8.2f}{unit} {f['sd_y']:>8.2f}{unit}")
        verdict[lab] = f

    print("\n판정")
    for lab, f in verdict.items():
        b, lo, hi = f["beta"], f["lo"], f["hi"]
        if lo <= 1 <= hi and not (lo <= 0 <= hi):
            v = "보정 안 됨 (β=1 을 포함하고 0 은 배제)"
        elif lo <= 0 <= hi and not (lo <= 1 <= hi):
            v = "보정 됨 (β=0 을 포함하고 1 은 배제)"
        elif lo <= 0 <= hi and lo <= 1 <= hi:
            v = "판별 불가 (구간이 0 과 1 을 모두 포함 — 검정력 부족)"
        else:
            v = "이례적 (구간이 0·1 을 모두 배제)"
        print(f"  {lab}: β={b:+.2f} → {v}")

    session_test()

    print("\n※ 이 검정은 야장 대조(audit_survey_table.py)와 독립이다. 야장은 "
          "성과표가\n   세션평균 원시값임을 직접 보였고, 이 검정은 그 결론을 "
          "상시관측 자료로 다시 확인한다.")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_일변화보정_적용여부.csv"
    w.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")
    return w


if __name__ == "__main__":
    main()
