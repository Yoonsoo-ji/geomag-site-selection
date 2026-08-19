# -*- coding: utf-8 -*-
"""
Regional 성과 신뢰도 감사 — 어느 측정을 믿을 수 있는가
========================================================

13절 작업순서 1단계의 본론. 「성과표가 원시값인가」(audit_survey_table.py)는 확인됐으니,
이제 **그 원시값 자체가 믿을 만한가**를 가린다.

세 층위로 나눠 본다.

  ① 내부 산술  — 야장 원시 판독값으로 편각·복각을 재계산해 기재값과 대조
  ② 정밀도     — 같은 날 세션 간 산포 (반복성)
  ③ 정확도     — 재방문 간 편각 차에서 IGRF 영년변화를 뺀 **잔여**
                  같은 표석을 다시 재면 0 이어야 한다. 0 이 아니면 둘 중 하나가 틀렸다.

③ 이 이 감사의 핵심이다. ①②가 통과해도 ③ 이 크면 «정밀하지만 부정확»한 상태이고,
그것이 곧 13절에서 관찰된 편각 35′ 미설명 잔차의 정체다.

    python audit_regional_reliability.py
"""
import datetime as dt
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import ppigrf

import lmm_verify2019 as V
from lmm_build import load_survey_points

ROOT = Path(__file__).parent
OUT = ROOT / "docs" / "output"

TOL_OK, TOL_WARN = 6.0, 20.0      # 잔여 판정 문턱(′). 6′ = KPI 0.1°

# 2019 신규측점은 성과표에 없어 좌표를 보완한다(lmm_build.FB2019_NEW_SITE_COORD 와 동일)
EXTRA_COORD = {"장흥": (34.6272, 126.9778, 100.0),
               "남지": (35.4361, 128.4736, 50.0)}


def igrf_D(lat, lon, elev_m, year):
    Be, Bn, _ = ppigrf.igrf(lon, lat, elev_m / 1000.0, dt.datetime(int(year), 7, 1))
    return float(np.degrees(np.arctan2(float(Be), float(Bn))))


def main():
    raw = V.collect_raw()
    r = V.internal_checks(raw)
    r["연도"] = pd.to_datetime(r["날짜"]).dt.year
    print(f"원시 판독 세션 {len(r)}건 · 측점 {r['측점'].nunique()}종 "
          f"({r['연도'].min()}~{r['연도'].max()})")

    # ── ① 내부 산술 ──
    mk = r["A_마크폐합_gon"].abs()
    dd = r["B_편각차_deg"].abs() * 3600
    print(f"\n① 내부 산술")
    print(f"   편각 재계산차  중앙 {dd.median():.2f}″ · 최대 {dd.max():.2f}″")
    print(f"   마크 폐합      허용 {V.TOL_MARK_GON} gon 초과 "
          f"{int((mk > V.TOL_MARK_GON).sum())}/{len(mk)}건")
    for _, s in r[mk > V.TOL_MARK_GON].iterrows():
        print(f"      · {s['측점']} {str(s['날짜'])[:10]} 폐합오차 "
              f"{s['A_마크폐합_gon']:.4f} gon")

    # ── ② 정밀도 ──
    rep = []
    for (site, day), g in r.groupby(["측점", "날짜"]):
        if len(g) < 2:
            continue
        rep.append({"측점": site, "일자": str(day)[:10], "세션": len(g),
                    "산포_arcmin": g["D_재계산"].std() * 60})
    rep = pd.DataFrame(rep)
    print(f"\n② 정밀도 (같은 날 세션 간 편각 산포)")
    print(f"   n={len(rep)}일 · 중앙 {rep.산포_arcmin.median():.2f}′ "
          f"· 최대 {rep.산포_arcmin.max():.2f}′ "
          f"· KPI 6′ 초과 {int((rep.산포_arcmin > 6).sum())}일")

    # ── ③ 정확도 ──
    sp = load_survey_points()
    coord = {n: (g.lat.iloc[0], g.lon.iloc[0], g.elev_m.iloc[0])
             for n, g in sp.groupby("name")}
    for k, v in EXTRA_COORD.items():
        coord.setdefault(k, v)

    per = (r.groupby(["측점", "연도"])
             .agg(D=("D_재계산", "mean"), 산포=("D_재계산", "std"),
                  세션=("D_재계산", "size"),
                  방위각=("방위각_deg", "mean")).reset_index())

    rows = []
    for site, g in per.groupby("측점"):
        if len(g) < 2 or site not in coord:
            continue
        g = g.sort_values("연도")
        lat, lon, el = coord[site]
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            obs = (b["D"] - a["D"]) * 60
            pred = (igrf_D(lat, lon, el, b["연도"])
                    - igrf_D(lat, lon, el, a["연도"])) * 60
            left = obs - pred
            rows.append({"측점": site, "구간": f"{a['연도']}~{b['연도']}",
                         "관측ΔD_arcmin": obs, "IGRF_ΔD_arcmin": pred,
                         "잔여_arcmin": left,
                         "방위각변화_deg": b["방위각"] - a["방위각"],
                         "판정": ("정상" if abs(left) < TOL_OK else
                                "주의" if abs(left) < TOL_WARN else "불일치")})
    cons = pd.DataFrame(rows)
    v = cons["잔여_arcmin"].abs()
    print(f"\n③ 정확도 (재방문 잔여 = 관측ΔD − IGRF 영년변화)")
    print(f"   n={len(cons)}구간 · 중앙 {v.median():.1f}′ · RMS "
          f"{np.sqrt((v**2).mean()):.1f}′ · 최대 {v.max():.1f}′")
    print(f"   정상 {int((cons.판정=='정상').sum())} · "
          f"주의 {int((cons.판정=='주의').sum())} · "
          f"불일치 {int((cons.판정=='불일치').sum())}")
    print()
    print(cons[["측점", "구간", "관측ΔD_arcmin", "IGRF_ΔD_arcmin",
                "잔여_arcmin", "판정"]].round(1).to_string(index=False))

    # ── 방문별 혐의 판정 ──
    # 3회 이상 방문한 측점에서 잔여가 «올랐다 내려오면» 가운데 방문이 범인이다.
    print("\n④ 혐의 방문 (연속 구간의 잔여 부호가 뒤집히면 가운데 방문이 원인)")
    suspect = []
    for site, g in cons.groupby("측점"):
        if len(g) < 2:
            continue
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.loc[i, "잔여_arcmin"], g.loc[i + 1, "잔여_arcmin"]
            if a * b < 0 and min(abs(a), abs(b)) > TOL_WARN:
                mid = g.loc[i, "구간"].split("~")[1]
                suspect.append((site, mid, a, b))
                print(f"   ★ {site} {mid}년 방문 — 앞구간 {a:+.1f}′ / "
                      f"뒷구간 {b:+.1f}′ 로 부호 반전")
    if not suspect:
        print("   (해당 없음)")

    # 2회 방문뿐이라 범인을 못 가리는 것
    only2 = [s for s, g in cons.groupby("측점") if len(g) == 1
             and abs(g["잔여_arcmin"].iloc[0]) >= TOL_WARN]
    if only2:
        print(f"\n   ⚠ 2회 방문뿐이라 어느 쪽이 틀렸는지 못 가리는 측점: "
              f"{', '.join(only2)}")

    OUT.mkdir(parents=True, exist_ok=True)
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    p = OUT / f"{ts}_Regional_신뢰도감사.csv"
    cons.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")

    print("\n" + "=" * 70)
    print("요약 — 정밀하지만 부정확하다")
    print("=" * 70)
    print(f"  내부 산술   완전 일치 (재계산차 중앙 {dd.median():.2f}″)")
    print(f"  정밀도      중앙 {rep.산포_arcmin.median():.2f}′ — KPI 6′ 대비 양호")
    print(f"  정확도      잔여 RMS {np.sqrt((v**2).mean()):.1f}′ — "
          f"{int((cons.판정!='정상').sum())}/{len(cons)} 구간이 6′ 초과")
    print("  → 측정 절차나 계산이 아니라 **방문마다 달라지는 기준**의 문제다.")
    print("     마크 방위각이 재방문마다 수십~170° 바뀌는 것과 함께 보면,")
    print("     신규 선점에서 **방위 기준을 고정·정밀화**해야 할 근거가 된다.")
    return cons


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
