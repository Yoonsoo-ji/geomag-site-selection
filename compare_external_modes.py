# -*- coding: utf-8 -*-
"""
External 적용 방식 A/B — end-to-end 재적합으로 판정
====================================================

13절 작업순서 2단계. 층이 결합돼 있어 순차 튜닝이 금지되므로, 모드마다
**처음부터 다시 적합**해 LOO 로 비교한다.

  none           : 미적용 (기준선)
  drop_storm     : 교란시(Kp>2) 성과를 표본에서 배제 — 값을 만지지 않는다
  subtract       : 전 세션 평균 외부장 편차를 뺀다
  subtract_quiet : 정온(Kp<=2) 세션 평균만으로 뺀다

⚠ 균일 V 근사가 남아 있다(CYG 한 곳의 변동을 전국 동일 적용). 2019 표본에서는
   이 근사가 편각을 악화시켰으므로, 개선이 나오더라도 표본 수를 함께 볼 것.

    python compare_external_modes.py
"""
import importlib
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd

MODES = ["none", "drop_storm", "subtract", "subtract_quiet", "subtract_D"]
# 보정량 자료원 — 청양 단독 vs 관측소 4소 공간보간(자문 권고)
SOURCES = ["cyg", "multi"]
QUALITY = ["none", "strict", "conservative"]   # Regional 신뢰도 필터
DEGREE = 1   # 생산 파이프라인이 LOO 로 채택하는 차수


_CRUSTAL = None
_CVEC = None


def crustal():
    """KIGAM 격자는 모드와 무관하므로 한 번만 읽는다."""
    global _CRUSTAL, _CVEC
    if _CRUSTAL is None:
        import lmm_build
        lons, lats, grid = lmm_build.load_kigam_grid()
        _CRUSTAL = lmm_build.CrustalGrid(lons, lats, grid)
        # 지각 벡터도 모드와 무관하므로 함께 한 번만 만든다
        bx, by, bz = lmm_build.crustal_vector(lons, lats, grid)
        _CVEC = lmm_build.CrustalVector(lons, lats, bx, by, bz)
    return _CRUSTAL


def run(mode, quality="none", source=None):
    import lmm_build
    importlib.reload(lmm_build)
    lmm_build.EXTERNAL_MODE = mode
    lmm_build.REGIONAL_FILTER = quality
    if source:
        lmm_build.EXTERNAL_SOURCE = source
        csv = (lmm_build.EXTERNAL_CSV_MULTI if source == "multi"
               else lmm_build.EXTERNAL_CSV_CYG)
        lmm_build.EXTERNAL_CSV = csv
    pts = lmm_build.load_all_points(include_2019=lmm_build.INCLUDE_2019)
    # main() 과 같은 순서: IGRF 잔차 -> 측점 집계 -> LOO
    res = lmm_build.igrf_residuals(pts)
    sites = lmm_build.aggregate_sites(pts, res)
    # 생산 설정과 같게 지각 벡터 기여를 붙인다 (없으면 열이 0 이라 종전 동작)
    crustal()
    sites = lmm_build.attach_crustal_di(
        sites, _CVEC if lmm_build.CRUSTAL_VECTOR else None)
    # 생산 파이프라인은 LOO 로 차수를 고르며 degree 1 이 채택된다.
    # 모드 비교도 같은 차수로 해야 본 모델 수치와 맞댈 수 있다.
    loo = lmm_build.loo_cv(sites, crustal(), degree=DEGREE)
    return dict(mode=mode, quality=quality, source=source or "cyg",
                n_obs=len(pts), n=len(sites),
                **{k: loo[k] for k in ("D", "I", "F")})


def main():
    rows = []
    for m in MODES:
        srcs = ["cyg"] if m in ("none", "drop_storm") else SOURCES
        for src in srcs:
            try:
                rows.append(run(m, source=src))
                print(f"  [{m}/{src}] 완료")
            except Exception as e:      # noqa: BLE001
                print(f"  [{m}/{src}] 실패: {type(e).__name__} {e}")
    df = pd.DataFrame(rows)
    if df.empty:
        return
    base = df[df["mode"] == "none"].iloc[0]
    print("\n" + "=" * 68)
    print(f"{'모드':16} {'자료원':>6} {'표본':>5} {'LOO D(°)':>10} "
          f"{'LOO I(°)':>10} {'LOO F(nT)':>11}")
    print("-" * 68)
    for _, r in df.iterrows():
        mark = ""
        if r["mode"] != "none":
            dD = r["D"] - base["D"]
            dF = r["F"] - base["F"]
            mark = f"   D {dD:+.4f} · F {dF:+.2f}"
        print(f"{r['mode']:16} {r['source']:>6} {int(r['n']):5} {r['D']:10.4f} "
              f"{r['I']:10.4f} {r['F']:11.2f}{mark}")
    print("=" * 68)
    best = df.loc[df["F"].idxmin()]
    print(f"\nF 최소 모드: {best['mode']}/{best['source']}  "
          f"(D {best['D']:.4f}° · F {best['F']:.2f} nT)")
    bd = df.loc[df["D"].idxmin()]
    print(f"D 최소 모드: {bd['mode']}/{bd['source']}  "
          f"(D {bd['D']:.4f}° · F {bd['F']:.2f} nT)")
    print("⚠ 개선폭이 잔차 수준(D 0.35° · F 60 nT)에 비해 작으면 잡음이다.")

    out = Path(__file__).parent / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    import datetime as dt
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_External_모드비교.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"[저장] {p}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
