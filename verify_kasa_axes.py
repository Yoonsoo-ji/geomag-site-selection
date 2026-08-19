# -*- coding: utf-8 -*-
"""
KASA 관측소 X·Y·Z·F 좌표축 규약 검증
=======================================

E2(4소 공간보간)의 개선폭이 0.5123° → 0.4987°, 곧 **약 0.8′** 뿐이다. 편차 벡터를
측점 좌표계로 옮길 때 쓰는 회전각이 조금만 틀려도 이 우위는 사라진다. 그래서
회전을 넣기 전에 **원자료의 축 규약 자체를 확인**한다.

판정 방법 — 각 관측소에서 관측값과 IGRF 를 직접 견준다.

    지리축이면   X ≈ X_IGRF,  Y ≈ Y_IGRF  (Y 는 음수 약 −4,400 nT)
    계기축이면   X ≈ H_IGRF,  Y ≈ 0       (센서가 자북에 정렬)
    벡터 정합    F(파일) ≈ √(X²+Y²+Z²)

    python verify_kasa_axes.py
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

import lmm_build as LB                      # noqa: E402
from lmm_external_multi import STATIONS, load_kasa   # noqa: E402
from lmm_cyg import fetch_range             # noqa: E402


def igrf_at(lat, lon, when):
    D, I, F, X, Y, Z = LB.igrf_dif(np.array([lat]), np.array([lon]),
                                   np.array([0.0]), when)
    return dict(D=float(D[0]), I=float(I[0]), F=float(F[0]),
                X=float(X[0]), Y=float(Y[0]), Z=float(Z[0]),
                H=float(np.hypot(X[0], Y[0])))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 92)
    print("KASA·CYG 좌표축 규약 검증 — 관측값 vs IGRF")
    print("=" * 92)

    rows = []
    for key, st in STATIONS.items():
        if st["src"] == "kasa":
            d = load_kasa(key)
        else:
            d = fetch_range(dt.date(2024, 6, 10), dt.date(2024, 6, 19),
                            quiet=True)
            d["time"] = pd.to_datetime(d["time"], utc=True)
        d = d.dropna(subset=["X", "Y", "Z"])
        if d.empty:
            print(f"  ! {st['name']}({key}) 자료 없음")
            continue
        # 정온 표본 — 자료 중앙 부근 하루
        mid = d.iloc[len(d) // 2]["time"]
        when = pd.Timestamp(mid).to_pydatetime().replace(tzinfo=None)
        seg = d[(d.time >= mid - pd.Timedelta("12h"))
                & (d.time <= mid + pd.Timedelta("12h"))]
        if len(seg) < 60:
            seg = d
        X, Y, Z = (float(seg.X.median()), float(seg.Y.median()),
                   float(seg.Z.median()))
        Fcol = float(seg["F"].median()) if "F" in seg else np.nan
        ig = igrf_at(st["lat"], st["lon"], when)

        H_obs = float(np.hypot(X, Y))
        D_obs = float(np.degrees(np.arctan2(Y, X)))
        F_vec = float(np.sqrt(X * X + Y * Y + Z * Z))

        # 규약 판정
        if abs(Y - ig["Y"]) < 500:
            conv = "지리축 (X=북, Y=동)"
        elif abs(Y) < 800 and abs(X - ig["H"]) < 800:
            conv = "계기축 (X≈H, Y≈0)"
        else:
            conv = "판정 보류"

        rows.append(dict(
            관측소=f"{st['name']}({key})", 표본일=str(when.date()),
            X=round(X, 1), Y=round(Y, 1), Z=round(Z, 1),
            X_IGRF=round(ig["X"], 1), Y_IGRF=round(ig["Y"], 1),
            H_IGRF=round(ig["H"], 1),
            D관측=round(D_obs, 2), D_IGRF=round(ig["D"], 2),
            회전각=round(ig["D"] - D_obs, 2),
            F파일=round(Fcol, 1) if Fcol == Fcol else np.nan,
            F벡터=round(F_vec, 1), F_IGRF=round(ig["F"], 1),
            규약=conv))

    df = pd.DataFrame(rows)
    print()
    print(df[["관측소", "표본일", "X", "Y", "Z", "X_IGRF", "Y_IGRF", "H_IGRF"]]
          .to_string(index=False))
    print()
    print(df[["관측소", "D관측", "D_IGRF", "회전각", "F파일", "F벡터", "F_IGRF",
              "규약"]].to_string(index=False))

    print("\n판정")
    for _, r in df.iterrows():
        dF = (r["F파일"] - r["F벡터"]) if r["F파일"] == r["F파일"] else np.nan
        note = ""
        if dF == dF:
            note = (f" · F 파일↔벡터 차 {dF:+.1f} nT"
                    + ("(정합)" if abs(dF) < 5 else " ⚠불일치"))
        print(f"  {r['관측소']:12} {r['규약']}{note}")

    print("\n해석")
    print("  · 「지리축」이면 회전을 넣으면 안 된다(현재 코드는 회전을 넣는다).")
    print("  · 「계기축」이면 회전이 필요하나, IGRF 로 추정하는 것보다 관측소의")
    print("    baseline orientation 을 쓰는 편이 낫다. 그 자료가 있는지 확인할 것.")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_KASA_좌표축_검증.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")
    return df


if __name__ == "__main__":
    main()
