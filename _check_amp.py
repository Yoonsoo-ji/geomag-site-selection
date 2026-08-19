# -*- coding: utf-8 -*-
"""관측소별 일변동 진폭 비교 — 제주 절대값 불일치가 「변동량」에도 영향을 주는가.

E2 는 절대값이 아니라 **편차**만 쓴다. 그래도 어느 관측소의 감도(스케일)가
다르면 편차 자체가 왜곡된다. 공통일의 일변동 진폭을 견줘 확인한다.
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from lmm_external_multi import STATIONS, load_kasa, KST


def daily_amp(d):
    d = d.dropna(subset=["X", "Y", "Z"]).copy()
    d["day"] = d.time.dt.tz_convert(KST).dt.date
    g = d.groupby("day").agg(X=("X", lambda v: v.max() - v.min()),
                             Y=("Y", lambda v: v.max() - v.min()),
                             Z=("Z", lambda v: v.max() - v.min()),
                             n=("X", "size"))
    return g[g.n > 1200]


amps = {}
for key, st in STATIONS.items():
    if st["src"] != "kasa":
        continue
    d = load_kasa(key)
    if not d.empty:
        amps[st["name"]] = daily_amp(d)

common = set.intersection(*[set(a.index) for a in amps.values()])
print(f"KASA 3소 공통일 {len(common)}일 ({min(common)} ~ {max(common)})")

rows = []
for name, a in amps.items():
    a = a[a.index.isin(common)]
    rows.append(dict(관측소=name, 일수=len(a),
                     X진폭=round(a.X.median(), 1),
                     Y진폭=round(a.Y.median(), 1),
                     Z진폭=round(a.Z.median(), 1)))
df = pd.DataFrame(rows)
print(df.to_string(index=False))

b = df[df["관측소"] == "이천"].iloc[0]
print("\n이천(지리축) 대비 비율")
for _, r in df.iterrows():
    print(f"  {r['관측소']}: X {r['X진폭']/b['X진폭']:.2f}  "
          f"Y {r['Y진폭']/b['Y진폭']:.2f}  Z {r['Z진폭']/b['Z진폭']:.2f}")
print("\n비율이 1 에서 크게 벗어나면 그 관측소의 감도(스케일)가 다르다는 뜻이고,")
print("편차만 쓰는 E2 에서도 보정량이 그만큼 왜곡된다.")
