# -*- coding: utf-8 -*-
"""야장 원본에서 F(자력) 관련 기록을 전수 훑는다 — 원격점/센서위치 재측정·시각·위치차."""
import warnings; warnings.filterwarnings("ignore")
import re, pandas as pd
from pathlib import Path
from collections import Counter
from lmm_fieldbook import ngii_files      # 지리원 야장 원본 순회

KEY = re.compile(r"자력|전자력|총자|F값|GSM|프로톤|오버하우저|양성자|센서|이격|거리|위치")
cnt = Counter(); samples = {}
nfile = 0
for _p in ngii_files():
    path = Path(_p)
    nfile += 1
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        continue
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh, header=None, dtype=object)
        except Exception:
            continue
        for r in range(df.shape[0]):
            for c in range(min(df.shape[1], 3)):
                v = df.iat[r, c]
                if isinstance(v, str) and KEY.search(v):
                    lab = v.strip()[:22]
                    cnt[lab] += 1
                    if lab not in samples:
                        row = [str(df.iat[r, k])[:14] for k in range(min(df.shape[1], 10))]
                        samples[lab] = (path.name[:28], sh[:12], r, row)
    if nfile >= 25:
        break
print(f"훑은 파일 {nfile}개\n")
print("자력·센서·위치 관련 라벨 빈도")
for lab, n in cnt.most_common(25):
    f, sh, r, row = samples[lab]
    print(f"  {n:4d}  {lab:24} | {f} [{sh}] r{r}")
    print(f"        {row}")
