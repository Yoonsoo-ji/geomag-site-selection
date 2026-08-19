# -*- coding: utf-8 -*-
"""
야장 F(전자력) 기록 전수 감사 — 관측시각과 센서 위치차
========================================================

2020 지구물리측량 연구보고서의 절차는 단순한 F 한 값이 아니다.

    ① 편각·복각 측정점에서 **5 m 이상 떨어져** 총자기장 측정
    ② 편복각 측정 후 **절대측정기 센서 위치에서 다시** 총자기장 측정
    ③ 두 위치의 차이를 보정

야장 「※ 전자력 측정 결과」 절에 `측정값`·`기본값`·`측정일시`·반복측정표가 있다.
이것이 위 두 지점의 값인지, 그리고 시각이 실제로 기입돼 있는지를 전수로 센다.

    python audit_fieldbook_F.py
"""
import datetime as dt
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
from lmm_fieldbook import ngii_files      # noqa: E402

LBL = {"head": "전자력 측정 결과", "meas": "자력측정", "time": "측정일시",
       "res": "결과", "avg": "전자력(평균)"}


def _s(v):
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v).strip()


def scan_sheet(df):
    """한 시트에서 F 절의 값들을 뽑는다."""
    out = []
    nrow, ncol = df.shape
    for r in range(nrow):
        c0 = _s(df.iat[r, 0])
        if LBL["head"] not in c0:
            continue
        rec = {"row": r}
        for rr in range(r, min(r + 8, nrow)):
            lab = _s(df.iat[rr, 0])
            vals = [_s(df.iat[rr, c]) for c in range(1, min(ncol, 12))]
            vals = [v for v in vals if v]
            if LBL["meas"] in lab:
                nums = [v for v in vals if re.fullmatch(r"\d{4,6}(\.\d+)?", v)]
                if len(nums) >= 1:
                    rec["측정값"] = float(nums[0])
                if len(nums) >= 2:
                    rec["기본값"] = float(nums[1])
            elif LBL["time"] in lab:
                rec["시각원문"] = " | ".join(vals[:4])
            elif lab == LBL["res"]:
                nums = [v for v in vals if re.fullmatch(r"\d{4,6}(\.\d+)?", v)]
                if nums:
                    rec["결과"] = float(nums[0])
        out.append(rec)
    return out


def has_real_time(txt):
    """측정일시 칸에 실제 시각이 있는가(00:00 · 날짜만은 결측으로 본다)."""
    if not txt:
        return False
    if re.search(r"\b([01]?\d|2[0-3]):[0-5]\d", txt):
        return not re.search(r"\b00:00(:00)?\b", txt) or \
               len(set(re.findall(r"\b([01]?\d|2[0-3]):[0-5]\d", txt))) > 1
    return False


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 78)
    print("야장 F(전자력) 기록 감사 — 시각·두 지점 측정 여부")
    print("=" * 78)

    rows = []
    files = [Path(p) for p in ngii_files()]
    print(f"대상 야장 {len(files)}개")
    for i, f in enumerate(files, 1):
        try:
            xl = pd.ExcelFile(f)
        except Exception:
            continue
        for sh in xl.sheet_names:
            try:
                df = xl.parse(sh, header=None, dtype=object)
            except Exception:
                continue
            for rec in scan_sheet(df):
                rec.update(파일=f.name, 시트=sh)
                rows.append(rec)
        if i % 40 == 0:
            print(f"  … {i}/{len(files)}")

    if not rows:
        print("F 절을 찾지 못했다.")
        return

    d = pd.DataFrame(rows)
    for c in ("측정값", "기본값", "결과"):
        if c not in d:
            d[c] = np.nan
    if "시각원문" not in d:
        d["시각원문"] = ""
    d["시각원문"] = d["시각원문"].fillna("")
    d["시각있음"] = d["시각원문"].map(has_real_time)
    d["위치차_nT"] = d["측정값"] - d["기본값"]

    print(f"\nF 절 {len(d)}건 발견 (파일 {d['파일'].nunique()}개)")
    print(f"  측정값 기재     {int(d['측정값'].notna().sum())}건")
    print(f"  기본값 기재     {int(d['기본값'].notna().sum())}건")
    both = d["측정값"].notna() & d["기본값"].notna()
    diff = both & (d["위치차_nT"].abs() > 0.05)
    print(f"  두 값 모두 기재 {int(both.sum())}건 · "
          f"**서로 다른 값** {int(diff.sum())}건")
    print(f"  측정일시 실기재 {int(d['시각있음'].sum())}건")

    if diff.any():
        print("\n두 값이 다른 사례 (위치차 후보)")
        show = d[diff][["파일", "시트", "측정값", "기본값", "위치차_nT"]].head(15)
        print(show.to_string(index=False))
        print(f"\n  위치차 |Δ| 중앙 {d.loc[diff, '위치차_nT'].abs().median():.1f} nT · "
              f"최대 {d.loc[diff, '위치차_nT'].abs().max():.1f} nT")
    else:
        print("\n⚠ 측정값과 기본값이 다른 사례가 없다 —")
        print("  이 야장 양식의 「기본값」은 두 번째 지점의 측정값이 아니라")
        print("  같은 값을 옮겨 적은 칸으로 보인다. 즉 **센서 위치차 보정에 필요한")
        print("  두 지점 F 는 야장에 남아 있지 않다.**")

    if not d["시각있음"].any():
        print("\n⚠ F 측정일시가 실제 시각으로 기입된 사례가 없다 —")
        print("  현재 F 는 세션 전체 시간창으로 대신 보정할 수밖에 없다.")

    out = ROOT / "docs" / "output"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{dt.datetime.now():%Y%m%d_%H%M%S}_야장_F기록_감사.csv"
    d.to_csv(p, index=False, encoding="utf-8-sig")
    print(f"\n[저장] {p}")
    return d


if __name__ == "__main__":
    main()
