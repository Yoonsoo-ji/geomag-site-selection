# -*- coding: utf-8 -*-
r"""카드 ② 결과표의 「P0→10 m」을 «마지막으로 잰 거리»로 바꾼다 (2026-09-21)

    python patch_card_lastdist.py --into "docs/output/…_구배자동계산.xlsx"

⚠️ 종전 수식은 `C26`(10 m 줄)만 봤다. 오호처럼 동 6 m · 남 4 m 로 끝난 방향은 그 칸이
통째로 비어 P0→끝 변화량을 못 냈다(발주자 지적: 「최대 몇 m 까지 갔는지 확인해서 항목을
변경해야 한다」). 이제 **마지막으로 값이 있는 거리**를 찾아 그 거리까지의 ΔF 와 평균구배를
낸다 — 방향마다 간 거리가 달라도 각각 맞게 나온다.

`LOOKUP(2,1/(범위<>""),범위)` 는 마지막 비어 있지 않은 칸을 집는 옛 관용식이라 엑셀 판을
가리지 않는다. 거리(m)는 같은 방법으로 행 번호에서 얻는다(17행 = 1 m … 26행 = 10 m).
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
DIR_ROWS = {31: "C", 32: "E", 33: "G", 34: "I"}      # 동·서·남·북 결과 행 ↔ 평균 F 열
P0_ROW, FIRST_M, LAST_M = 16, 17, 26
NOTE_ROW = 37


def patch(target: Path, out: Path) -> list[str]:
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    done: list[str] = []
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        for ws in wb.Worksheets:
            if ws.Range("A29").Text != "② 방향별 수평구배 결과 (자동 계산)":
                continue
            ws.Range("F30").Value = "P0→끝\nΔF (nT)"
            ws.Range("G30").Value = "P0→끝\n평균구배\n(nT/m)"
            for r, col in DIR_ROWS.items():
                rng = f"{col}{FIRST_M}:{col}{LAST_M}"
                last = f'LOOKUP(2,1/({rng}<>""),{rng})'
                dist = f'LOOKUP(2,1/({rng}<>""),ROW({rng})-{P0_ROW})'
                ws.Range(f"F{r}").Formula = (
                    f'=IF(OR(COUNT({col}{P0_ROW})=0,COUNT({rng})=0),"",'
                    f'{last}-{col}{P0_ROW})')
                ws.Range(f"G{r}").Formula = (
                    f'=IF(F{r}="","",ROUND(ABS(F{r})/{dist},2))')
            note = ws.Range(f"A{NOTE_ROW}").Value or ""
            tail = ("「P0→끝」은 그 방향에서 «마지막으로 잰 거리»까지의 차이와 그 평균입니다 — "
                    "6 m 까지 갔으면 6 m, 10 m 까지 갔으면 10 m 가 기준입니다.")
            if tail not in note:
                ws.Range(f"A{NOTE_ROW}").Value = (note.rstrip() + " " + tail).strip()
            done.append(ws.Name)
        wb.Worksheets(3).Activate()
        wb.Save()
        wb.Close(True)
    finally:
        xl.Quit()
        pythoncom.CoUninitialize()
    return done


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--into", required=True)
    ap.add_argument("--tag", default="구배자동계산_전체50")
    a = ap.parse_args()
    target = Path(a.into).resolve()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    done = patch(target, out)
    print(f"결과표 고친 카드 {len(done)}장: {done[:3]} …")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
