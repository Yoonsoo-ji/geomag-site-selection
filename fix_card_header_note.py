# -*- coding: utf-8 -*-
r"""카드 머리글이 50장 전부 「S01 태백 신규점」이었다 · 안내문 한 줄이 잘렸다 (2026-09-21)

    python fix_card_header_note.py --into "docs/output/…_카드자동연결_전체50.xlsx"

① **머리글** — 3차 검토 M7 로 「여러 쪽으로 나뉘면 2쪽부터 지점을 알 수 없다」를 고치며
   오른쪽 머리글에 Site ID 를 넣었는데, 50장을 한꺼번에 낼 때 그 값이 복제돼 **어느 카드를
   인쇄해도 2~4쪽 머리에 「S01 태백 신규점」이 찍혔다.** 카드 제목(A1)에서 그 카드의
   이름을 읽어 각자 제 이름을 갖게 한다.

② **2. 수평 자기구배 안내문(A13)** — 세 줄이 필요한데 행 높이가 두 줄치라 마지막 줄
   (「시각은 HHMM 숫자 4자리로…」)이 인쇄에서 잘렸다. 높이를 늘린다.
"""
from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
NOTE_ROW, NOTE_H = 13, 52.0


def site_label(title: str) -> str:
    """'[S19] 오호 신규점 — 지자기 …' → 'S19 오호 신규점'"""
    head = title.split("—")[0].strip()
    return re.sub(r"^\[(.+?)\]\s*", r"\1 ", head).strip()


def fix(target: Path, out: Path) -> list[tuple[str, str]]:
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    done: list[tuple[str, str]] = []
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            lab = site_label(str(ws.Range("A1").Value or ""))
            if lab:
                ws.PageSetup.RightHeader = lab
            ws.Rows(NOTE_ROW).RowHeight = NOTE_H
            done.append((ws.Name, lab))
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
    ap.add_argument("--tag", default="머리글정정_전체50")
    a = ap.parse_args()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    done = fix(Path(a.into).resolve(), out)
    print(f"고친 카드 {len(done)}장 — 예: {done[:3]}")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
