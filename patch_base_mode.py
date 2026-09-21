# -*- coding: utf-8 -*-
r"""카드 안내문을 Base 모드로 고친다 (2026-09-21 발주자 확정)

    python patch_base_mode.py --file "docs/output/…_오호입력_전체50.xlsx"

오호 원문이 **3열(시각·총자력·품질)** — 매뉴얼 부록 B.2 의 Base 모드 형식이고,
발주자가 「base 모드로 잰 것」이라 확인했다. 카드 안내문은 내가 매뉴얼 5.3.1·5.3.4 를
읽고 적은 **Mobile 모드 + X/Y 좌표**라 실제와 어긋나 있었다.

⚠️ 모드가 바뀌면 **파일 형식이 바뀐다** — Mobile 은 6열(시각·라인·스테이션·원시·환원·품질)
이라 ⑨ 의 파서(2번째 토막을 F 로 읽음)가 통째로 어긋난다. 즉 이 안내문은 문구가 아니라
**계산의 전제**다.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

OLD = "자력계는 Mobile 모드, 좌표는 X/Y(P0 = 0,0 · 동·북 + · 서·남 −)로 둡니다."
NEW = ("자력계는 Base 모드로 두고 한 자리에서 자동 반복 판독을 받습니다. 측점을 옮길 때마다 "
       "기록을 끊어 새 블록이 생기게 하고, 거리는 블록 순서로만 남으니 순서를 건너뛰지 마세요.")
NOTE_ROW, NOTE_H = 13, 64.0


def patch(f: Path) -> list[str]:
    import pythoncom
    import win32com.client as win32

    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    xl.ScreenUpdating = False
    done: list[str] = []
    try:
        wb = xl.Workbooks.Open(str(f), UpdateLinks=0)
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            t = str(ws.Range(f"A{NOTE_ROW}").Value or "")
            if OLD in t:
                ws.Range(f"A{NOTE_ROW}").Value = t.replace(OLD, NEW)
                ws.Rows(NOTE_ROW).RowHeight = NOTE_H
                done.append(ws.Name)
        wb.Worksheets("⑨ 구배 자동계산").Activate()
        wb.Save()
        wb.Close(True)
    finally:
        xl.ScreenUpdating = True
        xl.Quit()
        pythoncom.CoUninitialize()
    return done


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    a = ap.parse_args()
    f = Path(a.file).resolve()
    bak = Path(__file__).parent / "_base_mode_backup.xlsx"
    shutil.copy(f, bak)
    done = patch(f)
    print(f"안내문 고친 카드 {len(done)}장 · 되돌릴 사본 {bak.name}")


if __name__ == "__main__":
    main()
