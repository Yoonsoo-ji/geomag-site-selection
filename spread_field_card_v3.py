# -*- coding: utf-8 -*-
r"""3판 「5. 방위표지」 구역을 50개 카드 전부에 펴는 패치 (2026-09-21)

    python spread_field_card_v3.py --into "docs/output/20260918_112540_…_태백시험.xlsx"

9-18 작업에서 **01_태백 한 장만** 3판으로 바뀌어 있었다 — 안내문 줄을 빼고, 「측정 방법」을
드롭다운으로 바꾸고, 자리표시 문구(`[ 도로(추가요망) ]` 등)를 지운 것이다. 나머지 49장은
9-18 오전 판(현장개선) 그대로라 두 판이 한 파일 안에 섞여 있었다.

⚠️ **이 구역은 지점과 무관한 칸뿐이다**(표지 위도·경도·측정 방법·시준거리·비고). 그래서
태백의 `A54:I57` 을 엑셀 범위 복사로 그대로 얹는다 — 서식·병합·드롭다운이 함께 간다.
행 높이는 범위 복사로 따라오지 않으므로 따로 맞춘다.

⚠️ **행을 넣거나 지우지 않는다.** 카드 수식·인쇄 행바꿈(38·58·106)이 행 번호에 물려 있어
한 줄만 밀려도 3쪽 구성이 흐트러진다. 두 판의 행 번호가 같은 것을 확인하고 값만 바꾼다.

⚠️ 9-17 2판에 있던 **관측일자 메모(주석)와 숨긴 열 그룹(1·2·＋ 표시)** 은 9-18 판에서 이미
빠져 있다 — 이 스크립트는 남아 있는지 확인만 하고, 있으면 지운다.
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
SRC_SHEET = "01_태백 신규점"
BLOCK = "A54:I57"                    # 5. 방위표지 — 안내문 줄 ~ 방위표지2
ROW_H = {54: 6.0, 55: 25.05, 56: 27.0, 57: 27.0}
# ⚠️ 태백 드롭다운의 「후대GNSS」는 「휴대GNSS」 오타로 보여 바로잡는다(H50 목록과 같은 표기).
METHOD_LIST = "RTK-GNSS,네트워크RTK,정적GNSS,휴대GNSS,기존 성과 이용,보도로 측정"
XL_VALIDATE_LIST, XL_BETWEEN = 3, 1


def is_card(ws) -> bool:
    return ws.Range("A48").Text.startswith("4. 기준점 좌표")


def patch(target: Path, out: Path) -> dict:
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    report = {"패치": [], "건너뜀": [], "메모 지움": [], "열그룹 푼 시트": []}
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        src = wb.Worksheets(SRC_SHEET)
        for ws in wb.Worksheets:
            if not is_card(ws):
                report["건너뜀"].append(ws.Name)
                continue
            if ws.Name != SRC_SHEET:
                src.Range(BLOCK).Copy(ws.Range("A54"))
                for r, h in ROW_H.items():
                    ws.Rows(r).RowHeight = h
                dv = ws.Range("F56:F57").Validation
                dv.Delete()
                dv.Add(Type=XL_VALIDATE_LIST, AlertStyle=1, Operator=XL_BETWEEN,
                       Formula1=METHOD_LIST)
                dv.IgnoreBlank = True
                dv.InCellDropdown = True
                report["패치"].append(ws.Name)
            # 2판 잔재 — 메모(주석)와 숨긴 열 그룹
            if ws.Range("A7").Comment is not None:
                ws.Range("A7").Comment.Delete()
                report["메모 지움"].append(ws.Name)
            if ws.Columns("K").OutlineLevel > 1 or ws.Columns("AG").OutlineLevel > 1:
                ws.Cells.ClearOutline()
                report["열그룹 푼 시트"].append(ws.Name)
        # 태백 드롭다운도 같은 목록으로 맞춘다(오타 정정)
        dv = src.Range("F56:F57").Validation
        dv.Delete()
        dv.Add(Type=XL_VALIDATE_LIST, AlertStyle=1, Operator=XL_BETWEEN, Formula1=METHOD_LIST)
        dv.IgnoreBlank = True
        dv.InCellDropdown = True
        xl.CutCopyMode = False
        wb.Worksheets(SRC_SHEET).Activate()
        wb.Save()
        wb.Close(True)
    finally:
        xl.Quit()
        pythoncom.CoUninitialize()
    return report


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--into", required=True, help="3판 야장 xlsx (원본은 건드리지 않음)")
    ap.add_argument("--tag", default="전체50")
    a = ap.parse_args()
    target = Path(a.into).resolve()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    rep = patch(target, out)
    for k, v in rep.items():
        print(f"{k}: {len(v)}건 {v[:4]}{' …' if len(v) > 4 else ''}")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
