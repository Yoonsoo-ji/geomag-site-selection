# -*- coding: utf-8 -*-
r"""원문을 지점마다 «그 자리에» 보관한다 — 붙여넣기가 곧 저장 (2026-09-21)

    python add_raw_archive.py --into "docs/output/…_카드자동연결_전체50.xlsx"

발주자 요청: ① 원문(txt)도 이 엑셀에 남기고 ② 지점을 바꾸면 새로 붙여 넣을 수 있고
③ 「저장」을 눌러 그 지점 자료가 보관되게.

⚠️ **버튼은 매크로가 있어야 하는데 이 PC 는 COM 으로 VBA 를 넣지 못한다**
(엑셀 「VBA 프로젝트 개체 모델에 대한 액세스를 신뢰함」이 꺼져 있다 — 보안 설정이라
임의로 켜지 않는다). 그래서 **저장 단계 자체를 없애는 쪽**으로 풀었다.

    「⑩ 원문 보관」 시트에 지점마다 다섯 칸(동·서·남·북·수직)을 미리 만들어 두고,
    ⑨ 계산 시트는 «고른 지점의 칸»을 읽는다.

붙여 넣는 자리가 곧 보관 자리이므로 **붙여 넣는 순간 저장**되고, 지점을 바꿔 돌아와도
원문이 그대로 있어 결과가 다시 살아난다. 카드 값을 굳이 고정하지 않아도 잃지 않는다.

⑨ 의 A9:E708 은 이제 **보관 시트를 읽는 수식**이라 직접 붙여 넣는 칸이 아니다(회색).
「원문 칸으로 가기」 하이퍼링크가 그 지점 자리로 데려다 준다.

같이 고치는 것:
  · 카드 ② 비고 — 방향마다 「읽음/제외 · 확인」이 자동으로 들어간다(발주자 지적)
  · 카드 3. 수직 비고 — 「읽은 시각 · N회 평균 · M회 제외」
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
G, GS = "'⑨ 구배 자동계산'", "⑨ 구배 자동계산"
A, AS_ = "'⑩ 원문 보관'", "⑩ 원문 보관"
BLOCK = 703                  # 지점 한 칸: 이름 1 + 머리글 1 + 자료 700 + 여백 1
PASTE_TOP, PASTE_N = 9, 700
COLS = ["동", "서", "남", "북", "수직(삼각대)"]
START = "$CM$6"              # 고른 지점의 자료 시작행 (숨긴 칸)

GUIDE = (
    "① 오른쪽 위에서 **지점을 고르고**, 바로 아래 「원문 칸으로 가기」를 누르면 그 지점의 "
    "다섯 칸(동·서·남·북·수직)으로 갑니다. 거기에 그 방향 txt 를 «통째로» 붙여 넣으세요"
    "(메모장에서 Ctrl+A → Ctrl+C → 칸 하나만 고르고 Ctrl+V).   "
    "② 붙여 넣는 자리가 곧 보관 자리라 **따로 저장할 필요가 없습니다** — 지점을 바꿨다가 "
    "돌아와도 원문이 그대로 있고 결과도 다시 나옵니다.   "
    "③ 오른쪽 초록 표가 바로 채워지고, 고른 지점 카드에도 저절로 들어갑니다 — 블록을 "
    "시각순으로 줄 세워 P0(전)·1 m·…·P0(후) 로 배정하고, 0.00 읽음과 품질 뒷자리 0, "
    "중앙값에서 ±5 nT 벗어난 값을 빼고 평균을 냅니다.   "
    "④ 아래 A9:E708 은 보관 칸을 비춰 보는 곳이라 여기에 붙여 넣지 않습니다."
)
GUIDE4 = ('=IF($B$4="","지점을 고르면 그 지점 원문 칸과 카드가 이어집니다 — 먼저 골라 주세요",'
          '"「"&$B$4&"」 원문은 ⑩ 보관 시트에 남고, 결과는 같은 이름 카드에 자동으로 들어갑니다")')
LINK = ('=IF($B$4="","▶ 먼저 지점을 골라 주세요",'
        f'HYPERLINK("#{A}!A"&{START},'
        '"▶ 「"&$B$4&"」 원문 칸으로 가기 — 거기에 붙여 넣으면 이 표가 채워지고 그대로 보관됩니다"))')


def build(target: Path, out: Path) -> tuple[int, int]:
    import pythoncom
    import win32com.client as win32
    from win32com.client import constants as _c  # noqa: F401

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    xl.ScreenUpdating = False
    cards = 0
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        g = wb.Worksheets(GS)
        sites = [g.Cells(9 + i, 91).Value for i in range(50)]

        # ── ⑩ 원문 보관 ────────────────────────────────────────
        for s in wb.Worksheets:
            if s.Name == AS_:
                s.Delete()
        arc = wb.Worksheets.Add(After=g)
        arc.Name = AS_
        arc.Tab.Color = 0xD8B48C
        for i, col in enumerate("ABCDE", 1):
            arc.Columns(col).ColumnWidth = 38
        arc.Columns("F").ColumnWidth = 22
        arc.Range("A1:F1").Merge()
        arc.Range("A1").Value = ("⑩ 원문 보관 — 지점마다 다섯 칸. 여기에 붙여 넣는 것이 곧 저장입니다"
                                 " (⑨ 시트가 고른 지점의 칸을 읽어 계산합니다)")
        arc.Range("A1").Font.Bold = True
        arc.Range("A1").Font.Size = 12
        arc.Range("A1").Interior.Color = 0x4F2F1A
        arc.Range("A1").Font.Color = 0xFFFFFF
        arc.Rows(1).RowHeight = 24
        for i, name in enumerate(sites, 1):
            h = 2 + (i - 1) * BLOCK
            arc.Range(f"A{h}:F{h}").Merge()
            arc.Range(f"A{h}").Value = f"{name}  —  아래 다섯 칸에 그 방향 txt 를 통째로 붙여 넣습니다"
            arc.Range(f"A{h}").Font.Bold = True
            arc.Range(f"A{h}").Interior.Color = 0xE8D5B7
            arc.Rows(h).RowHeight = 21
            for j, lab in enumerate(COLS):
                c = arc.Cells(h + 1, 1 + j)
                c.Value = f"{lab} 원문"
                c.Font.Bold = True
                c.Font.Color = 0xFFFFFF
                c.Interior.Color = 0x6B4423
                c.HorizontalAlignment = -4108
            arc.Cells(h + 1, 6).Formula = f'=HYPERLINK("#{G}!A1","◂ ⑨ 계산 시트로")'
            arc.Range(f"A{h + 2}:E{h + 4}").Interior.Color = 0xCCF2FF   # 붙일 자리 세 줄만 노랗게
        arc.Rows(3).Select()
        arc.Range("A1").Select()

        # ── ⑨ — 보관 칸을 읽는다 ───────────────────────────────
        g.Range("CM6").Formula = (f'=IFERROR(MATCH($B$4,$CM$9:$CM$58,0)*{BLOCK}-{BLOCK - 4},0)')
        idx = f"{START}+ROW()-{PASTE_TOP}"
        pick = f"INDEX({A}!$A:$E,{idx},COLUMN())"
        g.Range(f"A{PASTE_TOP}:E{PASTE_TOP + PASTE_N - 1}").Formula = (
            f'=IF({START}=0,"",IF({pick}="","",{pick}))')
        for j, lab in enumerate(COLS):
            c = g.Cells(8, 1 + j)
            c.Value = f"{lab} 원문 (⑩ 보관에서 읽음)"
        g.Range(f"A{PASTE_TOP}:E{PASTE_TOP + 2}").Interior.Color = 0xE8E8E8
        g.Range("A2").Value = GUIDE
        g.Range("C4").Formula = GUIDE4
        g.Range("A6:O6").Merge()
        g.Range("A6").Formula = LINK
        g.Range("A6").Font.Bold = True
        g.Range("A6").Font.Size = 11
        g.Range("A6").Interior.Color = 0xCCF2FF
        g.Rows(6).RowHeight = 22

        # ── 카드 비고 ─────────────────────────────────────────
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            gate = f'{G}!$B$4<>"{ws.Name}"'
            ws.Range("H31:I34").Formula = (
                f'=IF({gate},"",IF({G}!$M30="","",{G}!$M30&"  ·  "&{G}!$O30))')
            ws.Range("G42:G45").Formula = (
                f'=IF({gate},"",IF({G}!$J23="","",{G}!$J23&"  ·  "&{G}!$L23))')
            ws.Range("A14").Formula = (
                '="① 방향별 관측자료"&IF(' + f'{G}!$B$4="{ws.Name}"'
                + ',"    ←  ⑨ 시트에서 자동 입력 중 · 원문은 ⑩ 보관에 남습니다","")')
            cards += 1
        wb.Worksheets(GS).Activate()
        g.Range("B4").Select()
        wb.Save()
        wb.Close(True)
    finally:
        xl.ScreenUpdating = True
        xl.Quit()
        pythoncom.CoUninitialize()
    return len(sites), cards


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--into", required=True)
    ap.add_argument("--tag", default="원문보관_전체50")
    a = ap.parse_args()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    n, cards = build(Path(a.into).resolve(), out)
    print(f"보관 칸 {n}지점 · 비고 붙인 카드 {cards}장")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
