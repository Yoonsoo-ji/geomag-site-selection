# -*- coding: utf-8 -*-
r"""카드 ② 결과표에 「잰 거리」 열을 넣는다 (2026-09-21)

    python add_measured_dist_col.py --into "docs/output/…_원문보관_전체50.xlsx"

발주자 지적: 「P0→끝」의 «끝»이 어디인지 표에 나와야 한다. 방향 이름에 붙여 두었더니
빈 카드에서는 보이지 않아, ⑨ 미리보기처럼 **전용 열**로 뺀다.

    방향 | 잰 거리 | P0 전후차 | F 변화폭 | 최대 구간구배 | 발생구간 | P0→끝 ΔF | P0→끝 평균구배 | 비고
      A       B          C           D            E             F           G             H            I

기존 B~H 를 한 칸씩 오른쪽으로 옮기고 B 를 새로 만든다. ⚠️ **수식을 «복사»하면 안 된다**
— 표 안에서 서로를 가리키는 참조(발생구간→최대구배, 평균구배→ΔF)가 함께 밀려 엉뚱한
칸을 가리키게 된다. 서식만 붙여넣고 수식은 다시 적는다.
"""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"
XL_FORMATS = -4122
XL_CENTER = -4108
ROWS = {31: ("동", "C", "O", "P"), 32: ("서", "E", "T", "U"),
        33: ("남", "G", "Y", "Z"), 34: ("북", "I", "AD", "AE")}
HEAD = ["방향", "잰 거리", "P0 전후차\n(nT)", "관측구간\nF 변화폭 (nT)",
        "최대 구간\n수평구배\n(nT/m)", "발생구간", "P0→끝\nΔF (nT)",
        "P0→끝\n평균구배\n(nT/m)", "비고"]
TAIL = ("「잰 거리」는 그 방향에서 마지막으로 값을 적은 거리이고, 「P0→끝」 두 칸은 "
        "그 거리까지의 차이와 평균입니다.")


def formulas(r: int, gsheet: str, name: str) -> dict[str, str]:
    d, f, g1, g2 = ROWS[r]
    rng = f"{f}17:{f}26"
    last = f'LOOKUP(2,1/({rng}<>""),{rng})'
    dist = f'LOOKUP(2,1/({rng}<>""),ROW({rng})-16)'
    q = f"'{gsheet}'"
    return {
        "A": d,
        "B": f'=IF(COUNT({rng})=0,"",{dist}&" m")',
        "C": f'=IF(COUNT({f}16,{f}27)<2,"",{f}27-{f}16)',
        "D": f'=IF(COUNT({f}16:{f}26)=0,"",MAX({f}16:{f}26)-MIN({f}16:{f}26))',
        "E": f'=IF(COUNT({g1}17:{g1}26)=0,"",MAX({g1}17:{g1}26))',
        "F": f'=IF(E{r}="","",INDEX({g2}17:{g2}26,MATCH(E{r},{g1}17:{g1}26,0)))',
        "G": f'=IF(OR(COUNT({f}16)=0,COUNT({rng})=0),"",{last}-{f}16)',
        "H": f'=IF(G{r}="","",ROUND(ABS(G{r})/{dist},2))',
        "I": (f'=IF({q}!$B$4<>"{name}","",IF({q}!$M{r - 1}="","",'
              f'{q}!$M{r - 1}&"  ·  "&{q}!$O{r - 1}))'),
    }


def patch(target: Path, out: Path) -> list[str]:
    import pythoncom
    import win32com.client as win32

    shutil.copy(target, out)
    pythoncom.CoInitialize()
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    xl.DisplayAlerts = False
    xl.AskToUpdateLinks = False
    xl.ScreenUpdating = False
    done: list[str] = []
    try:
        wb = xl.Workbooks.Open(str(out), UpdateLinks=0)
        gname = "⑨ 구배 자동계산"
        for ws in wb.Worksheets:
            if ws.Range("A15").Text != "거리" or not ws.Range("A29").Text.startswith("②"):
                continue
            for a in ("H30:I30", "H31:I34"):
                if ws.Range(a).MergeCells:
                    ws.Range(a).UnMerge()
            ws.Range("B30:H34").Copy()
            ws.Range("C30:I34").PasteSpecial(Paste=XL_FORMATS)
            xl.CutCopyMode = False
            # ⚠️ 이 엑셀은 COM 으로 영어 `NumberFormat` 을 거부한다 — 로컬 이름을 쓴다
            ws.Range("B30:B34").NumberFormatLocal = "G/표준"
            ws.Range("B31:B34").HorizontalAlignment = XL_CENTER
            for i, h in enumerate(HEAD):
                ws.Cells(30, 1 + i).Value = h
            for r in range(31, 35):
                for col, f in formulas(r, gname, ws.Name).items():
                    cell = ws.Range(f"{col}{r}")
                    if f.startswith("="):
                        cell.Formula = f
                    else:
                        cell.Value = f
            note = ws.Range("A37").Value or ""
            head = note.split("「P0→끝」")[0].rstrip()
            ws.Range("A37").Value = (head + " " + TAIL).strip()
            done.append(ws.Name)
        wb.Worksheets(gname).Activate()
        wb.Worksheets(gname).Range("B4").Select()
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
    ap.add_argument("--into", required=True)
    ap.add_argument("--tag", default="잰거리열_전체50")
    a = ap.parse_args()
    out = OUT_DIR / f"{dt.datetime.now():%Y%m%d_%H%M%S}_시험탐사_현장야장_3판_{a.tag}.xlsx"
    done = patch(Path(a.into).resolve(), out)
    print(f"결과표 고친 카드 {len(done)}장: {done[:3]} …")
    print(f"[저장] {out}  ({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
