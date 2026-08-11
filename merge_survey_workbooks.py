# -*- coding: utf-8 -*-
"""
회신 조사서 워크북 병합 → 배포 양식(총괄 목록 + 카드 103장) 재구성
====================================================================

11개 본부가 회신한 조사서(각 = 총괄 목록 + 채운 카드 시트들)를 원본 배포 양식
`20260707_조사카드_v2_총괄_전체103건.xlsx` 과 **동일한 구조**로 한 파일에 합친다.

  · 총괄 목록: 원본 템플릿의 총괄 목록(103행)을 그대로 복사
  · 카드 시트: 각 회신본의 카드를 DS 번호순으로 복사 —
    값·서식(글꼴·채움·테두리·정렬·표시형식)·병합셀·행높이·열너비·
    임베드 사진·드롭다운(데이터 검증)·하이퍼링크까지 보존

openpyxl 은 워크북 간 시트 복사를 지원하지 않으므로 셀 단위로 충실 복사한다.

사용:
    python merge_survey_workbooks.py            # 전체 11파일
    python merge_survey_workbooks.py --test     # 소형 파일만(1·2·11) 미리보기
"""
import argparse
import copy
import io
import re
import sys
from datetime import datetime
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.drawing.image import Image as XLImage

ROOT = Path(__file__).parent
SRC_DIR = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구\20260811_사전 현장 조사서 취합"
               r"\2.조사서_부산울산")
TEMPLATE = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구\20260707_지자기 사전 현장조사 엑셀"
                r"\20260707_조사카드_v2_총괄_전체103건.xlsx")
OUT_DIR = ROOT / "docs" / "output"


def copy_sheet(src_ws, dst_wb, title):
    """워크북 간 시트 충실 복사."""
    dst = dst_wb.create_sheet(title)

    # 열 너비 / 숨김
    for k, dim in src_ws.column_dimensions.items():
        d = dst.column_dimensions[k]
        if dim.width is not None:
            d.width = dim.width
        d.hidden = dim.hidden
    # 행 높이
    for k, dim in src_ws.row_dimensions.items():
        if dim.height is not None:
            dst.row_dimensions[k].height = dim.height

    # 셀 값 + 스타일 + 하이퍼링크
    for row in src_ws.iter_rows():
        for c in row:
            if c.value is None and not c.has_style and not c.hyperlink:
                continue
            nc = dst.cell(row=c.row, column=c.column, value=c.value)
            if c.has_style:
                nc.font = copy.copy(c.font)
                nc.fill = copy.copy(c.fill)
                nc.border = copy.copy(c.border)
                nc.alignment = copy.copy(c.alignment)
                nc.protection = copy.copy(c.protection)
                nc.number_format = c.number_format
            if c.hyperlink:
                nc.hyperlink = copy.copy(c.hyperlink)

    # 병합 셀
    for mc in list(src_ws.merged_cells.ranges):
        dst.merge_cells(str(mc))

    # 임베드 이미지 (지도·사진)
    for img in getattr(src_ws, "_images", []):
        try:
            ni = XLImage(io.BytesIO(img._data()))
            ni.anchor = img.anchor
            dst.add_image(ni)
        except Exception as e:   # noqa: BLE001
            print(f"    ! 이미지 복사 실패({title}): {e}", file=sys.stderr)

    # 데이터 검증(드롭다운)
    try:
        for dv in src_ws.data_validations.dataValidation:
            ndv = copy.copy(dv)
            dst.add_data_validation(ndv)
            ndv.sqref = dv.sqref
    except Exception as e:   # noqa: BLE001
        print(f"    ! 데이터검증 복사 실패({title}): {e}", file=sys.stderr)

    # 시트 뷰 속성
    if src_ws.freeze_panes:
        dst.freeze_panes = src_ws.freeze_panes
    if src_ws.auto_filter.ref:
        dst.auto_filter.ref = src_ws.auto_filter.ref
    dst.sheet_view.showGridLines = src_ws.sheet_view.showGridLines
    if src_ws.sheet_properties.tabColor:
        dst.sheet_properties.tabColor = src_ws.sheet_properties.tabColor
    return dst


def card_index(name):
    m = re.match(r"(\d+)_", name)
    return int(m.group(1)) if m else 9999


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--dir", default=str(SRC_DIR))
    a = ap.parse_args()

    files = sorted(Path(a.dir).glob("*.xlsx"),
                   key=lambda p: int(re.match(r"(\d+)", p.name).group(1))
                   if re.match(r"(\d+)", p.name) else 999)
    if a.test:
        files = [f for f in files if re.match(r"(\d+)", f.name)
                 and int(re.match(r"(\d+)", f.name).group(1)) in (1, 2, 11)]

    out = Workbook()
    out.remove(out.active)

    # ── 총괄 목록: 템플릿에서 복사 ──
    print("총괄 목록 복사(템플릿)...")
    tpl = load_workbook(TEMPLATE)
    copy_sheet(tpl["총괄 목록"], out, "총괄 목록")
    tpl.close()

    # ── 카드 시트: 파일별로 열고 닫으며 복사(메모리 보호), 이후 DS 순 재정렬 ──
    total = 0
    for fi, f in enumerate(files, 1):
        print(f"[{fi}/{len(files)}] {f.name} 여는 중...")
        wb = load_workbook(f)
        names = [s for s in wb.sheetnames if s != "총괄 목록"]
        for name in names:
            copy_sheet(wb[name], out, name)
            total += 1
        wb.close()
        del wb
        print(f"    카드 {len(names)}장 복사 (누적 {total})")

    # DS 번호순 정렬 (총괄 목록은 맨 앞 고정)
    head = out["총괄 목록"]
    rest = sorted((ws for ws in out.worksheets if ws is not head),
                  key=lambda ws: card_index(ws.title))
    out._sheets = [head] + rest
    print(f"카드 {total}장 병합 완료")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = "미리보기" if a.test else f"전체{total}건"
    dst = OUT_DIR / f"{ts}_현장조사_취합_총괄_{tag}.xlsx"
    dst.parent.mkdir(parents=True, exist_ok=True)
    print("저장 중(대용량, 수 분 소요 가능)...")
    out.save(dst)
    mb = dst.stat().st_size / 1e6
    print(f"\n저장: {dst}  ({mb:.1f} MB, 시트 {len(out.sheetnames)}개)")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
