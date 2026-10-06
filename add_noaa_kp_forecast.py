"""Add live NOAA 3-day and 27-day Kp forecasts to the standard fieldbook.

The workbook remains macro-free. Excel legacy web queries refresh the two NOAA
plain-text products on open or through Data > Refresh All. KST observations in
sheet ⑧ are converted to UTC by subtracting nine hours before lookup.
"""

from __future__ import annotations

import argparse
import shutil
from datetime import datetime
from pathlib import Path

import fitz
import pythoncom
import win32com.client


XL_OPEN_XML_WORKBOOK = 51
XL_OVERWRITE_CELLS = 0
XL_WEB_SELECTION_ENTIRE_PAGE = 1
XL_WEB_FORMATTING_NONE = 1
XL_CENTER = -4108
XL_LEFT = -4131
XL_RIGHT = -4152
XL_BOTTOM = -4107
XL_CONTINUOUS = 1
XL_THIN = 2
XL_COLUMN_CLUSTERED = 51
XL_CATEGORY = 1
XL_VALUE = 2

NOAA_3DAY = "https://services.swpc.noaa.gov/text/3-day-forecast.txt"
NOAA_27DAY = "https://services.swpc.noaa.gov/text/27-day-outlook.txt"
OBS_SHEET = "⑧ 관측조건·변동관측"
KP_SHEET = "⑧a Kp 예보"


def rgb(hex_color: str) -> int:
    value = hex_color.lstrip("#")
    r, g, b = int(value[:2], 16), int(value[2:4], 16), int(value[4:], 16)
    return r + (g << 8) + (b << 16)


def add_web_query(sheet, cell: str, name: str, url: str):
    query = sheet.QueryTables.Add(Connection=f"URL;{url}", Destination=sheet.Range(cell))
    query.Name = name
    query.WebSelectionType = XL_WEB_SELECTION_ENTIRE_PAGE
    query.WebFormatting = XL_WEB_FORMATTING_NONE
    query.RefreshStyle = XL_OVERWRITE_CELLS
    query.RefreshOnFileOpen = True
    query.BackgroundQuery = False
    query.SaveData = True
    query.AdjustColumnWidth = False
    query.PreserveFormatting = True
    query.EnableRefresh = True
    query.Refresh(False)
    return query


def style_header(cell_range, fill: str = "#16324F") -> None:
    cell_range.Interior.Color = rgb(fill)
    cell_range.Font.Color = rgb("#FFFFFF")
    cell_range.Font.Bold = True
    cell_range.HorizontalAlignment = XL_CENTER
    cell_range.VerticalAlignment = XL_CENTER
    cell_range.WrapText = True
    cell_range.Borders(XL_BOTTOM).LineStyle = XL_CONTINUOUS
    cell_range.Borders(XL_BOTTOM).Weight = XL_THIN
    cell_range.Borders(XL_BOTTOM).Color = rgb("#A8B6C7")


def style_label(cell_range) -> None:
    cell_range.Interior.Color = rgb("#EAF1F8")
    cell_range.Font.Color = rgb("#16324F")
    cell_range.Font.Bold = True
    cell_range.VerticalAlignment = XL_CENTER


def configure_kp_sheet(sheet) -> None:
    sheet.Cells.Font.Name = "Aptos"
    sheet.Cells.Font.Size = 10
    sheet.Activate()
    sheet.Application.ActiveWindow.DisplayGridlines = False

    sheet.Range("A1").Value = "NOAA Kp 예보 조회"
    sheet.Range("A1").Font.Size = 15
    sheet.Range("A1").Font.Bold = True
    sheet.Range("A1").Font.Color = rgb("#16324F")
    sheet.Range("A1:S1").Borders(XL_BOTTOM).LineStyle = XL_CONTINUOUS
    sheet.Range("A1:S1").Borders(XL_BOTTOM).Color = rgb("#1D4ED8")
    sheet.Range("A1:S1").Borders(XL_BOTTOM).Weight = XL_THIN

    sheet.Range("A2").Value = (
        "⑧ 시트 C열에 탐사일시를 한국시간(KST)으로 입력하면 E·F열이 자동 계산됩니다. "
        "인터넷 연결 후 ‘데이터 > 모두 새로 고침’을 실행하십시오."
    )
    sheet.Range("A2:S2").Merge()
    sheet.Range("A2").Font.Color = rgb("#44546A")
    sheet.Range("A2").WrapText = True
    sheet.Rows(2).RowHeight = 33

    sheet.Range("A4:A9").Value = [
        ["조회 날짜·시각 (KST)"],
        ["UTC 환산"],
        ["3일 예보 Kp (3시간)"],
        ["27일 예보 Kp (UTC 일 최대)"],
        ["사용 예보"],
        ["운영 원칙"],
    ]
    style_label(sheet.Range("A4:A9"))
    sheet.Range("B4").Interior.Color = rgb("#FFF2CC")
    sheet.Range("B4").Value = datetime.now().replace(second=0, microsecond=0)
    sheet.Range("B4").NumberFormat = "yyyy-mm-dd hh:mm"
    sheet.Range("B4").Validation.Delete()
    sheet.Range("B4").Validation.Add(Type=4, AlertStyle=1, Operator=1, Formula1="2020-01-01", Formula2="2100-12-31")
    sheet.Range("B4").Validation.InputTitle = "한국시간 입력"
    sheet.Range("B4").Validation.InputMessage = "예: 2026-10-14 13:24"
    sheet.Range("B5").Formula = '=IF(B4="","",B4-9/24)'
    sheet.Range("B5").NumberFormat = "yyyy-mm-dd hh:mm"
    sheet.Range("B6").Formula = '=IF(B4="","",IF(AND(B5>=$A$12,B5<$B$35),INDEX($E$12:$E$35,MATCH(B5,$A$12:$A$35,1)),"범위 밖"))'
    sheet.Range("B7").Formula = '=IF(B4="","",IFERROR(INDEX($I$12:$I$38,MATCH(INT(B5),$H$12:$H$38,0)),"범위 밖"))'
    sheet.Range("B8").Formula = '=IF(B4="","",IF(ISNUMBER(B6),"3일 예보",IF(ISNUMBER(B7),"27일 예보","예보범위 밖")))'
    sheet.Range("B9").Value = "Kp·예보는 기록용입니다. 배제 기준은 아직 확정되지 않았습니다."
    sheet.Range("B6:B7").NumberFormat = "0.00"
    sheet.Range("B4:B9").VerticalAlignment = XL_CENTER
    sheet.Range("B9").WrapText = True

    sheet.Range("D4:D10").Value = [
        ["3일 예보 발행 (UTC)"],
        ["27일 예보 발행 (UTC)"],
        ["새로 고침"],
        ["NOAA 3일 원문"],
        ["NOAA 27일 원문"],
        ["시간 기준"],
        ["인터넷 미연결"],
    ]
    style_label(sheet.Range("D4:D10"))
    sheet.Range("E4").Formula = '=IF($U$3="","미수신",MID($U$3,9,24))'
    sheet.Range("E5").Formula = '=IF($W$3="","미수신",MID($W$3,9,24))'
    sheet.Range("E6").Value = "데이터 > 모두 새로 고침"
    sheet.Hyperlinks.Add(Anchor=sheet.Range("E7"), Address=NOAA_3DAY, TextToDisplay="3-day-forecast.txt")
    sheet.Hyperlinks.Add(Anchor=sheet.Range("E8"), Address=NOAA_27DAY, TextToDisplay="27-day-outlook.txt")
    sheet.Range("E9").Value = "KST = UTC+9 · 조회 시 KST에서 9시간을 뺍니다."
    sheet.Range("E10").Value = "저장된 마지막 예보가 남을 수 있으므로 발행시각을 확인합니다."
    sheet.Range("E9:E10").WrapText = True

    sheet.Range("A11:E11").Value = [["UTC 시작", "UTC 종료", "KST 시작", "KST 종료", "Kp 예보"]]
    style_header(sheet.Range("A11:E11"))
    sheet.Range("H11:K11").Value = [["UTC 날짜", "일 최대 Kp", "A index", "F10.7"]]
    style_header(sheet.Range("H11:K11"))

    months = '{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"}'
    base_date = f'DATE(VALUE(MID($U$3,10,4)),MATCH(MID($U$3,15,3),{months},0),VALUE(MID($U$3,19,2)))'
    sheet.Range("A12").Formula = f'=IF($U$3="","",{base_date}+INT((ROW()-12)/8)+MOD(ROW()-12,8)*3/24)'
    sheet.Range("A12:A35").FillDown()
    sheet.Range("B12").Formula = '=IF(A12="","",A12+3/24)'
    sheet.Range("B12:B35").FillDown()
    sheet.Range("C12").Formula = '=IF(A12="","",A12+9/24)'
    sheet.Range("C12:C35").FillDown()
    sheet.Range("D12").Formula = '=IF(B12="","",B12+9/24)'
    sheet.Range("D12:D35").FillDown()
    utc_slots = 'CHOOSE(MOD(ROW()-12,8)+1,"00-03UT*","03-06UT*","06-09UT*","09-12UT*","12-15UT*","15-18UT*","18-21UT*","21-00UT*")'
    sheet.Range("E12").Formula = (
        f'=IFERROR(VALUE(MID(INDEX($U$2:$U$70,MATCH({utc_slots},$U$2:$U$70,0)),'
        '15+13*INT((ROW()-12)/8),4)),"")'
    )
    sheet.Range("E12:E35").FillDown()
    sheet.Range("A12:D35").NumberFormat = "yyyy-mm-dd hh:mm"
    sheet.Range("E12:E35").NumberFormat = "0.00"

    date27 = f'DATE(VALUE(LEFT($W13,4)),MATCH(MID($W13,6,3),{months},0),VALUE(MID($W13,10,2)))'
    sheet.Range("H12").Formula = f'=IF($W13="","",{date27})'
    sheet.Range("H12:H38").FillDown()
    sheet.Range("I12").Formula = '=IFERROR(VALUE(MID($W13,42,2)),"")'
    sheet.Range("I12:I38").FillDown()
    sheet.Range("J12").Formula = '=IFERROR(VALUE(MID($W13,30,3)),"")'
    sheet.Range("J12:J38").FillDown()
    sheet.Range("K12").Formula = '=IFERROR(VALUE(MID($W13,17,3)),"")'
    sheet.Range("K12:K38").FillDown()
    sheet.Range("H12:H38").NumberFormat = "yyyy-mm-dd"
    sheet.Range("I12:K38").NumberFormat = "0.00"

    sheet.Range("A12:E35").HorizontalAlignment = XL_CENTER
    sheet.Range("H12:K38").HorizontalAlignment = XL_CENTER
    sheet.Range("A4:A9,D4:D10").HorizontalAlignment = XL_LEFT
    sheet.Range("B4:B9,E4:E10").HorizontalAlignment = XL_LEFT
    sheet.Range("A4:E10").Borders(XL_BOTTOM).Color = rgb("#D9E2F3")

    widths = {"A": 23, "B": 31, "C": 19, "D": 22, "E": 32, "F": 3, "G": 3, "H": 16, "I": 14, "J": 12, "K": 12}
    for col, width in widths.items():
        sheet.Columns(col).ColumnWidth = width
    sheet.Rows("11:38").RowHeight = 20
    sheet.Rows("4:10").RowHeight = 23
    sheet.Rows(9).RowHeight = 37
    sheet.Rows(10).RowHeight = 37

    sheet.Range("L37:S38").Merge()
    sheet.Range("L37").Value = "노란색 KST 날짜·시각을 바꾸면 그래프도 해당 날짜로 갱신됩니다. ‘데이터 > 모두 새로 고침’ 후 최신 NOAA 예보를 확인하십시오. Kp는 기록 항목이며 배제 기준은 아직 확정되지 않았습니다."
    sheet.Range("L37").WrapText = True
    sheet.Range("L37").Font.Color = rgb("#44546A")
    sheet.Range("U1").Value = "NOAA 3-day raw"
    sheet.Range("W1").Value = "NOAA 27-day raw"
    sheet.Range("Y1:Z1").Value = [["선택일 KST 시각", "Kp"]]
    sheet.Range("Y2").Formula = '=IF($B$4="","",TEXT(INT($B$4)+(ROW()-2)*3/24,"hh:mm"))'
    sheet.Range("Y2:Y9").FillDown()
    sheet.Range("Z2").Formula = (
        '=IF($B$4="",NA(),IF(AND(INT($B$4)+(ROW()-2)*3/24-9/24>=$A$12,'
        'INT($B$4)+(ROW()-2)*3/24-9/24<$B$35),'
        'INDEX($E$12:$E$35,MATCH(INT($B$4)+(ROW()-2)*3/24-9/24,$A$12:$A$35,1)),NA()))'
    )
    sheet.Range("Z2:Z9").FillDown()
    sheet.Range("AB1:AC1").Value = [["선택일 전후 UTC 날짜", "일 최대 Kp"]]
    window_start = 'MAX($H$12,MIN(INT($B$4-9/24)-3,$H$38-6))'
    sheet.Range("AB2").Formula = f'=IF($B$4="","",TEXT({window_start}+ROW()-2,"m/d"))'
    sheet.Range("AB2:AB8").FillDown()
    sheet.Range("AC2").Formula = (
        f'=IF($B$4="",NA(),IFERROR(INDEX($I$12:$I$38,MATCH({window_start}+ROW()-2,$H$12:$H$38,0)),NA()))'
    )
    sheet.Range("AC2:AC8").FillDown()
    sheet.Range("Z2:Z9,AC2:AC8").NumberFormat = "0.00"
    sheet.Columns("U:AC").Hidden = True


def add_forecast_charts(sheet) -> None:
    """Add editable date-driven column charts bound to helper cells."""
    charts = (
        ("kp_3day_chart", "선택 날짜 Kp 예보 · 3시간 간격 (KST)", "Y2:Y9", "Z2:Z9", "#1D4ED8", "L3"),
        ("kp_27day_chart", "선택 날짜 전후 7일 · 일 최대 Kp (UTC)", "AB2:AB8", "AC2:AC8", "#E8531F", "L20"),
    )
    for name, title, x_range, y_range, color, anchor in charts:
        left = sheet.Range(anchor).Left
        top = sheet.Range(anchor).Top
        # 인쇄·PDF 범위(A:S) 안에 차트 전체가 들어오게 폭을 맞춘다.
        # 고정 500 pt는 오른쪽 마지막 막대를 S열 경계 밖으로 잘랐다.
        obj = sheet.ChartObjects().Add(left, top, sheet.Range("L:S").Width, 245)
        obj.Name = name
        chart = obj.Chart
        chart.ChartType = XL_COLUMN_CLUSTERED
        # 차트 원본은 사용자가 보지 않도록 숨긴 보조 열에 있다. 숨김 셀도
        # 그리지 않도록 하는 Excel 기본값을 해제해야 막대가 표시된다.
        chart.PlotVisibleOnly = False
        series = chart.SeriesCollection().NewSeries()
        series.Name = "Kp"
        series.XValues = sheet.Range(x_range)
        series.Values = sheet.Range(y_range)
        series.Format.Fill.ForeColor.RGB = rgb(color)
        series.Format.Line.ForeColor.RGB = rgb(color)
        chart.ChartGroups(1).GapWidth = 60
        chart.HasTitle = True
        chart.ChartTitle.Text = title
        chart.HasLegend = False
        chart.ChartArea.Format.Line.ForeColor.RGB = rgb("#D9E2F3")
        chart.PlotArea.Format.Fill.ForeColor.RGB = rgb("#FFFFFF")
        category = chart.Axes(XL_CATEGORY)
        category.HasTitle = True
        category.AxisTitle.Text = "시각 (KST)" if "3day" in name else "날짜 (UTC)"
        value = chart.Axes(XL_VALUE)
        value.HasTitle = True
        value.AxisTitle.Text = "Kp 지수"
        value.MinimumScale = 0
        value.MaximumScale = 9
        value.MajorUnit = 1
        value.HasMajorGridlines = True
        value.MajorGridlines.Format.Line.ForeColor.RGB = rgb("#D9E2F3")


def connect_observation_sheet(sheet) -> None:
    sheet.Range("C4").Value = "관측일시 (KST)"
    sheet.Range("E4").Value = "Kp 예보값 (자동)"
    sheet.Range("F4").Value = "예보 출처·UTC 구간"

    kp = "'⑧a Kp 예보'"
    utc = "C5-9/24"
    in_3day = f"AND({utc}>={kp}!$A$12,{utc}<{kp}!$B$35)"
    formula_e = (
        f'=IF(C5="","",IF({in_3day},'
        f'INDEX({kp}!$E$12:$E$35,MATCH({utc},{kp}!$A$12:$A$35,1)),'
        f'IFERROR(INDEX({kp}!$I$12:$I$38,MATCH(INT({utc}),{kp}!$H$12:$H$38,0)),"예보범위 밖")))'
    )
    formula_f = (
        f'=IF(C5="","",IF({in_3day},'
        f'"3일 예보 · UTC "&TEXT(INDEX({kp}!$A$12:$A$35,MATCH({utc},{kp}!$A$12:$A$35,1)),"yyyy-mm-dd hh:mm")&'
        f'"–"&TEXT(INDEX({kp}!$B$12:$B$35,MATCH({utc},{kp}!$A$12:$A$35,1)),"hh:mm"),'
        f'IF(COUNTIF({kp}!$H$12:$H$38,INT({utc}))>0,"27일 예보 · UTC 일 최대","예보범위 밖")))'
    )
    sheet.Range("E5").Formula = formula_e
    sheet.Range("E5:E54").FillDown()
    sheet.Range("F5").Formula = formula_f
    sheet.Range("F5:F54").FillDown()
    sheet.Range("E5:E54").NumberFormat = "0.00"
    sheet.Columns("C").ColumnWidth = max(sheet.Columns("C").ColumnWidth, 19)
    sheet.Columns("E").ColumnWidth = 16
    sheet.Columns("F").ColumnWidth = 31
    sheet.Range("E5:F54").Font.Color = rgb("#1D4ED8")
    sheet.Range("E5:F54").VerticalAlignment = XL_CENTER
    sheet.Range("F5:F54").WrapText = False


def connect_field_sheets(book) -> list[str]:
    """Connect each site card's observation date/time to the forecast sheet."""
    connected: list[str] = []
    kp = f"'{KP_SHEET}'"
    local_dt = "B11+D11"
    utc = f"{local_dt}-9/24"
    in_3day = f"AND({utc}>={kp}!$A$12,{utc}<{kp}!$B$35)"
    kp_formula = (
        f'=IF(OR(B11="",D11=""),"",IF({in_3day},'
        f'INDEX({kp}!$E$12:$E$35,MATCH({utc},{kp}!$A$12:$A$35,1)),'
        f'IFERROR(INDEX({kp}!$I$12:$I$38,MATCH(INT({utc}),{kp}!$H$12:$H$38,0)),"예보범위 밖")))'
    )
    source_formula = (
        f'=IF(OR(B11="",D11=""),"",IF({in_3day},'
        f'"3일 · UTC "&TEXT(INDEX({kp}!$A$12:$A$35,MATCH({utc},{kp}!$A$12:$A$35,1)),"m/d hh:mm")&'
        f'"–"&TEXT(INDEX({kp}!$B$12:$B$35,MATCH({utc},{kp}!$A$12:$A$35,1)),"hh:mm"),'
        f'IF(COUNTIF({kp}!$H$12:$H$38,INT({utc}))>0,"27일 · UTC 일 최대","예보범위 밖")))'
    )
    for sheet in book.Worksheets:
        if sheet.Name == KP_SHEET:
            continue
        if sheet.Range("A11").Text != "관측일자" or sheet.Range("A12").Text != "Kp 지수":
            continue
        sheet.Range("B12").Formula = kp_formula
        sheet.Range("B12").NumberFormat = "0.00"
        sheet.Range("B12").Font.Color = rgb("#1D4ED8")
        sheet.Range("D12").Formula = source_formula
        sheet.Range("D12").Font.Color = rgb("#1D4ED8")
        sheet.Range("D12").Font.Size = 8
        sheet.Range("D12").WrapText = True
        connected.append(sheet.Name)
    return connected


def export_range_png(sheet, address: str, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rng = sheet.Range(address)
    pdf_path = output_path.with_suffix(".pdf")
    setup = sheet.PageSetup
    old_print_area = setup.PrintArea
    old_orientation = setup.Orientation
    old_zoom = setup.Zoom
    old_fit_wide = setup.FitToPagesWide
    old_fit_tall = setup.FitToPagesTall
    try:
        setup.PrintArea = rng.Address
        setup.Orientation = 2
        setup.Zoom = False
        setup.FitToPagesWide = 1
        setup.FitToPagesTall = 1
        sheet.ExportAsFixedFormat(0, str(pdf_path), 0, True, False)
        with fitz.open(pdf_path) as document:
            pixmap = document[0].get_pixmap(matrix=fitz.Matrix(1.8, 1.8), alpha=False)
            pixmap.save(output_path)
    finally:
        setup.PrintArea = old_print_area
        setup.Orientation = old_orientation
        setup.Zoom = old_zoom
        setup.FitToPagesWide = old_fit_wide
        setup.FitToPagesTall = old_fit_tall
        if pdf_path.exists():
            pdf_path.unlink()


def build(input_path: Path, output_path: Path, preview_dir: Path | None) -> None:
    if output_path.exists():
        raise FileExistsError(f"Output already exists: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, output_path)

    pythoncom.CoInitialize()
    excel = win32com.client.DispatchEx("Excel.Application")
    excel.Visible = False
    excel.DisplayAlerts = False
    book = None
    try:
        book = excel.Workbooks.Open(str(output_path), UpdateLinks=0, ReadOnly=False)
        active_name = book.ActiveSheet.Name
        sheet_names = [sheet.Name for sheet in book.Worksheets]
        standard_mode = OBS_SHEET in sheet_names
        anchor = book.Worksheets(OBS_SHEET) if standard_mode else book.Worksheets(1)
        kp_sheet = book.Worksheets.Add(After=anchor)
        kp_sheet.Name = KP_SHEET
        kp_sheet.Tab.Color = rgb("#1D4ED8")

        configure_kp_sheet(kp_sheet)
        add_web_query(kp_sheet, "U2", "NOAA_3day_forecast", NOAA_3DAY)
        add_web_query(kp_sheet, "W2", "NOAA_27day_outlook", NOAA_27DAY)
        if standard_mode:
            observation = book.Worksheets(OBS_SHEET)
            connect_observation_sheet(observation)
        else:
            connect_field_sheets(book)

        # 새 수식 범위를 먼저 계산한 뒤 차트를 생성해야 Excel이 첫 저장 전에
        # 계열 캐시와 축을 채운다. 차트를 먼저 만들면 PDF/PNG 내보내기에서
        # 빈 차트 프레임만 나오는 경우가 있다.
        excel.CalculateFullRebuild()
        add_forecast_charts(kp_sheet)
        for index in range(1, kp_sheet.ChartObjects().Count + 1):
            kp_sheet.ChartObjects(index).Chart.Refresh()
        if preview_dir is not None:
            export_range_png(kp_sheet, "A1:S38", preview_dir / "kp_forecast.png")
            if standard_mode:
                export_range_png(observation, "A1:F16", preview_dir / "observation_kp.png")
            else:
                export_range_png(book.Worksheets(3), "A1:I18", preview_dir / "field_card_kp.png")
        book.Worksheets(active_name).Activate()
        book.Save()
    finally:
        if book is not None:
            book.Close(SaveChanges=False)
        excel.Quit()
        pythoncom.CoUninitialize()


def verify(output_path: Path) -> dict[str, object]:
    pythoncom.CoInitialize()
    excel = win32com.client.DispatchEx("Excel.Application")
    excel.Visible = False
    excel.DisplayAlerts = False
    book = None
    try:
        book = excel.Workbooks.Open(str(output_path), UpdateLinks=0, ReadOnly=True)
        kp_sheet = book.Worksheets(KP_SHEET)
        sheet_names = [sheet.Name for sheet in book.Worksheets]
        standard_mode = OBS_SHEET in sheet_names
        observation = book.Worksheets(OBS_SHEET) if standard_mode else None
        field_sheets = [
            sheet
            for sheet in book.Worksheets
            if sheet.Range("A11").Text == "관측일자" and sheet.Range("A12").Text == "Kp 지수"
        ]
        test_time = datetime(2026, 9, 14, 13, 0)
        excel_serial = (test_time - datetime(1899, 12, 30)).total_seconds() / 86400
        kp_sheet.Range("B4").Value = excel_serial
        if standard_mode:
            observation.Range("C5").Value = excel_serial
        else:
            field_sheets[0].Range("B11").Value = int(excel_serial)
            field_sheets[0].Range("D11").Value = excel_serial - int(excel_serial)
        excel.CalculateFull()
        result = {
            "mode": "standard" if standard_mode else "field",
            "sheets": book.Worksheets.Count,
            "connections": book.Connections.Count,
            "queries": kp_sheet.QueryTables.Count,
            "charts": kp_sheet.ChartObjects().Count,
            "chart_types": [
                kp_sheet.ChartObjects(index).Chart.ChartType
                for index in range(1, kp_sheet.ChartObjects().Count + 1)
            ],
            "chart_titles": [
                kp_sheet.ChartObjects(index).Chart.ChartTitle.Text
                for index in range(1, kp_sheet.ChartObjects().Count + 1)
            ],
            "connected_field_sheets": len(field_sheets),
            "refresh_on_open": [bool(kp_sheet.QueryTables(index).RefreshOnFileOpen) for index in range(1, kp_sheet.QueryTables.Count + 1)],
            "saved_query_data": [bool(kp_sheet.QueryTables(index).SaveData) for index in range(1, kp_sheet.QueryTables.Count + 1)],
            "issued_3day": kp_sheet.Range("E4").Text,
            "issued_27day": kp_sheet.Range("E5").Text,
            "test_kst": kp_sheet.Range("B4").Text,
            "test_utc": kp_sheet.Range("B5").Text,
            "quick_3day_kp": kp_sheet.Range("B6").Value,
            "quick_27day_kp": kp_sheet.Range("B7").Value,
            "selected_3day_labels": [row[0] for row in kp_sheet.Range("Y2:Y9").Value],
            "selected_3day_values": [row[0] for row in kp_sheet.Range("Z2:Z9").Value],
            "selected_27day_labels": [row[0] for row in kp_sheet.Range("AB2:AB8").Value],
            "selected_27day_values": [row[0] for row in kp_sheet.Range("AC2:AC8").Value],
        }
        if standard_mode:
            result["observation_kp"] = observation.Range("E5").Value
            result["observation_source"] = observation.Range("F5").Text
        else:
            result["field_kp"] = field_sheets[0].Range("B12").Value
            result["field_source"] = field_sheets[0].Range("D12").Text
            result["field_kp_formula"] = field_sheets[0].Range("B12").Formula
        first_day_label = kp_sheet.Range("AB2").Text
        kp_sheet.Range("B4").Value = excel_serial + 5
        excel.CalculateFull()
        result["date_change_updates_chart"] = kp_sheet.Range("AB2").Text != first_day_label
        fallback_time = kp_sheet.Range("D35").Value2
        kp_sheet.Range("B4").Value = fallback_time
        if standard_mode:
            observation.Range("C5").Value = fallback_time
        else:
            field_sheets[0].Range("B11").Value = int(fallback_time)
            field_sheets[0].Range("D11").Value = fallback_time - int(fallback_time)
        excel.CalculateFull()
        result["fallback_source"] = kp_sheet.Range("B8").Text
        if standard_mode:
            result["fallback_observation_kp"] = observation.Range("E5").Value
            result["fallback_observation_source"] = observation.Range("F5").Text
        else:
            result["fallback_field_kp"] = field_sheets[0].Range("B12").Value
            result["fallback_field_source"] = field_sheets[0].Range("D12").Text

        kp_sheet.Range("B4").ClearContents()
        if standard_mode:
            observation.Range("C5").ClearContents()
        else:
            field_sheets[0].Range("B11").ClearContents()
            field_sheets[0].Range("D11").ClearContents()
        excel.CalculateFull()
        result["blank_quick_result"] = kp_sheet.Range("B6").Text
        result["blank_observation_result"] = (
            observation.Range("E5").Text if standard_mode else field_sheets[0].Range("B12").Text
        )
        error_tokens = ("#REF!", "#DIV/0!", "#VALUE!", "#NAME?", "#N/A", "#NUM!", "#NULL!", "#SPILL!", "#CALC!")
        errors: list[str] = []
        checks = [(kp_sheet, "A1:S38")]
        if standard_mode:
            checks.append((observation, "C4:F54"))
        else:
            checks.extend((sheet, "B11:D12") for sheet in field_sheets)
        for sheet, address in checks:
            values = sheet.Range(address).Value
            if not isinstance(values, tuple):
                values = ((values,),)
            for r_idx, row in enumerate(values, start=1):
                for c_idx, value in enumerate(row, start=1):
                    if isinstance(value, str) and any(token in value for token in error_tokens):
                        errors.append(f"{sheet.Name}:{address}:{r_idx},{c_idx}:{value}")
        result["formula_errors"] = errors
        return result
    finally:
        if book is not None:
            book.Close(SaveChanges=False)
        excel.Quit()
        pythoncom.CoUninitialize()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--preview-dir", type=Path)
    args = parser.parse_args()
    build(args.input.resolve(), args.output.resolve(), args.preview_dir.resolve() if args.preview_dir else None)
    for key, value in verify(args.output.resolve()).items():
        print(f"{key}={ascii(value)}")


if __name__ == "__main__":
    main()
