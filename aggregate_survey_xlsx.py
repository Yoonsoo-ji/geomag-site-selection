# -*- coding: utf-8 -*-
"""
현장 조사서(엑셀 카드) → 취합 + 선점 검토 + 웹 표출용 정리
==========================================================

각 본부가 회신한 조사서 워크북(총괄 목록 + 카드 시트들)을 모두 읽어

  1) [현장조사 취합]  한 카드 = 한 행. 도상선점 정보 + 현장 조사 결과 전부.
  2) [선점 검토]      취합 결과를 근거로 선점 가능성 등급·결론·핵심 교란요인 정리
                      (웹 페이지 표출을 염두에 둔 열 구성).
  3) [검토 요약]      본부별·판정별 분포, 주요 교란요인 통계.

카드 시트는 make_survey_card_v2.py 로 생성된 고정 레이아웃이므로 셀 위치가 일정하다.
자동판정(F열·집계·종합판정)은 캐시된 수식값 대신 존재여부에서 규칙으로 재계산해
파일이 Excel 에서 평가되지 않았어도 안전하게 판정한다.

사용:
    python aggregate_survey_xlsx.py                        # 기본 폴더 전체
    python aggregate_survey_xlsx.py --dir <폴더> --out <파일.xlsx>
"""
import argparse
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).parent
# 회신 조사서가 놓인 베이스 폴더(하위 폴더명은 자주 바뀌므로 재귀 탐색한다)
DEF_BASE = Path(r"D:\LX_yoons\2026_research\2026_지자기 연구\20260811_사전 현장 조사서 취합")
DEF_DIR = DEF_BASE   # 하위호환 별칭
OUT_DIR = ROOT / "docs" / "output"


def survey_files(base=None):
    """베이스 폴더 아래에서 'N.조사서_지역.xlsx' 회신본을 재귀 탐색(폴더명 변경에 견고)."""
    base = Path(base) if base else DEF_BASE
    fs = [p for p in base.rglob("*.조사서_*.xlsx") if not p.name.startswith("~$")]
    return sorted(fs, key=lambda p: int(re.match(r"(\d+)", p.name).group(1))
                  if re.match(r"(\d+)", p.name) else 999)

# (카드 라벨, 취합 열 이름)
DIST_ITEMS = [("차량 영향", "차량영향"), ("전력시설", "전력시설"), ("통신시설", "통신시설"),
              ("금속류", "금속류"), ("매설물", "매설물"), ("건축물", "건축물"),
              ("자연환경", "자연환경")]

# ── 검토 단계 수동 재판정 (도상선점↔현장 기준점 좌표 1km 이상 불일치 5건) ──
# {관리번호: (종합판정, 재판정 사유)}. 원본 카드는 불변, 취합·검토·총괄·웹에만 반영.
_SHADOW = "도상 후보지 통신 음영지역으로 확인 어려운 지역 → 부적합"
JUDGMENT_OVERRIDE = {
    "DS-002": ("부적합", _SHADOW),
    "DS-041": ("부적합", _SHADOW),
    "DS-061": ("부적합", _SHADOW),
    "DS-011": ("부적합", "현장 이설·좌표 불일치 → 부적합"),
    "DS-087": ("부적합", "상공장애로 취득 불가 → 부적합"),
}


# ── 카드 시트 파싱 ────────────────────────────────────────────────────────────
def _norm(v):
    return re.split(r"[\n(]", str(v))[0].strip() if v is not None else ""


def label_rows(ws, col):
    m = {}
    for r in range(1, ws.max_row + 1):
        k = _norm(ws.cell(r, col).value)
        if k:
            m.setdefault(k, r)
    return m


def _s(v):
    return "" if v is None else str(v).strip()


def dms2dd(v):
    """도분초 문자열 → 십진도. 형식 무관(대시·°'\"·공백). 앞 숫자 3개를 도·분·초로 본다."""
    if v is None:
        return None
    s = str(v).strip()
    if not s or s == "-":
        return None
    nums = re.findall(r"\d*\.?\d+", s)   # '.6614'(앞자리 0 생략)도 초로 인식
    if len(nums) < 3:
        return None
    d, mi, se = float(nums[0]), float(nums[1]), float(nums[2])
    return d + mi / 60 + se / 3600


def _fmt_date(v):
    if v is None or v == "":
        return ""
    s = str(v).strip()
    if re.fullmatch(r"20\d{6}", s):
        return f"{s[:4]}-{s[4:6]}-{s[6:8]}"
    if hasattr(v, "strftime"):
        return v.strftime("%Y-%m-%d")
    return s


def parse_card(ws):
    A = label_rows(ws, 1)
    D = label_rows(ws, 4)
    d = {}

    def bA(name):  # A라벨 → B값
        r = A.get(name)
        return _s(ws.cell(r, 2).value) if r else ""

    def eD(name):  # D라벨 → E값
        r = D.get(name)
        return _s(ws.cell(r, 5).value) if r else ""

    d["관리번호"] = bA("관리번호")
    d["후보지명"] = eD("후보지명")
    d["도엽번호"] = bA("도엽번호")
    d["도엽명"] = eD("도엽명")
    d["위도"] = bA("위도")          # '위도(십진)' → norm '위도'
    d["경도"] = eD("경도")          # '경도(십진)' → norm '경도'
    d["표고"] = bA("표고").replace(" m", "").replace("m", "").strip()
    d["유형"] = eD("후보지 유형")
    d["지번주소"] = bA("지번 주소")
    d["도로명주소"] = bA("도로명 주소")
    d["소유권"] = bA("소유권 조사")
    d["관할본부"] = eD("관할 본부")
    d["조사대상"] = bA("조사대상")
    d["조사일"] = _fmt_date(ws.cell(D["조사일"], 5).value) if "조사일" in D else ""
    d["조사자"] = bA("조사자")
    d["기상"] = eD("기상상태")
    d["연계기존점"] = bA("연계 기존점")
    d["접근경로"] = bA("접근 경로").replace("\n", " ")
    d["차량진입"] = bA("차량 진입")
    d["유의사항"] = bA("유의사항")
    d["도상간섭"] = bA("도상 간섭요인").replace("\n", " ")

    # ── 자기교란 8항목: D=존재여부, E=실측이격 ──
    dist = {}
    for lab, key in DIST_ITEMS:
        r = A.get(lab)
        if r:
            dist[key] = (_s(ws.cell(r, 4).value), _s(ws.cell(r, 5).value))
    d["자기교란"] = dist
    rb = A.get("방위표지")
    d["방위표지"] = _s(ws.cell(rb, 4).value) if rb else ""   # 가능/불가

    # ── 종합판정 규칙 재계산 (캐시 수식값에 의존하지 않음) ──
    hits, missing = [], 0
    for lab, key in DIST_ITEMS:
        ex, gap = dist.get(key, ("", ""))
        if ex == "":
            missing += 1
        elif ex == "있음":
            hits.append((key, gap))
    if d["방위표지"] == "":
        missing += 1
    bang_bad = d["방위표지"] == "불가"
    nbad = len(hits) + (1 if bang_bad else 0)
    if missing > 0:
        verdict = "조사 미완료"
    elif nbad == 0:
        verdict = "적합"
    elif nbad <= 2:
        verdict = "조건부 적합"
    else:
        verdict = "부적합"
    d["종합판정"] = verdict
    d["부적합수"] = nbad
    d["교란hits"] = hits
    d["방위불가"] = bang_bad
    d["판정의견"] = bA("판정 의견")

    # ── 검토 단계 수동 재판정(원본 카드는 불변, 취합·검토·총괄·웹에만 반영) ──
    if d["관리번호"] in JUDGMENT_OVERRIDE:
        verdict_ov, reason = JUDGMENT_OVERRIDE[d["관리번호"]]
        d["종합판정"] = verdict_ov
        d["재판정"] = reason
        note = d["판정의견"].strip()
        d["판정의견"] = (note + "  ※재판정: " + reason) if note else reason
    else:
        d["재판정"] = ""

    # ── 방위표지 좌표 상세 (③) ──
    # A맵의 '위도'는 A5(위도(십진))와 충돌하므로 '경도'(A34) 앵커에서 offset 으로 읽는다.
    r_lon = A.get("경도")   # 34
    detail = {}
    d["기준점ll"] = d["표지1ll"] = d["표지2ll"] = None
    if r_lon:
        def coord(col):   # col 2=기준점(B), 4=표지1(D), 6=표지2(F)
            return {"경도": _s(ws.cell(r_lon, col).value),
                    "위도": _s(ws.cell(r_lon + 1, col).value),
                    "방위각": _s(ws.cell(r_lon + 3, col).value),
                    "거리": _s(ws.cell(r_lon + 4, col).value)}
        detail = {"기준점": coord(2), "표지1": coord(4), "표지2": coord(6)}

        def ll(col):   # (경도, 위도) 십진 or None. 경도/위도 칸 스왑 입력 자동 교정.
            a = dms2dd(ws.cell(r_lon, col).value)      # '경도' 라벨 칸
            b = dms2dd(ws.cell(r_lon + 1, col).value)  # '위도' 라벨 칸
            if a is None or b is None:
                return None
            # 한반도 범위(경도 124~132 · 위도 33~40)로 판별 — 본부별 칸 뒤바뀜 대응
            if 124 <= a <= 132 and 33 <= b <= 40:
                return (a, b)
            if 124 <= b <= 132 and 33 <= a <= 40:
                return (b, a)
            return None   # 범위 밖이면 무효
        d["기준점ll"] = ll(2)
        d["표지1ll"] = ll(4)
        d["표지2ll"] = ll(6)
    d["방위표지상세"] = detail
    az1 = detail.get("표지1", {}).get("방위각", "")
    d["방위표지좌표"] = "입력" if az1 and az1 != "-" else "미입력"
    return d


# ── 카드 임베드 사진 추출 (중심·동·서·남·북·교란원) ──────────────────────────
def _photo_slot(row, col):
    """anchor (row,col) 0-index → 사진 슬롯. 지도(col>=7)·해당없음은 None."""
    if col >= 7:
        return None
    if col in (0, 1):   # A·B 열
        if 28 <= row <= 31:
            return "교란원"
        if 43 <= row <= 45:
            return "중심"
        if 46 <= row <= 48:
            return "동"
        if 49 <= row <= 52:
            return "남"
    if col in (3, 4):   # D·E 열
        if 46 <= row <= 48:
            return "서"
        if 49 <= row <= 52:
            return "북"
    return None


def extract_photos(ws, code, outdir, maxpx=300, quality=75):
    """카드 시트의 방위/중심 사진을 다운스케일 JPEG 로 저장. {슬롯: 파일명} 반환."""
    import io as _io
    from PIL import Image
    outdir.mkdir(parents=True, exist_ok=True)
    got = {}
    for im in getattr(ws, "_images", []):
        try:
            a = im.anchor._from
        except AttributeError:
            continue
        slot = _photo_slot(a.row, a.col)
        if not slot or slot in got:
            continue
        try:
            img = Image.open(_io.BytesIO(im._data())).convert("RGB")
        except Exception:   # noqa: BLE001
            continue
        img.thumbnail((maxpx, maxpx))
        fname = f"{code}_{slot}.jpg"
        img.save(outdir / fname, "JPEG", quality=quality, optimize=True)
        got[slot] = fname
    return got


def parse_workbook(path, photo_dir=None):
    wb = load_workbook(path, data_only=True)
    out = []
    for name in wb.sheetnames:
        if name == "총괄 목록":
            continue
        ws = wb[name]
        d = parse_card(ws)
        if d.get("관리번호"):
            d["_src"] = path.name
            if photo_dir is not None:
                d["사진"] = extract_photos(ws, d["관리번호"], Path(photo_dir))
            out.append(d)
    wb.close()
    return out


# ── 검토(선점 가능성) 파생 ────────────────────────────────────────────────────
def key_disturb(d):
    parts = []
    for key, gap in d["교란hits"]:
        parts.append(f"{key}({gap}m)" if gap and gap != "-" else key)
    if d["방위불가"]:
        parts.append("방위표지 불가")
    return ", ".join(parts)


def review(d):
    """(등급, 선점 검토 결론, 검토 의견) — 웹 표출용."""
    v = d["종합판정"]
    kd = key_disturb(d)
    op = d["판정의견"].replace("\n", " ").strip()
    if v == "조사 미완료":
        return "미완료", "재조사 필요 — 조사 항목 누락", op or "조사 항목 미입력"
    if d["방위불가"]:
        grade = "C"
        concl = "부적합 — 방위표지 설치 불가로 대체 후보지 검토"
    elif v == "적합":
        grade, concl = "A", "선점 가능 — 자기구배 조사 진행"
    elif v == "조건부 적합":
        grade, concl = "B", f"조건부 가능 — 현장 재확인 후 결정 ({kd})"
    else:
        grade, concl = "C", f"부적합 — 대체 후보지 검토 ({kd})"
    note = op if op else ("주변 자기교란 요소 없음" if v == "적합" else kd)
    return grade, concl, note


# ── 스타일 ────────────────────────────────────────────────────────────────────
F_TITLE = Font(name="맑은 고딕", size=13, bold=True, color="FFFFFF")
F_HDR = Font(name="맑은 고딕", size=9, bold=True, color="FFFFFF")
F_VAL = Font(name="맑은 고딕", size=9)
F_BOLD = Font(name="맑은 고딕", size=9, bold=True)
FILL_TITLE = PatternFill("solid", fgColor="1F3864")
FILL_HDR = PatternFill("solid", fgColor="4472C4")
FILL_BAD = PatternFill("solid", fgColor="FDE0DC")
FILL_COND = PatternFill("solid", fgColor="FFF2CC")
FILL_OK = PatternFill("solid", fgColor="E2EFDA")
FILL_GRAY = PatternFill("solid", fgColor="EDEDED")
GRADE_FILL = {"A": FILL_OK, "B": FILL_COND, "C": FILL_BAD, "미완료": FILL_GRAY}
AL_L = Alignment(horizontal="left", vertical="center", wrap_text=True)
AL_C = Alignment(horizontal="center", vertical="center", wrap_text=True)
_thin = Side(style="thin", color="BBBBBB")
BORDER = Border(left=_thin, right=_thin, top=_thin, bottom=_thin)


def _title(ws, ncol, text):
    ws.merge_cells(f"A1:{get_column_letter(ncol)}1")
    c = ws["A1"]
    c.value = text
    c.font = F_TITLE
    c.fill = FILL_TITLE
    c.alignment = AL_L
    ws.row_dimensions[1].height = 26


def _header(ws, cols, row=2):
    for j, h in enumerate(cols, 1):
        c = ws.cell(row, j, h)
        c.font = F_HDR
        c.fill = FILL_HDR
        c.alignment = AL_C
        c.border = BORDER


def dcell(item):
    if not item:
        return ""
    ex, gap = item
    if ex == "있음":
        return f"있음 ({gap}m)" if gap and gap != "-" else "있음"
    return ex


# ── 시트 1: 현장조사 취합 ────────────────────────────────────────────────────
FULL_COLS = ["관할 본부", "관리번호", "후보지명", "도엽번호", "도엽명", "지번 주소",
             "위도", "경도", "표고(m)", "유형", "연계 기존점", "소유권",
             "조사일", "조사자", "기상",
             "차량영향", "전력시설", "통신시설", "금속류", "매설물", "건축물", "자연환경",
             "방위표지", "방위표지좌표", "부적합 건수", "종합 판정", "판정 의견",
             "접근 경로", "유의사항", "도상 간섭요인", "회신 파일"]


def sheet_full(wb, recs):
    ws = wb.create_sheet("현장조사 취합")
    n = len(FULL_COLS)
    _title(ws, n, f"지자기 도상선점 — 현장 조사서 취합 ({len(recs)}건)   ·   {datetime.now():%Y-%m-%d}")
    _header(ws, FULL_COLS)
    ds = FULL_COLS.index("차량영향") + 1
    de = FULL_COLS.index("자연환경") + 1
    bcol = FULL_COLS.index("방위표지") + 1
    jcol = FULL_COLS.index("종합 판정") + 1
    for ri, d in enumerate(recs, 3):
        dist = d["자기교란"]
        row = [d["관할본부"], d["관리번호"], d["후보지명"], d["도엽번호"], d["도엽명"],
               d["지번주소"], d["위도"], d["경도"], d["표고"], d["유형"], d["연계기존점"],
               d["소유권"], d["조사일"], d["조사자"], d["기상"],
               dcell(dist.get("차량영향")), dcell(dist.get("전력시설")), dcell(dist.get("통신시설")),
               dcell(dist.get("금속류")), dcell(dist.get("매설물")), dcell(dist.get("건축물")),
               dcell(dist.get("자연환경")), d["방위표지"], d["방위표지좌표"], d["부적합수"],
               d["종합판정"], d["판정의견"].replace("\n", " "),
               d["접근경로"], d["유의사항"], d["도상간섭"], d["_src"]]
        for j, v in enumerate(row, 1):
            c = ws.cell(ri, j, v)
            c.font = F_VAL
            c.border = BORDER
            c.alignment = AL_L if j in (6, 27, 28, 29, 30) else AL_C
            if ds <= j <= de and isinstance(v, str) and v.startswith("있음"):
                c.fill = FILL_BAD
            if j == bcol and v == "불가":
                c.fill = FILL_BAD
            if j == jcol:
                c.fill = {"적합": FILL_OK, "조건부 적합": FILL_COND,
                          "부적합": FILL_BAD}.get(v, FILL_GRAY)
    widths = [14, 8, 13, 11, 8, 24, 9, 9, 7, 8, 9, 8, 11, 12, 6,
              9, 8, 8, 8, 8, 8, 8, 8, 9, 8, 10, 30, 24, 18, 30, 16]
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "D3"
    ws.auto_filter.ref = f"A2:{get_column_letter(n)}{len(recs) + 2}"


# ── 시트 2: 선점 검토 (웹 표출용) ────────────────────────────────────────────
WEB_COLS = ["등급", "관할 본부", "관리번호", "후보지명", "위도", "경도", "표고(m)",
            "종합 판정", "핵심 교란요인", "방위표지", "선점 검토 결론", "검토 의견", "조사일"]


def sheet_web(wb, recs):
    ws = wb.create_sheet("선점 검토")
    n = len(WEB_COLS)
    _title(ws, n, "지자기 도상선점 — 선점 위치 검토 (웹 표출용)")
    # 등급 범례
    ws.merge_cells(f"A2:{get_column_letter(n)}2")
    lg = ws["A2"]
    lg.value = ("등급  A = 선점 가능(자기구배 조사)   ·   B = 조건부(현장 재확인)   ·   "
                "C = 부적합(대체 후보지 검토)   ·   방위표지 '불가'는 C 로 분류")
    lg.font = Font(name="맑은 고딕", size=8.5, color="555555")
    lg.alignment = AL_L
    _header(ws, WEB_COLS, row=3)
    order = {"A": 0, "B": 1, "C": 2, "미완료": 3}
    rr = sorted(recs, key=lambda d: (order.get(review(d)[0], 9), d["관할본부"], d["관리번호"]))
    for ri, d in enumerate(rr, 4):
        grade, concl, note = review(d)
        row = [grade, d["관할본부"], d["관리번호"], d["후보지명"], d["위도"], d["경도"],
               d["표고"], d["종합판정"], key_disturb(d) or "-", d["방위표지"],
               concl, note, d["조사일"]]
        for j, v in enumerate(row, 1):
            c = ws.cell(ri, j, v)
            c.font = F_BOLD if j == 1 else F_VAL
            c.border = BORDER
            c.alignment = AL_L if j in (9, 11, 12) else AL_C
            if j == 1:
                c.fill = GRADE_FILL.get(grade, FILL_GRAY)
            if j == 8:
                c.fill = {"적합": FILL_OK, "조건부 적합": FILL_COND,
                          "부적합": FILL_BAD}.get(v, FILL_GRAY)
            if j == 10 and v == "불가":
                c.fill = FILL_BAD
    widths = [6, 14, 8, 14, 9, 9, 7, 10, 26, 8, 34, 40, 11]
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "A4"
    ws.auto_filter.ref = f"A3:{get_column_letter(n)}{len(rr) + 3}"


# ── 시트 3: 검토 요약 ─────────────────────────────────────────────────────────
def sheet_summary(wb, recs):
    ws = wb.create_sheet("검토 요약", 0)
    _title(ws, 6, "지자기 도상선점 현장조사 — 검토 요약")
    grades = Counter(review(d)[0] for d in recs)
    verd = Counter(d["종합판정"] for d in recs)
    bang = sum(1 for d in recs if d["방위불가"])
    dist_cnt = Counter()
    for d in recs:
        for key, _ in d["교란hits"]:
            dist_cnt[key] += 1

    r = 3
    def line(label, val, fill=None, bold=False):
        nonlocal r
        a = ws.cell(r, 1, label)
        a.font = F_BOLD
        a.alignment = AL_L
        ws.merge_cells(f"A{r}:C{r}")
        b = ws.cell(r, 4, val)
        b.font = F_BOLD if bold else F_VAL
        b.alignment = AL_L
        ws.merge_cells(f"D{r}:F{r}")
        if fill:
            for cc in (a, b):
                cc.fill = fill
        r += 1

    ws.cell(r, 1, "■ 전체 현황").font = F_BOLD
    r += 1
    line("총 조사 후보지", f"{len(recs)} 건", bold=True)
    line("선점 가능 (등급 A)", f"{grades['A']} 건", FILL_OK)
    line("조건부 (등급 B)", f"{grades['B']} 건", FILL_COND)
    line("부적합·재검토 (등급 C)", f"{grades['C']} 건", FILL_BAD)
    if grades.get("미완료"):
        line("조사 미완료", f"{grades['미완료']} 건", FILL_GRAY)
    r += 1
    ws.cell(r, 1, "■ 종합 판정 분포 (조사자 카드 기준)").font = F_BOLD
    r += 1
    for k in ("적합", "조건부 적합", "부적합", "조사 미완료"):
        if verd.get(k):
            line(k, f"{verd[k]} 건")
    r += 1
    ws.cell(r, 1, "■ 특이사항").font = F_BOLD
    r += 1
    line("방위표지 설치 불가", f"{bang} 건", FILL_BAD if bang else None)
    r += 1
    ws.cell(r, 1, "■ 주요 교란요인 (있음 건수)").font = F_BOLD
    r += 1
    for key, c in dist_cnt.most_common():
        line(key, f"{c} 건")
    r += 1
    ws.cell(r, 1, "■ 본부별 등급 분포").font = F_BOLD
    r += 1
    hqs = sorted({d["관할본부"] for d in recs})
    hdr = ["관할 본부", "A", "B", "C", "미완료", "계"]
    for j, h in enumerate(hdr, 1):
        c = ws.cell(r, j, h)
        c.font = F_HDR
        c.fill = FILL_HDR
        c.alignment = AL_C
        c.border = BORDER
    r += 1
    for hq in hqs:
        sub = [d for d in recs if d["관할본부"] == hq]
        g = Counter(review(d)[0] for d in sub)
        vals = [hq, g["A"], g["B"], g["C"], g.get("미완료", 0), len(sub)]
        for j, v in enumerate(vals, 1):
            c = ws.cell(r, j, v)
            c.font = F_VAL
            c.alignment = AL_L if j == 1 else AL_C
            c.border = BORDER
        r += 1
    for j, w in enumerate([16, 6, 6, 6, 8, 6], 1):
        ws.column_dimensions[get_column_letter(j)].width = w


# ── 일괄취합(플랫) 워크북 — 분석·등급 없이 모든 원자료를 한 시트에 ──────────
def _action(v):
    return {"적합": "자기구배 조사 진행", "조건부 적합": "추가 자기조사 후 재판정",
            "부적합": "대체 후보지 검토"}.get(v, "② 조사 완료 후")


FLAT_COLS = ["관할 본부", "관리번호", "후보지명", "도엽번호", "도엽명", "위도", "경도",
             "표고(m)", "유형", "지번 주소", "도로명 주소", "연계 기존점",
             "소유권", "조사대상", "조사일", "조사자", "기상",
             "차량영향 존재", "차량영향 이격", "전력시설 존재", "전력시설 이격",
             "통신시설 존재", "통신시설 이격", "금속류 존재", "금속류 이격",
             "매설물 존재", "매설물 이격", "건축물 존재", "건축물 이격",
             "자연환경 존재", "자연환경 이격", "방위표지 가능",
             "부적합 건수", "종합 판정", "판정 의견", "향후 조치",
             "방위1 경도", "방위1 위도", "방위1 방위각", "방위1 거리(m)",
             "방위2 경도", "방위2 위도", "방위2 방위각", "방위2 거리(m)",
             "접근 경로", "차량 진입", "유의사항", "도상 간섭요인", "회신 파일"]


def build_flat_workbook(recs, out_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "일괄취합"
    n = len(FLAT_COLS)
    _title(ws, n, f"지자기 도상선점 — 현장 조사서 일괄취합 ({len(recs)}건)   ·   {datetime.now():%Y-%m-%d}")
    _header(ws, FLAT_COLS)
    for ri, d in enumerate(recs, 3):
        dist = d["자기교란"]

        def g(key):
            ex, gap = dist.get(key, ("", ""))
            return ex, gap
        b1 = d.get("방위표지상세", {}).get("표지1", {})
        b2 = d.get("방위표지상세", {}).get("표지2", {})
        row = [d["관할본부"], d["관리번호"], d["후보지명"], d["도엽번호"], d["도엽명"],
               d["위도"], d["경도"], d["표고"], d["유형"], d["지번주소"], d["도로명주소"],
               d["연계기존점"], d["소유권"], d["조사대상"], d["조사일"], d["조사자"], d["기상"],
               *g("차량영향"), *g("전력시설"), *g("통신시설"), *g("금속류"),
               *g("매설물"), *g("건축물"), *g("자연환경"), d["방위표지"],
               d["부적합수"], d["종합판정"], d["판정의견"].replace("\n", " "), _action(d["종합판정"]),
               b1.get("경도", ""), b1.get("위도", ""), b1.get("방위각", ""), b1.get("거리", ""),
               b2.get("경도", ""), b2.get("위도", ""), b2.get("방위각", ""), b2.get("거리", ""),
               d["접근경로"], d["차량진입"], d["유의사항"], d["도상간섭"], d["_src"]]
        for j, v in enumerate(row, 1):
            c = ws.cell(ri, j, v)
            c.font = F_VAL
            c.border = BORDER
            c.alignment = AL_L if FLAT_COLS[j - 1] in (
                "지번 주소", "도로명 주소", "판정 의견", "접근 경로", "유의사항", "도상 간섭요인") else AL_C
    ws.freeze_panes = "C3"
    ws.auto_filter.ref = f"A2:{get_column_letter(n)}{len(recs) + 2}"
    for j, col in enumerate(FLAT_COLS, 1):
        w = 30 if col in ("판정 의견", "도상 간섭요인") else \
            24 if col in ("지번 주소", "접근 경로") else \
            18 if col == "유의사항" else \
            13 if "방위" in col or col in ("조사일", "조사자", "관할 본부", "후보지명") else 9
        ws.column_dimensions[get_column_letter(j)].width = w
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(DEF_BASE))
    ap.add_argument("--out")
    a = ap.parse_args()

    files = survey_files(a.dir)
    if not files:
        print(f"조사서 xlsx 를 찾지 못했습니다: {a.dir}")
        return

    recs = []
    for f in files:
        cards = parse_workbook(f)
        print(f"  {f.name:28} 카드 {len(cards):3}건")
        recs += cards
    recs.sort(key=lambda d: (d["관할본부"] or "", d["관리번호"] or ""))
    print(f"총 {len(recs)}건 취합")

    wb = Workbook()
    wb.remove(wb.active)
    sheet_summary(wb, recs)
    sheet_web(wb, recs)
    sheet_full(wb, recs)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path(a.out) if a.out else OUT_DIR / f"{ts}_현장조사_취합검토_{len(recs)}건.xlsx"
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out)
    print(f"\n저장(취합·검토): {out}")

    flat = OUT_DIR / f"{ts}_현장조사_일괄취합_{len(recs)}건.xlsx"
    build_flat_workbook(recs, flat)
    print(f"저장(일괄취합): {flat}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
