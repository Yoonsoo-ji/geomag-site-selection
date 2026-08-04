# -*- coding: utf-8 -*-
"""
KASA 우주환경 빅데이터 플랫폼 지자기 관측자료 수집기.

    python kasa_geomag.py --probe                 # 연도별 커버리지 조사 -> 리포트
    python kasa_geomag.py --station three --start 2022-06-14 --end 2022-06-30
    python kasa_geomag.py --station boh   --start 2020-01-01 --end 2020-12-31

우주환경 빅데이터 플랫폼(spaceweather.kasa.go.kr)의 데이터셋 다운로드 경로를 쓴다.
Open API(/api/magnetism/*)는 최근 ~1일 구간만 주므로 과거자료에는 데이터셋 경로를 쓴다.

인증키·로그인 불필요. 흐름:
  POST dataSetSampleListExcelDownAjax.do (form) -> {filePath, fileName}
  GET  common/fileDownload.do?filePath=..&fileName=..  -> XLSX (확장자만 .csv)

주의
  · 반환 파일은 이름만 .csv 이고 실제로는 XLSX(ZIP). inlineStr + dimension 태그가
    깨져 있어 openpyxl read_only 로는 못 읽는다. sheet1.xml 을 직접 파싱한다.
  · 요청당 기간 상한: 3관측소(제주·강릉·이천)=5일, 보현산(천문연)=1개월.
  · 이천(ICH)은 관측기 점검으로 2023-10 이후 신규 데이터 중단.
  · downFile.csv 는 세션에 묶인 임시파일 → POST 직후 같은 세션으로 GET.
  · 정부 서버이므로 요청 사이 지연(REQ_DELAY)을 둔다.
"""
import argparse
import io
import re
import sys
import time
import urllib.parse
import zipfile
import datetime as dt
from pathlib import Path

import requests

BASE = "https://spaceweather.kasa.go.kr"
HERE = Path(__file__).parent
OUT_DIR = HERE / "data" / "kasa"
REQ_DELAY = 0.8          # 초, 요청 간 간격
RETRY = 3
MISSING = "-99999"

# ── 관측소 레지스트리 ──────────────────────────────────────────────────────
# region_map: 파일 내 지역코드 -> 한글명 (3관측소 세트 전용)
STATIONS = {
    "three": {
        "tbl": "TB_MI_MNT_AVG_MAG",
        "inst": "KSWC",
        "dataNm": "국내 지자기 관측데이터(제주,강릉,이천)",
        "max_days": 5,
        "kind": "xyzf",                 # REGN_CD, OBSR_YMD, X, Y, Z, F, ...품질
        "region_map": {"JJ": "제주", "GN": "강릉", "ICH": "이천"},
    },
    "boh": {
        "tbl": "TB_DG_L0_KOR_KASI_BOH_GEOMAG_1MNT",
        "inst": "KASI",
        "dataNm": "국내천문연보현산지자기1분",
        "max_days": 31,
        "kind": "hdz",                  # OBSR_YMD, H, D, Z, ...  (규약 검증 필요)
        "region_map": None,
    },
    "boh_k": {
        "tbl": "TB_DG_L0_KOR_KASI_GEOMAG_BOH_K",
        "inst": "KASI",
        "dataNm": "국내천문연보현산K지수",
        "max_days": 31,
        "kind": "raw",
        "region_map": None,
    },
}


def new_session():
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Mozilla/5.0 (research data collection; LX geomag)",
        "X-Requested-With": "XMLHttpRequest",
        "Referer": BASE + "/openpotal/datasetInfo/dataSetList.do",
    })
    s.get(BASE + "/openpotal/datasetInfo/dataSetList.do", timeout=30)
    return s


def _fetch_xlsx(s, st, start, end):
    """단일 요청(<= max_days). sheet1.xml 행 리스트 반환. 데이터 없으면 []."""
    form = {
        "obsrType": "", "imgNm": "",
        "tableEngNm": st["tbl"], "tblEngNm": st["tbl"],
        "pvsnInstCd": st["inst"], "dataNm": st["dataNm"], "cnncManageNo": "stdb",
        "startDate": start, "endDate": end,
    }
    for attempt in range(RETRY):
        try:
            r = s.post(BASE + "/openpotal/datasetInfo/dataSetSampleListExcelDownAjax.do",
                       data=form, timeout=90)
            j = r.json()
            if j.get("resultCode") != "00":
                return []                       # 01 = 해당 기간 데이터 없음
            fp, fn = j["filePath"], j["fileName"]
            dl = (BASE + "/common/fileDownload.do?filePath="
                  + urllib.parse.quote(urllib.parse.quote(fp + fn, safe=""))
                  + "&fileName=" + urllib.parse.quote(urllib.parse.quote(fn, safe="")))
            content = s.get(dl, timeout=180).content
            xml = zipfile.ZipFile(io.BytesIO(content)).read(
                "xl/worksheets/sheet1.xml").decode("utf-8", "replace")
            rows = [_parse_row(m.group(1))
                    for m in re.finditer(r"<row[^>]*>(.*?)</row>", xml, re.S)]
            # 헤더(한글 컬럼명이 든 행) 제거: 첫 셀이 날짜형/지역코드가 아니면 헤더
            data = [rr for rr in rows
                    if rr and (re.match(r"\d{4}-\d\d-\d\d", rr[0] or "")
                               or (rr[0] in ("JJ", "GN", "ICH")))]
            return data
        except Exception as e:                  # noqa: BLE001
            if attempt == RETRY - 1:
                print(f"    ! 실패 {start}~{end}: {e}", file=sys.stderr)
                return []
            time.sleep(2.0)
    return []


def _col_idx(ref):
    """셀 ref('B12')의 열문자를 0-기반 인덱스로. 빈 셀로 인한 열 밀림 방지."""
    letters = re.match(r"([A-Z]+)", ref).group(1)
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - 64)
    return idx - 1


def _parse_row(inner):
    """<row> 안의 셀들을 열 위치 기준으로 배열한다(빈 셀은 '' 로 채움)."""
    cells = {}
    maxc = -1
    for cm in re.finditer(r'<c r="([A-Z]+\d+)"[^>]*>(.*?)</c>', inner, re.S):
        ci = _col_idx(cm.group(1))
        tm = re.search(r"<t[^>]*>(.*?)</t>", cm.group(2), re.S)
        cells[ci] = tm.group(1) if tm else ""
        maxc = max(maxc, ci)
    return [cells.get(i, "") for i in range(maxc + 1)]


def _daterange_chunks(start, end, max_days):
    d0 = dt.date.fromisoformat(start)
    d1 = dt.date.fromisoformat(end)
    step = dt.timedelta(days=max_days - 1)      # 경계 포함이므로 -1
    cur = d0
    while cur <= d1:
        seg_end = min(cur + step, d1)
        yield cur.isoformat(), seg_end.isoformat()
        cur = seg_end + dt.timedelta(days=1)


def _valid_counts(st, rows):
    """유효(비결측) 레코드 수를 반환. 3관측소는 지역별 dict."""
    if st["region_map"]:
        cnt = {k: 0 for k in st["region_map"]}
        for r in rows:
            reg = r[0]
            if reg in cnt and len(r) > 5 and r[2] not in (MISSING, "", None):
                cnt[reg] += 1
        return cnt
    # 보현산: 2열(H)이 유효한 행 수
    n = sum(1 for r in rows if len(r) > 3 and r[1] not in (MISSING, "", None))
    return {"보현산": n}


# ── 커버리지 조사 ─────────────────────────────────────────────────────────
def probe():
    s = new_session()
    years = list(range(2018, 2027))
    # 각 연도 대표 표본일 (야장이 있는 5·9·10월 근처 + 연중 기준)
    print("KASA 지자기 관측자료 — 연도별 커버리지 조사")
    print("(각 연도 표본 2일치 요청, 유효 분(minute) 레코드 수)\n")
    for key in ("three", "boh"):
        st = STATIONS[key]
        label = "제주·강릉·이천" if key == "three" else "보현산(1분)"
        print(f"■ {label}  [{st['tbl']}]")
        for y in years:
            # 표본: 6/15~6/16 (없으면 데이터 없음으로 집계)
            a = f"{y}-06-15"
            b = f"{y}-06-16"
            rows = _fetch_xlsx(s, st, a, b)
            time.sleep(REQ_DELAY)
            vc = _valid_counts(st, rows)
            tot = sum(vc.values())
            mark = "―" if tot == 0 else " ".join(f"{k}:{v}" for k, v in vc.items())
            print(f"   {y}  {mark}")
        print()
    # 2019 야장 시기 정밀 확인 (5·9·10월)
    print("■ 2019 야장 시기 정밀 확인 (제주·강릉·이천 / 보현산)")
    for a, b, tag in [("2019-05-14", "2019-05-17", "미원·이원 시기"),
                      ("2019-09-18", "2019-09-19", "성산 시기"),
                      ("2019-10-21", "2019-10-23", "장흥 시기")]:
        r3 = _fetch_xlsx(s, STATIONS["three"], a, b); time.sleep(REQ_DELAY)
        rb = _fetch_xlsx(s, STATIONS["boh"], a, b); time.sleep(REQ_DELAY)
        c3 = sum(_valid_counts(STATIONS["three"], r3).values())
        cb = sum(_valid_counts(STATIONS["boh"], rb).values())
        print(f"   {a}~{b} ({tag})  3관측소:{c3}  보현산:{cb}")


# ── 실제 수집 ─────────────────────────────────────────────────────────────
# 한글 Excel(CP949)에서 열어도 안 깨지도록 UTF-8 BOM(utf-8-sig) 로 저장한다.
# 파일은 관측소별로 나눈다 — 3관측소를 한 파일에 합치면 연 단위가 약 158만 행이라
# Excel 행 한도(1,048,576)를 넘긴다. 관측소별이면 연 최대 ~52.7만 행으로 안전.
ENC = "utf-8-sig"


def collect(key, start, end):
    st = STATIONS[key]
    s = new_session()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gaps = []

    # 원본(API/데이터셋)은 최신→과거 역순으로 내려준다. 청크별로 버퍼에 모은 뒤
    # 시각 오름차순으로 정렬해 저장한다(2019-01-01 00:00 → 12-31 23:59).
    if st["region_map"]:
        header = "station,time_kst,X_nT,Y_nT,Z_nT,F_nT\n"
        buf = {code: [] for code in st["region_map"]}    # code -> [(time, line)]
        for a, b in _daterange_chunks(start, end, st["max_days"]):
            rows = _fetch_xlsx(s, st, a, b)
            time.sleep(REQ_DELAY)
            vc = _valid_counts(st, rows)
            if sum(vc.values()) == 0:
                gaps.append((a, b))
                print(f"  · {a}~{b}  데이터 없음")
                continue
            for r in rows:
                code = r[0]
                if code not in buf or len(r) < 6:
                    continue
                reg = st["region_map"][code]
                buf[code].append((r[1], f"{reg},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]}\n"))
            print(f"  · {a}~{b}  " + " ".join(f"{k}:{v}" for k, v in vc.items()))
        out = []
        for code, items in buf.items():
            items.sort(key=lambda t: t[0])               # 시각 오름차순
            p = OUT_DIR / f"kasa_{code}_{start}_{end}.csv"
            with open(p, "w", encoding=ENC, newline="") as fh:
                fh.write(header)
                fh.writelines(line for _, line in items)
            print(f"저장: {p.name}  ({len(items)} rows)")
            out.append(p)
    else:
        out = OUT_DIR / f"kasa_BOH_{start}_{end}.csv"
        buf, real = [], 0
        for a, b in _daterange_chunks(start, end, st["max_days"]):
            rows = _fetch_xlsx(s, st, a, b)
            time.sleep(REQ_DELAY)
            got = 0
            for r in rows:
                if len(r) < 4 or not re.match(r"\d{4}-\d\d-\d\d", r[0]):
                    continue
                # 상수 placeholder(H=D=Z) 는 그대로 기록(무효는 분석 단계에서 판별)
                buf.append((r[0], f"보현산,{r[0]},{r[1]},{r[2]},{r[3]}\n"))
                got += 1
                if not (r[1] == r[2] == r[3]):
                    real += 1
            if got == 0:
                gaps.append((a, b))
            print(f"  · {a}~{b}  기록 {got} (실측 {real})")
        buf.sort(key=lambda t: t[0])
        with open(out, "w", encoding=ENC, newline="") as fh:
            fh.write("station,time_kst,H,D,Z\n")
            fh.writelines(line for _, line in buf)
        print(f"저장: {out.name}  ({len(buf)} rows)")

    if gaps:
        print(f"결측 구간 {len(gaps)}건:")
        for a, b in gaps:
            print(f"   - {a} ~ {b}")
    return out, gaps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--station", choices=list(STATIONS))
    ap.add_argument("--start")
    ap.add_argument("--end")
    args = ap.parse_args()
    if args.probe:
        probe()
    elif args.station and args.start and args.end:
        collect(args.station, args.start, args.end)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
