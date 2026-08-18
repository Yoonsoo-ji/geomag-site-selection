# -*- coding: utf-8 -*-
"""
지자기 야장(野帳) 파서 — 관측 일시 추출
========================================

외부장(External) 보정은 "측정한 그 순간"의 변동을 빼는 작업이므로
관측 일시가 반드시 필요하다. 성과표(지자기측량 성과정리)에는 연도만 있어
야장에서 일시를 복원해야 한다.

야장 양식(MAG-01H)
------------------
세로로 관측 블록이 반복되고, 가로로는 3열 간격으로 관측회차가 배치된다.

    행 1~2   점번호 · 지자기측점명 · 관측일자 · 도엽명
    행 4     일자      (관측회차별 날짜)
    행 14    관측시간   (편각 세션 시작)
    행 16~21 E D / W D / E U / W U 각각의 시각과 측정값
    행 28~36 S U / N D / S D / N U (복각)
    행 38    자력 F
    행 40~41 편각 · 복각 결과

행 위치가 파일마다 1~2행씩 밀리므로 고정 좌표 대신 라벨을 탐색한다.

실행:
    python lmm_fieldbook.py
"""

import datetime as dt
import glob
import re
import sys
from pathlib import Path

import pandas as pd

BASE = Path(__file__).parent
FIELDBOOK_DIR = BASE / "docs" / "data"

# 국토지리정보원 야장 원본 트리 (2026-08 확보). 연도·측점별 폴더로 정리돼 있고
# 2022~2025 분에도 **분 단위 관측시각이 기록돼 있다** — External 층 봉쇄를 푸는 자료.
# 파일명이 제각각이라(「야장」이 빠진 것도 있다) 확장자로 훑고 내용으로 판별한다.
NGII_FIELDBOOK_DIR = Path(
    r"D:\LX_yoons\2026_research\2026_지자기 연구"
    r"\20260811_지리원 야장 자료\지구물리측량\지자기측량")

# 야장 파일 목록(파일명이 제각각이라 패턴으로 수집)
PATTERNS = ["*야장*.xlsx"]

# NGII 트리에서 제외할 것 — 점조서·성과 대장은 야장이 아니고 수십 MB 라 파싱이 느리다
NGII_SKIP_WORDS = ("조서", "점조서", "성과목록", "기본측량성과")
NGII_MAX_BYTES = 3_000_000

# 무시할 시트(양식·참고자료)
SKIP_SHEETS = {"야장양식", "간편야장양식", "스마트전자지도web자료",
               "2010연구성과보고서결과", "극좌표변환"}


def _norm(x) -> str:
    return re.sub(r"\s+", "", str(x))


def _is_time(v) -> bool:
    """실제 관측시각인지 판정. 00:00:00 은 미기입으로 본다."""
    if isinstance(v, dt.time):
        return not (v.hour == 0 and v.minute == 0 and v.second == 0)
    return False


def _is_date(v):
    if isinstance(v, (dt.datetime, pd.Timestamp)):
        return pd.Timestamp(v).normalize()
    return None


def parse_sheet(df: pd.DataFrame):
    """
    한 시트에서 (측점명, 날짜, 시각목록, F) 세션들을 추출.

    블록 경계는 '일    자' 라벨 행으로 잡고, 그 아래 다음 블록 전까지를
    한 블록으로 본다. 블록 안에서 열별로 날짜·시각을 묶는다.
    """
    nrow, ncol = df.shape
    labels = {r: _norm(df.iat[r, 0]) for r in range(nrow)}

    # 측점명: '지자기측점명' 또는 '도엽(측점)' 라벨 바로 아래 행
    site = None
    for r in range(min(6, nrow)):
        row = [_norm(df.iat[r, c]) for c in range(min(ncol, 12))]
        for key in ("지자기측점명", "도엽(측점)"):
            if key in row and r + 1 < nrow:
                site = df.iat[r + 1, row.index(key)]
                break
        if site is not None:
            break

    # 한 시트에 '일자' 표와 '측정일자' 표가 연달아 오므로 둘 다 경계로 삼는다.
    # (경계를 놓치면 다음 표의 시각까지 앞 세션에 섞여 들어간다)
    BLOCK_LABELS = ("일자", "측정일자", "관측일자")
    block_rows = [r for r in range(nrow) if labels[r] in BLOCK_LABELS]

    if not block_rows:
        # 성산 양식: 행1이 머리글, 행2가 값 (일자 라벨 행이 없음)
        for r in range(min(4, nrow)):
            row = [_norm(df.iat[r, c]) for c in range(min(ncol, 30))]
            if "관측일자" in row and r + 1 < nrow:
                block_rows = [r + 1]
                break
    if not block_rows:
        return []

    sessions = []
    for bi, br in enumerate(block_rows):
        end = block_rows[bi + 1] if bi + 1 < len(block_rows) else nrow

        # 열별 날짜
        col_date = {}
        for c in range(1, ncol):
            d = _is_date(df.iat[br, c])
            if d is not None:
                col_date[c] = d
        if not col_date:
            continue

        # 시각 수집은 '자력'(F) 행에서 끝난다 — DI 관측 시각은 모두 그 앞에 있다
        f_rows = [r for r in range(br, end) if labels[r] == "자력"]
        time_end = f_rows[0] if f_rows else end

        # 편각/복각 구간을 나눈다. 야장에 따라 두 구간 시각이 어긋난 사례가 있어
        # (장흥 2019-10-22 6회차: 복각 13:10대 < 편각 14:03대) 따로 봐야 한다.
        inc_rows = [r for r in range(br, time_end) if "복각" in labels[r]]
        inc_start = inc_rows[0] if inc_rows else time_end

        col_times = {c: [] for c in col_date}
        col_dtimes = {c: [] for c in col_date}
        col_itimes = {c: [] for c in col_date}
        for r in range(br, time_end):
            for c in range(1, ncol):
                if _is_time(df.iat[r, c]):
                    # 날짜가 있는 열 중 가장 가까운 왼쪽 열에 귀속
                    owner = max([k for k in col_date if k <= c], default=None)
                    if owner is None:
                        continue
                    col_times[owner].append(df.iat[r, c])
                    (col_dtimes if r < inc_start else col_itimes)[owner].append(
                        df.iat[r, c])

        # 자력(F)
        col_F = {}
        if f_rows:
            for c in col_date:
                v = df.iat[f_rows[0], c]
                if isinstance(v, (int, float)) and pd.notna(v) and v > 1000:
                    col_F[c] = float(v)

        for c, d in col_date.items():
            ts = sorted(set(col_times[c]))
            dt_ = sorted(set(col_dtimes[c]))
            it_ = sorted(set(col_itimes[c]))

            # 편각·복각 구간이 30분 넘게 떨어지면 기입 오류 의심
            flag = ""
            if dt_ and it_:
                gap = (dt.datetime.combine(d.date(), it_[0])
                       - dt.datetime.combine(d.date(), dt_[-1])).total_seconds() / 60
                if gap < -5:
                    flag = f"복각시각이 편각보다 {abs(gap):.0f}분 빠름"
                elif gap > 30:
                    flag = f"편각-복각 간격 {gap:.0f}분"

            sessions.append({
                "측점": site,
                "날짜": d.date(),
                "관측수": len(ts),
                "시작": ts[0] if ts else None,
                "종료": ts[-1] if ts else None,
                "편각구간": f"{dt_[0]:%H:%M}~{dt_[-1]:%H:%M}" if dt_ else None,
                "복각구간": f"{it_[0]:%H:%M}~{it_[-1]:%H:%M}" if it_ else None,
                "이상": flag,
                "F_nT": col_F.get(c),
            })
    return sessions


def parse_file(path: Path):
    out = []
    try:
        xl = pd.ExcelFile(path)
    except Exception as e:
        print(f"  [열기 실패] {path.name}: {e}")
        return out

    for s in xl.sheet_names:
        if s in SKIP_SHEETS:
            continue
        try:
            df = xl.parse(s, header=None)
        except Exception:
            continue
        if df.shape[0] < 10:
            continue
        for rec in parse_sheet(df):
            rec["파일"] = path.name
            rec["시트"] = s
            out.append(rec)
    return out


def ngii_files():
    """지리원 야장 트리에서 야장 후보 파일을 모은다(연도 폴더 하위 재귀)."""
    if not NGII_FIELDBOOK_DIR.exists():
        return []
    out = []
    for p in NGII_FIELDBOOK_DIR.rglob("*.xls*"):
        if p.name.startswith("~$"):
            continue
        if any(w in p.name for w in NGII_SKIP_WORDS):
            continue
        try:
            if p.stat().st_size > NGII_MAX_BYTES:
                continue
        except OSError:
            continue
        out.append(str(p))
    return sorted(out)


def collect():
    files = []
    for p in PATTERNS:
        files += glob.glob(str(FIELDBOOK_DIR / p))
    files += ngii_files()
    files = sorted(set(files))

    rows = []
    for f in files:
        rows += parse_file(Path(f))

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # 같은 측점·날짜·시간대는 시트 중복이므로 하나로
    df = df.dropna(subset=["측점", "날짜"])
    df = df.sort_values(["측점", "날짜", "시작"])
    df = df.drop_duplicates(subset=["측점", "날짜", "시작", "종료"], keep="first")

    # 시각 없는 행은, 같은 측점·날짜에 시각 있는 세션이 이미 있으면 중복이므로 제거
    timed = set(zip(df.loc[df["시작"].notna(), "측점"],
                    df.loc[df["시작"].notna(), "날짜"]))
    keep = df["시작"].notna() | ~df.apply(
        lambda r: (r["측점"], r["날짜"]) in timed, axis=1)
    df = df[keep]
    return df.reset_index(drop=True)


SESSION_CSV = BASE / "docs" / "data" / "fieldbook_sessions.csv"


def save_sessions(df):
    """세션표를 CSV 로 남긴다 — 원자료가 아니라 일시·측점만 담은 파생 자료."""
    cols = [c for c in ("측점", "날짜", "시작", "종료", "편각구간", "복각구간",
                        "관측수", "F_nT", "이상", "파일", "시트") if c in df.columns]
    out = df[cols].copy()
    out["날짜"] = out["날짜"].astype(str)
    SESSION_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(SESSION_CSV, index=False, encoding="utf-8-sig")
    print(f"\n세션표 저장: {SESSION_CSV}  ({len(out)}행)")


def main():
    print("=" * 74)
    print("지자기 야장 관측일시 추출")
    print("=" * 74)

    df = collect()
    if df.empty:
        print("추출된 세션 없음")
        return

    have = df[df["시작"].notna()]
    lack = df[df["시작"].isna()]

    print(f"\n총 세션 {len(df)}건  |  시각 있음 {len(have)}건  |  시각 없음 {len(lack)}건\n")

    show = df.copy()
    for col in ("편각구간", "복각구간"):
        show[col] = show[col].fillna("—")
    show["F_nT"] = show["F_nT"].apply(lambda v: f"{v:.1f}" if pd.notna(v) else "—")
    show["이상"] = show["이상"].fillna("")
    print(show[["측점", "날짜", "편각구간", "복각구간", "관측수", "F_nT", "이상"]]
          .to_string(index=False))

    print("\n--- 측점별 요약 ---")
    for site, g in df.groupby("측점"):
        ok = g["시작"].notna().sum()
        days = sorted(set(str(d) for d in g["날짜"]))
        print(f"  {site:6s} {len(g)}세션 (시각 {ok}건)  일자: {', '.join(days)}")

    if len(lack):
        print(f"\n[주의] 시각 미기입 {len(lack)}건 — 해당 세션은 외부장 보정 불가")
        for _, r in lack.iterrows():
            print(f"   · {r['측점']} {r['날짜']} ({r['파일']})")

    # 연도별 요약 — 성과표(연도만 기재)와 대조하기 위한 것
    df = df.copy()
    df["연도"] = df["날짜"].apply(lambda d: pd.Timestamp(d).year)
    print("\n--- 연도별 ---")
    for y, g in df.groupby("연도"):
        print(f"  {y}  측점 {g['측점'].nunique():2}종 · 세션 {len(g):3}건 "
              f"· 시각 {g['시작'].notna().sum():3}건")

    save_sessions(df)
    return df


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
