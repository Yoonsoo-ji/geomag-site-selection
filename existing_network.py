# -*- coding: utf-8 -*-
"""
기존 1등 지자기점 관측망 30점 — 단일 출처
==========================================

`make_survey_map.py`(survey_review.html)와 `geomag_site_selection.py`
(입지 선정 지도의 밀도 설계)가 **같은 관측망**을 봐야 한다. 두 곳에 명단을
두면 조용히 갈라진다(CLAUDE.md 12-G 의 별칭 사고와 같은 함정).

⚠️ **`docs/data/existing_sites.geojson` 의 33개를 그대로 쓰면 안 된다.**
그 파일은 성과표 두 개(22~25 · '10~'19)를 도엽명으로 병합한 **자료 보유
목록**이지 관측망이 아니다. 실제로 옛 이름 3건이 따로 실려 있는데, 그중
둘은 **좌표가 다르다** — 즉 폐지된 옛 위치가 현 측점과 나란히 앉는다:

| 옛 이름(연도) | 현 이름(연도) | 거리 |
|---|---|---|
| 언양 (2022) | 양산 (2022) | 9 m — 같은 표석, 단순 중복 |
| **경주 (2015)** | **영천 (2024)** | **23.2 km** |
| **임계 (2016)** | **삼척 (2024)** | **31.1 km** |

경주·임계를 관측망에 포함하면 **경북·강원에 있지도 않은 측점이 생겨**
공백이 실제보다 작게 나온다. 밀도 설계에서 하필 그 두 권역이 쟁점이므로
결과를 직접 왜곡한다. `SUPERSEDED` 로 명시해 배제한다.
"""
from __future__ import annotations

from pathlib import Path

# 기존 지자기 측정점 30점 (2026-08-19 확정) — 1등 지자기점 관측망
EXISTING_NETWORK = [
    "남양", "춘양", "청송", "상주", "미원", "서산", "청양", "이원", "부안",
    "와도", "순창", "강화", "함양", "가야", "영천", "언양", "남지", "거제",
    "남해", "순천", "장흥", "조도", "포천", "성산", "여주", "화천", "설악",
    "봉평", "삼척", "제천",
]

# 기존 30점 중 **선점 대상 15점**.
# 나머지 15점은 저수지(제방)에 설치된 점으로 판정되어 선점에서 뺀다 —
# 지반 거동·수위 변동으로 표석 지속성이 확보되지 않는다.
EXISTING_TARGET = [
    "상주", "미원", "부안", "강화", "언양", "남지", "거제", "조도",
    "포천", "성산", "여주", "화천", "봉평", "삼척", "제천",
]

# 확정 이름이 existing_sites.geojson 과 직접 일치하므로 별칭 치환이 필요 없다.
EXISTING_ALIAS: dict[str, str] = {}

# 옛 이름 → 현 이름. 성과표·야장에는 옛 이름으로 남아 있으나 관측망에서는
# 현 이름 쪽 좌표만 쓴다(`lmm_build.SITE_ALIAS` 와 같은 대응).
SUPERSEDED = {"양산": "언양", "경주": "영천", "임계": "삼척"}

# ── 좌표 정본 — 지자기 점조서 기준 30점 (2026-08-27 확정) ──────────────
#
# `data/geomag_network_30.csv`. 지점명·지점코드·도엽번호·
# 위경도(도분초)·소재지를 담은 30행이며, 지점코드(1-1 ~ 1-30)가 관측망 자체의
# 일련번호다.
#
# ⚠️ **「원장」·「정본」이라고 부르지 말 것** (2026-08-27 검토 지적). 원본 스캔·
# 공문·추출 로그가 저장소에 없어 원 문서와의 문자 단위 대조가 불가능하다.
# 확인된 것은 ① 내부 정합성(30행·30이름·좌표중복 0·도엽 포함 28/30) ②
# `'10~'19년 지자기점 관측현황(최종).xls`(표제 «지자기점 현황관리»)와
# **점의번호 30/30 · 도엽명 30/30 일치** — 즉 **같은 문서계열의 다른 판본**이라는
# 것까지다. 출처 3항목은 `geomag_network_30.meta.json` 에 공란으로 남겨 뒀다.
#
# ⚠️ **`'10~'19년 지자기점 관측현황(최종).xls` 의 좌표에는 전사 오류가 있다.**
# 대조표와 견주니 **주소는 같은데 좌표만 다른** 세 건이 나왔다 — 이설이 아니라
# 옮겨 적으면서 틀린 것이다:
#
#   | 지점 | 차이 | 내용 |
#   |---|---|---|
#   | **서산** | **52.0 km** | 남양 좌표가 통째로 복사돼 있었다(중복). 참값은 서산시 잠홍저수지 |
#   | **남양** | 3.6 km | 경도 분·초가 틀림(126°54′46″ vs 대조표 126°57′10.9″) |
#   | **와도** | 1.45 km | 경도 초 단위 오차 |
#
# 남양·서산이 같은 값이었던 것은 **두 행이 함께 오염**된 결과다. 서산은 남양의
# (그것도 이미 틀린) 값을 물려받았다. 대조표를 넣으면 둘 다 제자리로 간다.
#
# ⚠️ 반면 **22~25 성과표에 있는 15점은 대조표보다 최신**이다. 대조표와 크게
# 어긋나는 세 건은 **주소부터 다르므로 실제 이설**이다 — 대조표로 덮으면 안 된다:
#
#   포천 10.3 km (동두천시 보산동 → 연천군 청산면) · 여주 7.0 km (여주시 대신면
#   → 양평군 개군면) · 양산 21.3 km (미호리 산89-2 → 산59-2)
#
# → **규칙: 22~25 성과표에 있으면 성과표 좌표, 없으면 대조표 좌표.**
REGISTER_CSV = "data/geomag_network_30.csv"

# ── 성과표보다 **점조서를 우선**하는 측점 ────────────────────────────
#
# 기본 규칙은 「22~25 성과표에 있으면 성과표 좌표」다. 그 성과표가 최신 실측이기
# 때문인데, 아래 둘은 예외다(2026-08-27 검증).
#
#   여주 — 성과표 경도가 **1′ 크게 전사**됐다. 점조서 127°34′23.9581″ 대
#          성과표 127°35′23.96″ 로 **초까지 같고 분만 다르다.** 점조서 기재
#          표고 177.03 m 에 대해 DEM 조회가 점조서 좌표 171 m(차 6 m) ·
#          성과표 좌표 284 m(차 107 m) 를 주어 점조서가 맞음이 확인됐다.
#          점조서 경로도 「주읍리 마을회관 산수유꽃길」(개군면)로 주소와 맞는다.
#
#   언양 — 양산이 **언양으로 개명·대체**된 점이다(1-23-1, 미호리 산3,
#          2024.10.16 매설). 22~25 성과표(2022)는 그 이전이라 옛 산59-2
#          좌표를 담고 있다. 현재 표석은 점조서 값이다.
FORCE_REGISTER = {"여주", "언양"}


def load_register(root=None):
    """
    지자기점 관측망 30점 좌표 정본(점조서 기준)을 DataFrame 으로.
    컬럼: 지점명·지점코드·도엽번호·
    위도·경도(십진도)·소재지.
    """
    import pandas as pd

    base = Path(root) if root else Path(__file__).parent
    df = pd.read_csv(base / REGISTER_CSV)
    df["위도"] = df.위도_도 + df.위도_분 / 60 + df.위도_초 / 3600
    df["경도"] = df.경도_도 + df.경도_분 / 60 + df.경도_초 / 3600
    return df[["지점명", "지점코드", "도엽번호", "위도", "경도", "소재지"]]


def apply_register(df, name_col="도엽명", lat_col="위도", lon_col="경도",
                   only_names=None, root=None):
    """
    측점표의 좌표를 **대조표 값**으로 교정한다.

    `only_names` 를 주면 그 이름들만 손댄다 — **22~25 성과표에서 온 행은
    최신 실측이므로 제외**해야 한다(포천·여주·양산의 실제 이설을 되돌리게 된다).

    반환은 `(교정된 df, [(이름, 이동거리_km), ...])`.
    """
    import numpy as np
    import pandas as pd

    reg = load_register(root).set_index("지점명")
    out = df.copy()
    moved = []
    for i, row in out.iterrows():
        nm = str(row[name_col]).strip()
        if nm not in reg.index:
            continue
        if only_names is not None and nm not in only_names:
            continue
        la0 = pd.to_numeric(row[lat_col], errors="coerce")
        lo0 = pd.to_numeric(row[lon_col], errors="coerce")
        la1, lo1 = float(reg.at[nm, "위도"]), float(reg.at[nm, "경도"])
        if np.isfinite(la0) and np.isfinite(lo0):
            m = np.radians((la0 + la1) / 2)
            d = float(np.hypot((lo1 - lo0) * 111.320 * np.cos(m),
                               (la1 - la0) * 110.574))
            if d > 0.005:                     # 5 m 이상만 기록
                moved.append((nm, round(d, 3)))
        out.at[i, lat_col] = la1
        out.at[i, lon_col] = lo1
    return out, moved


def drop_duplicate_coords(df, name_col: str = "도엽명",
                          lat_col: str = "위도", lon_col: str = "경도",
                          decimals: int = 5):
    """
    좌표가 겹치는 측점을 하나만 남긴다 — 공간 계산용 관측망.

    ⚠️ **행 수를 그대로 측점 수로 쓰면 안 된다.** 담당면적(보로노이) 계산에서
    같은 좌표의 두 점 중 하나는 담당 셀을 하나도 못 받아 조용히 사라지는데,
    정규화 분모에는 계속 세어져 밀도를 과대평가한다.

    대조표 좌표를 적용한 뒤에는 중복이 없어야 한다 — 걸리면 새 결함이므로
    호출부에서 반드시 알릴 것.
    """
    import pandas as pd

    d = df.copy()
    d["_la"] = pd.to_numeric(d[lat_col], errors="coerce").round(decimals)
    d["_lo"] = pd.to_numeric(d[lon_col], errors="coerce").round(decimals)
    keep = d.drop_duplicates(["_la", "_lo"], keep="first")
    dropped = [str(n).strip() for n in
               d.loc[~d.index.isin(keep.index), name_col].tolist()]
    keep = keep.drop(columns=["_la", "_lo"])
    return keep, dropped


def network_names(target_only: bool = False) -> list[str]:
    """관측망 명단. `target_only=True` 면 선점 대상 15점만."""
    return list(EXISTING_TARGET if target_only else EXISTING_NETWORK)


def select_rows(df, name_col: str = "도엽명", target_only: bool = False):
    """
    측점 DataFrame 에서 관측망 행만 남긴다.

    `geomag_site_selection.load_existing_sites()` 의 반환값(도엽명·위도·경도)을
    그대로 받도록 만들었다. 명단에 있으나 표에 없는 이름은 조용히 빠지므로
    호출부에서 반환 길이를 확인할 것.
    """
    names = set(network_names(target_only))
    return df[df[name_col].astype(str).str.strip().isin(names)]


def missing_from(df, name_col: str = "도엽명", target_only: bool = False) -> list[str]:
    """명단에 있으나 표에 없는 이름 — 자료 누락 점검용."""
    have = set(df[name_col].astype(str).str.strip())
    return [n for n in network_names(target_only) if n not in have]


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8")
    print(f"관측망 {len(EXISTING_NETWORK)}점 · 선점 대상 {len(EXISTING_TARGET)}점")
    print(f"옛 이름 배제 {len(SUPERSEDED)}건: "
          + " · ".join(f"{k}→{v}" for k, v in SUPERSEDED.items()))
    reg = load_register()
    print(f"\n좌표 정본 {len(reg)}점 — 지자기 점조서 기준 ({REGISTER_CSV})")
    dup = reg.round({"위도": 5, "경도": 5}).duplicated(["위도", "경도"]).sum()
    print(f"  좌표 중복 {int(dup)}건")
    miss = [n for n in EXISTING_NETWORK if n not in set(reg.지점명)]
    extra = [n for n in reg.지점명 if n not in set(EXISTING_NETWORK)]
    print(f"  명단에만 있음: {miss or '없음'}")
    print(f"  대조표에만 있음: {extra or '없음'}  (옛 이름 — SUPERSEDED 대응)")
    rest = [n for n in EXISTING_NETWORK if n not in EXISTING_TARGET]
    print(f"\n선점 대상 15: {' · '.join(EXISTING_TARGET)}")
    print(f"저수지 제외 {len(rest)}: {' · '.join(rest)}")
