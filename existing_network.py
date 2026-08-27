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

# 기존 지자기 측정점 30점 (2026-08-19 확정) — 1등 지자기점 관측망
EXISTING_NETWORK = [
    "남양", "춘양", "청송", "상주", "미원", "서산", "청양", "이원", "부안",
    "와도", "순창", "강화", "함양", "가야", "영천", "양산", "남지", "거제",
    "남해", "순천", "장흥", "조도", "포천", "성산", "여주", "화천", "설악",
    "봉평", "삼척", "제천",
]

# 기존 30점 중 **선점 대상 15점**.
# 나머지 15점은 저수지(제방)에 설치된 점으로 판정되어 선점에서 뺀다 —
# 지반 거동·수위 변동으로 표석 지속성이 확보되지 않는다.
EXISTING_TARGET = [
    "상주", "미원", "부안", "강화", "양산", "남지", "거제", "조도",
    "포천", "성산", "여주", "화천", "봉평", "삼척", "제천",
]

# 확정 이름이 existing_sites.geojson 과 직접 일치하므로 별칭 치환이 필요 없다.
EXISTING_ALIAS: dict[str, str] = {}

# 옛 이름 → 현 이름. 성과표·야장에는 옛 이름으로 남아 있으나 관측망에서는
# 현 이름 쪽 좌표만 쓴다(`lmm_build.SITE_ALIAS` 와 같은 대응).
SUPERSEDED = {"언양": "양산", "경주": "영천", "임계": "삼척"}

# ── 좌표 결함 — 원본 성과표의 오류 (2026-08-27 확인) ─────────────────────
#
# ⚠️ **`'10~'19년 지자기점 관측현황(최종).xls` 에서 남양과 서산이 같은 좌표를
# 쓴다** — 둘 다 `126° 54′ 46″ / 37° 06′ 53″`. 전 파일에서 유일한 중복이다.
# 소재지는 서로 다르다(남양 = 화성 동탄면 장지리저수지 · 서산 = 서산 잠홍동
# 잠홍저수지, 실거리 약 70 km).
#
# **어느 쪽이 옳은지는 도엽으로 판별된다.** 측점명은 1:50,000 도엽명이고,
# 문제의 좌표는 **「남양」 도엽(NJ52-7-28) 안**에 떨어진다. 즉 남양이 맞고
# **서산 행이 남양 좌표를 복사해 온 것**이다.
#
# 서산의 참 좌표는 이 저장소의 어느 자료에도 없다(파생 파일들은 전부 같은
# 오류를 물려받았다). **지리원 원장에서 확인해 채워 넣어야 한다.**
#
# ⚠️ 그때까지 서산은 **공간 계산에서 빠진다** — 관측망은 명목 30점이나
# 공간적으로는 **29점**이다. 밀도 정규화에 30 을 쓰면 없는 측점을 세는 것이다.
BAD_COORDS = {
    "서산": ("남양 좌표와 동일(126°54′46″/37°06′53″) — 그 좌표는 남양 도엽 "
             "안이므로 서산 쪽이 오기. 지리원 원장 확인 필요"),
}


def drop_duplicate_coords(df, name_col: str = "도엽명",
                          lat_col: str = "위도", lon_col: str = "경도",
                          decimals: int = 5):
    """
    좌표가 겹치는 측점을 하나만 남긴다 — 공간 계산용 관측망.

    ⚠️ **행 수를 그대로 측점 수로 쓰면 안 된다.** 담당면적(보로노이) 계산에서
    같은 좌표의 두 점 중 하나는 담당 셀을 하나도 못 받아 조용히 사라지는데,
    정규화 분모에는 계속 세어져 밀도를 과대평가한다.

    `BAD_COORDS` 에 사유가 적힌 이름을 우선 버리고, 없으면 뒤에 오는 행을
    버린다. 반환은 `(정리된 df, 버린 이름 목록)`.
    """
    import pandas as pd  # 지연 import — 이 모듈은 명단만으로도 쓰인다

    d = df.copy()
    d["_la"] = pd.to_numeric(d[lat_col], errors="coerce").round(decimals)
    d["_lo"] = pd.to_numeric(d[lon_col], errors="coerce").round(decimals)
    # 결함이 알려진 이름을 뒤로 보내 drop_duplicates(keep="first") 에서 탈락시킨다
    d["_bad"] = d[name_col].astype(str).str.strip().isin(BAD_COORDS).astype(int)
    d = d.sort_values("_bad", kind="stable")
    keep = d.drop_duplicates(["_la", "_lo"], keep="first")
    dropped = [str(n).strip() for n in
               d.loc[~d.index.isin(keep.index), name_col].tolist()]
    keep = keep.sort_index().drop(columns=["_la", "_lo", "_bad"])
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
    print(f"\n좌표 결함 {len(BAD_COORDS)}건 — 공간 계산에서 제외:")
    for k, why in BAD_COORDS.items():
        print(f"  {k}: {why}")
    print(f"  → 명목 {len(EXISTING_NETWORK)}점 / 공간 "
          f"{len(EXISTING_NETWORK) - len(BAD_COORDS)}점")
    rest = [n for n in EXISTING_NETWORK if n not in EXISTING_TARGET]
    print(f"\n선점 대상 15: {' · '.join(EXISTING_TARGET)}")
    print(f"저수지 제외 {len(rest)}: {' · '.join(rest)}")
