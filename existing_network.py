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
    rest = [n for n in EXISTING_NETWORK if n not in EXISTING_TARGET]
    print(f"\n선점 대상 15: {' · '.join(EXISTING_TARGET)}")
    print(f"저수지 제외 {len(rest)}: {' · '.join(rest)}")
