# -*- coding: utf-8 -*-
"""
현장 야장 기입 안내서 — docx · **hwpx** 두 벌
================================================

    python create_field_guide.py            # 둘 다
    python create_field_guide.py --hwpx     # 한글만
    python create_field_guide.py --docx     # 워드만

`make_trial_field_card.py` 가 내는 현장 카드를 «무엇을 어떻게 채우는가» 설명하는
인쇄물이다. 카드와 함께 들고 다닐 수 있게 A4 몇 장으로 끝낸다.

## 내용은 한 벌, 렌더러는 둘

`content()` 가 블록 목록을 내고 `render_docx()` · `render_hwpx()` 가 각각 그린다.
**본문을 두 번 적지 않는다** — 두 파일이 갈라지면 현장에서 어느 쪽을 믿어야 할지
알 수 없다.

⚠️ 수치와 규격도 여기서 다시 적지 않는다. 측선 거리·참고값·발주자 결정은 전부
`trial_survey_spec.py` 에서 읽는다.

## hwpx 는 한컴 COM 을 거친다

이 PC 에 한글이 깔려 있다(HWPFrame.HwpObject 10.0). 다만 **docx 는 못 연다**
(`Open()` 이 어떤 format 문자열로도 False). HTML 은 열리므로 **HTML 을 만들어
한글로 열고 hwpx 로 저장**한다.

| 함정 | 대응 |
|---|---|
| UTF-8 HTML 을 주면 본문이 통째로 깨진다 | **CP949 로 저장**하고 `charset=euc-kr` 을 명시 |
| CP949 에 없는 글자가 있다 | `—`(U+2014) → `―`, `−`(U+2212) → `-` 로 치환(`_cp949()`) |
| 보안 모듈 미등록이면 파일 접근마다 대화상자 | `RegisterModule("FilePathCheckDLL", "FilePathCheckerModule")` |

말투는 카드와 같이 **설명하는 투**로 쓴다 — 「~한다」로 끊지 않는다.
"""
from __future__ import annotations

import datetime as dt
import html as _html
import os
import sys
from pathlib import Path

import trial_survey_spec as SP

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "docs" / "output"

NAVY_HEX = "1F3864"
RED_HEX = "B2182B"
GREY_HEX = "555555"
CALLOUT_BG = "FDE9E7"


# ══════════════════════════════════════════════════════════════
def content():
    """안내서 블록 목록. 렌더러가 이것만 보고 그린다."""
    nH = len(SP.H_DIRECTIONS) * (len(SP.H_OFFSETS_M) + 2)
    nV = len(SP.V_HEIGHTS_CM) + 2
    step = SP.V_HEIGHTS_CM[1] - SP.V_HEIGHTS_CM[0]
    offs = ", ".join(str(x) for x in SP.H_OFFSETS_M)

    B = []
    A = B.append

    A(("cover", "지자기 시험 탐사 현장 야장", "기입 안내서",
       f"{dt.date.today():%Y년 %m월}   ·   선점 검토 대상 50지점"))
    A(("p", "이 안내서는 현장 야장 카드를 어떻게 채우는지 설명합니다. 카드는 점마다 "
            "한 장씩 있고, 한 장을 다 채우면 그 지점의 기록이 끝납니다. 카드와 함께 "
            "들고 다니시면 됩니다."))
    A(("callout",
       "이 탐사는 법정 선점이 아니고, 얼마 이하여야 합격이라는 기준도 아직 "
       "없습니다. 오히려 그 기준을 만들 자료를 모으는 것이 이번 일의 목적입니다. "
       "그래서 값이 크게 나오더라도 현장에서 그 자리를 접지 마시고, 나온 대로 "
       "적어 주시는 것이 가장 중요합니다."))

    A(("h1", "0. 먼저 알아 둘 것"))
    A(("h2", "카드에 나오는 말"))
    A(("table", ["말", "뜻"],
       [["P0", "중심점. 그 지점의 한가운데이고, 모든 거리를 여기서 잽니다."],
        ["Site ID", "지점 번호(S01 …). 카드마다 미리 찍혀 있습니다."],
        ["Visit ID", "방문 번호. 그 지점에 몇 번째 가는지입니다(V1, V2 …). "
                     "현장에서 직접 적으십니다."],
        ["Location ID", "표석 자리 번호. 중심점을 옮기면 새 번호가 됩니다."],
        ["Set ID", "한 방향 또는 수직 한 벌의 묶음 번호(S01-H동, S01-V …)."],
        ["차수", "그 Set 을 몇 번째로 재는지. 처음이 1 이고 다시 잴 때마다 "
                "1씩 올립니다."],
        ["정·반 시준", "같은 표지를 망원경을 뒤집어 두 번 보는 것입니다."]],
       [18, 82]))

    A(("h2", "시각은 이렇게 적습니다"))
    A(("callout",
       "시각은 «한국 표준시(KST)» 로, 날짜와 시·분·초까지 적습니다. "
       "예: 2026-10-14 13:24:05. 그리고 도착 시각 한 번이 아니라 "
       "«F 를 잰 모든 줄»에 그때의 시각을 적으셔야 합니다. P0 앞·뒤도 마찬가지고, "
       "수직 측정의 높이별 줄도 그렇습니다. 이것 하나가 빠지면 그 줄은 계산에서 "
       "빠집니다."))

    A(("h2", "이 카드가 원본입니다"))
    A(("p", "종이에 적은 것이 원본이고, 장비에 저장된 자료는 나중에 따로 불러옵니다. "
            "그래서 «장비 원시파일명»과 «레코드 범위»를 1 구역에 적어 주셔야 합니다. "
            "철수할 때 몰아서 적으면 어느 구간이었는지 기억이 흐려지니, 한 지점을 "
            "시작할 때 파일명을 적고 끝낼 때 레코드 범위를 적는 편이 낫습니다."))
    A(("callout",
       "종이 카드에는 자동으로 계산되는 칸이 없습니다. 「P0 전후차」나 「폐합차」 "
       "같은 초록 칸은 «사무실에서 채우는 자리»입니다. 현장에서 확인하고 싶으시면 "
       "두 값을 직접 빼 보시면 됩니다. (같은 카드를 엑셀로 여시면 그 칸은 자동으로 "
       "계산됩니다.)"))

    A(("pagebreak",))

    A(("h1", "1. 하루의 흐름"))
    A(("p", "한 지점에서 하는 일을 순서대로 적으면 이렇습니다. 옆에 적은 시간은 "
            "«아직 재어 보지 않은 가늠치»입니다. 접근·설치·기상 대기·재측정은 "
            "빼고 센 값이라 실제로는 더 걸릴 수 있습니다. 이번 탐사로 실제 시간을 "
            "재어 볼 참이니, 시간에 쫓겨 기록을 빠뜨리는 편이 훨씬 손해입니다."))
    A(("table", ["순서", "무엇을", "카드 구역", "대략 걸리는 시간"],
       [["1", "몸에 지닌 자기성 물품 정리, 차량·전자기기 확인", "0", "5분"],
        ["2", "도착 기록, Kp·예보 확인, 장비 준비", "1", "10분"],
        ["3", "GNSS 로 기준점 좌표 실측", "2", "10분"],
        ["4", "수평 자기구배 — 동·서·남·북 네 방향", "3", "1시간 남짓"],
        ["5", "수직 자기구배 — 중심점에서 높이별", "4", "20분"],
        ["6", "방위표지 시준", "5", "20분"],
        ["7", "사진 촬영", "6", "15분"],
        ["8", "중심점 결정, 소견 작성, 완전성 점검", "7", "10분"]],
       [8, 44, 12, 20]))
    A(("p", f"수평은 한 방향에 {len(SP.H_OFFSETS_M)}개 지점({offs} m)을 재고, 방향이 "
            f"바뀔 때마다 중심점으로 돌아와 한 번 더 잽니다. 네 방향이니 {nH}번, "
            f"수직까지 하면 한 지점에서 {nH + nV}번쯤 재게 됩니다. 한 자리에서 한 "
            f"번씩만 재시면 됩니다 — 기기가 안에서 여러 번 재어 평균을 내 줍니다."))

    A(("h1", "2. 이것만은 꼭 — 빠지면 되살릴 수 없는 것"))
    A(("p", "아래 항목들은 나중에 사무실에서 어떤 방법으로도 복원할 수 없습니다. "
            "실제로 예전 야장 68개를 열어 보니 총자력 측정시각이 한 건도 적혀 있지 "
            "않았고, 그 때문에 시간에 따른 보정을 아예 할 수 없게 된 일이 있었습니다."))
    A(("table", ["항목", "어디에", "왜 없으면 안 되는가"],
       [["측정 시각 (시·분·초)", "3·4 구역",
         "이 시각으로 재는 동안 자기장이 흘러간 만큼을 빼냅니다. 분까지만 "
         "적으면 계산이 되지 않습니다."],
        ["P0 앞·뒤 두 번", "3 구역",
         "한 방향을 다 재고 중심점으로 돌아와 한 번 더 재야 그 사이를 이을 수 "
         "있습니다. 한쪽이 비면 그 방향 전체가 계산되지 않습니다."],
        ["차수", "3·4 구역",
         "다시 잰 것인지 처음 잰 것인지를 가르는 유일한 표시입니다. 안 적으면 "
         "1차와 2차가 섞여 둘 다 못 쓰게 됩니다."],
        ["장비 원시파일명·레코드 범위", "1 구역",
         "종이 기록과 장비에 저장된 원본을 잇는 유일한 끈입니다."],
        ["실제 센서 높이", "4 구역",
         "목표 높이가 아니라 센서 가운데가 실제로 몇 cm 였는지를 적어야 나중에 "
         "같은 조건으로 다시 재어 볼 수 있습니다."],
        ["표지 좌표와 취득방법", "5 구역",
         "폐합차가 작아도 참방위각이 틀릴 수 있습니다. 같은 점을 다시 찾아가 "
         "재었을 때 편각이 33.7분이나 어긋난 일이 있었는데, 그 원인이 방문할 "
         "때마다 달라진 참방위각이었습니다."]],
       [22, 12, 66]))

    A(("pagebreak",))

    A(("h1", "3. 구역별로 무엇을 적는가"))

    A(("h2", "0 구역 — 측정 전 확인"))
    A(("p", "측정을 시작하기 전에 몸에 지닌 자기성 물품부터 걷어내 주세요. 시계나 펜, "
            "금속 혁대 같은 것들입니다. 센서 가까이 쇠붙이가 있으면 그 영향이 값에 "
            "섞여 들어오는데, 나중에 보면 그것이 사람 때문인지 그 자리의 성질인지 "
            "가려낼 방법이 없습니다."))
    A(("p", "차량은 측정 범위(중심점에서 10 m) 밖으로 충분히 물려 두시고, 무전기나 "
            "노트북 같은 전자기기는 꺼 주세요. 얼마나 떨어뜨려야 하는지는 아직 "
            "정해지지 않았으니, 여유 있게 두시고 어느 정도 떨어뜨렸는지 적어 "
            "주시면 됩니다."))
    A(("note", "사정상 빼지 못한 것이 있다면 무엇을 왜 못 뺐는지 적어 주시면 됩니다. "
               "적어 두시면 나중에 그 방향 값을 해석할 때 참고할 수 있습니다."))

    A(("h2", "1 구역 — 도착 · 관측 조건 · 원시자료 연결"))
    A(("p", "관측일자와 도착 시각, 관측자, 기상을 적습니다. Kp 지수와 우주기상 예보도 "
            "확인해서 적어 주시는데, 값이 크게 나왔다고 해서 그 자료를 빼지는 마세요. "
            "몇 이상이면 걸러야 하는지가 아직 정해지지 않았고, 사실 그 기준을 만들려고 "
            "이 탐사를 하는 것이기 때문입니다."))
    A(("p", "야간 관측은 하지 않습니다."))
    A(("callout",
       "다시 재는 것은 «값이 커서»가 아니라 «재는 절차가 깨져서» 하는 것입니다. "
       "예를 들어 한 방향을 재는 동안 자기장이 크게 흔들렸다면 시간 보정이 믿을 만하지 "
       "않으니 다시 잽니다. 그러나 값 자체가 크게 나온 것은 다시 잴 이유가 되지 "
       "않습니다. 한 방향을 시작했으면 값이 어떻게 나오든 «끝까지» 재 주세요. "
       "큰 값만 골라 다시 재거나 빼면 나중에 기준을 정할 때 자료가 한쪽으로 "
       "치우칩니다."))
    A(("callout",
       "장비 원시파일명과 레코드 범위를 꼭 적어 주세요. 장비에 저장된 원본은 나중에 "
       "따로 불러올 예정인데, 그때 어느 파일의 어느 구간이 이 카드에 해당하는지 "
       "알아낼 방법은 여기 적힌 정보밖에 없습니다."))

    A(("h2", "2 구역 — 기준점 좌표"))
    A(("p", "GNSS 로 실측한 위도·경도를 «십진도»로 적습니다. 위도를 먼저, 경도를 "
            "나중에 쓰시고 소수점 여섯 자리까지 적어 주세요(예: 37.020705 / "
            "128.983930). 장비 화면이 도·분·초로 나온다면 십진도로 바꿔 주는 설정이 "
            "있는지 먼저 확인해 주시고, 어려우면 화면 그대로 적으시고 그 사실을 "
            "옆에 적어 주세요. 예전에 좌표 형식이 세 가지로 섞여 들어와 애를 먹은 "
            "적이 있습니다."))
    A(("p", "사전좌표와 얼마나 차이가 나는지는 사무실에서 계산합니다. 많이 옮기셨다면 "
            "어떤 사정이었는지 소견에 적어 주시면 됩니다. 도상에서 고른 자리라 실제로 "
            "가 보면 진입이 어렵거나 나무가 우거진 경우가 있습니다."))

    A(("h2", "3 구역 — 수평 자기구배"))
    A(("p", "가장 시간이 많이 걸리는 부분입니다. 순서는 이렇습니다."))
    A(("steps", [
        "중심점(P0)에서 한 번 재고 시각과 값을 적습니다.",
        f"동쪽으로 {SP.H_OFFSETS_M[0]} m, {SP.H_OFFSETS_M[1]} m, 그다음부터는 "
        f"1 m 씩 {SP.H_OFFSETS_M[-1]} m 까지 옮겨 가며 잽니다.",
        "다시 중심점으로 돌아와 한 번 더 잽니다.",
        "서·남·북 방향도 같은 방식으로 반복합니다."]))
    A(("p", "중심점을 앞뒤로 두 번 재는 이유는, 그 사이에 자기장이 저절로 변한 만큼을 "
            "빼내기 위해서입니다. 자력계가 한 대뿐이라 위치를 옮기는 동안 시간도 함께 "
            "흘러가는데, 앞뒤 두 값을 직선으로 이으면 각 측점을 잰 순간의 중심점 값을 "
            "추정할 수 있습니다. 이것을 빼지 않으면 시간이 흘러간 몫까지 그 자리의 "
            "값인 것처럼 계산됩니다."))
    A(("p", "한 방향을 마치면 P0 앞뒤 두 값을 직접 빼 보시면 좋습니다. 이 차이가 수십 "
            "nT 로 벌어졌다면 재는 동안 자기장이 꽤 흔들렸다는 뜻이라 시간 보정이 "
            "믿을 만하지 않습니다. 이럴 때는 차수를 올려 그 방향을 다시 재 주시고, "
            "여건이 안 되면 소견에 적어 주세요. 다시 재실 때도 앞서 적은 값은 "
            "지우지 마시고 그대로 두시면 됩니다."))

    A(("h2", "4 구역 — 수직 자기구배"))
    A(("p", f"중심점에서 센서 높이만 바꿔 가며 잽니다. 지상 {SP.V_HEIGHTS_CM[0]} cm "
            f"부터 {SP.V_HEIGHTS_CM[-1]} cm 까지 {step} cm 간격이고, 시작할 때와 끝날 "
            f"때 기준높이에서 한 번씩 더 재 주시면 됩니다. 수평과 같은 이유입니다."))
    A(("p", "여기서 「기준높이」는 앞뒤로 두 번 재는 그 높이를 말합니다. 어느 높이로 "
            "할지는 정해져 있지 않으니 재기 편한 높이 하나를 고르시고(보통 가운데쯤인 "
            "100 cm 를 씁니다), 카드 4 구역 머리의 「기준높이」 칸에 그 값을 적어 "
            "주세요. 한 지점 안에서는 같은 높이를 쓰셔야 합니다."))
    A(("p", "높이는 목표값이 아니라 센서 가운데가 실제로 몇 cm 였는지를 적어 주세요. "
            "그리고 어디를 0 으로 보았는지(지면인지 표석 윗면인지)도 함께 적어 주셔야 "
            "나중에 같은 조건을 재현할 수 있습니다."))

    A(("h2", "5 구역 — 방위표지 시준"))
    A(("p", "기준점에서 방위표지를 바라본 방향을 잽니다. 진북을 기준으로 북쪽이 0도, "
            "시계 방향입니다. 방향은 언제나 기준점에서 표지 쪽이며, 반대로 적으면 그 "
            "방문의 편각이 통째로 틀어집니다."))
    A(("p", "각도는 곤(gon)으로 적습니다. 한 바퀴가 400 이고 반대쪽이 200 입니다. "
            "정·반 시준의 폐합차는 카드가 계산해 주는데, 0.02 gon(약 65초)을 넘으면 "
            "다시 시준해 주시는 편이 좋습니다."))
    A(("p", "카드 5 구역 윗줄에 표지마다 «좌표»와 «그 좌표를 어떻게 얻었는지»를 적는 "
            "칸이 있습니다. GNSS 로 재셨으면 어떤 방식이었는지(RTK 인지 휴대용인지)와 "
            "정확도도 함께 적어 주세요. 그리고 예전에 쓰던 참방위각이 있으면 그 값도 "
            "옮겨 적어 주시면, 사무실에서 이번 값과 견주어 볼 수 있습니다."))
    A(("callout",
       "폐합차가 작다고 해서 참방위각이 맞다는 뜻은 아닙니다. 같은 표지를 쓰는지가 "
       "더 중요합니다. 재방문이라면 예전과 «같은 표지»를 시준하셨는지 확인하시고, "
       "다른 표지를 쓰셨다면 그 사실과 이유를 반드시 적어 주세요. 예전에 이것이 "
       "기록되지 않아, 방문마다 기준이 달라졌는데도 같은 자료인 것처럼 다뤄진 일이 "
       "있었습니다. 사진 칸에 표지 사진을 남겨 주시면 확인에 큰 도움이 됩니다."))

    A(("h2", "6 구역 — 현장 사진"))
    A(("p", "아홉 칸이 있습니다. 중심점 전경, 측정 모습, 방위표지 1·2, 동·서·남·북 "
            "네 방향 전경, 그리고 자기교란 요소입니다. 방위별 전경은 중심점에 서서 그 "
            "방향을 보고, 측선과 10 m 끝점이 함께 담기도록 찍어 주세요. 나중에 어느 "
            "방향 값이 크게 나왔을 때 그 사진에서 원인을 찾게 됩니다."))
    A(("p", "사진마다 원본 파일명과 촬영시각을 적어 주세요. 촬영시각까지 적는 이유는, "
            "카메라가 여럿이면 파일명이 겹치거나 옮기는 과정에서 이름이 바뀌는 일이 "
            "있어 파일명만으로는 어느 사진인지 되짚기 어렵기 때문입니다."))

    A(("h2", "7 구역 — 중심점 결정과 소견"))
    A(("p", "후보지 안에서 더 나은 자리를 찾는 것도 이 조사의 목적 가운데 하나입니다. "
            "한 방향이 유난히 크게 나왔다면 옮기는 것을 검토해 보셔도 됩니다. 다만 "
            "«처음 자리에서 잰 것을 먼저 끝내신 뒤»에 옮기세요. 중간에 그만두고 "
            "옮기면 그 자리가 어땠는지가 남지 않습니다."))
    A(("p", "옮기실 때는 이 카드를 고쳐 쓰지 마시고 새 Location ID 로 카드를 한 장 더 "
            "쓰시면 됩니다. 새 카드는 좌표부터 사진까지 «위치에 딸린 구역을 모두» 다시 "
            "채우셔야 합니다 — 2·3·4·5·6 구역이 그렇습니다. 자리가 바뀌면 방위표지와의 "
            "관계도, 주변 사진도 달라지기 때문입니다. 앞 카드는 그대로 두시고, 새 "
            "카드 소견에 「S01 카드에서 옮김」처럼 관계만 적어 주세요."))
    A(("p", "마지막 소견 칸에는 현장에서 느낀 것을 자유롭게 적어 주세요. 값만으로는 "
            "알 수 없는 것들 — 주변에 무엇이 있었는지, 접근이 어땠는지, 다시 온다면 "
            "무엇을 준비해야 하는지 — 이 나중에 큰 도움이 됩니다."))

    A(("pagebreak",))

    A(("h1", "4. 이럴 때는 어떻게 하나"))
    A(("table", ["상황", "이렇게 하시면 됩니다"],
       [["재는 도중 값이 이상하다",
         "그 방향을 처음부터 다시 재시고 차수를 2 로 올려 주세요. 앞서 적은 값은 "
         "지우지 마시고 그대로 두시면 됩니다. 사무실에서 차수로 갈라 보기 때문에 "
         "1차와 2차가 섞이지 않습니다."],
        ["P0 전후차가 크게 벌어졌다",
         "재는 동안 자기장이 흔들렸다는 뜻입니다. 여유가 되면 차수를 올려 다시 "
         "재시고, 어려우면 소견에 적어 주세요."],
        ["중간에 중단했다",
         "중단한 시각과 사유를 적어 두시고, 다시 시작하실 때 차수를 올려 주세요."],
        ["값이 참고값보다 크게 나온다",
         f"그대로 적어 주세요. 반경 10 m 안에서 {SP.IAGA_RANGE_NT:.0f} nT, 구배 "
         f"{SP.EURO_GRAD_NT_PER_M:.0f} nT/m 라는 국제 권고가 있긴 하지만 국내 "
         f"기준은 아직 없습니다. 현장에서 그 자리를 접지 마세요."],
        ["방위표지를 바꿔야 했다",
         "바꾸신 것 자체는 문제가 되지 않습니다. 다만 새 표지의 좌표와 예전 값, "
         "바꾼 이유를 반드시 적어 주세요."],
        ["기준점이 사전좌표와 많이 다르다",
         "실측한 좌표를 적으시고 왜 옮기셨는지 소견에 남겨 주세요. 도상에서 고른 "
         "자리라 진입이 어려운 경우가 있습니다."],
        ["Variometer 를 못 돌렸다",
         "「미운용」으로 표시해 주세요. 값을 못 쓰는 것은 아니지만, 시간에 따른 "
         "변화를 검증할 수 없다는 표시가 남습니다."]],
       [24, 76]))

    A(("h2", "빈칸을 남기지 않는 법"))
    A(("p", "종이에서는 빈칸이 「깜빡한 것」인지 「해당이 없는 것」인지 구분되지 "
            "않습니다. 그래서 값을 못 적을 때는 그냥 비워 두지 마시고 아래처럼 "
            "적어 주세요. 한 글자라도 적혀 있으면 사무실에서 판단할 수 있습니다."))
    A(("table", ["이런 경우", "이렇게 적습니다"],
       [["장애물이 있어 그 지점을 못 쟀다", "측정불가 + 사유(예: 「측정불가 — 바위」)"],
        ["그 지점에는 해당하지 않는다", "해당없음"],
        ["장비 오류로 값이 못 미덥다", "무효 + 사유. 값은 지우지 말고 그대로 두세요"],
        ["방위표지가 하나뿐이다", "표지2 줄에 「해당없음」"],
        ["사진을 못 찍었다", "그 칸에 「미촬영 — 사유」"]],
       [34, 66]))
    A(("callout",
       "잘못 적으셨을 때는 지우거나 덧쓰지 마시고, «가로줄 하나로 긋고 옆에 다시» "
       "적어 주세요. 그리고 왜 고쳤는지 짧게 남겨 주시면 됩니다. 원래 무엇이라고 "
       "적혀 있었는지가 남아 있어야 나중에 되짚을 수 있습니다."))

    A(("h1", "5. 철수하기 전에"))
    A(("p", "카드 3 구역 아래에 「완전성 점검」 줄이 있습니다. 다섯 가지를 스스로 "
            "확인하시고 「예 / 아니오」를 적어 주세요. 하나라도 「아니오」가 있으면 "
            "현장을 떠나기 전에 채우시는 편이 좋습니다. 돌아온 뒤에는 채울 방법이 "
            "없습니다."))
    A(("table", ["점검 항목", "확인할 것"],
       [["필수행 다 채움", "빈 칸으로 남은 측점이 없는지"],
        ["시각 시:분:초", "분까지만 적은 곳이 없는지"],
        ["P0 전·후 있음", "네 방향 모두 앞뒤로 중심점을 쟀는지"],
        ["측정시각 P0 사이", "측점 시각이 P0 두 시각 사이에 들어오는지"],
        ["Variometer 연결", "파일명을 적었는지 (미운용이면 그렇게 표시)"]],
       [28, 72]))
    A(("p", "그리고 사진 아홉 칸의 파일명과 촬영시각, 장비 원시파일명과 레코드 범위를 "
            "적으셨는지 한 번 더 봐 주세요. 이 둘이 빠지면 나중에 자료를 이을 수 "
            "없습니다."))

    A(("h1", "6. 안전"))
    A(("p", "산지 현장이고 혼자 다니실 때도 있습니다. 자세한 절차는 기관의 안전 "
            "지침을 따르시되, 이 일에 한해서는 아래 세 가지만 먼저 기억해 주세요."))
    A(("table", ["", "내용"],
       [["작업 중지가 먼저", "날씨가 나빠지거나 몸이 좋지 않으시면 측정을 "
                          "멈추십시오. 한 지점을 못 끝내는 것보다 다치는 것이 "
                          "훨씬 큰 손실입니다. 못 끝낸 이유를 소견에 적어 "
                          "주시면 됩니다."],
        ["출발 전 알리기", "어디로 가서 언제 돌아올지 사무실에 알리고 나가 "
                        "주세요. 통신이 안 되는 곳이 많습니다."],
        ["출입 확인", "사유지나 임도는 미리 허락을 받아야 하는 곳이 있습니다. "
                    "현장에서 제지를 받으면 다투지 마시고 물러난 뒤 알려 "
                    "주세요."]],
       [22, 78]))
    A(("note", "비상연락처와 적용 안전지침 문서는 배포 시 이 자리에 채워 넣습니다."))

    A(("h1", "7. 참고 — 이번 탐사의 성격"))
    A(("p", "야장을 채우다 보면 「이건 왜 이렇게 하지」 싶은 곳이 있을 텐데, 대개 "
            "아래 네 가지 결정에서 나옵니다."))
    # ⚠️ 「정해진 것」은 SP.DECISIONS 원문 그대로 — 갈라지지 않게. 「현장에서는」만
    #    이 안내서의 말투로 따로 쓴다. 원문을 마침표로 자르면 「2.6 h」가 끊긴다.
    FIELD_MEANING = {
        "법정 성격": "값이 크게 나와도 현장에서 그 자리를 접지 마시고 그대로 적어 "
                   "주세요. 합격선이 아직 없습니다.",
        "관측 조건": "밤에 나가실 필요는 없습니다. Kp 는 적어 두기만 하시고 그 값으로 "
                   "자료를 빼지는 마세요.",
        "현장 원본 매체": "이 카드가 원본입니다. 장비 원시파일명과 레코드 범위를 적어 "
                       "두셔야 나중에 기기 자료와 이을 수 있습니다.",
        "작업량": "한 지점에 몇 시간이 걸리는지는 이번 탐사로 재어 보는 중입니다. "
                "서두르시다 기록을 빠뜨리는 편이 더 손해입니다.",
    }
    A(("table", ["항목", "정해진 것", "현장에서는"],
       [[k, v.replace("**", ""), FIELD_MEANING.get(k, "")]
        for k, v, _ in SP.DECISIONS],
       [18, 34, 48]))
    A(("note", "측정 설계와 참고값은 IAGA 「Guide for Magnetic Repeat Station "
               "Surveys」(1996) 제4.2절과 유럽 반복관측점 권고를 따랐습니다. "
               "자세한 근거는 중간보고 자료를 보시면 됩니다."))
    return B


# ══════════════════════════════════════════════════════════════ docx
def render_docx(blocks, path):
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt, RGBColor

    NAVY = RGBColor(0x1F, 0x38, 0x64)
    RED = RGBColor(0xB2, 0x18, 0x2B)
    GREY = RGBColor(0x55, 0x55, 0x55)
    WIDE = 15.8                       # 판면 폭(cm)

    def shade(cell, hexc):
        el = OxmlElement("w:shd")
        el.set(qn("w:fill"), hexc)
        cell._tc.get_or_add_tcPr().append(el)

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "맑은 고딕"
    st._element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")
    st.font.size = Pt(10.5)
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.35
    for lv, size in ((0, 16), (1, 13), (2, 11.5)):
        hs = doc.styles[f"Heading {lv + 1}"]
        hs.font.name = "맑은 고딕"
        hs._element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")
        hs.font.size = Pt(size)
        hs.font.color.rgb = NAVY
    for s in doc.sections:
        s.left_margin = s.right_margin = Cm(2.0)
        s.top_margin = s.bottom_margin = Cm(1.8)

    def para(text, size=10.5, color=None, indent=0):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(6)
        if indent:
            p.paragraph_format.left_indent = Cm(indent)
        run = p.add_run(text)
        run.font.size = Pt(size)
        if color is not None:
            run.font.color.rgb = color

    for b in blocks:
        kind = b[0]
        if kind == "cover":
            for txt, size, col, bold in ((b[1], 22, NAVY, True),
                                         (b[2], 15, GREY, False),
                                         (b[3], 10, GREY, False)):
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                r = p.add_run(txt)
                r.font.size = Pt(size)
                r.font.bold = bold
                r.font.color.rgb = col
            doc.add_paragraph()
        elif kind == "h1":
            doc.add_heading(b[1], level=1)
        elif kind == "h2":
            doc.add_heading(b[1], level=2)
        elif kind == "p":
            para(b[1])
        elif kind == "note":
            para(b[1], size=10, color=GREY)
        elif kind == "steps":
            for i, t in enumerate(b[1], 1):
                para(f"{'①②③④⑤⑥⑦⑧⑨'[i-1]} {t}", indent=0.5)
        elif kind == "callout":
            t = doc.add_table(rows=1, cols=1)
            t.style = "Table Grid"
            c = t.rows[0].cells[0]
            c.text = ""
            r = c.paragraphs[0].add_run(b[1])
            r.font.size = Pt(9.5)
            r.font.color.rgb = RED
            shade(c, CALLOUT_BG)
            doc.add_paragraph().paragraph_format.space_after = Pt(4)
        elif kind == "table":
            heads, rows, pct = b[1], b[2], b[3]
            t = doc.add_table(rows=1, cols=len(heads))
            t.style = "Table Grid"
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for i, h in enumerate(heads):
                c = t.rows[0].cells[i]
                c.text = ""
                r = c.paragraphs[0].add_run(h)
                r.font.bold = True
                r.font.size = Pt(9.5)
                r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
                c.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
                shade(c, NAVY_HEX)
            for row in rows:
                cells = t.add_row().cells
                for i, v in enumerate(row):
                    cells[i].text = ""
                    r = cells[i].paragraphs[0].add_run(str(v))
                    r.font.size = Pt(9.5)
                    if i == 0:
                        r.font.bold = True
            for rr in t.rows:
                for i, w in enumerate(pct):
                    rr.cells[i].width = Cm(WIDE * w / 100)
            doc.add_paragraph().paragraph_format.space_after = Pt(4)
        elif kind == "pagebreak":
            doc.add_page_break()
    doc.save(path)
    return path


# ══════════════════════════════════════════════════════════════ hwpx
def _cp949(t):
    """CP949 에 없는 글자를 바꾼다 — 안 바꾸면 한글이 못 읽는다."""
    return t.replace("—", "―").replace("−", "-")


def render_html(blocks, path):
    """한컴에 넘길 중간 HTML. **CP949 로 저장**한다."""
    e = lambda t: _html.escape(_cp949(str(t)))
    out = ['<html><head><meta http-equiv="Content-Type" '
           'content="text/html; charset=euc-kr"><style>',
           'body{font-family:"맑은 고딕";font-size:10.5pt;line-height:1.5}',
           f'h1{{font-size:13pt;color:#{NAVY_HEX};margin:16pt 0 6pt}}',
           f'h2{{font-size:11.5pt;color:#{NAVY_HEX};margin:12pt 0 4pt}}',
           'p{margin:0 0 6pt}',
           'table{border-collapse:collapse;width:100%;margin:4pt 0 10pt}',
           'td,th{border:1px solid #999;padding:3pt 5pt;font-size:9.5pt;'
           'vertical-align:top}',
           f'th{{background:#{NAVY_HEX};color:#ffffff;text-align:center}}',
           '</style></head><body>']
    for b in blocks:
        k = b[0]
        if k == "cover":
            out.append(f'<p style="text-align:center;font-size:22pt;font-weight:bold;'
                       f'color:#{NAVY_HEX}">{e(b[1])}</p>')
            out.append(f'<p style="text-align:center;font-size:15pt;'
                       f'color:#{GREY_HEX}">{e(b[2])}</p>')
            out.append(f'<p style="text-align:center;font-size:10pt;'
                       f'color:#{GREY_HEX}">{e(b[3])}</p><p>&nbsp;</p>')
        elif k == "h1":
            out.append(f"<h1>{e(b[1])}</h1>")
        elif k == "h2":
            out.append(f"<h2>{e(b[1])}</h2>")
        elif k == "p":
            out.append(f"<p>{e(b[1])}</p>")
        elif k == "note":
            out.append(f'<p style="font-size:10pt;color:#{GREY_HEX}">{e(b[1])}</p>')
        elif k == "steps":
            for i, t in enumerate(b[1], 1):
                out.append(f'<p style="margin-left:14pt">'
                           f'{"①②③④⑤⑥⑦⑧⑨"[i-1]} {e(t)}</p>')
        elif k == "callout":
            out.append(f'<table><tr><td style="background:#{CALLOUT_BG};'
                       f'color:#{RED_HEX}">{e(b[1])}</td></tr></table>')
        elif k == "table":
            heads, rows, pct = b[1], b[2], b[3]
            cols = "".join(f'<col style="width:{w}%">' for w in pct)
            th = "".join(f"<th>{e(h)}</th>" for h in heads)
            body = "".join(
                "<tr>" + "".join(
                    f'<td{" style=font-weight:bold" if i == 0 else ""}>{e(v)}</td>'
                    for i, v in enumerate(r)) + "</tr>" for r in rows)
            out.append(f"<table>{cols}<tr>{th}</tr>{body}</table>")
        elif k == "pagebreak":
            out.append('<p style="page-break-before:always">&nbsp;</p>')
    out.append("</body></html>")
    Path(path).write_text("\n".join(out), encoding="cp949", errors="replace")
    return path


def render_hwpx(blocks, path):
    """HTML 을 한글로 열어 hwpx 로 저장한다."""
    import win32com.client as win32

    tmp_html = Path(os.environ["TEMP"]) / "_field_guide.html"
    render_html(blocks, tmp_html)
    if os.path.exists(path):
        os.remove(path)
    hwp = win32.gencache.EnsureDispatch("HWPFrame.HwpObject")
    try:
        hwp.RegisterModule("FilePathCheckDLL", "FilePathCheckerModule")
    except Exception:
        pass
    if not hwp.Open(str(tmp_html), "HTML", ""):
        hwp.Quit()
        raise RuntimeError(f"한글이 HTML 을 열지 못했다: {tmp_html}")
    ok = hwp.SaveAs(str(path), "HWPX", "")
    text = hwp.GetTextFile("TEXT", "")
    hwp.Quit()
    if not ok or not os.path.exists(path):
        raise RuntimeError("hwpx 저장 실패")
    return path, len(text)


# ══════════════════════════════════════════════════════════════
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    want = set(a for a in sys.argv[1:] if a.startswith("--"))
    both = not want
    blocks = content()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stamp = f"{dt.datetime.now():%Y%m%d_%H%M%S}"
    base = "시험탐사_현장야장_기입안내"

    if both or "--docx" in want:
        p = render_docx(blocks, OUT_DIR / f"{stamp}_{base}.docx")
        print(f"[저장] {p.name}  ({p.stat().st_size/1000:.0f} KB)")
    if both or "--hwpx" in want:
        p = OUT_DIR / f"{stamp}_{base}.hwpx"
        p, n = render_hwpx(blocks, p)
        print(f"[저장] {Path(p).name}  ({Path(p).stat().st_size/1000:.0f} KB) "
              f"· 본문 {n:,}자")
    return 0


if __name__ == "__main__":
    sys.exit(main())
