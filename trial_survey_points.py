# -*- coding: utf-8 -*-
"""
시험 탐사 대상 50점 마스터 — 단일 출처
======================================

`export_selection_50.py` 가 내는 세 갈래(등급 A 5 · B 도엽대표 29 · 기존 관측망
선점대상 16)를 **야장이 쓰기 좋은 한 벌의 레코드**로 만든다. 야장 생성기가
좌표·방위표지·예측 구배를 직접 계산하지 않게 하려는 것이다 — 계산이 두 곳에
있으면 조용히 갈라진다.

⚠️ **등급·명단은 여기서 다시 정하지 않는다.** `export_selection_50.gather_*` 를
그대로 호출하므로 지도·엑셀과 언제나 같은 50점이다.

레코드 한 건이 담는 것:

    구분 · 번호 · 지점명 · 도엽 · 기준점 좌표 · 표고 · 소재지
    방위표지 1·2 좌표 / 방위각 / 거리      ← 현장에서 재확인할 기준값
    예측 |∇ΔT| · s5                      ← 시험 탐사가 검증할 값
    최근접 상시관측소 · 거리               ← 외부장 보정 설계용
    기존점이면 직전 성과(편각·복각·총자력·최종관측)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent


def _f(v):
    # ⚠️ `AnomalyGradient.grad()` 는 0차원 ndarray 를 준다 — `float()` 에 바로
    #    넣으면 NumPy 2 에서 DeprecationWarning 이 뜨므로 `.item()` 을 먼저 쓴다.
    try:
        if isinstance(v, np.ndarray):
            v = v.item()
        x = float(v)
        return None if x != x else x
    except (TypeError, ValueError):
        return None


def _km(la1, lo1, la2, lo2):
    m = np.radians((la1 + la2) / 2)
    return float(np.hypot((lo2 - lo1) * 111.320 * np.cos(m),
                          (la2 - la1) * 110.574))


def nearest_observatory(lat, lon):
    """(관측소명, 거리 km) — 외부장 보정의 공간보간 설계에 쓴다.

    ⚠️ 이 거리가 곧 「균일 V 근사가 성립하는가」의 척도다. 이 프로젝트는
    76~229 km 에서 그 근사가 깨짐을 실측했다(CLAUDE.md 12절).
    """
    from lmm_external_multi import STATIONS
    best = min((_km(lat, lon, s["lat"], s["lon"]), s["name"])
               for s in STATIONS.values())
    return best[1], best[0]


def load_points(with_gradient: bool = True):
    """시험 탐사 50점 레코드 목록."""
    import export_selection_50 as X

    sv, merged = X.gather_survey()
    ex = X.gather_existing(merged)

    ag = None
    if with_gradient:
        import anomaly_gradient as AG
        ag = AG.AnomalyGradient()

    def grad_of(lat, lon):
        if ag is None:
            return None, None
        import anomaly_gradient as AG
        g = _f(ag.grad(lat, lon))
        return g, (_f(AG.representativeness_score(g)) if g is not None else None)

    out = []
    for d in sv:
        lat, lon = float(d["위도"]), float(d["경도"])
        det = d.get("방위표지상세") or {}
        g, s5 = grad_of(lat, lon)
        obs, obs_km = nearest_observatory(lat, lon)
        rec = {
            "구분": d["_구분"], "번호": d["관리번호"], "지점명": d["후보지명"],
            "관할본부": d["관할본부"], "도엽번호": d["도엽번호"], "도엽명": d["도엽명"],
            "위도": lat, "경도": lon, "표고": _f(d["표고"]),
            "소재지": d["지번주소"] or d["도로명주소"],
            "선정근거": d["_사유"], "검토결론": d["_결론"],
            "연계기존점": d.get("연계기존점", ""),
            "예측구배": g, "s5": s5,
            "최근접관측소": obs, "관측소거리": obs_km,
            "종전성과": None,
        }
        for i in (1, 2):
            m = det.get(f"표지{i}") or {}
            ll = d.get(f"표지{i}ll")
            rec[f"표지{i}_위도"] = _f(ll[1]) if ll else None
            rec[f"표지{i}_경도"] = _f(ll[0]) if ll else None
            rec[f"표지{i}_방위각"] = m.get("방위각", "")
            rec[f"표지{i}_거리"] = _f(m.get("거리"))
        out.append(rec)

    for r in ex:
        lat, lon = float(r["위도"]), float(r["경도"])
        g, s5 = grad_of(lat, lon)
        obs, obs_km = nearest_observatory(lat, lon)
        out.append({
            "구분": "기존점", "번호": r["지점코드"], "지점명": r["지점명"],
            "관할본부": "", "도엽번호": r["도엽번호"], "도엽명": r["지점명"],
            "위도": lat, "경도": lon, "표고": r["표고"], "소재지": r["소재지"],
            "선정근거": "기존 관측망 · 표석 지속성 확보",
            "검토결론": "재관측 — 방위표지 재확인 후 성과 갱신",
            "연계기존점": "",
            "예측구배": g, "s5": s5,
            "최근접관측소": obs, "관측소거리": obs_km,
            # 기존점은 직전 성과가 있으므로 시험 탐사 결과와 직접 대조할 수 있다.
            "종전성과": {"편각": r["편각"], "복각": r["복각"], "총자력": r["총자력"],
                        "최종관측": r["최종관측"], "관측장비": r["관측장비"]},
            "표지1_위도": None, "표지1_경도": None, "표지1_방위각": "", "표지1_거리": None,
            "표지2_위도": None, "표지2_경도": None, "표지2_방위각": "", "표지2_거리": None,
            "답사": r.get("답사 관리번호", ""),
        })
    return out


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    pts = load_points()
    from collections import Counter
    print(f"시험 탐사 대상 {len(pts)}점 — {dict(Counter(p['구분'] for p in pts))}")
    gap = [p["지점명"] for p in pts if p["예측구배"] is None]
    print(f"항공자력 자료공백 {len(gap)}점 — {', '.join(gap)}")
    print("  → 이 점들은 예측 구배가 없어 **현장 측정이 유일한 근거**다.")
    d = np.array([p["관측소거리"] for p in pts])
    print(f"최근접 상시관측소: 중앙 {np.median(d):.0f} km · 최대 {d.max():.0f} km")
