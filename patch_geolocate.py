# -*- coding: utf-8 -*-
"""
「📍 내 위치」 버튼을 이미 만들어진 survey_review.html 에 끼운다.

현장에서 휴대폰으로 지도를 열었을 때 자기 자리를 보려면 버튼이 필요한데,
`docs/survey_review.html` 은 Folium 이 통째로 찍어 낸 파일이라 생성기를
다시 돌리지 않으면 바뀌지 않는다. 생성기(`make_survey_map.py`)를 돌리려면
현장조사 회신본이나 취합본이 있어야 하므로, 그 전까지 쓰는 경로다.

⚠️ **버튼의 HTML·JS 는 여기에 적지 않는다.** `make_survey_map.geolocate_ui()`
를 그대로 불러 쓰므로 생성기와 배포본이 갈라지지 않는다. 생성기를 다시
돌리면 같은 버튼이 처음부터 들어간다.

**멱등**이다 — 이미 들어 있으면 손대지 않는다.

    python patch_geolocate.py [--check]
"""
import sys
from pathlib import Path

ROOT = Path(__file__).parent
HTML = ROOT / "docs" / "survey_review.html"

OPEN, CLOSE = "<!-- GEO_UI -->", "<!-- /GEO_UI -->"


def patch(check_only=False):
    """버튼 블록을 넣거나, 이미 있으면 **최신 내용으로 갈아 끼운다**.

    ⚠️ 「이미 있으면 건너뛴다」로 두면 안내 문구를 고쳐도 배포본이 옛 문구를
    계속 쓴다(실제로 아이폰 안내가 틀린 채 남아 있었다). 그래서 마커 사이를
    통째로 바꾼다 — `build_sv_viewer.py` 의 `SV_CARD` 와 같은 방식이다.
    """
    if not HTML.exists():
        print(f"[건너뜀] {HTML.name} 없음")
        return 1
    s = HTML.read_text(encoding="utf-8")

    from make_survey_map import geolocate_ui
    block = geolocate_ui()

    i, j = s.find(OPEN), s.find(CLOSE)
    if i != -1 and j != -1:
        cur = s[i:j + len(CLOSE)]
        if cur == block.strip():
            print("■ 내 위치 버튼 — 이미 최신")
            return 0
        s = s[:i] + block.strip() + s[j + len(CLOSE):]
        print(f"■ 내 위치 버튼 — 갱신 ({len(block):,} 바이트)")
    elif "id='geoWrap'" in s:
        # 마커를 붙이기 «전»에 끼운 판본 — 구조로 찾아 갈아 끼운다.
        # 블록은 <style>…</style> + <div id='geoWrap'>…</div> + <script>…</script>
        # 하나뿐이라 style 시작 뒤 첫 </script> 가 블록의 끝이다.
        a = s.find("<style>\n#geoWrap{")
        b = s.find("</script>", a)
        if a == -1 or b == -1:
            print("■ 내 위치 버튼 — 옛 판본을 찾지 못함(손으로 지울 것)")
            return 1
        s = s[:a] + block.strip() + s[b + len("</script>"):]
        print(f"■ 내 위치 버튼 — 옛 판본 교체 ({len(block):,} 바이트)")
    else:
        # Folium 산출물에는 </body> 가 없을 수 있다 — 그때는 </html> 앞에 넣는다.
        for tag in ("</body>", "</html>"):
            k = s.rfind(tag)
            if k != -1:
                s = s[:k] + block + "\n" + s[k:]
                break
        else:
            s = s + block
        print(f"■ 내 위치 버튼 — 추가 ({len(block):,} 바이트)")
    if not check_only:
        HTML.write_text(s, encoding="utf-8")
        print(f"    [저장] {HTML.name}")
    return 1


def main():
    rc = patch("--check" in sys.argv)
    if "--check" in sys.argv:
        print(f"\n[점검] 추가 필요: {rc}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
