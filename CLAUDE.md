# 지구자기장 측정 입지 선정 프로젝트 — Claude 작업 규칙

## 출력 파일 저장 경로

**모든 작업 산출물은 `docs/output/` 폴더에 저장한다.**

- Word 문서(.docx), CSV, 보고서 등 생성되는 파일은 반드시 `docs/output/` 경로 사용
- `output/` (루트 하위) 에 저장하지 않는다 — GitHub Pages 배포 폴더인 `docs/` 기준이므로 커밋 대상과 일치시켜야 함
- 스크립트 내 절대경로 하드코딩 금지. 대신 `Path(__file__).parent / "docs" / "output"` 패턴 사용

## 파일명 규칙

- 날짜/시간이 포함된 파일을 저장할 때는 반드시 `YYYYMMDD_HHMMSS` 접두사 사용
  - 예: `20260514_102055_site_selection_methodology_geology.docx`

## Git 커밋 규칙

- 커밋 대상: `docs/` 폴더 전체 (HTML, GeoJSON, output 산출물 포함)
- `output/` 루트 폴더는 커밋하지 않는다
- push는 `git push origin HEAD:main` 사용 (브랜치명 불일치 방지)

## 주요 스크립트

| 스크립트 | 역할 |
|---|---|
| `geomag_site_selection.py` | 메인 분석 + Folium 지도 생성 |
| `build_docs.py` | `output/` → `docs/` 빌드 (GeoJSON 단순화 포함) |
| `create_methodology_doc.py` | 기술 방법론 Word 문서 생성 → `docs/output/` 저장 |

## Python 환경

- 인터프리터: `C:\Users\YOONS\anaconda3\python.exe`
- 실행 예시: `& "C:\Users\YOONS\anaconda3\python.exe" create_methodology_doc.py`
