# -*- coding: utf-8 -*-
"""
관측소 영년변화 뷰어용 로컬 서버.

  · docs/ 정적 파일 서빙(lmm.html 등)
  · /api/minute?ym=YYYY-MM[&stn=GN,JJ,ICH,CYG]  → 해당 월 분(分) 단위 자료 JSON

lmm.html 의 '분 단위(월별)' 기능이 이 API 를 호출한다. 서버 없이(순수 오프라인)에는
일 중앙값 뷰어만 동작하고, 분 단위는 안내 메시지로 폴백한다.

    python sv_server.py           # http://127.0.0.1:8899/lmm.html
    python sv_server.py --port 9000

분 단위 원자료는 크므로(월·관측소당 ~43,200행) 웹에 내장하지 않고 요청 시에만 읽는다.
KASA 연도파일과 CYG 일자파일은 로드 후 메모리 캐시한다.
"""
import argparse
import calendar
import csv
import glob
import io
import json
import re
import datetime as dt
from functools import lru_cache
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).parent
DOCS = ROOT / "docs"
KA = ROOT / "data" / "kasa"
CY = ROOT / "data" / "cyg"

STN = {  # code -> (name, color)
    "CYG": ("청양", "#1f6feb"), "GN": ("강릉", "#6b7280"),
    "JJ": ("제주", "#0E8A6B"), "ICH": ("이천", "#c86a2b"),
}
BAND = {"X": (20000, 40000), "Y": (-9000, 3000), "Z": (28000, 46000), "F": (40000, 60000)}


def _clean(v, comp):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    lo, hi = BAND[comp]
    return round(x) if lo <= x <= hi else None


@lru_cache(maxsize=32)
def _kasa_year(code, year):
    """KASA 연도파일 → { 'YYYY-MM-DD HH:MM:00': (X,Y,Z,F) }."""
    p = KA / f"kasa_{code}_{year}.csv"
    out = {}
    if not p.exists():
        return out
    with open(p, encoding="utf-8-sig") as fh:
        rd = csv.reader(fh)
        next(rd, None)
        for r in rd:
            if len(r) < 6:
                continue
            out[r[1]] = (_clean(r[2], "X"), _clean(r[3], "Y"),
                         _clean(r[4], "Z"), _clean(r[5], "F"))
    return out


@lru_cache(maxsize=400)
def _cyg_day(day_iso):
    """CYG 하루 파일(UTC) → { 'YYYY-MM-DD HH:MM:00'(KST): (X,Y,Z,F) }."""
    ymd = day_iso.replace("-", "")
    matches = glob.glob(str(CY / "**" / f"cyg_{ymd}.csv"), recursive=True)
    out = {}
    if not matches:
        return out
    with open(matches[0], encoding="utf-8") as fh:
        rd = csv.reader(fh)
        hdr = next(rd, None)
        for r in rd:
            if len(r) < 5:
                continue
            # time 은 UTC ISO → KST(+9h), 초 버림
            t = r[0][:19].replace("T", " ")
            try:
                tt = dt.datetime.fromisoformat(t) + dt.timedelta(hours=9)
            except ValueError:
                continue
            key = tt.strftime("%Y-%m-%d %H:%M:00")
            out[key] = (_clean(r[1], "X"), _clean(r[2], "Y"),
                        _clean(r[3], "Z"), _clean(r[4], "F"))
    return out


MAX_MONTHS = 12   # 범위 상한(전송·메모리 보호)


def _month_iter(sy, sm, ey, em):
    y, m = sy, sm
    while (y, m) <= (ey, em):
        yield y, m
        m += 1
        if m > 12:
            m = 1
            y += 1


def minute_range(start_ym, end_ym, codes):
    """[시작 월 ~ 종료 월] 범위의 분 단위 자료를 1분 격자에 정렬해 반환."""
    sy, sm = map(int, start_ym.split("-"))
    ey, em = map(int, end_ym.split("-"))
    months = list(_month_iter(sy, sm, ey, em))
    t0 = dt.datetime(sy, sm, 1)
    end_last = calendar.monthrange(ey, em)[1]
    t_end = dt.datetime(ey, em, end_last) + dt.timedelta(days=1)   # 배타적 끝
    n = int((t_end - t0).total_seconds() // 60)
    years = sorted({y for y, _ in months})
    days_iso = [f"{y:04d}-{m:02d}-{d:02d}"
                for (y, m) in months for d in range(1, calendar.monthrange(y, m)[1] + 1)]

    out = {"start": start_ym, "end": end_ym,
           "t0": t0.strftime("%Y-%m-%dT%H:%M:%S"), "step": 60, "n": n, "stn": {}}
    for code in codes:
        if code not in STN:
            continue
        if code == "CYG":
            src = {}
            for diso in days_iso:
                src.update(_cyg_day(diso))
        else:
            src = {}
            for y in years:
                src.update(_kasa_year(code, str(y)))
        X = [None] * n
        Y = [None] * n
        Z = [None] * n
        F = [None] * n
        for key, vals in src.items():
            try:
                tt = dt.datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                continue
            idx = int((tt - t0).total_seconds() // 60)
            if 0 <= idx < n:
                X[idx], Y[idx], Z[idx], F[idx] = vals
        name, color = STN[code]
        out["stn"][code] = {"name": name, "c": color, "X": X, "Y": Y, "Z": Z, "F": F}
    return out


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=str(DOCS), **k)

    def log_message(self, *a):
        pass  # 조용히

    def do_GET(self):
        if self.path.startswith("/api/minute"):
            return self._api_minute()
        return super().do_GET()

    def do_POST(self):
        if self.path.startswith("/api/chat"):
            return self._api_chat()
        self.send_error(404)

    def do_OPTIONS(self):   # CORS preflight
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def _api_chat(self):
        try:
            n = int(self.headers.get("Content-Length", 0))
            body = json.loads(self.rfile.read(n).decode("utf-8")) if n else {}
        except Exception as e:            # noqa: BLE001
            return self._send_json({"answer": f"요청 파싱 오류: {e}", "actions": [],
                                    "tier": "data-only"}, 400)
        try:
            import globe_agent
            out = globe_agent.answer(body.get("message", ""))
        except Exception as e:            # noqa: BLE001
            out = {"answer": f"에이전트 오류: {e}", "actions": [], "tier": "data-only"}
        return self._send_json(out)

    def _send_json(self, obj, code=200):
        body = json.dumps(obj, separators=(",", ":")).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def _api_minute(self):
        q = dict(re.findall(r"([^?&=]+)=([^&]*)", self.path))
        # 단일 월(ym) 또는 범위(start,end) 지원. ym 만 오면 start=end=ym.
        start = q.get("start") or q.get("ym", "")
        end = q.get("end") or start
        if not (re.fullmatch(r"\d{4}-\d{2}", start) and re.fullmatch(r"\d{4}-\d{2}", end)):
            return self._send_json({"ok": False, "err": "start/end=YYYY-MM 필요"}, 400)
        sy, sm = map(int, start.split("-"))
        ey, em = map(int, end.split("-"))
        if (ey, em) < (sy, sm):
            return self._send_json({"ok": False, "err": "종료 월이 시작 월보다 빠릅니다"}, 400)
        nmonths = (ey - sy) * 12 + (em - sm) + 1
        if nmonths > MAX_MONTHS:
            return self._send_json(
                {"ok": False, "err": f"최대 {MAX_MONTHS}개월까지 조회 가능(요청 {nmonths}개월)"}, 400)
        codes = [c for c in q.get("stn", "CYG,GN,JJ,ICH").split(",") if c in STN]
        try:
            data = minute_range(start, end, codes)
        except Exception as e:                       # noqa: BLE001
            return self._send_json({"ok": False, "err": str(e)}, 500)
        data["ok"] = True
        return self._send_json(data)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8899)
    ap.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"관측소 뷰어 서버: http://{a.host}:{a.port}/lmm.html")
    print("  분 단위 API: /api/minute?ym=2022-06&stn=GN,JJ,ICH,CYG")
    print("  (Ctrl+C 로 종료)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
