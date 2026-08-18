# -*- coding: utf-8 -*-
"""발표자료용 화면캡처 수집기 — 페이지가 보낸 dataURL 을 PNG 로 저장한다.

헤드리스 브라우저는 WebGL 이 막혀 지구본 캡처가 안 되므로, 실제 브라우저에서
`window.__shot()` 로 뽑은 dataURL 을 여기로 POST 받는다.

    python _shot_server.py          # http://127.0.0.1:8902/save
"""
import base64, json
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

OUT = Path(__file__).parent / "docs" / "output" / "_web_shots"


class H(BaseHTTPRequestHandler):
    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "content-type")

    def do_OPTIONS(self):
        self.send_response(204); self._cors(); self.end_headers()

    def do_POST(self):
        n = int(self.headers.get("content-length", 0))
        body = json.loads(self.rfile.read(n) or b"{}")
        name = "".join(c for c in body.get("name", "shot") if c.isalnum() or c in "._-")
        data = body.get("data", "").split(",", 1)[-1]
        OUT.mkdir(parents=True, exist_ok=True)
        p = OUT / f"{name}.png"
        p.write_bytes(base64.b64decode(data))
        print(f"saved {p} ({p.stat().st_size:,} bytes)", flush=True)
        self.send_response(200); self._cors()
        self.send_header("content-type", "application/json"); self.end_headers()
        self.wfile.write(json.dumps({"ok": True, "path": str(p)}).encode())

    def log_message(self, *a):
        pass


if __name__ == "__main__":
    print("listening on http://127.0.0.1:8902/save", flush=True)
    HTTPServer(("127.0.0.1", 8902), H).serve_forever()
