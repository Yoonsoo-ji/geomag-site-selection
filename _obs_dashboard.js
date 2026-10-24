/* 기존 측정점 자료 점검 대시보드 — 화면 동작
 *
 * 설계: 지도·목록·상세가 「지금 고른 측점」 하나를 공유한다(연결된 선택).
 * 어디서 고르든 나머지가 따라온다.
 *
 * ⚠ 외부 라이브러리를 쓰지 않는다 — 단일 파일로 열려야 하고, 점검 PC 가
 *   인터넷에 닿지 않을 수 있다. 지도·그래프는 캔버스로 직접 그린다.
 */
var D = JSON.parse(document.getElementById('DATA').textContent);
var SEL = null, TAB = '시계열', MODE = 'res', SCALE = 'auto';
var SORT = {k: '측점', a: true};
var FLT = {q: '', f: ''};
var tip = document.getElementById('tip');

/* ── 공통 ───────────────────────────────────────────────── */
function esc(s) {
  return String(s === null || s === undefined ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function num(v, nd, unit) {
  if (v === null || v === undefined || v !== v) return '-';
  var s = Number(v).toLocaleString('ko-KR',
    {minimumFractionDigits: nd, maximumFractionDigits: nd});
  return s + (unit || '');
}
function byStation(arr, s) {
  return arr.filter(function (r) { return r['측점'] === s; });
}
function showTip(e, html) {
  tip.innerHTML = html;
  tip.style.display = 'block';
  var w = tip.offsetWidth, h = tip.offsetHeight;
  var x = e.clientX + 14, y = e.clientY + 14;
  if (x + w > innerWidth - 8) x = e.clientX - w - 14;
  if (y + h > innerHeight - 8) y = e.clientY - h - 14;
  tip.style.left = x + 'px';
  tip.style.top = y + 'px';
}
function hideTip() { tip.style.display = 'none'; }
addEventListener('scroll', hideTip, true);
addEventListener('resize', function () { hideTip(); drawMap(); render(); });

/* ── 지도 ───────────────────────────────────────────────── */
var MAP = {x0: 125.0, x1: 130.1, y0: 33.0, y1: 38.8, pad: 10};
function proj(lon, lat, W, H) {
  /* 위도에 따라 경도 폭이 줄어든다 — 그대로 늘이면 남북이 찌그러진다 */
  var mid = (MAP.y0 + MAP.y1) / 2, k = Math.cos(mid * Math.PI / 180);
  var sx = (W - 2 * MAP.pad) / ((MAP.x1 - MAP.x0) * k);
  var sy = (H - 2 * MAP.pad) / (MAP.y1 - MAP.y0);
  var s = Math.min(sx, sy);
  var ox = MAP.pad + ((W - 2 * MAP.pad) - (MAP.x1 - MAP.x0) * k * s) / 2;
  var oy = MAP.pad + ((H - 2 * MAP.pad) - (MAP.y1 - MAP.y0) * s) / 2;
  return [ox + (lon - MAP.x0) * k * s, oy + (MAP.y1 - lat) * s];
}
function stColor(st) {
  if (st['확인필요'] > 0) return D['색_심각도']['확인필요'];
  if (st['플래그'] > 0) return D['색_심각도']['주의'];
  return '#238B45';
}
function drawMap() {
  var cv = document.getElementById('map');
  var wrap = document.getElementById('mapwrap');
  var W = Math.max(240, wrap.clientWidth || 300);
  var H = Math.round(W * 1.32);
  var r = devicePixelRatio || 1;
  cv.width = W * r; cv.height = H * r;
  cv.style.width = W + 'px'; cv.style.height = H + 'px';
  var c = cv.getContext('2d');
  c.setTransform(r, 0, 0, r, 0, 0);
  c.clearRect(0, 0, W, H);
  c.fillStyle = '#F7FAFD'; c.fillRect(0, 0, W, H);

  c.strokeStyle = '#9FB3C8'; c.lineWidth = 0.8;
  (D['해안선'] || []).forEach(function (ring) {
    c.beginPath();
    ring.forEach(function (p, i) {
      var q = proj(p[0], p[1], W, H);
      if (i === 0) c.moveTo(q[0], q[1]); else c.lineTo(q[0], q[1]);
    });
    c.closePath(); c.fillStyle = '#FFFFFF'; c.fill(); c.stroke();
  });

  MAP.hits = [];
  D['측점'].forEach(function (st) {
    if (st['위도'] === null || st['경도'] === null) return;
    var q = proj(st['경도'], st['위도'], W, H);
    var on = (SEL === st['측점']);
    /* 플래그가 많을수록 크게 — 눈이 먼저 가는 곳이 먼저 볼 곳이다 */
    var rad = 3.6 + Math.min(4.6, (st['플래그'] || 0) * 0.3);
    MAP.hits.push({x: q[0], y: q[1], r: rad + 5, s: st['측점']});
    if (on) {
      c.beginPath(); c.arc(q[0], q[1], rad + 5, 0, 7);
      c.strokeStyle = '#16324F'; c.lineWidth = 2.2; c.stroke();
    }
    c.beginPath(); c.arc(q[0], q[1], rad, 0, 7);
    c.fillStyle = stColor(st); c.fill();
    c.strokeStyle = '#fff'; c.lineWidth = 1.1; c.stroke();
    c.fillStyle = on ? '#16324F' : '#4B5563';
    c.font = (on ? 'bold ' : '') + '10px "맑은 고딕",sans-serif';
    c.fillText(st['측점'], q[0] + rad + 3, q[1] + 3.5);
  });
}
(function () {
  var cv = document.getElementById('map');
  function hit(e) {
    var b = cv.getBoundingClientRect();
    var x = e.clientX - b.left, y = e.clientY - b.top, best = null, bd = 1e9;
    (MAP.hits || []).forEach(function (h) {
      var d = (h.x - x) * (h.x - x) + (h.y - y) * (h.y - y);
      if (d < h.r * h.r && d < bd) { bd = d; best = h; }
    });
    return best;
  }
  cv.addEventListener('click', function (e) {
    var h = hit(e); if (h) select(h.s);
  });
  cv.addEventListener('mousemove', function (e) {
    var h = hit(e);
    cv.style.cursor = h ? 'pointer' : 'default';
    if (!h) { hideTip(); return; }
    var st = D['측점'].filter(function (s) { return s['측점'] === h.s; })[0];
    showTip(e, '<b>' + esc(st['측점']) + '</b><br>' +
      '방문 ' + st['방문'] + ' · ' + (st['최초'] || '') + '~' + (st['최종'] || '') +
      '<br>편각 잔차 폭 ' + num(st['rD_폭'], 1, '분') +
      '<br>확인 필요 ' + (st['확인필요'] || 0) + ' · 플래그 ' + (st['플래그'] || 0));
  });
  cv.addEventListener('mouseleave', hideTip);
})();

/* ── 측점 목록 ──────────────────────────────────────────── */
var LCOLS = [
  {k: '측점', t: '측점'},
  {k: '방문', t: '방문', n: 1},
  {k: '최초', t: '최초', n: 1},
  {k: '최종', t: '최종', n: 1},
  {k: 'rD_폭', t: '편각폭', n: 1, nd: 0},
  {k: '확인필요', t: '확인', n: 1},
  {k: '플래그', t: '전체', n: 1}
];
function listRows() {
  var q = FLT.q.trim(), f = FLT.f;
  var rows = D['측점'].filter(function (s) {
    if (q && s['측점'].indexOf(q) < 0) return false;
    if (f === 'flag' && !(s['확인필요'] > 0)) return false;
    if (f === 'lmm' && s['LMM보유'] !== '있음') return false;
    if (f === 'nolmm' && s['LMM보유'] === '있음') return false;
    return true;
  });
  var k = SORT.k, a = SORT.a ? 1 : -1;
  rows.sort(function (x, y) {
    var u = x[k], v = y[k];
    if (u === null || u === undefined) return 1;
    if (v === null || v === undefined) return -1;
    if (typeof u === 'string') return a * u.localeCompare(v, 'ko');
    return a * (u - v);
  });
  return rows;
}
function drawList() {
  var rows = listRows();
  document.getElementById('listcount').textContent = rows.length + ' / ' +
    D['측점'].length;
  var h = '<thead><tr>' + LCOLS.map(function (c) {
    return '<th class="' + (c.n ? 'n' : '') + '" data-k="' + c.k + '">' +
      c.t + (SORT.k === c.k ? (SORT.a ? ' ▲' : ' ▼') : '') + '</th>';
  }).join('') + '</tr></thead><tbody>';
  rows.forEach(function (s) {
    h += '<tr data-s="' + esc(s['측점']) + '"' +
      (SEL === s['측점'] ? ' class="sel"' : '') + '>' +
      '<td><span class="dot" style="background:' + stColor(s) + '"></span>' +
      esc(s['측점']) + '</td>' +
      '<td class="n">' + s['방문'] + '</td>' +
      '<td class="n">' + (s['최초'] || '-') + '</td>' +
      '<td class="n">' + (s['최종'] || '-') + '</td>' +
      '<td class="n">' + num(s['rD_폭'], 0) + '</td>' +
      '<td class="n"' + (s['확인필요'] > 0 ? ' style="color:#B2182B;font-weight:700"' : '') +
      '>' + (s['확인필요'] || '') + '</td>' +
      '<td class="n muted">' + (s['플래그'] || '') + '</td></tr>';
  });
  document.getElementById('list').innerHTML = h + '</tbody>';
}
document.getElementById('list').addEventListener('click', function (e) {
  var th = e.target.closest('th');
  if (th && th.dataset.k) {
    if (SORT.k === th.dataset.k) SORT.a = !SORT.a;
    else { SORT.k = th.dataset.k; SORT.a = (th.dataset.k === '측점'); }
    drawList(); return;
  }
  var tr = e.target.closest('tr[data-s]');
  if (tr) select(tr.dataset.s);
});
document.getElementById('q').addEventListener('input', function (e) {
  FLT.q = e.target.value; drawList();
});
document.getElementById('filt').addEventListener('change', function (e) {
  FLT.f = e.target.value; drawList();
});

/* ── 시계열 그래프 ──────────────────────────────────────── */
var COMPS = [
  {k: 'rD', v: 'D', t: '편각 D', u: '분', vu: '°', nd: 1, vnd: 4},
  {k: 'rI', v: 'I', t: '복각 I', u: '분', vu: '°', nd: 1, vnd: 4},
  {k: 'rF', v: 'F', t: '총자력 F', u: ' nT', vu: ' nT', nd: 0, vnd: 1}
];
function drawChart(cv, rows, comp) {
  var W = cv.clientWidth || 600, H = 150, r = devicePixelRatio || 1;
  cv.width = W * r; cv.height = H * r;
  cv.style.width = '100%'; cv.style.height = H + 'px';
  var c = cv.getContext('2d');
  c.setTransform(r, 0, 0, r, 0, 0);
  c.clearRect(0, 0, W, H);
  var L = 54, R = 10, T = 16, B = 24;
  var key = MODE === 'res' ? comp.k : comp.v;
  var pts = rows.filter(function (d) {
    return d[key] !== null && d[key] !== undefined && d['연도'] !== null;
  });
  c.fillStyle = '#374151';
  c.font = 'bold 11px "맑은 고딕",sans-serif';
  c.fillText(comp.t + (MODE === 'res' ? ' — IGRF 잔차 (' + comp.u.trim() +
    ')' : ' — 실측값'), L, 11);
  if (!pts.length) {
    c.fillStyle = '#9CA3AF'; c.font = '11px "맑은 고딕",sans-serif';
    c.fillText('자료 없음', L, H / 2); return;
  }
  var ys = pts.map(function (d) { return d[key]; });
  /* ⚠ 이상값 하나가 눈금을 독점하면 나머지가 평평한 선이 되어 아무것도
   *   읽히지 않는다(거제 2014 의 복각 98° 가 실제로 그랬다). 기본은
   *   **중앙값 ± 4×MAD** 로 자르고, 벗어난 점은 **가장자리에 세모**로
   *   그려 둔다 — 지우지 않으므로 몇 개가 어느 쪽으로 벗어났는지 보인다.
   *   「전체 눈금」 단추로 원래 범위를 본다. */
  var y0, y1, clipped = 0;
  var srt = ys.slice().sort(function (a, b) { return a - b; });
  var med = srt[Math.floor(srt.length / 2)];
  var mad = srt.map(function (v) { return Math.abs(v - med); })
    .sort(function (a, b) { return a - b; })[Math.floor(srt.length / 2)];
  if (SCALE === 'auto' && srt.length >= 4 && mad > 0) {
    y0 = med - 4 * mad; y1 = med + 4 * mad;
    ys.forEach(function (v) { if (v < y0 || v > y1) clipped++; });
    if (!clipped) { y0 = Math.min.apply(null, ys); y1 = Math.max.apply(null, ys); }
  } else {
    y0 = Math.min.apply(null, ys); y1 = Math.max.apply(null, ys);
  }
  if (MODE === 'res') { y0 = Math.min(y0, 0); y1 = Math.max(y1, 0); }
  var pad = (y1 - y0) * 0.15 || Math.abs(y1 || 1) * 0.1 || 1;
  y0 -= pad; y1 += pad;
  var xs = pts.map(function (d) { return d['연도']; });
  var x0 = Math.min.apply(null, xs) - 0.6, x1 = Math.max.apply(null, xs) + 0.6;
  function px(y) { return H - B - (y - y0) / (y1 - y0) * (H - T - B); }
  function qx(x) { return L + (x - x0) / (x1 - x0) * (W - L - R); }

  c.strokeStyle = '#E8EBEF'; c.lineWidth = 1;
  c.font = '10px "맑은 고딕",sans-serif'; c.fillStyle = '#6B7280';
  for (var i = 0; i <= 3; i++) {
    var yv = y0 + (y1 - y0) * i / 3, yy = px(yv);
    c.beginPath(); c.moveTo(L, yy); c.lineTo(W - R, yy); c.stroke();
    c.textAlign = 'right';
    c.fillText(num(yv, comp.nd), L - 5, yy + 3);
  }
  c.textAlign = 'center';
  var span = Math.max(1, Math.round((x1 - x0) / 8));
  for (var y = Math.ceil(x0); y <= x1; y++) {
    if (y % span) continue;
    c.fillText(y, qx(y), H - 8);
  }
  if (MODE === 'res' && y0 < 0 && y1 > 0) {
    /* 0 선 — 같은 표석이면 여기에 붙어 있어야 한다 */
    c.strokeStyle = '#B2182B'; c.lineWidth = 1.2; c.setLineDash([4, 3]);
    c.beginPath(); c.moveTo(L, px(0)); c.lineTo(W - R, px(0)); c.stroke();
    c.setLineDash([]);
  }
  /* 출처마다 선을 잇는다 — 누가 무엇을 적었는지가 값만큼 중요하다 */
  var bySrc = {};
  pts.forEach(function (d) { (bySrc[d['출처']] = bySrc[d['출처']] || []).push(d); });
  Object.keys(bySrc).forEach(function (s) {
    var g = bySrc[s].slice().sort(function (a, b) { return a['연도'] - b['연도']; });
    c.strokeStyle = D['색_출처'][s] || '#888'; c.lineWidth = 1.4;
    c.beginPath();
    var started = false;
    g.forEach(function (d) {
      var v = d[key];
      if (v < y0 || v > y1) { started = false; return; }
      var X = qx(d['연도']), Y = px(v);
      if (!started) { c.moveTo(X, Y); started = true; } else c.lineTo(X, Y);
    });
    c.stroke();
  });
  cv._hits = [];
  pts.forEach(function (d) {
    var v = d[key], out = (v < y0 || v > y1);
    var X = qx(d['연도']);
    var Y = out ? (v > y1 ? T + 6 : H - B - 6) : px(v);
    var bad = (d['D판정'] === '불가' || d[comp.v + '판정'] === '불가');
    c.fillStyle = D['색_출처'][d['출처']] || '#888';
    c.strokeStyle = bad || out ? '#B2182B' : '#fff';
    c.lineWidth = (bad || out) ? 2 : 1.1;
    if (out) {                      /* 눈금 밖 — 가장자리에 세모 */
      var s = v > y1 ? -1 : 1;
      c.beginPath();
      c.moveTo(X, Y + s * 5.5); c.lineTo(X - 5, Y - s * 3.5);
      c.lineTo(X + 5, Y - s * 3.5); c.closePath();
      c.fill(); c.stroke();
    } else {
      c.beginPath(); c.arc(X, Y, bad ? 4.6 : 3.4, 0, 7);
      c.fill(); c.stroke();
    }
    cv._hits.push({x: X, y: Y, d: d, out: out});
  });
  if (clipped) {
    c.fillStyle = '#B2182B';
    c.font = 'bold 10px "맑은 고딕",sans-serif';
    c.textAlign = 'right';
    c.fillText('눈금 밖 ' + clipped + '점 ▲▼', W - R, 11);
    c.textAlign = 'left';
  }
  cv.onmousemove = function (e) {
    var b = cv.getBoundingClientRect(), best = null, bd = 169;
    (cv._hits || []).forEach(function (h) {
      var dd = Math.pow(h.x - (e.clientX - b.left), 2) +
        Math.pow(h.y - (e.clientY - b.top), 2);
      if (dd < bd) { bd = dd; best = h; }
    });
    if (!best) { hideTip(); cv.style.cursor = 'default'; return; }
    cv.style.cursor = 'pointer';
    var d = best.d;
    showTip(e, '<b>' + d['연도'] + ' · ' + esc(d['출처']) + '</b><br>' +
      '일자 ' + esc(d['일자']) + ' · 세션 ' + d['세션'] + '<br>' +
      '실측 ' + num(d[comp.v], comp.vnd, comp.vu) + '<br>' +
      '잔차 ' + num(d[comp.k], comp.nd, comp.u) +
      (d[comp.v + '산포'] !== null ? '<br>같은 날 산포 ' +
        num(d[comp.v + '산포'], 2, comp.u) : '') +
      (d[comp.v + '판정'] ? '<br>판정 ' + esc(d[comp.v + '판정']) : '') +
      (best.out ? '<br><b style="color:#FF9B9B">눈금 밖 — 「전체 눈금」으로 봅니다</b>' : ''));
  };
  cv.onmouseleave = hideTip;
}

/* ── 상세 ───────────────────────────────────────────────── */
var TABS = ['시계열', '방문 이력', '확인 필요', '원본 위치', '세션 원자료'];
function select(s) {
  SEL = s; drawMap(); drawList(); render();
  var el = document.querySelector('tr.sel');
  if (el && el.scrollIntoView) el.scrollIntoView({block: 'nearest'});
}
function render() {
  var body = document.getElementById('body');
  document.getElementById('tabs').innerHTML = TABS.map(function (t) {
    return '<div class="tab' + (TAB === t ? ' on' : '') + '" data-t="' + t +
      '">' + t + '</div>';
  }).join('');
  if (!SEL) {
    document.getElementById('title').firstChild.nodeValue = '측점을 고르세요';
    document.getElementById('subtitle').textContent = '';
    body.innerHTML = '<div class="note">왼쪽 지도나 목록에서 측점을 누르면 ' +
      '그 측점의 연도별 자료·플래그·원본 위치가 여기에 열립니다. ' +
      '머리의 <b>확인 필요</b> 숫자를 누르면 전 측점 플래그를 한 번에 봅니다.</div>' +
      '<h3>지도 색</h3><div class="legend">' +
      '<span><span class="dot" style="background:#B2182B"></span>확인 필요 있음</span>' +
      '<span><span class="dot" style="background:#E08214"></span>주의·참고만</span>' +
      '<span><span class="dot" style="background:#238B45"></span>플래그 없음</span>' +
      '<span class="muted">원이 클수록 플래그가 많습니다</span></div>';
    return;
  }
  var st = D['측점'].filter(function (x) { return x['측점'] === SEL; })[0];
  var vs = byStation(D['방문'], SEL);
  var fs = byStation(D['플래그'], SEL);
  document.getElementById('title').firstChild.nodeValue = SEL;
  document.getElementById('subtitle').innerHTML =
    (st['소재지'] ? esc(st['소재지']) + ' · ' : '') +
    (st['점의부호'] ? esc(st['점의부호']) + ' · ' : '') +
    num(st['위도'], 5) + ', ' + num(st['경도'], 5) +
    ' · 방문 ' + st['방문'] + ' · ' + (st['최초'] || '') + '~' + (st['최종'] || '');

  if (TAB === '시계열') body.innerHTML = viewChart(st, vs, fs);
  else if (TAB === '방문 이력') body.innerHTML = viewVisits(vs);
  else if (TAB === '확인 필요') body.innerHTML = viewFlags(fs, false);
  else if (TAB === '원본 위치') body.innerHTML = viewSrc(vs);
  else body.innerHTML = viewSessions(byStation(D['세션'], SEL));

  if (TAB === '시계열') {
    COMPS.forEach(function (c, i) {
      drawChart(document.getElementById('ch' + i), vs, c);
    });
  }
}
document.getElementById('tabs').addEventListener('click', function (e) {
  var t = e.target.closest('.tab');
  if (t) { TAB = t.dataset.t; render(); }
});
function viewChart(st, vs, fs) {
  var nf = fs.filter(function (f) { return f['심각도'] === '확인필요'; }).length;
  var h = '';
  if (nf) {
    h += '<div class="note bad"><b>확인 필요 ' + nf + '건</b> — ' +
      '「확인 필요」 탭에서 봅니다.</div>';
  }
  h += '<div class="tools"><b>표시</b>' +
    '<button id="mres" class="' + (MODE === 'res' ? 'on' : '') + '">IGRF 잔차</button>' +
    '<button id="mval" class="' + (MODE === 'val' ? 'on' : '') + '">실측값</button>' +
    '<b style="margin-left:10px">눈금</b>' +
    '<button id="sauto" class="' + (SCALE === 'auto' ? 'on' : '') +
    '">이상값 접기</button>' +
    '<button id="sall" class="' + (SCALE === 'all' ? 'on' : '') +
    '">전체</button>' +
    '<span class="muted" style="margin-left:6px">' +
    (MODE === 'res'
      ? '잔차는 영년변화가 빠진 양입니다. 같은 표석이면 붉은 0선에 붙어 있어야 합니다.'
      : '실측값은 영년변화로 저절로 기울어 있습니다.') + '</span></div>';
  h += '<div class="legend">' + Object.keys(D['색_출처']).map(function (s) {
    return '<span><span class="dot" style="background:' + D['색_출처'][s] +
      '"></span>' + esc(s) + '</span>';
  }).join('') + '</div>';
  h += COMPS.map(function (c, i) {
    return '<canvas class="chart" id="ch' + i + '"></canvas>';
  }).join('');
  h += '<h3>이 측점의 변동</h3><table><tbody>' +
    tr('편각 잔차 중앙 / 폭', num(st['rD_중앙'], 1, '분') + ' / ' +
      num(st['rD_폭'], 1, '분'), st['rD_폭'] > 60) +
    tr('총자력 잔차 중앙 / 폭', num(st['rF_중앙'], 0, ' nT') + ' / ' +
      num(st['rF_폭'], 0, ' nT'), st['rF_폭'] > 1000) +
    tr('같은 날 세션 산포(중앙)', num(st['같은날산포'], 2, '분'), false) +
    tr('출처 수 / 모델 투입 방문', st['출처수'] + ' / ' + st['투입D'], false) +
    tr('현재 모델 보유', st['LMM보유'], false) +
    '</tbody></table>';
  return h;
}
function tr(k, v, bad) {
  return '<tr><td style="color:#4B5563">' + k + '</td><td class="' +
    (bad ? 'warn' : '') + '">' + v + '</td></tr>';
}
function viewVisits(vs) {
  if (!vs.length) return '<div class="note">방문 자료가 없습니다.</div>';
  var h = '<div class="note">같은 해를 여러 출처가 적은 경우 <b>모두</b> 싣습니다 — ' +
    '점검에서는 어느 문서가 무엇을 적었는지가 값만큼 중요합니다.</div>' +
    '<table><thead><tr><th>연도</th><th>일자</th><th>출처</th>' +
    '<th class="n">세션</th><th class="n">편각</th><th class="n">복각</th>' +
    '<th class="n">총자력</th><th class="n">rD</th><th class="n">rI</th>' +
    '<th class="n">rF</th><th>판정 D/I/F</th><th>비고</th></tr></thead><tbody>';
  vs.slice().sort(function (a, b) { return (a['연도'] || 0) - (b['연도'] || 0); })
    .forEach(function (d) {
      var note = [];
      if (d['좌표신뢰'] && d['좌표신뢰'].indexOf('이설') >= 0) note.push('이설 이전');
      if (d['일자정밀도'] === '연') note.push('연도만');
      if (d['일자정밀도'] === '월') note.push('월까지');
      h += '<tr><td class="n">' + (d['연도'] || '-') + '</td><td>' +
        esc(d['일자']) + '</td><td><span class="dot" style="background:' +
        (D['색_출처'][d['출처']] || '#888') + '"></span>' + esc(d['출처']) +
        '</td><td class="n">' + d['세션'] + '</td>' +
        '<td class="n">' + num(d['D'], 4) + '</td>' +
        '<td class="n">' + num(d['I'], 4) + '</td>' +
        '<td class="n">' + num(d['F'], 1) + '</td>' +
        '<td class="n' + (Math.abs(d['rD'] || 0) > 60 ? ' warn' : '') + '">' +
        num(d['rD'], 1) + '</td>' +
        '<td class="n' + (Math.abs(d['rI'] || 0) > 60 ? ' warn' : '') + '">' +
        num(d['rI'], 1) + '</td>' +
        '<td class="n' + (Math.abs(d['rF'] || 0) > 400 ? ' warn' : '') + '">' +
        num(d['rF'], 0) + '</td>' +
        '<td class="muted">' + esc(d['D판정'].slice(0, 2)) + '/' +
        esc(d['I판정'].slice(0, 2)) + '/' + esc(d['F판정'].slice(0, 2)) +
        '</td><td class="muted">' + esc(note.join(' · ')) + '</td></tr>';
    });
  return h + '</tbody></table><div class="legend"><span class="muted">' +
    'rD·rI 는 분, rF 는 nT. 붉은 값은 한계(편각·복각 60분 · 총자력 400 nT)를 ' +
    '넘은 것입니다.</span></div>';
}
function viewFlags(fs, showSite) {
  if (!fs.length) return '<div class="note">확인이 필요한 항목이 없습니다.</div>';
  /* ⚠ 머리글·칸·<col> 을 세 곳에 따로 적으면 조건이 붙는 순간 어긋난다.
   *   실제로 「측점」 칸을 showSite 가 거짓일 때 **아예 안 내보내면서**
   *   <col> 은 8개 그대로 둬서 열이 한 칸씩 밀렸고, 「종류」가 46px 칸에
   *   들어가 글자가 석 줄로 접혔다.
   *   → 열 정의를 **한 배열**에 두고 머리글·칸·<col> 을 모두 거기서 만든다. */
  var COL = [
    {t: '심각도', w: 62, f: function (f) {
      return '<span class="sev" style="background:' +
        (D['색_심각도'][f['심각도']] || '#888') + '">' + esc(f['심각도']) +
        '</span>'; }},
    {t: '측점', w: 60, on: showSite, f: function (f) {
      return '<a href="#" onclick="select(\'' + esc(f['측점']) +
        '\');return false">' + esc(f['측점']) + '</a>'; }},
    {t: '연도', w: 48, cls: 'n', f: function (f) { return f['연도'] || '-'; }},
    {t: '종류', w: 132, f: function (f) { return esc(f['종류']); }},
    {t: '내용', w: 218, f: function (f) {
      return esc(f['내용']) + (f['상세']
        ? '<br><span class="muted">' + esc(f['상세']) + '</span>' : ''); }},
    {t: '왜 보나', min: 210, cls: 'muted', f: function (f) {
      return esc(f['왜']); }},
    {t: '출처', w: 136, cls: 'muted', f: function (f) {
      return esc(f['출처']); }},
    {t: '원본 파일', w: 272, cls: 'path', f: function (f) {
      return f['원본파일'] ? esc(f['원본파일']) : '-'; }}
  ].filter(function (c) { return c.on !== false; });

  var h = '<table class="fx"><colgroup>' + COL.map(function (c) {
    return '<col style="' + (c.w ? 'width:' + c.w + 'px'
      : 'min-width:' + c.min + 'px') + '">';
  }).join('') + '</colgroup><thead><tr>' + COL.map(function (c) {
    return '<th' + (c.cls === 'n' ? ' class="n"' : '') + '>' + c.t + '</th>';
  }).join('') + '</tr></thead><tbody>';
  fs.forEach(function (f) {
    h += '<tr>' + COL.map(function (c) {
      return '<td' + (c.cls ? ' class="' + c.cls + '"' : '') + '>' +
        c.f(f) + '</td>';
    }).join('') + '</tr>';
  });
  return h + '</tbody></table>';
}
function viewSrc(vs) {
  var seen = {}, rows = [];
  vs.forEach(function (d) {
    var k = d['출처'] + '|' + d['원본파일'];
    if (seen[k]) { seen[k].years.push(d['연도']); return; }
    seen[k] = {s: d['출처'], f: d['원본파일'], sh: d['원본시트'],
               years: [d['연도']]};
    rows.push(seen[k]);
  });
  var h = '<div class="note">각 값이 <b>어느 파일의 어느 시트</b>에서 왔는지입니다. ' +
    '경로를 복사해 탐색기 주소창에 붙여 넣으면 그 파일이 열립니다.</div>' +
    '<table><thead><tr><th>출처</th><th>연도</th><th>시트</th>' +
    '<th>원본 파일 (전체 경로)</th></tr></thead><tbody>';
  rows.forEach(function (r) {
    var ys = r.years.filter(function (y) { return y; })
      .sort(function (a, b) { return a - b; });
    h += '<tr><td><span class="dot" style="background:' +
      (D['색_출처'][r.s] || '#888') + '"></span>' + esc(r.s) + '</td>' +
      '<td class="muted">' + esc(ys.join(', ')) + '</td>' +
      '<td class="muted">' + esc(r.sh || '-') + '</td>' +
      '<td class="path">' + (r.f ? esc(r.f) : '-') + '</td></tr>';
  });
  return h + '</tbody></table>';
}
function viewSessions(ss) {
  if (!ss.length) return '<div class="note">세션 자료가 없습니다.</div>';
  var h = '<div class="note">회차별 원시 판독값입니다 — 가장 원본에 가깝습니다. ' +
    '방문 대표값은 여기서 중앙값으로 묶은 것입니다.</div>' +
    '<table><thead><tr><th>연도</th><th>일자</th><th>출처</th>' +
    '<th class="n">회차</th><th class="n">편각</th><th class="n">복각</th>' +
    '<th class="n">총자력</th><th class="n">rD</th><th>판정</th>' +
    '</tr></thead><tbody>';
  ss.forEach(function (d) {
    var ok = (d['판정'] === '정상/정상/정상');
    h += '<tr><td class="n">' + (d['연도'] || '-') + '</td><td>' +
      esc(d['일자']) + '</td><td class="muted">' + esc(d['출처']) + '</td>' +
      '<td class="n">' + (d['회차'] === null ? '-' : d['회차']) + '</td>' +
      '<td class="n">' + num(d['D'], 4) + '</td>' +
      '<td class="n">' + num(d['I'], 4) + '</td>' +
      '<td class="n">' + num(d['F'], 1) + '</td>' +
      '<td class="n">' + num(d['rD'], 1) + '</td>' +
      '<td class="' + (ok ? 'muted' : 'warn') + '">' + esc(d['판정']) +
      '</td></tr>';
  });
  return h + '</tbody></table>';
}

/* ── 전 측점 플래그 ─────────────────────────────────────── */
function allFlags(sev) {
  SEL = null; drawMap(); drawList();
  var fs = D['플래그'].filter(function (f) {
    return !sev || f['심각도'] === sev;
  });
  document.getElementById('tabs').innerHTML = '';
  document.getElementById('title').firstChild.nodeValue =
    '확인 필요 — 전 측점' + (sev ? ' (' + sev + ')' : '');
  document.getElementById('subtitle').textContent = fs.length + '건';
  var by = {};
  fs.forEach(function (f) { by[f['종류']] = (by[f['종류']] || 0) + 1; });
  var h = '<div class="legend">' + Object.keys(by).sort(function (a, b) {
    return by[b] - by[a];
  }).map(function (k) {
    return '<span><b>' + esc(k) + '</b> ' + by[k] + '</span>';
  }).join('') + '</div>';
  h += '<div class="note">측점 이름을 누르면 그 측점으로 갑니다. ' +
    '자동으로 버리지 않습니다 — 이 목록은 <b>확인할 것</b>이지 ' +
    '<b>배제할 것</b>이 아닙니다.</div>';
  document.getElementById('body').innerHTML = h + viewFlags(fs, true);
}

/* ── 시작 ───────────────────────────────────────────────── */
document.addEventListener('click', function (e) {
  if (e.target.id === 'mres') { MODE = 'res'; render(); }
  if (e.target.id === 'mval') { MODE = 'val'; render(); }
  if (e.target.id === 'sauto') { SCALE = 'auto'; render(); }
  if (e.target.id === 'sall') { SCALE = 'all'; render(); }
});
document.getElementById('maplegend').innerHTML =
  '<span><span class="dot" style="background:#B2182B"></span>확인 필요</span>' +
  '<span><span class="dot" style="background:#E08214"></span>주의·참고</span>' +
  '<span><span class="dot" style="background:#238B45"></span>없음</span>';
drawMap(); drawList(); render();
