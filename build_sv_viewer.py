# -*- coding: utf-8 -*-
"""
docs/lmm.html 에 '관측소 영년변화 시계열' 대화형 뷰어를 내장한다.

  · 일 중앙값(장기): 청양·강릉·제주·이천 X·Y·Z·F 를 일 중앙값으로 압축해 JSON 내장
    (오프라인 단일 파일 유지). 성분 탭·연도 버튼·기간 검색·관측소 토글.
  · 분 단위(월별): sv_server.py 의 /api/minute 에서 해당 월 분자료를 on-demand 로 fetch.
    서버가 없으면(순수 오프라인) 안내 메시지로 폴백.

  · 멱등: 마커(SV_CARD / SV_JS)로 재삽입·교체. CYG 수집 갱신 시 재실행만.

    python build_sv_viewer.py

주의: JS 문자열 안에 "\n" 을 쓰면 파일에 실제 개행으로 들어가 문자열이 끊긴다.
      라벨은 개행 없이 조립할 것.
"""
import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"C:\LG_gram_backup_users\LX\2026_geomag")
HTML = ROOT / "docs" / "lmm.html"
KA = ROOT / "data" / "kasa"
CY = ROOT / "data" / "cyg"

T0 = pd.Timestamp("2019-01-01")
T1 = pd.Timestamp("2025-12-31")
DAYS = pd.date_range(T0, T1, freq="D")
BAND = {"X": (20000, 40000), "Y": (-9000, 3000), "Z": (28000, 46000), "F": (40000, 60000)}


def _daily(series):
    d = series.resample("D").median().reindex(DAYS)
    return [None if pd.isna(v) else int(round(v)) for v in d.values]


def load_kasa(code):
    fr = [pd.read_csv(f, encoding="utf-8-sig") for f in
          sorted(glob.glob(str(KA / f"kasa_{code}_20*.csv")))]
    d = pd.concat(fr, ignore_index=True)
    t = pd.to_datetime(d["time_kst"], errors="coerce")
    out = {}
    for comp, col in [("X", "X_nT"), ("Y", "Y_nT"), ("Z", "Z_nT"), ("F", "F_nT")]:
        v = pd.to_numeric(d[col], errors="coerce")
        lo, hi = BAND[comp]
        v = v.where((v >= lo) & (v <= hi))
        out[comp] = _daily(pd.Series(v.values, index=t).dropna().sort_index())
    return out


def load_cyg():
    fr = []
    for f in glob.glob(str(CY / "**" / "*.csv"), recursive=True):
        try:
            fr.append(pd.read_csv(f, usecols=["time", "X", "Y", "Z", "F"]))
        except Exception:
            pass
    d = pd.concat(fr, ignore_index=True).drop_duplicates("time")
    t = (pd.to_datetime(d["time"], utc=True, errors="coerce")
         .dt.tz_convert("Asia/Seoul").dt.tz_localize(None))
    out = {}
    for comp in ("X", "Y", "Z", "F"):
        v = pd.to_numeric(d[comp], errors="coerce")
        out[comp] = _daily(pd.Series(v.values, index=t).dropna().sort_index())
    return out


STN = [("CYG", "청양", "#1f6feb"), ("GN", "강릉", "#6b7280"),
       ("JJ", "제주", "#0E8A6B"), ("ICH", "이천", "#c86a2b")]

data = {"t0": T0.strftime("%Y-%m-%d"), "step": 86400, "n": len(DAYS), "stn": {}}
for code, name, color in STN:
    comps = load_cyg() if code == "CYG" else load_kasa(code)
    data["stn"][code] = {"name": name, "c": color, **comps}
    nF = sum(1 for v in comps["F"] if v is not None)
    print(f"  {name}({code}): F 유효일 {nF}/{len(DAYS)}")

SVDATA = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
print(f"SVDATA JSON: {len(SVDATA)//1024} KB")

# ── 카드 HTML ─────────────────────────────────────────────────────────────
CARD = """  <!-- SV_CARD -->
  <div class="card" id="sv-card">
    <div class="maphead">
      <div>
        <div class="mtitle">관측소 영년변화 시계열</div>
        <div class="msub">청양·강릉·제주·이천 · 검색 가능 (일 중앙값 장기 / 분 단위 월별)</div>
      </div>
      <div class="seg" id="svcomp"></div>
    </div>
    <div class="seg" id="svmode" style="margin-bottom:12px"></div>
    <div id="svdaily">
      <div class="seg" id="svyear" style="margin-bottom:10px"></div>
    </div>
    <div id="svminute" style="display:none">
      <div class="grid" style="grid-template-columns:repeat(auto-fit,minmax(130px,1fr));\
align-items:end;margin-bottom:8px">
        <div><label>시작 월</label><input id="svmon0" type="month" min="2019-01" max="2025-12"\
 value="2022-06"></div>
        <div><label>종료 월 (최대 12개월)</label><input id="svmon1" type="month" min="2019-01"\
 max="2025-12" value="2022-06"></div>
        <div><button id="svload" style="margin-top:0">분 단위 불러오기</button></div>
      </div>
      <div class="lgnote" id="svmstat" style="margin-bottom:8px"></div>
    </div>
    <div class="grid" style="grid-template-columns:repeat(auto-fit,minmax(130px,1fr));\
margin-bottom:10px">
      <div><label>시작일 (구간 좁히기)</label><input id="svstart" type="date"></div>
      <div><label>종료일</label><input id="svend" type="date"></div>
    </div>
    <div id="svstn" style="display:flex;flex-wrap:wrap;gap:14px;margin-bottom:10px;\
font-size:13px"></div>
    <div class="mapwrap" style="max-width:none">
      <canvas id="svchart" style="width:100%;display:block"></canvas>
      <div id="svtip" class="tip"></div>
    </div>
    <div class="lgnote" id="svnote"></div>
  </div>
  <!-- /SV_CARD -->
"""

# ── 뷰어 JS ───────────────────────────────────────────────────────────────
JS = r"""  <!-- SV_JS -->
  <script>window.SVDATA = __SVDATA__;</script>
  <script>
  (function(){
    var D = window.SVDATA; if(!D) return;
    var codes = Object.keys(D.stn);
    var ACT = D;                                   // 활성 소스(일=D, 분=fetch 결과)
    var state = { comp:"F", mode:"daily", sel:{}, i0:0, i1:D.n-1 };
    codes.forEach(function(c){ state.sel[c]=true; });
    var $ = function(id){ return document.getElementById(id); };

    function t0ms(src){ var s=src.t0; if(s.length<=10) s+="T00:00:00"; return Date.parse(s+"Z"); }
    function idxDate(i){ return new Date(t0ms(ACT)+i*ACT.step*1000); }
    function idxISO(i){ return idxDate(i).toISOString().slice(0,10); }
    function dateToIdx(str){ return Math.round((Date.parse(str+"T00:00:00Z")-t0ms(ACT))/(ACT.step*1000)); }
    function pad2(x){ return (x<10?"0":"")+x; }

    // 성분 탭
    var comps=[["X","X"],["Y","Y"],["Z","Z"],["F","총자력 F"]];
    $("svcomp").innerHTML = comps.map(function(c){
      return '<button data-c="'+c[0]+'" aria-pressed="'+(c[0]==state.comp)+'">'+c[1]+'</button>'; }).join("");
    $("svcomp").onclick=function(e){ var b=e.target.closest("button"); if(!b)return; state.comp=b.dataset.c;
      [].forEach.call(this.children,function(x){x.setAttribute("aria-pressed",x.dataset.c==state.comp);}); draw(); };

    // 모드 탭
    var modes=[["daily","일 중앙값 · 장기"],["minute","분 단위 · 월별"]];
    $("svmode").innerHTML = modes.map(function(m){
      return '<button data-m="'+m[0]+'" aria-pressed="'+(m[0]=="daily")+'">'+m[1]+'</button>'; }).join("");
    $("svmode").onclick=function(e){ var b=e.target.closest("button"); if(!b)return; var m=b.dataset.m;
      [].forEach.call(this.children,function(x){x.setAttribute("aria-pressed",x.dataset.m==m);});
      state.mode=m;
      $("svdaily").style.display = (m=="daily")?"":"none";
      $("svminute").style.display = (m=="minute")?"":"none";
      if(m=="daily"){ ACT=D; setActiveFull(); }
      else { if(ACT===D){ $("svmstat").innerHTML="월을 고르고 <b>분 단위 불러오기</b>를 누르세요."; } else { setActiveFull(); } }
      draw(); };

    // 연도 버튼(일 모드)
    var years=["2019","2020","2021","2022","2023","2024","2025","전체"];
    $("svyear").innerHTML = years.map(function(y){return '<button data-y="'+y+'">'+y+'</button>';}).join("");
    $("svyear").onclick=function(e){ var b=e.target.closest("button"); if(!b)return; var y=b.dataset.y;
      if(y=="전체"){ setRange(D.t0, idxISO(D.n-1)); } else { setRange(y+"-01-01", y+"-12-31"); } };

    // 관측소 토글
    $("svstn").innerHTML = codes.map(function(c){ var s=D.stn[c];
      return '<label style="display:inline-flex;align-items:center;gap:5px;cursor:pointer;color:var(--ink)">'
        +'<input type="checkbox" data-c="'+c+'" checked style="width:auto">'
        +'<span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:'+s.c+'"></span>'
        +s.name+'</label>'; }).join("");
    $("svstn").onchange=function(e){ var c=e.target.dataset.c; if(!c)return; state.sel[c]=e.target.checked; draw(); };

    // 날짜 입력(일 모드에서 사용)
    function setRange(a,b){ $("svstart").value=a; $("svend").value=b;
      state.i0=Math.max(0,dateToIdx(a)); state.i1=Math.min(ACT.n-1,dateToIdx(b));
      if(state.i1<state.i0){var t=state.i0;state.i0=state.i1;state.i1=t;} draw(); }
    function setActiveFull(){ state.i0=0; state.i1=ACT.n-1;
      $("svstart").min=$("svend").min=idxISO(0); $("svstart").max=$("svend").max=idxISO(ACT.n-1);
      $("svstart").value=idxISO(0); $("svend").value=idxISO(ACT.n-1); }
    $("svstart").onchange=$("svend").onchange=function(){
      state.i0=Math.max(0,dateToIdx($("svstart").value||idxISO(0)));
      state.i1=Math.min(ACT.n-1,dateToIdx($("svend").value||idxISO(ACT.n-1)));
      if(state.i1<state.i0){var t=state.i0;state.i0=state.i1;state.i1=t;} draw(); };

    // 분 단위 로드(월 범위)
    $("svload").onclick=function(){ var a=$("svmon0").value, b=$("svmon1").value;
      if(!a)return; if(!b)b=a;
      if(b<a){ var t=a;a=b;b=t; $("svmon0").value=a; $("svmon1").value=b; }
      var nm=(parseInt(b.slice(0,4))-parseInt(a.slice(0,4)))*12+(parseInt(b.slice(5,7))-parseInt(a.slice(5,7)))+1;
      if(nm>12){ $("svmstat").innerHTML='<span style="color:var(--warn)">최대 12개월까지 조회 가능</span> (요청 '+nm+'개월). 범위를 줄이세요.'; return; }
      $("svmstat").innerHTML="불러오는 중… ("+a+(a==b?"":" ~ "+b)+", "+nm+"개월)";
      fetch("/api/minute?start="+a+"&end="+b).then(function(r){return r.json();}).then(function(j){
        if(!j || !j.ok) throw new Error(j&&j.err||"실패");
        ACT=j; state.mode="minute";
        setActiveFull(); draw();
        var ok=[]; codes.forEach(function(c){ var arr=j.stn[c]; if(!arr)return;
          var k=arr.F.filter(function(v){return v!=null;}).length; if(k>0) ok.push(D.stn[c].name+" "+k.toLocaleString()+"분"); });
        $("svmstat").innerHTML=a+(a==b?"":" ~ "+b)+" 분 단위 로드됨 — "+(ok.join(" · ")||"유효 자료 없음");
      }).catch(function(err){
        $("svmstat").innerHTML='<span style="color:var(--warn)">분 단위는 로컬 서버가 필요합니다.</span> '
          +'터미널에서 <code>python sv_server.py</code> 실행 후 이 페이지를 서버 주소로 여세요. ('+err.message+')';
      }); };

    // ── 캔버스 차트 ──
    var cv=$("svchart"), tip=$("svtip"), pad={l:66,r:16,t:14,b:40};
    function css(v){ return getComputedStyle(document.body).getPropertyValue(v).trim(); }
    function ticks(i0,i1){
      // 가시 구간에 맞는 눈금 생성: [{i, label}]
      var spanD=(i1-i0)*ACT.step/86400, out=[], d0=idxDate(i0), d1=idxDate(i1);
      function push(dt,lab){ var i=Math.round((dt.getTime()-t0ms(ACT))/(ACT.step*1000));
        if(i>=i0&&i<=i1) out.push({i:i,label:lab}); }
      if(spanD>420){ for(var y=d0.getUTCFullYear();;y++){ var dt=new Date(Date.UTC(y,0,1));
        if(dt>d1)break; push(dt,""+y);} }
      else if(spanD>75){ var c=new Date(Date.UTC(d0.getUTCFullYear(),d0.getUTCMonth(),1));
        while(c<=d1){ push(c,(c.getUTCMonth()==0)?""+c.getUTCFullYear():(c.getUTCMonth()+1)+"월");
          c=new Date(Date.UTC(c.getUTCFullYear(),c.getUTCMonth()+1,1)); } }
      else if(spanD>5){ var st=Math.max(1,Math.ceil(spanD/12));
        var c2=new Date(Date.UTC(d0.getUTCFullYear(),d0.getUTCMonth(),d0.getUTCDate()));
        while(c2<=d1){ push(c2,(c2.getUTCMonth()+1)+"/"+c2.getUTCDate());
          c2=new Date(c2.getTime()+st*86400000); } }
      else { var hrs=(i1-i0)*ACT.step/3600, hst=hrs>18?6:(hrs>6?3:1);
        var c3=new Date(d0.getTime()); c3.setUTCMinutes(0,0,0); c3.setUTCHours(Math.ceil(c3.getUTCHours()/hst)*hst);
        while(c3<=d1){ push(c3, c3.getUTCDate()+"일 "+pad2(c3.getUTCHours())+"시");
          c3=new Date(c3.getTime()+hst*3600000); } }
      return out;
    }
    function draw(){
      var W=cv.clientWidth||cv.parentNode.clientWidth||760, H=340, dpr=window.devicePixelRatio||1;
      cv.width=W*dpr; cv.height=H*dpr; cv.style.height=H+"px";
      var g=cv.getContext("2d"); g.setTransform(dpr,0,0,dpr,0,0); g.clearRect(0,0,W,H);
      var muted=css("--muted")||"#5b6672", line=css("--line")||"#dfe3e8";
      var i0=state.i0,i1=state.i1,comp=state.comp, pw=W-pad.l-pad.r, ph=H-pad.t-pad.b;
      var lo=Infinity,hi=-Infinity;
      codes.forEach(function(c){ if(!state.sel[c]||!ACT.stn[c])return; var a=ACT.stn[c][comp];
        for(var i=i0;i<=i1;i++){var v=a[i]; if(v!=null){ if(v<lo)lo=v; if(v>hi)hi=v; }}});
      if(!isFinite(lo)){ g.fillStyle=muted; g.font="13px sans-serif";
        g.fillText("선택 구간에 자료 없음",pad.l,pad.t+ph/2); cv._st=null; return; }
      if(hi==lo){hi+=1;lo-=1;} var mg=(hi-lo)*0.08; lo-=mg; hi+=mg;
      var X=function(i){ return pad.l+(i-i0)/Math.max(1,(i1-i0))*pw; };
      var Y=function(v){ return pad.t+(hi-v)/(hi-lo)*ph; };
      g.strokeStyle=line; g.fillStyle=muted; g.font="11px sans-serif"; g.lineWidth=1;
      g.textAlign="right"; g.textBaseline="middle";
      for(var k=0;k<=4;k++){ var vv=lo+(hi-lo)*k/4, yy=Y(vv);
        g.globalAlpha=.6; g.beginPath(); g.moveTo(pad.l,yy); g.lineTo(W-pad.r,yy); g.stroke();
        g.globalAlpha=1; g.fillText(Math.round(vv).toLocaleString(),pad.l-8,yy); }
      g.textAlign="center"; g.textBaseline="top";
      ticks(i0,i1).forEach(function(t){ var xx=X(t.i);
        g.strokeStyle=line; g.globalAlpha=.5; g.beginPath(); g.moveTo(xx,pad.t); g.lineTo(xx,pad.t+ph); g.stroke();
        g.globalAlpha=1; g.fillStyle=muted; g.fillText(t.label,xx,pad.t+ph+6); });
      codes.forEach(function(c){ if(!state.sel[c]||!ACT.stn[c])return; var a=ACT.stn[c][comp];
        g.strokeStyle=D.stn[c].c; g.lineWidth=1.4; g.beginPath(); var pen=false;
        for(var i=i0;i<=i1;i++){ var v=a[i]; if(v==null){pen=false;continue;}
          var xx=X(i),yy=Y(v); if(!pen){g.moveTo(xx,yy);pen=true;} else g.lineTo(xx,yy); } g.stroke(); });
      cv._st={i0:i0,i1:i1,pw:pw};
    }

    cv.onmousemove=function(e){ var st=cv._st; if(!st){tip.style.opacity=0;return;}
      var r=cv.getBoundingClientRect(), mx=e.clientX-r.left;
      if(mx<pad.l||mx>cv.clientWidth-pad.r){tip.style.opacity=0;return;}
      var i=Math.round(st.i0+(mx-pad.l)/st.pw*(st.i1-st.i0)); i=Math.max(st.i0,Math.min(st.i1,i));
      var dtt=idxDate(i), lab=dtt.getUTCFullYear()+"-"+pad2(dtt.getUTCMonth()+1)+"-"+pad2(dtt.getUTCDate());
      if(ACT.step<86400) lab+=" "+pad2(dtt.getUTCHours())+":"+pad2(dtt.getUTCMinutes());
      var rows=lab, any=false;
      codes.forEach(function(c){ if(!state.sel[c]||!ACT.stn[c])return; var v=ACT.stn[c][state.comp][i];
        if(v!=null){ any=true; rows+='<br><span style="color:'+D.stn[c].c+'">●</span> '
          +D.stn[c].name+" "+v.toLocaleString()+" nT"; }});
      if(!any){tip.style.opacity=0;return;}
      tip.innerHTML=rows; tip.style.opacity=1;
      var tx=Math.min(mx+12, cv.clientWidth-tip.offsetWidth-4);
      tip.style.left=tx+"px"; tip.style.top=(e.clientY-r.top+12)+"px"; };
    cv.onmouseleave=function(){ tip.style.opacity=0; };

    $("svnote").innerHTML="Z·총자력 증가, Y 감소(편각 서편화)는 외핵 지오다이나모에 의한 <b>영년변화</b>다 — "
      +"지구 자기장이 매년 변하므로 반복측점의 주기적 재측정이 필요하다. '분 단위·월별'은 로컬 서버"
      +"(<code>python sv_server.py</code>) 실행 시 특정 월을 분 단위로 조회한다. 강릉·제주 X·Y 는 계기"
      +"(변화계) 좌표라 절대 편각 환산 금지(추세는 유효). 자료: 우주환경 빅데이터 플랫폼 · INTERMAGNET.";

    setActiveFull(); draw(); requestAnimationFrame(draw);
    window.addEventListener("resize", draw);
  })();
  </script>
  <!-- /SV_JS -->
"""

JS = JS.replace("__SVDATA__", SVDATA)

# ── 삽입/교체 (멱등) ────────────────────────────────────────────────────────
html = HTML.read_text(encoding="utf-8")
if "<!-- SV_CARD -->" in html:
    html = re.sub(r"  <!-- SV_CARD -->.*?  <!-- /SV_CARD -->\n", CARD, html, flags=re.S)
else:
    anchor = '  <div class="card">\n    <div class="note" id="acc"></div>'
    assert anchor in html, "카드 삽입 앵커 없음"
    html = html.replace(anchor, CARD + "\n" + anchor, 1)
if "<!-- SV_JS -->" in html:
    html = re.sub(r"  <!-- SV_JS -->.*?  <!-- /SV_JS -->\n", JS, html, flags=re.S)
else:
    html = html.replace("</body>", JS + "</body>", 1)

HTML.write_text(html, encoding="utf-8")
print(f"완료 · lmm.html 총 {len(html)//1024} KB")
