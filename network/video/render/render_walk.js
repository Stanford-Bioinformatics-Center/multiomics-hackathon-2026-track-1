// =====================================================================================================
// render/render_walk.js — the fly-through: the figure 17 page, node by node along the walk, one JPEG per frame
// =====================================================================================================
// Called by exvideo/render.py with a spec file (render_spec.json): the page, size, fps, song length, walk, per-node
// segments (label, persona, fact lines) and per-bar timing. For every frame (t = frame / fps) it
//   - sets the camera: the whole network during the intro; then a fly from node to node (ease in / out, with a
//     slight zoom-out mid-flight) and a slow zoom-in while the node's bars play; back to the whole walk at the end;
//   - on entering a node: focuses it (its neighbours stay bright), switches on the page's "arm-specific edges"
//     colouring (red = strong after endurance only, blue = resistance only) and paints the walked edges gold;
//   - updates the overlay: title, walk breadcrumbs, the node's persona and fact card, the current bar and the next;
//   - writes the frame as a JPEG. Frame-stepped (not screen-recorded), so no frames are dropped and reruns match.
// =====================================================================================================
const puppeteer = require('puppeteer');
const fs = require('fs');
const path = require('path');

const spec = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const wait = ms => new Promise(r => setTimeout(r, ms));

(async () => {
  const browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--hide-scrollbars'] });
  const page = await browser.newPage();
  await page.setViewport({ width: spec.width, height: spec.height, deviceScaleFactor: 1 });
  page.on('pageerror', e => console.error('page error:', e.message));
  await page.goto('file://' + spec.html, { waitUntil: 'load' });
  await wait(1500);

  // ---- stage set-up: the network full screen, arm-specific edges on, the overlay added ------------------------------
  const ok = await page.evaluate((S) => {
    const g = document.querySelector('[id^=graph]'); if (!g || !g.chart) return 'no network widget on the page';
    const net = g.chart, missing = S.walk.filter(n => net.body.data.nodes.get(n) === null);
    if (missing.length) return 'walk nodes not on the page: ' + missing.join(', ');
    document.body.style.cssText += ';margin:0;overflow:hidden;background:#fff';
    const host = g.closest('.html-widget') || g.parentElement;
    host.style.cssText += `;position:fixed;left:0;top:0;width:${S.width}px;height:${S.height}px;z-index:5;background:#fff`;
    [g, g.parentElement].forEach(el => { el.style.width = S.width + 'px'; el.style.height = S.height + 'px'; });
    net.setSize(S.width + 'px', S.height + 'px'); net.redraw();
    // hide the page's own titles, controls and legend: every element that is neither the network nor its ancestor
    document.querySelectorAll('body *').forEach(el => { if (!host.contains(el) && !el.contains(host)) el.style.visibility = 'hidden'; });
    // ... and everything next to the network at every level up to the page (the widget's own title / subtitle)
    for (let el = g; el && el !== document.body; el = el.parentElement)
      [...(el.parentElement ? el.parentElement.children : [])].forEach(c => { if (c !== el && c.id !== 'mv') c.style.display = 'none'; });
    net.redraw();
    const sp = document.querySelector('.hk-spec'); if (sp) { sp.checked = true; sp.dispatchEvent(new Event('change')); }
    net.fit({ animation: false });
    const k = S.height / 720, css = `
      #mv { position: fixed; inset: 0; z-index: 10; pointer-events: none; font-family: Helvetica, Arial, sans-serif; }
      #mv-top { position: absolute; left: ${24 * k}px; right: ${24 * k}px; top: ${16 * k}px; display: flex; justify-content: space-between; align-items: baseline; }
      #mv-brand { font-weight: 800; font-size: ${22 * k}px; letter-spacing: ${2 * k}px; color: #111; }
      #mv-brand span { font-weight: 400; color: #555; letter-spacing: 0; margin-left: ${10 * k}px; font-size: ${18 * k}px; }
      #mv-crumbs { font-size: ${16 * k}px; color: #999; } #mv-crumbs b { color: #111; } #mv-crumbs i { color: #C98A00; font-style: normal; font-weight: 700; }
      #mv-card { position: absolute; left: ${24 * k}px; top: ${64 * k}px; background: rgba(255,255,255,0.92); border-left: ${5 * k}px solid #C98A00;
                 padding: ${10 * k}px ${14 * k}px; max-width: ${430 * k}px; box-shadow: 0 ${2 * k}px ${10 * k}px rgba(0,0,0,0.12); }
      #mv-node { font-size: ${30 * k}px; font-weight: 800; color: #111; } #mv-persona { font-size: ${19 * k}px; color: #C98A00; font-style: italic; margin-bottom: ${6 * k}px; }
      #mv-facts div { font-size: ${14 * k}px; color: #333; line-height: 1.45; }
      #mv-bars { position: absolute; left: 0; right: 0; bottom: 0; background: rgba(12,12,12,0.80); padding: ${14 * k}px ${40 * k}px ${16 * k}px; text-align: center; }
      #mv-now { color: #fff; font-size: ${27 * k}px; font-weight: 700; line-height: 1.25; } #mv-next { color: rgba(255,255,255,0.55); font-size: ${18 * k}px; margin-top: ${6 * k}px; }
      .vis-navigation, .vis-button { display: none !important; }
      #mv-center { position: absolute; left: 50%; top: 36%; transform: translateX(-50%); text-align: center; background: rgba(255,255,255,0.88);
                   padding: ${16 * k}px ${36 * k}px; box-shadow: 0 ${2 * k}px ${14 * k}px rgba(0,0,0,0.15); white-space: nowrap; }
      #mv-center .t { font-size: ${48 * k}px; font-weight: 800; color: #111; text-shadow: 0 0 ${12 * k}px #fff; } #mv-center .s { font-size: ${22 * k}px; color: #555; margin-top: ${10 * k}px; }`;
    const st = document.createElement('style'); st.textContent = css; document.head.appendChild(st);
    const ov = document.createElement('div'); ov.id = 'mv';
    ov.innerHTML = `<div id="mv-top"><div id="mv-brand">TEAM 2-PAC<span></span></div><div id="mv-crumbs"></div></div>
      <div id="mv-card"><div id="mv-node"></div><div id="mv-persona"></div><div id="mv-facts"></div></div>
      <div id="mv-center"><div class="t"></div><div class="s"></div></div><div id="mv-bars"><div id="mv-now"></div><div id="mv-next"></div></div>`;
    document.body.appendChild(ov);
    ov.querySelector('#mv-brand span').textContent = S.title;
    // everything the frame function needs
    const P = net.getPositions(S.walk), home = { pos: net.getViewPosition(), scale: net.getScale() };
    const zoom = Math.max(home.scale * 3.0, 0.9);
    const key = (a, b) => a < b ? a + '|' + b : b + '|' + a;
    const walkEdge = {}; net.body.data.edges.get().forEach(e => { walkEdge[key(e.from, e.to)] = e.id; });
    window.mv = { S, net, P, home, zoom, key, walkEdge, seg: -2 };
    const lerp = (a, b, u) => a + (b - a) * u, ease = u => u < 0.5 ? 2 * u * u : 1 - Math.pow(-2 * u + 2, 2) / 2, clamp = u => Math.max(0, Math.min(1, u));
    // paint the walk's edges up to (and including) step i gold; focus node i
    function enter(i) {
      const M = window.mv; if (M.seg === i) return; M.seg = i;
      if (window.hkState && window.hkApply) { window.hkState.focus = i >= 0 && i < S.walk.length ? [S.walk[i]] : null; window.hkApply(); }
      const upto = i < 0 ? 0 : Math.min(i, S.walk.length - 1);
      const ups = []; for (let j = 1; j <= upto; j++) { const id = M.walkEdge[M.key(S.walk[j - 1], S.walk[j])]; if (id !== undefined) ups.push({ id, width: 10, hidden: false, color: { color: '#F2A900', highlight: '#F2A900', opacity: 1 } }); }
      if (i >= S.walk.length) for (let j = 1; j < S.walk.length; j++) { const id = M.walkEdge[M.key(S.walk[j - 1], S.walk[j])]; if (id !== undefined) ups.push({ id, width: 10, hidden: false, color: { color: '#F2A900', highlight: '#F2A900', opacity: 1 } }); }
      if (ups.length) M.net.body.data.edges.update(ups);
      if (i >= 0 && i < S.walk.length) M.net.selectNodes([S.walk[i]]); else M.net.unselectAll();
    }
    // the state of the video at time t (seconds)
    window.mvFrame = function (t) {
      const M = window.mv, segs = S.segments, n = S.walk.length, first = segs[0].start, last = segs[n - 1].end;
      let cam, i;
      if (t < first) { i = -1; const u = clamp(t / Math.max(first, 0.1)); cam = { pos: M.home.pos, scale: M.home.scale * (1 + 0.06 * u) }; }
      else if (t >= last) { i = n; const u = ease(clamp((t - last) / 1.5)), from = M.P[S.walk[n - 1]];
        cam = { pos: { x: lerp(from.x, M.home.pos.x, u), y: lerp(from.y, M.home.pos.y, u) }, scale: lerp(M.zoom, M.home.scale, u) }; }
      else { i = segs.findIndex(s => t >= s.start && t < s.end); if (i < 0) i = n - 1;
        const s = segs[i], to = M.P[S.walk[i]], from = i === 0 ? { x: M.home.pos.x, y: M.home.pos.y } : M.P[S.walk[i - 1]];
        const fromScale = i === 0 ? M.home.scale : M.zoom, trans = Math.min(1.6, 0.35 * (s.end - s.start));
        const u = ease(clamp((t - s.start) / trans)), after = clamp((t - s.start - trans) / Math.max(s.end - s.start - trans, 0.1));
        const dip = i === 0 ? 0 : 0.35 * Math.sin(Math.PI * u);
        cam = { pos: { x: lerp(from.x, to.x, u), y: lerp(from.y, to.y, u) }, scale: lerp(fromScale, M.zoom, u) * (1 - dip) * (1 + 0.10 * after) }; }
      enter(i);
      M.net.moveTo({ position: cam.pos, scale: cam.scale, animation: false });
      // overlay text
      const crumbs = S.walk.map((w, j) => j === i ? `<i>${w}</i>` : j < i ? `<b>${w}</b>` : w).join(' &rarr; ');
      document.getElementById('mv-crumbs').innerHTML = crumbs;
      const card = document.getElementById('mv-card'), center = document.getElementById('mv-center'), barsBox = document.getElementById('mv-bars');
      if (i >= 0 && i < n) { const sg = segs[i]; card.style.display = 'block'; center.style.display = 'none';
        document.getElementById('mv-node').textContent = sg.label; document.getElementById('mv-persona').textContent = sg.persona;
        document.getElementById('mv-facts').innerHTML = sg.facts.map(f => `<div>${f}</div>`).join('');
      } else { card.style.display = 'none'; center.style.display = 'block';
        center.querySelector('.t').textContent = i < 0 ? S.title : S.walk.join(' → ');
        center.querySelector('.s').textContent = i < 0 ? S.subtitle : 'TEAM 2-PAC · physical links, exercise-weighted'; }
      const b = S.bars.findIndex(x => t >= x.start && t < x.end);
      barsBox.style.display = b >= 0 ? 'block' : 'none';
      if (b >= 0) { const now = document.getElementById('mv-now'); now.textContent = S.bars[b].text; now.style.color = S.bars[b].header ? '#F2A900' : '#fff';
        document.getElementById('mv-next').textContent = b + 1 < S.bars.length ? S.bars[b + 1].text : ''; }
      return new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
    };
    return 'ok';
  }, spec);
  if (ok !== 'ok') { console.error('renderer set-up failed: ' + ok); await browser.close(); process.exit(1); }

  // ---- frames ------------------------------------------------------------------------------------------------------
  const nFrames = Math.ceil(spec.duration * spec.fps); const t0 = Date.now();
  for (let f = 0; f < nFrames; f++) {
    await page.evaluate(t => window.mvFrame(t), f / spec.fps);
    await page.screenshot({ path: path.join(spec.frames_dir, String(f).padStart(5, '0') + '.jpg'), type: 'jpeg', quality: 88, optimizeForSpeed: true });
    if (f % Math.max(1, Math.round(nFrames / 10)) === 0) console.log(`  frame ${f}/${nFrames} (${((Date.now() - t0) / 1000).toFixed(0)} s)`);
  }
  console.log(`  ${nFrames} frames in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
