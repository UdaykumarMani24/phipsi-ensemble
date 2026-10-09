/* =========================================================================
   Phi/Psi Ensemble Analyzer - page logic
   Developed by Dr. Udayakumar Mani and Dr. Senthilkumar Rathnasamy
   ========================================================================= */
(function () {
  'use strict';
  const C = window.PhiPsiCore;
  const $ = id => document.getElementById(id);
  const PHI = 'φ', PSI = 'ψ', DEGS = '°', RHO = 'ρ', PM = '±', DASH = '—';
  const FLEX_CAP = 60;

  /* ---------- remember the last inputs (optional; ignored if storage is blocked) ---------- */
  const FIELDS = ['pdbid', 'chain', 'uniprot', 'offset'];
  try {
    const saved = JSON.parse(localStorage.getItem('phipsi-inputs') || '{}');
    FIELDS.forEach(f => { if (saved[f] !== undefined) $(f).value = saved[f]; });
  } catch (e) { /* storage unavailable */ }
  function saveInputs() {
    try { localStorage.setItem('phipsi-inputs', JSON.stringify(Object.fromEntries(FIELDS.map(f => [f, $(f).value])))); }
    catch (e) { /* storage unavailable */ }
  }

  /* ---------- messages (always plain text, never HTML) ---------- */
  function clearMessages() { $('messages').textContent = ''; }
  function message(text, kind = 'warn') {
    const d = document.createElement('div');
    d.className = 'msg ' + kind; d.textContent = text;
    $('messages').appendChild(d);
  }

  /* ---------- downloads ---------- */
  async function getText(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error('HTTP ' + r.status);
    return r.text();
  }
  const readFile = f => f.text();

  async function fetchPlddt(acc) {
    acc = acc.trim().toUpperCase();
    if (!/^[A-Z0-9]{6,10}$/.test(acc)) { message('Invalid UniProt accession: ' + acc); return {}; }
    const urls = [], notes = [];
    try {
      let data = JSON.parse(await getText('https://alphafold.ebi.ac.uk/api/prediction/' + acc));
      if (data && !Array.isArray(data)) data = [data];
      if (data && data[0] && data[0].pdbUrl) {
        urls.push(data[0].pdbUrl);
        if (data.length > 1) message('AlphaFold split this protein into fragments; only the first is used.');
      }
    } catch (e) { notes.push('API: ' + e.message); }
    ['v6', 'v4'].forEach(v => urls.push('https://alphafold.ebi.ac.uk/files/AF-' + acc + '-F1-model_' + v + '.pdb'));
    for (const u of [...new Set(urls)]) {
      try {
        const map = C.parsePlddt(await getText(u));
        if (Object.keys(map).length) return map;
        notes.push(u.split('/').pop() + ': no CA atoms');
      } catch (e) { notes.push(u.split('/').pop() + ': ' + e.message); }
    }
    message('AlphaFold download failed for ' + acc + ' (' + notes.join(' | ') +
            '). Download the PDB file from alphafold.ebi.ac.uk/entry/' + acc + ' and upload it instead.');
    return {};
  }

  /* ---------- run the analysis ---------- */
  let state = null;

  async function run() {
    clearMessages(); saveInputs();
    const btn = $('runBtn'); btn.disabled = true; btn.textContent = 'Working…';
    try {
      const chain = $('chain').value.trim().slice(0, 1);
      const offset = parseInt($('offset').value, 10) || 0;
      const texts = [];
      const file = $('pdbfile').files[0];
      if (file) {
        if (file.size > 60 * 1024 * 1024) { message('File too large (max 60 MB).', 'err'); return; }
        texts.push([file.name, await readFile(file)]);
      } else if ($('pdbid').value.trim()) {
        let ids = $('pdbid').value.trim().toUpperCase().split(/[\s,;]+/).filter(Boolean);
        ids = [...new Set(ids)].slice(0, 30);
        for (const id of ids) {
          if (!/^[0-9A-Z]{4}$/.test(id)) { message('Skipped invalid ID: ' + id); continue; }
          try { texts.push([id, await getText('https://files.rcsb.org/download/' + id + '.pdb')]); }
          catch (e) { message('Could not download ' + id + ' (' + e.message + '). Very large entries have no .pdb file.'); }
        }
      } else { message('Enter one or more PDB IDs, or upload a PDB file.', 'err'); return; }

      const allModels = [], labels = [];
      let viewText = null, source = '';
      texts.forEach(([label, t], k) => {
        if (k === 0) { viewText = t; source = label; }
        const ms = C.parseModels(t, chain);
        if (!ms.length) { message('No backbone atoms found in ' + label + (chain ? ' (chain ' + chain + ')' : '') + '.'); return; }
        ms.forEach((m, i) => { allModels.push(C.modelTorsions(m)); labels.push(label + (ms.length > 1 ? '#' + (i + 1) : '')); });
      });
      if (texts.length > 1) source += ' + ' + (texts.length - 1) + ' more';
      if (!allModels.length) { message('No usable protein backbone found.', 'err'); return; }

      let plddt = {};
      const af = $('affile').files[0];
      if (af) {
        plddt = C.parsePlddt(await readFile(af));
        if (!Object.keys(plddt).length) message('Could not read pLDDT from the uploaded AlphaFold file (needs PDB format).');
      } else if ($('uniprot').value.trim()) {
        plddt = await fetchPlddt($('uniprot').value);
      }

      const rows = C.aggregate(allModels, plddt, offset);
      let stats = null;
      if (Object.keys(plddt).length) {
        if (!rows.some(r => r.plddt !== null)) {
          message('pLDDT loaded but no residue numbers matched. Check the offset (UniProt position = PDB number + offset).');
        } else {
          const pairs = rows.filter(r => r.flex !== null && r.plddt !== null);
          stats = C.spearman(pairs.map(r => r.flex), pairs.map(r => r.plddt));
        }
      }
      state = { rows, labels, NM: labels.length, viewText, source, stats, multi: texts.length > 1, selected: null, sorted: false };
      show();
    } finally {
      btn.disabled = false; btn.textContent = 'Analyze';
    }
  }

  $('form').addEventListener('submit', e => { e.preventDefault(); run(); });
  $('exampleBtn').addEventListener('click', () => {
    $('pdbid').value = '1D3Z'; $('chain').value = 'A'; $('uniprot').value = 'P62987'; $('offset').value = '0';
    $('pdbfile').value = ''; $('affile').value = ''; run();
  });

  /* ---------- colours ---------- */
  function rgb(t) {
    const A = [59, 91, 219], B = [255, 212, 59], D = [224, 49, 49];
    const [p, q, u] = t < .5 ? [A, B, t / .5] : [B, D, (t - .5) / .5];
    return [0, 1, 2].map(i => Math.round(p[i] + (q[i] - p[i]) * u));
  }
  const isNum = v => v !== null && v !== undefined;
  const cssColor = f => (isNum(f) ? 'rgb(' + rgb(Math.min(1, f / FLEX_CAP)).join(',') + ')' : '#999');
  const hexColor = f => { if (!isNum(f)) return 0x999999; const c = rgb(Math.min(1, f / FLEX_CAP)); return (c[0] << 16) | (c[1] << 8) | c[2]; };
  const f1 = v => (isNum(v) ? v.toFixed(1) : DASH);
  const canvasXY = (cv, e) => { const b = cv.getBoundingClientRect(); return [(e.clientX - b.left) * cv.width / b.width, (e.clientY - b.top) * cv.height / b.height]; };

  /* ---------- render results ---------- */
  let viewer = null, nv = 1, flexMap = {};

  function show() {
    const S = state, R = S.rows;
    $('results').hidden = false;
    flexMap = {}; R.forEach(r => { flexMap[r.chain + '|' + r.resSeq] = r.flex; });
    $('info').textContent = S.source + ' ' + DASH + ' ' + S.NM + ' model(s), ' + R.length + ' residues';

    let sum;
    if (S.NM < 2) sum = 'Only 1 model found: spread cannot be computed. Use an NMR entry or several PDB IDs.';
    else {
      const top = R.filter(r => isNum(r.flex)).sort((a, b) => b.flex - a.flex).slice(0, 6)
        .map(r => r.name + r.resSeq + ' (' + r.flex.toFixed(0) + DEGS + ')').join(', ');
      sum = 'Most flexible: ' + (top || 'n/a');
      const two = R.filter(r => r.twoState).map(r => r.name + r.resSeq);
      if (two.length) sum += '  |  Two-state residues: ' + two.join(', ');
      if (S.stats) sum += '  |  Spearman ' + RHO + '(flexibility, pLDDT) = ' + S.stats.rho.toFixed(2) + ' (n = ' + S.stats.n +
        ', permutation p ' + (S.stats.p < 0.001 ? '< 0.001' : '= ' + S.stats.p.toFixed(3)) + '; approximate, neighbouring residues are not independent)';
    }
    $('summary').textContent = sum;
    $('viewerTitle').textContent = '3D structure coloured by flexibility' + (S.NM > 1 ? '' : ' (single model)');
    $('viewerNote').textContent = S.multi ? 'Several entries are not superimposed, so the 3D view shows the first entry only.' : '';
    $('scWrap').hidden = !R.some(r => isNum(r.plddt));
    $('plLegend').hidden = $('scWrap').hidden;

    $('viewer').textContent = '';
    viewer = $3Dmol.createViewer($('viewer'), { backgroundColor: 'white' });
    const vm = viewer.addModels(S.viewText, 'pdb');
    nv = Array.isArray(vm) ? vm.length : 1;
    viewer.zoomTo();
    S.sorted = false; $('sortBtn').textContent = 'Sort by flexibility';
    // taller, larger-text profile on phones
    if (window.innerWidth < 700) { pc.width = 600; pc.height = 340; } else { pc.width = 1000; pc.height = 240; }
    renderAll();
    $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function renderAll() { render3D(); buildTable(); drawProfile(); drawRama(); drawScatter(); }

  const selSpec = r => { const s = { resi: r.resSeq }; if (r.chain.trim()) s.chain = r.chain; return s; };
  function render3D() {
    if (!viewer) return;
    const first = $('firstOnly').checked, S = state;
    viewer.removeAllLabels();
    viewer.setStyle({}, { cartoon: { colorfunc: a => hexColor(flexMap[a.chain + '|' + a.resi]), opacity: first ? 1 : 0.45 } });
    if (first) for (let i = 1; i < nv; i++) viewer.setStyle({ model: i }, {});
    if (S.selected !== null) {
      const r = S.rows[S.selected];
      viewer.addStyle(first ? Object.assign({ model: 0 }, selSpec(r)) : selSpec(r), { stick: { radius: 0.14, colorscheme: 'greenCarbon' } });
      const ca = viewer.selectedAtoms(Object.assign({ atom: 'CA', model: 0 }, selSpec(r)))[0];
      if (ca) viewer.addLabel(r.name + r.resSeq, { position: ca, fontSize: 13, backgroundColor: 'white', backgroundOpacity: .8, fontColor: 'black', inFront: true });
      viewer.zoomTo(selSpec(r), 500);
    }
    viewer.render();
  }
  $('firstOnly').addEventListener('change', render3D);

  function select(i) {
    const S = state; S.selected = i;
    const r = S.rows[i];
    $('info').textContent = r.name + ' ' + r.resSeq + r.icode + ' (chain ' + (r.chain.trim() || '-') + ')   n = ' + r.n +
      '   ' + PHI + ' = ' + f1(r.phiMean) + DEGS + ' ' + PM + ' ' + f1(r.phiSD) +
      '   ' + PSI + ' = ' + f1(r.psiMean) + DEGS + ' ' + PM + ' ' + f1(r.psiSD) +
      (r.twoState ? '   [two states: ' + r.twoState.sizes.join(' + ') + ' models]' : '');
    renderAll();
    const tr = document.querySelector('#tbl tbody tr[data-i="' + i + '"]');
    if (tr) tr.scrollIntoView({ block: 'nearest' });
  }

  /* ---------- table ---------- */
  function buildTable() {
    const S = state, tb = document.querySelector('#tbl tbody');
    const order = S.rows.map((_, i) => i);
    if (S.sorted) order.sort((a, b) => (S.rows[b].flex ?? -1) - (S.rows[a].flex ?? -1));
    tb.textContent = '';
    order.forEach(i => {
      const r = S.rows[i], tr = document.createElement('tr');
      tr.dataset.i = i;
      if (i === S.selected) tr.className = 'sel';
      const cells = [r.chain.trim() || '-', r.name + ' ' + r.resSeq + r.icode, r.n,
        f1(r.phiMean), f1(r.phiSD), f1(r.psiMean), f1(r.psiSD), f1(r.flex), f1(r.plddt)];
      cells.forEach((v, j) => {
        const td = document.createElement('td'); td.textContent = v;
        if (j === 7) { td.style.color = cssColor(r.flex); td.style.fontWeight = '700'; }
        tr.appendChild(td);
      });
      const td = document.createElement('td');
      if (r.twoState) { const s = document.createElement('span'); s.className = 'tag'; s.textContent = '2-state'; td.appendChild(s); }
      tr.appendChild(td);
      tr.addEventListener('click', () => select(i));
      tb.appendChild(tr);
    });
  }
  $('sortBtn').addEventListener('click', e => {
    state.sorted = !state.sorted; e.target.textContent = state.sorted ? 'Sort by sequence' : 'Sort by flexibility'; buildTable();
  });

  /* ---------- flexibility profile ---------- */
  const pc = $('profile'), px = pc.getContext('2d');
  const PL = 46, PT = 12, PB = 28;
  const PR = () => (state.rows.some(r => isNum(r.plddt)) ? 44 : 12);
  const yMax = () => Math.max(FLEX_CAP, ...state.rows.map(r => r.flex ?? 0));
  const barW = () => (pc.width - PL - PR()) / state.rows.length;
  function drawProfile() {
    const R = state.rows, H = pc.height - PT - PB, bw = barW(), ym = yMax(), pr = PR();
    px.clearRect(0, 0, pc.width, pc.height);
    px.font = '13px sans-serif'; px.fillStyle = '#555'; px.strokeStyle = '#ddd'; px.lineWidth = 1;
    [0, ym / 2, ym].forEach(v => {
      const y = PT + H - v / ym * H;
      px.beginPath(); px.moveTo(PL, y); px.lineTo(pc.width - pr, y); px.stroke();
      px.fillText(Math.round(v) + DEGS, 6, y + 4);
    });
    R.forEach((r, i) => {
      const h = isNum(r.flex) ? Math.max(1, r.flex / ym * H) : 2;
      px.fillStyle = cssColor(r.flex);
      px.fillRect(PL + i * bw, PT + H - h, Math.max(1, bw - 0.5), h);
      if (r.twoState) { px.fillStyle = '#e8590c'; px.beginPath(); px.arc(PL + (i + .5) * bw, PT + H - h - 6, 3.5, 0, 7); px.fill(); }
    });
    if (R.some(r => isNum(r.plddt))) {
      px.strokeStyle = '#111'; px.lineWidth = 2; px.beginPath();
      let started = false;
      R.forEach((r, i) => {
        if (!isNum(r.plddt)) { started = false; return; }
        const x = PL + (i + .5) * bw, y = PT + H - r.plddt / 100 * H;
        if (!started) { px.moveTo(x, y); started = true; } else px.lineTo(x, y);
      });
      px.stroke();
      px.fillStyle = '#111';
      [0, 50, 100].forEach(v => px.fillText(String(v), pc.width - pr + 6, PT + H - v / 100 * H + 4));
    }
    px.fillStyle = '#555';
    px.fillText(String(R[0].resSeq), PL, pc.height - 8);
    px.fillText(String(R[R.length - 1].resSeq), pc.width - pr - 28, pc.height - 8);
    if (state.selected !== null) {
      px.strokeStyle = '#000'; px.lineWidth = 2;
      px.strokeRect(PL + state.selected * bw - 1, PT, Math.max(3, bw) + 1, H);
    }
  }
  pc.addEventListener('click', e => {
    const [x] = canvasXY(pc, e), i = Math.floor((x - PL) / barW());
    if (i >= 0 && i < state.rows.length) select(i);
  });

  /* ---------- Ramachandran ---------- */
  const rc = $('rama'), rx = rc.getContext('2d');
  const M = 38, RW = rc.width - M - 12, RH = rc.height - M - 12;
  const X = a => M + (a + 180) / 360 * RW, Y = a => 12 + (180 - a) / 360 * RH;
  function drawRama() {
    const S = state;
    rx.clearRect(0, 0, rc.width, rc.height);
    rx.strokeStyle = '#ccc'; rx.lineWidth = 1; rx.fillStyle = '#555'; rx.font = '12px sans-serif';
    for (let a = -180; a <= 180; a += 90) {
      rx.beginPath(); rx.moveTo(X(a), 12); rx.lineTo(X(a), 12 + RH); rx.stroke();
      rx.beginPath(); rx.moveTo(M, Y(a)); rx.lineTo(M + RW, Y(a)); rx.stroke();
      rx.fillText(String(a), X(a) - 12, 12 + RH + 16); rx.fillText(String(a), 2, Y(a) + 4);
    }
    rx.strokeStyle = '#333'; rx.strokeRect(M, 12, RW, RH);
    rx.fillStyle = '#000'; rx.font = 'italic 15px serif';
    rx.fillText(PHI, M + RW / 2, rc.height - 2); rx.fillText(PSI, 4, 12 + RH / 2 - 10);
    rx.fillStyle = 'rgba(120,120,120,.35)';
    S.rows.forEach(r => { for (let m = 0; m < S.NM; m++) if (isNum(r.phi[m]) && isNum(r.psi[m])) rx.fillRect(X(r.phi[m]) - 1, Y(r.psi[m]) - 1, 2, 2); });
    if (S.selected === null) return;
    const r = S.rows[S.selected];
    for (let m = 0; m < S.NM; m++) {
      if (!isNum(r.phi[m]) || !isNum(r.psi[m])) continue;
      rx.beginPath(); rx.arc(X(r.phi[m]), Y(r.psi[m]), 4.5, 0, 7);
      rx.fillStyle = 'hsla(' + (m * 360 / Math.max(S.NM, 1)) + ',75%,50%,.85)'; rx.fill();
      rx.strokeStyle = '#fff'; rx.lineWidth = 1; rx.stroke();
    }
    const cross = (cx, cy, col) => { rx.strokeStyle = col; rx.lineWidth = 2.5; rx.beginPath(); rx.moveTo(cx - 8, cy); rx.lineTo(cx + 8, cy); rx.moveTo(cx, cy - 8); rx.lineTo(cx, cy + 8); rx.stroke(); };
    if (r.twoState) r.twoState.centres.forEach(c => cross(X(c[0]), Y(c[1]), '#e8590c'));
    else if (isNum(r.phiMean) && isNum(r.psiMean)) cross(X(r.phiMean), Y(r.psiMean), '#000');
  }
  rc.addEventListener('click', e => {
    const S = state, [mx, my] = canvasXY(rc, e);
    let best = -1, bd = 14;
    S.rows.forEach((r, i) => { for (let m = 0; m < S.NM; m++) {
      if (!isNum(r.phi[m]) || !isNum(r.psi[m])) continue;
      const d = Math.hypot(X(r.phi[m]) - mx, Y(r.psi[m]) - my);
      if (d < bd) { bd = d; best = i; }
    }});
    if (best >= 0) select(best);
  });

  /* ---------- flexibility vs pLDDT ---------- */
  const sc = $('scatter'), sx = sc.getContext('2d');
  const SM = 44, SW = sc.width - SM - 14, SH = sc.height - SM - 14;
  const SX = v => SM + Math.min(v, yMax()) / yMax() * SW, SY = v => 14 + (100 - v) / 100 * SH;
  function drawScatter() {
    const S = state;
    if (!S.rows.some(r => isNum(r.plddt))) return;
    const ym = yMax();
    sx.clearRect(0, 0, sc.width, sc.height);
    sx.strokeStyle = '#ddd'; sx.lineWidth = 1; sx.fillStyle = '#555'; sx.font = '12px sans-serif';
    [0, 50, 100].forEach(v => { sx.beginPath(); sx.moveTo(SM, SY(v)); sx.lineTo(SM + SW, SY(v)); sx.stroke(); sx.fillText(String(v), 10, SY(v) + 4); });
    [0, ym / 2, ym].forEach(v => { sx.beginPath(); sx.moveTo(SX(v), 14); sx.lineTo(SX(v), 14 + SH); sx.stroke(); sx.fillText(String(Math.round(v)), SX(v) - 8, 14 + SH + 16); });
    sx.strokeStyle = '#333'; sx.strokeRect(SM, 14, SW, SH);
    sx.fillStyle = '#000'; sx.fillText('Flexibility (' + DEGS + ')', SM + SW / 2 - 34, sc.height - 4); sx.save(); sx.translate(12, 14 + SH / 2 + 18); sx.rotate(-Math.PI / 2); sx.fillText('pLDDT', 0, 0); sx.restore();
    S.rows.forEach((r, i) => {
      if (!isNum(r.flex) || !isNum(r.plddt)) return;
      const sel = i === S.selected;
      sx.beginPath(); sx.arc(SX(r.flex), SY(r.plddt), sel ? 7 : 3.5, 0, 7);
      sx.fillStyle = cssColor(r.flex); sx.fill();
      sx.strokeStyle = sel ? '#000' : 'rgba(0,0,0,.25)'; sx.lineWidth = sel ? 2.5 : 1; sx.stroke();
    });
  }
  sc.addEventListener('click', e => {
    const [mx, my] = canvasXY(sc, e);
    let best = -1, bd = 14;
    state.rows.forEach((r, i) => {
      if (!isNum(r.flex) || !isNum(r.plddt)) return;
      const d = Math.hypot(SX(r.flex) - mx, SY(r.plddt) - my);
      if (d < bd) { bd = d; best = i; }
    });
    if (best >= 0) select(best);
  });

  /* ---------- CSV export ---------- */
  function download(name, text) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type: 'text/csv' }));
    a.download = name; document.body.appendChild(a); a.click(); a.remove();
  }
  const cell = v => (isNum(v) ? v : '');
  $('csv1').addEventListener('click', () => {
    let s = 'chain,resSeq,icode,resName,n_models,phi_mean,phi_circ_sd,psi_mean,psi_circ_sd,flexibility_deg,alphafold_plddt,two_state\n';
    state.rows.forEach(r => { s += [r.chain.trim(), r.resSeq, r.icode, r.name, r.n, cell(r.phiMean), cell(r.phiSD),
      cell(r.psiMean), cell(r.psiSD), cell(r.flex), cell(r.plddt), r.twoState ? 'yes' : ''].join(',') + '\n'; });
    download('phi_psi_summary.csv', s);
  });
  $('csv2').addEventListener('click', () => {
    let s = 'model,chain,resSeq,icode,resName,phi,psi\n';
    for (let m = 0; m < state.NM; m++) state.rows.forEach(r => {
      s += [state.labels[m], r.chain.trim(), r.resSeq, r.icode, r.name, cell(r.phi[m]), cell(r.psi[m])].join(',') + '\n';
    });
    download('phi_psi_per_model.csv', s);
  });

  /* ---------- installable app ---------- */
  let installEvt = null;
  window.addEventListener('beforeinstallprompt', e => { e.preventDefault(); installEvt = e; $('installBtn').hidden = false; });
  $('installBtn').addEventListener('click', async () => {
    if (!installEvt) return;
    installEvt.prompt(); await installEvt.userChoice; installEvt = null; $('installBtn').hidden = true;
  });
  if ('serviceWorker' in navigator && (location.protocol === 'https:' || location.hostname === 'localhost')) {
    navigator.serviceWorker.register('sw.js').catch(() => {});
  }
})();
