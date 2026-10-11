/* SeptaTyper front end. No framework. All server text goes through esc(). */
"use strict";
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = n => Number(n).toLocaleString("en-US");
let META = null;
const SP = {};                       // key -> display name
const SHORT = { ecoli: "E. coli", ssonnei: "S. sonnei", senterica: "S. enterica", kpneumoniae: "K. pneumoniae" };
const SITE = { FECAL: "Fecal", URINE: "Urine", RESPIRATORY: "Respiratory", ALL: "All sites / unknown" };
const POS = [["gyrA", "GyrA-83"], ["gyrA", "GyrA-87"], ["parC", "ParC-80"], ["parC", "ParC-84"], ["recN", "RecN-346"], ["qnr", "qnrB"], ["qnr", "qnrS"]];
const WILD = ["S", "D", "S", "E", "L", "0", "0"];
const LOADED = {};

async function api(url, opts) {
  let r, j = null;
  try { r = await fetch(url, opts); } catch (e) { throw new Error("Cannot reach the server. Check your connection."); }
  try { j = await r.json(); } catch (e) { /* not json */ }
  if (!r.ok) {
    const d = j && j.detail;
    throw new Error(typeof d === "string" ? d : Array.isArray(d) ? d.map(x => x.msg).join("; ") : `Request failed (${r.status})`);
  }
  return j;
}
const callName = c => c === "Resistant" ? "Resistant (Cres)" : c === "Susceptible" ? "Susceptible (Csus)" : c;
const badge = (c, sm) => `<span class="badge ${esc(c)} ${sm ? "sm" : ""}">${esc(sm ? c : callName(c))}</span>`;
const bar = (r, s) => { const t = r + s || 1, p = Math.round(100 * r / t); return `<div class="bar" title="${r} resistant, ${s} susceptible"><i class="r" style="width:${p}%"></i><i class="s" style="width:${100 - p}%"></i></div>`; };
const mutSep = sep => (sep || "").split("-").map((v, i) => i < 5 ? (v !== WILD[i] ? `<span class="mut">${esc(v)}</span>` : esc(v)) : (v === "1" ? `<span class="mut">1</span>` : esc(v))).join("-");
const opts = (arr, first) => (first ? `<option value="">${first}</option>` : "") + arr.map(([v, t]) => `<option value="${esc(v)}">${esc(t)}</option>`).join("");
function setSites(sel, sp, withAll) {
  const list = (META.habitats[sp] || []).filter(h => withAll || h !== "ALL");
  sel.innerHTML = opts(list.map(h => [h, SITE[h]]));
  if (list.includes("ALL") && !withAll) sel.value = list[0];
}

/* ------------------------------------------------------------------ router */
const PAGES = ["home", "predict", "browse", "analytics", "submit", "theory", "data", "about", "curate"];
function route() {
  const p = (location.hash.replace(/^#\/?/, "") || "home").split("?")[0];
  const page = PAGES.includes(p) ? p : "home";
  $$(".page").forEach(e => e.classList.toggle("on", e.id === "p-" + page));
  $$("#nav a").forEach(a => a.classList.toggle("on", a.dataset.p === page));
  window.scrollTo(0, 0);
  if (META && !LOADED[page]) { LOADED[page] = true; (INIT[page] || (() => {}))(); }
}
window.addEventListener("hashchange", route);

/* -------------------------------------------------------------------- home */
let demoState = 0, demoData = null;
const DEMO = [["S", "D", "S", "E", "L", "0", "0"], ["L", "N", "I", "E", "L", "0", "0"]];
function drawDemo(animate) {
  const cur = DEMO[demoState], prev = DEMO[1 - demoState];
  $("#demo-tiles").innerHTML = cur.map((v, i) => {
    const changed = v !== WILD[i];
    return `<div class="tile ${changed ? "changed" : ""} ${animate && v !== prev[i] ? "flip" : ""}"><div class="cell">${esc(v)}</div><small>${POS[i][1]}</small></div>`;
  }).join("");
  const sep = cur.join("-"), d = demoData && demoData[sep];
  const call = demoState ? "Resistant" : "Susceptible";
  $("#demo-badge").className = "badge " + call; $("#demo-badge").textContent = callName(call);
  $("#demo-note").textContent = d ? `${sep}: ${d.n_resistant} of ${d.n_isolates} reference isolates resistant` : sep;
}
const INIT = {};
INIT.home = async () => {
  const c = META.counts;
  $("#c-genomes").textContent = fmt(c.genomes); $("#c-species").textContent = c.species;
  $("#c-countries").textContent = c.countries; $("#c-septa").textContent = c.septatypes;
  try {
    const r = await api("/api/v1/rules?species=ecoli&habitat=ALL");
    demoData = Object.fromEntries(r.rows.map(x => [x.septatype, x]));
  } catch (e) { /* demo still works without counts */ }
  drawDemo(false);
  const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (!reduce) setTimeout(() => { demoState = 1; drawDemo(true); }, 900);
  $("#demo-toggle").onclick = () => { demoState = 1 - demoState; drawDemo(!reduce); };
  const rows = META.stats.filter(s => s.habitat !== "ALL");
  $("#home-stats").innerHTML = `<tr><th>Species</th><th>Site</th><th>Isolates</th><th>Resistant</th><th>Septatypes</th><th>Accuracy</th></tr>` +
    rows.map(s => `<tr><td><i>${esc(SP[s.species])}</i></td><td>${esc(SITE[s.habitat])}</td><td class="stat">${fmt(s.n_isolates)}</td><td class="stat">${fmt(s.n_resistant)}</td><td class="stat">${s.n_septatypes}</td><td class="stat">${s.accuracy}%</td></tr>`).join("");
};

/* ----------------------------------------------------------------- predict */
let PMODE = "file", PFILE = null;
function renderResult(r, extra) {
  const p = r.prediction, a = r.analysis;
  let h = `<div class="band result">`;
  if (!p) {
    h += `<h3>No call could be made</h3><p>A complete septatype could not be built from this input.</p>`;
  } else {
    const tier = p.tier === "database" ? `Majority rule from ${fmt(p.n_isolates)} reference isolates`
      : p.tier === "database_other_site" ? `Majority rule from ${fmt(p.n_isolates)} isolates at other sites`
      : "Physicochemical grammar (new septatype)";
    h += `<div class="verdict-row">${badge(p.call)}<span class="muted">${esc(tier)}</span></div>`;
    const codon = (a && a.recN_codon) || r.recn_codon;
    h += `<div class="sep">${p.septatype.split("-").map((v, i) => {
      const ch = v !== WILD[i] && !(i >= 5 && v === "0");
      return `<div class="tile ${ch ? "changed" : ""}"><div class="cell">${esc(v)}</div><small>${POS[i][1]}</small>${i === 4 && codon ? `<span class="codon">${esc(codon.toLowerCase())}</span>` : ""}</div>`;
    }).join("")}</div>`;
    if (r.species === "senterica" && (a && a.recN_syn === "ctg" || String(codon).toUpperCase() === "CTG"))
      h += `<div class="info">RecN-346 uses the ctg codon for leucine, a high-usage synonymous change flagged in <i>S. enterica</i>.</div>`;
    if (p.known) {
      h += `<p class="stat">${fmt(p.n_resistant)} resistant and ${fmt(p.n_susceptible)} susceptible among reference isolates with this septatype.</p>${bar(p.n_resistant, p.n_susceptible)}`;
    } else if (p.similar && p.similar.length) {
      h += `<p class="muted" style="margin-top:12px">Known septatypes with the same GyrA/ParC tetratype:</p><div class="scroll"><table><tr><th>Septatype</th><th>Isolates</th><th>Resistant</th><th>Susceptible</th></tr>${p.similar.map(x => `<tr><td class="mono">${esc(x.septatype)}</td><td>${x.n_isolates}</td><td>${x.n_resistant}</td><td>${x.n_susceptible}</td></tr>`).join("")}</table></div>`;
    }
    h += `<h3>Why</h3><p class="stat">Combined mean hydrophobicity &lt;H&gt; of GyrA-83, GyrA-87, ParC-80 and ParC-84: <b>${p.mean_h}</b> (Eisenberg scale, cutoff 0.0).</p><ul class="why">` +
      p.mechanism.map(m => `<li><span class="dot ${esc(m.effect)}"></span><span><b>${esc(m.position)} ${esc(m.residue)}</b>. ${esc(m.text)}</span></li>`).join("") + `</ul>`;
    h += `<details><summary>Grammar votes: ${esc(p.grammar.call)} (${esc(p.grammar.rule)})</summary><div class="scroll"><table><tr><th>Position</th><th>Residue</th><th>Votes Resistant</th><th>Reason</th></tr>${p.grammar.votes.map(v => `<tr><td>${esc(v.position)}</td><td class="mono">${esc(v.residue)}</td><td>${v.vote ? "Yes" : "No"}</td><td>${esc(v.why)}</td></tr>`).join("")}</table></div></details>`;
  }
  if (r.warnings && r.warnings.length) h += `<div class="note"><ul>${r.warnings.map(w => `<li>${esc(w)}</li>`).join("")}</ul></div>`;
  if (a && a.positions) {
    h += `<details><summary>Sequence details</summary><div class="scroll"><table><tr><th>Position</th><th>Codon</th><th>Residue</th><th>Reference</th></tr>` +
      Object.values(a.positions).map(x => `<tr><td>${esc(x.gene)}-${x.pos}</td><td class="mono">${esc((x.codon || "-").toLowerCase())}</td><td class="mono">${esc(x.aa || "-")}</td><td class="mono">${esc(x.ref_aa || "-")}</td></tr>`).join("") +
      `</table></div><div class="scroll" style="margin-top:10px"><table><tr><th>Gene</th><th>Identity</th><th>Coverage</th><th>Contig</th><th>Strand</th></tr>` +
      Object.entries(a.genes).map(([g, x]) => `<tr><td>${esc(g)}</td><td>${x.found ? x.identity + "%" : "not found"}</td><td>${x.found ? x.coverage + "%" : "-"}</td><td>${esc(x.contig || "-")}</td><td>${esc(x.strand || "-")}</td></tr>`).join("") +
      `</table></div><div class="scroll" style="margin-top:10px"><table><tr><th>qnr family</th><th>Status</th><th>Basis</th><th>Best hit</th></tr>` +
      Object.entries(a.qnr).map(([k, x]) => `<tr><td>${k === "qnrB2" ? "qnrB" : "qnrS"}</td><td>${x.present ? "Present" : "Absent"}</td><td>${esc(x.source)}</td><td>${esc(x.best_hit || "-")}</td></tr>`).join("") +
      `</table></div><p class="muted" style="font-size:.85rem;margin-top:8px">${a.input.records} sequence record(s), ${fmt(a.input.total_bp)} bp.</p></details>`;
  }
  if (p && p.novel) h += `<div class="info">This septatype is new for ${esc(SHORT[r.species])}. If you have the laboratory phenotype, <a href="#/submit">submit the isolate</a> so curators can add it.</div>`;
  h += (extra || "") + `<p class="muted" style="font-size:.8rem">Research use only. Confirm with laboratory testing.</p></div>`;
  return h;
}
INIT.predict = () => {
  $("#sp").innerHTML = opts(Object.entries(SP));
  const upd = () => {
    setSites($("#hb"), $("#sp").value, true);
    const ref = META.references[$("#sp").value];
    $("#refwarn").innerHTML = ref && !ref.ready ? `<div class="note">Sequence analysis is not available yet because reference sequences are missing on the server. Enter residues works now.</div>` : "";
  };
  $("#sp").onchange = upd; upd();
  $$("#pmode button").forEach(b => b.onclick = () => {
    PMODE = b.dataset.m;
    $$("#pmode button").forEach(x => x.classList.toggle("on", x === b));
    $("#m-file").style.display = PMODE === "file" ? "" : "none";
    $("#m-paste").style.display = PMODE === "paste" ? "" : "none";
    $("#m-gene").style.display = PMODE === "manual" ? "none" : "";
    $("#m-manual").style.display = PMODE === "manual" ? "" : "none";
  });
  const drop = $("#drop"), fi = $("#file");
  const setF = f => { PFILE = f || null; $("#fname").textContent = f ? `${f.name} (${(f.size / 1048576).toFixed(1)} MB)` : ".fa .fasta .fna .gz, up to " + META.max_upload_mb + " MB"; };
  drop.onclick = () => fi.click();
  drop.onkeydown = e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); fi.click(); } };
  fi.onchange = () => setF(fi.files[0]);
  ["dragover", "dragenter"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove("over"); }));
  drop.addEventListener("drop", e => setF(e.dataTransfer.files[0]));
  $("#ex").onclick = () => { ["L", "N", "I", "E", "L"].forEach((v, i) => $("#" + ["r83", "r87", "r80", "r84", "r346"][i]).value = v); $("#mqb").checked = false; $("#mqs").checked = false; };
  $("#pf").onsubmit = async e => {
    e.preventDefault();
    const err = $("#perr"), btn = $("#go"); err.innerHTML = "";
    btn.disabled = true; const old = btn.textContent; btn.textContent = "Analysing";
    $("#presult").innerHTML = `<div class="band muted">${PMODE === "manual" ? "Looking up the septatype" : "Reading the sequence. A full genome takes up to a minute."}</div>`;
    try {
      let res;
      const species = $("#sp").value, habitat = $("#hb").value;
      if (PMODE === "manual") {
        const g = id => $(id).value.trim().toUpperCase();
        res = await api("/api/v1/predict/manual", { method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ species, habitat, gyrA_83: g("#r83"), gyrA_87: g("#r87"), parC_80: g("#r80"), parC_84: g("#r84"), recN_346: g("#r346"), qnrB2: $("#mqb").checked, qnrS1: $("#mqs").checked, recN_codon: $("#mc").value }) });
      } else {
        const fd = new FormData();
        fd.append("species", species); fd.append("habitat", habitat); fd.append("qnrB2", $("#qb").value); fd.append("qnrS1", $("#qs").value);
        if (PMODE === "file") { if (!PFILE) throw new Error("Choose a FASTA file first."); if (PFILE.size > META.max_upload_mb * 1048576) throw new Error(`File is larger than ${META.max_upload_mb} MB.`); fd.append("file", PFILE); }
        else { if (!$("#seq").value.trim()) throw new Error("Paste a nucleotide sequence first."); fd.append("sequence", $("#seq").value); }
        res = await api("/api/v1/predict", { method: "POST", body: fd });
      }
      $("#presult").innerHTML = renderResult(res);
      if (innerWidth < 860) $("#presult").scrollIntoView({ behavior: "smooth" });
    } catch (ex) { $("#presult").innerHTML = ""; err.innerHTML = `<div class="err">${esc(ex.message)}</div>`; }
    finally { btn.disabled = false; btn.textContent = old; }
  };
};

/* ------------------------------------------------------------------ browse */
let BS = { page: 1, sort: "genome_id", desc: false }, BTOT = 0;
const PER = 25;
function bfilters() {
  const f = { species: $("#f-sp").value, habitat: $("#f-hb").value, st: $("#f-st").value, country: $("#f-country").value, septatype: $("#f-sep").value.trim().toUpperCase(),
    phenotype: $("#f-ph").value, year_from: $("#f-y1").value, year_to: $("#f-y2").value, q: $("#f-q").value.trim(), source: $("#f-src").value };
  if ($("#f-wins").checked) f.wins = "true";
  const p = new URLSearchParams(); Object.entries(f).forEach(([k, v]) => { if (v) p.set(k, v); });
  return p;
}
async function loadBrowse() {
  const p = bfilters(), qs = p.toString();
  $("#dl-csv").href = "/api/v1/download/csv?" + qs; $("#dl-json").href = "/api/v1/download/json?" + qs;
  p.set("page", BS.page); p.set("per_page", PER); p.set("sort", BS.sort); if (BS.desc) p.set("desc", "true");
  try {
    const d = await api("/api/v1/isolates?" + p); BTOT = d.total;
    $("#b-count").textContent = `${fmt(d.total)} isolates`;
    $("#b-page").textContent = `Page ${d.page} of ${Math.max(1, Math.ceil(d.total / PER))}`;
    $("#b-prev").disabled = d.page <= 1; $("#b-next").disabled = d.page * PER >= d.total;
    const th = (k, t) => `<th class="sort" data-k="${k}">${t}${BS.sort === k ? (BS.desc ? " (desc)" : " (asc)") : ""}</th>`;
    const bv = META.bvLoaded;
    $("#btable").innerHTML = `<tr>${th("genome_id", "Genome ID and name")}${th("species", "Species and site")}${th("st", "ST")}${th("country", "Country, year")}${th("septatype", "Septatype")}${th("lab", "Lab")}<th>BV-BRC</th><th>SeptaTyper</th></tr>` +
      (d.rows.length ? d.rows.map(r => `<tr>
        <td class="nw">${r.bvbrc_url ? `<a href="${esc(r.bvbrc_url)}" target="_blank" rel="noopener">${esc(r.genome_id)}</a>` : esc(r.genome_id)}<span class="alias">${esc(r.genome_name || r.isolate_name)}</span></td>
        <td class="nw"><i>${esc(SHORT[r.species])}</i><span class="alias">${esc(SITE[r.habitat])}${r.source === "community" ? ", community" : ""}</span></td>
        <td class="nw">${esc(r.st || "-")}</td>
        <td class="nw">${esc(r.country || "-")}<span class="alias">${esc(r.year_raw || "-")}</span></td>
        <td class="mono nw">${mutSep(r.septatype)}</td>
        <td>${badge(r.lab_phenotype, 1)}</td>
        <td>${r.bvbrc_call ? badge(r.bvbrc_call, 1) : `<span class="muted" title="Not loaded yet">-</span>`}</td>
        <td class="nw">${r.septatyper_call ? badge(r.septatyper_call, 1) : "-"} ${r.star ? `<span class="star" title="SeptaTyper matches the lab result and BV-BRC does not">&#9733;</span>` : ""}</td></tr>`).join("")
        : `<tr><td colspan="8" class="muted">No isolates match these filters.</td></tr>`);
    $$("#btable th.sort").forEach(h => h.onclick = () => { BS.desc = BS.sort === h.dataset.k ? !BS.desc : false; BS.sort = h.dataset.k; BS.page = 1; loadBrowse(); });
  } catch (e) { $("#btable").innerHTML = `<tr><td class="err">${esc(e.message)}</td></tr>`; }
}
INIT.browse = async () => {
  const f = META.facets;
  $("#f-sp").innerHTML = opts(Object.entries(SP), "All");
  $("#f-st").innerHTML = opts(f.st.slice(0, 60).map(x => [x, x]), "All");
  $("#f-country").innerHTML = opts(f.country.map(x => [x, x]), "All");
  try { const s = await api("/api/v1/stats/accuracy"); META.bvLoaded = s.benchmark.loaded; } catch (e) { META.bvLoaded = false; }
  if (META.bvLoaded) $("#f-wins").disabled = false;
  else $("#bv-note").innerHTML = `<div class="info">BV-BRC computational predictions have not been loaded on this server yet. The BV-BRC column fills in after the administrator runs <code>scripts/fetch_bvbrc.py</code>.</div>`;
  let t; const reload = () => { BS.page = 1; loadBrowse(); };
  $$("#p-browse select, #p-browse input[type=number], #f-wins").forEach(e => e.onchange = reload);
  $$("#p-browse input[type=search]").forEach(e => e.oninput = () => { clearTimeout(t); t = setTimeout(reload, 350); });
  $("#b-prev").onclick = () => { BS.page--; loadBrowse(); }; $("#b-next").onclick = () => { BS.page++; loadBrowse(); };
  $("#f-reset").onclick = () => { $$("#p-browse .filters select, #p-browse .filters input").forEach(e => e.value = ""); $("#f-wins").checked = false; reload(); };
  loadBrowse();
};

/* --------------------------------------------------------------- analytics */
async function loadBars() {
  const d = await api("/api/v1/stats/septatypes?species=" + $("#a-sp").value + "&top=12");
  const max = Math.max(...d.rows.map(r => Math.max(...Object.values(r.by_habitat))), 1);
  $("#a-bars").innerHTML = d.rows.map(r => `<div class="hb"><span class="mono">${mutSep(r.septatype)}</span><div class="bars">${Object.entries(r.by_habitat).map(([h, n]) => `<div class="b ${h}" style="width:${Math.max(1, 100 * n / max)}%" title="${SITE[h]}: ${n}"></div>`).join("")}</div></div>`).join("") +
    `<p class="muted" style="font-size:.85rem;margin-top:8px">Top 12 septatypes of ${fmt(d.total)} isolates. Bar length is the number of isolates at each site.</p>`;
}
let plotlyP = null;
function loadPlotly() {
  if (window.Plotly) return Promise.resolve();
  return plotlyP || (plotlyP = new Promise((ok, no) => { const s = document.createElement("script"); s.src = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"; s.onload = ok; s.onerror = () => no(new Error("map library")); document.head.appendChild(s); }));
}
async function loadMap() {
  const p = new URLSearchParams(); [["species", "#g-sp"], ["st", "#g-st"], ["tetratype", "#g-tt"], ["phenotype", "#g-ph"]].forEach(([k, s]) => { if ($(s).value) p.set(k, $(s).value); });
  const d = await api("/api/v1/stats/geography?" + p);
  $("#g-table").innerHTML = `<tr><th>Country</th><th>Isolates</th><th>Resistant</th></tr>` + d.rows.slice(0, 15).map(r => `<tr><td>${esc(r.country)}</td><td class="stat">${fmt(r.n)}</td><td class="stat">${fmt(r.res)}</td></tr>`).join("");
  try {
    await loadPlotly();
    const dark = matchMedia("(prefers-color-scheme: dark)").matches;
    Plotly.react("map", [{ type: "choropleth", locations: d.rows.map(r => r.iso3), z: d.rows.map(r => r.n), text: d.rows.map(r => `${r.country}: ${r.n} isolates, ${r.res} resistant`), hoverinfo: "text",
      colorscale: [[0, "#e8f0fa"], [1, "#0a2540"]], marker: { line: { color: "#8aa0b8", width: .4 } }, colorbar: { title: "Isolates", thickness: 12 } }],
      { margin: { t: 0, b: 0, l: 0, r: 0 }, paper_bgcolor: "rgba(0,0,0,0)", geo: { showframe: false, showcoastlines: false, projection: { type: "natural earth" }, bgcolor: "rgba(0,0,0,0)", showcountries: true, countrycolor: dark ? "#26364a" : "#d5dee9" } },
      { displayModeBar: false, responsive: true });
  } catch (e) { $("#map").innerHTML = `<p class="muted">The map library could not be loaded. The table below has the same counts.</p>`; }
}
INIT.analytics = () => {
  $("#a-sp").innerHTML = opts(Object.entries(SP)); $("#a-sp").onchange = loadBars; loadBars();
  $("#g-sp").innerHTML = opts(Object.entries(SP), "All");
  $("#g-st").innerHTML = opts(META.facets.st.slice(0, 15).map(x => [x, x]), "All");
  $("#g-tt").innerHTML = opts(META.facets.tetratype.slice(0, 12).map(x => [x, x]), "All");
  $$("#p-analytics .filters select").forEach(e => e.onchange = loadMap); loadMap();
};

/* ------------------------------------------------------------------ submit */
INIT.submit = () => {
  $("#s-sp").innerHTML = opts(Object.entries(SP));
  const upd = () => { setSites($("#s-hb"), $("#s-sp").value, false); }; $("#s-sp").onchange = upd; upd();
  $("#s-look").onclick = async () => {
    const id = $("#s-gid").value.trim(); if (!id) return;
    $("#s-name").textContent = "Looking up";
    try { const r = await api("/api/v1/bvbrc/" + encodeURIComponent(id)); $("#s-name").textContent = r.found ? "BV-BRC name: " + r.genome_name : r.reason; $("#s-name").dataset.name = r.found ? r.genome_name : ""; }
    catch (e) { $("#s-name").textContent = e.message; }
  };
  $("#sf").onsubmit = async e => {
    e.preventDefault();
    const btn = $("#s-go"), err = $("#s-err"); err.innerHTML = ""; btn.disabled = true; const old = btn.textContent; btn.textContent = "Analysing";
    try {
      const f = $("#s-file").files[0]; if (!f) throw new Error("Choose a genome FASTA file.");
      const fd = new FormData();
      [["species", "#s-sp"], ["habitat", "#s-hb"], ["genome_id", "#s-gid"], ["country", "#s-country"], ["year", "#s-year"], ["source_text", "#s-src"], ["mic", "#s-mic"], ["lab_phenotype", "#s-ph"], ["email", "#s-email"]].forEach(([k, s]) => fd.append(k, $(s).value.trim()));
      fd.append("genome_name", $("#s-name").dataset.name || ""); fd.append("file", f);
      const r = await api("/api/v1/submissions", { method: "POST", body: fd });
      const msg = { needs_review: "Needs review: the septatype is new or the prediction disagrees with your laboratory result.", pending: "Pending: the prediction matches your result or no result was given. A curator will check it." }[r.status] || r.status;
      $("#s-result").innerHTML = renderResult(r.result, `<div class="info"><b>Submission ${esc(r.job_id)}</b><br>${esc(msg)}${r.matches === null ? " Without a laboratory phenotype the isolate cannot be approved." : ""}</div>`);
      $("#sf").reset(); $("#s-name").textContent = ""; $("#s-name").dataset.name = ""; upd();
    } catch (ex) { err.innerHTML = `<div class="err">${esc(ex.message)}</div>`; }
    finally { btn.disabled = false; btn.textContent = old; }
  };
};

/* ------------------------------------------------------------------ theory */
function drawPocket(mech) {
  const st = k => { const m = mech.find(x => x.position === k); return m ? m.effect : "intact"; };
  const lost = mech.filter(m => m.effect === "lost").length, ctx = mech.filter(m => m.effect === "context").length;
  const col = e => e === "lost" ? "var(--res)" : e === "context" ? "var(--gold)" : "var(--sus)";
  const node = (x, y, k) => { const m = mech.find(z => z.position === k); return `<g><circle cx="${x}" cy="${y}" r="30" fill="var(--white)" stroke="${col(m.effect)}" stroke-width="3"/><text x="${x}" y="${y + 8}" text-anchor="middle" font-size="24" font-weight="600" font-family="var(--mono)" fill="var(--ink)">${esc(m.residue)}</text><text x="${x}" y="${y - 40}" text-anchor="middle" font-size="13" fill="var(--muted)">${k}</text></g>`; };
  const link = (x, y, k, kind) => { const e = st(k), c = col(e), cx = 320, cy = 140; const dash = e === "lost" ? "6 6" : kind === "H-bond" ? "2 5" : "0";
    return `<line x1="${x}" y1="${y}" x2="${cx}" y2="${cy}" stroke="${c}" stroke-width="3" stroke-dasharray="${dash}"/>` + (e === "lost" ? `<text x="${(x + cx) / 2}" y="${(y + cy) / 2 - 6}" text-anchor="middle" font-size="22" font-weight="700" fill="${c}">&#215;</text>` : ""); };
  const bridge = lost >= 2 ? ["var(--res)", "disrupted", "6 5"] : (lost === 1 || ctx) ? ["var(--gold)", "weakened", "0"] : ["var(--sus)", "intact", "0"];
  $("#pocket").innerHTML = `
    ${link(120, 70, "GyrA-83", "H-bond")}${link(120, 240, "GyrA-87", "charge")}${link(520, 70, "ParC-80", "H-bond")}${link(520, 240, "ParC-84", "charge")}
    <line x1="320" y1="140" x2="320" y2="262" stroke="var(--blue)" stroke-width="4"/>
    <circle cx="320" cy="140" r="34" fill="var(--white)" stroke="${bridge[0]}" stroke-width="4" stroke-dasharray="${bridge[2]}"/>
    <text x="320" y="137" text-anchor="middle" font-size="16" font-weight="600" fill="var(--ink)">Mg&#178;&#8314;</text>
    <text x="320" y="154" text-anchor="middle" font-size="11" fill="var(--muted)">+ water</text>
    <rect x="240" y="262" width="160" height="44" rx="6" fill="var(--blue)"/>
    <text x="320" y="289" text-anchor="middle" font-size="14" font-weight="600" fill="#fff">Ciprofloxacin</text>
    ${node(120, 70, "GyrA-83")}${node(120, 240, "GyrA-87")}${node(520, 70, "ParC-80")}${node(520, 240, "ParC-84")}
    <text x="320" y="20" text-anchor="middle" font-size="14" font-weight="600" fill="${bridge[0]}">Water-Mg&#178;&#8314; bridge ${bridge[1]}</text>`;
}
async function showTetra(t) {
  const aas = t.split("");
  if (aas.length !== 4 || aas.some(a => !/[ACDEFGHIKLMNPQRSTVWY]/.test(a))) return;
  const sp = $("#t-sp").value;
  try {
    const r = await api("/api/v1/predict/manual", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ species: sp, habitat: (META.habitats[sp].includes("ALL") ? "ALL" : "FECAL"), gyrA_83: aas[0], gyrA_87: aas[1], parC_80: aas[2], parC_84: aas[3], recN_346: "L" }) });
    const p = r.prediction;
    drawPocket(p.mechanism);
    $("#t-kv").innerHTML = `<dt>Tetratype</dt><dd class="mono">${aas.join("-")}</dd><dt>Grammar call</dt><dd>${badge(p.grammar.call, 1)} ${esc(p.grammar.rule)}</dd><dt>Mean hydrophobicity &lt;H&gt;</dt><dd class="stat">${p.mean_h}</dd>`;
    $("#t-why").innerHTML = p.mechanism.map(m => `<li><span class="dot ${esc(m.effect)}"></span><span><b>${esc(m.position)} ${esc(m.residue)}</b>. ${esc(m.text)}</span></li>`).join("");
    $$("#t-pick button").forEach(b => b.classList.toggle("on", b.dataset.t === t));
  } catch (e) { $("#t-why").innerHTML = `<li class="err">${esc(e.message)}</li>`; }
}
INIT.theory = async () => {
  $("#t-sp").innerHTML = opts(Object.entries(SP));
  const pick = async () => {
    const r = await api(`/api/v1/rules?species=${$("#t-sp").value}&habitat=ALL`);
    const seen = new Map(); r.rows.forEach(x => { if (!seen.has(x.tetratype)) seen.set(x.tetratype, 0); seen.set(x.tetratype, seen.get(x.tetratype) + x.n_isolates); });
    const top = [...seen.entries()].sort((a, b) => b[1] - a[1]).slice(0, 8).map(x => x[0].replace(/-/g, ""));
    $("#t-pick").innerHTML = top.map(t => `<button type="button" data-t="${t}">${t.split("").join("-")}</button>`).join("");
    $$("#t-pick button").forEach(b => b.onclick = () => showTetra(b.dataset.t));
    if (top.length) showTetra(top[0]);
  };
  $("#t-sp").onchange = pick; pick();
  $("#t-custom").oninput = e => { const v = e.target.value.toUpperCase(); if (v.length === 4) showTetra(v); };
  try {
    const [s, aa] = await Promise.all([api("/api/v1/stats/accuracy"), api("/api/v1/aminoacids")]);
    const b = s.benchmark, rows = s.rows.filter(x => x.habitat !== "ALL");
    let h = `<h3 style="font-size:1rem">Majority rule and grammar against laboratory phenotype</h3><div class="scroll"><table><tr><th>Species</th><th>Site</th><th>Isolates</th><th>Majority rule</th><th>Grammar only</th></tr>${rows.map(x => `<tr><td><i>${esc(SHORT[x.species])}</i></td><td>${esc(SITE[x.habitat])}</td><td class="stat">${fmt(x.n_isolates)}</td><td class="stat">${x.accuracy}%</td><td class="stat">${x.grammar_accuracy}%</td></tr>`).join("")}</table></div><p class="muted" style="font-size:.85rem">Both are scored on the founding isolates the rules were built from. They describe the reference set, not independent performance. The grammar only calls septatypes missing from the database.</p>`;
    h += `<h3 style="font-size:1rem;margin-top:18px">Head to head with BV-BRC computational predictions</h3>`;
    if (b.loaded) h += `<p class="stat">${fmt(b.compared)} isolates have a BV-BRC computational prediction. SeptaTyper matches the lab result in ${fmt(b.septatyper_ok)} (${(100 * b.septatyper_ok / b.compared).toFixed(1)}%) and BV-BRC in ${fmt(b.bvbrc_ok)} (${(100 * b.bvbrc_ok / b.compared).toFixed(1)}%). SeptaTyper is right and BV-BRC wrong in ${fmt(b.wins)}. The reverse happens in ${fmt(b.losses)}.</p><div class="scroll"><table><tr><th>Genome</th><th>Septatype</th><th>Lab</th><th>BV-BRC</th><th>SeptaTyper</th></tr>${b.examples.map(x => `<tr><td>${esc(x.genome_id)}<span class="alias">${esc(x.isolate_name)}</span></td><td class="mono">${esc(x.septatype)}</td><td>${badge(x.lab_phenotype, 1)}</td><td>${badge(x.bvbrc_call, 1)}</td><td>${badge(x.septatyper_call, 1)}</td></tr>`).join("")}</table></div>`;
    else h += `<div class="info">BV-BRC predictions have not been loaded yet, so the head-to-head comparison is empty. Run <code>scripts/fetch_bvbrc.py</code> on the server to fill it.</div>`;
    h += `<h3 style="font-size:1rem;margin-top:18px">Susceptible tetratype S-D-S-E, split by qnr</h3><div class="scroll"><table><tr><th>Species</th><th>qnrB</th><th>qnrS</th><th>Isolates</th><th>Resistant</th></tr>${s.edgecases.wildtype_tetratype_by_qnr.map(x => `<tr><td><i>${esc(SHORT[x.species])}</i></td><td>${x.qnrB2 ? "yes" : "no"}</td><td>${x.qnrS1 ? "yes" : "no"}</td><td class="stat">${x.n}</td><td class="stat">${x.res} (${Math.round(100 * x.res / x.n)}%)</td></tr>`).join("")}</table></div><p class="muted" style="font-size:.85rem">A qnr gene on a susceptible baseline is usually resistant, but not always. A gene found in sequence may not be expressed.</p>`;
    $("#t-valid").innerHTML = h;
    $("#t-aa").innerHTML = `<tr><th>Residue</th><th>Name</th><th>Hydrophobicity</th><th>Class</th></tr>` + aa.rows.sort((a, b) => b.hydrophobicity - a.hydrophobicity).map(x => `<tr><td class="mono">${x.aa}</td><td>${esc(x.name)}</td><td class="stat">${x.hydrophobicity.toFixed(2)}</td><td>${esc(x.class)}</td></tr>`).join("");
  } catch (e) { $("#t-valid").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
};

/* -------------------------------------------------------------- data, about */
INIT.data = () => {
  const o = location.origin;
  $("#ex-curl").textContent = `# Resistant E. coli urine isolates, first page\ncurl "${o}/api/v1/isolates?species=ecoli&habitat=URINE&phenotype=Resistant"\n\n# Predict from a genome file\ncurl -F file=@genome.fasta -F species=kpneumoniae -F habitat=URINE ${o}/api/v1/predict`;
  $("#ex-py").textContent = `import requests\nbase = "${o}/api/v1"\nr = requests.post(f"{base}/predict", data={"species": "ecoli", "habitat": "URINE"},\n                  files={"file": open("genome.fasta", "rb")}).json()\nprint(r["septatype"], r["prediction"]["call"])`;
};
INIT.about = () => { $("#about-build").textContent = `Rules built ${META.build.built || ""} from ${META.build.source || ""} (${META.build.n_isolates || ""} isolates).`; };

/* ----------------------------------------------------------------- curation */
INIT.curate = () => {
  const load = async () => {
    $("#c-err").innerHTML = "";
    try {
      const d = await api("/api/v1/admin/queue?status=" + $("#c-status").value, { headers: { "X-Admin-Token": $("#c-token").value } });
      $("#c-table").innerHTML = `<tr><th>Submission</th><th>Isolate</th><th>Septatype</th><th>Prediction</th><th>Lab</th><th>Status</th><th>Decision</th></tr>` +
        (d.rows.length ? d.rows.map(r => `<tr><td class="nw">${esc(r.job_id)}<span class="alias">${esc(r.created)}</span></td>
          <td>${esc(r.genome_id || "no ID")}<span class="alias">${esc(SHORT[r.species])}, ${esc(SITE[r.habitat])}, ${esc(r.country || "-")} ${esc(r.year || "")}</span></td>
          <td class="mono nw">${mutSep(r.septatype)}${r.novel ? `<span class="alias">new septatype</span>` : ""}</td>
          <td>${r.call ? badge(r.call, 1) : "-"}<span class="alias">${esc(r.tier || "")}</span></td>
          <td>${r.lab_phenotype ? badge(r.lab_phenotype, 1) : "-"}${r.matches === 0 ? `<span class="alias">disagrees</span>` : ""}</td>
          <td>${esc(r.status)}<span class="alias">${esc(r.curator_note || "")}</span></td>
          <td>${["approved", "rejected"].includes(r.status) ? "" : `<input type="text" placeholder="note" data-n="${esc(r.job_id)}" style="margin-bottom:6px"><button class="btn solid" data-a="approve" data-j="${esc(r.job_id)}">Approve</button> <button class="btn line" data-a="reject" data-j="${esc(r.job_id)}">Reject</button>`}</td></tr>`).join("")
          : `<tr><td colspan="7" class="muted">The queue is empty.</td></tr>`);
      $$("#c-table button[data-a]").forEach(b => b.onclick = async () => {
        try {
          await api("/api/v1/admin/decide", { method: "POST", headers: { "Content-Type": "application/json", "X-Admin-Token": $("#c-token").value }, body: JSON.stringify({ job_id: b.dataset.j, action: b.dataset.a, note: ($(`[data-n="${b.dataset.j}"]`) || {}).value || "" }) });
          LOADED.browse = false; META.counts.genomes += b.dataset.a === "approve" ? 1 : 0; load();
        } catch (e) { $("#c-err").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
      });
    } catch (e) { $("#c-err").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
  };
  $("#c-load").onclick = load;
};

/* -------------------------------------------------------------------- start */
(async () => {
  try {
    META = await api("/api/v1/meta");
    Object.assign(SP, META.species);
    LOADED.home = false;
    route();
  } catch (e) {
    $("#main").innerHTML = `<div class="band"><div class="err">${esc(e.message)}</div></div>`;
  }
})();
