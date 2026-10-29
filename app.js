/* Ciprotyper front end. No framework. All server text goes through esc(). */
"use strict";
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = n => Number(n).toLocaleString("en-US");
let META = null;
const SP = {};                       // key -> display name
const SHORT = { ecoli: "E. coli", ssonnei: "S. sonnei", senterica: "S. enterica", kpneumoniae: "K. pneumoniae" };
const SITE = { ALL: "Unknown / not specified" };           // the rest comes from the server (META.sources)
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
const withN = arr => arr.map(([v, n]) => [v, `${v} (${fmt(n)})`]);
const dl = arr => arr.map(([v, n]) => `<option value="${esc(v)}" label="${fmt(n)} isolates"></option>`).join("");
const opts = (arr, first) => (first ? `<option value="">${first}</option>` : "") + arr.map(([v, t]) => `<option value="${esc(v)}">${esc(t)}</option>`).join("");
function setSites(sel, sp, withAll) {
  const keys = (META.habitats[sp] || []).filter(h => h !== "ALL");
  sel.innerHTML = opts((withAll ? ["ALL", ...keys] : keys).map(h => [h, SITE[h] || h]));
}
function useSources() { META.sources.forEach(s => { SITE[s.key] = s.label; SITE_LC[s.key] = s.short; }); }
async function refreshMeta() {          // approved submissions change the counts and lists
  try { const m = await api("/api/v1/meta"); META.facets = m.facets; META.counts = m.counts; META.stats = m.stats; META.sources = m.sources; useSources(); } catch (e) { /* keep what we have */ }
}
const PAL = ["#1d4e89", "#4d9de0", "#9cc4f5", "#00897b", "#8e6bbf", "#c98a00", "#6d7f95", "#d98aa0", "#4f9d69", "#a0522d", "#2a9d8f", "#7a7a7a"];
const hcol = h => PAL[Math.max(0, META.sources.findIndex(s => s.key === h)) % PAL.length];

/* ------------------------------------------------------------------ router */
const PAGES = ["home", "predict", "browse", "analytics", "submit", "theory", "data", "about", "curate", "calc"];
function route() {
  const p = (location.hash.replace(/^#\/?/, "") || "home").split("?")[0];
  const page = PAGES.includes(p) ? p : "home";
  $$(".page").forEach(e => e.classList.toggle("on", e.id === "p-" + page));
  $$("#nav a").forEach(a => a.classList.toggle("on", a.dataset.p === page));
  if (page === "calc") {
    $$("#nav a").forEach(a => a.classList.toggle("on", a.dataset.p === "theory"));
    if (META) INIT.calc();
    return;
  }
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
  $("#home-stats").innerHTML = `<tr><th>Species</th><th>Source</th><th>Isolates</th><th>Resistant</th><th>Septatypes</th><th>Accuracy</th></tr>` +
    rows.map(s => `<tr><td><i>${esc(SP[s.species])}</i></td><td>${esc(SITE[s.habitat])}</td><td class="stat">${fmt(s.n_isolates)}</td><td class="stat">${fmt(s.n_resistant)}</td><td class="stat">${s.n_septatypes}</td><td class="stat">${s.accuracy}%</td></tr>`).join("");
};

/* ----------------------------------------------------------------- predict */
let PMODE = "file", PFILE = null;

/* ------------------------------------------------------- exports and batch */
let LAST = null, BATCH = [], PFILES = [];
const csvCell = v => `"${String(v ?? "").replace(/"/g, '""')}"`;
const toCsv = rows => { const cols = Object.keys(rows[0]); return [cols.join(","), ...rows.map(r => cols.map(c => csvCell(r[c])).join(","))].join("\n"); };
function flat(res, name, error) {
  const p = (res && res.prediction) || {}, a = (res && res.analysis) || {}, pos = a.positions || {}, cd = k => (pos[k] || {}).codon || "";
  const bv = (res && res.bvbrc) || {};
  return { file: name, bvbrc_genome_id: bv.genome_id || "", bvbrc_genome_name: bv.genome_name || "", bvbrc_url: bv.url || "", reference_isolate: bv.in_reference ? "yes" : (bv.genome_id ? "no" : ""), species: res ? res.species_name : "", isolation_source: res ? res.habitat : "", septatype: (res && res.septatype) || "",
    call: error ? "Error" : (p.call || "No call"), basis: p.tier || "", n_isolates: p.n_isolates ?? "", n_resistant: p.n_resistant ?? "", n_susceptible: p.n_susceptible ?? "",
    sample_st: ((a.mlst || {}).st) || "", mlst_status: ((a.mlst || {}).status) || "",
    reference_sites: (p.reference_sites || []).map(x => `${x.habitat}:${x.n}`).join("; "),
    reference_countries: (p.reference_countries || []).map(x => `${x.country}:${x.n}`).join("; "),
    reference_sts: (p.reference_sts || []).map(x => `${x.st}:${x.n} (R ${x.n_resistant}, S ${x.n_susceptible})`).join("; "),
    H_GyrA_helix: ((p.criteria || {}).hydrophobicity || {}).gyr ? p.criteria.hydrophobicity.gyr.mean : "", H_ParC_helix: ((p.criteria || {}).hydrophobicity || {}).parc ? p.criteria.hydrophobicity.parc.mean : "", H_QRDR_helix: p.mean_h ?? "", gyrA83_codon: cd("gyrA83"), gyrA87_codon: cd("gyrA87"), parC80_codon: cd("parC80"), parC84_codon: cd("parC84"), recN346_codon: cd("recN346"),
    notes: error || ((res && res.warnings) || []).join(" | ") };
}
function download(name, text, type) {
  const u = URL.createObjectURL(new Blob([text], { type })), a = document.createElement("a");
  a.href = u; a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(u), 1000);
}
function printNow() {
  const ds = $$("details"), was = ds.map(d => d.open); ds.forEach(d => d.open = true);
  window.print(); ds.forEach((d, i) => d.open = was[i]);
}
async function maybeGzip(f) {      // smaller upload on a phone connection; the server accepts .gz
  if (!("CompressionStream" in window) || /\.gz$/i.test(f.name) || f.size < 500000) return f;
  try { const b = await new Response(f.stream().pipeThrough(new CompressionStream("gzip"))).blob(); return new File([b], f.name + ".gz", { type: "application/gzip" }); }
  catch (e) { return f; }
}
const TIER = { database: "Database", database_other_site: "Database (other source)", grammar: "Grammar (new)" };
function drawBatch(done, total, unit = "file") {
  const rows = BATCH.map((x, i) => {
    if (x.error) return `<tr><td>${esc(x.name)}</td><td colspan="5" class="muted">${esc(x.error)}</td></tr>`;
    const p = x.res.prediction;
    const bv = x.res.bvbrc, sub = bv ? `<span class="alias"><i>${esc(SHORT[x.res.species] || "")}</i> ${esc(bv.genome_name)}${bv.in_reference ? " (reference isolate)" : ""}</span>` : "";
    return `<tr class="clk" data-act="detail" data-i="${i}"><td>${esc(x.name)}${sub}</td><td class="mono nw">${p ? mutSep(x.res.septatype) : "-"}</td><td class="nw">${esc(((x.res.analysis || {}).mlst || {}).st || "-")}</td><td>${p ? badge(p.call, 1) : "No call"}</td><td>${p ? esc(TIER[p.tier] || p.tier) : "-"}</td><td class="stat">${p && p.known ? p.n_resistant + " of " + p.n_isolates + " resistant" : "-"}</td></tr>`;
  }).join("");
  const working = done < total;
  $("#presult").innerHTML = `<div class="band"><h3>Batch result</h3>${working ? `<p class="muted">Analysing ${unit} ${done + 1} of ${total}. Keep this page open.</p><div class="prog"><i style="width:${100 * done / total}%"></i></div>` : ""}
    <div class="scroll"><table><tr><th>${unit === "genome" ? "Genome" : "File"}</th><th>Septatype</th><th>ST</th><th>Call</th><th>Call from</th><th>Support</th></tr>${rows}</table></div>
    ${working ? "" : `<div class="actions no-print"><button class="btn solid" data-act="bcsv" type="button">Download CSV</button><button class="btn line" data-act="print" type="button">Print or save as PDF</button></div><p class="muted" style="font-size:.85rem">Tap a row for the full result.</p>`}
    <p class="print-only">Ciprotyper batch report, ${new Date().toLocaleDateString()}. Research use only.</p></div><div id="bdetail"></div>`;
}
async function runBatch(species, habitat) {
  if (!PFILES.length) throw new Error("Choose at least one FASTA file.");
  BATCH = [];
  for (let i = 0; i < PFILES.length; i++) {
    drawBatch(i, PFILES.length);
    const f = PFILES[i], fd = new FormData();
    fd.append("species", species); fd.append("habitat", habitat); fd.append("qnrB2", "auto"); fd.append("qnrS1", "auto"); fd.append("file", await maybeGzip(f));
    try { const res = await api("/api/v1/predict", { method: "POST", body: fd }); res._name = f.name; BATCH.push({ name: f.name, res }); }
    catch (e) { BATCH.push({ name: f.name, error: e.message }); }
  }
  drawBatch(PFILES.length, PFILES.length);
}
async function runBvbrc(fallback, habitat) {
  const ids = [...new Set($("#bv-ids").value.split(/[\s,;]+/).map(x => x.trim()).filter(Boolean))];
  if (!ids.length) throw new Error("Type at least one BV-BRC genome ID, for example 562.7710.");
  if (ids.length > 10) throw new Error(`Up to 10 IDs at a time. You typed ${ids.length}.`);
  BATCH = [];
  for (let i = 0; i < ids.length; i++) {
    drawBatch(i, ids.length, "genome");
    const id = ids[i];
    if (!/^\d{1,9}\.\d{1,6}$/.test(id)) { BATCH.push({ name: id, error: "Not a BV-BRC genome ID. An ID looks like 562.7710." }); continue; }
    try {
      const res = await api("/api/v1/predict/bvbrc", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ genome_id: id, species: "auto", fallback_species: fallback, habitat }) });
      res._name = `${id} ${res.bvbrc.genome_name}`; BATCH.push({ name: id, res });
    } catch (e) { BATCH.push({ name: id, error: e.message }); }
  }
  if (ids.length === 1 && BATCH[0].res) { $("#presult").innerHTML = renderResult(BATCH[0].res); return; }      // one ID: show the full result at once
  drawBatch(ids.length, ids.length, "genome");
}
document.addEventListener("click", e => {
  const b = e.target.closest("[data-act]"); if (!b) return;
  const act = b.dataset.act;
  if (act === "print") printNow();
  else if (act === "csv" && LAST) download("ciprotyper_result.csv", toCsv([flat(LAST, LAST._name || "result")]), "text/csv");
  else if (act === "json" && LAST) download("ciprotyper_result.json", JSON.stringify(LAST, null, 2), "application/json");
  else if (act === "bcsv") download("ciprotyper_batch.csv", toCsv(BATCH.map(x => flat(x.res, x.name, x.error))), "text/csv");
  else if (act === "detail") { const x = BATCH[+b.dataset.i]; if (x && x.res) { $("#bdetail").innerHTML = renderResult(x.res); $("#bdetail").scrollIntoView({ behavior: "smooth" }); } }
});

const SITE_LC = {};
const PROT = { gyrA: "GyrA", parC: "ParC", recN: "RecN" };
const f2 = v => (v < 0 ? "-" : "") + Math.abs(v).toFixed(2);
const SYNTAX = [["GyrA-83", "S Y T", "F L I W"], ["GyrA-87", "D N E", "N D G Y A E"], ["ParC-80", "S Y T", "I F W"], ["ParC-84", "D E", "E V D"]];
const f3 = v => (v < 0 ? "-" : "") + Math.abs(v).toFixed(3);
/* QRDR helix as a strip of letters; the four tetratype positions are marked */
const helixStrip = b => `<span class="helix">${b.values.map(v => `<span class="${v.sub ? "sub" : ""}">${esc(v.residue)}</span>`).join("")}</span>`;
/* the three lines of the <H> calculation for one sample: GyrA helix, ParC helix, combined */
const helixLines = (label, v) => `<p class="formula wrap">${label} GyrA helix ${v.gyr.start}-${v.gyr.end} ${helixStrip(v.gyr)} sum ${f2(v.gyr.sum)} / ${v.gyr.n} = <b>${f3(v.gyr.mean)}</b></p>
    <p class="formula wrap">${label} ParC helix ${v.parc.start}-${v.parc.end} ${helixStrip(v.parc)} sum ${f2(v.parc.sum)} / ${v.parc.n} = <b>${f3(v.parc.mean)}</b></p>
    <p class="formula wrap">${label} &lt;H&gt; = (${f3(v.gyr.mean)} + ${f3(v.parc.mean)}) / 2 = <b>${f3(v.mean)}</b></p>`;
function critTable(positions) {
  const chips = (list, side) => list.map(c => `<span class="rchip${c.hit ? " hit " + side : ""}${c.exception ? " exc" : ""}">${esc(c.residue)}${c.exception ? "<small>exception</small>" : ""}</span>`).join("");
  return `<div class="scroll"><table class="rtab"><tr><th>Position</th><th>Susceptible</th><th>Resistant</th></tr>${positions.map(x =>
    `<tr><td class="nw"><b>${esc(x.position)}</b><span class="alias">sample: ${esc(x.residue)}</span></td><td>${chips(x.susceptible, "s")}</td><td>${chips(x.resistant, "r")}</td></tr>`).join("")}</table></div>`;
}
/* number line: sample and reference <H> against the cutoff; extra = typical resistant and susceptible isolates */
function hChart(sample, ref, cutoff, extra) {
  const lo = 0.1, hi = 0.6, W = 560, pad = 28, H = extra ? 124 : 94;
  const x = v => pad + (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo) * (W - 2 * pad), lab = v => Math.max(64, Math.min(W - 64, x(v)));
  const col = sample > cutoff ? "var(--res)" : "var(--sus)";
  let ticks = "";
  if (extra) {
    const left = x(extra.sus) <= x(extra.res), tk = (v, t, c, anchor, dx) =>
      `<path d="M${x(v)} 62 l-6 11 h12 z" fill="${c}"/><text x="${x(v) + dx}" y="${H - 8}" text-anchor="${anchor}" font-size="11" fill="var(--muted)">${t} ${f3(v)}</text>`;
    ticks = tk(extra.sus, "typical susceptible isolate", "var(--sus)", left ? "end" : "start", left ? 8 : -8) + tk(extra.res, "typical resistant isolate", "var(--res)", left ? "start" : "end", left ? -8 : 8);
  }
  return `<svg viewBox="0 0 ${W} ${H}" class="hchart" role="img" aria-label="Sample hydrophobicity ${f3(sample)}, reference ${f3(ref)}, cutoff ${f3(cutoff)}">
    <line x1="${x(lo)}" x2="${x(cutoff)}" y1="48" y2="48" stroke="var(--sus-bg)" stroke-width="14"/><line x1="${x(cutoff)}" x2="${x(hi)}" y1="48" y2="48" stroke="var(--res-bg)" stroke-width="14"/>
    <text x="${x(lo)}" y="12" font-size="11" fill="var(--muted)">lower H</text><text x="${x(hi)}" y="12" font-size="11" text-anchor="end" fill="var(--muted)">higher H</text>
    <line x1="${x(cutoff)}" x2="${x(cutoff)}" y1="38" y2="76" stroke="var(--ink)" stroke-width="2" stroke-dasharray="4 3"/><text x="${x(cutoff)}" y="14" font-size="12" text-anchor="middle" fill="var(--ink)">cutoff ${f3(cutoff)}</text>
    <circle cx="${x(ref)}" cy="48" r="11" fill="var(--white)" stroke="var(--muted)" stroke-width="3"/><circle cx="${x(sample)}" cy="48" r="7" fill="${col}"/>
    <text x="${lab(sample)}" y="30" text-anchor="middle" font-size="13" font-weight="600" fill="${col}" class="hs">Sample ${f3(sample)}</text>
    <text x="${lab(ref)}" y="${extra ? 90 : 82}" text-anchor="middle" font-size="12" fill="var(--muted)">Reference ${f3(ref)}</text>${ticks}</svg>`;
}
function whyBlock(p) {
  const c = p.criteria; if (!c) return "";
  const k = c.checkpoints, r1 = k.residues, hy = c.hydrophobicity, rf = c.reference, dec = c.decided_by, res = c.positions.map(x => x.residue).join("");
  const intro = p.tier === "grammar" ? "Call from the three checkpoints, because this septatype is new." : `Call from the database majority rule${c.grammar_call !== p.call ? `. The three checkpoints alone would say ${esc(c.grammar_call)}` : ""}.`;
  const tag = on => on ? `<span class="tag">decides</span>` : "";
  const exc = (r1.exceptions || []).length ? ` (exception: the ParC pair S-K)` : "";
  const why1 = r1.override ? `${r1.pattern ? "All four residues are in the Susceptible column, but " : "The four residues are not all in one column. "}${esc(r1.rule)}: counted Resistant. This is an extra rule from the reference isolates, not one of the three checkpoints.`
    : r1.pattern ? `All four residues are in the ${r1.pattern === "resistant" ? "Resistant" : "Susceptible"} column${exc}.`
    : r1.resolved_by === "hydrophobicity" ? `The four residues are not all in one column, so the written rules do not decide. One position votes Resistant, so the QRDR helix mean &lt;H&gt; settles it.`
    : `The four residues are not all in one column, so the written rules do not decide. The residue votes settle it: ${r1.n_votes} of 4 positions vote Resistant.`;
  const hBlock = r1.resolved_by === "hydrophobicity" ? `${hChart(hy.mean, rf.mean, hy.cutoff)}
    ${helixLines("Sample", hy)}${helixLines("Reference", rf)}
    <p class="muted small" style="margin:4px 0">Every residue of the GyrA and ParC helices counts. Your residues at the four tetratype positions are marked red. Means are rounded to three decimals.</p>
    <p class="small"><a class="calc-link" href="#/calc?r=${esc(res)}&ref=${esc(rf.residues)}">How is this calculated?</a> &nbsp; <a class="calc-link" href="#/calc?r=${esc(res)}&ref=${esc(rf.residues)}&s=cutoff">Why a cutoff of ${f3(hy.cutoff)}?</a></p>` : "";
  const dim = x => x === "Resistant" ? "hit" : "dim";
  let h = `<h3>Why ${esc(p.call)}</h3><p class="muted small">${intro} Highlighted is what applies to this sample.</p>`;
  h += `<div class="crit ${dim(r1.call)}"><h4>Checkpoint 1: residues at GyrA-83, GyrA-87, ParC-80 and ParC-84 ${badge(r1.call, 1)} ${tag(["pattern", "votes", "hydrophobicity", "override"].includes(dec))}</h4>${critTable(c.positions)}
    <p style="margin:8px 0 0">${why1} <b>Checkpoint 1 says ${esc(r1.call)}.</b></p>${hBlock}</div>`;
  h += `<div class="crit ${dim(k.recn.call)}"><h4>Checkpoint 2: RecN-346 ${badge(k.recn.call, 1)} ${tag(dec === "checkpoint2")}</h4>
    <p style="margin:0">${esc(k.recn.text)}.</p><p class="muted small" style="margin:6px 0 0">Susceptible: leucine encoded by cta. Resistant: leucine encoded by ctg, or phenylalanine.</p></div>`;
  h += `<div class="crit ${dim(k.qnr.call)}"><h4>Checkpoint 3: PMQR genes ${badge(k.qnr.call, 1)} ${tag(dec === "checkpoint3")}</h4>
    <p style="margin:0">${esc(k.qnr.text.charAt(0).toUpperCase() + k.qnr.text.slice(1))}.</p><p class="muted small" style="margin:6px 0 0">Susceptible: both qnrB2 and qnrS1 absent. Resistant: either or both present.</p></div>`;
  h += `<p class="small" style="margin:10px 0 0"><b>Combined:</b> susceptible only if all three checkpoints are susceptible, resistant if any is. Here checkpoint 1 is ${esc(r1.call.toLowerCase())}, checkpoint 2 is ${esc(k.recn.call.toLowerCase())} and checkpoint 3 is ${esc(k.qnr.call.toLowerCase())}, so the three checkpoints give <b>${esc(c.grammar_call)}</b>.</p>`;
  return h;
}
function renderResult(r, extra) {
  const p = r.prediction, a = r.analysis;
  LAST = r;
  let h = `<div class="band result"><p class="print-only">Ciprotyper report, ${new Date().toLocaleDateString()}. ${esc(r.species_name)}, ${esc(SITE[r.habitat] || r.habitat)}. ${esc(r._name || "")}</p>`;
  if (r.bvbrc) h += `<p class="small muted" style="margin:0 0 8px">BV-BRC genome <a href="${esc(r.bvbrc.url)}" target="_blank" rel="noopener">${esc(r.bvbrc.genome_id)}</a>: ${esc(r.bvbrc.genome_name)}${r.bvbrc.country ? ", " + esc(r.bvbrc.country) : ""}${r.bvbrc.year ? ", " + esc(r.bvbrc.year) : ""}${r.bvbrc.isolation_source ? ". Source: " + esc(r.bvbrc.isolation_source) : ""}.</p>`;
  if (!p) {
    h += `<h3>No call could be made</h3><p>A complete septatype could not be built from this input.</p>`;
  } else {
    const sitesTxt = (p.reference_sites || []).map(x => `${SITE_LC[x.habitat] || String(x.habitat).toLowerCase()} ${x.n}`).join(", ");
    const cs = p.reference_countries || [], top = cs.slice(0, 3).map(x => `${esc(x.country)} ${x.n}`).join(", ");
    const countryTxt = !cs.length ? "" : cs.length === 1 ? `, all from ${esc(cs[0].country)}` : `, from ${top}${cs.length > 3 ? ` and ${cs.length - 3} other ${cs.length - 3 === 1 ? "country" : "countries"}` : ""}`;
    const basis = p.tier === "database" ? `This septatype is in <b>${fmt(p.n_isolates)}</b> reference isolates (${esc(sitesTxt)})${countryTxt}. Based on the majority rule, it is <b>${esc(p.call.toLowerCase())}</b>.`
      : p.tier === "database_other_site" ? `Not seen for this isolation source. It is in <b>${fmt(p.n_isolates)}</b> reference isolates from other sources (${esc(sitesTxt)})${countryTxt}. Based on the majority rule, it is <b>${esc(p.call.toLowerCase())}</b>.`
      : `New septatype, so the call comes from the grammar.`;
    h += `<div class="verdict-row">${badge(p.call)}</div><p class="basis">${basis}</p>`;
    const codon = (a && a.recN_codon) || r.recn_codon;
    h += `<div class="sep">${p.septatype.split("-").map((v, i) => {
      const ch = v !== WILD[i] && !(i >= 5 && v === "0");
      return `<div class="tile ${ch ? "changed" : ""}"><div class="cell">${esc(v)}</div><small>${POS[i][1]}</small>${i === 4 && codon ? `<span class="codon">${esc(codon.toLowerCase())}</span>` : ""}</div>`;
    }).join("")}</div>`;
    if (r.species === "senterica" && (a && a.recN_syn === "ctg" || String(codon).toUpperCase() === "CTG"))
      h += `<div class="info">RecN-346 uses the ctg codon for leucine, a high-usage synonymous change flagged in <i>S. enterica</i>.</div>`;
    const ml = a && a.mlst;
    if (ml) {
      h += `<p class="stline"><span class="muted">Sequence type</span> ${ml.st ? `<b class="mono">${esc(ml.st)}</b>` : `<span class="muted">${esc(ml.note || "not typed")}</span>`}</p>`;
      if (ml.alleles) h += `<details><summary>MLST alleles (${esc(ml.scheme_name || "")})</summary><div class="scroll"><table><tr>${Object.keys(ml.alleles).map(k => `<th>${esc(k)}</th>`).join("")}</tr><tr>${Object.values(ml.alleles).map(v => `<td class="mono">${esc(v)}</td>`).join("")}</tr></table></div></details>`;
    }
    if (p.known) {
      h += `<p class="stat">${fmt(p.n_resistant)} resistant, ${fmt(p.n_susceptible)} susceptible.</p>${bar(p.n_resistant, p.n_susceptible)}`;
      if (p.reference_sts && p.reference_sts.length) {
        const mine = ml && ml.st;
        h += `<h3>Those reference isolates, by ST</h3><div class="scroll"><table><tr><th>ST</th><th>Isolates</th><th>Resistant</th><th>Susceptible</th><th></th></tr>${p.reference_sts.map(x => `<tr><td class="mono">${esc(x.st)}</td><td class="stat">${x.n}</td><td class="stat">${x.n_resistant}</td><td class="stat">${x.n_susceptible}</td><td>${mine && x.st === mine ? "<b>your genome</b>" : ""}</td></tr>`).join("")}</table></div>`;
        if (mine && !p.reference_sts.some(x => x.st === mine)) h += `<p class="muted small">Your genome (${esc(mine)}) is not one of these STs.</p>`;
      }
      if (cs.length) h += `<details><summary>Those reference isolates, by country and source</summary><div class="scroll"><table><tr><th>Country</th><th>Isolates</th><th>Source and phenotype</th></tr>${cs.map(x => `<tr><td>${esc(x.country)}</td><td class="stat">${x.n}</td><td>${Object.entries(x.site_phenotypes || {}).sort((a, b) => b[1].n - a[1].n).map(([k, v]) => { const nm = String(SITE_LC[k] || k); return `<span class="nw">${esc(nm.charAt(0).toUpperCase() + nm.slice(1))} ${v.n} (R-${v.r}, S-${v.s})</span>`; }).join(", ")}</td></tr>`).join("")}</table></div></details>`;
    } else if (p.similar && p.similar.length) {
      h += `<p class="muted" style="margin-top:12px">Known septatypes with the same GyrA/ParC tetratype:</p><div class="scroll"><table><tr><th>Septatype</th><th>Isolates</th><th>Resistant</th><th>Susceptible</th></tr>${p.similar.map(x => `<tr><td class="mono">${esc(x.septatype)}</td><td>${x.n_isolates}</td><td>${x.n_resistant}</td><td>${x.n_susceptible}</td></tr>`).join("")}</table></div>`;
    }
    h += whyBlock(p);
    h += `<details><summary>Structural reading (schematic)</summary><ul class="why">${p.mechanism.map(m => `<li><span class="dot ${esc(m.effect)}"></span><span><b>${esc(m.position)} ${esc(m.residue)}</b>. ${esc(m.text)}</span></li>`).join("")}</ul></details>`;
  }
  if (r.warnings && r.warnings.length) h += `<div class="note"><ul>${r.warnings.map(w => `<li>${esc(w)}</li>`).join("")}</ul></div>`;
  if (a && a.positions) {
    const ref = a.reference || { genes: {}, qnr: {} };
    h += `<details><summary>Sequence details</summary>
      <p class="muted small">Your sequence is compared with the ${esc(ref.species || "")} reference genes (${Object.entries(ref.genes || {}).map(([g, ac]) => `${esc(g)} ${esc(ac)}`).join(", ")}).</p>
      <div class="scroll"><table><tr><th>Position</th><th>Sample codon</th><th>Sample residue</th><th>Reference codon</th><th>Reference residue</th></tr>` +
      Object.values(a.positions).map(x => `<tr><td class="nw">${PROT[x.gene] || esc(x.gene)}-${x.pos}</td><td class="mono">${esc((x.codon || "-").toLowerCase())}</td><td class="mono">${esc(x.aa || "-")}</td><td class="mono">${esc((x.ref_codon || "-").toLowerCase())}</td><td class="mono">${esc(x.ref_aa || "-")}</td></tr>`).join("") +
      `</table></div><div class="scroll" style="margin-top:10px"><table><tr><th>qnr gene</th><th>Result</th><th>Match to reference</th></tr>` +
      Object.entries(a.qnr).map(([k, x]) => `<tr><td>${k === "qnrB2" ? "qnrB" : "qnrS"}</td><td>${x.present ? "Present" : "Absent"}${x.source === "set by user" ? " (set by you)" : ""}</td><td>${x.identity != null ? `${x.identity}% identity, ${x.coverage}% coverage to ${esc(x.reference)}` : "-"}</td></tr>`).join("") +
      `</table></div><p class="muted small" style="margin-top:8px">${a.input.records} sequence record(s), ${fmt(a.input.total_bp)} bp.</p></details>`;
  }
  if (p && p.novel) h += `<div class="info">This septatype is new for ${esc(SHORT[r.species])}. If you have the laboratory phenotype, <a href="#/submit">submit the isolate</a> so curators can add it.</div>`;
  h += (extra || "") + `<div class="actions no-print"><button class="btn line" data-act="print" type="button">Print or save as PDF</button><button class="btn line" data-act="csv" type="button">Download CSV</button><button class="btn line" data-act="json" type="button">Download JSON</button></div><p class="muted small">Research use only. Confirm with laboratory testing.</p></div>`;
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
    $("#m-batch").style.display = PMODE === "batch" ? "" : "none";
    $("#m-bvbrc").style.display = PMODE === "bvbrc" ? "" : "none";
    $("#m-gene").style.display = (PMODE === "file" || PMODE === "paste") ? "" : "none";
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
  const bd = $("#bdrop"), bf = $("#bfile");
  const setB = fs => { PFILES = [...fs].slice(0, 10); $("#bnames").textContent = PFILES.length ? PFILES.map(f => f.name).join(", ") + (fs.length > 10 ? " (first 10 used)" : "") : "One genome per file. Results appear in one table you can download."; };
  bd.onclick = () => bf.click(); bd.onkeydown = e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); bf.click(); } };
  bf.onchange = () => setB(bf.files);
  ["dragover", "dragenter"].forEach(ev => bd.addEventListener(ev, e => { e.preventDefault(); bd.classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => bd.addEventListener(ev, e => { e.preventDefault(); bd.classList.remove("over"); }));
  bd.addEventListener("drop", e => setB(e.dataTransfer.files));
  $("#ex").onclick = () => { ["L", "N", "I", "E", "L"].forEach((v, i) => $("#" + ["r83", "r87", "r80", "r84", "r346"][i]).value = v); $("#mqb").checked = false; $("#mqs").checked = false; };
  $("#pf").onsubmit = async e => {
    e.preventDefault();
    const err = $("#perr"), btn = $("#go"); err.innerHTML = "";
    btn.disabled = true; const old = btn.textContent; btn.textContent = "Analysing";
    if (PMODE === "bvbrc") {
      try { await runBvbrc($("#sp").value, $("#hb").value); if (innerWidth < 860) $("#presult").scrollIntoView({ behavior: "smooth" }); }
      catch (ex) { $("#presult").innerHTML = ""; err.innerHTML = `<div class="err">${esc(ex.message)}</div>`; }
      finally { btn.disabled = false; btn.textContent = old; }
      return;
    }
    if (PMODE === "batch") {
      try { await runBatch($("#sp").value, $("#hb").value); if (innerWidth < 860) $("#presult").scrollIntoView({ behavior: "smooth" }); }
      catch (ex) { $("#presult").innerHTML = ""; err.innerHTML = `<div class="err">${esc(ex.message)}</div>`; }
      finally { btn.disabled = false; btn.textContent = old; }
      return;
    }
    $("#presult").innerHTML = `<div class="band muted">${PMODE === "manual" ? "Looking up the septatype" : "Compressing, uploading and reading the sequence. Usually a few seconds."}</div>`;
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
        if (PMODE === "file") { if (!PFILE) throw new Error("Choose a FASTA file first."); if (PFILE.size > META.max_upload_mb * 1048576) throw new Error(`File is larger than ${META.max_upload_mb} MB.`); fd.append("file", await maybeGzip(PFILE)); }
        else { if (!$("#seq").value.trim()) throw new Error("Paste a nucleotide sequence first."); fd.append("sequence", $("#seq").value); }
        res = await api("/api/v1/predict", { method: "POST", body: fd });
      }
      res._name = PMODE === "manual" ? "manual entry" : PMODE === "paste" ? "pasted sequence" : (PFILE && PFILE.name) || "";
      $("#presult").innerHTML = renderResult(res);
      if (innerWidth < 860) $("#presult").scrollIntoView({ behavior: "smooth" });
    } catch (ex) { $("#presult").innerHTML = ""; err.innerHTML = `<div class="err">${esc(ex.message)}</div>`; }
    finally { btn.disabled = false; btn.textContent = old; }
  };
};

/* ------------------------------------------------------------------ browse */
let BS = { page: 1, sort: "genome_id", desc: false }, BTOT = 0;
let PER = 100;
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
    const arrow = k => BS.sort === k ? (BS.desc ? " &#9660;" : " &#9650;") : "";
    const th = (k, t, rs = "") => `<th class="sort" data-k="${k}"${rs ? ` rowspan="2"` : ""}>${t}${arrow(k)}</th>`;
    const pubLink = u => u ? `<a href="${esc(u)}" target="_blank" rel="noopener">${esc(u.replace(/^https?:\/\//, ""))}</a>` : "-";
    $("#btable").innerHTML = `<thead><tr>${th("genome_id", "Genome ID", 1)}${th("species", "Species and source", 1)}${th("country", "Country", 1)}${th("septatype", "Septatype", 1)}${th("pred", "Ciprotyper predicted phenotype", 1)}<th class="grp" colspan="2">BV-BRC Phenotype Info</th>${th("submitter", "Submitter", 1)}${th("date", "Date", 1)}${th("affiliation", "Affiliation", 1)}<th rowspan="2">Publication</th></tr>
      <tr>${th("lab", "Lab")}${th("bv", "Computationally predicted")}</tr></thead><tbody>` +
      (d.rows.length ? d.rows.map(r => `<tr>
        <td class="nw">${r.bvbrc_url ? `<a href="${esc(r.bvbrc_url)}" target="_blank" rel="noopener">${esc(r.genome_id)}</a>` : esc(r.genome_id)}<span class="alias">${esc(r.genome_name || r.isolate_name)}</span></td>
        <td class="nw"><i>${esc(SHORT[r.species])}</i><span class="alias">${esc(SITE[r.habitat])}${r.source === "community" ? ", community" : ""}</span></td>
        <td class="nw">${esc(r.country || "-")}</td>
        <td class="mono nw">${mutSep(r.septatype)}</td>
        <td class="nw">${r.septatyper_call ? badge(r.septatyper_call, 1) : "-"} ${r.star ? `<span class="star" title="Ciprotyper matches the laboratory phenotype. The BV-BRC computational prediction does not.">&#9733;</span>` : ""}</td>
        <td>${badge(r.lab_phenotype, 1)}</td>
        <td>${r.bvbrc_call ? badge(r.bvbrc_call, 1) : `<span class="muted" title="Not loaded yet">-</span>`}</td>
        <td class="nw">${esc(r.submitter || "-")}</td>
        <td class="nw">${esc(r.date_display || "-")}</td>
        <td class="wrap">${esc(r.affiliation || "-")}</td>
        <td class="nw">${pubLink(r.publication)}</td></tr>`).join("")
        : `<tr><td colspan="11" class="muted">No isolates match these filters.</td></tr>`) + "</tbody>";
    const wrap = $("#bt-wrap"), tr1 = $("#btable thead tr");
    if (tr1) wrap.style.setProperty("--h1", tr1.offsetHeight + "px");
    wrap.scrollTop = 0;
    if (BS.jump) { BS.jump = false; wrap.scrollIntoView({ block: "start" }); }       // after a page or size change, bring the table back into view
    $$("#btable th.sort").forEach(h => h.onclick = () => { BS.desc = BS.sort === h.dataset.k ? !BS.desc : false; BS.sort = h.dataset.k; BS.page = 1; loadBrowse(); });
  } catch (e) { $("#btable").innerHTML = `<tr><td class="err">${esc(e.message)}</td></tr>`; }
}
async function bvNote() {
  let st; try { st = await api("/api/v1/bvbrc-status"); } catch (e) { st = { loaded: false }; }
  META.bvLoaded = st.loaded; $("#f-wins").disabled = !st.loaded;
  const n = $("#bv-note");
  if (st.loaded) { n.innerHTML = ""; return; }
  if (st.running) { n.innerHTML = `<div class="info">Fetching BV-BRC predictions: ${st.done} of ${st.total} genomes. This table updates when it finishes.</div>`; setTimeout(async () => { await bvNote(); if (META.bvLoaded) loadBrowse(); }, 6000); }
  else if (st.error) n.innerHTML = `<div class="info">The BV-BRC comparison is not available right now (${esc(st.error)}). It retries when the server restarts.</div>`;
  else n.innerHTML = `<div class="info">The BV-BRC comparison is being prepared and will appear here soon.</div>`;
}
INIT.browse = async () => {
  await refreshMeta();
  const f = META.facets;
  $("#f-hb").innerHTML = opts(f.habitat.map(([k, n]) => [k, `${SITE[k] || k} (${fmt(n)})`]), "All");
  $("#f-sp").innerHTML = opts(Object.entries(SP), "All");
  $("#dl-st").innerHTML = dl(f.st);
  $("#f-country").innerHTML = opts(withN(f.country), "All");
  bvNote();
  let t; const reload = () => { BS.page = 1; loadBrowse(); };
  $$("#p-browse select:not(#b-per), #p-browse input[type=number]:not(#b-per-custom), #f-wins, #f-st").forEach(e => e.onchange = reload);
  $$("#p-browse input[type=search]:not(#f-st)").forEach(e => e.oninput = () => { clearTimeout(t); t = setTimeout(reload, 350); });
  $("#b-prev").onclick = () => { BS.page--; BS.jump = true; loadBrowse(); }; $("#b-next").onclick = () => { BS.page++; BS.jump = true; loadBrowse(); };
  const setPer = v => { const n = Math.max(1, Math.min(1000, parseInt(v, 10) || 100)); if (n !== PER) { PER = n; BS.page = 1; BS.jump = true; loadBrowse(); } };
  $("#b-per").onchange = e => { const c = e.target.value === "custom"; if (!c) $("#b-per-custom").value = ""; $("#b-per-custom").style.display = c ? "" : "none"; if (c) $("#b-per-custom").focus(); else setPer(e.target.value); };
  $("#b-per-custom").onchange = e => { if (e.target.value) setPer(e.target.value); };
  $("#f-toggle").onclick = () => { const band = $("#f-toggle").closest(".band"), hid = band.classList.toggle("fhide"); $("#f-toggle").textContent = hid ? "Show filters" : "Hide filters"; };
  $("#f-reset").onclick = () => { $$("#p-browse .filters select, #p-browse .filters input").forEach(e => e.value = ""); $("#f-wins").checked = false; reload(); };
  loadBrowse();
};

/* --------------------------------------------------------------- analytics */
async function loadBars() {
  const d = await api("/api/v1/stats/septatypes?species=" + $("#a-sp").value + "&top=12");
  const max = Math.max(...d.rows.map(r => Math.max(...Object.values(r.by_habitat))), 1);
  const present = META.sources.map(s => s.key).filter(k => d.rows.some(r => k in r.by_habitat));
  $("#a-legend").innerHTML = present.map(k => `<span><i style="background:${hcol(k)}"></i>${esc(SITE[k] || k)}</span>`).join("");
  $("#a-bars").innerHTML = d.rows.map(r => `<div class="hb"><span class="mono">${mutSep(r.septatype)}</span><div class="bars">${Object.entries(r.by_habitat).map(([h, n]) => `<div class="b" style="width:${Math.max(1, 100 * n / max)}%;background:${hcol(h)}" title="${esc(SITE[h] || h)}: ${n}"></div>`).join("")}</div></div>`).join("") +
    `<p class="muted" style="font-size:.85rem;margin-top:8px">Top 12 septatypes of ${fmt(d.total)} isolates. Bar length is the number of isolates from each isolation source.</p>`;
}
let plotlyP = null;
function loadPlotly() {
  if (window.Plotly) return Promise.resolve();
  return plotlyP || (plotlyP = new Promise((ok, no) => { const s = document.createElement("script"); s.src = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"; s.onload = ok; s.onerror = () => no(new Error("map library")); document.head.appendChild(s); }));
}
async function loadMap() {
  const p = new URLSearchParams(); [["species", "#g-sp"], ["st", "#g-st"], ["septatype", "#g-sep"], ["phenotype", "#g-ph"]].forEach(([k, s]) => { if ($(s).value) p.set(k, $(s).value); });
  const d = await api("/api/v1/stats/geography?" + p);
  const tot = d.rows.reduce((a, r) => a + r.n, 0), res = d.rows.reduce((a, r) => a + r.res, 0);
  $("#g-sum").textContent = d.rows.length || d.no_country ? `${d.rows.length} countries, ${fmt(tot)} isolates with a country (${fmt(res)} resistant)${d.no_country ? `, plus ${d.no_country} with no country recorded: ${fmt(tot + d.no_country)} in all.` : "."}` : "No isolates match these filters.";
  $("#g-table").innerHTML = `<tr><th>Country</th><th>Isolates</th><th>Resistant</th><th>% resistant</th></tr>` + d.rows.map(r => `<tr><td>${esc(r.country)}</td><td class="stat">${fmt(r.n)}</td><td class="stat">${fmt(r.res)}</td><td class="stat">${Math.round(100 * r.res / r.n)}%</td></tr>`).join("") + (d.no_country ? `<tr><td class="muted">No country recorded</td><td class="stat">${fmt(d.no_country)}</td><td class="muted">-</td><td class="muted">-</td></tr>` : "");
  try {
    await loadPlotly();
    const dark = false;
    Plotly.react("map", [{ type: "choropleth", locations: d.rows.map(r => r.iso3), z: d.rows.map(r => r.n), text: d.rows.map(r => `${r.country}: ${r.n} isolates, ${r.res} resistant`), hoverinfo: "text",
      colorscale: [[0, "#e8f0fa"], [1, "#0a2540"]], marker: { line: { color: "#8aa0b8", width: .4 } }, colorbar: { title: "Isolates", thickness: 12 } }],
      { margin: { t: 0, b: 0, l: 0, r: 0 }, paper_bgcolor: "rgba(0,0,0,0)", geo: { showframe: false, showcoastlines: false, projection: { type: "natural earth" }, bgcolor: "rgba(0,0,0,0)", showcountries: true, countrycolor: dark ? "#26364a" : "#d5dee9" } },
      { displayModeBar: false, responsive: true });
  } catch (e) { $("#map").innerHTML = `<p class="muted">The map library could not be loaded. The table below has the same counts.</p>`; }
}
INIT.analytics = async () => {
  await refreshMeta();
  $("#a-sp").innerHTML = opts(Object.entries(SP)); $("#a-sp").onchange = loadBars; loadBars();
  $("#g-sp").innerHTML = opts(Object.entries(SP), "All");
  $("#dl-st2").innerHTML = dl(META.facets.st);
  $("#dl-sep").innerHTML = dl(META.facets.septatype);
  $$("#p-analytics .filters select, #g-st, #g-sep").forEach(e => e.onchange = loadMap); loadMap();
};

/* ------------------------------------------------------------------ submit */
INIT.submit = () => {
  const setDate = () => { $("#s-date").value = new Date().toLocaleDateString("en-GB", { day: "numeric", month: "long", year: "numeric" }); };
  setDate();
  $("#s-sp").innerHTML = opts(Object.entries(SP));
  $("#dl-country").innerHTML = META.countries.map(c => `<option value="${esc(c)}"></option>`).join("");
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
      [["species", "#s-sp"], ["habitat", "#s-hb"], ["genome_id", "#s-gid"], ["country", "#s-country"], ["year", "#s-year"], ["lab_phenotype", "#s-ph"], ["submitter", "#s-submitter"], ["affiliation", "#s-aff"], ["publication", "#s-pub"], ["email", "#s-email"]].forEach(([k, s]) => fd.append(k, $(s).value.trim()));
      fd.append("genome_name", $("#s-name").dataset.name || ""); fd.append("file", await maybeGzip(f));
      const r = await api("/api/v1/submissions", { method: "POST", body: fd });
      const msg = { needs_review: "Needs review: the septatype is new or the prediction disagrees with your laboratory result.", pending: "Pending: the prediction matches your result or no result was given. A curator will check it." }[r.status] || r.status;
      $("#s-result").innerHTML = renderResult(r.result, `<div class="info"><b>Submission ${esc(r.job_id)}</b><br>${esc(msg)}${r.matches === null ? " Without a laboratory phenotype the isolate cannot be approved." : ""}</div>`);
      $("#sf").reset(); setDate(); $("#s-name").textContent = ""; $("#s-name").dataset.name = ""; upd();
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
    $("#t-kv").innerHTML = `<dt>Tetratype</dt><dd class="mono">${aas.join("-")}</dd><dt>Grammar call</dt><dd>${badge(p.grammar.call, 1)} ${esc(p.grammar.rule)}</dd><dt>QRDR helix mean hydrophobicity &lt;H&gt;</dt><dd class="stat">${p.mean_h}</dd>`;
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
    let h = `<h3 style="font-size:1rem">Majority rule and grammar against laboratory phenotype</h3><div class="scroll"><table><tr><th>Species</th><th>Source</th><th>Isolates</th><th>Majority rule</th><th>Grammar only</th></tr>${rows.map(x => `<tr><td><i>${esc(SHORT[x.species])}</i></td><td>${esc(SITE[x.habitat])}</td><td class="stat">${fmt(x.n_isolates)}</td><td class="stat">${x.accuracy}%</td><td class="stat">${x.grammar_accuracy}%</td></tr>`).join("")}</table></div><p class="muted" style="font-size:.85rem">Both are scored on the founding isolates the rules were built from. They describe the reference set, not independent performance. The grammar only calls septatypes missing from the database.</p>`;
    h += `<h3 style="font-size:1rem;margin-top:18px">Head to head with BV-BRC computational predictions</h3>`;
    if (b.loaded) h += `<p class="stat">${fmt(b.compared)} isolates have a BV-BRC computational prediction. Ciprotyper matches the lab result in ${fmt(b.septatyper_ok)} (${(100 * b.septatyper_ok / b.compared).toFixed(1)}%) and BV-BRC in ${fmt(b.bvbrc_ok)} (${(100 * b.bvbrc_ok / b.compared).toFixed(1)}%). Ciprotyper matches the laboratory phenotype where BV-BRC does not in ${fmt(b.wins)} isolates. The reverse happens in ${fmt(b.losses)}.</p><div class="scroll"><table><tr><th>Genome</th><th>Septatype</th><th>Lab</th><th>BV-BRC</th><th>Ciprotyper</th></tr>${b.examples.map(x => `<tr><td>${esc(x.genome_id)}<span class="alias">${esc(x.isolate_name)}</span></td><td class="mono">${esc(x.septatype)}</td><td>${badge(x.lab_phenotype, 1)}</td><td>${badge(x.bvbrc_call, 1)}</td><td>${badge(x.septatyper_call, 1)}</td></tr>`).join("")}</table></div>`;
    else h += `<div class="info">BV-BRC predictions have not been loaded yet, so the head-to-head comparison is empty. They load automatically when the server can reach BV-BRC.</div>`;
    h += `<h3 style="font-size:1rem;margin-top:18px">Susceptible tetratype S-D-S-E, split by qnr</h3><div class="scroll"><table><tr><th>Species</th><th>qnrB</th><th>qnrS</th><th>Isolates</th><th>Resistant</th></tr>${s.edgecases.wildtype_tetratype_by_qnr.map(x => `<tr><td><i>${esc(SHORT[x.species])}</i></td><td>${x.qnrB2 ? "yes" : "no"}</td><td>${x.qnrS1 ? "yes" : "no"}</td><td class="stat">${x.n}</td><td class="stat">${x.res} (${Math.round(100 * x.res / x.n)}%)</td></tr>`).join("")}</table></div><p class="muted" style="font-size:.85rem">A qnr gene on a susceptible baseline is usually resistant, but not always. A gene found in sequence may not be expressed.</p>`;
    $("#t-valid").innerHTML = h;
    $("#t-aa").innerHTML = `<tr><th>Residue</th><th>Name</th><th>Hydrophobicity</th><th>Class</th></tr>` + aa.rows.sort((a, b) => b.hydrophobicity - a.hydrophobicity).map(x => `<tr><td class="mono">${x.aa}</td><td>${esc(x.name)}</td><td class="stat">${x.hydrophobicity.toFixed(2)}</td><td>${esc(x.class)}</td></tr>`).join("");
  } catch (e) { $("#t-valid").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
};

/* ------------------------------------------------------------------- calculation */
INIT.calc = async () => {
  const qs = new URLSearchParams(location.hash.split("?")[1] || "");
  const clean = v => (v || "").toUpperCase().replace(/[^A-Z]/g, "").slice(0, 4);
  const start = clean(qs.get("r")) || "LNIE", startRef = clean(qs.get("ref")) || "SDSE";
  let aa, cs;
  try { [aa, cs] = await Promise.all([api("/api/v1/aminoacids"), api("/api/v1/stats/cutoff")]); }
  catch (e) { $("#calc-body").innerHTML = `<div class="band"><div class="err">${esc(e.message)}</div></div>`; return; }
  const rows = aa.rows.slice().sort((a, b) => a.hydrophobicity - b.hydrophobicity);
  const wt = cs.wild_type, rs = cs.resistant_reference, cut = cs.cutoff;
  const ctab = () => `<div class="scroll"><table class="step"><tr><th class="num">Cutoff</th><th class="num">Whole grammar, all isolates</th><th class="num">&lt;H&gt; alone, all isolates</th><th class="num">&lt;H&gt; alone, where it decides</th></tr>${cs.grammar_all.map((x, i) => `<tr${x.cutoff === cut ? ' class="hit-row"' : ""}><td class="num">${x.cutoff.toFixed(x.cutoff === cut ? 3 : 2)}${x.cutoff === cut ? " (used by the site)" : x.cutoff === wt.mean ? " (wild type)" : ""}</td><td class="num">${x.accuracy}%</td><td class="num">${cs.all[i].accuracy}%</td><td class="num">${cs.single_vote[i].accuracy}%</td></tr>`).join("")}</table></div>`;
  $("#calc-body").innerHTML = `
  <div class="band"><h3>1. The hydrophobicity scale</h3>
    <p>Each amino acid has a hydrophobicity value H from the scale in the paper (the Eisenberg normalised consensus scale). Negative means hydrophilic, positive means hydrophobic.</p>
    <div class="scroll"><table class="step"><tr><th>Amino acid</th><th>Letter</th><th class="num">H</th><th>Class</th></tr>${rows.map(x => `<tr><td>${esc(x.name)}</td><td class="mono">${x.aa}</td><td class="num">${x.hydrophobicity.toFixed(2)}</td><td>${esc(x.class)}</td></tr>`).join("")}</table></div></div>
  <div class="band"><h3>2. Calculate the QRDR helix &lt;H&gt;</h3>
    <p>&lt;H&gt; is read over the whole quinolone resistance-determining region (QRDR) helix of both enzymes, not only over the four tetratype positions. The helices are fixed, and your four residues replace the marked positions:</p>
    <p class="formula">GyrA 82-93 &nbsp; D(S)AVY(D)TIVRMA &nbsp; 12 residues, S = GyrA-83, D = GyrA-87</p>
    <p class="formula">ParC 80-90 &nbsp; (S)ACY(E)AMVLMA &nbsp; 11 residues, S = ParC-80, E = ParC-84</p>
    <ol class="small"><li>Add the H values of the 12 GyrA helix residues and divide by 12. That is the GyrA helix mean.</li><li>Add the H values of the 11 ParC helix residues and divide by 11. That is the ParC helix mean.</li><li>&lt;H&gt; = (GyrA helix mean + ParC helix mean) / 2.</li></ol>
    <p class="muted small">Type four residues in the order GyrA-83, GyrA-87, ParC-80, ParC-84. Rows in amber are those four positions.</p>
    <div class="row2"><div><label for="calc-in">Sample residues</label><input type="text" id="calc-in" maxlength="4" value="${esc(start)}" class="mono-in"></div>
    <div><label for="calc-ref">Reference residues</label><input type="text" id="calc-ref" maxlength="4" value="${esc(startRef)}" class="mono-in"></div></div>
    <div id="calc-steps"></div></div>
  <div class="band" id="calc-cutoff"><h3>3. Where the cutoff comes from</h3>
    <p><b>The cutoff is the midpoint between two reference helices.</b> The susceptible wild-type residues ${esc(wt.residues.split("").join("-"))} give a GyrA helix mean of ${f3(wt.gyr.mean)} and a ParC helix mean of ${f3(wt.parc.mean)}, so &lt;H&gt; = <b>${f3(wt.mean)}</b>. The resistant reference residues ${esc(rs.residues.split("").join("-"))} (the commonest fully changed tetratype) give ${f3(rs.gyr.mean)} and ${f3(rs.parc.mean)}, so &lt;H&gt; = <b>${f3(rs.mean)}</b>. The midpoint is (${f3(wt.mean)} + ${f3(rs.mean)}) / 2 = <b>${f3(cut)}</b>. When the four residues are not all in one column of the table and exactly one position votes Resistant, the sample counts as Resistant only if its &lt;H&gt; is above this value. The helices are mostly oily residues, so both values are above zero. The cutoff is set from these two reference helices, not by searching the isolates for the best score.</p>
    <p><b>How much it matters.</b> &lt;H&gt; settles a call only where the four residues are not all in one column of the table and exactly one position votes Resistant. That is ${fmt(cs.n_single_vote)} of the ${fmt(cs.n_isolates)} reference isolates (${fmt(cs.single_vote_resistant)} resistant). The whole grammar scores ${cs.grammar_all.find(x => x.cutoff === wt.mean).accuracy}% with the wild-type value ${f3(wt.mean)} as cutoff and ${cs.grammar_all.find(x => x.cutoff === cut).accuracy}% with ${f3(cut)}.</p>
    <p><b>Agreement with the paper.</b> The paper reports &lt;H&gt; of 0.27 to 0.32 for the susceptible tetratypes (S-D-S-E, S-N-S-E, Y-D-S-E, S-Y-S-E, L-D-S-E) and 0.40 to 0.48 for the resistant ones (L-N-I-E, F-N-I-E, I-D-I-E, I-N-I-E, L-G-I-E, I-Y-I-E, F-A-I-E, L-N-I-V), with F-N-S-K an exception at 0.297. This calculation gives the same values (for example S-D-S-E ${f3(wt.mean)}, L-N-I-E ${f3(rs.mean)}). The cutoff ${f3(cut)} lies in the gap between the two groups, above L-D-S-E (0.322) and F-D-S-E (0.327).</p>
    <p><b>Check against the data.</b> Across the ${fmt(cs.n_isolates)} reference isolates the mean &lt;H&gt; is <b>${f3(cs.mean_resistant)}</b> for resistant and <b>${f3(cs.mean_susceptible)}</b> for susceptible isolates. The values run from ${f3(cs.range[0])} to ${f3(cs.range[1])}. The two groups differ on average, but only by a small step.</p>
    <p>Accuracy at each cutoff, scored on the reference isolates the rules were built from. The first column is the whole grammar (the three checkpoints, with &lt;H&gt; settling single votes). The other two use &lt;H&gt; alone. The ${fmt(cs.n_single_vote)} isolates where &lt;H&gt; decides (${fmt(cs.single_vote_resistant)} resistant) are the last column.</p>${ctab()}
    <p style="margin-top:12px">Where &lt;H&gt; decides it is a weak separator (best ${cs.best_single_vote.accuracy}% on this grid).${cs.h_decided_max !== null && cs.h_decided_max <= cut ? " At the site's cutoff every isolate it decides is called Susceptible, so on these data it works as a last tie-breaker toward Susceptible." : ""}</p>
    <p class="muted" style="font-size:.9rem;margin-top:8px">The grammar is only a fallback. When a septatype was seen in the reference isolates, the majority rule is used instead. These tables describe the reference set. They are not a validation.</p></div>
  <div class="band"><h3>4. The three checkpoints</h3>
    <p><b>Checkpoint 1.</b> Context-sensitive grammar of GyrA-83, GyrA-87, ParC-80 and ParC-84. The isolate is <b>susceptible</b> if all four residues are in the Susceptible column, and <b>resistant</b> if all four are in the Resistant column.</p>
    <div class="scroll"><table class="rtab"><tr><th>Position</th><th>Susceptible</th><th>Resistant</th></tr>${SYNTAX.map(([p, a, b]) => `<tr><td class="nw"><b>${p}</b></td><td>${a.split(" ").map(x => `<span class="rchip">${x}</span>`).join("")}</td><td>${b.split(" ").map(x => `<span class="rchip">${x}</span>`).join("")}</td></tr>`).join("")}</table></div>
    <p class="muted small" style="margin-top:6px">A residue can be in both columns (D at GyrA-87, E at ParC-84), so the four positions are read together, which is the context sensitivity. One exception counts as Resistant: the ParC pair S-K (S at ParC-80 together with K at ParC-84, the F-N-S-K tetratype). S alone is the susceptible S-E state, and K alone is not listed. In a result the exception appears in the table only when your sample carries the whole pair. The residues follow the paper's Table S6 and Fig. 6.</p>
    <p><b>Checkpoint 2.</b> RecN-346. Susceptible if it is leucine (L) encoded by cta. Resistant if it is leucine encoded by ctg, or phenylalanine (F). Other codons, and an unknown codon, give no resistant signal.</p>
    <p><b>Checkpoint 3.</b> PMQR genes. Susceptible if both qnrB2 and qnrS1 are absent. Resistant if either or both are present.</p>
    <p><b>Combined.</b> Susceptible if all checkpoints say susceptible. Resistant if any checkpoint says resistant.</p>
    <p class="muted" style="font-size:.9rem">The table gives no answer when the four residues are not all in one column (for example L-D-S-E, where L at GyrA-83 is in the Resistant column but S at ParC-80 is in the Susceptible column). There the residue votes decide (two or more positions outside the Susceptible column: Resistant; none: Susceptible), and one vote is settled by &lt;H&gt; as in section 3. This is an assumption added by the site. <i>S. sonnei</i> also keeps one extra rule from the reference isolates: a single change at GyrA-83 or GyrA-87 counts as Resistant.</p></div>`;
  const htab = b => `<div class="scroll"><table class="step"><tr><th>${esc(b.gene)} helix ${b.start}-${b.end}</th><th>Residue</th><th class="num">H</th></tr>${b.values.map(x => `<tr${x.sub ? ' class="hit-row"' : ""}><td>${esc(x.position)}</td><td class="mono">${esc(x.residue)}</td><td class="num">${f2(x.h)}</td></tr>`).join("")}<tr><td colspan="2">Sum of ${b.n}</td><td class="num">${f2(b.sum)}</td></tr><tr><td colspan="2"><b>${esc(b.gene)} mean = sum / ${b.n}</b></td><td class="num"><b>${f3(b.mean)}</b></td></tr></table></div>`;
  const side = (t, v) => `<div><h3 style="font-size:1rem">${t}</h3>${htab(v.gyr)}${htab(v.parc)}<p class="formula wrap">&lt;H&gt; = (${f3(v.gyr.mean)} + ${f3(v.parc.mean)}) / 2 = <b>${f3(v.mean)}</b></p></div>`;
  const draw = async () => {
    const r = $("#calc-in").value.toUpperCase(), rf = $("#calc-ref").value.toUpperCase(), el = $("#calc-steps");
    const bad = v => v.length !== 4 || /[^ACDEFGHIKLMNPQRSTVWY]/.test(v);
    if (bad(r) || bad(rf)) { el.innerHTML = `<p class="muted">Type four amino-acid letters in each box.</p>`; return; }
    try {
      const c = await api(`/api/v1/grammar?residues=${r}&reference=${rf}`), hy = c.hydrophobicity, rr = c.reference;
      el.innerHTML = `<div class="row2">${side("Sample", hy)}${side("Reference", rr)}</div>${hChart(hy.mean, rr.mean, hy.cutoff, { res: cs.mean_resistant, sus: cs.mean_susceptible })}
        <p>Sample &lt;H&gt; is ${hy.side === "above" ? "above" : "at or below"} the cutoff of ${f3(hy.cutoff)}, which points to <b>${hy.side === "above" ? "resistant" : "susceptible"}</b>. ${c.checkpoints.residues.resolved_by === "hydrophobicity" ? "Here it settles checkpoint 1, because the four residues are not all in one column and exactly one position votes Resistant." : "It does not settle this sample: " + (c.checkpoints.residues.pattern ? "the residue classes below do." : esc(c.checkpoints.residues.rule) + ".")}</p>
        <p><b>Checkpoint 1 says ${esc(c.checkpoints.residues.call)}.</b> ${c.checkpoints.residues.pattern ? `All four residues are in the ${c.checkpoints.residues.pattern === "resistant" ? "Resistant" : "Susceptible"} column.` : "The four residues are not all in one column."}</p>${critTable(c.positions)}`;
    } catch (e) { el.innerHTML = `<div class="err">${esc(e.message)}</div>`; }
  };
  $("#calc-in").oninput = draw; $("#calc-ref").oninput = draw;
  await draw();
  if (qs.get("s") === "cutoff") $("#calc-cutoff").scrollIntoView({ behavior: "smooth" }); else window.scrollTo(0, 0);
};

/* -------------------------------------------------------------- data, about */
INIT.data = () => {
  const o = location.origin;
  $("#ex-curl").textContent = `# Resistant E. coli urine isolates, first page\ncurl "${o}/api/v1/isolates?species=ecoli&habitat=URINE&phenotype=Resistant"\n\n# Predict from a genome file\ncurl -F file=@genome.fasta -F species=kpneumoniae -F habitat=URINE ${o}/api/v1/predict`;
  $("#ex-py").textContent = `import requests\nbase = "${o}/api/v1"\nr = requests.post(f"{base}/predict", data={"species": "ecoli", "habitat": "URINE"},\n                  files={"file": open("genome.fasta", "rb")}).json()\nprint(r["septatype"], r["analysis"]["mlst"]["st"], r["prediction"]["call"])`;
};
INIT.about = () => { $("#about-build").textContent = `Rules built ${META.build.built || ""} from ${META.build.source || ""} (${META.build.n_isolates || ""} isolates).`; };

/* ----------------------------------------------------------------- curation */
const tokenOk = () => {
  const t = $("#c-token").value.trim();
  if (!t) { $("#c-err").innerHTML = `<div class="err">Type the curator token first.</div>`; return null; }
  if (!/^[\x20-\x7E]+$/.test(t)) { $("#c-err").innerHTML = `<div class="err">The token has a character the browser cannot send. Use plain letters and numbers.</div>`; return null; }
  return t;
};
INIT.curate = () => {
  if (!META.curation_enabled) $("#c-err").innerHTML = `<div class="note">Curation is not switched on for this server. In Render, open the service, then Environment, add ADMIN_TOKEN with a value, save, and wait for the redeploy.</div>`;
  const TABS = [["needs_review", "Needs review"], ["pending", "Pending"], ["approved", "Approved"], ["rejected", "Rejected"], ["", "All"]];
  const LABEL = { needs_review: "Needs review", pending: "Pending", approved: "Approved", rejected: "Rejected" };
  let status = "needs_review";
  const load = async () => {
    $("#c-err").innerHTML = ""; const tk = tokenOk(); if (!tk) return;
    try {
      const d = await api("/api/v1/admin/queue?status=" + status, { headers: { "X-Admin-Token": tk } });
      const tabs = $("#c-tabs"); tabs.style.display = "";
      tabs.innerHTML = TABS.map(([k, t]) => `<button type="button" data-s="${k}" class="${k === status ? "on" : ""}">${t} (${d.counts[k || "all"] ?? 0})</button>`).join("");
      $$("#c-tabs button").forEach(b => b.onclick = () => { status = b.dataset.s; load(); });
      const actions = r => {
        const note = `<input type="text" placeholder="note (optional)" data-n="${esc(r.job_id)}" style="margin-bottom:6px">`;
        const btn = (a, t, cls) => `<button class="btn ${cls}" data-a="${a}" data-j="${esc(r.job_id)}" data-st="${esc(r.status)}">${t}</button>`;
        if (r.status === "approved") return `${note}${btn("reject", "Reject and remove", "line")}<span class="alias">In the database as ${esc(r.genome_id || "SUB-" + r.job_id)}</span>`;
        if (r.status === "rejected") return `${note}${btn("approve", "Approve instead", "solid")}`;
        return `${note}${btn("approve", "Approve", "solid")} ${btn("reject", "Reject", "line")}`;
      };
      $("#c-table").innerHTML = `<tr><th>Submission</th><th>Isolate</th><th>Septatype</th><th>ST</th><th>Prediction</th><th>Lab</th><th>Status</th><th>Decision</th></tr>` +
        (d.rows.length ? d.rows.map(r => `<tr><td class="nw">${esc(r.job_id)}<span class="alias">submitted ${esc(r.created)}</span></td>
          <td>${esc(r.genome_id || "no ID")}<span class="alias">${esc(SHORT[r.species])}, ${esc(SITE[r.habitat])}, ${esc(r.country || "-")} ${esc(r.year || "")}</span>
            <span class="alias">${esc(r.submitter || "-")}${r.affiliation ? ", " + esc(r.affiliation) : ""}${r.email ? " (" + esc(r.email) + ")" : ""}</span>${r.publication ? `<span class="alias"><a href="${esc(r.publication)}" target="_blank" rel="noopener">${esc(r.publication.replace(/^https?:\/\//, ""))}</a></span>` : ""}</td>
          <td class="mono nw">${mutSep(r.septatype)}${r.novel ? `<span class="alias">new septatype</span>` : ""}</td>
          <td class="nw">${esc(r.st || "-")}</td>
          <td>${r.call ? badge(r.call, 1) : "-"}<span class="alias">${esc(r.tier || "")}</span></td>
          <td>${r.lab_phenotype ? badge(r.lab_phenotype, 1) : "-"}${r.matches === 0 ? `<span class="alias">disagrees</span>` : ""}</td>
          <td>${esc(LABEL[r.status] || r.status)}${r.decided ? `<span class="alias">decided ${esc(r.decided)}</span>` : ""}${r.curator_note ? `<span class="alias">${esc(r.curator_note)}</span>` : ""}</td>
          <td>${actions(r)}</td></tr>`).join("")
          : `<tr><td colspan="8" class="muted">Nothing in this list.</td></tr>`);
      $$("#c-table button[data-a]").forEach(b => b.onclick = async () => {
        if (b.dataset.st === "approved" && !confirm("Reject this approved submission? Its isolate is removed from the repository and the rules.")) return;
        try {
          const res = await api("/api/v1/admin/decide", { method: "POST", headers: { "Content-Type": "application/json", "X-Admin-Token": $("#c-token").value.trim() },
            body: JSON.stringify({ job_id: b.dataset.j, action: b.dataset.a, note: ($(`[data-n="${b.dataset.j}"]`) || {}).value || "" }) });
          LOADED.browse = false; LOADED.home = false; LOADED.analytics = false;
          if (b.dataset.a === "approve") META.counts.genomes += 1; else if (res.revoked) META.counts.genomes -= 1;
          load();
        } catch (e) { $("#c-err").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
      });
    } catch (e) { $("#c-err").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
  };
  const pstatus = async () => {
    const tk = $("#c-token").value.trim(); if (!tk) return;
    try {
      const s = await api("/api/v1/admin/persist-status", { headers: { "X-Admin-Token": tk } });
      $("#c-pst").textContent = s.github_backup ? `Automatic GitHub backup is on${s.last_saved ? `, last saved ${s.last_saved} UTC` : ", nothing saved yet"}${s.error ? `. Last problem: ${s.error}` : ""}.` : "Automatic GitHub backup is off (set BACKUP_GITHUB_REPO and BACKUP_GITHUB_TOKEN on the server).";
    } catch (e) { /* the sign-in message is shown by the queue */ }
  };
  $("#c-load").onclick = () => { load(); pstatus(); };
  $("#c-bv").onclick = async () => {
    $("#c-err").innerHTML = ""; const tk = tokenOk(); if (!tk) return;
    try { const s = await api("/api/v1/admin/bvbrc/load", { method: "POST", headers: { "X-Admin-Token": tk } });
      $("#c-bvs").textContent = s.loaded ? "Already loaded." : "Started. Check the Browse page in a minute."; }
    catch (e) { $("#c-err").innerHTML = `<div class="err">${esc(e.message)}</div>`; }
  };
};

/* -------------------------------------------------------------------- start */
(async () => {
  try {
    META = await api("/api/v1/meta");
    Object.assign(SP, META.species); useSources();
    LOADED.home = false;
    route();
  } catch (e) {
    $("#main").innerHTML = `<div class="band"><div class="err">${esc(e.message)}</div></div>`;
  }
})();
