/* Browser-side orchestration. Scientific calculations stay in the Python workers. */
const state = {
  analysisGeneration: 0,
  refreshing: false,
  analysisId: null,
  detail: null,
  result: null,
  pairId: null,
  surface: null,
  stage: null,
  structureText: null,
  structureFormat: null,
  structureComponents: {},
  cartoonReps: {},
  meshComponents: {},
  meshReps: {},
  meshComponentChains: new Map(),
  meshGeometry: {},
  meshSourceData: new WeakMap(),
  meshSource: null,
  meshSelectionKey: null,
  meshColorKey: null,
  selectionReps: [],
  pocketComponent: null,
  pocketResult: null,
  pocketHighlighted: false,
  pocketHighlightCancelled: false,
  pocketTimer: null,
  pocketJobId: null,
  pocketGeneration: 0,
  selectedResidues: [],
  activeTab: "residue",
  surfaceScale: "6",
  surfaceOpacity: 0.55,
  surfacePadding: 2,
  surfaceSeparation: 0,
  visibleChains: {},
  contact: {cut: 5, type: "all", vals: true, sasa: true},
  viewerLoaded: false,
  curvatureLoaded: false,
  timer: null,
};

const $ = (id) => document.getElementById(id);
const esc = (value) => String(value ?? "").replace(/[&<>"']/g, (c) => ({
  "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
}[c]));
const aa3to1 = {
  ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
  GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
  PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V", MSE: "M",
  SEC: "C", PYL: "K",
};
const formalCharge = {ASP: -1, GLU: -1, LYS: 1, ARG: 1};
const residueClassGlyph = {polar: "●", hyd: "▲", arom: "◆", spec: "○"};
const residueClassName = {
  pos: "cationic", neg: "anionic", polar: "polar", hyd: "hydrophobic",
  arom: "aromatic", spec: "Gly / Pro / Cys",
};
const distanceBins = [3, 3.5, 4, 5, 6];
const distanceLabels = ["≤3.0", "3.0–3.5", "3.5–4.0", "4.0–5.0", "5.0–6.0"];
const curvatureColors = [
  [0, [142, 13, 19]], [0.25, [239, 106, 97]], [0.5, [244, 244, 241]],
  [0.75, [116, 173, 209]], [1, [33, 102, 172]],
];
const cartoonColors = ["#dfaa25", "#36a58b", "#8b70c7", "#d67439", "#4979b8", "#7c9b5c"];

async function getJSON(url, options) {
  const response = await fetch(url, options);
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(body.detail || body.error || response.statusText);
  return body;
}

function fmt(value, digits = 1) {
  return value == null || Number.isNaN(Number(value)) ? "n/a" : Number(value).toFixed(digits);
}

function showError(error) {
  closeDrawer();
  $("status-message").hidden = false;
  $("status-message").textContent = error.message || String(error);
}

function clearError() {
  $("status-message").hidden = true;
  $("status-message").textContent = "";
}

function residueKey(identity) {
  if (!identity) return "";
  return identity.key || `${identity.chain_id}:${identity.number}${identity.insertion_code || ""}`;
}

function residueLabel(identity) {
  return identity?.label || `${aa3to1[identity?.resname] || "X"}${identity?.number ?? ""}${identity?.insertion_code || ""}`;
}

function residueDisplay(identity) {
  return identity ? `${identity.chain_id}:${residueLabel(identity)}` : "n/a";
}

// Molecule name declared in the coordinate file, with the chain ID kept so
// identical copies (homodimers) stay distinguishable.
function chainLabel(chainId) {
  const chain = (state.result?.metadata?.chains || []).find((item) => item.id === chainId);
  const name = chain?.name;
  return name && name !== `Chain ${chainId}` ? `${name} (${chainId})` : `Chain ${chainId ?? "?"}`;
}

function setStatus(status) {
  $("analysis-status").textContent = status;
  $("analysis-status").className = `status-pill ${esc(status)}`;
}

const activeStatuses = ["queued", "running"];
const statusWords = {queued: "Queued", running: "Running", failed: "Failed", cancelled: "Cancelled", interrupted: "Interrupted"};

function isActiveAnalysis(analysis) {
  return activeStatuses.includes(analysis.status)
    || (analysis.jobs || []).some((job) => activeStatuses.includes(job.status));
}

function relativeTime(iso) {
  const date = new Date(iso);
  const today = new Date();
  const days = Math.round((new Date(today.toDateString()) - new Date(date.toDateString())) / 86400000);
  if (days === 0) return date.toLocaleTimeString([], {hour: "numeric", minute: "2-digit"});
  if (days === 1) return "Yesterday";
  return date.toLocaleDateString([], {day: "numeric", month: "short", ...(date.getFullYear() === today.getFullYear() ? {} : {year: "numeric"})});
}

// One segment per chain, sized by residue count and coloured like that
// chain's cartoon in the 3D view.
function chainStrip(chains) {
  if (!chains?.length) return '<span class="chain-strip pending" aria-hidden="true"><i style="flex-grow:1"></i></span>';
  const description = chains.map((chain) => `${chain.name && chain.name !== `Chain ${chain.id}` ? `${chain.name} (${chain.id})` : `Chain ${chain.id}`}, ${chain.residue_count} residues`).join("; ");
  const segments = chains.map((chain, index) => `<i style="flex-grow:${Math.max(1, Number(chain.residue_count) || 1)};background:${cartoonColors[index % cartoonColors.length]}"></i>`).join("");
  return `<span class="chain-strip" role="img" aria-label="${esc(description)}" title="${esc(description)}">${segments}</span>`;
}

function historyRow(analysis) {
  const busy = isActiveAnalysis(analysis);
  const status = busy ? (analysis.status === "queued" ? "queued" : "running") : analysis.status;
  const word = statusWords[status];
  const stamp = new Date(analysis.created_at).toLocaleString();
  const whenClass = word ? (busy ? "running" : "failed") : "";
  const when = word || relativeTime(analysis.created_at);
  const whenTitle = analysis.error ? `${word}: ${analysis.error}` : `${word ? `${word} · ` : ""}${stamp}`;
  const remove = busy
    ? '<span aria-hidden="true"></span>'
    : `<button class="history-delete" type="button" data-delete="${esc(analysis.id)}" title="Delete from history" aria-label="Delete ${esc(analysis.source_name)} from history">×</button>`;
  return `<div class="history-row ${analysis.id === state.analysisId ? "active" : ""}" data-id="${esc(analysis.id)}" tabindex="0"
      aria-label="Open ${esc(analysis.source_name)}, ${esc(word || "complete")}, ${esc(stamp)}">
    <span class="history-name" title="${esc(analysis.source_name)}">${esc(analysis.source_name)}</span>
    <span class="history-when ${whenClass}" title="${esc(busy ? `${whenTitle} · delete after the analysis finishes` : whenTitle)}">${esc(when)}</span>
    ${remove}
    ${chainStrip(analysis.chains)}
  </div>`;
}

async function refreshHistory() {
  try {
    const data = await getJSON("/api/analyses");
    const host = $("history");
    const count = data.analyses.length;
    $("history-count").textContent = count ? `(${count})` : "";
    $("clear-history").disabled = !count;
    if (!count) {
      host.innerHTML = '<p class="history-empty">Structures you open appear here.</p>';
      return;
    }
    host.innerHTML = data.analyses.map(historyRow).join("");
    host.querySelectorAll(".history-row").forEach((row) => {
      row.addEventListener("click", () => openAnalysis(row.dataset.id));
      row.addEventListener("keydown", (event) => {
        if (event.target !== row) return;
        if (event.key === "Enter" || event.key === " ") {
          event.preventDefault();
          openAnalysis(row.dataset.id);
        }
      });
    });
    host.querySelectorAll(".history-delete").forEach((button) => {
      button.addEventListener("click", (event) => {
        event.stopPropagation();
        deleteAnalysis(button.dataset.delete);
      });
    });
  } catch (error) {
    showError(error);
  }
}

function closeWorkspace() {
  resetState();
  state.analysisId = null;
  $("workspace").hidden = true;
  $("empty-workspace").hidden = false;
}

async function deleteAnalysis(id) {
  const row = $("history").querySelector(`.history-row[data-id="${CSS.escape(id)}"]`);
  const name = row?.querySelector(".history-name")?.textContent || "this analysis";
  if (!window.confirm(`Delete "${name}" and its saved results? This can't be undone.`)) return;
  clearError();
  try {
    const response = await fetch(`/api/analyses/${encodeURIComponent(id)}`, {method: "DELETE"});
    if (!response.ok && response.status !== 404) {
      const body = await response.json().catch(() => ({}));
      throw new Error(body.detail || "Delete failed");
    }
    if (state.analysisId === id) closeWorkspace();
  } catch (error) {
    showError(error);
  }
  refreshHistory();
}

async function clearHistory() {
  if (!window.confirm("Delete all recent analyses and their saved results? This can't be undone.")) return;
  clearError();
  try {
    const body = await getJSON("/api/analyses", {method: "DELETE"});
    if (body.deleted.includes(state.analysisId)) closeWorkspace();
    if (body.skipped.length) {
      showError(new Error(`${body.skipped.length} running analysis${body.skipped.length === 1 ? " was" : "es were"} kept. Delete ${body.skipped.length === 1 ? "it" : "them"} after the analysis finishes.`));
    }
  } catch (error) {
    showError(error);
  }
  refreshHistory();
}

// Above 1100px the menu sits beside the workspace and collapses to a rail.
// At or below it, the menu is an off-canvas drawer over the workspace.
const overlayQuery = window.matchMedia("(max-width: 1100px)");
let drawerReturnFocus = null;

function setDrawerCollapsed(collapsed) {
  const button = $("drawer-toggle");
  button.setAttribute("aria-expanded", String(!collapsed));
  button.textContent = collapsed ? "›" : "‹";
  button.title = collapsed ? "Expand menu" : "Collapse menu";
  button.setAttribute("aria-label", button.title);
  $("app-layout").classList.toggle("drawer-collapsed", collapsed);
  try { localStorage.setItem("drawerCollapsed", collapsed ? "1" : "0"); } catch (_) {}
}

function syncDrawerMode() {
  const button = $("drawer-toggle");
  const drawer = $("source-drawer");
  if (overlayQuery.matches) {
    const open = document.body.classList.contains("drawer-open-state");
    button.textContent = "‹";
    button.title = "Close menu";
    button.setAttribute("aria-label", "Close menu");
    button.setAttribute("aria-expanded", String(open));
    drawer.setAttribute("role", "dialog");
    drawer.toggleAttribute("aria-modal", open);
  } else {
    closeDrawer(false);
    drawer.removeAttribute("role");
    drawer.removeAttribute("aria-modal");
    setDrawerCollapsed($("app-layout").classList.contains("drawer-collapsed"));
  }
}

function openDrawer() {
  if (!overlayQuery.matches) {
    setDrawerCollapsed(false);
    return;
  }
  drawerReturnFocus = document.activeElement;
  document.body.classList.add("drawer-open-state");
  $("drawer-backdrop").hidden = false;
  $("main-column").inert = true;
  syncDrawerMode();
  $("pdb-id").focus();
}

function closeDrawer(restoreFocus = true) {
  if (!document.body.classList.contains("drawer-open-state")) return;
  const hadFocus = $("source-drawer").contains(document.activeElement);
  document.body.classList.remove("drawer-open-state");
  $("drawer-backdrop").hidden = true;
  $("main-column").inert = false;
  if (overlayQuery.matches) syncDrawerMode();
  if (restoreFocus && hadFocus && drawerReturnFocus?.isConnected) drawerReturnFocus.focus();
  drawerReturnFocus = null;
}

function setChosenFile(file) {
  $("file-label").innerHTML = file ? esc(file.name) : "Drop a .pdb or .cif file, or <u>browse</u>";
  $("upload-form").classList.toggle("has-file", Boolean(file));
  $("upload-submit").hidden = !file;
}

function bindDropZone() {
  const zone = $("upload-form");
  const input = $("structure-file");
  input.addEventListener("change", () => setChosenFile(input.files[0]));
  zone.addEventListener("dragover", (event) => {
    event.preventDefault();
    zone.classList.add("dragover");
  });
  zone.addEventListener("dragleave", (event) => {
    if (!zone.contains(event.relatedTarget)) zone.classList.remove("dragover");
  });
  zone.addEventListener("drop", (event) => {
    event.preventDefault();
    zone.classList.remove("dragover");
    const file = event.dataTransfer?.files?.[0];
    if (!file) return;
    const transfer = new DataTransfer();
    transfer.items.add(file);
    input.files = transfer.files;
    setChosenFile(file);
  });
}

async function submitUpload(event) {
  event.preventDefault();
  clearError();
  const input = $("structure-file");
  if (!input.files.length) return;
  const data = new FormData();
  data.append("file", input.files[0]);
  try {
    const response = await fetch("/api/analyses/upload", {method: "POST", body: data});
    const body = await response.json();
    if (!response.ok) throw new Error(body.detail || "Upload failed");
    await openAnalysis(body.analysis_id);
    refreshHistory();
  } catch (error) {
    showError(error);
  }
}

async function submitPdb(event) {
  event.preventDefault();
  clearError();
  try {
    const body = await getJSON("/api/analyses/pdb", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({
        pdb_id: $("pdb-id").value,
        assembly_id: $("assembly-id").value || null,
      }),
    });
    await openAnalysis(body.analysis_id);
    refreshHistory();
  } catch (error) {
    showError(error);
  }
}

async function inspectPdb() {
  const id = $("pdb-id").value.trim();
  const note = $("assembly-note");
  if (!/^[A-Za-z0-9]{4}$/.test(id)) {
    note.hidden = true;
    return;
  }
  try {
    const data = await getJSON(`/api/rcsb/${id}/metadata`);
    const ids = data.rcsb_entry_container_identifiers?.assembly_ids || [];
    note.textContent = ids.length
      ? `Assembl${ids.length === 1 ? "y" : "ies"} ${ids.join(", ")} available. Leave blank for the deposited coordinates.`
      : "No assemblies listed. The deposited coordinates will be used.";
  } catch (_) {
    note.textContent = "Couldn't check assemblies. Leave blank for the deposited coordinates.";
  }
  note.hidden = false;
}

function resetState() {
  state.analysisGeneration += 1;
  state.refreshing = false;
  if (state.timer) clearInterval(state.timer);
  state.timer = null;
  stopPocketPolling();
  cleanupViewer();
  state.detail = null;
  state.result = null;
  state.pairId = null;
  state.surface = null;
  state.structureText = null;
  state.structureFormat = null;
  state.pocketResult = null;
  state.pocketHighlighted = false;
  state.pocketHighlightCancelled = false;
  state.selectedResidues = [];
  state.curvatureLoaded = false;
  state.viewerLoaded = false;
  state.visibleChains = {};
  state.contact = {cut: 5, type: "all", vals: true, sasa: true};
  state.surfaceScale = "6";
  state.surfaceOpacity = 0.55;
  state.surfacePadding = 2;
  state.surfaceSeparation = 0;
  $("contact-cutoff").value = "5";
  $("contact-type").value = "all";
  $("contact-values").checked = true;
  $("contact-sasa").checked = true;
  $("surface-resolution").value = "6";
  $("surface-opacity").value = "0.55";
  $("surface-padding").value = "2";
  $("surface-separation").value = "0";
  $("protein-stats").textContent = "";
  $("surface-pair-stats").textContent = "";
  $("curvature-table").textContent = "";
  $("surface-note").textContent = "";
  $("surface-selection").textContent = "Select a residue in the surface, contact map, or interface table.";
  $("interfaces").innerHTML = '<p class="muted">Waiting for interface results…</p>';
  $("chain-toggles").textContent = "";
  hideContactTip();
  hideSurfaceTip();
  $("pocket-result").innerHTML = '<p class="muted">Choose a residue anchor and submit a pocket analysis.</p>';
  $("pocket-job-status").hidden = true;
  $("pocket-job-status").textContent = "";
  $("pair-workspace").hidden = true;
  $("pocket-panel").hidden = true;
  $("surface-status").textContent = "waiting";
  $("surface-status").className = "status-pill";
}

async function openAnalysis(id) {
  resetState();
  state.analysisId = id;
  $("workspace").hidden = false;
  $("empty-workspace").hidden = true;
  closeDrawer(false);
  $("history").querySelectorAll(".history-row").forEach((row) => {
    row.classList.toggle("active", row.dataset.id === id);
  });
  clearError();
  const generation = state.analysisGeneration;
  await refreshAnalysis();
  if (generation === state.analysisGeneration && !state.timer
    && !["complete", "failed", "cancelled", "interrupted"].includes(state.detail?.status)) {
    state.timer = setInterval(refreshAnalysis, 1200);
  }
}

async function refreshAnalysis() {
  if (!state.analysisId || state.refreshing) return;
  state.refreshing = true;
  const analysisId = state.analysisId;
  const generation = state.analysisGeneration;
  try {
    const detail = await getJSON(`/api/analyses/${analysisId}`);
    if (generation !== state.analysisGeneration) return;
    state.detail = detail;
    const analysis = state.detail;
    $("analysis-title").textContent = analysis.source_name;
    $("analysis-meta").textContent = `${analysis.source_kind} · ${analysis.source_format} · ${analysis.id}`;
    $("export-link").href = `/api/analyses/${state.analysisId}/export`;
    setStatus(analysis.status);
    if (analysis.error) {
      $("status-message").hidden = false;
      $("status-message").textContent = analysis.error;
    }
    if (analysis.result_path && !state.result) {
      const result = await getJSON(`/api/analyses/${analysisId}/result`);
      if (generation !== state.analysisGeneration) return;
      state.result = result;
      renderStructureMetadata();
      renderInterfaces();
      await loadStructureViewer();
      if (generation !== state.analysisGeneration) return;
    }
    renderJobs(analysis.jobs || []);
    if (["complete", "failed", "cancelled", "interrupted"].includes(analysis.status)) {
      if (state.timer) {
        clearInterval(state.timer);
        refreshHistory();
      }
      state.timer = null;
    }
  } catch (error) {
    if (generation === state.analysisGeneration) showError(error);
  } finally {
    if (generation === state.analysisGeneration) state.refreshing = false;
  }
}

function renderStructureMetadata() {
  const chains = state.result?.metadata?.chains || [];
  $("chain-count").textContent = `${chains.length} protein chain${chains.length === 1 ? "" : "s"}`;
  buildChainToggles(chains);
  renderSurface();
}

function renderJobs(jobs) {
  const core = [...jobs].reverse().find((job) => job.kind === "core");
  const curvature = [...jobs].reverse().find((job) => job.kind === "curvature");
  if (core && core.status === "running") {
    $("curvature").innerHTML = `<strong>Interface analysis running</strong><p class="muted">${esc(core.stage)} · ${Math.round(core.progress * 100)}% · calculating contacts and SASA…</p>`;
    setSurfaceStatus("waiting", "Surface queued");
  }
  if (core && core.status === "failed") {
    $("curvature").innerHTML = `<strong>Interface analysis failed</strong><p class="muted">${esc(core.error || "The interface job did not complete.")}</p>`;
    setSurfaceStatus("failed", "Surface unavailable");
  }
  if (curvature) {
    $("curvature").innerHTML = `<strong>Surface curvature</strong><p class="muted">${esc(curvature.stage)} · ${Math.round(curvature.progress * 100)}% · ${esc(curvature.status)}</p>`;
    if (curvature.status === "running" || curvature.status === "queued") {
      setSurfaceStatus("waiting", `Surface ${curvature.status}`);
    }
    if (curvature.result_path && !state.curvatureLoaded) loadCurvature();
    if (curvature.status === "failed") setSurfaceStatus("failed", "Surface unavailable");
  }
}

function setSurfaceStatus(status, text) {
  $("surface-status").textContent = text;
  $("surface-status").className = `status-pill ${status}`;
}

async function loadCurvature() {
  state.curvatureLoaded = true;
  const analysisId = state.analysisId;
  const generation = state.analysisGeneration;
  try {
    const surface = await getJSON(`/api/analyses/${analysisId}/curvature`);
    if (generation !== state.analysisGeneration) return;
    state.surface = surface;
    if (state.surface.status !== "complete" || !state.surface.report) {
      setSurfaceStatus("failed", state.surface.reason || "Surface unavailable");
      $("curvature").innerHTML = `<strong>Surface curvature unavailable</strong><p class="muted">${esc(state.surface.reason || "The optional surface job did not produce a mesh.")}</p>`;
      return;
    }
    setSurfaceStatus("complete", "Surface ready");
    renderSurfaceControls();
    renderContactMap();
    renderSurface();
  } catch (error) {
    if (generation !== state.analysisGeneration) return;
    state.curvatureLoaded = false;
    setSurfaceStatus("failed", "Surface unavailable");
    showError(error);
  }
}

function cleanupViewer() {
  clearPocketHighlight();
  for (const {component, representation} of state.selectionReps) {
    try { component.removeRepresentation(representation); } catch (_) {}
  }
  state.selectionReps = [];
  if (state.stage) {
    try { state.stage.removeAllComponents(); } catch (_) {}
    try { state.stage.dispose(); } catch (_) {}
  }
  state.stage = null;
  state.structureComponents = {};
  state.cartoonReps = {};
  state.meshComponents = {};
  state.meshReps = {};
  state.meshComponentChains = new Map();
  state.meshGeometry = {};
  state.meshSourceData = new WeakMap();
  state.meshSource = null;
  state.meshSelectionKey = null;
  state.meshColorKey = null;
}

function nglCompatibleMmcif(text) {
  const lines = String(text).split(/\r?\n/);
  const output = [];
  for (let i = 0; i < lines.length; i += 1) {
    if (lines[i].trim() !== "loop_") {
      output.push(lines[i]);
      continue;
    }
    let j = i + 1;
    const columns = [];
    while (j < lines.length) {
      const value = lines[j].trim();
      if (!value || value.startsWith("#")) { j += 1; continue; }
      if (!value.startsWith("_")) break;
      columns.push(value.split(/\s+/)[0]);
      j += 1;
    }
    if (columns.length && columns.every((column) => column.startsWith("_chem_comp."))) {
      while (j < lines.length && lines[j].trim() !== "#" && lines[j].trim() !== "stop_") j += 1;
      if (j < lines.length && lines[j].trim() === "#") j += 1;
      i = j - 1;
      continue;
    }
    output.push(lines[i]);
  }
  return output.join("\n");
}

async function loadStructureViewer() {
  if (state.viewerLoaded || !window.NGL || !state.result) return;
  state.viewerLoaded = true;
  const analysisId = state.analysisId;
  const generation = state.analysisGeneration;
  try {
    const response = await fetch(`/api/analyses/${analysisId}/structure`);
    if (!response.ok) throw new Error(response.statusText || "Coordinates could not be retrieved");
    const raw = await response.text();
    if (generation !== state.analysisGeneration) return;
    const format = state.result.source?.format || state.detail.source_format;
    state.structureFormat = format;
    state.structureText = format === "mmcif" ? nglCompatibleMmcif(raw) : raw;
    const stage = new NGL.Stage($("surface-view"), {backgroundColor: "#eef2f1", sampleLevel: 0});
    state.stage = stage;
    const extension = format === "mmcif" ? "cif" : "pdb";
    const chains = state.result.metadata?.chains || [];
    const blob = new Blob([state.structureText], {type: "text/plain"});
    await Promise.all(chains.map(async (chain, index) => {
      const component = await stage.loadFile(blob, {ext: extension, firstModelOnly: true, name: `${state.detail.source_name} chain ${chain.id}`});
      if (generation !== state.analysisGeneration) {
        stage.removeComponent(component);
        return;
      }
      const selection = `:${chain.id}`;
      state.structureComponents[chain.id] = component;
      state.cartoonReps[chain.id] = component.addRepresentation("cartoon", {
        sele: selection,
        color: cartoonColors[index % cartoonColors.length],
        opacity: 1,
        quality: "medium",
      });
    }));
    if (generation !== state.analysisGeneration) return;
    stage.autoView();
    bindStageSignals();
    renderSurface();
  } catch (error) {
    if (generation !== state.analysisGeneration) return;
    state.viewerLoaded = false;
    setSurfaceStatus("failed", "Viewer unavailable");
    $("surface-note").textContent = `The structure viewer could not load these coordinates: ${error.message || error}`;
  }
}

function renderInterfaces() {
  const pairs = state.result?.pairs || [];
  $("interfaces").innerHTML = pairs.length
    ? pairs.map((pair) => `<article class="interface-card ${state.pairId === pair.id ? "active" : ""}" data-pair="${esc(pair.id)}" tabindex="0">
        <h3>${esc(chainLabel(pair.chain_a))} × ${esc(chainLabel(pair.chain_b))}</h3>
        <div class="numbers"><span class="chip">BSA ${fmt(pair.buried_surface_area_A2, 0)} Å²</span><span class="chip">${pair.n_contact_residue_pairs} contact pairs</span><span class="chip">${pair.n_hbond_like} polar</span><span class="chip">${pair.n_salt_bridges} salt</span></div>
      </article>`).join("")
    : '<p class="muted">Only one protein chain was found; no pairwise interface exists.</p>';
  $("interfaces").querySelectorAll(".interface-card").forEach((card) => {
    card.addEventListener("click", () => selectPair(card.dataset.pair));
    card.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        selectPair(card.dataset.pair);
      }
    });
  });
  if (!pairs.length) {
    $("pair-workspace").hidden = true;
    $("pocket-panel").hidden = true;
    return;
  }
  if (!state.pairId || !pairs.some((pair) => pair.id === state.pairId)) selectPair(pairs[0].id);
}

function selectPair(id) {
  const pair = (state.result?.pairs || []).find((item) => item.id === id);
  if (!pair) return;
  state.pairId = id;
  state.activeTab = "residue";
  state.selectedResidues = [];
  state.pocketResult = null;
  state.pocketHighlightCancelled = false;
  stopPocketPolling();
  clearPocketHighlight();
  $("pocket-result").innerHTML = '<p class="muted">Choose a residue anchor and submit a pocket analysis.</p>';
  $("pocket-job-status").hidden = true;
  $("pocket-job-status").textContent = "";
  $("pair-workspace").hidden = false;
  $("pocket-panel").hidden = false;
  $("pair-title").textContent = `${chainLabel(pair.chain_a)} × ${chainLabel(pair.chain_b)}`;
  $("pair-summary").textContent = `${pair.chain_a} · ${pair.chain_b}`;
  renderInterfacesActiveOnly();
  renderPair(pair);
  renderTabs();
  renderContactMap();
  renderSurface();
}

function renderInterfacesActiveOnly() {
  $("interfaces").querySelectorAll(".interface-card").forEach((card) => {
    card.classList.toggle("active", card.dataset.pair === state.pairId);
  });
}

function renderPair(pair) {
  $("pair-stats").innerHTML = [
    ["BSA", `${fmt(pair.buried_surface_area_A2, 0)} Å²`],
    ["Interface residues", pair.n_interface_residues],
    ["H-bond-like", pair.n_hbond_like],
    ["Salt bridges", pair.n_salt_bridges],
    ["Like charge", pair.n_like_charges || 0],
  ].map(([label, value]) => `<div class="stat"><b>${esc(value)}</b><span>${esc(label)}</span></div>`).join("");

  $("residue-table").innerHTML = pair.residues.length
    ? pair.residues.map((residue) => `<tr class="focus-row" data-key="${esc(residueKey(residue))}" tabindex="0">
        <td><strong>${esc(residueDisplay(residue))}</strong></td><td>${fmt(residue.dSASA_A2)}</td><td>${residue.n_contacts}</td><td>${fmt(residue.min_dist_A, 2)} Å</td><td>${esc((residue.partners || []).join(", "))}</td>
      </tr>`).join("")
    : '<tr><td colspan="5" class="muted">No reported interface residues.</td></tr>';
  $("contact-table").innerHTML = pair.contacts.length
    ? pair.contacts.slice(0, 240).map((contact) => `<tr class="focus-row" data-keys="${esc([residueKey(contact.residue_a), residueKey(contact.residue_b)].join("|"))}" tabindex="0">
        <td>${esc(residueDisplay(contact.residue_a))}</td><td>${esc(residueDisplay(contact.residue_b))}</td><td>${fmt(contact.min_dist_A, 2)} Å</td><td>${contact.n_contacts}</td><td>${contact.hbond_like.length ? contact.hbond_like.length : "—"}</td><td>${contact.salt_bridges.length ? contact.salt_bridges.length : "—"}</td><td>${contact.like_charges?.length ? contact.like_charges.length : "—"}</td>
      </tr>`).join("")
    : '<tr><td colspan="7" class="muted">No contacts were found within 6 Å.</td></tr>';
  $("residue-table").querySelectorAll(".focus-row").forEach((row) => bindSelectionRow(row, [row.dataset.key]));
  $("contact-table").querySelectorAll(".focus-row").forEach((row) => bindSelectionRow(row, row.dataset.keys.split("|")));

  const pairChains = [pair.chain_a, pair.chain_b];
  const chainOptions = pairChains.map((chain) => `<option value="${esc(chain)}">${esc(chainLabel(chain))}</option>`).join("");
  $("pocket-target").innerHTML = chainOptions;
  $("pocket-partner").innerHTML = chainOptions;
  updateAnchors();
}

function bindSelectionRow(row, keys) {
  const choose = () => selectResidues(keys.filter(Boolean));
  row.addEventListener("click", choose);
  row.addEventListener("keydown", (event) => {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      choose();
    }
  });
}

function renderTabs() {
  const residueActive = state.activeTab === "residue";
  $("residue-tab").classList.toggle("active", residueActive);
  $("residue-tab").setAttribute("aria-selected", String(residueActive));
  $("contact-tab").classList.toggle("active", !residueActive);
  $("contact-tab").setAttribute("aria-selected", String(!residueActive));
  $("residue-panel").hidden = !residueActive;
  $("contact-panel").hidden = residueActive;
}

function chainResidues(chain) {
  return (state.result?.metadata?.chains || []).find((item) => item.id === chain)?.residues || [];
}

function surfaceResidueKey(chain, label) {
  const residue = chainResidues(chain).find((item) => item.label === label);
  return residue ? residueKey(residue) : `${chain}:${label.replace(/^[A-Z]/, "")}`;
}

function updateAnchors() {
  const target = $("pocket-target").value;
  $("pocket-anchor").innerHTML = chainResidues(target).map((residue) => `<option value="${esc(`${residue.number}:${residue.insertion_code || ""}`)}">${esc(residue.label)} · ${esc(residue.resname)}</option>`).join("");
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  if (pair) $("pocket-partner").value = target === pair.chain_a ? pair.chain_b : pair.chain_a;
}

function selectResidues(keys) {
  state.selectedResidues = [...new Set(keys)];
  renderSelectionState();
  renderSurface();
}

function renderSelectionState() {
  document.querySelectorAll(".focus-row[data-key], .focus-row[data-keys]").forEach((row) => {
    const keys = (row.dataset.keys || row.dataset.key || "").split("|");
    row.classList.toggle("selected", keys.some((key) => state.selectedResidues.includes(key)));
  });
  document.querySelectorAll("[data-residue-key]").forEach((element) => {
    element.classList.toggle("selected", state.selectedResidues.includes(element.dataset.residueKey));
  });
  document.querySelectorAll("[data-cell-key]").forEach((element) => {
    const keys = (element.dataset.cellKey || "").split("|");
    element.classList.toggle("selected", keys.some((key) => state.selectedResidues.includes(key)));
  });
}

function selectionParts(key) {
  const match = String(key).match(/^([^:]+):(\-?\d+)(.*)$/);
  if (!match) return null;
  return {chain: match[1], number: match[2], insertion: match[3] || ""};
}

function selectionString(key) {
  const part = selectionParts(key);
  return part ? `:${part.chain} and ${part.number}${part.insertion ? `^${part.insertion}` : ""}` : null;
}

function updateNglSelection() {
  for (const {component, representation} of state.selectionReps) {
    try { component.removeRepresentation(representation); } catch (_) {}
  }
  state.selectionReps = [];
  if (!state.stage) return;
  for (const key of state.selectedResidues) {
    const part = selectionParts(key);
    if (!part || !state.structureComponents[part.chain]) continue;
    try {
      const component = state.structureComponents[part.chain];
      const representation = component.addRepresentation("ball+stick", {
        sele: selectionString(key), color: "#d67439", scale: 0.85,
      });
      state.selectionReps.push({component, representation});
    } catch (_) {}
  }
}

function curvatureReport() {
  return state.surface?.report || null;
}

function curvatureScale() {
  return curvatureReport()?.scales?.[state.surfaceScale] || null;
}

function interfaceSurface(chain) {
  return curvatureScale()?.interfaces?.[state.pairId]?.[chain] || null;
}

function wholeSurface(chain) {
  return curvatureScale()?.proteins?.[chain] || null;
}

function fmtH(value) {
  return value == null ? "n/a" : `${value >= 0 ? "+" : ""}${Number(value).toFixed(3)} Å⁻¹`;
}

function fmtPercent(value) {
  return value == null ? "n/a" : `${(100 * Number(value)).toFixed(1)}%`;
}

function renderSurfaceControls() {
  const report = curvatureReport();
  if (!report) return;
  const scaleSelect = $("surface-resolution");
  scaleSelect.innerHTML = Object.keys(report.scales || {}).map((scale) => `<option value="${esc(scale)}">${esc(scale)} Å · ${scale === "6" ? "default" : scale === "4" ? "fine" : "broad"}</option>`).join("");
  scaleSelect.value = report.scales[state.surfaceScale] ? state.surfaceScale : Object.keys(report.scales)[0] || "6";
  state.surfaceScale = scaleSelect.value;
  buildChainToggles(state.result?.metadata?.chains || []);
}

function buildChainToggles(chains) {
  const host = $("chain-toggles");
  host.textContent = "";
  for (const chain of chains) {
    if (!(chain.id in state.visibleChains)) state.visibleChains[chain.id] = true;
    const label = document.createElement("label");
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = state.visibleChains[chain.id];
    checkbox.addEventListener("change", () => {
      state.visibleChains[chain.id] = checkbox.checked;
      renderSurface();
    });
    label.append(checkbox, ` ${chainLabel(chain.id)}`);
    host.appendChild(label);
  }
}

function interpolateColor(value) {
  let index = 0;
  while (index < curvatureColors.length - 2 && value > curvatureColors[index + 1][0]) index += 1;
  const [a, first] = curvatureColors[index];
  const [b, second] = curvatureColors[index + 1];
  const amount = (value - a) / (b - a);
  return first.map((item, component) => (item + amount * (second[component] - item)) / 255);
}

function surfaceMeshGeometry(mesh, selectedLabels) {
  let source = state.meshSourceData.get(mesh);
  if (!source) {
    const position = new Float32Array(mesh.position);
    const index = new Uint32Array(mesh.index);
    // Compute smooth normals on the original surface once. Splitting vertices
    // at a yellow patch boundary must not change the surface's lighting.
    const buffer = new NGL.MeshBuffer({position, index, color: new Float32Array(position.length)});
    const normal = new Float32Array(buffer.geometry.attributes.normal.array);
    buffer.dispose();
    source = {position, index, normal};
    state.meshSourceData.set(mesh, source);
  }
  const vertexCount = source.position.length / 3;
  const selected = mesh.residue.map((label) => selectedLabels.has(label));
  const selectedFaces = [];
  const copies = new Map();
  for (let face = 0; face < source.index.length; face += 3) {
    const vertices = source.index.subarray(face, face + 3);
    if (!vertices.some((vertex) => selected[vertex])) continue;
    selectedFaces.push(face);
    for (const vertex of vertices) {
      if (!copies.has(vertex)) copies.set(vertex, vertexCount + copies.size);
    }
  }
  const vertexMap = Uint32Array.from({length: vertexCount + copies.size}, (_, index) => index);
  if (!copies.size) return {...source, vertexMap, highlightedVertexStart: vertexCount};

  // Give every selected face its own yellow vertices. Neighboring unselected
  // faces retain their original vertices and curvature colors, without color
  // interpolation or a second transparent mesh overlapping the surface.
  const position = new Float32Array(vertexMap.length * 3);
  const normal = new Float32Array(vertexMap.length * 3);
  const index = new Uint32Array(source.index);
  position.set(source.position);
  normal.set(source.normal);
  for (const [original, copy] of copies) {
    position.set(source.position.subarray(original * 3, original * 3 + 3), copy * 3);
    normal.set(source.normal.subarray(original * 3, original * 3 + 3), copy * 3);
    vertexMap[copy] = original;
  }
  for (const face of selectedFaces) {
    for (let corner = 0; corner < 3; corner += 1) index[face + corner] = copies.get(source.index[face + corner]);
  }
  return {position, index, normal, vertexMap, highlightedVertexStart: vertexCount};
}

function meshColors(chain) {
  const mesh = state.surface?.meshes?.[chain];
  const report = curvatureReport();
  if (!mesh || !report) return new Float32Array();
  const values = mesh.h?.[state.surfaceScale] || [];
  const geometry = state.meshGeometry[chain];
  const colors = new Float32Array((geometry?.vertexMap.length || values.length) * 3);
  const limit = Number(report.color_limit_Ainv) || 1;
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  const partner = pair && pair.chain_a === chain ? pair.chain_b : pair && pair.chain_b === chain ? pair.chain_a : null;
  const distances = partner ? mesh.contact_distance_dA?.[partner] : null;
  const padding = Number(state.surfacePadding) * 10;
  const grey = [157 / 255, 164 / 255, 172 / 255];
  for (let index = 0; index < values.length; index += 1) {
    // The reference contact-map view colors only the two chains in the
    // selected interface. Other visible chains remain available as context,
    // but are rendered grey rather than being mistaken for interface shape.
    const inPadding = pair
      ? Boolean(distances) && distances[index] <= padding
      : true;
    const color = inPadding
      ? interpolateColor(Math.max(0, Math.min(1, (values[index] / limit + 1) / 2)))
      : grey;
    colors[index * 3] = color[0];
    colors[index * 3 + 1] = color[1];
    colors[index * 3 + 2] = color[2];
  }
  for (let index = geometry?.highlightedVertexStart ?? values.length; index < colors.length / 3; index += 1) {
    colors[index * 3] = 1;
    colors[index * 3 + 1] = 1;
    colors[index * 3 + 2] = 0;
  }
  return colors;
}

function installSurfaceMesh(chain) {
  const geometry = state.meshGeometry[chain];
  const shapeName = `surface-${encodeURIComponent(chain)}`;
  const shape = new NGL.Shape(shapeName);
  shape.addMesh(geometry.position, meshColors(chain), geometry.index, geometry.normal, `${chain} surface`);
  const component = state.stage.addComponentFromObject(shape);
  state.meshComponents[chain] = component;
  state.meshReps[chain] = component.addRepresentation("buffer", {
    opacity: state.surfaceOpacity,
    side: "double",
    metalness: 0,
    roughness: 1,
  });
  state.meshComponentChains.set(shapeName, chain);
}

function installSurfaceMeshes() {
  if (!state.stage || state.surface?.status !== "complete") return;
  const meshes = state.surface.meshes || {};
  const selectionKey = JSON.stringify([...state.selectedResidues].sort());
  if (state.meshSource === meshes && state.meshSelectionKey === selectionKey) return;
  for (const component of Object.values(state.meshComponents)) {
    try { state.stage.removeComponent(component); } catch (_) {}
  }
  state.meshComponents = {};
  state.meshReps = {};
  state.meshComponentChains = new Map();
  state.meshGeometry = {};
  state.meshSource = meshes;
  state.meshSelectionKey = selectionKey;
  state.meshColorKey = null;
  for (const [chain, mesh] of Object.entries(meshes)) {
    if (!mesh.faces?.length && !mesh.index?.length) continue;
    try {
      const labels = new Set(state.selectedResidues
        .filter((key) => selectionParts(key)?.chain === chain)
        .map(residueLabelFromKey));
      state.meshGeometry[chain] = surfaceMeshGeometry(mesh, labels);
      installSurfaceMesh(chain);
    } catch (error) {
      $("surface-note").textContent = `Surface mesh for chain ${chain} could not be displayed: ${error.message || error}`;
    }
  }
}

function updateMeshColors() {
  const key = JSON.stringify([state.pairId, state.surfaceScale, state.surfacePadding, state.meshSelectionKey]);
  if (state.meshColorKey === key) return;
  for (const [chain, component] of Object.entries(state.meshComponents)) {
    try {
      component.object.bufferList[0].setAttributes({color: meshColors(chain)});
    } catch (_) {
      // Some NGL builds do not expose bufferList; reinstalling keeps the
      // scientific result usable across the bundled viewer versions.
      if (!state.meshGeometry[chain]) continue;
      try { state.stage.removeComponent(component); } catch (_) {}
      delete state.meshComponents[chain];
      delete state.meshReps[chain];
      installSurfaceMesh(chain);
    }
  }
  state.meshColorKey = key;
}

function pairDirection() {
  return curvatureReport()?.pair_directions?.[state.pairId] || [0, 0, 0];
}

function chainShift(chain) {
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  if (!pair || !pairDirection().some((value) => value !== 0)) return [0, 0, 0];
  const side = chain === pair.chain_a ? -1 : chain === pair.chain_b ? 1 : 0;
  const amount = side * Number(state.surfaceSeparation) / 2;
  return pairDirection().map((value) => value * amount);
}

function applyPosition(component, chain) {
  if (!component?.setPosition) return;
  component.setPosition(chainShift(chain));
}

function renderSurfaceStats() {
  const report = curvatureReport();
  const scale = curvatureScale();
  const proteinHost = $("protein-stats");
  const pairHost = $("surface-pair-stats");
  proteinHost.textContent = "";
  pairHost.textContent = "";
  if (!report || !scale) return;
  for (const chain of Object.keys(scale.proteins || {})) {
    const summary = scale.proteins[chain];
    const card = document.createElement("div");
    card.className = "surface-stat-card";
    card.innerHTML = `<strong>${esc(chainLabel(chain))} · whole resolved surface</strong><span>Area ${fmt(summary.area_A2)} Å²</span><span>Mean H ${esc(fmtH(summary.mean_H_Ainv))}</span><span>Convex ${esc(fmtPercent(summary.convex_fraction))} · concave ${esc(fmtPercent(summary.concave_fraction))}</span>`;
    proteinHost.appendChild(card);
  }
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  const interfaceSummary = pair ? scale.interfaces?.[pair.id] : null;
  if (interfaceSummary) {
    for (const chain of [pair.chain_a, pair.chain_b]) {
      const summary = interfaceSummary[chain];
      if (!summary) continue;
      const card = document.createElement("div");
      card.className = "surface-stat-card interface-stat";
      card.innerHTML = `<strong>${esc(chainLabel(chain))} · interface face</strong><span>Patch area ${fmt(summary.area_A2)} Å²</span><span>Mean H ${esc(fmtH(summary.mean_H_Ainv))}</span><span>Convex ${esc(fmtPercent(summary.convex_fraction))} · concave ${esc(fmtPercent(summary.concave_fraction))}</span>`;
      pairHost.appendChild(card);
    }
    const sideA = interfaceSummary[pair.chain_a];
    const sideB = interfaceSummary[pair.chain_b];
    if (sideA && sideB) {
      const card = document.createElement("div");
      card.className = "surface-stat-card interface-stat";
      card.innerHTML = `<strong>Reference interface area</strong><span>${fmt(((sideA.reference_dSASA_A2 || 0) + (sideB.reference_dSASA_A2 || 0)) / 2, 0)} Å² from ΔSASA</span><span>Mesh patches: ${fmt(sideA.area_A2, 0)} / ${fmt(sideB.area_A2, 0)} Å²</span><span>Areas are computed independently.</span>`;
      pairHost.appendChild(card);
    }
  }
  $("surface-note").textContent = `Positive H is convex/outward. Colors are clipped at ±${fmt(report.color_limit_Ainv, 3)} Å⁻¹ for display; numeric values are not clipped. Contact padding changes only the colored display, while separation changes only viewing positions.`;
}

function renderCurvatureTable() {
  const tbody = $("curvature-table");
  tbody.textContent = "";
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  const scale = curvatureScale();
  if (!pair || !scale?.interfaces?.[pair.id]) return;
  const dsasa = new Map();
  for (const entry of [...(pair.contact_map?.rows || []), ...(pair.contact_map?.cols || [])]) dsasa.set(entry.key, entry.dSASA_A2);
  const records = [];
  for (const chain of [pair.chain_a, pair.chain_b]) {
    for (const item of scale.interfaces[pair.id][chain]?.residues || []) records.push({chain, item});
  }
  records.sort((a, b) => (b.item.area_A2 || 0) - (a.item.area_A2 || 0));
  for (const record of records) {
    const key = surfaceResidueKey(record.chain, record.item.label);
    const row = document.createElement("tr");
    row.className = "focus-row";
    row.tabIndex = 0;
    row.dataset.key = key;
    row.innerHTML = `<td>${esc(record.chain)}</td><td>${esc(record.item.label)}</td><td>${fmt(record.item.area_A2)} Å²</td><td>${esc(fmtH(record.item.mean_H_Ainv))}</td><td>${esc(fmtPercent(record.item.convex_fraction))}</td><td>${esc(fmtPercent(record.item.concave_fraction))}</td><td>${fmt(dsasa.get(key))}</td>`;
    bindSelectionRow(row, [key]);
    tbody.appendChild(row);
  }
}

function residueLabelFromKey(key) {
  const part = selectionParts(key);
  if (!part) return "";
  const chain = (state.result?.metadata?.chains || []).find((item) => item.id === part.chain);
  return chain?.residues?.find((residue) => String(residue.number) === part.number && String(residue.insertion_code || "") === part.insertion)?.label || "";
}

function updatePocketPosition() {
  if (!state.pocketComponent) return;
  const chain = state.pocketResult?.target_chain;
  applyPosition(state.pocketComponent, chain);
  state.pocketComponent.setVisibility(state.visibleChains[chain] !== false);
}

function renderSurface() {
  const report = curvatureReport();
  const opacityText = `${Math.round(Number(state.surfaceOpacity) * 100)}%`;
  const paddingText = `${Number(state.surfacePadding).toFixed(1)} Å`;
  const separationText = `${Number(state.surfaceSeparation).toFixed(1)} Å`;
  $("surface-opacity-value").value = opacityText;
  $("surface-opacity-value").textContent = opacityText;
  $("surface-padding-value").value = paddingText;
  $("surface-padding-value").textContent = paddingText;
  $("surface-separation-value").value = separationText;
  $("surface-separation-value").textContent = separationText;
  if (report) {
    const limit = Number(report.color_limit_Ainv) || 0;
    $("surface-color-min").textContent = `−${limit.toFixed(3)}`;
    $("surface-color-max").textContent = `+${limit.toFixed(3)}`;
    $("surface-color-label").textContent = `H (Å⁻¹) · ${state.surfaceScale} Å neighborhoods`;
  }
  renderSurfaceStats();
  renderCurvatureTable();
  selectedSummary();
  if (!state.stage) return;
  // Residue selections should work as soon as the coordinate viewer is ready;
  // they do not depend on the optional curvature job.
  if (!state.surface || state.surface.status !== "complete") {
    const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
    for (const chain of Object.keys(state.structureComponents)) {
      applyPosition(state.structureComponents[chain], chain);
      state.structureComponents[chain].setVisibility(state.visibleChains[chain] !== false);
      const isPair = pair && (chain === pair.chain_a || chain === pair.chain_b);
      state.cartoonReps[chain]?.setParameters({opacity: isPair ? 1 : 0.3});
    }
    updateNglSelection();
    if (state.pocketResult && !state.pocketHighlighted && !state.pocketHighlightCancelled) showPocketHighlight();
    updatePocketPosition();
    updatePocketButtons();
    state.stage.viewer?.requestRender();
    return;
  }
  installSurfaceMeshes();
  updateMeshColors();
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  for (const chain of Object.keys(state.structureComponents)) {
    applyPosition(state.structureComponents[chain], chain);
    applyPosition(state.meshComponents[chain], chain);
    state.structureComponents[chain].setVisibility(state.visibleChains[chain] !== false);
    if (state.meshComponents[chain]) state.meshComponents[chain].setVisibility(state.visibleChains[chain] !== false);
    const isPair = pair && (chain === pair.chain_a || chain === pair.chain_b);
    state.cartoonReps[chain]?.setParameters({opacity: isPair ? 1 : 0.3});
    state.meshReps[chain]?.setParameters({opacity: isPair ? state.surfaceOpacity : state.surfaceOpacity * 0.28});
  }
  updateNglSelection();
  if (state.pocketResult && !state.pocketHighlighted && !state.pocketHighlightCancelled) showPocketHighlight();
  updatePocketPosition();
  updatePocketButtons();
  state.stage.viewer?.requestRender();
}

function selectedSummary() {
  const detail = $("surface-selection");
  const scale = curvatureScale();
  if (!state.selectedResidues.length) {
    detail.textContent = "Select a residue in the surface, contact map, or interface table.";
    renderSelectionState();
    return;
  }
  if (!scale) {
    detail.textContent = `Selected: ${state.selectedResidues.map((key) => {
      const part = selectionParts(key);
      return `${part?.chain || "?"}:${residueLabelFromKey(key) || key}`;
    }).join(" × ")}`;
    renderSelectionState();
    return;
  }
  detail.textContent = state.selectedResidues.map((key) => {
    const part = selectionParts(key);
    const label = residueLabelFromKey(key);
    const whole = scale.proteins?.[part?.chain]?.residues?.find((item) => item.label === label);
    const patch = interfaceSurface(part?.chain)?.residues?.find((item) => item.label === label);
    const name = chainLabel(part?.chain);
    if (!whole) return `${name} ${label}: no surface vertices assigned`;
    return `${name} ${label}: whole H ${fmtH(whole.mean_H_Ainv)}, ${fmt(whole.area_A2)} Å²${patch ? `; interface H ${fmtH(patch.mean_H_Ainv)}, ${fmt(patch.area_A2)} Å²` : "; no selected interface patch"}`;
  }).join("  |  ");
  renderSelectionState();
}

function contactMapMatches(cell) {
  if (cell.min_dist_A > state.contact.cut) return false;
  switch (state.contact.type) {
    case "tight": return cell.n_tight > 0;
    case "hb": return cell.hbond_like.length > 0;
    case "sb": return cell.salt_bridges.length > 0;
    case "like": return cell.like_charges?.length > 0 && cell.salt_bridges.length === 0;
    default: return true;
  }
}

function contactMapSlice(data) {
  const cells = (data?.cells || []).filter(contactMapMatches);
  const rows = [...new Set(cells.map((cell) => cell.row))].sort((a, b) => a - b);
  const cols = [...new Set(cells.map((cell) => cell.col))].sort((a, b) => a - b);
  return {cells, rows, cols};
}

function distanceBin(distance) {
  for (let index = 0; index < distanceBins.length; index += 1) if (distance <= distanceBins[index]) return index;
  return distanceBins.length - 1;
}

function chargeClass(residue) {
  const charge = residue.charge ?? formalCharge[residue.resname] ?? 0;
  if (charge > 0) return "pos";
  if (charge < 0) return "neg";
  const hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET"];
  const aromatic = ["PHE", "TRP", "TYR"];
  const special = ["GLY", "PRO", "CYS"];
  if (aromatic.includes(residue.resname)) return "arom";
  if (hydrophobic.includes(residue.resname)) return "hyd";
  if (special.includes(residue.resname)) return "spec";
  return "polar";
}

function laneGlyph(residue) {
  const cls = chargeClass(residue);
  return cls === "pos" ? "+" : cls === "neg" ? "−" : residueClassGlyph[cls];
}

function chargeLabel(residue) {
  const charge = residue.charge ?? formalCharge[residue.resname] ?? 0;
  return charge > 0 ? "+1" : charge < 0 ? "−1" : "0";
}

function renderDistanceLegend() {
  const host = $("distance-legend");
  host.innerHTML = distanceLabels.map((label, index) => `<span><i class="map-bin b${index}"></i>${label}</span>`).join("");
}

function addGridElement(grid, className, text, column, row, attrs = {}) {
  const element = document.createElement("div");
  element.className = className;
  if (text != null) element.textContent = text;
  element.style.gridColumn = String(column);
  element.style.gridRow = String(row);
  Object.entries(attrs).forEach(([key, value]) => element.setAttribute(key, value));
  grid.appendChild(element);
  return element;
}

function contactAxisLayout(axis, indices) {
  const tracks = [];
  const gaps = [];
  let track = 4;
  const positions = indices.map((index, position) => {
    if (position > 0 && axis[index].number - axis[indices[position - 1]].number > 1) {
      tracks.push("8px");
      gaps.push(track++);
    }
    tracks.push("30px");
    return track++;
  });
  return {tracks, positions, gaps};
}

function decorateMapLabel(label, residue, column) {
  const summary = interfaceSurface(residue.chain_id)?.residues?.find((item) => item.label === residue.label)
    || wholeSurface(residue.chain_id)?.residues?.find((item) => item.label === residue.label);
  label.title = `${residueDisplay(residue)} · ${residue.resname}`;
  if (summary?.mean_H_Ainv != null) {
    const limit = Number(curvatureReport()?.color_limit_Ainv) || 1;
    const color = interpolateColor(Math.max(0, Math.min(1, (summary.mean_H_Ainv / limit + 1) / 2)));
    label.style.boxShadow = `inset ${column ? "0 -4px" : "4px 0"} 0 rgb(${color.map((value) => Math.round(value * 255)).join(",")})`;
    label.title += ` · mean H ${fmtH(summary.mean_H_Ainv)}`;
  }
}

function addSasaLane(grid, residue, column, row) {
  const sasa = addGridElement(grid, `map-sasa${residue.dSASA_A2 == null ? " none" : ""}`, "", column, row, {"aria-hidden": "true"});
  sasa.title = residue.dSASA_A2 == null ? "ΔSASA not computed at the reported interface cutoff" : `ΔSASA ${fmt(residue.dSASA_A2)} Å²`;
  if (residue.dSASA_A2 != null) sasa.style.opacity = Math.max(0.08, Math.min(1, residue.dSASA_A2 / 120));
}

function showContactTip(element, cell, data) {
  if (!cell) return;
  const row = data.rows[cell.row];
  const col = data.cols[cell.col];
  const tip = $("contact-tip");
  tip.innerHTML = `<strong>${esc(residueDisplay(row))} × ${esc(residueDisplay(col))}</strong><dl><dt>Minimum distance</dt><dd>${fmt(cell.min_dist_A, 2)} Å</dd><dt>Atom contacts</dt><dd>${cell.n_contacts}</dd><dt>Tight contacts</dt><dd>${cell.n_tight}</dd>${cell.hbond_like.length ? `<dt>H-bonds</dt><dd>${cell.hbond_like.length}</dd>` : ""}${cell.salt_bridges.length ? `<dt>Salt bridges</dt><dd>${cell.salt_bridges.length}</dd>` : ""}${cell.like_charges?.length ? `<dt>Like charge</dt><dd>${cell.like_charges.length}</dd>` : ""}</dl>`;
  tip.classList.add("on");
  tip.setAttribute("aria-hidden", "false");
  const bounds = element.getBoundingClientRect();
  const box = tip.getBoundingClientRect();
  const left = bounds.right + 12 + box.width > window.innerWidth ? bounds.left - box.width - 12 : bounds.right + 12;
  const top = Math.min(bounds.top, window.innerHeight - box.height - 8);
  tip.style.left = `${Math.max(8, left)}px`;
  tip.style.top = `${Math.max(8, top)}px`;
}

function hideContactTip() {
  $("contact-tip").classList.remove("on");
  $("contact-tip").setAttribute("aria-hidden", "true");
}

function renderContactMap() {
  const pair = state.result?.pairs?.find((item) => item.id === state.pairId);
  const data = pair?.contact_map;
  const grid = $("contact-grid");
  grid.textContent = "";
  if (!data) {
    $("contact-map-count").textContent = "No contact map data";
    $("contact-col-axis").textContent = "";
    $("contact-row-axis").textContent = "";
    $("contact-map-table").textContent = "";
    $("contact-map-note").textContent = "";
    return;
  }
  const slice = contactMapSlice(data);
  const rowCount = slice.rows.length;
  const colCount = slice.cols.length;
  $("contact-map-count").textContent = `${slice.cells.length} pairs · ${rowCount} × ${colCount} residues`;
  const columns = contactAxisLayout(data.cols, slice.cols);
  const rows = contactAxisLayout(data.rows, slice.rows);
  grid.style.gridTemplateColumns = `minmax(44px, max-content) 28px ${state.contact.sasa ? "12px" : "0px"} ${columns.tracks.join(" ") || "30px"}`;
  grid.style.gridTemplateRows = `82px 28px ${state.contact.sasa ? "12px" : "0px"} ${rows.tracks.join(" ") || "30px"}`;
  grid.setAttribute("aria-label", `${chainLabel(data.row_chain)} by ${chainLabel(data.col_chain)} residue contact map`);
  $("contact-col-axis").textContent = `${chainLabel(data.col_chain)} · columns →`;
  $("contact-row-axis").textContent = `${chainLabel(data.row_chain)} · rows ↓`;
  for (const row of rows.gaps) {
    const line = addGridElement(grid, "map-gap-line horizontal", "", 4, row, {"aria-hidden": "true"});
    line.style.gridColumn = `4 / span ${Math.max(1, columns.tracks.length)}`;
  }
  for (const column of columns.gaps) {
    const line = addGridElement(grid, "map-gap-line vertical", "", column, 4, {"aria-hidden": "true"});
    line.style.gridRow = `4 / span ${Math.max(1, rows.tracks.length)}`;
  }
  const cellByPosition = new Map(slice.cells.map((cell) => [`${cell.row}:${cell.col}`, cell]));
  slice.cols.forEach((columnIndex, index) => {
    const residue = data.cols[columnIndex];
    const column = columns.positions[index];
    const label = addGridElement(grid, "map-col-label", residue.label, column, 1, {"data-residue-key": residueKey(residue), tabindex: "0"});
    decorateMapLabel(label, residue, true);
    bindMapAxisLabel(label, residue);
    const lane = addGridElement(grid, `map-lane q-${residue.charge > 0 ? "pos" : residue.charge < 0 ? "neg" : "0"}`, laneGlyph(residue), column, 2, {"aria-hidden": "true"});
    lane.title = `${residueDisplay(residue)} · ${residueClassName[chargeClass(residue)]}`;
    if (state.contact.sasa) addSasaLane(grid, residue, column, 3);
  });
  slice.rows.forEach((rowIndex, index) => {
    const residue = data.rows[rowIndex];
    const row = rows.positions[index];
    const label = addGridElement(grid, "map-row-label", residue.label, 1, row, {"data-residue-key": residueKey(residue), tabindex: "0"});
    decorateMapLabel(label, residue, false);
    bindMapAxisLabel(label, residue);
    const lane = addGridElement(grid, `map-lane q-${residue.charge > 0 ? "pos" : residue.charge < 0 ? "neg" : "0"}`, laneGlyph(residue), 2, row, {"aria-hidden": "true"});
    lane.title = `${residueDisplay(residue)} · ${residueClassName[chargeClass(residue)]}`;
    if (state.contact.sasa) addSasaLane(grid, residue, 3, row);
    slice.cols.forEach((columnIndex, column) => {
      const cell = cellByPosition.get(`${rowIndex}:${columnIndex}`);
      const element = addGridElement(grid, `map-cell${cell ? ` hit b${distanceBin(cell.min_dist_A)}` : ""}`, cell && state.contact.vals ? cell.min_dist_A.toFixed(2) : "", columns.positions[column], row, {"role": "gridcell", tabindex: cell ? "0" : "-1"});
      if (cell) {
        element.dataset.cellKey = `${residueKey(residue)}|${residueKey(data.cols[columnIndex])}`;
        element.dataset.cell = `${cell.row}:${cell.col}`;
        element.setAttribute("aria-label", `${residueDisplay(residue)} to ${residueDisplay(data.cols[columnIndex])}: ${fmt(cell.min_dist_A, 2)} angstrom`);
        if (cell.salt_bridges.length || cell.like_charges?.length || cell.hbond_like.length) {
          const badge = document.createElement("span");
          badge.className = `map-badge ${cell.salt_bridges.length ? "salt" : cell.like_charges?.length ? "like" : "hbond"}`;
          badge.textContent = cell.salt_bridges.length ? "+" : cell.like_charges?.length ? "×" : "H";
          element.appendChild(badge);
        }
        element.addEventListener("pointerover", () => showContactTip(element, cell, data));
        element.addEventListener("focus", () => showContactTip(element, cell, data));
        element.addEventListener("pointerout", hideContactTip);
        element.addEventListener("blur", hideContactTip);
        element.addEventListener("click", () => selectResidues([residueKey(residue), residueKey(data.cols[columnIndex])]));
        element.addEventListener("keydown", (event) => {
          if (event.key === "Enter" || event.key === " ") {
            event.preventDefault();
            element.click();
          }
        });
      }
    });
  });
  renderDistanceLegend();
  $("contact-map-note").textContent = `Cells retain all atom contacts through ${data.cutoff_A.toFixed(1)} Å. Wider gaps mark sequence gaps. The reported interface cutoff is ${data.reported_cutoff_A.toFixed(1)} Å; hollow ΔSASA lanes mark residues outside that report. Salt and like-charge contacts use charged functional-group atoms.`;
  renderContactMapTable(slice, data);
  renderSelectionState();
}

function bindMapAxisLabel(element, residue) {
  const choose = () => selectResidues([residueKey(residue)]);
  element.addEventListener("click", choose);
  element.addEventListener("keydown", (event) => {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      choose();
    }
  });
}

function renderContactMapTable(slice, data) {
  const tbody = $("contact-map-table");
  tbody.textContent = "";
  for (const cell of slice.cells) {
    const row = data.rows[cell.row];
    const col = data.cols[cell.col];
    const tr = document.createElement("tr");
    tr.className = "focus-row";
    tr.tabIndex = 0;
    tr.dataset.keys = `${residueKey(row)}|${residueKey(col)}`;
    tr.innerHTML = `<td>${esc(residueDisplay(row))}</td><td>${esc(residueClassName[chargeClass(row)])}</td><td>${chargeLabel(row)}</td><td>${fmt(row.dSASA_A2)}</td><td>${esc(residueDisplay(col))}</td><td>${esc(residueClassName[chargeClass(col)])}</td><td>${chargeLabel(col)}</td><td>${fmt(col.dSASA_A2)}</td><td>${fmt(cell.min_dist_A, 2)} Å</td><td>${cell.n_contacts}</td><td>${cell.n_tight}</td><td>${cell.hbond_like.length || "—"}</td><td>${cell.salt_bridges.length || "—"}</td><td>${cell.like_charges?.length || "—"}</td>`;
    bindSelectionRow(tr, tr.dataset.keys.split("|"));
    tbody.appendChild(tr);
  }
  if (!slice.cells.length) tbody.innerHTML = '<tr><td colspan="14" class="muted">No contacts match the selected filters.</td></tr>';
}

function renderPocketCard(item, label) {
  if (!item) return `<div class="pocket-card"><h3>${esc(label)} · no pocket detected</h3><p class="muted">The selected anchor region does not contain a qualifying component.</p></div>`;
  const rows = [
    ["Class", item.class], ["Volume", `${fmt(item.volume_A3)} Å³`], ["Depth", item.depth_A == null ? "n/a (closed)" : `${fmt(item.depth_A)} Å`],
    ["Enclosure", `${fmt(item.enclosure_fraction * 100)}%`], ["Opening area", `${fmt(item.opening_area_A2)} Å²`],
    ["Lining SASA", `${fmt(item.lining_sasa_A2)} Å²`], ["Solvent exposure", `${fmt(item.solvent_exposure_A2)} Å²`], ["Lining residues", item.lining_residues.length],
  ];
  return `<div class="pocket-card"><h3>${esc(label)} <span class="chip">${esc(item.class)}</span></h3>${rows.map(([name, value]) => `<div class="metric"><span>${esc(name)}</span><b>${esc(value)}</b></div>`).join("")}<details><summary>Show lining residues</summary><p class="muted">${item.lining_residues.map((residue) => esc(residueDisplay(residue))).join(", ") || "none"}</p></details></div>`;
}

function clearPocketHighlight() {
  if (state.pocketComponent && state.stage) {
    try { state.stage.removeComponent(state.pocketComponent); } catch (_) {}
  }
  state.pocketComponent = null;
  state.pocketHighlighted = false;
}

function stopPocketPolling() {
  if (state.pocketTimer) clearInterval(state.pocketTimer);
  state.pocketTimer = null;
  state.pocketJobId = null;
  state.pocketGeneration += 1;
}

function showPocketHighlight() {
  const item = state.pocketResult?.bound?.primary;
  if (!item?.mesh?.faces?.length || !state.stage || !window.NGL) return;
  clearPocketHighlight();
  try {
    const shape = new NGL.Shape("pocket-highlight");
    const colors = new Float32Array(item.mesh.vertices.length * 3);
    for (let index = 0; index < item.mesh.vertices.length; index += 1) {
      colors[index * 3] = 0.95;
      colors[index * 3 + 1] = 0.25;
      colors[index * 3 + 2] = 0.08;
    }
    shape.addMesh(new Float32Array(item.mesh.vertices.flat()), colors, new Uint32Array(item.mesh.faces.flat()), undefined, "bound pocket");
    state.pocketComponent = state.stage.addComponentFromObject(shape);
    state.pocketComponent.userChain = state.pocketResult.target_chain;
    state.pocketComponent.addRepresentation("buffer", {opacity: 0.2, side: "double", metalness: 0, roughness: 1});
    state.pocketHighlighted = true;
    updatePocketPosition();
    state.stage.viewer?.requestRender();
  } catch (error) {
    $("surface-note").textContent = `The pocket mesh could not be displayed: ${error.message || error}`;
  }
  updatePocketButtons();
}

function renderPocket(result) {
  state.pocketResult = result;
  state.pocketHighlightCancelled = false;
  clearPocketHighlight();
  const comparison = result.comparison || {};
  $("pocket-result").innerHTML = `<div class="pocket-cards">${renderPocketCard(result.free?.primary, "Free")}${renderPocketCard(result.bound?.primary, "Bound")}</div><div class="delta"><strong>Free → bound</strong><p class="muted">State change: ${esc(comparison.state_change || "not comparable")} · volume Δ ${comparison.delta_volume_A3 == null ? "n/a" : `${fmt(comparison.delta_volume_A3)} Å³`} · opening Δ ${comparison.delta_opening_area_A2 == null ? "n/a" : `${fmt(comparison.delta_opening_area_A2)} Å²`} · exposure Δ ${comparison.delta_solvent_exposure_A2 == null ? "n/a" : `${fmt(comparison.delta_solvent_exposure_A2)} Å²`}</p></div><div class="pocket-actions"><button id="show-pocket-highlight" class="secondary" type="button" ${result.bound?.primary?.mesh?.faces?.length ? "" : "disabled"}>Highlight bound pocket</button><button id="cancel-pocket-highlight" class="quiet" type="button" disabled>Cancel highlight</button></div><details class="raw"><summary>Calculation details</summary><pre>${esc(JSON.stringify(result, null, 2))}</pre></details>`;
  if (result.bound?.primary?.mesh?.faces?.length) showPocketHighlight();
  updatePocketButtons();
}

function updatePocketButtons() {
  const show = $("show-pocket-highlight");
  const cancel = $("cancel-pocket-highlight");
  if (!show || !cancel) return;
  show.disabled = !state.pocketResult?.bound?.primary?.mesh?.faces?.length || state.pocketHighlighted;
  cancel.disabled = !state.pocketHighlighted;
}

async function submitPocket(event) {
  event.preventDefault();
  const raw = $("pocket-anchor").value.split(":");
  const request = {
    target_chain: $("pocket-target").value,
    partner_chain: $("pocket-partner").value,
    anchor_residue: {number: Number(raw[0]), insertion_code: raw.slice(1).join(":")},
    radius_A: Number($("pocket-radius").value),
    grid_A: 0.6,
  };
  stopPocketPolling();
  clearPocketHighlight();
  state.pocketResult = null;
  state.pocketHighlightCancelled = false;
  $("pocket-result").innerHTML = '<p class="muted">Measuring the selected pocket…</p>';
  $("pocket-job-status").hidden = false;
  $("pocket-job-status").textContent = "queued";
  const analysisId = state.analysisId;
  const analysisGeneration = state.analysisGeneration;
  const pocketGeneration = state.pocketGeneration;
  try {
    const body = await getJSON(`/api/analyses/${analysisId}/pockets`, {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify(request),
    });
    if (analysisGeneration !== state.analysisGeneration || pocketGeneration !== state.pocketGeneration) return;
    pollPocket(body.job_id);
  } catch (error) {
    if (analysisGeneration !== state.analysisGeneration || pocketGeneration !== state.pocketGeneration) return;
    showError(error);
    $("pocket-job-status").textContent = "failed";
  }
}

function pollPocket(jobId) {
  stopPocketPolling();
  state.pocketJobId = jobId;
  const analysisId = state.analysisId;
  const pairId = state.pairId;
  const analysisGeneration = state.analysisGeneration;
  const pocketGeneration = state.pocketGeneration;
  let polling = false;
  const isCurrent = () => state.pocketJobId === jobId
    && state.analysisGeneration === analysisGeneration && state.pocketGeneration === pocketGeneration;
  state.pocketTimer = setInterval(async () => {
    if (polling || !isCurrent()) return;
    polling = true;
    try {
      const job = await getJSON(`/api/jobs/${jobId}`);
      if (!isCurrent()) return;
      $("pocket-job-status").textContent = `${job.status} · ${Math.round(job.progress * 100)}%`;
      if (job.status === "complete") {
        clearInterval(state.pocketTimer);
        state.pocketTimer = null;
        const result = await getJSON(`/api/jobs/${jobId}/result`);
        if (!isCurrent()) return;
        stopPocketPolling();
        if (state.analysisId === analysisId && state.pairId === pairId) renderPocket(result);
      }
      if (["failed", "interrupted", "cancelled"].includes(job.status)) stopPocketPolling();
    } catch (error) {
      if (!isCurrent()) return;
      stopPocketPolling();
      showError(error);
    } finally {
      polling = false;
    }
  }, 700);
}

function focusSurfaceHit(pick) {
  if (pick?.type !== "mesh" || !pick.component) return null;
  const name = pick.component.name;
  const chain = state.meshComponentChains.get(name);
  if (!chain) return null;
  const mesh = state.surface?.meshes?.[chain];
  const vertexMap = state.meshGeometry[chain]?.vertexMap;
  if (!Number.isInteger(pick.pid) || !vertexMap || pick.pid < 0 || pick.pid >= vertexMap.length) return null;
  const vertex = vertexMap[pick.pid];
  if (!mesh || !Number.isInteger(vertex) || vertex < 0 || vertex >= mesh.residue.length) return null;
  return {chain, label: mesh.residue[vertex], vertex, h: mesh.h?.[state.surfaceScale]?.[vertex]};
}

function showSurfaceTip(hit, event) {
  const tip = $("surface-tip");
  const name = chainLabel(hit.chain);
  tip.textContent = `${name} · ${hit.label} · local H ${fmtH(hit.h)}`;
  tip.classList.add("on");
  tip.setAttribute("aria-hidden", "false");
  if (event) {
    const box = tip.getBoundingClientRect();
    tip.style.left = `${Math.max(8, Math.min(event.clientX + 14, window.innerWidth - box.width - 8))}px`;
    tip.style.top = `${Math.max(8, Math.min(event.clientY + 14, window.innerHeight - box.height - 8))}px`;
  }
}

function hideSurfaceTip() {
  $("surface-tip").classList.remove("on");
  $("surface-tip").setAttribute("aria-hidden", "true");
}

function bindStageSignals() {
  if (!state.stage) return;
  state.stage.signals.hovered.add((pick) => {
    const hit = focusSurfaceHit(pick);
    if (hit) {
      const canvas = pick?.canvasPosition;
      const bounds = $("surface-view").getBoundingClientRect();
      showSurfaceTip(hit, canvas ? {clientX: bounds.left + canvas.x, clientY: bounds.bottom - canvas.y} : null);
    }
    else hideSurfaceTip();
  });
  state.stage.signals.clicked.add((pick) => {
    const atom = pick?.atom;
    if (atom && atom.chainname && aa3to1[atom.resname]) {
      selectResidues([`${atom.chainname}:${atom.resno}${atom.inscode || ""}`]);
      return;
    }
    const hit = focusSurfaceHit(pick);
    if (hit) selectResidues([surfaceResidueKey(hit.chain, hit.label)]);
  });
}

function initSurfaceEvents() {
  $("surface-opacity").addEventListener("input", (event) => {
    state.surfaceOpacity = Number(event.target.value);
    renderSurface();
  });
  $("surface-padding").addEventListener("input", (event) => {
    state.surfacePadding = Number(event.target.value);
    renderSurface();
  });
  $("surface-separation").addEventListener("input", (event) => {
    state.surfaceSeparation = Number(event.target.value);
    renderSurface();
  });
  $("surface-resolution").addEventListener("change", (event) => {
    state.surfaceScale = event.target.value;
    renderContactMap();
    renderSurface();
  });
  $("contact-cutoff").addEventListener("change", (event) => {
    state.contact.cut = Number(event.target.value);
    renderContactMap();
  });
  $("contact-type").addEventListener("change", (event) => {
    state.contact.type = event.target.value;
    renderContactMap();
  });
  $("contact-values").addEventListener("change", (event) => {
    state.contact.vals = event.target.checked;
    renderContactMap();
  });
  $("contact-sasa").addEventListener("change", (event) => {
    state.contact.sasa = event.target.checked;
    renderContactMap();
  });
}

$("upload-form").addEventListener("submit", submitUpload);
$("pdb-form").addEventListener("submit", submitPdb);
$("pocket-form").addEventListener("submit", submitPocket);
$("pocket-target").addEventListener("change", updateAnchors);
bindDropZone();
$("pdb-id").addEventListener("blur", inspectPdb);
$("refresh-history").addEventListener("click", refreshHistory);
$("clear-history").addEventListener("click", clearHistory);
$("drawer-toggle").addEventListener("click", () => {
  if (overlayQuery.matches) closeDrawer();
  else setDrawerCollapsed(!$("app-layout").classList.contains("drawer-collapsed"));
});
document.querySelectorAll(".drawer-open").forEach((button) => button.addEventListener("click", openDrawer));
$("drawer-backdrop").addEventListener("click", () => closeDrawer());
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape" && document.body.classList.contains("drawer-open-state")) closeDrawer();
});
try { if (localStorage.getItem("drawerCollapsed") === "1") $("app-layout").classList.add("drawer-collapsed"); } catch (_) {}
overlayQuery.addEventListener("change", syncDrawerMode);
syncDrawerMode();
$("interfaces-toggle").addEventListener("click", () => {
  const button = $("interfaces-toggle");
  const expanded = button.getAttribute("aria-expanded") === "true";
  button.setAttribute("aria-expanded", String(!expanded));
  button.textContent = expanded ? "+" : "−";
  button.title = expanded ? "Expand interface results" : "Collapse interface results";
  $("interfaces-body").hidden = expanded;
  $("interfaces-sidebar").classList.toggle("collapsed", expanded);
});
$("residue-tab").addEventListener("click", () => { state.activeTab = "residue"; renderTabs(); });
$("contact-tab").addEventListener("click", () => { state.activeTab = "contact"; renderTabs(); });
$("pocket-result").addEventListener("click", (event) => {
  if (event.target.id === "show-pocket-highlight") {
    state.pocketHighlightCancelled = false;
    showPocketHighlight();
    updatePocketButtons();
  }
  if (event.target.id === "cancel-pocket-highlight") {
    clearPocketHighlight();
    state.pocketHighlightCancelled = true;
    updatePocketButtons();
  }
});

initSurfaceEvents();
renderDistanceLegend();
refreshHistory();
