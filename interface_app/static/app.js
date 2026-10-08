/* Browser-side orchestration. Scientific calculations stay in the Python workers. */
const state = {
  analysisGeneration: 0,
  refreshing: false,
  analysisId: null,
  detail: null,
  result: null,
  pairId: null,
  surface: null,
  surfaceUnavailable: false,
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
  surfaceViewTab: "interactive",
  faceRotationCoupled: true,
  faceView: null,
  faceGeneration: 0,
  faceCache: new Map(),
  faceLoading: false,
  faceError: null,
  selectionReps: [],
  pocketComponent: null,
  pocketResult: null,
  pocketHighlighted: false,
  pocketHighlightCancelled: false,
  pocketTimer: null,
  pocketJobId: null,
  pocketGeneration: 0,
  selectedResidues: [],
  highlightUnmatched: [],
  highlightMessage: "",
  highlightError: false,
  pocketMode: "anchor",
  pocketRadii: {anchor: 8, residues: 6},
  activeTab: "residue",
  surfaceScale: "6",
  surfaceOpacity: 0.55,
  surfacePadding: 2,
  surfaceSeparation: 0,
  visibleChains: {},
  contact: {cut: 5, type: "all", vals: true, sasa: true},
  viewerLoaded: false,
  curvatureLoaded: false,
  curvatureJobId: null,
  electrostatics: null,
  electrostaticsLoaded: false,
  electrostaticsJobId: null,
  electrostaticsStatus: "waiting",
  electrostaticsRetryable: false,
  electrostaticsRevision: 0,
  surfaceMode: "convexity",
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
  [0, [217, 95, 138]], [0.5, [250, 250, 250]], [1, [27, 158, 119]],
];
const convexitySurfaceColors = curvatureColors;
const electrostaticsColors = [
  [0, [33, 102, 172]], [0.5, [247, 247, 247]], [1, [178, 24, 43]],
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
  state.surfaceUnavailable = false;
  state.structureText = null;
  state.structureFormat = null;
  state.pocketResult = null;
  state.pocketHighlighted = false;
  state.pocketHighlightCancelled = false;
  state.selectedResidues = [];
  state.highlightUnmatched = [];
  state.highlightMessage = "";
  state.highlightError = false;
  state.pocketRadii = {anchor: 8, residues: 6};
  $("highlight-chain").textContent = "";
  $("highlight-residues").value = "";
  $("pocket-residues").value = "";
  setPocketMode("anchor", true);
  renderHighlightFeedback();
  state.curvatureLoaded = false;
  state.curvatureJobId = null;
  state.electrostatics = null;
  state.electrostaticsLoaded = false;
  state.electrostaticsJobId = null;
  state.electrostaticsStatus = "waiting";
  state.electrostaticsRetryable = false;
  state.electrostaticsRevision = 0;
  state.surfaceMode = "convexity";
  state.viewerLoaded = false;
  state.visibleChains = {};
  state.contact = {cut: 5, type: "all", vals: true, sasa: true};
  state.surfaceScale = "6";
  state.surfaceOpacity = 0.55;
  state.surfacePadding = 2;
  state.surfaceSeparation = 0;
  state.faceCache.clear();
  setSurfaceViewTab("interactive");
  $("contact-cutoff").value = "5";
  $("contact-type").value = "all";
  $("contact-values").checked = true;
  $("contact-sasa").checked = true;
  $("surface-resolution").value = "6";
  $("surface-opacity").value = "0.55";
  $("surface-padding").value = "2";
  $("surface-separation").value = "0";
  setElectrostaticsStatus("Preparing optional APBS potential…", "waiting", false);
  $("protein-stats").textContent = "";
  $("surface-pair-stats").textContent = "";
  $("curvature-table").textContent = "";
  $("surface-note").textContent = "";
  $("surface-selection").textContent = "Select a residue in the surface, contact map, or interface table.";
  $("interfaces").innerHTML = '<p class="muted">Waiting for interface results…</p>';
  $("chain-toggles").textContent = "";
  hideContactTip();
  hideSurfaceTip();
  $("pocket-result").innerHTML = '<p class="muted">Choose a residue anchor or residue list and submit a pocket analysis.</p>';
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
    const activeJobs = (analysis.jobs || []).some((job) => activeStatuses.includes(job.status));
    if (["complete", "failed", "cancelled", "interrupted"].includes(analysis.status) && !activeJobs) {
      if (state.timer) {
        clearInterval(state.timer);
        refreshHistory();
      }
      state.timer = null;
    } else if (!state.timer) {
      state.timer = setInterval(refreshAnalysis, 1200);
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
  renderHighlightChains();
  renderSurface();
}

function electrostaticsReady() {
  const data = state.electrostatics;
  const limit = Number(data?.report?.color_limit_kT_e);
  const chains = Object.keys(state.surface?.meshes || {});
  return data?.status === "complete" && Number.isFinite(limit) && limit > 0
    && chains.length > 0 && chains.every(chain =>
      data.meshes?.[chain]?.potential_kT_e?.length === state.surface.meshes[chain].position.length / 3)
    && (!data.curvature_job_id || data.curvature_job_id === state.curvatureJobId);
}

function electrostaticsJobRetryable() {
  const job = state.detail?.jobs?.find(item => item.id === state.electrostaticsJobId);
  return job && ["failed", "cancelled", "interrupted"].includes(job.status);
}

function setElectrostaticsStatus(text, status = "waiting", retryable = false) {
  state.electrostaticsStatus = status;
  state.electrostaticsRetryable = retryable;
  const host = $("electrostatics-status");
  const retry = $("electrostatics-retry");
  if (host) {
    host.textContent = text.length > 240 ? `${text.slice(0, 237)}…` : text;
    host.title = text;
    host.className = `surface-mode-status ${status}`;
  }
  if (retry) {
    retry.hidden = !retryable;
    retry.textContent = electrostaticsJobRetryable() ? "Retry" : "Reload";
  }
  if (renderSurfaceModeControls()) {
    hideSurfaceTip();
    renderSurface();
    if (state.faceView?.ready) updateFaceSurfaces(state.faceView);
  }
}

function renderSurfaceModeControls() {
  const convexity = $("surface-mode-convexity");
  const electro = $("surface-mode-electrostatics");
  if (!convexity || !electro) return;
  const ready = electrostaticsReady();
  const fellBack = state.surfaceMode === "electrostatics" && !ready;
  if (fellBack) state.surfaceMode = "convexity";
  convexity.classList.toggle("active", state.surfaceMode === "convexity");
  electro.classList.toggle("active", state.surfaceMode === "electrostatics");
  convexity.setAttribute("aria-pressed", String(state.surfaceMode === "convexity"));
  electro.setAttribute("aria-pressed", String(state.surfaceMode === "electrostatics"));
  electro.disabled = !ready;
  return fellBack;
}

function setSurfaceMode(mode) {
  if (mode !== "convexity" && mode !== "electrostatics") return;
  if (mode === "electrostatics" && !electrostaticsReady()) return;
  state.surfaceMode = mode;
  hideSurfaceTip();
  renderSurfaceModeControls();
  renderSurface();
  if (state.faceView?.ready) updateFaceSurfaces(state.faceView);
}

function renderElectrostaticsResultStatus() {
  const data = state.electrostatics;
  if (!data) return;
  if (data.status !== "complete") {
    setElectrostaticsStatus(data.reason || "Electrostatics is unavailable.", "failed", true);
  } else if (state.surface?.status !== "complete") {
    setElectrostaticsStatus("Waiting for surface mesh…", "waiting", false);
  } else if (!electrostaticsReady()) {
    setElectrostaticsStatus("The electrostatic potentials do not match this surface. Analyze the coordinates again to generate matching results.", "failed", true);
  } else {
    setElectrostaticsStatus("Ready · APBS potential", "complete", false);
  }
}

async function loadElectrostatics() {
  if (state.electrostaticsLoaded || !state.analysisId) return;
  state.electrostaticsLoaded = true;
  const analysisId = state.analysisId;
  const generation = state.analysisGeneration;
  const jobId = state.electrostaticsJobId;
  try {
    const result = await getJSON(`/api/analyses/${encodeURIComponent(analysisId)}/electrostatics`);
    if (generation !== state.analysisGeneration || jobId !== state.electrostaticsJobId) return;
    const limit = Number(result.report?.color_limit_kT_e);
    if (result.status === "complete" && (!Number.isFinite(limit) || limit <= 0
      || !Object.values(result.meshes || {}).length
      || Object.values(result.meshes).some(mesh => !Array.isArray(mesh.potential_kT_e)
        || mesh.potential_kT_e.some(value => !Number.isFinite(value))))) {
      throw new Error("The electrostatics result has invalid surface potentials.");
    }
    state.electrostatics = result;
    state.electrostaticsRevision += 1;
    renderElectrostaticsResultStatus();
    renderSurface();
  } catch (error) {
    if (generation !== state.analysisGeneration || jobId !== state.electrostaticsJobId) return;
    state.electrostaticsLoaded = false;
    setElectrostaticsStatus(error.message || "Electrostatics is unavailable.", "failed", true);
  }
}

async function retryElectrostatics() {
  const jobId = state.electrostaticsJobId;
  if (!jobId) return;
  if (!electrostaticsJobRetryable()) {
    state.electrostaticsLoaded = false;
    await loadElectrostatics();
    return;
  }
  const generation = state.analysisGeneration;
  $("electrostatics-retry").disabled = true;
  try {
    await getJSON(`/api/jobs/${encodeURIComponent(jobId)}/retry`, {method: "POST"});
    if (generation !== state.analysisGeneration) return;
    state.electrostatics = null;
    state.electrostaticsLoaded = false;
    setElectrostaticsStatus("Retrying APBS potential…", "waiting", false);
    await refreshAnalysis();
  } catch (error) {
    if (generation === state.analysisGeneration) showError(error);
  } finally {
    if (generation === state.analysisGeneration) $("electrostatics-retry").disabled = false;
  }
}

function renderJobs(jobs) {
  const core = [...jobs].reverse().find((job) => job.kind === "core");
  const curvature = [...jobs].reverse().find((job) => job.kind === "curvature");
  const electrostatics = [...jobs].reverse().find((job) => job.kind === "electrostatics");
  if (core && core.status === "running") {
    $("curvature").innerHTML = `<strong>Interface analysis running</strong><p class="muted">${esc(core.stage)} · ${Math.round(core.progress * 100)}% · calculating contacts and SASA…</p>`;
    setSurfaceStatus("waiting", "Surface queued");
  }
  if (core && core.status === "failed") {
    $("curvature").innerHTML = `<strong>Interface analysis failed</strong><p class="muted">${esc(core.error || "The interface job did not complete.")}</p>`;
    setSurfaceStatus("failed", "Surface unavailable");
  }
  if (curvature) {
    if (state.curvatureJobId !== curvature.id) {
      state.curvatureJobId = curvature.id;
      state.curvatureLoaded = false;
    }
    $("curvature").innerHTML = `<strong>Surface curvature</strong><p class="muted">${esc(curvature.stage)} · ${Math.round(curvature.progress * 100)}% · ${esc(curvature.status)}</p>`;
    if (curvature.status === "running" || curvature.status === "queued") {
      setSurfaceStatus("waiting", `Surface ${curvature.status}`);
    }
    if (curvature.result_path && !state.curvatureLoaded) loadCurvature();
    if (curvature.status === "failed") setSurfaceStatus("failed", "Surface unavailable");
  }
  if (!electrostatics) {
    state.electrostaticsJobId = null;
    setElectrostaticsStatus(activeStatuses.includes(state.detail?.status)
      ? "Waiting for surface calculation…"
      : "This saved analysis has no electrostatics. Analyze the coordinates again to generate it.", "waiting", false);
  } else {
    if (state.electrostaticsJobId !== electrostatics.id) {
      state.electrostaticsLoaded = false;
      state.electrostatics = null;
    }
    state.electrostaticsJobId = electrostatics.id;
    if (electrostatics.status === "queued" || electrostatics.status === "running") {
      setElectrostaticsStatus(`${electrostatics.stage} · ${Math.round(electrostatics.progress * 100)}%`, "waiting", false);
    }
    if (electrostatics.result_path && !state.electrostaticsLoaded) loadElectrostatics();
    if (electrostatics.status === "failed") {
      setElectrostaticsStatus(electrostatics.error || "Electrostatics is unavailable.", "failed", true);
    }
    if (["interrupted", "cancelled"].includes(electrostatics.status)) {
      setElectrostaticsStatus(electrostatics.error || `Electrostatics ${electrostatics.status}.`, "failed", true);
    }
  }
}

function setSurfaceStatus(status, text) {
  $("surface-status").textContent = text;
  $("surface-status").className = `status-pill ${status}`;
  state.surfaceUnavailable = status === "failed";
  if (state.faceView?.ready) updateFaceSurfaces(state.faceView);
}

async function loadCurvature() {
  state.curvatureLoaded = true;
  const analysisId = state.analysisId;
  const generation = state.analysisGeneration;
  const jobId = state.curvatureJobId;
  try {
    const surface = await getJSON(`/api/analyses/${analysisId}/curvature`);
    if (generation !== state.analysisGeneration || jobId !== state.curvatureJobId) return;
    state.surface = surface;
    if (state.surface.status !== "complete" || !state.surface.report) {
      setSurfaceStatus("failed", state.surface.reason || "Surface unavailable");
      $("curvature").innerHTML = `<strong>Surface curvature unavailable</strong><p class="muted">${esc(state.surface.reason || "The optional surface job did not produce a mesh.")}</p>`;
      renderHighlightFeedback();
      return;
    }
    setSurfaceStatus("complete", "Surface ready");
    renderElectrostaticsResultStatus();
    renderSurfaceControls();
    renderContactMap();
    renderSurface();
  } catch (error) {
    if (generation !== state.analysisGeneration || jobId !== state.curvatureJobId) return;
    state.curvatureLoaded = false;
    setSurfaceStatus("failed", "Surface unavailable");
    renderHighlightFeedback();
    showError(error);
  }
}

function cleanupViewer() {
  cleanupFaceView();
  clearPocketHighlight();
  for (const {component, representation} of state.selectionReps) {
    try { component.removeRepresentation(representation); } catch (_) {}
  }
  state.selectionReps = [];
  if (state.stage) disposeNglStage(state.stage);
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
  cleanupFaceView();
  state.pairId = id;
  state.activeTab = "residue";
  state.selectedResidues = [];
  state.highlightUnmatched = [];
  state.highlightError = false;
  state.highlightMessage = $("highlight-residues").value.trim() ? "Input has not been applied to this interface." : "";
  renderHighlightChains(pair.chain_a);
  setPocketFeedback("");
  state.pocketResult = null;
  state.pocketHighlightCancelled = false;
  stopPocketPolling();
  clearPocketHighlight();
  $("pocket-result").innerHTML = '<p class="muted">Choose a residue anchor or residue list and submit a pocket analysis.</p>';
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

function parseResidueSpec(text, defaultChain) {
  const chains = new Map((state.result?.metadata?.chains || []).map(chain => [chain.id, chain.residues || []]));
  const keys = new Set();
  const unmatched = [];
  const single = /^(-?\d+)([A-Za-z]?)$/;
  const range = /^(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)$/;
  let chain = defaultChain;
  const findIndex = (residues, spec) => {
    const match = spec.match(single);
    if (!match || !Number.isSafeInteger(Number(match[1]))) return -1;
    return residues.findIndex(item => item.number === Number(match[1]) && (item.insertion_code || "") === match[2]);
  };
  for (const token of String(text || "").split(/[,\s]+/).filter(Boolean)) {
    let spec = token;
    const colon = token.indexOf(":");
    if (colon !== -1) {
      chain = token.slice(0, colon);
      spec = token.slice(colon + 1);
    }
    const residues = chains.get(chain);
    if (!residues) {
      unmatched.push(token);
      continue;
    }
    // Also accept a spaced prefix, such as "A: 45, 46".
    if (!spec && colon !== -1) continue;
    const span = spec.match(range);
    if (span) {
      const start = findIndex(residues, span[1]);
      const end = findIndex(residues, span[2]);
      if (start === -1 || end < start) unmatched.push(token);
      else residues.slice(start, end + 1).forEach(item => keys.add(residueKey(item)));
    } else {
      const index = findIndex(residues, spec);
      if (index === -1) unmatched.push(token);
      else keys.add(residueKey(residues[index]));
    }
  }
  return {keys: [...keys], unmatched};
}

function renderHighlightChains(defaultChain) {
  const select = $("highlight-chain");
  const previous = select.value;
  const chains = state.result?.metadata?.chains || [];
  select.innerHTML = chains.map(chain => `<option value="${esc(chain.id)}">${esc(chainLabel(chain.id))}</option>`).join("");
  select.value = chains.some(chain => chain.id === defaultChain) ? defaultChain
    : chains.some(chain => chain.id === previous) ? previous : chains[0]?.id || "";
}

function setHighlightPanelCollapsed(collapsed) {
  $("surface-stage").classList.toggle("side-collapsed", collapsed);
  $("residue-highlight-body").hidden = collapsed;
  const button = $("residue-highlight-toggle");
  const label = collapsed ? "Expand residue highlights" : "Collapse residue highlights";
  button.setAttribute("aria-expanded", String(!collapsed));
  button.setAttribute("aria-label", label);
  button.title = label;
  button.textContent = collapsed ? "+" : "−";
  try { localStorage.setItem("highlightPanelCollapsed", collapsed ? "1" : "0"); } catch (_) {}
  requestAnimationFrame(() => state.stage?.handleResize());
}

function renderHighlightFeedback() {
  const parts = [];
  if (state.highlightMessage) parts.push(state.highlightMessage);
  if (state.selectedResidues.length) {
    parts.push(state.selectedResidues.length + " residue(s) selected.");
    if (state.surface?.status === "complete" && curvatureScale()) {
      const missing = state.selectedResidues.filter(key => {
        const part = selectionParts(key);
        return !wholeSurface(part?.chain)?.residues?.some(item => item.label === residueLabelFromKey(key));
      });
      if (missing.length) parts.push("No surface vertices assigned: " + missing.join(", ") + ".");
    } else {
      parts.push("Surface patches are unavailable; atom highlighting is available when the coordinate viewer is ready.");
    }
  } else if (!state.highlightMessage) {
    parts.push("Enter residue numbers to highlight their surface patches.");
  }
  if (state.highlightUnmatched.length) parts.push("Unmatched input: " + state.highlightUnmatched.join(", ") + ".");
  $("highlight-feedback").textContent = parts.join(" ");
  $("highlight-feedback").classList.toggle("input-error", state.highlightError || !!state.highlightUnmatched.length);
}

function applyResidueHighlight(event) {
  event?.preventDefault();
  const {keys, unmatched} = parseResidueSpec($("highlight-residues").value, $("highlight-chain").value);
  if (keys.length) {
    selectResidues(keys);
    state.highlightUnmatched = unmatched;
  } else {
    state.highlightUnmatched = unmatched;
    state.highlightError = true;
    state.highlightMessage = unmatched.length ? "No residues matched; the current selection was retained." : "Enter at least one residue number.";
  }
  renderHighlightFeedback();
}

function initResidueHighlightPanel() {
  $("residue-highlight-form").addEventListener("submit", applyResidueHighlight);
  $("residue-highlight-toggle").addEventListener("click", () => {
    setHighlightPanelCollapsed($("residue-highlight-toggle").getAttribute("aria-expanded") === "true");
  });
  $("highlight-clear").addEventListener("click", () => {
    selectResidues([]);
    state.highlightMessage = "Selection cleared; input retained.";
    renderHighlightFeedback();
  });
  for (const id of ["highlight-residues", "highlight-chain"]) {
    $(id).addEventListener(id === "highlight-chain" ? "change" : "input", () => {
      state.highlightUnmatched = [];
      state.highlightError = false;
      state.highlightMessage = "Edited input has not been applied.";
      renderHighlightFeedback();
    });
  }
  let collapsed = false;
  try { collapsed = localStorage.getItem("highlightPanelCollapsed") === "1"; } catch (_) {}
  setHighlightPanelCollapsed(collapsed);
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
  setPocketFeedback("");
}

function selectResidues(keys) {
  state.selectedResidues = [...new Set(keys)];
  state.highlightUnmatched = [];
  state.highlightMessage = "";
  state.highlightError = false;
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

function interpolateColor(value, palette = curvatureColors) {
  let index = 0;
  while (index < palette.length - 2 && value > palette[index + 1][0]) index += 1;
  const [a, first] = palette[index];
  const [b, second] = palette[index + 1];
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

function meshColors(chain, geometry = state.meshGeometry[chain]) {
  const mesh = state.surface?.meshes?.[chain];
  const report = curvatureReport();
  if (!mesh || !report) return new Float32Array();
  const electroMode = state.surfaceMode === "electrostatics";
  const electroMesh = state.electrostatics?.meshes?.[chain];
  const values = electroMode
    ? (electroMesh?.potential_kT_e || [])
    : (mesh.h?.[state.surfaceScale] || []);
  if (electroMode && !electrostaticsReady()) return new Float32Array();
  const colors = new Float32Array((geometry?.vertexMap.length || values.length) * 3);
  const limit = electroMode
    ? Number(state.electrostatics.report?.color_limit_kT_e) || 1
    : Number(report.color_limit_Ainv) || 1;
  const palette = electroMode ? electrostaticsColors : convexitySurfaceColors;
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
      ? interpolateColor(Math.max(0, Math.min(1, (values[index] / limit + 1) / 2)), palette)
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

function installSurfacePicking(shape, geometry) {
  // NGL interpolates vertex picking colors across each triangle. Distinct
  // vertex IDs can therefore decode to an unrelated vertex anywhere in the
  // mesh. Use a separate picking geometry with a constant ID per triangle;
  // leave the indexed display mesh and its smooth normals untouched.
  const position = new Float32Array(geometry.index.length * 3);
  const primitiveId = new Float32Array(geometry.index.length);
  for (let corner = 0; corner < geometry.index.length; corner += 1) {
    const vertex = geometry.index[corner];
    position.set(geometry.position.subarray(vertex * 3, vertex * 3 + 3), corner * 3);
    primitiveId[corner] = Math.floor(corner / 3);
  }
  const picking = new NGL.MeshBuffer({position, primitiveId,
    color: new Float32Array(0), normal: new Float32Array(0)});
  const buffer = shape.bufferList[0];
  const getPickingMesh = buffer.getPickingMesh.bind(buffer);
  buffer.getPickingMesh = () => {
    const mesh = getPickingMesh();
    mesh.geometry = picking.geometry;
    return mesh;
  };
  const dispose = buffer.dispose.bind(buffer);
  buffer.dispose = () => { picking.dispose(); dispose(); };
}

function installSurfaceMesh(chain) {
  const geometry = state.meshGeometry[chain];
  const shapeName = `surface-${encodeURIComponent(chain)}`;
  const shape = new NGL.Shape(shapeName);
  shape.addMesh(geometry.position, meshColors(chain), geometry.index, geometry.normal, `${chain} surface`);
  installSurfacePicking(shape, geometry);
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
  const key = JSON.stringify([
    state.pairId, state.surfaceScale, state.surfacePadding, state.surfaceMode,
    state.electrostaticsRevision, state.electrostatics?.version, state.meshSelectionKey,
  ]);
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
  $("surface-note").textContent = state.surfaceMode === "electrostatics" && electrostaticsReady()
    ? `Blue is negative potential; red is positive potential. APBS/PDB2PQR uses AMBER charges at pH 7.4 on the full protein assembly. Colors are clipped at ±${fmt(state.electrostatics.report.color_limit_kT_e, 3)} kT/e for display; numeric values are not clipped. The curvature resolution control applies to the curvature summaries. Contact padding and separation affect only the display.`
    : `Pink is concave/inward; green is convex/outward (positive H). Colors are clipped at ±${fmt(report.color_limit_Ainv, 3)} Å⁻¹ for display; numeric values are not clipped. Contact padding changes only the colored display, while separation changes only viewing positions.`;
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
  renderSurfaceModeControls();
  const electroMode = state.surfaceMode === "electrostatics" && electrostaticsReady();
  const colorRamp = $("surface-color-ramp");
  colorRamp?.classList.toggle("electrostatics", electroMode);
  colorRamp?.classList.toggle("convexity", !electroMode);
  if (electroMode) {
    const limit = Number(state.electrostatics.report?.color_limit_kT_e) || 0;
    $("surface-color-min").textContent = `Negative −${limit.toFixed(3)}`;
    $("surface-color-max").textContent = `Positive +${limit.toFixed(3)}`;
    $("surface-color-label").textContent = "Electrostatic potential φ (kT/e)";
  } else if (report) {
    const limit = Number(report.color_limit_Ainv) || 0;
    $("surface-color-min").textContent = `Concave −${limit.toFixed(3)}`;
    $("surface-color-max").textContent = `Convex +${limit.toFixed(3)}`;
    $("surface-color-label").textContent = `H (Å⁻¹) · ${state.surfaceScale} Å neighborhoods`;
  }
  renderSurfaceStats();
  renderCurvatureTable();
  selectedSummary();
  renderHighlightFeedback();
  if (state.surfaceViewTab === "faces") ensureFaceView();
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

function setPocketFeedback(message, error = false) {
  const feedback = $("pocket-feedback");
  feedback.textContent = message;
  feedback.hidden = !message;
  feedback.classList.toggle("input-error", error);
}

function setPocketMode(mode, reset = false) {
  if (!["anchor", "residues"].includes(mode)) return;
  const radius = Number($("pocket-radius").value);
  if (!reset && Number.isFinite(radius) && radius >= 4 && radius <= 12) {
    state.pocketRadii[state.pocketMode] = radius;
  }
  state.pocketMode = mode;
  $("pocket-mode-anchor").checked = mode === "anchor";
  $("pocket-mode-residues").checked = mode === "residues";
  $("pocket-anchor-field").hidden = mode !== "anchor";
  $("pocket-anchor").disabled = mode !== "anchor";
  $("pocket-residue-field").hidden = mode !== "residues";
  $("pocket-residues").disabled = mode !== "residues";
  $("pocket-form").classList.toggle("residue-mode", mode === "residues");
  $("pocket-radius-label").textContent = mode === "residues" ? "Distance from residues" : "Radius";
  $("pocket-radius").value = state.pocketRadii[mode];
  setPocketFeedback("");
}

function pocketResidueRequest() {
  const {keys, unmatched} = parseResidueSpec($("pocket-residues").value, $("pocket-target").value);
  if (unmatched.length) throw new Error("Unknown or invalid pocket residues: " + unmatched.join(", ") + ".");
  if (!keys.length) throw new Error("Enter at least one residue to define the pocket.");
  if (keys.length > 60) throw new Error("Pocket definitions are limited to 60 residues.");
  const allowed = new Set([$("pocket-target").value, $("pocket-partner").value]);
  const outside = keys.filter(key => !allowed.has(selectionParts(key)?.chain));
  if (outside.length) throw new Error("Pocket residues must belong to the target or partner chain: " + outside.join(", ") + ".");
  return keys.map(key => {
    const part = selectionParts(key);
    return {chain_id: part.chain, number: Number(part.number), insertion_code: part.insertion};
  });
}

function usePocketSelection() {
  if (!state.selectedResidues.length) {
    setPocketFeedback("Select or highlight residues before using the current selection.", true);
    return;
  }
  setPocketMode("residues");
  $("pocket-residues").value = state.selectedResidues.join(", ");
  try {
    const residues = pocketResidueRequest();
    setPocketFeedback("Copied " + residues.length + " selected residue(s).");
  } catch (error) {
    setPocketFeedback(error.message, true);
  }
}

function renderPocketCard(item, label) {
  if (!item) return `<div class="pocket-card"><h3>${esc(label)} · no pocket detected</h3><p class="muted">The selected region does not contain a qualifying component.</p></div>`;
  const rows = [
    ["Class", item.class], ["Volume", `${fmt(item.volume_A3)} Å³`], ["Depth", item.depth_A == null ? "n/a (closed)" : `${fmt(item.depth_A)} Å`],
    ["Enclosure", `${fmt(item.enclosure_fraction * 100)}%`], ["Opening area", `${fmt(item.opening_area_A2)} Å²`],
    ["Lining SASA", `${fmt(item.lining_sasa_A2)} Å²`], ["Solvent exposure", `${fmt(item.solvent_exposure_A2)} Å²`], ["Lining residues", item.lining_residues.length],
  ];
  if (item.defined_residues_present != null) {
    rows.push(["Defined residues lining", item.defined_residues_present
      ? item.defined_residues_lining.length + " / " + item.defined_residues_present + " (" + fmtPercent(item.defined_residue_coverage) + ")"
      : "n/a (no defined residues present)"]);
  }
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
  const definition = result.mode === "residues"
    ? "Defined residues: " + (result.pocket_residues || []).map(residueDisplay).join(", ") + " · distance " + fmt(result.radius_A) + " Å"
    : "Anchor: " + residueDisplay(result.anchor_residue) + " · radius " + fmt(result.radius_A) + " Å";
  $("pocket-result").innerHTML = `<p class="pocket-definition muted">${esc(definition)}</p><div class="pocket-cards">${renderPocketCard(result.free?.primary, "Free")}${renderPocketCard(result.bound?.primary, "Bound")}</div><div class="delta"><strong>Free → bound</strong><p class="muted">State change: ${esc(comparison.state_change || "not comparable")} · volume Δ ${comparison.delta_volume_A3 == null ? "n/a" : `${fmt(comparison.delta_volume_A3)} Å³`} · opening Δ ${comparison.delta_opening_area_A2 == null ? "n/a" : `${fmt(comparison.delta_opening_area_A2)} Å²`} · exposure Δ ${comparison.delta_solvent_exposure_A2 == null ? "n/a" : `${fmt(comparison.delta_solvent_exposure_A2)} Å²`}</p></div><div class="pocket-actions"><button id="show-pocket-highlight" class="secondary" type="button" ${result.bound?.primary?.mesh?.faces?.length ? "" : "disabled"}>Highlight bound pocket</button><button id="cancel-pocket-highlight" class="quiet" type="button" disabled>Cancel highlight</button></div><details class="raw"><summary>Calculation details</summary><pre>${esc(JSON.stringify(result, null, 2))}</pre></details>`;
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
  const request = {
    target_chain: $("pocket-target").value,
    partner_chain: $("pocket-partner").value,
    radius_A: Number($("pocket-radius").value),
    grid_A: 0.6,
  };
  try {
    if (!request.target_chain || !request.partner_chain || request.target_chain === request.partner_chain) {
      throw new Error("Choose different target and partner chains.");
    }
    if (!Number.isFinite(request.radius_A) || request.radius_A < 4 || request.radius_A > 12) {
      throw new Error("Pocket search distance must be between 4 and 12 Å.");
    }
    if (state.pocketMode === "residues") {
      request.pocket_residues = pocketResidueRequest();
    } else {
      const raw = $("pocket-anchor").value.split(":");
      if (!raw[0]) throw new Error("Choose an anchor residue.");
      request.anchor_residue = {number: Number(raw[0]), insertion_code: raw.slice(1).join(":")};
    }
  } catch (error) {
    setPocketFeedback(error.message, true);
    return;
  }
  setPocketFeedback("");
  clearError();
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

function surfaceHitVertex(pick, mesh, geometry) {
  if (!pick?.component || !mesh || !geometry || !Number.isInteger(pick.pid)
    || pick.pid < 0 || pick.pid * 3 + 2 >= geometry.index.length) return null;
  const cursor = pick.canvasPosition;
  if (!cursor || !Number.isFinite(cursor.x) || !Number.isFinite(cursor.y)) return null;
  const controls = pick.component.stage.viewerControls;
  const point = new NGL.Vector3();
  let closest = null;
  let minDistance = Infinity;
  for (let corner = 0; corner < 3; corner += 1) {
    const displayed = geometry.index[pick.pid * 3 + corner];
    const vertex = geometry.vertexMap[displayed];
    if (vertex >= mesh.residue.length) continue;
    point.fromArray(geometry.position, displayed * 3);
    if (pick.instance) point.applyMatrix4(pick.instance.matrix);
    point.applyMatrix4(pick.component.matrix);
    const screen = controls.getPositionOnCanvas(point);
    const distance = (screen.x - cursor.x) ** 2 + (screen.y - cursor.y) ** 2;
    if (distance < minDistance) { closest = vertex; minDistance = distance; }
  }
  return closest;
}

function focusSurfaceHit(pick) {
  if (pick?.type !== "mesh" || !pick.component) return null;
  const name = pick.component.name;
  const chain = state.meshComponentChains.get(name);
  if (!chain || pick.component !== state.meshComponents[chain]) return null;
  const mesh = state.surface?.meshes?.[chain];
  const vertex = surfaceHitVertex(pick, mesh, state.meshGeometry[chain]);
  if (vertex === null) return null;
  return surfaceHitDetail(chain, mesh, vertex);
}

function surfaceHitDetail(chain, mesh, vertex) {
  return {
    chain,
    label: mesh.residue[vertex],
    vertex,
    h: mesh.h?.[state.surfaceScale]?.[vertex],
    potential: state.electrostatics?.meshes?.[chain]?.potential_kT_e?.[vertex],
  };
}

function showSurfaceTip(hit, event) {
  const tip = $("surface-tip");
  const name = chainLabel(hit.chain);
  tip.textContent = state.surfaceMode === "electrostatics" && electrostaticsReady()
    ? `${name} · ${hit.label} · potential ${hit.potential == null ? "n/a" : `${hit.potential >= 0 ? "+" : ""}${Number(hit.potential).toFixed(3)} kT/e`}`
    : `${name} · ${hit.label} · local H ${fmtH(hit.h)}`;
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

function setSurfaceViewTab(tab) {
  state.surfaceViewTab = tab;
  const faces = tab === "faces";
  clearFaceHover();
  hideSurfaceTip();
  for (const [name, selected] of [["interactive", !faces], ["faces", faces]]) {
    const button = $(`surface-${name}-tab`);
    button.classList.toggle("active", selected);
    button.setAttribute("aria-selected", String(selected));
    button.tabIndex = selected ? 0 : -1;
    $(`surface-${name}-panel`).hidden = !selected;
  }
  $("surface-separation-control").hidden = faces;
  $("chain-toggles").hidden = faces;
  if (!faces) $("surface-view-help").textContent = "Drag to rotate, scroll to zoom, or select a residue in the surface, contact map, or tables.";
  renderFaceRotationControls();
  if (faces) {
    ensureFaceView();
    requestAnimationFrame(() => fitFaceViews(false));
  } else {
    requestAnimationFrame(() => state.stage?.handleResize());
  }
}

function faceViewMessage(message, retry = false) {
  $("face-view-status").textContent = message;
  $("face-view-retry").hidden = !retry;
}

function renderFaceRotationControls() {
  const coupled = state.faceRotationCoupled;
  $("face-rotation-coupled").checked = coupled;
  for (const panel of state.faceView?.panels || []) {
    panel.host.setAttribute("aria-label", `Binding face of ${chainLabel(panel.chain)} with ${coupled ? "mirrored coupled" : "independent"} rotation`);
  }
  if (state.surfaceViewTab === "faces") {
    const action = coupled ? "rotate both binding faces with mirrored motion" : "rotate that binding face independently";
    $("surface-view-help").textContent = `Drag either partner to ${action}; scroll or use + / − to zoom. Table selections appear here; hover a residue to find its closest partner.`;
  }
}

function setFaceRotationCoupled(coupled) {
  const changed = state.faceRotationCoupled !== coupled;
  state.faceRotationCoupled = coupled;
  const face = state.faceView;
  if (changed && coupled && face?.ready) applyFaceRotation(face);
  renderFaceRotationControls();
  if (face?.ready) updateFaceSurfaces(face);
}

function disposeNglStage(stage) {
  // This bundled NGL version leaves observation and animation loops running
  // after Stage.dispose(). Stop them and detach input behaviors first.
  if (stage.mouseObserver) {
    stage.mouseObserver.setParameters({hoverTimeout: -1});
    stage.mouseObserver._listen = () => {};
  }
  stage.viewer.animate = () => {};
  stage.viewer.render = () => {};
  stage.viewer.requestRender = () => {};
  for (const name of ["pickingBehavior", "mouseBehavior", "keyBehavior", "animationBehavior", "mouseObserver"]) {
    stage[name]?.dispose();
  }
  // KeyBehavior.dispose in this NGL build removes the wrong down/up callbacks.
  const keyboard = stage.keyBehavior;
  if (keyboard) {
    keyboard.domElement.removeEventListener("keydown", keyboard._onKeydown);
    keyboard.domElement.removeEventListener("keyup", keyboard._onKeyup);
  }
  stage.signals.hovered.removeAll();
  stage.signals.clicked.removeAll();
  try { stage.removeAllComponents(); } catch (_) {}
  try { stage.dispose(); } catch (_) {}
  stage.tooltip?.remove();
  stage.viewer.wrapper?.remove();
}

function cleanupFaceView() {
  state.faceGeneration += 1;
  state.faceLoading = false;
  state.faceError = null;
  const face = state.faceView;
  state.faceView = null;
  if (face) {
    face.observer?.disconnect();
    if (face.hoverFrame) cancelAnimationFrame(face.hoverFrame);
    if (face.resizeFrame) cancelAnimationFrame(face.resizeFrame);
    for (const panel of face.panels) {
      for (const [event, handler] of panel.handlers) panel.host.removeEventListener(event, handler);
      if (panel.rotationChanged) panel.stage.viewerControls.signals.changed.remove(panel.rotationChanged);
      if (panel.hoverRefreshFrame) cancelAnimationFrame(panel.hoverRefreshFrame);
      disposeNglStage(panel.stage);
      panel.host.replaceChildren();
    }
  }
  $("face-view-grid").hidden = true;
  $("face-hover-readout").textContent = "Hover a residue to highlight its closest partner.";
  faceViewMessage("Choose an interface to view its binding faces.");
}

function faceTransform(structure, frame) {
  const normal = new NGL.Vector3(...frame.normal);
  const up = new NGL.Vector3(...frame.up);
  const right = new NGL.Vector3().crossVectors(normal, up);
  // NGL's camera is on -Z. Map the outward normal to -Z and up to +Y,
  // using a proper rotation so neither molecule is reflected.
  const matrix = new NGL.Matrix4().set(
    right.x, right.y, right.z, 0,
    up.x, up.y, up.z, 0,
    -normal.x, -normal.y, -normal.z, 0,
    0, 0, 0, 1,
  );
  const bounds = new NGL.Box3().makeEmpty();
  const point = new NGL.Vector3();
  structure.eachAtom(atom => bounds.expandByPoint(point.set(atom.x, atom.y, atom.z).applyMatrix4(matrix)));
  if (bounds.isEmpty()) throw new Error("The partner chain has no displayable atoms.");
  const center = bounds.getCenter(new NGL.Vector3());
  matrix.setPosition(center.clone().negate());
  // Reserve room for the molecular surface before it arrives, keeping zoom
  // stable when cartoons gain their curvature mesh.
  bounds.translate(center.negate()).expandByScalar(3);
  return {matrix, bounds};
}

function createFacePanel(face, chain, side) {
  const host = $(`face-view-${side}`);
  const stage = new NGL.Stage(host, {backgroundColor: "#eef2f1", sampleLevel: 0,
    cameraType: "orthographic", tooltip: false});
  const panel = {stage, host, chain, side, handlers: [], hoverInside: false,
    lastRotation: stage.viewerControls.rotation.clone(),
    hoverRefreshFrame: null,
    mesh: null, meshComponent: null, meshRep: null, meshGeometry: null,
    meshSelectionKey: null, highlightKey: null, colorKey: null, lastWidth: null, lastHeight: null};
  // Register immediately so a later build failure can dispose every stage.
  face.panels.push(panel);
  stage.mouseControls.clear();
  stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
  stage.mouseControls.add("drag-ctrl-right", NGL.MouseActions.zRotateDrag);
  stage.mouseControls.add("scroll", NGL.MouseActions.zoomScroll);
  stage.keyControls.clear();
  stage.keyControls.disabled = true;
  panel.rotationChanged = () => syncFaceRotation(face, panel);
  stage.viewerControls.signals.changed.add(panel.rotationChanged);
  const view = state.structureComponents[chain].structure.getView(new NGL.Selection(`:${chain}`));
  const frame = face.geometry.chains[chain];
  const transform = faceTransform(view, frame);
  panel.transform = transform.matrix;
  panel.bounds = transform.bounds;
  // StructureViews own their listeners, not the original parsed coordinates.
  // Disposing this tab must leave Interactive 3D usable.
  panel.component = stage.addComponentFromObject(view, {name: `binding-face-${chain}`});
  panel.component.setTransform(panel.transform);
  const index = (state.result.metadata?.chains || []).findIndex(item => item.id === chain);
  panel.component.addRepresentation("cartoon", {color: cartoonColors[Math.max(0, index) % cartoonColors.length],
    quality: "medium", opacity: 1});
  panel.highlightRep = panel.component.addRepresentation("ball+stick", {
    sele: "none", color: "#ffff00", scale: 0.85,
  });
  $(`face-chain-${side}-label`).textContent = chainLabel(chain);
  host.setAttribute("aria-label", `Binding face of ${chainLabel(chain)} with ${state.faceRotationCoupled ? "mirrored coupled" : "independent"} rotation`);
  $(`face-plane-${side}`).textContent = `${frame.binding_residue_keys.length} binding residues · ${frame.estimated ? "estimated direction" : "fitted binding plane"}`;
  const enter = () => { panel.hoverInside = true; };
  const leave = () => { panel.hoverInside = false; queueFaceHover(face, null); hideSurfaceTip(); };
  for (const [event, handler] of [["mouseenter", enter], ["mouseleave", leave]]) {
    host.addEventListener(event, handler);
    panel.handlers.push([event, handler]);
  }
  stage.signals.hovered.add(pick => {
    if (state.faceView !== face || state.surfaceViewTab !== "faces" || !panel.hoverInside) return;
    queueFaceHover(face, facePickedResidue(panel, pick));
    const vertex = pick?.type === "mesh" && pick.component === panel.meshComponent
      ? surfaceHitVertex(pick, panel.mesh, panel.meshGeometry) : null;
    if (vertex !== null && pick.canvasPosition) {
      const bounds = host.getBoundingClientRect();
      showSurfaceTip(surfaceHitDetail(panel.chain, panel.mesh, vertex), {
        clientX: bounds.left + pick.canvasPosition.x,
        clientY: bounds.bottom - pick.canvasPosition.y,
      });
    } else hideSurfaceTip();
  });
  return panel;
}

function mirroredFaceRotation(rotation) {
  // Corresponding sites appear at opposite X positions in the two outward
  // binding views. Conjugate by this reflection to keep their height and
  // depth aligned: X tilt is shared, while Y turns and Z roll are reversed.
  // This composes correctly for mixed rotations and is its own inverse.
  return new NGL.Quaternion(rotation.x, -rotation.y, -rotation.z, rotation.w);
}

function syncFaceRotation(face, source) {
  if (state.faceView !== face || !face.ready) return;
  const sourceRotation = source.stage.viewerControls.rotation;
  // Zoom also emits changed. Record only rotation changes, including while
  // uncoupled, so re-coupling uses the last pane the user actually rotated.
  if (source.lastRotation.equals(sourceRotation)) return;
  source.lastRotation.copy(sourceRotation);
  if (face.syncingRotation) return;
  face.rotation.copy(source.side === "b" ? mirroredFaceRotation(sourceRotation) : sourceRotation);
  if (state.faceRotationCoupled) applyFaceRotation(face, source);
}

function applyFaceRotation(face, source = null) {
  face.syncingRotation = true;
  try {
    for (const panel of face.panels) {
      if (panel !== source) panel.stage.viewerControls.rotate(
        panel.side === "b" ? mirroredFaceRotation(face.rotation) : face.rotation);
    }
  } finally {
    face.syncingRotation = false;
  }
}

async function ensureFaceView() {
  if (state.surfaceViewTab !== "faces" || !state.result) return;
  const pair = state.result.pairs.find(item => item.id === state.pairId);
  if (!pair) {
    faceViewMessage("Choose an interface to view its binding faces.");
    return;
  }
  if (state.faceView?.pairId === pair.id && state.faceView.ready) {
    updateFaceSurfaces(state.faceView);
    return;
  }
  if (state.faceLoading || state.faceError) return;
  if (!state.structureComponents[pair.chain_a] || !state.structureComponents[pair.chain_b]) {
    faceViewMessage("Loading partner coordinates…");
    if (!state.viewerLoaded) loadStructureViewer();
    return;
  }
  const generation = state.faceGeneration;
  const analysisId = state.analysisId;
  state.faceLoading = true;
  faceViewMessage("Calculating binding planes and closest residues…");
  try {
    let geometry = state.faceCache.get(pair.id);
    if (!geometry) {
      geometry = await getJSON(`/api/analyses/${encodeURIComponent(analysisId)}/pairs/${encodeURIComponent(pair.id)}/face-view`);
      if (generation !== state.faceGeneration || analysisId !== state.analysisId) return;
      state.faceCache.set(pair.id, geometry);
    }
    if (state.surfaceViewTab !== "faces") return;
    if (!geometry.available) {
      state.faceError = geometry.reason || "Binding faces are unavailable for this interface.";
      faceViewMessage(state.faceError);
      return;
    }
    const face = {pairId: pair.id, geometry, panels: [], ready: false,
      rotation: new NGL.Quaternion(), syncingRotation: false,
      hoverKey: null, hoverKeys: [], pendingHover: null, hoverFrame: null, resizeFrame: null};
    state.faceView = face;
    $("face-view-grid").hidden = false;
    createFacePanel(face, pair.chain_a, "a");
    createFacePanel(face, pair.chain_b, "b");
    face.ready = true;
    updateFaceSurfaces(face);
    fitFaceViews(true);
    face.observer = new ResizeObserver(() => {
      if (face.resizeFrame) cancelAnimationFrame(face.resizeFrame);
      face.resizeFrame = requestAnimationFrame(() => {
        face.resizeFrame = null;
        if (state.faceView === face) fitFaceViews(false);
      });
    });
    face.panels.forEach(panel => face.observer.observe(panel.host));
  } catch (error) {
    if (generation !== state.faceGeneration || analysisId !== state.analysisId) return;
    cleanupFaceView();
    state.faceError = error.message || String(error);
    faceViewMessage(`Binding faces could not load: ${state.faceError}`, true);
  } finally {
    if (generation === state.faceGeneration) state.faceLoading = false;
  }
}

function fitFaceViews(force) {
  const face = state.faceView;
  if (!face?.ready || state.surfaceViewTab !== "faces") return;
  const sizes = face.panels.map(panel => panel.host.getBoundingClientRect());
  if (sizes.some(size => size.width <= 0 || size.height <= 0)) return;
  const changed = face.panels.some((panel, index) =>
    panel.lastWidth !== sizes[index].width || panel.lastHeight !== sizes[index].height);
  if (!force && !changed) return;
  let unitsPerPixel = 0;
  for (let i = 0; i < face.panels.length; i += 1) {
    const panel = face.panels[i];
    panel.stage.handleResize();
    const rotation = new NGL.Matrix4().makeRotationFromQuaternion(panel.stage.viewerControls.rotation);
    const extent = panel.bounds.clone().applyMatrix4(rotation).getSize(new NGL.Vector3());
    unitsPerPixel = Math.max(unitsPerPixel, extent.x / sizes[i].width, extent.y / sizes[i].height);
  }
  for (let i = 0; i < face.panels.length; i += 1) {
    const panel = face.panels[i];
    panel.lastWidth = sizes[i].width;
    panel.lastHeight = sizes[i].height;
    const fov = panel.stage.viewer.perspectiveCamera.fov * Math.PI / 180;
    const distance = sizes[i].height * unitsPerPixel * 1.1 / (2 * Math.tan(fov / 2));
    panel.stage.viewerControls.distance(-Math.max(1, distance));
    panel.stage.viewer.requestRender();
  }
}

function facePickedResidue(panel, pick) {
  let key = null;
  let atom = pick?.atom;
  if (!atom && pick?.bond) {
    // NGL's closestBondAtom helper projects untransformed coordinates. Apply
    // the fixed component frame before deciding which end was hovered.
    const atoms = [pick.bond.atom1, pick.bond.atom2];
    const distance = candidate => {
      const point = candidate.positionToVector3(new NGL.Vector3());
      if (pick.instance) point.applyMatrix4(pick.instance.matrix);
      point.applyMatrix4(pick.component.matrix);
      return panel.stage.viewerControls.getPositionOnCanvas(point).distanceTo(pick.canvasPosition);
    };
    atom = pick.canvasPosition && distance(atoms[1]) < distance(atoms[0]) ? atoms[1] : atoms[0];
  }
  if (atom && atom.chainname === panel.chain && aa3to1[atom.resname]) {
    key = `${atom.chainname}:${atom.resno}${atom.inscode || ""}`;
  } else if (pick?.type === "mesh" && pick.component === panel.meshComponent) {
    const vertex = surfaceHitVertex(pick, panel.mesh, panel.meshGeometry);
    if (vertex !== null) key = surfaceResidueKey(panel.chain, panel.mesh.residue[vertex]);
  }
  return state.faceView?.geometry.nearest[panel.chain]?.[key] ? key : null;
}

function queueFaceHover(face, key) {
  face.pendingHover = key;
  if (face.hoverFrame) return;
  face.hoverFrame = requestAnimationFrame(() => {
    face.hoverFrame = null;
    if (state.faceView === face && state.surfaceViewTab === "faces") applyFaceHover(face, face.pendingHover);
  });
}

function clearFaceHover() {
  const face = state.faceView;
  if (!face) return;
  if (face.hoverFrame) cancelAnimationFrame(face.hoverFrame);
  face.hoverFrame = null;
  face.pendingHover = null;
  face.panels.forEach(panel => { panel.hoverInside = false; });
  applyFaceHover(face, null);
}

function applyFaceHover(face, key) {
  if (face.hoverKey === key) return;
  face.hoverKey = key;
  const chain = selectionParts(key)?.chain;
  const match = chain && face.geometry.nearest[chain]?.[key];
  face.hoverKeys = match ? [key, match.partner_key] : [];
  updateFaceSurfaces(face);
  renderFaceHoverReadout(face);
}

function renderFaceHoverReadout(face) {
  const readout = $("face-hover-readout");
  const chain = selectionParts(face.hoverKey)?.chain;
  const match = chain && face.geometry.nearest[chain]?.[face.hoverKey];
  if (!match) {
    readout.textContent = "Hover a residue to highlight its closest partner.";
    return;
  }
  const keys = face.panels.map(panel => face.hoverKeys.find(key => selectionParts(key)?.chain === panel.chain));
  const labels = keys.map(key => `${selectionParts(key).chain}:${residueLabelFromKey(key)}`);
  const missing = face.panels.filter((panel, index) => panel.mesh && !panel.meshLabels.has(residueLabelFromKey(keys[index])));
  readout.textContent = `${labels.join(" ↔ ")} · ${fmt(match.distance_A, 2)} Å (minimum heavy-atom distance)`
    + (missing.length ? ` · No surface vertices for ${missing.map(panel => panel.chain).join(", ")}; atoms highlighted.` : "");
}

function installFaceMesh(panel) {
  const shape = new NGL.Shape(`face-surface-${encodeURIComponent(panel.chain)}`);
  const geometry = panel.meshGeometry;
  shape.addMesh(geometry.position, meshColors(panel.chain, geometry), geometry.index, geometry.normal, `${panel.chain} binding face`);
  installSurfacePicking(shape, geometry);
  panel.meshComponent = panel.stage.addComponentFromObject(shape);
  panel.meshComponent.setTransform(panel.transform);
  panel.meshRep = panel.meshComponent.addRepresentation("buffer", {
    opacity: state.surfaceOpacity, side: "double", metalness: 0, roughness: 1,
  });
}

function updateFaceSurfaces(face) {
  if (!face.ready) return;
  for (const panel of face.panels) {
    const keys = [...new Set([...state.selectedResidues, ...face.hoverKeys])]
      .filter(key => selectionParts(key)?.chain === panel.chain).sort();
    const highlightKey = JSON.stringify(keys);
    if (panel.highlightKey !== highlightKey) {
      panel.highlightRep.setSelection(keys.map(key => `(${selectionString(key)})`).join(" or ") || "none");
      panel.highlightKey = highlightKey;
    }
    const mesh = state.surface?.status === "complete" ? state.surface.meshes?.[panel.chain] : null;
    const labels = new Set(keys.map(residueLabelFromKey).filter(Boolean));
    const selectionKey = JSON.stringify([...labels].sort());
    const changed = panel.mesh !== mesh || panel.meshSelectionKey !== selectionKey;
    if (changed) {
      if (panel.mesh !== mesh) panel.meshLabels = new Set(mesh?.residue || []);
      if (panel.meshComponent) panel.stage.removeComponent(panel.meshComponent);
      panel.meshComponent = null;
      panel.meshRep = null;
      panel.meshGeometry = null;
      panel.mesh = mesh;
      panel.meshSelectionKey = selectionKey;
      panel.colorKey = null;
      if (mesh?.index?.length) {
        panel.meshGeometry = surfaceMeshGeometry(mesh, labels);
        installFaceMesh(panel);
        const component = panel.meshComponent;
        panel.stage.tasks.onZeroOnce(() => {
          // NGL picks once after mouse movement. A pick during a mesh rebuild
          // can miss the surface, so refresh it once the replacement is ready.
          if (state.faceView !== face || panel.meshComponent !== component) return;
          if (panel.hoverRefreshFrame) cancelAnimationFrame(panel.hoverRefreshFrame);
          panel.hoverRefreshFrame = requestAnimationFrame(() => {
            panel.hoverRefreshFrame = null;
            if (state.faceView === face && state.surfaceViewTab === "faces"
              && panel.meshComponent === component && panel.hoverInside
              && panel.stage.mouseObserver.overElement) panel.stage.mouseObserver.hovering = false;
          });
        });
      }
    }
    const colorKey = JSON.stringify([
      state.surfaceScale, state.surfacePadding, state.surfaceMode,
      state.electrostaticsRevision, state.electrostatics?.version, selectionKey,
    ]);
    if (panel.meshComponent && panel.colorKey !== colorKey) {
      panel.meshComponent.object.bufferList[0].setAttributes({color: meshColors(panel.chain, panel.meshGeometry)});
      panel.colorKey = colorKey;
    }
    panel.meshRep?.setParameters({opacity: state.surfaceOpacity});
    panel.stage.viewer.requestRender();
  }
  const meshCount = face.panels.filter(panel => panel.meshComponent).length;
  const rotationMode = state.faceRotationCoupled ? "Mirrored rotation" : "Independent rotation";
  faceViewMessage(meshCount === 2
    ? `${rotationMode} · independent zoom · table selections highlighted · closest residues measured in the original complex.`
    : state.surface || state.surfaceUnavailable
      ? `${rotationMode} · surfaces unavailable; selected atoms are highlighted; hover the cartoons to find closest residues.`
      : `${rotationMode} · selected atoms are highlighted; hover the cartoons while partner surfaces load.`);
  renderFaceHoverReadout(face);
}

function initFaceViewEvents() {
  $("face-rotation-coupled").addEventListener("change", event => setFaceRotationCoupled(event.target.checked));
  renderFaceRotationControls();
  for (const [name, tab] of [["interactive", "interactive"], ["faces", "faces"]]) {
    const button = $(`surface-${name}-tab`);
    button.addEventListener("click", () => setSurfaceViewTab(tab));
    button.addEventListener("keydown", event => {
      if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
      event.preventDefault();
      const next = event.key === "Home" ? "interactive" : event.key === "End" ? "faces"
        : tab === "faces" ? "interactive" : "faces";
      setSurfaceViewTab(next);
      $(`surface-${next}-tab`).focus();
    });
  }
  $("face-view-retry").addEventListener("click", () => {
    state.faceCache.delete(state.pairId);
    cleanupFaceView();
    ensureFaceView();
  });
  document.querySelectorAll("[data-face-zoom]").forEach(button => {
    button.addEventListener("click", () => {
      const panel = state.faceView?.panels.find(item => item.side === button.dataset.faceZoom);
      panel?.stage.viewerControls.zoom(button.dataset.zoom === "in" ? 0.15 : -0.18);
    });
  });
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
  $("surface-mode-convexity").addEventListener("click", () => setSurfaceMode("convexity"));
  $("surface-mode-electrostatics").addEventListener("click", () => setSurfaceMode("electrostatics"));
  $("electrostatics-retry").addEventListener("click", retryElectrostatics);
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
$("pocket-partner").addEventListener("change", () => setPocketFeedback(""));
$("pocket-residues").addEventListener("input", () => setPocketFeedback(""));
$("pocket-use-selection").addEventListener("click", usePocketSelection);
document.querySelectorAll('input[name="pocket-mode"]').forEach(input => {
  input.addEventListener("change", () => setPocketMode(input.value));
});
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
initFaceViewEvents();
initResidueHighlightPanel();
setPocketMode("anchor");
renderDistanceLegend();
refreshHistory();
