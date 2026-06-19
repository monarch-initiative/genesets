const state = {
  bundles: [],
  bundleId: null,
  facets: null,
  rows: [],
  selectedQueryId: null,
  activeTab: "diffs",
};

const els = {
  bundleSelect: document.querySelector("#bundleSelect"),
  refreshButton: document.querySelector("#refreshButton"),
  bundleSummary: document.querySelector("#bundleSummary"),
  searchInput: document.querySelector("#searchInput"),
  sourceFilter: document.querySelector("#sourceFilter"),
  stratumFilter: document.querySelector("#stratumFilter"),
  sortSelect: document.querySelector("#sortSelect"),
  queryCount: document.querySelector("#queryCount"),
  statusText: document.querySelector("#statusText"),
  queryTableBody: document.querySelector("#queryTableBody"),
  emptyState: document.querySelector("#emptyState"),
  detailContent: document.querySelector("#detailContent"),
  detailSource: document.querySelector("#detailSource"),
  detailTitle: document.querySelector("#detailTitle"),
  detailDescription: document.querySelector("#detailDescription"),
  detailMetrics: document.querySelector("#detailMetrics"),
  diffClassFilter: document.querySelector("#diffClassFilter"),
  diffCount: document.querySelector("#diffCount"),
  diffTableBody: document.querySelector("#diffTableBody"),
  oldResultsBody: document.querySelector("#oldResultsBody"),
  newResultsBody: document.querySelector("#newResultsBody"),
  geneList: document.querySelector("#geneList"),
};

function fmtInt(value) {
  if (value === null || value === undefined) return "0";
  return Number(value).toLocaleString();
}

function fmtMaybeInt(value) {
  if (value === null || value === undefined) return "";
  return Number(value).toLocaleString();
}

function fmtFloat(value, digits = 2) {
  if (value === null || value === undefined) return "";
  return Number(value).toFixed(digits);
}

function fmtP(value) {
  if (value === null || value === undefined) return "";
  const n = Number(value);
  if (!Number.isFinite(n)) return "";
  if (n === 0) return "0";
  if (n < 0.001) return n.toExponential(2);
  return n.toPrecision(3);
}

function classBadge(value) {
  const label = String(value || "").replaceAll("_", " ");
  let cls = "";
  if (value === "lost_significant") cls = "lost";
  if (value === "gained_significant") cls = "gained";
  if (value === "shared_significant") cls = "shared";
  return `<span class="badge ${cls}">${escapeHtml(label || "unknown")}</span>`;
}

function changeNote(term) {
  if (!term || !term.target_id) return "";
  if (term.delta !== null && term.delta !== undefined) {
    return `delta ${fmtFloat(term.delta, 2)}`;
  }
  return String(term.class || "").replaceAll("_", " ");
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) {
    const body = await response.text();
    throw new Error(`${response.status} ${response.statusText}: ${body}`);
  }
  return await response.json();
}

function currentBundle() {
  return state.bundles.find((bundle) => bundle.id === state.bundleId);
}

function renderBundleOptions() {
  els.bundleSelect.innerHTML = state.bundles
    .map((bundle) => {
      const label = bundle.analysis || bundle.id;
      return `<option value="${escapeHtml(bundle.id)}">${escapeHtml(label)}</option>`;
    })
    .join("");
  els.bundleSelect.value = state.bundleId;
}

function renderBundleSummary() {
  const bundle = currentBundle();
  if (!bundle) {
    els.bundleSummary.innerHTML = "";
    return;
  }
  const counts = bundle.result_counts || {};
  const timings = bundle.timings || {};
  const items = [
    ["Query sets", bundle.selected_count],
    ["A targets", counts.old?.rows],
    ["B targets", counts.new?.rows],
    ["Lost / gained", `${fmtInt(counts.diff?.lost_significant)} / ${fmtInt(counts.diff?.gained_significant)}`],
    ["Runtime", timings.total_seconds ? `${fmtFloat(timings.total_seconds, 1)}s` : ""],
  ];
  els.bundleSummary.innerHTML = items
    .map(
      ([label, value]) => `
        <div class="summary-item">
          <span class="summary-label">${escapeHtml(label)}</span>
          <span class="summary-value">${escapeHtml(value ?? "")}</span>
        </div>
      `,
    )
    .join("");
}

async function loadFacets() {
  state.facets = await fetchJson(`/api/bundles/${encodeURIComponent(state.bundleId)}/facets`);
  renderFacetSelect(els.sourceFilter, "All source families", state.facets.source_class || []);
  renderFacetSelect(els.stratumFilter, "All strata", state.facets.stratum || []);
}

function renderFacetSelect(select, emptyLabel, rows) {
  const current = select.value;
  select.innerHTML = `<option value="">${escapeHtml(emptyLabel)}</option>`;
  for (const row of rows) {
    const option = document.createElement("option");
    option.value = row.value;
    option.textContent = `${row.value} (${fmtInt(row.count)})`;
    select.appendChild(option);
  }
  if ([...select.options].some((option) => option.value === current)) {
    select.value = current;
  }
}

let queryTimer = null;

function scheduleLoadQuerySets() {
  clearTimeout(queryTimer);
  queryTimer = setTimeout(loadQuerySets, 180);
}

async function loadQuerySets() {
  if (!state.bundleId) return;
  els.statusText.textContent = "Loading";
  const params = new URLSearchParams({
    sort: els.sortSelect.value,
    order: els.sortSelect.value === "name" ? "asc" : "desc",
    limit: "150",
  });
  if (els.searchInput.value.trim()) params.set("search", els.searchInput.value.trim());
  if (els.sourceFilter.value) params.set("source_class", els.sourceFilter.value);
  if (els.stratumFilter.value) params.set("stratum", els.stratumFilter.value);
  const data = await fetchJson(`/api/bundles/${encodeURIComponent(state.bundleId)}/query-sets?${params}`);
  state.rows = data.rows || [];
  els.queryCount.textContent = `${fmtInt(data.total)} gene sets`;
  els.statusText.textContent = data.total > data.rows.length ? `Showing ${fmtInt(data.rows.length)}` : "";
  renderQueryTable();
  const visibleIds = new Set(state.rows.map((row) => row.id));
  if (state.rows.length && (!state.selectedQueryId || !visibleIds.has(state.selectedQueryId))) {
    selectQuery(state.rows[0].id);
  } else if (state.selectedQueryId) {
    highlightSelectedRow();
  } else {
    els.emptyState.classList.remove("hidden");
    els.detailContent.classList.add("hidden");
  }
}

function renderQueryTable() {
  els.queryTableBody.innerHTML = state.rows
    .map((row) => {
      const top = row.top_crossing_term || row.top_changed_term || {};
      const source = row.source_class || "unknown";
      return `
        <tr data-query-id="${escapeHtml(row.id)}">
          <td>
            <div class="gene-title">${escapeHtml(row.name || row.id)}</div>
            <div class="subtext">${escapeHtml(row.id)}</div>
          </td>
          <td>
            <span class="badge">${escapeHtml(source.replace(/^msigdb_/, ""))}</span>
            <div class="subtext">${escapeHtml(row.stratum || "")}</div>
            <div class="subtext">${fmtInt(row.gene_count)} genes</div>
          </td>
          <td class="num">${fmtInt(row.lost_specific_count)}<div class="subtext">${fmtInt(row.lost_count)} all</div></td>
          <td class="num">${fmtInt(row.gained_specific_count)}<div class="subtext">${fmtInt(row.gained_count)} all</div></td>
          <td>
            <div class="term-title">${escapeHtml(top.target_name || "")}</div>
            <div class="subtext">${escapeHtml(top.target_id || "")} ${escapeHtml(changeNote(top))}</div>
          </td>
        </tr>
      `;
    })
    .join("");
  for (const tr of els.queryTableBody.querySelectorAll("tr")) {
    tr.addEventListener("click", () => selectQuery(tr.dataset.queryId));
  }
  highlightSelectedRow();
}

function highlightSelectedRow() {
  for (const tr of els.queryTableBody.querySelectorAll("tr")) {
    tr.classList.toggle("selected", tr.dataset.queryId === state.selectedQueryId);
  }
}

async function selectQuery(queryId) {
  if (!queryId) return;
  state.selectedQueryId = queryId;
  highlightSelectedRow();
  const params = new URLSearchParams({ query_id: queryId });
  const detail = await fetchJson(`/api/bundles/${encodeURIComponent(state.bundleId)}/query-set?${params}`);
  renderDetail(detail);
  await Promise.all([loadDiffs(), loadResults("old"), loadResults("new")]);
}

function renderDetail(detail) {
  els.emptyState.classList.add("hidden");
  els.detailContent.classList.remove("hidden");
  els.detailSource.textContent = `${detail.source_class || "unknown"} / ${detail.stratum || "unknown"}`;
  els.detailTitle.textContent = detail.name || detail.id;
  els.detailDescription.textContent = detail.description || detail.id;
  const metrics = [
    ["Genes", detail.gene_count],
    ["Lost", detail.lost_count],
    ["Specific lost", detail.lost_specific_count],
    ["Gained", detail.gained_count],
    ["Shared", detail.shared_count],
  ];
  els.detailMetrics.innerHTML = metrics
    .map(
      ([label, value]) => `
        <div class="metric">
          <span class="metric-value">${fmtInt(value)}</span>
          <span class="metric-label">${escapeHtml(label)}</span>
        </div>
      `,
    )
    .join("");
  const genes = detail.genes || [];
  els.geneList.innerHTML = genes.length
    ? genes.map((gene) => `<span class="gene-chip">${escapeHtml(gene)}</span>`).join("")
    : `<span class="muted">No query genes were available in the report bundle.</span>`;
}

async function loadDiffs() {
  if (!state.selectedQueryId) return;
  const params = new URLSearchParams({
    query_id: state.selectedQueryId,
    limit: "300",
  });
  if (els.diffClassFilter.value) params.set("class_filter", els.diffClassFilter.value);
  const data = await fetchJson(`/api/bundles/${encodeURIComponent(state.bundleId)}/query-set/diffs?${params}`);
  els.diffCount.textContent = `${fmtInt(data.total)} rows`;
  els.diffTableBody.innerHTML = (data.rows || [])
    .map((row) => {
      const targetSize =
        row.class === "gained_significant" ? row.right_target_size : row.left_target_size;
      return `
        <tr>
          <td>${classBadge(row.class)}</td>
          <td>
            <div class="term-title">${escapeHtml(row.target_name || "")}</div>
            <div class="subtext">${escapeHtml(row.target_id || "")}</div>
          </td>
          <td class="num">${fmtP(row.left_p_adjust)}</td>
          <td class="num">${fmtP(row.right_p_adjust)}</td>
          <td class="num">${fmtFloat(row.delta_neg_log10_p_adjust, 2)}</td>
          <td class="num">${fmtMaybeInt(targetSize)}</td>
        </tr>
      `;
    })
    .join("");
}

async function loadResults(run) {
  if (!state.selectedQueryId) return;
  const params = new URLSearchParams({
    query_id: state.selectedQueryId,
    run,
    limit: "200",
  });
  const data = await fetchJson(`/api/bundles/${encodeURIComponent(state.bundleId)}/query-set/results?${params}`);
  const body = run === "old" ? els.oldResultsBody : els.newResultsBody;
  body.innerHTML = (data.rows || [])
    .map(
      (row) => `
        <tr>
          <td>
            <div class="term-title">${escapeHtml(row.target_name || "")}</div>
            <div class="subtext">${escapeHtml(row.target_id || "")}</div>
          </td>
          <td class="num">${fmtMaybeInt(row.overlap)}</td>
          <td class="num">${fmtMaybeInt(row.target_size)}</td>
          <td class="num">${fmtP(row.p_adjust_bonferroni)}</td>
        </tr>
      `,
    )
    .join("");
}

function activateTab(tabName) {
  state.activeTab = tabName;
  for (const tab of document.querySelectorAll(".tab")) {
    tab.classList.toggle("active", tab.dataset.tab === tabName);
  }
  for (const panel of document.querySelectorAll(".tab-panel")) {
    panel.classList.toggle("active", panel.id === `${tabName}Panel`);
  }
}

async function loadBundles() {
  const data = await fetchJson("/api/bundles");
  state.bundles = data.bundles || [];
  if (!state.bundles.length) throw new Error("No report bundles were loaded.");
  state.bundleId = state.bundleId || state.bundles[0].id;
  renderBundleOptions();
  renderBundleSummary();
  await loadFacets();
  await loadQuerySets();
}

function bindEvents() {
  els.bundleSelect.addEventListener("change", async () => {
    state.bundleId = els.bundleSelect.value;
    state.selectedQueryId = null;
    renderBundleSummary();
    await loadFacets();
    await loadQuerySets();
  });
  els.refreshButton.addEventListener("click", () => loadBundles().catch(showError));
  els.searchInput.addEventListener("input", scheduleLoadQuerySets);
  els.sourceFilter.addEventListener("change", loadQuerySets);
  els.stratumFilter.addEventListener("change", loadQuerySets);
  els.sortSelect.addEventListener("change", loadQuerySets);
  els.diffClassFilter.addEventListener("change", loadDiffs);
  for (const tab of document.querySelectorAll(".tab")) {
    tab.addEventListener("click", () => activateTab(tab.dataset.tab));
  }
}

function showError(error) {
  console.error(error);
  els.statusText.textContent = "Error";
  els.queryTableBody.innerHTML = `
    <tr>
      <td colspan="5">
        <div class="term-title">Could not load explorer data</div>
        <div class="subtext">${escapeHtml(error.message || error)}</div>
      </td>
    </tr>
  `;
}

bindEvents();
loadBundles().catch(showError);
