// ============================================================
// volcano.js
// Purpose: Interactive volcano plot with gene search,
//          NS toggle, and highlight functionality
// ============================================================

const DATA_PATH = "../data/processed/";

// ── CSV PARSER ───────────────────────────────────────────────
function parseCSV(text) {
  const lines   = text.trim().split("\n").filter(l => l.trim() !== "");
  const headers = lines[0].split(",").map(h => h.trim().replace(/"/g, ""));
  return lines.slice(1).map(line => {
    const values = line.split(",").map(v => v.trim().replace(/"/g, ""));
    const obj = {};
    headers.forEach((h, i) => {
      obj[h] = isNaN(values[i]) ? values[i] : parseFloat(values[i]);
    });
    return obj;
  });
}

// ── RENDER VOLCANO ───────────────────────────────────────────
async function renderVolcano(region, fcCutoff, showNS = true, searchGene = "") {

  const file = region === "EC"
    ? "deg_Entorhinal_Cortex.csv"
    : "deg_Superior_Frontal_Gyrus.csv";

  const regionLabel = region === "EC"
    ? "Entorhinal Cortex"
    : "Superior Frontal Gyrus";

  const response = await fetch(DATA_PATH + file);
  const text     = await response.text();
  const data     = parseCSV(text);

  // ── CATEGORIZE ────────────────────────────────────────────
  const up       = { x: [], y: [], text: [], genes: [] };
  const down     = { x: [], y: [], text: [], genes: [] };
  const ns       = { x: [], y: [], text: [], genes: [] };
  const searched = { x: [], y: [], text: [], genes: [] };

  data.forEach(row => {
    const logp = -Math.log10(Math.max(row["adj.P.Val"], 1e-20));
    const fc   = row["logFC"];
    const gene = row["gene"];

    const hoverText = [
      `<b>${gene}</b>`,
      `logFC: ${fc.toFixed(3)}`,
      `adj.P: ${row["adj.P.Val"].toExponential(2)}`,
      `-log10(P): ${logp.toFixed(2)}`
    ].join("<br>");

    // If gene matches search — put in special highlighted category
    if (searchGene && gene.toUpperCase() === searchGene.toUpperCase()) {
      searched.x.push(fc); searched.y.push(logp);
      searched.text.push(hoverText); searched.genes.push(gene);
    } else if (row["adj.P.Val"] < 0.05 && fc > fcCutoff) {
      up.x.push(fc); up.y.push(logp); up.text.push(hoverText); up.genes.push(gene);
    } else if (row["adj.P.Val"] < 0.05 && fc < -fcCutoff) {
      down.x.push(fc); down.y.push(logp); down.text.push(hoverText); down.genes.push(gene);
    } else {
      ns.x.push(fc); ns.y.push(logp); ns.text.push(hoverText); ns.genes.push(gene);
    }
  });

  // ── COUNT DEGs for subtitle ───────────────────────────────
  const upCount   = up.x.length;
  const downCount = down.x.length;

  // ── TRACES ────────────────────────────────────────────────
  const traces = [];

  // NS trace — only if showNS is true
  if (showNS) {
    traces.push({
      type: "scatter", mode: "markers", name: `NS (n=${ns.x.length})`,
      x: ns.x, y: ns.y, text: ns.text,
      hovertemplate: "%{text}<extra></extra>",
      marker: { color: "#cbd5e1", size: 3, opacity: 0.25, line: { width: 0 } }
    });
  }

  // Downregulated
  traces.push({
    type: "scatter", mode: "markers",
    name: `Downregulated (n=${downCount})`,
    x: down.x, y: down.y, text: down.text,
    hovertemplate: "%{text}<extra></extra>",
    marker: { color: "#3b82f6", size: 5, opacity: 0.7, line: { width: 0 } }
  });

  // Upregulated
  traces.push({
    type: "scatter", mode: "markers",
    name: `Upregulated (n=${upCount})`,
    x: up.x, y: up.y, text: up.text,
    hovertemplate: "%{text}<extra></extra>",
    marker: { color: "#ef4444", size: 5, opacity: 0.7, line: { width: 0 } }
  });

  // Searched gene — highlighted in yellow
  if (searched.x.length > 0) {
    traces.push({
      type: "scatter", mode: "markers+text",
      name: `${searchGene} (highlighted)`,
      x: searched.x, y: searched.y,
      text: searched.genes,
      textposition: "top center",
      hovertemplate: "%{text}<extra></extra>",
      marker: {
        color: "#f59e0b", size: 14, opacity: 1,
        line: { color: "#92400e", width: 2 },
        symbol: "star"
      }
    });
  }

  // ── THRESHOLD LINES ───────────────────────────────────────
  const shapes = [
    {
      type: "line", x0: -12, x1: 12,
      y0: -Math.log10(0.05), y1: -Math.log10(0.05),
      line: { color: "#94a3b8", width: 1, dash: "dash" }
    },
    {
      type: "line", x0: fcCutoff, x1: fcCutoff, y0: 0, y1: 22,
      line: { color: "#94a3b8", width: 1, dash: "dash" }
    },
    {
      type: "line", x0: -fcCutoff, x1: -fcCutoff, y0: 0, y1: 22,
      line: { color: "#94a3b8", width: 1, dash: "dash" }
    }
  ];

  // ── LAYOUT ────────────────────────────────────────────────
  const layout = {
    title: {
      text: `<b>${regionLabel}</b> — AD vs normal<br><sup style="color:#94a3b8">↑ ${upCount} upregulated · ↓ ${downCount} downregulated</sup>`,
      font: { size: 13, color: "#0f172a" }
    },
    xaxis: {
      title: "Log₂ Fold Change",
      zeroline: false, gridcolor: "#f1f5f9",
      range: [-12, 12]
    },
    yaxis: {
      title: "-Log₁₀ (adj. P-value)",
      gridcolor: "#f1f5f9"
    },
    shapes,
    legend: { orientation: "h", yanchor: "bottom", y: 1.02, xanchor: "right", x: 1 },
    plot_bgcolor: "#ffffff", paper_bgcolor: "#ffffff",
    margin: { t: 80, r: 20, b: 60, l: 60 },
    hovermode: "closest"
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"],
    toImageButtonOptions: { format: "png", filename: `volcano_${region}` }
  };

  Plotly.newPlot("volcano-plot", traces, layout, config);
}