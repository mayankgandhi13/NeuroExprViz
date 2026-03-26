// ============================================================
// summary.js
// Purpose: DEG summary bar chart + GSEA bubble chart
// ============================================================

// ── DEG SUMMARY BAR CHART ────────────────────────────────────
// Shows upregulated and downregulated DEG counts
// for both EC and SFG side by side

async function renderSummary() {

  // Load both region CSVs simultaneously using Promise.all
  // This is faster than loading them one after the other
  const [ecResp, sfgResp] = await Promise.all([
    fetch(DATA_PATH + "deg_Entorhinal_Cortex.csv"),
    fetch(DATA_PATH + "deg_Superior_Frontal_Gyrus.csv")
  ]);

  const ecText  = await ecResp.text();
  const sfgText = await sfgResp.text();
  const ec      = parseCSV(ecText);
  const sfg     = parseCSV(sfgText);

  // ── COUNT DEGs ─────────────────────────────────────────────
  // Count up and down DEGs for each region
  // using the same thresholds as R pipeline

  function countDEGs(data, fcCutoff = 1) {
    let up   = 0;
    let down = 0;
    data.forEach(row => {
      if (row["adj.P.Val"] < 0.05 && row["logFC"] >  fcCutoff) up++;
      if (row["adj.P.Val"] < 0.05 && row["logFC"] < -fcCutoff) down++;
    });
    return { up, down };
  }

  const ecCounts  = countDEGs(ec);
  const sfgCounts = countDEGs(sfg);

  // ── BUILD TRACES ────────────────────────────────────────────
  // Two traces — one for upregulated, one for downregulated
  // Each trace has two bars — one per region

  const upTrace = {
    type: "bar",
    name: "Upregulated",
    x: ["Entorhinal Cortex", "Superior Frontal Gyrus"],
    y: [ecCounts.up, sfgCounts.up],
    marker: {
      color: "#ef4444",
      opacity: 0.85
    },
    // Show count on top of each bar
    text: [ecCounts.up, sfgCounts.up],
    textposition: "outside",
    textfont: { size: 11, color: "#0f172a" }
  };

  const downTrace = {
    type: "bar",
    name: "Downregulated",
    x: ["Entorhinal Cortex", "Superior Frontal Gyrus"],
    y: [ecCounts.down, sfgCounts.down],
    marker: {
      color: "#3b82f6",
      opacity: 0.85
    },
    text: [ecCounts.down, sfgCounts.down],
    textposition: "outside",
    textfont: { size: 11, color: "#0f172a" }
  };

  const layout = {
    title: {
      text: "<b>DEG Counts</b> per Brain Region",
      font: { size: 13, color: "#0f172a" }
    },
    barmode: "group",          // grouped bars side by side
    xaxis: {
      tickfont: { size: 10 },
      gridcolor: "#f1f5f9"
    },
    yaxis: {
      title: "Number of DEGs",
      gridcolor: "#f1f5f9"
    },
    legend: {
      orientation: "h",
      yanchor: "bottom",
      y: 1.02,
      xanchor: "right",
      x: 1
    },
    plot_bgcolor:  "#ffffff",
    paper_bgcolor: "#ffffff",
    margin: { t: 60, r: 20, b: 80, l: 60 }
  };

  const config = {
    responsive: true,
    displayModeBar: false      // hide toolbar for this simple chart
  };

  Plotly.newPlot("summary-plot", [upTrace, downTrace], layout, config);
}

// ── GSEA BUBBLE CHART ─────────────────────────────────────────
// Shows enriched pathways as bubbles
// X axis  = NES (Normalized Enrichment Score)
//           positive NES = pathway upregulated in AD
//           negative NES = pathway downregulated in AD
// Y axis  = -log10(padj) — higher = more significant
// Size    = number of leading edge genes
// Color   = direction (red = activated, blue = suppressed)

async function renderGSEA(region) {

  const file = region === "EC"
    ? "gsea_EC.csv"
    : "gsea_SFG.csv";

  const regionLabel = region === "EC"
    ? "Entorhinal Cortex"
    : "Superior Frontal Gyrus";

  const response = await fetch(DATA_PATH + file);
  const text     = await response.text();
  const data     = parseCSV(text);

  // Filter to significant pathways only
  const sig = data.filter(row => row["padj"] < 0.05);

  // Clean pathway names — remove "HALLMARK_" prefix for readability
  // e.g. "HALLMARK_INFLAMMATORY_RESPONSE" → "Inflammatory Response"
  function cleanName(name) {
    return name
      .replace("HALLMARK_", "")
      .replace(/_/g, " ")
      .toLowerCase()
      .replace(/\b\w/g, c => c.toUpperCase()); // Title Case
  }

  // Separate into activated (positive NES) and suppressed (negative NES)
  const activated  = sig.filter(row => row["NES"] > 0);
  const suppressed = sig.filter(row => row["NES"] < 0);

  // Build traces for activated and suppressed separately
  function buildTrace(rows, name, color) {
    return {
      type: "scatter",
      mode: "markers",
      name: name,
      x: rows.map(r => r["NES"]),
      y: rows.map(r => -Math.log10(r["padj"])),
      text: rows.map(r => {
        const edges = typeof r["leadingEdge"] === "string"
          ? r["leadingEdge"].split(";").length
          : 0;
        return [
          `<b>${cleanName(r["pathway"])}</b>`,
          `NES: ${r["NES"].toFixed(3)}`,
          `padj: ${r["padj"].toExponential(2)}`,
          `Leading edge genes: ${edges}`
        ].join("<br>");
      }),
      hovertemplate: "%{text}<extra></extra>",
      marker: {
        color: color,
        // Bubble size based on leading edge gene count
        size: rows.map(r => {
          const edges = typeof r["leadingEdge"] === "string"
            ? r["leadingEdge"].split(";").length
            : 10;
          // Scale size — min 8, max 30
          return Math.min(30, Math.max(8, edges / 5));
        }),
        opacity: 0.75,
        line: { color: "#ffffff", width: 1 }
      }
    };
  }

  const traces = [
    buildTrace(activated,  "Activated in AD",  "#ef4444"),
    buildTrace(suppressed, "Suppressed in AD", "#3b82f6")
  ];

  // Add vertical line at NES = 0
  const shapes = [{
    type: "line",
    x0: 0, x1: 0,
    y0: 0, y1: 20,
    line: { color: "#94a3b8", width: 1, dash: "dot" }
  }];

  const layout = {
    title: {
      text: `<b>GSEA Hallmark Pathways</b> — ${regionLabel}`,
      font: { size: 13, color: "#0f172a" }
    },
    xaxis: {
      title: "Normalized Enrichment Score (NES)",
      zeroline: false,
      gridcolor: "#f1f5f9"
    },
    yaxis: {
      title: "-Log₁₀ (padj)",
      gridcolor: "#f1f5f9"
    },
    shapes: shapes,
    legend: {
      orientation: "h",
      yanchor: "bottom",
      y: 1.02,
      xanchor: "right",
      x: 1
    },
    plot_bgcolor:  "#ffffff",
    paper_bgcolor: "#ffffff",
    margin: { t: 60, r: 20, b: 60, l: 70 },
    hovermode: "closest"
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"],
    toImageButtonOptions: {
      format: "png",
      filename: `gsea_${region}`
    }
  };

  Plotly.newPlot("gsea-plot", traces, layout, config);
}