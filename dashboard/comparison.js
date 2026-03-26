// ============================================================
// comparison.js
// Purpose: Renders the JS version of each plot inside
//          the tool comparison tabs
// ============================================================

// Renders the Plotly.js version of a plot inside the
// comparison tab — called when user switches to JS tab
async function renderComparisonPlot(plotType) {

  const containerId = "comparison-volcano";

  if (plotType === "volcano") {
    const response = await fetch(DATA_PATH + "deg_Entorhinal_Cortex.csv");
    const data     = parseCSV(await response.text());

    const up   = { x: [], y: [], text: [] };
    const down = { x: [], y: [], text: [] };
    const ns   = { x: [], y: [], text: [] };

    data.forEach(row => {
      const logp = -Math.log10(Math.max(row["adj.P.Val"], 1e-20));
      const fc   = row["logFC"];
      const gene = row["gene"];
      const ht   = `<b>${gene}</b><br>logFC: ${fc.toFixed(2)}<br>adj.P: ${row["adj.P.Val"].toExponential(2)}`;

      if (row["adj.P.Val"] < 0.05 && fc > 1)       { up.x.push(fc);   up.y.push(logp);   up.text.push(ht); }
      else if (row["adj.P.Val"] < 0.05 && fc < -1)  { down.x.push(fc); down.y.push(logp); down.text.push(ht); }
      else                                            { ns.x.push(fc);   ns.y.push(logp);   ns.text.push(ht); }
    });

    const traces = [
      { type: "scatter", mode: "markers", name: "NS",
        x: ns.x, y: ns.y, text: ns.text, hovertemplate: "%{text}<extra></extra>",
        marker: { color: "#cbd5e1", size: 3, opacity: 0.25 } },
      { type: "scatter", mode: "markers", name: "Downregulated",
        x: down.x, y: down.y, text: down.text, hovertemplate: "%{text}<extra></extra>",
        marker: { color: "#3b82f6", size: 5, opacity: 0.7 } },
      { type: "scatter", mode: "markers", name: "Upregulated",
        x: up.x, y: up.y, text: up.text, hovertemplate: "%{text}<extra></extra>",
        marker: { color: "#ef4444", size: 5, opacity: 0.7 } }
    ];

    const layout = {
      title: { text: "<b>Entorhinal Cortex</b> — JS / Plotly.js", font: { size: 13 } },
      xaxis: { title: "Log₂ Fold Change", gridcolor: "#f1f5f9", zeroline: false },
      yaxis: { title: "-Log₁₀ (adj. P)", gridcolor: "#f1f5f9" },
      shapes: [
        { type: "line", x0: -12, x1: 12, y0: -Math.log10(0.05), y1: -Math.log10(0.05), line: { color: "#94a3b8", width: 1, dash: "dash" } },
        { type: "line", x0: 1, x1: 1, y0: 0, y1: 22, line: { color: "#94a3b8", width: 1, dash: "dash" } },
        { type: "line", x0: -1, x1: -1, y0: 0, y1: 22, line: { color: "#94a3b8", width: 1, dash: "dash" } }
      ],
      plot_bgcolor: "#ffffff", paper_bgcolor: "#f8fafc",
      margin: { t: 50, r: 20, b: 50, l: 60 },
      hovermode: "closest",
      legend: { orientation: "h", y: 1.1 }
    };

    Plotly.newPlot(containerId, traces, layout, { responsive: true, displayModeBar: false });

  } else if (plotType === "pca") {
    // Simple placeholder for PCA in JS tab
    Plotly.newPlot(containerId, [{
      type: "scatter", mode: "markers",
      x: [1, 2, 3], y: [1, 2, 3],
      marker: { color: "#3b82f6" }
    }], {
      title: { text: "PCA — coming soon in JS" },
      plot_bgcolor: "#ffffff", paper_bgcolor: "#f8fafc"
    }, { responsive: true, displayModeBar: false });

  } else {
    // For heatmap and enrichment — just clear the container
    document.getElementById(containerId).innerHTML =
      `<div style="display:flex;align-items:center;justify-content:center;height:400px;color:#94a3b8;font-family:monospace;font-size:12px">
        Interactive ${plotType} — see main dashboard above ↑
      </div>`;
  }
}