"use client";
import React, { useState } from "react";
import { BarChart3, Play, Download, Maximize2, X } from "lucide-react";

interface AnalysisViewProps {
  normalizedCsvBlob: Blob | null;
  geneIdCol: string;
  apiBaseUrl: string;
}

export default function AnalysisView({ normalizedCsvBlob, geneIdCol, apiBaseUrl }: AnalysisViewProps) {
  const [method, setMethod] = useState<"spearman" | "pearson">("spearman");
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [expandedImg, setExpandedImg] = useState<string | null>(null);

  const runAnalysis = async () => {
    if (!normalizedCsvBlob) return;
    setIsLoading(true);
    try {
      const formData = new FormData();
      formData.append("file", normalizedCsvBlob, "normalized_data.csv");
      formData.append("gene_id_col", geneIdCol);
      formData.append("method", method);

      const res = await fetch(`${apiBaseUrl}/api/analyze`, {
        method: "POST",
        body: formData,
      });

      if (!res.ok) {
        const errText = await res.text();
        throw new Error(errText);
      }

      const data = await res.json();
      setResults(data);
    } catch (err: any) {
      console.error("Analysis error:", err);
      alert(`Analysis Failed:\n${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  const downloadCorrelationMatrix = () => {
    if (!results?.correlation?.matrix || !results?.correlation?.columns) return;
    const cols = results.correlation.columns;
    const matrix = results.correlation.matrix;
    const header = ["", ...cols].join(",");
    const rows = matrix.map((row: number[], i: number) =>
      [cols[i], ...row.map((v: number) => v.toFixed(4))].join(",")
    );
    const csv = [header, ...rows].join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = `correlation_matrix_${method}.csv`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <div className="mt-10">
      {/* Fullscreen Overlay */}
      {expandedImg && (
        <div
          className="fixed inset-0 z-50 bg-black/80 backdrop-blur-sm flex items-center justify-center p-8 cursor-pointer"
          onClick={() => setExpandedImg(null)}
        >
          <button
            className="absolute top-6 right-6 text-white/70 hover:text-white transition-colors"
            onClick={() => setExpandedImg(null)}
          >
            <X className="w-8 h-8" />
          </button>
          <img
            src={`data:image/png;base64,${expandedImg}`}
            alt="Expanded plot"
            className="max-w-full max-h-full object-contain rounded-xl shadow-2xl"
          />
        </div>
      )}

      <div className="glass-panel rounded-3xl p-6 md:p-8 shadow-2xl relative overflow-hidden">
        <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-rose-500 to-orange-500" />

        <h2 className="text-2xl md:text-3xl font-bold flex items-center gap-3 mb-6">
          <BarChart3 className="text-rose-400 w-7 h-7" />
          Data Analysis & Visualizations
        </h2>

        {/* Controls */}
        <div className="flex flex-wrap items-end gap-4 mb-8">
          <div>
            <label className="block text-sm font-semibold text-slate-400 mb-2 uppercase tracking-wide">
              Correlation Method
            </label>
            <select
              className="bg-slate-900 border-2 border-slate-700 rounded-xl p-3 text-slate-200 font-medium focus:outline-none focus:border-rose-500 transition-all cursor-pointer min-w-[200px]"
              value={method}
              onChange={(e) => setMethod(e.target.value as any)}
            >
              <option value="spearman">spearman</option>
              <option value="pearson">pearson</option>
            </select>
          </div>

          <button
            onClick={runAnalysis}
            disabled={isLoading}
            className={`px-6 py-3 rounded-xl flex items-center gap-2 font-bold text-base transition-all ${
              isLoading
                ? "bg-slate-800 text-slate-500 cursor-not-allowed"
                : "bg-gradient-to-r from-rose-600 to-orange-600 hover:from-rose-500 hover:to-orange-500 text-white shadow-lg shadow-rose-900/30 hover:scale-105"
            }`}
          >
            {isLoading ? (
              <>
                <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
                Analyzing...
              </>
            ) : (
              <>
                <Play className="w-4 h-4 fill-current" />
                Run Analysis
              </>
            )}
          </button>
        </div>

        {/* Results */}
        {results && (
          <div className="space-y-10 animate-in fade-in duration-500">
            <hr className="border-slate-700/50" />

            {/* 1. Correlation Heatmap */}
            {results.correlation && (
              <div>
                <h3 className="text-xl font-bold text-slate-200 mb-4">1. Correlation Heatmap</h3>
                <div className="relative bg-white rounded-xl overflow-hidden inline-block">
                  <button
                    className="absolute top-2 right-2 z-10 p-1.5 bg-slate-800/70 hover:bg-slate-700 rounded-lg text-white/70 hover:text-white transition-all"
                    onClick={() => setExpandedImg(results.correlation.image)}
                    title="Expand"
                  >
                    <Maximize2 className="w-4 h-4" />
                  </button>
                  <img
                    src={`data:image/png;base64,${results.correlation.image}`}
                    alt="Correlation Heatmap"
                    className="max-w-full"
                    style={{ maxHeight: "500px" }}
                  />
                </div>
                <div className="mt-4">
                  <button
                    onClick={downloadCorrelationMatrix}
                    className="px-4 py-2.5 rounded-xl border border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white transition-all flex items-center gap-2"
                  >
                    <Download className="w-4 h-4" />
                    Download Correlation Matrix
                  </button>
                </div>
              </div>
            )}

            {/* 2. PCA */}
            {results.pca && (
              <div>
                <hr className="border-slate-700/50 mb-10" />
                <h3 className="text-xl font-bold text-slate-200 mb-4">2. Principal Component Analysis (PCA)</h3>
                <div className="relative bg-white rounded-xl overflow-hidden inline-block">
                  <button
                    className="absolute top-2 right-2 z-10 p-1.5 bg-slate-800/70 hover:bg-slate-700 rounded-lg text-white/70 hover:text-white transition-all"
                    onClick={() => setExpandedImg(results.pca)}
                    title="Expand"
                  >
                    <Maximize2 className="w-4 h-4" />
                  </button>
                  <img
                    src={`data:image/png;base64,${results.pca}`}
                    alt="PCA Plot"
                    className="max-w-full"
                    style={{ maxHeight: "500px" }}
                  />
                </div>
              </div>
            )}

            {/* 3. Distribution */}
            {results.distribution && (
              <div>
                <hr className="border-slate-700/50 mb-10" />
                <h3 className="text-xl font-bold text-slate-200 mb-4">3. Distribution Comparisons</h3>
                <div className="relative bg-white rounded-xl overflow-hidden inline-block">
                  <button
                    className="absolute top-2 right-2 z-10 p-1.5 bg-slate-800/70 hover:bg-slate-700 rounded-lg text-white/70 hover:text-white transition-all"
                    onClick={() => setExpandedImg(results.distribution)}
                    title="Expand"
                  >
                    <Maximize2 className="w-4 h-4" />
                  </button>
                  <img
                    src={`data:image/png;base64,${results.distribution}`}
                    alt="Distribution Comparison"
                    className="max-w-full"
                    style={{ maxHeight: "600px" }}
                  />
                </div>
              </div>
            )}

            {/* No results message */}
            {!results.correlation && !results.pca && !results.distribution && (
              <div className="text-center text-slate-500 py-12">
                <BarChart3 className="w-12 h-12 mx-auto mb-3 opacity-40" />
                <p>Not enough sample columns to generate analysis plots.</p>
                <p className="text-sm mt-1">Need at least 2 numeric sample columns.</p>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}
