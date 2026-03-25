import React, { useState } from "react";
import { CheckCircle, Download, RotateCcw, Table2 } from "lucide-react";

export default function ResultsView({ data, geneIdCol, onReset }: any) {
  const [activeTab, setActiveTab] = useState<"preview" | "conversion">("preview");

  const columns = data && data.length > 0 ? Object.keys(data[0]) : [];
  
  const handleDownload = () => {
    if (!data || data.length === 0) return;
    const csv = [
      columns.map(c => `"${c}"`).join(","),
      ...data.map((row: any) => columns.map(c => `"${row[c] !== null && row[c] !== undefined ? row[c] : ""}"`).join(","))
    ].join("\n");

    const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.setAttribute("download", "omicsforge_normalized.csv");
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <div className="glass-panel rounded-3xl p-6 md:p-8 shadow-2xl relative overflow-hidden">
      <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-emerald-400 to-cyan-500" />
      
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-6 gap-4">
        <div>
          <h2 className="text-2xl md:text-3xl font-bold flex items-center gap-3">
            <CheckCircle className="text-emerald-400 w-7 h-7" />
            Normalisation Complete
          </h2>
          <p className="text-slate-400 mt-2 text-sm md:text-base">Processed <span className="text-cyan-400 font-semibold">{data.length}</span> valid genes successfully.</p>
        </div>

        <div className="flex gap-3">
          <button 
            onClick={onReset}
            className="px-4 py-2.5 rounded-xl border border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white transition-all flex items-center gap-2 group"
          >
            <RotateCcw className="w-4 h-4 group-hover:-rotate-180 transition-transform duration-500" /> Start Over
          </button>
          <button 
            onClick={handleDownload}
            className="px-5 py-2.5 rounded-xl bg-gradient-to-r from-cyan-600 to-blue-600 hover:from-cyan-500 hover:to-blue-500 text-white transition-all flex items-center gap-2 shadow-lg shadow-cyan-900/30 hover:scale-105"
          >
            <Download className="w-4 h-4" /> Export CSV
          </button>
        </div>
      </div>

      <div className="flex gap-2 border-b border-slate-700/50 mb-6">
        <button 
          onClick={() => setActiveTab("preview")}
          className={`pb-3 px-4 font-medium transition-colors border-b-2 ${activeTab === "preview" ? "border-cyan-400 text-cyan-400" : "border-transparent text-slate-400 hover:text-slate-300"}`}
        >
          <span className="flex items-center gap-2"><Table2 className="w-4 h-4" /> Data Preview</span>
        </button>
      </div>

      {/* Upgraded SOTA Table Wrapper */}
      <div className="bg-slate-900/40 rounded-2xl overflow-hidden border border-slate-700/50 shadow-inner">
        <div className="overflow-x-auto overflow-y-auto max-h-[500px]">
          <table className="w-full text-sm text-left text-slate-300 whitespace-nowrap">
            <thead className="bg-slate-800/90 backdrop-blur-md sticky top-0 z-10 shadow-sm border-b border-slate-700">
              <tr>
                {columns.slice(0, 10).map((c) => (
                  <th key={c} className="px-5 py-4 font-semibold text-slate-200 tracking-wide border-r border-slate-700/30 last:border-0">{c}</th>
                ))}
                {columns.length > 10 && <th className="px-5 py-4 font-semibold text-slate-400">...</th>}
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-800/50">
              {data.slice(0, 50).map((row: any, i: number) => (
                <tr key={i} className={`transition-colors hover:bg-slate-800/60 ${i % 2 === 0 ? "bg-slate-900/20" : "bg-transparent"}`}>
                  {columns.slice(0, 10).map((c, colIdx) => (
                    <td 
                      key={c} 
                      className={`px-5 py-3 max-w-[180px] truncate border-r border-slate-800/30 last:border-0 ${colIdx !== 0 ? 'font-mono text-slate-300' : 'font-medium text-slate-200'}`} 
                      title={String(row[c])}
                    >
                      {typeof row[c] === 'number' ? row[c].toFixed(3) : row[c]}
                    </td>
                  ))}
                  {columns.length > 10 && <td className="px-5 py-3 text-slate-500">...</td>}
                </tr>
              ))}
            </tbody>
          </table>
          {data.length === 0 && (
            <div className="p-12 text-center text-slate-500 flex flex-col items-center">
              <Table2 className="w-12 h-12 mb-3 text-slate-600 opacity-50" />
              <p>No renderable data found.</p>
            </div>
          )}
        </div>
        {data.length > 50 && (
          <div className="bg-slate-800/40 p-3 text-center text-xs text-slate-400 border-t border-slate-700/50">
            Showing <span className="font-semibold text-slate-300">50</span> of <span className="font-semibold text-slate-300">{data.length}</span> rows. Please export the CSV to view the entire dataset.
          </div>
        )}
      </div>
    </div>
  );
}
