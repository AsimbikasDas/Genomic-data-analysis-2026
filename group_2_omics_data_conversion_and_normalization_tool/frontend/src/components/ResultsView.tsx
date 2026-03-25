import React, { useState } from "react";
import { CheckCircle, Download, RotateCcw, Table2, Calculator, Play } from "lucide-react";

export default function ResultsView({ data, geneIdCol, onReset }: any) {
  const [activeTab, setActiveTab] = useState<"preview" | "conversion">("preview");

  // Deep clone data to allow adding new calculated columns client-side without disrupting parent state
  const [localData, setLocalData] = useState(() => data.map((d: any) => ({ ...d })));
  
  const [convType, setConvType] = useState<"RPKM_TO_TPM" | "TPM_TO_RPKM">("RPKM_TO_TPM");
  const [selectedCols, setSelectedCols] = useState<string[]>([]);
  
  const allColumns = localData && localData.length > 0 ? Object.keys(localData[0]) : [];
  
  // Automatically display the core identifiers and any dynamically generated TPM/RPKM columns
  const columns = allColumns.filter((c: string) => 
    c === geneIdCol || 
    c === 'gene_length_bp' || 
    c.startsWith('TPM_') || 
    c.startsWith('RPKM_')
  );

  const availableColsForConv = columns.filter((c: string) => 
    convType === "RPKM_TO_TPM" ? c.startsWith("RPKM_") : c.startsWith("TPM_")
  );

  const runConversion = () => {
    if (selectedCols.length === 0) return;
    
    const newData = [...localData];
    
    selectedCols.forEach((col: string) => {
      const isTPMtoRPKM = convType === "TPM_TO_RPKM";
      const newColName = isTPMtoRPKM ? col.replace("TPM", "RPKM_conv") : col.replace("RPKM", "TPM_conv");
      
      // RPKM to TPM can be derived statically from ratio sums natively.
      if (!isTPMtoRPKM) {
         let sum = 0;
         newData.forEach(row => { sum += (parseFloat(row[col]) || 0); });
         if (sum === 0) return; 
         
         newData.forEach((row) => {
            const val = parseFloat(row[col]) || 0;
            row[newColName] = (val / sum) * 1000000;
         });
      } else {
         // TPM backwards to RPKM is mathematically impossible without scaling against original Raw Counts.
         const rawColName = col.replace("TPM_", "");
         if (newData[0] && newData[0].hasOwnProperty(rawColName) && newData[0].hasOwnProperty("gene_length_bp")) {
            let totalReads = 0;
            newData.forEach(row => { totalReads += (parseFloat(row[rawColName]) || 0); });
            if (totalReads === 0) return;
            
            newData.forEach(row => {
               const rawCounts = parseFloat(row[rawColName]) || 0;
               const len = parseFloat(row["gene_length_bp"]);
               if (len && len > 0) {
                 // Genuine mathematical RPKM computation re-constructed from ground truth depths
                 row[newColName] = (rawCounts * 1e9) / (len * totalReads);
               } else {
                 row[newColName] = 0;
               }
            });
         } else {
            alert(`Biological Data Error:\n\nMathematically, TPM cannot be unscaled back to RPKM without the original raw read sequencing depth. TPM permanently standardizes samples to 1,000,000 reads.\n\nCould not reconstruct from raw column: '${rawColName}'.`);
            return;
         }
      }
    });

    setLocalData(newData);
    setActiveTab("preview");
    setSelectedCols([]); // Clean up bounds
  };

  const handleDownload = () => {
    if (!localData || localData.length === 0) return;
    const csv = [
      columns.map((c: string) => `"${c}"`).join(","),
      ...localData.map((row: any) => columns.map((c: string) => `"${row[c] !== null && row[c] !== undefined ? row[c] : ""}"`).join(","))
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
    <div className="glass-panel rounded-3xl p-6 md:p-8 shadow-2xl relative overflow-hidden animate-in fade-in zoom-in duration-500">
      <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-emerald-400 to-cyan-500" />
      
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-6 gap-4">
        <div>
          <h2 className="text-2xl md:text-3xl font-bold flex items-center gap-3">
            <CheckCircle className="text-emerald-400 w-7 h-7" />
            Normalisation Complete
          </h2>
          <p className="text-slate-400 mt-2 text-sm md:text-base">Successfully processed <span className="text-cyan-400 font-semibold">{localData.length}</span> genes.</p>
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

      {/* Tabs */}
      <div className="flex gap-2 border-b border-slate-700/50 mb-6">
        <button 
          onClick={() => setActiveTab("preview")}
          className={`pb-3 px-4 font-medium transition-colors border-b-2 ${activeTab === "preview" ? "border-cyan-400 text-cyan-400" : "border-transparent text-slate-400 hover:text-slate-300"}`}
        >
          <span className="flex items-center gap-2"><Table2 className="w-4 h-4" /> Data Preview</span>
        </button>
        <button 
          onClick={() => setActiveTab("conversion")}
          className={`pb-3 px-4 font-medium transition-colors border-b-2 ${activeTab === "conversion" ? "border-purple-400 text-purple-400" : "border-transparent text-slate-400 hover:text-slate-300"}`}
        >
          <span className="flex items-center gap-2"><Calculator className="w-4 h-4" /> Convert Metrics</span>
        </button>
      </div>

      {/* Dynamic Tab Panes */}
      {activeTab === "preview" ? (
        <div className="bg-slate-900/40 rounded-2xl overflow-hidden border border-slate-700/50 shadow-inner">
          <div className="overflow-x-auto overflow-y-auto max-h-[500px] custom-scrollbar">
            <table className="w-full text-sm text-left text-slate-300 whitespace-nowrap">
              <thead className="bg-slate-800/90 backdrop-blur-md sticky top-0 z-10 shadow-sm border-b border-slate-700">
                <tr>
                  {columns.map((c: string) => (
                    <th key={c} className="px-5 py-4 font-semibold text-slate-200 tracking-wide border-r border-slate-700/30 last:border-0">{c}</th>
                  ))}
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-800/50">
                {localData.slice(0, 50).map((row: any, i: number) => (
                  <tr key={i} className={`transition-colors hover:bg-slate-800/60 ${i % 2 === 0 ? "bg-slate-900/20" : "bg-transparent"}`}>
                    {columns.map((c: string, colIdx: number) => (
                      <td 
                        key={c} 
                        className={`px-5 py-3 max-w-[180px] truncate border-r border-slate-800/30 last:border-0 ${colIdx !== 0 ? 'font-mono text-slate-300' : 'font-medium text-slate-200'}`} 
                        title={String(row[c])}
                      >
                        {typeof row[c] === 'number' ? row[c].toFixed(3) : row[c]}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
            {localData.length === 0 && (
              <div className="p-12 text-center text-slate-500 flex flex-col items-center">
                <Table2 className="w-12 h-12 mb-3 text-slate-600 opacity-50" />
                <p>No renderable data found.</p>
              </div>
            )}
          </div>
          {localData.length > 50 && (
            <div className="bg-slate-800/40 p-3 text-center text-xs text-slate-400 border-t border-slate-700/50">
              Showing <span className="font-semibold text-slate-300">50</span> of <span className="font-semibold text-slate-300">{localData.length}</span> rows. Export CSV to view all.
            </div>
          )}
        </div>
      ) : (
        <div className="bg-slate-900/40 rounded-2xl border border-slate-700/50 p-8 shadow-inner animate-in fade-in slide-in-from-bottom-2 duration-300">
          <div className="max-w-2xl">
            <h3 className="text-xl font-bold text-slate-200 mb-6 flex items-center gap-3">
               <Calculator className="text-purple-400 w-6 h-6" /> Internal Conversion Engine
            </h3>
            
            <div className="space-y-6">
              <div>
                 <label className="block text-sm font-semibold text-slate-400 mb-3 uppercase tracking-wide">Target Path</label>
                 <select 
                   className="w-full bg-slate-900 border-2 border-slate-700 rounded-xl p-3.5 text-slate-200 font-medium focus:outline-none focus:border-purple-500 transition-all cursor-pointer"
                   value={convType}
                   onChange={(e) => {
                     setConvType(e.target.value as any);
                     setSelectedCols([]);
                   }}
                 >
                   <option value="RPKM_TO_TPM">Convert existing RPKM ➔ TPM</option>
                   <option value="TPM_TO_RPKM">Convert existing TPM ➔ RPKM</option>
                 </select>
              </div>

              <div>
                 <label className="block text-sm font-semibold text-slate-400 mb-3 uppercase tracking-wide">Columns to Convert</label>
                 <div className="bg-slate-800/50 border border-slate-700 rounded-xl p-4 max-h-60 overflow-y-auto custom-scrollbar">
                    {availableColsForConv.length === 0 ? (
                       <p className="text-amber-500/80 text-sm p-2 flex items-center gap-2"><Table2 className="w-4 h-4"/> No valid columns available in the table for this path.</p>
                    ) : (
                       availableColsForConv.map((col: string) => (
                         <label key={col} className="flex items-center p-3 hover:bg-slate-700/40 rounded-lg cursor-pointer transition-colors border border-transparent hover:border-slate-600">
                           <input 
                             type="checkbox" 
                             className="form-checkbox h-5 w-5 text-purple-500 rounded bg-slate-900 border-slate-600 focus:ring-0 focus:ring-offset-0"
                             checked={selectedCols.includes(col)}
                             onChange={(e) => {
                               if (e.target.checked) setSelectedCols([...selectedCols, col]);
                               else setSelectedCols(selectedCols.filter(c => c !== col));
                             }}
                           />
                           <span className="ml-3 font-mono text-slate-300 text-sm">{col}</span>
                         </label>
                       ))
                    )}
                 </div>
              </div>

              <button 
                onClick={runConversion}
                disabled={selectedCols.length === 0}
                className={`py-3.5 px-6 rounded-xl flex items-center justify-center font-bold text-lg transition-all w-full md:w-auto ${
                  selectedCols.length === 0 
                  ? "bg-slate-800 text-slate-500 border border-slate-700 cursor-not-allowed" 
                  : "bg-gradient-to-r from-purple-600 to-indigo-600 hover:from-purple-500 hover:to-indigo-500 text-white shadow-lg shadow-purple-900/30 border border-purple-400/30"
                }`}
              >
                <Play className="w-5 h-5 mr-3 fill-current" /> Execute Mathematics
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
