"use client";
import React, { useState, useRef } from "react";
import { motion } from "framer-motion";
import { UploadCloud, FileType, CheckCircle2, ChevronRight, Settings2 } from "lucide-react";

export default function UploadArea({ 
  onUpload, 
  preview, 
  onNormalize, 
  isProcessing,
  selectedGeneCol,
  setSelectedGeneCol
}: any) {
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [isDragging, setIsDragging] = useState(false);
  
  const [doTpm, setDoTpm] = useState(true);
  const [doRpkm, setDoRpkm] = useState(true);

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = () => setIsDragging(false);

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      onUpload(e.dataTransfer.files[0]);
    }
  };

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      onUpload(e.target.files[0]);
    }
  };

  return (
    <div className="w-full flex flex-col md:flex-row gap-6 lg:gap-8">
      {/* Upload Zone */}
      <div className="flex-1 glass-panel rounded-[2rem] p-8 relative overflow-hidden group shadow-xl">
        <div className="absolute inset-0 bg-gradient-to-br from-cyan-500/5 to-purple-500/5 opacity-0 group-hover:opacity-100 transition-opacity duration-700 pointer-events-none" />
        
        {!preview ? (
          <div 
            className={`flex flex-col items-center justify-center h-[320px] border-2 border-dashed rounded-3xl transition-all duration-300 ease-out cursor-pointer ${isDragging ? "border-cyan-400 bg-cyan-400/10 scale-[0.98]" : "border-slate-600 hover:border-cyan-500/50 hover:bg-slate-800/30"}`}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
            onClick={() => fileInputRef.current?.click()}
          >
            <input 
              type="file" 
              accept=".csv" 
              className="hidden" 
              ref={fileInputRef} 
              onChange={handleFileChange} 
            />
            <motion.div 
              animate={{ y: [0, -8, 0] }} 
              transition={{ repeat: Infinity, duration: 2.5, ease: "easeInOut" }}
              className="bg-slate-800/50 p-5 rounded-full mb-5 shadow-inner border border-slate-700"
            >
              <UploadCloud className={`w-12 h-12 transition-colors duration-300 ${isDragging ? "text-cyan-400" : "text-slate-400 group-hover:text-cyan-300"}`} />
            </motion.div>
            <h3 className="text-2xl font-bold mb-2 text-slate-200">Drag & Drop your CSV</h3>
            <p className="text-slate-500 font-medium">or click here to browse your computer</p>
          </div>
        ) : (
          <div className="h-[320px] flex flex-col justify-center animate-in fade-in zoom-in duration-500">
            <div className="flex items-center gap-5 mb-6 bg-slate-800/40 p-4 rounded-2xl border border-slate-700/50">
              <div className="p-4 bg-emerald-500/10 rounded-xl text-emerald-400 shadow-inner">
                <FileType className="w-8 h-8" />
              </div>
              <div>
                <h3 className="text-xl font-bold text-slate-100">Dataset Loaded Ready</h3>
                <p className="text-slate-400 flex items-center gap-2 text-sm mt-1">
                  <CheckCircle2 className="w-4 h-4 text-emerald-500" /> Structure parsed successfully
                </p>
              </div>
            </div>
            
            <div className="bg-slate-900/50 rounded-2xl p-5 border border-slate-700/50 overflow-x-auto custom-scrollbar">
              <p className="text-xs text-slate-500 uppercase tracking-widest mb-3 font-bold">Data Preview (First 3 rows)</p>
              <table className="w-full text-sm text-left text-slate-300 whitespace-nowrap">
                <thead>
                  <tr>
                    {preview.columns.slice(0, 6).map((c: string) => (
                      <th key={c} className="px-4 py-2 font-medium text-slate-400 border-b border-slate-700/50">{c}</th>
                    ))}
                    {preview.columns.length > 6 && <th className="px-4 py-2 font-medium text-slate-400 border-b border-slate-700/50">...</th>}
                  </tr>
                </thead>
                <tbody>
                  {preview.preview.slice(0, 3).map((row: any, i: number) => (
                    <tr key={i} className="hover:bg-slate-800/40 transition-colors">
                      {preview.columns.slice(0, 6).map((c: string) => (
                        <td key={c} className="px-4 py-3 border-b border-slate-800/50 truncate max-w-[140px]">{row[c]}</td>
                      ))}
                      {preview.columns.length > 6 && <td className="px-4 py-3 border-b border-slate-800/50 text-slate-500">...</td>}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>

      {/* Configuration Zone */}
      <div className={`md:w-96 glass-panel rounded-[2rem] p-8 shadow-xl transition-all duration-700 ${!preview ? "opacity-30 pointer-events-none grayscale translate-x-4" : "opacity-100 translate-x-0"}`}>
        <div className="flex items-center gap-3 mb-8">
          <div className="p-2 bg-purple-500/10 rounded-lg">
            <Settings2 className="w-6 h-6 text-purple-400" />
          </div>
          <h3 className="text-2xl font-bold text-slate-100">Parameters</h3>
        </div>

        <div className="space-y-7">
          <div>
            <label className="block text-sm font-semibold text-slate-400 mb-3 uppercase tracking-wide">Gene ID Target</label>
            <div className="relative">
              <select 
                className="w-full bg-slate-900/80 border-2 border-slate-700 rounded-xl p-3.5 text-slate-200 font-medium focus:outline-none focus:border-cyan-500 transition-all appearance-none cursor-pointer"
                value={selectedGeneCol}
                onChange={(e) => setSelectedGeneCol(e.target.value)}
              >
                <option value="" disabled>Select index column...</option>
                {preview?.columns.map((c: string) => (
                  <option key={c} value={c}>{c}</option>
                ))}
              </select>
              <ChevronRight className="absolute right-4 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-500 rotate-90 pointer-events-none" />
            </div>
            {!selectedGeneCol && preview?.gene_id_col === null && (
              <p className="text-amber-500 flex items-center gap-2 text-xs mt-3 font-medium bg-amber-500/10 p-2 rounded-lg">
                <span className="animate-pulse">⚠️</span> Could not auto-detect Ensembl column.
              </p>
            )}
          </div>

          <div>
            <label className="block text-sm font-semibold text-slate-400 mb-3 uppercase tracking-wide">Execution Methods</label>
            <div className="space-y-3">
              <label className={`flex items-center p-4 border-2 rounded-xl cursor-pointer transition-all duration-300 ${doTpm ? "border-cyan-500/50 bg-cyan-500/5" : "border-slate-700 hover:border-slate-600 bg-slate-900/50"}`}>
                <input type="checkbox" className="form-checkbox h-5 w-5 text-cyan-500 rounded focus:ring-0 focus:ring-offset-0 bg-transparent border-slate-500" checked={doTpm} onChange={(e) => setDoTpm(e.target.checked)} />
                <span className={`ml-3 font-semibold transition-colors ${doTpm ? "text-cyan-300" : "text-slate-400"}`}>TPM Protocol</span>
              </label>
              <label className={`flex items-center p-4 border-2 rounded-xl cursor-pointer transition-all duration-300 ${doRpkm ? "border-purple-500/50 bg-purple-500/5" : "border-slate-700 hover:border-slate-600 bg-slate-900/50"}`}>
                <input type="checkbox" className="form-checkbox h-5 w-5 text-purple-500 rounded focus:ring-0 focus:ring-offset-0 bg-transparent border-slate-500" checked={doRpkm} onChange={(e) => setDoRpkm(e.target.checked)} />
                <span className={`ml-3 font-semibold transition-colors ${doRpkm ? "text-purple-300" : "text-slate-400"}`}>RPKM Protocol</span>
              </label>
            </div>
          </div>

          <motion.button 
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.97 }}
            disabled={!selectedGeneCol || (!doTpm && !doRpkm) || isProcessing}
            onClick={() => onNormalize({ gene_id_col: selectedGeneCol, do_tpm: doTpm, do_rpkm: doRpkm })}
            className={`w-full py-4.5 rounded-xl flex items-center justify-center font-bold text-lg transition-all duration-300 ${
              (!selectedGeneCol || (!doTpm && !doRpkm)) 
                ? "bg-slate-800 text-slate-500 border border-slate-700 cursor-not-allowed" 
                : "bg-gradient-to-r from-cyan-600 to-blue-600 hover:from-cyan-500 hover:to-blue-500 text-white shadow-xl shadow-cyan-900/30 border border-cyan-400/30"
            }`}
          >
            {isProcessing ? (
              <>
                <div className="w-5 h-5 border-2 border-white/30 border-t-white rounded-full animate-spin mr-3" />
                Processing Dataset...
              </>
            ) : (
              <>
                Ignite Pipeline
                <ChevronRight className="w-5 h-5 ml-2" />
              </>
            )}
          </motion.button>
        </div>
      </div>
    </div>
  );
}
