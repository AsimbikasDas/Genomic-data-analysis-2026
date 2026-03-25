"use client";

import { useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Dna } from "lucide-react";
import UploadArea from "../components/UploadArea";
import ResultsView from "../components/ResultsView";

export default function Home() {
  const [file, setFile] = useState<File | null>(null);
  const [dataPreview, setDataPreview] = useState<any>(null);
  const [isProcessing, setIsProcessing] = useState(false);
  const [normalizedData, setNormalizedData] = useState<any>(null);
  const [geneIdCol, setGeneIdCol] = useState<string>("");

  let API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";
  if (API_BASE_URL && !API_BASE_URL.startsWith("http")) {
    API_BASE_URL = "https://" + API_BASE_URL; // Defends against missed https:// in Vercel config
  }

  const handleUpload = async (uploadedFile: File) => {
    setFile(uploadedFile);
    const formData = new FormData();
    formData.append("file", uploadedFile);
    try {
      const res = await fetch(`${API_BASE_URL}/api/preview`, {
        method: "POST",
        body: formData,
      });
      const data = await res.json();
      setDataPreview(data);
      if (data.gene_id_col) {
        setGeneIdCol(data.gene_id_col);
      }
    } catch (err) {
      console.error(err);
      alert("Failed to read file.");
    }
  };

  const handleNormalize = async (config: { gene_id_col: string, do_tpm: boolean, do_rpkm: boolean }) => {
    if (!file) return;
    setIsProcessing(true);
    const formData = new FormData();
    formData.append("file", file);
    formData.append("gene_id_col", config.gene_id_col);
    formData.append("do_tpm", config.do_tpm.toString());
    formData.append("do_rpkm", config.do_rpkm.toString());

    try {
      const res = await fetch(`${API_BASE_URL}/api/normalize`, {
        method: "POST",
        body: formData,
      });
      
      if (!res.ok) {
        const text = await res.text();
        throw new Error(`Server Error (${res.status}): ${text}`);
      }
      
      const data = await res.json();
      
      if (!data || !data.result) {
        throw new Error("Invalid payload returned from the backend pipeline.");
      }
      
      if (data.result.length === 0) {
        alert("Warning: Normalisation resulted in an empty dataset! This usually means NONE of your Gene IDs mapped successfully against the Exonic Length Database.");
      }

      setNormalizedData(data.result);
      setGeneIdCol(config.gene_id_col);
    } catch (err: any) {
      console.error(err);
      alert(`Normalization Process Failed:\n\n${err.message}`);
    } finally {
      setIsProcessing(false);
    }
  };

  const resetAll = () => {
    setNormalizedData(null);
    setFile(null);
    setDataPreview(null);
  };

  return (
    <main className="min-h-screen p-8 text-slate-200 font-sans">
      <div className="dna-bg" />

      <div className="max-w-6xl mx-auto pt-10">
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          className="text-center mb-12"
        >
          <div className="inline-flex items-center justify-center p-4 rounded-2xl glass-panel mb-6 border-cyan-500/30">
            <Dna className="w-10 h-10 text-cyan-400" />
          </div>
          <h1 className="text-6xl font-extrabold bg-clip-text text-transparent bg-gradient-to-r from-cyan-400 via-blue-400 to-purple-500 tracking-tight pb-2">
            OmicsForge
          </h1>
          <p className="mt-4 text-slate-400 text-xl max-w-2xl mx-auto font-light">
            Ultra-fast RNA normalisation powered by <span className="text-cyan-300 font-medium">Pandas & NumPy</span>.
            Instant TPM and RPKM calculations in a beautiful interface.
          </p>
        </motion.div>

        <AnimatePresence mode="wait">
          {!normalizedData ? (
            <motion.div
              key="upload"
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              exit={{ opacity: 0, scale: 0.95, filter: "blur(10px)" }}
              transition={{ duration: 0.4, type: "spring", bounce: 0.3 }}
            >
              <UploadArea
                onUpload={handleUpload}
                preview={dataPreview}
                onNormalize={handleNormalize}
                isProcessing={isProcessing}
                selectedGeneCol={geneIdCol}
                setSelectedGeneCol={setGeneIdCol}
              />
            </motion.div>
          ) : (
            <motion.div
              key="results"
              initial={{ opacity: 0, y: 40 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.5, delay: 0.1, type: "spring" }}
            >
              <ResultsView
                data={normalizedData}
                geneIdCol={geneIdCol}
                onReset={resetAll}
              />
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </main>
  );
}
