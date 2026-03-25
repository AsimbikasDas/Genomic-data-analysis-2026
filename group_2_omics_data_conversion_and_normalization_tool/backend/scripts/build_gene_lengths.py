import pandas as pd
from collections import defaultdict

def build_exonic_gene_lengths(gtf_path):
    gene_exons = defaultdict(list)

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")

            # Only exon entries
            if parts[2] != "exon":
                continue

            start, end = int(parts[3]), int(parts[4])

            info = parts[8]
            try:
                gene_id = info.split('gene_id "')[1].split('"')[0]
            except:
                continue

            gene_exons[gene_id].append((start, end))

    # Merge overlapping exons and calculate length
    gene_lengths = {}

    for gene, exons in gene_exons.items():
        exons = sorted(exons)
        merged = []

        for start, end in exons:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)

        length = sum(e[1] - e[0] + 1 for e in merged)
        gene_lengths[gene] = length

    return gene_lengths


if __name__ == "__main__":
    gtf_file = "Homo_sapiens.GRCh38.115.gtf" 

    print("Building gene length database...")
    gene_lengths = build_exonic_gene_lengths(gtf_file)

    df = pd.DataFrame({
        "gene_id": list(gene_lengths.keys()),
        "gene_length_bp": list(gene_lengths.values())
    })

    df.to_csv("gene_lengths_exonic.csv", index=False)

    print(f"Done Saved {len(df)} genes")
