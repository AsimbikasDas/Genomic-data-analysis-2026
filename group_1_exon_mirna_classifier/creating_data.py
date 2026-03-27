import random
import pandas as pd
from Bio import SeqIO

def approach_1(fasta_file, gff_file, fixed_length=22):
    print("Loading FASTA file...")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    mirna_coords = []
    exon_coords = []
    
    print("Parsing GFF file...")
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
                
            chrom = parts[0]
            feature_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            if feature_type == 'mirna':

                mirna_len = end - start + 1
                if mirna_len > fixed_length:
                    mid = (start + end) // 2
                    new_start = mid - (fixed_length // 2)
                    new_end = new_start + fixed_length - 1
                elif mirna_len < fixed_length:
                    new_start = start
                    new_end = end+fixed_length-mirna_len
                else:
                    new_start = start
                    new_end = end
                mirna_coords.append((chrom, new_start, new_end, strand))
            elif feature_type == 'exon':
                if (end - start + 1) >= fixed_length:
                    exon_coords.append((chrom, start, end, strand))

    def extract_seqs(coords, label, enforce_length=False):
        extracted = []
        for chrom, start, end, strand in coords:
            if chrom in fasta_dict:

                seq_obj = fasta_dict[chrom].seq[start-1:end]
                
                if strand == '-':
                    seq_obj = seq_obj.reverse_complement()
                
                seq_str = str(seq_obj).upper()
                
                if enforce_length and len(seq_str) >= fixed_length:
                    max_start = len(seq_str) - fixed_length
                    slice_start = random.randint(0, max_start)
                    seq_str = seq_str[slice_start:slice_start + fixed_length]
                    
                if 'N' not in seq_str:
                    extracted.append({'sequence': seq_str, 'label': label})
                    
        return extracted

    print(f"Extracting sequences for {len(mirna_coords)} miRNAs...")
    mirna_data = extract_seqs(mirna_coords, 'miRNA', enforce_length=False)
    
    
    print(f"Extracting sequences for {len(exon_coords)} exons...")
    exon_data = extract_seqs(exon_coords, 'exon', enforce_length=True)
    
    num_mirnas = len(mirna_data)
    target_exons = num_mirnas * 4
    
    if len(exon_data) >= target_exons:
        exon_data = random.sample(exon_data, target_exons)
        print(f"Sampled {target_exons} exons to match the 4:1 ratio target.")
    else:
        print(f"Warning: Only found {len(exon_data)} valid exons.")
        
    combined_data = mirna_data + exon_data
    df = pd.DataFrame(combined_data)
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    print(f"Final dataset created with {len(df)} total sequences.")
    print(df['label'].value_counts())
    
    return df

def approach_2(fasta_file, gff_file, target_length=100):
    print("Loading FASTA file...")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    mirna_coords = []
    exon_coords = []
    
    print("Parsing GFF file...")
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
                
            chrom = parts[0]
            feature_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            if chrom not in fasta_dict:
                continue
                
            chrom_len = len(fasta_dict[chrom].seq)
            
            if feature_type == 'mirna':
                mid = (start + end) // 2
                new_start = mid - (target_length // 2)
                new_end = new_start + target_length - 1
                
                if new_start >= 1 and new_end <= chrom_len:
                    mirna_coords.append((chrom, new_start, new_end, strand))
                    
            elif feature_type == 'exon':
                if (end - start + 1) >= target_length:
                    max_allowed_start = end - target_length + 1
                    new_start = random.randint(start, max_allowed_start)
                    new_end = new_start + target_length - 1
                    
                    exon_coords.append((chrom, new_start, new_end, strand))

    def extract_seqs(coords, label):
        extracted = []
        for chrom, start, end, strand in coords:
            seq_obj = fasta_dict[chrom].seq[start-1:end]
            
            if strand == '-':
                seq_obj = seq_obj.reverse_complement()
            
            seq_str = str(seq_obj).upper()
            
            if len(seq_str) == target_length and 'N' not in seq_str:
                extracted.append({'sequence': seq_str, 'label': label})
                
        return extracted

    print(f"Extracting {target_length}nt sequences for {len(mirna_coords)} miRNAs...")
    mirna_data = extract_seqs(mirna_coords, 'miRNA')
    
    print(f"Extracting {target_length}nt sequences for {len(exon_coords)} exons...")
    exon_data = extract_seqs(exon_coords, 'exon')
    
    num_mirnas = len(mirna_data)
    target_exons = num_mirnas * 4
    
    if len(exon_data) >= target_exons:
        exon_data = random.sample(exon_data, target_exons)
        print(f"Sampled {target_exons} exons to match the 4:1 ratio target.")
    else:
        print(f"Warning: Only found {len(exon_data)} valid exons (needed {target_exons}).")
        
    combined_data = mirna_data + exon_data
    df = pd.DataFrame(combined_data)
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    print(f"Final dataset created with {len(df)} total sequences.")
    print(df['label'].value_counts())
    
    return df
