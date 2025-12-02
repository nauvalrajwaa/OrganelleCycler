# /scripts/assess_assemblies.py

import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

# --- CONFIG & INPUTS ---
if 'snakemake' in globals():
    # Input files
    in_flye = snakemake.input.flye
    in_raven = snakemake.input.raven
    in_canu = snakemake.input.canu  # AKTIF: Input dari Canu
    
    # Reference file (Dinamis: Plastome atau Mito tergantung Snakefile)
    ref_file = snakemake.input.ref
    
    # Output files
    out_fasta = snakemake.output.best
    out_report = snakemake.output.report
    
    # Parameters (Diambil dari 'params' Snakefile agar dinamis)
    min_len = snakemake.params.min_len
    max_len = snakemake.params.max_len
else:
    sys.exit("Script ini harus dijalankan via Snakemake.")

# --- FUNGSI BANTUAN ---

def get_blast_metrics(query_fasta, subject_fasta):
    """
    Menjalankan BLASTN untuk mendapatkan Coverage dan Identity.
    Output fmt 6: qcovs (Query Coverage), pident (Percentage Identity)
    """
    out_temp = "temp_blast.tsv"
    cmd = [
        "blastn", 
        "-query", query_fasta, 
        "-subject", subject_fasta, 
        "-outfmt", "6 qcovs pident", 
        "-max_target_seqs", "1"
    ]
    
    try:
        subprocess.run(cmd, stdout=open(out_temp, "w"), stderr=subprocess.DEVNULL, check=True)
        
        # Jika hasil kosong (tidak ada hit)
        if os.path.getsize(out_temp) == 0:
            return 0, 0
            
        df = pd.read_csv(out_temp, sep="\t", header=None, names=["qcovs", "pident"])
        
        # Ambil nilai maksimum (single best hit strategy)
        cov = df["qcovs"].max()
        ident = df["pident"].mean() # Rata-rata identity dari alignment blocks
        
        return cov, ident
        
    except Exception as e:
        return 0, 0
    finally:
        if os.path.exists(out_temp): os.remove(out_temp)

def check_overlap_circularity(seq_obj, overlap_len=1000):
    """
    Cek manual: Apakah ujung depan dan ujung belakang sama?
    (Raven/Canu kadang menghasilkan output linear tapi sebenarnya ujungnya overlap)
    """
    seq = str(seq_obj.seq)
    if len(seq) < overlap_len * 2: return False
    
    start = seq[:overlap_len]
    end = seq[-overlap_len:]
    
    # Simple string matching (exact match)
    if start == end: 
        return True
    return False

# --- LOGIKA UTAMA ---

def evaluate_assemblies():
    candidates = []
    
    # List tools yang akan dinilai
    # Format: (Nama_Tool, File_Path)
    tools_to_check = [
        ("Flye", in_flye),
        ("Raven", in_raven),
        ("Canu", in_canu) # AKTIF: Canu dimasukkan ke penilaian
    ]
    
    print(f"--- MEMULAI PENILAIAN KANDIDAT (Target Size: {min_len}-{max_len} bp) ---")
    
    for tool_name, fasta_path in tools_to_check:
        if not os.path.exists(fasta_path):
            continue
            
        for record in SeqIO.parse(fasta_path, "fasta"):
            length = len(record.seq)
            
            # 1. HARD FILTER: Ukuran (Dinamis sesuai Target)
            if not (min_len <= length <= max_len):
                # Skip contig yang terlalu kecil/besar
                continue
            
            # 2. CEK SIRKULARITAS
            is_circular = False
            # Flye memberi label 'circular=true'
            if "circular=true" in record.description:
                is_circular = True
            # Cek manual overlap (untuk Raven/Canu)
            elif check_overlap_circularity(record):
                is_circular = True
                
            # 3. CEK BLAST (Kualitas vs Reference yang sesuai)
            # Tulis contig tunggal ke temp file untuk di-blast
            temp_ctg = f"temp_{tool_name}.fasta"
            SeqIO.write(record, temp_ctg, "fasta")
            
            cov, ident = get_blast_metrics(temp_ctg, ref_file)
            
            if os.path.exists(temp_ctg): os.remove(temp_ctg)
            
            # 4. HITUNG SKOR
            # Rumus Skor:
            # - Base Score: Coverage + Identity (Max 200)
            # - Circular Bonus: +500 (Sangat berharga)
            # - Penalty: Jika identity < 80% (mungkin salah spesies/junk), skor dikurangi
            
            score = cov + ident
            if is_circular:
                score += 500
            
            if ident < 80:
                score -= 1000 # Hukuman berat untuk sequence ngawur
            
            candidates.append({
                "tool": tool_name,
                "id": record.id,
                "len": length,
                "circ": is_circular,
                "cov": cov,
                "ident": ident,
                "score": score,
                "record": record
            })

    # --- PENENTUAN PEMENANG ---
    
    # Urutkan berdasarkan Score Tertinggi
    candidates.sort(key=lambda x: x["score"], reverse=True)
    
    report_lines = ["Tool | ContigID | Length | Circ | Cov% | Ident% | SCORE"]
    report_lines.append("-" * 70)
    
    best_record = None
    
    for c in candidates:
        line = f"{c['tool']} | {c['id']} | {c['len']} | {c['circ']} | {c['cov']:.1f} | {c['ident']:.1f} | {c['score']:.1f}"
        report_lines.append(line)
    
    if candidates:
        winner = candidates[0]
        # Validasi akhir: Skor harus positif (berarti bukan junk)
        if winner["score"] > 0:
            best_record = winner["record"]
            # Rename header biar rapi
            best_record.description = f"source={winner['tool']} length={winner['len']} circular={winner['circ']} coverage={winner['cov']:.2f} identity={winner['ident']:.2f}"
            best_record.id = "Target_Consensus"
            report_lines.append(f"\nWINNER: {winner['tool']} (Score: {winner['score']})")
        else:
            report_lines.append("\nNO WINNER: Top candidate score is too low (Bad Identity).")
    else:
        report_lines.append("\nNO CANDIDATES PASSED LENGTH FILTER.")

    # Tulis Output
    with open(out_report, "w") as f:
        f.write("\n".join(report_lines))
        
    if best_record:
        SeqIO.write(best_record, out_fasta, "fasta")
    else:
        # Jika gagal total, buat file kosong agar snakemake tidak error,
        # tapi user harus cek report.
        open(out_fasta, 'a').close()

if __name__ == "__main__":
    evaluate_assemblies()