# /scripts/assess_assemblies.py

import sys
import os
import pandas as pd
from Bio import SeqIO
import subprocess

# --- CONFIG & INPUTS ---
if 'snakemake' in globals():
    in_flye = snakemake.input.flye
    in_raven = snakemake.input.raven
    in_canu = snakemake.input.canu
    ref_file = snakemake.input.ref
    out_fasta = snakemake.output.best
    out_report = snakemake.output.report
    min_len = snakemake.params.min_len
    max_len = snakemake.params.max_len
else:
    sys.exit("Script ini harus dijalankan via Snakemake.")

# --- FUNGSI BANTUAN ---

def get_blast_metrics(query_fasta, subject_fasta):
    """
    BLAST semua contig sekaligus ke referensi.
    Kita cari TOTAL Coverage (bukan max) dan Rata-rata Identity.
    """
    out_temp = "temp_blast.tsv"
    # qcovs = Query Coverage per subject
    cmd = ["blastn", "-query", query_fasta, "-subject", subject_fasta, 
           "-outfmt", "6 qcovs pident length", "-max_target_seqs", "5"]
    
    try:
        subprocess.run(cmd, stdout=open(out_temp, "w"), stderr=subprocess.DEVNULL, check=True)
        if os.path.getsize(out_temp) == 0: return 0, 0
        
        df = pd.read_csv(out_temp, sep="\t", header=None, names=["qcovs", "pident", "len"])
        
        # 1. Total Coverage: Asumsi sederhana, ambil max coverage per alignment blok
        # (Ini pendekatan kasar tapi cepat untuk ranking)
        cov = df["qcovs"].max() # Maksimal referensi yang tercover oleh satu contig terbaik
        if cov > 100: cov = 100
        
        # 2. Average Identity: Rata-rata kemiripan
        ident = df["pident"].mean()
        
        return cov, ident
    except:
        return 0, 0
    finally:
        if os.path.exists(out_temp): os.remove(out_temp)

def is_circular(record):
    if "circular=true" in record.description: return True
    # Cek overlap manual sederhana
    s = str(record.seq)
    if len(s) > 2000 and s[:1000] == s[-1000:]: return True
    return False

# --- LOGIKA UTAMA ---

def evaluate_assemblies():
    # Target Ukuran Biologis (Titik Tengah Range)
    target_min = float(min_len)
    target_max = float(max_len)
    target_ideal = (target_min + target_max) / 2
    
    tools = [("Flye", in_flye), ("Raven", in_raven), ("Canu", in_canu)]
    scoreboard = []

    print(f"--- RANKING MODE (Target: {int(target_min)}-{int(target_max)} bp) ---")

    for tool_name, fasta_path in tools:
        # Jika file tidak ada (misal Canu gagal total), beri nilai 0
        if not os.path.exists(fasta_path):
            scoreboard.append({"tool": tool_name, "score": -999, "note": "FAILED/MISSING"})
            continue
            
        # 1. BACA & FILTER NOISE
        valid_records = []
        total_len = 0
        has_circle = False
        
        for rec in SeqIO.parse(fasta_path, "fasta"):
            # Kita tetap buang sampah mikro (<1kb) agar tidak merusak perhitungan rata-rata
            if len(rec.seq) < 1000: continue
            
            valid_records.append(rec)
            total_len += len(rec.seq)
            if is_circular(rec): has_circle = True
            
        # Jika tidak ada contig valid
        if total_len == 0:
            scoreboard.append({"tool": tool_name, "score": -999, "note": "EMPTY"})
            continue

        # 2. HITUNG SKOR
        # A. BLAST Score (Quality)
        cov, ident = get_blast_metrics(fasta_path, ref_file)
        # Bobot Coverage x2 (Penting LENGKAP), Identity x1 (Penting BENAR)
        blast_score = (cov * 2) + ident
        
        # B. Size Score (Proximity to Target)
        # Semakin dekat ke target_ideal, semakin kecil penaltinya.
        dist = abs(total_len - target_ideal)
        # Rumus: Penalti persentase selisih * 50
        size_penalty = (dist / target_ideal) * 50
        
        # C. Structure Bonus
        circ_bonus = 100 if has_circle else 0
        
        # TOTAL SCORE
        final_score = blast_score - size_penalty + circ_bonus
        
        # Simpan Data
        scoreboard.append({
            "tool": tool_name,
            "score": final_score,
            "total_len": total_len,
            "cov": cov,
            "ident": ident,
            "circ": has_circle,
            "records": valid_records,
            "note": "OK"
        })
        
        print(f"> {tool_name}: Len={total_len} Cov={cov:.1f}% Circ={has_circle} -> SCORE={final_score:.1f}")

    # --- PENENTUAN RANKING ---
    # Sort dari skor tertinggi ke terendah
    scoreboard.sort(key=lambda x: x["score"], reverse=True)
    
    # --- OUTPUT REPORT ---
    lines = ["RANK | TOOL | TotalLen | Circ | Cov% | Ident% | SCORE"]
    lines.append("-" * 65)
    
    # Loop untuk menulis ranking 1-2-3
    for i, data in enumerate(scoreboard):
        rank = i + 1
        if data["score"] == -999:
            lines.append(f"FAIL | {data['tool']} | - | - | - | - | FAILED")
        else:
            lines.append(f"  {rank}  | {data['tool']} | {data['total_len']} | {data['circ']} | {data['cov']:.1f} | {data['ident']:.1f} | {data['score']:.1f}")

    # --- OUTPUT FASTA JUARA 1 ---
    # Kita ambil Juara 1 (Index 0), siapapun dia, sejelek apapun dia.
    winner = scoreboard[0]
    
    if winner["score"] != -999:
        lines.append(f"\nWINNER: {winner['tool']} (Length: {winner['total_len']} bp)")
        
        # Tulis semua contig milik pemenang ke file output
        output_recs = []
        for i, rec in enumerate(winner["records"]):
            rec.id = f"ctg_{i+1}_{winner['tool']}_{rec.id}"
            rec.description = f"[Rank 1 Assembly] score={winner['score']:.1f}"
            output_recs.append(rec)
            
        SeqIO.write(output_recs, out_fasta, "fasta")
    else:
        lines.append("\nNO WINNER: All tools failed.")
        open(out_fasta, 'a').close()
        
    # Simpan Report
    with open(out_report, "w") as f:
        f.write("\n".join(lines))

if __name__ == "__main__":
    evaluate_assemblies()