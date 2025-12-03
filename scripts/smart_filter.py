# /scripts/smart_filter.py

import sys
import os
import subprocess
from Bio import SeqIO

# --- CONFIG & INPUTS ---
if 'snakemake' in globals():
    input_fastas = snakemake.input.fastas 
    input_gfas = snakemake.input.gfas 
    
    ref_plastome = snakemake.input.ref_p
    ref_mito = snakemake.input.ref_m
    
    output_blacklist = snakemake.output.blacklist
    output_filtered_gfas = snakemake.output.filtered_gfas 
    log_file = snakemake.output.log
    
    config = snakemake.config
    target_mode = config.get("target_assembly", "PLASTOME")
    
    if target_mode == "MITO":
        valid_min = config.get("mito_min_len", 200000)
        valid_max = config.get("mito_max_len", 2000000) 
    else:
        valid_min = config.get("plastome_min", 1000) # Diturunkan ke 1kb agar fragmen kecil terdeteksi
        valid_max = config.get("plastome_max", 500000)
else:
    sys.exit("Run via Snakemake only.")

# --- FUNGSI BANTUAN ---

def parse_gfa_circularity(gfa_file):
    circular_contigs = set()
    if os.path.exists(gfa_file):
        with open(gfa_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[0] == 'L' and parts[1] == parts[3]:
                    circular_contigs.add(parts[1])
    return circular_contigs

def run_minimap(query, ref, output_paf):
    # Log error minimap dihilangkan dari output standar agar bersih
    cmd = f"minimap2 -x asm5 -t 4 {ref} {query} > {output_paf} 2> /dev/null"
    subprocess.run(cmd, shell=True, check=False)

def get_coverage_score(paf_file):
    scores = {}
    if not os.path.exists(paf_file): return scores
    with open(paf_file, 'r') as f:
        for line in f:
            p = line.split('\t')
            if len(p) < 10: continue
            scores[p[0]] = scores.get(p[0], 0) + int(p[9])
    return scores

def filter_and_write_gfa(input_gfa, output_gfa, kept_ids, max_output=5):
    if not os.path.exists(input_gfa): return
    
    # 1. Ambil panjang sequence
    seq_lengths = []
    with open(input_gfa, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                if parts[1] in kept_ids:
                    seq_lengths.append((len(parts[2]), parts[1]))
    
    # 2. Sort & Filter Top N
    seq_lengths.sort(key=lambda x: x[0], reverse=True)
    final_ids = set([x[1] for x in seq_lengths[:max_output]])

    # 3. Rewrite GFA
    with open(input_gfa, 'r') as fin, open(output_gfa, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                if parts[1] in final_ids: fout.write(line)
            elif parts[0] == 'L':
                if parts[1] in final_ids and parts[3] in final_ids: fout.write(line)

# --- LOGIKA UTAMA ---

def main():
    pid = os.getpid()
    tools = ["Raven", "Flye"]
    
    # Header Log
    all_logs = [f"MODE: {target_mode} | Tool | ID | Len | Circ | Target% | Contam% | DECISION"]
    all_blacklist = []
    kept_ids_map = {"Raven": set(), "Flye": set()}

    for i, tool in enumerate(tools):
        in_fasta = input_fastas[i]
        in_gfa = input_gfas[i]
        
        paf_p = f"temp_p_{tool}_{pid}.paf"
        paf_m = f"temp_m_{tool}_{pid}.paf"
        
        try:
            circular_ids = parse_gfa_circularity(in_gfa)
            run_minimap(in_fasta, ref_plastome, paf_p)
            run_minimap(in_fasta, ref_mito, paf_m)
            
            s_plast = get_coverage_score(paf_p)
            s_mito = get_coverage_score(paf_m)
            
            for record in SeqIO.parse(in_fasta, "fasta"):
                cid = record.id
                length = len(record.seq)
                is_circ = cid in circular_ids
                
                # Hitung % Identity
                t_cov = (s_plast.get(cid, 0) / length) * 100
                c_cov = (s_mito.get(cid, 0) / length) * 100
                
                # Swap logic jika mode MITO
                if target_mode == "MITO":
                    t_cov, c_cov = c_cov, t_cov # Swap variables for easy logic
                    t_lbl, c_lbl = "Mitochondria", "Plastome"
                    is_jumbo_bad = False
                else:
                    t_lbl, c_lbl = "Plastome", "Mitochondria"
                    is_jumbo_bad = True

                # --- DECISION MATRIX (OPTIMIZED) ---
                decision = "UNKNOWN"
                action = "KEEP"

                # --- VERSI LONGGAR (LOOSE) UNTUK PENYELAMATAN DATA ---

                # 1. GOLDEN (Tetap)
                if (is_circ and valid_min <= length <= valid_max):
                    decision = f"{t_lbl} (Circular)"
                    action = "KEEP"

                # 2. SHARED MTPT (Diperluas)
                # Dulu >50%, sekarang >30% sudah kita anggap jembatan
                elif t_cov > 30 and c_cov > 30:
                    decision = "SHARED MTPT (Bridge)"
                    action = "KEEP"

                # 3. HIGH IDENTITY (Diturunkan)
                # Dulu >30%, sekarang >10% saja (asalkan lebih besar dari contam)
                elif t_cov > 10 and t_cov > c_cov:
                    decision = f"{t_lbl} (Potential)"
                    action = "KEEP"

                # 4. KONTAMINASI (Dipersempit)
                # Dulu >10%, sekarang >20% baru kita blacklist.
                # Biarkan mito yang "agak mirip" lolos dulu, nanti Flye yang urus.
                elif c_cov > 20 and c_cov > t_cov:
                    decision = f"{c_lbl} Contamination"
                    action = "BLACKLIST"

                # Sisanya (Low signal, dll) biarkan default (UNKNOWN -> KEEP)
                # Hapus atau komentar bagian "Low Signal Blacklist"
                # elif t_cov < 5 and c_cov < 5: ... (JANGAN DIAKTIFKAN)

                # Logging
                all_logs.append(f"{tool} | {cid} | {length} | {is_circ} | {t_cov:.1f} | {c_cov:.1f} | {decision} -> {action}")
                
                if action == "BLACKLIST":
                    all_blacklist.append(record)
                else:
                    kept_ids_map[tool].add(cid)
            
            filter_and_write_gfa(in_gfa, output_filtered_gfas[i], kept_ids_map[tool])

        finally:
            if os.path.exists(paf_p): os.remove(paf_p)
            if os.path.exists(paf_m): os.remove(paf_m)

    SeqIO.write(all_blacklist, output_blacklist, "fasta")
    with open(log_file, "w") as f:
        f.write("\n".join(all_logs))

if __name__ == "__main__":
    main()