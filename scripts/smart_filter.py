# /scripts/smart_filter.py

import sys
import os
import subprocess
from Bio import SeqIO

# --- CONFIG & INPUTS ---
if 'snakemake' in globals():
    # Input list (Urutan harus konsisten: Raven dulu, baru Flye)
    input_fastas = snakemake.input.fastas 
    input_gfas = snakemake.input.gfas 
    
    ref_plastome = snakemake.input.ref_p
    ref_mito = snakemake.input.ref_m
    
    # Output
    output_blacklist = snakemake.output.blacklist
    output_filtered_gfas = snakemake.output.filtered_gfas # [filtered_raven.gfa, filtered_flye.gfa]
    log_file = snakemake.output.log
    
    # Config
    config = snakemake.config
    target_mode = config.get("target_assembly", "PLASTOME")
    
    if target_mode == "MITO":
        valid_min = config.get("mito_min_len", 200000)
        valid_max = config.get("mito_max_len", 2000000) 
    else:
        valid_min = config.get("plastome_min", 100000)
        valid_max = config.get("plastome_max", 200000)
else:
    sys.exit("Run via Snakemake only.")

# --- FUNGSI BANTUAN ---

def parse_gfa_circularity(gfa_file):
    """Mendeteksi contig sirkular"""
    circular_contigs = set()
    if os.path.exists(gfa_file):
        with open(gfa_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[0] == 'L' and parts[1] == parts[3]:
                    circular_contigs.add(parts[1])
    return circular_contigs

def run_minimap(query, ref, output_paf):
    # Hapus 2> /dev/null agar log terlihat
    cmd = f"minimap2 -x asm5 -t 4 {ref} {query} > {output_paf}"
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
    """
    Membaca GFA asli, hanya menulis node yang ada di kept_ids.
    Dibatasi maksimal 'max_output' contig terpanjang.
    """
    if not os.path.exists(input_gfa): return

    # 1. Identifikasi contig mana yang akan ditulis (Top N Longest)
    # Kita perlu baca S-lines dulu untuk tahu panjangnya
    seq_lengths = []
    with open(input_gfa, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                seg_id = parts[1]
                # Sequence ada di index 2
                seq = parts[2]
                if seg_id in kept_ids:
                    seq_lengths.append((len(seq), seg_id))
    
    # Urutkan dari yang terpanjang, ambil Top N
    seq_lengths.sort(key=lambda x: x[0], reverse=True)
    final_allowed_ids = set([x[1] for x in seq_lengths[:max_output]])

    # 2. Tulis ulang GFA
    with open(input_gfa, 'r') as fin, open(output_gfa, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            # Handle Segment (S)
            if parts[0] == 'S':
                if parts[1] in final_allowed_ids:
                    fout.write(line)
            # Handle Link (L)
            elif parts[0] == 'L':
                # Tulis link hanya jika KEDUA node ada di final list
                if parts[1] in final_allowed_ids and parts[3] in final_allowed_ids:
                    fout.write(line)

# --- LOGIKA UTAMA ---

def main():
    pid = os.getpid()
    
    # Kita proses per-tool agar scoring adil & GFA terpisah
    # Input list diasumsikan: [0]=Raven, [1]=Flye
    tools = ["Raven", "Flye"]
    
    # Header Log (Format Baru)
    all_logs = [f"MODE: {target_mode} | Tool | ID | Len | Circ | Target% | Contam% | DECISION"]
    all_blacklist = []
    
    # Set global untuk menyimpan ID yang lolos (Keep) per tool
    # Format: {"Raven": {id1, id2}, "Flye": {id3}}
    kept_ids_map = {"Raven": set(), "Flye": set()}

    for i, tool in enumerate(tools):
        in_fasta = input_fastas[i]
        in_gfa = input_gfas[i]
        
        # Temp files
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
                
                p_cov = (s_plast.get(cid, 0) / length) * 100
                m_cov = (s_mito.get(cid, 0) / length) * 100
                
                decision = "UNKNOWN"
                action = "KEEP" # Default
                
                # Logic Switcher
                if target_mode == "MITO":
                    t_cov, c_cov = m_cov, p_cov
                    t_lbl, c_lbl = "Mitochondria", "Plastome"
                    is_jumbo_bad = False
                else:
                    t_cov, c_cov = p_cov, m_cov
                    t_lbl, c_lbl = "Plastome", "Mitochondria"
                    is_jumbo_bad = True
                
                # MATRIX DECISION
                if (is_circ and valid_min <= length <= valid_max):
                    decision = f"{t_lbl} (Circular)"
                    action = "KEEP"
                elif t_cov > 30 and t_cov > c_cov:
                    decision = f"{t_lbl} (High Identity)"
                    action = "KEEP"
                elif c_cov > 20 and c_cov > t_cov:
                    decision = f"{c_lbl} Contamination"
                    action = "BLACKLIST"
                elif is_jumbo_bad and length > 250000 and t_cov < 10:
                    decision = "Jumbo Junk"
                    action = "BLACKLIST"
                elif length < 20000: 
                    decision = "Small Fragment"
                    action = "BLACKLIST"
                
                # Simpan Log
                all_logs.append(f"{tool} | {cid} | {length} | {is_circ} | {t_cov:.1f} | {c_cov:.1f} | {decision} -> {action}")
                
                if action == "BLACKLIST":
                    all_blacklist.append(record)
                else:
                    # Simpan ID yang lolos untuk filtering GFA
                    kept_ids_map[tool].add(cid)
            
            # --- FILTER GFA LANGSUNG ---
            # Menulis GFA baru yang hanya berisi graph KEEP (Max 5 terpanjang)
            filter_and_write_gfa(in_gfa, output_filtered_gfas[i], kept_ids_map[tool], max_output=5)

        finally:
            if os.path.exists(paf_p): os.remove(paf_p)
            if os.path.exists(paf_m): os.remove(paf_m)

    # Output Final
    SeqIO.write(all_blacklist, output_blacklist, "fasta")
    with open(log_file, "w") as f:
        f.write("\n".join(all_logs))

if __name__ == "__main__":
    main()