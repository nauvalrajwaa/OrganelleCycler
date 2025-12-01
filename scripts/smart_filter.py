import sys
import os
import subprocess
from Bio import SeqIO

# --- CONFIG & INPUTS ---
# Fallback logic agar script bisa ditest tanpa snakemake
if 'snakemake' in globals():
    input_fasta = snakemake.input.fasta
    input_gfa = snakemake.input.gfa
    ref_plastome = snakemake.input.ref_p
    ref_mito = snakemake.input.ref_m
    output_blacklist = snakemake.output.blacklist
    log_file = snakemake.output.log
    min_len = snakemake.config["plastome_min"]
    max_len = snakemake.config["plastome_max"]
else:
    print("Run via Snakemake only.")
    sys.exit(1)

# --- FUNGSI BANTUAN ---

def parse_gfa_circularity(gfa_file):
    """Mendeteksi contig sirkular (Self-loop) dari GFA Raven"""
    circular_contigs = set()
    try:
        with open(gfa_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # Format Link: L seg1 orient seg2 orient ...
                if parts[0] == 'L':
                    seg1 = parts[1]
                    seg2 = parts[3]
                    if seg1 == seg2: # Self Loop -> Circular
                        circular_contigs.add(seg1)
    except Exception as e:
        print(f"Warning parsing GFA: {e}")
    return circular_contigs

def run_minimap(query, ref, output_paf):
    """Jalankan minimap2 untuk cek coverage"""
    # -x asm5: assembly-to-assembly alignment
    cmd = f"minimap2 -x asm5 -t 4 {ref} {query} > {output_paf} 2> /dev/null"
    subprocess.run(cmd, shell=True, check=False)

def get_coverage_score(paf_file):
    """Hitung berapa base yg match dari PAF"""
    scores = {}
    if not os.path.exists(paf_file): return scores
    with open(paf_file, 'r') as f:
        for line in f:
            p = line.split('\t')
            ctg_id = p[0]
            # Kolom 10 = Number of residue matches
            matches = int(p[9]) 
            scores[ctg_id] = scores.get(ctg_id, 0) + matches
    return scores

# --- LOGIKA UTAMA ---

def main():
    # 1. Analisis Graph
    circular_ids = parse_gfa_circularity(input_gfa)
    
    # 2. Mapping Identity (Minimap2)
    paf_plastid = "temp_plastid.paf"
    paf_mito = "temp_mito.paf"
    
    run_minimap(input_fasta, ref_plastome, paf_plastid)
    run_minimap(input_fasta, ref_mito, paf_mito)
    
    plastid_matches = get_coverage_score(paf_plastid)
    mito_matches = get_coverage_score(paf_mito)
    
    blacklist_records = []
    logs = ["ContigID | Length | Circular | PlastidScore | MitoScore | DECISION"]
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        cid = record.id
        length = len(record.seq)
        is_circ = cid in circular_ids
        
        # Hitung % Identity (Coverage terhadap panjang contig sendiri)
        p_cov = (plastid_matches.get(cid, 0) / length) * 100
        m_cov = (mito_matches.get(cid, 0) / length) * 100
        
        decision = "UNKNOWN"
        action = "KEEP" # Default aman
        
        # --- ATURAN MAIN (DECISION MATRIX) ---
        
        # A. PLASTOME PASTI (KEEP)
        # Circular & Size range wajar ATAU Linear tapi match plastid tinggi
        if (is_circ and min_len <= length <= max_len):
            decision = "Plastome (Circular)"
            action = "KEEP"
        elif p_cov > 30 and p_cov > m_cov:
            decision = "Plastome (High Identity)"
            action = "KEEP"
            
        # B. MITOKONDRIA / SAMPAH (BLACKLIST)
        # Match mito dominan
        elif m_cov > 20 and m_cov > p_cov:
            decision = "Mitochondria"
            action = "BLACKLIST"
            
        # C. UKURAN JUMBO TAPI BUKAN PLASTID
        # Misal > 200kb, tidak sirkular, plastid cov rendah
        elif length > 250000 and p_cov < 10:
            decision = "Jumbo Junk"
            action = "BLACKLIST"
            
        # D. FRAGMEN KECIL
        elif length < 20000: # Di bawah 20kb
            decision = "Small Fragment"
            action = "BLACKLIST"
            
        # Logging
        logs.append(f"{cid} | {length} | {is_circ} | {p_cov:.1f}% | {m_cov:.1f}% | {decision} -> {action}")
        
        if action == "BLACKLIST":
            blacklist_records.append(record)
            
    # Output Files
    SeqIO.write(blacklist_records, output_blacklist, "fasta")
    
    with open(log_file, "w") as f:
        f.write("\n".join(logs))

    # Cleanup Temp
    if os.path.exists(paf_plastid): os.remove(paf_plastid)
    if os.path.exists(paf_mito): os.remove(paf_mito)

if __name__ == "__main__":
    main()
