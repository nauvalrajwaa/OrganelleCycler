import subprocess
import sys
import os
import random
import argparse
import shutil

# ==========================================
# KONFIGURASI TARGET ORGANEL (DEFAULT)
# ==========================================
# Nilai ini akan ditimpa oleh argumen dari command line (Config)
ESTIMATED_GENOME_SIZE = 150000 
MIN_COVERAGE = 20    
MAX_COVERAGE = 200   

# ==========================================
# BAGIAN 1: DEPENDENCIES & UTILS
# ==========================================

def check_dependencies():
    """Memastikan tools eksternal wajib terinstall."""
    required_tools = ["minimap2", "samtools"]
    for tool in required_tools:
        if subprocess.call(f"which {tool}", shell=True, stdout=subprocess.DEVNULL) != 0:
            print(f"[ERROR] Tool '{tool}' tidak ditemukan. Harap install terlebih dahulu.")
            sys.exit(1)

def get_fastq_stats(fastq_file):
    """
    Menghitung total base pair (bp) dan jumlah reads dalam file FASTQ.
    Adaptasi logika matematika round_statistics.py
    """
    total_bp = 0
    read_count = 0
    try:
        with open(fastq_file, 'r') as f:
            while True:
                # Skip header
                if not f.readline(): break
                # Baca Sequence
                seq = f.readline().strip()
                # Skip + dan Qual
                f.readline()
                f.readline()
                
                total_bp += len(seq)
                read_count += 1
    except Exception:
        return 0, 0
    return total_bp, read_count

# ==========================================
# BAGIAN 2: LOGIKA FILTERING (seq_parser.py)
# ==========================================

def calculate_avg_qual(quality_str):
    if not quality_str: return 0
    total_score = sum(ord(c) - 33 for c in quality_str)
    return total_score / len(quality_str)

def filter_reads(input_fastq, filtered_fastq, min_len=1000, min_qual=7):
    print(f"[PRE-PROCESS] Memfilter reads mentah...")
    print(f"   -> Kriteria: Panjang >= {min_len} bp, Avg Quality >= {min_qual}")
    
    kept_reads = 0
    try:
        with open(input_fastq, 'r') as fin, open(filtered_fastq, 'w') as fout:
            while True:
                header = fin.readline()
                if not header: break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                
                if len(seq.strip()) < min_len: continue
                
                # Hitung kualitas hanya jika lolos panjang (optimasi)
                if calculate_avg_qual(qual.strip()) >= min_qual:
                    fout.write(header + seq + plus + qual)
                    kept_reads += 1
                    
    except IOError as e:
        print(f"[ERROR] I/O Error: {e}")
        return None

    print(f"   -> Selesai. Reads valid tersimpan: {kept_reads}")
    if kept_reads == 0: return None
    return filtered_fastq

# ==========================================
# BAGIAN 3: SUBSAMPLING (Adaptasi Logic Coverage Control)
# ==========================================

def subsample_reads(input_fastq, output_fastq, target_bp):
    """
    Mengurangi jumlah reads jika data terlalu besar (Overload).
    Mencegah Flye crash karena memori penuh.
    """
    print(f"[OPTIMIZER] Melakukan Subsampling (Downsampling)...")
    print(f"   -> Target: ~{target_bp/1_000_000:.2f} Mb data")
    
    total_bp, total_reads = get_fastq_stats(input_fastq)
    if total_bp == 0: return False
    
    # Probabilitas ambil read
    keep_ratio = target_bp / total_bp
    
    kept_count = 0
    with open(input_fastq, 'r') as fin, open(output_fastq, 'w') as fout:
        while True:
            header = fin.readline()
            if not header: break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            
            if random.random() < keep_ratio:
                fout.write(header + seq + plus + qual)
                kept_count += 1
                
    print(f"   -> Subsampling selesai. Menyimpan {kept_count} reads.")
    return True

# ==========================================
# BAGIAN 4: PIPELINE UTAMA
# ==========================================

def run_recruitment(seed_fasta, total_reads_input, output_dir, threads=8, min_len=1000, min_qual=7):
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    # 1. FILTERING
    clean_input = os.path.join(output_dir, "clean_input.fastq")
    # Update: Pass min_len and min_qual arguments
    if not filter_reads(total_reads_input, clean_input, min_len, min_qual):
        return None
        
    # 2. MAPPING (MINIMAP2)
    # Gunakan clean_input agar lebih cepat
    raw_recruited = os.path.join(output_dir, "raw_recruited.fastq")
    print(f"[RECRUITER] Memulai Minimap2 (Seed-based)...")
    
    cmd = (
        f"minimap2 -ax map-ont -Y -t {threads} --secondary=no {seed_fasta} {clean_input} "
        f"| samtools view -b -F 4 -q 10 -@ 4 "
        f"| samtools fastq -@ 4 > {raw_recruited}"
    )
    
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        print("[ERROR] Mapping gagal.")
        return None

    # 3. STATISTIK & SUBSAMPLING (Adaptasi round_statistics.py)
    final_fastq = os.path.join(output_dir, "final_recruited.fastq")
    
    total_bp, count = get_fastq_stats(raw_recruited)
    coverage = total_bp / ESTIMATED_GENOME_SIZE
    
    print(f"[STATS] Total Data Terekrut: {total_bp/1000:.1f} kb ({count} reads)")
    print(f"[STATS] Estimasi Coverage Organel: {coverage:.2f}x")
    
    if coverage < MIN_COVERAGE:
        print(f"[WARNING] Coverage sangat rendah (<{MIN_COVERAGE}x). Assembly mungkin gagal/terputus.")
        # Tetap lanjut, tapi user diberi tahu
        os.rename(raw_recruited, final_fastq)
        
    elif coverage > MAX_COVERAGE:
        print(f"[INFO] Coverage terlalu tinggi (>{MAX_COVERAGE}x). Melakukan subsampling agar Flye optimal.")
        # Hitung target BP (misal: 200x * 150kb = 30MB)
        target_bp = MAX_COVERAGE * ESTIMATED_GENOME_SIZE
        subsample_reads(raw_recruited, final_fastq, target_bp)
        
    else:
        print(f"[INFO] Coverage optimal. Siap untuk assembly.")
        os.rename(raw_recruited, final_fastq)
    
    # Bersihkan file temporary (kecuali hasil akhir)
    if os.path.exists(clean_input): os.remove(clean_input)
    if os.path.exists(raw_recruited): os.remove(raw_recruited)
    
    return final_fastq

# ==========================================
# MAIN EXECUTION (ARGPARSE UPDATE)
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organelle Reads Recruiter (GetOrganelle Logic)")
    
    # Required Arguments
    parser.add_argument("--bait", required=True, help="Path to Hybrid Bait Fasta")
    parser.add_argument("--reads", required=True, help="Path to Raw Reads Fastq")
    parser.add_argument("--output", required=True, help="Path to Output Concentrated Fastq")
    
    # Optional Arguments (Configurable via Snakemake)
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--est_size", type=int, default=150000, help="Estimated Genome Size")
    parser.add_argument("--max_cov", type=int, default=200, help="Target Max Coverage for Subsampling")
    parser.add_argument("--min_len", type=int, default=1000, help="Min Read Length Filter")
    parser.add_argument("--min_qual", type=int, default=7, help="Min Read Quality Filter")

    args = parser.parse_args()

    # Update Global Variables based on Arguments
    ESTIMATED_GENOME_SIZE = args.est_size
    MAX_COVERAGE = args.max_cov
    # Min coverage hardcoded atau bisa ditambah argumen jika perlu, default 20 aman.

    # Check Tools
    check_dependencies()

    # Prepare Temp Directory
    # Menggunakan path output untuk membuat folder temp disebelahnya
    temp_dir = os.path.dirname(args.output) + f"/temp_recruitment_{os.path.basename(args.output)}"
    
    print(f"[START] Recruitment pipeline started...")
    print(f"[CONFIG] Size: {args.est_size}bp | MaxCov: {args.max_cov}x | Threads: {args.threads}")

    # Run Pipeline
    final_out = run_recruitment(
        seed_fasta=args.bait, 
        total_reads_input=args.reads, 
        output_dir=temp_dir, 
        threads=args.threads,
        min_len=args.min_len,
        min_qual=args.min_qual
    )
    
    # Finalize
    if final_out and os.path.exists(final_out):
        # Move final result to destination
        shutil.move(final_out, args.output)
        print(f"[SUCCESS] Result saved to: {args.output}")
        # Cleanup directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    else:
        print("[FAIL] Pipeline failed to produce output.")
        sys.exit(1)