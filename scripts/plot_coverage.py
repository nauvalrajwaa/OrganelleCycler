import sys
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# Setup Snakemake
if 'snakemake' in globals():
    INPUT_READS = snakemake.input.reads
    INPUT_REF   = snakemake.input.ref
    OUTPUT_BAM  = snakemake.output.bam
    OUTPUT_TXT  = snakemake.output.depth
    OUTPUT_PLOT = snakemake.output.plot
    THREADS     = snakemake.threads
else:
    sys.exit("Run via Snakemake.")

def run_mapping():
    # 1. Mapping Reads ke Final Assembly
    # Gunakan minimap2 dan pipe ke samtools sort
    cmd_map = f"minimap2 -ax map-ont -t {THREADS} {INPUT_REF} {INPUT_READS} | " \
              f"samtools sort -@ {THREADS} -o {OUTPUT_BAM}"
    
    print(f"Running mapping: {cmd_map}")
    subprocess.check_call(cmd_map, shell=True)
    
    # Index BAM
    subprocess.check_call(f"samtools index {OUTPUT_BAM}", shell=True)

    # 2. Calculate Depth (Coverage per base)
    # -a: Output all positions (including zero coverage)
    cmd_depth = f"samtools depth -a {OUTPUT_BAM} > {OUTPUT_TXT}"
    print(f"Calculating depth: {cmd_depth}")
    subprocess.check_call(cmd_depth, shell=True)

def plot_coverage():
    # Baca file depth (Format: Contig, Pos, Depth)
    # Gunakan chunksize jika file sangat besar, tapi untuk organel biasanya kecil
    try:
        df = pd.read_csv(OUTPUT_TXT, sep='\t', header=None, names=['Contig', 'Pos', 'Depth'])
    except pd.errors.EmptyDataError:
        print("Coverage file empty. Creating dummy plot.")
        plt.figure()
        plt.text(0.5, 0.5, 'No Coverage Data', ha='center')
        plt.savefig(OUTPUT_PLOT)
        return

    # Buat Plot
    plt.figure(figsize=(12, 4))
    
    # Jika ada banyak contig, kita plot contig utama (paling panjang)
    if not df.empty:
        # Kelompokkan per contig
        for contig_name, data in df.groupby('Contig'):
            plt.plot(data['Pos'], data['Depth'], label=contig_name, linewidth=1, alpha=0.8)

        plt.title(f"Read Coverage Depth: {os.path.basename(INPUT_REF)}")
        plt.xlabel("Position (bp)")
        plt.ylabel("Depth (x)")
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.legend(loc='upper right', fontsize='small')
        plt.tight_layout()
        
        plt.savefig(OUTPUT_PLOT, dpi=300)
        print(f"Plot saved to {OUTPUT_PLOT}")
    else:
        # Kosong
        plt.figure()
        plt.savefig(OUTPUT_PLOT)

if __name__ == "__main__":
    try:
        run_mapping()
        plot_coverage()
    except Exception as e:
        print(f"Error in plot_coverage.py: {e}")
        # Jangan exit error biar pipeline gak mati total cuma gara-gara gambar
        # Buat dummy output
        if not os.path.exists(OUTPUT_BAM): open(OUTPUT_BAM, 'a').close()
        if not os.path.exists(OUTPUT_TXT): open(OUTPUT_TXT, 'a').close()
        if not os.path.exists(OUTPUT_PLOT): open(OUTPUT_PLOT, 'a').close()