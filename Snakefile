configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# SNAKEMAKE PIPELINE: ORGANELLE ASSEMBLY (PLASTOME/MITOME)
# Version: 2.0 (Smart Switching)
# =============================================================================

# 1. Tentukan Mode Operasi (PLASTOME / MITOME)
# Default ke PLASTOME jika tidak ditentukan lewat command line (--config mode=...)
mode = config.get("mode", "PLASTOME").upper()

print("="*50)
print(f"   🚀 RUNNING PIPELINE IN MODE: [ {mode} ]")
print("="*50)

# 2. Setup Variabel Dinamis Berdasarkan Mode
# Kita sesuaikan dengan struktur baru di config.yaml (nested paths)

if mode == "PLASTOME":
    # Akses ke dalam config["paths"]["plastome"]
    # Perhatikan struktur 'specific' untuk referensi
    target_ref  = config["paths"]["plastome"]["specific"]["fasta"]
    target_gbk  = config["paths"]["plastome"]["specific"]["gbk"]
    
    # Metadata (Seed, Label, Size) ada di level plastome (bukan specific)
    active_seed = config["paths"]["plastome"]["seed"]
    active_label= config["paths"]["plastome"]["label"]
    target_size = config["paths"]["plastome"]["size"]
    
    # Path output final (opsional, untuk rule all)
    final_assembly = "results/assembly_plastome/final_circular.fasta"

elif mode == "MITOME":
    # Akses ke dalam config["paths"]["mitome"]
    target_ref  = config["paths"]["mitome"]["specific"]["fasta"]
    target_gbk  = config["paths"]["mitome"]["specific"]["gbk"]
    
    active_seed = config["paths"]["mitome"]["seed"]
    active_label= config["paths"]["mitome"]["label"]
    target_size = config["paths"]["mitome"]["size"]
    
    final_assembly = "results/assembly_mitome/final_circular.fasta"

else:
    raise ValueError(f"Mode tidak dikenal: {mode}. Gunakan PLASTOME atau MITOME.")

# Simpan variabel ke config global agar bisa diakses oleh rule lain jika perlu
config["active_seed"] = active_seed
config["active_label"] = active_label
config["target_ref"] = target_ref
config["target_gbk"] = target_gbk  # <-- Penting untuk rule extract_genes

# --- DEBUG INFO ---
print(f"[INFO] Target Reference    : {config['target_ref']}")
print(f"[INFO] Cleaner Gene Label  : {config['active_label']}")
print(f"[INFO] GetOrganelle Seed   : {config['active_seed']}")
print("-" * 50)

# =============================================================================
# 2. SETUP GLOBAL
# =============================================================================
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)
ASSEMBLERS = [x.strip() for x in config.get("loop_mitohifi", "flye").split(",")]

# =============================================================================
# 3. TARGET DEFINITION (Output yang diinginkan)
# =============================================================================
TARGETS = []

# --- A. VISUALISASI INTERMEDIATE (ALWAYS RUN) ---
# Memetakan hasil assembly Flye/Raven yang sudah dibersihkan (Step 02) ke referensi
TARGETS.extend(expand("results/{sample}/09_viz/{assembler}_alignment_map.html", 
                      sample=samples.index, assembler=ASSEMBLERS))

# --- B. HASIL AKHIR ---
# Kita selalu mengincar hasil 'Final Rescue'.
# Jika rescue dimatikan di config, logic internal 04_rescue.smk akan handle fallback-nya.
TARGETS.extend(expand("results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta", 
                      sample=samples.index))

# --- C. VISUALISASI FINAL ---
# Memetakan hasil akhir (Final Rescue) ke referensi
TARGETS.extend(expand("results/{sample}/09_viz/final_alignment_map.html", 
                      sample=samples.index))


rule all:
    input:
        final_assembly

# =============================================================================
# 4. INCLUDE MODULES
# =============================================================================
# Urutan include tidak terlalu berpengaruh pada eksekusi, tapi memudahkan pembacaan.

include: "rules/00_prep.smk"       # Download Ref & GetOrganelle
include: "rules/01_rough.smk"      # Hybrid Baiting & Rough Assembly (Flye/Raven Raw)
include: "rules/02_assembly.smk"   # Graph Cleaning & Solving
include: "rules/03_qc.smk"         # (Opsional) QC intermediate
include: "rules/04_rescue.smk"     # MitoHiFi & Final Selection
include: "rules/05_visualization.smk" # Plotting