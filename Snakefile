import pandas as pd
import os

# =============================================================================
# SNAKEMAKE PIPELINE: ORGANELLE ASSEMBLY (PLASTOME/MITOME)
# Version: 2.1 (Resilient Config)
# =============================================================================

configfile: "config/config.yaml"

# -----------------------------------------------------------------------------
# 1. SETUP MODE & VARIABEL DINAMIS
# -----------------------------------------------------------------------------

# Baca Mode dari Config (gunakan 'run_mode' sesuai YAML Anda)
mode = config.get("run_mode", "PLASTOME").upper()

print("="*50)
print(f"   🚀 RUNNING PIPELINE IN MODE: [ {mode} ]")
print("="*50)

# Default Metadata (Agar script tidak error jika 'seed' tidak ada di config)
DEFAULTS = {
    "PLASTOME": {
        "seed": "embplant_pt.fasta", 
        "label": "embplant_pt.fasta", 
        "size": "150k"
    },
    "MITOME": {
        "seed": "embplant_mt.fasta", 
        "label": "embplant_mt.fasta", 
        "size": "400k"
    }
}

# Logika Penentuan Path & Variabel
if mode == "PLASTOME":
    # 1. Path Referensi (Sesuai struktur config Anda)
    target_ref = config["paths"]["plastome"]["specific"]["fasta"]
    target_gbk = config["paths"]["plastome"]["specific"]["gbk"]
    
    # 2. Metadata (Coba ambil dari config, jika kosong pakai DEFAULTS)
    # Ini mencegah KeyError 'seed'
    plastome_conf = config["paths"]["plastome"]
    active_seed  = plastome_conf.get("seed", DEFAULTS["PLASTOME"]["seed"])
    active_label = plastome_conf.get("label", DEFAULTS["PLASTOME"]["label"])
    target_size  = plastome_conf.get("size", DEFAULTS["PLASTOME"]["size"])
    
    # 3. Output Final
    final_assembly = "results/assembly_plastome/final_circular.fasta"

elif mode == "MITOME":
    # 1. Path Referensi
    target_ref = config["paths"]["mitome"]["specific"]["fasta"]
    target_gbk = config["paths"]["mitome"]["specific"]["gbk"]
    
    # 2. Metadata
    mitome_conf = config["paths"]["mitome"]
    active_seed  = mitome_conf.get("seed", DEFAULTS["MITOME"]["seed"])
    active_label = mitome_conf.get("label", DEFAULTS["MITOME"]["label"])
    target_size  = mitome_conf.get("size", DEFAULTS["MITOME"]["size"])
    
    # 3. Output Final
    final_assembly = "results/assembly_mitome/final_circular.fasta"

else:
    raise ValueError(f"Mode tidak dikenal: {mode}. Gunakan PLASTOME atau MITOME di config.")

# Simpan variabel ke config global agar bisa diakses oleh rule lain (rules/*.smk)
config["active_seed"]  = active_seed
config["active_label"] = active_label
config["target_ref"]   = target_ref
config["target_gbk"]   = target_gbk
config["target_size"]  = target_size

# --- DEBUG INFO ---
print(f"[INFO] Target Reference    : {target_ref}")
print(f"[INFO] Reference GenBank   : {target_gbk}")
print(f"[INFO] GetOrganelle Seed   : {active_seed}")
print("-" * 50)

# =============================================================================
# 2. SETUP DATA & ASSEMBLERS
# =============================================================================
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)

# Mengambil list assembler dari config (misal: "flye, raven")
ASSEMBLERS = [x.strip() for x in config.get("loop_mitohifi", "flye").split(",")]

# =============================================================================
# 3. TARGET DEFINITION (Output yang diinginkan)
# =============================================================================
TARGETS = []

# --- A. VISUALISASI INTERMEDIATE (ALWAYS RUN) ---
# Memetakan hasil assembly Flye/Raven yang sudah dibersihkan ke referensi
TARGETS.extend(expand("results/{sample}/09_viz/{assembler}_alignment_map.html", 
                      sample=samples.index, assembler=ASSEMBLERS))

# --- B. HASIL AKHIR (Rescue/Circularization) ---
# Kita targetkan hasil akhir per sample
TARGETS.extend(expand("results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta", 
                      sample=samples.index))

# --- C. VISUALISASI FINAL ---
# Memetakan hasil akhir (Final Rescue) ke referensi
TARGETS.extend(expand("results/{sample}/09_viz/final_alignment_map.html", 
                      sample=samples.index))

# =============================================================================
# 4. PIPELINE ENTRY POINT
# =============================================================================

rule all:
    input:
        # Menjalankan pipeline untuk semua sampel yang ada di samples.tsv
        TARGETS

# =============================================================================
# 5. INCLUDE MODULES
# =============================================================================

include: "rules/00_prep.smk"       # Download Ref & GetOrganelle
include: "rules/01_rough.smk"      # Hybrid Baiting & Rough Assembly
include: "rules/02_assembly.smk"   # Graph Cleaning & Solving
include: "rules/03_qc.smk"         # QC intermediate
include: "rules/04_rescue.smk"     # MitoHiFi & Final Selection
include: "rules/05_visualization.smk" # Plotting