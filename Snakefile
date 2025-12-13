configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# SNAKEMAKE PIPELINE: ORGANELLE ASSEMBLY (PLASTOME/MITOME)
# Version: 2.0 (Smart Switching)
# =============================================================================

# --- 1. SETUP KONFIGURASI DINAMIS ---
# Bagian ini membaca 'run_mode' dari config.yaml dan menentukan
# path mana yang menjadi TARGET, BLACKLIST, SEED, dan LABEL.

MODE = config["run_mode"].upper()

print(f"\n{'='*50}")
print(f"   🚀 RUNNING PIPELINE IN MODE: [ {MODE} ]")
print(f"{'='*50}\n")

if MODE == "PLASTOME":
    # --- MODE PLASTOME ---
    # Target utama adalah Plastome. 
    config["target_ref"]      = config["paths"]["plastome"]["fasta"]
    config["target_gbk"]      = config["paths"]["plastome"]["gbk"]
    config["blacklist_ref"]   = config["paths"]["mitome"]["fasta"]
    config["target_est_size"] = config["paths"]["plastome"]["size"]
    
    # Seed untuk GetOrganelle (Initial Recruitment)
    config["active_seed"]     = config["paths"]["plastome"]["seed"]
    
    # Label untuk Graph Cleaner (Gene Database)
    config["active_label"]    = config["paths"]["plastome"]["label"]
    
    # Label internal
    config["target_assembly"] = "PLASTOME" 

elif MODE == "MITOME":
    # --- MODE MITOME ---
    # Target utama adalah Mitome.
    config["target_ref"]      = config["paths"]["mitome"]["fasta"]
    config["target_gbk"]      = config["paths"]["mitome"]["gbk"]
    config["blacklist_ref"]   = config["paths"]["plastome"]["fasta"]
    config["target_est_size"] = config["paths"]["mitome"]["size"]
    
    # Seed untuk GetOrganelle (Initial Recruitment)
    config["active_seed"]     = config["paths"]["mitome"]["seed"]
    
    # Label untuk Graph Cleaner (Gene Database)
    config["active_label"]    = config["paths"]["mitome"]["label"]
    
    # Label internal
    config["target_assembly"] = "MITOME"

else:
    # Error Handling
    raise ValueError(f"❌ FATAL ERROR: Mode '{MODE}' tidak dikenali! Harap set 'run_mode' di config.yaml menjadi 'PLASTOME' atau 'MITOME'.")

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
    input: TARGETS

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