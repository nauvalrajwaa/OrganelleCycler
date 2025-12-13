import pandas as pd
import os

# =============================================================================
# SNAKEMAKE PIPELINE: ORGANELLE ASSEMBLY
# Version: 2.8 (Fixed: Deep Injection for 00_prep.smk)
# =============================================================================

configfile: "config/config.yaml"

# 1. SETUP MODE
# -----------------------------------------------------------------------------
mode = config.get("run_mode", "PLASTOME").upper()

print(f"🚀 RUNNING PIPELINE IN MODE: [ {mode} ]")

# 2. VARIABEL DINAMIS & DEFAULTS
# -----------------------------------------------------------------------------
DEFAULTS = {
    "PLASTOME": { "seed": "embplant_pt.fasta", "label": "embplant_pt.fasta", "size": "150k", "int_size": 150000 },
    "MITOME":   { "seed": "embplant_mt.fasta", "label": "embplant_mt.fasta", "size": "400k", "int_size": 400000 }
}

if mode == "PLASTOME":
    target_ref = config["paths"]["plastome"]["specific"]["fasta"]
    target_gbk = config["paths"]["plastome"]["specific"]["gbk"]
    blacklist_ref = config["paths"]["mitome"]["specific"]["fasta"]
    
    p_conf = config["paths"]["plastome"]
    active_seed  = p_conf.get("seed", DEFAULTS["PLASTOME"]["seed"])
    active_label = p_conf.get("label", DEFAULTS["PLASTOME"]["label"])
    target_size  = p_conf.get("size", DEFAULTS["PLASTOME"]["size"])
    default_int_size = DEFAULTS["PLASTOME"]["int_size"]

elif mode == "MITOME":
    target_ref = config["paths"]["mitome"]["specific"]["fasta"]
    target_gbk = config["paths"]["mitome"]["specific"]["gbk"]
    blacklist_ref = config["paths"]["plastome"]["specific"]["fasta"]
    
    m_conf = config["paths"]["mitome"]
    active_seed  = m_conf.get("seed", DEFAULTS["MITOME"]["seed"])
    active_label = m_conf.get("label", DEFAULTS["MITOME"]["label"])
    target_size  = m_conf.get("size", DEFAULTS["MITOME"]["size"])
    default_int_size = DEFAULTS["MITOME"]["int_size"]

else:
    raise ValueError(f"Mode {mode} tidak dikenal.")

# -----------------------------------------------------------------------------
# [CRITICAL FIX] MEMANIPULASI STRUKTUR CONFIG
# -----------------------------------------------------------------------------

# 1. Inject Global variables (untuk rules lain)
config["active_seed"]   = active_seed
config["active_label"]  = active_label
config["target_ref"]    = target_ref
config["target_gbk"]    = target_gbk
config["blacklist_ref"] = blacklist_ref
config["ref_gb"]        = target_gbk

# 2. FIX KHUSUS UNTUK 00_prep.smk (Wajib ada di dalam ['recruitment'])
# Pastikan section 'recruitment' ada dulu
if "recruitment" not in config:
    config["recruitment"] = {}

# Cek apakah 'est_genome_size' sudah ada di dalam 'recruitment'
# Jika TIDAK ADA, kita paksa masukkan defaultnya ke SANA.
if "est_genome_size" not in config["recruitment"]:
    config["recruitment"]["est_genome_size"] = default_int_size
    print(f"[WARNING] 'est_genome_size' disisipkan otomatis ke config['recruitment']: {default_int_size}")

# Jembatan untuk 03_qc.smk (yang mungkin pakai nama variabel lain)
config["target_est_size"] = config["recruitment"]["est_genome_size"]

# Debug Info
# =============================================================================
# CONSOLE LOGGING (STATUS REPORT)
# =============================================================================
print("\n" + "="*70)
print(f"🚀  ORGANELLE CYCLER PIPELINE REPORT")
print("="*70)
print(f"{'🔹 RUN MODE':<25} : [ {mode} ]")
print(f"{'🔹 EST. GENOME SIZE':<25} : {config['recruitment']['est_genome_size']} bp")
print("-" * 70)
print("📂 REFERENCE FILES:")
print(f"   - Target Ref (Fasta)   : {target_ref}")
print(f"   - Target Ref (GBK)     : {target_gbk}")
print(f"   - Blacklist (Negative) : {blacklist_ref}")
print(f"   - Seed Database        : {active_seed}")
print(f"   - Label Database       : {active_label}")
print("-" * 70)
print("⚙️  CONFIGURATION:")
print(f"   - Active Assemblers    : {', '.join(ASSEMBLERS)}")
print(f"   - Total Samples        : {len(samples)} sample(s)")
print(f"   - Samples List         : {', '.join(samples.index[:5])}{'...' if len(samples) > 5 else ''}")
print("="*70 + "\n")

# 3. SETUP TARGETS
# -----------------------------------------------------------------------------
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)
ASSEMBLERS = [x.strip() for x in config.get("loop_mitohifi", "flye").split(",")]

TARGETS = []
TARGETS.extend(expand("results/{sample}/09_viz/{assembler}_alignment_map.html", 
                      sample=samples.index, assembler=ASSEMBLERS))
TARGETS.extend(expand("results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta", 
                      sample=samples.index))
TARGETS.extend(expand("results/{sample}/09_viz/final_alignment_map.html", 
                      sample=samples.index))

rule all:
    input: TARGETS

# 4. MODULES
# -----------------------------------------------------------------------------
include: "rules/00_prep.smk"
include: "rules/01_rough.smk"
include: "rules/02_assembly.smk"
include: "rules/03_qc.smk"
include: "rules/04_rescue.smk"
include: "rules/05_visualization.smk"