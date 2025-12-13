import pandas as pd
import os

# =============================================================================
# SNAKEMAKE PIPELINE: ORGANELLE ASSEMBLY
# Version: 2.3 (Clean & Simple)
# =============================================================================

configfile: "config/config.yaml"

# 1. SETUP MODE
# -----------------------------------------------------------------------------
mode = config.get("run_mode", "PLASTOME").upper()

print(f"🚀 RUNNING PIPELINE IN MODE: [ {mode} ]")

# 2. VARIABEL DINAMIS (PATH)
# -----------------------------------------------------------------------------
# Default Values (Agar tidak crash jika config kosong)
DEFAULTS = {
    "PLASTOME": { "seed": "embplant_pt.fasta", "label": "embplant_pt.fasta", "size": "150k" },
    "MITOME":   { "seed": "embplant_mt.fasta", "label": "embplant_mt.fasta", "size": "400k" }
}

if mode == "PLASTOME":
    target_ref = config["paths"]["plastome"]["specific"]["fasta"]
    target_gbk = config["paths"]["plastome"]["specific"]["gbk"]
    
    # Ambil metadata aman
    p_conf = config["paths"]["plastome"]
    active_seed  = p_conf.get("seed", DEFAULTS["PLASTOME"]["seed"])
    active_label = p_conf.get("label", DEFAULTS["PLASTOME"]["label"])
    target_size  = p_conf.get("size", DEFAULTS["PLASTOME"]["size"])
    
    final_assembly = "results/assembly_plastome/final_circular.fasta"

elif mode == "MITOME":
    target_ref = config["paths"]["mitome"]["specific"]["fasta"]
    target_gbk = config["paths"]["mitome"]["specific"]["gbk"]
    
    m_conf = config["paths"]["mitome"]
    active_seed  = m_conf.get("seed", DEFAULTS["MITOME"]["seed"])
    active_label = m_conf.get("label", DEFAULTS["MITOME"]["label"])
    target_size  = m_conf.get("size", DEFAULTS["MITOME"]["size"])
    
    final_assembly = "results/assembly_mitome/final_circular.fasta"

else:
    raise ValueError(f"Mode {mode} tidak dikenal.")

# Simpan ke config global agar rule lain bisa akses
config["active_seed"]  = active_seed
config["active_label"] = active_label
config["target_ref"]   = target_ref
config["target_gbk"]   = target_gbk

# 3. SETUP TARGETS
# -----------------------------------------------------------------------------
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)
ASSEMBLERS = [x.strip() for x in config.get("loop_mitohifi", "flye").split(",")]

TARGETS = []
# Visualisasi
TARGETS.extend(expand("results/{sample}/09_viz/{assembler}_alignment_map.html", 
                      sample=samples.index, assembler=ASSEMBLERS))
# Hasil Akhir
TARGETS.extend(expand("results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta", 
                      sample=samples.index))
# Visualisasi Akhir
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