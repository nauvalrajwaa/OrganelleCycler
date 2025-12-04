configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# 1. SETUP GLOBAL
# =============================================================================
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)

# =============================================================================
# 2. HELPER FUNCTIONS (Global)
# =============================================================================
def get_target_ref(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_out"]
    else: return config["reference_out"]

def get_min_len(wildcards):
    if config.get("target_assembly") == "MITO": return config.get("mito_min_len", 200000)
    else: return config.get("plastome_min", 80000)

def get_max_len(wildcards):
    if config.get("target_assembly") == "MITO": return config.get("mito_max_len", 2000000)
    else: return config.get("plastome_max", 180000)

def get_rescue_ref_gbk(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_gbk"]
    else: return config["reference_gbk"]

def get_rescue_ref(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_out"]
    else: return config["reference_out"]

# =============================================================================
# 3. TARGET DEFINITION
# =============================================================================
TARGETS = []

if config.get("target_assembly") in ["PLASTOME", "MITO", "BOTH"]:
    TARGETS.extend(expand("results/{sample}/02_blacklist/filtering_report.html", sample=samples.index))
    TARGETS.extend(expand("results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html", sample=samples.index))

    if config.get("run_circularization_rescue") == "YES":
        TARGETS.extend(expand("results/{sample}/08_rescue/final_circularized.fasta", sample=samples.index))
        TARGETS.extend(expand("results/{sample}/08_rescue/final_mitogenome.gb", sample=samples.index))
    else:
        TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))

rule all:
    input: TARGETS

# =============================================================================
# 4. INCLUDE MODULES
# =============================================================================
include: "rules/00_prep.smk"
include: "rules/01_rough.smk"
include: "rules/02_assembly.smk"
include: "rules/03_qc.smk"
include: "rules/04_rescue.smk"