# PHASE 0 & 1: PRE-PROCESSING & FETCHING

# --- 1. FETCH EXPANDED (Untuk Master Bait) ---
rule fetch_plastome_expanded:
    output: 
        fasta = config["reference_expanded_out"],
        gbk   = config["reference_expanded_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = True  # <--- SAKLAR ON
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

rule fetch_mito_expanded:
    output: 
        fasta = config["mito_reference_expanded_out"],
        gbk   = config["mito_reference_expanded_gbk"]
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"],
        expand_lineage = True  # <--- SAKLAR ON
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# --- 2. FETCH SPECIFIC (Untuk Smart Filter & Assessment) ---
rule fetch_plastome_specific:
    output: 
        fasta = config["reference_out"],
        gbk   = config["reference_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = False # <--- SAKLAR OFF (Hanya Tebu)
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

rule fetch_mito_specific:
    output: 
        fasta = config["mito_reference_out"],
        gbk   = config["mito_reference_gbk"]
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"],
        expand_lineage = False # <--- SAKLAR OFF (Hanya Tebu)
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# --- 3. CREATE MASTER BAIT ---
# Menggabungkan hasil yang EXPANDED
rule create_master_bait:
    input:
        ref_p = config["reference_expanded_out"], # Pakai Expanded
        ref_m = config["mito_reference_expanded_out"] # Pakai Expanded
    output:
        bait = "resources/MASTER_BAIT_combined.fasta"
    shell:
        "cat {input.ref_p} {input.ref_m} > {output.bait}"