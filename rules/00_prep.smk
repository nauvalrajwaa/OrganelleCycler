# =============================================================================
# PHASE 0 & 1: PRE-PROCESSING & FETCHING
# File: rules/00_prep.smk
# =============================================================================

# --- BAGIAN 1: FETCHING REFERENSI (4 Rules) ---

# 1.A. Fetch Plastome EXPANDED (Untuk Master Bait) -> Saklar ON
rule fetch_plastome_expanded:
    output: 
        fasta = config["reference_expanded_out"],
        gbk   = config["reference_expanded_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = True  # <--- SAKLAR ON
    conda: "../envs/blast_biopython.yaml" # Gunakan .. jika file env ada di root/envs
    script: "../scripts/fetch_organelle_ref.py" # Gunakan .. jika script ada di root/scripts

# 1.B. Fetch Mito EXPANDED (Untuk Master Bait) -> Saklar ON
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

# 1.C. Fetch Plastome SPECIFIC (Untuk Smart Filter & Assessment) -> Saklar OFF
rule fetch_plastome_specific:
    output: 
        fasta = config["reference_out"],
        gbk   = config["reference_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = False # <--- SAKLAR OFF
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# 1.D. Fetch Mito SPECIFIC (Untuk Smart Filter) -> Saklar OFF
rule fetch_mito_specific:
    output: 
        fasta = config["mito_reference_out"],
        gbk   = config["mito_reference_gbk"]
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"],
        expand_lineage = False # <--- SAKLAR OFF
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"


# --- BAGIAN 2: PRE-FILTERING (Membuat Bahan Baku) ---

# 2.A. Create Master Bait (Gabungan dari yang EXPANDED)
rule create_master_bait:
    input:
        ref_p = config["reference_expanded_out"],     # Input dari Rule 1.A
        ref_m = config["mito_reference_expanded_out"] # Input dari Rule 1.B
    output:
        bait = "resources/MASTER_BAIT_combined.fasta"
    shell:
        "cat {input.ref_p} {input.ref_m} > {output.bait}"

# 2.B. The Big Sieve (Pukat Harimau 80GB -> 1GB)
rule pre_filter_huge_data:
    input:
        # Ambil path 80GB dari TSV
        raw_reads = lambda wildcards: samples.loc[wildcards.sample, "reads_path"],
        # Ambil Master Bait dari langkah sebelumnya
        bait      = "resources/MASTER_BAIT_combined.fasta"
    output:
        clean_fastq = "resources/{sample}_organelle_concentrated.fastq"
    threads: 64
    conda: "../envs/minimap2.yaml"
    shell:
        """
        # Minimap2 dengan parameter LOOSE (Pukat Harimau)
        minimap2 -ax map-ont \
            -k12 -A1 -B2 -O2 -E1 -s40 \
            -t {threads} \
            {input.bait} {input.raw_reads} \
        | samtools view -@ 4 -b -F 4 - \
        | samtools fastq -@ 4 - > {output.clean_fastq}
        """