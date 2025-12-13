# =============================================================================
# PHASE 0 & 1: PRE-PROCESSING & FETCHING
# File: rules/00_prep.smk
# =============================================================================

# --- BAGIAN 1: FETCHING REFERENSI (4 Rules) ---
# Tujuannya: Mengisi "Gudang" (Storage Paths) yang didefinisikan di config.
# Kita butuh versi EXPANDED (untuk Pukat Harimau) dan SPECIFIC (untuk Target/Blacklist).

# 1.A. Fetch Plastome EXPANDED (Untuk Bahan Pukat Harimau)
rule fetch_plastome_expanded:
    output: 
        fasta = config["plastome_expanded_fasta"], # <-- Update: Storage Path
        gbk   = config["plastome_expanded_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = True   # SAKLAR ON: Cari kerabat luas
    resources: ncbi_connection=1 
    conda: "../envs/blast_biopython.yaml" 
    script: "../scripts/fetch_organelle_ref.py" 

# 1.B. Fetch Mito EXPANDED (Untuk Bahan Pukat Harimau)
rule fetch_mito_expanded:
    output: 
        fasta = config["mitome_expanded_fasta"],   # <-- Update: Storage Path
        gbk   = config["mitome_expanded_gbk"]
    params:
        search_term = config["mito_search_term"], 
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"],
        expand_lineage = True   # SAKLAR ON: Cari kerabat luas
    resources: ncbi_connection=1 
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# 1.C. Fetch Plastome SPECIFIC (Untuk Disimpan di Slot Plastome)
rule fetch_plastome_specific:
    output: 
        fasta = config["plastome_specific_fasta"], # <-- Update: Storage Path
        gbk   = config["plastome_specific_gbk"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"],
        expand_lineage = False  # SAKLAR OFF: Spesifik spesies ini
    resources: ncbi_connection=1 
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# 1.D. Fetch Mito SPECIFIC (Untuk Disimpan di Slot Mitome)
rule fetch_mito_specific:
    output: 
        fasta = config["mitome_specific_fasta"],   # <-- Update: Storage Path
        gbk   = config["mitome_specific_gbk"]
    params:
        search_term = config["mito_search_term"], 
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"],
        expand_lineage = False  # SAKLAR OFF: Spesifik spesies ini
    resources: ncbi_connection=1 
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# --- BAGIAN 2: HYBRID RECRUITMENT (ADAPTASI GETORGANELLE) ---

# Rule 2.A. Create Hybrid Bait (Gabungan Seed Database + Target Spesifik Kita)
rule create_hybrid_bait:
    input:
        # UPDATE: Menggabungkan Path Dir + Nama File dari Config
        seed_db = os.path.join(config["seed_db_dir"], config["active_seed_file"]),
        
        # Referensi spesifik (Expanded Version)
        ref_p   = config["plastome_expanded_fasta"],
        ref_m   = config["mitome_expanded_fasta"]
    output:
        hybrid_bait = "resources/HYBRID_BAIT.fasta"
    shell:
        """
        cat {input.seed_db} {input.ref_p} {input.ref_m} > {output.hybrid_bait}
        """

# 2.B. The Smart Recruiter (Menggantikan Big Sieve Minimap2 Biasa)
# Menggunakan script python 'recruiter.py' untuk:
# 1. Filter Kualitas (Pre-map)
# 2. Mapping ke Hybrid Bait
# 3. Normalisasi Coverage (Subsampling jika > 200x)
rule smart_recruitment:
    input:
        raw_reads   = lambda wildcards: samples.loc[wildcards.sample, "reads_path"],
        hybrid_bait = "resources/HYBRID_BAIT.fasta",
        script      = "scripts/recruiter.py" # Dependency script
    output:
        # Output ini SEKARANG sudah bersih dan ternormalisasi
        concentrated = "resources/{sample}_organelle_concentrated.fastq" 
    params:
        # Mengirim parameter dari config ke script python
        out_dir = "resources/temp_recruitment_{sample}", # Folder temp
        est_size = config["recruitment"]["est_genome_size"],
        min_cov  = config["recruitment"]["min_coverage"],
        max_cov  = config["recruitment"]["max_coverage"],
        min_len  = config["recruitment"]["min_read_len"],
        min_qual = config["recruitment"]["min_read_qual"]
    threads: 32 # Minimap2 di dalam script butuh thread
    conda: "../envs/minimap2.yaml" # Script butuh samtools & minimap2
    shell:
        """
        # Menjalankan script python secara shell command
        # Usage: python recruiter.py <bait> <reads> <out_file> <threads> <est_size> <max_cov> ...
        
        # Kita perlu sedikit modifikasi cara panggil scriptnya agar sesuai dengan argumen sys.argv
        # Asumsi recruiter.py dimodifikasi sedikit untuk menerima Argumen via CLI argparse atau sys.argv urut
        
        python {input.script} \
            --bait {input.hybrid_bait} \
            --reads {input.raw_reads} \
            --output {output.concentrated} \
            --threads {threads} \
            --est_size {params.est_size} \
            --max_cov {params.max_cov} \
            --min_len {params.min_len} \
            --min_qual {params.min_qual}
        
        # Bersihkan folder temp jika ada
        rm -rf {params.out_dir}
        """