# ==============================================================================
# PHASE 1: PREPARATION & FETCHING
# File: rules/00_prep.smk
# ==============================================================================

# --- Helper Python: Ambil Nama Spesies & Genus ---
SPECIES = config["project"]["species_name"]       # Contoh: "Saccharum officinarum"
GENUS   = SPECIES.split()[0]                      # Ambil kata pertama: "Saccharum"

# ------------------------------------------------------------------------------
# BAGIAN 1: FETCH SPECIFIC (Selalu Jalan)
# ------------------------------------------------------------------------------

rule fetch_plastome_specific:
    output:
        fasta = config["paths"]["plastome"]["specific"]["fasta"],
        gbk   = config["paths"]["plastome"]["specific"]["gbk"]
    params:
        # OTOMATIS: "Saccharum officinarum[orgn] AND chloroplast[filter]..."
        search_term = f'"{SPECIES}"[orgn] AND chloroplast[filter] AND complete genome[title]',
        min_len     = config["settings"]["filter_params"]["p_min"],
        max_len     = config["settings"]["filter_params"]["p_max"],
        expand_lineage = False 
    resources: ncbi_connection=1
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

rule fetch_mito_specific:
    output:
        fasta = config["paths"]["mitome"]["specific"]["fasta"],
        gbk   = config["paths"]["mitome"]["specific"]["gbk"]
    params:
        # OTOMATIS: "Saccharum officinarum[orgn] AND mitochondrion[filter]..."
        search_term = f'"{SPECIES}"[orgn] AND mitochondrion[filter] AND complete genome[title]',
        min_len     = config["settings"]["filter_params"]["m_min"],
        max_len     = config["settings"]["filter_params"]["m_max"],
        expand_lineage = False
    resources: ncbi_connection=1
    conda: "../envs/blast_biopython.yaml"
    script: "../scripts/fetch_organelle_ref.py"

# ------------------------------------------------------------------------------
# BAGIAN 2: FETCH EXPANDED (Kondisional)
# ------------------------------------------------------------------------------

if config["settings"]["use_expanded_search"]:

    rule fetch_plastome_expanded:
        output:
            fasta = config["paths"]["plastome"]["expanded"]["fasta"],
            gbk   = config["paths"]["plastome"]["expanded"]["gbk"]
        params:
            # OTOMATIS: Pakai GENUS saja -> "Saccharum[orgn] AND ..."
            search_term = f'"{GENUS}"[orgn] AND chloroplast[filter] AND complete genome[title]',
            min_len     = config["settings"]["filter_params"]["p_min"],
            max_len     = config["settings"]["filter_params"]["p_max"],
            expand_lineage = True
        resources: ncbi_connection=1
        conda: "../envs/blast_biopython.yaml"
        script: "../scripts/fetch_organelle_ref.py"

    rule fetch_mito_expanded:
        output:
            fasta = config["paths"]["mitome"]["expanded"]["fasta"],
            gbk   = config["paths"]["mitome"]["expanded"]["gbk"]
        params:
            # OTOMATIS: Pakai GENUS saja -> "Saccharum[orgn] AND ..."
            search_term = f'"{GENUS}"[orgn] AND mitochondrion[filter] AND complete genome[title]',
            min_len     = config["settings"]["filter_params"]["m_min"],
            max_len     = config["settings"]["filter_params"]["m_max"],
            expand_lineage = True
        resources: ncbi_connection=1
        conda: "../envs/blast_biopython.yaml"
        script: "../scripts/fetch_organelle_ref.py"

# --- BAGIAN 2: HYBRID RECRUITMENT (ADAPTASI GETORGANELLE) ---
import os

# Rule 2.A. Create Hybrid Bait (Gabungan Seed Database + Target + Blacklist)
# Tujuannya: Membuat pancingan "Super" yang menangkap semua DNA organel (baik target maupun kontaminan)
# agar tidak ada yang lolos di tahap awal. Pemisahan dilakukan nanti.

rule create_hybrid_bait:
    input:
        # 1. Seed Database (Seed GetOrganelle yang aktif)
        # Variabel 'active_seed' sudah diset otomatis di Snakefile (embplant_pt atau embplant_mt)
        seed_db = os.path.join(config["seed_db_dir"], config["active_seed"]),
        
        # 2. Reference Target (HANYA Target, bukan keduanya)
        # Variabel 'target_ref' sudah diset otomatis di Snakefile sesuai MODE
        target_ref = config["target_ref"]
    output:
        hybrid_bait = "resources/HYBRID_BAIT.fasta"
    shell:
        """
        # Gabungkan Seed Database + Target Reference Genome
        # Ini membuat baiting fokus hanya pada organel yang sedang dikerjakan.
        cat {input.seed_db} {input.target_ref} > {output.hybrid_bait}
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