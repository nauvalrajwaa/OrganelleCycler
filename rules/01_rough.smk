# ==============================================================================
# PHASE 2: PRE-QC ROUGH ASSEMBLY
# ==============================================================================

# --- 0. PREPARATION: Extract & Prepare Cleaner Bait ---

# RULE BARU: Ekstrak Gen (CDS/tRNA/rRNA) dari Target Reference GenBank (.gb)
# Ini memastikan kita punya database gen spesifik dari spesies target.
rule extract_target_genes:
    input:
        gbk = config["target_gbk"]  # File .gb yang didownload di prep
    output:
        genes_fasta = "resources/target_extracted_genes.fasta"
    conda: "../envs/python_utils.yaml"  # Pastikan env ini ada biopython
    shell:
        "python scripts/extract_genes.py {input.gbk} {output.genes_fasta}"

# UPDATE RULE: Membuat 'Super Bait' untuk Cleaner
# Gabungan: Label DB (Gen Umum) + Target Extracted Genes (Gen Spesifik)
rule create_cleaner_bait:
    input:
        # 1. Gene Database Umum (misal: embplant_pt.fasta)
        label_db = os.path.join(config["label_db_dir"], config["active_label"]),
        
        # 2. Gene Spesifik dari Target (Hasil ekstraksi di atas)
        target_genes = "resources/target_extracted_genes.fasta"
    output:
        bait_db = "resources/CLEANER_BAIT_TARGET.fasta"
    shell:
        """
        cat {input.label_db} {input.target_genes} > {output.bait_db}
        """

# --- 1. ASSEMBLERS ---

rule rough_assembly_flye:
    input:
        reads = "resources/{sample}_organelle_concentrated.fastq" 
    output:
        fasta = "results/{sample}/01_rough/flye/assembly.fasta",
        gfa   = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    log: "logs/{sample}/01_rough_flye.log"
    params: outdir = "results/{sample}/01_rough/flye"
    threads: config["threads"]
    conda: "../envs/flye.yaml"
    shell: 
        "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --iterations 0 2> {log}"

rule rough_assembly_raven:
    input:
        reads = "resources/{sample}_organelle_concentrated.fastq"
    output:
        fasta = "results/{sample}/01_rough/raven.fasta",
        gfa   = "results/{sample}/01_rough/raven.gfa"
    log: "logs/{sample}/01_rough_raven.log"
    threads: config["threads"]
    conda: "../envs/raven.yaml"
    shell: 
        "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta} 2> {log}"

# ==============================================================================
# PIPELINE: FLYE (Clean -> Solve -> Polish Medaka)
# ==============================================================================

rule clean_rough_flye:
    input:
        gfa = "results/{sample}/01_rough/flye/assembly_graph.gfa",
        # UPDATE: Menggunakan Bait Gabungan (Gene + Genome) yang baru dibuat
        ref = "resources/CLEANER_BAIT_TARGET.fasta"
    output:
        gfa_clean = "results/{sample}/01_rough/flye/assembly_graph.cleaned.gfa"
    params:
        temp_dir = "results/{sample}/01_rough/flye/temp_cleaner"
    log: "logs/{sample}/01_rough_clean_flye.log"
    conda: "../envs/python_utils.yaml"
    shell:
        """
        # Hapus temp dir jika ada
        rm -rf {params.temp_dir}

        python scripts/cleaner.py \
            --gfa {input.gfa} \
            --ref {input.ref} \
            --out {output.gfa_clean} \
            --temp {params.temp_dir} > {log} 2>&1
        """

rule solve_rough_flye:
    input:
        gfa_clean = "results/{sample}/01_rough/flye/assembly_graph.cleaned.gfa"
    output:
        fasta_circ = "results/{sample}/01_rough/flye/assembly_circular.fasta"
    log: "logs/{sample}/01_rough_solve_flye.log"
    conda: "../envs/python_utils.yaml"
    shell:
        "python scripts/solver.py {input.gfa_clean} {output.fasta_circ} > {log} 2>&1"

rule polish_rough_flye:
    input:
        # Polishing dilakukan pada FASTA hasil Solver
        assembly = "results/{sample}/01_rough/flye/assembly_circular.fasta",
        reads    = "resources/{sample}_organelle_concentrated.fastq"
    output:
        polished = "results/{sample}/01_rough/flye/assembly_circular_polished.fasta"
    log: "logs/{sample}/01_rough_polish_flye.log"
    threads: config["threads"]
    conda: "../envs/medaka.yaml"
    params:
        model = config.get("medaka_model", "r10.4.1_e8.2_400bps_sup@v5.2.0"),
        outdir = "results/{sample}/01_rough/flye/medaka_temp"
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            touch {output.polished}
            exit 0
        fi

        # Medaka output ke folder, kita pindahkan hasilnya
        medaka_consensus -i {input.reads} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model} > {log} 2>&1
        
        # Pindahkan consensus ke output file
        mv {params.outdir}/consensus.fasta {output.polished}
        
        # Bersihkan folder temp medaka
        rm -rf {params.outdir}
        """

# ==============================================================================
# PIPELINE: RAVEN (Clean -> Solve -> Polish Medaka)
# ==============================================================================

rule clean_rough_raven:
    input:
        gfa = "results/{sample}/01_rough/raven.gfa",
        # UPDATE: Menggunakan Bait Gabungan (Gene + Genome)
        ref = "resources/CLEANER_BAIT_TARGET.fasta"
    output:
        gfa_clean = "results/{sample}/01_rough/raven.cleaned.gfa"
    params:
        temp_dir = "results/{sample}/01_rough/raven_temp_cleaner"
    log: "logs/{sample}/01_rough_clean_raven.log"
    conda: "../envs/python_utils.yaml"
    shell:
        """
        rm -rf {params.temp_dir}

        python scripts/cleaner.py \
            --gfa {input.gfa} \
            --ref {input.ref} \
            --out {output.gfa_clean} \
            --temp {params.temp_dir} > {log} 2>&1
        """

rule solve_rough_raven:
    input:
        gfa_clean = "results/{sample}/01_rough/raven.cleaned.gfa"
    output:
        fasta_circ = "results/{sample}/01_rough/raven_circular.fasta"
    log: "logs/{sample}/01_rough_solve_raven.log"
    conda: "../envs/python_utils.yaml"
    shell:
        "python scripts/solver.py {input.gfa_clean} {output.fasta_circ} > {log} 2>&1"

rule polish_rough_raven:
    input:
        assembly = "results/{sample}/01_rough/raven_circular.fasta",
        reads    = "resources/{sample}_organelle_concentrated.fastq"
    output:
        polished = "results/{sample}/01_rough/raven_circular_polished.fasta"
    log: "logs/{sample}/01_rough_polish_raven.log"
    threads: config["threads"]
    conda: "../envs/medaka.yaml"
    params:
        model = config.get("medaka_model", "r10.4.1_e8.2_400bps_sup@v5.2.0"),
        outdir = "results/{sample}/01_rough/raven/medaka_temp"
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            touch {output.polished}
            exit 0
        fi

        medaka_consensus -i {input.reads} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model} > {log} 2>&1
        
        mv {params.outdir}/consensus.fasta {output.polished}
        rm -rf {params.outdir}
        """

# ==============================================================================
# VISUALIZATION (RAW GFA vs CLEANED GFA)
# ==============================================================================

rule plot_rough_graphs:
    input:
        raven_gfa   = "results/{sample}/01_rough/raven.gfa",
        flye_gfa    = "results/{sample}/01_rough/flye/assembly_graph.gfa",
        raven_clean = "results/{sample}/01_rough/raven.cleaned.gfa",
        flye_clean  = "results/{sample}/01_rough/flye/assembly_graph.cleaned.gfa"
    output:
        raven_png       = "results/{sample}/01_rough/viz/raven_graph.png",
        flye_png        = "results/{sample}/01_rough/viz/flye_graph.png",
        raven_clean_png = "results/{sample}/01_rough/viz/raven_graph_cleaned.png",
        flye_clean_png  = "results/{sample}/01_rough/viz/flye_graph_cleaned.png"
    log: "logs/{sample}/01_plot_rough.log"
    conda: "../envs/visualization.yaml"
    params: 
        opts = "--height 1000 --width 1000 --nodewidth 15 --fontsize 20 --names --lengths --colour depth"
    shell:
        """
        (
        # Raw Plots
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts}; then touch {output.raven_png}; fi
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts}; then touch {output.flye_png}; fi

        # Cleaned Plots (Cek file ada isinya tidak)
        if [ -s {input.raven_clean} ]; then
            Bandage image {input.raven_clean} {output.raven_clean_png} {params.opts}
        else
            touch {output.raven_clean_png}
        fi

        if [ -s {input.flye_clean} ]; then
            Bandage image {input.flye_clean} {output.flye_clean_png} {params.opts}
        else
            touch {output.flye_clean_png}
        fi
        ) 2> {log}
        """