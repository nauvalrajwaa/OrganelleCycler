# =============================================================================
# PHASE 3, 4, 5: RECRUITMENT, ASSEMBLY, POLISHING
# File: rules/02_assembly.smk
# Strategy: 
# 1. Filter Reads (Buang Blacklist, Ambil Target)
# 2. Assembly (Flye/Raven) pakai Filtered Reads
# 3. Clean -> Solve -> Polish (Medaka)
# =============================================================================

# --- STEP 1: FILTERING DENGAN STATISTIK ---

rule filter_reads_dual:
    input:
        reads     = "resources/{sample}_organelle_concentrated.fastq",
        blacklist = config["blacklist_ref"], 
        target    = config["target_ref"]
    output:
        temp_neg    = temp("results/{sample}/03_temp/filtered_neg.fastq"),
        final_reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        stats       = "results/{sample}/03_filtered_reads/filtering_stats.txt"
    log: "logs/{sample}/03_filter_reads.log"
    threads: config["threads"]
    conda: "../envs/minimap2.yaml"
    shell:
        """
        # 1. Hitung Total Awal (Baris / 4 = Jumlah Reads)
        TOTAL_LINES=$(wc -l < {input.reads})
        TOTAL_READS=$((TOTAL_LINES / 4))

        # 2. Blacklist Filtering (Buang yang nempel Blacklist)
        # -f 4 : Keep Unmapped (artinya yang TIDAK nempel ke blacklist disimpan)
        (minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - \
        | samtools fastq - > {output.temp_neg}) 2> {log}
        
        # Hitung sisa setelah blacklist
        NEG_LINES=$(wc -l < {output.temp_neg})
        NEG_READS=$((NEG_LINES / 4))

        # 3. Target Recruitment (Ambil yang nempel Target)
        # -F 4 : Discard Unmapped (artinya HANYA yang nempel ke target disimpan)
        (minimap2 -ax map-ont -t {threads} {input.target} {output.temp_neg} \
        | samtools view -b -F 4 - \
        | samtools fastq - > {output.final_reads}) 2>> {log}
        
        # Hitung Final
        FINAL_LINES=$(wc -l < {output.final_reads})
        FINAL_READS=$((FINAL_LINES / 4))

        # Simpan Statistik ke File
        echo "Total_Raw_Reads: $TOTAL_READS" > {output.stats}
        echo "After_Blacklist: $NEG_READS" >> {output.stats}
        echo "Final_Recruited: $FINAL_READS" >> {output.stats}
        """

# --- STEP 1.5: GENERATE FILTERING REPORT (HTML) ---

rule generate_filtering_report:
    input:
        stats = "results/{sample}/03_filtered_reads/filtering_stats.txt"
    output:
        html = "results/{sample}/03_filtered_reads/FILTERING_REPORT.html"
    conda: "../envs/python_utils.yaml" # Pastikan env ini punya pandas/jinja2 jika script butuh
    script: "../scripts/report_filtering.py"

# --- STEP 2: ASSEMBLY (Input: Filtered Reads) ---

rule final_assembly_flye:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        gfa   = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        fasta = "results/{sample}/04_assemblies/flye/assembly.fasta"
    log: "logs/{sample}/04_flye_final.log"
    params: outdir = "results/{sample}/04_assemblies/flye"
    threads: config["threads"]
    conda: "../envs/flye.yaml"
    shell:
        # Gunakan --meta agar Flye tidak membuang plasmid/coverage aneh
        "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --meta --iterations 0 2> {log}"

rule final_assembly_raven:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        gfa   = "results/{sample}/04_assemblies/raven/assembly_graph.gfa",
        fasta = "results/{sample}/04_assemblies/raven/assembly.fasta"
    log: "logs/{sample}/04_raven_final.log"
    threads: config["threads"]
    conda: "../envs/raven.yaml"
    shell:
        "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta} 2> {log}"

# --- STEP 3: CLEANING & SOLVING (Sama seperti Rules 01) ---

# --- FLYE PATH ---
rule clean_final_flye:
    input:
        gfa = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        # UPDATE: Gunakan 'CLEANER_BAIT_TARGET' (Gene + Genome) agar konsisten dengan 01_rough
        ref = "resources/CLEANER_BAIT_TARGET.fasta"
    output:
        gfa_clean = "results/{sample}/04_assemblies/flye/assembly_graph.cleaned.gfa"
    params: temp_dir = "results/{sample}/04_assemblies/flye/temp_cleaner"
    log: "logs/{sample}/04_clean_flye.log"
    conda: "../envs/python_utils.yaml"
    shell:
        """
        rm -rf {params.temp_dir}
        python scripts/cleaner.py --gfa {input.gfa} --ref {input.ref} --out {output.gfa_clean} --temp {params.temp_dir} > {log} 2>&1
        """

rule solve_final_flye:
    input: "results/{sample}/04_assemblies/flye/assembly_graph.cleaned.gfa"
    output: "results/{sample}/04_assemblies/flye/assembly_circular.fasta"
    log: "logs/{sample}/04_solve_flye.log"
    conda: "../envs/python_utils.yaml"
    shell:
        "python scripts/solver.py {input} {output} > {log} 2>&1"

# --- RAVEN PATH ---
rule clean_final_raven:
    input:
        gfa = "results/{sample}/04_assemblies/raven/assembly_graph.gfa",
        # UPDATE: Gunakan 'CLEANER_BAIT_TARGET' (Gene + Genome)
        ref = "resources/CLEANER_BAIT_TARGET.fasta"
    output:
        gfa_clean = "results/{sample}/04_assemblies/raven/assembly_graph.cleaned.gfa"
    params: temp_dir = "results/{sample}/04_assemblies/raven/temp_cleaner"
    log: "logs/{sample}/04_clean_raven.log"
    conda: "../envs/python_utils.yaml"
    shell:
        """
        rm -rf {params.temp_dir}
        python scripts/cleaner.py --gfa {input.gfa} --ref {input.ref} --out {output.gfa_clean} --temp {params.temp_dir} > {log} 2>&1
        """

rule solve_final_raven:
    input: "results/{sample}/04_assemblies/raven/assembly_graph.cleaned.gfa"
    output: "results/{sample}/04_assemblies/raven/assembly_circular.fasta"
    log: "logs/{sample}/04_solve_raven.log"
    conda: "../envs/python_utils.yaml"
    shell:
        "python scripts/solver.py {input} {output} > {log} 2>&1"

# --- STEP 4: FINAL POLISHING (Medaka) ---
# Menggunakan Reads yang SUDAH DIFILTER untuk akurasi maksimal

rule polish_final_flye:
    input:
        assembly = "results/{sample}/04_assemblies/flye/assembly_circular.fasta",
        reads    = "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        final = "results/{sample}/05_final_assembly/flye_polished.fasta"
    log: "logs/{sample}/05_polish_flye.log"
    threads: config["threads"]
    conda: "../envs/medaka.yaml"
    params:
        model = config.get("medaka_model", "r10.4.1_e8.2_400bps_sup@v5.2.0"),
        outdir = "results/{sample}/05_final_assembly/flye_medaka_temp"
    shell:
        """
        if [ ! -s {input.assembly} ]; then touch {output.final}; exit 0; fi
        
        medaka_consensus -i {input.reads} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model} > {log} 2>&1
        mv {params.outdir}/consensus.fasta {output.final}
        rm -rf {params.outdir}
        """

rule polish_final_raven:
    input:
        assembly = "results/{sample}/04_assemblies/raven/assembly_circular.fasta",
        reads    = "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        final = "results/{sample}/05_final_assembly/raven_polished.fasta"
    log: "logs/{sample}/05_polish_raven.log"
    threads: config["threads"]
    conda: "../envs/medaka.yaml"
    params:
        model = config.get("medaka_model", "r10.4.1_e8.2_400bps_sup@v5.2.0"),
        outdir = "results/{sample}/05_final_assembly/raven_medaka_temp"
    shell:
        """
        if [ ! -s {input.assembly} ]; then touch {output.final}; exit 0; fi
        
        medaka_consensus -i {input.reads} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model} > {log} 2>&1
        mv {params.outdir}/consensus.fasta {output.final}
        rm -rf {params.outdir}
        """

# --- STEP 5: VISUALIZATION (Comparison) ---

rule plot_final_graphs:
    input:
        gfa_flye = "results/{sample}/04_assemblies/flye/assembly_graph.cleaned.gfa",
        gfa_raven = "results/{sample}/04_assemblies/raven/assembly_graph.cleaned.gfa"
    output:
        png_flye = "results/{sample}/05_final_assembly/viz/flye_cleaned_graph.png",
        png_raven = "results/{sample}/05_final_assembly/viz/raven_cleaned_graph.png"
    log: "logs/{sample}/05_plot_final.log"
    conda: "../envs/visualization.yaml"
    params: opts = "--height 1000 --width 1000 --nodewidth 15 --fontsize 20 --names --lengths --colour depth"
    shell:
        """
        (
        if [ -s {input.gfa_flye} ]; then Bandage image {input.gfa_flye} {output.png_flye} {params.opts}; else touch {output.png_flye}; fi
        if [ -s {input.gfa_raven} ]; then Bandage image {input.gfa_raven} {output.png_raven} {params.opts}; else touch {output.png_raven}; fi
        ) 2> {log}
        """