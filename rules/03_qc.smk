# ==============================================================================
# PHASE 6: SELECTION & QC
# File: rules/03_qc.smk
# Tujuan: Memilih assembly terbaik antara Flye vs Raven & QC mendalam
# ==============================================================================

# --- 1. SELEKSI KANDIDAT TERBAIK ---
rule select_best_candidate:
    input:
        flye  = "results/{sample}/05_final_assembly/flye_polished.fasta",
        raven = "results/{sample}/05_final_assembly/raven_polished.fasta",
        ref   = config["target_ref"]
    output:
        best   = "results/{sample}/07_best_candidate/best_candidate.fasta",
        report = "results/{sample}/07_best_candidate/selection_report.txt"
    log: "logs/{sample}/06_selection.log"
    params:
        # Mengambil range ukuran yang diharapkan dari config (misal: 150000 +- 20%)
        target_size = config["target_est_size"], 
        tolerance   = 0.2 
    conda: "../envs/python_utils.yaml" # Butuh Biopython
    script: "../scripts/select_best.py"

# --- 2. QUAST QC (Kualitas Assembly vs Referensi) ---
rule run_quast:
    input:
        assembly = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref      = config["target_ref"]
    output:
        report_txt = "results/{sample}/07_best_candidate/quast/report.txt",
        report_pdf = "results/{sample}/07_best_candidate/quast/report.pdf"
    log: "logs/{sample}/07_quast.log"
    params:
        outdir = "results/{sample}/07_best_candidate/quast"
    conda: "../envs/quast.yaml"
    shell:
        """
        # Jalankan QUAST jika assembly tidak kosong
        if [ -s {input.assembly} ]; then
            quast.py {input.assembly} -r {input.ref} -o {params.outdir} --threads 1 > {log} 2>&1
        else
            touch {output.report_txt} {output.report_pdf}
        fi
        """

# --- 3. VISUALISASI DOTPLOT (MUMmer) ---
rule plot_dotplot:
    input:
        query = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref   = config["target_ref"]
    output:
        plot = "results/{sample}/07_best_candidate/viz/dotplot.png"
    params:
        prefix = "results/{sample}/07_best_candidate/viz/nucmer"
    log: "logs/{sample}/07_dotplot.log"
    conda: "../envs/qc_plot.yaml"
    shell:
        """
        if [ -s {input.query} ]; then
            nucmer --prefix={params.prefix} {input.ref} {input.query} 2> {log}
            # --fat = force alignment even if too many matches
            # --layout = arrange layout for plot
            # --filter = filter for best hits
            mummerplot --png --fat --layout --filter -p {params.prefix} {params.prefix}.delta 2>> {log}
            mv {params.prefix}.png {output.plot}
        else
            touch {output.plot}
        fi
        """

# --- 4. VISUALISASI GRAPH (Cleaned GFA) ---
# Mengambil graph sebelum dipoles (karena polisher outputnya FASTA)
rule plot_final_graphs_qc:
    input:
        flye_gfa  = "results/{sample}/04_assemblies/flye/assembly_graph.cleaned.gfa",
        raven_gfa = "results/{sample}/04_assemblies/raven/assembly_graph.cleaned.gfa"
    output:
        flye_png  = "results/{sample}/07_best_candidate/viz/flye_graph.png",
        raven_png = "results/{sample}/07_best_candidate/viz/raven_graph.png"
    log: "logs/{sample}/07_plot_graphs.log"
    conda: "../envs/visualization.yaml"
    params: 
        opts = "--height 1000 --width 1000 --nodewidth 20 --fontsize 24 --names --lengths --colour depth"
    shell:
        """
        (
        if [ -s {input.flye_gfa} ]; then Bandage image {input.flye_gfa} {output.flye_png} {params.opts}; else touch {output.flye_png}; fi
        if [ -s {input.raven_gfa} ]; then Bandage image {input.raven_gfa} {output.raven_png} {params.opts}; else touch {output.raven_png}; fi
        ) 2> {log}
        """

# --- 5. COVERAGE PLOT (Reads Mapping back to Best Assembly) ---
rule map_reads_final:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        ref   = "results/{sample}/07_best_candidate/best_assembly.fasta"
    output:
        bam = "results/{sample}/07_best_candidate/qc/mapping.bam",
        bai = "results/{sample}/07_best_candidate/qc/mapping.bam.bai"
    log: "logs/{sample}/07_map_final.log"
    threads: config["threads"]
    conda: "../envs/minimap2.yaml"
    shell:
        """
        if [ -s {input.ref} ]; then
            minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
        else
            touch {output.bam} {output.bai}
        fi
        """

rule plot_coverage:
    input:
        bam = "results/{sample}/07_best_candidate/qc/mapping.bam",
        bai = "results/{sample}/07_best_candidate/qc/mapping.bam.bai"
    output:
        plot = "results/{sample}/07_best_candidate/viz/coverage_plot.png",
        stats = "results/{sample}/07_best_candidate/qc/coverage_stats.txt"
    log: "logs/{sample}/07_plot_coverage.log"
    conda: "../envs/python_utils.yaml" # Menggunakan script python matplotlib
    script: "../scripts/plot_coverage.py"

# --- 6. FINAL REPORT ---
rule generate_final_report:
    input:
        sel_report = "results/{sample}/07_best_candidate/selection_report.txt",
        quast_txt  = "results/{sample}/07_best_candidate/quast/report.txt",
        flye_img   = "results/{sample}/07_best_candidate/viz/flye_graph.png",
        raven_img  = "results/{sample}/07_best_candidate/viz/raven_graph.png",
        dotplot    = "results/{sample}/07_best_candidate/viz/dotplot.png",
        covplot    = "results/{sample}/07_best_candidate/viz/coverage_plot.png"
    output:
        html = "results/{sample}/FINAL_REPORT_{sample}.html"
    log: "logs/{sample}/08_report.log"
    script: "../scripts/generate_final_report.py"