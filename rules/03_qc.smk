# PHASE 6 & 7: JUDGING & REPORTING

rule select_best_candidate:
    input:
        flye  = "results/{sample}/04_assemblies/flye/assembly.fasta",
        raven = "results/{sample}/05_polished/raven_polished.fasta",
        canu  = "results/{sample}/05_polished/canu_polished.fasta",
        ref   = get_target_ref # Primary Ref Only
    output:
        best   = "results/{sample}/07_best_candidate/best_assembly.fasta",
        report = "results/{sample}/07_best_candidate/selection_report.txt"
    log: "logs/{sample}/06_selection.log"
    params:
        min_len = get_min_len,
        max_len = get_max_len
    conda: "envs/blast_biopython.yaml"
    script: "scripts/assess_assemblies.py"

rule plot_final_graphs:
    input:
        flye_gfa  = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        raven_gfa = "results/{sample}/04_assemblies/raven/assembly_graph.gfa",
        canu_gfa  = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    output:
        flye_png  = "results/{sample}/07_best_candidate/viz/flye_final.png",
        raven_png = "results/{sample}/07_best_candidate/viz/raven_final.png",
        canu_png  = "results/{sample}/07_best_candidate/viz/canu_final.png"
    log: "logs/{sample}/07_plot_final.log"
    conda: "envs/visualization.yaml"
    params: 
        opts = "--height 1000 --width 1000 --nodewidth 20 --fontsize 24 --names --lengths --colour depth"
    shell:
        """
        (
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then touch {output.flye_png}; fi
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then touch {output.raven_png}; fi
        if ! Bandage image {input.canu_gfa} {output.canu_png} {params.opts} 2>/dev/null; then touch {output.canu_png}; fi
        ) 2> {log}
        """

rule plot_dotplot:
    input:
        query = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref   = get_target_ref
    output:
        plot = "results/{sample}/07_best_candidate/qc/dotplot.png"
    params:
        prefix = "results/{sample}/07_best_candidate/qc/nucmer"
    conda: "envs/qc_plot.yaml"
    shell:
        """
        if [ -s {input.query} ]; then
            nucmer --prefix={params.prefix} {input.ref} {input.query}
            mummerplot --png --fat --layout --filter -p {params.prefix} {params.prefix}.delta
            mv {params.prefix}.png {output.plot}
        else
            touch {output.plot}
        fi
        """

rule plot_coverage:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        ref   = "results/{sample}/07_best_candidate/best_assembly.fasta"
    output:
        bam   = "results/{sample}/07_best_candidate/qc/mapping.bam",
        depth = "results/{sample}/07_best_candidate/qc/coverage.txt",
        plot  = "results/{sample}/07_best_candidate/qc/coverage_plot.png"
    threads: config["threads"]
    conda: "envs/qc_plot.yaml"
    script: "scripts/plot_coverage.py"

rule generate_final_report:
    input:
        report = "results/{sample}/07_best_candidate/selection_report.txt",
        plots  = ["results/{sample}/07_best_candidate/viz/flye_final.png",
                  "results/{sample}/07_best_candidate/viz/raven_final.png",
                  "results/{sample}/07_best_candidate/viz/canu_final.png"],
        dotplot = "results/{sample}/07_best_candidate/qc/dotplot.png",
        covplot = "results/{sample}/07_best_candidate/qc/coverage_plot.png"
    output:
        html   = "results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html"
    script: "scripts/generate_final_report.py"