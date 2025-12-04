# PHASE 2: ROUGH ASSEMBLY & FILTERING

rule rough_assembly_flye:
    input:
        reads = "resources/{sample}_organelle_concentrated.fastq" 
    output:
        fasta = "results/{sample}/01_rough/flye/assembly.fasta",
        gfa   = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    log: "logs/{sample}/01_rough_flye.log"
    params: outdir = "results/{sample}/01_rough/flye"
    threads: config["threads"]
    conda: "envs/flye.yaml"
    shell: "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --iterations 0 2> {log}"

rule rough_assembly_raven:
    input:
        reads = "resources/{sample}_organelle_concentrated.fastq"
    output:
        fasta = "results/{sample}/01_rough/raven.fasta",
        gfa   = "results/{sample}/01_rough/raven.gfa"
    log: "logs/{sample}/01_rough_raven.log"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell: "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta} 2> {log}"

rule plot_rough_graphs:
    input:
        raven_gfa = "results/{sample}/01_rough/raven.gfa",
        flye_gfa  = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    output:
        raven_png = "results/{sample}/01_rough/viz/raven_graph.png",
        flye_png  = "results/{sample}/01_rough/viz/flye_graph.png"
    log: "logs/{sample}/01_plot_rough.log"
    conda: "envs/visualization.yaml"
    params: 
        opts = "--height 1000 --width 1000 --nodewidth 15 --fontsize 20 --names --lengths --colour depth"
    shell:
        """
        (
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts}; then touch {output.raven_png}; fi
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts}; then touch {output.flye_png}; fi
        ) 2> {log}
        """

rule create_smart_blacklist:
    input:
        fastas = ["results/{sample}/01_rough/raven.fasta", 
                  "results/{sample}/01_rough/flye/assembly.fasta"],
        gfas   = ["results/{sample}/01_rough/raven.gfa", 
                  "results/{sample}/01_rough/flye/assembly_graph.gfa"],
        ref_p = config["reference_out"],
        ref_m = config["mito_reference_out"]
    output:
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        log       = "results/{sample}/02_blacklist/filter_logic.log",
        filtered_gfas = ["results/{sample}/02_blacklist/raven_filtered.gfa",
                         "results/{sample}/02_blacklist/flye_filtered.gfa"]
    log: "logs/{sample}/02_smart_filter.log"
    conda: "envs/blast_biopython.yaml"
    script: "scripts/smart_filter.py"

rule plot_filtered_graphs:
    input:
        raven_gfa = "results/{sample}/02_blacklist/raven_filtered.gfa",
        flye_gfa  = "results/{sample}/02_blacklist/flye_filtered.gfa"
    output:
        raven_png = "results/{sample}/02_blacklist/viz/raven_filtered.png",
        flye_png  = "results/{sample}/02_blacklist/viz/flye_filtered.png"
    log: "logs/{sample}/02_plot_filtered.log"
    conda: "envs/visualization.yaml"
    params: opts = "--height 800 --width 800 --nodewidth 15 --fontsize 20 --names --lengths --colour depth"
    shell:
        """
        (
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then touch {output.raven_png}; fi
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then touch {output.flye_png}; fi
        ) 2> {log}
        """

rule generate_filter_report:
    input:
        log   = "results/{sample}/02_blacklist/filter_logic.log",
        plots = ["results/{sample}/02_blacklist/viz/raven_filtered.png",
                 "results/{sample}/02_blacklist/viz/flye_filtered.png"]
    output:
        html  = "results/{sample}/02_blacklist/filtering_report.html"
    script: "scripts/generate_report.py"