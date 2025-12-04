# PHASE 3, 4, 5: RECRUITMENT, ASSEMBLY, POLISHING

rule filter_reads_dual:
    input:
        reads = "resources/{sample}_organelle_concentrated.fastq", 
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        bait_ref  = config["super_bait_path"] # Pake Super Bait!
    output:
        temp_filtered = temp("results/{sample}/03_temp/filtered_neg.fastq"),
        final_reads   = "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    log: "logs/{sample}/03_filter_reads.log"
    threads: config["threads"]
    conda: "envs/minimap2.yaml"
    shell:
        """
        # 1. Filter Negatif (Blacklist - Mito)
        (minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - | samtools fastq - > {output.temp_filtered}) 2> {log}

        # 2. Baiting Target (Loose Sensitivity)
        (minimap2 -ax map-ont -k12 -A1 -B2 -O2 -E1 -s 40 -t {threads} \
         {input.bait_ref} {output.temp_filtered} \
        | samtools view -b -F 4 - | samtools fastq - > {output.final_reads}) 2>> {log}
        """

rule assemble_flye:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/flye/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        info  = "results/{sample}/04_assemblies/flye/assembly_info.txt"
    log: "logs/{sample}/04_flye.log"
    threads: config["threads"]
    params:
        size = config["target_est_size"],
        outdir = "results/{sample}/04_assemblies/flye"
    conda: "envs/flye.yaml"
    shell:
        "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --genome-size {params.size} --meta --iterations 2 2> {log}"

rule assemble_raven_final:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/raven/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/raven/assembly_graph.gfa"
    log: "logs/{sample}/04_raven.log"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell:
        "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta} 2> {log}"

rule assemble_canu:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/canu/assembly.contigs.fasta",
        gfa   = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    log: "logs/{sample}/04_canu.log"
    threads: config["threads"]
    params:
        size = config["target_est_size"],
        dir = "results/{sample}/04_assemblies/canu"
    conda: "envs/canu.yaml"
    shell:
        """
        if canu -p assembly -d {params.dir} genomeSize={params.size} \
             -nanopore {input} maxThreads={threads} useGrid=false \
             stopOnLowCoverage=5 minInputCoverage=5 > {log} 2>&1; then
            if [ -s {output.fasta} ]; then
                if [ ! -f {output.gfa} ]; then
                    awk '/^>/ {{printf("\\nS\\t%s\\t",substr($1,2))}} !/^>/ {{printf("%s",$0)}} END {{printf("\\n")}}' {output.fasta} | sed '/^$/d' > {output.gfa}
                fi
            else
                touch {output.fasta} {output.gfa}
            fi
        else
            touch {output.fasta} {output.gfa}
        fi
        """

rule polish_raven:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        draft = "results/{sample}/04_assemblies/raven/assembly_graph.gfa"
    output: "results/{sample}/05_polished/raven_polished.fasta"
    log: "logs/{sample}/05_polish_raven.log"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        """
        if minipolish -t {threads} {input.reads} {input.draft} > {output}.temp_gfa 2> {log}; then
            awk '/^S/{{print ">"$2"\\n"$3}}' {output}.temp_gfa | fold > {output}
            rm {output}.temp_gfa
        else
            awk '/^S/{{print ">"$2"\\n"$3}}' {input.draft} | fold > {output}
        fi
        """

rule polish_canu:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        draft = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    output: "results/{sample}/05_polished/canu_polished.fasta"
    log: "logs/{sample}/05_polish_canu.log"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        """
        if minipolish -t {threads} {input.reads} {input.draft} > {output}.temp_gfa 2> {log}; then
            awk '/^S/{{print ">"$2"\\n"$3}}' {output}.temp_gfa | fold > {output}
            rm {output}.temp_gfa
        else
            awk '/^S/{{print ">"$2"\\n"$3}}' {input.draft} | fold > {output}
        fi
        """