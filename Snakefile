configfile: "config/config.yaml"
import pandas as pd

# 1. Baca daftar sampel
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)

# 2. Target Akhir (All Rule)
# Pipeline selesai jika semua sampel sudah punya 'best_assembly.fasta'
rule all:
    input:
        expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index)

# =============================================================================
# PHASE 1: PREPARATION (Fetch References)
# =============================================================================

# Download Referensi Plastome (Otomatis via Python)
rule fetch_plastome_ref:
    output:
        config["reference_out"]
    conda:
        "envs/blast_biopython.yaml"
    script:
        "scripts/fetch_refs.py"

# Download Referensi Mitokondria (Otomatis via Python)
rule fetch_mito_ref:
    output:
        config["mito_reference_out"]
    conda:
        "envs/blast_biopython.yaml"
    script:
        "scripts/fetch_mito_refs.py"

# =============================================================================
# PHASE 2: ROUGH ASSEMBLY & SMART IDENTIFICATION
# =============================================================================

# Step 1: Raven Rough Assembly (Untuk dapat Graph & Contig Besar)
rule rough_assembly_raven:
    input:
        reads = lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/assembly.fasta",
        gfa   = "results/{sample}/01_rough/assembly_graph.gfa"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell:
        """
        raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input.reads} > {output.fasta}
        """

# Step 2: Smart Filtering (Python Script: Graph + Identity Check)
rule create_smart_blacklist:
    input:
        fasta = "results/{sample}/01_rough/assembly.fasta",
        gfa   = "results/{sample}/01_rough/assembly_graph.gfa",
        ref_p = config["reference_out"],      # Dependency ke Phase 1
        ref_m = config["mito_reference_out"]  # Dependency ke Phase 1
    output:
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        log       = "results/{sample}/02_blacklist/filter_logic.log"
    conda: "envs/blast_biopython.yaml"
    script:
        "scripts/smart_filter.py"

# =============================================================================
# PHASE 3: DUAL FILTERING (Cleaning Reads)
# =============================================================================

# Step 3: Negative Filter (Buang Mito) -> Positive Filter (Baiting Plastome)
rule filter_reads_dual:
    input:
        reads = lambda wildcards: samples.loc[wildcards.sample, "reads_path"],
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        ref_p = config["reference_out"]
    output:
        temp_no_mito = temp("results/{sample}/03_temp/no_mito.fastq"),
        final_reads  = "results/{sample}/03_filtered_reads/recruited_plastome.fastq"
    threads: config["threads"]
    conda: "envs/minimap2.yaml"
    shell:
        """
        # 1. Map ke Blacklist -> Ambil Unmapped (-f 4)
        minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - | samtools fastq - > {output.temp_no_mito}

        # 2. Map ke Plastome Ref -> Ambil Mapped (-F 4) (Baiting)
        minimap2 -ax map-ont -t {threads} {input.ref_p} {output.temp_no_mito} \
        | samtools view -b -F 4 - | samtools fastq - > {output.final_reads}
        """

# =============================================================================
# PHASE 4: THE COMPETITION (Multi-Assembler)
# =============================================================================

# Assembler 1: Flye (Specialist Circular)
rule assemble_flye:
    input:
        "results/{sample}/03_filtered_reads/recruited_plastome.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/flye/assembly.fasta",
        info  = "results/{sample}/04_assemblies/flye/assembly_info.txt"
    threads: config["threads"]
    params:
        size = config["plastome_est_size"],
        outdir = "results/{sample}/04_assemblies/flye"
    conda: "envs/flye.yaml"
    shell:
        """
        flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} \
             --genome-size {params.size} --meta --iterations 2
        """

# Assembler 2: Raven (Aggressive)
rule assemble_raven_final:
    input:
        "results/{sample}/03_filtered_reads/recruited_plastome.fastq"
    output:
        "results/{sample}/04_assemblies/raven/assembly.fasta"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell:
        """
        raven --threads {threads} {input} > {output}
        """

# Assembler 3: Canu (High Accuracy - Optional but included)
rule assemble_canu:
    input:
        "results/{sample}/03_filtered_reads/recruited_plastome.fastq"
    output:
        "results/{sample}/04_assemblies/canu/assembly.contigs.fasta"
    threads: config["threads"]
    params:
        size = config["plastome_est_size"],
        dir = "results/{sample}/04_assemblies/canu"
    conda: "envs/canu.yaml"
    shell:
        """
        # Canu perlu parameter khusus agar tidak terlalu lambat di genom kecil
        canu -p assembly -d {params.dir} genomeSize={params.size} \
             -nanopore {input} maxThreads={threads} useGrid=false \
             stopOnLowCoverage=5 minInputCoverage=5
        """

# =============================================================================
# PHASE 5: EQUALIZATION (Polishing)
# =============================================================================

# Poles Raven (karena output raw-nya kasar)
rule polish_raven:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_plastome.fastq",
        draft = "results/{sample}/04_assemblies/raven/assembly.fasta"
    output:
        "results/{sample}/05_polished/raven_polished.fasta"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        "minipolish -t {threads} {input.reads} {input.draft} > {output}"

# Poles Canu (output canu bernama assembly.contigs.fasta)
rule polish_canu:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_plastome.fastq",
        draft = "results/{sample}/04_assemblies/canu/assembly.contigs.fasta"
    output:
        "results/{sample}/05_polished/canu_polished.fasta"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        "minipolish -t {threads} {input.reads} {input.draft} > {output}"

# =============================================================================
# PHASE 6: THE JUDGE (Final Selection)
# =============================================================================

rule select_best_candidate:
    input:
        flye  = "results/{sample}/04_assemblies/flye/assembly.fasta",
        raven = "results/{sample}/05_polished/raven_polished.fasta",
        canu  = "results/{sample}/05_polished/canu_polished.fasta",
        ref   = config["reference_out"]
    output:
        best   = "results/{sample}/07_best_candidate/best_assembly.fasta",
        report = "results/{sample}/07_best_candidate/selection_report.txt"
    conda: "envs/blast_biopython.yaml"
    script:
        "scripts/assess_assemblies.py"
