configfile: "config/config.yaml"
import pandas as pd

# 1. Baca daftar sampel
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)

# 2. Logic Target Akhir (Dynamic Rule All)
# Pipeline akan menentukan output berdasarkan config['target_assembly']
TARGETS = []

# Jika mode PLASTOME atau BOTH, tambahkan hasil akhir plastome assembly
if config["target_assembly"] in ["PLASTOME", "BOTH"]:
    TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))

# Jika mode MITO atau BOTH, tambahkan hasil akhir mito 
# (Untuk saat ini kita targetkan referensinya dulu karena pipeline mito belum full)
if config["target_assembly"] in ["MITO", "BOTH"]:
    TARGETS.append(config["mito_reference_out"])

rule all:
    input: TARGETS

# =============================================================================
# PHASE 1: UNIVERSAL FETCHING (Code Baru)
# =============================================================================

# Download Referensi Plastome (Menggunakan Script Universal)
rule fetch_plastome_ref:
    output:
        config["reference_out"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

# Download Referensi Mitokondria (Menggunakan Script Universal yang SAMA)
rule fetch_mito_ref:
    output:
        config["mito_reference_out"]
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

# =============================================================================
# PHASE 2: ROUGH ASSEMBLY & SMART IDENTIFICATION
# =============================================================================

# Step 1: Raven Rough Assembly
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

# Step 2: Smart Filtering (Graph + Identity Check)
rule create_smart_blacklist:
    input:
        fasta = "results/{sample}/01_rough/assembly.fasta",
        gfa   = "results/{sample}/01_rough/assembly_graph.gfa",
        ref_p = config["reference_out"],
        ref_m = config["mito_reference_out"]
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

# Assembler 1: Flye
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

# Assembler 2: Raven
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

# Assembler 3: Canu
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
        canu -p assembly -d {params.dir} genomeSize={params.size} \
             -nanopore {input} maxThreads={threads} useGrid=false \
             stopOnLowCoverage=5 minInputCoverage=5
        """

# =============================================================================
# PHASE 5: EQUALIZATION (Polishing)
# =============================================================================

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