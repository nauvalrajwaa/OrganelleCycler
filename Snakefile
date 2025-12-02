configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# 1. LOAD SAMPLES (WAJIB PALING ATAS SETELAH IMPORT)
# =============================================================================
# Kita harus baca file ini duluan supaya variabel 'samples.index' tersedia
samples = pd.read_csv(config["samples_file"], sep="\t").set_index("sample_id", drop=False)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def get_target_ref(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_out"]
    else: return config["reference_out"]

def get_min_len(wildcards):
    if config.get("target_assembly") == "MITO": return config.get("mito_min_len", 200000)
    else: return config.get("plastome_min", 100000)

def get_max_len(wildcards):
    if config.get("target_assembly") == "MITO": return config.get("mito_max_len", 2000000)
    else: return config.get("plastome_max", 200000)

# =============================================================================
# 3. DEFINISI TARGET AKHIR (RULE ALL)
# =============================================================================
TARGETS = []

# Cek mode assembly
if config.get("target_assembly") in ["PLASTOME", "MITO", "BOTH"]:
    # 1. Final Assembly Fasta (Hasil Utama)
    TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))
    
    # 2. Smart Filter Report (Laporan Tengah - HTML)
    TARGETS.extend(expand("results/{sample}/02_blacklist/filtering_report.html", sample=samples.index))
    
    # 3. FINAL REPORT (Laporan Akhir - HTML) <-- BARU
    TARGETS.extend(expand("results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html", sample=samples.index))

rule all:
    input: TARGETS

TARGETS = []

# Tentukan output akhir berdasarkan target
# Apapun targetnya (PLASTOME atau MITO), kita ingin file final 'best_assembly.fasta'
if config.get("target_assembly") in ["PLASTOME", "MITO", "BOTH"]:
    # 1. Final Assembly
    TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))
    # 2. Smart Filter Report (HTML) - Wajib ada agar rule generate report dijalankan
    TARGETS.extend(expand("results/{sample}/02_blacklist/filtering_report.html", sample=samples.index))

rule all:
    input: TARGETS

# =============================================================================
# PHASE 1: UNIVERSAL FETCHING
# =============================================================================

# Download Referensi Plastome
rule fetch_plastome_ref:
    output:
        config["reference_out"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

# Download Referensi Mitokondria
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
# PHASE 2: ROUGH ASSEMBLY & FILTERING (VISUALISASI HASIL FILTER)
# =============================================================================

# 1. Rough Raven
rule rough_assembly_raven:
    input: lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/raven.fasta",
        gfa   = "results/{sample}/01_rough/raven.gfa"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell: "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input.reads} > {output.fasta}"

# 2. Rough Flye
rule rough_assembly_flye:
    input: lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/flye/assembly.fasta",
        gfa   = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    params: outdir = "results/{sample}/01_rough/flye"
    threads: config["threads"]
    conda: "envs/flye.yaml"
    shell: "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --iterations 0"

# 3. Smart Filtering (Menghasilkan GFA Bersih)
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
        # Output Baru: GFA yang sudah difilter (hanya target)
        filtered_gfas = ["results/{sample}/02_blacklist/raven_filtered.gfa",
                         "results/{sample}/02_blacklist/flye_filtered.gfa"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/smart_filter.py"

# 4. Visualisasi Graph (SEKARANG MENGGUNAKAN GFA TERFILTER)
rule plot_filtered_graphs:
    input:
        raven_gfa = "results/{sample}/02_blacklist/raven_filtered.gfa",
        flye_gfa  = "results/{sample}/02_blacklist/flye_filtered.gfa"
    output:
        raven_png = "results/{sample}/02_blacklist/viz/raven_filtered.png",
        flye_png  = "results/{sample}/02_blacklist/viz/flye_filtered.png"
    conda: "envs/visualization.yaml"
    params: opts = "--height 800 --width 800 --nodewidth 15 --fontsize 20"
    shell:
        """
        Bandage image {input.raven_gfa} {output.raven_png} {params.opts}
        Bandage image {input.flye_gfa} {output.flye_png} {params.opts}
        """

# 5. Generate Report (Menggunakan Gambar Filtered)
rule generate_filter_report:
    input:
        log   = "results/{sample}/02_blacklist/filter_logic.log",
        plots = ["results/{sample}/02_blacklist/viz/raven_filtered.png",
                 "results/{sample}/02_blacklist/viz/flye_filtered.png"]
    output:
        html  = "results/{sample}/02_blacklist/filtering_report.html"
    script: "scripts/generate_report.py"

# =============================================================================
# PHASE 3: DUAL FILTERING (Cleaning Reads)
# =============================================================================

# Step 3: Negative Filter (Buang Blacklist) -> Positive Filter (Ambil Target)
rule filter_reads_dual:
    input:
        reads = lambda wildcards: samples.loc[wildcards.sample, "reads_path"],
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        # DINAMIS: Gunakan referensi yang sesuai target untuk memancing (Baiting)
        bait_ref = get_target_ref 
    output:
        temp_filtered = temp("results/{sample}/03_temp/filtered_neg.fastq"),
        final_reads   = "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    threads: config["threads"]
    conda: "envs/minimap2.yaml"
    shell:
        """
        # 1. Negative Filter: Map ke Blacklist -> Ambil yang TIDAK nempel (Unmapped -f 4)
        minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - | samtools fastq - > {output.temp_filtered}

        # 2. Positive Filter: Map ke Target Ref -> Ambil yang NEMPEL (Mapped -F 4) 

[Image of pairwise sequence alignment]

        minimap2 -ax map-ont -t {threads} {input.bait_ref} {output.temp_filtered} \
        | samtools view -b -F 4 - | samtools fastq - > {output.final_reads}
        """

# =============================================================================
# PHASE 4: THE COMPETITION (Updated Raven for GFA)
# =============================================================================

# Assembler 1: Flye (Output GFA sudah ada default)
rule assemble_flye:
    input:
        "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/flye/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/flye/assembly_graph.gfa", # Pastikan ini terambil
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

# Assembler 2: Raven (Updated: Add --graphical-fragment-assembly)
rule assemble_raven_final:
    input:
        "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/raven/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/raven/assembly_graph.gfa" # NEW OUTPUT
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell:
        # Tambahkan flag GFA
        "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}"

# Assembler 3: Canu (GFA is auto-generated named *.contigs.gfa)
rule assemble_canu:
    input:
        "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/canu/assembly.contigs.fasta",
        # Canu menaruh GFA di folder output dengan nama prefix.contigs.gfa
        # Kita definisikan pathnya agar snakemake tahu
        gfa   = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
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

# ... (Phase 5 Polishing TETAP SAMA) ...
# Note: Polishing tidak mengubah struktur graph secara drastis, 
# jadi kita akan visualisasikan GFA dari Phase 4 untuk report.

# =============================================================================
# PHASE 6: THE JUDGE
# =============================================================================

rule select_best_candidate:
    input:
        flye  = "results/{sample}/04_assemblies/flye/assembly.fasta",
        raven = "results/{sample}/05_polished/raven_polished.fasta",
        canu  = "results/{sample}/05_polished/canu_polished.fasta",
        ref   = get_target_ref  
    output:
        best   = "results/{sample}/07_best_candidate/best_assembly.fasta",
        report = "results/{sample}/07_best_candidate/selection_report.txt"
    params:
        min_len = get_min_len,
        max_len = get_max_len
    conda: "envs/blast_biopython.yaml"
    script:
        "scripts/assess_assemblies.py"

# =============================================================================
# PHASE 7: FINAL VISUALIZATION & REPORT (NEW!)
# =============================================================================

# 1. Plot 3 Finalis
rule plot_final_graphs:
    input:
        flye_gfa  = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        raven_gfa = "results/{sample}/04_assemblies/raven/assembly_graph.gfa",
        canu_gfa  = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    output:
        flye_png  = "results/{sample}/07_best_candidate/viz/flye_final.png",
        raven_png = "results/{sample}/07_best_candidate/viz/raven_final.png",
        canu_png  = "results/{sample}/07_best_candidate/viz/canu_final.png"
    conda: "envs/visualization.yaml"
    params:
        opts = "--height 800 --width 800 --nodewidth 20 --fontsize 24"
    shell:
        """
        Bandage image {input.flye_gfa} {output.flye_png} {params.opts}
        Bandage image {input.raven_gfa} {output.raven_png} {params.opts}
        Bandage image {input.canu_gfa} {output.canu_png} {params.opts}
        """

# 2. Generate Final HTML Report
rule generate_final_report:
    input:
        report = "results/{sample}/07_best_candidate/selection_report.txt",
        plots  = ["results/{sample}/07_best_candidate/viz/flye_final.png",
                  "results/{sample}/07_best_candidate/viz/raven_final.png",
                  "results/{sample}/07_best_candidate/viz/canu_final.png"]
    output:
        html   = "results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html"
    script:
        "scripts/generate_final_report.py"