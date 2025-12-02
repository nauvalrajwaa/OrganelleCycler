configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# 1. LOAD SAMPLES (WAJIB PALING ATAS)
# =============================================================================
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

if config.get("target_assembly") in ["PLASTOME", "MITO", "BOTH"]:
    # 1. Final Assembly
    TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))
    # 2. Smart Filter Report
    TARGETS.extend(expand("results/{sample}/02_blacklist/filtering_report.html", sample=samples.index))
    # 3. Final Report
    TARGETS.extend(expand("results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html", sample=samples.index))

rule all:
    input: TARGETS

# =============================================================================
# PHASE 1: UNIVERSAL FETCHING
# =============================================================================

rule fetch_plastome_ref:
    output: config["reference_out"]
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

rule fetch_mito_ref:
    output: config["mito_reference_out"]
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

# =============================================================================
# PHASE 2: ROUGH ASSEMBLY & FILTERING
# =============================================================================

# 1. Rough Raven
rule rough_assembly_raven:
    input: lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/raven.fasta",
        gfa   = "results/{sample}/01_rough/raven.gfa"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell: "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}"

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

# 3. Visualisasi Graph Kasar
rule plot_rough_graphs:
    input:
        raven_gfa = "results/{sample}/01_rough/raven.gfa",
        flye_gfa  = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    output:
        raven_png = "results/{sample}/01_rough/viz/raven_graph.png",
        flye_png  = "results/{sample}/01_rough/viz/flye_graph.png"
    conda: "envs/visualization.yaml"
    params: opts = "--height 800 --width 800 --nodewidth 15 --fontsize 20"
    shell:
        """
        Bandage image {input.raven_gfa} {output.raven_png} {params.opts}
        Bandage image {input.flye_gfa} {output.flye_png} {params.opts}
        """

# 4. Smart Filtering (LOGIC ONLY)
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
    conda: "envs/blast_biopython.yaml"
    script: "scripts/smart_filter.py"

# 5. Plot Filtered Graphs (FAIL-SAFE VERSION)
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
        # Coba plot Raven. Jika error (misal file kosong), buat file png kosong agar lanjut.
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then
            echo "[WARNING] Raven GFA invalid or empty. Creating dummy image."
            touch {output.raven_png}
        fi

        # Coba plot Flye. Jika error, buat file png kosong.
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then
            echo "[WARNING] Flye GFA invalid or empty. Creating dummy image."
            touch {output.flye_png}
        fi
        """

# 6. Generate Report (VISUALIZATION)
rule generate_filter_report:
    input:
        log   = "results/{sample}/02_blacklist/filter_logic.log",
        plots = ["results/{sample}/02_blacklist/viz/raven_filtered.png",
                 "results/{sample}/02_blacklist/viz/flye_filtered.png"]
    output:
        html  = "results/{sample}/02_blacklist/filtering_report.html"
    script: "scripts/generate_report.py"

# =============================================================================
# PHASE 3: DUAL FILTERING
# =============================================================================

rule filter_reads_dual:
    input:
        reads = lambda wildcards: samples.loc[wildcards.sample, "reads_path"],
        blacklist = "results/{sample}/02_blacklist/smart_blacklist.fasta",
        bait_ref = get_target_ref 
    output:
        temp_filtered = temp("results/{sample}/03_temp/filtered_neg.fastq"),
        final_reads   = "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    threads: config["threads"]
    conda: "envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - | samtools fastq - > {output.temp_filtered}

        minimap2 -ax map-ont -t {threads} {input.bait_ref} {output.temp_filtered} \
        | samtools view -b -F 4 - | samtools fastq - > {output.final_reads}
        """

# =============================================================================
# PHASE 4: THE COMPETITION (ASSEMBLY)
# =============================================================================

rule assemble_flye:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/flye/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/flye/assembly_graph.gfa",
        info  = "results/{sample}/04_assemblies/flye/assembly_info.txt"
    threads: config["threads"]
    params:
        size = config["target_est_size"],
        outdir = "results/{sample}/04_assemblies/flye"
    conda: "envs/flye.yaml"
    shell:
        "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --genome-size {params.size} --meta --iterations 2"

rule assemble_raven_final:
    input: "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/raven/assembly.fasta",
        gfa   = "results/{sample}/04_assemblies/raven/assembly_graph.gfa"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell:
        "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}"

# Assembler 3: Canu (Fixed AWK Syntax)
rule assemble_canu:
    input:
        "results/{sample}/03_filtered_reads/recruited_reads.fastq"
    output:
        fasta = "results/{sample}/04_assemblies/canu/assembly.contigs.fasta",
        gfa   = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    threads: config["threads"]
    params:
        size = config["target_est_size"],
        dir = "results/{sample}/04_assemblies/canu"
    conda: "envs/canu.yaml"
    shell:
        """
        # 1. Jalankan Canu
        # (Tambahkan '|| true' agar jika Canu exit code aneh tapi output ada, dia tetap lanjut)
        canu -p assembly -d {params.dir} genomeSize={params.size} \
             -nanopore {input} maxThreads={threads} useGrid=false \
             stopOnLowCoverage=5 minInputCoverage=5

        # 2. Fail-Safe GFA Generation
        if [ ! -f {output.gfa} ]; then
            echo "[INFO] Canu did not produce GFA. Generating one from FASTA..."
            
            # AWK command dengan escaping yang benar untuk Snakemake:
            # Perhatikan {{ }} untuk setiap kurung kurawal AWK
            awk '/^>/ {{printf("\\nS\\t%s\\t",substr($1,2))}} !/^>/ {{printf("%s",$0)}} END {{printf("\\n")}}' {output.fasta} | sed '/^$/d' > {output.gfa}
        fi
        """

# =============================================================================
# PHASE 5: EQUALIZATION (POLISHING) - FIX INPUT/OUTPUT FORMAT
# =============================================================================

rule polish_raven:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        # FIX: Gunakan GFA sebagai input, bukan FASTA
        draft = "results/{sample}/04_assemblies/raven/assembly_graph.gfa"
    output:
        "results/{sample}/05_polished/raven_polished.fasta"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        """
        # 1. Jalankan Minipolish (Outputnya GFA)
        # 2. Pipe (|) ke AWK untuk convert GFA jadi FASTA
        #    (Ambil baris 'S', cetak '>Nama', lalu cetak 'Sequence')
        minipolish -t {threads} {input.reads} {input.draft} \
        | awk '/^S/{{print ">"$2"\\n"$3}}' | fold > {output}
        """

rule polish_canu:
    input:
        reads = "results/{sample}/03_filtered_reads/recruited_reads.fastq",
        # FIX: Gunakan GFA (yang kita generate manual tadi) sebagai input
        draft = "results/{sample}/04_assemblies/canu/assembly.contigs.gfa"
    output:
        "results/{sample}/05_polished/canu_polished.fasta"
    threads: config["threads"]
    conda: "envs/minipolish.yaml"
    shell:
        """
        minipolish -t {threads} {input.reads} {input.draft} \
        | awk '/^S/{{print ">"$2"\\n"$3}}' | fold > {output}
        """

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
    script: "scripts/assess_assemblies.py"

# =============================================================================
# PHASE 7: FINAL VISUALIZATION & REPORT
# =============================================================================

# 1. Plot 3 Finalis (FAIL-SAFE VERSION)
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
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then
            touch {output.flye_png}
        fi
        
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then
            touch {output.raven_png}
        fi
        
        if ! Bandage image {input.canu_gfa} {output.canu_png} {params.opts} 2>/dev/null; then
            touch {output.canu_png}
        fi
        """

rule generate_final_report:
    input:
        report = "results/{sample}/07_best_candidate/selection_report.txt",
        plots  = ["results/{sample}/07_best_candidate/viz/flye_final.png",
                  "results/{sample}/07_best_candidate/viz/raven_final.png",
                  "results/{sample}/07_best_candidate/viz/canu_final.png"]
    output:
        html   = "results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html"
    script: "scripts/generate_final_report.py"