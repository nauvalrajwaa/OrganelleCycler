configfile: "config/config.yaml"
import pandas as pd

# =============================================================================
# 1. LOAD SAMPLES
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

# Helper untuk mengambil path GBK referensi
def get_rescue_ref_gbk(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_gbk"]
    else: return config["reference_gbk"]

# =============================================================================
# 3. DEFINISI TARGET AKHIR (CONDITIONAL)
# =============================================================================
TARGETS = []

# Base Targets (Report selalu dibuat)
if config.get("target_assembly") in ["PLASTOME", "MITO", "BOTH"]:
    TARGETS.extend(expand("results/{sample}/02_blacklist/filtering_report.html", sample=samples.index))
    TARGETS.extend(expand("results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html", sample=samples.index))

    # LOGIC CABANG:
    if config.get("run_circularization_rescue") == "YES":
        # Jika YES, kita minta hasil yang sudah di-rescue (MitoHiFi)
        TARGETS.extend(expand("results/{sample}/08_rescue/final_circularized.fasta", sample=samples.index))
    else:
        # Jika NO, kita cukup sampai best candidate saja
        TARGETS.extend(expand("results/{sample}/07_best_candidate/best_assembly.fasta", sample=samples.index))

rule all:
    input: TARGETS

# =============================================================================
# PHASE 1: UNIVERSAL FETCHING
# =============================================================================

rule fetch_plastome_ref:
    output: 
        fasta = config["reference_out"],
        gbk   = config["reference_gbk"] # <-- Output Baru
    log: "logs/universal/fetch_plastome.log"
    params:
        search_term = config["plastome_search_term"],
        min_len = config["plastome_min_len"],
        max_len = config["plastome_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

rule fetch_mito_ref:
    output: 
        fasta = config["mito_reference_out"],
        gbk   = config["mito_reference_gbk"] # <-- Output Baru
    log: "logs/universal/fetch_mito.log"
    params:
        search_term = config["mito_search_term"],
        min_len = config["mito_min_len"],
        max_len = config["mito_max_len"]
    conda: "envs/blast_biopython.yaml"
    script: "scripts/fetch_organelle_ref.py"

# =============================================================================
# PHASE 2: ROUGH ASSEMBLY & FILTERING
# =============================================================================

rule rough_assembly_raven:
    input: lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/raven.fasta",
        gfa   = "results/{sample}/01_rough/raven.gfa"
    log: "logs/{sample}/01_rough_raven.log"
    threads: config["threads"]
    conda: "envs/raven.yaml"
    shell: "raven --threads {threads} --graphical-fragment-assembly {output.gfa} {input} > {output.fasta} 2> {log}"

rule rough_assembly_flye:
    input: lambda wildcards: samples.loc[wildcards.sample, "reads_path"]
    output:
        fasta = "results/{sample}/01_rough/flye/assembly.fasta",
        gfa   = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    log: "logs/{sample}/01_rough_flye.log"
    params: outdir = "results/{sample}/01_rough/flye"
    threads: config["threads"]
    conda: "envs/flye.yaml"
    shell: "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} --iterations 0 2> {log}"

rule plot_rough_graphs:
    input:
        raven_gfa = "results/{sample}/01_rough/raven.gfa",
        flye_gfa  = "results/{sample}/01_rough/flye/assembly_graph.gfa"
    output:
        raven_png = "results/{sample}/01_rough/viz/raven_graph.png",
        flye_png  = "results/{sample}/01_rough/viz/flye_graph.png"
    log: "logs/{sample}/01_plot_rough.log"
    conda: "envs/visualization.yaml"
    # UPDATED: Added --names --lengths --depth
    params: opts = "--height 800 --width 800 --nodewidth 15 --fontsize 20 --names --lengths --depth"
    shell:
        """
        (Bandage image {input.raven_gfa} {output.raven_png} {params.opts}
        Bandage image {input.flye_gfa} {output.flye_png} {params.opts}) 2> {log}
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
    params: opts = "--height 800 --width 800 --nodewidth 15 --fontsize 20 --names --lengths --depth"
    shell:
        """
        (
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then
            touch {output.raven_png}
        fi

        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then
            touch {output.flye_png}
        fi
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
    log: "logs/{sample}/03_filter_reads.log"
    threads: config["threads"]
    conda: "envs/minimap2.yaml"
    shell:
        """
        # 1. Filter Negatif (Blacklist) - Biarkan standar
        (minimap2 -ax map-ont -t {threads} {input.blacklist} {input.reads} \
        | samtools view -b -f 4 - | samtools fastq - > {output.temp_filtered}) 2> {log}

        # 2. BAITING TARGET (RESCUE MODE / SUPER LOOSE)
        # Parameter Penjelasan:
        # -k12 : Gunakan k-mer pendek (default 15). Membantu nempel di daerah banyak error.
        # -A1 -B2 : Skor match=1, mismatch=2 (Default B=4). Hukuman mismatch diperkecil.
        # -O2 -E1 : Gap open/extension penalty diturunkan. Agar reads bolong-bolong tetap diambil.
        # -s 40 : Minimal peak score diturunkan.
        
        (minimap2 -ax map-ont -k12 -A1 -B2 -O2 -E1 -s 40 -t {threads} \
         {input.bait_ref} {output.temp_filtered} \
        | samtools view -b -F 4 - | samtools fastq - > {output.final_reads}) 2>> {log}
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

# Assembler 3: Canu (FAIL-SAFE VERSION)
rule assemble_canu:
    input:
        "results/{sample}/03_filtered_reads/recruited_reads.fastq"
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
        # 1. Jalankan Canu
        # Gunakan '|| true' agar jika Canu error, script bash tidak langsung mati
        if canu -p assembly -d {params.dir} genomeSize={params.size} \
             -nanopore {input} maxThreads={threads} useGrid=false \
             stopOnLowCoverage=5 minInputCoverage=5 > {log} 2>&1; then
            
            # --- JIKA CANU SUKSES ---
            
            # Cek apakah output FASTA terbentuk dan ada isinya
            if [ -s {output.fasta} ]; then
                # Cek GFA. Jika tidak ada, buat manual dari FASTA
                if [ ! -f {output.gfa} ]; then
                    echo "[INFO] Generating GFA from FASTA..." >> {log}
                    awk '/^>/ {{printf("\\nS\\t%s\\t",substr($1,2))}} !/^>/ {{printf("%s",$0)}} END {{printf("\\n")}}' {output.fasta} | sed '/^$/d' > {output.gfa}
                fi
            else
                # Canu sukses run tapi tidak menghasilkan contig (kosong)
                echo "[WARNING] Canu finished but produced no contigs." >> {log}
                touch {output.fasta} {output.gfa}
            fi

        else
            # --- JIKA CANU GAGAL/CRASH ---
            echo "[WARNING] Canu crashed or failed. Creating dummy files to continue pipeline." >> {log}
            touch {output.fasta} {output.gfa}
        fi
        """

# =============================================================================
# PHASE 5: EQUALIZATION (POLISHING) - ROBUST FAIL-SAFE VERSION
# =============================================================================

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
        # 1. Coba jalankan Minipolish ke file temporary
        if minipolish -t {threads} {input.reads} {input.draft} > {output}.temp_gfa 2> {log}; then
            # JIKA SUKSES: Convert GFA Polished ke FASTA
            awk '/^S/{{print ">"$2"\\n"$3}}' {output}.temp_gfa | fold > {output}
            rm {output}.temp_gfa
        else
            # JIKA GAGAL: Gunakan Draft Asli (Unpolished) sebagai fallback
            echo "[WARNING] Minipolish failed. Using unpolished draft as fallback." >> {log}
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
        # Logic yang sama untuk Canu
        if minipolish -t {threads} {input.reads} {input.draft} > {output}.temp_gfa 2> {log}; then
            awk '/^S/{{print ">"$2"\\n"$3}}' {output}.temp_gfa | fold > {output}
            rm {output}.temp_gfa
        else
            echo "[WARNING] Minipolish failed. Using unpolished draft as fallback." >> {log}
            awk '/^S/{{print ">"$2"\\n"$3}}' {input.draft} | fold > {output}
        fi
        """

# =============================================================================
# PHASE 6 & 7: JUDGING & REPORTING
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
    params: opts = "--height 800 --width 800 --nodewidth 20 --fontsize 24 --names --lengths --depth"
    shell:
        """
        (
        if ! Bandage image {input.flye_gfa} {output.flye_png} {params.opts} 2>/dev/null; then touch {output.flye_png}; fi
        if ! Bandage image {input.raven_gfa} {output.raven_png} {params.opts} 2>/dev/null; then touch {output.raven_png}; fi
        if ! Bandage image {input.canu_gfa} {output.canu_png} {params.opts} 2>/dev/null; then touch {output.canu_png}; fi
        ) 2> {log}
        """

# Update bagian ini di Snakefile Anda:

rule generate_final_report:
    input:
        report = "results/{sample}/07_best_candidate/selection_report.txt",
        plots  = ["results/{sample}/07_best_candidate/viz/flye_final.png",
                  "results/{sample}/07_best_candidate/viz/raven_final.png",
                  "results/{sample}/07_best_candidate/viz/canu_final.png"],
        # INPUT BARU: QC PLOTS
        dotplot = "results/{sample}/07_best_candidate/qc/dotplot.png",
        covplot = "results/{sample}/07_best_candidate/qc/coverage_plot.png"
    output:
        html   = "results/{sample}/07_best_candidate/FINAL_ASSEMBLY_REPORT.html"
    script: "scripts/generate_final_report.py"
    
# =============================================================================
# PHASE 8: THE RESCUE (CIRCULARIZATION) - HYBRID INSTALL
# =============================================================================

def get_rescue_ref(wildcards):
    if config.get("target_assembly") == "MITO": return config["mito_reference_out"]
    else: return config["reference_out"]

rule run_rescue_mitohifi:
    input:
        contigs = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref     = get_rescue_ref,     # FASTA (-f)
        ref_gb  = get_rescue_ref_gbk  # GBK (-g) <-- INPUT BARU
    output:
        final = "results/{sample}/08_rescue/final_circularized.fasta",
        gbk   = "results/{sample}/08_rescue/final_mitogenome.gb"
    params:
        outdir = "results/{sample}/08_rescue",
        tool_dir = "resources/tools",
        organism_type = "plant",
        perc_match    = "80",
        # Genetic Code 11 (Bacterial/Plant Plastid)
        # Jika Mito, ganti jadi 1 (Standard)
        gencode       = "11" 
    threads: config["threads"]
    conda: "envs/mitohifi.yaml"
    shell:
        """
        # Clone MitoHiFi (sama seperti sebelumnya)
        if [ ! -d "{params.tool_dir}/MitoHiFi" ]; then
            mkdir -p {params.tool_dir}
            git clone https://github.com/marcelauliano/MitoHiFi.git {params.tool_dir}/MitoHiFi
        fi
        
        rm -rf {params.outdir}
        
        # Jalankan dengan parameter lengkap (-g dan -o)
        python {params.tool_dir}/MitoHiFi/src/mitohifi.py \
            -c {input.contigs} \
            -f {input.ref} \
            -g {input.ref_gb} \
            -a {params.organism_type} \
            -p {params.perc_match} \
            -o {params.gencode} \
            -t {threads} \
            -o {params.outdir} 
        
        # Finalize output (sama seperti sebelumnya)
        if [ -f {params.outdir}/final_mitogenome.fasta ]; then
            cp {params.outdir}/final_mitogenome.fasta {output.final}
            cp {params.outdir}/final_mitogenome.gb {output.gbk}
        else
            echo "[WARNING] Rescue failed. Copying linear input."
            cp {input.contigs} {output.final}
            touch {output.gbk}
        fi
        """