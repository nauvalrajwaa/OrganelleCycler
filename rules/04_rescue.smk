# PHASE 8: THE RESCUE (CIRCULARIZATION) - MULTI ASSEMBLER LOOP
# Menjalankan MitoHiFi pada Canu, Raven, dan Flye secara paralel.

ruleorder: flag_best_rescue > run_rescue_mitohifi

# --- HELPER FUNCTION ---
def get_draft_by_assembler(wildcards):
    tool = wildcards.assembler.lower()
    if tool == "flye":
        return "results/{sample}/04_assemblies/flye/assembly.fasta"
    elif tool == "raven":
        return "results/{sample}/05_polished/raven_polished.fasta"
    elif tool == "canu":
        return "results/{sample}/05_polished/canu_polished.fasta"
    else:
        return "results/{sample}/07_best_candidate/best_assembly.fasta"

# --- RULE UTAMA: MITOHIFI LOOP ---
rule run_rescue_mitohifi:
    input:
        contigs   = get_draft_by_assembler,
        ref_fasta = get_rescue_ref,     
        ref_gb    = get_rescue_ref_gbk 
    output:
        final_fasta = "results/{sample}/08_rescue/{assembler}/final_circularized.fasta",
        final_gbk   = "results/{sample}/08_rescue/{assembler}/final_mitogenome.gb"
    params:
        out_dir_base = "results/{sample}/08_rescue/{assembler}",
        gcode        = "11" if config["target_assembly"] == "PLASTOME" else "1",
        organism     = "plant",
        perc_match   = "50", # Tetap longgar agar aman
        species      = config.get("species_name", "Saccharum")
    threads: config["threads"]
    
    container: "docker://ghcr.io/marcelauliano/mitohifi:master"
    
    shell:
        r"""
        # 1. Bersihkan Folder Output
        rm -rf {params.out_dir_base}
        mkdir -p {params.out_dir_base}
        
        echo "[INFO] Processing MitoHiFi for: {wildcards.assembler}"

        # 2. FILTER REFERENSI (Single Sequence Strict Mode)
        awk -v sp="{params.species}" '
            BEGIN {{ found=0 }}
            /^>/ {{ if (found) exit; if ($0 ~ sp) found=1; }}
            found {{ print }}         
        ' {input.ref_fasta} > {params.out_dir_base}/input_ref.fasta

        if [ ! -s {params.out_dir_base}/input_ref.fasta ]; then
            awk '/^>/ {{if (seq_count++ > 0) exit}} {{print}}' {input.ref_fasta} > {params.out_dir_base}/input_ref.fasta
        fi

        awk -v sp="{params.species}" -v RS="//" '$0 ~ sp {{ print $0 "//"; exit }}' {input.ref_gb} > {params.out_dir_base}/input_ref.gb
        if [ ! -s {params.out_dir_base}/input_ref.gb ]; then
            awk '{{print}} /^\/\/$/ {{exit}}' {input.ref_gb} > {params.out_dir_base}/input_ref.gb
        fi

        # 3. FILTER INPUT CONTIGS (INITIATE FORCE MODE)
        # Masalah: MitoHiFi sering mengambil contig kecil (4kb-8kb) daripada yang besar.
        # Solusi: Kita filter menggunakan Python One-Liner. Hanya ambil contig > 20.000 bp.
        
        # Copy dulu file asli ke temp
        cp {input.contigs} {params.out_dir_base}/raw_contigs.fasta
        
        # Masuk folder
        cd {params.out_dir_base}

        echo "[INFO] Filtering contigs < 20kb..."
        
        # Python Script untuk membuang contig kecil
        python3 -c "import sys, Bio.SeqIO; Bio.SeqIO.write([r for r in Bio.SeqIO.parse('raw_contigs.fasta', 'fasta') if len(r.seq) > 20000], 'input_contigs.fasta', 'fasta')"
        
        # Cek apakah hasil filter kosong? (Jaga-jaga kalau Canu gagal total)
        if [ ! -s input_contigs.fasta ]; then
             echo "[WARNING] No contigs > 20kb found! Reverting to original input."
             cp raw_contigs.fasta input_contigs.fasta
        else
             COUNT=$(grep -c "^>" input_contigs.fasta)
             echo "[INFO] Filter success. Keeping $COUNT large contigs for MitoHiFi."
        fi

        # 4. JALANKAN MITOHIFI
        # Parameter -p sudah diset 50 di params rule ini
        
        mitohifi.py \
            -c input_contigs.fasta \
            -f input_ref.fasta \
            -g input_ref.gb \
            -t {threads} \
            -o {params.gcode} \
            -a {params.organism} \
            -p {params.perc_match} 
        
        # 5. RENAME & SAVE OUTPUT (FIXED LOGIC)
        TARGET_FASTA=$(basename {output.final_fasta})
        TARGET_GB=$(basename {output.final_gbk})
        
        FOUND_FASTA=$(find . -name "final_mitogenome.fasta" | head -n 1)
        FOUND_GB=$(find . -name "final_mitogenome.gb" | head -n 1)

        if [ ! -z "$FOUND_FASTA" ]; then
            echo "[INFO] Found FASTA: $FOUND_FASTA"
            mv -f "$FOUND_FASTA" "$TARGET_FASTA"
        else
            echo "[WARNING] MitoHiFi FASTA not found. Fallback."
            cp input_contigs.fasta "$TARGET_FASTA"
        fi

        if [ ! -z "$FOUND_GB" ]; then
            if [ ! "$FOUND_GB" -ef "$TARGET_GB" ]; then mv -f "$FOUND_GB" "$TARGET_GB"; fi
        else
            touch "$TARGET_GB"
        fi
        
        # Cleanup
        """

# --- RULE FINAL: FLAG BEST CANDIDATE ---
rule flag_best_rescue:
    input:
        report = "results/{sample}/07_best_candidate/selection_report.txt",
        all_rescues = expand("results/{{sample}}/08_rescue/{assembler}/final_circularized.fasta", 
                             assembler=[x.strip() for x in config.get("loop_mitohifi", "flye").split(",")])
    output:
        best_fasta = "results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta",
        best_gbk   = "results/{sample}/08_rescue/FINAL_BEST/final_mitogenome.gb",
        info       = "results/{sample}/08_rescue/FINAL_BEST/source_info.txt"
    shell:
        r"""
        WINNER=$(grep "WINNER:" {input.report} | awk '{{print $2}}' | tr '[:upper:]' '[:lower:]')
        echo "Detected Winner: $WINNER"
        
        SOURCE_DIR="results/{wildcards.sample}/08_rescue/$WINNER"
        
        if [ -d "$SOURCE_DIR" ]; then
            cp -f $SOURCE_DIR/final_circularized.fasta {output.best_fasta}
            cp -f $SOURCE_DIR/final_mitogenome.gb {output.best_gbk}
            echo "Best assembly: $WINNER" > {output.info}
            echo "Source path: $SOURCE_DIR" >> {output.info}
        else
            touch {output.best_fasta} {output.best_gbk}
            echo "Error copying best candidate" > {output.info}
        fi
        """