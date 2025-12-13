# =============================================================================
# PHASE 4: THE RESCUE (CIRCULARIZATION) - MITOHIFI LOOP
# File: rules/04_rescue.smk
# Strategy: 
# 1. Jalankan MitoHiFi pada Flye dan Raven (Parallel).
# 2. Ambil hasil terbaik berdasarkan siapa pemenang di Step 02 (select_best).
# =============================================================================

# Config default untuk looping (Hanya Flye dan Raven)
ASSEMBLERS = ["flye", "raven"]

ruleorder: finalize_best_rescue > run_mitohifi_loop

# --- HELPER FUNCTION ---
def get_draft_by_assembler(wildcards):
    """Mengambil Raw Assembly dari Step 02"""
    tool = wildcards.assembler.lower()
    if tool == "flye":
        return "results/{sample}/04_assemblies/flye/assembly.fasta"
    elif tool == "raven":
        return "results/{sample}/04_assemblies/raven/assembly.fasta"
    else:
        # Fallback safe
        return "results/{sample}/04_assemblies/flye/assembly.fasta"

# --- RULE 1: MITOHIFI LOOP (FLYE & RAVEN) ---
rule run_mitohifi_loop:
    input:
        contigs   = get_draft_by_assembler,
        ref_fasta = config["target_ref"],      
        ref_gb    = config["ref_gb"] # Pastikan ini ada di config.yaml
    output:
        final_fasta = "results/{sample}/08_rescue/{assembler}/final_mitogenome.fasta",
        final_gbk   = "results/{sample}/08_rescue/{assembler}/final_mitogenome.gb"
    params:
        out_dir_base = "results/{sample}/08_rescue/{assembler}",
        species      = config.get("species_name", "Plant"),
        # Parameter MitoHiFi
        perc_match   = 50, 
        min_len      = 10000, # Filter contig input (10kb)
        organism     = "plant",
        gcode        = "11" if config.get("target_assembly", "PLASTOME") == "PLASTOME" else "1"
    threads: config["threads"]
    container: "docker://ghcr.io/marcelauliano/mitohifi:master"
    log: "logs/{sample}/08_rescue_{assembler}.log"
    shell:
        r"""
        # 1. Setup Folder
        rm -rf {params.out_dir_base}
        mkdir -p {params.out_dir_base}
        
        echo "[INFO] Running MitoHiFi for: {wildcards.assembler}" > {log}

        # 2. FILTER REFERENSI (Ambil seq pertama saja agar aman)
        awk '/^>/ {{if (seq_count++ > 0) exit}} {{print}}' {input.ref_fasta} > {params.out_dir_base}/input_ref.fasta
        awk '{{print}} /^\/\/$/ {{exit}}' {input.ref_gb} > {params.out_dir_base}/input_ref.gb

        # 3. FILTER INPUT CONTIGS (Hanya ambil yang cukup panjang)
        # Menggunakan Python inline yang lebih aman
        python3 -c "import sys, Bio.SeqIO; 
recs = [r for r in Bio.SeqIO.parse('{input.contigs}', 'fasta') if len(r.seq) > {params.min_len}];
Bio.SeqIO.write(recs, '{params.out_dir_base}/filtered_contigs.fasta', 'fasta')" >> {log} 2>&1

        cd {params.out_dir_base}

        # Cek apakah ada contig tersisa
        if [ ! -s filtered_contigs.fasta ]; then
            echo "[WARNING] No contigs > {params.min_len}bp found. Using original." >> ../../../{log}
            cp ../../../{input.contigs} filtered_contigs.fasta
        fi

        # 4. JALANKAN MITOHIFI
        # find_mito_reference.py sudah deprecated di versi baru, langsung mitohifi.py
        mitohifi.py \
            -c filtered_contigs.fasta \
            -f input_ref.fasta \
            -g input_ref.gb \
            -t {threads} \
            -o {params.gcode} \
            -a {params.organism} \
            -p {params.perc_match} >> ../../../{log} 2>&1 || true

        # 5. RENAMING OUTPUT (Standardisasi nama file output)
        # MitoHiFi outputnya agak random namanya, kita cari yang 'final_mitogenome'
        
        FOUND_FASTA=$(find . -maxdepth 1 -name "*final_mitogenome.fasta" | head -n 1)
        FOUND_GB=$(find . -maxdepth 1 -name "*final_mitogenome.gb" | head -n 1)

        if [ ! -z "$FOUND_FASTA" ]; then
            mv "$FOUND_FASTA" final_mitogenome.fasta
            mv "$FOUND_GB" final_mitogenome.gb
            echo "[SUCCESS] MitoHiFi finished." >> ../../../{log}
        else
            echo "[FAIL] MitoHiFi failed to circularize. Creating dummy files." >> ../../../{log}
            # Buat file kosong/dummy agar rule tidak error, tapi user tahu ini gagal
            cp filtered_contigs.fasta final_mitogenome.fasta
            touch final_mitogenome.gb
        fi
        """

# --- RULE 2: FINALIZE (PICK WINNER) ---
rule finalize_best_rescue:
    input:
        report      = "results/{sample}/07_best_candidate/selection_report.txt",
        # Kita paksa rule di atas jalan untuk kedua tool
        all_rescues = expand("results/{{sample}}/08_rescue/{assembler}/final_mitogenome.fasta", 
                             assembler=ASSEMBLERS)
    output:
        best_fasta = "results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta",
        best_gbk   = "results/{sample}/08_rescue/FINAL_BEST/final_mitogenome.gb",
        info       = "results/{sample}/08_rescue/FINAL_BEST/source_info.txt"
    shell:
        r"""
        # 1. Baca siapa pemenang dari Step 02
        WINNER=$(grep "WINNER:" {input.report} | awk '{{print $2}}' | tr '[:upper:]' '[:lower:]')
        
        # Default ke flye jika logic error
        if [ -z "$WINNER" ] || [ "$WINNER" == "none" ]; then WINNER="flye"; fi
        
        echo "Selected Winner based on QC: $WINNER"
        
        SOURCE_DIR="results/{wildcards.sample}/08_rescue/$WINNER"
        
        # 2. Copy Hasil
        mkdir -p $(dirname {output.best_fasta})
        
        if [ -s "$SOURCE_DIR/final_mitogenome.gb" ]; then
            # Jika MitoHiFi sukses menghasilkan GenBank (tanda sirkularisasi berhasil)
            cp $SOURCE_DIR/final_mitogenome.fasta {output.best_fasta}
            cp $SOURCE_DIR/final_mitogenome.gb {output.best_gbk}
            echo "Source: $WINNER (MitoHiFi Success)" > {output.info}
        else
            # Jika MitoHiFi gagal pada pemenang, coba cek tool satunya (Fallback)
            OTHER="flye"
            if [ "$WINNER" == "flye" ]; then OTHER="raven"; fi
            
            ALT_DIR="results/{wildcards.sample}/08_rescue/$OTHER"
            
            if [ -s "$ALT_DIR/final_mitogenome.gb" ]; then
                 echo "Winner ($WINNER) failed rescue, but ($OTHER) succeeded. Switching."
                 cp $ALT_DIR/final_mitogenome.fasta {output.best_fasta}
                 cp $ALT_DIR/final_mitogenome.gb {output.best_gbk}
                 echo "Source: $OTHER (Fallback Rescue)" > {output.info}
            else
                 echo "Both failed rescue. Keeping linear draft of $WINNER."
                 cp $SOURCE_DIR/final_mitogenome.fasta {output.best_fasta}
                 touch {output.best_gbk} # Empty GBK
                 echo "Source: $WINNER (Rescue Failed)" > {output.info}
            fi
        fi
        """