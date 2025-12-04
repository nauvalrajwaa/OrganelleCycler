# PHASE 8: THE RESCUE (CIRCULARIZATION)
# Menggunakan Docker Container resmi MitoHiFi

rule run_rescue_mitohifi:
    input:
        contigs   = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref_fasta = get_rescue_ref,     
        ref_gb    = get_rescue_ref_gbk 
    output:
        # PERBAIKAN: Kembalikan nama ke 'final_circularized.fasta' agar cocok dengan Snakefile
        final_fasta = "results/{sample}/08_rescue/final_circularized.fasta",
        final_gbk   = "results/{sample}/08_rescue/final_mitogenome.gb"
    params:
        out_dir_base = "results/{sample}/08_rescue",
        gcode        = "11" if config["target_assembly"] == "PLASTOME" else "1",
        organism     = "plant",
        perc_match   = "80"
    threads: config["threads"]
    
    # Gunakan Docker via Singularity/Apptainer
    container: "docker://ghcr.io/marcelauliano/mitohifi:master"
    
    shell:
        """
        # 1. Persiapan Folder
        rm -rf {params.out_dir_base}
        mkdir -p {params.out_dir_base}

        # 2. Pindah ke folder kerja & Copy Input
        # (Wajib copy karena Docker kadang susah baca path relative ../../)
        cd {params.out_dir_base}
        
        cp ../../../{input.contigs} ./input_contigs.fasta
        cp ../../../{input.ref_fasta} ./input_ref.fasta
        cp ../../../{input.ref_gb} ./input_ref.gb

        # 3. Jalankan MitoHiFi
        # Output bawaan MitoHiFi selalu bernama: final_mitogenome.fasta
        mitohifi.py \
            -c input_contigs.fasta \
            -f input_ref.fasta \
            -g input_ref.gb \
            -t {threads} \
            -o {params.gcode} \
            -a {params.organism} \
            -p {params.perc_match} \
            --debug
        
        # 4. Rename Hasil ke Target Snakemake
        # Snakemake minta 'final_circularized.fasta', MitoHiFi kasih 'final_mitogenome.fasta'
        
        if [ -f final_mitogenome.fasta ]; then
            echo "[INFO] MitoHiFi Success."
            mv final_mitogenome.fasta ../../../{output.final_fasta}
            mv final_mitogenome.gb ../../../{output.final_gbk}
        else
            echo "[WARNING] MitoHiFi Failed. Fallback to linear input."
            cp input_contigs.fasta ../../../{output.final_fasta}
            touch ../../../{output.final_gbk}
        fi
        
        # Cleanup file temporary
        rm input_contigs.fasta input_ref.fasta input_ref.gb
        """