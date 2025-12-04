# PHASE 8: THE RESCUE (CIRCULARIZATION)
# Menggunakan Docker Container resmi MitoHiFi untuk stabilitas maksimal

rule run_rescue_mitohifi:
    input:
        contigs   = "results/{sample}/07_best_candidate/best_assembly.fasta",
        ref_fasta = get_rescue_ref,     
        ref_gb    = get_rescue_ref_gbk 
    output:
        final_fasta = "results/{sample}/08_rescue/final_mitogenome.fasta",
        final_gbk   = "results/{sample}/08_rescue/final_mitogenome.gb"
    params:
        # MitoHiFi membuat folder berdasarkan parameter '-c' (input fasta).
        # Kita perlu trick sedikit agar outputnya masuk ke folder yang kita mau.
        out_dir_base = "results/{sample}/08_rescue",
        # Genetic Code: 11 (Plastid/Bacterial), 1 (Mito Standard)
        gcode   = "11" if config["target_assembly"] == "PLASTOME" else "1",
        organism = "plant",
        perc_match = "80"
    threads: config["threads"]
    
    # --- DOCKER CONTAINER DIRECTIVE ---
    # Snakemake akan otomatis pull image ini dan menjalankannya via Apptainer/Singularity
    container: "docker://ghcr.io/marcelauliano/mitohifi:master"
    
    shell:
        """
        # 1. Bersihkan output folder lama agar bersih
        rm -rf {params.out_dir_base}
        mkdir -p {params.out_dir_base}

        # 2. EKSEKUSI MITOHIFI (DOCKERIZED)
        # Catatan: Di dalam container, script 'mitohifi.py' sudah ada di PATH global.
        # Kita tidak perlu path absolut resource/tools lagi.
        
        # Kita masuk ke folder output dulu agar MitoHiFi menulis di sana
        cd {params.out_dir_base}

        # Perhatian: MitoHiFi butuh path input absolut atau relative dari tempat script dijalankan.
        # Karena kita 'cd', kita harus menyesuaikan path input (naik 3 level: ../../../)
        # Atau lebih aman pakai path absolut snakemake: {input.contigs} (biasanya relative dari root).
        # Trik paling aman: Copy input ke folder kerja saat ini.
        
        cp ../../../{input.contigs} ./input_contigs.fasta
        cp ../../../{input.ref_fasta} ./input_ref.fasta
        cp ../../../{input.ref_gb} ./input_ref.gb

        mitohifi.py \
            -c input_contigs.fasta \
            -f input_ref.fasta \
            -g input_ref.gb \
            -t {threads} \
            -o {params.gcode} \
            -a {params.organism} \
            -p {params.perc_match} \
            --debug
        
        # 3. Rename/Cleanup Output agar sesuai target Snakemake
        # MitoHiFi biasanya menghasilkan output dengan nama: final_mitogenome.fasta
        
        # Cek apakah sukses
        if [ -f final_mitogenome.fasta ]; then
            echo "[INFO] MitoHiFi Success."
            # Tidak perlu cp karena kita sudah di folder tujuan, file sudah bernama final_mitogenome.fasta
        else
            echo "[WARNING] MitoHiFi Failed to circularize. Fallback to linear input."
            cp input_contigs.fasta final_mitogenome.fasta
            touch final_mitogenome.gb
        fi
        
        # Hapus file temporary copy tadi
        rm input_contigs.fasta input_ref.fasta input_ref.gb
        """