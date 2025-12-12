# PHASE 9: VISUALISASI ALIGNMENT
# Memetakan contig assembly ke referensi untuk melihat coverage & overhang.

# --- HELPER FUNCTION ---
def get_viz_assembly_input(wildcards):
    tool = wildcards.assembler.lower()
    if tool == "flye":
        # Flye ambil dari raw output
        return "results/{sample}/04_assemblies/flye/assembly.fasta"
    elif tool == "raven":
        return "results/{sample}/05_polished/raven_polished.fasta"
    elif tool == "canu":
        return "results/{sample}/05_polished/canu_polished.fasta"
    else:
        return "results/{sample}/07_best_candidate/best_assembly.fasta"

rule visualize_contigs:
    input:
        assembly = get_viz_assembly_input,
        # Referensi utama
        reference = config.get("reference_out", "resources/ref_tebu_plastome.fasta")
    output:
        html = "results/{sample}/09_viz/{assembler}_alignment_map.html"
    
    # GUNAKAN ENV BARU DI SINI
    conda:
        "../envs/viz_env.yaml"
        
    script:
        "../scripts/visualize_contig_placement.py"