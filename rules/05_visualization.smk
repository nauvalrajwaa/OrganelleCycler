# =============================================================================
# PHASE 5: VISUALIZATION
# File: rules/05_visualization.smk
# Tugas: Memetakan contig ke referensi untuk melihat Synteny & Coverage
# =============================================================================

# --- HELPER FUNCTION ---
def get_viz_assembly_input(wildcards):
    tool = wildcards.assembler.lower()
    
    # 1. FINAL RESULT (Hasil Rescue/MitoHiFi)
    if tool == "final":
        return "results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta"
    
    # 2. INTERMEDIATE RESULTS (Hasil Cleaning & Solving Step 02)
    elif tool == "flye":
        return "results/{sample}/04_assemblies/flye/assembly_circular.fasta"
    elif tool == "raven":
        return "results/{sample}/04_assemblies/raven/assembly_circular.fasta"
        
    # Default Fallback
    return "results/{sample}/08_rescue/FINAL_BEST/final_circularized.fasta"

rule visualize_contigs:
    input:
        assembly  = get_viz_assembly_input,
        reference = config["target_ref"]  # Menggunakan referensi target (Plastome/Mito)
    output:
        html = "results/{sample}/09_viz/{assembler}_alignment_map.html"
    
    # Pastikan env ini punya: python, biopython, plotly, pandas, minimap2
    conda: "../envs/visualization.yaml"
    
    script:
        "../scripts/visualize_contig_placement.py"