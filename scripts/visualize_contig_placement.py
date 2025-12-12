import sys
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import base64
from io import BytesIO, StringIO

# --- INPUTS ---
assembly_fasta = snakemake.input.assembly
reference_fasta = snakemake.input.reference
output_html = snakemake.output.html
sample_name = snakemake.wildcards.sample
assembler_name = snakemake.wildcards.assembler

# --- 1. JALANKAN BLAST (TABULAR) UNTUK GRAFIK & LAPORAN GAP ---
def run_blast_tabular(query, subject):
    fmt = "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
    cmd = [
        "blastn", 
        "-query", query, 
        "-subject", subject, 
        "-outfmt", fmt, 
        "-max_target_seqs", "10",
        "-evalue", "1e-10"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return pd.DataFrame()
    
    col_names = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore"]
    if not result.stdout:
        return pd.DataFrame(columns=col_names)
        
    return pd.read_csv(StringIO(result.stdout), sep="\t", names=col_names)

# --- 2. JALANKAN BLAST (ANCHORED) UNTUK VISUALISASI TEKS ---
def run_blast_pairwise(query, subject):
    """
    Menggunakan outfmt 1 (Query-anchored).
    Kecocokan ditandai dengan TITIK (.), Perbedaan ditandai HURUF.
    Ini memudahkan melihat GAP dan MISMATCH.
    """
    cmd = [
        "blastn", 
        "-query", query, 
        "-subject", subject, 
        "-outfmt", "1",  # <-- GANTI KE FORMAT 1 (Anchored)
        "-max_target_seqs", "1", 
        "-evalue", "1e-10"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return f"Error running BLAST: {result.stderr}"
    
    return result.stdout

# --- 3. ANALISIS GAP (UNALIGNED REGIONS) ---
def analyze_unaligned_regions(df, contig_lengths):
    """
    Menganalisis bagian contig mana yang TIDAK ter-align (Gap di ujung atau tengah).
    """
    report = []
    
    # Kelompokkan berdasarkan contig
    for qseqid, group in df.groupby("qseqid"):
        ctg_len = contig_lengths.get(qseqid, 0)
        
        # Urutkan berdasarkan posisi start di query
        sorted_hits = group.sort_values("qstart")
        
        # Cek Gap Awal (Ujung Kiri)
        first_start = sorted_hits.iloc[0]['qstart']
        if first_start > 1:
            report.append({
                "contig": qseqid,
                "type": "Unaligned Start",
                "start": 1,
                "end": first_start - 1,
                "size": first_start - 1
            })
            
        # Cek Gap Akhir (Ujung Kanan)
        last_end = sorted_hits.iloc[-1]['qend']
        if last_end < ctg_len:
            report.append({
                "contig": qseqid,
                "type": "Unaligned End",
                "start": last_end + 1,
                "end": ctg_len,
                "size": ctg_len - last_end
            })
            
        # Cek Gap Tengah (Internal Gaps antar HSPs)
        prev_end = sorted_hits.iloc[0]['qend']
        for i in range(1, len(sorted_hits)):
            curr_start = sorted_hits.iloc[i]['qstart']
            if curr_start > prev_end + 1:
                report.append({
                    "contig": qseqid,
                    "type": "Internal Gap",
                    "start": prev_end + 1,
                    "end": curr_start - 1,
                    "size": (curr_start - 1) - (prev_end + 1)
                })
            prev_end = max(prev_end, sorted_hits.iloc[i]['qend'])

    return pd.DataFrame(report)

# --- 4. GENERATE PLOTS ---
def fig_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight', dpi=100)
    buf.seek(0)
    plt.close(fig)
    return base64.b64encode(buf.read()).decode('utf-8')

def plot_linear_alignment(df, ref_len, contig_lengths):
    fig, ax = plt.subplots(figsize=(15, len(contig_lengths) * 0.4 + 2))
    ax.plot([0, ref_len], [0, 0], color="black", linewidth=3, label="Reference")
    ax.text(ref_len/2, 0.5, f"Reference ({ref_len/1000:.1f} kb)", ha='center', va='bottom', fontsize=12, fontweight='bold')

    sorted_contigs = sorted(contig_lengths.items(), key=lambda x: x[1], reverse=True)
    y_pos = -1
    labels = []
    yticks = []

    for ctg_id, ctg_len in sorted_contigs:
        ctg_hits = df[df['qseqid'] == ctg_id]
        if not ctg_hits.empty:
            for _, hit in ctg_hits.iterrows():
                start = min(hit['sstart'], hit['send'])
                end = max(hit['sstart'], hit['send'])
                width = end - start
                color = "green" if hit['pident'] >= 95 else "orange" if hit['pident'] >= 80 else "red"
                rect = patches.Rectangle((start, y_pos - 0.4), width, 0.8, linewidth=0, edgecolor=None, facecolor=color, alpha=0.8)
                ax.add_patch(rect)
        labels.append(f"{ctg_id} ({ctg_len/1000:.1f}kb)")
        yticks.append(y_pos)
        y_pos -= 1

    ax.set_yticks(yticks)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Reference Genome Position (bp)")
    ax.set_title(f"Alignment Map: {assembler_name}", fontsize=14)
    ax.set_ylim(y_pos, 2)
    ax.set_xlim(0, ref_len)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    return fig_to_base64(fig)

def plot_dotplot(df, best_contig_id, ref_len, ctg_len):
    fig, ax = plt.subplots(figsize=(10, 10))
    hits = df[df['qseqid'] == best_contig_id]
    for _, hit in hits.iterrows():
        x = [hit['sstart'], hit['send']]
        y = [hit['qstart'], hit['qend']]
        color = 'blue' if (hit['sstart'] < hit['send']) == (hit['qstart'] < hit['qend']) else 'red'
        ax.plot(x, y, color=color, linewidth=2, alpha=0.7)
    ax.set_xlim(0, ref_len)
    ax.set_ylim(0, ctg_len)
    ax.set_xlabel(f"Reference ({ref_len} bp)")
    ax.set_ylabel(f"Contig ({ctg_len} bp)")
    ax.set_title(f"Dotplot: {best_contig_id}", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.3)
    return fig_to_base64(fig)

# --- MAIN EXECUTION ---

# 1. Parsing Fasta
try:
    ref_record = list(SeqIO.parse(reference_fasta, "fasta"))[0]
    ref_len = len(ref_record.seq)
    contig_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(assembly_fasta, "fasta")}
except Exception as e:
    print(f"Error parsing FASTA: {e}")
    sys.exit(1)

# 2. Run BLAST
df = run_blast_tabular(assembly_fasta, reference_fasta)
pairwise_text = run_blast_pairwise(assembly_fasta, reference_fasta)

# 3. Analisis Gap
gap_report_df = analyze_unaligned_regions(df, contig_lengths)

# 4. Build HTML
html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Alignment: {sample_name} ({assembler_name})</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; color: #333; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; border-bottom: 1px solid #eee; padding-bottom: 5px;}}
        .section {{ background: #fff; padding: 20px; margin-bottom: 30px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .stats {{ background: #e8f6f3; padding: 15px; border-left: 5px solid #1abc9c; margin-bottom: 20px; }}
        img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}

        /* TEXT ALIGNMENT BOX */
        .seq-box {{
            background-color: #f8f9fa;
            border: 1px solid #ccc;
            padding: 15px;
            overflow-x: auto;
            font-family: 'Courier New', Courier, monospace; 
            font-size: 13px;
            white-space: pre; 
            max-height: 600px;
            overflow-y: auto;
            line-height: 1.2;
        }}
        .legend {{ font-size: 0.9em; color: #666; margin-bottom: 10px; font-style: italic; }}
    </style>
</head>
<body>
    <h1>Alignment Report: {sample_name}</h1>
    
    <div class="stats">
        <p><strong>Assembler:</strong> {assembler_name}</p>
        <p><strong>Reference Length:</strong> {ref_len:,} bp</p>
        <p><strong>Aligned Hits:</strong> {len(df)} regions</p>
    </div>
"""

if not df.empty:
    img_linear = plot_linear_alignment(df, ref_len, contig_lengths)
    
    # Pilih contig terbaik
    best_contig_id = max(contig_lengths, key=contig_lengths.get)
    best_len = contig_lengths[best_contig_id]
    img_dot = plot_dotplot(df, best_contig_id, ref_len, best_len)
    
    # 1. GAP REPORT TABLE
    html_content += """
    <div class="section">
        <h2>1. Missing / Unaligned Regions Report</h2>
        <p>Tabel ini menunjukkan bagian dari Contig yang <strong>TIDAK</strong> ditemukan pada referensi (Gap).</p>
    """
    if not gap_report_df.empty:
        html_content += gap_report_df.to_html(index=False, classes="table")
    else:
        html_content += "<p><strong>Perfect Coverage!</strong> Tidak ada gap signifikan pada contig (Start-to-End aligned).</p>"
    html_content += "</div>"

    # 2. GRAFIK
    html_content += f"""
    <div class="section">
        <h2>2. Visual Alignment Map</h2>
        <img src="data:image/png;base64,{img_linear}">
    </div>

    <div class="section">
        <h2>3. Structural Dotplot ({best_contig_id})</h2>
        <img src="data:image/png;base64,{img_dot}">
    </div>
    """
    
    # 3. TEXT ALIGNMENT
    html_content += f"""
    <div class="section">
        <h2>4. Sequence Alignment (Anchored View)</h2>
        <p class="legend">
            <strong>Cara Membaca:</strong><br>
            - Baris Atas: Query (Contig Anda)<br>
            - Baris Bawah: Reference<br>
            - Tanda <strong>Titik (.)</strong> = MATCH (Sama persis)<br>
            - Tanda <strong>Huruf (A/T/C/G)</strong> pada baris bawah = MISMATCH (Beda)<br>
            - Tanda <strong>Strip (-)</strong> = GAP (Indel)<br>
            Fokuslah pada bagian yang <strong>bukan titik</strong> untuk melihat perbedaannya.
        </p>
        <div class="seq-box">
{pairwise_text}
        </div>
    </div>
    """

else:
    html_content += "<h2 style='color:red'>No Alignment Found!</h2>"

html_content += "</body></html>"

with open(output_html, "w") as f:
    f.write(html_content)