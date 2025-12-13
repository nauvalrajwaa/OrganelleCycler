import os
import sys

# --- INPUT HANDLING ---
if 'snakemake' in globals():
    input_metric = snakemake.input.sel_metric
    input_quast  = snakemake.input.quast_txt
    
    # Images
    plot_flye    = snakemake.input.flye_img
    plot_raven   = snakemake.input.raven_img
    plot_dot     = snakemake.input.dotplot
    plot_cov     = snakemake.input.covplot
    
    output_html  = snakemake.output.html
    sample_name  = snakemake.wildcards.sample
else:
    sys.exit("Run via Snakemake only.")

def get_rel_path(target_path):
    """Agar HTML bisa dibuka di mana saja (relative path)"""
    if target_path and os.path.exists(target_path):
        return os.path.relpath(target_path, os.path.dirname(output_html))
    return None

def read_file_content(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return f.read()
    return "File not found."

def parse_selection_metric(filepath):
    """Parsing output dari select_best.py"""
    winner = "Unknown"
    reason = ""
    rows = []
    
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("WINNER:"):
                    winner = line.split(":")[1].strip()
                elif line.startswith("REASON:"):
                    reason = line.split(":")[1].strip()
                elif "|" in line: 
                    # Format: TOOL | LENGTH | CIRCULAR | DIFF
                    parts = [p.strip() for p in line.split("|")]
                    if len(parts) >= 4:
                        rows.append(parts)
    return winner, reason, rows

def main():
    winner_tool, winner_reason, stats_rows = parse_selection_metric(input_metric)
    quast_content = read_file_content(input_quast)
    
    # Get Relative Paths for Images
    p_flye = get_rel_path(plot_flye)
    p_raven = get_rel_path(plot_raven)
    p_dot = get_rel_path(plot_dot)
    p_cov = get_rel_path(plot_cov)

    # --- HTML GENERATION ---
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Assembly Report: {sample_name}</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: #f4f4f9; color: #333; margin: 0; padding: 20px; }}
            .container {{ max-width: 1000px; margin: auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 5px 15px rgba(0,0,0,0.05); }}
            
            h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #2980b9; margin-top: 40px; }}
            
            /* Winner Box */
            .winner-box {{ background: linear-gradient(135deg, #6dd5ed, #2193b0); color: white; padding: 20px; border-radius: 8px; text-align: center; margin-bottom: 30px; }}
            .winner-box h1 {{ border: none; color: white; margin: 0; font-size: 2.5em; }}
            .winner-reason {{ font-style: italic; opacity: 0.9; margin-top: 5px; }}

            /* Table */
            table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
            th, td {{ padding: 12px; border-bottom: 1px solid #ddd; text-align: left; }}
            th {{ background-color: #f8f9fa; color: #2c3e50; }}
            
            /* Grid for Images */
            .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 20px; }}
            .card {{ border: 1px solid #eee; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
            .card-header {{ background: #f8f9fa; padding: 10px; font-weight: bold; text-align: center; color: #555; }}
            .card img {{ width: 100%; height: auto; display: block; }}
            .card-desc {{ padding: 10px; font-size: 0.85em; color: #777; background: #fff; }}

            /* QUAST Box */
            .code-box {{ background: #2d3436; color: #dfe6e9; padding: 15px; border-radius: 5px; overflow-x: auto; font-family: monospace; font-size: 0.9em; }}
        </style>
    </head>
    <body>
    <div class="container">
        
        <div style="text-align:center; margin-bottom:30px;">
            <p style="margin:0; color:#7f8c8d;">Organelle Assembly Pipeline</p>
            <h1 style="border:none; margin-top:5px;">Final Report: {sample_name}</h1>
        </div>

        <div class="winner-box">
            <p style="margin:0;">Best Candidate Selected:</p>
            <h1>🏆 {winner_tool}</h1>
            <div class="winner-reason">Decision: {winner_reason}</div>
        </div>

        <h2>1. Candidate Comparison</h2>
        <table>
            <thead>
                <tr>
                    <th>Assembler</th>
                    <th>Length (bp)</th>
                    <th>Circular?</th>
                    <th>Diff to Target</th>
                </tr>
            </thead>
            <tbody>
    """
    
    # Isi Tabel dari select_best.py
    for row in stats_rows:
        tool = row[0]
        length = row[1]
        circ = row[2]
        diff = row[3]
        
        # Highlight pemenang di tabel
        style = "font-weight:bold; color:#27ae60;" if tool == winner_tool else ""
        html += f"<tr style='{style}'><td>{tool}</td><td>{length}</td><td>{circ}</td><td>{diff}</td></tr>"

    html += """
            </tbody>
        </table>

        <h2>2. Assembly Graphs (Cleaned)</h2>
        <div class="grid">
            <div class="card">
                <div class="card-header">Flye Graph</div>
                <div><img src="{}" alt="Flye Graph"></div>
            </div>
            <div class="card">
                <div class="card-header">Raven Graph</div>
                <div><img src="{}" alt="Raven Graph"></div>
            </div>
        </div>
    """.format(p_flye, p_raven)

    html += """
        <h2>3. Final Validation (QC)</h2>
        <div class="grid">
            <div class="card">
                <div class="card-header">Synteny Dotplot (vs Reference)</div>
                <div><img src="{}" alt="Dotplot"></div>
                <div class="card-desc">
                    <b>Check for:</b> A clear diagonal line. An 'X' pattern indicates Inverted Repeats (common in plastomes). Breaks indicate structural errors.
                </div>
            </div>
            <div class="card">
                <div class="card-header">Read Coverage Depth</div>
                <div><img src="{}" alt="Coverage"></div>
                <div class="card-desc">
                    <b>Check for:</b> Uniform coverage (blue area). Zero coverage means the assembly has gaps unsupported by reads.
                </div>
            </div>
        </div>
    """.format(p_dot, p_cov)

    html += f"""
        <h2>4. QUAST Metrics</h2>
        <div class="code-box">
            <pre>{quast_content}</pre>
        </div>
        
    </div>
    </body>
    </html>
    """

    with open(output_html, "w") as f:
        f.write(html)

if __name__ == "__main__":
    main()