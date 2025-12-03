# /scripts/generate_final_report.py

import os
import sys

# --- INPUT HANDLING ---
if 'snakemake' in globals():
    input_report = snakemake.input.report
    input_plots  = snakemake.input.plots  # [flye.png, raven.png, canu.png]
    
    # Input QC Baru (Gunakan getattr agar aman jika file belum terbentuk)
    input_dotplot = getattr(snakemake.input, 'dotplot', None)
    input_covplot = getattr(snakemake.input, 'covplot', None)
    
    output_html  = snakemake.output.html
    sample_name  = snakemake.wildcards.sample
else:
    sys.exit("Run via Snakemake only.")

def get_rel_path(target_path):
    """Mengubah absolute path menjadi relative path untuk HTML"""
    if target_path and os.path.exists(target_path):
        return os.path.relpath(target_path, os.path.dirname(output_html))
    return None

def parse_selection_report(report_path):
    """
    Membaca report txt format baru (RANKING MODE).
    Format: RANK | TOOL | TotalLen | Circ | Cov% | Ident% | SCORE
    """
    data = []
    winner_info = "No Winner"
    
    with open(report_path, 'r') as f:
        lines = f.readlines()
        
    # 1. Cari Winner
    for line in lines:
        if line.startswith("WINNER:"):
            winner_info = line.strip().replace("WINNER:", "").strip()
            
    # 2. Cari Data Tabel
    start_reading = False
    for line in lines:
        # Deteksi Header Baru
        if "RANK | TOOL" in line:
            start_reading = True
            continue
        if "-----" in line:
            continue
            
        if start_reading:
            if not line.strip(): break
            if "WINNER" in line or "NO WINNER" in line or "NO CANDIDATES" in line: break
            
            parts = [p.strip() for p in line.split('|')]
            # Pastikan baris memiliki kolom yang cukup (minimal 7 kolom)
            if len(parts) >= 7:
                data.append({
                    "rank": parts[0],
                    "tool": parts[1],
                    "len": parts[2],
                    "circ": parts[3],
                    "cov": parts[4],
                    "ident": parts[5],
                    "score": parts[6]
                })
    return data, winner_info

def main():
    table_data, winner_text = parse_selection_report(input_report)
    
    # Paths gambar QC
    dotplot_path = get_rel_path(input_dotplot)
    covplot_path = get_rel_path(input_covplot)
    
    # HTML Styling
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Organelle Assembly Report - {sample_name}</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <style>
            :root {{
                --primary: #2c3e50;
                --accent: #3498db;
                --success: #27ae60;
                --bg: #f4f7f6;
                --card-bg: #ffffff;
            }}
            body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: var(--bg); padding: 20px; color: #333; line-height: 1.6; }}
            .container {{ max-width: 1200px; margin: auto; }}
            
            /* Header */
            .header {{ text-align: center; margin-bottom: 40px; padding: 20px; background: white; border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }}
            .header h1 {{ margin: 0; color: var(--primary); font-size: 2.2em; }}
            .header p {{ color: #7f8c8d; font-size: 1.2em; margin-top: 10px; }}
            
            /* Winner Banner */
            .winner-box {{ background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); color: white; padding: 30px; border-radius: 12px; text-align: center; margin-bottom: 40px; box-shadow: 0 10px 20px rgba(56, 239, 125, 0.3); transform: translateY(-5px); }}
            .winner-box h2 {{ margin: 0; font-size: 2.5em; text-shadow: 0 2px 4px rgba(0,0,0,0.2); }}
            .winner-box p {{ margin-top: 10px; font-size: 1.1em; opacity: 0.9; }}
            
            /* Section Titles */
            h2 {{ color: var(--primary); border-bottom: 3px solid var(--accent); display: inline-block; padding-bottom: 5px; margin-top: 40px; }}
            
            /* Cards & Grid */
            .grid-container {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 30px; margin-bottom: 40px; }}
            .card {{ background: var(--card-bg); border-radius: 12px; overflow: hidden; box-shadow: 0 4px 6px rgba(0,0,0,0.05); transition: transform 0.3s; }}
            .card:hover {{ transform: translateY(-5px); box-shadow: 0 10px 20px rgba(0,0,0,0.1); }}
            
            .card-header {{ background: #f8f9fa; padding: 15px; border-bottom: 1px solid #eee; font-weight: 700; color: var(--primary); text-align: center; letter-spacing: 0.5px; }}
            .card-body {{ padding: 15px; text-align: center; }}
            .card img {{ max-width: 100%; height: auto; border-radius: 4px; }}
            .card-desc {{ font-size: 0.9em; color: #666; margin-top: 10px; padding: 0 10px; }}

            /* Table */
            .table-wrapper {{ overflow-x: auto; background: white; border-radius: 12px; padding: 20px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }}
            table {{ width: 100%; border-collapse: collapse; }}
            th {{ background: var(--primary); color: white; padding: 15px; text-align: left; font-weight: 600; }}
            td {{ padding: 15px; border-bottom: 1px solid #eee; }}
            tr:hover {{ background-color: #f1f1f1; }}
            
            .badge {{ padding: 5px 10px; border-radius: 20px; font-size: 0.85em; font-weight: bold; }}
            .rank-1 {{ background: #ffd700; color: #B37400; }}
            .rank-2 {{ background: #C0C0C0; color: #555; }}
            .rank-3 {{ background: #CD7F32; color: white; }}
            .score {{ color: var(--accent); font-weight: bold; font-size: 1.1em; }}
            
        </style>
    </head>
    <body>
    <div class="container">
        <div class="header">
            <h1>🧬 Final Assembly Report</h1>
            <p>Sample ID: <strong>{sample_name}</strong></p>
        </div>
        
        <div class="winner-box">
            <h2>🏆 {winner_text}</h2>
            <p>Best Candidate Consensus</p>
        </div>
        
        <h2>1. Assembly Graph Topologies</h2>
        <div class="grid-container">
    """
    
    # --- LOOP GAMBAR ASSEMBLY (Bandage) ---
    tool_map = {"flye": "Flye (Graph)", "raven": "Raven (Graph)", "canu": "Canu (Graph)"}
    
    for plot in input_plots:
        filename = os.path.basename(plot).lower()
        tool_name = "Unknown"
        for key, val in tool_map.items():
            if key in filename:
                tool_name = val
                break
        
        rel_path = get_rel_path(plot)
        if rel_path:
            html += f"""
            <div class="card">
                <div class="card-header">{tool_name}</div>
                <div class="card-body">
                    <img src="{rel_path}" alt="{tool_name}">
                </div>
            </div>
            """
        else:
            html += f"""<div class="card"><div class="card-body">Image not found: {filename}</div></div>"""
            
    html += "</div>"
    
    # --- SECTION 2: QUALITY CONTROL (BARU) ---
    html += """
        <h2>2. Quality Control & Validation</h2>
        <div class="grid-container">
    """
    
    # Card 1: Dotplot
    if dotplot_path:
        html += f"""
        <div class="card">
            <div class="card-header">Synteny Dotplot (Assembly vs Reference)</div>
            <div class="card-body">
                <img src="{dotplot_path}" alt="Dotplot">
                <div class="card-desc">
                    <ul style="text-align:left; font-size:0.85em;">
                        <li><strong>Diagonal Line:</strong> Good alignment / synteny.</li>
                        <li><strong>X-Shape:</strong> Indicates Inverted Repeats (IR) present (Good for Plastome).</li>
                        <li><strong>Breaks/Shifts:</strong> Structural rearrangements or misassembly.</li>
                    </ul>
                </div>
            </div>
        </div>
        """
    else:
        html += """<div class="card"><div class="card-body">Dotplot not available (QC failed?)</div></div>"""
        
    # Card 2: Coverage
    if covplot_path:
        html += f"""
        <div class="card">
            <div class="card-header">Read Coverage Distribution</div>
            <div class="card-body">
                <img src="{covplot_path}" alt="Coverage Plot">
                <div class="card-desc">
                    Checks if the assembly is supported by reads uniformly. Sharp drops to zero indicate gaps.
                </div>
            </div>
        </div>
        """
    else:
        html += """<div class="card"><div class="card-body">Coverage plot not available</div></div>"""
        
    html += "</div>" # End QC Grid
    
    # --- SECTION 3: METRICS TABLE ---
    html += """
    <h2>3. Comparative Metrics (Ranked)</h2>
    <div class="table-wrapper">
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Tool</th>
                <th>Total Length (bp)</th>
                <th>Circular?</th>
                <th>Coverage %</th>
                <th>Identity %</th>
                <th>Score</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for row in table_data:
        # Style Ranking
        rank_display = row['rank']
        if row['rank'] == '1':
            rank_display = f'<span class="badge rank-1">RANK 1</span>'
        elif row['rank'] == '2':
            rank_display = f'<span class="badge rank-2">RANK 2</span>'
        elif row['rank'] == '3':
            rank_display = f'<span class="badge rank-3">RANK 3</span>'
            
        html += f"""
        <tr>
            <td>{rank_display}</td>
            <td><b>{row['tool']}</b></td>
            <td>{row['len']}</td>
            <td>{row['circ']}</td>
            <td>{row['cov']}</td>
            <td>{row['ident']}</td>
            <td class="score">{row['score']}</td>
        </tr>
        """
        
    html += """
        </tbody>
    </table>
    </div>
    
    <div style="margin-top: 30px; padding: 20px; background: #e8f4fd; border-left: 5px solid #3498db; border-radius: 4px; color: #444;">
        <strong>💡 Interpretation Guide:</strong>
        <ul style="margin-top:10px; padding-left:20px;">
            <li><strong>Total Length:</strong> Sum of all valid contigs. Should be close to target (Plastome ~150kb).</li>
            <li><strong>Score:</strong> Composite metric based on (Coverage x 2) + Identity - Length Penalty.</li>
            <li><strong>QC Dotplot:</strong> Look for 'X' pattern (Inverted Repeats) to confirm complete plastome structure.</li>
        </ul>
    </div>
    
    </div>
    </body>
    </html>
    """
    
    with open(output_html, "w") as f:
        f.write(html)

if __name__ == "__main__":
    main()