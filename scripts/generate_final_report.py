import os
import sys

if 'snakemake' in globals():
    input_report = snakemake.input.report # selection_report.txt dari assess_assemblies.py
    input_plots  = snakemake.input.plots  # [flye.png, raven.png, canu.png]
    output_html  = snakemake.output.html
    sample_name  = snakemake.wildcards.sample
else:
    sys.exit("Run via Snakemake only.")

def parse_selection_report(report_path):
    """Membaca file txt report dan mengubahnya jadi list of dicts"""
    data = []
    winner_info = "No Winner"
    
    with open(report_path, 'r') as f:
        lines = f.readlines()
        
    # Cari baris winner
    for line in lines:
        if line.startswith("WINNER:"):
            winner_info = line.strip()
            
    # Cari data tabel (Skip header dan separator)
    # Format di txt: Tool | ContigID | Length | Circ | Cov% | Ident% | SCORE
    start_reading = False
    for line in lines:
        if "Tool | ContigID" in line:
            start_reading = True
            continue
        if "-------" in line:
            continue
        if start_reading:
            if not line.strip(): break # Stop di baris kosong
            if "WINNER" in line or "NO WINNER" in line or "NO CANDIDATES" in line: break
            
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 7:
                data.append({
                    "tool": parts[0],
                    "id": parts[1],
                    "len": parts[2],
                    "circ": parts[3],
                    "cov": parts[4],
                    "ident": parts[5],
                    "score": parts[6]
                })
    return data, winner_info

def main():
    table_data, winner_text = parse_selection_report(input_report)
    
    # HTML Styling
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Final Assembly Report - {sample_name}</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, sans-serif; background: #f0f2f5; padding: 20px; color: #333; }}
            .container {{ max-width: 1200px; margin: auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 5px 15px rgba(0,0,0,0.1); }}
            
            /* Header Section */
            .header {{ text-align: center; border-bottom: 2px solid #eee; padding-bottom: 20px; margin-bottom: 30px; }}
            .header h1 {{ margin: 0; color: #2c3e50; }}
            .header p {{ color: #7f8c8d; font-size: 1.1em; }}
            
            /* Winner Banner */
            .winner-box {{ background: linear-gradient(135deg, #28a745 0%, #218838 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; margin-bottom: 40px; box-shadow: 0 4px 6px rgba(40, 167, 69, 0.3); }}
            .winner-box h2 {{ margin: 0; font-size: 2em; }}
            
            /* Grid Layout for Graphs */
            .grid-container {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 25px; margin-bottom: 40px; }}
            .card {{ background: white; border: 1px solid #e1e4e8; border-radius: 8px; overflow: hidden; transition: transform 0.2s; }}
            .card:hover {{ transform: translateY(-5px); box-shadow: 0 5px 15px rgba(0,0,0,0.1); }}
            .card-header {{ background: #f8f9fa; padding: 15px; border-bottom: 1px solid #e1e4e8; font-weight: bold; text-align: center; text-transform: uppercase; letter-spacing: 1px; }}
            .card-body {{ padding: 10px; text-align: center; }}
            .card img {{ max-width: 100%; height: auto; max-height: 300px; }}
            
            /* Table Styling */
            table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
            th {{ background: #34495e; color: white; padding: 12px; text-align: left; }}
            td {{ padding: 12px; border-bottom: 1px solid #eee; }}
            tr:hover {{ background-color: #f8f9fa; }}
            .score {{ font-weight: bold; color: #2980b9; }}
            .bool-true {{ color: green; font-weight: bold; }}
            .bool-false {{ color: #ccc; }}
            
        </style>
    </head>
    <body>
    <div class="container">
        <div class="header">
            <h1>Final Assembly Competition Report</h1>
            <p>Sample: {sample_name}</p>
        </div>
        
        <div class="winner-box">
            <h2>🏆 {winner_text}</h2>
        </div>
        
        <h2>1. Candidate Visualization</h2>
        <div class="grid-container">
    """
    
    # --- Bagian Gambar ---
    # Mapping nama file ke judul yang cantik
    tool_map = {"flye": "Flye (Assembly)", "raven": "Raven (Polished)", "canu": "Canu (Polished)"}
    
    for plot in input_plots:
        filename = os.path.basename(plot).lower()
        tool_name = "Unknown"
        for key, val in tool_map.items():
            if key in filename:
                tool_name = val
                break
        
        rel_path = os.path.relpath(plot, os.path.dirname(output_html))
        
        html += f"""
        <div class="card">
            <div class="card-header">{tool_name}</div>
            <div class="card-body">
                <img src="{rel_path}" alt="{tool_name}">
            </div>
        </div>
        """
    html += "</div>"
    
    # --- Bagian Tabel ---
    html += """
    <h2>2. Detailed Metrics</h2>
    <table>
        <thead>
            <tr>
                <th>Tool</th>
                <th>Contig ID</th>
                <th>Length (bp)</th>
                <th>Circular?</th>
                <th>Coverage %</th>
                <th>Identity %</th>
                <th>Final Score</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for row in table_data:
        circ_class = "bool-true" if "True" in row['circ'] or "1" in row['circ'] else "bool-false"
        html += f"""
        <tr>
            <td><b>{row['tool']}</b></td>
            <td>{row['id']}</td>
            <td>{row['len']}</td>
            <td class="{circ_class}">{row['circ']}</td>
            <td>{row['cov']}</td>
            <td>{row['ident']}</td>
            <td class="score">{row['score']}</td>
        </tr>
        """
        
    html += """
        </tbody>
    </table>
    
    <div style="margin-top: 30px; padding: 15px; background: #e8f4fd; border-left: 4px solid #3498db; border-radius: 4px;">
        <strong>Note:</strong> Score is calculated based on Identity + Coverage + Circularity Bonus (+500). High scores indicate likely complete organelle genomes.
    </div>
    
    </div>
    </body>
    </html>
    """
    
    with open(output_html, "w") as f:
        f.write(html)

if __name__ == "__main__":
    main()