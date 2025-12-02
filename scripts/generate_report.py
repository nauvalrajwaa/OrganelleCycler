# /scripts/generate_report.py

import os
import sys

if 'snakemake' in globals():
    input_log = snakemake.input.log
    input_plots = snakemake.input.plots
    output_html = snakemake.output.html
else:
    sys.exit("Run via Snakemake only.")

def main():
    with open(input_log, 'r') as f:
        log_lines = f.readlines()
        
    header = log_lines[0] if log_lines else ""
    
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Smart Filter Report</title>
        <style>
            body { font-family: sans-serif; padding: 20px; background: #f9f9f9; color: #333; }
            .container { max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 10px rgba(0,0,0,0.1); }
            h1 { border-bottom: 2px solid #007bff; padding-bottom: 15px; color: #007bff; }
            h2 { margin-top: 30px; color: #555; border-left: 5px solid #007bff; padding-left: 10px; }
            .img-grid { display: flex; gap: 20px; margin-bottom: 30px; flex-wrap: wrap; }
            .img-box { flex: 1; min-width: 300px; text-align: center; border: 1px solid #eee; padding: 15px; border-radius: 8px; background: #fff; }
            .img-box h3 { margin-top: 0; color: #666; }
            img { max-width: 100%; height: auto; border-radius: 4px; border: 1px solid #ddd; }
            table { width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; }
            th, td { padding: 12px 15px; border-bottom: 1px solid #ddd; text-align: left; }
            th { background: #333; color: white; text-transform: uppercase; letter-spacing: 1px; }
            .keep { color: #28a745; font-weight: bold; background-color: #e6ffed; padding: 5px 10px; border-radius: 4px; }
            .blacklist { color: #dc3545; font-weight: bold; background-color: #ffe6e6; padding: 5px 10px; border-radius: 4px; }
        </style>
    </head>
    <body>
    <div class="container">
        <h1>Smart Filter Analysis Report (Filtered Result)</h1>
        
        <h2>1. Filtered Graph Visualization (Max 5 Top Contigs)</h2>
        <div class='img-grid'>
    """
    
    for plot in input_plots:
        if os.path.exists(plot):
            rel_path = os.path.relpath(plot, os.path.dirname(output_html))
            name = os.path.basename(plot).replace("_filtered.png", "").upper()
            html += f"<div class='img-box'><h3>{name} FILTERED</h3><img src='{rel_path}'></div>"
            
    html += f"</div><h2>2. Filtering Decisions ({header.split('|')[0].strip()})</h2>"
    
    html += """
    <table>
        <thead>
            <tr>
                <th>Tool</th>
                <th>Contig ID</th>
                <th>Length</th>
                <th>Circular?</th>
                <th>Target Score</th>
                <th>Contam Score</th>
                <th>Decision</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for line in log_lines[1:]:
        parts = line.strip().split(" | ")
        if len(parts) < 7: continue
        
        tool, cid, length, circ, score_t, score_c, decision = parts
        css_class = "blacklist" if "BLACKLIST" in decision else "keep"
        
        html += f"""
        <tr>
            <td><b>{tool}</b></td>
            <td>{cid}</td>
            <td>{length}</td>
            <td>{circ}</td>
            <td>{score_t}</td>
            <td>{score_c}</td>
            <td><span class="{css_class}">{decision}</span></td>
        </tr>
        """
        
    html += "</tbody></table></div></body></html>"
    
    with open(output_html, "w") as f:
        f.write(html)

if __name__ == "__main__":
    main()