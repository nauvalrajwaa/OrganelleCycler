import sys
import os

# --- INPUT HANDLING ---
if 'snakemake' in globals():
    input_stats = snakemake.input.stats
    output_html = snakemake.output.html
    sample_name = snakemake.wildcards.sample
else:
    sys.exit("Run via Snakemake only.")

def main():
    # 1. Baca Data Statistik
    stats = {}
    with open(input_stats, 'r') as f:
        for line in f:
            key, val = line.strip().split(': ')
            stats[key] = int(val)
    
    total = stats.get('Total_Raw_Reads', 0)
    after_bl = stats.get('After_Blacklist', 0)
    final = stats.get('Final_Recruited', 0)
    
    # Hitung Delta (Yang dibuang)
    removed_blacklist = total - after_bl
    removed_nontarget = after_bl - final
    
    # Hitung Persentase (handle division by zero)
    def calc_pct(val, tot):
        return (val / tot * 100) if tot > 0 else 0

    pct_final = calc_pct(final, total)
    pct_bl = calc_pct(removed_blacklist, total)
    pct_nt = calc_pct(removed_nontarget, total)

    # 2. Generate HTML
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Filtering Report: {sample_name}</title>
        <style>
            body {{ font-family: 'Segoe UI', sans-serif; background: #f4f6f8; color: #333; padding: 20px; }}
            .container {{ max-width: 800px; margin: auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 4px 20px rgba(0,0,0,0.05); }}
            h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            
            .summary-box {{ display: flex; gap: 20px; margin-bottom: 40px; }}
            .stat-card {{ flex: 1; background: #f8f9fa; padding: 20px; border-radius: 8px; text-align: center; border: 1px solid #eee; }}
            .stat-num {{ font-size: 2em; font-weight: bold; color: #2c3e50; }}
            .stat-label {{ color: #7f8c8d; font-size: 0.9em; text-transform: uppercase; letter-spacing: 1px; }}
            
            /* Bar Chart CSS */
            .chart-container {{ margin-top: 30px; }}
            .bar-group {{ margin-bottom: 25px; }}
            .bar-label {{ font-weight: bold; margin-bottom: 8px; display: flex; justify-content: space-between; }}
            .progress {{ height: 25px; background-color: #e9ecef; border-radius: 15px; overflow: hidden; }}
            .bar {{ height: 100%; line-height: 25px; color: white; text-align: center; font-size: 0.85em; font-weight: bold; transition: width 1s ease-in-out; }}
            
            .bg-keep {{ background: #27ae60; }}     /* Hijau */
            .bg-bl {{ background: #c0392b; }}       /* Merah */
            .bg-trash {{ background: #f39c12; }}    /* Oranye */
            
            .legend {{ margin-top: 30px; font-size: 0.9em; color: #666; background: #fffbe6; padding: 15px; border-left: 4px solid #f39c12; }}
        </style>
    </head>
    <body>
    <div class="container">
        <h1>🔍 Reads Filtering Report</h1>
        <p>Sample: <strong>{sample_name}</strong></p>

        <div class="summary-box">
            <div class="stat-card">
                <div class="stat-num">{total:,}</div>
                <div class="stat-label">Total Input Reads</div>
            </div>
            <div class="stat-card" style="border-bottom: 4px solid #27ae60;">
                <div class="stat-num">{final:,}</div>
                <div class="stat-label">Final Recruited</div>
            </div>
        </div>

        <h2>Where did the reads go?</h2>
        
        <div class="chart-container">
            <div class="bar-group">
                <div class="bar-label">
                    <span>🚫 Removed by Blacklist (Nuclear/Mito)</span>
                    <span>{removed_blacklist:,} reads ({pct_bl:.1f}%)</span>
                </div>
                <div class="progress">
                    <div class="bar bg-bl" style="width: {pct_bl}%"></div>
                </div>
            </div>

            <div class="bar-group">
                <div class="bar-label">
                    <span>🗑️ Removed (Non-Target / Junk)</span>
                    <span>{removed_nontarget:,} reads ({pct_nt:.1f}%)</span>
                </div>
                <div class="progress">
                    <div class="bar bg-trash" style="width: {pct_nt}%"></div>
                </div>
            </div>

            <div class="bar-group">
                <div class="bar-label">
                    <span>✅ Kept (Target Matches)</span>
                    <span>{final:,} reads ({pct_final:.1f}%)</span>
                </div>
                <div class="progress">
                    <div class="bar bg-keep" style="width: {pct_final}%"></div>
                </div>
            </div>
        </div>
        

        <div class="legend">
            <strong>Interpretation:</strong>
            <ul>
                <li><strong>Removed by Blacklist:</strong> DNA yang cocok dengan referensi blacklist (misal: genom inti inang). Ini dibuang duluan.</li>
                <li><strong>Removed (Non-Target):</strong> DNA yang lolos blacklist tapi TIDAK cocok dengan referensi target. Kemungkinan kontaminasi bakteri atau DNA sampah.</li>
                <li><strong>Kept:</strong> DNA yang spesifik menempel pada referensi organel target Anda. Ini yang dipakai untuk Assembly.</li>
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