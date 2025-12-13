import sys
import os
from Bio import SeqIO

# --- INPUT HANDLING ---
# Menggunakan getattr untuk keamanan jika salah satu input kosong/tidak ada
flye_file = getattr(snakemake.input, 'flye', None)
raven_file = getattr(snakemake.input, 'raven', None)

target_size = int(snakemake.params.target_size)
tolerance = float(snakemake.params.tolerance) # Misal 0.2

out_fasta = snakemake.output.best
out_report = snakemake.output.report

def get_stats(fasta_file):
    """
    Mengembalikan dictionary: {
        'len': panjang seq, 
        'circ': boolean (is circular?), 
        'rec': SeqRecord object
    }
    """
    if not fasta_file or not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
        return {'len': 0, 'circ': False, 'rec': None}

    try:
        recs = list(SeqIO.parse(fasta_file, "fasta"))
        if not recs: return {'len': 0, 'circ': False, 'rec': None}
        
        # Ambil contig terpanjang
        longest = max(recs, key=lambda x: len(x.seq))
        
        # Cek sirkularitas dari header (Output Solver biasanya memberi label 'circular=true' atau sejenisnya)
        # Atau cek manual ujung sequence (simple check)
        seq_str = str(longest.seq)
        is_circ = False
        
        # Cek header
        if "circular=true" in longest.description.lower():
            is_circ = True
        # Cek overlap ujung (jika solver belum menandai)
        elif len(seq_str) > 500 and seq_str[:50] == seq_str[-50:]:
            is_circ = True
            
        return {'len': len(seq_str), 'circ': is_circ, 'rec': longest}
    except:
        return {'len': 0, 'circ': False, 'rec': None}

def main():
    # 1. Ambil Statistik
    flye_stat = get_stats(flye_file)
    raven_stat = get_stats(raven_file)
    
    flye_diff = abs(flye_stat['len'] - target_size)
    raven_diff = abs(raven_stat['len'] - target_size)
    
    winner = "None"
    final_rec = None
    reason = "None"

    # 2. Logika Seleksi (Decision Tree)
    
    # Keduanya Kosong
    if flye_stat['len'] == 0 and raven_stat['len'] == 0:
        winner = "None"
        reason = "Both assemblies failed/empty."
        
    # Salah Satu Kosong
    elif flye_stat['len'] == 0:
        winner = "Raven"
        final_rec = raven_stat['rec']
        reason = "Flye failed."
    elif raven_stat['len'] == 0:
        winner = "Flye"
        final_rec = flye_stat['rec']
        reason = "Raven failed."
        
    # Keduanya Ada: Cek Sirkularitas Dulu!
    elif flye_stat['circ'] and not raven_stat['circ']:
        winner = "Flye"
        final_rec = flye_stat['rec']
        reason = "Flye is circular (Raven is linear)."
    elif raven_stat['circ'] and not flye_stat['circ']:
        winner = "Raven"
        final_rec = raven_stat['rec']
        reason = "Raven is circular (Flye is linear)."
        
    # Keduanya Sama Status Sirkularitasnya (Sama-sama bulat atau sama-sama lurus)
    # Pilih berdasarkan kedekatan dengan Target Size
    else:
        if flye_diff < raven_diff:
            winner = "Flye"
            final_rec = flye_stat['rec']
            reason = f"Both {'circular' if flye_stat['circ'] else 'linear'}, Flye closer to target."
        else:
            winner = "Raven"
            final_rec = raven_stat['rec']
            reason = f"Both {'circular' if raven_stat['circ'] else 'linear'}, Raven closer to target."

    # 3. Tulis Output FASTA
    if final_rec:
        # Bersihkan header
        final_rec.id = f"{winner}_Selected"
        final_rec.description = f"len={len(final_rec.seq)} circ={flye_stat['circ'] if winner=='Flye' else raven_stat['circ']}"
        with open(out_fasta, "w") as out:
            SeqIO.write(final_rec, out, "fasta")
    else:
        # Buat file kosong agar rule tidak error
        open(out_fasta, 'a').close()

    # 4. Tulis Report Text
    with open(out_report, "w") as f:
        f.write(f"WINNER: {winner}\n")
        f.write(f"REASON: {reason}\n")
        f.write("-" * 40 + "\n")
        f.write(f"Target Size: {target_size}\n")
        f.write(f"Flye : Len={flye_stat['len']} | Circ={flye_stat['circ']} | Diff={flye_diff}\n")
        f.write(f"Raven: Len={raven_stat['len']} | Circ={raven_stat['circ']} | Diff={raven_diff}\n")

if __name__ == "__main__":
    main()