import sys
import os
import subprocess
import shutil
import argparse
from collections import defaultdict

# ==========================================
# BAGIAN 1: STRUKTUR DATA GRAFIK
# ==========================================

class GraphNode:
    def __init__(self, name, length, coverage, sequence):
        self.name = name
        self.length = length
        self.coverage = coverage
        self.sequence = sequence
        
        # Status Node
        self.is_anchor = False      # Punya gen target (Label)?
        self.keep = False           # Apakah akan disimpan?
        self.blast_hit_len = 0      # Panjang alignment BLAST terbaik

class AssemblyGraph:
    def __init__(self):
        self.nodes = {}
        self.edges = defaultdict(set) # Adjacency list: node -> {neighbors}

    def add_node(self, name, length, coverage, sequence):
        self.nodes[name] = GraphNode(name, length, coverage, sequence)

    def add_edge(self, u, v):
        # Grafik DNA bersifat bidirectional, simpan kedua arah
        self.edges[u].add(v)
        self.edges[v].add(u)

    def remove_node(self, name):
        if name in self.nodes:
            del self.nodes[name]
        # Hapus semua edge yang mengarah ke node ini
        if name in self.edges:
            del self.edges[name]
        for u in self.edges:
            if name in self.edges[u]:
                self.edges[u].remove(name)

# ==========================================
# BAGIAN 2: PARSER GFA (FLYE)
# ==========================================

def parse_flye_gfa(gfa_path):
    print(f"[CLEANER] Membaca grafik GFA: {gfa_path}")
    graph = AssemblyGraph()
    
    try:
        with open(gfa_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if not parts: continue
                
                # Parsing Node (Segment)
                if parts[0] == 'S':
                    name = parts[1]
                    seq = parts[2]
                    length = len(seq)
                    
                    # Cari coverage Flye (format dp:f:123.45)
                    cov = 0.0
                    for tag in parts[3:]:
                        if tag.startswith('dp:f:'):
                            cov = float(tag.split(':')[-1])
                            break
                    
                    graph.add_node(name, length, cov, seq)

                # Parsing Edge (Link)
                elif parts[0] == 'L':
                    u, v = parts[1], parts[3]
                    # Pastikan kedua node ada
                    if u in graph.nodes and v in graph.nodes:
                        graph.add_edge(u, v)
    except Exception as e:
        print(f"[ERROR] Gagal membaca GFA: {e}")
        return None
                    
    print(f"   -> Total Node Awal: {len(graph.nodes)}")
    return graph

# ==========================================
# BAGIAN 3: IDENTIFIKASI LABEL (BLAST)
# ==========================================

def identify_anchors(graph, label_db, temp_dir):
    print("[CLEANER] Menjalankan BLASTN untuk identifikasi Label...")
    
    if not os.path.exists(temp_dir): os.makedirs(temp_dir)
    
    # 1. Tulis Node ke Fasta Sementara
    query_fasta = os.path.join(temp_dir, "nodes.fasta")
    with open(query_fasta, 'w') as f:
        for name, node in graph.nodes.items():
            f.write(f">{name}\n{node.sequence}\n")
            
    # 2. Jalankan BLASTN
    # -evalue 1e-5: Cukup longgar untuk ONT yang error-prone
    # -outfmt 6: Tabular output
    # -dust no: PENTING! Jangan filter low-complexity (AT-rich)
    blast_out = os.path.join(temp_dir, "blast_hits.txt")
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-subject", label_db,
        "-outfmt", "6 qseqid length pident", 
        "-evalue", "1e-5",
        "-dust", "no" 
    ]
    
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("[ERROR] Gagal menjalankan BLASTN. Pastikan ncbi-blast+ terinstall.")
        return False

    # 3. Parsing Hasil
    # Kita butuh set node yang punya hit signifikan
    anchors_found = 0
    with open(blast_out, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts: continue
            
            node_name = parts[0]
            aln_len = int(parts[1])
            
            if node_name in graph.nodes:
                # Filter tambahan: Panjang alignment minimal 50bp agar tidak random noise
                if aln_len > 50:
                    node = graph.nodes[node_name]
                    node.is_anchor = True
                    node.keep = True
                    # Simpan hit terbaik
                    if aln_len > node.blast_hit_len:
                        node.blast_hit_len = aln_len
                        
    # Hitung jumlah anchor unik
    anchors_found = sum(1 for n in graph.nodes.values() if n.is_anchor)
    print(f"   -> Ditemukan {anchors_found} Node Jangkar (Anchors).")
    
    return anchors_found > 0

# ==========================================
# BAGIAN 4: LOGIKA SLIMMING (INTI)
# ==========================================

def slim_graph_topology(graph):
    """
    Logika: Hitung rata-rata coverage Anchor, lalu selamatkan tetangga 
    yang terhubung dan memiliki coverage yang 'mirip'.
    """
    print("[CLEANER] Memulai proses Slimming & Neighbor Rescue...")
    
    # 1. Hitung Target Coverage (Baseline)
    anchor_covs = [n.coverage for n in graph.nodes.values() if n.is_anchor]
    
    if not anchor_covs:
        return graph # Tidak ada anchor, kembalikan apa adanya (atau error)
        
    avg_target_cov = sum(anchor_covs) / len(anchor_covs)
    print(f"   -> Target Coverage Organel: {avg_target_cov:.2f}x")
    
    # Ambang batas penyelamatan: 
    # Node tetangga harus punya coverage minimal 10% dari target.
    min_rescue_cov = avg_target_cov * 0.1
    
    # 2. Propagasi Penyelamatan (Neighbor Rescue Loop)
    iteration = 0
    changed = True
    
    while changed:
        changed = False
        iteration += 1
        rescued_in_this_round = 0
        
        for name, node in graph.nodes.items():
            # Jika sudah disimpan, skip
            if node.keep: continue
            
            # Cek Tetangga
            # Apakah ada tetangga yang SUDAH berstatus 'keep=True'?
            neighbors = graph.edges[name]
            connected_to_safe_zone = any(graph.nodes[n].keep for n in neighbors)
            
            if connected_to_safe_zone:
                # Cek Syarat Coverage
                if node.coverage >= min_rescue_cov:
                    node.keep = True
                    changed = True
                    rescued_in_this_round += 1
        
        if rescued_in_this_round > 0:
            print(f"      [Iterasi {iteration}] Menyelamatkan {rescued_in_this_round} node tetangga.")

    # 3. Hapus Sampah
    all_nodes = list(graph.nodes.keys())
    removed_count = 0
    
    for name in all_nodes:
        if not graph.nodes[name].keep:
            graph.remove_node(name)
            removed_count += 1
            
    print(f"   -> Pembersihan Selesai. Membuang {removed_count} node sampah.")
    print(f"   -> Sisa Node Valid: {len(graph.nodes)}")
    
    return graph

# ==========================================
# BAGIAN 5: OUTPUT WRITER
# ==========================================

def write_clean_gfa(graph, original_gfa_path, output_gfa_path):
    print(f"[CLEANER] Menyimpan grafik bersih ke: {output_gfa_path}")
    
    valid_names = set(graph.nodes.keys())
    
    with open(original_gfa_path, 'r') as fin, open(output_gfa_path, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            if not parts: continue
            
            # Tulis Node jika valid
            if parts[0] == 'S':
                if parts[1] in valid_names:
                    fout.write(line)
            
            # Tulis Edge jika kedua ujung valid
            elif parts[0] == 'L':
                if parts[1] in valid_names and parts[3] in valid_names:
                    fout.write(line)

# ==========================================
# BAGIAN 6: WRAPPER & MAIN (UPDATE ARGPARSE)
# ==========================================

def clean_graph(gfa_input, label_db, gfa_output, custom_temp_dir=None):
    """
    Fungsi utama yang dipanggil oleh main.
    UPDATE: Menambahkan handling custom_temp_dir untuk parallel processing.
    """
    # Gunakan temp dir unik (jika custom_temp_dir ada) atau generate dari output path
    if custom_temp_dir:
        temp_dir = custom_temp_dir
    else:
        # Fallback agar tidak bentrok jika dijalankan manual
        temp_dir = os.path.dirname(gfa_output) + "/temp_cleaner_" + os.path.basename(gfa_output)
    
    # 1. Load Graph
    g = parse_flye_gfa(gfa_input)
    if not g: return None
    
    # 2. Labeling (BLAST)
    has_anchors = identify_anchors(g, label_db, temp_dir)
    
    if has_anchors:
        # 3. Slimming (Logic tetangga)
        g = slim_graph_topology(g)
        
        # 4. Save
        write_clean_gfa(g, gfa_input, gfa_output)
        
        # Cleanup (Hanya hapus temp jika berhasil)
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            
        return gfa_output
    else:
        print("[WARNING] Tidak ada jangkar ditemukan. Grafik tidak dibersihkan.")
        if os.path.exists(temp_dir): shutil.rmtree(temp_dir)
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assembly Graph Cleaner (BLAST-based)")
    
    # Arguments
    parser.add_argument("--gfa", required=True, help="Input GFA from Flye")
    parser.add_argument("--ref", required=True, help="Reference Fasta for labeling (Hybrid Bait)")
    parser.add_argument("--out", required=True, help="Output Cleaned GFA")
    parser.add_argument("--temp", required=False, help="Temporary directory (optional)")

    args = parser.parse_args()
    
    # Run
    clean_graph(args.gfa, args.ref, args.out, args.temp)