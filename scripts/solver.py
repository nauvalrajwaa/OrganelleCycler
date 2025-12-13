import sys
import os
import shutil
from collections import defaultdict

# ==========================================
# BAGIAN 1: STRUKTUR DATA & UTILS
# ==========================================

# Tabel transisi untuk Reverse Complement DNA
TRANS_TABLE = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")

def reverse_complement(seq):
    """Membalik urutan DNA (A->T, G->C, dst)"""
    return seq.translate(TRANS_TABLE)[::-1]

class Node:
    def __init__(self, name, length, coverage, sequence):
        self.name = name
        self.length = length
        self.coverage = coverage
        self.sequence = sequence
        self.copy_number = 1   # Default: 1x (Single Copy)
        self.visited_count = 0 # Counter saat penelusuran

# ==========================================
# BAGIAN 2: LOGIKA PARSING (INPUT)
# ==========================================

def parse_flye_gfa(gfa_path):
    print(f"[SOLVER] [TRACE] Membaca grafik GFA: {gfa_path}")
    nodes = {}
    edges = defaultdict(list) 

    try:
        with open(gfa_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if not parts: continue

                # --- Parsing Node (Segment) ---
                if parts[0] == 'S':
                    name = parts[1]
                    seq = parts[2]
                    length = len(seq)
                    cov = 0.0
                    for tag in parts[3:]:
                        if tag.startswith('dp:f:'):
                            cov = float(tag.split(':')[-1])
                            break
                    
                    # Filter noise: Abaikan node sangat pendek (<500bp)
                    if length > 500:
                        nodes[name] = Node(name, length, cov, seq)

                # --- Parsing Edge (Link) ---
                elif parts[0] == 'L':
                    u, u_ori, v, v_ori = parts[1], parts[2], parts[3], parts[4]
                    if u in nodes and v in nodes:
                        edges[u][u_ori].append((v, v_ori))
                        rev_u = '-' if u_ori == '+' else '+'
                        rev_v = '-' if v_ori == '+' else '+'
                        edges[v][rev_v].append((u, rev_u))
    except Exception as e:
        print(f"[ERROR] Gagal membaca GFA: {e}")
        return None, None

    print(f"   -> [STATS] Grafik dimuat. Total Node: {len(nodes)}")
    return nodes, edges

# ==========================================
# BAGIAN 3: LOGIKA MATEMATIKA (VERBOSE TABLE)
# ==========================================

def estimate_multiplicity(nodes):
    print("\n[SOLVER] [LOGIC] Analisis Multiplisitas (Copy Number):")
    
    total_len = 0
    total_cov_len = 0
    
    for n in nodes.values():
        total_len += n.length
        total_cov_len += (n.length * n.coverage)
    
    if total_len == 0: return None
    avg_cov = total_cov_len / total_len
    
    print(f"   [STATS] Weighted Avg Coverage: {avg_cov:.2f}x")
    
    # Ambang batas 1.6x rata-rata biasanya indikator kuat untuk IR
    ir_threshold = avg_cov * 1.6
    
    # Cetak Tabel Keputusan
    print("   [TABLE] Detil Node & Keputusan:")
    print("   Name        | Len   | Cov    | Decision")
    print("   " + "-"*45)
    
    for name, node in nodes.items():
        decision = "SINGLE (1x)"
        if node.coverage >= ir_threshold:
            node.copy_number = 2
            decision = "REPEAT (2x) [IR]"
        else:
            node.copy_number = 1
            
        print(f"   {name:<11} | {node.length:<5} | {node.coverage:<6.1f} | {decision}")

    return avg_cov

# ==========================================
# BAGIAN 4: LOGIKA SOLVER (SIRKULAR & LINEAR WITH TRACE)
# ==========================================

def find_path(nodes, edges, mode="circular"):
    print(f"\n[SOLVER] [LOGIC] Memulai Pencarian Jalur -> Mode: {mode.upper()}")

    # 1. Cari Titik Start (Harus LSC: Panjang & Single Copy)
    start_node = None
    max_len = 0
    
    for name, n in nodes.items():
        if n.copy_number == 1 and n.length > max_len:
            max_len = n.length
            start_node = name

    if not start_node:
        print("   [ERROR] Tidak menemukan kandidat LSC (Start Node). Grafik mungkin fragmentasi.")
        return None
    
    print(f"   [START] Titik awal: {start_node} (+)")

    # Reset visit count
    for n in nodes.values(): n.visited_count = 0
    
    best_path = []
    max_path_len = 0

    # 2. Algoritma Backtracking (DFS) dengan Visual Trace
    def dfs(current_u, current_ori, path_stack, current_len_bp, depth):
        nonlocal best_path, max_path_len
        
        # Indentasi visual untuk log
        indent = "      " + "| " * depth
        
        node_obj = nodes[current_u]
        
        # Uncomment baris ini jika ingin melihat setiap langkah kunjungan (sangat verbose)
        # print(f"{indent}-> Visit: {current_u}({current_ori}) [{node_obj.visited_count}/{node_obj.copy_number}]")

        # A. Cek Kuota Kunjungan
        if node_obj.visited_count >= node_obj.copy_number:
            # print(f"{indent}   X Kuota Habis.")
            return False
        
        # B. Ambil Tiket & Catat Path
        node_obj.visited_count += 1
        seq_fragment = node_obj.sequence if current_ori == '+' else reverse_complement(node_obj.sequence)
        path_stack.append((current_u, current_ori, seq_fragment))
        current_len_bp += len(seq_fragment)
        
        # LOGIKA LINEAR: Update jika rekor terpecahkan
        if mode == "linear":
            if current_len_bp > max_path_len:
                max_path_len = current_len_bp
                best_path = list(path_stack)
                # print(f"{indent}   ! New Record Linear: {current_len_bp} bp")

        # C. Cek Kemenangan (Sirkularitas)
        if mode == "circular" and len(path_stack) > 1:
            start_u, start_ori, _ = path_stack[0]
            if current_u == start_u and current_ori == start_ori:
                print(f"{indent}   *** LOOP DETECTED! Sirkularitas Terkonfirmasi ***")
                best_path = list(path_stack)
                return True

        # D. Cari Tetangga (Next Step)
        found_next = False
        if current_u in edges and current_ori in edges[current_u]:
            neighbors = edges[current_u][current_ori]
            # Heuristic: Prioritaskan tetangga dengan coverage tertinggi
            neighbors.sort(key=lambda x: nodes[x[0]].coverage, reverse=True)
            
            for v_name, v_ori in neighbors:
                if dfs(v_name, v_ori, path_stack, current_len_bp, depth + 1):
                    if mode == "circular": return True 
                    found_next = True
        
        if not found_next and mode == "circular":
             pass # Dead end
        
        # E. Backtrack
        node_obj.visited_count -= 1
        path_stack.pop()
        return False

    # Jalankan DFS
    dummy_stack = []
    found = dfs(start_node, '+', dummy_stack, 0, 0)
    
    if mode == "circular":
        if found:
            print(f"   [RESULT] Sirkular: YES. Path Length: {len(best_path)} nodes.")
            best_path.pop() # Buang node terakhir (duplikat start)
            return best_path
        print("   [RESULT] Sirkular: NO (Gagal menemukan loop tertutup).")
        return None
    else:
        print(f"   [RESULT] Linear: YES. Max Length: {max_path_len} bp.")
        return best_path if len(best_path) > 0 else None

# ==========================================
# BAGIAN 5: OUTPUT GENERATOR
# ==========================================

def save_fasta(path_data, output_file, header_tag="Assembly"):
    print(f"[SOLVER] Menyimpan hasil ke {output_file}")
    
    full_seq = ""
    header_parts = []
    
    for name, ori, seq in path_data:
        full_seq += seq
        header_parts.append(f"{name}({ori})")
    
    header = f">{header_tag} path={'_'.join(header_parts)} length={len(full_seq)}"
    
    with open(output_file, 'w') as f:
        f.write(f"{header}\n{full_seq}\n")
    
    print("[SUCCESS] Selesai.")

# ==========================================
# BAGIAN 6: WRAPPER FUNCTION (INTEGRASI)
# ==========================================

def solve_graph(gfa_input, fasta_output):
    """
    Fungsi wrapper utama yang dipanggil oleh main_ont.py
    Mengadaptasi logika 'Fallback': Coba Sirkular dulu, jika gagal coba Linear.
    """
    # 1. Parse
    nodes_dict, edges_dict = parse_flye_gfa(gfa_input)
    if not nodes_dict:
        print("[ERROR] Grafik kosong atau tidak valid.")
        return False

    # 2. Analyze
    estimate_multiplicity(nodes_dict)
    
    # 3. Solve Strategy 1: CIRCULAR (Prioritas Utama)
    print("\n" + "-"*30)
    print(" STRATEGI 1: MENCARI LINGKARAN")
    print("-" * 30)
    path_result = find_path(nodes_dict, edges_dict, mode="circular")
    
    if path_result:
        print("\n[SOLVER] [DECISION] ==> BERHASIL SIRKULAR")
        save_fasta(path_result, fasta_output, header_tag="Circular_Organelle")
        return True
    
    # 4. Solve Strategy 2: LINEAR / SCAFFOLD (Fallback)
    print("\n" + "-"*30)
    print(" STRATEGI 2: FALLBACK KE LINEAR")
    print("-" * 30)
    print("[WARNING] Gagal membentuk lingkaran sempurna. Mencoba mengekstrak scaffold terpanjang...")
    
    path_result = find_path(nodes_dict, edges_dict, mode="linear")
    
    if path_result:
        print("\n[SOLVER] [DECISION] ==> BERHASIL LINEAR (SCAFFOLD)")
        # Beri nama output berbeda agar user tahu ini bukan sirkular
        base, ext = os.path.splitext(fasta_output)
        linear_output = f"{base}_linear_scaffold{ext}"
        save_fasta(path_result, linear_output, header_tag="Linear_Scaffold_Incomplete")
        
        # Copy ke output utama juga agar pipeline tidak crash, tapi user harus cek
        shutil.copyfile(linear_output, fasta_output)
        return True
    
    print("\n[SOLVER] [DECISION] ==> GAGAL TOTAL")
    return False

# ==========================================
# MAIN EXECUTION (STANDALONE)
# ==========================================

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ont_solver.py <input_graph.gfa> <output.fasta>")
        sys.exit(1)
        
    input_gfa = sys.argv[1]
    output_fasta = sys.argv[2]
    
    solve_graph(input_gfa, output_fasta)