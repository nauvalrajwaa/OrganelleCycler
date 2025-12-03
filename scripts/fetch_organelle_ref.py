import time
import sys
from Bio import Entrez

# --- SETUP DARI SNAKEMAKE ---
if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    
    # OUTPUTS (Wajib disesuaikan dengan rule Snakefile)
    OUT_FASTA = snakemake.output.fasta
    OUT_GBK   = snakemake.output.gbk
    
    # PARAMETERS
    ORGANELLE_QUERY = snakemake.params.search_term
    MIN_LEN = snakemake.params.min_len
    MAX_LEN = snakemake.params.max_len
    
    # TAXONOMY CONFIG
    TAXA_CONFIG = snakemake.config["query_taxa"]
else:
    sys.exit("Script ini wajib dijalankan via Snakemake.")

def search_and_fetch():
    print(f"--- [NCBI FETCH] Mode: Super Bait & GenBank ---")
    print(f"--- Query: {ORGANELLE_QUERY} ({MIN_LEN}-{MAX_LEN} bp) ---")
    
    unique_ids = set()
    
    # 1. SEARCHING PHASE (Mengumpulkan ID)
    for category, organism_list in TAXA_CONFIG.items():
        # Prioritas jumlah download
        if category.lower() == "primary": max_ret = 20
        elif category.lower() == "secondary": max_ret = 5
        else: max_ret = 2 
            
        print(f"\nScanning Category: {category} (Max per organism: {max_ret})")
        
        for organism in organism_list:
            # Query yang presisi
            term = (f'"{organism}"[Organism] AND ({ORGANELLE_QUERY}) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
            
            try:
                # Search ID saja dulu
                handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_ret, idtype="acc")
                record = Entrez.read(handle)
                handle.close()
                
                id_list = record["IdList"]
                if id_list:
                    print(f"   [+] {organism}: Found {len(id_list)} IDs -> {', '.join(id_list[:3])}...")
                    unique_ids.update(id_list)
                else:
                    print(f"   [-] {organism}: No matching reference found.")
                
                time.sleep(0.3) # Biar tidak diblokir NCBI
                
            except Exception as e:
                print(f"   [!] Error querying {organism}: {e}")

    # 2. FALLBACK PHASE (Jika kosong total)
    if not unique_ids:
        print("\n[WARNING] Tidak ada referensi spesifik ditemukan dari list config.")
        if "mitochondrion" in ORGANELLE_QUERY.lower():
            print("   -> Fallback: Menggunakan Zea mays Mitochondrion (NC_007982.1)")
            unique_ids.add("NC_007982.1")
        else:
             print("   -> Fallback: Menggunakan Saccharum officinarum Plastome (NC_006084.1)")
             unique_ids.add("NC_006084.1")

    # 3. DOWNLOADING PHASE (Batch Download)
    num_ids = len(unique_ids)
    print(f"\n--- Downloading {num_ids} sequences (FASTA & GBK) ---")
    
    id_list = list(unique_ids)
    
    try:
        # A. Download FASTA (Untuk Baiting & MitoHiFi -f)
        print(f"   > Fetching FASTA format...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        fasta_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_FASTA, "w") as out_f:
            out_f.write(fasta_data)
            
        # B. Download GENBANK (Untuk MitoHiFi -g Anotasi/Rotasi)
        print(f"   > Fetching GENBANK format...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
        gbk_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_GBK, "w") as out_g:
            out_g.write(gbk_data)
            
        print("   > Success! References saved.")
        
    except Exception as e:
        print(f"[FATAL ERROR] Download failed: {e}")
        # Hapus file output jika ada biar Snakemake tau ini gagal
        if os.path.exists(OUT_FASTA): os.remove(OUT_FASTA)
        if os.path.exists(OUT_GBK): os.remove(OUT_GBK)
        sys.exit(1)

if __name__ == "__main__":
    search_and_fetch()