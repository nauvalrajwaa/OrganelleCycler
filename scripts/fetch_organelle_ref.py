import time
import sys
from Bio import Entrez

# --- SETUP DARI SNAKEMAKE ---
if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    OUTPUT_FILE = snakemake.output[0]
    
    # Ambil parameter dinamis dari 'params' di Snakefile
    ORGANELLE_QUERY = snakemake.params.search_term
    MIN_LEN = snakemake.params.min_len
    MAX_LEN = snakemake.params.max_len
    
    # Ambil list taksonomi dari config
    TAXA_CONFIG = snakemake.config["query_taxa"]
else:
    sys.exit("Script ini wajib dijalankan via Snakemake.")

def search_and_fetch():
    print(f"--- [NCBI FETCH] Target: {OUTPUT_FILE} ---")
    print(f"--- Query Type: {ORGANELLE_QUERY} ---")
    
    unique_ids = set()
    
    # Loop melalui Primary, Secondary, Outgroup yang ada di Config
    for category, organism_list in TAXA_CONFIG.items():
        # Tentukan prioritas jumlah download
        if category.lower() == "primary": max_ret = 20
        elif category.lower() == "secondary": max_ret = 5
        else: max_ret = 2 
            
        print(f"\nScanning Category: {category} (Max: {max_ret})")
        
        for organism in organism_list:
            # Query Universal
            term = (f'"{organism}"[Organism] AND ({ORGANELLE_QUERY}) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
            
            try:
                handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_ret, idtype="acc")
                record = Entrez.read(handle)
                handle.close()
                
                id_list = record["IdList"]
                if id_list:
                    print(f"   [+] {organism}: {len(id_list)} found.")
                    unique_ids.update(id_list)
                else:
                    print(f"   [-] {organism}: None.")
                
                time.sleep(0.3) 
                
            except Exception as e:
                print(f"   [!] Error {organism}: {e}")

    if not unique_ids:
        print("WARNING: Tidak ada referensi spesifik ditemukan.")
        # Fallback sederhana agar pipeline tidak crash total
        if "mitochondrion" in ORGANELLE_QUERY:
            print("Mencoba fallback Mito Zea mays...")
            unique_ids.add("NC_007982.1")
        else:
             print("Mencoba fallback Plastome Saccharum...")
             unique_ids.add("NC_006084.1")

    print(f"\nDownloading {len(unique_ids)} sequences...")
    try:
        net_handle = Entrez.efetch(db="nucleotide", id=list(unique_ids), rettype="fasta", retmode="text")
        data = net_handle.read()
        net_handle.close()
        
        with open(OUTPUT_FILE, "w") as out_f:
            out_f.write(data)
        print("Download Selesai.")
        
    except Exception as e:
        print(f"Error download: {e}")
        sys.exit(1)

if __name__ == "__main__":
    search_and_fetch()