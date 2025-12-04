import time
import sys
import os
from Bio import Entrez
from urllib.error import HTTPError, URLError

# --- SETUP DARI SNAKEMAKE ---
if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    # Jika punya API key, pakai. Jika tidak, biarkan.
    if "api_key" in snakemake.config:
        Entrez.api_key = snakemake.config["api_key"]
    
    OUT_FASTA = snakemake.output.fasta
    OUT_GBK   = snakemake.output.gbk
    
    ORGANELLE_QUERY = snakemake.params.search_term
    MIN_LEN = snakemake.params.min_len
    MAX_LEN = snakemake.params.max_len
    TAXA_CONFIG = snakemake.config["query_taxa"]
    
    # Ambil parameter expand, default False
    EXPAND_LINEAGE = snakemake.params.get("expand_lineage", False)
else:
    sys.exit("Script ini wajib dijalankan via Snakemake.")

# --- KONFIGURASI NCBI AGAR TIDAK KENA BANNED ---
Entrez.max_tries = 5             # Coba ulang sampai 5 kali jika gagal
Entrez.sleep_between_tries = 10  # Tunggu 10 detik setiap gagal (Penting!)

def robust_entrez_search(term, retmax):
    """
    Wrapper search yang sangat sabar menghadapi server NCBI.
    """
    attempt = 0
    max_attempts = 5
    
    while attempt < max_attempts:
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax, idtype="acc")
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        
        except (HTTPError, URLError, RuntimeError) as e:
            # Cek jika error 429 (Too Many Requests)
            print(f"   [!] Connection warning (attempt {attempt+1}/{max_attempts}): {e}")
            print("       -> Sleeping for 20 seconds...")
            time.sleep(20) # Jeda panjang biar IP didinginkan
            attempt += 1
            
    print(f"   [X] Failed to search term after {max_attempts} attempts.")
    return []

def get_parent_taxa(organism):
    try:
        # Cari TaxID
        handle = Entrez.esearch(db="taxonomy", term=organism)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']: return []
        tax_id = record['IdList'][0]
        
        # Fetch Detail
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if "Lineage" in records[0]:
            full_lineage = records[0]["Lineage"].split("; ")
            return full_lineage[::-1]
            
    except Exception as e:
        print(f"   [!] Taxonomy lookup warning: {e}")
        return []
    return []

def search_and_fetch():
    mode_str = "Expanded (Kerabat)" if EXPAND_LINEAGE else "Specific (Target)"
    print(f"--- [NCBI FETCH] Mode: {mode_str} | Query: {ORGANELLE_QUERY} ---")
    
    unique_ids = set()
    
    # 1. SEARCHING
    for category, organism_list in TAXA_CONFIG.items():
        is_primary = (category.lower() == "primary")
        
        # Tentukan batas download
        if is_primary: max_ret = 20
        elif category.lower() == "secondary": max_ret = 5
        else: max_ret = 2 
            
        for organism in organism_list:
            # A. Cari Spesies Target
            term = (f'"{organism}"[Organism] AND ({ORGANELLE_QUERY}) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
            
            ids = robust_entrez_search(term, max_ret)
            
            if ids:
                print(f"   [+] {organism}: Found {len(ids)} IDs.")
                unique_ids.update(ids)
            else:
                print(f"   [-] {organism}: No direct match found.")
            
            # Wajib sleep antar request
            time.sleep(2) 

            # B. Lineage Expansion (Hanya jika Mode Expanded & Primary)
            if is_primary and EXPAND_LINEAGE:
                print(f"   [Lineage] Expanding search for: {organism}...")
                parents = get_parent_taxa(organism)
                
                # Ambil 2 level di atas (Genus, Family)
                target_parents = parents[:2] 
                
                for parent in target_parents:
                    parent_term = (f'"{parent}"[Organism] AND ({ORGANELLE_QUERY}) '
                                   f'AND (complete genome[Title] OR complete sequence[Title]) '
                                   f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
                    
                    parent_ids = robust_entrez_search(parent_term, 5) # Limit 5
                    
                    if parent_ids:
                        new_ids = [x for x in parent_ids if x not in unique_ids]
                        if new_ids:
                            print(f"       [+] Added {len(new_ids)} new IDs from {parent}")
                            unique_ids.update(new_ids)
                    
                    time.sleep(2) # Sleep antar parent

    # 2. FALLBACK MEKANISME (Jika pencarian gagal total)
    if not unique_ids:
        print("\n[WARNING] Tidak ada ID ditemukan. Menggunakan FALLBACK NCBI Standard.")
        if "mitochondrion" in ORGANELLE_QUERY.lower():
            # Zea mays mitochondrion (RefSeq standard)
            fallback_id = "NC_007982.1"
            print(f"   -> Fallback ID: {fallback_id}")
            unique_ids.add(fallback_id)
        else:
             # Saccharum officinarum plastid (RefSeq standard)
             fallback_id = "NC_006084.1"
             print(f"   -> Fallback ID: {fallback_id}")
             unique_ids.add(fallback_id)

    # 3. DOWNLOADING
    num_ids = len(unique_ids)
    print(f"\n--- Downloading {num_ids} sequences ---")
    
    id_list = list(unique_ids)
    
    try:
        # Download FASTA
        print("   > Downloading FASTA...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        fasta_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_FASTA, "w") as out_f:
            out_f.write(fasta_data)
            
        # Download GBK
        print("   > Downloading GBK...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
        gbk_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_GBK, "w") as out_g:
            out_g.write(gbk_data)
            
        print(f"   > Success! Output saved.")
        
    except Exception as e:
        print(f"[FATAL ERROR] Download failed: {e}")
        # Hapus file output jika ada biar Snakemake tau job ini gagal dan bisa diulang
        if os.path.exists(OUT_FASTA): os.remove(OUT_FASTA)
        if os.path.exists(OUT_GBK): os.remove(OUT_GBK)
        sys.exit(1)

if __name__ == "__main__":
    search_and_fetch()