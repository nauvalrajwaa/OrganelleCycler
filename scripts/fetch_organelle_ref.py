import time
import sys
import os
from Bio import Entrez

# --- SETUP DARI SNAKEMAKE ---
if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    
    # OUTPUTS
    OUT_FASTA = snakemake.output.fasta
    OUT_GBK   = snakemake.output.gbk
    
    # PARAMETERS
    ORGANELLE_QUERY = snakemake.params.search_term
    MIN_LEN = snakemake.params.min_len
    MAX_LEN = snakemake.params.max_len
    
    # TAXONOMY CONFIG
    TAXA_CONFIG = snakemake.config["query_taxa"]
    
    # --- PARAMETER BARU: SAKLAR LINEAGE ---
    # Default False jika tidak didefinisikan di Snakefile
    EXPAND_LINEAGE = snakemake.params.get("expand_lineage", False) 
else:
    sys.exit("Script ini wajib dijalankan via Snakemake.")

def get_parent_taxa(organism):
    """
    Mengambil daftar parent taxonomy (Genus, Family, dst) dari NCBI Taxonomy.
    """
    try:
        # 1. Cari TaxID
        handle = Entrez.esearch(db="taxonomy", term=organism)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            return []
            
        tax_id = record['IdList'][0]
        
        # 2. Fetch Detail Lineage
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # 3. Parse Lineage
        if "Lineage" in records[0]:
            full_lineage = records[0]["Lineage"].split("; ")
            return full_lineage[::-1] # Dibalik agar Genus/Family ada di index awal
            
    except Exception as e:
        print(f"   [!] Gagal mengambil lineage untuk {organism}: {e}")
        return []
    return []

def search_ncbi(term, max_ret):
    """Helper function untuk search"""
    try:
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_ret, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"   [!] Error searching term '{term}': {e}")
        return []

def search_and_fetch():
    mode_str = "Expanded (Primary + Lineage)" if EXPAND_LINEAGE else "Specific (Primary Only)"
    print(f"--- [NCBI FETCH] Mode: {mode_str} ---")
    print(f"--- Query: {ORGANELLE_QUERY} ({MIN_LEN}-{MAX_LEN} bp) ---")
    
    unique_ids = set()
    
    # 1. SEARCHING PHASE
    for category, organism_list in TAXA_CONFIG.items():
        is_primary = (category.lower() == "primary")
        
        # Config Max Return
        if is_primary: max_ret = 20
        elif category.lower() == "secondary": max_ret = 5
        else: max_ret = 2 
            
        print(f"\nScanning Category: {category} (Max: {max_ret})")
        
        for organism in organism_list:
            # A. Cari Spesies Target (Selalu dijalankan)
            term = (f'"{organism}"[Organism] AND ({ORGANELLE_QUERY}) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
            
            ids = search_ncbi(term, max_ret)
            if ids:
                print(f"   [+] {organism}: Found {len(ids)} IDs.")
                unique_ids.update(ids)
            else:
                print(f"   [-] {organism}: No direct match.")
            
            time.sleep(0.3) 

            # B. [KHUSUS PRIMARY & JIKA SAKLAR ON] Cari Lineage (Genus/Family) 
            if is_primary and EXPAND_LINEAGE:
                print(f"   [Lineage Search] Mengambil kerabat dekat untuk Primary: {organism}...")
                parents = get_parent_taxa(organism)
                
                # Ambil 2 level di atasnya saja
                target_parents = parents[:2] 
                
                for parent in target_parents:
                    print(f"       -> Searching parent taxon: {parent}")
                    parent_term = (f'"{parent}"[Organism] AND ({ORGANELLE_QUERY}) '
                                   f'AND (complete genome[Title] OR complete sequence[Title]) '
                                   f'AND {MIN_LEN}:{MAX_LEN}[Sequence Length]')
                    
                    parent_ids = search_ncbi(parent_term, max_ret=5)
                    
                    if parent_ids:
                        new_ids = [x for x in parent_ids if x not in unique_ids]
                        if new_ids:
                            print(f"       [+] Added {len(new_ids)} new IDs from {parent}")
                            unique_ids.update(new_ids)
                        else:
                            print(f"       [.] {parent} found hits but already in list.")
                    time.sleep(0.3)
            elif is_primary and not EXPAND_LINEAGE:
                 print(f"   [Lineage Search] Skipped (Mode Specific).")

    # 2. FALLBACK PHASE
    if not unique_ids:
        print("\n[WARNING] Tidak ada referensi spesifik ditemukan sama sekali.")
        if "mitochondrion" in ORGANELLE_QUERY.lower():
            print("   -> Fallback: Menggunakan Zea mays Mitochondrion (NC_007982.1)")
            unique_ids.add("NC_007982.1")
        else:
             print("   -> Fallback: Menggunakan Saccharum officinarum Plastome (NC_006084.1)")
             unique_ids.add("NC_006084.1")

    # 3. DOWNLOADING PHASE
    num_ids = len(unique_ids)
    print(f"\n--- Downloading {num_ids} unique sequences (FASTA & GBK) ---")
    
    if num_ids == 0:
        print("[ERROR] No IDs to download.")
        sys.exit(1)

    id_list = list(unique_ids)
    
    try:
        # A. Download FASTA
        print(f"   > Fetching FASTA format...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        fasta_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_FASTA, "w") as out_f:
            out_f.write(fasta_data)
            
        # B. Download GENBANK
        print(f"   > Fetching GENBANK format...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
        gbk_data = net_handle.read()
        net_handle.close()
        
        with open(OUT_GBK, "w") as out_g:
            out_g.write(gbk_data)
            
        print(f"   > Success! References saved to {OUT_FASTA} and {OUT_GBK}")
        
    except Exception as e:
        print(f"[FATAL ERROR] Download failed: {e}")
        if os.path.exists(OUT_FASTA): os.remove(OUT_FASTA)
        if os.path.exists(OUT_GBK): os.remove(OUT_GBK)
        sys.exit(1)

if __name__ == "__main__":
    search_and_fetch()