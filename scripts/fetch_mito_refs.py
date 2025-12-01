import time
from Bio import Entrez
import sys

if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    OUTPUT_FILE = snakemake.output[0]
else:
    Entrez.email = "example@email.com"
    OUTPUT_FILE = "ref_tebu_mito_manual.fasta"

# Target Mitokondria (Sama, karena kita cari kontaminan dari spesies yg sama)
targets = {
    "Primary": ["Saccharum", "Sorghum", "Zea mays"], # Genus saja cukup
}

def search_and_fetch():
    print(f"--- [FETCH-MITO] Target Output: {OUTPUT_FILE} ---")
    unique_ids = set()
    
    for category, organism_list in targets.items():
        for organism in organism_list:
            # Query spesifik MITOKONDRIA (Range size > 200kb)
            term = (f'"{organism}"[Organism] AND (mitochondrion[Title] OR mitochondrial[Title]) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND 200000:10000000[Sequence Length]')
            
            try:
                handle = Entrez.esearch(db="nucleotide", term=term, retmax=5, idtype="acc")
                record = Entrez.read(handle)
                handle.close()
                
                id_list = record["IdList"]
                if id_list:
                    print(f"   [+] {organism}: {len(id_list)} found.")
                    unique_ids.update(id_list)
                time.sleep(0.5)
                
            except Exception as e:
                print(f"   [!] Error {organism}: {e}")

    if not unique_ids:
        # Fallback: Jika gagal fetch mito tebu, ambil jagung/sorghum manual
        print("WARNING: Tidak ada ref mito spesifik. Mencoba Zea mays umum...")
        unique_ids.add("NC_007982.1") # Zea mays mitochondrion ref

    print(f"Downloading {len(unique_ids)} sequences...")
    try:
        net_handle = Entrez.efetch(db="nucleotide", id=list(unique_ids), rettype="fasta", retmode="text")
        data = net_handle.read()
        net_handle.close()
        
        with open(OUTPUT_FILE, "w") as out_f:
            out_f.write(data)
    except Exception as e:
        print(f"Error download: {e}")
        sys.exit(1)

if __name__ == "__main__":
    search_and_fetch()
