import time
from Bio import Entrez
import sys

# Setup komunikasi dengan Snakemake atau Manual
if 'snakemake' in globals():
    Entrez.email = snakemake.config["email"]
    OUTPUT_FILE = snakemake.output[0]
else:
    Entrez.email = "example@email.com" # Ganti untuk test manual
    OUTPUT_FILE = "ref_tebu_plastome_manual.fasta"

# Target Taksonomi (Bisa diedit jika project berubah tanaman)
targets = {
    "Primary": ["Saccharum officinarum", "Saccharum spontaneum"],
    "Secondary": ["Miscanthus", "Erianthus", "Sorghum bicolor"],
    "Outgroup": ["Zea mays"]
}

def search_and_fetch():
    print(f"--- [FETCH-PLASTOME] Target Output: {OUTPUT_FILE} ---")
    unique_ids = set()
    
    for category, organism_list in targets.items():
        max_ret = 20 if category == "Primary" else 5
        
        for organism in organism_list:
            # Query spesifik Plastome
            term = (f'"{organism}"[Organism] AND (chloroplast[Title] OR plastid[Title]) '
                    f'AND (complete genome[Title] OR complete sequence[Title]) '
                    f'AND 100000:250000[Sequence Length]')
            
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
                time.sleep(0.5) 
                
            except Exception as e:
                print(f"   [!] Error {organism}: {e}")

    if not unique_ids:
        print("FATAL: Tidak ada referensi ditemukan!")
        sys.exit(1)

    print(f"Downloading {len(unique_ids)} sequences...")
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
