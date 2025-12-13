import sys
from Bio import SeqIO

def extract_genes_from_gbk(gbk_path, output_fasta):
    """
    Mengekstrak CDS, tRNA, dan rRNA dari file GenBank ke format FASTA.
    """
    print(f"[EXTRACTOR] Membaca GenBank: {gbk_path}")
    count = 0
    
    with open(output_fasta, "w") as out_f:
        # Parse GenBank
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type in ["CDS", "tRNA", "rRNA"]:
                    # Ambil nama gen
                    gene_name = "unknown"
                    if "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                    elif "product" in feature.qualifiers:
                        gene_name = feature.qualifiers["product"][0]
                    
                    # Ambil sequence
                    try:
                        seq = feature.extract(record.seq)
                        # Tulis ke Fasta
                        out_f.write(f">{gene_name}_from_{record.id} type={feature.type}\n{seq}\n")
                        count += 1
                    except Exception as e:
                        print(f"[WARNING] Gagal mengekstrak fitur {gene_name}: {e}")
                        continue

    print(f"[EXTRACTOR] Berhasil mengekstrak {count} gen ke {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_genes.py <input.gb> <output.fasta>")
        sys.exit(1)
    
    gbk_input = sys.argv[1]
    fasta_output = sys.argv[2]
    
    extract_genes_from_gbk(gbk_input, fasta_output)