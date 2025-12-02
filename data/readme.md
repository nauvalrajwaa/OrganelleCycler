# 📂 Data Directory: Input Instructions

This directory serves as the storage location for your **Raw Sequencing Data**. The *OrganelleCycler* pipeline looks for input files here to begin the assembly process.

## ✅ Required Data Format

To ensure the pipeline runs smoothly, your data must meet the following criteria:

1.  **Sequencing Platform:** Oxford Nanopore Technologies (ONT).
    * *Note:* This pipeline is optimized for **Long Reads**. Short reads (Illumina) are not supported.
2.  **File Format:** FASTQ.
3.  **Compression:** Gzipped (`.fastq.gz`) is highly recommended to save space, though uncompressed (`.fastq`) is accepted.
4.  **Content:** Total Genomic DNA (reads can contain a mix of Nucleus, Mitochondria, and Plastid DNA). The pipeline will automatically filter and extract the specific organelle reads.

---

## 🔗 How to Connect Data to the Pipeline

Placing files here is not enough; you must register them in the configuration file.

1.  **Place your file:**
    e.g., `data/tebu_sample_01.fastq.gz`

2.  **Update `config/samples.tsv`:**
    Open the `samples.tsv` file in the main directory and ensure the path matches exactly.

    ```tsv
    sample_id       reads_path
    sample_01       data/tebu_sample_01.fastq.gz
    sample_02       data/tebu_sample_02.fastq.gz
    ```

---

## 🧬 Recommended Specifications

For the best assembly results (Circular & High Quality):

* **Basecalling:** High-accuracy basecalling (HAC or SUP) using Guppy or Dorado is preferred.
* **Read Length:** N50 > 5kb is ideal. Ultra-long reads help resolve the Inverted Repeats (IR) regions in plastomes.
* **Coverage:**
    * Organelle genomes usually have much higher coverage than the nuclear genome.
    * **Minimum:** ~200-500 MB of total raw data is usually sufficient to get >50x coverage for plastomes in most plant species.
    * **Maximum:** If you have >10 GB of data, consider subsampling to speed up the process, as excessive coverage provides diminishing returns for organelle assembly.

---

## 🧪 Testing (Dry Run)

If you do not have your own data yet and want to test the pipeline installation, you can download a small subset of public data (e.g., *Arabidopsis* or *Zea mays*) into this folder.

**Example Command to download test data:**
```bash
# Download a small subset of Zea mays reads
wget -O data/dryrun_test.fastq.gz [https://github.com/nauvalraj/test-data/raw/main/zea_mays_subset.fastq.gz](https://github.com/nauvalraj/test-data/raw/main/zea_mays_subset.fastq.gz)
````

*Remember to update `config/samples.tsv` to point to this test file before running.*

```
```