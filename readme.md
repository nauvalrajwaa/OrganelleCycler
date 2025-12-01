# OrganelleCycler: Smart De Novo Assembly Pipeline for Plant Organelles

**OrganelleCycler** is a robust, automated Snakemake pipeline designed for the *de novo* assembly of plant organelle genomes (Plastomes and Mitochondria) using **Oxford Nanopore** long reads. 

It specifically addresses common challenges in plant assembly, such as **MTPTs (Mitochondrial Plastid DNA insertions)** and high repetitiveness, by employing a "Smart Filtering" strategy and a multi-assembler competition logic.

---

## 🚀 Key Features

* **automated Reference Fetching:** Automatically downloads taxonomically related references (e.g., *Saccharum*, *Oryza*, *Zea*) from NCBI based on your configuration.
* **Smart Blacklisting:**
    * Uses a **Rough Assembly (Raven)** step to identify large mitochondrial contigs.
    * Filters based on **Graph Circularity**, **Identity**, and **Length** to distinguish true Plastomes from Mitokondria or Nuclear junk.
    * Prevents accidental data loss from MTPTs (Mitochondrial DNA containing Plastid insertions).
* **Dual-Filtering Strategy:**
    1.  **Negative Filter:** Removes reads mapping to the mitochondrial blacklist.
    2.  **Positive Filter (Baiting):** Recruits remaining reads that map to the fetched plastome reference.
* **Multi-Assembler Competition:** Runs **Flye**, **Raven**, and **Canu** in parallel on the cleaned reads.
* **Automated Judging:** A Python script evaluates all assemblies based on circularity, size, and BLAST identity to select the single best consensus sequence.

---

## 🛠️ Prerequisites

* **Linux/Unix** environment (or WSL on Windows).
* **Conda** (Miniconda or Mambaforge) installed.
* **Git** installed.

All bioinformatics tools (Raven, Flye, Canu, Minimap2, BLAST+, etc.) are managed automatically via Snakemake and Conda environments.

---

## 📦 Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/yourusername/OrganelleCycler.git](https://github.com/yourusername/OrganelleCycler.git)
    cd OrganelleCycler
    ```

2.  **Create the Snakemake environment:**
    ```bash
    conda create -n snakemake -c conda-forge -c bioconda snakemake
    conda activate snakemake
    ```

---

## ⚙️ Configuration

### 1. Setup Samples (`config/samples.tsv`)
Map your sample names to their raw Nanopore FASTQ files.
```tsv
sample_id       reads_path
tebu_sample01   data/raw_reads_01.fastq.gz
tebu_sample02   data/raw_reads_02.fastq.gz
````

### 2\. Configure Parameters (`config/config.yaml`)

Edit this file to define your target organism and assembly goals.

  * **`target_assembly`**: Choose `"PLASTOME"`, `"MITO"`, or `"BOTH"`.
  * **`query_taxa`**: Define the taxonomy for reference fetching (e.g., *Saccharum*, *Oryza*).
  * **`plastome_est_size`**: Estimated genome size (e.g., "150k").

-----

## ▶️ Usage

### 1\. Dry Run (Test Configuration)

Check if the pipeline is set up correctly without executing heavy jobs.

```bash
snakemake -n
```

### 2\. Run the Pipeline

Execute the workflow. The `--use-conda` flag is **mandatory** to handle tool dependencies automatically.

```bash
snakemake --cores 8 --use-conda
```

*(Replace `8` with the number of CPU threads you wish to use).*

-----

## 📂 Output Structure

After a successful run, results are stored in the `results/` directory:

```text
results/{sample_id}/
├── 01_rough/            # Initial rough assembly (Raven) & GFA graphs
├── 02_blacklist/        # Smart blacklist of mitochondrial contigs
├── 03_filtered_reads/   # Final clean reads used for assembly (recruited)
├── 04_assemblies/       # Draft assemblies from Flye, Raven, and Canu
├── 05_polished/         # Polished assemblies (using Medaka/Racon)
└── 07_best_candidate/   # 🏆 FINAL OUTPUT
    ├── best_assembly.fasta      # The chosen best plastome assembly
    └── selection_report.txt     # detailed score report of the competition
```

-----

## 🧠 Logic & Workflow

1.  **Reference Fetching:** Python script queries NCBI Entrez to download related genomes.
2.  **Rough Assembly:** **Raven** generates a raw graph to identify large junk contigs (mitochondria).
3.  **Smart Filtering:**
      * *Is it circular?* → Keep.
      * *Is it \>200kb & Mito-like?* → Blacklist.
4.  **Read Recruitment:** Raw reads are mapped against the blacklist (remove matches) and then mapped against the Plastome reference (keep matches).
5.  **Assembly & Polishing:** The recruited reads are assembled by **Flye**, **Raven**, and **Canu**. Resulting drafts are polished with **Minipolish**.
6.  **The Judge:** Assemblies are scored:
      * **+500 points** for circularity.
      * **Score** based on BLAST Identity & Coverage against the reference.
      * The highest scoring assembly is saved as `best_assembly.fasta`.

-----

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

-----

**Developed for Plant Organelle Assembly**
*Author: [Your Name]*

```
```