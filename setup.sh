#!/usr/bin/env bash
set -euo pipefail

# ===============================
# Setup working directory
# ===============================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p resources
mkdir -p dryrun

# ===============================
# Conda environment setup
# ===============================
if ! command -v conda &> /dev/null; then
    echo "❌ conda tidak ditemukan. Pastikan conda sudah terinstal."
    exit 1
fi

echo "▶ Membuat conda environment..."
conda env create -f environment.yml || echo "⚠ Environment mungkin sudah ada, lanjut..."

ENV_NAME=$(grep '^name:' environment.yml | awk '{print $2}')

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# ===============================
# Download GetOrganelleDB
# ===============================
echo "▶ Download GetOrganelleDB..."
curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz \
  | tar zx --strip-components=1 -C resources

echo "📂 GetOrganelleDB tersedia di: resources/"

# ===============================
# Download FASTQ (dry-run example)
# ===============================
echo "▶ Download example FASTQ (dry-run, Arabidopsis)..."

wget -qO- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR217/003/ERR2173373/ERR2173373.fastq.gz \
  | zcat \
  | head -n 400000 \
  | gzip > dryrun/arabidopsis_dryrun.fastq.gz

echo "📂 Dry-run FASTQ tersedia di: dryrun/arabidopsis_dryrun.fastq.gz"

# ===============================
# Notes
# ===============================
cat << EOF

====================================
✅ Setup selesai
====================================

📌 NOTES:
- File FASTQ di folder dryrun/ hanya berisi subset data
  (head 400000 baris) dan digunakan untuk:
  ▶ example run
  ▶ testing pipeline
  ▶ debugging workflow

- GetOrganelleDB berada di folder:
  ▶ resources/

====================================
▶ Cara menjalankan pipeline
====================================

snakemake --cores 64 --use-conda --use-singularity

====================================

🐍 Conda env aktif: ${ENV_NAME}
🚀 Siap menjalankan workflow
EOF
