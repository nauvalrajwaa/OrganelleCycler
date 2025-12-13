#!/usr/bin/env bash
set -euo pipefail

# Buat conda environment dari YAML
echo "▶ Membuat conda environment..."
conda env create -f environment.yml || echo "⚠ Environment mungkin sudah ada, lanjut..."

# Ambil nama environment dari YAML
ENV_NAME=$(grep '^name:' setup.yml | awk '{print $2}')

# Inisialisasi conda untuk bash
source "$(conda info --base)/etc/profile.d/conda.sh"

# Aktifkan environment
conda activate "$ENV_NAME"

# Download & extract GetOrganelleDB ke folder resources
echo "▶ Download GetOrganelleDB..."
curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz \
  | tar zx -C resources

echo "✅ Setup selesai."
echo "📂 Database tersedia di: resources/"
echo "🐍 Conda env aktif: $ENV_NAME"
