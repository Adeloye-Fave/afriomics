# AfriOmics 🌍
## Free, Open-Source Multi-Omics Microbiome Pipeline for African Research

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Free to Use](https://img.shields.io/badge/Free-Yes-brightgreen)](https://github.com/afriomics/pipeline)
[![Africa Centred](https://img.shields.io/badge/Context-Africa--Centred-gold)](https://github.com/afriomics/pipeline)

---

## What is AfriOmics?

AfriOmics is a complete, free, Africa-contextualised multi-omics microbiome pipeline integrating:

| Module | Omics Layer | Biological Question |
|--------|-------------|---------------------|
| 🧬 Metagenomics | Shotgun DNA sequencing | **Who is there?** (species composition) |
| 🔬 Metatranscriptomics | RNA sequencing | **What are they actively doing?** |
| ⚗️ Metabolomics | LC-MS/MS mass spectrometry | **What are they producing?** |
| 🔗 Integration | MOFA+ · DIABLO · Networks | **How do they interact?** |
| 🦠 Disease Modelling | Random Forest · ROC | **What predicts disease?** |

Results are benchmarked against curated African reference microbiome datasets (AWI-Gen, H3Africa) across six body sites, and linked to endemic African infectious disease models (malaria, HIV, TB, schistosomiasis).

---

## Quick Start

### Option A — No-Code GUI (Easiest)
```bash
# Install R and Shiny
Rscript -e "install.packages(c('shiny','bslib','ggplot2','dplyr','vegan','plotly','DT'))"

# Launch the GUI
Rscript -e "shiny::runApp('shiny/')"
# Open http://localhost:3838 in your browser
# Upload your MetaPhlAn4 profile and metadata → click Run Analysis
```

### Option B — Jupyter Notebooks (Step-by-step)
```bash
# Install dependencies
pip install pandas numpy matplotlib seaborn plotly scipy scikit-learn networkx pyvis

# Open notebooks in order:
jupyter notebook notebooks/01_metagenomics_walkthrough.ipynb
jupyter notebook notebooks/02_metatranscriptomics_walkthrough.ipynb
jupyter notebook notebooks/03_metabolomics_walkthrough.ipynb
jupyter notebook notebooks/04_integration_disease_modelling.ipynb
```

### Option C — Full Snakemake Pipeline (Advanced)
```bash
# Clone repository
git clone https://github.com/afriomics/pipeline
cd afriomics

# Install conda environments (one-time setup, ~30 min)
conda env create -f envs/environments.yaml

# Edit config
nano config/config_complete.yaml   # Set database paths
nano config/samples.tsv            # Add your sample information

# Dry run (see what will execute)
snakemake --use-conda --cores 16 -s workflow/Snakefile_master -n

# Run everything
snakemake --use-conda --cores 16 -s workflow/Snakefile_master

# Run on HPC (SLURM)
snakemake --use-conda --cores 64 --cluster "sbatch -c {threads} --mem={resources.max_memory_mb}M" \
    -s workflow/Snakefile_master
```

### Option D — Google Colab (No local installation)
```
Open any notebook in notebooks/ and upload to Google Colab.
Runtime → Run All → Upload your data when prompted.
```

---

## Pipeline Architecture

```
Raw Data Input
├── FASTQ reads (metagenomics)
├── FASTQ reads (metatranscriptomics)
├── mzML files (metabolomics)
└── Metadata TSV

          ▼
┌─────────────────────────────────────────────────────────────┐
│                    MODULE 1: METAGENOMICS                   │
│  FastQC → Trimmomatic → Bowtie2 (host removal)             │
│  → Kraken2 + Bracken (taxonomy)                            │
│  → MetaPhlAn4 (species profiles)                           │
│  → HUMAnN3 (functional potential)                          │
│  → QIIME2 (diversity)                                      │
│  → African reference comparison (AWI-Gen/H3Africa)         │
└─────────────────────────────────────────────────────────────┘
          ▼
┌─────────────────────────────────────────────────────────────┐
│               MODULE 2: METATRANSCRIPTOMICS                 │
│  SortMeRNA (rRNA removal) → STAR (alignment)               │
│  → featureCounts (quantification)                          │
│  → DESeq2 (differential expression)                        │
│  → HUMAnN3 RNA (active pathways)                          │
│  → DNA:RNA ratio (potential vs active)                     │
└─────────────────────────────────────────────────────────────┘
          ▼
┌─────────────────────────────────────────────────────────────┐
│                   MODULE 3: METABOLOMICS                    │
│  MZmine3 / XCMS (peak detection, alignment)                │
│  → GNPS (spectral library annotation)                      │
│  → SIRIUS + CSI:FingerID (structure prediction)            │
│  → MetaboAnalyst-compatible stats (PCA, volcano)           │
│  → KEGG/MetaCyc pathway enrichment                        │
│  → African diet metabolome context                         │
└─────────────────────────────────────────────────────────────┘
          ▼
┌─────────────────────────────────────────────────────────────┐
│                  MODULE 4: INTEGRATION                      │
│  Data harmonisation (CLR, log-norm, z-score)               │
│  → MOFA+ (unsupervised multi-omics factors)               │
│  → DIABLO / mixOmics (supervised integration)              │
│  → SparCC (taxa ↔ metabolite correlations)                │
│  → Multi-omics network (igraph + visNetwork)               │
│  → Cross-body-site axes (gut-oral, gut-vaginal, gut-lung)  │
│  → African population structure (PERMANOVA)                │
└─────────────────────────────────────────────────────────────┘
          ▼
┌─────────────────────────────────────────────────────────────┐
│                MODULE 5: DISEASE MODELLING                  │
│  Random Forest (cross-validated, multi-omics features)     │
│  → ROC curves per disease class                            │
│  → SHAP / permutation feature importance                   │
│  → Disease-specific multi-omics signatures                 │
│  → Interactive integrated HTML report                      │
└─────────────────────────────────────────────────────────────┘
          ▼
    Final Report: afriomics_integrated_report.html
```

---

## African Microbiome Context

AfriOmics is the first pipeline built specifically to contextualise microbiome results within African population diversity.

### Body Sites Supported
| Site | Key African Health Context |
|------|---------------------------|
| 🦠 Gut | Malaria, HIV, TB, schistosomiasis, cholera |
| 👄 Oral | Periodontal disease → systemic inflammation |
| 🩺 Vaginal | BV, HIV susceptibility, maternal health |
| 👃 Nasal/Respiratory | TB co-infections, respiratory infections |
| 🩹 Skin | Buruli ulcer, cutaneous leishmaniasis, NTDs |
| 🩸 Blood/Systemic | Malaria translocation, sepsis biomarkers |

### Reference Databases
- **AWI-Gen Cohort** — Pan-African genomic + microbiome data (West, East, Southern, Central Africa)
- **H3Africa** — Human Heredity and Health in Africa consortium
- **MicrobiomeDB** — African subset of curated microbiome studies
- **GNPS African Food Metabolome** — Metabolite library enriched for African dietary compounds
- **African Diet Metabolite Database** — Curated markers from fermented foods, traditional staples, medicinal plants

### African Population Variables Handled
- Geographic region (West, East, Southern, Central, North Africa)
- Diet type (traditional, mixed, western) — including fermented food markers
- Urbanisation level (rural, peri-urban, urban)
- Dietary staples (cassava, millet, sorghum, moringa, fermented cereals)

---

## Disease Models

| Disease | Omics Signatures | African Burden |
|---------|-----------------|----------------|
| 🦟 Malaria (P. falciparum) | ↓ butyrate, ↑ kynurenine, ↓ Faecalibacterium | >90% global cases in Africa |
| 🦠 HIV/AIDS | Gut barrier disruption, ↓ SCFA, ↑ LPS | >70% global burden in SSA |
| 🫁 Tuberculosis | Gut-lung axis, kynurenine pathway, dysbiosis | 30% global TB cases in Africa |
| 🪱 Schistosomiasis | Helminth-microbiome crosstalk, immune modulation | Endemic in 42 African countries |
| 💧 Cholera / Diarrhoeal | Colonisation resistance, SCFA, Vibrio interaction | WASH-linked, high burden |
| 🩺 Sepsis | Microbial translocation, blood metabolomics | High ICU mortality in Africa |

---

## All Tools Used (100% Free)

### Metagenomics
| Tool | Version | Purpose | License |
|------|---------|---------|---------|
| FastQC | 0.12.1 | Quality control | GPL |
| MultiQC | 1.21 | QC aggregation | GPL |
| Trimmomatic | 0.39 | Adapter trimming | GPL |
| Bowtie2 | 2.5.3 | Host removal | GPL |
| Kraken2 | 2.1.3 | Taxonomic classification | MIT |
| Bracken | 2.9 | Abundance re-estimation | GPL |
| MetaPhlAn4 | 4.1.0 | Species profiling | MIT |
| HUMAnN3 | 3.9 | Functional profiling | MIT |
| QIIME2 | 2024.2 | Diversity analysis | BSD |

### Metatranscriptomics
| Tool | Purpose | License |
|------|---------|---------|
| SortMeRNA | rRNA removal | LGPL |
| STAR | RNA alignment | MIT |
| featureCounts | Gene quantification | GPL |
| DESeq2 (R) | Differential expression | LGPL |
| clusterProfiler (R) | Pathway enrichment | Artistic-2.0 |

### Metabolomics
| Tool | Purpose | License |
|------|---------|---------|
| MZmine3 | Peak detection | GPL |
| XCMS (R) | Peak alignment | LGPL |
| GNPS | Spectral annotation | Free (web) |
| SIRIUS | Structure prediction | Academic free |
| MetaboAnalyst | Statistical analysis | Free (web/R) |
| MaAsLin2 (R) | Microbiome-metabolite association | MIT |

### Integration & Modelling
| Tool | Purpose | License |
|------|---------|---------|
| MOFA+ (R/Python) | Multi-omics factor analysis | GPL |
| mixOmics/DIABLO (R) | Supervised integration | GPL |
| SparCC (R) | Compositional correlations | BSD |
| igraph (R) | Network analysis | GPL |
| visNetwork (R) | Interactive network viz | MIT |
| randomForest (R) | Disease classification | GPL |
| Snakemake | Workflow management | MIT |
| R Shiny | No-code GUI | GPL |

---

## Repository Structure

```
afriomics/
├── workflow/
│   ├── Snakefile                         # Metagenomics module
│   ├── Snakefile_metatranscriptomics     # Metatranscriptomics module
│   ├── Snakefile_metabolomics            # Metabolomics module
│   ├── Snakefile_integration             # Integration + disease modelling
│   ├── Snakefile_master                  # Master orchestration
│   ├── scripts/
│   │   ├── diversity_analysis.R
│   │   ├── african_context.R
│   │   ├── deseq2_analysis.R
│   │   ├── mofa_analysis.R
│   │   ├── disease_model.R
│   │   ├── build_multiomics_network.R
│   │   ├── cross_body_site.R
│   │   ├── african_population_structure.R
│   │   └── ...
│   └── report/
│       ├── metagenomics_report.Rmd
│       ├── metatranscriptomics_report.Rmd
│       ├── metabolomics_report.Rmd
│       └── integrated_report.Rmd
│
├── config/
│   ├── config_complete.yaml              # Full pipeline configuration
│   ├── samples.tsv                       # Sample sheet template
│   └── mzmine3_batch.xml                 # MZmine3 batch settings
│
├── envs/
│   └── environments.yaml                 # All conda environments
│
├── shiny/
│   └── app.R                             # No-code Shiny GUI
│
├── notebooks/
│   ├── 01_metagenomics_walkthrough.ipynb
│   ├── 02_metatranscriptomics_walkthrough.ipynb
│   ├── 03_metabolomics_walkthrough.ipynb
│   └── 04_integration_disease_modelling.ipynb
│
├── resources/                            # Databases (download separately)
│   ├── indexes/                          # Human genome Bowtie2 index
│   ├── kraken2_db/                       # Kraken2 PlusPF database
│   ├── metaphlan4_db/                    # MetaPhlAn4 marker database
│   ├── humann3/                          # HUMAnN3 databases
│   ├── african_reference/                # AWI-Gen / H3Africa profiles
│   └── sortmerna_db/                     # SortMeRNA rRNA databases
│
├── data/
│   ├── raw/                              # Raw FASTQ files
│   └── metabolomics/                     # Raw mzML files
│
└── results/                              # All outputs (auto-generated)
    ├── taxonomy/
    ├── functional/
    ├── diversity/
    ├── african_context/
    ├── metatranscriptomics/
    ├── metabolomics/
    ├── integration/
    └── report/
```

---

## Sample Sheet Format

Create `config/samples.tsv` with the following columns:

```tsv
sample_id   r1                          r2                          body_site   disease   region          diet_type     urbanisation
AFRI_001    data/raw/AFRI_001_R1.fq.gz  data/raw/AFRI_001_R2.fq.gz  gut         healthy   west_africa     traditional   rural
AFRI_002    data/raw/AFRI_002_R1.fq.gz  data/raw/AFRI_002_R2.fq.gz  gut         malaria   east_africa     mixed         urban
AFRI_003    data/raw/AFRI_003_R1.fq.gz  data/raw/AFRI_003_R2.fq.gz  gut         HIV       southern_africa western       urban
```

**Required columns:** `sample_id`, `r1`, `r2`

**Recommended columns:** `body_site`, `disease`, `region`, `diet_type`, `urbanisation`

---

## Database Setup (One-Time)

```bash
# 1. Human genome (host removal) — 3.1 GB
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/\
     GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
bowtie2-build GRCh38.fa resources/indexes/human_GRCh38

# 2. Kraken2 PlusPF database — 8 GB (standard) or 70 GB (full)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
tar -xzf k2_standard_08gb_20240112.tar.gz -C resources/kraken2_db/

# 3. MetaPhlAn4 database — auto-downloads on first run
metaphlan --install --bowtie2db resources/metaphlan4_db

# 4. HUMAnN3 databases — ~50 GB total
humann_databases --download chocophlan full resources/humann3/nucleotide
humann_databases --download uniref uniref90_diamond resources/humann3/protein

# 5. SortMeRNA rRNA databases — 500 MB
wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.6/database.tar.gz
tar -xzf database.tar.gz -C resources/sortmerna_db/

# 6. African reference profiles (from AWI-Gen / H3Africa)
# Register at: https://h3africa.org/consortium/resources
# Download MetaPhlAn4 profiles and place in:
# resources/african_reference/awigen_metaphlan_profiles.tsv
```

---

## Computing Requirements

| Component | Minimum | Recommended | Cloud Option |
|-----------|---------|-------------|--------------|
| RAM | 8 GB | 64 GB | Google Colab Pro (25 GB) |
| CPU | 4 cores | 16+ cores | Any HPC cluster |
| Storage | 50 GB | 500 GB | Google Drive / institutional |
| OS | Linux/macOS | Linux | Ubuntu 22.04 |
| Python | 3.9+ | 3.10+ | Included in Colab |
| R | 4.2+ | 4.3+ | Included in conda env |

**Low-resource alternative:** Use Galaxy (usegalaxy.eu) — free, cloud-based, no installation needed. Most tools are pre-installed.

---

## Citation

If you use AfriOmics in your research, please cite:

```
AfriOmics: A free, open-source multi-omics pipeline for African microbiome research.
[Preprint in preparation]

Key tool citations:
- MetaPhlAn4: Blanco-Miguez et al. (2023) Nature Methods
- HUMAnN3: Beghini et al. (2021) eLife
- MOFA+: Argelaguet et al. (2020) Genome Biology
- DIABLO: Singh et al. (2019) Bioinformatics
- AWI-Gen: Asiki et al. (2013) European Journal of Epidemiology
```

---

## Contributing

AfriOmics is built for African scientists, by the global microbiome community.

- 🐛 **Bug reports:** Open a GitHub issue
- 🌍 **New African reference data:** Submit a pull request to `resources/african_reference/`
- 🔧 **New disease modules:** Fork and add to `workflow/scripts/`
- 🌐 **Translations:** Help translate documentation to French, Swahili, Hausa, or Portuguese

---

## License

MIT License — free to use, modify, and distribute.

---

*AfriOmics is dedicated to African scientists working under resource constraints. Science should be free.*
