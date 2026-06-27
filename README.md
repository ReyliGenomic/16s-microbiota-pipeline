# 16S rRNA Microbiota Analysis Pipeline — Direct Joining (DJ) Method

A reproducible, end-to-end bioinformatics pipeline for 16S rRNA amplicon sequencing analysis using the **Direct Joining (DJ) method** (Kim et al., 2025). This pipeline processes paired-end reads across multiple batches, performs denoising with DADA2 via QIIME 2, applies batch correction with ConQuR, and runs downstream ecological and differential abundance analysis in R using **microeco**.

---

## Pipeline Overview

```
Raw FASTQ (paired-end)
        │
        ▼
[1] Adapter Removal — fastp (Docker)
        │
        ▼
[2] DJ Concatenation — concatenate_dj.py (Python)
    Forward + NNNNNNNN + RevComp(Reverse)
        │
        ▼
[3] Spacer Removal — clear_strings.py (Python)
        │
        ▼
[4] Manifest Creation — create_manifest_dj.py (Python)
        │
        ▼
[5] QIIME 2: Import → DADA2 Denoising → Taxonomy (SILVA 138)
        │
        ▼
[6] Feature Table Merging (5 batches) + Phylogenetic Tree
        │
        ▼
[7] Batch Correction — ConQuR (R)
        │
        ▼
[8] Alpha Diversity — Rarefaction, Shannon, Chao1, KW, Wilcoxon, LM (microeco)
        │
        ▼
[9] Beta Diversity — PCoA, NMDS (Bray-Curtis / Aitchison) + PERMANOVA (vegan)
        │
        ▼
[10] Differential Abundance — LEfSe + Cladogram (microeco)
```

---

## 🛠️ Requirements

### Software
- Docker (for fastp and QIIME 2)
- Miniconda / conda
- Python 3.9+
- R 4.x+

### Python environment
```bash
conda create -n microbioma python=3.9 -y
conda activate microbioma
conda install -c conda-forge biopython -y
```

### R packages
```r
install.packages(c("ggplot2", "vegan", "compositions", "doParallel",
                   "magrittr", "stringr", "readr", "cowplot"))

# Bioconductor
BiocManager::install(c("phyloseq", "metagenomeSeq", "edgeR"))

# CRAN/GitHub
install.packages(c("qiime2R", "microeco", "file2meco", "ConQuR", "ape"))
```

### QIIME 2 (via Docker)
```bash
docker pull quay.io/qiime2/amplicon:2025.10
```

---

## 📁 Repository Structure

```
16s-microbiota-pipeline/
├── README.md
├── 1_preprocessing/
│   └── fastp_adapter_removal.sh          # Paired-end adapter trimming
├── 2_concatenation/
│   ├── concatenate_dj.py                 # DJ-method read concatenation
│   └── clear_strings.py                  # Remove NNNNNNNN spacers
├── 3_manifest/
│   └── create_manifest_dj.py             # QIIME 2 manifest generator
├── 4_qiime2/
│   └── qiime2_pipeline.sh                # Import → DADA2 → Taxonomy → Merge
├── 5_batch_correction/
│   └── conqur_batch_correction.R and     # ConQuR batch effect removal
|   └── permanova_analysis.R     
└── 6_downstream/
    └── Alpha_Beta_diversity.Rmd          # Alpha/Beta diversity + LEfSe
```

---

## 🚀 Usage

### Step 1 — Adapter Removal

Run fastp inside Docker to trim primers:

```bash
docker run -it --rm \
    --user $(id -u):$(id -g) \
    -v /path/to/your/data/:/data \
    staphb/fastp:latest /bin/bash
```

Execute `1_preprocessing/fastp_adapter_removal.sh` inside the container.
Key parameters: `--trim_front1 17` (R1), `--trim_front2 21` (R2), `--detect_adapter_for_pe`.

---

### Step 2 — Direct Joining Concatenation

```bash
conda activate microbioma
python 2_concatenation/concatenate_dj.py
```

Concatenates reads as: `Forward + NNNNNNNN (spacer) + RevComp(Reverse)` → `joined_reads/`

---

### Step 3 — Spacer Removal

```bash
python 2_concatenation/clear_strings.py
```

Removes `NNNNNNNN` sequence spacers and `IIIIIIII` quality spacers → `joined_reads_clean/`

---

### Step 4 — Manifest File

```bash
python 3_manifest/create_manifest_dj.py joined_reads_clean/
```

Generates `manifest.tsv` with Docker-compatible `/data/...` paths.

---

### Step 5 — QIIME 2 Pipeline

```bash
docker run -it --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd):/data -w /data \
  quay.io/qiime2/amplicon:2025.10 /bin/bash
```

Run `4_qiime2/qiime2_pipeline.sh`, which covers:

1. **Import** — `SingleEndFastqManifestPhred33V2`
2. **Denoising** — DADA2 (`--p-trunc-len 440`, `--p-max-ee 5.0`, 16 threads)
3. **Feature table merging** — across 5 batches
4. **Phylogenetic tree** — MAFFT + FastTree
5. **Taxonomic classification** — SILVA 138 99% NB classifier (16 jobs)
6. **Filtering** — removes mitochondria, chloroplasts; min-frequency 10, min-samples 5

---

### Step 6 — Batch Correction (ConQuR)

Open `5_batch_correction/conqur_batch_correction.R` in RStudio and run all code.

**Inputs:**
- `final_table_clean.qza` — filtered ASV table from QIIME 2
- `taxonomy.qza` — taxonomic classifications
- `sample-metadata.tsv` — sample metadata

**What it does:**
- Imports QIIME 2 artifacts via `qiime2R` into `phyloseq`
- Applies **ConQuR** using `Lot` as the batch variable
- Covariates: `Sx`, `agey`, `BMI`, `glu`, `SBP`, `chol`, `Class..SBP_CTRL`
- Evaluates correction with PCoA plots (Bray-Curtis and Aitchison) before/after
- Saves corrected tables: `feature_table_corrected_L3.rds`, `feature_table_corrected_L5.rds`

---

### Step 7 — Downstream Analysis (Alpha, Beta, LEfSe)

Open `6_downstream/Alpha_Beta_diversity.Rmd` in RStudio and run all chunks.

#### 7a. Data Import into microeco

```r
# Loads batch-corrected ASV table, taxonomy, metadata, phylogenetic tree
meco <- microtable$new(
  otu_table   = asv_table,
  sample_table = metadata,
  tax_table   = tax_fixed
)
```

Groups analyzed: `Normotenso` (Controlado / No controlado / Normotenso), `Sx` (M / F)

---

#### 7b. Alpha Diversity

- **Rarefaction depth:** 5,000 reads per sample
- **Metrics:** Shannon, Chao1, observed richness
- **Statistical tests:**
  - Kruskal-Wallis (`KW`) and Dunn post-hoc (`KW_dunn`)
  - Wilcoxon test per group
  - Linear model: `Shannon ~ Normotenso + Sx`
- **Outputs:** rarefaction curve, boxplots, error bar plots, linear model heatmap

```r
t2 <- trans_alpha$new(dataset = amplicon_16S_microtable_rarefy_5000, group = "Grupo")
t2$cal_diff(method = "KW")

t3 <- trans_alpha$new(dataset = amplicon_16S_microtable_rarefy_5000, group = "Normotenso")
t3$cal_diff(method = "wilcox", measure = "Shannon")
t3$plot_alpha(measure = "Shannon")
```

---

#### 7c. Beta Diversity

- **Distances:** Bray-Curtis (rarefied), Aitchison/CLR (non-rarefied)
- **Ordination:** PCoA and NMDS
- **Statistical tests:** PERMANOVA (`adonis2`), beta-dispersion (`betadisper`)

```r
# PCoA — Bray-Curtis
t1 <- trans_beta$new(dataset = amplicon_16S_microtable_rarefy_5000,
                     group = "Normotenso", measure = "bray")
t1$cal_ordination(method = "PCoA")
t1$plot_ordination(plot_color = "Normotenso", plot_shape = "Sx",
                   plot_type = c("point", "ellipse"))

# NMDS — Bray-Curtis
t1$cal_ordination(method = "NMDS")
```

**PERMANOVA results (Aitchison distance, 999 permutations):**

```
                 Df SumOfSqs      R2      F Pr(>F)
batchid           4    29175 0.07760 11.386  0.001 ***
agey              1     8248 0.02194 12.877  0.001 ***
Sx                1    10895 0.02898 17.009  0.001 ***
BMI               1    15598 0.04149 24.351  0.001 ***
glu               1     9434 0.02509 14.727  0.001 ***
SBP               1    24242 0.06448 37.845  0.001 ***
chol              1    10225 0.02720 15.962  0.001 ***
Class..SBP_CTRL   1    10364 0.02757 16.180  0.001 ***
Residual        391   250461 0.66620
Total           402   375956 1.00000
```

---

#### 7d. Differential Abundance — LEfSe

```r
tmp_microtable$filter_taxa(rel_abund = 0.0001)  # filter low-abundance taxa

t1 <- trans_diff$new(
  dataset = tmp_microtable,
  method  = "lefse",
  group   = "Normotenso",
  alpha   = 0.05,
  p_adjust_method = "none"
)

t1$plot_diff_bar(use_number = 1:30, group_order = c("No controlado", "Controlado", "Normotenso"))
t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5)
```

Groups compared: `Controlado` vs `No controlado` vs `Normotenso`  
Taxonomic level: **Genus**

---

## 📊 Methods Summary

| Step | Tool | Version / Notes |
|------|------|-----------------|
| Adapter trimming | fastp | 1.0.1, Docker |
| DJ concatenation | Custom Python | Biopython |
| Denoising | DADA2 via QIIME 2 | QIIME 2 2025.10, trunc-len 440 |
| Taxonomy | SILVA 138 NB classifier | 99% identity |
| Batch correction | ConQuR | Lot as batch variable |
| Alpha diversity | microeco | Shannon, Chao1; KW, Wilcoxon, LM |
| Beta diversity | microeco + vegan | Bray-Curtis, Aitchison; PCoA, NMDS, PERMANOVA |
| Differential abundance | LEfSe (microeco) | Genus level, α = 0.05 |
| Phylogenetic tree | MAFFT + FastTree | via QIIME 2 |

---

## 📖 References

- Kim, K.S., Noh, J., Kim, BS. et al. (2025). Refining microbiome diversity analysis by concatenating and integrating dual 16S rRNA amplicon reads. npj Biofilms Microbiomes 11, 57.
- Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.
- Bolyen, E. et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37, 852–857.
- Quast, C. et al. (2013). The SILVA ribososomal RNA gene database project. *Nucleic Acids Research*, 41, D590–D596.
- Gu, W. et al. (2022). ConQuR: Batch effects removal for microbiome data via regression. *Nature Communications*.
- Liu, C. et al. (2021). microeco: An R package for data mining in microbial community ecology. *FEMS Microbiology Ecology*, 97(2).

---

## 👤 Author

**Reyli Sanchez**  
Intern research — Statistical Genomics and Population Health Laboratory, Laboratorio Internacional de Investigación sobre el Genoma Humano, UNAM  
[LinkedIn](www.linkedin.com/in/reyli-sanchez-3b932a219) · [GitHub](https://github.com/ReyliGenomic)
