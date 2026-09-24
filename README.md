# Metagenomics to Pathogen Genomics Workshop
### From raw reads to quality-checked, taxonomically classified MAGs — using Galaxy

This tutorial walks through a complete metagenomic binning pipeline in [Galaxy](https://usegalaxy.eu/), from raw paired-end sequencing reads to final Metagenome-Assembled Genomes (MAGs) with quality scores and taxonomic classification.

Based on and adapted from the official [GTN Binning of metagenomic sequencing data tutorial](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/metagenomics-binning/tutorial.html).

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Getting Demo Data](#getting-demo-data)
3. [Quality Control](#1-quality-control-kneaddata)
4. [Assembly](#2-assembly)
5. [Read Mapping](#3-read-mapping-bowtie2)
6. [Binning](#4-binning)
   - [MetaBAT2](#41-metabat2)
   - [MaxBin2](#42-maxbin2)
   - [CONCOCT](#43-concoct)
7. [Bin Refinement](#5-bin-refinement)
8. [Quality Assessment](#6-quality-assessment-checkm)
9. [Taxonomic Classification](#7-taxonomic-classification-gtdb-tk)
10. [Troubleshooting](#troubleshooting--common-pitfalls)

---

## Prerequisites

- A [Galaxy](https://usegalaxy.eu/) account
- Paired-end fastq reads (forward + reverse)
- Basic familiarity with the Galaxy interface (uploading data, running tools, building collections)

---

## Getting Demo Data

For a workshop-scale demo (fast runtime, small download), use the official GTN test dataset:

```
https://zenodo.org/records/17661262/files/reads_forward.fastqsanger.gz
https://zenodo.org/records/17661262/files/reads_reverse.fastqsanger.gz
```

Upload both files to a new Galaxy history, then build a **paired dataset collection**:
1. Select both files in your history
2. Choose **"Build Dataset Pair"**
3. Assign forward/reverse correctly, name it, and create

> ⚠️ If you're using your own real data instead of the demo set, and it's large, **subsample first** to keep runtimes workshop-friendly:
> ```bash
> seqtk sample -s100 reads_1.fastq.gz 100000 | gzip > sub_R1.fastq.gz
> seqtk sample -s100 reads_2.fastq.gz 100000 | gzip > sub_R2.fastq.gz
> ```
> Use the **same seed** (`-s100`) for both files so read pairs stay synced. Verify both outputs have identical read counts before proceeding.

---

## 1. Quality Control (KneadData)

**Tool:** `KneadData`

| Parameter | Value |
|---|---|
| Input | your paired collection |

**Outputs to keep:** use the **"Paired output reads"** output — this is your final, QC'd, host-decontaminated reads. Ignore "Trimmed paired reads" and "Repeats removed paired reads" (intermediate steps).

---

## 2. Assembly

**Tool:** `MEGAHIT` (recommended for speed) or `metaSPAdes` (higher quality, slower)

| Parameter | Value |
|---|---|
| Input reads | KneadData "Paired output reads" |

Use the **Contigs** output (not Scaffolds) for all downstream binning steps — scaffolds contain gap characters (`N`s) that interfere with coverage calculation and clustering.

---

## 3. Read Mapping (Bowtie2)

**Tool:** `Bowtie2`

| Parameter | Value |
|---|---|
| Is this single or paired library | Paired-end |
| FASTQ Paired Dataset | your reads collection |
| Reference genome source | Use a genome from the history and build index |
| Select reference genome | your assembly Contigs |
| Set read groups information? | Do not set |
| Select analysis mode | 1: Default setting only |
| Save mapping statistics to history | Yes |

Then sort the output:

**Tool:** `Samtools sort`

| Parameter | Value |
|---|---|
| BAM file | Bowtie2 output |
| Primary sort key | coordinate |

---

## 4. Binning

### Shared step: Calculate contig depths

**Tool:** `Calculate contig depths` (jgi_summarize_bam_contig_depths)

| Parameter | Value |
|---|---|
| Mode to process BAM files | One by one |
| Sorted bam files | Samtools sort output |
| Select a reference genome? | No |

This depth file is reused by both MetaBAT2 and MaxBin2 below.

### 4.1 MetaBAT2

**Tool:** `MetaBAT2`

| Parameter | Value |
|---|---|
| Fasta file containing contigs | assembly Contigs |
| **Use a base coverage depth file?** | **Yes** → select the depth matrix from above |
| Minimum size of a contig for binning | 2500 (lower to 500–1000 for small/fragmented demo assemblies) |
| Minimum size of a bin as the output | 200000 (lower to ~10,000–25,000 for small demo assemblies) |

### 4.2 MaxBin2

**Tool:** `MaxBin2`

| Parameter | Value |
|---|---|
| Contig file | assembly Contigs |
| Assembly type used to generate contig(s) | Assembly of sample(s) one by one (individual assembly) |
| Input type | Abundances |
| **Abundance file** | the same depth matrix (**not** the assembly fasta — see [Troubleshooting](#troubleshooting--common-pitfalls)) |
| Outputs → all four toggles | Yes |

### 4.3 CONCOCT

CONCOCT needs its own multi-step chain, since it clusters cut-up contig fragments rather than whole contigs.

**Step 1 — Cut up contigs**

| Parameter | Value |
|---|---|
| Fasta contigs file | assembly Contigs |
| Chunk size | 10000 |
| Overlap size | 0 |
| **Concatenate final part to last contig?** | **Yes** ⚠️ (critical — see Troubleshooting) |
| Output bed file? | Yes |

**Step 2 — Generate the input coverage table**

| Parameter | Value |
|---|---|
| Contigs BEDFile | BED output from Step 1 |
| Type of assembly | Individual assembly: 1 run per BAM file |
| Sorted BAM file | Samtools sort output |

**Step 3 — Run CONCOCT**

| Parameter | Value |
|---|---|
| Coverage file | output of Step 2 |
| Composition file with sequences | cut-up fasta from Step 1 |
| Read length for coverage | your actual read length (e.g. 150 — a plain number, not a placeholder) |

**Step 4 — Merge cut clusters**

| Parameter | Value |
|---|---|
| Clusters generated by CONCOCT | the **"Clusters"** output of Step 3 (not "PCA transformed clusters") |

**Step 5 — Extract a fasta file**

| Parameter | Value |
|---|---|
| Original contig file | assembly Contigs (original, full-length — not cut-up) |
| CONCOCT clusters | merged clusters from Step 4 |

> ⚠️ CONCOCT's Gaussian clustering model often struggles on small/uneven-coverage demo datasets — expect many low-completeness bins. This is documented, expected behavior, not a sign of misconfiguration.

---

## 5. Bin Refinement

Combine the outputs of all three binners into one consensus, non-redundant bin set.

**First, convert each binner's fasta bins into a contig-to-bin table:**

**Tool:** `Converts genome bins in fasta format` — run once per binner, selecting each binner's final fasta bin output (MetaBAT2's "Bin sequences", MaxBin2's "Bins", CONCOCT's "Extract a fasta file" output).

### Option A: Binette

**Tool:** `Build list` — combine the three contig-to-bin tables into one list (Insert Dataset ×3, label each `Index`)

**Tool:** `Binette`

| Parameter | Value |
|---|---|
| Input contig table | Build list output |
| Input contig file | assembly Contigs |
| Database | cached database → CheckM2 diamond DB |
| Minimum completeness | 0 (demo data) |

### Option B: DAS_Tool (alternative to Binette)

**Tool:** `DAS Tool for genome-resolved metagenomics`

| Parameter | Value |
|---|---|
| Contig sequences | assembly Contigs |
| Bins (repeat ×3) | each binner's contig-to-bin table, with a label (`metabat2`, `maxbin2`, `concoct`) |

**Key outputs:**
- **"Bins"** → final refined MAG fasta files
- **"Quality and completeness estimates of input bin sets"** → per-binner quality before refinement (good for a before/after comparison)
- **"Summary of output bins"** → final bin stats

---

## 6. Quality Assessment (CheckM)

**Tool:** `CheckM lineage_wf`

| Parameter | Value |
|---|---|
| Bins | your final refined bins (Binette or DAS_Tool output) |

Gives completeness %, contamination %, and strain heterogeneity per bin.

---

## 7. Taxonomic Classification (GTDB-Tk)

Galaxy's GTDB-Tk tool requires a specific pre-cached database release that may not be available on your instance. If unavailable, run locally instead:

```bash
conda activate gtdbtk-conda

# Check if the reference database is already set up
echo $GTDBTK_DATA_PATH
ls $GTDBTK_DATA_PATH | head

# If empty, download it (large — 70-100+ GB, do this ahead of time)
download-db.sh

# Download your final bins from Galaxy first, then classify
gtdbtk classify_wf \
  --genome_dir path/to/your/bins/ \
  --out_dir gtdbtk_output/ \
  --extension fasta \
  --cpus 4

# Results
cat gtdbtk_output/gtdbtk.bac120.summary.tsv
```

> Low-completeness bins may only classify to a shallow taxonomic level (e.g. phylum rather than species) — this is expected with partial genomes.

---

## Troubleshooting / Common Pitfalls

These are real issues encountered while building this pipeline — listed here so you don't have to rediscover them.

| Symptom | Cause | Fix |
|---|---|---|
| MetaBAT2 runs "fine" but produces 0 bins | `Use a base coverage depth file?` left on default "No", or bin size threshold too high | Set to Yes + select depth file; lower "Minimum size of a bin" for small demo assemblies |
| MaxBin2: `Failed to get abundance information` | Assembly fasta accidentally used as the abundance file | Use the **depth matrix** output, not the contigs fasta |
| MaxBin2 job held/killed by scheduler | Exceeded requested memory (marker-gene search is resource-heavy) | Retry with a smaller dataset, or drop MaxBin2 if repeatedly failing on limited infrastructure |
| CONCOCT: `TypeError: Invalid value ... for dtype 'float64'` | Composition/coverage files mismatched, or wrong file (e.g. a BED file) plugged in as "Coverage file" | Regenerate the coverage table end-to-end using a matching BED + BAM; double-check you're selecting the actual coverage table, not the BED file |
| CONCOCT: `can only concatenate str (not "float") to str` | `Concatenate final part to last contig?` left on default "No", creating tiny leftover fragments with unreliable coverage | Set to **Yes** in "Cut up contigs" |
| CONCOCT bins mostly 0% completeness | Known limitation — CONCOCT's clustering struggles on small/low-coverage assemblies | Expected; not a configuration error. Document as a discussion point. |
| GTDB-Tk (Galaxy): "No options available... requires release 232" | No cached database on the Galaxy instance | Run locally via conda instead, or ask a Galaxy admin to install it |
| GTDB-Tk (local): `AttributeError: module 'numpy' has no attribute 'bool'` | numpy ≥1.24 removed the deprecated `np.bool` alias that GTDB-Tk 2.1.1 still uses | `pip install "numpy<1.24" --force-reinstall` in the gtdbtk conda environment |
| Paired collection auto-pairing fails | Filenames like `reads_forward`/`reads_reverse` don't match Galaxy's default `_1`/`_2` pattern | Use "configure auto-pairing" to add a custom pattern, or pair manually |

---

## Summary: Expected Results on the Demo Dataset

| Binner | Bins produced | Best completeness | Notes |
|---|---|---|---|
| MetaBAT2 | 1 | ~15.7% | Clean, no contamination |
| MaxBin2 | 2 | 10.5% / 6.9% | No contamination |
| CONCOCT | ~10 | ~15.7% (1 bin); rest ~0% | Struggles on small/uneven-coverage assemblies |
| Refined (Binette/DAS_Tool) | 2 | ~15.7% | Best bins combined across all three tools |

Recovering even 1-2 low-to-moderate completeness MAGs from a small demo dataset is a **successful, expected outcome** — real-world analyses with full sequencing depth recover more numerous and more complete genomes.

---

## References

- [GTN: Binning of metagenomic sequencing data](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/metagenomics-binning/tutorial.html)
- [GTN: Metagenomics Assembly tutorial](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/metagenomics-assembly/tutorial.html)
- [GTN: MAG building tutorial](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/mags-building/tutorial.html)
- Sieber, C.M.K., et al. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology*, 3(7), 836-843.
