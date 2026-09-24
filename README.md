# Exploring the Microbial COSMOS III
### From a single-end mock sample to quality-checked, taxonomically classified MAGs — using Galaxy

This tutorial walks through a complete metagenomic binning pipeline in [Galaxy](https://usegalaxy.eu/), starting from a **single-end** mock/demo sample, through quality control, taxonomic profiling, assembly, binning, refinement, quality assessment, and taxonomic classification of the final MAGs.

Adapted from the official [GTN Binning of metagenomic sequencing data tutorial](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/metagenomics-binning/tutorial.html), with parameters adjusted for single-end input.

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [About the Sample](#about-the-sample)
3. [Workflow at a Glance](#workflow-at-a-glance)
4. [Step 1 — FastQC](#step-1--fastqc)
5. [Step 2 — KneadData](#step-2--kneaddata-single-end-mode)
6. [Step 3 — MetaPHlAn](#step-3--metaphlan-on-the-original-raw-fastqgz)
7. [Step 4 — Assembly (MEGAHIT)](#step-4--megahit-assembly-single-end-mode)
8. [Step 5 — Read Mapping (Bowtie2)](#step-5--bowtie2-single-end-mode)
9. [Step 6 — Binning](#step-6--binning)
10. [Step 7 — Bin Refinement (DAS_Tool)](#step-7--bin-refinement-das_tool)
11. [Step 8 — Quality Assessment & Taxonomy](#step-8--quality-assessment--taxonomy)
12. [Step 9 — KBase: Phylogenetic Classification](#step-9--kbase-phylogenetic-classification)
13. [Step 10 — Pathogenwatch: Pathogen Identification](#step-10--pathogenwatch-pathogen-identification)


---

## Prerequisites

- A [Galaxy](https://usegalaxy.eu/) account
- A single-end fastq.gz sample (e.g. a mock/demo community)
- Basic familiarity with the Galaxy interface (uploading data, running tools)

**Before the workshop, please complete the following setup:**

1. **Create a new history** in Galaxy — this keeps all your workshop files organized in one place (click the **+** icon at the top of the History panel, or **Data → Histories → Create new → name it WORKSHOP_2026**).
2. **Upload your data** — drag and drop your `.fastq.gz` file into Galaxy, or use the **Upload Data** button, and wait until it turns green (finished) in your history.
3. **Run FastQC** — use the tool search bar on the left-hand tool panel, search for "**fastqc**", select the tool, and run it on your uploaded file (see [Step 1](#step-1--fastqc) below for full details).
4. **Run KneadData** — search for "**kneaddata**" in the tool search bar and run it on the same raw file (see [Step 2](#step-2--kneaddata-single-end-mode) below for full details).

Completing these two tool runs (FastQC and KneadData) ahead of time means we can dive straight into the more interesting parts of the pipeline during the workshop itself.

---

## About the Sample

Some mock/demo samples — particularly synthetically generated, reference-based mock communities — are **single-end only** and **cannot be converted to paired-end**. Pairing reflects a real physical sequencing process (two ends of the same DNA fragment being sequenced); a single tiled or single-end fastq file has no second read to pair with.

**A tell-tale sign of a synthetically tiled mock sample:** headers where each successive read starts exactly one base later than the previous one, combined with uniform maximum quality scores across every base. If you see this pattern, treat the file as single-end and do not attempt to force-pair or interleave it.

Because of this, this tutorial uses **single-end-specific settings** at every relevant step, and **does not use MetaWRAP** — MetaWRAP hard-requires a paired dataset collection and will not accept single-end reads under any workaround.

---

## Workflow at a Glance

```
raw_mock.fastq.gz
   │
   ├──> FastQC ──> check stats
   │
   ├──> KneadData (single-end) ──> cleaned reads
   │                                    │
   │                                    └──> MEGAHIT (single-end) ──> assembly (Contigs)
   │                                                                        │
   └──> MetaPHlAn (single-end, on the ORIGINAL raw fastq.gz)                │
        [taxonomic profile — independent branch, for comparison]           │
                                                                            ▼
                                                              Bowtie2 (single-end) ──> BAM
                                                                            │
                                                                     Samtools sort
                                                                            │
                                                          Calculate contig depths
                                                                            │
                                      ┌─────────────────────────────────────┼─────────────────────────────────────┐
                                      ▼                                     ▼                                     ▼
                                 MetaBAT2                                MaxBin2                              CONCOCT
                                      └─────────────────────────────────────┼─────────────────────────────────────┘
                                                                            ▼
                                                                        DAS_Tool
                                                                            │
                                                                            ▼
                                                                  CheckM  →  GTDB-Tk (local)
                                                                            │
                                                                            ▼
                                              KBase: Import Assembly → Build AssemblySet
                                                                            │
                                                                            ▼
                                                              KBase: GTDB-Tk Classify
                                                              (taxonomy + phylogenetic tree)
                                                                            │
                                                                            ▼
                                                                 Identify candidate pathogen
                                                                            │
                                                                            ▼
                                                            Pathogenwatch (AMR + typing)
```

---

## Step 1 — FastQC

Run on the **raw** fastq.gz file — no special settings needed. This establishes a baseline before any processing.

| Parameter | Value |
|---|---|
| Input | raw fastq.gz |

> A synthetically tiled mock sample will typically show **uniform maximum quality scores** and **no adapter content** — quite different from real sequencing data. Worth pointing out to attendees as a sign this is simulated, not instrument-generated, data.

---


## Step 2 — KneadData (single-end mode)

Run on the **raw fastq.gz**, even though the mock sample is already clean — this step is included for demonstration purposes, so attendees see the full standard pipeline.

| Parameter | Value |
|---|---|
| Input | **Single-end** — select the one raw fastq.gz file |
| Reference database | your host-decontamination DB (Human Genome) |

Expect minimal trimming/removal, since the sample is already clean — a fine teaching point in itself ("here's what KneadData reports when there's nothing to clean"). Use the main cleaned-reads output going forward for assembly and taxnomic profiling.

---

## Step 3 — MetaPHlAn (on the ORIGINAL raw fastq.gz)

Run MetaPHlAn on the **original, untouched raw fastq.gz** — this is a separate, parallel branch from the assembly path. MetaPHlAn profiles taxonomic composition directly from reads using marker genes, giving a quick community snapshot to later compare against what the assembly/binning pipeline recovers.

| Parameter | Value |
|---|---|
| Input | **Single-end** mode |
| Input file | the **original raw** fastq.gz (not KneadData's output) |

Keep this output aside — you'll compare it against the GTDB-Tk classification of your final MAGs in Step 8.

---

## Step 4 — MEGAHIT (assembly, single-end mode)

| Parameter | Value |
|---|---|
| Input reads | KneadData's cleaned single-end output |
| Library type | Single-end |

MEGAHIT natively supports single-end assembly (unlike metaSPAdes, which expects paired input). Use the **Contigs** output (not Scaffolds) for all downstream steps.

---

## Step 5 — Bowtie2 (single-end mode)

**Tool:** `Bowtie2`

| Parameter | Value |
|---|---|
| Is this single or paired library | **Single-end** |
| FASTQ file |  single-end reads |
| Reference genome source | Use a genome from the history and build index |
| Select reference genome | MEGAHIT assembly Contigs |
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

## Step 6 — Binning

### Shared step: Calculate contig depths

**Tool:** `Calculate contig depths` 

| Parameter | Value |
|---|---|
| Mode to process BAM files | One by one |
| Sorted bam files | Samtools sort output |
| Select a reference genome? | No |

This depth file is reused by MetaBAT2 and MaxBin2 below.

### 6.1 MetaBAT2

**Tool:** `MetaBAT2`

| Parameter | Value |
|---|---|
| Fasta file containing contigs | assembly Contigs |


### 6.2 MaxBin2

**Tool:** `MaxBin2`

| Parameter | Value |
|---|---|
| Contig file | assembly Contigs |
| Assembly type used to generate contig(s) | Assembly of sample(s) one by one (individual assembly) |
| Input type | Abundances |
| **Abundance file** | the same depth matrix (**not** the assembly fasta) |

### 6.3 CONCOCT

CONCOCT needs its own multi-step chain, since it clusters cut-up contig fragments rather than whole contigs.

**Step A — Cut up contigs**

| Parameter | Value |
|---|---|
| Fasta contigs file | assembly Contigs |
| Chunk size | 10000 |
| Overlap size | 0 |
| **Concatenate final part to last contig?** | **Yes** ⚠️ (critical — see Troubleshooting) |
| Output bed file? | Yes |

**Step B — Generate the input coverage table**

| Parameter | Value |
|---|---|
| Contigs BEDFile | BED output from Step A |
| Type of assembly | Individual assembly: 1 run per BAM file |
| Sorted BAM file | Samtools sort output |

**Step C — Run CONCOCT**

| Parameter | Value |
|---|---|
| Coverage file | output of Step B |
| Composition file with sequences | cut-up fasta from Step A |
| Read length for coverage | your actual read length (e.g. 150 — a plain number) |

**Step D — Merge cut clusters**

| Parameter | Value |
|---|---|
| Clusters generated by CONCOCT | the **"Clusters"** output of Step C (not "PCA transformed clusters") |

**Step E — Extract a fasta file**

| Parameter | Value |
|---|---|
| Original contig file | assembly Contigs by MEGAHIT (original, full-length — not cut-up) |
| CONCOCT clusters | merged clusters from Step D |

> ⚠️ CONCOCT's Gaussian clustering model often struggles on small/uneven-coverage demo datasets — expect many low-completeness bins. This is documented, expected behavior, not a sign of misconfiguration.

> ⚠️ **Do not use MetaWRAP** anywhere in this pipeline — it hard-requires a `collection_type="paired"` input and will not accept single-end reads under any workaround.

---

## Step 7 — Bin Refinement (DAS_Tool)

First, convert each binner's fasta bins into a contig-to-bin table:

**Tool:** `Converts genome bins in fasta format` — run once per binner, selecting each binner's final fasta bin output (MetaBAT2's "Bin sequences", MaxBin2's "Bins", CONCOCT's "Extract a fasta file" output).

Then run:

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

## Step 8 — Quality Assessment & Taxonomy

### CheckM

**Tool:** `CheckM lineage_wf`

| Parameter | Value |
|---|---|
| Bins | DAS_Tool "Bins" output |

Gives completeness %, contamination %, and strain heterogeneity per bin.

### GTDB-Tk

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

**Finally, compare methods:** put your **MetaPHlAn read-based taxonomic profile** (Step 4) side by side with the **GTDB-Tk classification of your assembled/binned MAGs** (this step) — a good discussion point on how read-based vs. assembly-based taxonomic methods can agree or diverge.

---

## Step 9 — KBase: MAGs Phylogenetic Classification

As an alternative (or complement) to running GTDB-Tk locally, you can run the same classification inside [KBase](https://www.kbase.us/) — a free, browser-based platform that also gives you a proper phylogenetic tree placement for your MAGs, not just a summary table.

### 9.1 — Upload your high-quality MAGs to KBase

1. Take the **ranked/high-quality bins** from your DAS_Tool output (Step 8) — typically the ones with the best completeness/contamination scores from Step 9's CheckM report.
2. In a KBase Narrative, go to **Upload** → **Staging Area**, and upload each MAG fasta file.
3. For each fasta file, run:

   **App:** `Import FASTA as Assembly from Staging`

   | Parameter | Value |
   |---|---|
   | Staging file | your MAG fasta file |
   | Assembly name | a clear name per MAG (e.g. `mag_1_assembly`) |
   | Type | draft isolate (or metagenome, depending on how you want it labeled) |
   | Min contig length | 0 (don't filter further — you've already refined these bins) |

Repeat for each MAG you want to classify.

### 9.2 — Build an AssemblySet

GTDB-Tk in KBase does **not** accept individual Assembly objects directly — they must first be grouped into an **AssemblySet** (this avoids running the app inefficiently, once per genome).

**App:** `Build AssemblySet`

| Parameter | Value |
|---|---|
| Assemblies | select all your imported MAG assemblies |
| Output AssemblySet name | e.g. `workshop_mags_assemblyset` |

### 9.3 — Run GTDB-Tk Classify

**App:** `GTDB-Tk Classify` (`kb_gtdbtk/run_kb_gtdbtk_classify_wf`)

| Parameter | Value |
|---|---|
| Input object | your AssemblySet from Step 10.2 |
| Reference data | keep the default (currently GTDB R07-RS207 / R08-RS214, class-level subtrees — lighter on memory than the full tree) |

**Output:** a taxonomic classification per MAG (domain → species, as far as confidently resolvable) plus a phylogenetic placement, viewable directly in the Narrative.

> Just like with local GTDB-Tk, low-completeness MAGs may only resolve to a shallow taxonomic level (e.g. phylum or genus rather than species).

### 9.4 — Identify the pathogen

Review the GTDB-Tk classification output for each MAG. If a MAG classifies to a genus/species with known pathogenic members (e.g. *Salmonella*, *Escherichia*, *Klebsiella*, *Mycobacterium*, *Streptococcus* etc.), that's your candidate for the next step. Cross-check the identification against what you'd expect from your **MetaPHlAn** read-based profile (Step 4) as a sanity check — the two methods should broadly agree.

---

## Step 10 — Pathogenwatch: Pathogen Identification

Once you've identified a MAG of interest as a likely pathogen, upload its fasta file to [Pathogenwatch](https://pathogen.watch/) for pathogen-specific genomic analysis (AMR gene detection, MLST typing, and species-specific typing schemes where available).

1. Download the specific MAG's fasta file from KBase (or directly from your Galaxy DAS_Tool "Bins" output).
2. Go to [pathogen.watch/upload](https://pathogen.watch/upload) and create/sign in to an account.
3. Select the correct **organism/species scheme** matching your GTDB-Tk classification (Pathogenwatch supports specific pathogens — e.g. *Streptococcus*, *Salmonella*, *E. coli*, *Klebsiella*, *M. tuberculosis*, *Neisseria*, and others — check their supported organism list, since unsupported organisms won't have a dedicated typing scheme).
4. Upload the fasta file.
5. Review the results: AMR gene predictions, sequence typing (MLST/cgMLST where supported), and clustering against Pathogenwatch's global genome collection.

> This final step ties the whole pipeline together: raw reads → assembly → binning → refined MAG → taxonomic identity → pathogen-specific genomic surveillance, mirroring a real genomic epidemiology workflow.

---




---

## References

- [GTN: Binning of metagenomic sequencing data](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/metagenomics-binning/tutorial.html)
- [GTN: Metagenomics Assembly tutorial](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/metagenomics-assembly/tutorial.html)
- [GTN: MAG building tutorial](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/mags-building/tutorial.html)
- Sieber, C.M.K., et al. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology*, 3(7), 836-843.
