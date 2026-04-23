# GoSeq — GO Enrichment Analysis Pipeline

A lightweight R pipeline for **Gene Ontology (GO) enrichment analysis** of differentially expressed genes (DEGs), built around the Bioconductor [`goseq`](https://bioconductor.org/packages/release/bioc/html/goseq.html) package. The pipeline corrects for gene-length bias, runs enrichment tests across all input files in batch, and produces both an annotated Excel results file and a publication-ready bubble plot for each comparison.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input File Format](#input-file-format)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
- [Supported Organisms](#supported-organisms)
- [Notes](#notes)

---

## Overview

The core workflow for each input file:

1. Reads a DEG list (gene IDs + binary expression flag).
2. Constructs a probability weighting function (PWF) via `nullp()`, correcting for gene-length bias.
3. Runs `goseq()` to compute over- and under-representation p-values for all GO terms.
4. Annotates results with the specific DEGs driving each enriched category.
5. Saves an Excel results file and a bubble plot (top 20 Biological Process terms, sorted by p-value).

---

## Project Structure

```
GoSeq/
├── Main script.r                  # Main analysis script
│
├── input/                         # DEG input files (tab-delimited or .xlsx)
│   └── expression_*.tab
│
├── output/                        # Generated automatically on first run
│   ├── *_GO_results.xlsx          # Enrichment results with DEG annotation
│   └── *_GO_bubble_plot.pdf       # Bubble plot (top 20 BP terms)
│
├── Go_annotations/                # GO term annotation files (per organism)
│   ├── ArabidopsisGeneID_GO.xlsx
│   ├── E_coliGeneID_GO.xlsx
│   ├── E_coliGeneID_GO_locusbased.xlsx
│   └── ITAG3.0_GeneID_GO.xlsx
│
└── Length/                        # Gene length files (per organism)
    ├── GeneID_length.xlsx
    ├── GeneID_Length_1.xlsx
    ├── GeneID_Length_Arabidopsis.xlsx
    ├── GeneID_Length_e_coli_locus.xlsx
    └── ITAG3.0_GeneID_Length.xlsx
```

---

## Requirements

- R ≥ 4.0
- Bioconductor ≥ 3.12

### R packages

| Package | Source |
|---|---|
| `goseq` | Bioconductor |
| `ggplot2` | CRAN |
| `readxl` | CRAN |
| `dplyr` | CRAN |
| `tidyr` | CRAN |
| `writexl` | CRAN |
| `BiocManager` | CRAN |

---

## Installation

All required packages are installed automatically when the script is first run. To install them manually:

```r
install.packages(c("ggplot2", "readxl", "dplyr", "tidyr", "writexl", "BiocManager"))
BiocManager::install("goseq")
```

---

## Input File Format

Input files go in the `input/` folder. Accepted formats: `.tab`, `.tabular`, or `.xlsx`.

The file must contain at least two columns:

| Column | Description |
|---|---|
| `GeneID` | Gene identifier matching those in the annotation and length files |
| `Expression` | Binary flag: `TRUE` (or `1`) = differentially expressed, `FALSE` (or `0`) = not |

**Example** (`expression_15minTOXIN vs 15minCONTROL_...tab`):

```
GeneID    Expression
b0001     TRUE
b0002     FALSE
b0003     TRUE
...
```

The pipeline processes **all** files in `input/` matching the pattern in one run.

---

## Configuration

At the top of `Main script.r`, set the paths for your organism:

```r
setwd("path/to/GoSeq")                                 # working directory
go_annotation_path <- "Go_annotations/E_coliGeneID_GO.xlsx"  # GO term annotations
gene_length_path   <- "Length/GeneID_length.xlsx"             # gene lengths
input_folder  <- "input"
output_folder <- "output"
```

Swap in the appropriate annotation and length files from `Go_annotations/` and `Length/` when analysing a different organism (see [Supported Organisms](#supported-organisms)).

---

## Usage

Open R or RStudio, then run:

```r
source("Main script.r")
```

The script will:
- Install any missing packages.
- Iterate over every file in `input/`.
- Print progress messages to the console.
- Write results to `output/`.

---

## Output

For each input file two outputs are generated in `output/`:

### 1. `*_GO_results.xlsx`

An Excel workbook with one row per GO category, containing:

| Column | Description |
|---|---|
| `category` | GO term ID |
| `over_represented_pvalue` | P-value for over-representation |
| `under_represented_pvalue` | P-value for under-representation |
| `numDEInCat` | Number of DEGs in this GO category |
| `numInCat` | Total genes in this GO category |
| `term` | GO term description |
| `ontology` | GO namespace (`BP`, `MF`, `CC`) |
| `GeneID` | Comma-separated list of DEGs driving the category |

### 2. `*_GO_bubble_plot.pdf`

A bubble plot of the **top 20 Biological Process (BP)** terms, sorted by p-value (ascending):

- **X-axis** — percentage of DEGs in the category
- **Y-axis** — GO category label
- **Bubble size** — count of DEGs in the category
- **Bubble colour** — p-value (grey → dark red = low → high p-value)

---

## Supported Organisms

The `Go_annotations/` and `Length/` folders include reference files for:

| Organism | Annotation file | Length file |
|---|---|---|
| *E. coli* (gene ID) | `E_coliGeneID_GO.xlsx` | `GeneID_length.xlsx` |
| *E. coli* (locus tag) | `E_coliGeneID_GO_locusbased.xlsx` | `GeneID_Length_e_coli_locus.xlsx` |
| *Arabidopsis thaliana* | `ArabidopsisGeneID_GO.xlsx` | `GeneID_Length_Arabidopsis.xlsx` |
| Tomato (ITAG 3.0) | `ITAG3.0_GeneID_GO.xlsx` | `ITAG3.0_GeneID_Length.xlsx` |

To add a new organism, provide an annotation file with columns `GeneID` and `GOTerm`, and a length file with columns `GeneID` and `Length`, then update the paths in the script.

---

## Notes

- The bubble plot filters for **Biological Process (BP)** ontology only. To include Molecular Function or Cellular Component, adjust the `filter(ontology == "BP")` line in `generate_bubble_plot()`.
- The `setwd()` path in the script is an absolute Windows path — update it to your local environment before running.
- Output files are named after the input file (sans extension), so running the same file twice will **overwrite** previous results.
