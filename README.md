# GOseq GO enrichment analysis for *Escherichia coli*

## Overview

This repository performs Gene Ontology (GO) enrichment analysis using the
Bioconductor [`goseq`](https://bioconductor.org/packages/release/bioc/html/goseq.html)
package. GOseq accounts for gene-length bias in RNA-seq differential-expression
analyses.

The pipeline processes every supported file in the `input/` directory and
produces:

- a complete Excel workbook containing GOseq results;
- Benjamini–Hochberg-adjusted P values;
- a separate worksheet containing significantly overrepresented GO terms;
- a Biological Process bubble plot containing up to the 20 most significant
  terms.

## Important: definition of the tested gene set

`Main script.r` does not calculate differential expression, determine
regulation direction, or require a `log2FC` column.

Each input file defines its tested gene set through the `Expression` column:

- `Expression = TRUE` indicates that the gene belongs to the gene set being
  tested;
- `Expression = FALSE` indicates that the gene belongs to the background gene
  universe but not to the tested set.

Combined, upregulated, and downregulated enrichment analyses can be performed
by providing separate input files in which the `Expression` column has already
been defined appropriately during the preceding differential-expression
analysis.

## Repository structure

```text
GoSeq/
├── Main script.r
├── input/
│   └── one or more .tab, .tabular, or .xlsx input files
├── Go_annotations/
│   └── E_coliGeneID_GO.xlsx
├── Length/
│   └── GeneID_length.xlsx
└── output/
```

Example input files may represent:

```text
expression_15min_combined.tab
expression_15min_upregulated.tab
expression_15min_downregulated.tab
expression_30min_combined.tab
expression_30min_upregulated.tab
expression_30min_downregulated.tab
```

The upregulated and downregulated input files are optional. They must be
prepared before running this pipeline.

File and directory names are case-sensitive on Linux. In particular, the
required gene-length file is:

```text
Length/GeneID_length.xlsx
```

## Requirements

The analysis requires R and the following packages:

- `goseq`
- `readxl`
- `ggplot2`
- `dplyr`
- `writexl`

Install the CRAN packages with:

```r
install.packages(c("readxl", "ggplot2", "dplyr", "writexl"))
```

Install GOseq through Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("goseq")
```

The script checks whether the required packages are available and stops with
an informative message if any package is missing.

## Input files

The script processes every `.tab`, `.tabular`, or `.xlsx` file placed in the
`input/` directory.

Each input file must contain the following columns:

| Column | Description |
|---|---|
| `GeneID` | Gene identifier matching the GO annotation and gene-length files |
| `Expression` | Logical indicator specifying whether the gene belongs to the tested gene set |

Example:

```text
GeneID    Expression
b0001     FALSE
b0002     TRUE
b0003     FALSE
```

The complete tested-gene universe should be retained in every input file.
Genes belonging to the tested gene set should have `Expression = TRUE`; all
other genes in the background universe should have `Expression = FALSE`.

Alternative column names accepted by the script include:

- `gene_id` for the gene identifier;
- `DE` or `is_DE` for the gene-set indicator.

### Combined analysis

For a combined DEG analysis:

- `Expression = TRUE` for all selected differentially expressed genes;
- `Expression = FALSE` for all other genes in the tested background.

### Upregulated analysis

For an upregulated analysis:

- `Expression = TRUE` only for selected upregulated DE genes;
- `Expression = FALSE` for all other genes.

### Downregulated analysis

For a downregulated analysis:

- `Expression = TRUE` only for selected downregulated DE genes;
- `Expression = FALSE` for all other genes.

The upregulated and downregulated classifications must be generated during the
preceding differential-expression analysis. `Main script.r` does not use
`log2FC` and does not assign regulation direction.

## Statistical analysis

GO enrichment is performed using the Bioconductor `goseq` package, with
gene-length bias estimated using `nullp()`.

For each input file, `goseq()` calculates over-representation and
under-representation P values for the annotated GO terms.

Benjamini–Hochberg multiple-testing correction is applied across all tested
GO terms—Biological Process, Molecular Function, and Cellular Component
combined—separately for each input file.

Raw over- and under-representation P values are adjusted separately.

GO terms are considered significantly overrepresented when:

```text
BH-adjusted P value ≤ 0.05
```

The adjusted values are reported as:

- `over_represented_p_adjusted`;
- `under_represented_p_adjusted`.

The significance threshold is defined in `Main script.r`:

```r
PADJ_CUTOFF <- 0.05
```

## Bubble plots

Bubble plots display Biological Process terms only.

The plotting procedure:

1. retains terms with `over_represented_p_adjusted <= 0.05`;
2. ranks them by adjusted P value;
3. displays up to the 20 most significant terms.

If fewer than 20 Biological Process terms pass the significance threshold,
all significant terms are displayed.

Plot elements represent:

- **X-axis:** percentage of genes in the GO category that belong to the tested
  gene set;
- **Y-axis:** GO identifier and term name;
- **Bubble size:** number of selected genes assigned to the GO term;
- **Bubble colour:** BH-adjusted over-representation P value (`P-adjusted`).

If no Biological Process term passes the adjusted P-value threshold, the plot
reports that no significant terms were identified.

## Configuration

By default, the pipeline uses:

```r
go_annotation_path <- file.path(
  base_folder,
  "Go_annotations",
  "E_coliGeneID_GO.xlsx"
)

gene_length_path <- file.path(
  base_folder,
  "Length",
  "GeneID_length.xlsx"
)
```

The script automatically creates the following output directory:

```text
output/
```

To analyse another organism, update the annotation and gene-length paths in
`Main script.r` and provide compatible files.

## Running the analysis

From the repository root, run:

```bash
Rscript "Main script.r"
```

Alternatively, in R or RStudio, set the repository root as the working
directory and run:

```r
source("Main script.r")
```

The script automatically processes every supported file in `input/`.

## Output files

All results are written to:

```text
output/
```

For every input file, the pipeline produces:

```text
<input-name>_GOseq_BH_p_adjusted_results.xlsx
<input-name>_GOseq_BH_p_adjusted_top20_significant_BP_bubble_plot.pdf
```

Each Excel workbook contains two worksheets:

| Worksheet | Contents |
|---|---|
| `All_GO_terms` | All GO terms tested by GOseq, including raw and BH-adjusted P values |
| `Significant_overrep` | Terms with BH-adjusted over-representation P values ≤ 0.05 |

The output columns include:

| Column | Description |
|---|---|
| `category` | GO term identifier |
| `over_represented_pvalue` | Raw over-representation P value |
| `under_represented_pvalue` | Raw under-representation P value |
| `numDEInCat` | Number of selected genes assigned to the GO category |
| `numInCat` | Total number of genes assigned to the GO category |
| `term` | GO term description |
| `ontology` | GO namespace: `BP`, `MF`, or `CC` |
| `DE_GeneID` | Selected genes assigned to the GO category |
| `over_represented_p_adjusted` | BH-adjusted over-representation P value |
| `under_represented_p_adjusted` | BH-adjusted under-representation P value |

The Excel workbooks contain results for all three GO ontologies. The PDF
bubble plots display significant Biological Process results only.

## Reproducibility notes

The analysis uses only files included in this repository. It does not require
files from a volcano-plot directory or any other external analysis folder.

For reproducible results:

- retain the complete tested-gene universe in every input file;
- define the tested gene set through the `Expression` column;
- prepare upregulated and downregulated gene sets before running this pipeline;
- use the supplied GO annotation and gene-length files;
- preserve case-sensitive directory and filenames;
- record the repository release or commit used for the analysis.

## Version 1.1.0

Version 1.1.0 adds:

- Benjamini–Hochberg multiple-testing correction;
- adjusted over- and under-representation P values;
- an adjusted P-value significance threshold of 0.05;
- processing of separately prepared combined, upregulated, or downregulated
  gene-set input files;
- Biological Process plots showing up to 20 significant terms;
- updated Excel workbooks;
- repository-relative file paths with no external analysis-folder dependency.

This release supersedes `x1.0.1` for the corrected GO enrichment analyses.
