# GOseq enrichment with Benjamini-Hochberg multiple-testing correction.
#
# For every file in input/ this script analyses all genes marked as
# differentially expressed in the Expression column.
#
# BH correction is applied across all GO terms (BP, MF and CC together) within
# each input comparison. The script writes a complete workbook and
# one Biological Process plot containing up to 20 significant BP terms.

required_packages <- c("readxl", "ggplot2", "dplyr", "writexl", "goseq")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install CRAN packages normally and install goseq with ",
    "BiocManager::install('goseq')."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(writexl)
  library(goseq)
})

# User-adjustable significance setting.
PADJ_CUTOFF <- 0.05
MAX_BP_TERMS_TO_PLOT <- 20

get_script_directory <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)
  if (length(file_argument) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  script_path <- sub("^--file=", "", file_argument[[1]])
  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}

base_folder <- get_script_directory()
input_folder <- file.path(base_folder, "input")
output_folder <- file.path(base_folder, "output")
go_annotation_path <- file.path(
  base_folder, "Go_annotations", "E_coliGeneID_GO.xlsx"
)
gene_length_path <- file.path(base_folder, "Length", "GeneID_length.xlsx")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

normalise_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_column <- function(data, accepted_names, description) {
  normalised_columns <- normalise_name(names(data))
  match_index <- match(normalise_name(accepted_names), normalised_columns)
  match_index <- match_index[!is.na(match_index)]
  if (length(match_index) == 0) {
    stop(
      "Could not find ", description, ". Available columns: ",
      paste(names(data), collapse = ", ")
    )
  }
  names(data)[match_index[[1]]]
}

as_expression_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  if (is.numeric(x)) return(replace(x != 0, is.na(x), FALSE))
  normalised <- tolower(trimws(as.character(x)))
  flag <- normalised %in% c("true", "t", "1", "yes", "y")
  replace(flag, is.na(flag), FALSE)
}

go_annotations <- readxl::read_excel(go_annotation_path)
gene_lengths <- readxl::read_excel(gene_length_path)

annotation_gene_column <- find_column(
  go_annotations, c("GeneID", "gene_id"), "the annotation gene-ID column"
)
annotation_term_column <- find_column(
  go_annotations, c("GOTerm", "GO", "GOID"), "the GO-term column"
)
length_gene_column <- find_column(
  gene_lengths, c("GeneID", "gene_id"), "the length-table gene-ID column"
)
length_value_column <- find_column(
  gene_lengths, c("Length", "gene_length"), "the gene-length column"
)

go_annotations_list <- split(
  as.character(go_annotations[[annotation_term_column]]),
  as.character(go_annotations[[annotation_gene_column]])
)
gene_lengths_vector <- setNames(
  suppressWarnings(as.numeric(gene_lengths[[length_value_column]])),
  as.character(gene_lengths[[length_gene_column]])
)

add_de_gene_ids <- function(go_results, selected_genes) {
  cat2gene <- stack(go_annotations_list)[2:1]
  names(cat2gene) <- c("GeneID", "category")

  de_cat2gene <- cat2gene %>%
    filter(GeneID %in% selected_genes) %>%
    distinct(category, GeneID) %>%
    group_by(category) %>%
    summarise(DE_GeneID = paste(sort(unique(GeneID)), collapse = ", "),
              .groups = "drop")

  go_results %>%
    left_join(de_cat2gene, by = "category") %>%
    mutate(
      over_represented_p_adjusted = p.adjust(
        over_represented_pvalue, method = "BH"
      ),
      under_represented_p_adjusted = p.adjust(
        under_represented_pvalue, method = "BH"
      )
    ) %>%
    arrange(over_represented_p_adjusted, over_represented_pvalue)
}

run_goseq_analysis <- function(de_flags) {
  if (!is.logical(de_flags) || is.null(names(de_flags))) {
    stop("The DE vector is not a named logical vector.")
  }

  selected_count <- sum(de_flags, na.rm = TRUE)
  if (selected_count == 0) {
    stop("No genes were selected for GOseq analysis.")
  }

  bias_vector <- gene_lengths_vector[names(de_flags)]
  missing_lengths <- names(de_flags)[is.na(bias_vector)]
  if (length(missing_lengths) > 0) {
    stop(
      "Missing gene lengths for ", length(missing_lengths), " gene(s), e.g. ",
      paste(head(missing_lengths, 10), collapse = ", ")
    )
  }

  message("  Selected ", selected_count, " DE genes")
  pwf <- goseq::nullp(de_flags, bias.data = unname(bias_vector))
  raw_results <- goseq::goseq(pwf, gene2cat = go_annotations_list)
  add_de_gene_ids(raw_results, names(de_flags)[de_flags])
}

save_bubble_plot <- function(results, analysis_title, file_stem) {
  # Keep the plotting style of the original 20240314-lo.r script. The only
  # analytical changes are BP-only filtering, BH-adjusted p-values, and keeping
  # up to the 20 most significant terms after applying the adjusted-p cutoff.
  plot_data <- results %>%
    filter(
      ontology == "BP",
      is.finite(over_represented_p_adjusted),
      over_represented_p_adjusted <= PADJ_CUTOFF,
      numInCat > 0
    ) %>%
    mutate(percent_DEG = (numDEInCat / numInCat) * 100) %>%
    arrange(over_represented_p_adjusted) %>%
    slice_head(n = MAX_BP_TERMS_TO_PLOT) %>%
    mutate(
      go_term_label = paste(category, term, sep = " - ")
    )

  plot_data$go_term_label <- factor(
    plot_data$go_term_label,
    levels = rev(plot_data$go_term_label)
  )

  if (nrow(plot_data) == 0) {
    plot_object <- ggplot() +
      annotate(
        "text", x = 0, y = 0,
        label = paste0("No BP terms passed adjusted p <= ", PADJ_CUTOFF)
      ) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = analysis_title) +
      theme_void()
  } else {
    plot_object <- ggplot(
      plot_data,
      aes(
        x = percent_DEG,
        y = go_term_label,
        size = numDEInCat,
        color = over_represented_p_adjusted
      )
    ) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(
        low = "slategrey",
        high = "firebrick4",
        name = "P-adjusted"
      ) +
      labs(
        title = analysis_title,
        x = "% of DEG",
        y = "GO term",
        size = "Count of DEGs"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
  }

  plot_path <- file.path(output_folder, paste0(file_stem, ".pdf"))
  ggsave(
    plot_path,
    plot = plot_object,
    device = "pdf",
    width = 12,
    height = 6
  )
}

save_analysis_outputs <- function(results, comparison_stem) {
  output_stem <- paste0(comparison_stem, "_GOseq_BH_p_adjusted")
  significant_results <- results %>%
    filter(over_represented_p_adjusted <= PADJ_CUTOFF)

  workbook_path <- file.path(output_folder, paste0(output_stem, "_results.xlsx"))
  writexl::write_xlsx(
    list(
      All_GO_terms = results,
      Significant_overrep = significant_results
    ),
    workbook_path
  )

  save_bubble_plot(
    results,
    "Significant Biological Processes (up to 20)",
    paste0(output_stem, "_top20_significant_BP_bubble_plot")
  )
}

input_files <- list.files(
  input_folder,
  pattern = "\\.(xlsx|tab|tabular)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(input_files) == 0) {
  stop("No .xlsx, .tab or .tabular files found in ", input_folder)
}

for (file_path in input_files) {
  message("Processing: ", basename(file_path))

  gene_expression <- if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
    readxl::read_excel(file_path)
  } else {
    read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  }

  input_gene_column <- find_column(
    gene_expression, c("GeneID", "gene_id"), "the input gene-ID column"
  )
  expression_column <- find_column(
    gene_expression, c("Expression", "DE", "is_DE"), "the DE flag column"
  )
  analysis_data <- gene_expression %>%
    transmute(
      GeneID = trimws(as.character(.data[[input_gene_column]])),
      is_DE = as_expression_flag(.data[[expression_column]])
    ) %>%
    filter(!is.na(GeneID), GeneID != "") %>%
    distinct(GeneID, .keep_all = TRUE)

  de_genes <- setNames(analysis_data$is_DE, analysis_data$GeneID)
  comparison_stem <- tools::file_path_sans_ext(basename(file_path))
  results <- run_goseq_analysis(de_genes)
  save_analysis_outputs(results, comparison_stem)

  message("Finished: ", basename(file_path))
}

message("All BH-adjusted GOseq outputs saved in: ", output_folder)
