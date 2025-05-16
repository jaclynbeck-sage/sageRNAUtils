#' Format FastQC Output
#'
#' Loads a folder of fastqc.zip files output by FastQC, combines data for all
#' samples, and cleans up the final data. This is a wrapper around
#' [fastqcr::qc_read_collection()].
#'
#' @details
#' This function uses [fastqcr::qc_read_collection()] to parse all of the zip
#' files and combine them into data.frames. Some of the data is then cleaned up
#' to fix column typing, naming, and/or missing Percentage columns. The final
#' output will be a list of data.frames corresponding to the modules requested
#' from [fastqcr::qc_read_collection()] (by default, all modules).
#'
#' **Note:** This function can take anywhere from several seconds to over an hour, depending
#' on how many files there are.
#'
#' @param files a character vector of file names (including absolute or relative
#'   file paths) to parse. Files should all end with `_fastqc.zip`, as output by
#'   FastQC.
#' @param sample_names a character vector the same length as `files` with the
#'   values to put in the `sample` column of each combined data frame. Sample
#'   names should be unique to account for the fact that there are potentially
#'   two read files per specimen.
#' @param ... additional arguments to pass to [fastqcr::qc_read_collection()],
#'   for example different values for `modules` or `verbose`.
#'
#' @return a named list of data.frames containing stats for all samples listed
#'   in `files`. There will be one data frame per FastQC "module". See the
#'   `modules` argument of [fastqcr::qc_read_collection()] for the list of data
#'   available.
#'
#' @export
#'
#' @seealso [fastqcr::qc_read_collection()], [fastqcr::qc_read()]
#'
#' @examples
#' \dontrun{
#' # Fastqc output for 2 samples (example1 and example2) that had both R1 and R2
#' # fastq files.
#' files <- c(
#'   "data/example1_1_fastqc.zip", "data/example1_2_fastqc.zip",
#'   "data/example2_1_fastqc.zip", "data/example2_2_fastqc.zip"
#' )
#' sample_names <- c("example1_R1", "example1_R2", "example2_R1", "example2_R2")
#'
#' fastqc_data <- load_fastqc_output(files, sample_names)
#'
#' # Only get certain modules
#' fastqc_data_short <- load_fastqc_output(
#'   files, sample_names,
#'   modules = c("Summary", "Basic Statistics", "Overrepresented sequences")
#' )
#' }
load_fastqc_output <- function(files, sample_names, ...) {
  message("Reading fastqc.zip files...")

  # Stop readr from printing output for every file it reads in
  old_readr_opt <- getOption("readr.show_col_types")
  options(readr.show_col_types = FALSE)

  # Read in all of the statistics from the zip file for every sample
  fq_stats <- fastqcr::qc_read_collection(
    files = files,
    sample_names = sample_names,
    ...
  )

  # Restore previous option
  options(readr.show_col_types = old_readr_opt)

  message("Cleaning data frames...")

  # The following data frames in fq_stats do not need editing:
  #   per_base_sequence_quality, per_tile_sequence_quality,
  #   per_base_sequence_content, per_base_n_content,
  #   sequence_duplication_levels, overrepresented_sequences, adapter_content,
  #   kmer_content

  # summary - convert from long to wide format
  if ("summary" %in% names(fq_stats)) {
    fq_stats$summary <- tidyr::pivot_wider(
      data = fq_stats$summary,
      names_from = "module",
      values_from = "status"
    )
  }

  # basic_statistics -- convert from long to wide format and fix typing on
  # numerical columns. Also add which read each file came from.
  if ("basic_statistics" %in% names(fq_stats)) {
    fq_stats$basic_statistics <- fq_stats$basic_statistics |>
      tidyr::pivot_wider(names_from = "Measure", values_from = "Value") |>
      dplyr::mutate(
        # Unknown if this is universal to all FastQC output or specific to the
        # test data sets
        read = dplyr::case_when(
          grepl("_1.gz", Filename) ~ "R1",
          grepl("_2.gz", Filename) ~ "R2",
          .default = "Unknown"
        ),
        # These values are strings originally
        dplyr::across(
          c(`Total Sequences`, `Sequences flagged as poor quality`, `%GC`),
          as.numeric
        )
      ) |>
      dplyr::relocate(read, .after = Filename) |>
      as.data.frame()
  }

  # total_deduplicated_percentage -- The percentage column is erroneously named
  # "x[[module_names[i]]]"
  if ("total_deduplicated_percentage" %in% names(fq_stats)) {
    colnames(fq_stats$total_deduplicated_percentage)[2] <- "Percentage"
  }

  # Several data frames need a "Percentage" column added to normalize the
  # counts. They are all formatted the same way, so this function works on all
  # of them
  add_percentage <- function(df) {
    df |>
      dplyr::group_by(sample) |>
      dplyr::mutate(Percentage = Count / sum(Count, na.rm = TRUE) * 100)
  }

  # per_sequence_quality_scores
  if ("per_sequence_quality_scores" %in% names(fq_stats)) {
    fq_stats$per_sequence_quality_scores <- add_percentage(fq_stats$per_sequence_quality_scores)
  }

  # per_sequence_gc_content
  if ("per_sequence_gc_content" %in% names(fq_stats)) {
    fq_stats$per_sequence_gc_content <- add_percentage(fq_stats$per_sequence_gc_content)
  }

  # sequence_length_distribution
  if ("sequence_length_distribution" %in% names(fq_stats)) {
    fq_stats$sequence_length_distribution <- add_percentage(fq_stats$sequence_length_distribution)
  }

  return(fq_stats)
}


#' Find Outliers by PCA, Split by Group
#'
#' Finds outliers via PCA, but splits the data into groups and runs PCA outlier
#' detection on each group separately.
#'
#' This function splits the input expression matrix into groups of samples as
#' specified by `pca_group` and runs PCA outlier detection on each group
#' separately. This is recommended if the data has noticeable batch effects or
#' comes from separate tissues. See [find_pca_outliers] for more details about
#' input, filtering, and plotting.
#'
#' @param pca_group either a vector or a single character string describing the
#'   groups of samples. If `metadata` is provided, `pca_group` should be the
#'   name of the column in `metadata` that provides the group values. If
#'   `metadata` is missing or `NULL`, `pca_group` should be a vector the same
#'   length as the number of columns in `data`, where each entry is the group
#'   assignment of the corresponding sample in `data`.
#' @inheritParams find_pca_outliers
#'
#' @return a named list with the following items:
#'
#' \item{group_results}{a list containing the return values from [find_pca_outliers()] for each group}
#' \item{outliers}{a character vector of sample names that were marked as outliers in any group}
#'
#' Each item in `group_results` can be used to plot the results with [plot_pca_outliers()].
#' @export
#'
#' @seealso [find_pca_outliers()], [simple_lognorm()],
#'   [plot_pca_outliers()], [get_gc_content_biomart()], [get_gc_content_gtf()]
#'
#' @examples
#' \dontrun{
#' test
#' }
find_pca_outliers_by_group <- function(data, pca_group,
                                       n_sds = 4,
                                       metadata = NULL,
                                       sample_colname = "specimenID",
                                       gene_info = NULL) {
  metadata$pca_group <- metadata[, pca_group]
  pca_valid <- FALSE

  results <- lapply(unique(metadata$pca_group), function(grp) {
    meta_group <- subset(metadata, pca_group == grp)
    data_group <- data[, meta_group[, sample_colname]]

    find_pca_outliers(meta_group, data_group,
                      n_sds = n_sds,
                      sample_colname = sample_colname,
                      gene_info = gene_info)
  })

  names(results) <- unique(metadata$pca_group)

  return(list(
    group_results = results,
    outliers = unlist(lapply(results, "[[", "outliers"))
  ))
}


#' Find Outliers by PCA
#'
#' Runs PCA on normalized data and detects outliers that are a certain number of
#' standard deviations away from the mean on either PC1 or PC2.
#'
#' @details
#' This function runs `prcomp` on the provided data and marks samples as
#' outliers if they fall outside an ellipse defined by `radius1 = n_sds *
#' sd(PC1)` and `radius2 = n_sds * sd(PC2)`.
#'
#' **Gene filtering**
#'
#' The PCA can run using either protein-coding autosomal genes only, or all
#' genes in the expression matrix. Either set of genes is then filtered to genes
#' that are expressed in at least 80% of samples and have a variance > 0.001, to
#' remove genes that are all 0 or have extremely low variance.
#'
#' To run using protein-coding genes, a `gene_info` data.frame must be supplied
#' that contains the following columns: `ensembl_gene_id`, `gene_biotype`,
#' `chromosome_name`. The values in `ensembl_gene_id` should contain values in
#' the rownames of `data`, and all protein-coding genes should have a
#' `gene_biotype` value of "protein_coding". This data frame can be obtained
#' with [get_gc_content_biomart()] or [get_gc_content_gtf()], or built with
#' other tools to match the required structure.
#'
#' **Plotting**
#'
#' This function returns all the information required to easily plot the results
#' with [plot_pca_outliers()]. If `metadata` is provided, the result will also
#' include all columns that exist in `metadata`, enabling the use of those
#' columns with [ggplot2::aes] specifications like color and shape.
#'
#' @param data a matrix, data.frame, or matrix-like object of expression data
#'   where rows are genes and columns are samples. Data should be normalized and
#'   on the log or log2 scale, for example as returned by [simple_lognorm()].
#' @param n_sds (optional) samples will be labeled as outliers if they are
#'   outside the ellipse defined by `radius1 = n_sds * sd(PC1)` and `radius2 =
#'   n_sds * sd(PC2)`. Defaults to 4.
#' @param metadata (optional) a data.frame where rows are samples and columns
#'   are covariates. If supplied, it must contain, at minimum, a column with
#'   sample IDs, and any additional columns will be merged into the result to
#'   make plotting simpler. If omitted or `NULL`, the result will contain only
#'   the columns returned in the `$x` matrix from `prcomp`. See details for more
#'   information.
#' @param sample_colname (optional) if `metadata` is provided, a character
#'   string naming the column in `metadata` that contains sample labels. Ignored
#'   if `metadata` is `NULL`. Defaults to "specimenID".
#' @param gene_info (optional) a data.frame where rows are genes and columns are
#'   variables that provide metadata about each gene. If provided, the PCA will
#'   run using only protein-coding autosomal genes, which are extracted from
#'   `gene_info`. If `NULL` or not provided, the PCA will use all genes in
#'   `data`. Defaults to `NULL`. See details for more information.
#'
#' @return a named list with the following items:
#'
#' \item{pca_df}{the `$x` data.frame returned by `prcomp`, plus any columns from
#' `metadata` if it was supplied}
#' \item{pc1_threshold, pc2_threshold}{the calculated radii of the ellipse for
#'  PC1 and PC2}
#' \item{outliers}{a character vector of sample names that were marked as outliers}
#'
#' The first three items in the list can be used to plot the results with
#' [plot_pca_outliers()].
#' @export
#'
#' @seealso [find_pca_outliers_by_group()], [simple_lognorm()],
#'   [plot_pca_outliers()], [get_gc_content_biomart()], [get_gc_content_gtf()]
#'
#' @examples
#' \dontrun{
#' test
#' }
find_pca_outliers <- function(data,
                              n_sds = 4,
                              metadata = NULL,
                              sample_colname = "specimenID",
                              gene_info = NULL) {
  # If gene biotype information is provided, use protein-coding autosomal genes
  # only. Otherwise, use all genes in the data.
  if (!is.null(gene_info) & "gene_biotype" %in% colnames(gene_info) &
      "chromosome_name" %in% colnames(gene_info)) {
    genes_use <- subset(gene_info, gene_biotype == "protein_coding" &
                          !grepl("(X|Y|M|K|G)", chromosome_name)) |>
      dplyr::pull(ensembl_gene_id)
  } else {
    genes_use <- rownames(data)
  }

  data <- data[genes_use, ]

  # Remove genes that are mostly 0's, which may be 0 or a negative number in
  # log2-scale. For PCA, restrict to genes expressed in >= 80% of samples
  genes_keep <- rowSums(data > min(data)) >= 0.8 * ncol(data)
  data <- data[genes_keep, ]

  # Remove genes with very low variance
  variance <- matrixStats::rowVars(data)
  genes_keep <- names(variance)[variance > 0.001]
  data <- data[genes_keep, ]

  pca_res <- stats::prcomp(t(data), center = TRUE, scale. = TRUE)

  if (!is.null(metadata)) {
    pca_df <- merge(metadata, pca_res$x,
                    by.x = sample_colname, by.y = "row.names")
  } else {
    sample_colname <- "sample"
    pca_df <- pca_res$x
    pca_df$sample <- rownames(pca_df)
  }

  pc1_thresh <- stats::sd(pca_df$PC1) * n_sds
  pc2_thresh <- stats::sd(pca_df$PC2) * n_sds

  in_ellipse <- (pca_df$PC1^2 / pc1_thresh^2) + (pca_df$PC2^2 / pc2_thresh^2)

  pca_df$is_outlier <- in_ellipse > 1
  outliers <- pca_df[, sample_colname][pca_df$is_outlier]

  return(list(
    pca_df = pca_df,
    pc1_threshold = pc1_thresh,
    pc2_threshold = pc2_thresh,
    outliers = outliers
  ))
}


#' Plot PCA Outlier Detection Results
#'
#' Plots PC1 vs PC2 and marks outlier samples as found by [find_pca_outliers()].
#' It also draws the ellipse that was used to determine outliers.
#'
#' @param pca_df a data.frame where rows are samples and columns are PCs as
#'   output by [prcomp()], plus any additional covariates added by
#'   [find_pca_outliers()].
#' @param pc1_threshold the outlier threshold for PC1, which is calculated by
#'   [find_pca_outliers()] as `n_sds * sd(PC1)`.
#' @param pc2_threshold the outlier threshold for PC2, which is calculated by
#'   [find_pca_outliers()] as `n_sds * sd(PC2)`.
#' @param color (optional) the name of the column in `pca_df` to use as the
#'   "color" in the [ggplot2::aes()] specification. If `color` is `NULL` or is
#'   not a valid column name, all points will be black instead. Defaults to
#'   `is_outlier`.
#' @param ... (optional) additional [ggplot2::aes] specifications to pass to the
#'   plot, in the format <aes_spec> = "column_name_as_string". For example,
#'   `shape = "diagnosis"` or `size = "quantity"`.
#'
#' @return a [ggplot2::ggplot] object with the built plot. The plot is also
#'   printed out.
#' @export
#'
#' @examples
#' \dontrun{
#' test
#' }
plot_pca_outliers <- function(pca_df, pc1_threshold, pc2_threshold,
                              color = "is_outlier", ...) {
  # Formula for an ellipse is (x^2 / a^2) + (y^2 / b^2) = 1, so
  # y = +/- sqrt((1 - x^2 / a^2) * b^2)
  ellipse_points <- function(axis_A, axis_B) {
    x <- seq(from = -axis_A, to = axis_A, length.out = 1000)
    y <- sqrt((1 - x^2 / axis_A^2) * axis_B^2)

    data.frame(x = c(x, rev(x)), y = c(y, rev(-y)))
  }

  ellipse_df <- ellipse_points(pc1_threshold, pc2_threshold)

  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (X %in% colnames(pca_df)) {
      return(rlang::sym(X))
    } else {
      message(paste0("\"", X, "\" is not a valid column in `pca_df`. ",
                     "It will not be used in the plot."))
    }
  })

  if (!is.null(color)) {
    if (color %in% colnames(pca_df)) {
      aes_opts$color = ifelse(is.character(color), rlang::sym(color), color)
    } else {
      message(paste0("\"", color, "\" is not a valid column in `pca_df`. ",
                     "No color aes will be used."))
    }
  }

  plt <- ggplot2::ggplot(pca_df,
                         ggplot2::aes(x = PC1, y = PC2, !!!aes_opts)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_path(data = ellipse_df,
                       ggplot2::aes(x = x, y = y),
                       linetype = "dotdash",
                       color = "blue",
                       inherit.aes = FALSE) +
    ggplot2::theme_bw()

  print(plt)

  return(plt)
}


#' Find Mismatches Between Reported Sex and Y Chromosome Expression
#'
#' Estimates the sex of each sample based on expression of several genes on the
#' Y chromosome, and compares estimated sex to the reported sex to find
#' potential mismatches.
#'
#' @details
#' This function uses the mean expression of four Y-chromosome marker genes to
#' estimate the sex of each sample. The four genes are RPS4Y1, EIF1AY, DDX3Y,
#' and KDM5D, which were determined to be good markers for sex in [Staedtler,
#' et. al (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23829492).
#'
#' **Note:** This function currently only works for human genes.
#'
#' @param metadata a data.frame containing at least two columns, one for sample
#'   IDs and one for the reported sex of each sample. The reported sex should be
#'   either "male" or "female", all lowercase.
#' @param data a matrix, data.frame, or matrix-like object of expression data
#'   where rows are genes and columns are samples. The rownames of data should
#'   either be gene symbols or Ensembl IDs. Data should be normalized and on the
#'   log or log2 scale, for example as returned by [simple_lognorm()].
#' @param sample_colname (optional) a character string naming the column in
#'   `metadata` that contains sample labels. Defaults to "specimenID".
#' @param sex_colname (optional) a character string naming the column in
#'   `metadata` that contains sex labels. Defaults to "sex".
#' @param y_expr_threshold (optional) a numeric threshold that splits samples
#'   based on mean expression of the Y-chromosome marker genes. Samples with
#'   mean expression above the threshold will be labeled as "male" while samples
#'   with mean expression below the threshold will be labeled as "female".
#'   Defaults to 2.0, though this value may need to be be calibrated to the
#'   dataset.
#'
#' @return a named list with the following items:
#'
#' \item{sex_check_df}{a data.frame containing the reported sex, estimated sex,
#'   and expression values of the sex marker genes for each sample. It also
#'   contains a `sex_valid` column denoting whether the reported sex matches the
#'   estimated sex}
#' \item{mismatches}{a character vector of samples whose estimated sex does not
#'   match the reported sex}
#'
#' The `sex_check_df` data frame can be used with [plot_sex_mismatch_results()]
#' to view the results.
#' @export
#'
#' @seealso [plot_sex_mismatch_results()], [simple_lognorm()]
#'
#' @examples
#' \dontrun{
#' test = 1
#' }
find_sex_mismatches <- function(metadata, data,
                                sample_colname = "specimenID",
                                sex_colname = "sex",
                                y_expr_threshold = 2.0) {
  y_genes <- c(ENSG00000129824 = "RPS4Y1",
               ENSG00000198692 = "EIF1AY",
               ENSG00000067048 = "DDX3Y",
               ENSG00000012817 = "KDM5D")
  sex_genes <- c(ENSG00000229807 = "XIST", y_genes)

  # Strip any version numbers off of the gene names if they are Ensembl IDs
  # before checking that all necessary genes exist
  rownames(data) <- stringr::str_replace(rownames(data), "\\.[0-9]+", "")

  if (!all(names(sex_genes) %in% rownames(data)) &&
      !all(sex_genes %in% rownames(data))) {
    stop(paste0("Data is missing the required genes. Ensure that these genes ",
                "are present: ", paste(sex_genes, collapse = ", "), " (",
                paste(names(sex_genes), collapse = ", "), ")."))
  }

  # If genes are Ensembl IDs, convert to gene symbols
  if (any(grepl("ENSG00", rownames(data)))) {
    data <- data[names(sex_genes), ]
    rownames(data) <- sex_genes
  }

  sex_check <- metadata |>
    dplyr::select(dplyr::all_of(c(sample_colname, sex_colname))) |>
    merge(
      t(data[sex_genes, ]),
      by.x = sample_colname, by.y = "row.names"
    ) |>
    dplyr::mutate(
      mean_Y = rowMeans(dplyr::across(dplyr::all_of(y_genes))),
      reported_sex = .data[[sex_colname]],
      estimated_sex = dplyr::case_when(mean_Y > y_expr_threshold ~ "male",
                                       .default = "female"),
      sex_valid = (reported_sex == estimated_sex)
    ) |>
    # Put FALSE last so they are plotted on top of other dots
    dplyr::arrange(dplyr::desc(sex_valid)) |>
    as.data.frame()

  mismatches <- sex_check[, sample_colname][!sex_check$sex_valid]

  return(list(
    sex_check_df = sex_check,
    y_expr_threshold = y_expr_threshold,
    mismatches = mismatches
  ))
}


#' Plot Sex Mismatch Results
#'
#' This function plots the expression of XIST vs the mean expression of four
#' Y-chromosome genes (RPS4Y1, EIF1AY, DDX3Y, KDM5D) for each sample and
#' highlights mismatches between reported and estimated sex.
#'
#' @param sex_check_df a data.frame containing `reported_sex`, `estimated_sex`,
#'   and expression values for `XIST` and `mean_Y`, as output by
#'   [find_sex_mismatches()]. It should also contain a `sex_valid` column
#'   denoting whether the reported sex matches the estimated sex.
#' @param y_expr_threshold (optional) a numeric threshold that splits samples
#'   based on mean expression of the Y-chromosome marker genes. The threshold
#'   will be drawn as a line across the plots. Defaults to 2.0.
#'
#' @return
#' a list of two [ggplot2::ggplot] plots, where the first plot is colored by
#' reported sex and the second plot is colored by mismatch status. Each plot
#' is also printed out.
#' @export
#'
#' @seealso [find_sex_mismatches()]
#'
#' @examples
#' \dontrun{
#' test
#' }
plot_sex_mismatch_results <- function(sex_check_df, y_expr_threshold = 2.0) {
  plt1 <- ggplot2::ggplot(sex_check_df,
                          ggplot2::aes(x = XIST, y = mean_Y, color = reported_sex)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  status_colors <- c("TRUE" = "black", "FALSE" = "red")

  plt2 <- ggplot2::ggplot(sex_check_df,
                          ggplot2::aes(x = XIST, y = mean_Y, color = sex_valid)) +
    ggplot2::geom_point(size = ifelse(sex_check_df$sex_valid, 0.5, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::geom_hline(yintercept = y_expr_threshold, linetype = "dotdash",
                        color = "red", alpha = 0.8)

  print(plt1)
  print(plt2)
  return(list(plt1, plt2))
}
