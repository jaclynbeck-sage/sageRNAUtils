#' Find Outliers by the IQR Method
#'
#' Labels points as outliers if they are below `Q1 - 1.5 * IQR` and/or above
#' `Q3 + 1.5 * IQR`, with the option to modify the multiplier.
#'
#' Note that `NA` values in `data` are ignored in the calculation of `Q1`, `Q3`,
#' and `IQR`, and all `NA` data points will be marked as outliers in the final
#' result.
#'
#' @param data a numeric vector of data points
#' @param tail (optional) either "both", "lower", or "upper". If "lower", only
#'   points below `Q1 - IQR_mult * IQR` will be marked as outliers. If "upper",
#'   only points above `Q3 + IQR_mult * IQR` will be marked as outliers. If
#'   "both", points in both tails will be marked as outliers. Defaults to
#'   "both".
#' @param IQR_mult (optional) a multiplier to apply to the IQR. Defaults to 1.5.
#'
#' @return a vector of TRUE/FALSE values, where TRUE indicates that the point is
#'   an outlier. Any `NA` values in `data` will be marked as TRUE.
#' @export
#'
#' @seealso [is_outlier_SD()]
#'
#' @examples
#' data <- rnorm(n = 1000)
#' outliers <- is_outlier_IQR(data = data, tail = "upper")
#' print(data[outliers])
is_outlier_IQR <- function(data, tail = "both", IQR_mult = 1.5) {
  iqr <- stats::IQR(data, na.rm = TRUE) * IQR_mult
  q1 <- stats::quantile(data, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(data, 0.75, na.rm = TRUE)

  switch(
    tail,
    "both" = (data < q1 - iqr) | (data > q3 + iqr),
    "upper" = data > q3 + iqr,
    "lower" = data < q1 - iqr
  ) | is.na(data)
}


#' Find Outliers using Mean and Standard Deviation
#'
#' Labels points as outliers if they are below `mean - n_sds * sd` and/or above
#' `mean + n_sds * sd`, where `n_sds` is a configurable multiplier.
#'
#' Note that `NA` values in `data` are ignored in the calculation of `mean` and
#' `sd`, and all `NA` data points will be marked as outliers in the final
#' result.
#'
#' @param data a numeric vector of data points
#' @param tail (optional) either "both", "lower", or "upper". If "lower", only
#'   points below `mean - n_sds * SD` will be marked as outliers. If "upper",
#'   only points above `mean + n_sds * SD` will be marked as outliers. If
#'   "both", points in both tails will be marked as outliers. Defaults to
#'   "both".
#' @param n_sds (optional) the number of standard deviations away from the mean
#'   that a point needs to be in order to be marked as an outlier. Defaults to
#'   4.
#'
#' @return a vector of TRUE/FALSE values, where TRUE indicates that the point is
#'   an outlier. Any `NA` values in `data` will be marked as TRUE.
#' @export
#'
#' @seealso [is_outlier_IQR()]
#'
#' @examples
#' data <- rnorm(n = 1000)
#' outliers <- is_outlier_SD(data = data, n_sds = 3)
#' print(data[outliers])
is_outlier_SD <- function(data, tail = "both", n_sds = 4) {
  mean_d <- mean(data, na.rm = TRUE)
  stdev <- n_sds * stats::sd(data, na.rm = TRUE)

  switch(
    tail,
    "both" = (data < mean_d - stdev) | (data > mean_d + stdev),
    "upper" = data > mean_d + stdev,
    "lower" = data < mean_d - stdev
  ) | is.na(data)
}


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
          dplyr::all_of(c("Total Sequences", "Sequences flagged as poor quality", "%GC")),
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


#' Load MultiQC JSON file
#'
#' Loads stats from a MultiQC JSON file and merges them all into one data frame.
#'
#' NOTE: This function has not been thoroughly tested.
#'
#' @param files a single path to a file or a vector or list of paths, which can
#'   be relative or absolute paths
#'
#' @return a data.frame with fields from the following metrics: picard, rsem,
#' samtools, samtools idx, cutadapt. If multiple JSON files were supplied, all
#' data will be concatenated into one data.frame.
#' @export
#'
#' @examples
#' \dontrun{
#' # One file
#' stats_df <- load_multiqc_json("data/multiqc_data.json")
#'
#' # Multiple files
#' stats_df <- load_multiqc_json(
#'   c("data/multiqc_data_1.json", "data/multiqc_data_2.json", "data/multiqc_data_3.json")
#' )
#' }
load_multiqc_json <- function(files) {
  data <- lapply(files, function(js_file) {
    js_txt <- readr::read_file(js_file)
    js_txt <- stringr::str_replace_all(js_txt, "NaN", stringr::str_escape('"NaN"'))
    js <- jsonlite::fromJSON(js_txt, simplifyVector = FALSE)

    # Remove general stats -- it's a mix of sample-level and fastqc read-level
    # data that doesn't merge together
    keep <- !grepl("multiqc_general_stats", names(js$report_saved_raw_data))
    js$report_saved_raw_data <- js$report_saved_raw_data[keep]

    lapply(js$report_saved_raw_data, function(item) {
      do.call(rbind, item)
    })
  })

  # Check for duplicate samples
  samps <- lapply(data, "[[", "multiqc_rsem") |>
    lapply(rownames) |>
    unlist()

  stopifnot(all(table(samps) == 1))

  # Merge the lists. Using sapply instead of lapply so it keeps the names
  data <- sapply(names(data[[1]]), function(item_name) {
    do.call(rbind, lapply(data, "[[", item_name))
  })

  # Helper function to convert items from "data", which are matrices, to
  # data frames with only the specified columns. Rownames are moved to a new
  # "specimenID" column and the specified prefix is added to the other column
  # names
  reformat_stats <- function(df, prefix, cols_keep) {
    df |>
      as.data.frame() |>
      dplyr::select({{cols_keep}}) |>
      # Add prefix to column names
      dplyr::rename_with(~ paste0(prefix, "_", .x), dplyr::everything()) |>
      # All columns are actually named lists, unlist them
      dplyr::mutate(dplyr::across(dplyr::everything(), ~unlist(.x))) |>
      # Make specimenID column
      tibble::rownames_to_column("specimenID") |>
      dplyr::mutate(specimenID = make.names(specimenID))
  }

  # TODO decide on:
  # multiqc_samtools_flagstat: ??
  # multiqc_fastqc: %GC, total_deduplicated_percentage <-- identical to fastqc data, but phred and base content not available

  picard_stats <- data$multiqc_picard_dups |>
    reformat_stats(prefix = "picard", cols_keep = PERCENT_DUPLICATION) |>
    dplyr::mutate(picard_PERCENT_DUPLICATION = picard_PERCENT_DUPLICATION * 100)

  rsem_stats <- data$multiqc_rsem |>
    reformat_stats(prefix = "rsem",
                   cols_keep = c(alignable_percent, Unique, Total)) |>
    dplyr::mutate(
      rsem_uniquely_aligned_percent = rsem_Unique / rsem_Total * 100
    ) |>
    dplyr::select(-rsem_Unique, -rsem_Total)

  samtools_stats <- data$multiqc_samtools_stats |>
    reformat_stats(
      prefix = "samtools",
      cols_keep = c(average_quality, insert_size_average, reads_mapped_percent,
                    reads_duplicated_percent, reads_MQ0_percent,
                    reads_QC_failed_percent)
    )

  # Each item in idxstats is a list with 2 numbers. The first number is the
  # number of reads mapped to the chromosome, the second number is the length
  # of the chromosome. We want the first number for X and Y
  xy_stats <- apply(data$multiqc_samtools_idxstats, 1, function(row) {
    c(chrX = row[["chrX"]][[1]], chrY = row[["chrY"]][[1]])
  }) |>
    t() |>
    reformat_stats(prefix = "samtools", cols_keep = c(chrX, chrY)) |>
    dplyr::mutate(total = samtools_chrX + samtools_chrY,
                  samtools_percent_mapped_X = samtools_chrX / total * 100,
                  samtools_percent_mapped_Y = samtools_chrY / total * 100) |>
    dplyr::select(-total, -samtools_chrX, -samtools_chrY)

  # This matrix may have separate rows for read 1 and read 2, which need to be
  # combined into a single row.
  cutadapt_stats <- data$multiqc_cutadapt |>
    reformat_stats(prefix = "cutadapt", cols_keep = percent_trimmed)

  if (all(grepl("_(1|2)", cutadapt_stats$specimenID))) {
    cutadapt_stats <- cutadapt_stats |>
      dplyr::mutate(read = ifelse(grepl("_1$", specimenID), "R1", "R2"),
                    specimenID = stringr::str_replace(specimenID, "_(1|2)$", "")) |>
      tidyr::pivot_wider(names_from = "read",
                         values_from = "cutadapt_percent_trimmed") |>
      dplyr::rename(cutadapt_percent_trimmed_R1 = R1,
                    cutadapt_percent_trimmed_R2 = R2) |>
      dplyr::mutate(
        cutadapt_mean_percent_trimmed = (cutadapt_percent_trimmed_R1 +
                                           cutadapt_percent_trimmed_R2) / 2
    )
  }

  technical_stats <- purrr::reduce(
    list(picard_stats, rsem_stats, samtools_stats, xy_stats, cutadapt_stats),
    dplyr::full_join
  )

  if (any(is.na(technical_stats))) {
    warning("NA values exist in technical stats data frame.")
  }

  return(technical_stats)

  # Possible QC metrics:
  # FastQC -
  #   median PHRED
  #   base content at each position
  # RSeqC -
  #   reads explained by first strand
  #   percent mapped reads
  #   number of reads mapped to genes
  #   percent reads mapped to junctions
  #   percent reads mapped to genes
  #   insert size inner distance
  #   percent known junctions
  #   80/20 ratio of gene body coverage
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
#' @param min_group_size (optional) the minimum number of samples that must be
#'   in a group in order to run outlier detection on that group. Defaults to 10.
#' @inheritParams find_pca_outliers
#'
#' @return a named list with the following items:
#'
#' \item{group_results}{a list containing the return values from [find_pca_outliers()] for each group}
#' \item{outliers}{a character vector of sample names that were marked as outliers in any group}
#'
#' Each item in `group_results` can be used to plot the results with [plot_pca_outliers()].
#' If a group was too small, it will not have an entry in `group_results`.
#' @export
#'
#' @seealso [find_pca_outliers()], [simple_log2norm()],
#'   [plot_pca_outliers()], [get_gc_content_biomart()], [get_gc_content_gtf()]
#'
#' @examples
#' \dontrun{
#' }
find_pca_outliers_by_group <- function(data, pca_group,
                                       n_sds = 4,
                                       metadata = NULL,
                                       sample_colname = "specimenID",
                                       gene_info = NULL,
                                       min_group_size = 10) {
  if (length(pca_group) == 1) {
    if (is.null(metadata)) {
      stop("`pca_group` is a single value but no metadata was supplied.")
    }
    if (!(pca_group %in% colnames(metadata))) {
      stop(paste0("\"", pca_group, "\" is not a valid column in `metadata`."))
    }
    pca_group <- metadata[, pca_group]
  } else if (length(pca_group) != ncol(data)) {
    stop("The number of samples in `pca_group` does not match the number of samples in `data`.")
  }

  if (!is.null(metadata)) {
    if (!(sample_colname %in% colnames(metadata))) {
      stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
    }
    if (!all(colnames(data) %in% metadata[, sample_colname])) {
      stop("`metadata` is missing samples that are present in `data`.")
    }
    # Make sure data and metadata are in the same sample order
    data <- data[, metadata[, sample_colname]]
  }

  results <- lapply(unique(pca_group), function(grp) {
    if (sum(pca_group == grp) < min_group_size) {
      message(paste0("Group '", grp, "' is too small to run outlier detection. Skipping..."))
      return(NULL)
    }

    data_group <- data[, pca_group == grp]

    if (is.null(metadata)) {
      meta_group <- NULL
    } else {
      meta_group <- metadata[pca_group == grp, ]
    }

    find_pca_outliers(data_group,
                      n_sds = n_sds,
                      metadata = meta_group,
                      sample_colname = sample_colname,
                      gene_info = gene_info)
  })

  names(results) <- unique(pca_group)

  # Remove any NULL entries from groups that were too small
  results <- results[lengths(results) > 0]

  return(list(
    group_results = results,
    outliers = lapply(results, "[[", "outliers") |> unlist() |> as.vector()
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
#' @param print_plot (optional) whether to print the plot out before returning,
#'   or not. Defaults to TRUE.
#' @param color (optional) the name of the column in `pca_df` to use as the
#'   "color" in the [ggplot2::aes()] specification. If `color` is `NULL` or is
#'   not a valid column name, all points will be black instead. Defaults to
#'   `is_outlier`.
#' @param ... (optional) additional [ggplot2::aes] specifications to pass to the
#'   plot, in the format <aes_spec> = "column_name_as_string". For example,
#'   `shape = "diagnosis"` or `size = "quantity"`.
#'
#' @return a [ggplot2::ggplot] object with the built plot.
#'
#' @examples
#' \dontrun{
#' }
plot_pca_outliers_ellipse <- function(pca_df, pc1_threshold, pc2_threshold,
                                      print_plot = TRUE,
                                      color = "is_outlier", ...) {
  # TODO add ggrepel labels and make outlier points larger. Find a way to get
  # title information
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

  if (print_plot) {
    print(plt)
  }

  return(plt)
}


#' Find Outliers by PCA
#'
#' Runs PCA on normalized data and detects outliers that are a certain number of
#' standard deviations away from the mean on any PC from 1:n_pcs.
#'
#' @details
#' This function runs `prcomp` on the provided data and marks samples as
#' outliers if they fall outside +/- n_sds * sd(PC) for any PC from 1:n_pcs. By
#' default, n_pcs = 3.
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
#'   on the log or log2 scale, for example as returned by [simple_log2norm()].
#' @param n_sds (optional) samples will be labeled as outliers if they are
#'   outside `n_sds` standard deviations from any PC. Defaults to 4.
#' @param n_pcs (optional) how many principal components to use for outlier
#'   detection. Typical values are 2 or 3. Defaults to 3.
#' @param metadata (optional) a data.frame where rows are samples and columns
#'   are covariates. If supplied, it must contain, at minimum, a column with
#'   sample IDs, and any additional columns will be merged into the result to
#'   make plotting simpler. If omitted or `NULL`, the result will contain only
#'   the columns returned in the `$x` matrix from `prcomp`. See details for more
#'   information.
#' @param sample_colname (optional) if `metadata` is provided, a character
#'   string naming the column in `metadata` that contains sample labels. Ignored
#'   if `metadata` is `NULL`. Defaults to "specimenID".
#' @param cutoff_method (optional) either `default` or `ellipse`. If `default`,
#'   samples are labeled as outliers if they are > `n_sds * sd(PC)` or < `-n_sds
#'   * sd(PC)` for any PC up to `n_pcs`. If `ellipse`, samples are labeled as
#'   outliers if they are outside an ellipse (or ellipsoid, if `n_pcs` > 2)
#'   defined by `radius1 = n_sds * sd(PC1)`, `radius2 = n_sds * sd(PC2)`, ...,
#'   `radiusN = n_sds * sd(PC_N)`, up to `n_pcs`.
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
#' \item{thresholds}{a named numeric vector of length `n_pcs` with the
#' thresholds used as cutoffs}
#' \item{outliers}{a character vector of sample names that were marked as outliers}
#'
#' The first two items in the list can be used to plot the results with
#' [plot_pca_outliers()].
#' @export
#'
#' @seealso [find_pca_outliers_by_group()], [simple_log2norm()],
#'   [plot_pca_outliers()], [get_gc_content_biomart()], [get_gc_content_gtf()]
#'
#' @examples
#' \dontrun{
#' }
find_pca_outliers <- function(data,
                              n_sds = 4,
                              n_pcs = 3,
                              metadata = NULL,
                              sample_colname = "specimenID",
                              cutoff_method = "default",
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

  data <- data[intersect(genes_use, rownames(data)), ]

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
    if (!(sample_colname %in% colnames(metadata))) {
      stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
    }

    if (!all(colnames(data) %in% metadata[, sample_colname])) {
      stop("`metadata` is missing samples that are present in `data`.")
    }

    pca_df <- merge(metadata, pca_res$x,
                    by.x = sample_colname, by.y = "row.names")
  } else {
    sample_colname <- "sample"
    pca_df <- as.data.frame(pca_res$x)
    pca_df$sample <- rownames(pca_df)
  }

  pc_names <- paste0("PC", 1:n_pcs)

  thresholds <- sapply(pc_names, function(pc) {
    stats::sd(pca_df[, pc]) * n_sds
  })

  if (cutoff_method == "default") {
    # Check that each PC value is within its threshold, for each PC separately
    in_bounds <- lapply(pc_names, function(pc) {
      abs(pca_df[, pc]) < thresholds[pc]
    }) |>
      as.data.frame() |>
      rowSums() == n_pcs # All values should be true across the row

  } else if (cutoff_method == "ellipse") {
    # Check that each point is inside an ellipse (or ellipsoid, if n_pcs > 2),
    # where each radius is one of the thresholds.
    # The formula for an ellipse/ellipsoid is (x^2 / a^2) + (y^2 / b^2) + ... = 1,
    # where (x, y, ...) are PCs and (a, b, ...) are radii defined by n_sds * sd(PC).
    in_bounds <- lapply(pc_names, function(pc) {
      pca_df[, pc]^2 / thresholds[pc]^2
    }) |>
      as.data.frame() |>
      rowSums() <= 1
  }

  pca_df$is_outlier <- !in_bounds
  outliers <- pca_df[, sample_colname][pca_df$is_outlier]

  return(list(
    pca_df = pca_df,
    thresholds = thresholds,
    outliers = outliers
  ))
}


#' Plot PCA Outlier Detection Results
#'
#' Plots PC1 vs PC2 and marks outlier samples as found by [find_pca_outliers()].
#' TODO: If more than two PCs were used, this function will plot a grid of PCx vs PCy,
#' up to 5 PCs.
#'
#' @param pca_df a data.frame where rows are samples and columns are PCs as
#'   output by [prcomp()], plus any additional covariates added by
#'   [find_pca_outliers()].
#' @param thresholds a named vector of outlier thresholds for each PC used in
#'   outlier detection, calculated in [find_pca_outliers()] as `n_sds * sd(PC)`.
#' @param print_plot (optional) whether to print the plot out before returning,
#'   or not. Defaults to TRUE.
#' @param color (optional) the name of the column in `pca_df` to use as the
#'   "color" in the [ggplot2::aes()] specification. If `color` is `NULL` or is
#'   not a valid column name, all points will be black instead. Defaults to
#'   `is_outlier`.
#' @param ... (optional) additional [ggplot2::aes] specifications to pass to the
#'   plot, in the format <aes_spec> = "column_name_as_string". For example,
#'   `shape = "diagnosis"` or `size = "quantity"`.
#'
#' @return a [ggplot2::ggplot] object with the built plot.
#' @export
#'
#' @examples
#' \dontrun{
#' }
plot_pca_outliers <- function(pca_df, thresholds,
                              print_plot = TRUE,
                              color = "is_outlier", ...) {
  # TODO add ggrepel labels and make outlier points larger. Find a way to get
  # title information
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
    ggplot2::geom_rect(ggplot2::aes(xmin = -PC1, xmax = PC1,
                                    ymin = -PC2, ymax = PC2),
                       data = as.data.frame(t(thresholds)),
                       linetype = "dotdash",
                       color = "blue",
                       fill = NA,
                       inherit.aes = FALSE) +
    ggplot2::theme_bw()

  if (print_plot) {
    print(plt)
  }

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
#'   log or log2 scale, for example as returned by [simple_log2norm()].
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
#' @seealso [plot_sex_mismatch_results()], [simple_log2norm()]
#'
#' @examples
#' \dontrun{
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

  if (!(sample_colname %in% colnames(metadata))) {
    stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
  }

  if (!(sex_colname %in% colnames(metadata))) {
    stop(paste0("\"", sex_colname, "\" is not a valid column in `metadata`."))
  }

  if (!all(colnames(data) %in% metadata[, sample_colname])) {
    stop("`metadata` is missing samples that are present in `data`.")
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
#' @param print_plot (optional) whether to print the two plots out before
#'   returning, or not. Defaults to TRUE.
#'
#' @return
#' a list of two [ggplot2::ggplot] plots, where the first plot is colored by
#' reported sex and the second plot is colored by mismatch status.
#' @export
#'
#' @seealso [find_sex_mismatches()]
#'
#' @examples
#' \dontrun{
#' }
plot_sex_mismatch_results <- function(sex_check_df, y_expr_threshold = 2.0,
                                      print_plot = TRUE) {
  # Put FALSE last so they are plotted on top of other dots
  sex_check_df <- dplyr::arrange(sex_check_df, dplyr::desc(sex_valid))

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

  if (print_plot) {
    print(plt1)
    print(plt2)
  }
  return(list(plt1, plt2))
}
