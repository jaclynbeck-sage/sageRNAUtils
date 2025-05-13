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
        across(
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
