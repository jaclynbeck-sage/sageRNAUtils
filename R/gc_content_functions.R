#' Get Gene GC Content from BioMart Data
#'
#' Queries BioMart for exon coordinates and sequences for each gene and
#' calculates GC content from the results.
#'
#' @details
#' This function is a wrapper around [calculate_gc_content()], which calculates
#' GC content based on either the entire gene, including both introns and exons,
#' or just exon sequences only. Coordinates and sequences are obtained from
#' BioMart, which may or may not have data on all genes present in the counts
#' matrix.
#'
#' Previous versions of Ensembl can be used instead of the latest version by
#' adding `version = <Ensembl version>` as an argument. GC content will be the
#' most accurate if the Ensembl version corresponding to the alignment reference
#' is used.
#'
#' This function obtains the following information from BioMart for each gene:
#' `ensembl_gene_id`, `ensembl_exon_id`, `exon_chrom_start`, `exon_chrom_end`,
#' `strand`, `start_position`, `end_position`, `external_gene_name`,
#' `chromosome_name`, `gene_biotype`.
#'
#' To add additional columns from Biomart, add them to the
#' `additional_attributes` argument. Attribute names must match names in
#' [biomaRt::listAttributes()].
#'
#' @param gene_ids a character vector of Ensembl gene IDs to get GC content for.
#' @param include_introns (optional) if `FALSE`, GC content will be calculated
#'   using exon sequences only. If `TRUE`, content will be calculated using the
#'   entire gene sequence, including introns. Defaults to `FALSE`.
#' @param mart (optional) a [`biomaRt::Mart`][biomaRt::Mart-class] object that
#'   is connected to a specific BioMart database and dataset. If `NULL`, this
#'   function will connect to BioMart using the `dataset` and `version`
#'   arguments instead.
#' @param dataset (optional) the name of the dataset to use from BioMart, if
#'   `mart` is `NULL`. Defaults to "hsapiens_gene_ensembl" for human genes. Both
#'   "hsapiens_gene_ensembl" and "mmusculus_gene_ensembl" have been tested with
#'   this function, but other species may work as well.
#' @param ensembl_version (optional) the Ensembl version to use (e.g. 110 or
#'   113), if `mart` is `NULL` and data should be fetched from an older archive.
#'   If `NULL`, the latest version of Ensembl will be queried.
#' @param additional_attributes (optional) a character vector of attributes /
#'   columns to fetch from BioMart in addition to the default columns returned
#'   by this function (see details). If `NULL`, only the default columns will be
#'   used.
#' @param verbose (optional) whether to print status messages
#' @param ... (optional) additional arguments to [biomaRt::useEnsembl()], for
#'   example an Ensembl mirror could be specified if BioMart is having
#'   connection issues.
#'
#' @return a data.frame where rows are genes and columns are variables. It will
#'   contain, at minimum, columns for `ensembl_gene_id`, `external_gene_name`,
#'   `gene_biotype`, `chromosome_name`, `gene_length`, and `percent_gc_content`.
#'   It will also contain columns listed in the `additional_attributes`
#'   argument.
#' @export
#'
#' @seealso [get_gc_content_gtf], [calculate_gc_content]
#'
#' @examples
#' \dontrun{
#' # Access version 110 instead of the latest Ensembl, and retrieve gene
#' # description during the query
#' gc_content <- get_gc_content_biomart(
#'   gene_ids = c("ENSG00000215388", "ENSG00000253771", "ENSG00000269584"),
#'   ensembl_version = 110,
#'   additional_attributes = c("description")
#' )
#'
#' # Optionally pass in a pre-made biomaRt object
#' mart <- biomaRt::useEnsembl(
#'   biomart = "ensembl",
#'   dataset = "hsapiens_gene_ensembl"
#' )
#' gc_content <- get_gc_content_biomart(
#'   gene_ids = c("ENSG00000215388", "ENSG00000253771", "ENSG00000269584"),
#'   mart = mart
#' )
#' }
get_gc_content_biomart <- function(gene_ids,
                                   include_introns = FALSE,
                                   mart = NULL,
                                   dataset = "hsapiens_gene_ensembl",
                                   ensembl_version = NULL,
                                   additional_attributes = NULL,
                                   verbose = TRUE,
                                   ...) {
  if (is.null(mart)) {
    if (verbose) {
      message("Connecting to BioMart...")
    }
    mart <- biomaRt::useEnsembl(biomart = "ensembl",
                                dataset = dataset,
                                version = ensembl_version,
                                ...)
  }

  if (include_introns) {
    attributes_list <- c(
      "ensembl_gene_id", "strand", "start_position", "end_position",
      "external_gene_name", "chromosome_name", "gene_biotype"
      # percentage_gene_gc_content?
    )
  } else {
    attributes_list <- c(
      "ensembl_gene_id", "ensembl_exon_id", "exon_chrom_start",
      "exon_chrom_end", "strand", "start_position", "end_position",
      "external_gene_name", "chromosome_name", "gene_biotype"
    )
  }

  attributes_list <- c(attributes_list, additional_attributes)

  if (verbose) {
    message(paste("Querying BioMart for", length(gene_ids), "genes..."))
  }

  gene_info <- biomaRt::getBM(
    attributes = attributes_list,
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  )

  # Fix strand to be + and - to work with GenomicRanges
  strand_map <- c("1" = "+", "-1" = "-")
  gene_info$strand <- strand_map[as.character(gene_info$strand)]

  gene_info <- gene_info |>
    # "start_position" is gene start. If the gene is on the reverse strand,
    # coordinates need to be reversed since the sequence returned by Biomart
    # is the reverse complement instead of the forward strand.
    dplyr::mutate(
      start = dplyr::case_match(strand,
                                "+" ~ exon_chrom_start - start_position + 1,
                                "-" ~ end_position - exon_chrom_end + 1),
      end = dplyr::case_match(strand,
                              "+" ~ exon_chrom_end - start_position + 1,
                              "-" ~ end_position - exon_chrom_start + 1)
    )

  if (verbose) {
    message(paste("BioMart returned",
                  length(unique(gene_info$ensembl_gene_id)),
                  "genes. Fetching sequences..."))
  }
  seq_results <- biomaRt::getSequence(id = unique(gene_info$ensembl_gene_id),
                                      type = "ensembl_gene_id",
                                      seqType = "gene_exon_intron",
                                      mart = mart)

  sequences <- Biostrings::DNAStringSet(seq_results$gene_exon_intron)
  names(sequences) <- seq_results$ensembl_gene_id

  grange <- GenomicRanges::makeGRangesFromDataFrame(
    gene_info,
    keep.extra.columns = TRUE,
    seqnames.field = "ensembl_gene_id"
  )
  grange$ensembl_gene_id <- as.character(GenomicRanges::seqnames(grange))

  if (verbose) {
    message("Calculating length and GC content for each gene...")
  }
  gene_data <- calculate_gc_content(grange, sequences,
                                    additional_attributes = additional_attributes)

  return(gene_data)
}


#' Get Gene GC Content from GTF and FASTA reference files
#'
#' Uses GTF and FASTA reference files to calculate gene length and GC content.
#'
#' @details
#' This function is a wrapper around [calculate_gc_content()], which calculates
#' GC content based on either the entire gene, including both introns and exons,
#' or just exon sequences only. Coordinates and sequences are obtained from
#' reference GTF and FASTA files, which ideally should match the files that were
#' used to align the RNA seq data.
#'
#' This function has been tested with human genome data from Gencode, although
#' it should work for other species and sources as well.
#'
#' This method is likely to be more accurate than [get_gc_content_biomart] in
#' cases where no Ensembl version exactly corresponds to the reference, or where
#' the genome is custom-built.
#'
#' @param gtf_file the full file path of the reference GTF file. URLs are also
#'   acceptable. The GTF file can be compressed (.gtf.gz) or not (.gtf). The GTF
#'   file must have columns named "gene_id", "gene_name", and "gene_type".
#' @param fasta_file the full file path of the reference FASTA file. URLs are
#'   also acceptable. The FASTA file can be compressed (.fa.gz) or not (.fa).
#' @param additional_attributes (optional) a vector of additional column names
#'   from the GTF file to retain in the result beyond what is returned by default
#'   in this function. If `NULL`, only the default columns are retained in the
#'   result.
#' @inheritParams get_gc_content_biomart
#'
#' @return a data.frame where rows are genes and columns are variables. It will
#'   contain columns for `ensembl_gene_id`, `external_gene_name`,
#'   `gene_biotype`, `chromosome_name`, `gene_length`, and `percent_gc_content`.
#' @export
#'
#' @seealso [get_gc_content_biomart], [calculate_gc_content]
#'
#' @examples
#' \dontrun{
#' # Using local files
#' gc_content <- get_gc_content_gtf(
#'   gtf_file = "data/gencode.v43.primary_assembly.annotation.gtf.gz",
#'   fasta_file = "data/GRCH38.primary_assembly.genome.fa.gz"
#' )
#'
#' # Using URLS
#' gc_content <- get_gc_content_gtf(
#'   gtf_file = paste0(
#'     "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/",
#'     "gencode.v43.primary_assembly.basic.annotation.gtf.gz"
#'   ),
#'   fasta_file = paste0(
#'     "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/",
#'     "GRCh38.primary_assembly.genome.fa.gz"
#'   )
#' )
#' }
get_gc_content_gtf <- function(gtf_file, fasta_file,
                               include_introns = FALSE,
                               additional_attributes = NULL,
                               verbose = TRUE) {
  if (verbose) {
    message("Importing GTF data...")
  }

  feature <- ifelse(include_introns, "gene", "exon")
  gtf <- rtracklayer::import(gtf_file,
                             format = "gtf",
                             feature.type = feature)

  # Some renames to match Biomart naming conventions
  gtf$ensembl_gene_id <- gtf$gene_id
  gtf$chromosome_name <- as.character(GenomicRanges::seqnames(gtf))
  gtf$hgnc_symbol <- gtf$gene_name
  gtf$gene_biotype <- gtf$gene_type

  if (verbose) {
    message("Importing FASTA data...")
  }

  fasta <- rtracklayer::import(fasta_file,
                               format = "fasta",
                               type = "DNA")
  # Change chromosome names from things like "chr2 2" to "chr2" to match
  # GTF file
  names(fasta) <- stringr::str_replace(names(fasta), " .*", "")

  if (verbose) {
    message("Calculating length and GC content for each gene...")
  }

  gene_data <- calculate_gc_content(gtf, fasta)

  return(gene_data)
}


#' Calculate Gene GC Content
#'
#' Calculate gene length and GC content of a set of genes from a GRanges object
#' and a DNAStringSet. It is highly recommended that you run either
#' [get_gc_content_gtf()] or [get_gc_content_biomart()] instead of calling this
#' function directly.
#'
#' @details
#'
#' **GTF/FASTA input**
#'
#' If the original input was imported from GTF and FASTA files, each `DNAString`
#' in `sequences` will be the sequence of an entire chromosome, and the
#' coordinates in `ranges_obj` should be relative to the start of the
#' chromosome, which is the default in GTF files.
#'
#' **BioMart input**
#'
#' If the original input was from Biomart, each `DNAString` in `sequences`
#' will be the sequence of a single gene, and the coordinates in `ranges_obj`
#' should be relative to the start of the gene. Note that for genes on the
#' reverse strand, Biomart returns the gene sequence as the reverse complement
#' of the reference strand, such that position 1 on the reference strand is
#' position `<gene_length>` on the reverse strand, and position `<gene_length>`
#' of the reference strand is position 1 of the reverse strand. However, the
#' coordinates given by Biomart are always relative to the reference strand even
#' if it returns a reverse complement sequence. This means that one of two
#' adjustments need to be made *prior* to calling this function.
#'
#' Either:
#'   1. `sequences` needs to be updated so that the sequence for every gene on
#'   the negative/reverse strand is reverse-complemented, or
#'   2. coordinates for exons on the reverse strand need to be adjusted to be
#'   relative to the start of reverse complement sequence, e.g. if BioMart
#'   returns coordinates `start = 3` and `end = 10`, the start position is now
#'   the *end* coordinate relative to the reverse complement strand, and should
#'   be adjusted as position <gene_length - 3 + 1>. The original end position is
#'   now the start coordinate and is adjusted as <gene_length - 10 + 1>.
#'
#' [get_gc_content_biomart] does this adjustment for you.
#'
#' No such adjustment is necessary for data from GTF and FASTA files, as the
#' sequences are always on the reference strand.
#'
#' **Other input**
#'
#' If input is being constructed by some other means, ensure that the
#' coordinates in `ranges_obj` are relative to the start of the corresponding
#' sequence(s) in `sequences`. `names(sequences)` must match the values in
#' `seqnames(ranges_obj)`.
#'
#' @param ranges_obj a [GenomicRanges::GRanges] object containing the start and
#'   end coordinates of each exon for each gene. This object must contain either
#'   a `gene_id` or `ensembl_gene_id` metadata column. It should also have an
#'   `external_gene_name` column. Values in `seqnames(ranges_obj)` must match
#'   the values in `names(sequences)`, and coordinates should be relative to the
#'   corresponding sequence in `sequences`. It is assumed that coordinates start
#'   at 1.
#' @param sequences a [Biostrings::DNAStringSet] object where each item is a
#'   `Biostrings::DNAString` containing the sequence for a gene or chromosome.
#'   `names(sequences)` must match the values in `seqnames(ranges_obj)`.
#' @param additional_attributes (optional) a vector of additional column names
#'   in `ranges_obj` to retain in the result beyond what is returned by default
#'   in this function. If `NULL`, only the default columns are retained in the
#'   result.
#'
#' @return a data.frame where rows are genes and columns are variables. It will
#'   contain, at a minimum, columns for `ensembl_gene_id`, `external_gene_name`,
#'   `gene_biotype`, `chromosome_name`, `gene_length`, and `percent_gc_content`.
#'   It will also contain columns listed in the `additional_attributes`
#'   argument if it is not `NULL`.
#' @export
#'
#' @seealso [get_gc_content_biomart], [get_gc_content_gtf]
calculate_gc_content <- function(ranges_obj,
                                 sequences,
                                 additional_attributes = NULL) {
  # "genes" will be a GRangesList where each list item is a GRanges object for a
  # specific gene. The GRanges object will have one or more ranges in it
  # corresponding to the de-duplicated ranges for all exonic regions in that
  # gene.
  id_field <- ifelse(
    "ensembl_gene_id" %in% colnames(GenomicRanges::mcols(ranges_obj)),
    "ensembl_gene_id", "gene_id"
  )

  genes <- GenomicRanges::split(
    ranges_obj,
    f = GenomicRanges::mcols(ranges_obj)[[id_field]]
  ) |>
    GenomicRanges::reduce()

  gene_length <- sum(GenomicRanges::width(genes))

  # Creates a DNAStringSetList where each list item is a set of sequences
  # (DNAStringSet) for a specific gene, corresponding to the exonic regions in
  # that gene
  exon_seqs <- BSgenome::getSeq(sequences, genes)

  gc <- sapply(exon_seqs, function(gene) {
    # Count all bases across all exonic sequences for this single gene
    gc_counts <- colSums(BSgenome::alphabetFrequency(gene))

    # Proportion of C and G bases. Some values in the sequence might be "N" and
    # we don't count "N" in the total so we can't use the calculated gene length
    # in the denominator
    sum(gc_counts[c("C", "G")]) / sum(gc_counts[c("A", "C", "G", "T")])
  })

  # Collapse to gene level
  fields_keep <- c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                   "chromosome_name", additional_attributes)

  fields_keep <- intersect(fields_keep, colnames(GenomicRanges::mcols(ranges_obj)))

  # Get data frame with one row per gene and add the calculated length and GC
  # content.
  gene_final <- ranges_obj |>
    as.data.frame() |>
    dplyr::select(dplyr::all_of(fields_keep)) |>
    dplyr::distinct() |>
    # Make sure gene length and gc are in the same order as the IDs in this data frame
    dplyr::mutate(
      gene_length = gene_length[.data[[id_field]]],
      percent_gc_content = gc[.data[[id_field]]],
      # Some gene symbols might have been assigned the Ensembl ID because there's
      # no symbol for that gene. This fixes those to be blank.
      external_gene_name = dplyr::case_when(
        external_gene_name == ensembl_gene_id ~ "",
        .default = external_gene_name
      )
    ) |>
    # Sort for neatness and easier diff in case of updates
    dplyr::arrange(ensembl_gene_id)

  return(gene_final)
}
