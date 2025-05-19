# PCA outlier functions --------------------------------------------------------

matrix_for_pca_no_outliers <- function() {
  set.seed(101)

  # 5 samples with 20 genes expressed at ~1, 20 at ~5, and 20 at ~10
  matrix(c(stats::rnorm(20*20, mean = 1, sd = 0.5),
           stats::rnorm(20*20, mean = 5, sd = 1),
           stats::rnorm(20*20, mean = 10, sd = 2)),
         ncol = 20, byrow = TRUE,
         dimnames = list(paste0("gene", 1:60), paste0("sample", 1:20)))
}

matrix_for_pca_with_outliers <- function() {
  set.seed(101)

  cbind(matrix_for_pca_no_outliers(),
        matrix(stats::rnorm(60*2, mean = 15, sd = 5), ncol = 2,
               dimnames = list(paste0("gene", 1:60), c("sample_out1", "sample_out2"))))
}

matrix_for_pca_with_groups <- function() {
  set.seed(101)

  data1 <- matrix_for_pca_with_outliers()
  data2 <- matrix_for_pca_no_outliers() * 2 # One group with no outliers
  data3 <- matrix_for_pca_with_outliers()

  # Rename the samples by group
  colnames(data1) <- paste0("group1_", colnames(data1))
  colnames(data2) <- paste0("group2_", colnames(data2))

  # Move the outliers in data 3 and rename the samples
  data3 <- data3[, c(21:22, 1:20)]
  colnames(data3) <- c("group3_sample_out3", "group3_sample_out4", paste0("group3_sample", 3:22))

  cbind(data1, data2, data3)
}


# Sex mismatch functions -------------------------------------------------------

matrix_for_sex_mismatches_with_symbols <- function() {
  set.seed(202)

  male <- matrix(c(rnorm(10, mean = 0, sd = 1), # low XIST
                   rnorm(4 * 10, mean = 10, sd = 2), # high Y genes
                   rnorm(10 * 10, mean = 5, sd = 3)), # some extra genes
                 ncol = 10,
                 byrow = TRUE,
                 dimnames = list(c("XIST", "RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", paste0("gene", 1:10)),
                                 paste0("sample_M", 1:10)))

  female <- matrix(c(rnorm(10, mean = 10, sd = 2), # high XIST
                     rnorm(4 * 10, mean = 0, sd = 1), # low Y genes
                     rnorm(10 * 10, mean = 5, sd = 3)), # some extra genes
                   ncol = 10,
                   byrow = TRUE,
                   dimnames = list(c("XIST", "RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", paste0("gene", 1:10)),
                                   paste0("sample_F", 1:10)))

  cbind(male, female)
}

matrix_for_sex_mismatches_with_ensembl_ids <- function() {
  set.seed(202)

  data <- matrix_for_sex_mismatches_with_symbols()
  rownames(data) <- c("ENSG00000229807", "ENSG00000129824", "ENSG00000198692",
                      "ENSG00000067048", "ENSG00000012817",
                      paste0("ENSG000000000", 1:10))
  data
}
