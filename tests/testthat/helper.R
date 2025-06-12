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


# regression functions ---------------------------------------------------------

df_with_unusable_columns <- function() {
  set.seed(303)

  data.frame(
    # Unique values in non-numeric (character) field
    id_field = paste("sample", 1:10),

    # Numeric, logical, character, and factor columns that should not get removed
    age = runif(10, min = 50, max = 100),
    treatment = c(rep(TRUE, 5), rep(FALSE, 5)),
    sex = c(rep("male", 5), rep("female", 5)),
    dose = factor(c(rep(50, 5), rep(100, 5))),

    # Column with some NA values which should not get removed
    weight = c(runif(8, min = 50, max = 100), NA, NA),

    # Numeric, logical, character, and factor columns with same values
    num_1 = 7,
    log_1 = TRUE,
    char_1 = "Yes",
    factor_1 = factor(rep("group1", 10)),

    # Column with all NA values, which should get removed
    na_1 = rep(NA, 10),

    # Factor column with unique values, which should get removed
    factor_2 = factor(paste("group", 1:10))
  )
}


df_with_correlated_vars <- function() {
  num_base <- c(1:5, 5:1)

  data.frame(
    # highly-correlated numerics (R^2 ~0.56)
    num_1 = num_base,
    num_2 = num_base + c(6:10, 5:1),

    # Below 0.5 R^2 with num_1 and num_2 (0 and 0.33)
    num_3 = 1:10,

    # highly-correlated categoricals. They are also highly correlated with num_3
    cat_1 = c(rep("b1", 3), rep("b2", 3), rep("b3", 4)),
    cat_2 = c(rep("c1", 3), rep("c2", 4), rep("c3", 2), "c4"),

    # low correlation with cat_1 and cat_2
    cat_3 = c(rep(c("d1", "d2"), 4), "d2", "d2")
  )
}

fake_na_vars <- function(df) {
  na <- as.data.frame(df)[1, ]
  na[1, ] <- 0
  na
}

r2_matrix <- function() {
  matrix(c(1, 0.8, 0.2, 0.4,
           0.8, 1, 0.1, 0.1,
           0.2, 0.1, 1, 0.6,
           0.4, 0.1, 0.6, 1),
         nrow = 4,
         dimnames = list(paste0("v", 1:4), paste0("v", 1:4)))
}

r2_matrix_multiway <- function() {
  matrix(c(1, 0.8, 0.2, 0.4, 0.1,
           0.8, 1, 0.1, 0.1, 0.1,
           0.2, 0.1, 1, 0.6, 0.7,
           0.4, 0.1, 0.6, 1, 0.1,
           0.1, 0.1, 0.7, 0.1, 1),
         nrow = 5,
         dimnames = list(paste0("v", 1:5), paste0("v", 1:5)))
}
