# remove_unusable_covariates ---------------------------------------------------

test_that("remove_unusable_covariates removes NA and same-value columns", {
  df <- df_with_unusable_columns()
  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_s3_class(new_df, "data.frame")
  expect_identical(colnames(new_df), c("age", "treatment", "sex", "dose", "weight"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_unusable_covariates removes columns with same value after NA removal", {
  df <- df_with_unusable_columns()
  df$sex <- c(rep("male", 8), NA, NA)

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_identical(colnames(new_df), c("age", "treatment", "dose", "weight"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_unusable_covariates handles factors with unused levels", {
  df <- df_with_unusable_columns()

  # Make a factor with 2 levels, then subset to only one value
  df$factor_1 <- factor(c(rep("group1", 7), rep("group2", 3)))
  df <- subset(df, factor_1 == "group1")

  expect_identical(levels(df$factor_1), c("group1", "group2"))

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_identical(colnames(new_df), c("age", "treatment", "sex", "dose", "weight"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_unusable_covariates removes logical column with unique values", {
  df <- df_with_unusable_columns()[1:2, ]
  df$log_1 <- c(TRUE, FALSE)

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_identical(colnames(new_df), c("age", "weight"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_unusable_covariates returns empty data.frame when all columns are removed", {
  df <- df_with_unusable_columns()[, c("num_1", "log_1", "char_1", "factor_1", "na_1")]

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_s3_class(new_df, "data.frame")
  expect_equal(ncol(new_df), 0)
})

test_that("remove_unusable_covariates returns one-column data.frame when only one column remains", {
  df <- df_with_unusable_columns()[, c("age", "num_1", "log_1")]

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_s3_class(new_df, "data.frame")
  expect_equal(ncol(new_df), 1)
  expect_identical(new_df$age, df$age)
})

test_that("remove_unusable_covariates keeps always_keep columns", {
  df <- df_with_unusable_columns()

  # Single column
  new_df1 <- remove_unusable_covariates(df, always_keep = "id_field", verbose = FALSE)

  expect_identical(colnames(new_df1), c("id_field", "age", "treatment", "sex", "dose", "weight"))
  expect_identical(new_df1, df[, colnames(new_df1)])

  # Multiple columns, and listed out of order
  new_df2 <- remove_unusable_covariates(df,
                                        always_keep = c("char_1", "id_field", "num_1"),
                                        verbose = FALSE)

  expect_identical(colnames(new_df2),
                   c("id_field", "age", "treatment", "sex", "dose", "weight", "num_1", "char_1"))
  expect_identical(new_df2, df[, colnames(new_df2)])
})

test_that("remove_unusable_covariates retains row names", {
  df <- df_with_unusable_columns()
  rownames(df) <- df$id_field

  new_df <- remove_unusable_covariates(df, verbose = FALSE)

  expect_identical(rownames(new_df), rownames(df))
})

test_that("remove_unusable_covariates works on matrices", {
  df <- df_with_unusable_columns()
  mat <- df[, c("age", "weight", "num_1", "na_1")] |> as.matrix()
  rownames(mat) <- df$id_field

  new_df <- remove_unusable_covariates(mat, verbose = FALSE)

  expect_s3_class(new_df, "data.frame")
  expect_identical(colnames(new_df), c("age", "weight"))
  expect_identical(new_df, as.data.frame(mat[, colnames(new_df)]))
  expect_identical(rownames(new_df), rownames(mat))
})

test_that("remove_unusable_covariates prints removed variables when verbose is TRUE", {
  df <- df_with_unusable_columns()
  df <- df[, c("age", "weight", "num_1", "na_1")]

  expect_output(remove_unusable_covariates(df, verbose = TRUE),
                regexp = "Removing 2 columns: num_1, na_1")
})


# .make_symmetrical ------------------------------------------------------------

test_that(".make_symmetrical creates a symmetrical matrix", {
  set.seed(404)
  mat <- matrix(runif(6), nrow = 2, ncol = 3,
                dimnames = list(c("r1", "r2"), c("c1", "c2", "c3")))

  new_mat <- .make_symmetrical(mat)

  expect_identical(nrow(new_mat), ncol(new_mat))
  expect_identical(rownames(new_mat), colnames(new_mat))
  expect_identical(sum(is.na(new_mat)), as.integer(13)) # 3x3 + 2x2
  expect_true(isSymmetric(new_mat))
  expect_identical(mat, new_mat[rownames(mat), colnames(mat)])
})

test_that(".make_symmetrical tests for unique row and column names", {
  set.seed(404)
  mat <- matrix(runif(6), nrow = 2, ncol = 3,
                dimnames = list(c("r1", "r2"), c("c1", "c2", "r1")))

  expect_error(.make_symmetrical(mat),
               regexp = "row names of `mat`")
})


# .get_removals ----------------------------------------------------------------

test_that(".get_removals removes correlated values", {
  r2 <- r2_matrix()
  removed <- .get_removals(r2, fake_na_vars(r2))

  # "v1" has the highest mean correlation between v1 and v2, "v4" has the
  # highest mean correlation between v3 and v4
  expect_identical(removed, c("v1", "v4"))
})

test_that(".get_removals handles multi-way correlations", {
  r2 <- r2_matrix_multiway()
  removed <- .get_removals(r2, fake_na_vars(r2))

  # v3 is correlated with both v4 and v5, so it gets removed instead of v4
  expect_identical(removed, c("v1", "v3"))
})

test_that(".get_removals uses different R2 threshold", {
  r2 <- r2_matrix()
  removed <- .get_removals(r2, fake_na_vars(r2), R2_threshold = 0.7)

  # Only v1 vs v2 has a high enough correlation to trigger removal
  expect_identical(removed, "v1")
})

test_that(".get_removals prefers variables with more NAs", {
  r2 <- r2_matrix()
  na_vars <- fake_na_vars(r2)
  na_vars["v1"] <- 1
  na_vars["v2"] <- 3
  na_vars["v3"] <- 2

  removed <- .get_removals(r2, na_vars)

  # v2 and v3 get removed due to more NAs, even though they have lower mean
  # correlation than v1 and v4
  expect_identical(removed, c("v2", "v3"))
})

test_that(".get_removals returns highest correlated var when there are equal NAs", {
  r2 <- r2_matrix()
  na_vars <- fake_na_vars(r2)
  na_vars["v1"] <- 3
  na_vars["v2"] <- 3

  removed <- .get_removals(r2, na_vars)

  # v1 and v2 have equal NAs so v1 should get removed as the highest-correlated variable
  expect_identical(removed, c("v1", "v4"))
})

test_that(".get_removals retains variables in `removed`", {
  r2 <- r2_matrix()
  removed <- .get_removals(r2, fake_na_vars(r2), removed = "new_var")

  # "new_var" should be in the vector even though it's not in r2
  expect_identical(removed, c("new_var", "v1", "v4"))
})

test_that(".get_removals returns null on low correlation", {
  r2 <- r2_matrix()
  r2[r2 >= 0.5] <- 0.1

  removed <- .get_removals(r2, fake_na_vars(r2))
  expect_null(removed)
})

test_that(".get_removals returns previous `removed` value on low correlation", {
  r2 <- r2_matrix()
  r2[r2 >= 0.5] <- 0.1

  removed <- .get_removals(r2, fake_na_vars(r2), removed = c("new_var"))
  expect_identical(removed, "new_var")
})

test_that(".get_removals never removes `always_keep` variables", {
  r2 <- r2_matrix()
  removed1 <- .get_removals(r2, fake_na_vars(r2), always_keep = "v1")

  # v2 should get removed instead of v1
  expect_identical(removed1, c("v2", "v4"))

  removed2 <- .get_removals(r2, fake_na_vars(r2), always_keep = c("v1", "v2"))

  # v2 shouldn't get removed either
  expect_identical(removed2, "v4")
})

test_that(".get_removals checks for square and symmetrical matrix", {
  r2 <- r2_matrix()

  expect_error(.get_removals(r2[1:2, ], fake_na_vars(r2)),
               regexp = "`r2_mat` is not square or symmetrical")

  r2_2 <- r2
  colnames(r2_2) <- paste0("c", 1:4)
  expect_error(.get_removals(r2_2, fake_na_vars(r2_2)),
               regexp = "`r2_mat` is not square or symmetrical")

  r2_3 <- r2
  r2_3[1, 2] <- 0.5
  expect_error(.get_removals(r2_3, fake_na_vars(r2_3)),
               regexp = "`r2_mat` is not square or symmetrical")
})


# remove_correlated_covariates -------------------------------------------------

test_that("remove_correlated_covariates removes correlated numeric variables", {
  df <- df_with_correlated_vars()[, c("num_1", "num_2", "num_3")]

  new_df <- remove_correlated_covariates(df, verbose = FALSE)

  # In the example I constructed, num_2 is the most highly correlated with other
  # variables
  expect_identical(colnames(new_df), c("num_1", "num_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates removes correlated categorical variables", {
  df <- df_with_correlated_vars()[, c("cat_1", "cat_2", "cat_3")]

  new_df <- remove_correlated_covariates(df, verbose = FALSE)

  # In the example I constructed, cat_2 is the most highly correlated with the
  # other variables
  expect_identical(colnames(new_df), c("cat_1", "cat_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates removes correlated categorical and numeric variables", {
  df <- df_with_correlated_vars()

  new_df <- remove_correlated_covariates(df, verbose = FALSE)

  # In the example I constructed, cat_1 and cat_2 are both highly correlated
  # with all 3 numeric variables, while cat_3 isn't really correlated with
  # any of them. num_2 gets removed from num_1 vs num_2, then cat_2 from cat_1
  # vs cat_2, then cat_1 from num_3 vs cat_1
  expect_identical(colnames(new_df), c("num_1", "num_3", "cat_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates ignores the `id_col` variable", {
  df <- df_with_correlated_vars()[, c("cat_1", "cat_2", "cat_3")]

  new_df <- remove_correlated_covariates(df, id_cols = "cat_2", verbose = FALSE)

  # Nothing is removed because cat_2 is not compared to anything, and cat_1 and
  # cat_3 aren't very correlated
  expect_identical(colnames(new_df), colnames(df))
  expect_identical(new_df, df)
})

test_that("remove_correlated_covariates handles a vector of `id_cols`", {
  df <- df_with_correlated_vars()[, c("cat_1", "cat_2", "cat_3")]
  df$cat_4 <- df$cat_2

  new_df <- remove_correlated_covariates(df, id_cols = c("cat_2", "cat_4"), verbose = FALSE)

  # Nothing is removed
  expect_identical(new_df, df)
})

test_that("remove_correlated_covariates keeps the `always_keep` column", {
  df <- df_with_correlated_vars()[, c("num_1", "num_2", "num_3")]

  new_df <- remove_correlated_covariates(df, always_keep = "num_2", verbose = FALSE)

  # num_1 gets removed instead of num_2
  expect_identical(colnames(new_df), c("num_2", "num_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates handles a vector of `always_keep` variables", {
  df <- df_with_correlated_vars()[, c("num_1", "num_2", "num_3")]

  new_df <- remove_correlated_covariates(df, always_keep = c("num_1", "num_2"), verbose = FALSE)

  # nothing gets removed
  expect_identical(new_df, df)
})

test_that("remove_correlated_covariates does not test mixed effect column vs numeric column", {
  df <- df_with_correlated_vars()

  new_df <- remove_correlated_covariates(df, mixed_effects = "cat_1", verbose = FALSE)

  # cat_1 would get removed if it wasn't a mixed effect, but it shouldn't get
  # removed here
  expect_identical(colnames(new_df), c("num_1", "num_3", "cat_1", "cat_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates tests mixed effect column vs categorical columns", {
  df <- df_with_correlated_vars()

  new_df <- remove_correlated_covariates(df, mixed_effects = "cat_2", verbose = FALSE)

  # cat_2 is still removed
  expect_identical(colnames(new_df), c("num_1", "num_3", "cat_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates handles vector of mixed effects", {
  df <- df_with_correlated_vars()

  new_df <- remove_correlated_covariates(df, mixed_effects = c("cat_1", "cat_2"), verbose = FALSE)

  # cat_2 is still removed, cat_1 should stay
  expect_identical(colnames(new_df), c("num_1", "num_3", "cat_1", "cat_3"))
  expect_identical(new_df, df[, colnames(new_df)])
})

test_that("remove_correlated_covariates uses different R^2 threshold", {
  df <- df_with_correlated_vars()

  new_df1 <- remove_correlated_covariates(df, R2_threshold = 0.99, verbose = FALSE)

  # nothing should get removed
  expect_identical(new_df1, df)

  new_df2 <- remove_correlated_covariates(df, R2_threshold = 0.1, verbose = FALSE)

  # Everything but num_1 and num_3, which are completely uncorrelated, should get removed
  expect_identical(colnames(new_df2), c("num_1", "num_3"))
  expect_identical(new_df2, df[, colnames(new_df2)])
})

test_that("remove_correlated_covariates prints removed columns when verbose is TRUE", {
  df <- df_with_correlated_vars()

  expect_output(remove_correlated_covariates(df, verbose = TRUE),
                regexp = "Removing 3 columns: num_2, cat_2, cat_1")
})

