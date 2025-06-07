#' Remove Covariates that are Not Usable for Regression
#'
#' Given a data frame of metadata, this function removes columns that are either
#' all `NA` values, all the same value, or are not numeric and have a unique
#' value for each sample. Columns that have all the same value after removal of
#' `NA` values are also removed.
#'
#' @param data a data.frame or data.frame-like object of covariate data where
#'   rows are samples and columns are covariates
#' @param always_keep (optional) a single string or a character vector of
#'   variables to never remove even if they meet the criteria for removal, for
#'   example a column of IDs. Defaults to an empty vector (removes everything
#'   that meets the criteria).
#' @param verbose (optional) whether to print out which variables have been
#'   removed before returning. Defaults to `TRUE`.
#'
#' @return a data.frame that is a copy of `data` but without columns that meet
#'   the removal criteria
#' @export
#' @importFrom dplyr select summarize across where any_of everything
#'
#' @examples
#' df <- data.frame(
#'   id = 1:10,
#'   valid = rep(TRUE, 10),
#'   group = c(rep(1, 5), rep(2, 5))
#' )
#' df <- remove_unusable_covariates(df, always_keep = c("id"))
remove_unusable_covariates <- function(data, always_keep = c(), verbose = TRUE) {
  data <- as.data.frame(data)

  # Columns that have only one unique value or are all NA
  same_value <- data |>
    select(-any_of(always_keep)) |>
    summarize(across(everything(), ~ length(na.omit(unique(.x))) <= 1))

  # Columns that have a unique value for each sample
  all_unique <- data |>
    select(-any_of(always_keep), -where(is.numeric)) |>
    summarize(across(everything(), ~ length(unique(.x)) == nrow(data)))

  to_remove <- c(colnames(same_value)[same_value == TRUE],
                 colnames(all_unique)[all_unique == TRUE]) |>
    unique()

  if (verbose) {
    print(paste("Removing", length(to_remove), "columns:",
                paste(to_remove, collapse = ", ")))
  }

  data |> select(-any_of(to_remove))
}


#' Remove Highly Correlated Covariates
#'
#' This function examines the correlation between all given covariates and
#' removes at least one covariate from all pairs of highly-correlated covariates
#' in order to reduce the total correlation between variables in a regression
#' formula. [remove_unusable_covariates()] should be run before using this
#' function, to remove variables that are all `NA` or all the same value.
#'
#' @param id_cols (optional) a string or character vector of columns that should
#'   be treated like IDs and never compared to other covariates or considered
#'   for removal. Defaults to an empty vector.
#' @param mixed_effects (optional) a string or character vector of columns that
#'   will be mixed effects in a regression model. These variables are compared
#'   to other categoricals for potential removal but are not compared with
#'   numeric variables. Defaults to an empty vector.
#' @param always_keep (optional) a string or character vector of columns that
#'   should always be kept even if they are highly correlated with another
#'   variable. In that case, the other variable is removed instead, as long as
#'   it is not also in `always_keep`. Defaults to an empty vector.
#' @param R2_threshold (optional) the threshold for R^2, where variables with an
#'   R^2 value > `R2_threshold` are considered for removal. Defaults to 0.5.
#' @inheritParams remove_unusable_covariates
#'
#' @return a data.frame that is a copy of `data` but without columns that meet
#'   the removal criteria
#' @export
#' @importFrom dplyr select summarize across where any_of everything
#'
#' @examples
#' \dontrun{
#' }
remove_correlated_covariates <- function(data,
                                         id_cols = c(),
                                         mixed_effects = c(),
                                         always_keep = c(),
                                         R2_threshold = 0.5,
                                         verbose = TRUE) {
  data <- as.data.frame(data)

  # Number of NA entries in each column
  na_vars <- data |>
    summarize(across(everything(), ~sum(is.na(.x))))

  data_num <- data |> select(where(is.numeric), -any_of(id_cols))
  data_cat <- data |> select(where(is.character),
                             where(is.factor),
                             where(is.logical),
                             -any_of(id_cols))

  to_remove <- c()

  # Correlation between numerical values
  if (ncol(data_num) > 1) {
    r2_mat <- stats::cor(data_num, use = "na.or.complete")^2
    to_remove <- .get_removals(r2_mat, na_vars, R2_threshold, to_remove, always_keep)
  }

  # Cramer's V between categorical values
  if (ncol(data_cat) > 1) {
    cv_mat <- sapply(colnames(data_cat), function(col_name1) {
      sapply(colnames(data_cat), function(col_name2) {
        table(data_cat[, col_name1], data_cat[, col_name2]) |>
          DescTools::CramerV()
      })
    })
    cv_mat <- cv_mat^2 # Similar to R^2

    to_remove <- .get_removals(cv_mat, na_vars, R2_threshold, to_remove, always_keep)
  }

  # R^2 between categorical and numerical values via linear fit, excluding
  # any mixed-effect variables and any variables already getting removed
  data_num <- data_num |> select(-any_of(to_remove))
  data_cat <- data_cat |> select(-any_of(to_remove), -any_of(mixed_effects))

  if (ncol(data_num) > 1 && ncol(data_cat) > 1) {
    lm_mat <- sapply(colnames(data_cat), function(col_cat) {
      sapply(colnames(data_num), function(col_num) {
        fit <- stats::lm(data[, col_num] ~ data[, col_cat])
        r2 <- summary(fit)$r.squared
      })
    })

    to_remove <- .get_removals(lm_mat, na_vars, R2_threshold, to_remove, always_keep)
  }

  if (verbose) {
    print(paste("Removing", length(to_remove), "columns:",
                paste(to_remove, collapse = ", ")))
  }

  data |> select(-any_of(to_remove))
}

#' Get Which Variables to Remove Based on Correlation
#'
#' This is a helper function for [remove_correlated_covariates()] to decide
#' which variable in each highly-correlated pair of covariates should get
#' removed.
#'
#' @param cor_mat a matrix of correlation (or other association-like) values,
#'   which must have row and column names. This matrix does not need to be
#'   square.
#' @param na_vars a one-row data.frame where the columns are the covariates and
#'   the values are the number of `NA` values in each column.
#' @param removed (optional) a string or character vector of variables that have
#'   already been marked for removal. In the case where one covariate in a pair
#'   of highly-correlated variables is in `removed` already, the other variable
#'   will not be removed even if it would otherwise have met the criteria for
#'   removal. Defaults to an empty vector.
#' @inheritParams remove_correlated_covariates
#'
#' @return a character vector with the names of the columns that should be
#'   removed. May also be an empty vector.
#'
#' @examples
#' \dontrun{
#' }
.get_removals <- function(cor_mat, na_vars, R2_threshold = 0.5,
                         removed = c(), always_keep = c()) {
  if (nrow(cor_mat) == ncol(cor_mat) &&
      all(rownames(cor_mat) == colnames(cor_mat))) {
    diag(cor_mat) <- NA # Remove self-correlation
  }

  r2_melt <- cor_mat

  if (nrow(cor_mat) == ncol(cor_mat) &&
      all(rownames(cor_mat) == colnames(cor_mat))) {
    r2_melt[upper.tri(r2_melt, diag = TRUE)] <- 0  # Avoids picking up both (a vs b) and (b vs a)
  }

  r2_melt <- r2_melt |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(cols = -var1,
                        names_to = "var2", values_to = "value") |>
    subset(value >= R2_threshold) |>  # pairs with R^2 > 0.5 (cor ~ 0.7) only
    dplyr::arrange(dplyr::desc(value))

  if (nrow(r2_melt) == 0) {
    return(removed)
  }

  for (R in 1:nrow(r2_melt)) {
    vars <- as.character(c(r2_melt$var1[R], r2_melt$var2[R]))

    if (any(vars %in% removed)) {
      # No need to re-check if we're already removing one of these variables
      next
    } else if (any(vars %in% always_keep)) {
      # If either variable is in the always_keep list, only remove the variable
      # that isn't in always_keep. If both are in the list, neither one gets
      # removed.
      removed <- c(removed, setdiff(vars, always_keep))
    } else if (any(na_vars[1, vars] > 0) && (na_vars[1, vars[1]] != na_vars[1, vars[2]])) {
      # If either variable has any NA values, and they don't have the same number
      # of NAs, remove the one with the most NAs
      removed <- c(removed,
                   vars[which.max(na_vars[1, vars])])
    } else {
      # No NAs, no variables are in `always_keep` or `removed`
      cur_vars <- setdiff(rownames(cor_mat), removed)
      mean_r2 <- rowMeans(cor_mat[cur_vars, cur_vars], na.rm = TRUE)

      # Remove the variable with the largest mean R^2 with the other remaining variables
      removed <- c(removed,
                   vars[which.max(mean_r2[vars])])
    }
  }

  return(unique(removed))
}


