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
#' Correlation is calculated using [variancePartition::canCorPairs()].
#'
#' @param id_cols (optional) a string or character vector of columns that should
#'   be treated like IDs and never compared to other covariates or considered
#'   for removal. Defaults to an empty vector.
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
                                         always_keep = c(),
                                         R2_threshold = 0.5,
                                         verbose = TRUE) {
  data <- as.data.frame(data)

  # Number of NA entries in each column
  na_vars <- data |>
    summarize(across(everything(), ~sum(is.na(.x))))

  cols_analyze <- setdiff(colnames(data), id_cols)

  form <- paste("~", paste(cols_analyze, collapse = " + "))

  cor_mat <- variancePartition::canCorPairs(form, data, showWarnings = FALSE)
  r2_mat <- cor_mat^2

  to_remove <- .get_removals(r2_mat, na_vars, R2_threshold, always_keep)

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
#' @param r2_mat a matrix of R^2 values (or other positive association-like
#'   values), which must have row and column names. This matrix must be square.
#' @param na_vars a one-row data.frame where the columns are the covariates and
#'   the values are the number of `NA` values in each column.
#' @inheritParams remove_correlated_covariates
#'
#' @return a character vector with the names of the columns that should be
#'   removed. May also be an empty vector.
.get_removals <- function(r2_mat, na_vars, R2_threshold = 0.5, always_keep = c()) {
  if (!(nrow(r2_mat) == ncol(r2_mat)) ||
      !(all(rownames(r2_mat) == colnames(r2_mat))) ||
      !isSymmetric(r2_mat)) {
    stop("`r2_mat` is not square or symmetrical")
  }

  # Remove self-correlation
  diag(r2_mat) <- NA

  r2_melt <- r2_mat

  # Avoids picking up both (a vs b) and (b vs a)
  r2_melt[upper.tri(r2_melt, diag = TRUE)] <- 0

  r2_melt <- r2_melt |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(cols = -var1,
                        names_to = "var2", values_to = "value") |>
    subset(value >= R2_threshold) |>  # pairs with R^2 > 0.5 (cor ~ 0.7) only
    dplyr::arrange(dplyr::desc(value))

  if (nrow(r2_melt) == 0) {
    return(c())
  }

  removed <- c()

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
      cur_vars <- setdiff(rownames(r2_mat), removed)
      mean_r2 <- rowMeans(r2_mat[cur_vars, cur_vars], na.rm = TRUE)

      # Remove the variable with the largest mean R^2 with the other remaining variables
      removed <- c(removed,
                   vars[which.max(mean_r2[vars])])
    }
  }

  return(unique(removed))
}
