#' Generalized correlation
#'
#' Calculates correlation or correlation-like values between numeric variables,
#' categorical variables, and combinations of both. Categorical vs categorical
#' is calculated with Cramer's V, and categorical vs numeric is calculated as
#' the square root of the R^2 value from a linear model of `numeric ~ categorical`.
#'
#' @param data a data.frame, matrix, or object coercible to a data.frame, where
#'   rows are observations/samples and columns are variables.
#' @param exclude_cols (optional) a string or character vector containing the
#'   names of any columns in `data` to exclude from this calculation, for
#'   example an ID column. Defaults to `NULL`.
#'
#' @return an NxN symmetric matrix, where N = the number of columns of `data`
#'   minus any `exclude_cols`. Values will be the correlation or Cramer's V
#'   between each variable. Note that values for categorical vs categorical and
#'   categorical vs numerical will always be positive, but numeric vs numeric
#'   correlation may be positive or negative.
#' @export
#'
#' @examples
#' df <- data.frame(var1 = rnorm(10),
#'                  var2 = c(rnorm(5, mean = 0), rnorm(5, mean = 10)), # Correlated with group
#'                  var3 = c(rep("group1", 5), rep("group2", 5)),
#'                  var4 = c(rep("dose1", 3), rep("dose2", 3), rep("dose3", 4)))
#' generalized_correlation(df)
generalized_correlation <- function(data, exclude_cols = NULL) {
  data <- as.data.frame(data)

  numerical_vars <- dplyr::select(data, where(is.numeric), -any_of(exclude_cols))
  categorical_vars <- dplyr::select(data, !where(is.numeric), -any_of(exclude_cols)) |>
    dplyr::mutate(across(where(is.factor), droplevels)) # in case some levels have no data

  num_cor_mat <- cat_cv_mat <- lm_mat <- matrix()

  # Calculate correlation between numeric variables
  if (ncol(numerical_vars) > 0) {
    num_cor_mat <- stats::cor(numerical_vars, use = "na.or.complete")
  }

  # Calculate CramerV between categorical variables
  if (ncol(categorical_vars) > 0) {
    cat_cv_mat <- .categorical_cor(categorical_vars)
  }

  # Calculate correlation between pairs of numeric + categorical variables
  if (ncol(numerical_vars) > 0 && ncol(categorical_vars) > 0) {
    lm_mat <- .categorical_vs_numeric_cor(numerical_vars, categorical_vars)

    # Make lm_mat symmetrical and insert correlation/cv values in the empty
    # parts of the matrix. This works even if num_cor_mat or cat_cv_mat are empty.
    lm_mat <- .make_symmetrical(lm_mat)

    cor_vars <- intersect(rownames(lm_mat), rownames(num_cor_mat))
    cv_vars <- intersect(rownames(lm_mat), rownames(cat_cv_mat))

    lm_mat[cor_vars, cor_vars] <- num_cor_mat[cor_vars, cor_vars]
    lm_mat[cv_vars, cv_vars] <- cat_cv_mat[cv_vars, cv_vars]

    return(lm_mat)

  } else if (ncol(numerical_vars) > 0 && ncol(categorical_vars) == 0) {
    return(num_cor_mat) # No categorical variables

  } else {
    return(cat_cv_mat) # No numeric variables
  }
}


.categorical_cor <- function(categorical_df) {
  # For each column col_name1
  sapply(colnames(categorical_df), function(col_name1) {
    # Compare col_name1 against each column col_name2. col_name1 will get
    # compared to itself here, which should give a value of 1
    sapply(colnames(categorical_df), function(col_name2) {
      # Make contingency table and run CramerV
      table(categorical_df[, col_name1], categorical_df[, col_name2]) |>
        DescTools::CramerV()
    })
  })
}


.categorical_vs_numeric_cor <- function(numerical_df, categorical_df) {
  # For each categorical column col_cat
  lm_mat <- sapply(colnames(categorical_df), function(col_cat) {
    # Compare it to each numerical column col_num by fitting a linear model
    # and getting the model's R^2 value
    sapply(colnames(numerical_df), function(col_num) {
      fit <- stats::lm(numerical_df[, col_num] ~ categorical_df[, col_cat])
      r2 <- summary(fit)$r.squared
      sqrt(r2)
    })
  })
}


.make_symmetrical <- function(mat) {
  new_names <- c(rownames(mat), colnames(mat))

  if (length(unique(new_names)) != length(new_names)) {
    stop("row names of `mat` are not distinct from col names")
  }

  new_mat <- matrix(NA,
                    nrow = nrow(mat) + ncol(mat),
                    ncol = nrow(mat) + ncol(mat),
                    dimnames = list(new_names, new_names))

  # Insert mat at the right rows and columns
  new_mat[rownames(mat), colnames(mat)] <- mat

  # Insert the transpose at the right rows and columns
  new_mat[colnames(mat), rownames(mat)] <- t(mat)

  return(new_mat)
}
