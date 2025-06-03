#' Convert String to Numerical Seed
#'
#' Converts a string or character vector to an integer that can be used to set a
#' random seed.
#'
#' This function is especially useful for things like running a process under
#' different conditions and ensuring that each run/condition gets its own unique
#' but consistent seed.
#'
#' @param string a string or a character vector or list. If the input is a
#'   character vector/list, a seed will be generated for each item in the vector
#'   separately
#'
#' @return an integer, if `string` is a single string, or a vector of integers
#'   the same length as `string`, if `string` is a vector or list
#' @export
#'
#' @examples
#' # Single string
#' seed <- string_to_seed("Step 1")
#' set.seed(seed)
#'
#' # Vector of strings
#' seeds <- string_to_seed(c("Step 1", "Step 2", "Step 3"))
#' set.seed(seeds[1])
string_to_seed <- function(string) {
  if (is.null(string) || length(string) == 0 || any(nchar(string) == 0)) {
    stop("`string` must not be NULL or contain an emtpy string")
  }

  if (!is.character(string) && !(is.list(string) && all(sapply(string, is.character)))) {
    stop("`string` must be a string or character vector")
  }

  # This works for both single strings and vectors
  sapply(string, utf8ToInt, USE.NAMES = FALSE) |>
    colSums()
}
