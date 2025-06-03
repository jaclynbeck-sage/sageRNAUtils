# string_to_seed ---------------------------------------------------------------

test_that("string_to_seed converts a single string to an integer", {
  seed <- string_to_seed("1234")

  expect_identical(seed, sum(49, 50, 51, 52))
})

test_that("string_to_seed converts a vector of strings to separate integers", {
  seeds <- string_to_seed(c("test1", "test2", "test3"))

  expect_length(seeds, 3)
  expect_identical(seeds, c(497, 498, 499))
})

test_that("string_to_seed works on lists of strings", {
  seeds <- string_to_seed(list("test1", "test2", "test3"))

  expect_length(seeds, 3)
  expect_identical(seeds, c(497, 498, 499))
})

test_that("string_to_seed checks for non-character single input", {
  expect_error(string_to_seed(12),
               regexp = "must be a string")

  expect_error(string_to_seed(TRUE),
               regexp = "must be a string")
})

test_that("string_to_seed checks for non-character vector input", {
  expect_error(string_to_seed(c(12, 34, 56)),
               regexp = "must be a string")

  expect_error(string_to_seed(c(TRUE, FALSE)),
               regexp = "must be a string")
})

test_that("string_to_seed checks for NULL single input", {
  expect_error(string_to_seed(NULL),
               regexp = "must not be NULL")
})

test_that("string_to_seed checks for empty vector input", {
  expect_error(string_to_seed(c()),
               regexp = "must not be NULL")
})

test_that("string_to_seed checks for empty string single input", {
  expect_error(string_to_seed(""),
               regexp = "must not be NULL")
})

test_that("string_to_seed checks for empty string in vector input", {
  expect_error(string_to_seed(c("test1", "", "test2")),
               regexp = "must not be NULL")
})
