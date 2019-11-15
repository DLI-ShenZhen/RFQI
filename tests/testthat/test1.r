library(testthat)
library(RFQI)

test_that("convert_num_to_char is character", {
  expect_equal(convert_num_to_char("1",prefix = "CP",n=4), "CP0001")
})

test_that("convert_char_to_num is numeric", {
  expect_equal(convert_char_to_num("FT0001", prefix="FT"), 1)
})