test_that("automatic core detection distinguishes batch and single-sample work", {
  local_mocked_bindings(detectCores=function(logical=TRUE)if(logical)16L else 8L,
    .package="parallel")
  expect_equal(GCMS_threads(NULL),8L)
  expect_equal(GCMS_threads(NULL,logical=TRUE),16L)
  expect_equal(GCMS_threads(4L,logical=TRUE),4L)
  expect_error(GCMS_threads(0L,logical=TRUE),"positive integer")
})
