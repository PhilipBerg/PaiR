test_that("prnn returns dataframe", {
  data <- data.frame(a = letters, b = runif(length(letters)), c = runif(length(letters)))
  expect_equal(names(data), names(prnn(data, 'a', F)))
  expect_equal(typeof(data), typeof(prnn(data, 'a', F)))
  expect_equal("list", typeof(prnn(data, 'a', T)))
})
