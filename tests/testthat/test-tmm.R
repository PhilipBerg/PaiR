test_that("tmm returns a dataframe", {
  data1 <- data.frame(a = letters, b = runif(length(letters)), c = runif(length(letters)))
  expect_equal(names(data1), names(tmm(data1, "a", load_info = F)))
  expect_equal(typeof(data1), typeof(tmm(data1, "a", load_info = F)))
  expect_equal("list", typeof(tmm(data1, "a", load_info = T)))
  expect_equal(2L, length(tmm(data1, "a", load_info = T)))
  expect_equal(ncol(data1), length(tmm(data1, "a", load_info = F)))
  data2 <- data.frame(a1 = letters, a2 = LETTERS, b = runif(length(letters)), c = runif(length(letters)))
  expect_equal(names(data2), names(tmm(data2, "a1", load_info = F)))
  expect_equal(typeof(data2), typeof(tmm(data2, "a1", load_info = F)))
  expect_equal("list", typeof(tmm(data2, "a1", load_info = T)))
  expect_equal(2L, length(tmm(data2, "a1", load_info = T)))
  expect_equal(ncol(data2), length(tmm(data2, "a1", load_info = F)))
  data3 <- data.frame(
    a1 = letters,
    a2 = LETTERS,
    b = runif(length(letters)),
    c = runif(length(letters)),
    d = runif(length(letters))
  )
  expect_equal(ncol(data3), ncol(tmm(data3, "a1", target = c("b", "c"), load_info = F)))
  expect_equal(ncol(data3), ncol(tmm(data3, "a1", target = c("b", "c"), load_info = F, reference_sample = "d")))
})
