
library(testthat)
library(MSEtool)


context("Surplus production in assessment mode")

test_that("SP assess model", {
  expect_s4_class(SP(Data = swordfish), "Assessment")
  res <- expect_s4_class(SP(Data = swordfish, start = list(dep = 0.95)), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, UMSY = seq(0.05, 0.5, 0.005),
                            MSY = seq(0.1, 2, 0.01) * 1e4, figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


test_that("SP_SS assess model", {
  res <- expect_s4_class(SP_SS(Data = swordfish, start = list(dep = 0.95, tau = 0.1)), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, UMSY = seq(0.05, 0.5, 0.1),
                            MSY = seq(0.1, 2, 0.1) * 1e4, figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})
