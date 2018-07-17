
library(testthat)
library(MSEtool)

context("Delay-Difference in assessment mode")
Red_snapper@sigmaR <- 0.3

test_that("DD_TMB assess model", {
  expect_s4_class(DD_TMB(Data = SimulatedData), "Assessment")
  expect_s4_class(DD_TMB(Data = SimulatedData, SR = "Ricker"), "Assessment")

  res <- expect_s4_class(DD_TMB(Data = Red_snapper), "Assessment")
  expect_s4_class(DD_TMB(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, R0 = seq(0.75, 1.25, 0.025), h = seq(0.95, 1, 0.0025), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})

test_that("DD_SS assess model", {
  res <- expect_s4_class(DD_SS(Data = Red_snapper), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, SR = "Ricker"), "Assessment")

  expect_s4_class(DD_SS(Data = Red_snapper), "Assessment")
  expect_s4_class(DD_SS(Data = Red_snapper, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, R0 = seq(0.75, 1.25, 0.025), h = seq(0.9, 0.99, 0.01), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


