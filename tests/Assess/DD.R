
library(testthat)
library(MSEtool)

context("Delay-Difference in assessment mode")

test_that("DD_TMB assess model", {
  expect_s4_class(DD_TMB(Data = Simulation_1), "Assessment")
  expect_s4_class(DD_TMB(Data = Simulation_1, SR = "Ricker"), "Assessment")

  res <- expect_s4_class(DD_TMB(Data = sim_snapper), "Assessment")
  expect_s4_class(DD_TMB(Data = sim_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, UMSY = seq(0.05, 0.15, 0.01),
                            MSY = seq(1, 3, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})

test_that("DD_SS assess model", {
  res <- expect_s4_class(DD_SS(Data = SimulatedData), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = SimulatedData, SR = "Ricker"), "Assessment")

  expect_s4_class(DD_SS(Data = Simulation_1, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = Simulation_1, SR = "Ricker"), "Assessment")

  expect_s4_class(DD_SS(Data = sim_snapper), "Assessment")
  expect_s4_class(DD_SS(Data = sim_snapper, integrate = TRUE), "Assessment")
  expect_s4_class(DD_SS(Data = sim_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, UMSY = seq(0.05, 0.15, 0.01),
                            MSY = seq(1, 3, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


