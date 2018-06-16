
library(testthat)
library(MSEtool)

context("SCA in assessment mode")

test_that("SCA assess model", {
  expect_s4_class(SCA(Data = Simulation_1), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, integrate = TRUE, fix_tau = F), "Assessment")
  expect_s4_class(SCA(Data = Simulation_1, SR = "Ricker"), "Assessment")

  res <- expect_s4_class(SCA(Data = Red_snapper), "Assessment")
  expect_s4_class(SCA(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, meanR = seq(1.5, 2.5, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})

test_that("SCA assess model fix_h", {

  res <- expect_s4_class(SCA(Data = Simulation_1, fix_h = TRUE), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, meanR = seq(1.5, 2.5, 0.1), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


test_that("SCA2 assess model", {
  res <- expect_s4_class(SCA2(Data = Simulation_1), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, integrate = TRUE, fix_tau = F), "Assessment")
  expect_s4_class(SCA2(Data = Simulation_1, SR = "Ricker"), "Assessment")
  expect_s4_class(SCA2(Data = Red_snapper, SR = "Ricker"), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, UMSY = seq(0.1, 0.2, 0.01), MSY = seq(50, 150, 10), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})



test_that("SCA3 assess model", {
  res <- expect_s4_class(SCA3(Data = Simulation_1), "Assessment")
  expect_s4_class(SCA3(Data = Simulation_1, integrate = TRUE, fix_tau = F), "Assessment")
  expect_s4_class(SCA3(Data = Simulation_1, SR = "Ricker"), "Assessment")
  expect_s4_class(SCA3(Data = Red_snapper, SR = "Ricker"), "Assessment")


  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, R0 = seq(600, 700, 10), h = seq(0.7, 0.99, 0.025), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})



test_that("SCA3 assess model fix_h", {
  res <- expect_s4_class(SCA3(Data = Simulation_1, fix_h = TRUE), "Assessment")

  expect_equivalent(plot(res, save_figure = FALSE), invisible())

  pro <- profile_likelihood(res, R0 = seq(600, 700, 10), figure = FALSE)
  expect_type(pro, "list")
  expect_true(is.data.frame(pro))

  expect_type(retrospective(res, 5, save_figure = FALSE), "list")

})


