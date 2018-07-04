
library(testthat)
library(MSEtool)

context("Delay-Difference in MSE")

test_that("DD_TMB MP", {
  DD_MSY <- expect_s3_class(make_MP(DD_TMB, HCR_MSY, diagnostic = TRUE), "MP")
  DD_6020 <- expect_s3_class(make_MP(DD_TMB, HCR60_10, diagnostic = TRUE), "MP")

  myMSE <- expect_s4_class(runMSE(MPs = c("DD_MSY", "DD_6020")), "MSE")
  expect_gt(length(ls(DLMenv)), 0L)
})

test_that("DD_SS MP", {
  DDSS_MSY <- expect_s3_class(make_MP(DD_SS, HCR_MSY, diagnostic = TRUE), "MP")
  DDSS_6020 <- expect_s3_class(make_MP(DD_SS, HCR60_10, diagnostic = TRUE), "MP")

  myMSE <- expect_s4_class(runMSE(MPs = c("DDSS_MSY", "DDSS_6020")), "MSE")
  expect_gt(length(ls(DLMenv)), 0L)
})
