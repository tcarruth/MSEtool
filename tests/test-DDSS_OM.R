

# Test that default configurations are very, very robust when running MSEs.
library(testthat)
library(MSEtool)
library(DLMextra)

OMs <- avail("OM")

# Make hist objects once
#Hist <- lapply(OMs, function(x) {message(x); return(runMSE(get(x), Hist = TRUE))})
#Data_from_Hist <- lapply(Hist, getElement, "Data")
#save(Data_from_Hist, file = "tests_results/Data_from_Hist.RData")

# Run test
context("Run DD_SS from Hist object")

DDSS_MP <- make_MP(DD_SS, HCR_MSY)
load("tests_results/Data_from_Hist.RData")

out_csv <- list()

for(i in 1:length(OMs)) {

  #test_that(OMs[i], {
    message(paste("OM:", OMs[i]))
    my_OM <- get(OMs[i])

    Data <- Data_from_Hist[[i]]

    res <- lapply(1:my_OM@nsim, DDSS_MP, Data = Data)
    #res2 <- lapply(1:my_OM@nsim, DD_TMB, Data = Data)
    expect_true(all(vapply(res, inherits, logical(1), "Rec")))

    TAC <- vapply(res, function(x) x@TAC, numeric(1))
    expect_type(TAC, "double")

    #hist(TAC, main = OMs[i])

    is_na <- sum(is.na(TAC))
    message(paste0(is_na, " NAs out of ", my_OM@nsim, " (", round(100 * is_na/my_OM@nsim, 1), " %) for ", OMs[i], "\n"))

    out_csv[[i]] <- TAC

  #})

}

names(out_csv) <- OMs

save(out_csv, file = "tests_results/DDSS_OM.RData")
