

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
context("Run SCA from Hist object")

SCA_MP <- make_MP(SCA, HCR_MSY)
load("tests_results/Data_from_Hist.RData")

out_csv <- list()

sfInit(parallel = TRUE, cpus = 8)
sfLibrary(MSEtool)
sfLibrary(DLMextra)

sfExport(list = c("SCA_MP"))

for(i in 1:length(OMs)) {

  #test_that(OMs[i], {
    message(paste("OM:", OMs[i]))
    my_OM <- get(OMs[i])

    Data <- Data_from_Hist[[i]]
    sfExport(list = c("my_OM", "Data"))
    start_time <- proc.time()
    res <- sfLapply(1:my_OM@nsim, SCA_MP, Data = Data)
	  expect_true(all(vapply(res, inherits, logical(1), "Rec")))

	  message(paste("Timing:", round((proc.time() - start_time)[3], 2), "seconds"))
    TAC <- vapply(res, getElement, numeric(1), "TAC")
	  expect_type(TAC, "double")

	  is_na <- sum(is.na(TAC))
    message(paste0(is_na, " NAs out of ", my_OM@nsim, " (", round(100 * is_na/my_OM@nsim, 1), " %) for ", OMs[i], "\n"))

	  if(is_na > 0) message(paste("Check indices:", paste(which(is.na(TAC)), collapse = " "), "\n"))

    out_csv[[i]] <- TAC


  #})

}
sfStop()
names(out_csv) <- OMs

save(out_csv, file = "SCA_OM.RData")
