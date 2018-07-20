

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
context("Run Assess models from Hist object")

load("tests_results/Data_from_Hist.RData")


# Remove years if only one age-class (bug in runMSE)
remove_odd_CAA <- function(Data) {
  for(i in 1:dim(Data@CAA)[1]) {
    for(j in 1:dim(Data@CAA)[2]) {
      if(sum(Data@CAA[i, j, ] != 0) == 1) Data@CAA[i, j, ] <- NA
    }
  }
  return(Data)
}

Data_from_Hist <- lapply(Data_from_Hist, remove_odd_CAA)

sfInit(parallel = TRUE, cpus = 8)
sfLibrary(MSEtool)
sfLibrary(DLMextra)


# Other OMs need to be evaluated separately
for(i in c(1:2, 4, 6, 8:12, 14:18, 20:24)) {

  test_that(OMs[i], {
    message(paste("OM:", OMs[i]))
    my_OM <- get(OMs[i])

    Data <- Data_from_Hist[[i]]
    xvec <- 1:my_OM2@nsim

    sfExport(list = c("my_OM", "Data", "xvec"))
    start_time <- proc.time()

    # Main call: change desired Assess model as necessary here.
    res <- sfLapply(xvec, SCA, Data = Data)

	  expect_true(all(vapply(res, inherits, logical(1), "Assessment")))

	  message(paste("Timing:", round((proc.time() - start_time)[3], 2), "seconds"))

    conv <- vapply(res, function(x) is.character(x@SD), logical(1))
	  is_na <- sum(conv)
    message(paste0(is_na, " didn't converge out of ", my_OM@nsim, " (", round(100 * is_na/my_OM@nsim, 1), " %) for ", OMs[i], "\n"))
	  if(is_na > 0) message(paste("Check indices:", paste(which(conv), collapse = " "), "\n"))

  })

}

sfStop()

