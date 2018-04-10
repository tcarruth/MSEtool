# Internal DLMtool functions that are also needed for MSEtool
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))
mconv <- function (m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))

optimize_TMB_model <- function(obj) {
  # Use hessian for fixed-effects models
  if(is.null(obj$env$random)) {
    opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr, obj$he), error = function(e) as.character(e))
  } else {
    opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr), error = function(e) as.character(e))
  }
  return(opt)
}

get_sdreport <- function(obj, opt) {
  if(is.character(opt)) {
    res <- "Model did not converge with nlminb(). Did not run TMB::sdreport()."
  }
  else {
    res <- tryCatch(sdreport(obj, getReportCovariance = FALSE), error = function(e) as.character(e))
  }
  if(inherits(res, "sdreport") && any(eigen(res$cov.fixed)$value < 0)) {
    res <- "Estimated covariance matrix was not positive definite."
  }
  return(res)
}
