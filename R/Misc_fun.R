# Internal DLMtool functions that are also needed for MSEtool
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))


# Test convergence called in MP
check_convergence <- function() {
  obj <- get("obj", envir = parent.frame())
  SD <- try(sdreport(obj))
  if(inherits(SD, "try-error")) {
    warning("Convergence issue: could not obtain covariance matrix.")
  } else {
    if(any(eigen(SD$cov.fixed)$value < 0)) warning("Convergence issue: Hessian is not positive definite.")
  }
  return(invisible())
}

# Called during profile_likelihood. Add for state-space models.
profile_elapsed_time <- function() {
  output <- mget(c("Assessment", "nll", "i"), envir = parent.frame())
  MP <- output$Assessment@MP
  N <- length(output$nll)
  if(output$i %in% floor(seq(0, 1, 0.1) * N)) {
    message(paste0("Likelihood profile for ", paste0(MP), ": ", i, "/", N, " (",
                   i/N, "%) completed."))
  }
  invisible()
}
