# Internal DLMtool functions that are also needed for MSEtool
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))
#mconv <- function (m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))

optimize_TMB_model <- function(obj, control = list()) {
  # Use hessian for fixed-effects models
  if(is.null(obj$env$random)) hess <- obj$he else hess <- NULL
  opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr, hess,
                         control = control), error = function(e) as.character(e))
  return(opt)
}

get_sdreport <- function(obj, opt) {
  if(is.character(opt)) {
    res <- "nlminb() optimization returned an error. Could not run TMB::sdreport()."
  } else {
    res <- tryCatch(sdreport(obj, par.fixed = opt$par, getReportCovariance = FALSE),
                    error = function(e) as.character(e))
  }
  if(inherits(res, "sdreport") && !res$pdHess) {
    res <- "Estimated covariance matrix was not positive definite."
  }
  return(res)
}


# Start random effects estimation only after the first year in which
# data from the criterion_vector (e.g. index, CAA) is available
random_map <- function(criterion_vector) {
  ind <- which(!is.na(criterion_vector))[1]
  map <- ifelse(1:length(criterion_vector) <= ind, NA, 1)
  nrandom <- sum(!is.na(map))
  map[!is.na(map)] <- 1:nrandom
  map <- map[2:length(map)]
  return(factor(map))
}

# If there are fewer years of CAA/CAL than Year, add NAs to matrix
expand_comp_matrix <- function(Data, comp_type = c("CAA", "CAL")) {
  comp_type <- match.arg(comp_type)
  ny <- length(Data@Year)

  comp <- slot(Data, comp_type)
  dim_comp <- dim(comp)
  ny_comp <- dim_comp[2]
  if(ny_comp < ny) {
    newcomp <- array(NA, dim = c(1, ny, dim_comp[3]))
    ind_new <- ny - ny_comp + 1
    newcomp[ , ind_new:ny, ] <- comp
    slot(Data, comp_type) <- newcomp
  }
  if(ny_comp > ny) {
    ind_new <- ny_comp - ny + 1
    newcomp <- comp[, ind_new:ny, ]
    slot(Data, comp_type) <- newcomp
  }
  return(Data)
}

