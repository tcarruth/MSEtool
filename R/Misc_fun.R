# Internal DLMtool functions that are also needed for MSEtool
iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))
mconv <- function (m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))

logit <- function(p) log(p/(1 - p))
ilogit <- function(x) 1/(1 + exp(-x))
ilogitm <- function(x) exp(x)/apply(exp(x), 1, sum)

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
    if(is.null(obj$env$random)) hess <- obj$he(opt$par) else hess <- NULL
    res <- tryCatch(sdreport(obj, par.fixed = opt$par, hessian.fixed = hess,
                             getReportCovariance = FALSE), error = function(e) as.character(e))
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

# var_div - report variables which are divided by the catch rescale
# var_mult - report variables which are multiplied by the catch rescale
# var_trans - transformed variables which need to be rescaled
# fun_trans - the function for rescaling the transformed variables (usually either "*" or "/")
# fun_fixed - the transformation from the output variable to the estimated variable indicated in var_trans (e.g. log, logit, NULL)
rescale_report <- function(var_div, var_mult, var_trans = NULL, fun_trans = NULL, fun_fixed = NULL) {
  output <- mget(c("report", "rescale", "SD"), envir = parent.frame(), ifnotfound = list(NULL))
  report <- output$report

  if(!is.null(var_div)) report[var_div] <- lapply(report[var_div], "/", output$rescale)
  if(!is.null(var_mult)) report[var_mult] <- lapply(report[var_mult], "*", output$rescale)
  assign("report", report, envir = parent.frame())

  if(!is.null(output$SD) && !is.character(output$SD)) {
    SD <- output$SD
    if(!is.null(var_trans)) {
      for(i in 1:length(var_trans)) {
        var_trans2 <- var_trans[i]
        fun_trans2 <- fun_trans[i]
        fun_fixed2 <- fun_fixed[i]

        ind <- pmatch(var_trans2, names(SD$value))
        SD$value[ind] <- do.call(match.fun(fun_trans2), list(SD$value[ind], output$rescale))
        SD$sd[ind] <- do.call(match.fun(fun_trans2), list(SD$sd[ind], output$rescale))

        if(!is.na(fun_fixed2)) {
          fixed_name <- paste0(fun_fixed2, "_", var_trans2)
          ind_fixed <- pmatch(fixed_name, names(SD$par.fixed))
          SD$par.fixed[ind_fixed] <- do.call(match.fun(fun_fixed2), list(SD$value[var_trans2]))
        }
      }
    }
    assign("SD", SD, envir = parent.frame())
  }
  invisible()
}
