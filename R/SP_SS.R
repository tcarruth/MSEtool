#' State-space surplus production model with UMSY and MSY as leading parameters
#'
#' A state-space surplus production model that estimates annual deviates in biomass.
#' The model uses time-series of catches and a relative abundance index and coded in TMB.
#' The model is conditioned on catch and estimates a predicted index.
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param report Indicates whether report will be produced for Data object.
#' @param n The exponent of the production function, which is fixed in the model.
#' Typically fixed at n = 2, where the biomass at MSY is half of carrying capacity.
#' @param B1frac The biomass to carrying capacity ratio in the first year of the model.
#' A value of 1 assumes virgin conditions.
#' @param Assessment Indicates whether variables and assessment results will be produced from MP.
#' @return If \code{Assessment = TRUE}, an object of class \code{Assessment}.
#' Otherwise, an object of class \code{Rec}.
#' @note The model is parameterized with UMSY and MSY as leading parameters.
#' Virgin conditions are assumed in the first year of the time series and
#' the production function is assumed to be symmetric.
#' @author Q. Huynh
#' @references Fletcher, Pella and Tomlinson
#' @export SP_SS
#' @seealso \code{\link{SP}}
#' @import TMB
#' @importFrom stats nlminb
#' @importFrom mvtnorm rmvnorm
#' @useDynLib MSEtool
SP_SS <- function(x, Data, reps = 100, n = 2, B1frac = 1, Assessment = FALSE) {
  dependencies = "Data@Cat, Data@Ind, Data@CV_Ind"
  yind <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  Catch <- Data@Cat[x, yind]
  if(any(is.na(Catch))) stop('Model is conditioned on complete catch time series, but there is missing catch.')
  Index <- Data@Ind[x, yind]
  Index[is.na(Index)] <- -1
  n_y <- length(Catch)
  AvC <- mean(Catch, na.rm = TRUE)
  if(!is.na(Data@Mort[x])) UMSYstart <- 1 - exp(-0.5*Data@Mort[x])
  else UMSYstart <- 0.2
  sigmaI <- sdconv(1, Data@CV_Ind[x])

  data <- list(model = "SP_SS", Catch = Catch, Index = Index, n_y = n_y)
  params <- list(logit_UMSY = log(UMSYstart/(1 - UMSYstart)),
                 log_MSY = log(3 * AvC), log_B1frac = log(B1frac), log_n = log(n),
                 log_sigma = log(sigmaI), log_tau = log(0.2),
                 log_B_dev = rep(0, n_y - 1))
  info <- list(Year = Year, data = data, params = params)

  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = list(log_B1frac = factor(NA), log_n = factor(NA),
                              log_sigma = factor(NA)),
                   random = "log_B_dev", DLL = "MSEtool", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  if(reps == 1) TAC <- obj$report()$TAC
  if(reps > 1) {
    SD <- sdreport(obj, getReportCovariance = FALSE)
    samps <- rmvnorm(reps, opt$par, round(SD$cov.fixed, 4))
    TAC <- rep(NA, reps)
    for (i in 1:reps) {
      params.new <- list(logit_UMSY = samps[i, 1], log_MSY = samps[i, 2],
                         log_B1frac = log(1), log_n = log(2), log_sigma = log(sigmaI),
                         log_tau = samps[i, 3], log_B_dev = obj$env$last.par.best[obj$env$random])
      obj.samp <- MakeADFun(data = info$data, parameters = params.new, DLL = "MSEtool")
      TAC[i] <- obj.samp$report()$TAC
    }
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)

  if(Assessment) {
    if(!exists("SD")) SD <- sdreport(obj, getReportCovariance = FALSE)
    return(return_Assessment())
  } else {
    return(Rec)
  }
}
class(SP_SS) <- "Assess"

