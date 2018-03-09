#' State-space surplus production model with UMSY and MSY as leading parameters
#'
#' A state-space surplus production model that estimates annual deviates in biomass.
#' The model uses time-series of catches and a relative abundance index and coded in TMB.
#' The model is conditioned on catch and estimates a predicted index.
#'
#' @param Data A data-limited methods data object
#' @param n The exponent of the production function, which is fixed in the model.
#' Typically fixed at n = 2, where the biomass at MSY is half of carrying capacity.
#' @param B1frac The biomass to carrying capacity ratio in the first year of the model,
#' which is fixed. A value of 1 assumes virgin conditions in the first year.
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @note The model is parameterized with UMSY and MSY as leading parameters.
#' Virgin conditions are assumed in the first year of the time series and
#' the production function is assumed to be symmetric.
#' @author Q. Huynh
#' @references
#' Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson system.
#' Fishery Bulletin 76:515:521.
#'
#' Pella, J. J. and Tomlinson, P. K. 1969. A generalized stock production model.
#' Inter-Am. Trop. Tuna Comm., Bull. 13:419-496.
#' @export SP_SS
#' @seealso \code{\link{SP}}
#' @import TMB
#' @importFrom stats nlminb
#' @importFrom mvtnorm rmvnorm
#' @useDynLib MSEtool
SP_SS <- function(Data, n = 2, B1frac = 1) {
  dependencies = "Data@Cat, Data@Ind, Data@CV_Ind"
  x <- 1
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
  info <- list(Year = Year, data = data, params = params, sigma = sigmaI)

  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = list(log_B1frac = factor(NA), log_n = factor(NA),
                              log_sigma = factor(NA)),
                   random = "log_B_dev", DLL = "MSEtool", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  SD <- sdreport(obj)

  Assessment <- return_Assessment()
  return(Assessment)
}
class(SP_SS) <- "Assess"

