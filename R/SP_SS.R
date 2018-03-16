#' State-space surplus production model with UMSY and MSY as leading parameters
#'
#' A state-space surplus production model that estimates annual deviates in biomass.
#' The model uses time-series of catches and a relative abundance index and coded in TMB.
#' The model is conditioned on catch and estimates a predicted index.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data A data-limited methods data object
#' @param n The exponent of the production function (fixed in the model).
#' Typically fixed at \code{n = 2}, where the biomass at MSY is half of carrying capacity.
#' @param Binit_frac The ratio of biomass to carrying capacity in the first year
#' (fixed in the model). A value of 1 assumes virgin conditions.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @note The model uses the Fletcher (1978) formulation and is parameterized with UMSY and MSY as
#' leading parameters. The default conditions assume virgin conditions in the first year of the time series
#' and a symmetric production function.
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
SP_SS <- function(x = 1, Data, n = 2, Binit_frac = 1, ...) {
  dependencies = "Data@Cat, Data@Ind, Data@CV_Ind"
  yind <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist))) stop('Model is conditioned on complete catch time series, but there is missing catch.')
  I_hist <- Data@Ind[x, yind]
  I_hist[I_hist < 0] <- NA
  ny <- length(C_hist)
  AvC <- mean(C_hist, na.rm = TRUE)
  if(!is.na(Data@Mort[x])) UMSYstart <- 1 - exp(-0.5 * Data@Mort[x])
  else UMSYstart <- 0.2
  sigmaI <- sdconv(1, Data@CV_Ind[x])

  data <- list(model = "SP_SS", C_hist = C_hist, I_hist = I_hist, ny = ny)
  params <- list(logit_UMSY = log(UMSYstart/(1 - UMSYstart)),
                 log_MSY = log(3 * AvC), log_Binit_frac = log(Binit_frac), log_n = log(n),
                 log_sigma = log(sigmaI), log_tau = log(0.2),
                 log_B_dev = rep(0, ny - 1))
  info <- list(Year = Year, data = data, params = params, sigma = sigmaI)

  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = list(log_Binit_frac = factor(NA), log_n = factor(NA),
                              log_sigma = factor(NA)),
                   random = "log_B_dev", DLL = "MSEtool", silent = TRUE)
  opt <- optimize_TMB_model(obj)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || opt$convergence != 0 || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SP_SS",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
    warning("Model did not properly converge. Check TMB objects (slots names: opt and SD).")
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearrandom <- Year[2]:max(Year)

    Assessment <- new("Assessment", Model = "SP_SS",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY,
                      B0 = report$K, U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$K, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Index = structure(report$Ipred, names = Year),
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Biomass deviations",
                      NLL = structure(c(opt$objective, report$nll_comp), names = c("Total", "Index", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[7],
                      SE_B_BMSY_final = SD$sd[8], SE_B_B0_final = SD$sd[9],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }

  return(Assessment)
}
class(SP_SS) <- "Assess"

