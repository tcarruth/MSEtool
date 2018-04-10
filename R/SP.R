#' Surplus production model with UMSY and MSY as leading parameters
#'
#' A surplus production model that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB.
#' The model is conditioned on catch and estimates a predicted index. The state-space version
#' estimates annual deviates in biomass.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data A data-limited methods data object
#' @param n The exponent of the production function (fixed in the model).
#' Typically fixed at \code{n = 2}, where the biomass at MSY is half of carrying capacity.
#' @param Binit_frac The ratio of biomass to carrying capacity in the first year
#' (fixed in the model). A value of 1 assumes virgin conditions.
#' @param start Optional list of starting values. See details.
#' @param silent (TRUE/FALSE) Whether tracing information is provided by TMB during optimization.
#' Used to aid convergence, help with diagnostics.
#' @param ... Additional arguments (not currently used).
#' @details
#' To provide starting values for the model, a named list can be provided for \code{UMSY} and
#' \code{MSY} via the start argument (see example).
#'
#' For \code{SP_SS}, a start value can also be provided for \code{tau}, the standard deviation
#' of the biomass deviates. The index standard deviation is fixed based on the value of \code{Data@@CV_Ind}.
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
#' @export SP
#' @seealso \link{SP_production}
#' @describeIn SP Fixed effects model
#' @examples
#' data(swordfish)
#'
#' #### Observation-error surplus production model
#' res <- SP(1, swordfish)
#'
#' # Provide starting values
#' start <- list(UMSY = 0.1, MSY = 1e5)
#' res <- SP(1, swordfish, start = start)
#'
#' # Assume B/K = 0.95 in first year of model
#' res <- SP(1, swordfish, Binit_frac = 0.95)
#'
#' #### State-space version
#' res <- SP_SS(1, swordfish)
#'
#' @import TMB
#' @importFrom stats nlminb
#' @useDynLib MSEtool
SP <- function(x = 1, Data, n = 2, Binit_frac = 1, start = NULL, silent = TRUE, ...) {
  dependencies = "Data@Cat, Data@Ind"
  ystart <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- ystart:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist))) stop('Model is conditioned on complete catch time series, but there is missing catch.')
  I_hist <- Data@Ind[x, yind]
  I_hist[I_hist < 0] <- NA
  ny <- length(C_hist)

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- log(start$UMSY[1]/(1 - start$UMSY[1]))
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- ifelse(is.na(Data@Mort[x]), 0.2, 1 - exp(-0.5 * Data@Mort[x]))
    params$logit_UMSY <- log(UMSY_start/(1 - UMSY_start))
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist, na.rm = TRUE)
    params$log_MSY <- log(3 * AvC)
  }

  data <- list(model = "SP", C_hist = C_hist, I_hist = I_hist, ny = ny)
  params$log_Binit_frac <- log(Binit_frac)
  params$log_n <- log(n)
  info <- list(Year = Year, data = data, params = params)

  map <- list(log_Binit_frac = factor(NA), log_n = factor(NA))
  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, DLL = "MSEtool", silent = silent)
  opt <- optimize_TMB_model(obj)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SP",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
    warning("Model did not properly converge. Check TMB objects (slots names: opt and SD).")
  } else {
    Yearplusone <- c(Year, max(Year) + 1)

    Assessment <- new("Assessment", Model = "SP",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY, VBMSY = report$BMSY,
                      B0 = report$K, VB0 = report$K, U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$K, names = Yearplusone),
                      VB_VBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      VB_VB0 = structure(report$B/report$K, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Index = structure(report$Ipred, names = Year),
                      NLL = structure(c(opt$objective, report$nll), names = c("Total", "Index")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[7],
                      SE_B_BMSY_final = SD$sd[8], SE_B_B0_final = SD$sd[9],
                      SE_VB_VBMSY_final = SD$sd[8], SE_VB_VB0_final = SD$sd[9],
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }

  return(Assessment)
}
class(SP) <- "Assess"


#' @describeIn SP State-space version
#' @export SP_SS
#' @import TMB
#' @importFrom stats nlminb
#' @useDynLib MSEtool
SP_SS <- function(x = 1, Data, n = 2, Binit_frac = 1, start = NULL, silent = TRUE, ...) {
  dependencies = "Data@Cat, Data@Ind, Data@CV_Ind"
  yind <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist))) stop('Model is conditioned on complete catch time series, but there is missing catch.')
  I_hist <- Data@Ind[x, yind]
  I_hist[I_hist < 0] <- NA
  ny <- length(C_hist)
  sigmaI <- sdconv(1, Data@CV_Ind[x])

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- log(start$UMSY[1]/(1 - start$UMSY[1]))
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- ifelse(is.na(Data@Mort[x]), 0.2, 1 - exp(-0.5 * Data@Mort[x]))
    params$logit_UMSY <- log(UMSY_start/(1 - UMSY_start))
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist, na.rm = TRUE)
    params$log_MSY <- log(3 * AvC)
  }
  if(is.null(params$log_tau)) params$log_tau <- log(0.3)

  data <- list(model = "SP_SS", C_hist = C_hist, I_hist = I_hist, ny = ny)
  params$log_Binit_frac <- log(Binit_frac)
  params$log_n <- log(n)
  params$log_sigma <- log(sigmaI)
  params$log_B_dev = rep(0, ny - 1)

  info <- list(Year = Year, data = data, params = params, sigma = sigmaI)
  map <- list(log_Binit_frac = factor(NA), log_n = factor(NA), log_sigma = factor(NA))
  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, random = "log_B_dev", DLL = "MSEtool", silent = silent)
  opt <- optimize_TMB_model(obj)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SP_SS",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
    warning("Model did not properly converge. Check TMB objects (slots names: opt and SD).")
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearrandom <- Year[2]:max(Year)

    Assessment <- new("Assessment", Model = "SP_SS",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY, VBMSY = report$BMSY,
                      B0 = report$K, VB0 = report$K, U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$K, names = Yearplusone),
                      VB_VBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      VB_VB0 = structure(report$B/report$K, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Index = structure(report$Ipred, names = Year),
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Biomass deviations",
                      NLL = structure(c(opt$objective, report$nll_comp), names = c("Total", "Index", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[7],
                      SE_B_BMSY_final = SD$sd[8], SE_B_B0_final = SD$sd[9],
                      SE_VB_VBMSY_final = SD$sd[8], SE_VB_VB0_final = SD$sd[9],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SP_SS) <- "Assess"
