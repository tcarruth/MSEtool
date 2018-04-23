#' Delay - Difference Stock Assessment in TMB with UMSY and MSY as leading parameters
#'
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB. The model
#' is conditioned on effort and estimates predicted catch. In the state-space version,
#' recruitment deviations from the stock-recruit relationship are estimated as a
#' random effect variable.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 when running an assessment.
#' @param Data An object of class \linkS4class{Data}.#'
#' @param start Optional list of starting values. See details.
#' @param silent (TRUE/FALSE) Passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @details
#' To provide starting values for the model, a named list can be provided for \code{UMSY},
#' \code{MSY}, and \code{q} via the \code{start} argument (see example).
#'
#' For \code{DD_SS}, a start value can also be provided for \code{tau}, the standard
#' deviation of the recruitment variability. The catch standard deviation is fixed based
#' on the value of \code{Data@@CV_Cat}.
#' @note Similar to many other assessment
#' models, the model depends on assumptions such as stationary productivity and
#' proportionality between the abundance index and real abundance.
#' Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers & Z. Siders. Zach Siders coded the TMB function.
#' @references
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#'
#' Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York.
#' @describeIn DD_TMB Observation-error only model
#' @import TMB
#' @importFrom stats nlminb
#' @examples
#' data(sim_snapper)
#'
#' #### Observation-error delay difference model
#' res <- DD_TMB(1, sim_snapper)
#'
#' ### State-space version
#' res <- DD_SS(1, sim_snapper)
#'
#' # Provide starting values
#' start <- list(UMSY = 0.05, MSY = 4, q = 0.3, tau = 0.5)
#' res <- DD_TMB(1, sim_snapper, start = start)
#'
#' summary(res@@SD) # Look at parameter estimates
#' @useDynLib MSEtool
#' @export
DD_TMB <- function(x, Data, start = NULL, silent = TRUE, control = list(), ...) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x,   ]))[1]
  yind <- ystart:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  E_hist <- C_hist/I_hist
  if(any(is.na(E_hist))) stop("Missing values in catch and index in Data object.")
  E_hist <- E_hist/mean(E_hist)
  ny <- length(C_hist)
  k <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k[k > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho <- (wa[k + 2] - Winf)/(wa[k + 1] - Winf)
  Alpha <- Winf * (1 - Rho)
  S0 <- exp(-Data@Mort[x])  # get So survival rate
  wk <- wa[k]

  data <- list(model = "DD", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist)

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- log(start$UMSY[1]/(1 - start$UMSY[1]))
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
    params$logit_UMSY <- log(UMSY_start/(1 - UMSY_start))
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist, na.rm = TRUE)
    params$log_MSY <- log(3 * AvC)
  }
  if(is.null(params$log_q)) params$log_q <- log(Data@Mort[x])

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist)
  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   DLL = "MSEtool", silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_TMB",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearplusk <- c(Year, max(Year) + 1:k)

    Assessment <- new("Assessment", Model = "DD_TMB",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY,
                      SSBMSY = report$BMSY, VBMSY = report$BMSY,
                      B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$B0, VB0 = report$B0, h = report$h,
                      U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$B0, names = Yearplusone),
                      SSB = structure(report$B, names = Yearplusone),
                      SSB_SSBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                      VB = structure(report$B, names = Yearplusone),
                      VB_VBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      VB_VB0 = structure(report$B/report$B0, names = Yearplusone),
                      R = structure(report$R, names = Yearplusk),
                      N = structure(report$N, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Catch = structure(report$Cpred, names = Year),
                      NLL = structure(c(opt$objective, report$nll), names = c("Total", "Catch")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2],
                      SE_U_UMSY_final = SD$sd[5], SE_B_BMSY_final = SD$sd[6],
                      SE_B_B0_final = SD$sd[7], SE_SSB_SSBMSY_final = SD$sd[6],
                      SE_SSB_SSB0_final = SD$sd[7], SE_VB_VBMSY_final = SD$sd[6],
                      SE_VB_VB0_final = SD$sd[7], info = info, obj = obj, opt = opt,
                      SD = SD, TMB_report = report, dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(DD_TMB) <- "Assess"


#' @describeIn DD_TMB State-Space version of Delay-Difference model
#' @export
DD_SS <- function(x, Data, start = NULL, silent = TRUE, control = list(), ...) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@CV_Cat, Data@Ind"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x,   ]))[1]
  yind <- ystart:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  E_hist <- C_hist/I_hist
  if(any(is.na(E_hist))) stop("Missing values in catch and index in Data object.")
  E_hist <- E_hist/mean(E_hist)
  ny <- length(C_hist)
  k <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k[k > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho <- (wa[k + 2] - Winf)/(wa[k + 1] - Winf)
  Alpha <- Winf * (1 - Rho)
  S0 <- exp(-Data@Mort[x])  # get So survival rate
  wk <- wa[k]
  data <- list(model = "DD_SS", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist)

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- log(start$UMSY[1]/(1 - start$UMSY[1]))
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
    params$logit_UMSY <- log(UMSY_start/(1 - UMSY_start))
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist, na.rm = TRUE)
    params$log_MSY <- log(3 * AvC)
  }
  if(is.null(params$log_q)) params$log_q <- log(Data@Mort[x])
  if(is.null(params$log_tau)) params$log_tau <- log(0.3)

  sigmaC <- max(0.05, sdconv(1, Data@CV_Cat[x]))

  params$log_sigma = log(sigmaC)
  params$log_rec_dev = rep(0, ny - k)
  info <- list(Year = Year, data = data, params = params, sigma = sigmaC,
               I_hist = I_hist, control = control)

  obj <- MakeADFun(data = info$data, parameters = info$params, random = "log_rec_dev",
                   map = list(log_sigma = factor(NA)), checkParameterOrder = FALSE,
                   DLL = "MSEtool", silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_SS",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearplusk <- c(Year, max(Year) + 1:k)
    Yearrandom <- seq(Year[1] + k, max(Year))

    Assessment <- new("Assessment", Model = "DD_SS",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY,
                      SSBMSY = report$BMSY, VBMSY = report$BMSY,
                      B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$B0, VB0 = report$B0, h = report$h,
                      U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$B0, names = Yearplusone),
                      SSB = structure(report$B, names = Yearplusone),
                      SSB_SSBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                      VB = structure(report$B, names = Yearplusone),
                      VB_VBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      VB_VB0 = structure(report$B/report$B0, names = Yearplusone),
                      R = structure(report$R, names = Yearplusk),
                      N = structure(report$N, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Catch = structure(report$Cpred, names = Year),
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$jnll_comp), names = c("Total", "Catch", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[6],
                      SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      SE_SSB_SSBMSY_final = SD$sd[7], SE_SSB_SSB0_final = SD$sd[8],
                      SE_VB_VBMSY_final = SD$sd[7], SE_VB_VB0_final = SD$sd[8],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(DD_SS) <- "Assess"


