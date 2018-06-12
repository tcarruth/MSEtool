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
#' @param Data An object of class \linkS4class{Data}.
#' @param SR Stock-recruit function (either \code{BH} for Beverton-Holt or \code{Ricker}).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional list of starting values. See details.
#' @param fix_U_equilibrium Logical, whether the equilibrium harvest rate prior to the first year of the model is
#' estimated. If TRUE, U_equilibruim is fixed to value provided in start (if provided), otherwise, equal to zero
#' (assumes virgin conditions).
#' @param fix_sigma Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Cat}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, equal to 1.
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a state-space variable).
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}} via \code{\link[TMB]{MakeADFun}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @details
#' To provide starting values for \code{DD_TMB}, a named list can be provided for \code{UMSY},
#' \code{MSY}, and \code{q} via the \code{start} argument (see example).
#'
#' For \code{DD_SS}, additional start values can be provided for and \code{sigma} and \code{tau}, the standard
#' deviation of the catch and recruitment variability, respectively.
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
#' #### Observation-error delay difference model
#' res <- DD_TMB(Data = sim_snapper)
#'
#'
#' # Provide starting values
#' start <- list(UMSY = 0.05, MSY = 4)
#' res <- DD_TMB(Data = sim_snapper, start = start)
#'
#' summary(res@@SD) # Look at parameter estimates
#'
#' \dontrun{
#' ### State-space version
#' res <- DD_SS(Data = sim_snapper)
#'
#' # Plot and save figures
#' plot(res)
#' }
#' @useDynLib MSEtool
#' @export
DD_TMB <- function(x = 1, Data, SR = c("BH", "Ricker"), rescale = "mean1", start = NULL,
                   fix_U_equilibrium = TRUE, silent = TRUE, control = list(iter.max = 1e6, eval.max = 1e6), ...) {
  SR <- match.arg(SR)
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

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "DD", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist * rescale, SR_type = SR)
  LH <- list(LAA = la, WAA = wa, maxage = Data@MaxAge, A50 = k)

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- logit(start$UMSY[1])
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equililbrium
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
    params$logit_UMSY <- logit(UMSY_start)
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist * rescale)
    params$log_MSY <- log(3 * AvC)
  }
  if(is.null(params$log_q)) params$log_q <- log(1)
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist, LH = LH, rescale = rescale, control = control)

  map <- list()
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)

  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, DLL = "MSEtool", silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_TMB",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {

    if(rescale != 1) {
      vars_div <- c("B0", "B", "Cpred", "BMSY", "MSY", "N0", "N", "R", "R0")
      vars_mult <- c("Brec")
      var_trans <- c("MSY")
      fun_trans <- c("/")
      fun_fixed <- c("log")
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    }
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
                      SE_UMSY = SD$sd[names(SD$value) == "UMSY"], SE_MSY = SD$sd[names(SD$value) == "MSY"],
                      SE_U_UMSY_final = SD$sd[names(SD$value) == "U_UMSY_final"],
                      SE_B_BMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_B_B0_final = SD$sd[names(SD$value) == "B_B0_final"],
                      SE_SSB_SSBMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_SSB_SSB0_final = SD$sd[names(SD$value) == "B_B0_final"],
                      SE_VB_VBMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_VB_VB0_final = SD$sd[names(SD$value) == "B_B0_final"], info = info, obj = obj, opt = opt,
                      SD = SD, TMB_report = report, dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(DD_TMB) <- "Assess"


#' @describeIn DD_TMB State-Space version of Delay-Difference model
#' @export
DD_SS <- function(x = 1, Data, SR = c("BH", "Ricker"), rescale = "mean1", start = NULL,
                  fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                  integrate = FALSE, silent = TRUE, control = list(iter.max = 1e6, eval.max = 1e6),
                  inner.control = list(), ...) {
  SR <- match.arg(SR)
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@CV_Cat, Data@Ind"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x, ]))[1]
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

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "DD_SS", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist * rescale, SR_type = SR)
  LH <- list(LAA = la, WAA = wa, maxage = Data@MaxAge, A50 = k)

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- logit(start$UMSY[1])
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equililbrium
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma[1])
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
    params$logit_UMSY <- logit(UMSY_start)
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist * rescale)
    params$log_MSY <- log(3 * AvC)
  }
  if(is.null(params$log_q)) params$log_q <- log(1)
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0
  if(is.null(params$log_sigma)) {
    sigmaC <- max(0.05, sdconv(1, Data@CV_Cat[x]))
    params$log_sigma <- log(sigmaC)
  }
  if(is.null(params$log_tau)) params$log_tau <- log(1)
  params$log_rec_dev = rep(0, ny - k)

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist, LH = LH,
               rescale = rescale, control = control, inner.control = inner.control)

  map <- list()
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)

  random <- NULL
  if(integrate) random <- "log_rec_dev"

  obj <- MakeADFun(data = info$data, parameters = info$params, random = random,
                   map = map, checkParameterOrder = FALSE,
                   DLL = "MSEtool", inner.control = inner.control, silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_SS",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {

    if(rescale != 1) {
      vars_div <- c("B0", "B", "Cpred", "BMSY", "MSY", "N0", "N", "R", "R0")
      vars_mult <- c("Brec")
      var_trans <- c("MSY")
      fun_trans <- c("/")
      fun_fixed <- c("log")
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    }

    Yearplusone <- c(Year, max(Year) + 1)
    Yearplusk <- c(Year, max(Year) + 1:k)
    YearDev <- seq(Year[1] + k, max(Year))

    if(integrate) {
      Dev <- SD$par.random
      SE_Dev <- sqrt(SD$diag.cov.random)
    } else {
      Dev <- SD$par.fixed[names(SD$par.fixed) == "log_rec_dev"]
      SE_Dev <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])
    }

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
                      Dev = structure(Dev, names = YearDev),
                      Dev_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$jnll_comp), names = c("Total", "Catch", "Dev")),
                      SE_UMSY = SD$sd[names(SD$value) == "UMSY"], SE_MSY = SD$sd[names(SD$value) == "MSY"],
                      SE_U_UMSY_final = SD$sd[names(SD$value) == "U_UMSY_final"],
                      SE_B_BMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_B_B0_final = SD$sd[names(SD$value) == "B_B0_final"],
                      SE_SSB_SSBMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_SSB_SSB0_final = SD$sd[names(SD$value) == "B_B0_final"],
                      SE_VB_VBMSY_final = SD$sd[names(SD$value) == "B_BMSY_final"],
                      SE_VB_VB0_final = SD$sd[names(SD$value) == "B_B0_final"],
                      SE_Dev = structure(SE_Dev, names = YearDev),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(DD_SS) <- "Assess"


