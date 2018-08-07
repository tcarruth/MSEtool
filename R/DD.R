#' Delay - Difference Stock Assessment in TMB
#'
#' A simple delay-difference assessment model using a
#' time-series of catches and a relative abundance index and coded in TMB. The model
#' is conditioned on effort and estimates predicted catch. In the state-space version,
#' recruitment deviations from the stock-recruit relationship are estimated.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 when running an assessment.
#' @param Data An object of class \linkS4class{Data}.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional list of starting values. See details.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the assessment model.
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
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}} via \code{\link[TMB]{MakeADFun}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @details
#' To provide starting values for \code{DD_TMB}, a named list can be provided for \code{R0} (virgin recruitment),
#' \code{h} (steepness), and \code{q} (catchability coefficient) via the \code{start} argument (see example).
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
#' @section Required Data:
#' \itemize{
#' \item \code{DD_TMB}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' \item \code{DD_SS}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' }
#' @section Optional Data:
#' \itemize{
#' \item \code{DD_TMB}: steep
#' \item \code{DD_SS}: steep, CV_Cat
#' }
#' @import TMB
#' @importFrom stats nlminb
#' @examples
#' \donttest{
#' #### Observation-error delay difference model
#' res <- DD_TMB(Data = DLMtool::Red_snapper)
#'
#' # Provide starting values
#' start <- list(R0 = 1, h = 0.95)
#' res <- DD_TMB(Data = DLMtool::Red_snapper, start = start)
#'
#' summary(res@@SD) # Parameter estimates
#'
#' ### State-space version
#' ### Set recruitment variability SD = 0.3 (since fix_tau = TRUE)
#' res <- DD_SS(Data = Red_snapper, start = list(tau = 0.3))
#' }
#' @seealso \link{plot,Assessment,ANY-method} \link{summary,Assessment-method} \link{make_MP}
#' @useDynLib MSEtool
#' @export
DD_TMB <- function(x = 1, Data, SR = c("BH", "Ricker"), rescale = "mean1", start = NULL, fix_h = FALSE,
                   fix_U_equilibrium = TRUE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                   control = list(iter.max = 5e3, eval.max = 1e4), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  SR <- match.arg(SR)
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x])
  a50V <- max(a50V, 1)
  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x,   ]))[1]
    yind <- ystart:length(Data@Cat[x, ])
  }
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
    if(!is.null(start$R0) && is.numeric(start$R0)) params$log_R0 <- log(start$R0[1] * rescale)
    if(!is.null(start$h) && is.numeric(start$h)) {
      if(SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equililbrium
  }
  if(is.null(params$log_R0)) {
    params$log_R0 <- ifelse(is.null(Data@OM$N0[x]),
                            log(mean(data$C_hist)) + 4, 1.1 * rescale * Data@OM$N0[x] * (1 - exp(-Data@Mort[x])))
  }
  if(is.null(params$transformed_h)) {
    h_start <- ifelse(is.na(Data@steep[x]), 0.9, Data@steep[x])
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$log_q)) params$log_q <- log(1)
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist, LH = LH, rescale = rescale, control = control)

  map <- list()
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)

  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, DLL = "MSEtool", silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  if(rescale != 1) {
    vars_div <- c("B0", "B", "Cpred", "N0", "N", "R0", "R")
    vars_mult <- c("Brec")
    var_trans <- c("R0")
    fun_trans <- c("/")
    fun_fixed <- c("log")
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }
  Yearplusone <- c(Year, max(Year) + 1)
  Yearplusk <- c(Year, max(Year) + 1:k)

  nll_report <- ifelse(is.character(opt), report$nll, opt$objective)
  Assessment <- new("Assessment", Model = "DD_TMB", Name = Data@Name, conv = !is.character(SD),
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$B0, VB0 = report$B0, h = report$h,
                    U = structure(report$U, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$B, names = Yearplusone),
                    SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                    VB = structure(report$B, names = Yearplusone),
                    VB_VB0 = structure(report$B/report$B0, names = Yearplusone),
                    R = structure(report$R, names = Yearplusk),
                    N = structure(report$N, names = Yearplusone),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Cpred/(C_hist/I_hist), names = Year),
                    NLL = structure(c(nll_report, report$nll - report$penalty, report$penalty),
                                    names = c("Total", "Catch", "Penalty")),
                    info = info, obj = obj, opt = opt,
                    SD = SD, TMB_report = report, dependencies = dependencies)

  if(!is.character(opt) && !is.character(SD)) {
    ref_pt <- get_MSY_DD(info$data, report$Arec, report$Brec)
    report <- c(report, ref_pt)

    Assessment@UMSY <- report$UMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- Assessment@SSBMSY <- Assessment@VBMSY <- report$BMSY
    Assessment@U_UMSY <- structure(report$U/report$UMSY, names = Year)
    Assessment@B_BMSY <- Assessment@SSB_SSBMSY <- Assessment@VB_VBMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(DD_TMB) <- "Assess"


#' @describeIn DD_TMB State-Space version of Delay-Difference model
#' @useDynLib MSEtool
#' @export
DD_SS <- function(x = 1, Data, SR = c("BH", "Ricker"), rescale = "mean1", start = NULL,
                  fix_h = FALSE, fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                  integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                  control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  SR <- match.arg(SR)
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x, ]))[1]
    yind <- ystart:length(Data@Cat[x, ])
  }
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
    if(!is.null(start$R0) && is.numeric(start$R0)) params$log_R0 <- log(start$R0[1] * rescale)
    if(!is.null(start$h) && is.numeric(start$h)) {
      if(SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if(!is.null(start$q) && is.numeric(start$q)) params$log_q <- log(start$q[1])
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equililbrium
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma[1])
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  if(is.null(params$log_R0)) {
    params$log_R0 <- ifelse(is.null(Data@OM$N0[x]),
                            log(mean(data$C_hist)) + 4, 1.1 * rescale * Data@OM$N0[x] * (1 - exp(-Data@Mort[x])))
  }
  if(is.null(params$transformed_h)) {
    h_start <- ifelse(is.na(Data@steep[x]), 0.9, Data@steep[x])
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$log_q)) params$log_q <- log(1)
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0
  if(is.null(params$log_sigma)) {
    sigmaC <- max(0.05, sdconv(1, Data@CV_Cat[x]))
    params$log_sigma <- log(sigmaC)
  }
  if(is.null(params$log_tau)) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  params$log_rec_dev = rep(0, ny - k)

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist, LH = LH,
               rescale = rescale, control = control, inner.control = inner.control)

  map <- list()
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)

  random <- NULL
  if(integrate) random <- "log_rec_dev"

  obj <- MakeADFun(data = info$data, parameters = info$params, random = random,
                   map = map, checkParameterOrder = FALSE,
                   DLL = "MSEtool", inner.control = inner.control, silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  if(rescale != 1) {
    vars_div <- c("B0", "B", "Cpred", "N0", "N", "R0", "R")
    vars_mult <- c("Brec")
    var_trans <- c("R0")
    fun_trans <- c("/")
    fun_fixed <- c("log")
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }

  Yearplusone <- c(Year, max(Year) + 1)
  Yearplusk <- c(Year, max(Year) + 1:k)
  YearDev <- seq(Year[1] + k, max(Year))

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "DD_SS", Name = Data@Name, conv = !is.character(SD),
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$B0, VB0 = report$B0, h = report$h,
                    U = structure(report$U, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$B, names = Yearplusone),
                    SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                    VB = structure(report$B, names = Yearplusone),
                    VB_VB0 = structure(report$B/report$B0, names = Yearplusone),
                    R = structure(report$R, names = Yearplusk),
                    N = structure(report$N, names = Yearplusone),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Cpred/(C_hist/I_hist), names = Year),
                    Dev = structure(report$log_rec_dev, names = YearDev),
                    Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty),
                                    names = c("Total", "Catch", "Dev", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(!is.character(opt) && !is.character(SD)) {
    ref_pt <- get_MSY_DD(info$data, report$Arec, report$Brec)
    report <- c(report, ref_pt)

    if(integrate) {
      SE_Dev <- sqrt(SD$diag.cov.random)
    } else {
      SE_Dev <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])
    }

    Assessment@UMSY <- report$UMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- Assessment@SSBMSY <- Assessment@VBMSY <- report$BMSY
    Assessment@U_UMSY <- structure(report$U/report$UMSY, names = Year)
    Assessment@B_BMSY <- Assessment@SSB_SSBMSY <- Assessment@VB_VBMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@SE_Dev <- structure(SE_Dev, names = YearDev)
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(DD_SS) <- "Assess"


get_MSY_DD <- function(TMB_data, Arec, Brec) {
  S0 <- TMB_data$S0
  Alpha <- TMB_data$Alpha
  Rho <- TMB_data$Rho
  wk <- TMB_data$wk
  SR <- TMB_data$SR_type

  solveMSY <- function(x) {
    U <- ilogit(x)
    SS <- S0 * (1 - U)
    Spr <- (SS * Alpha/(1 - SS) + wk)/(1 - Rho * SS)
    if(SR == "BH") Req <- (Arec * Spr * (1 - U) - 1)/(Brec * Spr * (1 - U))
    if(SR == "Ricker") Req <- log(Arec * Spr * (1 - U))/(Brec * Spr * (1 - U))
    Beq <- Spr * Req
    Yield <- U * Beq
    return(-1 * Yield)
  }

  opt2 <- optimize(solveMSY, interval = c(-50, 6))
  UMSY <- ilogit(opt2$minimum)
  MSY <- -1 * opt2$objective
  BMSY <- MSY/UMSY
  return(list(UMSY = UMSY, MSY = MSY, BMSY = BMSY))
}
