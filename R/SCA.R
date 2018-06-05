#' Statistical catch-at-age (SCA) model
#'
#' A generic statistical catch-at-age model (single fleet, single season) that uses catch, index, and catch-at-age composition
#' data. An annual harvest rate is calculated (assuming a pulse fishery) as described in Forrest et al. (2008).
#' There are two parameterizations for estimation of recruitment deviations, the stock-recruit relationship,
#' and reference points (see functions section below).
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param SR Stock-recruit function (either \code{BH} for Beverton-Holt or \code{Ricker}).
#' @param vulnerability Whether estimated vulnerability is logistic or dome (double-normal).
#' See details for parameterization.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix. See details.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the DLMtool operating model),
#' vulnerable biomass (VB), or spawning stock biomass (SSB).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional list of starting values. See details.
#' @param fix_U_equilibrium Logical, whether the equilibrium harvest rate prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{U_equilibruim} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes virgin conditions).
#' @param fix_sigma Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, equal to 1.
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a state-space variable).
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param control A named list of agruments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Other arguments to be passed.
#' @return An object of class \linkS4class{Assessment}.
#' @details
#' For the SCA model, the basic data inputs are catch (by weight), index
#' (by weight/biomass), and catch-at-age matrix (by numbers). Catches are
#' assumed to be known perfectly (the harvest rate in a given year is the ratio of the observed
#' catch to the vulnerable biomass at the beginning of the year), while the observation error
#' of the index is fixed based on the value of \code{CV_Ind} in the Data object.
#'
#' The annual sample sizes of the catch-at-age matrix is provided to the model (used in the
#' likelihood for catch-at-age, assuming a multinomial distribution),
#' and is manipulated via argument \code{CAA_multiplier}. This argument is
#' interpreted in two different ways depending on the value provided.
#' If \code{CAA_multiplier > 1}, then this value will cap the annual sample sizes
#' to that number. If \code{CAA_multiplier <= 1}, then all the annual samples sizes
#' will be re-scaled by that number. By default, sample sizes are capped at 50.
#'
#' Vulnerability can be specified to be either logistic or dome. If logistic, then the model
#' vector \code{vul_par} is of length 2, containing the ages of 50\% and 95\% vulnerability,
#' respectively. The age of 95\% vulnerability is an offset, i.e.,
#' \code{vul_95 = vul_par[1] + exp(vul_par[2])}.
#'
#' With dome vulnerability, a double normal parameterization is used, where \code{vul_par}
#' is an estimated vector of length 4:
#' \itemize{
#' \item \code{vulpar[1]}: \code{log(sd_asc)}, where sd_asc of the normal distribution function for the ascending limb
#' \item \code{vulpar[2]}: \code{mu_asc}, mean of the normal distribution function for the ascending limb
#' \item \code{vulpar[3]}: \code{mu_des}, mean of the normal distribution function for the descending limb.
#' This is parameterized as an offset, i.e., \code{mu_desc = mu_asc + exp(vulpar[3])} to ensure
#' \code{mu_desc > mu_asc}.
#' \item \code{vulpar[4]}: \code{log(sd_des)}, where sd_des of the normal distribution function for the descending limb.
#' }
#'
#' For \code{start}, a named list of starting values of estimates for:
#' \itemize{
#' \item \code{UMSY} Only for \code{SCA2}.
#' \item \code{MSY} Only for \code{SCA2}.
#' \item \code{meanR} Mean recruitment, only for \code{SCA}.
#' \item \code{U_equilibrium}. Harvest rate prior to the first year of model, e.g. zero means virgin conditions.
#' \item \code{vul_par} (length 2 vector for logistic or length 4 for dome, see above).
#' \item \code{sigma} Standard deviation of index.
#' \item \code{tau} Standard deviation of recruitment deviations.
#' }
#' @author Q. Huynh
#' @references
#' Cadigan, N.G. 2016. A state-space stock assessment model for northern cod, including under-reported catches and
#' variable natural mortality rates. Canadian Journal of Fisheries and Aquatic Science 72:296-308.
#'
#' Forrest, R.E., Martell, S.J.D., Melnychuk, M.C., and Walters, C.J. 2008.
#' An age-structured model with leading management parameters, incorporating
#' age-specific selectivity and maturity. Canadian Journal of Fisheries and Aquatic
#' Science 65:286-296.
#' @describeIn SCA The mean recruitment in the time series is estimated and recruitment deviations around this mean are estimated
#' (as penalized parameters similar to Cadigan 2016). This generally runs quickly. MSY reference points are estimated after the
#' assessment run.
#' @examples
#' res <- SCA(Data = DLMtool::SimulatedData)
#'
#' res <- SCA2(Data = DLMtool::Simulation_1)
#' @export
SCA <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"),
                CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1",
                start = NULL, fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                integrate = FALSE, silent = TRUE, control = list(iter.max = 1e6, eval.max = 1e6),
                inner.control = list(), ...) {
  dependencies = "Data@Cat, Data@Ind, Data@CAA, Data@Mort, Data@wla, Data@wlb, Data@vbLinf, Data@vbK, Data@vbt0, Data@L50, Data@L95, Data@MaxAge"
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
  yind <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  I_hist <- Data@Ind[x, yind]
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
  CAA_hist <- Data@CAA[x, yind, ]

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)

  n_y <- length(C_hist)
  max_age <- Data@MaxAge
  M <- rep(Data@Mort[x], max_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(1:max_age) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))
  LH <- list(LAA = La, WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA", C_hist = C_hist * rescale, I_hist = I_hist, CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M, weight = Wa, mat = mat_age,
               vul_type = vulnerability, I_type = I_type, est_rec_dev = rep(1L, length(CAA_n_nominal)))

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$meanR) && is.numeric(start$meanR)) params$log_meanR <- log(start$meanR)
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) params$vul_par <- start$vul_par
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau)
  }
  if(is.null(params$log_meanR)) params$log_meanR <- log(mean(C_hist * rescale)) + 1
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0
  if(is.null(params$vul_par)) {
    CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
    if(vulnerability == "logistic") {
      params$vul_par <- c(CAA_mode-1, log(1)) # 50 and log(95%-offset) vulnerability respectively
    }
    if(vulnerability == "dome") {
      params$vul_par <- c(log(1), CAA_mode, log(0.5), log(5)) # double normal: logsd(ascending), mean(asc), mean(desc)-logoffset, sd(desc)
    }
  }
  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.05, sdconv(1, Data@CV_Ind[x]))
    params$log_sigma <- log(sigmaI)
  }
  if(is.null(params$log_tau)) params$log_tau <- log(1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Data@Year, data = data, params = params, LH = LH, SR = SR, control = control,
               inner.control = inner.control, rescale = rescale)

  map <- list()
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)

  random <- NULL
  if(integrate) random <- "log_rec_dev"

  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)
  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {
    if(fix_U_equilibrium && params$U_equilibrium == 0) {
      SSB0 <- report$E[1]
      R0 <- report$R[1]
    } else SSB0 <- R0 <- NULL
    refpt <- get_refpt(SSB = report$E[1:(length(report$E) - 1)], rec = report$R[2:length(report$R)],
                       SSB0 = SSB0, R0 = R0, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    report <- c(report, refpt)
    if(rescale != 1) {
      vars_div <- c("meanR", "B", "E", "CAApred", "CN", "N", "VB", "R", "MSY", "VBMSY",
                    "RMSY", "BMSY", "EMSY", "VB0", "R0", "B0", "E0", "N0")
      vars_mult <- "Brec"
      var_trans <- c("meanR", "q")
      fun_trans <- c("/", "*")
      fun_fixed <- c("log", NA)
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    }

    Yearplusone <- c(Year, max(Year) + 1)
    YearDev <- Year

    if(integrate) {
      Dev <- SD$par.random
      SE_Dev <- sqrt(SD$diag.cov.random)
    } else {
      Dev <- SD$par.fixed[names(SD$par.fixed) == "log_rec_dev"]
      SE_Dev <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])
    }
    Assessment <- new("Assessment", Model = "SCA", UMSY = report$UMSY,
                      MSY = report$MSY, BMSY = report$BMSY, SSBMSY = report$EMSY,
                      VBMSY = report$VBMSY, B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$E0, VB0 = report$VB0, h = report$h,
                      U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$B0, names = Yearplusone),
                      SSB = structure(report$E, names = Yearplusone),
                      SSB_SSBMSY = structure(report$E/report$EMSY, names = Yearplusone),
                      SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                      VB = structure(report$VB, names = Yearplusone),
                      VB_VBMSY = structure(report$VB/report$VBMSY, names = Yearplusone),
                      VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                      R = structure(report$R, names = Yearplusone),
                      N = structure(rowSums(report$N), names = Yearplusone),
                      N_at_age = report$N,
                      Selectivity = matrix(report$vul, nrow = length(Year),
                                           ncol = max_age, byrow = TRUE),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Obs_C_at_age = CAA_hist,
                      Index = structure(report$Ipred, names = Year),
                      C_at_age = report$CAApred,
                      Dev = structure(Dev, names = YearDev),
                      Dev_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$nll_comp),
                                      names = c("Total", "Index", "CAA", "Dev")),
                      SE_Dev = structure(SE_Dev, names = YearDev),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SCA) <- "Assess"




get_refpt <- function(SSB, rec, SSB0 = NULL, R0 = NULL, M, weight, mat, vul, SR = c("BH", "Ricker")) {
  SR <- match.arg(SR)
  maxage <- length(M)

  solve_SR_par <- function(x, type = c("ab", "steepness")) {
    if(type == "steepness") {
      transformed_h <- x
      if(SR == "BH") {
        h <- 0.2 + 0.8 * 1/(1 + exp(-transformed_h))
        recpred <- ((0.8 * R0 * h * SSB)/(0.2 * SSBpR*R0*(1-h)+(h-0.2)*SSB))
      }
      if(SR == "Ricker") {
        h <- 0.2 + exp(transformed_h)
        recpred <- SSB/SSBpR * (5*h)^(1.25 * (1 - SSB/SSB0))
      }
    }

    if(type == "ab") {
      Arec <- exp(x[1])
      Brec <- exp(x[2])
      if(SR == "BH") recpred <- Arec * SSB / (1 + Brec * SSB)
      if(SR == "Ricker") recpred <- Arec * SSB * exp(-Brec * SSB)
    }
    sigmaR <- sqrt(sum((log(rec/recpred))^2)/length(recpred))
    nLL <- -sum(dnorm(log(rec/recpred), 0, sigmaR, log = TRUE))
    return(nLL)
  }

  if(is.null(SSB0) || is.null(R0)) {
    opt <- nlminb(log(c(10, 10)), solve_SR_par, type = "ab")
    Arec <- exp(opt$par[1])
    Brec <- exp(opt$par[2])
  } else {
    SSBpR <- SSB0/R0
    opt <- optimize(solve_SR_par, interval = c(-6, 6), type = "steepness")$minimum # steepness
    if(SR == "BH") {
      h <- 0.2 + 0.8/(1 + exp(-opt))
      Arec <- 4*h/(1-h)/SSBpR
      Brec <- (5*h-1)/(1-h)/SSB0
    }
    if(SR == "Ricker") {
      h <- 0.2 + exp(opt)
      Arec <- 1/SSBpR * (5*h)^1.25
      Brec <- 1.25 * log(5*h) / SSB0
    }
  }

  # virgin reference points
  surv0 <- exp(-M)
  NPR0 <- c(1, cumprod(surv0[1:(maxage-1)]))
  EPR0 <- sum(NPR0 * mat * weight)
  if(SR == "BH") {
    R0 <- (Arec * EPR0 - 1)/(Brec * EPR0)
    if(is.null(SSB0) || is.null(R0)) {
      CR <- Arec * EPR0
      h <- CR/(4 + CR)
    }
  }
  if(SR == "Ricker") {
    R0 <- log(Arec * EPR0)/(Brec * EPR0)
    if(is.null(SSB0) || is.null(R0)) h <- exp(0.8 * Brec * EPR0)
  }
  N0 <- R0 * sum(NPR0)
  E0 <- R0 * EPR0
  VB0 <- R0 * sum(NPR0 * weight * vul)
  B0 <- R0 * sum(NPR0 * weight)

  solveMSY <- function(logit_U) {
    U <- 1/(1 + exp(-logit_U))
    surv <- exp(-M) * (1 - vul * U)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <- sum(NPR * mat * weight)
    if(SR == "BH") Req <- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * U * NPR
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  opt2 <- optimize(solveMSY, interval = c(-6, 6))
  UMSY <- 1/(1 + exp(-opt2$minimum))
  MSY <- -1 * opt2$objective
  VBMSY <- MSY/UMSY

  surv_UMSY <- exp(-M) * (1 - vul * UMSY)
  NPR_UMSY <- c(1, cumprod(surv_UMSY[1:(maxage-1)]))
  NPR_UMSY[maxage] <- NPR_UMSY[maxage]/(1 - surv_UMSY[maxage])

  RMSY <- VBMSY/sum(vul * NPR_UMSY * weight)
  BMSY <- RMSY * sum(NPR_UMSY * weight)
  EMSY <- RMSY * sum(NPR_UMSY * weight * mat)
  return(list(h = h, Arec = Arec, Brec = Brec, UMSY = UMSY, MSY = MSY, VBMSY = VBMSY,
              RMSY = RMSY, BMSY = BMSY, EMSY = EMSY, VB0 = VB0, R0 = R0, B0 = B0, E0 = E0, N0 = N0))
}



SCA_jacobian <- function(obj, par_fixed, nll_comp, ...) {

  f <- function(x) {
    obj$report(x)[[nll_comp]]
  }

  jacobian(f, par_fixed, ...)
}
