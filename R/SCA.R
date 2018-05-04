#' Statistical catch-at-age (SCA) model
#'
#' A generic statistical catch-at-age model that uses catch, index, and catch-at-age composition
#' data. An annual harvest rate is calculated (assuming a pulse fishery) in lieu of an instantaneous
#' fishing mortality rate from continuous fishing as described in Forrest et al. (2008). Recruitment
#' deviations are estimated beginning in the year when age composition data are available. There are
#' several parameterizations for estimation of recruitment deviations, the stock-recruit relationship,
#' and reference points (see functions section below).
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data.
#' @param U_begin Indicates how to handle the equilibrium harvest rate prior to the first
#' year of model. \code{"virgin"} fixes \code{U = 0} prior to year 1 while \code{"est"} will
#' estimate it.
#' @param vulnerability Whether estimated vulnerability is logistic or dome (double-normal).
#' See details for parameterization.
#' @param SR Stock-recruit function (currently only supports Beverton-Holt).
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix. See details.
#' @param start Optional list of starting values. See details.
#' @param tau The standard deviation of the recruitment deviations from the estimated mean recruitment
#' or stock-recruit relationship (by default, euqal to one).
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood of
#' the recruitment deviations (thus, treating it as a state-space variable).
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param control A named list of agruments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Other arguments to be passed.
#' @return An object of class \linkS4class{Assessment}.
#' @note Convergence problems may occur if perfect data are used (e.g., \code{Data@@CV_Ind}
#' is very low).
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
#' Recruitment deviations are estimated beginning in the year when age composition are available.
#'
#' Vulnerability can be specified to be either logistic or dome. If logistic, then the model
#' vector \code{vul_par} is of length 2, containing the ages of 50\% and 95\% vulnerability,
#' respectively. The age of 95\% vulnerability is an offset, i.e.,
#' \code{vul_95 = vul_par[1] + exp(vul_par[2])}.
#'
#' With dome vulnerability, a double normal parameterization is used, where \code{vul_par}
#' is an estimated vector of length 4:
#' \itemize{
#' \item \code{vulpar[1]} \code{log(sd_asc)}, where sd_asc of the normal distribution function for the ascending limb
#' \item \code{vulpar[2]} \code{mu_asc}, mean of the normal distribution function for the ascending limb
#' \item \code{vulpar[3]} \code{mu_des}, mean of the normal distribution function for the descending limb.
#' This is parameterized as an offset, i.e., \code{mu_desc = mu_asc + exp(vulpar[3])} to ensure
#' \code{mu_desc > mu_asc}.
#' \item \code{vulpar[4]} \code{log(sd_des)}, where sd_des of the normal distribution function for the descending limb.
#' }
#'
#' For \code{start}, a named list of starting values of estimates for:
#' \itemize{
#' \item \code{UMSY} Only for \code{SCA2}.
#' \item \code{MSY} Only for \code{SCA2}.
#' \item \code{meanR} Mean recruitment, only for \code{SCA}.
#' \item \code{U_equilibrium} (typically 0). See \code{U_begin} argument.
#' \item \code{vul_par} (length 2 vector for logistic or length 4 for dome, see above).
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
#' assessment run, assuming that the recruitment in the first year of the model is the virgin recruitment.
#' @export
SCA <- function(x = 1, Data, U_begin = c("virgin", "est"), vulnerability = c("logistic", "dome"),
                SR = c("BH", "Ricker"), CAA_multiplier = 50, start = NULL,
                tau = 1, integrate = FALSE, silent = TRUE, control = list(eval.max = 1e3), inner.control = list(), ...) {
  dependencies = ""
  U_begin <- match.arg(U_begin)
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
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

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$meanR) && is.numeric(start$meanR)) params$log_meanR <- log(start$meanR)
	  if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
	  if(!is.null(start$vulnerability) && is.numeric(start$vulnerability)) params$vul_par <- start$vul_par
  }
  if(is.null(params$log_meanR)) params$log_meanR <- 2e3
	if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0
  if(is.null(params$vul_par)) {
    CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
    if(vulnerability == "logistic") {
      params$vul_par <- c(CAA_mode-1, log(1)) # 50 and log(95%-offset) vulnerability respectively
    }
    if(vulnerability == "dome") {
      # double normal: logsd(ascending), mean(asc), mean(desc)-logoffset, sd(desc)
      params$vul_par <- c(log(1), CAA_mode, log(0.5), log(5))
    }
  }
  sigma_I <- sdconv(1, max(0.05, Data@CV_Ind[x])) # can't fit to perfect data
  params$log_sigma <- log(sigma_I)
  params$log_tau <- log(tau)
  params$log_rec_dev <- rep(0, n_y)

  data <- list(model = "SCA", C_hist = C_hist, I_hist = I_hist, CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M, weight = Wa, mat = mat_age,
               vul_type = vulnerability, est_rec_dev = rep(1L, length(CAA_n_nominal)))
  info <- list(Year = Data@Year, data = data, params = params, control = control)

  map <- list(log_tau = factor(NA), log_sigma = factor(NA))
  if(U_begin == "virgin") map$U_equilibrium <- factor(NA)
  if(any(is.na(CAA_n_nominal))) {
    map$log_rec_dev <- random_map(CAA_n_nominal)
  }
  random <- NULL
  if(integrate) random <- "log_rec_dev"

  browser()

  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)
  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  else {
    refpt <- get_refpt(SSB = report$E[1:(length(report$E) - 1)], rec = report$R[2:length(report$R)],
                       SSB0 = report$E[1], R0 = report$R[1], M = M, weight = Wa, mat = mat_age, vul = report$vul,
                       SR = SR)
    report <- c(report, refpt)
    report$B0 <- report$B[1]
    report$R0 <- report$R[1]
    report$N0 <- report$N[1]
    report$E0 <- report$E[1]
    report$VB0 <- report$VB[1]

    Yearplusone <- c(Year, max(Year) + 1)
    Yearrandom <- seq(Year[2], max(Year))
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
                      Random = ifelse(integrate, structure(report$log_rec_dev, names = Yearrandom), numeric(0)),
                      Random_type = ifelse(integrate, "log-Recruitment deviations", character(0)),
                      NLL = structure(c(opt$objective, report$nll_comp),
                                      names = c("Total", "Index", "CAA", "log-Recruitment deviations")),
                      #SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[6],
                      #SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      #SE_SSB_SSBMSY_final = SD$sd[9], SE_SSB_SSB0_final = SD$sd[10],
                      #SE_VB_VBMSY_final = SD$sd[11], SE_VB_VB0_final = SD$sd[12],
                      #SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SCA) <- "Assess"
environment(SCA) <- asNamespace("MSEtool")




get_refpt <- function(SSB, rec, SSB0, R0, M, weight, mat, vul, SR = c("BH", "Ricker")) {
  SR <- match.arg(SR)
  SSBpR <- SSB0/R0

  solve_steep <- function(transformed_h) {
    if(SR == "BH") {
      h <- 0.2 + 0.8 * 1/(1 + exp(-transformed_h))
      recpred <- ((0.8 * R0 * h * SSB)/(0.2 * SSBpR*R0*(1-h)+(h-0.2)*SSB))
    }
    if(SR == "Ricker") {
      h <- 0.2 + exp(transformed_h)
      recpred <- SSB/SSBpR * (5*h)^(1.25 * (1 - SSB/SSB0))
    }
    sigmaR <- sqrt(sum((log(rec/recpred))^2)/length(recpred))
    nLL <- -sum(dnorm(log(rec/recpred), -0.5 * sigmaR^2, sigmaR, log = TRUE))
    return(nLL)
  }

  opt <- optimize(solve_steep, interval = c(-6, 6))$minimum
  if(SR == "BH") {
    h <- 0.2 + 1/(1+exp(-opt)) * 0.8
    Arec <- 4*h/(1-h)/SSBpR
    Brec <- (5*h-1)/(1-h)/SSB0
  }
  if(SR == "Ricker") {
    h <- 0.2 + exp(opt)
    Arec <- 1/SSBpR * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / SSB0
  }

  maxage <- length(M)

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
              RMSY = RMSY, BMSY = BMSY, EMSY = EMSY))
}





#SRpar <- function(SSB, rec, SSB0 = NA, R0 = NA, plot = FALSE, param = c("ab", "steepness"), SR = c("BH", "Ricker")) {
#  param <- match.arg(param)
#  SR <- match.arg(SR)
#
#  SSBpR <- SSB0/R0
#  solve_steep <- function(semilogit_h) {
#    h <- 0.2 + 0.8 * 1/(1 + exp(-semilogit_h))
#    recpred <- ((0.8 * R0 * h * SSB)/(0.2 * SSBpR*R0*(1-h)+(h-0.2)*SSB))
#    sigmaR <- sqrt(sum((log(rec/recpred))^2)/length(recpred))
#    nLL <- -sum(dnorm(log(rec/recpred), 0, sigmaR, log = TRUE))
#    return(nLL)
#  }
#  solve_ab <- function(pars) {
#    a <- exp(pars[1])
#    b <- exp(pars[2])
#    recpred <- a * SSB /(1 + b * SSB)
#    sigmaR <- sqrt(sum((log(rec/recpred))^2)/length(recpred))
#    nLL <- -sum(dnorm(log(rec/recpred), 0, sigmaR, log = TRUE))
#    return(nLL)
#  }
#
#  if(param == "steepness") {
#    opt <- optimize(solve_steep, interval = c(-6, 6))$minimum
#    if(plot) getBH(c(opt, log(R0)), SSB, rec, SSBpR, mode=2, plot=plot)
#    h <- 0.2 + 1/(1+exp(-opt)) * 0.8
#    return(h)
#  }
#  if(param == "ab") {
#    opt <- optim(c(0, 0), solve_ab, method = "BFGS")$par
#    return(exp(opt))
#  }
#}


