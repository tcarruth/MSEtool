#' Statistical catch-at-age (SCA) model
#'
#' A generic statistical catch-at-age model (single fleet, single season) that uses catch, index, and catch-at-age composition
#' data. An annual harvest rate is calculated (assuming a pulse fishery) as described in Forrest et al. (2008).
#' There are two parameterizations for estimation of recruitment deviations, the stock-recruit relationship,
#' and reference points (see functions section below).
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param vulnerability Whether estimated vulnerability is \code{"logistic"} or \code{"dome"} (double-normal).
#' See details for parameterization.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix. See details.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the DLMtool operating model),
#' vulnerable biomass (VB), or spawning stock biomass (SSB).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional list of starting values. See details.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the model for \code{SCA}. This only affects
#' calculation of reference points for \code{SCA2}.
#' @param fix_U_equilibrium Logical, whether the equilibrium harvest rate prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{U_equilibrium} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes virgin conditions).
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@sigmaR}.
#' @param common_dev Typically, a numeric for the number of most recent years in which a common recruitment deviation will
#' be estimated (in \code{SCA2}, uninformative years will have a recruitment closer to the mean, which can be very misleading,
#' especially near the end of the time series). By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @param early_dev Character string describing the years for which recruitment deviations are estimated in \code{SCA}. By default, \code{"comp_onegen"}
#' rec devs are estimated one full generation prior to the first year when catch-at-age (CAA) data are available. With \code{"comp"}, rec devs are
#' estimated starting in the first year with CAA. With \code{"all"}, rec devs start at the beginning of the model.
#' @param late_dev Typically, a numeric for the number of most recent years in which recruitment deviations will
#' not be estimated in \code{SCA} (recruitment in these years will be based on the mean predicted by stock-recruit relationship).
#' By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
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
#' For the statistical catch-at-age model, the basic data inputs are catch (by weight), index
#' (by weight/biomass), and catch-at-age matrix (by numbers). Catches are
#' assumed to be known perfectly (the harvest rate in a given year is the ratio of the observed
#' catch to the vulnerable biomass at the beginning of the year).
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
#' \item \code{R0} Virgin recruitment, only for \code{SCA}.
#' \item \code{h} Steepness, only for \code{SCA}. If not provided, the value in \code{Data@@steep} is used.
#' \item \code{meanR} Mean recruitment, only for \code{SCA2}.
#' \item \code{U_equilibrium}. Harvest rate prior to the first year of model, e.g. zero means virgin conditions.
#' \item \code{vul_par} (length 2 vector for logistic or length 4 for dome, see above).
#' \item \code{sigma} Standard deviation of index.
#' \item \code{tau} Standard deviation of recruitment deviations. If not provided, the value in \code{Data@@sigmaR} is used.
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
#' @examples
#' \dontrun{
#' res <- SCA(Data = DLMtool::SimulatedData)
#' res2 <- SCA2(Data = DLMtool::Simulation_1)
#' }
#' @describeIn SCA The parameterization with R0 and steepness as leading parameters.
#' @export
SCA <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"),
                CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1",
                start = NULL, fix_h = FALSE, fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                silent = TRUE, control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  dependencies = ""
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
  early_dev <- match.arg(early_dev)
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

  if(early_dev == "all") {
    est_early_rec_dev <- rep(1, max_age-1)
    est_rec_dev <- rep(1, n_y)
  }
  if(early_dev == "comp") {
    est_early_rec_dev <- rep(NA, max_age-1)
    ind1 <- which(!is.na(CAA_n_nominal))[1]
    est_rec_dev <- ifelse(c(1:n_y) < ind1, NA, 1)
  }
  if(early_dev == "comp_onegen") {
    ind1 <- which(!is.na(CAA_n_nominal))[1] - max_age
    if(ind1 < 0) {
      early_start <- max_age + ind1
      est_early_rec_dev <- ifelse(c(1:(max_age-1)) < early_start, NA, 1)
      est_rec_dev <- rep(1, n_y)
    } else {
      est_early_rec_dev <- rep(NA, max_age-1)
      est_rec_dev <- ifelse(c(1:n_y) < ind1, NA, 1)
    }
  }
  if(is.character(late_dev) && late_dev == "comp50") {
    CAA_all <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE))
    CAA_mode <- which.max(CAA_all)[1]
    comp50_ind <- which(CAA_all[1:CAA_mode] <= 0.5)
    comp50_ind <- comp50_ind[length(comp50_ind)]
    late_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  if(is.numeric(late_dev) && late_dev > 0) {
    if(late_dev > length(est_rec_dev)) late_dev <- length(est_rec_dev)
    ind_late <- (length(est_rec_dev) - late_dev + 1):length(est_rec_dev)
    est_rec_dev[ind_late] <- NA
  }

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA", C_hist = C_hist * rescale, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = vulnerability, I_type = I_type,
               SR_type = SR, est_early_rec_dev = est_early_rec_dev, est_rec_dev = est_rec_dev)

  # Starting values
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
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) params$vul_par <- start$vul_par
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau)
  }
  if(is.null(params$log_R0)) {
    params$log_R0 <- ifelse(is.null(Data@OM$N0[x]) | is.na(Data@OM$N0[x]),
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
  if(is.null(params$log_tau)) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Data@Year, data = data, params = params, LH = LH, control = control,
               inner.control = inner.control, rescale = rescale)

  map <- list()
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(any(is.na(est_early_rec_dev))) {
    n_est <- sum(!is.na(est_early_rec_dev))
    if(n_est == 0) map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
    else {
      est_early_rec_dev[!is.na(est_early_rec_dev)] <- 1:n_est
      map$log_early_rec_dev <- factor(est_early_rec_dev)
    }
  }
  if(any(is.na(est_rec_dev))) {
    n_est <- sum(!is.na(est_rec_dev))
    est_rec_dev[!is.na(est_rec_dev)] <- 1:n_est
    map$log_rec_dev <- factor(est_rec_dev)
  }

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev")

  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
  mod <- optimize_TMB_model(obj, control)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  if(rescale != 1) {
    vars_div <- c("B", "E", "CAApred", "CN", "N", "VB", "R", "R_early", "VB0", "R0", "B0", "E0", "N0")
    vars_mult <- "Brec"
    var_trans <- c("R0", "q")
    fun_trans <- c("/", "*")
    fun_fixed <- c("log", NA)
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }

  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  else {
    ref_pt <- get_MSY(Arec = report$Arec, Brec = report$Brec, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    report <- c(report, ref_pt)
    #get_MSY_opt(opt$par[1:2], M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    #covar <- SD$cov.fixed[1:2, 1:2]
    #grad_MSY <- jacobian(get_MSY_opt, x = opt$par[1:2], M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    #covar_MSY <- grad_MSY %*% covar %*% t(grad_MSY)
    Yearplusone <- c(Year, max(Year) + 1)
    YearEarly <- (Year[1] - max_age + 1):(Year[1] - 1)
    YearDev <- c(YearEarly, Year)
    YearR <- c(YearDev, max(YearDev) + 1)
    R <- c(rev(report$R_early), report$R)

    Dev <- c(rev(report$log_early_rec_dev), report$log_rec_dev)
    Dev_out <- structure(Dev, names = YearDev)

    if(integrate) {
      if(!all(is.na(est_early_rec_dev))) SE_Early <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_early_rec_dev"])
      SE_Main <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_rec_dev"])
    } else {
      if(!all(is.na(est_early_rec_dev))) SE_Early <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_early_rec_dev"])
      SE_Main <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])
    }

    SE_Early2 <- est_early_rec_dev
    if(!all(is.na(est_early_rec_dev))) {
      SE_Early2[!is.na(SE_Early2)] <- SE_Early
    }
    SE_Main2 <- est_rec_dev
    SE_Main2[!is.na(SE_Main2)] <- SE_Main

    SE_Dev <- structure(c(rev(SE_Early2), SE_Main2), names = YearDev)
    SE_Dev[is.na(SE_Dev)] <- 0

    first_non_zero <- which(Dev != 0)[1]
    if(first_non_zero > 1) {
      Dev_out <- Dev_out[-c(1:(first_non_zero - 1))]
      SE_Dev <- SE_Dev[-c(1:(first_non_zero - 1))]
    }

    Assessment <- new("Assessment", Model = "SCA", UMSY = report$UMSY,
                      MSY = report$MSY, BMSY = report$BMSY, SSBMSY = report$EMSY,
                      VBMSY = report$VBMSY, B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$E0, VB0 = report$VB0,
                      h = report$h, U = structure(report$U, names = Year),
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
                      R = structure(R, names = YearR),
                      N = structure(rowSums(report$N), names = Yearplusone),
                      N_at_age = report$N,
                      Selectivity = matrix(report$vul, nrow = length(Year),
                                           ncol = max_age, byrow = TRUE),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Obs_C_at_age = CAA_hist,
                      Index = structure(report$Ipred, names = Year),
                      C_at_age = report$CAApred,
                      Dev = Dev_out,
                      Dev_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$nll_comp),
                                      names = c("Total", "Index", "CAA", "Dev")),
                      #SE_UMSY = sqrt(covar_MSY[1,1]), SE_MSY = sqrt(covar_MSY[2,2]), #SE_U_UMSY_final = SD$sd[6],
                      #SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      #SE_SSB_SSBMSY_final = SD$sd[9], SE_SSB_SSB0_final = SD$sd[10],
                      #SE_VB_VBMSY_final = SD$sd[11], SE_VB_VB0_final = SD$sd[12],
                      SE_Dev = SE_Dev,
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SCA) <- "Assess"




get_MSY_opt <- function(x, M, weight, mat, vul, SR = c("BH", "Ricker"), vul_type = NULL) {
  SR <- match.arg(SR)
  maxage <- length(M)

  R0 <- exp(x[1])
  surv0 <- exp(-M)
  NPR0 <- c(1, cumprod(surv0[1:(maxage-1)]))
  NPR0[maxage] <- NPR0[maxage]/(1 - surv0[maxage])
  E0 <- R0 * sum(NPR0 * weight * mat)
  EPR0 <- E0/R0
  if(SR == "BH") {
    h <- 0.2 + 0.8/(1 + exp(-x[2]))
    Arec <- 4*h/(1-h)/EPR0
    Brec <- (5*h-1)/(1-h)/E0
  }
  if(SR == "Ricker") {
    h <- 0.2 + exp(x[2])
    Arec <- 1/EPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / E0
  }
  #vul generate vul parameters here

  solveMSY <- function(logit_U) {
    U <- ilogit(logit_U)
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
  UMSY <- ilogit(opt2$minimum)
  MSY <- -1 * as.numeric(opt2$objective)

  return(c(UMSY = UMSY, MSY = MSY))
}


vul_fn <- function(vul_par, maxage, type = c("logistic", "dome")) {
  age <- 1:maxage

  if(type == "logistic") {
    a50 <- exp(vul_par[1])
    a95 <- a50 + exp(vul_par[2])
    vul = 1/(1 + exp(-log(19) * (age - a50)/(a95 - a50)))
  }
  if(type == "dome") {
    sd_asc <- exp(vul_par[1])
    mu_asc <- exp(vul_par[2])
    mu_des <- mu_asc + exp(vul_par[3])
    sd_des <- exp(vul_par[4])

    denom_asc <- dnorm(mu_asc, mu_asc, sd_asc)
    denom_des <- dnorm(mu_des, mu_des, sd_des)

    vul <- rep(NA, maxage)
    for(i in age) {
      if(i <= mu_asc) {
        vul[i] <- dnorm(i, mu_asc, sd_asc)
      } else if(i <= mu_des) {
        vul[i] <- 1
      } else {
        vul[i] <- dnorm(i, mu_des, sd_des)
      }
    }
  }
  return(vul)
}
