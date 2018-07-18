#' @describeIn SCA The mean recruitment in the time series is estimated and recruitment deviations around this mean are estimated
#' as penalized parameters (similar to Cadigan 2016). This version is generally very fast and robust. Virgin and MSY reference points
#' are estimated after the assessment run.
#' @export
SCA2 <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"),
                 CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1",
                 start = NULL, fix_h = FALSE, fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                 common_dev = "comp50", integrate = FALSE, silent = TRUE, control = list(iter.max = 5e3, eval.max = 1e4),
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
  est_rec_dev <- rep(1L, length(CAA_n_nominal))
  est_early_rec_dev <- rep(1L, max_age - 1)

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA2", C_hist = C_hist * rescale, I_hist = I_hist, CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M, weight = Wa, mat = mat_age,
               vul_type = vulnerability, I_type = I_type, est_early_rec_dev = est_early_rec_dev,
               est_rec_dev = est_rec_dev)

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$meanR) && is.numeric(start$meanR)) params$log_meanR <- log(start$meanR)
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) params$vul_par <- start$vul_par
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau)
  }
  if(is.null(params$log_meanR)) params$log_meanR <- log(mean(C_hist * rescale)) + 2
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
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Data@Year, data = data, params = params, LH = LH, SR = SR, control = control,
               inner.control = inner.control, rescale = rescale)

  map <- list()
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)

  if(is.character(common_dev) && common_dev == "comp50") {
    CAA_all <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE)) #CAA with max = 1
    CAA_mode <- which.max(CAA_all)[1] # Find max
    comp50_ind <- which(CAA_all[1:CAA_mode] <= 0.5)[1]
    comp50_ind <- comp50_ind[length(comp50_ind)]
    common_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  map_log_rec_dev <- 1:length(params$log_rec_dev)
  if(is.numeric(common_dev) && !is.na(common_dev) && common_dev > 0) {
    ind <- (length(map_log_rec_dev) - common_dev + 1):length(map_log_rec_dev)
    map_log_rec_dev[ind] <- map_log_rec_dev[ind[1]-1]
    map$log_rec_dev <- factor(map_log_rec_dev)
  }

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev")

  obj <- MakeADFun(data = info$data, parameters = info$params, checkParameterOrder = FALSE,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
  mod <- optimize_TMB_model(obj, control)
  opt <- mod[[1]]
  SD <- mod[[2]]

  report <- obj$report(obj$env$last.par.best)
  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA2", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  } else {
    info$h <- ifelse(fix_h, Data@steep[x], NA)
    refpt <- get_refpt2(SSB = report$E[1:(length(report$E) - 1)], rec = report$R[2:length(report$R)],
                        SSBPR0 = report$EPR0, NPR0 = report$NPR_virgin, weight = Wa,
                        mat = mat_age, M = M, vul = report$vul, SR = SR, fix_h = fix_h, h = info$h)
    #refpt <- get_refpt(SSB = report$E[1:(length(report$E) - 1)], rec = report$R[2:length(report$R)],
    #                   SSB0 = SSB0, R0 = R0, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    report <- c(report, refpt)

    if(rescale != 1) {
      vars_div <- c("meanR", "B", "E", "CAApred", "CN", "N", "VB", "R", "R_early", "MSY", "VBMSY",
                    "RMSY", "BMSY", "EMSY", "VB0", "R0", "B0", "E0", "N0")
      vars_mult <- "Brec"
      var_trans <- c("meanR", "q")
      fun_trans <- c("/", "*")
      fun_fixed <- c("log", NA)
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    }

    Yearplusone <- c(Year, max(Year) + 1)
    YearEarly <- (Year[1] - max_age + 1):(Year[1] - 1)
    YearDev <- c(YearEarly, Year)
    YearR <- c(YearDev, max(YearDev) + 1)
    R <- c(rev(report$R_early), report$R)

    Dev <- c(rev(report$log_early_rec_dev), report$log_rec_dev)
    if(integrate) {
      SE_Early <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_early_rec_dev"])
      SE_Main <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_rec_dev"])[map_log_rec_dev]
    } else {
      SE_Early <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_early_rec_dev"])
      SE_Main <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])[map_log_rec_dev]
    }
    SE_Dev <- c(rev(SE_Early), SE_Main)
    Assessment <- new("Assessment", Model = "SCA2", UMSY = report$UMSY,
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
                      R = structure(R, names = YearR),
                      N = structure(rowSums(report$N), names = Yearplusone),
                      N_at_age = report$N,
                      Selectivity = matrix(report$vul, nrow = length(Year),
                                           ncol = max_age, byrow = TRUE),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Obs_C_at_age = CAA_hist,
                      Catch = structure(colSums(t(report$CAApred) * Wa), names = Year),
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
class(SCA2) <- "Assess"


get_refpt2 <- function(SSB, rec, SSBPR0, NPR0, weight, mat, M, vul, SR, fix_h, h) {
  maxage <- length(M)
  SSB0 <- recpred <- sigmaR <- NULL

  solve_SR_par <- function(x, h = NULL) {
    R0 <- exp(x[1])
    SSB0 <<- R0 * SSBPR0
    if(!fix_h) {
      if(SR == "BH") h <- 0.2 + 0.8 * ilogit(x[2])
      if(SR == "Ricker") h <- 0.2 + exp(x[2])
    }

    if(SR == "BH") recpred <<- (0.8 * R0 * h * SSB)/(0.2 * SSBPR0 * R0 *(1-h)+(h-0.2)*SSB)
    if(SR == "Ricker") recpred <<- SSB/SSBPR0 * (5*h)^(1.25 * (1 - SSB/SSB0))
    sigmaR <<- sqrt(sum((log(rec/recpred))^2)/length(recpred))
    nLL <- -sum(dnorm(log(rec/recpred), 0, sigmaR, log = TRUE))
    return(nLL)
  }

  if(fix_h) {
    opt <- optimize(solve_SR_par, interval = c(-10, 10), h = h)$minimum # steepness
    R0 <- exp(opt)
  } else {
    opt <- nlminb(c(10, 10), solve_SR_par)
    R0 <- exp(opt$par[1])
    if(SR == "BH") h <- 0.2 + 0.8 * ilogit(opt$par[2])
    if(SR == "Ricker") h <- 0.2 + exp(opt$par[2])
  }
  if(SR == "BH") {
    Arec <- 4*h/(1-h)/SSBPR0
    Brec <- (5*h-1)/(1-h)/SSB0
  }
  if(SR == "Ricker") {
    Arec <- 1/SSBPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / SSB0
  }

  # virgin reference points
  N0 <- R0 * sum(NPR0)
  E0 <- SSB0
  VB0 <- R0 * sum(NPR0 * weight * vul)
  B0 <- R0 * sum(NPR0 * weight)

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
