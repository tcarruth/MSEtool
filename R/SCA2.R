#' @describeIn SCA The parameterization by Forrest et al. (2008) is used where UMSY and MSY are leading parameters.
#' Recruitment deviations are only estimated beginning in the year when age composition are available.
#' @export
SCA2 <- function(x = 1, Data, U_begin = c("virgin", "est"), vulnerability = c("logistic", "dome"),
                 SR = c("BH", "Ricker"), CAA_multiplier = 50, start = NULL, tau = 1, integrate = TRUE,
                 silent = TRUE, control = list(eval.max = 1e3), inner.control = list(), ...) {
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
    if(!is.null(start$UMSY) && is.numeric(start$UMSY)) params$logit_UMSY <- log(start$UMSY[1]/(1 - start$UMSY[1]))
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$log_MSY <- log(start$MSY[1])
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
    if(!is.null(start$vulnerability) && is.numeric(start$vulnerability)) params$vul_par <- start$vul_par
  }
  if(is.null(params$logit_UMSY)) {
    UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
    params$logit_UMSY <- log(UMSY_start/(1 - UMSY_start))
  }
  if(is.null(params$log_MSY)) {
    AvC <- mean(C_hist, na.rm = TRUE)
    params$log_MSY <- log(3 * AvC)
  }
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
  params$log_rec_dev <- rep(0, n_y - 1)

  CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
  if(vulnerability == "logistic") {
    vul_par <- c(CAA_mode-1, log(1)) # 50 and log(95%-offset) vulnerability respectively
  }
  if(vulnerability == "dome") {
    # double normal: logsd(ascending), mean(asc), mean(desc)-logoffset, sd(desc)
    vul_par <- c(log(1), CAA_mode, log(0.5), log(5))
  }

  data <- list(model = "SCA", C_hist = C_hist, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = vulnerability,
               est_rec_dev = as.integer(random_map(CAA_n_nominal)))

  info <- list(Year = Data@Year, data = data, params = params, control = control)

  map <- list()
  if(U_begin == "virgin") map$U_equilibrium <- factor(NA)
  if(any(is.na(CAA_n_nominal))) {
    map$log_rec_dev <- random_map(CAA_n_nominal)
  }
  random <- NULL
  if(integrate) random <- "log_rec_dev"
  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
  opt <- optimize_TMB_model(obj, control)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)
  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA2", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearrandom <- seq(Year[2], max(Year))
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
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$nll_comp),
                                      names = c("Total", "Index", "CAA", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[6],
                      SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      SE_SSB_SSBMSY_final = SD$sd[9], SE_SSB_SSB0_final = SD$sd[10],
                      SE_VB_VBMSY_final = SD$sd[11], SE_VB_VB0_final = SD$sd[12],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SCA2) <- "Assess"

