


#' Stock-reduction analysis (SRA) for conditioning operating models
#'
#' Intended for conditioning operating models for data-limited stocks. From a historical time series of total catch, and potentially
#' age/length compositions and multiple indices of abundance, the SRA returns a range of values for depletion, selectivity,
#' unfished recruitment (R0), historical fishing effort, and recruitment deviations. This is done by sampling life history parameters
#' from an operating model provided by the user and fitting to the data in a statistical catch-at-age model. This function is intended
#' to generate potential depletion scenarios that could supported from sparse data.
#'
#' @param OM An object of class \linkS4class{OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
#' steepness, maturity parameters (L50 and L50_95), as well as catch/index uncertainty (Cobs, Iobs).
#' @param Chist A vector of historical catch, should be of length OM@@nyears. Ideally, the first year of the catch series represents
#' unfished conditions.
#' @param Index A matrix of historical indices of abundances, with rows indexing years and columns indexing fleets.
#' @param CAA Age compositions with nyears rows and OM@@maxage columns.
#' @param CAL Length compositions with nyears rows and columns index the length bin.
#' @param ML A vector of mean length observations (length OM@@nyears).
#' @param length_bin A vector for the midpoints of the length bins in \code{CAL}. All length bin widths should be equal in size.
#' @param I_type A character vector of length \code{nrow(Index)} to indicate the type of biomass for which the index follows. Either \code{"B"} for
#' total biomass, \code{"SSB"} for spawning biomass, or \code{"VB"} for vulnerable biomass. If \code{NULL}, "B" is used.
#' @param selectivity Whether to use logistic or dome selectivity in the stock-reduction analysis.
#' @param C_eq Equilibrium catch prior to the first year of the operating model. Zero implies unfished conditions in year one. Otherwise, this is used
#' would estimate depletion in year one.
#' @param cores Integer for the number of cores for the stock-reduction analysis.
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters (FALSE) or random effects (TRUE).
#' @param figure Logical, whether to plot diagnostic figures (histograms of estimated depletion and unfished recruitment, SRA outputs, model fits, etc.).
#'
#' @details One of indices, age compositions, or length compositions should be provided in addition to the historical catch.
#' Selectivity is fixed if no age or length compositions are provided.
#'
#' @note If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the SRA as well.
#'
#' @export
SRA_scope <- function(OM, Chist, Index = NULL, CAA = NULL, CAL = NULL, ML = NULL, length_bin = NULL, I_type = NULL,
                      selectivity = c("logistic", "dome"), C_eq = 0, cores = 1L, integrate = FALSE, figure = TRUE) {
  selectivity <- match.arg(selectivity)

  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- length(Chist)

  # Match number of historical years of catch to OM
  if(OM@nyears != nyears) {
    message("OM@nyears will be updated to length(Chist): ", nyears)
    OM@nyears <- nyears
  }

  # Sample life history parameters
  set.seed(OM@seed)
  StockPars <- SampleStockPars(OM, msg = FALSE)
  ObsPars <- SampleObsPars(OM)
  FleetPars <- SampleFleetPars(OM, msg = FALSE)

  # Interpolate missing catches
  if(any(is.na(Chist))) {
    message("One or more of the historical annual catch observations is missing. Linear interpolation has been used to fill these data.")
    Chistold <- Chist
    Chist <- approx(Chist)$y
    cond <- !is.na(Chistold)
    Chist[(1:nyears)[cond]] <- Chistold[(1:nyears)[cond]]
    print(data.frame(`Catches entered` = Chistold, `Catches interpolated` = Chist))
  }

  # Indices
  if(!is.null(Index)) {
    if(is.vector(Index)) {
      if(length(Index) != nyears) stop("Length of Index vector does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
      Index <- matrix(Index, ncol = 1)
    } else if(is.matrix(Index)) {
      if(nrow(Index) != nyears) stop("Number of rows of Index matrix does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    # Match index to I_type
    if(is.null(I_type)) I_type <- rep("VB", ncol(Index))
    if(length(I_type) != ncol(Index)) {
      stop("Length of I_type needs to be ", ncol(Index), call. = FALSE)
    }

    I_type_check <- match(I_type, c("VB", "SSB", "B"))
    if(any(is.na(I_type_check))) stop("I_type vector needs to be entries of either: \"SSB\", \"VB\", or \"B\".", call. = FALSE)
  } else {
    Index <- matrix(NA, ncol = 1, nrow = nyears)
    I_type <- "B"
  }

  # No comp data
  if(is.null(CAA) && is.null(CAL)) {
    message("No length or age compositions were provided. Selectivity is fixed to values from OM.")
  }
  fix_sel <- is.null(CAA) && is.null(CAL)

  # Process age comps
  if(!is.null(CAA)) {
    if(dim(CAA)[1] != nyears) {
      stop("Number of CAA rows (", dim(CAA)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAA)[2] != OM@maxage) {
      message("Number of CAA columns (", dim(CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no CAA for ages greater than ", dim(CAA)[2], " and filling with NA")
      addages <- OM@maxage - dim(CAA)[2]
      CAA2 <- matrix(NA, nrow = nrow(CAA), ncol = addages)
      CAA <- cbind(CAA, CAA2)
    }
  } else {
    CAA <- matrix(NA, nrow = nyears, ncol = maxage)
  }

  # Process length comps
  if(!is.null(CAL)) {
    if(is.null(length_bin)) {
      stop("You must specify the argument length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(CAL)[1] != nyears) {
      stop("Number of CAL rows (", dim(CAL)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAL)[2] != length(length_bin)) {
      stop("Number of CAL columns (", dim(CAL)[2], ") does not equal length(length_bin) (", length(length_bin), ").", call. = FALSE)
    }
  } else {
    CAL <- matrix(NA, nrow = nyears, ncol = length(StockPars$CAL_binsmid))
    length_bin <- StockPars$CAL_binsmid
  }
  weight_at_length <- StockPars$a * length_bin ^ StockPars$b

  # Process mean lengths
  if(!is.null(ML)) {
    if(length(ML) != nyears) stop("Mean length vector (ML) must be of length ", nyears, ".", call. = FALSE)
  } else ML <- rep(NA, nyears)

  # Catch matrix
  set.seed(OM@seed)
  Csd <- ObsPars$Csd
  Cat <- trlnorm(nsim * nyears, 1, Csd) * rep(Chist, each = nsim)
  Cat_new <- array(Cat, c(nsim, nyears))

  # Fit model
  message("Fitting model (", nsim, " simulations) ...")

  if(cores > 1) {
    DLMtool::setup(as.integer(cores))
    on.exit(snowfall::sfStop())
    res <- snowfall::sfLapply(1:nsim, SRA_scope_est, Catch = Cat_new, Index = Index, CAA = CAA, CAL = CAL, ML = ML, length_bin = length_bin,
                              wt_at_len = weight_at_length, I_type = I_type, C_eq = C_eq, selectivity = selectivity, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"),
                              integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  } else {
    res <- lapply(1:nsim, SRA_scope_est, Catch = Cat_new, Index = Index, CAA = CAA, CAL = CAL, ML = ML, length_bin = length_bin,
                  wt_at_len = weight_at_length, I_type = I_type, C_eq = C_eq, selectivity = selectivity, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"),
                  integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  }

  conv <- vapply(res, getElement, logical(1), name = "conv")
  message(sum(conv), " out of ", nsim , " iterations converged (", 100*sum(conv)/nsim, "%).\n")
  if(sum(conv) < nsim) {
    no_conv_ind <- !conv
    no_conv <- (conv)

    message("For non-converged iterations, values were re-sampled from converged iterations.\n")
  }

  ### R0
  OM@cpars$R0 <- vapply(res, getElement, numeric(1), name = "R0")
  message("Range of unfished recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))

  ### Depletion and init D
  OM@cpars$initD <- vapply(res, function(x) x$E[1]/x$E0, numeric(1))
  message("Range of initial depletion (OM@cpars$initD): ", paste(round(range(OM@cpars$initD), 2), collapse = " - "))

  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0, numeric(1))
  message("Range of depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "))

  ### Selectivity
  vul_par <- do.call(rbind, lapply(res, getElement, "vul_par"))
  vul <- do.call(rbind, lapply(res, getElement, "vul"))

  if(!fix_sel) {
    OM@isRel <- FALSE

    LFS <- MSEtool:::ilogit(vul_par[, 1]) * 0.75 * max(length_bin) # Actually L95 for logistics
    L50 <- LFS - exp(vul_par[, 2])

    if(selectivity == "logistic") {
      L5 <- 2 * L50 - LFS

      OM@cpars$L5 <- L5
      message("Range of OM@cpars$L5: ", paste(round(range(OM@cpars$L5), 2), collapse = " - "))

      OM@cpars$LFS <- LFS
      message("Range of OM@cpars$LFS: ", paste(round(range(OM@cpars$LFS), 2), collapse = " - "))

      OM@cpars$Vmaxlen <- rep(1, nsim)
      message("With logistic selectivity, setting OM@cpars$Vmaxlen = 1")

    } else {
      sd_asc <- sqrt((L50 - LFS)^2/log(4))
      L5 <- sd_asc * sqrt(-log(0.05)/0.5) + LFS

      OM@cpars$L5 <- L5
      message("Range of OM@cpars$L5: ", paste(round(range(OM@cpars$L5), 2), collapse = " - "))

      OM@cpars$LFS <- LFS
      message("Range of OM@cpars$LFS: ", paste(round(range(OM@cpars$LFS), 2), collapse = " - "))

      OM@cpars$Vmaxlen <- ilogit(vul_par[, 4])
      message("Range of OM@cpars$Vmaxlen: ", paste(round(range(OM@cpars$Vmaxlen), 2), collapse = " - "))
    }
  }

  ### Find
  U <- do.call(rbind, lapply(res, getElement, "U"))
  Find <- -log(1 - U)
  OM@cpars$Find <- Find
  message("Historical F trends set in OM@cpars$Find.")

  ### Rec devs
  log_rec_dev <- do.call(rbind, lapply(res, getElement, "log_rec_dev"))
  log_early_rec_dev <- do.call(rbind, lapply(res, getElement, "log_early_rec_dev"))
  message("Historical recruitment trends set in OM@cpars$Perr_y.")

  make_Perr <- function(x) exp(x$log_rec_dev - 0.5 * x$tau^2)
  Perr <- do.call(rbind, lapply(res, make_Perr))
  OM@cpars$Perr_y <- StockPars$Perr_y
  OM@cpars$Perr_y[, 1:(OM@maxage - 1)] <- 1
  OM@cpars$Perr_y[, OM@maxage:(OM@maxage + nyears - 1)] <- Perr

  ### Catch fits
  make_catch <- function(x, res, StockPars) {
    rowSums(res[[x]]$CAApred * t(StockPars$Wt_age[x, , 1:nyears]))
  }
  Cpred <- do.call(rbind, lapply(1:nsim, make_catch, res = res, StockPars = StockPars))

  ### Index fits
  Ipred <- array(do.call(cbind, lapply(res, getElement, "Ipred")), dim = c(nyears, ncol(Index), nsim))
  Ipred <- aperm(Ipred, c(3, 1, 2))
  if(all(is.na(Ipred))) {
    Ipred2 <- do.call(rbind, lapply(res, function(x) x$B[1:(length(x$B)-1)]))
    Ipred <- array(NA, dim = c(nsim, nyears, 1))
    Ipred[, , 1] <- Ipred2
  }

  ### Mean length fits
  MLpred <- do.call(rbind, lapply(res, getElement, "mlen_pred"))

  ### Mean age fits
  MApred <- do.call(cbind, lapply(res, function(x) x$CAApred %*% c(1:ncol(x$CAApred))/x$CN))

  ### Generate figures
  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2,3), mar = c(5, 4, 1, 1), oma = c(0, 0, 3, 0))

    # R0 histogram
    hist(OM@cpars$R0, main = "", xlab = expression(R[0]))

    # Depletion histograms
    hist(OM@cpars$initD, main = "", xlab = "Initial depletion")
    hist(OM@cpars$D, main = "", xlab = "Depletion")

    # selectivity
    matplot(matrix(length_bin, ncol = nsim, nrow = length(length_bin)), t(vul), typ = "l", col = "black",
            xlab = "Length", ylab = "Selectivity")
    abline(h = 0, col = "grey")

    # Find
    matplot(t(Find), typ = "l", col = "black", xlab = "Year", ylab = "F")
    abline(h = 0, col = "grey")

    # Perr
    matplot(t(Perr), typ = "l", col = "black", xlab = "Year", ylab = "Recruitment deviations", ylim = c(0, 1.1 * max(Perr)))
    abline(h = 0, col = "grey")

    mtext("Distributions of OM parameters", outer = TRUE, side = 3)


    # Sampled catches
    matplot(t(Cpred), typ = "l", col = "black", xlab = "Year", ylab = "Catch")
    lines(Chist, col = "red", lwd = 3)
    abline(h = 0, col = "grey")

    # Index fits
    for(i in 1:ncol(Index)) {
      matplot(t(Ipred[, , i]), typ = "l", col = "black", xlab = "Year", ylab = paste("Index", i))
      if(!all(is.na(Index))) lines(Index[, i], col = "red", lwd = 3)
      abline(h = 0, col = "grey")
    }

    # ML fits
    matplot(t(MLpred), typ = "l", col = "black", xlab = "Year", ylab = "Mean length")
    if(!all(is.na(CAL))) {
      lines(CAL %*% length_bin/rowSums(CAL, na.rm = TRUE), col = "red", lwd = 3)
    } else if(!all(is.na(ML))) lines(ML, col = "red", lwd = 3)

    # MA fits
    matplot(MApred, typ = "l", col = "black", xlab = "Year", ylab = "Mean age")
    if(!all(is.na(CAA))) lines(CAA %*% c(1:ncol(CAA))/rowSums(CAA, na.rm = TRUE), col = "red", lwd = 3)

    mtext("Observed (red) and predicted data (black) \nfrom operating model conditioning", outer = TRUE, side = 3)
  }

  return(OM)
}

#### Parameters
# R0
# sel
# rec_devs

#### Data list
# Cat vector by year
# Index matrix by year and fleet
# Index_type vector by fleet
# CAA matrix by year, age
# CAL matrix by year, length bin
# mean length vector by year

#### Life history
# M array by nsim, age, year
# mat_age array by nsim, age
# Length at age array by nsim, age
# Weight at age arraty by nsim, length
# CV_LAA vector by nsim
# h steepness

#### Other settings
# SR type
# selectivity type
# maxage integer
SRA_scope_est <- function(x, SRA_scope_est, Catch, Index = NULL, CAA = NULL, CAL = NULL, ML = NULL, length_bin,
                          I_type, C_eq = 0, selectivity = c("logistic", "dome"), SR_type = c("BH", "Ricker"), n = c(30, 30), StockPars, ObsPars,
                          FleetPars, wt_at_len, integrate = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3)) {
  selectivity <- match.arg(selectivity)
  SR_type <- match.arg(SR_type)
  I_type2 <- as.integer(ifelse(I_type == "B", 0, ifelse(I_type == "VB", 1, 2)))

  nyears <- length(Catch[x, ])
  max_age <- ncol(CAA)

  TMB_data <- list(model = "SRA_scope", C_hist = Catch[x, ], C_eq = C_eq, I_hist = Index, nfleet = ncol(Index),
                   CAA_hist = CAA, CAA_n = pmin(rowSums(CAA, na.rm = TRUE), n[1]),
                   CAL_hist = CAL, CAL_n = pmin(rowSums(CAL, na.rm = TRUE), n[2]), length_bin = length_bin, mlen = ML, n_y = nyears,
                   max_age = ncol(CAA), M = t(StockPars$M_ageArray[x, , 1:nyears]), len_age = t(StockPars$Len_age[x, , 1:nyears]),
                   CV_LAA = StockPars$LenCV[x], wt_at_len = wt_at_len, mat = t(StockPars$Mat_age[x, , 1:nyears]),
                   vul_type = selectivity, I_type = I_type2, SR_type = SR_type, CAA_dist = "multinomial",
                   est_C_eq = as.integer(C_eq > 0), est_early_rec_dev = rep(NA, max_age - 1), est_rec_dev = rep(1, nyears))

  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)

  L95 <- FleetPars$LFS[nyears, x]
  L50 <- mean(c(FleetPars$LFS[nyears, x], FleetPars$L5[nyears, x]))
  if(selectivity == "logistic") {
    vul_par <- c(logit(L95/max(length_bin)/0.75), log(L95 - L50))
  } else {
    vul_par <- c(logit(L95/max(length_bin)/0.75), log(L95 - L50), -20, logit(FleetPars$Vmaxlen[nyears, x]))
  }

  TMB_params <- list(log_R0 = 10, transformed_h = transformed_h, log_U_equilibrium = -10, vul_par = vul_par,
                     log_sigma_I = log(rep(sdconv(1, ObsPars$Isd[x]), ncol(Index))),
                     log_sigma_mlen = log(0.2), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = rep(0, max_age - 1), log_rec_dev = rep(0, nyears))

  map <- list()
  map$transformed_h <- map$log_tau <- factor(NA)
  map$log_sigma_I <- factor(rep(NA, ncol(Index)))
  if(C_eq == 0) map$log_U_equilibrium <- factor(NA)
  if(selectivity == "dome") map$vul_par <- factor(c(1, 2, NA, 3))
  if(all(is.na(ML))) map$log_sigma_mlen <- factor(NA)
  if(all(is.na(CAA)) && all(is.na(CAL))) {
    map$vul_par <- factor(rep(NA, length(TMB_params$vul_par)))
    map$log_sigma_mlen <- factor(NA)
  }
  map$log_early_rec_dev <- factor(rep(NA, max_age - 1))

  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev") else random <- NULL
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, map = map, random = random, inner.control = inner.control,
                   DLL = "MSEtool", silent = TRUE)
  mod <- optimize_TMB_model(obj, control, restart = 1)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  return(c(report, list(conv = !is.character(opt) && SD$pdHess)))
}

