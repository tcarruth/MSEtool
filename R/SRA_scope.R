


#' Stock-reduction analysis (SRA) for conditioning operating models
#'
#' Intended for conditioning operating models for data-limited stocks. From a historical time series of total catch, and potentially
#' age/length compositions and multiple indices of abundance, the SRA returns a range of values for depletion, selectivity,
#' unfished recruitment (R0), historical fishing effort, and recruitment deviations for the operating model. This is done by sampling life history parameters
#' provided by the user and fitting to the data in a statistical catch-at-age model (with the predicted catch equal to the observed catch).
#' This function is intended to generate a range of potential depletion scenarios that could be supported from sparse data. A full catch time series
#' is needed but missing data (as NAs) are allowed for all other data types.
#'
#' @param OM An object of class \linkS4class{OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
#' steepness, maturity parameters (L50 and L50_95), standard deviation of recruitment variability (Perr), as well as index uncertainty (Iobs).
#' @param Chist A vector of historical catch, should be of length OM@@nyears. If there are multiple fleets: a matrix of OM@@nyears rows and nfleet columns.
#' Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).
#' @param Index A vector of values of an index (of length OM@@nyears). If there are multiple surveys: a matrix of historical indices of abundances, with rows
#' indexing years and columns indexing surveys.
#' @param I_sd A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
#' If NULL, this function will use values from OM@@Iobs.
#' @param CAA Age composition matrix with nyears rows and OM@@maxage columns. If multiple fleets: an array with dimension: nyears, OM@@maxage, and nfleets.
#' @param CAL Length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
#' length bins, and nfleets.
#' @param ML A vector of mean length observations (length OM@@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets. Generally, should not
#' be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.
#' @param length_bin A vector for the midpoints of the length bins for \code{CAL}. All length bin widths should be equal in size.
#' @param C_eq A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
#' Zero implies unfished conditions in year one. Otherwise, this is used would estimate depletion in the first year of the data.
#' @param ML_sd The standard deviation (normal distribution) of the observed mean lengths. If there are multiple fleets, a vector of length nfleet.
#' If \code{NULL}, default value is 0.1.
#' @param selectivity A character vector of length nfleet to indicate \code{"logistic"} or \code{"dome"} selectivity for each fleet in \code{Chist}.
#' @param I_type A character vector (length nsurvey) to indicate the type of biomass for which each index follows. Either \code{"B"} for
#' total biomass, or \code{"SSB"} for spawning biomass. If \code{NULL}, "B" is used. Use numbers if the index corresponds to a fleet in \code{Chist}.
#' @param LWT A named list of likelihood weights for the SRA model. See details.
#' @param ESS A numeric vector of length two for the maximum effective samples size of the age and length compositions, respectively for the
#' multinomial likelihood function. The annual sample size of an age or length composition sample is the minimum of ESS or the number of observations.
#' @param cores Integer for the number of CPU cores for the stock reduction analysis.
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters (FALSE) or random effects (TRUE).
#' @param figure Logical, whether to plot diagnostic figures (histograms of estimated depletion and unfished recruitment, SRA outputs, model fits, etc.).
#' @param Year A vector of years for the historical period, used for plotting.
#' @param report Logical, whether to return all assessment output. See value section below.
#' @param report_list The list of assessment output returned by \code{SRA_scope} when \code{report = TRUE}.
#' @return
#' A named list containing the following:
#' \enumerate{
#' \item "OM" - an updated operating model with depletion, F, selectivity, and recruitment deviations from the SRA fits.
#' \item "output" - A list of output, e.g. spawning biomass and predicted catch at age, from the SRA fits:
#'
#' \itemize{
#' \item SSB - A matrix of \code{OM@@nsim} rows and \code{OM@@nyears+1} columns for estimated spawning biomass
#' \item N - An array of dimension \code{c(nsim, nyears+1, maxage)} for estimated abundance by simulation, year, and age.
#' \item CAA - An array of dimension \code{c(nsim, nyears+1, maxage, and nfleet)} for estimated catch at age by simulation, year, age, and fleet.
#' \item CAL - An array of dimension \code{c(nsim, nyears+1, maxage, and nfleet)} for estimated catch at length by simulation, year, length bin, and fleet.
#' \item conv - A logical vector of length nsim indicating convergence of the SRA in the i-th simulation.
#' }
#'
#' \item "report" - If \code{report = TRUE}, a list of length \code{OM@@nsim} containing all the assessment output from each model fit, otherwise returns NULL.
#' Less organized than "output".
#' }
#' If \code{figure = TRUE}, a set of diagnostic plots for the fits to the SRA for each simulation as well histograms of operating model parameters,
#' e.g., depletion.
#'
#' @details
#' \code{plot_SRA_scope} generates the plots from the SRA scope function.
#'
#' One of indices, age compositions, or length compositions should be provided in addition to the historical catch.
#' Selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.
#'
#' \code{LWT} is a named list containing the likelihood weights (values > 0) with the possible options:
#' \itemize{
#' \item Chist: A vector of length nfleet.
#' \item Index: A vector of length nsurvey.
#' \item CAA, CAL, ML, C_eq: A vector of length nfleet for each.
#' }
#' By default, all likelihood weights are equal to one if not specified by the user. Likelihoods for CAA and CAL can also be adjusted by changing the
#' multinomial sample size. See argument \code{ESS}.
#'
#' @note If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the SRA as well.
#' Time-varying life history parameters create ambiguity in the calculation and interpretation of depletion and reference points in \link[DLMtool]{runMSE}.
#' See section D.5 of \code{DLMtool::userguide()}.
#'
#' Here, the initial depletion (OM@@cpars$initD) is calculated based on unfished spawning biomass using growth and M in the first year. If growth/M varies,
#' then this reference point may no longer be relevant at the end of the historical period and in the projection period.
#'
#' The easiest way to turn off time-varying growth/M is by setting: \code{OM@@Msd <- OM@@Linfsd <- OM@@Ksd <- c(0, 0)}.
#'
#' @author Q. Huynh
#'
#' @export
SRA_scope <- function(OM, Chist, Index = NULL, I_sd = NULL, CAA = NULL, CAL = NULL, ML = NULL, length_bin = NULL, C_eq = 0, ML_sd = NULL,
                      selectivity = "logistic", I_type = NULL, LWT = list(), ESS = c(30, 30), cores = 1L,
                      integrate = FALSE, figure = TRUE, Year = NULL, report = FALSE) {

  # Convert single fleet inputs to multiple fleet, e.g. matrices to arrays
  if(!is.matrix(Chist)) Chist <- matrix(Chist, ncol = 1)

  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- nrow(Chist)
  nfleet <- ncol(Chist)
  message(nfleet, " fleet(s) detected.")


  # Match number of historical years of catch to OM
  if(OM@nyears != nyears) {
    message("OM@nyears will be updated to length(Chist): ", nyears)
    OM@nyears <- nyears
  }

  # Interpolate missing catches
  if(any(is.na(Chist))) {
    stop("One or more of the historical annual catch observations is missing. Suggestion: use linear interpolation to fill these data.")
  }

  # Indices
  if(!is.null(Index)) {
    if(is.vector(Index)) {
      if(length(Index) != nyears) stop("Length of Index vector does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
      Index <- matrix(Index, ncol = 1)
    } else if(is.matrix(Index)) {
      if(nrow(Index) != nyears) stop("Number of rows of Index matrix does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    nsurvey <- ncol(Index)

    # Match index to I_type
    if(is.null(I_type)) I_type <- rep("B", ncol(Index))
    if(length(I_type) != ncol(Index)) {
      stop("Length of I_type needs to be ", ncol(Index), call. = FALSE)
    }

    I_type_check <- match(I_type, c("B", "SSB", 1:ncol(Chist)))
    if(any(is.na(I_type_check))) stop("I_type vector needs to be entries of either: \"SSB\", \"B\", or 1 - ", ncol(Chist), ".", call. = FALSE)
  } else {
    nsurvey <- 0
    Index <- matrix(NA, ncol = 1, nrow = nyears)
    I_type <- "B"
  }
  message(nsurvey, " survey(s) detected.")

  I_type2 <- suppressWarnings(as.numeric(I_type))
  I_type2[I_type == "B"] <- -1
  I_type2[I_type == "SSB"] <- -2

  if(!is.null(I_sd)) {
    if(is.vector(I_sd)) {
      if(length(I_sd) != nyears) stop("Length of I_sd vector does not equal nyears (", nyears, ").", call. = FALSE)
      I_sd <- matrix(I_sd, ncol = 1)
    } else if(is.matrix(I_sd)) {
      if(nrow(I_sd) != nyears) stop("Number of rows of I_sd matrix does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(I_sd) != nsurvey) stop("Number of columns of I_sd matrix does not equal nsurvey (", nsurvey, ").", call. = FALSE)
    } else I_sd <- NULL
  }

  # No comp data
  if(is.null(CAA) && is.null(CAL)) {
    fix_sel <- TRUE
    message("No length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }

  # Selectivity
  if(!fix_sel) {
    sel_test <- match(selectivity, c("logistic", "dome"))
    if(any(is.na(sel_test))) stop("selectivity vector should be either \"logistic\" or \"dome\".", call. = FALSE)

    if(length(selectivity) < nfleet) stop("selectivity vector should be of length nfleet (", nfleet, ").", call. = FALSE)
  }
  sel <- ifelse(selectivity == "logistic", 1L, 0L)

  # Process age comps
  if(!is.null(CAA)) {

    if(is.matrix(CAA)) {
      CAA2 <- array(NA, c(dim(CAA), 1))
      CAA2[,,1] <- CAA
      CAA <- CAA2
    }

    if(dim(CAA)[1] != nyears) {
      stop("Number of CAA rows (", dim(CAA)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAA)[2] < OM@maxage) {
      message("Number of CAA columns (", dim(CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no CAA for ages greater than ", dim(CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(CAA)[2]
      CAA_new <- array(0, c(nyears, OM@maxage, nfleet))
      CAA_new[, 1:dim(CAA)[2], ] <- CAA
      CAA <- CAA_new
    }
    if(dim(CAA)[2] > OM@maxage) {
      OM@maxage <- maxage <- dim(CAA)[2]
      message("Increasing OM@maxage to ", OM@maxage, ".")
    }
    if(dim(CAA)[3] != nfleet) {
      stop("Number of CAA slices (", dim(CAA)[3], ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    CAA <- array(0, c(nyears, maxage, nfleet))
  }




  # Sample life history, selectivity, and obs parameters
  set.seed(OM@seed)
  StockPars <- SampleStockPars(OM, msg = FALSE)
  ObsPars <- SampleObsPars(OM)
  FleetPars <- SampleFleetPars(OM, msg = FALSE)




  # Process length comps
  if(!is.null(CAL)) {

    if(is.matrix(CAL)) {
      CAL2 <- array(NA, c(dim(CAL), 1))
      CAL2[,,1] <- CAL
      CAL <- CAL2
    }
    if(is.null(length_bin)) {
      stop("You must specify the argument length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(CAL)[1] != nyears) {
      stop("Number of CAL rows (", dim(CAL)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAL)[2] != length(length_bin)) {
      stop("Number of CAL columns (", dim(CAL)[2], ") does not equal length(length_bin) (", length(length_bin), ").", call. = FALSE)
    }
    if(dim(CAL)[3] != nfleet) {
      stop("Number of CAL slices (", dim(CAA)[3], ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    CAL <- array(0, c(nyears, length(StockPars$CAL_binsmid), nfleet))
    length_bin <- StockPars$CAL_binsmid
  }
  weight_at_length <- StockPars$a * length_bin ^ StockPars$b

  # Process mean lengths
  if(!is.null(ML)) {
    if(is.vector(ML)) {
      if(length(ML) != nyears) stop("Mean length vector (ML) must be of length ", nyears, ".", call. = FALSE)
      ML <- matrix(ML, ncol = 1)
    }
    if(nrow(ML) != nyears) stop("Number of ML rows (", nrow(ML), ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    if(ncol(ML) != nfleet) stop("Number of ML columns (", ncol(ML), ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)

    if(is.null(ML_sd)) {
      ML_sd <- rep(0.1, nfleet)
    } else if(length(ML_sd) == 1) ML_sd <- rep(ML_sd, nfleet)
    if(length(ML_sd) != nfleet) stop("Mean length SD vector (ML_sd) must be of length ", nfleet, ".", call. = FALSE)

  } else {
    ML <- matrix(0, nrow = nyears, ncol = nfleet)
    ML_sd <- rep(0.1, nfleet)
  }

  # Ceq
  if(length(C_eq) == 1) C_eq <- rep(C_eq, nfleet)
  if(length(C_eq) < nfleet) stop("C_eq needs to be of length nfleet (", nfleet, ").", call. = FALSE)

  # Likelihood weights
  if(is.null(LWT$Chist)) {
    LWT$Chist <- rep(1, nfleet)
  } else if(length(LWT$Chist) == 1 && nfleet > 1) {
    LWT$Chist <- rep(LWT$Chist, nfleet)
  }
  if(length(LWT$Chist) != nfleet) stop("LWT$Chist should be a vector of length ", nfleet, ".")

  if(is.null(LWT$Index)) {
    LWT$Index <- rep(1, max(1, nsurvey))
  } else if(length(LWT$Index) == 1 && nsurvey > 1) {
    LWT$Index <- rep(LWT$Index, nsurvey)
  }
  if(length(LWT$Index) != max(1, nsurvey)) stop("LWT$Index should be a vector of length ", nsurvey, ".")

  if(is.null(LWT$CAA)) {
    LWT$CAA <- rep(1, nfleet)
  } else if(length(LWT$CAA) == 1 && nfleet > 1) {
    LWT$CAA <- rep(LWT$CAA, nfleet)
  }
  if(length(LWT$CAA) != nfleet) stop("LWT$CAA should be a vector of length ", nfleet, ".")

  if(is.null(LWT$CAL)) {
    LWT$CAL <- rep(1, nfleet)
  } else if(length(LWT$CAL) == 1 && nfleet > 1) {
    LWT$CAL <- rep(LWT$CAL, nfleet)
  }
  if(length(LWT$CAL) != nfleet) stop("LWT$CAL should be a vector of length ", nfleet, ".")

  if(is.null(LWT$ML)) {
    LWT$ML <- rep(1, nfleet)
  } else if(length(LWT$ML) == 1 && nfleet > 1) {
    LWT$ML <- rep(LWT$ML, nfleet)
  }
  if(length(LWT$ML) != nfleet) stop("LWT$ML should be a vector of length ", nfleet, ".")

  if(is.null(LWT$C_eq)) {
    LWT$C_eq <- rep(1, max(1, nfleet))
  } else if(length(LWT$C_eq) == 1 && nfleet > 1) {
    LWT$C_eq <- rep(LWT$C_eq, nfleet)
  }
  if(length(LWT$C_eq) != nfleet) stop("LWT$C_eq should be a vector of length ", nfleet, ".")

  # SR
  SR <- ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker")
  message(SR, " stock-recruitment relationship used.")

  # Fit model
  message("Fitting model (", nsim, " simulations) ...")

  if(cores > 1) {
    DLMtool::setup(as.integer(cores))
    on.exit(snowfall::sfStop())
    mod <- snowfall::sfLapply(1:nsim, SRA_scope_est3, Catch = Chist, Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML,
                              ML_sd = ML_sd, length_bin = length_bin,
                              wt_at_len = weight_at_length, I_type = I_type2, C_eq = C_eq, selectivity = sel,
                              fix_selectivity = fix_sel, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, ESS = ESS,
                              integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  } else {
    mod <- lapply(1:nsim, SRA_scope_est3, Catch = Chist, Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML,
                  ML_sd = ML_sd, length_bin = length_bin,
                  wt_at_len = weight_at_length, I_type = I_type2, C_eq = C_eq, selectivity = sel,
                  fix_selectivity = fix_sel, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, ESS = ESS,
                  integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  }
  #assign('mod', mod, envir = globalenv())
  res <- lapply(mod, getElement, "report")
  conv <- vapply(res, getElement, logical(1), name = "conv")
  message(sum(conv), " out of ", nsim , " iterations converged (", 100*sum(conv)/nsim, "%).\n")
  #if(sum(conv) < nsim) {
  #  no_conv_ind <- !conv
  #  #no_conv <- (conv)
  #  message("For non-converged iterations, values were re-sampled from converged iterations.\n")
  #}

  ### R0
  OM@cpars$R0 <- vapply(res, getElement, numeric(1), name = "R0")
  message("Range of unfished recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))

  ### Depletion and init D
  OM@cpars$initD <- vapply(res, function(x) x$E[1]/x$E0_year1, numeric(1))
  message("Range of initial spawning depletion (OM@cpars$initD): ", paste(round(range(OM@cpars$initD), 2), collapse = " - "))

  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0, numeric(1))
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")

  ### Selectivity and F
  ### Find
  get_F_at_age <- function(report) {
    N <- report$N
    max_age <- ncol(N)
    nyears <- nrow(N) - 1

    surv_matrix <- N[2:(nyears+1), 2:max_age]/N[1:nyears, 1:(max_age-1)]
    surv_plusgroup <- N[2:(nyears+1), max_age]/(N[1:nyears, max_age-1] + N[1:nyears, max_age])

    Z_matrix <- -log(surv_matrix)
    Z_matrix[, max_age-1] <- -log(surv_plusgroup)

    Z_matrix <- cbind(Z_matrix, -log(surv_plusgroup))

    F_matrix <- Z_matrix - report$M
    F_matrix[F_matrix < 0] <- 0

    return(F_matrix)
  }
  F_matrix <- lapply(res, get_F_at_age)
  apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))
  Find <- do.call(rbind, apical_F)

  V <- Map("/", e1 = F_matrix, e2 = apical_F)
  expand_V_matrix <- function(x) {
    y <- matrix(x[nyears, ], proyears, maxage, byrow = TRUE)
    rbind(x, y)
  }
  V <- lapply(V, expand_V_matrix)
  V2 <- array(unlist(V), c(nyears+proyears, maxage, nsim))
  V2 <- aperm(V2, c(3, 2, 1))

  if(nfleet > 1) {
    OM@cpars$V <- V2
    OM@cpars$Find <- Find
    message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
    message("Selectivity during projection period is equal to that in most recent historical year.")
  } else {
    if(!fix_sel) {
      OM@isRel <- FALSE

      vul_par <- do.call(cbind, lapply(res, getElement, "vul_par"))
      LFS <- ilogit(vul_par[1, ]) * 0.75 * max(length_bin) # Actually L95 for logistics
      L50 <- LFS - exp(vul_par[2, ])

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

        OM@cpars$Vmaxlen <- ilogit(vul_par[4, ])
        message("Range of OM@cpars$Vmaxlen: ", paste(round(range(OM@cpars$Vmaxlen), 2), collapse = " - "))
      }
    }

    F_vector <- do.call(cbind, lapply(res, getElement, "F"))
    OM@cpars$Find <- t(F_vector)

    message("Historical F set in OM@cpars$Find.")
  }

  Eff <- apply(OM@cpars$Find, 2, range)
  OM@EffLower <- Eff[1, ]
  OM@EffUpper <- Eff[2, ]
  if(length(OM@EffYears) != nyears) OM@EffYears <- 1:nyears
  message("Historical effort trends set in OM@EffLower and OM@EffUpper.\n")

  ### Rec devs
  OM@cpars$Perr <- StockPars$procsd
  make_Perr <- function(x) exp(x$log_rec_dev - 0.5 * x$tau^2)
  Perr <- do.call(rbind, lapply(res, make_Perr))

  make_early_Perr <- function(x) exp(x$log_early_rec_dev - 0.5 * x$tau^2)
  early_Perr <- do.call(rbind, lapply(res, make_early_Perr))

  OM@cpars$Perr_y <- StockPars$Perr_y
  OM@cpars$Perr_y[, 1:(OM@maxage - 1)] <- early_Perr
  OM@cpars$Perr_y[, OM@maxage:(OM@maxage + nyears - 1)] <- Perr

  log_rec_dev <- do.call(rbind, lapply(res, getElement, "log_rec_dev"))
  log_early_rec_dev <- do.call(rbind, lapply(res, getElement, "log_early_rec_dev"))

  OM@cpars$AC <- apply(log_rec_dev, 1, function(x) acf(x, lag.max = 1, plot = FALSE)$acf[2])
  OM@AC <- range(OM@cpars$AC)

  pro_Perr_y <- matrix(rnorm(proyears * nsim, rep(StockPars$procmu, proyears), rep(StockPars$procsd, proyears)),
                       c(nsim, proyears))
  for(y in 2:proyears) pro_Perr_y[, y] <- OM@cpars$AC * pro_Perr_y[, y - 1] + pro_Perr_y[, y] * sqrt(1 - OM@cpars$AC^2)
  OM@cpars$Perr_y[, (OM@maxage+nyears):ncol(OM@cpars$Perr_y)] <- exp(pro_Perr_y)

  message("Recruitment standard deviations set in OM@cpars$Perr.")
  message("Historical recruitment trends set in OM@cpars$Perr_y.")
  message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
  message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")

  ### Assign OM variables that were used in the SRA to the output
  OM@cpars$Len_age <- StockPars$Len_age
  OM@cpars$Linf <- StockPars$Linf
  OM@cpars$K <- StockPars$K
  OM@cpars$t0 <- StockPars$t0
  OM@cpars$LenCV <- StockPars$LenCV

  OM@cpars$Mat_age <- StockPars$Mat_age
  OM@cpars$L50 <- StockPars$L50
  OM@cpars$L95 <- StockPars$L95

  OM@cpars$M_ageArray <- StockPars$M_ageArray

  OM@cpars$h <- StockPars$hs

  if(fix_sel) {
    OM@cpars$L5 <- if(is.matrix(FleetPars$L5)) FleetPars$L5[nyears, ] else FleetPars$L5
    OM@cpars$LFS <- if(is.matrix(FleetPars$LFS)) FleetPars$LFS[nyears, ] else FleetPars$LFS
    OM@cpars$Vmaxlen <- if(is.matrix(FleetPars$Vmaxlen)) FleetPars$Vmaxlen[nyears, ] else FleetPars$Vmaxlen
    V <- FleetPars$V
    maxV <- apply(FleetPars$V, c(1, 3), max)
    for(i in 1:maxage) V[,i,] <- V[,i,]/maxV
    OM@cpars$V <- V
  }

  OM@cpars$Iobs <- ObsPars$Iobs

  message("Growth, maturity, natural mortality, and steepness values from SRA are set in OM@cpars.\n")

  ### Output list
  E <- do.call(rbind, lapply(res, getElement, "E"))

  #N <- lapply(res, getElement, "N_full")
  #N_length <- lapply(N, function(x) rowSums(aperm(x, c(1, 3, 2)), dims = 2))
  N_age <- lapply(res, getElement, "N")
  N_age2 <- array(unlist(N_age), c(nyears+1, maxage, nsim))

  CAA_pred <- lapply(res, getElement, "CAApred")
  CAA_pred2 <- array(unlist(CAA_pred), c(nyears, maxage, nfleet, nsim))

  CAL_pred <- lapply(res, getElement, "CALpred")
  CAL_pred2 <- array(unlist(CAL_pred), c(nyears, length(length_bin), nfleet, nsim))

  output <- list(SSB = E, N = aperm(N_age2, c(3, 1, 2)),
                 CAA = aperm(CAA_pred2, c(4, 1, 2, 3)), CAL = aperm(CAL_pred2, c(4, 1, 2, 3)), conv = conv)

  ### Generate figures
  if(figure) plot_SRA_scope(OM, Chist, Index, CAA, CAL, ML, res, OM@EffYears)

  message("Complete.")
  return(list(OM = OM, output = output, report = if(report) res else NULL))
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
SRA_scope_est3 <- function(x, SRA_scope_est, Catch, Index = NULL, I_sd = NULL, CAA = NULL, CAL = NULL, ML = NULL, ML_sd = NULL, length_bin,
                           I_type, C_eq = 0, SR_type = c("BH", "Ricker"), LWT = list(), ESS = c(30, 30),
                           StockPars, ObsPars, FleetPars, wt_at_len, integrate = FALSE, selectivity, fix_selectivity = TRUE,
                           control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3)) {
  SR_type <- match.arg(SR_type)

  nyears <- nrow(Catch)
  max_age <- dim(CAA)[2]
  nfleet <- ncol(Catch)
  nsurvey <- ncol(Index)

  if(is.null(I_sd)) I_sd <- matrix(sdconv(1, ObsPars$Isd[x]), nyears, nsurvey)

  CAA_n <- t(rowSums(aperm(CAA, c(3, 1, 2)), dims = 2, na.rm = TRUE))
  CAL_n <- t(rowSums(aperm(CAL, c(3, 1, 2)), dims = 2, na.rm = TRUE))

  rescale <- 1/mean(Catch, na.rm = TRUE)

  LWT_C <- matrix(c(LWT$Chist, LWT$CAA, LWT$CAL, LWT$ML, LWT$C_eq), nrow = nfleet, ncol = 5)

  TMB_data <- list(model = "SRA_scope",
                   C_hist = Catch * rescale, C_eq = C_eq * rescale, I_hist = Index, sigma_I = I_sd,
                   CAA_hist = CAA, CAA_n = pmin(CAA_n, ESS[1]),
                   CAL_hist = CAL, CAL_n = pmin(CAL_n, ESS[2]), length_bin = length_bin, mlen = ML,
                   n_y = nyears, max_age = ncol(CAA), nfleet = nfleet, nsurvey = nsurvey,
                   M = t(StockPars$M_ageArray[x, , 1:nyears]), len_age = t(StockPars$Len_age[x, , 1:nyears]),
                   CV_LAA = StockPars$LenCV[x], wt_at_len = wt_at_len, mat = t(StockPars$Mat_age[x, , 1:nyears]),
                   vul_type = selectivity, I_type = I_type, SR_type = SR_type, LWT_C = LWT_C, LWT_Index = LWT$Index,
                   est_early_rec_dev = rep(NA, max_age - 1), est_rec_dev = c(rep(1, nyears-1), NA))

  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)

  L95 <- FleetPars$LFS[nyears, x]
  L50 <- mean(c(FleetPars$LFS[nyears, x], FleetPars$L5[nyears, x]))

  vul_par <- matrix(0, 4, nfleet)
  map_vul_par <- matrix(ifelse(fix_selectivity, NA, 0), 4, nfleet)

  for(ff in 1:nfleet) {
    if(selectivity[ff]) { #logistic
      vul_par[1:2, ff] <- c(logit(L95/max(length_bin)/0.75), log(L95 - L50))
      map_vul_par[3:4, ff] <- NA
    } else {
      vul_par[, ff] <- c(logit(L95/max(length_bin)/0.75), log(L95 - L50), -20, logit(FleetPars$Vmaxlen[nyears, x]))
      map_vul_par[3, ff] <- NA
    }
  }

  if(!fix_selectivity) {
    map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))
  }

  TMB_params <- list(log_R0 = 3, transformed_h = transformed_h, vul_par = vul_par,
                     log_F = matrix(log(0.05), nyears, nfleet), log_F_equilibrium = rep(log(0.05), nfleet),
                     log_sigma_mlen = log(ML_sd), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = rep(0, max_age - 1), log_rec_dev = rep(0, nyears))

  map <- list()
  map$transformed_h <- map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  if(any(C_eq == 0)) {
    map_log_F_equilibrium <- rep(NA, nfleet)
    map_log_F_equilibrium[C_eq > 0] <- 1:sum(C_eq > 0)
    map$log_F_equilibrium <- factor(map_log_F_equilibrium)
  }

  map$log_sigma_mlen <- factor(rep(NA, nfleet))
  map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
  map$log_rec_dev <- factor(c(1:(nyears-1), NA))

  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev") else random <- NULL

  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, map = map, random = random, inner.control = inner.control,
                   DLL = "MSEtool", silent = TRUE)

  mod <- optimize_TMB_model(obj, control, restart = 1)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  report$vul <- do.call(cbind, report$vul)

  vars_div <- c("B", "E", "Cat", "C_eq_pred", "CAApred", "CALpred", "CN", "Cpred", "N", "N_full", "VB",
                "R", "R_early", "R_eq", "VB0", "R0", "B0", "E0", "N0", "E0_year1")
  vars_mult <- c("Brec", "q")
  var_trans <- c("R0", "q")
  fun_trans <- c("/", "*")
  fun_fixed <- c("log", NA)
  rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)

  return(list(obj = obj, opt = opt, SD = SD, report = c(report, list(conv = !is.character(opt) && SD$pdHess))))
}

#' @rdname SRA_scope
#' @export
plot_SRA_scope <- function(OM, Chist, Index = matrix(NA, 0, 0), CAA = NA, CAL = NA, ML = NA, report_list, Year = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if(is.vector(Chist)) Chist <- matrix(Chist, ncol = 1)
  if(is.vector(Index)) Index <- matrix(Index, ncol = 1)
  if(is.matrix(CAA)) CAA <- array(CAA, c(dim(CAA), 1))
  if(is.matrix(CAL)) CAL <- array(CAL, c(dim(CAL), 1))
  if(is.vector(ML)) ML <- matrix(ML, ncol = 1)

  nsim <- OM@nsim
  maxage <- OM@maxage
  nyears <- OM@nyears
  if(is.null(Year)) {
    if(length(OM@EffYears) == nyears) Year <- OM@EffYears else Year <- 1:nyears
  }
  Year_matrix <- matrix(Year, ncol = nsim, nrow = nyears)
  Yearplusone_matrix <- matrix(c(Year, max(Year) + 1), ncol = nsim, nrow = nyears+1)
  nfleet <- ncol(Chist)
  nsurvey <- ncol(Index)
  length_bin <- report_list[[1]]$length_bin

  ###### First figure OM summary
  par(mfrow = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 3, 0))

  # R0 histogram
  hist(OM@cpars$R0, main = "", xlab = expression(R[0]))

  # Depletion histograms
  hist(OM@cpars$initD, main = "", xlab = "Initial depletion")
  hist(OM@cpars$D, main = "", xlab = "Depletion")

  # Perr
  Perr <- OM@cpars$Perr_y[, maxage:(maxage+nyears-1), drop = FALSE]
  matplot(Year_matrix, t(Perr), type = "l", col = "black", xlab = "Year", ylab = "Recruitment deviations", ylim = c(0, 1.1 * max(Perr)))
  abline(h = 0, col = "grey")

  # Find
  matplot(Year_matrix, t(OM@cpars$Find), type = "l", col = "black", xlab = "Year", ylab = "Apical F")
  abline(h = 0, col = "grey")

  # Selectivity
  if(nfleet == 1) {
    vul <- do.call(cbind, lapply(report_list, getElement, "vul"))
    matplot(matrix(length_bin, ncol = nsim, nrow = length(length_bin)), vul, type = "l", col = "black",
            xlab = "Length", ylab = "Selectivity", ylim = c(0, 1.1))
    abline(h = 0, col = "grey")
  } else {
    matplot(matrix(1:maxage, ncol = nsim, nrow = maxage), t(OM@cpars$V[, , nyears]), type = "l", col = "black",
            xlab = "Age", ylab = "Selectivity (last historical year)", ylim = c(0, 1.1))
    abline(h = 0, col = "grey")
  }

  mtext("Operating model parameter summary", outer = TRUE, side = 3)

  ###### Next ff figures by fleet
  for(ff in 1:nfleet) {
    par(mfrow = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 3, 0))
    # Selectivity
    vul_ff <- do.call(cbind, lapply(report_list, function(x) x$vul[, ff]))
    matplot(matrix(length_bin, ncol = nsim, nrow = length(length_bin)), vul_ff, type = "l", col = "black",
            xlab = "Length", ylab = paste("Selectivity of Fleet", ff))
    abline(h = 0, col = "grey")

    # Find
    FM <- do.call(cbind, lapply(report_list, function(x) x$F[, ff]))
    matplot(Year_matrix, FM, type = "l", col = "black", xlab = "Year", ylab = paste("Fishing Mortality of Fleet", ff))
    abline(h = 0, col = "grey")

    # Sampled catches
    Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ff]))
    matplot(Year_matrix, Cpred, type = "l", col = "black", xlab = "Year", ylab = paste("Catch of Fleet", ff))
    lines(Year, Chist[, ff], col = "red", lwd = 3)
    abline(h = 0, col = "grey")

    # ML fits
    MLpred <- do.call(cbind, lapply(report_list, function(x) x$mlen_pred[, ff]))
    matplot(Year_matrix, MLpred, type = "l", col = "black", xlab = "Year", ylab = "Mean length")
    if(!all(is.na(CAL))) {
      lines(Year, CAL[, , ff] %*% length_bin/rowSums(CAL[, , ff], na.rm = TRUE), col = "red", lwd = 3, typ = "o", pch = 16)
    } else if(!all(is.na(ML))) lines(Year, ML, col = "red", lwd = 3, typ = "o", pch = 16)

    # MA fits
    MApred <- do.call(cbind, lapply(report_list, function(x) x$CAApred[, , ff] %*% 1:maxage/x$CN[, ff]))
    matplot(Year_matrix, MApred, type = "l", col = "black", xlab = "Year", ylab = "Mean age")
    if(!all(is.na(CAA))) {
      lines(Year, CAA[,,ff] %*% c(1:ncol(CAA[,,ff]))/rowSums(CAA[,,ff], na.rm = TRUE), col = "red", lwd = 3, typ = "o", pch = 16)
    }

    mtext(paste0("Fleet ", ff, ": observed (red) and predicted data (black) \nfrom SRA"), outer = TRUE, side = 3)
  }

  if(nsurvey > 0) {
    if(nsurvey > 1) par(mfrow = c(2, min(ceiling(0.5 * nsurvey), 3))) else par(mfrow = c(1, 1))

    for(sur in 1:nsurvey) {
      Ipred <- do.call(cbind, lapply(report_list, function(x) x$Ipred[, sur]))
      matplot(Year_matrix, Ipred, type = "l", col = "black", xlab = "Year", ylab = paste("Index #", sur))
      lines(Year, Index[, sur], col = "red", lwd = 3, typ = "o", pch = 16)
      abline(h = 0, col = "grey")
    }
  } else {
    par(mfrow = c(1, 1))
  }

  E <- do.call(cbind, lapply(report_list, getElement, "E"))
  matplot(Yearplusone_matrix, E, type = "l", col = "black", xlab = "Year", ylab = "Spawning biomass")
  abline(h = 0, col = "grey")

  invisible()
}
