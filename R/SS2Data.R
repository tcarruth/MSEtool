
#' Reads data Stock Synthesis file structure into a Data object using package r4ss
#'
#' @description A function that uses the file location of a fitted SS3 model including input files to population
#' the various slots of an Data object.
#' @param SSdir A folder with Stock Synthesis input and output files in it
#' @param Name The name for the Data object
#' @param Common_Name Character string for the common name of the stock.
#' @param Species Scientific name of the species
#' @param Region Geographic region of the stock or fishery.
#' @param min_age_M Currently, the Data object supports a single value of M for all ages. The argument selects the
#' minimum age for calculating the mean of age-dependent M from the SS assessment.
#' @param comp_fleet A vector of indices corresponding to fleets in the assessment over which to aggregate the composition
#' (catch-at-length and catch-at-age) data. By default, characer string \code{"all"} will aggregate across all fleets.
#' @param comp_season Integer, for seasonal models, the season for which the value of the index will be used. By default, \code{"mean"}
#' will take the average across seasons.
#' @param comp_partition Integer vector for selecting length/age observations that are retained (2), discarded (1), or both (0). By default, \code{"all"}
#' sums over all available partitions.
#' @param comp_gender Integer vector for selecting length/age observations that are female (1), male (2), or both (0), or both scaled to sum to one (3).
#' By default, \code{"all"} sums over all gender codes.
#' @param index_fleet Integer for selecting the fleet of the index to put in the Data object. By default, \code{"SSB"}
#' will use the relative trend in spawning stock biomass as estimated in the model as the index.
#' @param index_season Integer, for seasonal models, the season for which the value of the index will be used. By default, \code{"mean"}
#' will take the average across seasons.
#' @param ... Arguments to pass to \link[r4ss]{SS_output}
#' @return An object of class Data.
#' @note Currently supports the latest version of r4ss on CRAN (v.1.24). Function may be incompatible with newer versions of r4ss on Github.
#' @author T. Carruthers
#' @export
#' @seealso \link{SS2OM}
#' @importFrom r4ss SS_output
SS2Data <- function(SSdir, Name = NULL, Common_Name = "", Species = "", Region = "",
                    min_age_M = 1, comp_fleet = "all", comp_season = "sum", comp_partition = "all", comp_gender = "all",
                    index_fleet = "SSB", index_season = "mean", ...) {

  dots <- list(dir = SSdir, ...)
  if(!any(names(dots) == "covar")) dots$covar <- FALSE
  if(!any(names(dots) == "forecast")) dots$forecast <- FALSE
  if(!any(names(dots) == "ncols")) dots$ncols <- 1e3
  if(!any(names(dots) == "printstats")) dots$printstats <- FALSE
  if(!any(names(dots) == "verbose")) dots$verbose <- FALSE
  if(!any(names(dots) == "warn")) dots$warn <- FALSE

  message(paste("-- Using function SS_output of package r4ss version", packageVersion("r4ss"), "to extract data from SS file structure --"))
  message(paste("Reading directory:", SSdir))
  replist <- do.call(SS_output, dots)
  message("-- End of r4ss operations --")

  season_as_years <- FALSE
  if(replist$nseasons == 1 && replist$seasduration < 1) {
    message(paste("Season-as-years detected in SS model. There is one season in the year with duration of", replist$seasduration, "year."))
    season_as_years <- TRUE
    nseas <- 1/replist$seasduration
    message("DLMtool operates on annual basis. Since the SS model is seasonal, we need to aggregate over seasons.")
  } else {
    nseas <- replist$nseasons
    if(nseas > 1) {
      message("DLMtool operating model is an annual model. Since the SS model is seasonal, we need to aggregate over seasons.")
    }
  }

  # Create Data object
  Data <- new("Data", stock = "MSE")
  Data@Common_Name <- Common_Name
  Data@Species <- Species
  Data@Region <- Region

  Data@MPs <- NA
  Data@TAC <- array(NA, dim = c(1, 1, 1))
  if(is.null(Name)) {
    Data@Name <- SSdir
  } else Data@Name <- Name
  Data@MPeff <- 1

  mainyrs <- replist$startyr:replist$endyr
  if(season_as_years) {
    nyears <- ceiling(length(mainyrs)/nseas)
    Data@Year <- 1:nyears

    seas1_yind_full <- expand.grid(nseas = 1:nseas, nyears = 1:nyears)
    seas1_yind <- which(seas1_yind_full$nseas == 1)
  } else {
    nyears <- length(mainyrs)
    Data@Year <- mainyrs
  }
  Data@LHYear <- Data@Year[length(Data@Year)]
  message(paste("Detected", nyears, "years in the assessment model."))
  message(paste0("First year: ", Data@Year[1], ", Last year: ", Data@Year[length(Data@Year)]))

  ##### Life history
  #### Growth --------------------------------------
  growdat <- getGpars(replist)      # Age-specific parameters in endyr

  # Max age
  Data@MaxAge <- maxage <- ceiling(nrow(growdat)/ifelse(season_as_years, nseas, 1))

  seas1_aind_full <- expand.grid(nseas = 1:nseas, age = 1:maxage)
  seas1_aind <- which(seas1_aind_full$nseas == 1)

  GP <- replist$Growth_Parameters   # Some growth parameters (presumably in endyr)
  if(nrow(GP)>1) {
    message(paste(nrow(GP),"different rows of growth parameters were reported by r4ss:"))
    print(GP)
    message("Only the first row of values will be used.\n")
  }

  #### Length at age --------------------------------------
  Len_age <- growdat$Len_Mid

  Data@vbLinf <- GP$Linf[1]
  t0 <- GP$A_a_L0[1]
  #t0[t0 > 1] <- 0
  Data@vbt0 <- t0
  muK <- GP$K[1]
  if(muK <= 0) { #Estimate K from Len_age if K < 0 (e.g., age-varying K with negative deviations in K).
    message("Negative K value was detected. Attempting to re-estimate K based on mean length-at-age...")
    get_K <- function(K, Lens, Linf, t0, ages) sum((Lens - (Linf * (1 - exp(-K * (ages - t0)))))^2)
    muK <- optimize(get_K, c(0, 2), Lens = Len_age, Linf = GP$Linf[1], t0 = t0, ages = 1:maxage)$minimum
  }
  Data@vbK <- muK

  message(paste0("Von Bertalanffy parameters: Linf = ", Data@vbLinf, ", K = ", muK, ", t0 = ", t0))

  LenCV <- GP$CVmax[1]
  if(LenCV > 1) LenCV <- LenCV/Data@vbLinf
  message(paste0("Data@LenCV = ", Data@LenCV))

  #### Weight
  Data@wla <- GP$WtLen1[1]
  Data@wlb <- GP$WtLen2[1]

  message(paste0("Length-weight parameters: a = ", Data@wla, ", b = ", Data@wlb))

  #### Maturity --------------------------------------
  if(min(growdat$Len_Mat < 1)) {                    # Condition to check for length-based maturity
    Mat <- growdat$Len_Mat/max(growdat$Len_Mat)
  } else {                                          # Use age-based maturity
    Mat <- growdat$Age_Mat/max(growdat$Age_Mat)
  }
  if(season_as_years) Mat <- Mat[seas1_aind]

  # Currently using linear interpolation of mat vs len, is robust and very close to true logistic model predictions
  Data@L50 <- LinInterp(Mat, Len_age, 0.5+1e-6)
  Data@L95 <- LinInterp(Mat, Len_age, 0.95)

  message(paste0("Lengths at 50% and 95% maturity: ", paste(Data@L50, Data@L95, collapse = " ")))

  #### M --------------------------------------
  M <- growdat$M

  if(length(unique(M)) > 1) {
    if(season_as_years) {
      seasonal_min_age <- min_age_M * nseas
      Data@Mort <- mean(M[growdat$Age >= seasonal_min_age])
    } else {
      Data@Mort <- mean(M[growdat$Age >= min_age_M])
    }
    message("Age-dependent natural mortality detected, but only a single value of M is currently supported for the Data object.")
    message(paste0("Using mean from ages >= ", min_age_M, " to set Data@Mort = ", Data@Mort, "."))
  } else {
    Data@Mort <- unique(M)
    message(paste0("Natural mortality Data@Mort = ", Data@Mort))
  }

  #### Composition data -------------------------
  ## Internal function -- don't move --
  get_comps <- function(dbase, comp_fleet, type = c("length", "age")) {
    type <- match.arg(type)
    dbase_ind <- match(dbase$Yr, mainyrs) # Match years
    dbase <- dbase[!is.na(dbase_ind), ]
    dbase$Obs2 <- dbase$Obs * dbase$N # Expand comp proportions to numbers

    comp_list <- split(dbase, dbase$Fleet) # List by fleet
    comp_mat <- lapply(comp_list, acast, formula = list("Yr", "Bin"), fun.aggregate = sum, value.var = "Obs2", fill = 0) # Convert to matrix

    comp_ind <- match(comp_fleet, as.numeric(names(comp_mat)))
    comp_mat <- comp_mat[comp_ind] # Subset fleets
    expand_matrix <- function(x) {
      if(type == "length") {
        ncol_dim <- length(replist$lbins)
        res <- matrix(NA, nrow = length(mainyrs), ncol = ncol_dim)
        res_ind <- match(as.numeric(rownames(x)), mainyrs)
        res_ind2 <- match(as.numeric(colnames(x)), replist$lbins)
      }
      if(type == "age") {
        ncol_dim <- maxage
        res <- matrix(NA, nrow = length(mainyrs), ncol = ncol_dim)
        res_ind <- match(as.numeric(rownames(x)), mainyrs)
        res_ind2 <- 1:ncol(x)
      }
      res[res_ind, res_ind2] <- x
      return(res)
    }
    comp_mat2 <- lapply(comp_mat, expand_matrix) # Expand matrices to full years (by mainyrs)
    comp_all <- do.call(rbind, comp_mat2)

    comp_res <- aggregate(comp_all, list(Yr = rep(1:length(mainyrs), length(comp_fleet))), sum, na.rm = TRUE) # Sum across fleets

    if(season_as_years) {
      if(is.numeric(comp_season) && comp_season <= nseas) {
        comp_res <- comp_res[seas1_yind_full$nseas == comp_season, ]
        message(paste("Using season", comp_season, "for", type, "comps."))
      }
      if(is.character(comp_season) && comp_season == "sum") {
        comp_res_list <- split(comp_res, seas1_yind_full$nyears)
        comp_res_list <- lapply(comp_res_list, colSums, na.rm = TRUE)
        comp_res <- do.call(rbind, comp_res_list)
      }
      message(paste(type, "comps summed across seasons."))
    }

    comp_res[comp_res == 0] <- NA
    return(comp_res[, -1])
  }

  #### CAA
  if(nrow(replist$agedbase) > 0) {
    if(is.character(comp_fleet) && comp_fleet == "all") {
      comp_fleet_age <- unique(replist$agedbase$Fleet)
    } else {
      comp_fleet_age <- comp_fleet
    }

    if(comp_partition != "all") {
      part_subset <- replist$agedbase$Part %in% comp_partition
    } else part_subset <- rep(TRUE, nrow(replist$agedbase))
    if(comp_gender != "all") {
      gender_subset <- replist$agedbase$Gender %in% comp_gender
    } else gender_subset <- rep(TRUE, nrow(replist$agedbase))

    agedbase <- replist$agedbase[part_subset & gender_subset, ]
    if(!season_as_years && nseas > 1) {
      if(is.numeric(comp_season) && comp_season <= nseas) {
        agedbase <- agedbase[agedbase$Seas == comp_season, ]
        message(paste("Using season", comp_season, "for age comps."))
      }
    }
    message(paste0("Using age comps partition codes: ", paste(unique(agedbase$Part), collapse = " ")))
    message(paste0("Using age comps gender codes: ", paste(unique(agedbase$Gender), collapse = " ")))

    if(any(!is.na(match(agedbase$Fleet, comp_fleet_age)))) {
      CAA <- get_comps(agedbase, comp_fleet = comp_fleet_age, type = "age")
      Data@CAA <- array(as.matrix(CAA), c(1, nyears, maxage))
      message(paste0("Collected age comps from Fleets: \n",
                     paste0(paste(comp_fleet_age, replist$FleetNames[comp_fleet_age], collapse = "\n"))))
    } else {
      message(paste0("Could not find any age comps from Fleets: \n",
                     paste0(paste(comp_fleet_age, replist$FleetNames[comp_fleet_age], collapse = "\n"))))
    }
  } else {
    message("No age comps found in SS assessment.")
  }

  #### CAL
  if(nrow(replist$lendbase) > 0) {
    if(is.character(comp_fleet) && comp_fleet == "all") {
      comp_fleet_length <- unique(replist$lendbase$Fleet)
    } else {
      comp_fleet_length <- comp_fleet
    }

    if(comp_partition != "all") {
      part_subset <- replist$lendbase$Part %in% comp_partition
    } else part_subset <- rep(TRUE, nrow(replist$lendbase))
    if(comp_gender != "all") {
      gender_subset <- replist$lendbase$Gender %in% comp_gender
    } else gender_subset <- rep(TRUE, nrow(replist$lendbase))

    lendbase <- replist$lendbase[part_subset & gender_subset, ]
    if(!season_as_years && nseas > 1) {
      if(is.numeric(comp_season) && comp_season <= nseas) {
        lendbase <- lendbase[lendbase$Seas == comp_season, ]
        message(paste("Using season", comp_season, "for length comps."))
      }
    }
    message(paste0("Using length comps partition codes: ", paste(unique(lendbase$Part), collapse = " ")))
    message(paste0("Using length comps gender codes: ", paste(unique(lendbase$Gender), collapse = " ")))

    if(any(!is.na(match(lendbase$Fleet, comp_fleet_length)))) {
      CAL <- get_comps(lendbase, comp_fleet_length, type = "length")
      CAL <- as.matrix(CAL)
      Data@CAL <- array(CAL, c(1, nyears, ncol(CAL)))

      message(paste0("Collected length comps from Fleets: \n",
                     paste0(paste(comp_fleet_length, replist$FleetNames[comp_fleet_length], collapse = "\n"))))

      CAL_bins <- replist$lbins # add one more length bin
      width_bin <- CAL_bins[length(CAL_bins)] - CAL_bins[length(CAL_bins)-1]
      plus_one <- CAL_bins[length(CAL_bins)] + width_bin

      Data@CAL_bins <- c(CAL_bins, plus_one)

      ML <- rowSums(CAL * rep(CAL_bins, each = nyears))/rowSums(CAL)
      ML[ML <= 0] <- NA
      Data@ML <- matrix(ML, nrow = 1)

      lcpos <- apply(CAL, 1, function(x) if(all(is.na(x))) return(NA) else which.max(x))
      Data@Lc <- matrix(CAL_bins[lcpos], nrow = 1)
      Data@Lc[Data@Lc <= 0] <- NA

      Data@Lbar <- matrix(NA, nrow = 1, ncol = ncol(Data@Lc))
      for(i in 1:ncol(Data@Lbar)) {
        if(!is.na(lcpos[i])) {
          Data@Lbar[1, i] <- weighted.mean(x = CAL_bins[lcpos[i]:length(replist$lbins)],
                                           w = CAL[i, lcpos[i]:length(replist$lbins)])
        }
      }

    } else {
      message(paste0("Could not find any length comps from Fleets: \n",
                     paste0(paste(comp_fleet_length, replist$FleetNames[comp_fleet_length], collapse = "\n"))))
    }
  } else {
    message("No length comps found in SS assessment.")
  }

  #### Catch -------------------------
  message("Adding total catch in weight across all fleets...")
  cat_yr_ind <- !is.na(match(replist$timeseries$Yr, mainyrs))
  ts <- replist$timeseries[cat_yr_ind, ]

  cat_col <- grepl("obs_cat", colnames(ts))
  cat <- ts[, cat_col, drop = FALSE]

  is_weight <- replist$catch_units[replist$IsFishFleet] == 1

  cat_weight <- cat[, is_weight]
  cat_numbers <- cat[, !is_weight]
  if(ncol(cat_numbers) > 0) {
    fleet_in_numbers <- which(replist$catch_units == 2)
    message(paste0("Catch in numbers was detected for Fleets: \n", paste(paste(fleet_in_numbers, replist$FleetNames[fleet_in_numbers]), collapse = "\n")))
    message("Will use estimated assessment output for catch in weight (dead biomass) for these Fleets.")

    cat_col2 <- grepl("dead\\(B", colnames(ts))
    cat2 <- ts[, cat_col2]
    cat_numbers <- cat2[, !is_weight]
  }
  total_catch <- aggregate(rowSums(cbind(cat_weight, cat_numbers)), list(Yr = ts$Yr), sum, na.rm = TRUE)

  tc_ind <- match(total_catch$Yr, mainyrs)
  total_catch_vec <- total_catch$x[tc_ind]
  if(season_as_years) {
    total_catch2 <- aggregate(total_catch_vec, list(Yr = seas1_yind_full$nyears), sum, na.rm = TRUE)
    total_catch_vec <- total_catch2$x
    message("Summing catch across seasons.")
  }
  Data@Cat <- matrix(total_catch_vec, nrow = 1)
  message(paste0(sum(!is.na(Data@Cat[1, ])), " years of catch in Data object."))
  message(paste0("Default CV of Catch, Data@CV_Cat = ", Data@CV_Cat))

  Data@AvC <- mean(total_catch_vec)
  message(paste0("Mean catch, Data@AvC = ", round(Data@AvC, 2)))

  #### Index -------------------------
  # (function argument to select index)
  if(is.numeric(index_fleet) && index_fleet %in% unique(replist$cpue$Fleet)) {
    cpue <- split(replist$cpue, replist$cpue$Fleet)
    cpue_ind <- match(index_fleet, as.numeric(names(cpue)))

    cpue_subset <- cpue[[cpue_ind]]

    if(!season_as_years && nseas > 1 && !is.null(cpue_subset)) {
      if(is.numeric(index_season) && index_season <= nseas) {
        cpue_subset <- cpue_subset[cpue_subset$Seas == index_season, ]
        message(paste("Using season", index_season, "for index."))
      }
      if(is.character(index_season) && index_season == "mean") {
        cpue_subset <- aggregate(cpue_subset$Obs, by = list(Yr = cpue_subset$Yr), mean, na.rm = TRUE)
        names(cpue_subset) <- c("Yr", "Obs")
        message("Taking mean of index across seasons.")
      }
    }

    cpue_ind2 <- match(mainyrs, cpue_subset$Yr)
    Ind <- cpue_subset$Obs[cpue_ind2]

    Data@CV_Ind <- mean(cpue_subset$SE/cpue_subset$Obs)
    message(paste0("CV of Index (as mean of the ratios of SE and observed values), Data@CV_Ind = ", round(Data@CV_Ind, 2)))
  }
  if(is.character(index_fleet) && index_fleet == "SSB") {
    SR_ind <- match(mainyrs, replist$recruit$year)
    SSB <- replist$recruit$spawn_bio[SR_ind]
    Ind <- SSB/mean(SSB)

    message(paste0("Default CV of Index, Data@CV_Ind = ", round(Data@CV_Ind, 2)))
  }

  if(exists("Ind")) {
    if(season_as_years) {
      if(is.character(index_season) && index_season == "mean") {
        Ind2 <- aggregate(Ind, by = list(Yr = seas1_yind_full$nyears), mean, na.rm = TRUE)$x
        Ind2[is.nan(Ind2)] <- NA
        message("Taking mean of index across seasons.")
      }
      if(is.numeric(index_season) && index_season <= nseas) {
        Ind2 <- Ind[seas1_yind_full$nseas == index_season]
        message(paste("Using season", index_season, "for index."))
      }
      Ind <- Ind2
    }

    if(is.na(Ind[length(Ind)])) warning("No index value for most recent year.")
    message(paste0(sum(!is.na(Ind)), " years of index values from Fleet ", index_fleet, " (", replist$FleetNames[index_fleet], ") in Data object."))
    Data@Ind <- matrix(Ind, nrow = 1)

  } else {
    message(paste("No index found for Fleet", index_fleet, replist$FleetNames[index_fleet]))
    Data@Ind <- matrix(NA, ncol = nyears, nrow = 1)
  }


  #### Recruitment
  rec_ind <- match(mainyrs, replist$recruit$year)
  rec <- replist$recruit$pred_recr[rec_ind]

  if(season_as_years) {
    rec2 <- aggregate(rec, by = list(Yr = seas1_yind_full$nyears), mean, na.rm = TRUE)$x
    rec <- rec2
    message("Summing recruitment across seasons.")
  }

  Data@Rec <- matrix(rec/mean(rec), nrow = 1)
  message("Relative recruitment strength to Data@Rec obtained from assessment.")

  #### Depletion
  SSB <- replist$recruit$spawn_bio[rec_ind]
  Data@Dt <- SSB[length(SSB)]/SSB[1]
  message(paste("Depletion since year", mainyrs[1], "(Data@Dt) =", round(Data@Dt, 2)))

  Data@Dep <- replist$current_depletion
  message(paste("Depletion from unfished conditions (Data@Dep) =", round(Data@Dep, 2)))

  Data@t <- length(Data@Year)

  #### Reference points ----------------------
  Data@Cref <- replist$derived_quants$Value[replist$derived_quants$LABEL == "TotYield_MSY"] * ifelse(season_as_years, nseas, 1)
  message(paste("Reference catch set to MSY, Data@Cref =", Data@Cref))

  FMSY <- replist$derived_quants$Value[replist$derived_quants$LABEL == "Fstd_MSY"] * ifelse(season_as_years, nseas, 1)
  Data@FMSY_M <- FMSY/Data@Mort
  message(paste0("FMSY = ", FMSY, ", Data@FMSY_M = ", Data@FMSY_M))

  Data@Bref <- replist$derived_quants$Value[replist$derived_quants$LABEL == "SSB_MSY"]
  message(paste("Reference biomass set to SSB at MSY, Data@Bref =", Data@Bref))

  SSB0 <- replist$derived_quants$Value[replist$derived_quants$LABEL == "SSB_Unfished"]
  Data@BMSY_B0 <- Data@Bref/SSB0
  message(paste("Data@BMSY_B0 =", Data@BMSY_B0))

  if(index_fleet == "SSB") {
    Data@Iref <- Data@Ind[1] * Data@Bref/SSB[1]
    message(paste("Data@Iref =", Data@Iref))
  }

  Data@Units <- ""
  OFLs <- replist$derived_quants[grepl("OFLCatch", replist$derived_quants$LABEL), ]
  if(season_as_years) {
    OFL_terminal <- sum(OFLs$Value[1:nseas])
  } else OFL_terminal <- OFLs$Value[1]

  if(is.na(OFL_terminal)) {
    Data@Ref <- Data@Cat[1, ncol(Data@Cat)]
    Data@Ref_type <- paste("Catch in Year", Data@Year[length(Data@Year)])
    message("No OFL detected from SS Assessment.")
  } else {
    Data@Ref <- OFL_terminal
    if(season_as_years) {
      Yr_OFL <- vapply(OFLs$LABEL[1:nseas], function(x) strsplit(x, "_")[[1]][2], character(1))
    } else {
      Yr_OFL <- strsplit(OFLs$LABEL[1], "_")[[1]][2]
    }
    if(length(Yr_OFL) == 1){
      Data@Ref_type <- paste("OFL in Year", Yr_OFL, "from SS Assessment")
    } else {
      Data@Ref_type <- paste("OFL in Years", Yr_OFL[1], "-", Yr_OFL[length(Yr_OFL)] , "from SS Assessment")
    }
  }

  message(paste0("Data@Ref = ", Data@Ref, " (", Data@Ref_type, ")"))

  #### Steepness --------------------------------------
  steep <- replist$parameters[grepl("steep", rownames(replist$parameters)), ]
  if(nrow(steep) == 1) {
    Data@steep <- steep$Value
  } else {
    SR_ind <- match(mainyrs, replist$recruit$year)
    SSB <- replist$recruit$spawn_bio[SR_ind]
    rec <- replist$recruit$pred_recr[SR_ind]

    res <- try(SSB0 <- replist$derived_quants[replist$derived_quants$LABEL == "SPB_Virgin", 2], silent = TRUE)
    if(inherits(res, "try-error")) SSB0 <- SSB[1]

    res2 <- try(R0 <- replist$derived_quants[replist$derived_quants$LABEL == "Recr_Virgin", 2], silent = TRUE)
    if(inherits(res, "try-error")) {
      surv <- c(1, exp(-cumsum(M[1:(maxage-1)])))
      Wt_age <- growdat$Wt_Mid
      if(season_as_years) Wt_age <- Wt_age[seas1_aind]
      SpR0 <- sum(Wt_age * Mat * surv)
      R0 <- SSB0/SpR0
    } else {
      SpR0 <- SSB0/(R0 * ifelse(season_as_years, nseas, 1))
    }
    if(replist$SRRtype == 3 || replist$SRRtype == 6) SR <- "BH"
    if(replist$SRRtype == 2) SR <- "Ricker"
    Data@steep <- mean(SRopt(100, SSB, rec, SpR0, plot = FALSE, type = SR), na.rm = TRUE)
  }
  message(paste0("Steepness = ", Data@steep))

  SpAbun_ind <- match(replist$endyr+1, replist$recruit$year)
  Data@SpAbun <- replist$recruit$spawn_bio[SpAbun_ind]

  #### Selectivity
  # Get F-at-age in terminal year, then obtain LFC and LFS
  ages <- growdat$Age
  cols <- match(ages, names(replist$Z_at_age))
  rows <- match(mainyrs, replist$Z_at_age$Year)

  Z_at_age <- replist$Z_at_age[rows, ]
  M_at_age <- replist$M_at_age[rows, ]

  rows2 <- Z_at_age$Gender == 1 & Z_at_age$Bio_Pattern == 1
  F_at_age <- t(Z_at_age[rows2, cols] - M_at_age[rows2, cols])
  F_at_age[nrow(F_at_age), ] <- F_at_age[nrow(F_at_age) - 1, ] # assume F at maxage = F at maxage-1

  if(ncol(F_at_age) == nyears - 1) { # Typically because forecast is off
    F_at_age_terminal <- F_at_age[, ncol(F_at_age)]
    F_at_age <- cbind(F_at_age, F_at_age_terminal)
  }
  if(ncol(F_at_age) == nyears && all(is.na(F_at_age[, ncol(F_at_age)]))) {
    F_at_age[, ncol(F_at_age)] <- F_at_age[, ncol(F_at_age)-1]
  }

  F_at_age[F_at_age < 1e-8] <- 1e-8

  if(season_as_years) {
    Ftab <- expand.grid(Age = 1:dim(F_at_age)[1], Yr = 1:dim(F_at_age)[2])
    Ftab$F_at_age <- as.vector(F_at_age)

    # Mean F across aggregated age (groups of nseas), then sum F across aggregated years (groups of nseas)
    sumF <- aggregate(Ftab[, 3], list(Age = seas1_aind_full[1:nrow(F_at_age), 2][Ftab[, 1]], Yr = Ftab[, 2]), mean, na.rm = TRUE)
    sumF <- aggregate(sumF[, 3], list(Age = sumF[, 1], Yr = rep(seas1_yind_full[1:ncol(F_at_age), 2], each = maxage)), sum, na.rm = TRUE)

    F_at_age <- matrix(sumF[, 3], nrow = maxage)
  }

  V_terminal <- F_at_age[, ncol(F_at_age)]/max(F_at_age[, ncol(F_at_age)])
  Data@LFC <- LinInterp(V_terminal, Len_age, 0.05, ascending = TRUE, zeroint = TRUE)
  Data@LFS <- Len_age[which.min((exp(V_terminal)-exp(1.05))^2 * 1:length(V_terminal))]
  message(paste0("Data@LFC = ", Data@LFC, ", Data@LFS = ", Data@LFS))

  Data@Log <- Data@Misc <- list(note = paste("This Data object was created by the SS2Data function from MSEtool version", packageVersion("MSEtool")))

  message("\nImport was successful.\n")

  Data

}

