
summary_cDD <- function(Assessment) {
  assign_Assessment_slots()

  if(conv) current_status <- c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  input_parameters <- data.frame(Value = as.numeric(c(h, info$data$M, info$data$Kappa, info$data$k, info$data$wk, info$data$Winf)),
                                 Description = c("Stock-recruit steepness", "Natural mortality", "Weight coefficient",
                                                 "Age of knife-edge selectivity", "Weight at age k", "Asymptotic weight"),
                                 stringsAsFactors = FALSE)
  rownames(input_parameters) <- c("h", "M", "Kappa", "k", "w_k", "Winf")
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  if(conv) derived <- c(B0, N0, MSY, FMSY, BMSY)
  else derived <- rep(NA, 5)
  derived <- data.frame(Value = derived,
                        Description = c("Unfished biomass", "Unfished abundance", "Maximum sustainable yield (MSY)",
                                        "Fishing mortality at MSY", "Biomass at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("B0", "N0", "MSY", "FMSY", "BMSY")

  if(!is.character(SD)) {
    model_estimates <- summary(SD)
    model_estimates <- model_estimates[is.na(model_estimates[, 2]) || model_estimates[, 2] > 0, ]
  } else {
    model_estimates <- SD
  }

  output <- list(model = "Continuous Delay-Differential",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_cDD <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary(Assessment)
    html_report(plot.dir, model = "Continuous Delay Differential", current_status = index.report$current_status,
                input_parameters = index.report$input_parameters, derived_quantities = index.report$derived_quantities,
                model_estimates = index.report$model_estimates, name = Name, report_type = "Index")
  }

  age <- 1:info$LH$maxage

  plot_generic_at_age(age, info$LH$LAA, label = 'Mean Length-at-age')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "lifehistory_1_length_at_age.png"))
    plot_generic_at_age(age, info$LH$LAA, label = 'Mean Length-at-age')
    dev.off()
    lh.file.caption <- c("lifehistory_1_length_at_age.png",
                         paste("Mean Length-at-age from Data object."))
  }

  plot_generic_at_age(age, info$LH$WAA, label = 'Mean Weight-at-age')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "lifehistory_2_mean_weight_at_age.png"))
    plot_generic_at_age(age, info$LH$WAA, label = 'Mean Weight-at-age')
    dev.off()
    lh.file.caption <- rbind(lh.file.caption,
                             c("lifehistory_2_mean_weight_at_age.png",
                               "Mean Weight at age from Data object."))
  }

  plot(info$LH$LAA, info$LH$WAA, typ = 'o', xlab = 'Length', ylab = 'Weight')
  abline(h = 0, col = 'grey')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "lifehistory_3_length_weight.png"))
    plot(info$LH$LAA, info$LH$WAA, typ = 'o', xlab = 'Length', ylab = 'Weight')
    abline(h = 0, col = 'grey')
    dev.off()
    lh.file.caption <- rbind(lh.file.caption,
                             c("lifehistory_3_length_weight.png",
                               "Length-weight relationship from Data object."))
  }

  k <- info$data$k
  sel <- ifelse(age < k, 0, 1)
  plot_ogive(age, sel, label = "Maturity")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "lifehistory_4_maturity.png"))
    plot_ogive(age, sel, label = "Maturity")
    dev.off()
    lh.file.caption <- rbind(lh.file.caption,
                             c("lifehistory_4_maturity.png", "Assumed knife-edge maturity at the age corresponding to the length of 50% maturity."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "Continuous Delay Differential", captions = lh.file.caption,
                name = Name, report_type = "Life_History")
  }

  Year <- info$Year

  plot_timeseries(as.numeric(names(Obs_Catch)), Obs_Catch, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_catch.png"))
    plot_timeseries(as.numeric(names(Obs_Catch)), Obs_Catch, label = "Catch")
    dev.off()
    data.file.caption <- c("data_catch.png", "Catch time series")
  }

  plot_timeseries(as.numeric(names(Obs_Index)), Obs_Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(as.numeric(names(Obs_Index)), Obs_Index, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption, c("data_index.png", "Index time series."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "Continuous Delay Differential", captions = data.file.caption,
                name = Name, report_type = "Data")
  }

  if(conv) {
    ind <- names(SD$par.fixed) == "log_R0"
    plot_lognormalvar(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), label = expression(Unfished~~recruitment~~(R[0])), logtransform = TRUE)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_R0.png"))
      plot_lognormalvar(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), label = expression(Unfished~~recruitment~~(R[0])), logtransform = TRUE)
      dev.off()
      assess.file.caption <- c("assessment_R0.png", "Estimate of R0, distribution based on
                               normal approximation of estimated covariance matrix.")
    }

    if(!"transformed_h" %in% names(obj$env$map)) {
      ind <- names(SD$par.fixed) == "transformed_h"
      plot_steepness(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), is_transform = TRUE, SR = info$data$SR_type)
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "assessment_h.png"))
        plot_steepness(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), is_transform = TRUE, SR = info$data$SR_type)
        dev.off()
        assess.file.caption <- rbind(assess.file.caption,
                                     c("assessment_h.png", "Estimate of steepness, distribution based on normal
                                       approximation of estimated covariance matrix."))
      }
    }
  }

  plot_ogive(age, sel)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_selectivity.png"))
    plot_ogive(age, sel)
    dev.off()
    if(conv) assess.file.caption <- rbind(assess.file.caption,
                                          c("assessment_selectivity.png", "Assumed knife-edge selectivity at the age corresponding to the length of 50% maturity."))
    else assess.file.caption <- c("assessment_selectivity.png",
                                  "Assumed knife-edge selectivity at the age corresponding to the length of 50% maturity.")
  }

  plot_timeseries(Year, Obs_Index, Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index.png"))
    plot_timeseries(Year, Obs_Index, Index, label = "Index")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index.png", "Observed (black) and predicted (red) index."))
  }

  plot_residuals(Year, log(Obs_Index/Index), label = "log(Index) Residual")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index_residual.png"))
    plot_residuals(Year, log(Obs_Index/Index), label = "log(Index) Residual")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index_residual.png", "Index residuals in log-space."))
  }

  qqnorm(log(Obs_Index/Index), main = "")
  qqline(log(Obs_Index/Index))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index_qqplot.png"))
    qqnorm(log(Obs_Index/Index), main = "")
    qqline(log(Obs_Index/Index))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index_qqplot.png", "QQ-plot of index residuals in log-space."))
  }

  plot_timeseries(Year, Obs_Catch, Catch, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch.png"))
    plot_timeseries(Year, Obs_Catch, Catch, label = "Catch")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch.png", "Observed (black) and predicted (red) catch. Predicted catch should match observed in this model."))
  }

  ny <- info$data$ny
  SSB <- B[1:ny]
  Arec <- TMB_report$Arec
  Brec <- TMB_report$Brec
  if(info$data$SR_type == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB)
  if(info$data$SR_type == "Ricker") expectedR <- Arec * SSB * exp(-Brec * SSB)

  plot_SR(SSB, expectedR, R0, B0)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_stock_recruit.png"))
    plot_SR(SSB, expectedR, R0, B0)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_stock_recruit.png", "Stock-recruitment relationship."))
  }

  plot_timeseries(as.numeric(names(B)), B, label = "Biomass")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_biomass.png"))
    plot_timeseries(as.numeric(names(B)), B, label = "Biomass")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_biomass.png", "Time series of biomass."))
  }

  if(conv) {
    plot_timeseries(as.numeric(names(B_BMSY)), B_BMSY, label = expression(B/B[MSY]))
    abline(h = 1, lty = 2)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_B_BMSY.png"))
      plot_timeseries(as.numeric(names(B_BMSY)), B_BMSY, label = expression(B/B[MSY]))
      abline(h = 1, lty = 2)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_B_BMSY.png", "Time series of B/BMSY."))
    }
  }

  plot_timeseries(as.numeric(names(B_B0)), B_B0, label = expression(B/B[0]))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_B_B0.png"))
    plot_timeseries(as.numeric(names(B_B0)), B_B0, label = expression(B/B[0]))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_B_B0.png", "Time series of biomass depletion."))
  }

  plot_timeseries(as.numeric(names(R)), R, label = "Recruitment")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_recruitment.png"))
    plot_timeseries(as.numeric(names(R)), R, label = "Recruitment")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_recruitment.png", "Time series of recruitment."))
  }

  plot_timeseries(as.numeric(names(N)), N, label = "Population Abundance (N)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_abundance.png"))
    plot_timeseries(as.numeric(names(N)), N, label = "Population Abundance (N)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_abundance.png", "Time series of abundance."))
  }

  plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing Mortality (F)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_F.png"))
    plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing Mortality (F)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_F.png", "Time series of fishing mortality."))
  }

  if(conv) {
    plot_timeseries(as.numeric(names(F_FMSY)), F_FMSY, label = expression(F/F[MSY]))
    abline(h = 1, lty = 2)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_F_FMSY.png"))
      plot_timeseries(as.numeric(names(F_FMSY)), F_FMSY, label = expression(F/F[MSY]))
      abline(h = 1, lty = 2)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_F_FMSY.png", "Time series of F/FMSY."))
    }

    plot_Kobe(B_BMSY, F_FMSY, ylab = expression(F/F[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_Kobe.png"))
      plot_Kobe(B_BMSY, F_FMSY, ylab = expression(F/F[MSY]))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_Kobe.png", "Kobe plot trajectory of stock."))
    }

    plot_yield_cDD(info$data, TMB_report, FMSY, MSY, xaxis = "F")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_yield_curve_F.png"))
      plot_yield_cDD(info$data, TMB_report, FMSY, MSY, xaxis = "F")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_yield_curve_F.png", "Yield plot relative to F."))
    }

    plot_yield_cDD(info$data, TMB_report, FMSY, MSY, xaxis = "Depletion")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_yield_curve_B_B0.png"))
      plot_yield_cDD(info$data, TMB_report, FMSY, MSY, xaxis = "Depletion")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_yield_curve_B_B0.png", "Yield plot relative to depletion."))
    }
  }

  plot_surplus_production(B, B0, Obs_Catch)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_surplus_production.png"))
    plot_surplus_production(B, B0, Obs_Catch)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_surplus_production.png", "Surplus production relative to depletion."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "Continuous Delay Differential (State-Space)", captions = assess.file.caption,
                name = Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }

  return(invisible())
}


#' @importFrom reshape2 acast
profile_likelihood_cDD <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was not found. See help file.")
  if(!is.null(dots$R0)) R0 <- dots$R0 else {
    R0 <- Assessment@R0
    profile_par <- "h"
  }
  if(!is.null(dots$h)) h <- dots$h else {
    h <- Assessment@h
    profile_par <- "R0"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(R0 = R0, h = h)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$log_R0 <- log(profile_grid[i, 1] * Assessment@info$rescale)
    if(Assessment@info$data$SR_type == "BH") {
      params$transformed_h <- logit((profile_grid[i, 2] - 0.2)/0.8)
    } else {
      params$transformed_h <- log(profile_grid[i, 2] - 0.2)
    }

    if(length(Assessment@opt$par) == 1) { # R0 is the only estimated parameter
      if(!joint_profile && profile_par == "R0") {
        nll <- Assessment@obj$fn(params$log_R0)
      } else {
        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
        opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
        if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      }
    } else if(length(Assessment@opt$par) == 2 && all(names(Assessment@opt$par) == c("log_R0", "transformed_h"))) {
      if(joint_profile) {
        nll <- Assessment@obj$fn(c(params$log_R0, params$transformed_h))
      } else {
        if(profile_par == "R0") map$log_R0 <- factor(NA) else map$transformed_h <- factor(NA)
        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
        opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
        if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      }
    } else if(length(Assessment@opt$par) == 3 && all(names(Assessment@opt$par) == c("log_R0", "transformed_h", "F_equilibrium"))) {
      map$log_R0 <- map$transformed_h <- factor(NA)
      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    } else { #R0, F are estimated parameters
      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    }
    if(!exists("nll")) nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(figure) {
    R0.MLE <- Assessment@R0
    h.MLE <- Assessment@h
    if(joint_profile) {
      z.mat <- acast(profile_grid, list("h", "R0"), value.var = "nll")
      contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]), nlevels = 20)
      points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
    } else {
      if(profile_par == "R0") xlab <- expression(R[0]) else xlab <- "Steepness"
      plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
    }

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      if(joint_profile) {
        z.mat <- acast(profile_grid, list("h", "R0"), value.var = "nll")
        contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]), nlevels = 20)
        points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
        msg <- "Joint profile likelihood of R0 and h. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate."
      } else {
        if(profile_par == "R0") xlab <- expression(R[0]) else xlab <- "Steepness"
        plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
        msg <- paste0("Profile likelihood of ", profile_par, ". Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      }
      dev.off()

      profile.file.caption <- c("profile_likelihood.png", msg)

      html_report(plot.dir, model = "Continuous Delay Differential", captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))

    }

  }
  return(profile_grid)
}


#' @importFrom gplots rich.colors
retrospective_cDD <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()
  data <- info$data
  ny <- data$ny
  k <- data$k

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k
  Year <- c(Year, moreRecruitYears)
  C_hist <- data$C_hist
  I_hist <- data$I_hist

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, B/BMSY, B/B0, N, R, F, F/FMSY
  retro_ts <- array(NA, dim = c(nyr+1, ny + k, 8))
  retro_est <- array(NA, dim = c(nyr+1, dim(summary(SD))))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    ny_ret <- ny - i
    C_hist_ret <- C_hist[1:ny_ret]
    I_hist_ret <- I_hist[1:ny_ret]
    data$ny <- ny_ret
    data$C_hist <- C_hist_ret
    data$I_hist <- I_hist_ret

    obj2 <- MakeADFun(data = data, parameters = info$params, map = obj$env$map, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- get_MSY_cDD(info$data, report$Arec, report$Brec)
      report <- c(report, ref_pt)

      if(rescale != 1) {
        vars_div <- c("B0", "B", "Cpred", "BMSY", "MSY", "N0", "N", "R", "R0")
        vars_mult <- c("Brec")
        var_trans <- c("R0")
        fun_trans <- c("/")
        fun_fixed <- c("log")
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      B <- c(report$B, rep(NA, k - 1 + i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/B0
      R <- c(report$R, rep(NA, i))
      N <- c(report$N, rep(NA, k - 1 + i))
      FMort <- c(report$F, rep(NA, k + i))
      F_FMSY <- FMort/report$FMSY

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, R, N, FMort, F_FMSY)
      if(!is.character(SD)) retro_est[i+1, , ] <- summary(SD)

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }

  }

  fix_h <- "transformed_h" %in% names(obj$env$map)
  if(fix_h) est_ind <- 4 else est_ind <- 4:5
  est_lab <- "R0 estimate"
  if(!fix_h) est_lab <- c(est_lab, "Steepness estimate")

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1], retro_est[, , 1][, est_ind, drop = FALSE],
                                 ts_lab = c("Biomass", "B_BMSY", "B_B0", "Recruitment", "Abundance", "F", "F_FMSY"),
                                 est_lab = est_lab)

  if(figure) {
    plot_retro_cDD(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                   nyr_label = 0:nyr, color = rich.colors(nyr+1), fix_h, data$SR_type)
  }

  return(Mohn_rho)
}

plot_retro_cDD <- function(retro_ts, retro_est, save_figure = FALSE, save_dir = tempdir(), nyr_label, color, fix_h, SR) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]), "Recruitment",
                "Population Abundance (N)", "Fishing Mortality (F)", expression(F/F[MSY]))
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "cDD"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(retro_ts[, , i+1], na.rm = TRUE)
    plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i],
         ylim = c(0, 1.1 * y.max), col = color[1])
    for(j in 2:length(nyr_label)) {
      lines(Year, retro_ts[j, , i+1], col = color[j])
    }
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    abline(h = 0, col = 'grey')
    if(i %in% c(2, 7)) abline(h = 1, lty = 2)

    if(save_figure) {
      create_png(filename = file.path(plot.dir, paste0("retrospective_", i, ".png")))
      plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i],
           ylim = c(0, 1.1 * y.max), col = color[1])
      for(j in 2:length(nyr_label)) {
        lines(Year, retro_ts[j, , i+1], col = color[j])
      }
      legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
             title = "Years removed:")
      abline(h = 0, col = 'grey')
      if(i %in% c(2, 7)) abline(h = 1, lty = 2)
      dev.off()
    }
  }

  plot_lognormalvar(retro_est[, 1, 1], retro_est[, 1, 2], label = expression(hat(R)[0]),
                    logtransform = TRUE, color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
         title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 1, ".png")))
    plot_lognormalvar(retro_est[, 1, 1], retro_est[, 1, 2], label = expression(hat(R)[0]),
                      logtransform = TRUE, color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    dev.off()
  }

  if(!fix_h) {
    plot_steepness(retro_est[, 2, 1], retro_est[, 2, 2], is_transform = TRUE, SR = SR, color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 2, ".png")))
      plot_steepness(retro_est[, 2, 1], retro_est[, 2, 2], is_transform = TRUE, SR = SR, color = color)
      legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
             title = "Years removed:")
      dev.off()
    }
  }

  if(save_figure) {
    ret.file.caption <- data.frame(x1 = paste0("retrospective_", c(1:(n_tsplots+1)), ".png"),
                                   x2 = paste0("Retrospective pattern in ",
                                               c("biomass", "B/BMSY", "biomass depletion", "recruitment",
                                                 "abundance", "F", "F/FMSY", "R0 estimate"), "."),
                                   stringsAsFactors = FALSE)
    if(!fix_h) {
      ret.file.caption <- rbind(ret.file.caption,
                                c(paste0("retrospective_", n_tsplots+2, ".png"), "Retrospective pattern in steepness estimate."))
    }

    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Continuous Delay Differential", captions = ret.file.caption,
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}


plot_yield_cDD <- function(data, report, fmsy, msy, F.vector = seq(0, 2.5 * fmsy, length.out = 1e2), xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  Kappa <- data$Kappa
  M <- data$M
  Winf <- data$Winf
  wk <- data$wk
  SR_type <- data$SR_type

  Arec <- report$Arec
  Brec <- report$Brec
  BMSY <- report$BMSY

  Z <- F.vector + M
  BPR <- (wk + Kappa * Winf/Z)/(Z+Kappa)
  if(SR_type == "BH") R <- (Arec * BPR - 1)/Brec/BPR
  if(SR_type == "Ricker") R <- log(Arec * BPR)/Brec/BPR

  Biomass <- BPR * R
  Yield <- F.vector * Biomass
  ind <- R >= 0

  if(xaxis == "F") {
    plot(F.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality", ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Biomass", ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/report$B0, Yield[ind], typ = 'l', xlab = expression(B/B[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/report$B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/report$B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible()
}
