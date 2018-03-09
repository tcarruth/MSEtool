
summary_DD_SS <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(B_BMSY[length(B_BMSY)], U_UMSY[length(U_UMSY)]))
  rownames(current_status) <- c("B/BMSY", "U/UMSY")

  input_parameters <- data.frame(Value = c(unlist(info$data[c(2,3,4,6,7)]), info$sigma),
                                 Description = c("Unfished survival = exp(-M)", "alpha = Winf * (1-rho)",
                                                 "rho = (W_k+2 - Winf)/(W_k+1 - Winf)",
                                                 "Age of knife-edge selectivity",
                                                 "Weight at age k", "Catch observation error (log-space)"),
                                 stringsAsFactors = FALSE)
  rownames(input_parameters) <- c("S0", "alpha", "rho", "k", "w_k", "sigma")

  derived <- data.frame(Value = c(h, B0, R0, N0, BMSY, TMB_report$Spr0_DD, TMB_report$Spr_DD),
                        Description = c("Stock-recruit steepness", "Virgin biomass", "Virgin recruitment",
                                        "Virgin abundance", "Biomass at MSY", "Virgin biomass-per-recruit",
                                        "Biomass-per-recruit at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("h", "B0", "R0", "N0", "BMSY", "BPR0", "BPR_UMSY")

  model_estimates <- summary(SD)
  model_estimates <- model_estimates[model_estimates[, 2] > 0, ]

  output <- list(model = "Delay Difference (State-Space)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_DD_SS <- function(Assessment, save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary_DD_SS(Assessment)
    html_report(plot.dir, model = "Delay Difference (State-Space)",
                current_status = index.report$current_status,
                input_parameters = index.report$input_parameters,
                model_estimates = index.report$model_estimates,
                derived_quantities = index.report$derived_quantities,
                name = Data@Name, report_type = "Index")
  }

  lh.file.caption <- plot_life_history(Data, save_figure = save_figure, save_dir = save_dir, Model = Model)

  if(save_figure) {
    html_report(plot.dir, model = "Delay Difference (State-Space)",
                captions = lh.file.caption, name = Data@Name, report_type = "Life_History")
  }

  Year <- info$Year
  C_hist <- info$data$C_hist

  plot_timeseries(Year, C_hist, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_catch.png"))
    plot_timeseries(Year, C_hist, label = "Catch")
    dev.off()
    data.file.caption <- c("data_catch.png", "Catch time series")
  }

  if(!is.na(Data@CV_Cat[1])) {
    plot_timeseries(Year, C_hist, obs_CV = Data@CV_Cat, label = "Catch")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_catch_with_CI.png"))
      plot_timeseries(Year, C_hist, obs_CV = Data@CV_Cat, label = "Catch")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_catch_with_CI.png", "Catch time series with 95% confidence interval."))
    }
  }

  plot_timeseries(Year, info$I_hist, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(Year, info$I_hist, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_index.png", "Index time series."))
  }

  if(!is.na(Data@CV_Ind[1])) {
    plot_timeseries(Year, info$I_hist, obs_CV = Data@CV_Ind, label = "Index")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_index_with_CI.png"))
      plot_timeseries(Year, info$I_hist, obs_CV = Data@CV_Ind, label = "Index")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_index_with_CI.png", "Index time series with 95% confidence interval."))
    }
  }

  if(save_figure) {
    html_report(plot.dir, model = "Delay Difference (State-Space)",
                captions = data.file.caption, name = Data@Name, report_type = "Data")
  }

  logit.umsy <- as.numeric(obj$env$last.par.best[1])
  logit.umsy.sd <- sqrt(diag(SD$cov.fixed)[1])

  plot_betavar(logit.umsy, logit.umsy.sd, logit = TRUE, label = expression(hat(U)[MSY]))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_UMSYestimate.png"))
    plot_betavar(logit.umsy, logit.umsy.sd, logit = TRUE, label = expression(hat(U)[MSY]))
    dev.off()
    assess.file.caption <- c("assessment_UMSYestimate.png", "Estimate of UMSY, distribution based on normal approximation of estimated covariance matrix.")
  }

  log.msy <- as.numeric(obj$env$last.par.best[2])
  log.msy.sd <- sqrt(diag(SD$cov.fixed)[2])

  plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_MSYestimate.png"))
    plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_MSYestimate.png", "Estimate of MSY, distribution based on normal approximation of estimated covariance matrix."))
  }

  k_DD <- info$data$k_DD
  age <- 1:Data@MaxAge
  sel <- ifelse(age < k_DD, 0, 1)

  plot_ogive(age, sel)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_selectivity.png"))
    plot_ogive(age, sel)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_selectivity.png", "Assumed knife-edge selectivity at the age corresponding to the length of 50% maturity."))
  }

  plot_timeseries(Year, C_hist, Catch, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch.png"))
    plot_timeseries(Year, C_hist, Catch, label = "Catch")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch.png", "Observed (black) and predicted (red) catch."))
  }

  plot_residuals(Year, log(C_hist/Catch), label = "log(Catch) Residual")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch_residual.png"))
    plot_residuals(Year, log(C_hist/Catch), label = "log(Catch) Residual")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch_residual.png", "Catch residuals in log-space."))
  }

  qqnorm(log(C_hist/Catch), main = "")
  qqline(log(C_hist/Catch))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch_qqplot.png"))
    qqnorm(log(C_hist/Catch), main = "")
    qqline(log(C_hist/Catch))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch_qqplot.png", "QQ-plot of catch residuals in log-space."))
  }

  ny_DD <- info$data$ny_DD
  SSB <- B[1:ny_DD] - Catch # B*(1-u)
  Arec <- TMB_report$Arec_DD
  Brec <- TMB_report$Brec_DD
  expectedR <- Arec * SSB / (1 + Brec * SSB)

  first_recruit_year <- k_DD + 1
  last_recruit_year <- length(info$Year) + k_DD
  ind_recruit <- first_recruit_year:last_recruit_year
  rec_dev <- R[ind_recruit]

  plot_SR(SSB, expectedR, R0, B0, rec_dev)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_stock_recruit.png"))
    plot_SR(SSB, expectedR, R0, B0, rec_dev)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_stock_recruit.png", "Stock-recruitment relationship."))
  }

  if(max(rec_dev) > 3 * max(expectedR)) {
    y_zoom <- 3
    plot_SR(SSB, expectedR, R0, B0, rec_dev, y_zoom = y_zoom)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_stock_recruit_zoomed.png"))
      plot_SR(SSB, expectedR, R0, B0, rec_dev, y_zoom = y_zoom)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_stock_recruit_zoomed.png", "Stock-recruitment relationship (zoomed in)."))
    }
  } else y_zoom <- NULL

  plot_SR(SSB, expectedR, R0, B0, rec_dev, trajectory = TRUE, y_zoom = y_zoom)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_stock_recruit_trajectory.png"))
    plot_SR(SSB, expectedR, R0, B0, rec_dev, trajectory = TRUE, y_zoom = y_zoom)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_stock_recruit_trajectory.png", "Stock-recruitment relationship (trajectory plot)."))

  }

  plot_timeseries(as.numeric(names(B)), B, label = "Biomass")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_biomass.png"))
    plot_timeseries(as.numeric(names(B)), B, label = "Biomass")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_biomass.png", "Time series of biomass."))
  }

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

  plot_residuals(as.numeric(names(Random)), Random, label = "Recruitment deviations")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_rec_devs.png"))
    plot_residuals(as.numeric(names(Random)), Random, label = "Recruitment deviations")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_rec_devs.png", "Time series of recruitment deviations."))
  }

  plot_residuals(as.numeric(names(Random)), Random, Random_SE,
                 label = "Recruitment deviations")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_rec_devs_with_CI.png"))
    plot_residuals(as.numeric(names(Random)), Random, Random_SE,
                   label = "Recruitment deviations")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_rec_devs_with_CI.png", "Time series of recruitment deviations with 95% confidence intervals."))
  }

  plot_timeseries(as.numeric(names(N)), N, label = "Population Abundance (N)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_abundance.png"))
    plot_timeseries(as.numeric(names(N)), N, label = "Population Abundance (N)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_abundance.png", "Time series of abundance."))
  }

  plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_exploitation.png"))
    plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_exploitation.png", "Time series of exploitation rate."))
  }

  plot_timeseries(as.numeric(names(U_UMSY)), U_UMSY, label = expression(U/U[MSY]))
  abline(h = 1, lty = 2)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_U_UMSY.png"))
    plot_timeseries(as.numeric(names(U_UMSY)), U_UMSY, label = expression(U/U[MSY]))
    abline(h = 1, lty = 2)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_U_UMSY.png", "Time series of U/UMSY."))
  }

  plot_Kobe(B_BMSY, U_UMSY)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_Kobe.png"))
    plot_Kobe(B_BMSY, U_UMSY)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_Kobe.png", "Kobe plot trajectory of stock."))
  }

  plot_yield_DD(info$data, TMB_report, UMSY, MSY, xaxis = "U")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_U.png"))
    plot_yield_DD(info$data, TMB_report, UMSY, MSY, xaxis = "U")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_U.png", "Yield plot relative to exploitation."))
  }

  plot_yield_DD(info$data, TMB_report, UMSY, MSY, xaxis = "Depletion")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_B_B0.png"))
    plot_yield_DD(info$data, TMB_report, UMSY, MSY, xaxis = "Depletion")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_B_B0.png", "Yield plot relative to depletion."))
  }

  plot_surplus_production(B, B0, C_hist)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_surplus_production.png"))
    plot_surplus_production(B, B0, C_hist)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_surplus_production.png", "Surplus production relative to depletion."))
  }

	if(save_figure) {
	  html_report(plot.dir, model = "Delay Difference (State-Space)",
	              captions = assess.file.caption, name = Data@Name, report_type = "Assessment")
	  browseURL(file.path(plot.dir, "Assessment.html"))
  }
  return(invisible())
}


#' @importFrom reshape2 acast
profile_likelihood_DD_SS <- function(Assessment, figure = TRUE, save_figure = TRUE,
                                     save_dir = getwd(), ...) {
  dots <- list(...)
  if(!"UMSY" %in% names(dots)) stop("Sequence of UMSY was not found. See help file.")
  if(!"MSY" %in% names(dots)) stop("Sequence of MSY was not found. See help file.")
  UMSY.MLE <- Assessment@UMSY
  MSY.MLE <- Assessment@MSY

  profile.grid <- expand.grid(UMSY = UMSY, MSY = MSY)
  nll <- rep(NA, nrow(profile.grid))
  for(i in 1:nrow(profile.grid)) {
    params <- list(logit_UMSY_DD = log(profile.grid[i, 1]/(1-profile.grid[i, 1])),
                   log_MSY_DD = log(profile.grid[i, 2]),
                   log_q_DD = Assessment@obj$env$last.par.best[3],
                   log_sigma_DD = Assessment@info$params$log_sigma_DD,
                   log_tau_DD = Assessment@obj$env$last.par.best[4],
                   log_rec_dev = Assessment@info$params$log_rec_dev)
    obj <- MakeADFun(data = Assessment@info$data, parameters = params,
                     map = list(logit_UMSY_DD = factor(NA), log_MSY_DD = factor(NA),
                                log_sigma_DD = factor(NA)),
                     random = "log_rec_dev",
                     DLL = "MSEtool", silent = TRUE)
    opt <- suppressWarnings(try(nlminb(obj$par, obj$fn, obj$gr), silent = TRUE))
    if(!inherits(opt, "try-error") && opt$convergence == 0) {
      if(opt$convergence == 0) nll[i] <- opt$objective
    }
  }
  profile.grid$nll <- nll - min(nll, na.rm = TRUE)
  if(figure) {
    z.mat <- acast(profile.grid, UMSY ~ MSY, value.var = "nll")
    contour(x = UMSY, y = MSY, z = z.mat, xlab = expression(U[MSY]), ylab = "MSY",
            nlevels = 20)
    MLE <- as.numeric(Assessment@obj$env$last.par.best) # Max. likelihood est.
    UMSY.MLE <- 1/(1 + exp(-MLE[1]))
    MSY.MLE <- exp(MLE[2])
    points(UMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      contour(x = UMSY, y = MSY, z = z.mat, xlab = expression(U[MSY]), ylab = "MSY",
              nlevels = 20)
      UMSY.MLE <- 1/(1 + exp(-MLE[1]))
      MSY.MLE <- exp(MLE[2])
      points(UMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
      dev.off()
	    profile.file.caption <- c("profile_likelihood.png",
	                              "Joint profile likelihood of UMSY and MSY. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      html_report(plot.dir, model = "Delay Difference (State-Space)",
                  captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Data@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_DD_SS <- function(Assessment, nyr, figure = TRUE,
                                save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()
  data <- info$data
  ny_DD <- data$ny_DD
  k_DD <- data$k_DD

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k_DD
  Year <- c(Year, moreRecruitYears)
  C_hist <- data$C_hist
  E_hist <- data$E_hist
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, N, R, U, relU, relB, dep, log_rec_dev
  retro_ts <- array(NA, dim = c(nyr+1, ny_DD + k_DD, 9))
  summSD <- summary(SD)
  summSD <- summSD[rownames(summSD) != "log_rec_dev", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(summSD)))

  for(i in 0:nyr) {
    ny_DD_ret <- ny_DD - i
    C_hist_ret <- C_hist[1:ny_DD_ret]
    E_hist_ret <- E_hist[1:ny_DD_ret]
    data$ny_DD <- ny_DD_ret
    data$C_hist <- C_hist_ret
    data$E_hist <- E_hist_ret
    params$log_rec_dev <- rep(0, ny_DD_ret - k_DD)

    obj <- MakeADFun(data = data, parameters = params,
                     map = list(log_sigma_DD = factor(NA)),
                     random = "log_rec_dev", DLL = "MSEtool", silent = TRUE)
    opt <- suppressWarnings(try(nlminb(obj$par, obj$fn, obj$gr), silent = TRUE))

    if(!inherits(opt, "try-error") && opt$convergence == 0) {
      B <- c(obj$report()$B_DD, rep(NA, k_DD - 1 + i))
      relB <- c(obj$report()$relB_DD, rep(NA, k_DD - 1 + i))
      dep <- B/obj$report()$Bo_DD
      R <- c(obj$report()$R_DD, rep(NA, i))
      N <- c(obj$report()$N_DD, rep(NA, k_DD - 1 + i))
      U <- c(obj$report()$U_DD, rep(NA, k_DD + i))
      relU <- c(obj$report()$relU_DD, rep(NA, k_DD + i))
      log_rec_dev <- c(obj$report()$log_rec_dev, rep(NA, 2 * k_DD + i))

      retro_ts[i+1, , ] <- cbind(Year, B, relB, dep, R, N, U, relU, log_rec_dev)
      summSD <- summary(sdreport(obj))
      retro_est[i+1, , ] <- summSD[rownames(summSD) != "log_rec_dev", ]

    } else {
      warning(paste("Non-convergence when", i, "years of data were removed."))
    }
  }
  if(figure) {
    plot_retro_DD_SS(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                      nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }
  # Need to write legend
  legend <- NULL
  return(invisible(list(legend = legend, retro_ts = retro_ts, retro_est = retro_est)))
}


plot_retro_DD_SS <- function(retro_ts, retro_est, save_figure = FALSE,
                              save_dir = getwd(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]), "Recruitment",
                "Population Abundance (N)", "Exploitation rate (U)",
                expression(U/U[MSY]), "Recruitment deviations")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "DD_SS"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(retro_ts[, , i+1], na.rm = TRUE)
    if(i < n_tsplots) {
      ylim = c(0, 1.1 * y.max)
    } else ylim = c(-y.max, y.max)
    plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i],
         ylim = ylim, col = color[1])
    for(j in 2:length(nyr_label)) {
      lines(Year, retro_ts[j, , i+1], col = color[j])
    }
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    if(i != 8) abline(h = 0, col = 'grey')
    if(i %in% c(2, 7)) abline(h = 1, lty = 2)
    if(i == 8) abline(h = 0, lty = 2)

    if(save_figure) {
      create_png(filename = file.path(plot.dir, paste0("retrospective_", i, ".png")))
      plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i],
           ylim = ylim, col = color[1])
      for(j in 2:length(nyr_label)) {
        lines(Year, retro_ts[j, , i+1], col = color[j])
      }
      legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
             title = "Years removed:")
      if(i != 8) abline(h = 0, col = 'grey')
      if(i %in% c(2, 7)) abline(h = 1, lty = 2)
      if(i == 8) abline(h = 0, lty = 2)
      dev.off()
    }
  }

  plot_betavar(retro_est[, 1, 1], retro_est[, 1, 2], logit = TRUE,
               label = expression(hat(U)[MSY]), color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
         title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 1, ".png")))
    plot_betavar(retro_est[, 1, 1], retro_est[, 1, 2], logit = TRUE,
                 label = expression(hat(U)[MSY]), color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    dev.off()
  }

  plot_lognormalvar(retro_est[, 2, 1], retro_est[, 2, 2], logtransform = TRUE,
                    label = expression(widehat(MSY)), color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
         title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 2, ".png")))
    plot_lognormalvar(retro_est[, 2, 1], retro_est[, 2, 2], logtransform = TRUE,
                      label = expression(widehat(MSY)), color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
           title = "Years removed:")
    dev.off()
  }

  if(save_figure) {
    ret.file.caption <- data.frame(x1 = paste0("retrospective_", c(1:(n_tsplots+2)), ".png"),
                                   x2 = paste0("Retrospective pattern in ",
                                               c("biomass", "B/BMSY", "biomass depletion", "recruitment",
                                                 "abundance", "exploitation", "U/UMSY", "recruitment deviations", "UMSY estimate", "MSY estimate"), "."))
    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Delay Difference (State-Space)", captions = ret.file.caption,
                name = Assessment@Data@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}
