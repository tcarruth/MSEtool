
summary_SCA3 <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  M <- info$data$M[1]
  maxage <- info$data$max_age
  Linf <- info$LH$Linf
  K <- info$LH$K
  t0 <- info$LH$t0
  A50 <- info$LH$A50
  A95 <- info$LH$A95
  Winf <- info$LH$a * Linf ^ info$LH$b

  Value <- c(h, M, maxage, Linf, K, t0, Winf, A50, A95)
  Description = c("Stock-recruit steepness", "Natural mortality", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  rownam <- c("h", "M", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  Value = c(VB0, SSB0, MSY, UMSY, VBMSY, SSBMSY)
  Description = c("Virgin vulnerable biomass",
                  "Virgin spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Harvest Rate at MSY",
                  "Vulnerable biomass at MSY", "SSB at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "UMSY", "VBMSY", "SSBMSY")

  model_estimates <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev", ]
  model_estimates <- model_estimates[model_estimates[, 2] > 0, ]
  dev_estimates <- cbind(Dev, SE_Dev)
  rownames(dev_estimates) <- paste0("log_rec_dev_", names(Dev))

  output <- list(model = "Statistical Catch-at-Age (SCA3)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = rbind(model_estimates, dev_estimates))
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_SCA3 <- function(Assessment, save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary(Assessment)
    html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)",
                current_status = index.report$current_status,
                input_parameters = index.report$input_parameters,
                model_estimates = index.report$model_estimates,
                derived_quantities = index.report$derived_quantities,
                name = Data@Name, report_type = "Index")
  }

  age <- 1:info$data$max_age

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

  plot_ogive(age, info$data$mat, label = "Maturity")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "lifehistory_4_maturity.png"))
    plot_ogive(age, info$data$mat, label = "Maturity")
    dev.off()
    lh.file.caption <- rbind(lh.file.caption,
                             c("lifehistory_4_maturity.png", "Maturity at age."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "Statistical Catch-at-Age",
                captions = lh.file.caption, name = Data@Name, report_type = "Life_History")
  }

  Year <- info$Year

  plot_timeseries(Year, Obs_Catch, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_catch.png"))
    plot_timeseries(Year, Obs_Catch, label = "Catch")
    dev.off()
    data.file.caption <- c("data_catch.png", "Catch time series")
  }

  if(!is.na(Data@CV_Cat[1]) && sdconv(1, Data@CV_Cat[1]) > 0.01) {
    plot_timeseries(Year, Obs_Catch, obs_CV = Data@CV_Cat[1], label = "Catch")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_catch_with_CI.png"))
      plot_timeseries(Year, Obs_Catch, obs_CV = Data@CV_Cat[1], label = "Catch")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_catch_with_CI.png", "Catch time series with 95% confidence interval."))
    }
  }

  plot_timeseries(Year, Obs_Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(Year, Obs_Index, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_index.png", "Index time series."))
  }

  if(!is.na(Data@CV_Ind[1]) && sdconv(1, Data@CV_Ind[1]) > 0.01) {
    plot_timeseries(Year, Obs_Index, obs_CV = Data@CV_Ind[1], label = "Index")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_index_with_CI.png"))
      plot_timeseries(Year, Obs_Index, obs_CV = Data@CV_Ind[1], label = "Index")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_index_with_CI.png", "Index time series with 95% confidence interval."))
    }
  }

  ind_valid <- rowSums(Obs_C_at_age, na.rm = TRUE) > 0
  Year2 <- Year[ind_valid]
  Obs_CAA <- Obs_C_at_age[ind_valid, ]
  plot_composition(Year2, Obs_CAA, plot_type = 'bubble_data', data_type = 'age')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_age_comps_bubble.png"))
    plot_composition(Year2, Obs_CAA, plot_type = 'bubble_data', data_type = 'age')
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_age_comps_bubble.png", "Age composition bubble plot."))
  }

  plot_composition(Year2, Obs_CAA, plot_type = 'annual', data_type = 'age')
  if(save_figure) {
    nplots <- ceiling(length(Year2)/16)
    for(i in 1:nplots) {
      ind <- (16*(i-1)+1):(16*i)
      if(i == nplots) ind <- (16*(i-1)+1):length(Year2)

      create_png(filename = file.path(plot.dir, paste0("data_age_comps_", i, ".png")))
      plot_composition(Year2, Obs_CAA, plot_type = 'annual', data_type = 'age', ind = ind)
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c(paste0("data_age_comps_", i, ".png"), paste0("Annual age compositions (", i, "/", nplots, ")")))
    }
  }

  if(save_figure) {
    html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)",
                captions = data.file.caption, name = Data@Name, report_type = "Data")
  }

  ind <- names(SD$par.fixed) == "log_R0"
  plot_lognormalvar(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), label = expression(Virgin~~recruitment~~(R[0])), logtransform = TRUE)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_R0.png"))
    plot_lognormalvar(SD$par.fixed[ind], sqrt(diag(SD$cov.fixed)[ind]), label = expression(Virgin~~recruitment~~(R[0])), logtransform = TRUE)
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

  plot_ogive(age, Selectivity[nrow(Selectivity), ])
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_selectivity.png"))
    plot_ogive(age, Selectivity[nrow(Selectivity), ])
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_selectivity.png", "Estimated selectivity at age."))
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


  Fit_CAA <- C_at_age[ind_valid, ]
  plot_composition(Year2, Obs_CAA, Fit_CAA, plot_type = 'bubble_residuals', data_type = 'age', bubble_adj = 35)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assess_age_comps_bubble_resids.png"))
    plot_composition(Year2, Obs_CAA, Fit_CAA, plot_type = 'bubble_residuals', data_type = 'age', bubble_adj = 35)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                               c("assess_age_comps_bubble_resids.png", "Age composition bubble plot of residuals (black are negative, white are positive)."))
  }

  plot_composition(Year2, Obs_CAA, Fit_CAA, N = info$data$CAA_n[ind_valid], plot_type = 'annual', data_type = 'age')
  if(save_figure) {
    nplots <- ceiling(length(Year2)/16)
    for(i in 1:nplots) {
      ind <- (16*(i-1)+1):(16*i)
      if(i == nplots) ind <- (16*(i-1)+1):length(Year2)

      create_png(filename = file.path(plot.dir, paste0("assess_age_comps_", i, ".png")))
      plot_composition(Year2, Obs_CAA, Fit_CAA, plot_type = 'annual', data_type = 'age', N = info$data$CAA_n[ind_valid], ind = ind)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c(paste0("assess_age_comps_", i, ".png"), paste0("Annual observed (black) and predicted (red) age compositions (",
                                                                                    i, "/", nplots, ")")))
    }
  }

  plot_composition(Year, Obs_C_at_age, C_at_age, plot_type = 'mean', data_type = 'age')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_mean_age.png"))
    plot_composition(Year, Obs_C_at_age, C_at_age, plot_type = 'mean', data_type = 'age')
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_mean_age.png", "Observed (black) and predicted (red) mean age of the composition data."))
  }

  Arec <- TMB_report$Arec
  Brec <- TMB_report$Brec
  SSB_plot <- SSB[1:(length(SSB)-1)]
  if(info$data$SR == "BH") expectedR <- Arec * SSB_plot / (1 + Brec * SSB_plot)
  if(info$data$SR == "Ricker") expectedR <- Arec * SSB_plot * exp(-Brec * SSB_plot)
  estR <- R[as.numeric(names(R)) > Year[1]]

  plot_SR(SSB_plot, expectedR, R0, SSB0, estR)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_stock_recruit.png"))
    plot_SR(SSB_plot, expectedR, R0, SSB0, estR)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_stock_recruit.png", "Stock-recruitment relationship."))
  }

  if(max(estR) > 3 * max(expectedR)) {
    y_zoom <- 3
    plot_SR(SSB_plot, expectedR, R0, SSB0, estR, y_zoom = y_zoom)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_stock_recruit_zoomed.png"))
      plot_SR(SSB_plot, expectedR, R0, SSB0, estR, y_zoom = y_zoom)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_stock_recruit_zoomed.png", "Stock-recruitment relationship (zoomed in)."))
    }
  } else y_zoom <- NULL

  plot_SR(SSB_plot, expectedR, R0, SSB0, estR, trajectory = TRUE, y_zoom = y_zoom)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_stock_recruit_trajectory.png"))
    plot_SR(SSB_plot, expectedR, R0, SSB0, estR, trajectory = TRUE, y_zoom = y_zoom)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_stock_recruit_trajectory.png", "Stock-recruitment relationship (trajectory plot)."))

  }

  plot_timeseries(as.numeric(names(SSB)), SSB, label = "Spawning Stock Biomass (SSB)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_spawning_biomass.png"))
    plot_timeseries(as.numeric(names(SSB)), SSB, label = "Spawning Stock Biomass (SSB)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_spawning_biomass.png", "Time series of spawning stock biomass."))
  }

  plot_timeseries(as.numeric(names(SSB_SSBMSY)), SSB_SSBMSY, label = expression(SSB/SSB[MSY]))
  abline(h = 1, lty = 2)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_SSB_SSBMSY.png"))
    plot_timeseries(as.numeric(names(SSB_SSBMSY)), SSB_SSBMSY, label = expression(SSB/SSB[MSY]))
    abline(h = 1, lty = 2)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_SSB_SSBMSY.png", "Time series of SSB/SSBMSY."))
  }

  plot_timeseries(as.numeric(names(SSB_SSB0)), SSB_SSB0, label = expression(SSB/SSB[0]))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_SSB_SSB0.png"))
    plot_timeseries(as.numeric(names(SSB_SSB0)), SSB_SSB0, label = expression(SSB/SSB[0]))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_SSB_SSB0.png", "Time series of spawning stock biomass depletion."))
  }

  plot_timeseries(as.numeric(names(R)), R, obs_ind_blue = as.numeric(names(R)) < Year[1], label = "Recruitment")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_recruitment.png"))
    plot_timeseries(as.numeric(names(R)), R, obs_ind_blue = as.numeric(names(R)) < Year[1], label = "Recruitment")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_recruitment.png", "Time series of recruitment (recruitment prior to the first year of the model in blue)."))
  }

  if(as.numeric(names(Dev))[1] < Year[1]) {
    obs_ind_blue <- as.numeric(names(Dev)) < Year[1]
    msg <- c("Time series of recruitment deviations  (deviations prior to the first year of the model in blue).",
             "Time series of recruitment deviations with 95% confidence intervals (deviations prior to the first year of the model in blue).")

  } else {
    obs_ind_blue <- NULL
    msg <- c("Time series of recruitment deviations.",
             "Time series of recruitment deviations with 95% confidence intervals.")
  }
  plot_residuals(as.numeric(names(Dev)), Dev, res_ind_blue = obs_ind_blue, label = Dev_type)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_rec_devs.png"))
    plot_residuals(as.numeric(names(Dev)), Dev, res_ind_blue = obs_ind_blue, label = Dev_type)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_rec_devs.png", msg[1]))
  }

  plot_residuals(as.numeric(names(Dev)), Dev, SE_Dev, res_ind_blue = obs_ind_blue, label = Dev_type)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_rec_devs_with_CI.png"))
    plot_residuals(as.numeric(names(Dev)), Dev, SE_Dev, res_ind_blue = obs_ind_blue, label = Dev_type)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_rec_devs_with_CI.png", msg[2]))
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

  plot_yield_SCA(info$data, TMB_report, UMSY, MSY, xaxis = "U", SR = info$data$SR_type)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_U.png"))
    plot_yield_SCA(info$data, TMB_report, UMSY, MSY, xaxis = "U", SR = info$data$SR_type)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_U.png", "Yield plot relative to exploitation."))
  }

  plot_yield_SCA(info$data, TMB_report, UMSY, MSY, xaxis = "Depletion", SR = info$data$SR_type)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_SSB_SSB0.png"))
    plot_yield_SCA(info$data, TMB_report, UMSY, MSY, xaxis = "Depletion", SR = info$data$SR_type)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_SSB_SSB0.png", "Yield plot relative to spawning depletion."))
  }

  plot_surplus_production(B, B0, Obs_Catch)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_surplus_production.png"))
    plot_surplus_production(B, B0, Obs_Catch)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_surplus_production.png", "Surplus production relative to depletion (total biomass)."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)",
                captions = assess.file.caption, name = Data@Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }
  return(invisible())
  }


#' @importFrom reshape2 acast
profile_likelihood_SCA3 <- function(Assessment, figure = TRUE, save_figure = TRUE,
                                    save_dir = getwd(), ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots)) stop("Sequence of R0 was not found. See help file.")
  if(!"transformed_h" %in% names(Assessment@obj$env$map) && !"h" %in% names(dots)) {
    stop("Sequence of h was not found. See help file.")
  }
  R0 <- dots$R0
  if(!"transformed_h" %in% names(Assessment@obj$env$map)) h <- dots$h else h <- Assessment@h

  profile.grid <- expand.grid(R0 = R0, h = h)
  nll <- rep(NA, nrow(profile.grid))
  params <- Assessment@info$params
  random <- Assessment@obj$env$random
  map <- Assessment@obj$env$map
  map$log_R0 <- map$transformed_h <- factor(NA)
  if(Assessment@info$data$SR_type == "BH") {
    transformed_h <- logit((profile.grid$h - 0.2)/0.8)
  } else {
    transformed_h <- log(profile.grid$h - 0.2)
  }
  for(i in 1:nrow(profile.grid)) {
    params$log_R0 <- log(profile.grid$R0[i] * Assessment@info$rescale)
    params$transformed_h <- transformed_h[i]
    obj <- MakeADFun(data = Assessment@info$data, parameters = params,
                     map = map, random = random, inner.control = Assessment@info$inner.control,
                     DLL = "MSEtool", silent = TRUE)
    opt <- optimize_TMB_model(obj, Assessment@info$control)[[1]]

    if(!is.character(opt)) nll[i] <- opt$objective
  }
  profile.grid$nll <- nll - Assessment@opt$objective
  if(figure) {
    if(length(h) > 1) {
      z.mat <- acast(profile.grid, h ~ R0, value.var = "nll")
      contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]),
              nlevels = 20)

      h.MLE <- Assessment@h
      R0.MLE <- Assessment@R0
      points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
      if(save_figure) {
        Model <- Assessment@Model
        prepare_to_save_figure()

        create_png(file.path(plot.dir, "profile_likelihood.png"))
        contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]),
                nlevels = 20)
        points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
        dev.off()
        profile.file.caption <- c("profile_likelihood.png",
                                  "Joint profile likelihood of h and R0. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
        html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)",
                    captions = matrix(profile.file.caption, nrow = 1),
                    name = Assessment@Data@Name, report_type = "Profile_Likelihood")
        browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
      }
    } else {
      plot(profile.grid$R0, nll, typ = 'o', pch = 16, xlab = expression(R[0]), ylab = "Change in negative log-likelihood")
      abline(v = Assessment@SD$value[names(Assessment@SD$value) == "R0"], lty = 2)

      if(save_figure) {
        Model <- Assessment@Model
        prepare_to_save_figure()

        create_png(file.path(plot.dir, "profile_likelihood.png"))
        plot(profile.grid$R0, nll, typ = 'o', pch = 16, xlab = expression(R[0]), ylab = "Change in negative log-likelihood")
        abline(v = Assessment@SD$value[names(Assessment@SD$value) == "R0"], lty = 2)
        dev.off()
        profile.file.caption <- c("profile_likelihood.png",
                                  "Profile likelihood of R0. Vertical, dashed line indicates maximum likelihood estimate.")
        html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)",
                    captions = matrix(profile.file.caption, nrow = 1),
                    name = Assessment@Data@Name, report_type = "Profile_Likelihood")
        browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
      }

    }

  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_SCA3 <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()
  data <- info$data
  n_y <- data$n_y

  Year <- c(info$Year, max(info$Year) + 1)
  C_hist <- data$C_hist
  I_hist <- data$I_hist
  CAA_hist <- data$CAA_hist
  CAA_n <- data$CAA_n
  est_rec_dev <- data$est_rec_dev
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, SSB, SSB_SSBMSY, SSB_SSB0, N, R, U, U_UMSY, log_rec_dev
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 9))
  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    n_y_ret <- n_y - i
    data$n_y <- n_y_ret
    data$C_hist <- C_hist[1:n_y_ret]
    data$I_hist <- I_hist[1:n_y_ret]
    data$CAA_hist <- CAA_hist[1:n_y_ret, ]
    data$CAA_n <- CAA_n[1:n_y_ret]
    data$est_rec_dev <- est_rec_dev[1:n_y_ret]
    params$log_rec_dev <- rep(0, n_y_ret)

    map <- obj$env$map
    if(any(names(map) == "log_rec_dev")) {
      new_map <- as.numeric(map$log_rec_dev) - i
      map$log_rec_dev <- factor(new_map[new_map > 0])
    }

    obj2 <- MakeADFun(data = data, parameters = params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- get_MSY(Arec = report$Arec, Brec = report$Brec, M = data$M, weight = data$weight, mat = data$mat,
                        vul = report$vul, SR = data$SR_type)
      report <- c(report, ref_pt)
      if(info$rescale != 1) {
        vars_div <- c("meanR", "B", "E", "CAApred", "CN", "N", "VB",
                      "R", "MSY", "VBMSY", "RMSY", "BMSY", "EMSY", "VB0", "R0",
                      "B0", "E0", "N0")
        vars_mult <- "Brec"
        var_trans <- c("R0", "q")
        fun_trans <- c("/", "*")
        fun_fixed <- c("log", NA)
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      SSB <- c(report$E, rep(NA, i))
      SSB_SSBMSY <- SSB/report$EMSY
      SSB_SSB0 <- SSB/report$E0
      R <- c(report$R, rep(NA, i))
      N <- c(rowSums(report$N), rep(NA, i))
      U <- c(report$U, rep(NA, i + 1))
      U_UMSY <- U/report$UMSY
      log_rec_dev <- c(report$log_rec_dev, rep(NA, i + 1))

      retro_ts[i+1, , ] <- cbind(Year, SSB, SSB_SSBMSY, SSB_SSB0, R, N, U, U_UMSY, log_rec_dev)
      retro_est[i+1, , ] <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev", ]

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }
  }
  if(figure) {
    fix_h <- "transformed_h" %in% names(obj$env$map)
    plot_retro_SCA3(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                    nyr_label = 0:nyr, color = rich.colors(nyr+1), fix_h, SR = data$SR_type)
  }
  # Need to write legend?
  legend <- NULL
  return(invisible(list(legend = legend, retro_ts = retro_ts, retro_est = retro_est)))
}


plot_retro_SCA3 <- function(retro_ts, retro_est, save_figure = FALSE,
                            save_dir = getwd(), nyr_label, color, fix_h, SR) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Spawning Stock Biomass", expression(SSB/SSB[MSY]), expression(SSB/SSB[0]), "Recruitment",
                "Population Abundance (N)", "Exploitation rate (U)",
                expression(U/U[MSY]), "Recruitment deviations")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "SCA3"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(abs(retro_ts[, , i+1]), na.rm = TRUE)
    if(i < n_tsplots) {
      ylim <- c(0, 1.1 * y.max)
    } else ylim <- c(-y.max, y.max)
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
                                               c("spawning stock biomass", "SSB/SSBMSY", "spawning depletion", "recruitment",
                                                 "abundance", "exploitation", "U/UMSY", "recruitment deviations",
                                                 "R0 estimate"), "."), stringsAsFactors = FALSE)
    if(!fix_h) {
      ret.file.caption <- rbind(ret.file.caption,
                                c(paste0("retrospective_", n_tsplots+2, ".png"), "Retrospective pattern in steepness estimate."))
    }

    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Statistical Catch-at-Age (SCA3)", captions = ret.file.caption,
                name = Assessment@Data@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}


