summary_spict <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  input_parameters <- data.frame()

  derived <- data.frame(Value = c(SD$value["r"], TMB_report$Bmsy, TMB_report$Fmsy, exp(SD$value["logalpha"]), exp(SD$value["logbeta"])),
                        Description = c("Intrinsic rate of population increase", "Fishing mortality at MSY", "Biomass at MSY",
                                        "Ratio of index and biomass standard deviations", "Ratio of catch and F standard deviations"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("r", "K", "BMSY", "alpha", "beta")

  model_estimates <- summary(SD, "fixed")
  SD_exp <- model_estimates
  SD_exp[, 1] <- exp(SD_exp[, 1])
  SD_exp[, 2] <- SD_exp[, 1] * SD_exp[, 2]
  rownames(SD_exp) <- c("MSY", "K", "q", "n", "sd_b", "sd_f", "sd_i", "sd_c")
  model_estimates <- rbind(model_estimates, SD_exp)

  output <- list(model = "spict", current_status = current_status,
                 input_parameters = input_parameters, derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_spict <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary_spict(Assessment)
    html_report(plot.dir, model = "SPiCT",
                current_status = index.report$current_status,
                derived_quantities = index.report$derived_quantities,
                model_estimates = index.report$model_estimates,
                name = Name, report_type = "Index")
  }

  Year <- as.numeric(names(Index))

  plot_timeseries(Year, Obs_Catch, label = "Catch")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_catch.png"))
    plot_timeseries(Year, Obs_Catch, label = "Catch")
    dev.off()
    data.file.caption <- c("data_catch.png", "Catch time series")
  }

  plot_timeseries(Year, Obs_Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(Year, Obs_Index, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_index.png", "Index time series."))
  }

  if(save_figure) {
    html_report(plot.dir, model = "SPiCT",
                captions = data.file.caption, name = Name, report_type = "Data")
  }

  if(conv) {

    plot_normalvar(FMSY, SE_FMSY, label = expression(hat(F)[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_FMSYestimate.png"))
      plot_normalvar(FMSY, SE_FMSY, label = expression(hat(F)[MSY]))
      dev.off()
      assess.file.caption <- c("assessment_FMSYestimate.png", "Estimate of FMSY, distribution based on normal approximation of estimated covariance matrix.")
    }

    msy.ind <- names(SD$par.fixed) == "logm"
    log.msy <- SD$par.fixed[msy.ind]
    log.msy.sd <- sqrt(diag(SD$cov.fixed)[msy.ind])

    plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_MSYestimate.png"))
      plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_MSYestimate.png", "Estimate of MSY, distribution based on normal approximation of estimated covariance matrix."))
    }

    Fy <- names(F_FMSY)[length(F_FMSY)]
    plot_normalvar(F_FMSY[length(F_FMSY)], SE_F_FMSY_final, label = bquote(F[.(Fy)]/F[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_F_FMSYestimate.png"))
      plot_normalvar(F_FMSY[length(F_FMSY)], SE_F_FMSY_final, label = bquote(F[.(Fy)]/F[MSY]))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_F_FMSYestimate.png",
                                     paste0("Estimate of F/FMSY in ", Fy, ", distribution based on
                                            normal approximation of estimated covariance matrix.")))
    }

    By <- names(B_BMSY)[length(B_BMSY)]
    plot_normalvar(B_BMSY[length(B_BMSY)], SE_B_BMSY_final, label = bquote(B[.(By)]/B[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_B_BMSYestimate.png"))
      plot_normalvar(B_BMSY[length(B_BMSY)], SE_B_BMSY_final, label = bquote(widehat(B[.(By)]/B[MSY])))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_B_BMSYestimate.png",
                                     paste0("Estimate of B/BMSY in ", By, ", distribution based on
                                            normal approximation of estimated covariance matrix.")))
    }

    #By <- names(B_B0)[length(B_B0)]
    #plot_normalvar(B_B0[length(B_B0)], SE_B_B0_final, label = bquote(B[.(By)]/B[0]))
    #if(save_figure) {
    #  create_png(filename = file.path(plot.dir, "assessment_B_B0estimate.png"))
    #  plot_normalvar(B_B0[length(B_B0)], SE_B_B0_final, label = bquote(widehat(B[.(By)]/B[0])))
    #  dev.off()
    #  assess.file.caption <- rbind(assess.file.caption,
    #                               c("assessment_B_B0estimate.png",
    #                                 paste0("Estimate of B/B0 in ", By, ", distribution based on
    #                                        normal approximation of estimated covariance matrix.")))
    #}
  }

  plot_timeseries(Year, Obs_Index, Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index.png"))
    plot_timeseries(Year, Obs_Index, Index, label = "Index")
    dev.off()
    if(conv) assess.file.caption <- rbind(assess.file.caption,
                                          c("assessment_index.png", "Observed (black) and predicted (red) index."))
    else assess.file.caption <- c("assessment_index.png", "Observed (black) and predicted (red) index.")
  }

  plot_residuals(Year, log(Obs_Index/Index), label = "log(Index) Residual")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index_residual.png"))
    plot_residuals(Year, log(Obs_Index/Index), label = "log(Index) Residual")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index_residual.png", "Index residuals in log-space. One-step-ahead residuals available in spict package."))
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
                                 c("assessment_catch.png", "Observed (black) and predicted (red) catch."))
  }

  plot_residuals(Year, log(Obs_Catch/Catch), label = "log(Catch) Residual")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch_residual.png"))
    plot_residuals(Year, log(Obs_Catch/Catch), label = "log(Catch) Residual")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch_residual.png", "Catch residuals in log-space. One-step-ahead residuals available in spict package."))
  }

  qqnorm(log(Obs_Catch/Catch), main = "")
  qqline(log(Obs_Catch/Catch))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_catch_qqplot.png"))
    qqnorm(log(Obs_Catch/Catch), main = "")
    qqline(log(Obs_Catch/Catch))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_catch_qqplot.png", "QQ-plot of catch residuals in log-space."))
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

  plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing Mortality F")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_F.png"))
    plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing Mortality F")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_F.png", "Time series of exploitation rate."))
  }

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

  pcurve_report <- list(K = B0, BMSY = BMSY, n = exp(SD$par.fixed["logn"]))
  plot_yield_SP(pcurve_report, FMSY, MSY, xaxis = "F")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_F.png"))
    plot_yield_SP(pcurve_report, FMSY, MSY, xaxis = "F")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_F.png", "Yield plot relative to exploitation."))
  }

  plot_yield_SP(pcurve_report, FMSY, MSY, xaxis = "Depletion")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_B_B0.png"))
    plot_yield_SP(pcurve_report, FMSY, MSY, xaxis = "Depletion")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_B_B0.png", "Yield plot relative to depletion."))
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
    html_report(plot.dir, model = "Surplus Production",
                captions = assess.file.caption, name = Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }
  invisible()
}



profile_likelihood_spict <- function(Assessment, figure = TRUE, save_figure = FALSE, save_dir = tempdir(), ...) {
  stop("Profiling currently not supported for spict in MSEtool.", call. = FALSE)
}



retrospective_spict <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = tempdir()) {
  stop("Retrospective analysis is currently not supported for spict in MSEtool, but is available in the spict package.", call. = FALSE)
}


