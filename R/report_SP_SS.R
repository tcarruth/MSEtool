summary_SP_SS <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  Value <- numeric(0)
  Description <- character(0)
  rownam <- character(0)
  if("log_dep" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$dep)
    Description <- c(Description, "Initial depletion")
    rownam <- c(rownam, "dep")
  }
  if("log_n" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$n)
    Description <- c(Description, "Production exponent")
    rownam <- c(rownam, "n")
  }
  if("log_sigma" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$sigma)
    Description <- c(Description, "Index SD (log-space)")
    rownam <- c(rownam, "sigma")
  }
  if("log_tau" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$tau)
    Description <- c(Description, "Biomass deviation SD (log-space)")
    rownam <- c(rownam, "tau")
  }
  if(length(Value) == 0) input_parameters <- data.frame() else {
    input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
    rownames(input_parameters) <- rownam
  }

  derived <- data.frame(Value = c(TMB_report$r, TMB_report$K, TMB_report$BMSY),
                        Description = c("Intrinsic rate of population increase", "Carrying capacity",
                                        "Biomass at MSY"), stringsAsFactors = FALSE)
  rownames(derived) <- c("r", "K", "BMSY")

  if(!is.character(SD)) {
    if(is.null(obj$env$random)) {
      model_estimates <- summary(SD)[rownames(summary(SD)) != "log_B_dev", ]
      dev_estimates <- summary(SD)[rownames(summary(SD)) == "log_B_dev", ]
    } else {
      model_estimates <- rbind(summary(SD, "fixed"), summary(SD, "report"))
      dev_estimates <- summary(SD, "random")
    }

    rownames(dev_estimates) <- paste0(rownames(dev_estimates), "_", names(Dev)[Dev != 0])
    model_estimates <- model_estimates[is.na(model_estimates[, 2]) || model_estimates[, 2] > 0, ]
    model_estimates <- rbind(model_estimates, dev_estimates)
  } else {
    model_estimates <- SD
  }

  output <- list(model = "Surplus Production (State-Space)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_SP_SS <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary_SP_SS(Assessment)
    html_report(plot.dir, model = "Surplus Production (State-Space)",
                current_status = index.report$current_status,
                derived_quantities = index.report$derived_quantities,
                model_estimates = index.report$model_estimates,
                name = Name, report_type = "Index")
  }

  Year <- info$Year

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
    html_report(plot.dir, model = "Surplus Production (State-Space)",
                captions = data.file.caption, name = Name, report_type = "Data")
  }

  if(conv) {
    Fmsy.ind <- names(SD$par.fixed) == "log_FMSY"
    log.Fmsy <- SD$par.fixed[Fmsy.ind]
    log.Fmsy.sd <- sqrt(diag(SD$cov.fixed)[Fmsy.ind])

    plot_lognormalvar(log.Fmsy, log.Fmsy.sd, logtransform = TRUE, label = expression(hat(F)[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_FMSYestimate.png"))
      plot_lognormalvar(log.Fmsy, log.Fmsy.sd, logtransform = TRUE, label = expression(hat(F)[MSY]))
      dev.off()
      assess.file.caption <- c("assessment_FMSYestimate.png", "Estimate of FMSY, distribution based on normal approximation of estimated covariance matrix.")
    }

    msy.ind <- names(SD$par.fixed) == "log_MSY"
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

    By <- names(B_B0)[length(B_B0)]
    plot_normalvar(B_B0[length(B_B0)], SE_B_B0_final, label = bquote(B[.(By)]/B[0]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_B_B0estimate.png"))
      plot_normalvar(B_B0[length(B_B0)], SE_B_B0_final, label = bquote(widehat(B[.(By)]/B[0])))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_B_B0estimate.png",
                                     paste0("Estimate of B/B0 in ", By, ", distribution based on
                                          normal approximation of estimated covariance matrix.")))
    }
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
                                 c("assessment_catch.png", "Observed (black) and predicted (red) catch. In this model, observed catch should match predicted catch."))
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

  plot_residuals(as.numeric(names(Dev)), Dev, label = Dev_type)
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_Biomass_devs.png"))
    plot_residuals(as.numeric(names(Dev)), Dev, label = Dev_type)
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_Biomass_devs.png", "Time series of biomass deviations."))
  }

  if(conv) {
    plot_residuals(as.numeric(names(Dev)), Dev, SE_Dev, label = Dev_type)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_Biomass_devs_with_CI.png"))
      plot_residuals(as.numeric(names(Dev)), Dev, SE_Dev, label = Dev_type)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_Biomass_devs_with_CI.png", "Time series of biomass deviations with 95% confidence intervals."))
    }
  }


  plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing mortality F")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_F.png"))
    plot_timeseries(as.numeric(names(FMort)), FMort, label = "Fishing mortality F")
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

  plot_yield_SP(TMB_report, FMSY, MSY, xaxis = "F")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_F.png"))
    plot_yield_SP(TMB_report, FMSY, MSY, xaxis = "F")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_F.png", "Yield plot relative to exploitation."))
  }

  plot_yield_SP(TMB_report, FMSY, MSY, xaxis = "Depletion")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_B_B0.png"))
    plot_yield_SP(TMB_report, FMSY, MSY, xaxis = "Depletion")
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
    html_report(plot.dir, model = "Surplus Production (State-Space)",
                captions = assess.file.caption, name = Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }
  return(invisible())
}



#' @importFrom reshape2 acast
profile_likelihood_SP_SS <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  dots <- list(...)
  if(!"FMSY" %in% names(dots) && !"MSY" %in% names(dots)) stop("Sequence of neither FMSY nor MSY was found. See help file.")
  if(!is.null(dots$FMSY)) FMSY <- dots$FMSY else {
    FMSY <- Assessment@FMSY
    profile_par <- "MSY"
  }
  if(!is.null(dots$MSY)) MSY <- dots$MSY else {
    MSY <- Assessment@MSY
    profile_par <- "FMSY"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(FMSY = FMSY, MSY = MSY)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$log_FMSY <- log(profile_grid[i, 1])
    params$log_MSY <- log(profile_grid[i, 2] * Assessment@info$rescale)

    if(joint_profile) {
	  map$log_MSY <- map$log_FMSY <- factor(NA)
	} else { 
	  if(profile_par == "MSY") map$log_MSY <- factor(NA) else map$log_FMSY <- factor(NA)
    }
	obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, inner.control = Assessment@info$inner.control,
                      random = Assessment@obj$env$random, DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]

    if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(figure) {
    FMSY.MLE <- Assessment@FMSY
    MSY.MLE <- Assessment@MSY
    if(joint_profile) {
      z.mat <- acast(profile_grid, list("FMSY", "MSY"), value.var = "nll")
      contour(x = FMSY, y = MSY, z = z.mat, xlab = expression(F[MSY]), ylab = "MSY", nlevels = 20)
      points(FMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
    } else {
      if(profile_par == "FMSY") xlab <- expression(F[MSY]) else xlab <- "MSY"
      plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
    }

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      if(joint_profile) {
        z.mat <- acast(profile_grid, list("FMSY", "MSY"), value.var = "nll")
        contour(x = FMSY, y = MSY, z = z.mat, xlab = expression(F[MSY]), ylab = "MSY", nlevels = 20)
        points(FMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
		msg <- "Joint profile likelihood of FMSY and MSY. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate."
      } else {
        if(profile_par == "FMSY") xlab <- expression(F[MSY]) else xlab <- "MSY"
        plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)		
		msg <- paste0("Profile likelihood of ", profile_par, ". Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      }
      dev.off()

      profile.file.caption <- c("profile_likelihood.png", msg)

      html_report(plot.dir, model = "Surplus Production (State-Space)", captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile_grid)
}




#' @importFrom gplots rich.colors
retrospective_SP_SS <- function(Assessment, nyr, figure = TRUE,
                                save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  data <- info$data
  ny <- data$ny

  Year <- info$Year
  Year <- c(Year, max(Year) + 1)
  C_hist <- data$C_hist
  I_hist <- data$I_hist
  est_B_dev <- data$est_B_dev
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, F, relF, relB, log_B_dev
  retro_ts <- array(NA, dim = c(nyr+1, ny + 1, 7))
  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_B_dev", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    ny_ret <- ny - i
    data$ny <- ny_ret
    data$Catch <- C_hist[1:ny_ret]
    data$Index <- I_hist[1:ny_ret]
    data$est_B_dev <- est_B_dev[1:ny_ret]
    params$log_B_dev <- rep(0, ny_ret)

    map <- obj$env$map
    if("log_B_dev" %in% names(map)) {
      new_map <- as.numeric(map$log_B_dev) - i
      map$log_B_dev <- factor(new_map[new_map > 0])
    }

    obj2 <- MakeADFun(data = data, parameters = params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      if(rescale != 1) {
        vars_div <- c("B", "BMSY", "K", "MSY", "Cpred", "SP")
        vars_mult <- NULL
        var_trans <- c("MSY", "K", "q")
        fun_trans <- c("/", "/", "*")
        fun_fixed <- c("log", NA, NA)
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }
      B <- c(report$B, rep(NA, i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/report$K

      FMort <- c(report$F, rep(NA, 1 + i))
      F_FMSY <- FMort/report$FMSY
      log_B_dev <- c(report$log_B_dev, rep(NA, 1 + i))

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, FMort, F_FMSY, log_B_dev)
      retro_est[i+1, , ] <- summary(SD)[rownames(summary(SD)) != "log_B_dev", ]

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }

  }

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1], retro_est[, 3:4, 1],
                                 ts_lab = c("Biomass", "B_BMSY", "B_B0", "F", "F_FMSY", "Biomass deviations"),
                                 est_lab = c("FMSY estimate", "MSY estimate"))

  if(figure) {
    plot_retro_SP_SS(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir, nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }

  return(Mohn_rho)
}

plot_retro_SP_SS <- function(retro_ts, retro_est, save_figure = FALSE, save_dir = tempdir(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]),
                "Fishing mortality (F)", expression(F/F[MSY]), "Biomass deviations")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "SP_SS"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(retro_ts[, , i+1], na.rm = TRUE)
    if(i < n_tsplots) {
      ylim <- c(0, 1.1 * y.max)
    } else ylim <- c(-y.max, y.max)
    plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i], ylim = ylim, col = color[1])
    for(j in 2:length(nyr_label)) lines(Year, retro_ts[j, , i+1], col = color[j])
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
    abline(h = 0, col = 'grey')
    if(i %in% c(2, 5)) abline(h = 1, lty = 2)

    if(save_figure) {
      create_png(filename = file.path(plot.dir, paste0("retrospective_", i, ".png")))
      plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i], ylim = ylim, col = color[1])
      for(j in 2:length(nyr_label)) lines(Year, retro_ts[j, , i+1], col = color[j])
      legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
      abline(h = 0, col = 'grey')
      if(i %in% c(2, 5)) abline(h = 1, lty = 2)
      dev.off()
    }
  }

  plot_lognormalvar(retro_est[, 1, 1], retro_est[, 1, 2], logtransform = TRUE, label = expression(hat(F)[MSY]), color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 1, ".png")))
    plot_lognormalvar(retro_est[, 1, 1], retro_est[, 1, 2], logtransform = TRUE, label = expression(hat(F)[MSY]), color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
    dev.off()
  }

  plot_lognormalvar(retro_est[, 2, 1], retro_est[, 2, 2], logtransform = TRUE, label = expression(widehat(MSY)), color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 2, ".png")))
    plot_lognormalvar(retro_est[, 2, 1], retro_est[, 2, 2], logtransform = TRUE, label = expression(widehat(MSY)), color = color)
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
    dev.off()
  }

  if(save_figure) {
    ret.file.caption <- data.frame(x1 = paste0("retrospective_", c(1:(n_tsplots+2)), ".png"),
                                   x2 = paste0("Retrospective pattern in ",
                                               c("biomass", "B/BMSY", "biomass depletion", "F", "F/FMSY", "Biomass deviations",
                                                 "FMSY estimate", "MSY estimate"), "."))
    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Surplus Production (State-Space)", captions = ret.file.caption, name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}



