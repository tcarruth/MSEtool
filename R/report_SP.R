summary_SP <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  input_parameters <- data.frame()

  derived <- data.frame(Value = c(TMB_report$r, TMB_report$K, TMB_report$BMSY),
                        Description = c("Intrinsic rate of population increase", "Carrying capacity",
                                        "Biomass at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("r", "K", "BMSY")

  model_estimates <- summary(SD)
  model_estimates <- model_estimates[model_estimates[, 2] > 0, ]

  output <- list(model = "Delay Difference",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)

}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_SP <- function(Assessment, save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary_SP(Assessment)
    html_report(plot.dir, model = "Surplus Production",
                current_status = index.report$current_status,
                derived_quantities = index.report$derived_quantities,
                model_estimates = index.report$model_estimates,
                name = Data@Name, report_type = "Index")
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

  I_hist <- info$data$I_hist

  plot_timeseries(Year, I_hist, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(Year, I_hist, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_index.png", "Index time series."))
  }

  if(!is.na(Data@CV_Cat[1])) {
    plot_timeseries(Year, I_hist, obs_CV = Data@CV_Ind, label = "Index")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_index_with_CI.png"))
      plot_timeseries(Year, I_hist, obs_CV = Data@CV_Ind, label = "Index")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_index_with_CI.png", "Index time series with 95% confidence interval."))
    }
  }

  if(save_figure) {
    html_report(plot.dir, model = "Surplus Production",
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

  Uy <- names(U_UMSY)[length(U_UMSY)]
  plot_normalvar(U_UMSY[length(U_UMSY)], SE_U_UMSY_final, label = bquote(U[.(Uy)]/U[MSY]))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_U_UMSYestimate.png"))
    plot_normalvar(U_UMSY[length(U_UMSY)], SE_U_UMSY_final, label = bquote(widehat(U[.(Uy)]/U[MSY])))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_U_UMSYestimate.png",
                                   paste0("Estimate of U/UMSY in ", Uy, ", distribution based on
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

  plot_timeseries(Year, I_hist, Index, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index.png"))
    plot_timeseries(Year, I_hist, Index, label = "Index")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index.png", "Observed (black) and predicted (red) index."))
  }

  plot_residuals(Year, log(I_hist/Index), label = "log(Index) Residual")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index_residual.png"))
    plot_residuals(Year, log(I_hist/Index), label = "log(Index) Residual")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index_residual.png", "Index residuals in log-space."))
  }

  qqnorm(log(I_hist/Index), main = "")
  qqline(log(I_hist/Index))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_index_qqplot.png"))
    qqnorm(log(I_hist/Index), main = "")
    qqline(log(I_hist/Index))
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_index_qqplot.png", "QQ-plot of index residuals in log-space."))
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

  plot_yield_SP(TMB_report, UMSY, MSY, xaxis = "U")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_U.png"))
    plot_yield_SP(TMB_report, UMSY, MSY, xaxis = "U")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_yield_curve_U.png", "Yield plot relative to exploitation."))
  }

  plot_yield_SP(TMB_report, UMSY, MSY, xaxis = "Depletion")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_yield_curve_B_B0.png"))
    plot_yield_SP(TMB_report, UMSY, MSY, xaxis = "Depletion")
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
    html_report(plot.dir, model = "Surplus Production",
                captions = assess.file.caption, name = Data@Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }
  invisible()
}


#' @importFrom reshape2 acast
profile_likelihood_SP <- function(Assessment, figure = TRUE, save_figure = TRUE,
                                  save_dir = getwd(), ...) {
  dots <- list(...)
  if(!"UMSY" %in% names(dots)) stop("Sequence of UMSY was not found. See help file.")
  if(!"MSY" %in% names(dots)) stop("Sequence of MSY was not found. See help file.")
  UMSY <- dots$UMSY
  MSY <- dots$MSY

  profile.grid <- expand.grid(UMSY = UMSY, MSY = MSY)
  nll <- rep(NA, nrow(profile.grid))
  for(i in 1:nrow(profile.grid)) {
    logit_UMSY <- log(profile.grid[i, 1]/(1-profile.grid[i, 1]))
    log_MSY <- log(profile.grid[i, 2])
    nll[i] <- Assessment@obj$fn(x = c(logit_UMSY, log_MSY))
  }
  profile.grid$nll <- nll - min(nll, na.rm = TRUE)
  if(figure) {
    z.mat <- acast(profile.grid, UMSY ~ MSY, value.var = "nll")
    contour(x = UMSY, y = MSY, z = z.mat, xlab = expression(U[MSY]), ylab = "MSY",
            nlevels = 20)

    UMSY.MLE <- Assessment@UMSY
    MSY.MLE <- Assessment@MSY
    points(UMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      contour(x = UMSY, y = MSY, z = z.mat, xlab = expression(U[MSY]), ylab = "MSY",
              nlevels = 20)
      points(UMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
      dev.off()
      profile.file.caption <- c("profile_likelihood.png",
                                "Joint profile likelihood of UMSY and MSY. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      html_report(plot.dir, model = "Surplus Production",
                  captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Data@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile.grid)
}




#' @importFrom gplots rich.colors
retrospective_SP <- function(Assessment, nyr, figure = TRUE,
                             save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()

  data <- info$data
  ny <- data$ny

  Year <- info$Year
  Year <- c(Year, max(Year) + 1)
  C_hist <- data$C_hist
  I_hist <- data$I_hist
  params <- as.list(obj$env$last.par.best)
  params$log_Binit_frac <- info$params$log_Binit_frac
  params$log_n <- info$params$log_n
  map <- obj$env$map

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, U, relU, relB, dep
  retro_ts <- array(NA, dim = c(nyr + 1, ny + 1, 6))
  retro_est <- array(NA, dim = c(nyr + 1, dim(summary(SD))))

  for(i in 0:nyr) {
    ny_ret <- ny - i
    data$ny <- ny_ret
    data$C_hist <- C_hist[1:ny_ret]
    data$I_hist <- I_hist[1:ny_ret]

    obj <- MakeADFun(data = data, parameters = params, map = map,
                     DLL = "MSEtool", silent = TRUE)
    opt <- optimize_TMB_model(obj)
    SD <- get_sdreport(obj, opt)

    if(!is.character(opt) && opt$convergence == 0 && !is.character(SD)) {
      report <- obj$report(obj$env$last.par.best)
	  B <- c(report$B, rep(NA, i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/report$K

      U <- c(report$U, rep(NA, 1 + i))
      U_UMSY <- U/report$UMSY

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, U, U_UMSY)
      retro_est[i+1, , ] <- summary(SD)

    } else {
      warning(paste("Non-convergence when", i, "years of data were removed."))
    }

  }
  if(figure) {
    plot_retro_SP(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                  nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }
  # Need to write legend
  legend <- NULL
  return(invisible(list(legend = legend, retro_ts = retro_ts, retro_est = retro_est)))
}

plot_retro_SP <- function(retro_ts, retro_est, save_figure = FALSE,
                              save_dir = getwd(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]),
                "Exploitation rate (U)", expression(U/U[MSY]))
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "SP"
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
    if(i %in% c(2, 5)) abline(h = 1, lty = 2)

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
      if(i %in% c(2, 5)) abline(h = 1, lty = 2)
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
                                               c("biomass", "B/BMSY", "biomass depletion",
                                                 "exploitation", "U/UMSY", "UMSY estimate", "MSY estimate"), "."))

    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Surplus Production", captions = ret.file.caption,
                name = Assessment@Data@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}


plot_yield_SP <- function(report, umsy, msy, BKratio = seq(0, 1, 0.01),
                          xaxis = c("U", "Biomass", "Depletion")) {
  K <- report$K
  gamma.par <- report$gamma
  n <- report$n
  BMSY <- report$BMSY
  Yield <- gamma.par * msy * (BKratio - BKratio^n)
  Biomass <- BKratio * K
  u.vector <- Yield/Biomass

  if(xaxis == "U") {
    plot(u.vector, Yield, typ = 'l', xlab = "Exploitation rate (U)", ylab = "Equilibrium yield")
    segments(x0 = umsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass, Yield, typ = 'l', xlab = "Biomass", ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(BKratio, Yield, typ = 'l',
         xlab = expression(B/B[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/K, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/K, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible()
}
