
#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_report_SP <- function(Assessment, figure = TRUE, save_figure = FALSE,
                                   save_dir = getwd()) {
  assign_Assessment_var()

  Year <- info$Year
  n_y <- info$data$n_y
  Catch <- info$data$Catch
  Index <- info$data$Index
  Index[Index <= 0] <- NA

  logit.umsy <- as.numeric(obj$env$last.par.best[1])
  logit.umsy.sd <- sqrt(diag(SD$cov.fixed)[1])
  umsy <- 1/(1 + exp(-logit.umsy))
  log.msy <- as.numeric(obj$env$last.par.best[2])
  log.msy.sd <- sqrt(diag(SD$cov.fixed)[2])
  msy <- exp(log.msy)

  input_data <- data.frame(Year = Year, Catch = Catch, Index = Index)
  input_parameters <- data.frame()
  Ipred <- data.frame(Year = Year, Value = report$Ipred,
                      Description = "Predicted Index")
  B <- data.frame(Year = c(Year, max(Year) + 1), Value = report$Biomass,
                  Description = "Biomass")
  U <- data.frame(Year = Year, Value = report$U, Description = "Exploitation rate (U)")
  relU <- data.frame(Year = Year, Value = report$relU,
                     Description = "U/UMSY")
  relB <- data.frame(Year = c(Year, max(Year) + 1), Value = report$relB,
                     Description = "B/BMSY")
  dep <- data.frame(Year = c(Year, max(Year) + 1), Value = report$Biomass/report$K,
                    Description = "B/B0")
  model_output <- rbind(B, U, relU, relB, dep)

  derived <- data.frame(Value = c(report$r, report$K, report$BMSY),
                        Description = c("Intrinsic rate of population increase", "Carrying capacity",
                                        "Biomass at MSY"))
  rownames(derived) <- c("r", "K", "BMSY")

  output <- list(model = "Surplus Production",
                 input_data = input_data, input_parameters = input_parameters,
                 model_output = model_output, derived_quantities = derived,
                 model_estimates = summary(SD))

  if(figure) {
    if(save_figure) {
      prepare_to_save_figure()
      html_report(plot.dir, model = "Surplus Production",
                  input_parameters = NULL,
                  model_estimates = output$model_estimates,
                  derived_quantities = output$derived_quantities,
                  report_type = "Index")
    }

    plot_timeseries(Year, Catch, label = "Catch")
	  if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_1a_catch.png"))
      plot_timeseries(Year, Catch, label = "Catch")
      dev.off()
      data.file.caption <- c("data_1a_catch.png", "Catch time series")
    }

    if(!is.na(Data@CV_Cat)) {
      plot_timeseries(Year, Catch, obs_CV = Data@CV_Cat, label = "Catch")
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "data_1b_catch_with_CV.png"))
        plot_timeseries(Year, Catch, obs_CV = Data@CV_Cat, label = "Catch")
        dev.off()
        data.file.caption <- rbind(data.file.caption,
                                   c("data_1b_catch_with_CV.png", "Catch time series with 95% confidence interval."))
      }
    }

    plot_timeseries(Year, Index, label = "Index")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "data_2a_index.png"))
      plot_timeseries(Year, Index, label = "Index")
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c("data_2a_index.png", "Index time series."))

    }

    if(!is.na(Data@CV_Cat)) {
      plot_timeseries(Year, Index, obs_CV = Data@CV_Ind, label = "Index")
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "data_2b_index_with_CV.png"))
        plot_timeseries(Year, Index, obs_CV = Data@CV_Ind, label = "Index")
        dev.off()
        data.file.caption <- rbind(data.file.caption,
                                   c("data_2b_index_with_CV.png", "Index time series with 95% confidence interval."))
      }
    }

    if(save_figure) {
      html_report(plot.dir, model = "Surplus Production",
                  captions = data.file.caption, report_type = "Data")
    }

    plot_betavar(logit.umsy, logit.umsy.sd, logit = TRUE, label = expression(hat(U)[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_1a_UMSYestimate.png"))
      plot_betavar(logit.umsy, logit.umsy.sd, logit = TRUE, label = expression(hat(U)[MSY]))
      dev.off()
      assess.file.caption <- c("assessment_1a_UMSYestimate.png", "Estimate of UMSY, distribution based on normal approximation of estimated covariance matrix.")
    }

    plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_1b_MSYestimate.png"))
      plot_lognormalvar(log.msy, log.msy.sd, logtransform = TRUE, label = expression(widehat(MSY)))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_1b_MSYestimate.png", "Estimate of MSY, distribution based on normal approximation of estimated covariance matrix."))
    }

    if(length(Rec@TAC) > 1) {
      plot_TAC(Rec@TAC)
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "assessment_2_TAC.png"))
        plot_TAC(Rec@TAC)
        dev.off()
        assess.file.caption <- rbind(assess.file.caption,
                                     c("assessment_2_TAC.png", "TAC recommendation, distribution based on resampling the covariance matrix."))
      }
    }

    plot_timeseries(Year, Index, report$Ipred, label = "Index")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_4a_index.png"))
      plot_timeseries(Year, Index, report$Ipred, label = "Index")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_4a_index.png", "Observed (black) and predicted (red) index."))

    }

    plot_residuals(Year, log(Index/report$Ipred), label = "log(Index) Residual")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_4b_index_residual.png"))
      plot_residuals(Year, log(Index/report$Ipred), label = "log(Index) Residual")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_4b_index_residual.png", "Index residuals in log-space."))
    }

    qqnorm(log(Index/report$Ipred), main = "")
    qqline(log(Index/report$Ipred))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_4c_index_qqplot.png"))
      qqnorm(log(Index/report$Ipred), main = "Q-Q Plot of log(Index) Residual")
      qqline(log(Index/report$Ipred))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_4c_index_qqplot.png", "QQ-plot of index residuals in log-space."))

    }

    plot_timeseries(B$Year, B$Value, label = "Biomass")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_6a_biomass.png"))
      plot_timeseries(B$Year, B$Value, label = "Biomass")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_6a_biomass.png", "Time series of biomass."))
    }

    plot_timeseries(relB$Year, relB$Value, label = expression(B/B[MSY]))
    abline(h = 1, lty = 2)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_6b_B_BMSY.png"))
      plot_timeseries(relB$Year, relB$Value, label = expression(B/B[MSY]))
      abline(h = 1, lty = 2)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_6b_B_BMSY.png", "Time series of B/BMSY."))
    }

    plot_timeseries(dep$Year, dep$Value, label = expression(B/B[0]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_6c_B_B0.png"))
      plot_timeseries(dep$Year, dep$Value, label = expression(B/B[0]))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_6c_B_B0.png", "Time series of biomass depletion."))
    }

    plot_timeseries(Year, U$Value, label = "Exploitation rate (U)")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_9a_exploitation.png"))
      plot_timeseries(Year, U$Value, label = "Exploitation rate (U)")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_9a_exploitation.png", "Time series of exploitation rate."))
    }

    plot_timeseries(Year, relU$Value, label = expression(U/U[MSY]))
    abline(h = 1, lty = 2)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_9b_U_UMSY.png"))
      plot_timeseries(Year, relU$Value, label = expression(U/U[MSY]))
      abline(h = 1, lty = 2)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_9b_U_UMSY.png", "Time series of U/UMSY."))
    }

    plot_Kobe(relB$Value, relU$Value)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_10_Kobe.png"))
      plot_Kobe(relB$Value, relU$Value)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_10_Kobe.png", "Kobe plot trajectory of stock."))
    }

    plot_yield_SP(report, umsy, msy, xaxis = "U")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_11a_yield_curve.png"))
      plot_yield_SP(report, umsy, msy, xaxis = "U")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_11a_yield_curve.png", "Yield plot relative to exploitation."))
    }

    plot_yield_SP(report, umsy, msy, xaxis = "Depletion")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_11b_yield_curve.png"))
      plot_yield_SP(report, umsy, msy, xaxis = "Depletion")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_11b_yield_curve.png", "Yield plot relative to depletion."))
    }

    plot_surplus_production(report$Biomass, report$K, Catch)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_12_surplus_production.png"))
      plot_surplus_production(report$Biomass, report$K, Catch)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_12_surplus_production.png", "Surplus production relative to depletion."))
    }

	  if(save_figure) {
	    html_report(plot.dir, model = "Surplus Production",
	                captions = assess.file.caption, report_type = "Assessment")
	    browseURL(file.path(plot.dir, "Assessment.html"))
    }
  }
  return(output)
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
  MLE <- as.numeric(Assessment@obj$env$last.par.best) # Max. likelihood est.
  UMSY.MLE <- 1/(1 + exp(-MLE[1]))
  MSY.MLE <- exp(MLE[2])

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
    points(UMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)

    if(save_figure) {
      MP <- Assessment@MP
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
                  report_type = "Profile_likelihood")
      browseURL(file.path(plot.dir, "Profile_likelihood.html"))
    }
  }
  return(profile.grid)
}




#' @importFrom gplots rich.colors
retrospective_SP <- function(Assessment, nyr, figure = TRUE,
                             save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_var()

  data <- info$data
  n_y <- data$n_y

  Year <- info$Year
  Year <- c(Year, max(Year) + 1)
  Catch <- data$Catch
  Index <- data$Index
  params <- as.list(obj$env$last.par.best)
  params$log_B1frac <- info$params$log_B1frac
  params$log_n <- info$params$log_n
  map <- obj$env$map

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, U, relU, relB, dep
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 6))
  class(SD) <- "sdreport"
  retro_est <- array(NA, dim = c(nyr+1, dim(summary(SD))))

  for(i in 0:nyr) {
    n_y_ret <- n_y - i
    Catch <- Catch[1:n_y_ret]
    Index <- Index[1:n_y_ret]
    data$n_y <- n_y_ret
    data$Catch <- Catch
    data$Index <- Index

    obj <- MakeADFun(data = data, parameters = params, map = map,
                     DLL = "MSEtool", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)

    if(opt$convergence == 0) {
      B <- c(obj$report()$Biomass, rep(NA, i))
      relB <- c(obj$report()$relB, rep(NA, i))
      dep <- B/obj$report()$K

      U <- c(obj$report()$U, rep(NA, 1 + i))
      relU <- c(obj$report()$relU, rep(NA, 1 + i))

      SD <- sdreport(obj)

      retro_ts[i+1, , ] <- cbind(Year, B, relB, dep, U, relU)
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
  return(list(legend = legend, retro_ts = retro_ts, retro_est = retro_est))
}

plot_retro_SP <- function(retro_ts, retro_est, save_figure = FALSE,
                              save_dir = getwd(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]),
                "Exploitation rate (U)", expression(U/U[MSY]))
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    MP <- "SP"
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
    html_report(plot.dir, model = "Surplus Production",
                captions = ret.file.caption, report_type = "Retrospective")
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
