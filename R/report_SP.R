summary_SP <- function(Assessment) {
  assign_Assessment_slots()

  current_status <- data.frame(Value = c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  if(length(obj$env$map) == 0) input_parameters <- data.frame()
  else {
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
    input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
    rownames(input_parameters) <- rownam
  }

  derived <- data.frame(Value = c(TMB_report$r, TMB_report$K, TMB_report$BMSY),
                        Description = c("Intrinsic rate of population increase", "Carrying capacity",
                                        "Biomass at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("r", "K", "BMSY")

  if(!is.character(SD)) {
    model_estimates <- summary(SD)
    model_estimates <- model_estimates[is.na(model_estimates[, 2]) || model_estimates[, 2] > 0, ]
  } else {
    model_estimates <- SD
  }

  output <- list(model = "Surplus Production", current_status = current_status,
                 input_parameters = input_parameters, derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)

}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_SP <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary_SP(Assessment)
    html_report(plot.dir, model = "Surplus Production",
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
    html_report(plot.dir, model = "Surplus Production",
                captions = data.file.caption, name = Name, report_type = "Data")
  }

  if(conv) {
    umsy.ind <- names(SD$par.fixed) == "logit_UMSY"
    logit.umsy <- SD$par.fixed[umsy.ind]
    logit.umsy.sd <- sqrt(diag(SD$cov.fixed)[umsy.ind])

    plot_betavar(logit.umsy, logit.umsy.sd, is_logit = TRUE, label = expression(hat(U)[MSY]))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_UMSYestimate.png"))
      plot_betavar(logit.umsy, logit.umsy.sd, is_logit = TRUE, label = expression(hat(U)[MSY]))
      dev.off()
      assess.file.caption <- c("assessment_UMSYestimate.png", "Estimate of UMSY, distribution based on normal approximation of estimated covariance matrix.")
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


#' @importFrom reshape2 acast
profile_likelihood_SP <- function(Assessment, figure = TRUE, save_figure = FALSE,
                                  save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"UMSY" %in% names(dots)) stop("Sequence of UMSY was not found. See help file.")
  if(!"MSY" %in% names(dots)) stop("Sequence of MSY was not found. See help file.")
  UMSY <- dots$UMSY
  MSY <- dots$MSY

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile.grid <- expand.grid(UMSY = UMSY, MSY = MSY)
  nll <- rep(NA, nrow(profile.grid))
  for(i in 1:nrow(profile.grid)) {
    params$logit_UMSY <- logit(profile.grid[i, 1])
    params$log_MSY <- log(profile.grid[i, 2] * Assessment@info$rescale)
    if(length(Assessment@obj$par) == 2) {
      nll[i] <- Assessment@obj$fn(x = c(params$logit_UMSY, params$log_MSY))
    } else { # More than 2 parameters
      map$logit_UMSY <- map$log_MSY <- factor(NA)
      obj2 <- MakeADFun(data = Assessment@info$data, parameters = Assessment@info$params,
                        map = map, DLL = "MSEtool", silent = TRUE)

      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll[i] <- opt2$objective
    }
  }
  profile.grid$nll <- nll - Assessment@opt$objective
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
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile.grid)
}




#' @importFrom gplots rich.colors
retrospective_SP <- function(Assessment, nyr, figure = TRUE,
                             save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  data <- info$data
  ny <- data$ny

  Year <- info$Year
  Year <- c(Year, max(Year) + 1)
  C_hist <- data$C_hist
  I_hist <- data$I_hist
  params <- info$params
  map <- obj$env$map

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, U, relU, relB, dep
  retro_ts <- array(NA, dim = c(nyr + 1, ny + 1, 6))
  retro_est <- array(NA, dim = c(nyr + 1, dim(summary(SD))))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    ny_ret <- ny - i
    data$ny <- ny_ret
    data$C_hist <- C_hist[1:ny_ret]
    data$I_hist <- I_hist[1:ny_ret]

    obj2 <- MakeADFun(data = data, parameters = params, map = map,
                      DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      if(rescale != 1) {
        vars_div <- c("B", "BMSY", "SP", "K", "MSY")
        vars_mult <- NULL
        var_trans <- c("MSY", "K", "q")
        fun_trans <- c("/", "/", "*")
        fun_fixed <- c("log", NA, NA)
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }
	    B <- c(report$B, rep(NA, i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/report$K

      U <- c(report$U, rep(NA, 1 + i))
      U_UMSY <- U/report$UMSY

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, U, U_UMSY)
      retro_est[i+1, , ] <- summary(SD)

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }

  }

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1], retro_est[, 3:4, 1],
                                 ts_lab = c("Biomass", "B_BMSY", "B_B0", "U", "U_UMSY"),
                                 est_lab = c("UMSY estimate", "MSY estimate"))

  if(figure) {
    plot_retro_SP(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                  nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }

  return(Mohn_rho)
}

plot_retro_SP <- function(retro_ts, retro_est, save_figure = FALSE,
                          save_dir = tempdir(), nyr_label, color) {
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

  plot_betavar(retro_est[, 1, 1], retro_est[, 1, 2], is_logit = TRUE,
               label = expression(hat(U)[MSY]), color = color)
  legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n",
         title = "Years removed:")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("retrospective_", n_tsplots + 1, ".png")))
    plot_betavar(retro_est[, 1, 1], retro_est[, 1, 2], is_logit = TRUE,
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
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}


plot_yield_SP <- function(report, umsy, msy, BKratio = seq(0, 1, 0.01),
                          xaxis = c("U", "Biomass", "Depletion"), relative_yaxis = FALSE) {
  K <- report$K
  n <- report$n
  BMSY <- report$BMSY

  if(n == 1) {
    Yield <- ifelse(BKratio == 0, 0, -exp(1) * msy * BKratio * log(BKratio))
  } else {
    gamma.par <- n^(n/(n-1))/n-1
    Yield <- gamma.par * msy * (BKratio - BKratio^n)
  }

  Biomass <- BKratio * K
  u.vector <- Yield/Biomass

  if(relative_yaxis) {
    Yield <- Yield/max(Yield)
    ylab <- "Relative Equilibrium Yield"
  } else ylab <- "Equilibrium Yield"

  if(xaxis == "U") {
    plot(u.vector, Yield, typ = 'l', xlab = "Exploitation rate (U)", ylab = ylab)
    segments(x0 = umsy, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass, Yield, typ = 'l', xlab = "Biomass", ylab = ylab)
    segments(x0 = BMSY, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(BKratio, Yield, typ = 'l', xlab = expression(B/B[0]), ylab = ylab)
    segments(x0 = BMSY/K, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = BMSY/K, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible()
}

#' Find the production parameter based on depletion that produces MSY
#'
#' For surplus production models, this function returns the production exponent n corresponding
#' to BMSY/K (Fletcher 1978).
#'
#' @param depletion The hypothesized depletion that produces MSY.
#' @param figure Local, plots figure of production function as a function of depletion (B/K)
#'
#' @author Q. Huynh
#' @references
#' Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson system. Fishery Bulletin 76:515:521.
#' @note May be useful for parameterizing \code{n} in \link{SP} and \link{SP_SS}.
#' @examples SP_production(0.5)
#' @return The production function exponent n (numeric).
#' @importFrom stats uniroot
#' @examples
#' SP_production(0.5)
#' @seealso \link{SP} \link{SP_SS}
#' @export SP_production
SP_production <- function(depletion, figure = TRUE) {

  if(length(depletion) > 1) {
    depletion <- depletion[1]
    message(paste("Function is not vectorized. Depletion value of", depletion, "is used."))
  }
  if(depletion <= 0 || depletion >= 1) stop(paste("Proposed depletion =", depletion, "but value must be between 0 and 1."))

  calc_depletion <- function(n) {
    depletion_MSY <- if(n==1) 1/exp(1) else n^(1/(1-n))
    return(depletion_MSY)
  }
  n_solver <- function(x) calc_depletion(x) - depletion
  get_n <- uniroot(f = n_solver, interval = c(0, 1e3))
  n_answer <- round(get_n$root, 3)

  if(figure) {
    #n_term <- if(n_answer == 1) exp(1) else n_answer^(n_answer/(n_answer-1))
    umsy <- 0.1
    msy <- umsy * depletion
    plot_yield_SP(report = list(n = n_answer, BMSY = depletion, K = 1), umsy = umsy,
                  msy = msy, xaxis = "Depletion", relative_yaxis = TRUE)
    title(paste0("Production exponent n = ", n_answer))
  }
  return(n_answer)
}
