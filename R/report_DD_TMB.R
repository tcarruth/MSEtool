
summary_DD_TMB <- function(Assessment) {
  assign_Assessment_slots()

  if(conv) current_status <- c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  input_parameters <- data.frame(Value = as.numeric(c(h, unlist(info$data[c(2,3,4,6,7)]))),
                                 Description = c("Stock-recruit steepness", "Unfished survival = exp(-M)", "alpha = Winf * (1-rho)",
                                                 "rho = (W_k+2 - Winf)/(W_k+1 - Winf)", "Age of knife-edge selectivity",
                                                 "Weight at age k"),
                                 stringsAsFactors = FALSE)
  rownames(input_parameters) <- c("h", "S0", "alpha", "rho", "k", "w_k")
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  if(conv) derived <- c(B0, N0, MSY, UMSY, BMSY)
  else derived <- rep(NA, 5)
  derived <- data.frame(Value = derived,
                        Description = c("Virgin biomass", "Virgin abundance", "Maximum sustainable yield (MSY)",
                                        "Harvest Rate at MSY", "Biomass at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("B0", "N0", "MSY", "UMSY", "BMSY")

  if(conv) {
    model_estimates <- summary(SD)
    model_estimates <- model_estimates[model_estimates[, 2] > 0, ]
  } else {
    model_estimates <- SD
  }

  output <- list(model = "Delay-Difference",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
generate_plots_DD_TMB <- function(Assessment, save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary(Assessment)
    html_report(plot.dir, model = "Delay Difference",
                current_status = index.report$current_status,
                input_parameters = index.report$input_parameters,
                derived_quantities = index.report$derived_quantities,
                model_estimates = index.report$model_estimates,
                name = Name, report_type = "Index")
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
    html_report(plot.dir, model = "Delay Difference", captions = lh.file.caption,
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

  #if(!is.na(Data@CV_Cat[1]) && sdconv(1, Data@CV_Cat[1]) > 0.01) {
  #  plot_timeseries(as.numeric(names(Obs_Catch)), Obs_Catch, obs_CV = Data@CV_Cat[1], label = "Catch")
  #  if(save_figure) {
  #    create_png(filename = file.path(plot.dir, "data_catch_with_CI.png"))
  #    plot_timeseries(as.numeric(names(Obs_Catch)), Obs_Catch, obs_CV = Data@CV_Cat[1], label = "Catch")
  #    dev.off()
  #    data.file.caption <- rbind(data.file.caption,
  #                               c("data_catch_with_CI.png", "Catch time series with 95% confidence interval."))
  #  }
  #}

  plot_timeseries(Year, info$I_hist, label = "Index")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_index.png"))
    plot_timeseries(Year, info$I_hist, label = "Index")
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_index.png", "Index time series."))
  }

  #if(!is.na(Data@CV_Ind[1]) && sdconv(1, Data@CV_Ind[1]) > 0.01) {
  #  plot_timeseries(Year, info$I_hist, obs_CV = Data@CV_Ind[1], label = "Index")
  #  if(save_figure) {
  #    create_png(filename = file.path(plot.dir, "data_index_with_CI.png"))
  #    plot_timeseries(Year, info$I_hist, obs_CV = Data@CV_Ind[1], label = "Index")
  #    dev.off()
  #    data.file.caption <- rbind(data.file.caption,
  #                               c("data_index_with_CI.png", "Index time series with 95% confidence interval."))
  #  }
  #}

  if(save_figure) {
    html_report(plot.dir, model = "Delay Difference", captions = data.file.caption,
                name = Name, report_type = "Data")
  }

  if(conv) {
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
                                 c("assessment_catch_residual.png", "Catch residuals in log-space."))
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

  ny <- info$data$ny
  SSB <- B[1:ny] - Catch # B*(1-u)
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

  plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_exploitation.png"))
    plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_exploitation.png", "Time series of exploitation rate."))
  }

  if(conv) {
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
    html_report(plot.dir, model = "Delay Difference", captions = assess.file.caption,
                name = Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }

  return(invisible())
}


#' @importFrom reshape2 acast
profile_likelihood_DD_TMB <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = getwd(), ...) {
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
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params,
                      map = map, DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
    if(!is.character(opt2)) nll[i] <- opt2$objective
  }
  profile.grid$nll <- nll - Assessment@opt$objective

  if(figure) {
    if(length(h) > 1) {
      z.mat <- acast(profile.grid, list("h", "R0"), value.var = "nll")
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
      }
    }
    html_report(plot.dir, model = "Delay Difference",
                captions = matrix(profile.file.caption, nrow = 1),
                name = Assessment@Name, report_type = "Profile_Likelihood")
    browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_DD_TMB <- function(Assessment, nyr, figure = TRUE,
                                 save_figure = FALSE, save_dir = getwd()) {
  assign_Assessment_slots()
  data <- info$data
  ny <- data$ny
  k <- data$k

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k
  Year <- c(Year, moreRecruitYears)
  C_hist <- data$C_hist
  E_hist <- data$E_hist

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, B/BMSY, B/B0, N, R, U, U/UMSY
  retro_ts <- array(NA, dim = c(nyr+1, ny + k, 8))
  retro_est <- array(NA, dim = c(nyr+1, dim(summary(SD))))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    ny_ret <- ny - i
    C_hist_ret <- C_hist[1:ny_ret]
    E_hist_ret <- E_hist[1:ny_ret]
    data$ny <- ny_ret
    data$C_hist <- C_hist_ret
    data$E_hist <- E_hist_ret

    obj2 <- MakeADFun(data = data, parameters = info$params, map = obj$env$map, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- get_MSY_DD(info$data, report$Arec, report$Brec)
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
      U <- c(report$U, rep(NA, k + i))
      U_UMSY <- U/report$UMSY

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, R, N, U, U_UMSY)
      if(!is.character(SD)) retro_est[i+1, , ] <- summary(SD)

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }

  }
  if(figure) {
    fix_h <- "transformed_h" %in% names(obj$env$map)
    plot_retro_DD_TMB(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                      nyr_label = 0:nyr, color = rich.colors(nyr+1), fix_h, data$SR_type)
  }
  # Need to write legend
  legend <- NULL
  return(invisible(list(legend = legend, retro_ts = retro_ts, retro_est = retro_est)))
}

plot_retro_DD_TMB <- function(retro_ts, retro_est, save_figure = FALSE,
                              save_dir = getwd(), nyr_label, color, fix_h, SR) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]), "Recruitment",
                "Population Abundance (N)", "Exploitation rate (U)",
                expression(U/U[MSY]))
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "DD_TMB"
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
                                                 "abundance", "exploitation", "U/UMSY", "R0 estimate"), "."),
                                   stringsAsFactors = FALSE)
    if(!fix_h) {
      ret.file.caption <- rbind(ret.file.caption,
                                c(paste0("retrospective_", n_tsplots+2, ".png"), "Retrospective pattern in steepness estimate."))
    }

    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Delay Difference", captions = ret.file.caption,
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}


plot_yield_DD <- function(data, report, umsy, msy, u.vector = seq(0, 1, 0.01),
                          xaxis = c("U", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  S0 <- data$S0
  Alpha <- data$Alpha
  wk <- data$wk
  Rho <- data$Rho
  SR_type  <- data$SR_type

  Arec <- report$Arec
  Brec <- report$Brec
  BMSY <- report$BMSY

  Surv <- S0 * (1 - u.vector)

  BPR <- (Surv * Alpha/(1 - Surv) + wk)/(1 - Rho * Surv)
  if(SR_type == "BH") R <- (Arec * BPR * (1 - u.vector)- 1)/(Brec * BPR * (1 - u.vector))
  if(SR_type == "Ricker") R <- log(Arec * BPR * (1 - u.vector))/(Brec * BPR * (1 - u.vector))

  Biomass <- BPR * R
  Yield <- u.vector * BPR * R
  ind <- R >= 0

  if(xaxis == "U") {
    plot(u.vector[ind], Yield[ind], typ = 'l', xlab = "Exploitation rate (U)",
         ylab = "Equilibrium yield")
    segments(x0 = umsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/report$B0, Yield[ind], typ = 'l',
         xlab = expression(B/B[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/report$B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/report$B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible()
}
