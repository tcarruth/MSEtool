
summary_VPA <- function(Assessment) {
  assign_Assessment_slots()

  if(conv) current_status <- c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- rep(NA, 3)
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  Value <- c(h, info$data$M[1], min(info$ages), max(info$ages), info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95)
  Description = c("Stock-recruit steepness", "Natural mortality", "Minimum age (minus-group)", "Maximum age (plus-group)", "Asymptotic length",
                  "Growth coefficient", "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  rownam <- c("h", "M", "minage", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam

  if(conv) Value <- c(VB0, SSB0, MSY, FMSY, VBMSY, SSBMSY)
  else Value <- rep(NA, 6)

  Description <- c("Virgin vulnerable biomass",
                  "Virgin spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Fishing mortality at MSY",
                  "Vulnerable biomass at MSY", "SSB at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "UMSY", "VBMSY", "SSBMSY")

  if(!is.character(SD)) {
    model_estimates <- summary(SD)
    model_estimates <- model_estimates[is.na(model_estimates[, 2]) || model_estimates[, 2] > 0, ]
  } else model_estimates <- SD

  output <- list(model = "Virtual Population Analysis (VPA)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

#' @import grDevices
#' @importFrom stats qqnorm qqline
#' @importFrom graphics persp
generate_plots_VPA <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  if(save_figure) {
    prepare_to_save_figure()
    index.report <- summary(Assessment)
    html_report(plot.dir, model = "Virtual Population Analysis (VPA)",
                current_status = index.report$current_status,
                input_parameters = index.report$input_parameters,
                model_estimates = index.report$model_estimates,
                derived_quantities = index.report$derived_quantities,
                name = Name, report_type = "Index")
  }

  age <- info$ages

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
    html_report(plot.dir, model = "Virtual Population Analysis (VPA)",
                captions = lh.file.caption, name = Name, report_type = "Life_History")
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

  plot_composition(Year, Obs_C_at_age, ages = age, plot_type = 'bubble_data')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "data_age_comps_bubble.png"))
    plot_composition(Year, Obs_C_at_age, ages = age, plot_type = 'bubble_data')
    dev.off()
    data.file.caption <- rbind(data.file.caption,
                               c("data_age_comps_bubble.png", "Catch-at-age bubble plot."))
  }

  plot_composition(Year, Obs_C_at_age, plot_type = 'annual', ages = age, annual_yscale = "raw", annual_ylab = "Catch-at-age")
  if(save_figure) {
    nplots <- ceiling(length(Year)/16)
    for(i in 1:nplots) {
      ind <- (16*(i-1)+1):(16*i)
      if(i == nplots) ind <- (16*(i-1)+1):length(Year)

      create_png(filename = file.path(plot.dir, paste0("data_age_comps_", i, ".png")))
      plot_composition(Year, Obs_C_at_age, plot_type = 'annual', ages = age, annual_yscale = "raw", annual_ylab = "Catch-at-age", ind = ind)
      dev.off()
      data.file.caption <- rbind(data.file.caption,
                                 c(paste0("data_age_comps_", i, ".png"), paste0("Annual catch-at-age (", i, "/", nplots, ")")))
    }
  }

  if(save_figure) {
    html_report(plot.dir, model = "Virtual Population Analysis (VPA)",
                captions = data.file.caption, name = Name, report_type = "Data")
  }


  plot_timeseries(as.numeric(names(FMort)), FMort, label = "Apical Fishing Mortality (F)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_F.png"))
    plot_timeseries(as.numeric(names(FMort)), FMort, label = "Apical Fishing Mortality (F)")
    dev.off()
    assess.file.caption <- c("assessment_F.png", "Time series of apical fishing mortality.")
  }

  plot_composition(Year, Selectivity, plot_type = 'annual', ages = age, annual_yscale = "raw", annual_ylab = "Selectivity")
  if(save_figure) {
    nplots <- ceiling(length(Year)/16)
    for(i in 1:nplots) {
      ind <- (16*(i-1)+1):(16*i)
      if(i == nplots) ind <- (16*(i-1)+1):length(Year)

      create_png(filename = file.path(plot.dir, paste0("assessment_selectivity_", i, ".png")))
      plot_composition(Year, Selectivity, plot_type = 'annual', ages = age, annual_yscale = "raw", annual_ylab = "Selectivity", ind = ind)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                 c(paste0("assessment_selectivity_", i, ".png"), paste0("Annual selectivity (", i, "/", nplots, ")")))
    }
  }

  persp(Year, age, Selectivity, expand = 0.35, ticktype = 'detailed', phi = 25,
        theta = 45, xlab = "Year", ylab = "Age", zlab = "Selectivity")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, paste0("assessment_selectivity_persp.png")))
    persp(Year, age, Selectivity, expand = 0.35, ticktype = 'detailed', phi = 25,
          theta = 45, xlab = "Year", ylab = "Age", zlab = "Selectivity")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c(paste0("assessment_selectivity_persp.png"), "Perspective plot of selectivity."))
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

    plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_exploitation.png"))
      plot_timeseries(as.numeric(names(U)), U, label = "Exploitation rate (U)")
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_exploitation.png", "Time series of exploitation rate (ratio of catch and vulnerable biomass)."))
    }

    plot_timeseries(as.numeric(names(U_UMSY)), U_UMSY, label = expression(U/U[MSY]))
    abline(h = 1, lty = 2)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_U_UMSY.png"))
      plot_timeseries(as.numeric(names(U_UMSY)), U_UMSY, label = expression(U/U[MSY]))
      abline(h = 1, lty = 2)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_U_UMSY.png", "Time series of U/UMSY, where UMSY = 1 - exp(-FMSY)."))
    }

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

  plot_composition(Year, Obs_C_at_age, C_at_age, N = rowSums(Obs_C_at_age), ages = age, plot_type = 'annual')
  if(save_figure) {
    nplots <- ceiling(length(Year)/16)
    for(i in 1:nplots) {
      ind <- (16*(i-1)+1):(16*i)
      if(i == nplots) ind <- (16*(i-1)+1):length(Year)

      create_png(filename = file.path(plot.dir, paste0("assess_age_comps_", i, ".png")))
      plot_composition(Year, Obs_C_at_age, C_at_age, N = rowSums(Obs_C_at_age), ages = age, plot_type = 'annual', ind = ind)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c(paste0("assess_age_comps_", i, ".png"), paste0("Annual observed (black) and predicted (red) catch-at-age (",
                                                                                    i, "/", nplots, ")")))
    }
  }

  plot_composition(Year, Obs_C_at_age, C_at_age, plot_type = 'mean')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_mean_age.png"))
    plot_composition(Year, Obs_C_at_age, C_at_age, plot_type = 'mean')
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_mean_age.png", "Observed (black) and predicted (red) mean age of the catch-at-age."))
  }

  if(conv) {

    Arec <- TMB_report$Arec
    Brec <- TMB_report$Brec
    SSB_plot <- SSB[1:(length(SSB)-1)]
    if(info$SR == "BH") expectedR <- Arec * SSB_plot / (1 + Brec * SSB_plot)
    if(info$SR == "Ricker") expectedR <- Arec * SSB_plot * exp(-Brec * SSB_plot)
    estR <- R[as.numeric(names(R)) > Year[1]]

    plot_SR(SSB_plot, expectedR, R0, SSB0, estR, ylab = paste0("Recruitment (age-", min(age), ")"))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_stock_recruit.png"))
      plot_SR(SSB_plot, expectedR, R0, SSB0, estR, ylab = paste0("Recruitment (age-", min(age), ")"))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_stock_recruit.png", "Stock-recruitment relationship."))
    }

    if(max(estR) > 3 * max(expectedR)) {
      y_zoom <- 3
      plot_SR(SSB_plot, expectedR, R0, SSB0, estR, y_zoom = y_zoom, ylab = paste0("Recruitment (age-", min(age), ")"))
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "assessment_stock_recruit_zoomed.png"))
        plot_SR(SSB_plot, expectedR, R0, SSB0, estR, y_zoom = y_zoom, ylab = paste0("Recruitment (age-", min(age), ")"))
        dev.off()
        assess.file.caption <- rbind(assess.file.caption,
                                     c("assessment_stock_recruit_zoomed.png", "Stock-recruitment relationship (zoomed in)."))
      }
    } else y_zoom <- NULL

    plot_SR(SSB_plot, expectedR, R0, SSB0, estR, trajectory = TRUE, y_zoom = y_zoom, ylab = paste0("Recruitment (age-", min(age), ")"))
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_stock_recruit_trajectory.png"))
      plot_SR(SSB_plot, expectedR, R0, SSB0, estR, trajectory = TRUE, y_zoom = y_zoom, ylab = paste0("Recruitment (age-", min(age), ")"))
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_stock_recruit_trajectory.png", "Stock-recruitment relationship (trajectory plot)."))

    }

  }

  plot_timeseries(as.numeric(names(SSB)), SSB, label = "Spawning Stock Biomass (SSB)")
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_spawning_biomass.png"))
    plot_timeseries(as.numeric(names(SSB)), SSB, label = "Spawning Stock Biomass (SSB)")
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                                 c("assessment_spawning_biomass.png", "Time series of spawning stock biomass."))
  }

  if(conv) {

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
  }


  plot_timeseries(as.numeric(names(R)), R, label = paste0("Recruitment (age-", min(age), ")"))
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_recruitment.png"))
    plot_timeseries(as.numeric(names(R)), R, label = paste0("Recruitment (age-", min(age), ")"))
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

  plot_composition(c(Year, max(Year) + 1), N_at_age, ages = age, plot_type = 'bubble_data')
  if(save_figure) {
    create_png(filename = file.path(plot.dir, "assessment_abundance_at_age_bubble.png"))
    plot_composition(c(Year, max(Year) + 1), N_at_age, ages = age, plot_type = 'bubble_data')
    dev.off()
    assess.file.caption <- rbind(assess.file.caption,
                               c("assessment_abundance_at_age_bubble.png", "Abundance at age bubble plot."))
  }

  if(conv) {

    plot_Kobe(B_BMSY, U_UMSY)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_Kobe.png"))
      plot_Kobe(B_BMSY, U_UMSY)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_Kobe.png", "Kobe plot trajectory of stock."))
    }

    plot_yield_VPA(info$data, TMB_report, info$vul_refpt, FMSY, MSY, xaxis = "F", SR = info$SR)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_yield_curve_F.png"))
      plot_yield_VPA(info$data, TMB_report, info$vul_refpt, FMSY, MSY, xaxis = "F", SR = info$SR)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_yield_curve_F.png", "Yield plot relative to fishing mortality."))
    }

    plot_yield_VPA(info$data, TMB_report, info$vul_refpt, FMSY, MSY, xaxis = "Depletion", SR = info$SR)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_yield_curve_SSB_SSB0.png"))
      plot_yield_VPA(info$data, TMB_report, info$vul_refpt, FMSY, MSY, xaxis = "Depletion", SR = info$SR)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_yield_curve_SSB_SSB0.png", "Yield plot relative to spawning depletion."))
    }

    plot_surplus_production(B, B0, Catch)
    if(save_figure) {
      create_png(filename = file.path(plot.dir, "assessment_surplus_production.png"))
      plot_surplus_production(B, B0, Catch)
      dev.off()
      assess.file.caption <- rbind(assess.file.caption,
                                   c("assessment_surplus_production.png", "Surplus production relative to depletion (total biomass)."))
    }
  }

  if(save_figure) {
    html_report(plot.dir, model = "Virtual Population Analysis (VPA)",
                captions = assess.file.caption, name = Name, report_type = "Assessment")
    browseURL(file.path(plot.dir, "Assessment.html"))
  }

  return(invisible())
}





plot_yield_VPA <- function(data, report, vul, fmsy, msy, f.vector = seq(0, 2, 0.01), SR, xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)

  M <- data$M
  mat <- data$mat
  weight <- data$weight
  maxage <- data$max_age

  BMSY <- report$SSBMSY
  B0 <- report$SSB0

  Arec <- report$Arec
  Brec <- report$Brec

  BPR <- Req <- NA
  solveMSY <- function(logF) {
    Fmort <- exp(logF)
    surv <- exp(-vul * Fmort - M)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    BPR <<- sum(NPR * mat * weight)
    if(SR == "BH") Req <<- (Arec * BPR - 1)/(Brec * BPR)
    if(SR == "Ricker") Req <<- log(Arec * BPR)/(Brec * BPR)
    CPR <- vul * Fmort/(vul * Fmort + M) * NPR * (1 - exp(-vul * Fmort - M))
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  Biomass <- Yield <- R <- rep(NA, length(f.vector))
  for(i in 1:length(f.vector)) {
    Yield[i] <- -1 * solveMSY(log(f.vector[i]))
    R[i] <- Req
    Biomass[i] <- BPR * Req
  }

  ind <- R >= 0

  if(xaxis == "F") {
    plot(f.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality (F)",
         ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Spawning Stock Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/B0, Yield[ind], typ = 'l',
         xlab = expression(SSB/SSB[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible()
}


#' @importFrom reshape2 acast
profile_likelihood_VPA <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"F_term" %in% names(dots)) stop("Sequence of F_term was not found. See help file.")
  F_term <- dots$F_term

  nll <- rep(NA, length(F_term))
  params <- Assessment@info$params
  map <- Assessment@obj$env$map
  map$logF_term <- factor(NA)
  for(i in 1:length(F_term)) {
    params$logF_term <- log(F_term[i])
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]

    if(!is.character(opt2)) nll[i] <- opt2$objective
  }
  profile.grid <- data.frame(F_term = F_term, nll = nll - Assessment@opt$objective)
  if(figure) {
    plot(F_term, profile.grid$nll, typ = 'o', pch = 16, xlab = "Terminal F", ylab = "Change in negative log-likelihood")
    abline(v = Assessment@SD$value[names(Assessment@SD$value) == "F_term"], lty = 2)

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      plot(F_term, profile.grid$nll, typ = 'o', pch = 16, xlab = "Terminal F", ylab = "Change in negative log-likelihood")
      abline(v = Assessment@SD$value[names(Assessment@SD$value) == "F_term"], lty = 2)
      dev.off()
      profile.file.caption <- c("profile_likelihood.png",
                                "Profile likelihood of terminal F. Vertical, dashed line indicates maximum likelihood estimate.")

      html_report(plot.dir, model = "Virtual Population Analysis (VPA)",
                  captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_VPA <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()
  data <- info$data
  n_y <- data$n_y

  Year <- c(info$Year, max(info$Year) + 1)
  I_hist <- data$I_hist
  CAA_hist <- data$CAA_hist
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, SSB, N, R, U, U_UMSY
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 6))
  retro_est <- array(NA, dim = c(nyr+1, dim(summary(SD))))

  SD <- NULL
  rescale <- info$rescale
  fix_h <- ifelse(is.null(info$h), FALSE, TRUE)

  for(i in 0:nyr) {
    n_y_ret <- n_y - i
    data$n_y <- n_y_ret
    data$I_hist <- I_hist[1:n_y_ret]
    data$CAA_hist <- CAA_hist[1:n_y_ret, ]

    map <- obj$env$map

    obj2 <- MakeADFun(data = data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {

      report <- obj2$report(obj2$env$last.par.best)
      if(rescale != 1) {
        vars_div <- c("B", "SSB", "VB", "N", "CAApred")
        vars_mult <- NULL
        var_trans <- "q"
        fun_trans <- "*"
        fun_fixed <- NA
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      report <- projection_VPA(report, info, data$n_Rpen)

      SSB <- c(report$SSB, rep(NA, i))
      R <- c(report$N[, 1], rep(NA, i))
      N <- c(rowSums(report$N), rep(NA, i))
      FMort <- c(apply(report$F, 1, max), rep(NA, i + 1))
      U <- c(colSums(t(report$CAApred) * data$weight)/report$VB[1:(length(report$VB)-1)], rep(NA, i + 1))

      retro_ts[i+1, , ] <- cbind(Year, SSB, R, N, FMort, U)
      retro_est[i+1, , ] <- summary(SD)

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }
  }

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1],
                                 ts_lab = c("SSB", "Recruitment", "Abundance", "Apical F", "U"))

  if(figure) {
    plot_retro_VPA(retro_ts, save_figure = save_figure, save_dir = save_dir,
                   nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }

  return(Mohn_rho)
}



plot_retro_VPA <- function(retro_ts, save_figure = FALSE, save_dir = tempdir(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Spawning Stock Biomass", "Recruitment",
                "Population Abundance (N)", "Apical Fishing Mortality (F)", "Exploitation rate (U)")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "VPA"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(abs(retro_ts[, , i+1]), na.rm = TRUE)
    ylim <- c(0, 1.1 * y.max)
    plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i], ylim = ylim, col = color[1])
    for(j in 2:length(nyr_label)) lines(Year, retro_ts[j, , i+1], col = color[j])
    legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
    abline(h = 0, col = 'grey')

    if(save_figure) {
      create_png(filename = file.path(plot.dir, paste0("retrospective_", i, ".png")))
      plot(Year, retro_ts[1, , i+1], typ = 'l', ylab = ts_label[i], ylim = ylim, col = color[1])
      for(j in 2:length(nyr_label)) lines(Year, retro_ts[j, , i+1], col = color[j])
      legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
      abline(h = 0, col = 'grey')
      dev.off()
    }
  }

  if(save_figure) {
    ret.file.caption <- data.frame(x1 = paste0("retrospective_", c(1:n_tsplots), ".png"),
                                   x2 = paste0("Retrospective pattern in ",
                                               c("spawning stock biomass", "recruitment",
                                                 "abundance", "apical fishing mortality", "exploitation"), "."))
    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Virtual Population Analysis (VPA)", captions = ret.file.caption,
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}
