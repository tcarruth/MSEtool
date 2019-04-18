
summary_VPA <- function(Assessment) {
  assign_Assessment_slots(Assessment)

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
  if(!info$fix_h) input_parameters <- input_parameters[-1, ]

  if(conv) Value <- c(VB0, SSB0, MSY, FMSY, VBMSY, SSBMSY, SSBMSY/SSB0)
  else Value <- rep(NA, 7)

  Description <- c("Unfished vulnerable biomass",
                  "Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Fishing mortality at MSY",
                  "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "UMSY", "VBMSY", "SSBMSY", "SSBMSY/SSB0")

  if(!is.character(SD)) {
    model_estimates <- summary(SD)
    model_estimates <- model_estimates[!is.na(model_estimates[, 2]) && model_estimates[, 2] > 0, ]
  } else model_estimates <- SD

  output <- list(model = "Virtual Population Analysis (VPA)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates)
  return(output)
}

rmd_VPA <- function(Assessment) {
  ss <- rmd_summary("Virtual Population Analysis (VPA)")

  # Life History
  age <- Assessment@info$ages
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, Assessment@info$data$mat, fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."))

  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"), rmd_data_timeseries("Index"),
                    rmd_data_age_comps("bubble", ages = vector2char(age)),
                    rmd_data_age_comps("annual", ages = vector2char(age), annual_yscale = "\"raw\"", annual_ylab = "\"Catch-at-age\""))

  # Assessment
  #### Pars and Fit
  assess_fit <- c("## Assessment {.tabset}\n### Estimates and Model Fit\n",
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Estimated selectivity at age."),
                  rmd_assess_fit("Index", "index"), rmd_assess_resid("Index"), rmd_assess_qq("Index", "index"),
                  rmd_fit_age_comps("annual", ages = vector2char(age), match = TRUE))

  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n", fig.cap = "apical fishing mortality"), rmd_F_FMSY(),
                 rmd_sel_annual(age), rmd_sel_persp(age),
                 rmd_U(fig.cap = "harvest rate (ratio of catch and vulnerable biomass)"),
                 rmd_U_UMSY(fig.cap = "U/UMSY, where UMSY = MSY/VBMSY"),
                 rmd_SSB(), rmd_SSB_SSBMSY(), rmd_SSB_SSB0(),
                 rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"),
                 rmd_R(), rmd_N(), rmd_N_at_age(), rmd_C_at_age(), rmd_C_mean_age())

  # Productivity
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  SSB <- Assessment@SSB[1:(length(Assessment@SSB)-1)]

  SR <- Assessment@info$SR
  if(SR == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB) else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }
  estR <- Assessment@R[as.numeric(names(Assessment@R)) > Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, ylab = paste0("Recruitment (age- ", min(age), ")"),
                           header = "### Productivity\n\n\n", conv_check = TRUE),
                    rmd_SR(SSB, expectedR, estR, ylab = paste0("Recruitment (age- ", min(age), ")"),
                           fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE, conv_check = TRUE),
                    rmd_yield_F("VPA"), rmd_yield_depletion("VPA"), rmd_sp(depletion = FALSE))

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}



#' @import grDevices
#' @importFrom stats qqnorm qqline
#' @importFrom graphics persp
generate_plots_VPA <- function(Assessment, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()

  plot_composition(Year, Selectivity, plot_type = 'annual', ages = age, annual_yscale = "raw", annual_ylab = "Selectivity")


  plot_composition(Year, Obs_C_at_age, C_at_age, ages = age, plot_type = 'annual')

  return(invisible())
}




plot_yield_VPA <- function(data, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion")) {
  plot_yield_SCA(data = data, report = report, fmsy = fmsy, msy = msy, xaxis = xaxis)
}


#' @importFrom reshape2 acast
profile_likelihood_VPA <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"F_term" %in% names(dots)) stop("Sequence of F_term was not found. See help file.")
  F_term <- dots$F_term

  params <- Assessment@info$params
  map <- Assessment@obj$env$map
  map$logF_term <- factor(NA)

  profile_fn <- function(i, Assessment, params, map) {
    params$logF_term <- log(F_term[i])
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
    if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:length(F_term), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid <- data.frame(F_term = F_term, nll = nll)

  if(figure) {
    plot(F_term, profile_grid$nll, typ = 'o', pch = 16, xlab = "Terminal F", ylab = "Change in negative log-likelihood")
    abline(v = Assessment@SD$value[names(Assessment@SD$value) == "F_term"], lty = 2)

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      plot(F_term, profile_grid$nll, typ = 'o', pch = 16, xlab = "Terminal F", ylab = "Change in negative log-likelihood")
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
  return(profile_grid)
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
        vars_div <- c("B", "E", "VB", "N", "CAApred")
        vars_mult <- NULL
        var_trans <- "q"
        fun_trans <- "*"
        fun_fixed <- NA
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      report <- projection_VPA(report, info, data$n_Rpen)

      SSB <- c(report$E, rep(NA, i))
      R <- c(report$N[, 1], rep(NA, i))
      N <- c(rowSums(report$N), rep(NA, i))
      FMort <- c(apply(report$F, 1, max), rep(NA, i + 1))

      Z_mat <- t(report$F) + data$M
      VB_mid <- t(report$N[-ncol(report$N), ]) * (1 - exp(-Z_mat))/Z_mat
      U <- c(colSums(t(report$CAApred) * data$weight)/colSums(VB_mid), rep(NA, i + 1))

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
