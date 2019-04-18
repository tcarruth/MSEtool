
summary_DD_TMB <- function(Assessment, state_space = FALSE) {
  assign_Assessment_slots()

  if(conv) current_status <- c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  Value <- c(unlist(info$data[c(2,3,4,6,7)]))
  Description <- c("Unfished survival = exp(-M)", "alpha = Winf * (1-rho)",
                  "rho = (W_k+2 - Winf)/(W_k+1 - Winf)",
                  "Age of knife-edge selectivity",
                  "Weight at age k")
  rownam <- c("S0", "alpha", "rho", "k", "w_k")
  if("log_sigma" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$sigma)
    Description <- c(Description, "Catch SD (log-space)")
    rownam <- c(rownam, "sigma")
  }
  if("log_tau" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$tau)
    Description <- c(Description, "log-Recruitment deviation SD")
    rownam <- c(rownam, "tau")
  }
  if("transformed_h" %in% names(obj$env$map)) {
    Value <- c(Value, h)
    Description <- c(Description, "Stock-recruit steepness")
    rownam <- c(rownam, "tau")
  }
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam

  if(conv) derived <- c(B0, N0, MSY, UMSY, BMSY, BMSY/B0)
  else derived <- rep(NA, 6)
  derived <- data.frame(Value = derived,
                        Description = c("Unfished biomass", "Unfished abundance", "Maximum sustainable yield (MSY)",
                                        "Harvest rate at MSY", "Biomass at MSY", "Depletion at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("B0", "N0", "MSY", "UMSY", "BMSY", "BMSY/B0")

  if(!is.character(SD)) {
    if(state_space) {
      if(is.null(obj$env$random)) {
        model_estimates <- summary(SD)[rownames(summary(SD)) != "log_rec_dev", ]
        dev_estimates <- summary(SD)[rownames(summary(SD)) == "log_rec_dev", ]
      } else {
        model_estimates <- rbind(summary(SD, "fixed"), summary(SD, "report"))
        dev_estimates <- summary(SD, "random")
      }
      rownames(dev_estimates) <- paste0(rownames(dev_estimates), "_", names(Dev))
      model_estimates <- rbind(model_estimates, dev_estimates)
      model_estimates <- model_estimates[!is.na(model_estimates[, 2]) && model_estimates[, 2] > 0, ]
    } else model_estimates <- summary(SD)

  } else {
    model_estimates <- SD
  }

  model_name <- "Delay-Difference"
  if(state_space) model_name <- paste(model_name, "(State-Space)")
  output <- list(model = model_name,
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates)
  return(output)
}


rmd_DD_TMB <- function(Assessment, state_space = FALSE) {
  if(state_space) {
    ss <- rmd_summary("Delay-Difference (State-Space)")
  } else ss <- rmd_summary("Delay-Difference")

  # Life History
  age <- 1:Assessment@info$LH$maxage
  k <- Assessment@info$data$k
  mat <- ifelse(age < k, 0, 1)
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, mat, fig.cap = "Assumed knife-edge maturity at age corresponding to length of 50% maturity."))

  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"), rmd_data_timeseries("Index"))

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h(),
                  rmd_sel(age, mat, fig.cap = "Knife-edge selectivity set to the age corresponding to the length of 50% maturity."),
                  rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"))

  if(state_space) {
    assess_fit2 <- c(rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type),
                     rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                  label = Assessment@Dev_type, conv_check = TRUE))
    assess_fit <- c(assess_fit, assess_fit2)
  }

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"),
                 rmd_R(), rmd_N())

  #### Productivity
  ny <- Assessment@info$data$ny
  SSB <- Assessment@SSB[1:ny]
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  if(Assessment@info$data$SR_type == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB) else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }

  first_recruit_year <- k + 1
  last_recruit_year <- length(Assessment@info$Year) + k
  ind_recruit <- first_recruit_year:last_recruit_year
  rec_dev <- Assessment@R[ind_recruit]

  productivity <- c(rmd_SR(SSB, expectedR, rec_dev, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, rec_dev, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("DD"), rmd_yield_depletion("DD"), rmd_sp())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))

}



#' @importFrom reshape2 acast
profile_likelihood_DD_TMB <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
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

        html_report(plot.dir, model = "Delay Difference", captions = matrix(profile.file.caption, nrow = 1),
                    name = Assessment@Name, report_type = "Profile_Likelihood")
        browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
      }
    }

  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_DD_TMB <- function(Assessment, nyr, figure = TRUE,
                                 save_figure = FALSE, save_dir = tempdir()) {
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

  fix_h <- "transformed_h" %in% names(obj$env$map)
  if(fix_h) est_ind <- 4 else est_ind <- 4:5
  est_lab <- "R0 estimate"
  if(!fix_h) est_lab <- c(est_lab, "Steepness estimate")

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1], retro_est[, , 1][, est_ind, drop = FALSE],
                                 ts_lab = c("Biomass", "B_BMSY", "B_B0", "Recruitment", "Abundance", "U", "U_UMSY"),
                                 est_lab = est_lab)

  if(figure) {
    plot_retro_DD_TMB(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                      nyr_label = 0:nyr, color = rich.colors(nyr+1), fix_h, data$SR_type)
  }

  return(Mohn_rho)
}

plot_retro_DD_TMB <- function(retro_ts, retro_est, save_figure = FALSE,
                              save_dir = tempdir(), nyr_label, color, fix_h, SR) {
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


plot_yield_DD <- function(data, report, umsy, msy, xaxis = c("U", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  u.vector <- seq(0, 1, 0.01)
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
  #if(SR_type == "BH") R <- (Arec * BPR * (1 - u.vector)- 1)/(Brec * BPR * (1 - u.vector))
  #if(SR_type == "Ricker") R <- log(Arec * BPR * (1 - u.vector))/(Brec * BPR * (1 - u.vector))
  if(SR_type == "BH") R <- (Arec * BPR - 1)/(Brec * BPR)
  if(SR_type == "Ricker") R <- log(Arec * BPR)/(Brec * BPR)

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
  invisible(data.frame(U = u.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/report$B0))
}
