
summary_DD_SS <- function(Assessment) summary_DD_TMB(Assessment, TRUE)

rmd_DD_SS <- function(Assessment) rmd_DD_TMB(Assessment, TRUE)



#' @importFrom reshape2 acast
profile_likelihood_DD_SS <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
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
    params$log_R0 = log(profile.grid$R0[i] * Assessment@info$rescale)
    params$transformed_h <- transformed_h[i]
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params,
                      map = map, random = random, inner.control = Assessment@info$inner.control,
                      DLL = "MSEtool", silent = TRUE)
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

        html_report(plot.dir, model = "Delay Difference (State-Space)",
                    captions = matrix(profile.file.caption, nrow = 1),
                    name = Assessment@Name, report_type = "Profile_Likelihood")
        browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
      }
    }
  }
  return(profile.grid)
}


#' @importFrom gplots rich.colors
retrospective_DD_SS <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()
  data <- info$data
  ny <- data$ny
  k <- data$k

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k
  Year <- c(Year, moreRecruitYears)
  C_hist <- data$C_hist
  E_hist <- data$E_hist
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, B, B_BMSY, B_B0, N, R, U, U_UMSY, log_rec_dev
  retro_ts <- array(NA, dim = c(nyr+1, ny + k, 9))
  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_rec_dev", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))

  SD <- NULL
  rescale <- info$rescale

  for(i in 0:nyr) {
    ny_ret <- ny - i
    C_hist_ret <- C_hist[1:ny_ret]
    E_hist_ret <- E_hist[1:ny_ret]
    data$ny <- ny_ret
    data$C_hist <- C_hist_ret
    data$E_hist <- E_hist_ret
    params$log_rec_dev <- rep(0, ny_ret - k)

    obj2 <- MakeADFun(data = data, parameters = params, map = obj$env$map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
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
      B_B0 <- B/report$B0
      R <- c(report$R, rep(NA, i))
      N <- c(report$N, rep(NA, k - 1 + i))
      U <- c(report$U, rep(NA, k + i))
      U_UMSY <- U/report$UMSY
      log_rec_dev <- c(report$log_rec_dev, rep(NA, 2 * k + i))

      retro_ts[i+1, , ] <- cbind(Year, B, B_BMSY, B_B0, R, N, U, U_UMSY, log_rec_dev)
      retro_est[i+1, , ] <- summary(SD)[rownames(summary(SD)) != "log_rec_dev", ]

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }
  }

  fix_h <- "transformed_h" %in% names(obj$env$map)
  if(fix_h) est_ind <- 5 else est_ind <- 5:6
  est_lab <- "R0 estimate"
  if(!fix_h) est_lab <- c(est_lab, "Steepness estimate")

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1], retro_est[, , 1][, est_ind, drop = FALSE],
                                 ts_lab = c("Biomass", "B_BMSY", "B_B0", "Recruitment", "Abundance", "U", "U_UMSY", "log_rec_dev"),
                                 est_lab = est_lab)

  if(figure) {
    plot_retro_DD_SS(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                      nyr_label = 0:nyr, color = rich.colors(nyr+1), fix_h, data$SR_type)
  }

  return(Mohn_rho)
}


plot_retro_DD_SS <- function(retro_ts, retro_est, save_figure = FALSE,
                             save_dir = tempdir(), nyr_label, color, fix_h, SR) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Biomass", expression(B/B[MSY]), expression(B/B[0]), "Recruitment",
                "Population Abundance (N)", "Exploitation rate (U)",
                expression(U/U[MSY]), "Recruitment deviations")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "DD_SS"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(retro_ts[, , i+1], na.rm = TRUE)
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
                                               c("biomass", "B/BMSY", "biomass depletion", "recruitment",
                                                 "abundance", "exploitation", "U/UMSY", "recruitment deviations", "R0 estimate"), "."),
                                   stringsAsFactors = FALSE)
    if(!fix_h) {
      ret.file.caption <- rbind(ret.file.caption,
                                c(paste0("retrospective_", n_tsplots+2, ".png"), "Retrospective pattern in steepness estimate."))
    }

    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Delay Difference (State-Space)", captions = ret.file.caption,
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}
