summary_SP_SS <- function(Assessment) summary_SP(Assessment, TRUE)

rmd_SP_SS <- function(Assessment) rmd_SP(Assessment, TRUE)




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



