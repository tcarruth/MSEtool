
summary_SCA2 <- function(Assessment) summary_SCA(Assessment, TRUE)

rmd_SCA2 <- function(Assessment) rmd_SCA(Assessment, TRUE)


#' @importFrom reshape2 acast
profile_likelihood_SCA2 <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"meanR" %in% names(dots)) stop("Sequence of meanR was not found. See help file.")
  meanR <- dots$meanR

  params <- Assessment@info$params
  map <- Assessment@obj$env$map
  map$log_meanR <- factor(NA)

  profile_fn <- function(i, Assessment, params, map) {
    params$log_meanR <- log(meanR[i] * Assessment@info$rescale)
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map,
                      random = Assessment@obj$env$random, inner.control = Assessment@info$inner.control,
                      DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
    if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:length(meanR), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid <- data.frame(meanR = meanR, nll = nll)

  if(figure) {
    plot(profile_grid$meanR, profile_grid$nll, typ = 'o', pch = 16, xlab = "Mean recruitment", ylab = "Change in negative log-likelihood")
    abline(v = Assessment@SD$value[names(Assessment@SD$value) == "meanR"], lty = 2)

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      plot(profile_grid$meanR, profile_grid$nll, typ = 'o', pch = 16, xlab = "Mean recruitment", ylab = "Change in negative log-likelihood")
      abline(v = Assessment@SD$value[names(Assessment@SD$value) == "meanR"], lty = 2)
      dev.off()
      profile.file.caption <- c("profile_likelihood.png",
                                "Profile likelihood of mean recruitment. Vertical, dashed line indicates maximum likelihood estimate.")

      html_report(plot.dir, model = "Statistical Catch-at-Age (SCA2)",
                  captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile_grid)
}


#' @importFrom gplots rich.colors
retrospective_SCA2 <- function(Assessment, nyr, figure = TRUE, save_figure = FALSE, save_dir = tempdir()) {
  assign_Assessment_slots()
  data <- info$data
  n_y <- data$n_y

  Year <- c(info$Year, max(info$Year) + 1)
  C_hist <- data$C_hist
  I_hist <- data$I_hist
  CAA_hist <- data$CAA_hist
  CAA_n <- data$CAA_n
  params <- info$params

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar Year, SSB, SSB_SSBMSY, SSB_SSB0, N, R, F, F_FMSY, log_rec_dev
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 9))
  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev" &
                             rownames(summary(SD)) != "logF", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))

  SD <- NULL
  rescale <- info$rescale
  fix_h <- ifelse(is.null(info$h), FALSE, TRUE)

  for(i in 0:nyr) {
    n_y_ret <- n_y - i
    data$n_y <- n_y_ret
    data$C_hist <- C_hist[1:n_y_ret]
    data$I_hist <- I_hist[1:n_y_ret]
    data$CAA_hist <- CAA_hist[1:n_y_ret, ]
    data$CAA_n <- CAA_n[1:n_y_ret]
    params$log_rec_dev <- rep(0, n_y_ret)

    map <- obj$env$map
    new_map <- as.numeric(map$log_rec_dev) - i
    map$log_rec_dev <- factor(new_map[new_map > 0])

    obj2 <- MakeADFun(data = data, parameters = params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      refpt <- SCA_refpt_calc(E = report$E[1:(length(report$E) - 1)], R = report$R[2:length(report$R)],
                              weight = data$weight, mat = data$mat, M = data$M, vul = report$vul, SR = info$SR, fix_h = fix_h, h = info$h)
      report <- c(report, refpt)
      if(info$rescale != 1) {
        vars_div <- c("meanR", "B", "E", "CAApred", "Cpred", "CN", "N", "VB",
                      "R", "MSY", "VBMSY", "RMSY", "BMSY", "EMSY", "VB0", "R0",
                      "B0", "E0", "N0")
        vars_mult <- "Brec"
        var_trans <- c("meanR", "q")
        fun_trans <- c("/", "*")
        fun_fixed <- c("log", NA)
        rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
      }

      SSB <- c(report$E, rep(NA, i))
      SSB_SSBMSY <- SSB/report$EMSY
      SSB_SSB0 <- SSB/report$E0
      R <- c(report$R, rep(NA, i))
      N <- c(rowSums(report$N), rep(NA, i))
      FMort <- c(report$F, rep(NA, i + 1))
      F_FMSY <- FMort/report$FMSY
      log_rec_dev <- c(report$log_rec_dev, rep(NA, i + 1))

      retro_ts[i+1, , ] <- cbind(Year, SSB, SSB_SSBMSY, SSB_SSB0, R, N, FMort, F_FMSY, log_rec_dev)
      retro_est[i+1, , ] <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev" &
                                          rownames(summary(SD)) != "logF", ]

    } else {
      message(paste("Non-convergence when", i, "years of data were removed."))
    }
  }

  Mohn_rho <- calculate_Mohn_rho(retro_ts[, , -1],
                                 ts_lab = c("SSB", "SSB_SSBMSY", "SSB_SSB0", "Recruitment", "Abundance", "F", "F_FMSY", "log_rec_dev"))

  if(figure) {
    plot_retro_SCA2(retro_ts, retro_est, save_figure = save_figure, save_dir = save_dir,
                    nyr_label = 0:nyr, color = rich.colors(nyr+1))
  }

  return(Mohn_rho)
}


plot_retro_SCA2 <- function(retro_ts, retro_est, save_figure = FALSE, save_dir = tempdir(), nyr_label, color) {
  n_tsplots <- dim(retro_ts)[3] - 1
  ts_label <- c("Spawning Stock Biomass", expression(SSB/SSB[MSY]), expression(SSB/SSB[0]), "Recruitment",
                "Population Abundance (N)", "Fishing Mortality", expression(F/F[MSY]), "Recruitment deviations")
  Year <- retro_ts[1, , 1]

  if(save_figure) {
    Model <- "SCA2"
    prepare_to_save_figure()
  }

  for(i in 1:n_tsplots) {
    y.max <- max(abs(retro_ts[, , i+1]), na.rm = TRUE)
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

  if(save_figure) {
    ret.file.caption <- data.frame(x1 = paste0("retrospective_", c(1:n_tsplots), ".png"),
                                   x2 = paste0("Retrospective pattern in ",
                                               c("spawning stock biomass", "SSB/SSBMSY", "spawning depletion", "recruitment",
                                                 "abundance", "fishing mortality", "F/FMSY", "recruitment deviations"), "."))
    Assessment <- get("Assessment", envir = parent.frame())
    html_report(plot.dir, model = "Statistical Catch-at-Age (SCA2)", captions = ret.file.caption,
                name = Assessment@Name, report_type = "Retrospective")
    browseURL(file.path(plot.dir, "Retrospective.html"))
  }

  invisible()
}

