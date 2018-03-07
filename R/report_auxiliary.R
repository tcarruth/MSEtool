# Call from inside MP
return_Assessment <- function() {
  Call.history <- sys.calls()
  lCall <- length(Call.history)
  Model <- as.character(Call.history[[lCall-1]])[1]

  output <- mget(c('Data', 'info', 'obj', 'opt', 'SD', 'dependencies'),
                 envir = parent.frame())
  report <- output$obj$report()

  if(Model == "DD_TMB") {
    Year <- output$info$Year
    Yearplusone <- c(Year, max(Year) + 1)
    k <- output$info$data$k_DD
    Yearplusk <- c(Year, (max(Year)+1):(max(Year)+k))
    Assessment <- new("Assessment", Model = Model,
                      MSY = report$MSY_DD, UMSY = report$UMSY_DD, BMSY = report$BMSY_DD,
                      B0 = report$Bo_DD, R0 = report$Ro_DD, N0 = report$No_DD,
                      SSB0 = report$Bo_DD, h = report$h,
                      U = structure(report$U_DD, names = Year),
                      U_UMSY = structure(report$relU_DD, names = Year),
                      B = structure(report$B_DD, names = Yearplusone),
                      B_BMSY = structure(report$relB_DD, names = Yearplusone),
                      B_B0 = structure(report$B_DD/report$Bo_DD, names = Yearplusone),
                      SSB = structure(report$B_DD, names = Yearplusone),
                      SSB_SSBMSY = structure(report$relB_DD, names = Yearplusone),
                      SSB_SSB0 = structure(report$B_DD/report$Bo_DD, names = Yearplusone),
                      N = structure(report$N_DD, names = Yearplusone),
                      R = structure(report$R_DD, names = Yearplusk),
                      Catch = structure(report$Cpred_DD, names = Year),
                      NLL = report$jnll,
                      info = output$info, obj = output$obj, opt = output$opt,
                      SD = output$SD, TMB_report = report,
                      dependencies = output$dependencies,
                      Data = output$Data)
  }

  if(Model == "DD_SS") {
    Year <- output$info$Year
    Yearplusone <- c(Year, max(Year) + 1)
    k <- output$info$data$k_DD
    Yearplusk <- c(Year, (max(Year)+1):(max(Year)+k))
    Yearrandom <- seq(Year[1] + k, max(Year))
    Assessment <- new("Assessment", Model = Model,
                      MSY = report$MSY_DD, UMSY = report$UMSY_DD, BMSY = report$BMSY_DD,
                      B0 = report$Bo_DD, R0 = report$Ro_DD, N0 = report$No_DD,
                      SSB0 = report$Bo_DD, h = report$h,
                      U = structure(report$U_DD, names = Year),
                      U_UMSY = structure(report$relU_DD, names = Year),
                      B = structure(report$B_DD, names = Yearplusone),
                      B_BMSY = structure(report$relB_DD, names = Yearplusone),
                      B_B0 = structure(report$B_DD/report$Bo_DD, names = Yearplusone),
                      SSB = structure(report$B_DD, names = Yearplusone),
                      SSB_SSBMSY = structure(report$relB_DD, names = Yearplusone),
                      SSB_SSB0 = structure(report$B_DD/report$Bo_DD, names = Yearplusone),
                      N = structure(report$N_DD, names = Yearplusone),
                      R = structure(report$R_DD, names = Yearplusk),
                      Catch = structure(report$Cpred_DD, names = Year),
                      Random = structure(output$SD$par.random, names = Yearrandom),
                      Random_SE = structure(sqrt(output$SD$diag.cov.random), names = Yearrandom),
                      NLL = report$jnll, NLL_Catch = report$jnll_comp[1],
                      NLL_Random = report$jnll_comp[2],
                      info = output$info, obj = output$obj, opt = output$opt,
                      SD = output$SD, TMB_report = report,
                      dependencies = output$dependencies,
                      Data = output$Data)
  }

  return(Assessment)
}


# Call from inside generate_plots() and summary.Assessment
assign_Assessment_slots <- function() {
  Assessment <- get("Assessment", envir = parent.frame())
  Nslots <- length(slotNames(Assessment))
  for(i in 1:Nslots) {
    assign(slotNames(Assessment)[i], slot(Assessment, slotNames(Assessment)[i]),
           envir = parent.frame())
  }
  invisible()
}

# Call from inside generate_plots(), profile_likelihood(), retrospective(),
prepare_to_save_figure <- function() {
  Model <- get("Model", envir = parent.frame())
  base.dir <- get("save_dir", envir = parent.frame()) # by default: getwd()
  Model.dir <- paste0("plots_", Model)
  plot.dir <- file.path(base.dir, Model.dir)

  if(!dir.exists(plot.dir)) {
    message(paste0("Creating directory: \n", plot.dir, "\n"))
    dir.create(plot.dir)
  }
  assign("plot.dir", plot.dir, envir = parent.frame())
  invisible()
}


create_png <- function(filename, units = "in", res = 400, height = 4, width = 6) {
  png(filename = filename, units = units, res = res, height = height, width = width)
  par(mar = c(5, 4, 1, 1))
  invisible()
}


#' Plots a lognormal variable
#'
#' Plots the probability distribution function of a lognormal variable from the
#' mean and standard deviation in either transformed (normal) or untransformed space.
#'
#' @param m A vector of means of the distribution.
#' @param sd A vector of standard deviations of the distribution.
#' @param label Name of the variable to be used as x-axis label.
#' @param logtransform Indicates whether the mean and standard deviation are in
#' transformed (normal) or untransformed space.
#' @param color A vector of colors.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution. This function can plot multiple curves when multiple means
#' and standard deviations are provided.
#' @author Q. Huynh
#' @export plot_lognormalvar
#' @seealso \code{\link{plot_betavar}} \code{\link{plot_steepness}}
#' @examples
#' mu <- 0.5
#' stddev <- 0.1
#' plot_lognormalvar(mu, stddev) # mean of plot should be 0.5
#'
#' #logtransformed parameters
#' mu <- 0
#' stddev <- 0.1
#' plot_lognormalvar(mu, stddev, logtransform = TRUE) # mean of plot should be 1
plot_lognormalvar <- function(m, sd, label = NULL, logtransform = FALSE, color = "black") {
  # plots life history parameters: Linf, K, t0, M, FMSY_M
  ncurve <- length(m)
  if(!logtransform) {
    true.m <- m
    if(all(m) < 0) m <- -1 * m # special case needed when t0 < 0
    mu <- DLMtool:::mconv(m, sd)
    sdlog <- sdconv(m, sd)
    support <- seq(0.001, max(m + 5*sdlog), length.out = 1e3)

    dist <- matrix(NA, nrow = length(support), ncol = length)
    for(i in 1:ncurve) dist[, i] <- dlnorm(support, mu[i], sdlog[i])
    dist[is.infinite(dist)] <- NA

    dist.max <- max(dist, na.rm = TRUE)
    tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
    tails <- which(!tails)
    ind.tails <- c(tails[1]:tails[length(tails)])

    support <- support[ind.tails]
    dist <- as.matrix(dist[ind.tails, ])

    if(all(true.m) < 0) {
      support <- -1 * support
      xlim_truncated <- range(pretty(support))
      plot(support, dist[, 1], typ = 'l', xlab = label,
           ylab = 'Probability density function', xlim = xlim_truncated,
           ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), color = color[1])
      if(ncurve > 1) {
        for(i in 2:ncurve) lines(support, dist[, i], color = color[i])
      }
    }
    if(all(true.m) > 0) {
      xlim_truncated <- range(pretty(support))
      plot(support, dist[, 1], typ = 'l', xlab = label,
           ylab = 'Probability density function', xlim = xlim_truncated,
           ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), color = color[1])
      if(ncurve > 1) {
        for(i in 2:ncurve) lines(support, dist[, i], color = color[i])
      }
    }
    abline(h = 0, col = 'grey')
    abline(v = true.m, lty = 2, col = color.vec)
  }

  if(logtransform) {
    #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
    #where f is the pdf of distribution, g(y) = exp(X) is the transformation
    #y is the lognormal variable, x is a normal variable
    support.norm <- seq(min(m - 5*sd, na.rm = TRUE), max(m+5*sd, na.rm = TRUE),
                        length.out = 1e3)
    support <- exp(support.norm)

    dist <- matrix(NA, nrow = length(support), ncol = length(m))
    for(i in 1:ncurve) dist[, i] <- dnorm(support.norm, m[i], sd[i])/abs(support)
    dist[is.infinite(dist)] <- NA

    dist.max <- max(dist, na.rm = TRUE)
    tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
    tails <- which(!tails)
    ind.tails <- c(tails[1]:tails[length(tails)])

    support <- support[ind.tails]
    dist <- as.matrix(dist[ind.tails, ])

    xlim_truncated <- range(pretty(support), finite = TRUE, na.rm = TRUE)

    plot(support, dist[, 1], typ = 'l', xlab = label, xlim = xlim_truncated,
         ylab = 'Probability density function',
         ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
    if(ncurve > 1) {
      for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
    }
    abline(h = 0, col = 'grey')
    abline(v = exp(m), lty = 2, col = color)
  }

  invisible()
}

#' Plots a beta variable
#'
#' Plots the probability distribution function of a beta variable from the
#' mean and standard deviation in either transformed (logit) or untransformed space.
#'
#' @param m A vector of means of the distribution.
#' @param sd A vector of standard deviations of the distribution.
#' @param label Name of the variable to be used as x-axis label.
#' @param logit Boolean that indicates whether the means and standard deviations are in
#' transformed (logit) or untransformed space.
#' @param color A vector of colors.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution. This function can plot multiple curves when multiple means
#' and standard deviations are provided.
#' @author Q. Huynh
#' @export plot_betavar
#' @seealso \code{\link{plot_lognormalvar}} \code{\link{plot_steepness}}
#' @examples
#' mu <- 0.5
#' stddev <- 0.1
#' plot_betavar(mu, stddev) # mean of plot should be 0.5
#'
#' #logit parameters
#' mu <- 0
#' stddev <- 0.1
#' plot_betavar(mu, stddev, logit = TRUE) # mean of plot should be 0.5
plot_betavar <- function(m, sd, label = NULL, logit = FALSE, color = "black") {
  # currently plots life history: BMSY_B0, depletion
  support <- seq(0.01, 0.99, length.out = 1e3)
  ncurve <- length(m)
  dist <- matrix(NA, nrow = length(support), ncol = ncurve)
  if(!logit) {
    a <- alphaconv(m, sd)
    b <- betaconv(m, sd)
    for(i in 1:ncurve) dist[, i] <- dbeta(support, a[i], b[i])
  }
  if(logit) {
    for(i in 1:ncurve) {
      #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
      #where f is the pdf of distribution, g(y) = 1/(1 + exp(-X)) is the transformation
      #y is the beta variable, x is a normal variable
      dist[, i] <- dnorm(log(support/(1-support)), m[i], sd[i])/abs(support * (1-support))
      m[i] <- 1/(1 + exp(-m[i]))
    }
  }
  dist[is.infinite(dist)] <- NA

  dist.max <- max(dist, na.rm = TRUE)
  tails <- apply(dist, 1, function(x) all(x < 0.001 * dist.max))
  tails <- which(!tails)
  ind.tails <- c(tails[1]:tails[length(tails)])

  support <- support[ind.tails]
  dist <- as.matrix(dist[ind.tails, ])
  xlim_truncated <- range(pretty(support), finite = TRUE)

  plot(support, dist[, 1], typ = 'l', xlab = label, ylab = 'Probability density function',
       xlim = xlim_truncated, ylim = c(0, 1.1 * max(dist, na.rm = TRUE)), col = color[1])
  if(ncurve > 1) {
    for(i in 2:ncurve) lines(support, dist[, i], col = color[i])
  }
  abline(h = 0, col = 'grey')
  abline(v = m, lty = 2, col = color)

  invisible()
}

#' Plots probability distribution function of stock-recruit steepness
#'
#' Plots the probability distribution function of steepness from the
#' mean and standard deviation.
#'
#' @param m The mean of the distribution.
#' @param sd The standard deviation of the distribution.
#' @return A plot of the probability distribution function. Vertical dotted line
#' indicates mean of distribution.
#' @note The function samples from a beta distribution with parameters alpha and beta
#' that are converted from the mean and standard deviation. Then, the distribution is
#' transformed from 0 - 1 to 0.2 - 1.
#' @author Q. Huynh
#' @export plot_steepness
#' @seealso \code{\link{plot_lognormalvar}} \code{\link{plot_betavar}}
#' @examples
#' data(Red_snapper)
#' mu <- Red_snapper@steep
#' stddev <- Red_snapper@steep * Red_snapper@CV_steep
#' plot_steepness(mu, stddev)
#'
plot_steepness <- function(m, sd) {
  #f_y(y) = f_x(g-1(y)) * abs(d/dy[g-1(y)])
  #where f is the pdf of distribution, g(y) = 0.8X + 0.2 is the transformation
  #y is steepness, x is a beta variable
  m.transformed <- (m - 0.2)/0.8
  a <- alphaconv(m = m.transformed, sd = sd/0.8)
  b <- betaconv(m = m.transformed, sd = sd/0.8)

  support <- seq(0.201, 0.999, 0.001)
  dist <- dbeta((support - 0.2)/0.8, a, b) / 0.8

  plot(support, dist, typ = 'l', xlab = 'Steepness (h)',
       ylab = 'Probability density function', xlim = c(0.2, 1),
       ylim = c(0, 1.1 * max(dist, na.rm = TRUE)))
  abline(h = 0, col = 'grey')
  abline(v = m, lty = 2)

  invisible()
}





#' Plot time series of data
#'
#' Plot time series of observed (with lognormally-distributed error bars) vs.
#' predicted data.
#'
#' @param Year A vector of years for the data.
#' @param obs A vector of observed data.
#' @param fit A vector of predicted data (e.g., from an assessment model).
#' @param obs_CV A vector of year-specific coefficient of variation in the observed data.
#' @param obs_CV_CI The confidence interval for the error bars based for \code{obs_CV}.
#' @param obs_upper A vector of year-specific upper bounds for the error bars of the observed data (in lieu of argument \code{obs_CV}).
#' @param obs_lower A vector of year-specific lower bounds for the error bars of the observed data (in lieu of argument \code{obs_CV}).
#' @param fit_linewidth Argument \code{lwd} for fitted line.
#' @param fit_color Color of fitted line.
#' @param label Character string that describes the data to label the y-axis.
#' @author Q. Huynh
#' @seealso \code{\link{plot_residuals}}
#' @examples
#' data(Red_snapper)
#' plot_timeseries(Red_snapper@Year, Red_snapper@Cat[1, ], obs_CV = Red_snapper@CV_Cat, label = "Catch")
#' @export plot_timeseries
plot_timeseries <- function(Year, obs, fit = NULL, obs_CV = NULL, obs_CV_CI = 0.95,
                            obs_upper = NULL, obs_lower = NULL, fit_linewidth = 3,
                            fit_color = "red", label = "Observed data") {
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))

  # Without CV interval
  if(is.null(obs_CV)) {
    y.max <- max(c(obs, fit), na.rm = TRUE)
    plot(Year, obs, typ = 'o', ylab = label, ylim = c(0, 1.1 * y.max))
  }

  # With CV interval
  if(!is.null(obs_CV) || (!is.null(obs_upper) & !is.null(obs_lower))) {
    sigma <- sdconv(1, obs_CV)
    if(is.null(obs_upper))
      obs_upper <- exp(log(obs) + qnorm(1-0.5*(1-obs_CV_CI)) * sigma)
    if(is.null(obs_lower))
      obs_lower <- exp(log(obs) + qnorm(0.5*(1-obs_CV_CI)) * sigma)
    y.max <- max(c(obs_lower, obs_upper, obs, fit), na.rm = TRUE)

    plot(Year, obs, typ = 'o', ylab = label, ylim = c(0, 1.1 * y.max))
    arrows(Year, obs_lower, Year, obs_upper, length = 0.025, angle = 90,
           code = 3, col = 'grey30')
    lines(Year, fit, lwd = 2)
  }
  if(!is.null(fit)) lines(Year, fit, lwd = fit_linewidth, col = fit_color)
  abline(h = 0, col = 'grey')

  invisible()
}

#' Plot residuals
#'
#' Plots figure of residuals (or any time series with predicted mean of zero).
#'
#' @param Year A vector of years for the data.
#' @param res A vector of residuals.
#' @param res_sd A vector of year specific standard deviation for \code{res}.
#' @param res_sd_CI The confidence interval for the error bars based for \code{res_sd}.
#' @param res_upper A vector of year-specific upper bounds for the error bars of the residual (in lieu of argument \code{res_CV}).
#' @param res_lower A vector of year-specific lower bounds for the error bars of the residual (in lieu of argument \code{res_CV}).
#' @param draw_zero Indicates whether a horizontal line should be drawn at zero.
#' @param zero_linetype Passes argument \code{lty} (e.g. solid line = 1, dotted = 2) to \code{draw_zero}.
#' @param label Character string that describes the data to label the y-axis.
#' @author Q. Huynh
#' @seealso \code{\link{plot_timeseries}}
#' @export plot_residuals
plot_residuals <- function(Year, res, res_sd = NULL, res_sd_CI = 0.95,
                           res_upper = NULL, res_lower = NULL, draw_zero = TRUE,
                           zero_linetype = 2, label = "Residual") {
  # Without sd interval
  if(is.null(res_sd)) {
    res.lim <- max(abs(res), na.rm = TRUE)
    plot(Year, res, typ = 'o', ylab = label, ylim = c(-res.lim, res.lim))
  }

  # With CV interval
  if(!is.null(res_sd) || (!is.null(res_upper) & !is.null(res_lower))) {
    if(is.null(res_upper))
      res_upper <- res + qnorm(1-0.5*(1-res_sd_CI)) * res_sd
    if(is.null(res_lower))
      res_lower <- res + qnorm(0.5*(1-res_sd_CI)) * res_sd
    res.lim <- max(abs(c(res_lower, res_upper, res)), na.rm = TRUE)

    plot(Year, res, typ = 'o', ylab = label, ylim = c(-res.lim, res.lim))
    arrows(Year, res_lower, Year, res_upper, length = 0.025, angle = 90,
           code = 3, col = 'grey30')
  }
  if(draw_zero) abline(h = 0, lty = zero_linetype)
  invisible()
}


#plot_ts(Red_snapper@Year, Red_snapper@Cat[1, ], Red_snapper@CV_Cat, "Catch")
#plot_ts(Red_snapper@Year, Red_snapper@Ind[1, ], Red_snapper@CV_Ind, "Abundance Index")
#nyrs <- dim(Red_snapper@CAL)[2]
#nbins <- dim(Red_snapper@CAL)[3]
#CAL_bins <- 0.5 * (Red_snapper@CAL_bins[1:nbins] + Red_snapper@CAL_bins[2:(nbins+1)])
#par(mar = c(5.1, 4.1, 2.1, 2.1), oma = c(0, 0, 0, 0))
#plot_composition(1:nyrs, obs = Red_snapper@CAL[1, , ],
#          plot_type = 'annual', data_type = 'length',
#          CAL_bins = CAL_bins)


#plot_composition(1:nyrs, obs = Red_snapper@CAL[1, , ],
#                 plot_type = 'bubble', data_type = 'length',
#                 CAL_bins = CAL_bins)

#nyrs <- dim(Red_snapper@CAA)[2]
#par(mar = c(5.1, 4.1, 2.1, 2.1), oma = c(0, 0, 0, 0))
#plot_composition(1:nyrs, obs = Red_snapper@CAA[1, , ],
#                 plot_type = 'mean', data_type = 'age')


#' Plot composition data
#'
#' Plots annual length or age composition data.
#'
#' @param Year A vector of years.
#' @param obs A matrix of either length or age composition data. For lengths, rows and columns
#' should index years and length bin, respectively. For ages, rows and columns should index
#' years and age, respectively.
#' @param fit A matrix of predicted length or age composition from an assessment model.
#' Same dimensions as obs.
#' @param plot_type Indicates which plots to create. Options include annual distributions,
#' bubble plot, and annual means.
#' @param data_type Indicates whether length or age data are being used.
#' @param CAL_bins A vector of lengths corresponding to the columns in \code{obs}
#' and \code{fit}. Ignored for age data.
#' @return Plots depending on \code{plot_type}.
#' @author Q. Huynh
#' @export plot_composition
plot_composition <- function(Year, obs, fit = NULL,
                             plot_type = c('annual', 'bubble', 'mean'),
                             data_type = c(NULL, 'length', 'age'), CAL_bins = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)

  plot_type <- match.arg(plot_type, several.ok = TRUE)
  data_type <- match.arg(data_type)
  if(is.null(data_type)) stop('Indicate in data_type whether age or length data are being considered.')
  if(data_type == 'length' & is.null(CAL_bins)) {
    stop('Need vector of length bins.')
  }

  if(data_type == 'length') {
    data_val <- CAL_bins
    data_lab <- "Length"
  }
  if(data_type == 'age') {
    MaxAge <- ncol(obs)
    data_val <- 1:MaxAge
    data_lab <- "Age"
  }

  # Annual comps (obs vs. fitted if available)
  if(plot_type %in% 'annual') {

    par(mfcol = c(4, 4), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))
    ylim <- c(0, 1.1)
    yaxp <- c(0, 1, 4)
    las <- 1
    type <- 'o'
    for(i in 1:length(Year)) {
      obs.vec <- obs[i, ]
      N.true <- sum(obs.vec)
      obs.vec <- obs.vec/N.true
      obs.vec <- obs.vec/max(obs.vec, na.rm = TRUE)
      N <- round(N.true, 0)

      if(i < (length(Year)-1)) {
        if(i %% 16 %in% c(1:4)) { # First column
          yaxt <- 's'

          # First three rows
          if(i %% 4 %in% c(1:3)) {
            xaxt <- 'n'
          } else {
            xaxt <- 's'
          }
        } else { # All other columns
          if(i %% 4 %in% c(1:3)) { # First three rows
            xaxt <- yaxt <- 'n'
          } else {
            xaxt <- 's'
          }
        }
        plot(data_val, obs.vec, typ = 'o', ylim = ylim, yaxp = yaxp,
             xaxt = xaxt, yaxt = yaxt, las = las)
        abline(h = 0, col = 'grey')
        legend('topright', legend = c(Year[i], paste0("N = ", N)), bty = 'n', xjust = 1)
      }

      if(i == length(Year)) {
        xaxt <- 's'
        if(i %% 16 %in% c(1:4)) yaxt <- 's' else yaxt <- 'n'
        plot(data_val, obs.vec, typ = 'o', ylim = ylim, yaxp = yaxp,
             xaxt = xaxt, yaxt = yaxt, las = las)
        abline(h = 0, col = 'grey')
        legend('topright', legend = c(Year[i], paste0("N = ", N)), bty = 'n', xjust = 1)
      }
      if(i %% 16 == 0 || i == length(Year)) {
        mtext(data_lab, side = 1, line = 3, outer = TRUE)
        mtext('Relative Frequency', side = 2, line = 3.5, outer = TRUE)
      }
    }

  }
  # Bubble plot (obs)
  if(plot_type %in% 'bubble') {
    par(mfcol = c(1, 1), oma = rep(0, 4), mar = c(5.1, 5.1, 2.1, 2.1))

    #x.pretty <- pretty(Year)
    radius <- 5 / max(obs, na.rm = TRUE)
    plot(NULL, NULL, typ = 'n', xlim = range(Year), xlab = "Year",
         ylim = c(0, max(data_val)), ylab = data_lab)
    for(i in 1:length(Year)) {
      #obs.vec <- obs[i, ]/sum(obs[i, ])
      for(j in 1:length(data_val)) {
        points(Year[i], data_val[j], pch = 1, cex = radius * obs[i, j])
      }
    }
  }
  # Bubble plot (residuals if applicable)

  # Mean length or age over time
  if(plot_type %in% 'mean') {
    mu <- numeric(length = length(Year))
    for(i in 1:length(mu)) mu[i] <- weighted.mean(data_val, obs[i, ], na.rm = TRUE)

    x.pretty <- pretty(Year)
    y.pretty <- pretty(mu)
    plot(Year, mu, xlim = range(x.pretty), ylim = range(y.pretty),
         ylab = paste0('Mean ', data_type), typ = 'o')
  }

  invisible()
}




#' Plot life history parameters from Data object
#'
#' Produces plots of growth parameters.
#'
#' @param Data An object of class Data.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @param Model Name of assessment model to save into appropriate sub-directory (optional).
#' @return Plots of length-at-age, weight-at-age, and weight-at-length
#' (if appropriate parameters are available).
#'
#' @author Q. Huynh
#' @seealso \code{\link{plot_steepness}}
#' \code{\link{plot_lognormalvar}} \code{\link{plot_betavar}}
#' @export plot_life_history
#' @examples plot_life_history(Red_snapper)
plot_life_history <- function(Data, save_figure = FALSE, save_dir = getwd(),
                              Model = NULL) {
  if(save_figure) prepare_to_save_figure()

  # Plot growth
  # Need: Linf, K, t0, Maxage
  # Optional: LenCV, a, b
  if(!is.na(Data@vbLinf) & !is.na(Data@vbK) & !is.na(Data@vbt0) & !is.na(Data@MaxAge)) {

    Linf <- Data@vbLinf
    K <- Data@vbK
    t0 <- Data@vbt0
    MaxAge <- Data@MaxAge
    if(!is.na(Data@LenCV)) LenCV <- Data@LenCV else LenCV <- NULL

    age <- 1:MaxAge
    Laa <- Linf * (1 - exp(-K * (age - t0)))
    if(!is.null(LenCV)) {
      sigmaLaa <- LenCV * Laa
      upper <- Laa + 2 * sigmaLaa
      lower <- Laa - 2 * sigmaLaa
      y.max <- max(upper)
    } else y.max <- max(Laa)
    plot(age, Laa, typ = 'n', ylim = c(0, 1.1 * y.max), xlab = 'Age', ylab = 'Length')
    if(!is.null(LenCV)) {
      polygon(c(age, rev(age)), c(upper, rev(lower)), col = 'grey80', border = NA)
      title.add <- " (with 95% probability interval)"
    } else title.add <- ""
    lines(age, Laa, typ = 'o')
    abline(h = 0, col = 'grey')

    if(save_figure) {
      create_png(filename = file.path(plot.dir, "lifehistory_1_length_at_age.png"))
      plot(age, Laa, typ = 'n', ylim = c(0, 1.1 * y.max), xlab = 'Age', ylab = 'Length')
      if(!is.null(LenCV)) {
        polygon(c(age, rev(age)), c(upper, rev(lower)), col = 'grey80', border = NA)
        title.add <- "(with 95% probability interval)"
      } else title.add <- ""
      lines(age, Laa, typ = 'o')
      abline(h = 0, col = 'grey')
      dev.off()
      lh.file.caption <- c("lifehistory_1_length_at_age.png",
                           paste("Length at age", title.add, "from Data object."))
    }

    if(!is.na(Data@wla) & !is.na(Data@wlb)) {
      a <- Data@wla
      b <- Data@wlb

      Waa <- a * Laa ^ b
      plot(age, Waa, typ = 'o', ylim = c(0, 1.1 * max(Waa)),
           xlab = 'Age', ylab = 'Weight')
      abline(h = 0, col = 'grey')
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "lifehistory_2_mean_weight_at_age.png"))
        plot(age, Waa, typ = 'o', ylim = c(0, 1.1 * max(Waa)),
             xlab = 'Age', ylab = 'Weight')
        abline(h = 0, col = 'grey')
        dev.off()
        lh.file.caption <- rbind(lh.file.caption,
                                 c("lifehistory_2_mean_weight_at_age.png",
                                   "Mean weight at age from Data object."))
      }

      plot(Laa, Waa, typ = 'o', xlab = 'Length', ylab = 'Weight')
      abline(h = 0, col = 'grey')
      if(save_figure) {
        create_png(filename = file.path(plot.dir, "lifehistory_3_length_weight.png"))
        plot(Laa, Waa, typ = 'o', xlab = 'Length', ylab = 'Weight')
        abline(h = 0, col = 'grey')
        dev.off()
        lh.file.caption <- rbind(lh.file.caption,
                                 c("lifehistory_2_mean_weight_at_age.png",
                                   "Length-weight relationship from Data object."))
      }
    }
  }

  if(save_figure) return(invisible(lh.file.caption))
  else return(invisible())
}



#' @importFrom graphics arrows
plot_surplus_production <- function(B, B0 = NULL, C, arrow_size = 0.07) {
  old.warning <- options()$warn
  on.exit(options(warn = old.warning))
  options(warn = -1)

  if(!is.null(B0)) {
    B <- B/B0
    xlab_label <- expression(B/B[0])
  } else {
    xlab_label <- "Biomass"
  }
  B_now <- B[1:(length(B)-1)]
  B_next <- B[2:length(B)]
  SP_now <- B_next - B_now + C

  xlim <- c(0, max(B))
  ylim <- c(min(0, min(SP_now)), max(SP_now))
  plot(B_now, SP_now, typ = 'n', xlab = xlab_label, xlim = xlim, ylim = ylim,
       ylab = "Surplus production")
  arrows(x0 = B_now, y0 = SP_now[1:(length(B)-1)], x1 = B_next, y1 = SP_now[2:length(B)],
         length = arrow_size)
  abline(h = 0, col = 'grey')

  invisible()
}


plot_Kobe <- function(biomass, exploit, arrow_size = 0.07, color = TRUE) {
  old.warning <- options()$warn
  on.exit(options(warn = old.warning))
  options(warn = -1)
  n.arrows <- length(exploit)
  if(length(biomass) > n.arrows) biomass <- biomass[1:n.arrows]

  x.max <- max(biomass)
  y.max <- max(exploit)
  plot(NULL, NULL, typ = 'n', xlab = expression(B/B[MSY]), ylab = expression(U/U[MSY]),
       xlim = c(0, 1.1 * x.max), ylim = c(0, 1.1 * y.max))
  if(color) {
    # Colors from https://www.rapidtables.com/web/color/html-color-codes.html
    green <- "#228B22"    #forestgreen
    yellow <- "#F0E68C"   #khaki
    red <- "#CD5C5C"      #indianred
    polygon(x = c(1, 1, 10*x.max, 10*x.max), y = c(1, -5, -5, 1), col = green, border = NA)
    polygon(x = c(-5, -5, 1, 1), y = c(1, -5, -5, 1), col = yellow, border = NA)
    polygon(x = c(1, 1, 10*x.max, 10*x.max), y = c(10*y.max, 1, 1, 10*y.max),
            col = yellow, border = NA)
    polygon(x = c(-5, -5, 1, 1), y = c(10*y.max, 1, 1, 10*y.max),
            col = red, border = NA)
    box()
  }
  arrows(x0 = biomass[1:(n.arrows-1)], y0 = exploit[1:(n.arrows-1)],
         x1 = biomass[2:n.arrows], y1 = exploit[2:n.arrows], length = arrow_size)
  abline(h = 0, col = 'grey')
  abline(v = 0, col = 'grey')
  if(!color) {
    abline(h = 1, lty = 2)
    abline(v = 1, lty = 2)
  }
  invisible()
}

plot_TAC <- function(TAC, median_line = TRUE) {
  TAC.dens <- density(TAC, na.rm = TRUE)
  plot(TAC.dens, xlab = "TAC", main = "")
  if(median_line) abline(v = median(TAC, na.rm = TRUE), lty = 2)
  abline(h = 0, col = "grey")
  invisible()
}


#' Plot stock-recruitment function
#'
#' Plot stock-recruitment (with recruitment deviations if estimated).
#'
#' @param Spawners A vector of the number of the spawners (x-axis).
#' @param expectedR A vector of the expected recruitment (from the
#' stock-recruit function) corresponding to values of \code{Spawners}.
#' @param R0 Virgin recruitment.
#' @param S0 Virgin spawners.
#' @param rec_dev If recruitment deviations are estimated, a vector of estimated recruitment
#' (in normal space) corresponding to values of \code{Spawners}.
#' @param trajectory Indicates whether arrows will be drawn showing the trajectory of
#' spawners and recruitment deviations over time.
#' @param y_zoom If recruitment deviations are plotted, the y-axis limit relative to
#' maximum expected recruitment \code{expectedR}. If \code{NULL}, all recruitments are plotted.
#' @author Q. Huynh
#' @return A stock-recruit plot
#' @export plot_SR
plot_SR <- function(Spawners, expectedR, R0, S0, rec_dev = NULL, trajectory = FALSE,
                    y_zoom = NULL) {
  if(is.null(rec_dev)) R.max <- 1.1 * max(expectedR)
  else {
    if(is.null(y_zoom)) R.max <- 1.1 * max(rec_dev)
    else R.max <- y_zoom * max(expectedR)
  }
  S.max <- 1.1 * max(Spawners)
  plot(Spawners, expectedR, typ = 'l', xlim = c(0, 1.05 * S.max), ylim = c(0, 1.1 * R.max),
       xlab = 'Spawners', ylab = 'Recruits')
  if(!trajectory) {
    if(is.null(rec_dev)) points(Spawners, expectedR)
    if(!is.null(rec_dev)) points(Spawners, rec_dev)
  }
  if(trajectory) {
    old.warning <- options()$warn
    on.exit(options(warn = old.warning), add = TRUE)
    options(warn = -1)

    n.arrows <- length(Spawners)

    arrows(x0 = Spawners[1:(n.arrows-1)], y0 = rec_dev[1:(n.arrows-1)],
           x1 = Spawners[2:n.arrows], y1 = rec_dev[2:n.arrows], length = 0.07)
  }
  points(S0, R0, col = 'red', pch = 16)
  abline(h = 0, col = 'grey')
  abline(v = 0, col = 'grey')
}

plot_ogive <- function(Age, ogive, label = "Selectivity") {
  plot(Age, ogive, ylab = label, typ = 'n', ylim = c(0, 1))
  abline(h = 0, col = 'grey')
  lines(Age, ogive, typ = 'o')

  invisible()
}
