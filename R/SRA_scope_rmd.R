

rmd_persp_plot <- function(x, y, z, xlab, ylab, zlab, phi, theta, expand, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
           paste0("persp(x = ", x, ", y = ", y, ", z = ", z, ", theta = ", theta, ", phi = ", phi, ", expand = ", expand, ", xlab = \"", xlab, "\",
                   ylab = \"", ylab, "\", zlab = \"", zlab, "\", ticktype = \"detailed\")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_matplot <- function(x, y, col, xlab, ylab, legend.lab = "Fleet", type = "l", lty = 1, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
           paste0("xx <- ", x, "; yy <- ", y),
           paste0("matplot(xx, yy, type = \"", type, "\", lty = ", lty, ", col = ", col,
                  ", ylim = c(0, 1.1 * max(yy, na.rm = TRUE)), xlab = \"", xlab, "\", ylab = \"", ylab, "\")"),
           "abline(h = 0, col = \"grey\")",
           paste0("if(ncol(xx) > 1) legend(\"topleft\", paste(\"", legend.lab, "\", 1:ncol(x)), text.col = ", col, ")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

# For SRA scope function
rmd_assess_fit2 <- function(year, obs, fit, fig.cap, label = fig.cap, match = FALSE) {
  fig.cap2 <- paste0("Observed (black) and predicted (red) ", fig.cap, ".")
  if(match) fig.cap2 <- paste(fig.cap2, "Predicted", fig.cap, "should match observed in this model.")

  c(paste0("```{r, fig.cap = \"", fig.cap2, "\"}"),
    paste0("plot_timeseries(", year, ", ", obs, ", ", fit, ", label = \"", label, "\")"),
    "```\n")
}

rmd_fit_comps <- function(year, obs, fit, type = c("bubble", "annual"), ages = "NULL", CAL_bins = "NULL", fig.cap) {
  type <- match.arg(type)
  if(type == "bubble") {
    arg <- paste0("\"bubble_data\", CAL_bins = ", CAL_bins, ", ages = ", ages)
  } else {
    arg <- paste0("\"annual\", CAL_bins = ", CAL_bins, ", ages = ", ages)
  }
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    paste0("ind_valid <- rowSums(", obs, ", na.rm = TRUE) > 0"),
    paste0("if(any(ind_valid)) plot_composition(", year, "[ind_valid], ", obs, "[ind_valid, ], ", fit, "[ind_valid, ], plot_type = ", arg, ")"),
    "```\n")
}

rmd_SRA_R0 <- function(fig.cap = "Histogram of R0 (unfished recruitment).") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$R0)) hist(OM@cpars$R0, main = \"\", xlab = expression(R[0]))",
    "```\n")
}




rmd_SRA_D <- function(fig.cap = "Histogram of historical depletion.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$D)) hist(OM@cpars$D, main = \"\", xlab = \"Depletion\")",
    "```\n")
}

rmd_SRA_Perr <- function(fig.cap = "Recruitment deviations among simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "Perr <- OM@cpars$Perr_y[, max_age:(max_age+nyears-1), drop = FALSE]",
    "matplot(Year_matrix, t(Perr), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment deviations\",",
    "        ylim = c(0, 1.1 * max(Perr)))",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_Find <- function(fig.cap = "Apical F from SRA model. These values may be subsequently re-scaled in the operating model in order to match the specified depletion") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "matplot(Year_matrix, t(OM@cpars$Find), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Apical F\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_sel <- function(fig.cap = "Operating model selectivity among simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "if(nfleet == 1) {",
    "  vul <- do.call(cbind, lapply(report_list, getElement, \"vul_len\"))",
    "  matplot(matrix(length_bin, ncol = nsim, nrow = length(length_bin)), vul, type = \"l\", col = \"black\",",
    "          xlab = \"Length\", ylab = \"Selectivity\", ylim = c(0, 1.1))",
    "} else {",
    "  matplot(matrix(age, ncol = nsim, nrow = max_age), t(OM@cpars$V[, , nyears]), type = \"l\", col = \"black\",",
    "          xlab = \"Age\", ylab = \"Selectivity (last historical year)\", ylim = c(0, 1.1))",
    "}",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_fleet_output <- function(ff) {
  if(ff == 1) header <- "## SRA output {.tabset}\n" else header <- NULL
  ans <- c(paste("### Fleet", ff, "\n"),
           paste0("```{r, fig.cap = \"Selectivity of fleet ", ff, ".\"}"),
           paste0("vul_ff <- do.call(cbind, lapply(report_list, function(x) x$vul_len[, ", ff, "]))"),
           "matplot(length_bin, vul_ff, type = \"l\", col = \"black\",",
           paste0("        xlab = \"Length\", ylab = \"Selectivity of Fleet ", ff, "\")"),
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Fishing Mortality of fleet ", ff, ".\"}"),
           paste0("FM <- do.call(cbind, lapply(report_list, function(x) x$F[, ", ff, "]))"),
           paste0("matplot(Year_matrix, FM, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Fishing Mortality of Fleet ", ff, "\")"),
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) catch from fleet ", ff, ".\"}"),
           paste0("if(any(data$Chist[, ", ff, "] > 0)) {"),
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]))"),
           paste0("  matplot(Year_matrix, Cpred, type = \"o\", pch = 1, col = \"black\", xlab = \"Year\", ylab = \"Catch of Fleet ", ff, "\")"),
           paste0("  lines(Year, data$Chist[, ", ff, "], col = \"red\", lwd = 3)"),
           "} else {",
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]/mean(x$Cpred[, ", ff, "])))"),
           paste0("  matplot(Year_matrix, Cpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Relative catch of Fleet ", ff, "\")"),
           "}",
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean ages from fleet ", ff, ".\"}"),
           paste0("MApred <- do.call(cbind, lapply(report_list, function(x) x$CAApred[, , ", ff, "] %*% age/x$CN[, ", ff, "]))"),
           "matplot(Year_matrix, MApred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean age\")",
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, (data$CAA[, , ", ff, "] %*% age)/rowSums(data$CAA[, , ", ff, "], na.rm = TRUE),",
           "  col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean lengths from fleet ", ff, ".\"}"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$mlen_pred[, ", ff, "]))"),
           "matplot(Year_matrix, MLpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean length\")",
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, (data$CAL[, , ", ff, "] %*% length_bin)/rowSums(data$CAL[, , ", ff, "], na.rm = TRUE),",
                  "  col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           paste0("} else if(any(data$ML[, ", ff, "] > 0, na.rm = TRUE)) lines(Year, data$ML[, ", ff, "], col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) age composition from fleet ", ff, ".\"}"),
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_SRA(Year, x@CAA[, , , ", ff, "], data$CAA[, , ", ff, "])"),
           "}",
           "```\n",
           paste0("```{r, fig.cap = \"Predicted age composition from fleet ", ff, ".\"}"),
           paste0("if(all(is.na(data$CAA[, , ", ff, "]))) {"),
           paste0("plot_composition_SRA(Year, x@CAA[, , , ", ff, "])"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) length composition from fleet ", ff, ".\"}"),
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("plot_composition_SRA(Year, x@CAL[, , , ", ff, "], data$CAL[, , ", ff, "], CAL_bins = data$length_bin)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Predicted length composition from fleet ", ff, ".\"}"),
           paste0("if(all(is.na(data$CAL[, , ", ff, "]))) {"),
           paste0("plot_composition_SRA(Year, x@CAL[, , , ", ff, "], data$CAL[, , ", ff, "], CAL_bins = data$length_bin)"),
           "}",
           "```\n")

  c(header, ans)
}

rmd_SRA_survey_output <- function(sur) {
  if(sur == 1) header <- "### Surveys\n" else header <- NULL
  ans <- c(paste0("```{r, fig.cap = \"Observed (red) and predicted (black) index values in survey ", sur, ".\"}"),
           paste0("Ipred <- do.call(cbind, lapply(report_list, function(x) x$Ipred[, ", sur, "]))"),
           paste0("matplot(Year_matrix, Ipred, type = \"l\", col = \"black\", ylim = c(0, 1.1 * max(c(Ipred, data$I_hist[, ", sur, "]), na.rm = TRUE)), xlab = \"Year\", ylab = \"Survey ", sur, "\")"),
           paste0("lines(Year, data$Index[, ", sur, "], col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "abline(h = 0, col = \"grey\")",
           "```\n")
  c(header, ans)
}

rmd_SRA_initD <- function(fig.cap = "Histogram of initial depletion among all simulations.") {
  c(paste0("```{r, fig.cap = \"", fig.cap, "\"}"),
    "initD <- vapply(report_list, function(x) x$E[1]/x$E0[1], numeric(1))",
    "hist(initD, main = \"\", xlab = \"Initial depletion\")",
    "```\n")
}

rmd_SRA_R_output <- function() {
  c("```{r, fig.cap = \"Estimated recruitment among all simulations.\"}",
    "R_out <- do.call(cbind, lapply(report_list, getElement, \"R\"))",
    "matplot(Yearplusone_matrix, R_out, ylim = c(0, 1.1 * max(R_out)), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_SSB_output <- function() {
  c("```{r, fig.cap = \"Estimated spawning biomass among all simulations.\"}",
    "E <- do.call(cbind, lapply(report_list, getElement, \"E\"))",
    "matplot(Yearplusone_matrix, E, ylim = c(0, 1.1 * max(E)), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Spawning biomass\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}


plot_composition_SRA <- function(Year, SRA, dat = NULL, CAL_bins = NULL, ages = NULL, annual_ylab = "Frequency",
                                 annual_yscale = c("proportions", "raw"), dat_linewidth = 3, dat_color = "red") {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfcol = c(4, 4), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))

  annual_yscale <- match.arg(annual_yscale)
  if(is.null(CAL_bins)) data_type <- "age" else data_type <- "length"

  #if(!is.null(dat) && !all(dim(dat) == dim(SRA))) stop("Dimensions of 'SRA' and 'dat' do not match.")

  if(data_type == 'length') {
    data_val <- CAL_bins
    data_lab <- "Length"
  }
  if(data_type == 'age') {
    data_val <- if(is.null(ages)) 1:dim(SRA)[3] else ages
    data_lab <- "Age"
  }

  # Annual comps (SRA vs. dat if available)
  # Dim of
  SRA_plot <- SRA
  dat_plot <- dat
  if(annual_yscale == "proportions") {
    for(i in 1:length(Year)) {
      SRA_plot[, i, ] <- SRA[, i, ]/max(SRA[, i, ])
      if(!is.null(dat)) dat_plot[i, ] <- dat[i, ]/max(dat[i, ])
    }
  }
  ylim <- c(0, 1.1 * max(SRA_plot, dat_plot, na.rm = TRUE))
  yaxp <- c(0, max(pretty(ylim, n = 4)), 4)
  if(max(SRA_plot, dat_plot, na.rm = TRUE) == 1) yaxp <- c(0, 1, 4)

  las <- 1

  for(i in 1:length(Year)) {

    if(i < length(Year)) {
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
      matplot(data_val, t(SRA_plot[, i, ]), typ = "l", col = "black", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = 'grey')
      if(!is.null(dat)) lines(data_val, dat_plot[i, ], lwd = dat_linewidth, col = dat_color)
      legend("topright", legend = Year[i], bty = "n", xjust = 1)
    }

    if(i == length(Year)) {
      xaxt <- 's'
      if(i %% 16 %in% c(1:4)) yaxt <- 's' else yaxt <- 'n'
      matplot(data_val, t(SRA_plot[, i, ]), typ = "l", col = "black", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = 'grey')
      if(!is.null(dat)) lines(data_val, dat_plot[i, ], lwd = dat_linewidth, col = dat_color)
      legend("topright", legend = Year[i], bty = "n", xjust = 1)
    }
    if(i %% 16 == 0 || i == length(Year)) {
      mtext(data_lab, side = 1, line = 3, outer = TRUE)
      mtext(annual_ylab, side = 2, line = 3.5, outer = TRUE)
    }
  }

  invisible()
}




