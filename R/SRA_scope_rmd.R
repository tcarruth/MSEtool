
#' @rdname SRA_scope
#' @importFrom rmarkdown render
#' @export
report_SRA_scope <- function(OM, report_list, sims = 1:OM@nsim, filename = "SRA_scope", dir = tempdir(), Year = NULL,
                             open_file = TRUE, quiet = TRUE, ...) {

  # Subset cpars
  OM <- Sub_cpars(OM, sims)
  report_list[[2]] <- report_list[[2]][sims]

  # Generate markdown report
  filename_html <- paste0(filename, ".html")
  filename_rmd <- paste0(filename, ".Rmd")
  if(!dir.exists(dir)) {
    message("Creating directory: \n", dir)
    dir.create(dir)
  }

  message("Writing markdown file: ", file.path(dir, filename_rmd))

  ####### Assign variables
  nsim <- OM@nsim

  SD <- report_list[[1]]$SD
  report <- report_list[[1]]$report
  data <- report_list[[1]]$obj$env$data
  data$C_hist <- report$C_hist
  data$E_hist <- report$E_hist

  conv <- report$conv
  length_bin <- report$length_bin

  max_age <- OM@maxage
  age <- 1:max_age
  nyears <- OM@nyears
  if(is.null(Year)) Year <- 1:nyears

  nfleet <- data$nfleet
  nsurvey <- data$nsurvey

  ####### Document header
  header <- rmd_head(paste("Operating model conditioning for", ifelse(nchar(OM@Name) > 0, OM@Name, substitute(OM))), FALSE)

  ####### Updated historical OM parameters
  Year_matrix <- matrix(Year, ncol = nsim, nrow = nyears)
  Yearplusone_matrix <- matrix(c(Year, max(Year) + 1), ncol = nsim, nrow = nyears+1)

  OM_update <- c("# Summary {.tabset}\n",
                 "## Updated historical operating model parameters\n", rmd_SRA_R0(),
                 rmd_SRA_initD(), rmd_SRA_D(), rmd_SRA_Perr(), rmd_SRA_Find(), rmd_SRA_sel())

  ####### Output from all simulations {.tabset}
  fleet_output <- lapply(1:nfleet, rmd_SRA_fleet_output)

  if(any(data$I_hist > 0)) {
    survey_output <- lapply(1:nsurvey, rmd_SRA_survey_output)
  } else survey_output <- NULL

  all_sims_output <- c(fleet_output, survey_output, "### Model output\n",
                       rmd_SRA_R_output(), rmd_SRA_SSB_output())

  ####### Fit to mean inputs from operating model
  # Generate summary table (parameter estimates)
  sumry <- c("## Fit to mean inputs from operating model {.tabset}\n",
             "### Model Estimates\n",
             "`r as.data.frame(summary(SD))`\n\n")

  # Life History section
  LH_varies_fn <- function(x) {
    n_unique <- apply(x, 2, function(y) length(unique(y)))
    any(n_unique > 1)
  }
  LAA <- rmd_at_age(age, data$len_age[nyears, ], header = "### Life History\n", fig.cap = "Length-at-age in last historical year.",
                    label = "Mean Length-at-age")
  if(LH_varies_fn(data$len_age)) {
    LAA_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data$len_age", xlab = "Year", ylab = "Age",
                                zlab = "Length-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual length-at-age.")
  } else LAA_persp <- NULL

  LW <- rmd_LW(data$length_bin, data$wt_at_len)

  mat <- rmd_mat(age, data$mat[nyears, ], fig.cap = "Maturity-at-age in last historical year.")
  if(LH_varies_fn(data$mat)) {
    mat_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data$mat", xlab = "Year", ylab = "Age",
                                zlab = "Maturity-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual maturity-at-age.")
  } else mat_persp <- NULL

  NatM <- rmd_at_age(age, data$M[nyears, ], fig.cap = "Natural mortality in last historical year.", label = "Natural mortality")
  if(LH_varies_fn(data$M)) {
    NatM_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data$M", xlab = "Year", ylab = "Age",
                                 zlab = "Natural mortality", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual M-at-age.")
  } else NatM_persp <- NULL

  LH_section <- c(LAA, LAA_persp, LW, mat, mat_persp, NatM, NatM_persp)

  # Data and fit section
  individual_matrix_fn <- function(i, obs, pred, fig.cap, label) {
    rmd_assess_fit2("Year", paste0(obs, "[, ", i, "]"), paste0(pred, "[, ", i, "]"),
                    fig.cap = paste(fig.cap, i), label = paste(label, i))
  }
  individual_array_fn <- function(i, obs, pred, comps = c("age", "length")) {
    comps <- match.arg(comps)
    obs2 <- paste(obs, "[, , ", i, "]")
    pred2 <- paste(pred, "[, , ", i, "]")
    fig.cap2 <- paste0("Observed (black) and predicted (red) ", comps, " composition from fleet ", i, ".")
    if(comps == "age") {
      rmd_fit_comps("Year", obs2, pred2, type = "annual", fig.cap = fig.cap2)
    } else rmd_fit_comps("Year", obs2, pred2, type = "annual", CAL_bins = "data$length_bin", fig.cap = fig.cap2)
  }

  if(any(data$C_hist > 0)) {
    C_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data$C_hist", col = "rich.colors(nfleet)",
                             xlab = "Year", ylab = "Catch", fig.cap = "Catch time series.", header = "### Data and Fit\n")

    if(data$condition == "effort" || ncol(data$C_hist) > 1) {
      C_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "data$C_hist", pred = "report$Cpred",
                        fig.cap = "catch from fleet", label = "Fleet")
    } else C_plots <- NULL
  } else C_matplot <- C_plots <- NULL

  if(!all(is.na(data$E_hist)) && any(data$E_hist > 0)) {
    E_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data$E_hist", col = "rich.colors(nfleet)",
                             xlab = "Year", ylab = "Effort", fig.cap = "Effort time series.")
  } else E_matplot <- NULL

  if(!all(is.na(report$Ipred))) {
    I_plots <- lapply(1:nsurvey, individual_matrix_fn, obs = "data$I_hist", pred = "report$Ipred",
                      fig.cap = "index from survey", label = "Survey")
  } else I_plots <- NULL

  if(any(data$CAA_hist > 0)) {
    CAA_plots <- lapply(1:nfleet, individual_array_fn, obs = "data$CAA_hist", pred = "report$CAApred", comps = "age")
  } else CAA_plots <- NULL

  if(any(data$CAL_hist > 0)) {
    CAL_plots <- lapply(1:nfleet, individual_array_fn, obs = "data$CAL_hist", pred = "report$CALpred", comps = "length")
  } else CAL_plots <- NULL

  if(any(data$mlen > 0)) {
    ML_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "data$mlen", pred = "report$Ipred",
                       fig.cap = "Mean Length from Fleet", label = "mean lengths from fleet")
  } else ML_plots <- NULL

  data_section <- c(C_matplot, E_matplot, C_plots, I_plots, CAA_plots, CAL_plots, ML_plots)

  # Model output
  sel_matplot <- rmd_matplot(x = "matrix(data$length_bin, nrow(report$vul), nfleet)", y = "report$vul", col = "rich.colors(nfleet)",
                             xlab = "Length", ylab = "Selectivity",
                             fig.cap = "Selectivity by fleet.", header = "### Output \n")

  F_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "report$F", col = "rich.colors(nfleet)",
                           xlab = "Year", ylab = "Fishing Mortality (F)",
                           fig.cap = "Time series of fishing mortality by fleet.")

  SSB <- structure(report$E, names = c(Year, max(Year) + 1))
  SSB0 <- structure(report$E0, names = Year)
  if(length(unique(SSB0)) > 1) {
    SSB_plot <- rmd_assess_timeseries("SSB0", "unfished spawning depletion (growth and/or M are time-varying)",
                                      "expression(SSB[0])")
  } else SSB_plot <- NULL

  SSB_SSB0 <- structure(report$E/c(report$E0, report$E0[length(report$E0)]), names = c(Year, max(Year) + 1))

  R <- structure(report$R, names = c(Year, max(Year) + 1))
  N_at_age <- report$N
  N <- structure(rowSums(N_at_age), names = c(Year, max(Year) + 1))

  log_rec_dev <- structure(report$log_rec_dev, names = Year)
  log_rec_dev_SE <- summary(SD)[, 2]
  log_rec_dev_SE <- c(log_rec_dev_SE[names(log_rec_dev_SE) == "log_rec_dev"], 0)

  N_bubble <- rmd_bubble("c(Year, max(Year)+1)", "N_at_age", fig.cap = "Predicted abundance-at-age.")

  CAA_all <- apply(report$CAApred, c(1, 2), sum)
  CAA_bubble <- rmd_bubble("Year", "CAA_all", fig.cap = "Predicted catch-at-age (summed over all fleets).")

  CAL_all <- apply(report$CALpred, c(1, 2), sum)
  CAL_bubble <- rmd_bubble("Year", "CAL_all", CAL_bins = "data$length_bin",
                           fig.cap = "Predicted catch-at-length (summed over all fleets).")

  ts_output <- c(sel_matplot, F_matplot, rmd_SSB(), SSB_plot, rmd_SSB_SSB0(FALSE), rmd_R(),
                 rmd_residual("log_rec_dev", fig.cap = "Time series of recruitment deviations.", label = "log-Recruitment deviations"),
                 rmd_residual("log_rec_dev", "log_rec_dev_SE", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                              label = "log-Recruitment deviations", conv_check = TRUE),
                 rmd_N(), N_bubble, CAA_bubble, CAL_bubble)


  rmd <- c(header, OM_update, all_sims_output, sumry, LH_section, data_section, ts_output, rmd_footer())
  if(is.list(rmd)) rmd <- do.call(c, rmd)

  write(rmd, file = file.path(dir, filename_rmd))

  # Rendering markdown file
  message("Rendering markdown file to HTML: ", file.path(dir, filename_html))

  output <- rmarkdown::render(file.path(dir, filename_rmd), "html_document", filename_html, dir,
                              output_options = list(df_print = "paged"), quiet = quiet, ...)
  message("Rendering complete.")

  if(open_file) browseURL(file.path(dir, filename_html))
  invisible(output)
}


rmd_persp_plot <- function(x, y, z, xlab, ylab, zlab, phi, theta, expand, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           paste0("persp(x = ", x, ", y = ", y, ", z = ", z, ", theta = ", theta, ", phi = ", phi, ", expand = ", expand, ", xlab = \"", xlab, "\",
                   ylab = \"", ylab, "\", zlab = \"", zlab, "\", ticktype = \"detailed\")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_matplot <- function(x, y, col, xlab, ylab, legend.lab = "Fleet", type = "l", lty = 1, fig.cap, header = NULL) {
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           paste0("x <- ", x, "; y <- ", y),
           paste0("matplot(x, y, type = \"", type, "\", lty = ", lty, ", col = ", col,
                  ", ylim = c(0, 1.1 * max(y, na.rm = TRUE)), xlab = \"", xlab, "\", ylab = \"", ylab, "\")"),
           "abline(h = 0, col = \"grey\")",
           paste0("if(ncol(x) > 1) legend(\"topleft\", paste(\"", legend.lab, "\", 1:ncol(x)), text.col = ", col, ")"),
           " ```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

# For SRA scope function
rmd_assess_fit2 <- function(year, obs, fit, fig.cap, label = fig.cap, match = FALSE) {
  fig.cap2 <- paste0("Observed (black) and predicted (red) ", fig.cap, ".")
  if(match) fig.cap2 <- paste(fig.cap2, "Predicted", fig.cap, "should match observed in this model.")

  c(paste0("```{r, fig.cap=\"", fig.cap2, "\"}"),
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
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    paste0("ind_valid <- rowSums(", obs, ", na.rm = TRUE) > 0"),
    paste0("if(any(ind_valid)) plot_composition(", year, "[ind_valid], ", obs, "[ind_valid, ], ", fit, "[ind_valid, ], plot_type = ", arg, ")"),
    "```\n")
}

rmd_SRA_R0 <- function(fig.cap = "Histogram of R0 (unfished recruitment).") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$R0)) hist(OM@cpars$R0, main = \"\", xlab = expression(R[0]))",
    "```\n")
}


rmd_SRA_initD <- function(fig.cap = "Histogram of initial depletion.") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$initD)) hist(OM@cpars$initD, main = \"\", xlab = \"Initial depletion\")",
    "```\n")
}

rmd_SRA_D <- function(fig.cap = "Histogram of historical depletion.") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(!is.null(OM@cpars$D)) hist(OM@cpars$D, main = \"\", xlab = \"Depletion\")",
    "```\n")
}

rmd_SRA_Perr <- function(fig.cap = "Recruitment deviations among simulations.") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "Perr <- OM@cpars$Perr_y[, max_age:(max_age+nyears-1), drop = FALSE]",
    "matplot(Year_matrix, t(Perr), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment deviations\",",
    "        ylim = c(0, 1.1 * max(Perr)))",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_Find <- function(fig.cap = "Apical F from model used as effort among simulations.") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "matplot(Year_matrix, t(OM@cpars$Find), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Apical F\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_sel <- function(fig.cap = "Operating model selectivity among simulations.") {
  c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
    "if(nfleet == 1) {",
    "  vul <- do.call(cbind, lapply(report_list[[2]], getElement, \"vul\"))",
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
  if(ff == 1) header <- "## Output from all simulations {.tabset}\n" else header <- NULL
  ans <- c(paste("### Fleet", ff, "\n"),
           paste0("```{r, fig.cap = \"Selectivity of fleet ", ff, ".\"}"),
           paste0("vul_ff <- do.call(cbind, lapply(report_list[[2]], function(x) x$vul[, ", ff, "]))"),
           "matplot(matrix(length_bin, ncol = nsim, nrow = length(length_bin)), vul_ff, type = \"l\", col = \"black\",",
           paste0("        xlab = \"Length\", ylab = \"Selectivity of Fleet ", ff, "\")"),
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Fishing Mortality of fleet ", ff, ".\"}"),
           paste0("FM <- do.call(cbind, lapply(report_list[[2]], function(x) x$F[, ", ff, "]))"),
           paste0("matplot(Year_matrix, FM, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Fishing Mortality of Fleet ", ff, "\")"),
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) catch of fleet ", ff, ".\"}"),
           paste0("if(any(data$C_hist[, ", ff, "] > 0)) {"),
           paste0("  Cpred <- do.call(cbind, lapply(report_list[[2]], function(x) x$Cpred[, ", ff, "]))"),
           paste0("  matplot(Year_matrix, Cpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Catch of Fleet ", ff, "\")"),
           paste0("  lines(Year, data$C_hist[, ", ff, "], col = \"red\", lwd = 3)"),
           "} else {",
           paste0("  Cpred <- do.call(cbind, lapply(report_list[[2]], function(x) x$Cpred[, ", ff, "]/mean(x$Cpred[, ", ff, "])))"),
           paste0("  matplot(Year_matrix, Cpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Relative catch of Fleet ", ff, "\")"),
           "}",
           "abline(h = 0, col = \"grey\")",
           "```\n",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean ages of fleet ", ff, ".\"}"),
           paste0("MApred <- do.call(cbind, lapply(report_list[[2]], function(x) x$CAApred[, , ", ff, "] %*% age/x$CN[, ", ff, "]))"),
           "matplot(Year_matrix, MApred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean age\")",
           paste0("if(any(data$CAA_hist[, , ", ff, "] > 0)) {"),
           paste0("  lines(Year, (data$CAA_hist[, , ", ff, "] %*% age)/rowSums(data$CAA_hist[, , ", ff, "], na.rm = TRUE),",
           "  col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "}",
           "```\n",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean lengths of fleet ", ff, ".\"}"),
           paste0("MLpred <- do.call(cbind, lapply(report_list[[2]], function(x) x$mlen_pred[, ", ff, "]))"),
           "matplot(Year_matrix, MLpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean length\")",
           paste0("if(any(data$CAL_hist[, , ", ff, "] > 0)) {"),
           paste0("  lines(Year, (data$CAL_hist[, , ", ff, "] %*% length_bin)/rowSums(data$CAL_hist[, , ", ff, "], na.rm = TRUE),",
                  "  col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           paste0("} else if(any(data$mlen[, ", ff, "] > 0)) lines(Year, data$mlen[, ", ff, "], col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "```\n")
  c(header, ans)
}

rmd_SRA_survey_output <- function(sur) {
  if(sur == 1) header <- "### Surveys\n" else header <- NULL
  ans <- c(paste0("```{r, fig.cap = \"Observed (red) and predicted (black) index values in survey ", sur, ".\"}"),
           paste0("Ipred <- do.call(cbind, lapply(report_list[[2]], function(x) x$Ipred[, ", sur, "]))"),
           paste0("matplot(Year_matrix, Ipred, type = \"l\", col = \"black\", ylim = c(0, 1.1 * max(c(Ipred, data$I_hist[, ", sur, "]), na.rm = TRUE)), xlab = \"Year\", ylab = \"Survey ", sur, "\")"),
           paste0("lines(Year, data$I_hist[, ", sur, "], col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "abline(h = 0, col = \"grey\")",
           "```\n")
  c(header, ans)
}

rmd_SRA_R_output <- function() {
  c("```{r, fig.cap = \"Estimated recruitment among all simulations.\"}",
    "R_out <- do.call(cbind, lapply(report_list[[2]], getElement, \"R\"))",
    "matplot(Yearplusone_matrix, R_out, ylim = c(0, 1.1 * max(R_out)), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Recruitment\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_SRA_SSB_output <- function() {
  c("```{r, fig.cap = \"Estimated spawning biomass among all simulations.\"}",
    "E <- do.call(cbind, lapply(report_list[[2]], getElement, \"E\"))",
    "matplot(Yearplusone_matrix, E, ylim = c(0, 1.1 * max(E)), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Spawning biomass\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}


#SRA_scope_rmd(OM_SYC_D_RecC[[1]], OM_SYC_D_RecC[[3]], dir = getwd())


#SRA_scope_rmd(pcod_no_comps_SRA[[1]], pcod_no_comps_SRA[[3]], dir = getwd(), Year = 1956:2018)

