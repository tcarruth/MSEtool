
#' Class-\code{SRA}
#'
#' An S4 class for the output from \link{SRA_scope}.
#'
#' @name SRA-class
#' @docType class
#'
#' @slot OM An updated operating model, class \linkS4class{OM}.
#' @slot SSB A matrix of estimated spawning biomass with \code{OM@@nsim} rows and \code{OM@@nyears+1} columns.
#' @slot NAA An array for the predicted numbers at age with dimension \code{OM@@nsim}, \code{OM@@nyears+1}, and \code{OM@@maxage}.
#' @slot CAA An array for the predicted catch at age with dimension \code{OM@@nsim}, \code{OM@@nyears}, \code{OM@@maxage}, and nfleet.
#' @slot CAL An array for the predicted catch at length with dimension \code{OM@@nsim}, \code{OM@@nyears}, length bins, and nfleet.
#' @slot conv A logical vector of length \code{OM@@nsim} indicating convergence of the SRA scoping model in the i-th simulation.
#' @slot Misc A list of length \code{OM@@nsim} with more output from the fitted SRA scoping model.
#' @slot mean_fit A list of output from fit to mean values of life history parameters in the operating model.
#' @slot data A list of the data inputs for the SRA scoping model.
#' @slot config A data frame describing configuration of the SRA scoping model.
#' @seealso \link{plot.SRA} \link{SRA_scope}
#' @author Q. Huynh
#' @export SRA
#' @exportClass SRA
SRA <- setClass("SRA", slots = c(OM = "ANY", SSB = "matrix", NAA = "array",
                                 CAA = "array", CAL = "array", conv = "logical", Misc = "list", mean_fit = "list",
                                 data = "list", config = "data.frame"))


#' @name plot.SRA
#' @aliases plot,SRA,missing-method
#' @title Plot SRA scope output
#' @description Produces HTML file (via markdown) figures of parameter estimates and output from an \linkS4class{Assessment} object.
#' Plots histograms of operating model parameters that are updated by the SRA scoping function, as well as diagnostic plots
#' for the fits to the SRA for each simulation.
#'
#' @param x An object of class \linkS4class{SRA} (output from \link{SRA_scope}).
#' @param compare Logical, if TRUE, the function will run \code{runMSE} to compare the historical period of the operating model
#' and the SRA model output.
#' @param filename Character string for the name of the markdown and HTML files.
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param sims A logical vector of length \code{x@@OM@@nsim} or a numeric vector indicating which simulations to keep.
#' @param Year Optional, a vector of years for the historical period for plotting.
#' @param open_file Logical, whether the HTML document is opened after it is rendered.
#' @param quiet Logical, whether to silence the markdown rendering function.
#' @param ... Other arguments to pass to \link[rmarkdown]{render}.
#' @return Returns invisibly the output from \link[rmarkdown]{render}.
#' @importFrom rmarkdown render
#' @seealso \linkS4class{SRA} \link{SRA_scope}
#' @exportMethod plot
setMethod("plot", signature(x = "SRA", y = "missing"),
          function(x, compare = TRUE, filename = "SRA_scope", dir = tempdir(), sims = 1:x@OM@nsim, Year = NULL, open_file = TRUE, quiet = TRUE, ...) {
            OM <- Sub_cpars(x@OM, sims)
            mean_fit <- x@mean_fit
            report_list <- x@Misc[sims]

            ####### Assign variables
            nsim <- OM@nsim

            data <- x@data

            max_age <- OM@maxage
            age <- 1:max_age
            nyears <- OM@nyears
            if(is.null(Year)) Year <- (OM@CurrentYr - nyears + 1):OM@CurrentYr

            nfleet <- x@data$nfleet
            nsurvey <- x@data$nsurvey
            length_bin <- x@data$length_bin

            ####### Document header
            header <- c("---",
                        "title: \"Operating model (OM) conditioning for `r ifelse(nchar(OM@Name) > 0, OM@Name, substitute(OM))`\"",
                        "subtitle: Output from Stock Reduction Analysis scoping function (SRA_scope)",
                        "date: \"`r Sys.Date()`\"",
                        "---",
                        "<style type=\"text/css\">",
                        "h1 { /* Header 1 */",
                        "  font-size: 24px;",
                        "}",
                        "</style>",
                        "",
                        "```{r setup, include = FALSE, echo = FALSE}",
                        "  knitr::opts_chunk$set(collapse = TRUE, echo = FALSE, message = FALSE,",
                        "  fig.width = 6, fig.height = 4.5, out.width = \"650px\", comment = \"#>\")",
                        "```\n")

            ####### Updated historical OM parameters
            Year_matrix <- matrix(Year, ncol = nsim, nrow = nyears)
            Yearplusone_matrix <- matrix(c(Year, max(Year) + 1), ncol = nsim, nrow = nyears+1)

            OM_update <- c("# Summary {.tabset}\n",
                           "## Updated historical OM parameters\n", rmd_SRA_R0(),
                           rmd_SRA_D(), rmd_SRA_Perr(), rmd_SRA_Find(), rmd_SRA_sel())

            ####### Output from all simulations {.tabset}
            fleet_output <- lapply(1:nfleet, rmd_SRA_fleet_output)

            if(any(data$Index > 0, na.rm = TRUE)) {
              survey_output <- lapply(1:nsurvey, rmd_SRA_survey_output)
            } else survey_output <- NULL

            all_sims_output <- c(fleet_output, survey_output, "### Model predictions\n",
                                 rmd_SRA_initD(), rmd_SRA_R_output(), rmd_SRA_SSB_output(), rmd_log_rec_dev())

            ####### Fit to mean inputs from operating model
            # Generate summary table (parameter estimates)

            if(length(x@mean_fit) > 0) {
              SD <- x@mean_fit$SD
              report <- x@mean_fit$report
              length_bin <- x@mean_fit$report$length_bin
              data_mean_fit <- x@mean_fit$obj$env$data

              conv <- report$conv

              sumry <- c("## Fit to mean parameters of the OM {.tabset}\n",
                         "### SRA Model Estimates\n",
                         "`r as.data.frame(summary(SD))`\n\n")

              # Life History section
              LH_varies_fn <- function(x) {
                n_unique <- apply(x, 2, function(y) length(unique(y)))
                any(n_unique > 1)
              }
              LAA <- rmd_at_age(age, data_mean_fit$len_age[nyears, ], header = "### Life History\n", fig.cap = "Length-at-age in last historical year.",
                                label = "Mean Length-at-age")
              if(LH_varies_fn(data_mean_fit$len_age)) {
                LAA_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$len_age[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Length-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual length-at-age.")
              } else LAA_persp <- NULL

              #LW <- rmd_LW(data_mean_fit$length_bin, data_mean_fit$wt_at_len)

              mat <- rmd_mat(age, data_mean_fit$mat[nyears, ], fig.cap = "Maturity-at-age in last historical year.")
              if(LH_varies_fn(data_mean_fit$mat)) {
                mat_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$mat[1:nyears, ]", xlab = "Year", ylab = "Age",
                                            zlab = "Maturity-at-age", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual maturity-at-age.")
              } else mat_persp <- NULL

              NatM <- rmd_at_age(age, data_mean_fit$M[nyears, ], fig.cap = "Natural mortality in last historical year.", label = "Natural mortality")
              if(LH_varies_fn(data_mean_fit$M)) {
                NatM_persp <- rmd_persp_plot(x = "Year", y = "age", z = "data_mean_fit$M[1:nyears, ]", xlab = "Year", ylab = "Age",
                                             zlab = "Natural mortality", phi = 35, theta = 45, expand = 0.55, fig.cap = "Annual M-at-age.")
              } else NatM_persp <- NULL

              LH_section <- c(LAA, LAA_persp, mat, mat_persp, NatM, NatM_persp)

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
                } else rmd_fit_comps("Year", obs2, pred2, type = "annual", CAL_bins = "data_mean_fit$length_bin", fig.cap = fig.cap2)
              }

              if(any(data_mean_fit$C_hist > 0, na.rm = TRUE)) {
                C_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data_mean_fit$C_hist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Catch", fig.cap = "Catch time series.", header = "### Data and Fit\n")

                if(data_mean_fit$condition == "effort" || ncol(data_mean_fit$C_hist) > 1) {
                  C_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "data_mean_fit$C_hist", pred = "report$Cpred",
                                    fig.cap = "catch from fleet", label = "Fleet")
                } else C_plots <- NULL
              } else C_matplot <- C_plots <- NULL

              if(any(data_mean_fit$E_hist > 0, na.rm = TRUE)) {
                E_matplot <- rmd_matplot(x = "matrix(Year, nyears, nfleet)", y = "data_mean_fit$E_hist", col = "rich.colors(nfleet)",
                                         xlab = "Year", ylab = "Effort", fig.cap = "Effort time series.")
              } else E_matplot <- NULL

              if(any(report$Ipred > 0, na.rm = TRUE)) {
                I_plots <- lapply(1:nsurvey, individual_matrix_fn, obs = "data_mean_fit$I_hist", pred = "report$Ipred",
                                  fig.cap = "index from survey", label = "Survey")
              } else I_plots <- NULL

              if(any(data_mean_fit$CAA_hist > 0, na.rm = TRUE)) {
                CAA_plots <- lapply(1:nfleet, individual_array_fn, obs = "data_mean_fit$CAA_hist", pred = "report$CAApred", comps = "age")
              } else CAA_plots <- NULL

              if(any(data_mean_fit$CAL_hist > 0, na.rm = TRUE)) {
                CAL_plots <- lapply(1:nfleet, individual_array_fn, obs = "data_mean_fit$CAL_hist", pred = "report$CALpred", comps = "length")
              } else CAL_plots <- NULL

              if(any(data_mean_fit$mlen > 0, na.rm = TRUE)) {
                ML_plots <- lapply(1:nfleet, individual_matrix_fn, obs = "data_mean_fit$mlen", pred = "report$mlen_pred",
                                   fig.cap = "mean lengths from fleet", label = "Mean Length from Fleet")
              } else ML_plots <- NULL

              data_section <- c(C_matplot, E_matplot, C_plots, I_plots, CAA_plots, CAL_plots, ML_plots)

              # Model output
              sel_matplot <- rmd_matplot(x = "matrix(data_mean_fit$length_bin, nrow(report$vul_len), nfleet)", y = "report$vul_len", col = "rich.colors(nfleet)",
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
              CAL_bubble <- rmd_bubble("Year", "CAL_all", CAL_bins = "data_mean_fit$length_bin",
                                       fig.cap = "Predicted catch-at-length (summed over all fleets).")

              ts_output <- c(sel_matplot, F_matplot, rmd_SSB(), SSB_plot, rmd_SSB_SSB0(FALSE), rmd_R(),
                             rmd_residual("log_rec_dev", fig.cap = "Time series of recruitment deviations.", label = "log-Recruitment deviations"),
                             rmd_residual("log_rec_dev", "log_rec_dev_SE", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                          label = "log-Recruitment deviations", conv_check = TRUE),
                             rmd_N(), N_bubble, CAA_bubble, CAL_bubble)

              mean_fit_rmd <- c(sumry, LH_section, data_section, ts_output)
            } else mean_fit_rmd <- c("## Fit to mean parameters of OM {.tabset}\n",
                                     "No model found. Re-run `SRA_scope()` with `mean_fit = TRUE`.\n\n")

            if(compare) {
              Hist <- runMSE(OM, Hist = TRUE, parallel = OM@nsim >= 48 & snowfall::sfIsRunning())

              compare_rmd <- c("## Updated OM {.tabset}\n",
                               "### OM historical period\n\n",
                               "```{r, fig.cap = \"Apical F from the operating model.\"}",
                               "Hist_F <- apply(Hist@AtAge$FM, c(1, 3), max)",
                               "matplot(Year, t(Hist_F), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM Apical F\", ylim = c(0, 1.1 * max(Hist_F)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$SSB), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM SSB\", ylim = c(0, 1.1 * max(Hist@TSdata$SSB)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) relative to MSY from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$SSB/Hist@Ref$SSBMSY), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = expression(OM~~SSB/SSB[MSY]), ylim = c(0, 1.1 * max(Hist@TSdata$SSB/Hist@Ref$SSBMSY)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning depletion from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$SSB/Hist@Ref$SSB0), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = expression(OM~~SSB/SSB[0]), ylim = c(0, 1.1 * max(Hist@TSdata$SSB/Hist@Ref$SSB0)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Recruitment from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$Rec), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM Recruitment\", ylim = c(0, 1.1 * max(Hist@TSdata$Rec)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Catch (total removals, including discards) from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$Removals), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM Catch\", ylim = c(0, 1.1 * max(Hist@TSdata$Catch)))",
                               "abline(h = 0, col = \"grey\")",
                               "```\n",
                               "",
                               "### OM/SRA Comparison\n\n",
                               "```{r, fig.cap = \"Difference in apical F between the OM and SRA. Positive values indicate higher F in the OM.\"}",
                               "matplot(Year, t(Hist_F - OM@cpars$Find), typ = \"n\", xlab = \"Year\", ylab = \"Difference in apical F\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t(Hist_F - OM@cpars$Find), col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in spawning biomass (SSB), relative to SSB0, between the OM and SRA, calculated as $(SSB^{OM}_y - SSB^{SRA}_y)/SSB^{OM}_0$. Positive values indicate higher SSB in the OM.\"}",
                               "matplot(Year, t((Hist@TSdata$SSB - x@SSB[sims, 1:OM@nyears, drop = FALSE])/Hist@Ref$SSB0), typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative SSB\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t((Hist@TSdata$SSB - x@SSB[sims, 1:OM@nyears, drop = FALSE])/Hist@Ref$SSB0), col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in recruitment (relative to R0) between the OM and SRA, calculated as $(R^{OM}_y - R^{SRA}_y)/R^{OM}_0$. Positive values indicate higher recruitment in the OM.\"}",
                               "matplot(Year, t(Hist@TSdata$Rec/OM@cpars$R0 - x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]/OM@cpars$R0), typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative recruitment\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, t(Hist@TSdata$Rec/OM@cpars$R0 - x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]/OM@cpars$R0),",
                               "         col = \"black\")",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in annual catch (relative to observed), calculated as $C^{OM}_y/C^{obs}_y - 1$. Positive values indicate higher catch in the OM. Catch in the OM is the total removals (both landings and discards).\"}",
                               "if(any(data$Chist > 0, na.rm = TRUE)) {",
                               "Catch_difference <- t(Hist@TSdata$Removals)/rowSums(data$Chist, na.rm = TRUE) - 1",
                               "Catch_difference[is.infinite(Catch_difference)] <- 0",
                               "matplot(Year, Catch_difference, typ = \"n\", xlab = \"Year\", ylab = \"Difference in relative catch\")",
                               "abline(h = 0, col = \"grey\")",
                               "matlines(Year, Catch_difference, col = \"black\")",
                               "}",
                               "```\n")

            } else compare_rmd <- c("## Updated OM\n",
                                    "Re-run `plot()` function with argument `compare = TRUE`.\n\n")

            rmd <- c(header, OM_update, all_sims_output, mean_fit_rmd, compare_rmd, rmd_footer())
            if(is.list(rmd)) rmd <- do.call(c, rmd)

            # Generate markdown report
            filename_html <- paste0(filename, ".html")
            filename_rmd <- paste0(filename, ".Rmd")
            if(!dir.exists(dir)) {
              message("Creating directory: \n", dir)
              dir.create(dir)
            }

            message("Writing markdown file: ", file.path(dir, filename_rmd))

            write(rmd, file = file.path(dir, filename_rmd))

            # Rendering markdown file
            message("Rendering markdown file to HTML: ", file.path(dir, filename_html))

            output <- rmarkdown::render(file.path(dir, filename_rmd), "html_document", filename_html, dir,
                                        output_options = list(df_print = "paged"), quiet = quiet, ...)
            message("Rendering complete.")

            if(open_file) browseURL(file.path(dir, filename_html))
            invisible(output)
          })



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
           "matplot(length_bin, vul_ff, type = \"l\", col = \"black\", xlab = \"Length\", ylab = \"Selectivity of Fleet ", ff, "\")",
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
           paste0("  Chist <- data$Chist[, ", ff, "]"),
           "  ylim <- c(0.9, 1.1) * range(c(Cpred, Chist), na.rm = TRUE)",
           paste0("  matplot(Year_matrix, Cpred, type = \"o\", pch = 1, col = \"black\", xlab = \"Year\", ylab = \"Catch of Fleet ", ff, "\", ylim = ylim)"),
           paste0("  lines(Year, Chist, col = \"red\", lwd = 3)"),
           "} else {",
           paste0("  Cpred <- do.call(cbind, lapply(report_list, function(x) x$Cpred[, ", ff, "]/mean(x$Cpred[, ", ff, "])))"),
           paste0("  matplot(Year_matrix, Cpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Predicted relative catch of Fleet ", ff, "\")"),
           "}",
           "abline(h = 0, col = \"grey\")",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean ages from fleet ", ff, ".\"}"),
           paste0("MApred <- do.call(cbind, lapply(report_list, function(x) x$CAApred[, , ", ff, "] %*% age/x$CN[, ", ff, "]))"),
           paste0("MAobs <- (data$CAA[, , ", ff, "] %*% age)/rowSums(data$CAA[, , ", ff, "], na.rm = TRUE)"),
           "ylim <- c(0.9, 1.1) * range(c(MApred, MAobs), na.rm = TRUE)",
           "matplot(Year_matrix, MApred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean age\", ylim = ylim)",
           paste0("if(any(data$CAA[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  lines(Year, MAobs, col = \"red\", lwd = 3, typ = \"o\", pch = 16)"),
           "}",
           "```\n",
           "",
           paste0("```{r, fig.cap = \"Observed (red) and predicted (black) mean lengths from fleet ", ff, ".\"}"),
           paste0("MLpred <- do.call(cbind, lapply(report_list, function(x) x$mlen_pred[, ", ff, "]))"),
           paste0("if(any(data$CAL[, , ", ff, "] > 0, na.rm = TRUE)) {"),
           paste0("  MLobs <- (data$CAL[, , ", ff, "] %*% length_bin)/rowSums(data$CAL[, , ", ff, "], na.rm = TRUE)"),
           paste0("} else if(any(data$ML[, ", ff, "] > 0, na.rm = TRUE)) MLobs <- data$ML[, ", ff, "] else MLobs <- NA"),
           "ylim <- c(0.9, 1.1) * range(c(MLpred, MLobs), na.rm = TRUE)",
           "matplot(Year_matrix, MLpred, type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Mean length\", ylim = ylim)",
           "if(!all(is.na(MLobs))) lines(Year, MLobs, col = \"red\", lwd = 3, typ = \"o\", pch = 16)",
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
           paste0("matplot(Year_matrix, Ipred, type = \"l\", col = \"black\", ylim = c(0, 1.1 * max(c(Ipred, data$Index[, ", sur, "]), na.rm = TRUE)), xlab = \"Year\", ylab = \"Survey ", sur, "\")"),
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
    "```\n",
    "",
    "```{r, fig.cap = \"Estimated spawning depletion among all simulations. Unfished spawning biomass is the value calculated from first year life history parameters.\"}",
    "E_E0 <- do.call(cbind, lapply(report_list, function(x) x$E/x$E0_SR))",
    "matplot(Yearplusone_matrix, E_E0, ylim = c(0, 1.1 * max(E_E0)), type = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Spawning depletion\")",
    "abline(h = 0, col = \"grey\")",
    "```\n")
}

rmd_log_rec_dev <- function() {
  c("```{r, fig.cap = \"Estimated recruitment deviations among all simulations.\"}",
    "log_rec_dev <- do.call(cbind, lapply(report_list, getElement, \"log_rec_dev\"))",
    "matplot(Year, log_rec_dev, type = \"n\", xlab = \"Year\", ylab = \"log-recruitment deviations\")",
    "abline(h = 0, col = \"grey\")",
    "matlines(Year, log_rec_dev, col = \"black\")",
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
      matplot(data_val, t(SRA_plot[, i, ]), type = "l", col = "black", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = 'grey')
      if(!is.null(dat)) lines(data_val, dat_plot[i, ], lwd = dat_linewidth, col = dat_color)
      legend("topright", legend = Year[i], bty = "n", xjust = 1)
    }

    if(i == length(Year)) {
      xaxt <- 's'
      if(i %% 16 %in% c(1:4)) yaxt <- 's' else yaxt <- 'n'
      matplot(data_val, t(SRA_plot[, i, ]), type = "l", col = "black", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
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





