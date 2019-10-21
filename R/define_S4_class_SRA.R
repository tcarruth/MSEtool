
#' Class-\code{SRA}
#'
#' An S4 class for the output from \link{SRA_scope}.
#'
#' @name SRA-class
#' @docType class
#'
#' @slot OM An updated operating model, class \linkS4class{OM}.
#' @slot SSB An matrix of estimated spawning biomass with \code{OM@@nsim} rows and \code{OM@@nyears+1} columns.
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
                                 rmd_SRA_initD(), rmd_SRA_R_output(), rmd_SRA_SSB_output())

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
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Spawning biomass (SSB) from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$SSB), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM SSB\", ylim = c(0, 1.1 * max(Hist@TSdata$SSB)))",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Recruitment from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$Rec), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM Recruitment\", ylim = c(0, 1.1 * max(Hist@TSdata$Rec)))",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Catch from the operating model.\"}",
                               "matplot(Year, t(Hist@TSdata$Catch), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"OM Catch\", ylim = c(0, 1.1 * max(Hist@TSdata$Catch)))",
                               "```\n",
                               "",
                               "### OM/SRA Comparison\n\n",
                               "```{r, fig.cap = \"Difference in apical F between the OM and SRA. Positive values indicate higher F in the OM.\"}",
                               "matplot(Year, t(Hist_F - OM@cpars$Find), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Difference in apical F\")",
                               "abline(h = 0, lwd = 3, lty = 2)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in spawning biomass (SSB) (relative to SSB0) between the OM and SRA. Positive values indicate higher SSB in the OM.\"}",
                               "matplot(Year, t((Hist@TSdata$SSB - x@SSB[sims, 1:OM@nyears, drop = FALSE])/Hist@Ref$SSB0), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Difference in relative SSB\")",
                               "abline(h = 0, lwd = 3, lty = 2)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in recruitment (relative to R0) between the OM and SRA. Positive values indicate higher recruitment in the OM.\"}",
                               "matplot(Year, t(Hist@TSdata$Rec/OM@cpars$R0 - x@NAA[, 1:OM@nyears, 1][sims, , drop = FALSE]/OM@cpars$R0), typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Difference in relative recruitment\")",
                               "abline(h = 0, lwd = 3, lty = 2)",
                               "```\n",
                               "",
                               "```{r, fig.cap = \"Difference in catch (relative to time series mean) between the OM and SRA. Positive values indicate higher catch in the OM.\"}",
                               "if(any(data$Chist > 0, na.rm = TRUE)) {",
                               "Catch_difference <- t(Hist@TSdata$Catch)/rowSums(data$Chist, na.rm = TRUE) - 1",
                               "Catch_difference[is.infinite(Catch_difference)] <- 0",
                               "matplot(Year, Catch_difference, typ = \"l\", col = \"black\", xlab = \"Year\", ylab = \"Difference in relative catch\")",
                               "abline(h = 0, lwd = 3, lty = 2)",
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


