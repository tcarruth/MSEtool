summary_spict <- function(Assessment) {
  assign_Assessment_slots(Assessment)

  input_parameters <- data.frame()

  if(conv) {
    current_status <- data.frame(Value = c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)],
                                           B_B0[length(B_B0)]))
    rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")
  } else current_status <- data.frame()

  if(!is.character(SD)) {
    derived <- data.frame(Value = c(SD$value["r"], TMB_report$Bmsy, TMB_report$Fmsy, exp(SD$value["logalpha"]), exp(SD$value["logbeta"]),
                                    BMSY/B0),
                          Description = c("Intrinsic rate of population increase", "Fishing mortality at MSY", "Biomass at MSY",
                                          "Ratio of index and biomass standard deviations", "Ratio of catch and F standard deviations",
                                          "Depletion at MSY"),
                          stringsAsFactors = FALSE)
    rownames(derived) <- c("r", "K", "BMSY", "alpha", "beta", "BMSY/B0")

    model_estimates <- summary(SD, "fixed")
    SD_exp <- model_estimates
    SD_exp[, 1] <- exp(SD_exp[, 1])
    SD_exp[, 2] <- SD_exp[, 1] * SD_exp[, 2]
    rownames(SD_exp) <- c("MSY", "K", "q", "n", "sd_b", "sd_f", "sd_i", "sd_c")
    model_estimates <- rbind(model_estimates, SD_exp)
  } else {
    current_status <- derived <- model_estimates <- data.frame()
  }

  output <- list(model = "SPiCT", current_status = current_status,
                 input_parameters = input_parameters, derived_quantities = derived,
                 model_estimates = model_estimates)
  return(output)
}

rmd_spict <- function(Assessment) {
  ss <- rmd_summary("SPiCT")

  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"), rmd_data_timeseries("Index"))

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_spict_FMSY(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_MSY("logm"),
                  rmd_F_FMSY_terminal(), rmd_B_BMSY_terminal(),
                  rmd_assess_fit("Index", "index"), rmd_assess_resid("Index"), rmd_assess_qq("Index", "index"),
                  rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"))

  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n"), rmd_F_FMSY(FALSE), rmd_B(), rmd_B_BMSY(FALSE),
                 rmd_B_B0(FALSE), rmd_Kobe("B_BMSY", xlab = "expression(B/B[MSY])", conv_check = FALSE))

  productivity <- c(rmd_spict_yield_header(), rmd_yield_F("SP", FALSE), rmd_yield_depletion("SP", FALSE), rmd_sp(FALSE))

  return(c(ss, data_section, assess_fit, ts_output, productivity))
}

rmd_spict_FMSY <- function(header = NULL) {
  fig.cap <- "Estimate of FMSY, distribution based on normal approximation of estimated covariance matrix."
  ans <- c(paste0("```{r, fig.cap=\"", fig.cap, "\"}"),
           "if(conv) plot_normalvar(FMSY, SE_FMSY, label = expression(hat(F)[MSY]))",
           "```\n")
  if(!is.null(header)) ans <- c(header, ans)
  return(ans)
}

rmd_spict_yield_header <- function(header = NULL) {
  c("### Productivity\n",
    "```{r}",
    "TMB_report <- c(TMB_report, list(K = B0, BMSY = BMSY, n = exp(SD$par.fixed[\"logn\"])))",
    "```\n")
}


profile_likelihood_spict <- function(Assessment, ...) {
  stop("Profiling currently not supported for spict in MSEtool.", call. = FALSE)
}



retrospective_spict <- function(Assessment, nyr) {
  stop("Retrospective analysis is currently not supported for spict in MSEtool, but is available in the spict package.", call. = FALSE)
}


