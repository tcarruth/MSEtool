
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


retrospective_SCA2 <- function(Assessment, nyr, figure = TRUE) retrospective_SCA(Assessment, nyr, figure, TRUE)
