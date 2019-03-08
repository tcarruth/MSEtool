#' What objects of this class are available
#'
#' Generic class finder
#'
#' Finds objects of the specified class in the global environment or in the
#' MSEtool and DLMtool packages. This function is an addendum to the \code{\link[DLMtool]{avail}}
#' function in DLMtool.
#'
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @param all_avail Logical. If TRUE, function will return all objects of class \code{classy} available to user.
#' If FALSE, returns only those objects included in MSEtool.
#' @author Q. Huynh
#' @examples
#' \donttest{
#' avail("Assess")
#' avail("HCR")
#' avail("Stock")
#' avail("MP")
#' avail("MP", all_avail = FALSE)
#' }
#' @export
avail <- function(classy, all_avail = TRUE) {
  temp <- try(class(classy), silent = TRUE)
  if(class(temp) == "try-error") classy <- deparse(substitute(classy))
  if(temp == "function") classy <- deparse(substitute(classy))

  temp <- ls("package:MSEtool")[vapply(ls("package:MSEtool"), getclass, logical(1), classy = classy)]

  if(all_avail) {
    temp_globalenv <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), getclass, logical(1), classy = classy)]
    temp <- c(temp, temp_globalenv)

    temp_DLMtool <- try(DLMtool::avail(classy), silent = TRUE)
    if(!inherits(temp_DLMtool, "try-error")) temp <- unique(c(temp, temp_DLMtool))
  }

  if(length(temp) < 1) stop("No objects of class '", classy, "' found", call. = FALSE)
  return(temp)
}

getclass <- function(x, classy) any(inherits(get(x), classy))

#' Get the MSEtool vignettes
#'
#' A convenient function to open a web browser with the MSEtool package vignettes
#' @examples
#' \dontrun{
#' MSEtool::userguide()
#' DLMtool::userguide()
#' }
#' @seealso \link[DLMtool]{userguide}
#' @export
userguide <- function() {
  message("For the DLMtool user guide, type in 'DLMtool::userguide()' to the console.")
  browseVignettes("MSEtool")
}
