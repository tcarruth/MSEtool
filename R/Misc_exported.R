#' What objects of this class are available
#'
#' Generic class finder
#'
#' Finds objects of the specified class in the global environment or in the
#' MSEtool and DLMtool packages. This function is an addendum to the \code{\link[DLMtool]{avail}}
#' function in DLMtool.
#'
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' avail("Assess")
#' avail("HCR")
#' avail("Stock")
#' }
#' @export
avail <- function(classy) {
  temp <- try(class(classy), silent = TRUE)
  if (class(temp) == "try-error") classy <- deparse(substitute(classy))
  if (temp == "function") classy <- deparse(substitute(classy))

  temp <- c(ls("package:MSEtool")[vapply(ls("package:MSEtool"), getclass, logical(1), classy = classy)],
            ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), getclass, logical(1), classy = classy)])
  temp_DLMtool <- try(DLMtool::avail(classy), silent = TRUE)

  if(!inherits(temp_DLMtool, "try-error")) temp <- unique(c(temp, temp_DLMtool))
  if(length(temp) < 1) stop("No objects of class '", classy, "' found", call. = FALSE)
  return(temp)
}

getclass <- function(x, classy) any(inherits(get(x), classy))
