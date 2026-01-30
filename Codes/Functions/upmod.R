#' Update `asreml` model
#'
#' Iteratively update `asreml` models until no component changes more than 1%.
#'
#' @param model An `asreml` object.
#'
#' @return An update and fully convergend `asreml` object.
#'
#' @author Saulo Chaves
#'
#' @export

up.mod = function (model) 
{
  if (model$converge) {
    repeat {
      if (any(stats::na.exclude(model$vparameters.pc) >= 1)) {
        model = suppressWarnings(asreml::update.asreml(model))
      }
      else {
        message("The model was updated and reached full convergence")
        break
      }
    }
  }
  else message("The model failed to converge")
  return(model)
}
