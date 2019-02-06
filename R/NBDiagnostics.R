# Fit the Negtive Binomial (NB) model
#
# Functions for fitting an NB model


#' Diagnostic functionality for NB models
#'
#' @param formula,data,...  arguments for the \code{glm.nb()} function
#' @param id_varname  variable name for subject ID
#' @param outcome_varname  variable name for the outcome count/rate
#' @param baseline_varname  variable name for the baseline count/rate
#' @param group_varname  variable name for group
#' @param addition  addition for log-transformation
#' @description Diagnostic functionality for NB models
#'
#' @return \code{nbdiagnostics()} returns an object of class \code{NBDiagnostics}
#' inheriting from \code{negbin}, \code{glm}, and \code{lm}.
#' @export
#' @import MASS
#'
#' @examples
#' ## simulate a dataset with outcome and baseline count
#' obs <- 200
#' s <- rgamma(obs, shape = 1/3, scale = 3)
#' y0 <- rpois(obs, 25 * s)
#' x <- rep(c(1, 0), each=obs/2)
#' y1 <- rpois(obs, 25 * s * exp(-0.7 * x))
#' x <- factor(x, levels = c(1, 0),
#'             labels = c('Intervention', 'Control'))
#' dat <- data.frame(id=1:obs, y0=y0, y1=y1, group=x)
#'
#' ## fit an NB model using nbdiagnostics() and specify
#' ## the varaibles names for the plots
#' mod_nb <- nbdiagnostics(
#'   y1 ~ group, data=dat, id_varname = "id",
#'   outcome_varname="y1", baseline_varname="y0",
#'   group_varname = "group")
#'
#' ## produce the baseline/outcome event plots
#' boeplot(mod_nb, diagnostic_stat = "cookd")
#' boeplot(mod_nb, diagnostic_stat = "leverage")
#' boeplot(mod_nb, diagnostic_stat = "anscombe_resid")
#' boeplot(mod_nb, diagnostic_stat = "dfbeta")
#'
#' ## covariate adjusted probability plots (Holling et al, 2016)
#' caprob_nb(20, mod_nb)
#' caprob_nb_poi(20, mod_nb, data=dat)

nbdiagnostics <- function(
  formula,
  data, ...,
  id_varname,
  outcome_varname,
  baseline_varname,
  group_varname,
  addition=0.5
) {
  col_vars = unique(
    c(id_varname, baseline_varname, outcome_varname, group_varname,
      all.vars(formula))
    )
  data <- na.omit(data[col_vars])
  nb_model <- glm.nb(formula=formula, data=data, ...)

  class(nb_model) <- c("NBDiagnostics", class(nb_model))
  nb_model$id_varname <- id_varname
  nb_model$outcome_varname <- outcome_varname
  nb_model$baseline_varname <- baseline_varname
  nb_model$group_varname <- group_varname
  nb_model$addition <- addition
  nb_model$data <- data

  return(nb_model)
}
