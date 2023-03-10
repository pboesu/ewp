#model functions

#' Exponentially weighted poisson regression model
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family choice of "ewp2" or "ewp3"
#' @param data a data frame containing the variables in the model.
#'
#' @return
#' @export
#'
ewp_reg <- function(formula, family = 'ewp3', data){
  mm <- model.matrix(formula, data)

}
