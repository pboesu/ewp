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
  mf <- model.frame(formula, data)
  X <- model.response(mf, "numeric")
  #get start values for lambda linpred from simple poisson regression
  start_values <- coef(glm(formula = formula, data = data, family = poisson))
  #add dispersion parameter start values
  start_values = c(start_values, beta1 = 0, beta2 = 0)
  print('start values are: \n')
  print(start_values)

  pllik3 <- function(par, mm, X){
    lambda = exp(mm %*% par[1:ncol(mm)])
    beta1 = unname(par['beta1'])
    beta2 = unname(par['beta2'])
    ll = numeric(nrow(mm))
    for (i in 1:nrow(mm)){
      ll[i] = log(dewp3(X[i],lambda[i],beta1,beta2))
    }
    return(-1*sum(ll))
  }

  resultp3 <- optim(par = start_values,
                    fn = pllik3, mm = mm, X = X,
                    method = 'BFGS',
                    hessian = TRUE,
                    control = list(trace = 1,REPORT=2))
  return(resultp3)
}


# ## covariances
# vc <- if(hessian) {
#   tryCatch(-solve(as.matrix(fit$hessian)),
#            error = function(e) {
#              warning(e$message, call = FALSE)
#              matrix(NA_real_, nrow = k + (dist == "negbin"), ncol = k + (dist == "negbin"))
#            })
# } else {
#   matrix(NA_real_, nrow = k + (dist == "negbin"), ncol = k + (dist == "negbin"))
# }
# if(dist == "negbin") {
#   SE.logtheta <- as.vector(sqrt(diag(vc)[k + 1]))
#   vc <- vc[-(k+1), -(k+1), drop = FALSE]
# } else {
#   SE.logtheta <- NULL
# }
# colnames(vc) <- rownames(vc) <- colnames(X)
