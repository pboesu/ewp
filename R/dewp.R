#likelihood functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

w_k = function(beta, k, lambda){
  exp(-beta*abs(k-lambda))
}#EWP_2 weights

W_inner = function(beta, k, lambda){
  exp(-lambda)*lambda^k*w_k(beta, k, lambda)/factorial(k)
}

W = function(beta, lambda, sum_limit = 30){
  sum(W_inner(beta,0:sum_limit,lambda))
}

w_k3 = function(beta1, beta2, k, lambda){
  ifelse(k <= lambda, exp(-beta1*(lambda-k)),exp(-beta2*(k-lambda)))
}#EWP_3 weights

W_inner3 = function(beta1, beta2, k, lambda){
  exp(-lambda)*lambda^k*w_k3(beta1, beta2, k, lambda)/factorial(k)
}


W3 = function(beta1, beta2, lambda, sum_limit = 30){
  sum(W_inner3(beta1, beta2, 0:sum_limit,lambda))
}


#' Probability mass function of the two-parameter EWP
#'
#' @param x vector of (positive integer) quantiles.
#' @param lambda centrality parameter
#' @param beta dispersion parameter
#'
#' @return
#' @export
#'
dewp2 <- function(x, lambda, beta){
  stopifnot(is.wholenumber(x))
  exp(-lambda)*lambda^x*w_k(beta, k=x, lambda)/(W(beta, lambda)*factorial(x))
}

#' Probability mass function of the three-parameter EWP
#'
#' @param x vector of (positive integer) quantiles.
#' @param lambda centrality parameter
#' @param beta1 lower-tail dispersion parameter
#' @param beta2 upper tail dispersion parameter
#'
#' @return
#' @export
#'
dewp3 <- function(x, lambda, beta1, beta2){
  stopifnot(is.wholenumber(x))
  exp(-lambda)*lambda^x*w_k3(beta1, beta2, k=x, lambda)/(W3(beta1, beta2, lambda)*factorial(x))
}
