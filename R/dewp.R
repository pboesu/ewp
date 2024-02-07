#likelihood functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

w_k = function(beta, k, lambda){
  exp(-beta*abs(k-lambda))
}#EWP_2 weights

W_inner = function(beta, k, lambda){
  exp(-lambda)*lambda^k*w_k(beta, k, lambda)/factorial(k)
}

W = function(beta, lambda, sum_limit = 30){#TODO:
  sum(W_inner(beta,0:sum_limit,lambda))
}

w_k3 = function(beta1, beta2, k, lambda){
  ifelse(k <= lambda, exp(-beta1*(lambda-k)),exp(-beta2*(k-lambda)))
}#EWP_3 weights

W_inner3 = function(beta1, beta2, k, lambda){
  exp(-lambda)*lambda^k*w_k3(beta1, beta2, k, lambda)/factorial(k)
}


W3 = function(beta1, beta2, lambda, sum_limit = max(x)*2.5){
  sum(W_inner3(beta1, beta2, 0:sum_limit,lambda))
}


#' Probability mass function of the two-parameter EWP
#'
#' @param x vector of (positive integer) quantiles.
#' @param lambda centrality parameter
#' @param beta dispersion parameter
#'
#' @return a vector of probabilities
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
#' @return a vector of probabilities
#' @export
#'
dewp3 <- function(x, lambda, beta1, beta2){
  stopifnot(is.wholenumber(x))
  exp(-lambda)*lambda^x*w_k3(beta1, beta2, k=x, lambda)/(W3(beta1, beta2, lambda)*factorial(x))
}

#' Random samples from the three-parameter EWP
#'
#' @param n number of observations
#' @param lambda centrality parameter
#' @param beta1 lower-tail dispersion parameter
#' @param beta2 upper tail dispersion parameter
#' @param sum_limit largest integer to evaluate the PMF sum for ()
#'
#' @return random deviates from the EWP_3 distribution
#' @importFrom stats runif
#' @export
#'
rewp3 <- function(n, lambda, beta1, beta2, sum_limit=30){#TODO:sum_limit should be a package option, and also there should be some way of automatically ensuring this is appropriate based on lambda and/or input data
    if(any(lambda >= sum_limit)) stop('sum_limit must be larger than lambda')
    probs <- vapply(1:sum_limit, dewp3_cpp, numeric(1), lambda, beta1, beta2, sum_limit)
    sample(x = 1:sum_limit, n, replace = T, prob = probs)
}
