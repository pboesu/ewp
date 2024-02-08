#######################################################################################################
####Original R version of function code
#######################################################################################################

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


W3 = function(beta1, beta2, lambda, x, sum_limit){
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

#######################################################################################################

microbenchmark::microbenchmark(r = dewp3(1,2,3,4),
                               cpp = dewp3_cpp(1,2,3,4),
                               times = 500)


library(ewp)
system.time(fit_null <- ewp_reg(eggs ~ 1, data = linnet))
print(fit_null)
summary(fit_null)


system.time(fit <- ewp_reg(eggs ~ cov1 + cov2, data = linnet))
cbind(fit$coefficients, fit$se)
print(fit)
summary(fit)
print(summary(fit))

plot(linnet$eggs,predict(fit))
abline(0,1, col = 'red')
plot(linnet$eggs,fit$residuals)

plot(linnet$eggs,predict(fit_null))
abline(0,1, col = 'red')
plot(linnet$eggs,fit_null$residuals)

system.time(fitNM <- ewp_reg(eggs ~ cov1 + cov2, data = linnet, method='Nelder-Mead'))
print(fitNM)

library(dplyr)
piedfly <- readRDS('../../2023_exponentially_weighted_poisson/data/piefl.rds') %>%
     select(max_num_eggs, min_first_egg, yearf,northing) %>%
     mutate(northing_scl = scale(northing)[,1],
            min_first_egg_scl = scale(min_first_egg)[,1])

m1<-glm(max_num_eggs~min_first_egg_scl + yearf + northing_scl, family = poisson, data = piedfly)
summary(m1)

piedfly_s <- piedfly %>%
  group_by(yearf) %>%
  slice_sample(n = 200)

system.time(m2 <- ewp_reg(max_num_eggs~min_first_egg_scl + yearf + northing_scl, data = piedfly_s, hessian = FALSE))#very slow - does not converge in 100 iterations / 42 minutes
system.time(m2b <- ewp_reg(max_num_eggs~min_first_egg_scl + yearf + northing_scl, data = piedfly_s, hessian = T))
system.time(m3 <- ewp_reg(max_num_eggs~min_first_egg + northing_scl, data = piedfly_s, hessian = FALSE))#hessian slowdown (unsurprisingly) caused by year factors - optimisation may not be stable for large sample sizes - maybe to do with parscale - maybe try Nelder-Mead instead?
system.time(m3b <- ewp_reg(max_num_eggs~min_first_egg + northing_scl, data = piedfly_s, hessian = T))
summary(m3b)


#simulate methods
fitp <- glm(eggs ~ cov1 + cov2, data = linnet, family=poisson)
fitqp <- glm(eggs ~ cov1 + cov2, data = linnet, family=quasipoisson)
fite <- ewp_reg(eggs ~ cov1 + cov2, data = linnet)
coef(fitp)
coef(fitqp)
coef(fite)

simp <- simulate(fitp, nsim = 10, seed = 1234)
str(simp)
simqp <- simulate(fitqp, nsim = 10, seed = 1234)
sime <- simulate(fite, nsim = 1)

sim_points <- function(sim, ...){
  simt <- lapply(sime, table)
  lapply(simt, function(x)points(as.numeric(names(x)),x,...))
}

hist(linnet$eggs)
sim_points(sime, col = 'red', pch = 16)
