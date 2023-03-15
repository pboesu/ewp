#model functions

#' Exponentially weighted poisson regression model
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family choice of "ewp2" or "ewp3"
#' @param data a data frame containing the variables in the model.
#' @param verbose logical, defaults to TRUE; print model fitting progress
#' @param hessian logical, defaults to TRUE; calculate Hessian?
#'
#' @return
#' @export
#'
ewp_reg <- function(formula, family = 'ewp3', data, verbose = TRUE, method = 'BFGS', hessian = TRUE){
  mm <- model.matrix(formula, data)
  mf <- model.frame(formula, data)
  X <- model.response(mf, "numeric")
  #get start values for lambda linpred from simple poisson regression
  start_values <- coef(glm(formula = formula, data = data, family = poisson))
  #add dispersion parameter start values
  start_values = c(start_values, beta1 = 0, beta2 = 0)
  if(verbose){
    cat('start values are: \n')
    print(start_values)
  }

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
                    method = method,
                    hessian = FALSE,
                    control = list(trace = verbose,REPORT=2*verbose,ndeps=rep(1e-5, ncol(mm)+2)))

  if(hessian){
    if(verbose) cat('\nCalculating Hessian. This may take a while.\n')
    resultp3$hessian <- optimHess(resultp3$par, fn = pllik3, mm = mm, X = X)
  }
  #estimate vcov
  vc = solve(resultp3$hessian)

  #TODO: rename coefficients to highlight lambda linpred? and/or separate mean and dispersion coefficients in output

  #output structure
  out <- list(
    coefficients = resultp3$par,
    vcov = vc,
    se = sqrt(diag(vc)),
    optim = resultp3,
    loglik = -resultp3$value,
    residuals = NA,
    fitted.values = NA,
    start = start_values,
    n = length(X),
    df.residual = length(X) - ncol(mm) - 2,
    converged = resultp3$convergence < 1,
    formula = formula,
    dist= 'ewp3'
  )
  class(out) <- "ewp"
  return(out)
}

#' Extract coefficients
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return
#' @export
#'
coef.ewp <- function(object, ...) {
  object$coefficients
}

#' Extract estimated variance-covariance matrix
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return
#' @export
#'
vcov.ewp <- function(object, ...) {
  object$vcov
}

#' Extract log likelihood
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return
#' @export
#'
logLik.ewp <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, nobs = object$n, class = "logLik")
}


#' Print ewp model object
#'
#' @param x ewp model object
#' @param digits digits to print
#' @param ... ignored
#'
#' @return
#' @export
#'
print.ewp <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    dist <- x$dist
    fixed <- FALSE

  #cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
      } else {
    cat(paste("Coefficients (", dist, " with log link on lambda):\n", sep = ""))
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    }

  invisible(x)
}

#' Model summary
#'
#' @param object ewp model fit
#' @param ... ignored
#'
#' @return
#' @export
#'
summary.ewp <- function(object,...)
{
  ## deviance residuals
  #object$residuals <- residuals(object, type = "deviance")

  ## compute z statistics
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  k <- length(cf)


  zstat <- cf/se
  pval <- 2*pnorm(-abs(zstat))
  cf <- cbind(cf, se, zstat, pval)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf

  ## number of iterations
  object$iterations <- tail(na.omit(object$optim$count), 1)

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.ewp"
  object
}

#' Print ewp model summary
#'
#' @param x ewp model summary
#' @param digits number of digits to print
#' @param ... additional arguments to printCoefmat()
#'
#' @return
#' @export
#'
print.summary.ewp <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  #cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {


      dist <- x$dist
      fixed <- FALSE
      npar_lambda <- nrow(x$coefficients) - 2


    cat("Deviance residuals:\n")
    #print(structure(quantile(x$residuals),
    #                names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)

    cat(paste("\nlambda coefficients (", dist, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients[1:npar_lambda, , drop = FALSE], digits = digits, ...)
    cat(paste("\ndispersion coefficients:\n"))
    printCoefmat(x$coefficients[(npar_lambda+1):(npar_lambda+2),], digits = digits, ...)

    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }

  invisible(x)
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
