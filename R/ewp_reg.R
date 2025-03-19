#model functions

#' Exponentially weighted Poisson regression model
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family choice of "ewp2" or "ewp3"
#' @param data a data frame containing the variables in the model.
#' @param verbose logical, defaults to TRUE; print model fitting progress
#' @param method string, passed to optim, defaults to 'BFGS'
#' @param hessian logical, defaults to TRUE; calculate Hessian?
#' @param autoscale logical, defaults to TRUE; automatically scale model parameters inside the optimisation routine based on initial estimates from a Poisson regression.
#' @param maxiter numeric, maximum number of iterations for optim
#' @param sum_limit numeric, defaults to 3*maximum count; upper limit for the sum used for the normalizing factor.
#' @param start_val list, defaults to fitting a Poisson regression; specify starting values
#'
#' @return an ewp model
#' @importFrom stats .getXlevels coef delete.response glm.fit model.frame model.matrix model.response na.omit na.pass optim optimHess poisson terms
#' @export
#'

ewp_reg <- function(formula, family = 'ewp3', data, verbose = TRUE, method = 'Nelder-Mead', hessian = TRUE, autoscale = TRUE, maxiter = 500, sum_limit = round(max(Y)*3), start_val=NULL){
  cl <- match.call()
  mt <- terms(formula, data = data)
  #if(missing(data)) data <- environment(formula)
  #mf <- match.call(expand.dots = FALSE)
  #m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  #mf <- mf[c(1, m)]
  #mf$drop.unused.levels <- TRUE

  ## call model.frame()
  #mf[[1]] <- as.name("model.frame")
  #mf <- eval(mf, parent.frame())
  mf <- model.frame(formula, data, drop.unused.levels = TRUE)

  ## model matrix, response
  mm <- model.matrix(formula, mf)
  Y <- model.response(mf, "numeric")

  if(any(Y >= 30)) ("Counts >= 30 detected. The likelihood estimation procedure is not currently set up to deal with this.")
  if(any(Y > 20)) warning("Counts > 20 detected. The likelihood estimation procedure is not currently set up to deal with counts in excess of 30. Results may be misleading if lambda >= 25 and beta2 < 1.")


  # set a warning message if the sum_limit < 15. Occurs for small values of max count (<5), or when sum_limit is set manually.
  if (sum_limit<15){
    warning("sum_limit < 15 detected. A sum_limit of 3 times the maximum count value is recommended, or of
            at least 15, in the case of a small maximum count.")
  }

  #get start values for lambda linpred from simple poisson regression
  if (is.null(start_val)){
    start_values <- coef(glm.fit(x = mm, y = Y, family = poisson()))
  } else{
    start_values <- start_val
  }

  #estimate relative effect sizes for optim - assumes dispersion parameter is approx 1!
  if (autoscale) {
     parscale_est = abs(c(start_values, beta1 = 1, beta2 = 1)/start_values[1])
  } else {
     parscale_est = rep(1, length(start_values)+2)
    }
  #add dispersion parameter start values
  start_values = c(start_values, beta1 = 0, beta2 = 0)
  if(verbose){
    cat('start values are: \n')
    print(start_values)
  }

  pllik3 <- function(par, mm, Y){
    lambda = exp(mm %*% par[1:ncol(mm)])
    beta1 = unname(par['beta1'])
    beta2 = unname(par['beta2'])
    #ll = numeric(nrow(mm))
    #for (i in 1:nrow(mm)){
    #  ll[i] = log(dewp3_cpp(Y[i],lambda[i],beta1,beta2))
    #}
    #return(-1*sum(ll))
    return(pllik3_part_cpp(Y, lambda, beta1, beta2, sum_limit))
  }

  resultp3 <- optim(par = start_values,
                    fn = pllik3, mm = mm, Y = Y,
                    method = method,
                    hessian = FALSE,
                    control = list(trace = verbose,
                                   REPORT=4*verbose,
                                   ndeps=rep(1e-5, ncol(mm)+2),
                                   parscale = parscale_est,
                                   maxit = maxiter))

  if(hessian){
    if(verbose) cat('\nCalculating Hessian. This may take a while.\n')
    resultp3$hessian <- optimHess(resultp3$par, fn = pllik3, mm = mm, Y = Y)
    #estimate vcov
    vc = solve(resultp3$hessian)
  } else {
    vc <- resultp3$hessian <- matrix(NA_real_, nrow = ncol(mm) + 2, ncol = ncol(mm) + 2)
  }



  ## fitted and residuals
  Yhat = exp(mm %*% resultp3$par[1:ncol(mm)])
  res <- (Y - Yhat)

  #output structure
  out <- list(
    coefficients = resultp3$par,
    vcov = vc,
    se = sqrt(diag(vc)),
    optim = resultp3,
    loglik = -resultp3$value,
    residuals = res,
    fitted.values = as.vector(Yhat),
    terms = mt,
    frame = data,
    call = cl,
    levels = .getXlevels(mt, mf),
    start = start_values,
    n = length(Y),
    df.residual = length(Y) - ncol(mm) - 2,
    converged = resultp3$convergence < 1,
    formula = formula,
    sum_limit=sum_limit,
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
#' @return a vector of coefficient values. Beware that the lambda parameters are on the log-link scale, whereas the betas are estimated using an identity link.
#' @importFrom stats coef
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
#' @return a matrix
#' @importFrom stats vcov
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
#' @return a numeric
#' @importFrom stats logLik
#' @export
#'
logLik.ewp <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, nobs = object$n, class = "logLik")
}

#' Extract fitted values
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return a vector of fitted values on the response scale
#' @importFrom stats fitted
#' @export
#'
fitted.ewp <- function(object, ...) {
  object$fitted.values
}


#' Print ewp model object
#'
#' @param x ewp model object
#' @param digits digits to print
#' @param ... ignored
#'
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
#' @importFrom stats pnorm
#' @importFrom utils tail
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
#' @return printout of the summary object
#' @importFrom stats printCoefmat
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


#' Predict from fitted model
#'
#' @param object ewp model object
#' @param newdata optional data.frame
#' @param type character; default="response", no other type implemented
#' @param na.action  defaults to na.pass()
#' @param ... ignored
#'
#' @return a vector of predictions
#' @importFrom stats predict
#' @importFrom stats weighted.mean
#' @export
#'
predict.ewp <- function(object, newdata, type = c("response"),
                              na.action = na.pass, ...)
{
  type <- match.arg(type)

  ## if no new data supplied
  if(missing(newdata)) {
    if(type != "response") {
      stop('Unknown prediction type')
    } else {
      x <- seq(0,object[["sum_limit"]], by= 1)

      pred_ewp <- vector()
      for (i in 1:length(object$fitted.values)){
        pmf_ewp <- dewp3(x, object$fitted.values[i], object$coefficients[["beta1"]],object$coefficients[["beta2"]])
        pred_ewp[i] <- weighted.mean(x,w=(pmf_ewp))
      }
      return(pred_ewp)
    }
  } else {
    mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms), mf)#, contrasts = object$contrasts)
    offset <- if(!is.null(off.num <- attr(object$terms, "offset")))
      eval(attr(object$terms, "variables")[[off.num + 1]], newdata)
    else if(!is.null(object$offset)) eval(object$call$offset, newdata)
    if(is.null(offset)) offset <- rep(0, NROW(X))
  }

  rval <- exp(X %*% object$coefficients[1:ncol(X)])[,1]

  x <- seq(0,object[["sum_limit"]], by= 1)

  pred_ewp <- vector()
  for (i in 1:nrow(newdata)){
    pmf_ewp <- dewp3(x, rval[i], object$coefficients[["beta1"]],object$coefficients[["beta2"]])
    pred_ewp[i] <- weighted.mean(x, w=(pmf_ewp))
  }
  return(pred_ewp)
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


#' simulate from fitted model
#'
#' @param object ewp model object
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param ... ignored
#'
#' @return a data frame with `nsim` columns.
#'
#' @importFrom stats simulate
#' @export
#'
simulate.ewp <- function(object, nsim=1, ...){
  ftd <- fitted(object)
  n <- length(ftd)
  ntot <- n * nsim

  ncoef <- length(coef(object))
  nm <- names(ftd)
  val <- matrix(NA_integer_, nrow = n, ncol = nsim)
  for (i in 1:nsim){
    val[,i] <- vapply(ftd, function(x)rewp3(n = 1, lambda = x, beta1 = coef(object)[ncoef-1], beta2 = coef(object)[ncoef]), numeric(1))
  }

  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  if (!is.null(nm))
    row.names(val) <- nm
  val
}


#' Estimate marginal means
#'
#' @param object ewp model object
#' @param cov character; covariate to find marginal mean for
#' @param ci  logical; whether or not to include confidence intervals, defaults to TRUE
#' @param nsamples  number of samples for use in obtaining the 95% confidence intervals, defaults to 250
#' @param ... ignored
#'
#' @return printout of the marginal means
#' @export
#'

mmean <- function(object,cov,ci=TRUE,nsamples=250,...){

  conts = list()

  coef_ewp = names(object$frame)

  n_fac = length(object$levels)

  facs = names(object$levels)

  cont_var = intersect(names(object$coefficients),coef_ewp)

  n_cont = length(cont_var)

  n_terms=n_fac+n_cont

  levs= vector(mode= "list", length=n_terms)

  cov_cho = intersect(coef_ewp,cov)



  if(is.factor(eval(parse(text=paste0("object$frame$",cov_cho))))==TRUE){

    facs = union(cov_cho,facs)
    full_term = union(facs,cont_var)


    levs[[1]] = levels(eval(parse(text=paste0("object$frame$",cov_cho))))

    if (n_fac>1){
      for (i in 2:n_fac){
        levs[[i]] = levels(eval(parse(text=paste0("object$frame$",facs[i]))))
      }
    }

    if (n_cont>0){
      for (i in (n_fac+1):n_terms){
        levs[[i]] = mean(eval(parse(text=paste0("object$frame$",full_term[i]))))
      }
    }

    RG = do.call("expand.grid",levs)

    names(RG) = full_term
    pred_mat = matrix(predict(object, newdata = RG), nrow = length(levs[[1]]))



    if (ci==TRUE){    ################# CI method

      resample_est <- mvtnorm::rmvnorm(nsamples, mean=object[["coefficients"]],sigma=object[["vcov"]])

      ### Create model frame

      mod_mat <- function(formula, data){
        mf <- model.frame(formula, data, drop.unused.levels = T)
        mm <- model.matrix(formula, mf)
        return(mm)
      }

      mm_und <- mod_mat(formula=delete.response(object$terms), data=RG)

      mmean_boot <-  vector(mode= "list", length=length(levs[[1]]))
      #Yhat <- vector()
      for (k in 1:nsamples){
        Yhat <- exp(mm_und %*% resample_est[k,][1:ncol(mm_und)])

        x <- seq(0,object[["sum_limit"]], by= 1)

        pred_ewp <- vector()

        for (j in 1:length(Yhat)){
          pmf_ewp <- dewp3(x, Yhat[j], resample_est[k,][ncol(resample_est)-1], resample_est[k,][ncol(resample_est)])
          pred_ewp[j] <- weighted.mean(x, w=(pmf_ewp))
        }

        pred_mat_boot = matrix(pred_ewp, nrow = length(levs[[1]]))

        mmean_out_boot = apply(pred_mat_boot, 1, mean)

        for (fac_lev in 1:length(levs[[1]])){
          mmean_boot[[fac_lev]][k] <- mmean_out_boot[fac_lev]
        }
      }

      ci_low <-  vector(mode= "numeric", length=length(levs[[1]]))
      ci_up <-  vector(mode= "numeric", length=length(levs[[1]]))
      for (fac_lev in 1:length(levs[[1]])){
        ci_low[[fac_lev]] <- quantile(mmean_boot[[fac_lev]], prob=0.025,na.rm=T)
        ci_up[[fac_lev]] <- quantile(mmean_boot[[fac_lev]], prob=0.975,na.rm=T)
      }
    }




    if (ci==TRUE){
      emtab = data.frame(levs[[1]],
                         mmean = apply(pred_mat, 1, mean),
                         lower.CL = ci_low,
                         upper.CL = ci_up)
    } else{
      emtab = data.frame(levs[[1]],
                         mmean = apply(pred_mat, 1, mean))
    }

    colnames(emtab)[1] = cov

    print(emtab,row.names = FALSE)

  } else{

    cont_var = union(cov_cho,cont_var)
    full_term = union(cont_var,facs)


    levs[[1]] = c(mean(eval(parse(text=paste0("object$frame$",cov_cho)))),mean(eval(parse(text=paste0("object$frame$",cov_cho))))+0.001)

    if (n_cont>1){
      for (i in 2:n_cont){
        levs[[i]] = mean(eval(parse(text=paste0("object$frame$",cont_var[i]))))
      }
    }

    if (n_fac>0){
      for (i in (n_cont+1):n_terms){
        levs[[i]] = levels(eval(parse(text=paste0("object$frame$",full_term[i]))))
      }
    }

    RG = do.call("expand.grid",levs)


    names(RG) = full_term

    pred_mat = matrix(predict(object, newdata = RG), nrow = length(levs[[1]]))

    if(n_fac>0){
      fac_mean = vector()
      for (i in 1:ncol(pred_mat)){
        fac_mean[i] = (pred_mat[2,i]-pred_mat[1,i])/0.001
      }

      mmean_out = mean(fac_mean)

    } else{
      mmean_out = (pred_mat[2,1]-pred_mat[1,1])/0.001
    }



    if (ci==TRUE){    ################# CI method

      resample_est <- mvtnorm::rmvnorm(nsamples, mean=object[["coefficients"]],sigma=object[["vcov"]])

      ### Create model frame

      mod_mat <- function(formula, data){
        mf <- model.frame(formula, data, drop.unused.levels = T)
        mm <- model.matrix(formula, mf)
        return(mm)
      }

      mm_und <- mod_mat(formula=delete.response(object$terms), data=RG)

      mmean_boot <-  vector(mode= "numeric", length=nsamples)
      #Yhat <- vector()
      for (k in 1:nsamples){
        Yhat <- exp(mm_und %*% resample_est[k,][1:ncol(mm_und)])

        x <- seq(0,object[["sum_limit"]], by= 1)

        pred_ewp <- vector()

        for (j in 1:length(Yhat)){
          pmf_ewp <- dewp3(x, Yhat[j], resample_est[k,][ncol(resample_est)-1], resample_est[k,][ncol(resample_est)])
          pred_ewp[j] <- weighted.mean(x, w=(pmf_ewp))
        }

        pred_mat_boot = matrix(pred_ewp, nrow = length(levs[[1]]))

        if(n_fac>0){
          fac_mean = vector()
          for (i in 1:ncol(pred_mat_boot)){
            fac_mean[i] = (pred_mat_boot[2,i]-pred_mat_boot[1,i])/0.001
          }

          mmean_out_boot = mean(fac_mean)

        } else{
          mmean_out_boot = (pred_mat_boot[2,1]-pred_mat_boot[1,1])/0.001
        }

        mmean_boot[k] <- mmean_out_boot
      }



      # }

      ci_low <- quantile(mmean_boot, prob=0.025,na.rm=T)
      ci_up <- quantile(mmean_boot, prob=0.975,na.rm=T)
    }

    if (ci==TRUE){
      emtab = data.frame(levs[[1]][1],
                         mmean = mmean_out,
                         lower.CL = ci_low,
                         upper.CL = ci_up)
    } else{
      emtab = data.frame(levs[[1]][1],
                         mmean = mmean_out)
    }


    colnames(emtab)[1] = cov

    print(emtab,row.names = FALSE)
  }
}

