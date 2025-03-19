library(ewp)
library(emmeans)
library(broom)
library(tidyverse)

linnet$cov3 <- c("ARC","BLU")
linnet$cov4 <- factor(round(linnet$cov2))
linnet$cov3 <- factor(linnet$cov3)

goon <- ewp::ewp_reg(eggs~cov1+cov3+cov4 ,data=linnet)


mean(linnet$cov1)

goon <- ewp::ewp_reg(eggs~cov3,data=linnet,maxiter=10000)

summary(goon)

avg_cov1 <- mean(linnet$cov1)

avg_cov2 <- mean(linnet$cov2)

new_dat_fitted <- predict(goon, newdata = expand.grid(cov1=c(avg_cov1,avg_cov1+0.001),cov3= unique(linnet$cov3)))

summary(goon)


object=goon


cov_frame <- goon$frame[,names(goon$frame)%in% names(goon$coefficients)]
cov_frame_new <- cov_frame %>% select(!paste(cov))
cov=c("cov3")
cov_frame


predict(goon, newdata=tibble(cov1=c(mean(select(goon$frame,paste(cov))[[1]]),mean(select(goon$frame,paste(cov))[[1]])+0.001),cov2=c(mean(cov_frame_new[,1]),mean(cov_frame_new[,1]))))

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






mmean(goon,"cov4",ci=TRUE)

ewp_ci(goon,data=linnet)






