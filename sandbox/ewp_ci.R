ewp_ci <- function(object, data, nsamples=250, year=F){
  resample_est <- mvtnorm::rmvnorm(nsamples, mean=object[["coefficients"]],sigma=object[["vcov"]])

  ### Create model frame

  mod_mat <- function(formula, data){
    mf <- model.frame(formula, data, drop.unused.levels = T)
    mm <- model.matrix(formula, mf)
    return(mm)
  }

  mm_und <- mod_mat(formula=object$formula, data=data)

  if (year==F){
    lambda_und <- vector()
    #Yhat <- vector()
    for (i in 1:nsamples){
      Yhat <- exp(mm_und %*% resample_est[i,][1:ncol(mm_und)])

      x <- seq(0,object[["sum_limit"]], by= 1)
      pred_ewp <- vector()
      for (j in 1:length(Yhat)){
        pmf_ewp <- dewp3(x, Yhat[j], resample_est[i,][ncol(resample_est)-1], resample_est[i,][ncol(resample_est)])
        pred_ewp[j] <- weighted.mean(x, w=(pmf_ewp))
      }
      lambda_und[i] <- mean(pred_ewp)
    }

    ci_low <- quantile(lambda_und, prob=0.025,na.rm=T)
    ci_up <- quantile(lambda_und, prob=0.975,na.rm=T)

    ci_out <- c(ci_low,ci_up)
  } else if (year==T){
    lambda_und <- list()
    #Yhat <- vector()
    for (i in 1:nsamples){
      Yhat <- exp(mm_und %*% resample_est[i,][1:ncol(mm_und)])

      x <- seq(0,object[["sum_limit"]], by= 1)
      pred_ewp <- vector()
      for (j in 1:length(Yhat)){
        pmf_ewp <- dewp3(x, Yhat[j], resample_est[i,][ncol(resample_est)-1], resample_est[i,][ncol(resample_est)])
        pred_ewp[j] <- weighted.mean(x, w=(pmf_ewp))
      }
      data$pred <- pred_ewp
      lambda_und[[i]] <- data %>% group_by(yearf) %>% summarise(m_pred=mean(pred))
    }


    lambda_join <- do.call(rbind,lambda_und)


    ci_low_year <- lambda_join %>% group_by(yearf) %>%
      summarize(low=quantile(m_pred,prob=0.025,na.rm=T))

    ci_up_year <- lambda_join %>% group_by(yearf) %>%
      summarize(up=quantile(m_pred,prob=0.975,na.rm=T))

    ci_year <- left_join(ci_low_year,ci_up_year,by='yearf')

    ci_out <- ci_year
  }

  return(ci_out)

}
