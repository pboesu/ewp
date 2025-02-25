library(ewp)
library(emmeans)
library(broom)
library(tidyverse)

linnet$cov3 <- c("ARC","BLU")
linnet$cov4 <- factor(round(linnet$cov2))

goon <- ewp::ewp_reg(eggs~cov1 ,data=linnet)

linnet$cov3 <- factor(linnet$cov3)

goon <- ewp::ewp_reg(eggs~cov1+cov3+cov4,data=linnet,maxiter=10000)

summary(goon)

avg_cov1 <- mean(linnet$cov1)

avg_cov2 <- mean(linnet$cov2)

new_dat_fitted <- predict(goon, newdata = expand.grid(cov1=c(avg_cov1,avg_cov1+0.001),cov3= unique(linnet$cov3)))

summary(goon)


object=goon


cov_frame <- goon$frame[,names(goon$frame)%in% names(goon$coefficients)]
cov_frame_new <- cov_frame %>% select(!paste(cov))
cov=c("cov1")
cov_frame


predict(goon, newdata=tibble(cov1=c(mean(select(goon$frame,paste(cov))[[1]]),mean(select(goon$frame,paste(cov))[[1]])+0.001),cov2=c(mean(cov_frame_new[,1]),mean(cov_frame_new[,1]))))

mmean <- function(object,cov,...){
  levs= list()
  conts = list()

  coef_ewp = names(object$frame)

  n_fac = length(object$levels)

  facs = names(object$levels)

  cont_var = intersect(names(object$coefficients),coef_ewp)

  n_cont = length(cont_var)

  n_terms=n_fac+n_cont

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


    emtab = data.frame(levs[[1]],
                       mmean = apply(pred_mat, 1, mean))

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


  emtab = data.frame(levs[[1]][1],
                     mmean = mmean_out)

  colnames(emtab)[1] = cov

  print(emtab,row.names = FALSE)
  }
}






mmean(goon,"cov1")

ewp_ci(goon,data=linnet)






