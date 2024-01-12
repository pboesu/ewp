#### Sensitivity analysis for sum limit

library(tidyverse)

## create a loop for applying the sum limit, testing both running times and coefficient estimates

fix_sum <- as.vector(round(seq(from = 8, to = 20, by=1)))

ewp_vary <- vector(mode = "list", length = length(fix_sum))
for (i in 1:length(fix_sum)){
  start.time <- proc.time()
  ewp_out <- ewp_reg(eggs~cov2,method='Nelder-Mead',data=linnet,sum_limit = fix_sum[i])
  end.time <- proc.time()
  ewp_vary[[i]]$coef_ewp <- coef(ewp_out)
  ewp_vary[[i]]$elapsed <- (end.time-start.time)[3] ##calculate elapsed time
  ewp_vary[[i]]$sum_limit <- fix_sum[i]
}

ewp_reg(eggs~cov2,method='Nelder-Mead',data=linnet,sum_limit = 30)
#save(ewp_vary,file="sensitivity_NM.rda")

time <- vector(mode = "numeric", length = length(fix_sum))
sum_lim <- vector(mode = "numeric", length = length(fix_sum))
for (i in 1:length(fix_sum)){
  time[i] <- as.numeric(ewp_vary[[i]][["elapsed"]])
  sum_lim[i] <- ewp_vary[[i]][["sum_limit"]]
}

time_df <- data.frame(time,sum_lim)


##plot time taken against limit of the sum of the normalising factor
ggplot(time_df, aes(x=sum_lim, y=time))+
  geom_point()+
  theme_bw(base_size = 15)


###########################################################################################
## Simulate count data

set.seed(1234)
eggs <- rewp3(n=5414,lambda=4.88,beta1=1.46,beta2=1,sum_limit=100)
eggs <- as.data.frame(eggs)
eggs$cov2 <- linnet$cov2


fix_sum <- as.vector(round(seq(from = max(eggs$eggs)+1, to = max(eggs$eggs)+20, by=1)))

ewp_vary <- vector(mode = "list", length = length(fix_sum))
for (i in 1:length(fix_sum)){
  start.time <- proc.time()
  ewp_out <- ewp_reg(eggs~cov2,method='Nelder-Mead',data=eggs,sum_limit = fix_sum[i])
  end.time <- proc.time()
  ewp_vary[[i]]$coef_ewp <- coef(ewp_out)
  ewp_vary[[i]]$elapsed <- (end.time-start.time)[3] ##calculate elapsed time
  ewp_vary[[i]]$sum_limit <- fix_sum[i]
}

ewp_reg(eggs~cov2,method='Nelder-Mead',data=eggs,sum_limit = 100)

save(ewp_vary,file="sensitivity_NM_lambda=4.88_b2=1.rda")

