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

#try_ewp <- ewp_reg(eggs~cov2,method='Nelder-Mead',data=linnet, sum_limit = 25)
#dewp3(linnet$eggs, lambda=4.88,beta1=1,beta2=3,sum_limit=10)

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
# simulate a few random sets. through n.iter

n.iter <- 100
set.seed(1234)
eggs <- vector(mode="list", length=n.iter)
fix_sum <- vector(mode="list", length=n.iter)
ewp_vary <- vector(mode="list", length=n.iter)
ewp_out <- vector(mode="list", length=n.iter)
true_out <- vector(mode="list", length=n.iter)
true_coef <- vector(mode="list", length=n.iter)

for (j in 1:n.iter){
eggs[[j]] <- rewp3(n=5414,lambda=4.88,beta1=1.46,beta2=1,sum_limit=100)
eggs[[j]] <- as.data.frame(eggs[[j]])
eggs[[j]] <- rename(eggs[[j]],'eggs'='eggs[[j]]')
eggs[[j]]$cov2 <- linnet$cov2

fix_sum[[j]] <- as.vector(round(seq(from = max(eggs[[j]]$eggs)+1, to = max(eggs[[j]]$eggs)+20, by=1)))
}


for (j in 1:n.iter){
  for (i in 1:length(fix_sum[[j]])){
  start.time <- proc.time()
  ewp_out[[j]][[i]] <- ewp_reg(eggs~cov2,method='Nelder-Mead',data=eggs[[j]],sum_limit = fix_sum[[j]][i])
  end.time <- proc.time()
  ewp_vary[[j]][[i]] <- coef(ewp_out[[j]][[i]])
  ewp_vary[[j]][[i]] <- as.data.frame(ewp_vary[[j]][[i]])
  ewp_vary[[j]][[i]] <- rename(ewp_vary[[j]][[i]],'coef_ewp'='ewp_vary[[j]][[i]]')
  ewp_vary[[j]][[i]]$sum_limit <- fix_sum[[j]][i]
  ewp_vary[[j]][[i]]$elapsed <- (end.time-start.time)[3] ##calculate elapsed time
  }
}

for (j in 1:n.iter){
  true_out[[j]] <- ewp_reg(eggs~cov2,method='Nelder-Mead',data=eggs[[j]],sum_limit = 100)
  true_coef[[j]] <- coef(true_out[[j]])
}


save(ewp_vary,file="sensitivity_NM_lambda=4.88_b2=1_100simdat.rda")
save(true_coef,file="true_coef_NM_lambda=4.88_b2=1_100simdat.rda")




##################################################################################################

## Plot the difference between the true value and the input value in the simulation
## load the appropriate output

### for a varying beta2

##set the true value
true_beta <- 1
## go through list
beta_2 <- vector(mode="list", length=n.iter)
for (i in 1:n.iter){
  for (j in 1:20){
    beta_2[[i]][[j]] <- (ewp_vary[[i]][[j]][["coef_ewp"]][4]-true_beta)/ewp_vary[[i]][[j]][["coef_ewp"]][4]
    beta_2[[i]][[j]] <- as.data.frame(beta_2[[i]][[j]])
    beta_2[[i]][[j]] <- rename(beta_2[[i]][[j]],'beta2'='beta_2[[i]][[j]]')
    beta_2[[i]][[j]]$sum_limit <- ewp_vary[[i]][[j]][["sum_limit"]][1]
  }
}

beta_2_1 <- bind_rows(beta_2)

beta_2_1$set <- rep(1:n.iter,each=20)
beta_2_1$input <- 1.5


beta_2_full <- bind_rows(beta_2_1,beta_2_1.5,beta_2_2.36,beta_2_4,beta_2_5)
beta_2_full$input <- as.factor(beta_2_full$input)
beta_2_full$set <- as.factor(beta_2_full$set)

beta_2_full$aftmax <- rep(1:20, length.out=length(beta_2_full$set))
beta_2_full$aftmax <- as.factor(beta_2_full$aftmax)


ggplot(beta_2_full, aes(x=aftmax, y=beta2))+
  geom_boxplot()+
  facet_grid(rows=vars(input))


################################################################################################
## For a varying lambda

true_beta <- 2.36
## go through list
beta_2 <- vector(mode="list", length=n.iter)
for (i in 1:n.iter){
  for (j in 1:length(fix_sum[[1]])){
    beta_2[[i]][[j]] <- (ewp_vary[[i]][[j]][["coef_ewp"]][4]-true_beta)/ewp_vary[[i]][[j]][["coef_ewp"]][4]
    beta_2[[i]][[j]] <- as.data.frame(beta_2[[i]][[j]])
    beta_2[[i]][[j]] <- rename(beta_2[[i]][[j]],'beta2'='beta_2[[i]][[j]]')
    beta_2[[i]][[j]]$sum_limit <- ewp_vary[[i]][[j]][["sum_limit"]][1]
  }
}


beta_2_l_4.88 <- bind_rows(beta_2)

beta_2_l_4.88$set <- rep(1:n.iter,each=20)
beta_2_l_4.88$input <- 4.88


beta_2_lambda_full <- bind_rows(beta_2_l_4.88,beta_2_l_8,beta_2_l_10,beta_2_l_12,beta_2_l_15,beta_2_l_18)
beta_2_lambda_full$input <- as.factor(beta_2_lambda_full$input)
beta_2_lambda_full$set <- as.factor(beta_2_lambda_full$set)

beta_2_lambda_full$aftmax <- rep(1:20, length.out=length(beta_2_lambda_full$set))
beta_2_lambda_full$aftmax <- as.factor(beta_2_lambda_full$aftmax)


ggplot(beta_2_lambda_full, aes(x=aftmax, y=beta2))+
  geom_boxplot()+
  facet_grid(rows=vars(input))


##########################################################################################################


################################################################################################
## Try comparing with 100 sum limit fit

## Have to load both the true coef and sensitivity saves

## go through list
beta_2 <- vector(mode="list", length=n.iter)
for (i in 1:n.iter){
  for (j in 1:20){
    beta_2[[i]][[j]] <- (ewp_vary[[i]][[j]][["coef_ewp"]][4]-true_coef[[i]][4])/ewp_vary[[i]][[j]][["coef_ewp"]][4]
    beta_2[[i]][[j]] <- as.data.frame(beta_2[[i]][[j]])
    beta_2[[i]][[j]] <- rename(beta_2[[i]][[j]],'beta2'='beta_2[[i]][[j]]')
    beta_2[[i]][[j]]$sum_limit <- ewp_vary[[i]][[j]][["sum_limit"]][1]
  }
}


beta_2_l_4.88 <- bind_rows(beta_2)

beta_2_l_4.88$set <- rep(1:n.iter,each=20)
beta_2_l_4.88$input <- 4.88


beta_2_lambda_full_100 <- bind_rows(beta_2_l_4.88,beta_2_l_8,beta_2_l_10,beta_2_l_12,beta_2_l_15,beta_2_l_18)

beta_2_lambda_full_100$rel_limit <- beta_2_lambda_full_100$sum_limit/beta_2_lambda_full_100$input
beta_2_lambda_full_100$aftmax <- rep(1:20, length.out=length(beta_2_lambda_full_100$set))
beta_2_lambda_full_100$max_count <- beta_2_lambda_full_100$sum_limit - beta_2_lambda_full_100$aftmax
beta_2_lambda_full_100$rel_to_count <- beta_2_lambda_full_100$sum_limit/beta_2_lambda_full_100$max_count


beta_2_lambda_full_100$input <- as.factor(beta_2_lambda_full_100$input)
beta_2_lambda_full_100$set <- as.factor(beta_2_lambda_full_100$set)
beta_2_lambda_full_100$sum_limit <- as.factor(beta_2_lambda_full_100$sum_limit)
beta_2_lambda_full_100$rel_limit <- as.character(beta_2_lambda_full_100$rel_limit)
beta_2_lambda_full_100$aftmax <- as.factor(beta_2_lambda_full_100$aftmax)
beta_2_lambda_full_100$rel_to_count <- as.character(beta_2_lambda_full_100$rel_to_count)

save(beta_2_lambda_full_100, file="varylambda_sensitivity_output_100it.rda")

ggplot(beta_2_lambda_full_100, aes(x=sum_limit, y=beta2))+
  geom_boxplot()+
  facet_grid(rows=vars(input))+
  theme_bw()

beta_2_lambda_full_100$rel_to_count <- as.numeric(beta_2_lambda_full_100$rel_to_count)
beta_2_lambda_full_100$rel_limit <- as.numeric(beta_2_lambda_full_100$rel_limit)

ggplot(beta_2_lambda_full_100, aes(x=rel_to_count, y=beta2,color=input))+
  geom_point(alpha=0.35)+
  xlab(expression(frac("Summation limit","Maximum count")))+ ylab(expression(frac(beta[2]-hat(beta)[2],beta[2])))+
  scale_color_discrete(name= expression("Input"~lambda))+
  theme_bw(base_size = 15)

ggplot(beta_2_lambda_full_100, aes(x=rel_limit, y=beta2,color=input))+
  geom_point(alpha=0.3)+
  #scale_x_discrete(breaks=seq(1,4,1))+
  theme_bw()


#####################################################################################
## Same for varying beta2


## Plot the difference between the true value and the input value in the simulation
## load the appropriate output

### for a varying beta2

## go through list
beta_2 <- vector(mode="list", length=n.iter)
for (i in 1:n.iter){
  for (j in 1:20){
    beta_2[[i]][[j]] <- (ewp_vary[[i]][[j]][["coef_ewp"]][4]-true_coef[[i]][4])/ewp_vary[[i]][[j]][["coef_ewp"]][4]
    beta_2[[i]][[j]] <- as.data.frame(beta_2[[i]][[j]])
    beta_2[[i]][[j]] <- rename(beta_2[[i]][[j]],'beta2'='beta_2[[i]][[j]]')
    beta_2[[i]][[j]]$sum_limit <- ewp_vary[[i]][[j]][["sum_limit"]][1]
  }
}

beta_2_1 <- bind_rows(beta_2)

beta_2_1$set <- rep(1:n.iter,each=20)
beta_2_1$input <- 1


beta_2_full_100 <- bind_rows(beta_2_1,beta_2_1.5,beta_2_2.36,beta_2_4,beta_2_5)

beta_2_full_100$rel_limit <- beta_2_full_100$sum_limit/4.88
beta_2_full_100$aftmax <- rep(1:20, length.out=length(beta_2_full_100$set))
beta_2_full_100$max_count <- beta_2_full_100$sum_limit - beta_2_full_100$aftmax
beta_2_full_100$rel_to_count <- beta_2_full_100$sum_limit/beta_2_full_100$max_count


beta_2_full_100$input <- as.factor(beta_2_full_100$input)
beta_2_full_100$set <- as.factor(beta_2_full_100$set)
beta_2_full_100$sum_limit <- as.factor(beta_2_full_100$sum_limit)
beta_2_full_100$rel_limit <- as.character(beta_2_full_100$rel_limit)
beta_2_full_100$aftmax <- as.factor(beta_2_full_100$aftmax)
beta_2_full_100$rel_to_count <- as.character(beta_2_full_100$rel_to_count)


save(beta_2_full_100, file="varybeta2_sensitivity_output_100it.rda")

ggplot(beta_2_full_100, aes(x=sum_limit, y=beta2))+
  geom_boxplot()+
  facet_grid(rows=vars(input))

beta_2_full_100$rel_to_count <- as.numeric(beta_2_full_100$rel_to_count)
beta_2_full_100$rel_limit <- as.numeric(beta_2_full_100$rel_limit)

ggplot(beta_2_full_100, aes(x=rel_to_count, y=beta2,color=input))+
  geom_point(alpha=0.3)+
  xlab(expression(frac("Summation limit","Maximum count")))+ ylab(expression(frac(beta[2]-hat(beta)[2],beta[2])))+
  scale_color_discrete(name= expression("Input"~beta[2]))+
  theme_bw(base_size = 15)

ggplot(beta_2_full_100, aes(x=rel_limit, y=beta2,color=input))+
  geom_point(alpha=0.3)+
  xlab(expression(frac("Sum limit",lambda)))+ ylab(expression(frac(beta[2]-hat(beta[2]),beta[2])))+
  theme_bw()
