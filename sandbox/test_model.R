system.time(fit_null <- ewp_reg(eggs ~ 1, data = linnet))
system.time(fit <- ewp_reg(eggs ~ cov1 + cov2, data = linnet))
cbind(fit$coefficients, fit$se)
print(fit)
summary(fit)
print(summary(fit))

system.time(fitNM <- ewp_reg(eggs ~ cov1 + cov2, data = linnet, method='Nelder-Mead'))

piedfly <- readRDS('../../2023_exponentially_weighted_poisson/data/piefl.rds') %>%
     select(max_num_eggs, min_first_egg, yearf,northing) %>%
     mutate(northing_scl = scale(northing)[,1])

m1<-glm(max_num_eggs~min_first_egg + yearf + northing_scl, family = poisson, data = piedfly)
summary(m1)

piedfly_s <- piedfly %>%
  group_by(yearf) %>%
  slice_sample(n = 200)

system.time(m2 <- ewp_reg(max_num_eggs~min_first_egg + yearf + northing_scl, data = piedfly_s))#very slow - converges within ~10 minutes, but hessian calculation does not wrap up for at least a further 60 minutes - aborted at 5000 secs
system.time(m3 <- ewp_reg(max_num_eggs~min_first_egg + northing_scl, data = piedfly_s))#hessian slowdown (unsurprisingly) caused by year factors - optimisation may not be stable for large sample sizes - maybe to do with parscale - maybe try Nelder-Mead instead?
