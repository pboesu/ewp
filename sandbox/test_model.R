system.time(fit <- ewp_reg(eggs ~ cov1 + cov2, data = linnet))
cbind(fit$coefficients, fit$se)
print(fit)
summary(fit)
print(summary(fit))


piedfly <- readRDS('../../2023_exponentially_weighted_poisson/data/piefl.rds') %>%
     select(max_num_eggs, min_first_egg, yearf,northing) %>%
     mutate(northing_scl = scale(northing)[,1])

m1<-glm(max_num_eggs~min_first_egg + yearf + northing_scl, family = poisson, data = piedfly)
summary(m1)

system.time(m2 <- ewp_reg(max_num_eggs~min_first_egg + yearf + northing_scl, data = piedfly))#very slow - converges within ~10 minutes, but hessian calculation does not wrap up for at least a further 40 minutes

