## code to prepare example datasets

#recreate Linnet data from Ridout & Besbeas 2004
linnet = data.frame(eggs = c(rep(1,18),
                             rep(2,35),
                             rep(3,210),
                             rep(4,1355),
                             rep(5,3492),
                             rep(6,299),
                             rep(7,5)
                             )
                    )
#add one noise only, one signal+noise covariate
linnet$cov1 <- rnorm(n = nrow(linnet))
linnet$cov2 <- linnet$eggs + rnorm(n = nrow(linnet))


usethis::use_data(linnet, overwrite = TRUE)

# #read in pied flycatcher data
# library(dplyr)
# piedfly <- readRDS('data-raw/piefl.rds') %>%
#   select(max_num_eggs, min_first_egg, yearf,northing) %>%
#   mutate(northing_scl = scale(northing)[,1])
#
#
# summary(piedfly)
#
# m1<-glm(max_num_eggs~min_first_egg + yearf + northing_scl, family = poisson, data = piedfly)
# summary(m1)
#
# usethis::use_data(piedfly, overwrite = TRUE)
