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
