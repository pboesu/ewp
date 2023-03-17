test_that("null model runthrough", {
  fit_null <- ewp_reg(eggs ~ 1, data = linnet)
  expect_s3_class(fit_null, "ewp")
  print(fit_null)
  summary(fit_null)
  #get estimates and compare against published estimates from Ridout & Besbeas (2004)
  coefs <- unname(coef(fit_null))
  expect_equal(exp(coefs[1]),4.88,tolerance=0.017)
  expect_equal(coefs[2],1.46,tolerance=0.056)
  expect_equal(coefs[3],2.36,tolerance=0.056)
  })

test_that("NA handling and dropped levels", {
  set.seed(1234)
  linnet2 <- linnet
  linnet2$eggs[sample(nrow(linnet2), size = 100)] <- NA
  linnet2$arbfac <- factor(sample(letters[1:10], size = nrow(linnet2), replace = TRUE))
  linnet2 <- subset(linnet2, arbfac %in% letters[1:5], drop = FALSE)
  fit_na <- ewp_reg(eggs ~ arbfac, data = linnet2)
  expect_s3_class(fit_na, "ewp")
})
