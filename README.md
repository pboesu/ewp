
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ewp

<!-- badges: start -->

[![R-CMD-check](https://github.com/pboesu/ewp/workflows/R-CMD-check/badge.svg)](https://github.com/pboesu/ewp/actions)
[![R-CMD-check](https://github.com/pboesu/ewp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pboesu/ewp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ewp is to provide a modelling interface for underdispersed
count data based on the exponentially weighted Poisson (EWP)
distribution described by Ridout & Besbeas (2004), allowing for
nest-level covariates on the location parameter
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
of the EWP. Currently only the three-parameter version of the
distribution (EWP_3) is implemented.

## Installation

You can install the development version of ewp directly from github like
so:

``` r
remotes::install_github('pboesu/ewp')
```

## Example

The package contains a reconstructed version of the linnet dataset from
Ridout & Besbeas (2004), which consists of 5414 clutch size records and
is augmented with two synthetic covariates, one that is random noise,
and one that is correlated with clutch size. The parameter estimates for
the intercept-only model presented in Ridout & Besbeas (2004) can be
reproduced like so:

``` r
library(ewp)
fit_null <- ewp_reg(eggs ~ 1, data = linnet)
#> start values are: 
#> (Intercept)       beta1       beta2 
#>    1.546823    0.000000    0.000000 
#> initial  value 9530.456809 
#> iter   2 value 7419.828863
#> iter   4 value 6892.193468
#> iter   6 value 5861.606330
#> iter   8 value 5372.176869
#> iter  10 value 5315.599815
#> iter  12 value 5305.705742
#> iter  14 value 5299.406589
#> iter  16 value 5299.362169
#> iter  16 value 5299.362148
#> iter  16 value 5299.362105
#> final  value 5299.362105 
#> converged
#> 
#> Calculating Hessian. This may take a while.
summary(fit_null)
#> Deviance residuals:
#> 
#> lambda coefficients (ewp3 with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.584657   0.003511   451.3   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> dispersion coefficients:
#>       Estimate Std. Error z value Pr(>|z|)    
#> beta1  1.46429    0.05588   26.20   <2e-16 ***
#> beta2  2.35675    0.05606   42.04   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Number of iterations in  optimization: 16 
#> Log-likelihood: -5299 on 3 Df
```

Note that the linear predictor on
![lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;lambda "lambda")
uses a log-link.

A model with nest-level covariates can be fitted by specifying a more
complex model formula - as in the base R `glm()`

``` r
fit <- ewp_reg(eggs ~ cov1 + cov2, data = linnet)
#> start values are: 
#>  (Intercept)         cov1         cov2        beta1        beta2 
#>  1.204279467 -0.001259307  0.071872488  0.000000000  0.000000000 
#> initial  value 9434.616462 
#> iter   2 value 6466.395926
#> iter   4 value 6387.476994
#> iter   6 value 5379.197851
#> iter   8 value 5353.658136
#> iter  10 value 4867.787688
#> iter  12 value 4649.639970
#> iter  14 value 4638.454116
#> iter  16 value 4514.116596
#> iter  18 value 4453.166131
#> iter  20 value 4425.890799
#> iter  22 value 4420.100521
#> iter  24 value 4420.096545
#> iter  26 value 4420.086482
#> iter  28 value 4420.080867
#> iter  30 value 4420.080707
#> final  value 4420.080600 
#> converged
#> 
#> Calculating Hessian. This may take a while.
summary(fit_null)
#> Deviance residuals:
#> 
#> lambda coefficients (ewp3 with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.584657   0.003511   451.3   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> dispersion coefficients:
#>       Estimate Std. Error z value Pr(>|z|)    
#> beta1  1.46429    0.05588   26.20   <2e-16 ***
#> beta2  2.35675    0.05606   42.04   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Number of iterations in  optimization: 16 
#> Log-likelihood: -5299 on 3 Df
```

:warning: **Note that the maximum likelihood optimisation procedure is
still experimental**

In particular:

-   Fitting is very slow (think minutes, not seconds), especially when
    estimating the Hessian matrix for more than a couple of parameters
    (think hours, not minutes)!
-   Estimates may not be stable results for models with many covariates
    and/or very large sample sizes (1000s)

:warning:
