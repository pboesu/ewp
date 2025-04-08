## Release summary
This is a re-submission to CRAN following feedback on the first submission.

We have addressed all feedback:

> Please do not start the description with "This package", package name,title or similar.

The package description has been updated.

> Please add \value to .Rd files regarding exported methods

\value Rd tags have been added to print.ewp and summary.ewp documentation

> Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN policies. ->  R/ewp_reg.R
> Please do not set a seed to a specific number within a function. -> R/ewp_reg.R

The simulate.ewp method has been updated to not modify .GlobalEnv, and to not require an input seed.


## Test environments
* Windows Server 2022 x64 (build 20348) R version 4.4.1 (2024-06-14 ucrt), R-devel (2024-06-24 r86823 ucrt)
* macOS Sonoma 14.5 (on github-actions) R version 4.4.1 (2024-06-14)
* Ubuntu 22.04.4 LTS (on github-actions) R version 4.3.3 (2024-02-29), R version 4.4.1 (2024-06-14), and R-devel (2024-06-24 r86823)


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* R CMD Check highlights the following possibly misspelled words in DESCRIPTION:
  Besbeas (15:200)
  Ridout (15:191)
  underdispersed (3:31, 15:73)
  These are all spelled correctly.
