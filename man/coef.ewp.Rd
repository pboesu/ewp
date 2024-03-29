% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ewp_reg.R
\name{coef.ewp}
\alias{coef.ewp}
\title{Extract coefficients}
\usage{
\method{coef}{ewp}(object, ...)
}
\arguments{
\item{object}{an object of class ewp}

\item{...}{ignored}
}
\value{
a vector of coefficient values. Beware that the lambda parameters are on the log-link scale, whereas the betas are estimated using an identity link.
}
\description{
Extract coefficients
}
