% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glcca.linear.R
\name{glcca.linear}
\alias{glcca.linear}
\title{Generalized Longitudial canonical correlation analysis}
\usage{
glcca.linear(
  X,
  varthresh = 0.95,
  projectthresh = 1,
  nboot = 1000,
  ncores = 12,
  seed.offset = 1234,
  verbose = FALSE
)
}
\arguments{
\item{projectthresh}{(default=1) threshold for the dimension reduction projection in lfpca.}

\item{verbose}{(default=FALSE) print all details}

\item{x}{object for input list(x, time, J, I, visit)}

\item{varthreshold}{(default=0.95) threshold to detemined the number of components of lpcs.}

\item{timeadjust}{(defult=FALSE)}

\item{method}{(defualt='Wilks') test statistic to be used. "Wilks","Hotelling", "Pillai", or "Roy".}
}
\value{
ccor

xcv_x0: Longitudinal Canonical vector for the intercept for x

xcv_x1: Longitudinal Canonical vector for the slope for x

xcv_y0: Longitudinal Canonical vector for the intercept for y

xcv_y1: Longitudinal Canonical vector for the slope for y
}
\description{
GLCCA handles more than three sets of data
Current version is implemented for only the linear trajectory.
}
\examples{
library(lcca)

}
