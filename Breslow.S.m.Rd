\name{Breslow.S.m}
\alias{Breslow.S.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimator of conditional survival function for the mth candidate model.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimate the conditional survival function for the mth candidate model by the Breslow type method.
}
\usage{
Breslow.S.m(n, delta, z, B_spline, m, beta_m, theta_m, covariate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{delta}{
status indicator vector
}
  \item{z}{
covariate matrix of train data
}
  \item{B_spline}{
B-spline basis matrix for the time-dependent coefficients
}
  \item{m}{
the index of candidate model
}
  \item{beta_m}{
estimator of time-independent coefficient vector in the mth candidate model
}
  \item{theta_m}{
estimator of B-spline coefficient vector in the mth candidate model
}
  \item{covariate}{
a new covariate vector corresponding to conditional survival function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The status indicator vector and the covariate matrix of dimension n*p need be arranged in increasing order of follow-up times.
}
\value{
%%  ~Describe the value returned
\item{Surv.m.covariate }{the Breslow estimator of the conditional survival function for the mth candidate model}
}

\note{
%%  ~~further notes~~
The dimension of covariate vector "covariate" should be the same as the number of columns of covariate matrix "z". When "Inf" occurs, it will be replaced with a finite value large enough for subsequent calculations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rm(list = ls())
library(survival)
library(timereg)
library(splines)
library(quadprog)

# Extract a train data set

data(sTRACE)
sTRACE=na.omit(sTRACE)
data=sTRACE[sTRACE$status==0 | sTRACE$status==9,c(2,4,5,6,7,9,8,3)]
data$status=data$status/9

# Estimate the candidate models and their weights

submodels=MAtdcCOX(data,c(5:8),test.plot=F,compute.S=T)
n=submodels$n
delta=submodels$delta.order
z=submodels$z.order
Bspline=submodels$Bspline
submodel_all=submodels$candidate_models
p=submodels$p
K=submodels$K_n

# For the first candidate model, estimate the conditional survival function given covariate vector "z[1,]" by the Breslow type method

Breslow_S_temp=Breslow.S.m(n,delta,z,Bspline,1,submodel_all[[1]][["coefficients"]][1:(p-1)],submodel_all[[1]][["coefficients"]][p:(p+K-1)],z[1,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
