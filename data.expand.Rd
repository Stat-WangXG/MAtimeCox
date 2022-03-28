\name{data.expand}
\alias{data.expand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Expand the original data set
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Expand the original data set in a repeated measurement format: the counting-process style.
}
\usage{
data.expand(delta2, time2, z2, bs7_2, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta2}{
status indicators in increasing order of follow-up time
}
  \item{time2}{
follow-up time in increasing order
}
  \item{z2}{
covariate matrix in increasing order of follow-up time
}
  \item{bs7_2}{
B-spline basis matrix for the time-dependent coefficients
}
  \item{K}{
number of B-spline basis functions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The data for each subject is presented as multiple rows or "observations", each of which applies to an interval of observation (start, stop]. A new variable will be defined to identify the records for each subject after expanding the original data set.
}
\value{
%%  ~Describe the value returned
\item{z_new }{the covariate matrix after expanding data}
\item{bs_new }{the B-spline basis matrix after expanding data}
\item{delta_new }{the status indicator vector after expanding data}
\item{time_start }{the left endpoint of each interval of observation}
\item{time_stop }{the right endpoint of each interval of observation}
\item{index }{a variable to identify the records for each subject}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Martinussen and Scheike, Dynamic Regression Models for Survival Data, Springer (2006).
}

\note{
%%  ~~further notes~~
It is assumed that there are no tied events, and if such are present one can  break the ties by adding random noise.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rm(list = ls())
library(timereg)
library(splines)

# Extract a train data set

data(sTRACE)
sTRACE=na.omit(sTRACE)
data=sTRACE[sTRACE$status==0 | sTRACE$status==9,c(2,4,5,6,7,9,8,3)]
data$status=data$status/9
n=nrow(data)
p=ncol(data)-2
K=9
z=as.matrix(data[,1:p])
time=data[,ncol(data)-1]
delta=data[,ncol(data)]
delta=delta[order(time)]
z=z[order(time),]
time=time[order(time)]
Bspline=bs(time,df=K,knots =quantile(time[delta==1],probs = seq(1:(K-4))/(K-3)),intercept=T,degree = 3)

# Expand the original data set in the counting-process style

data_new=data.expand(delta, time, z, Bspline,K)
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
