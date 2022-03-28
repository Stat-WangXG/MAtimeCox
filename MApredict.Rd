\name{MApredict}
\alias{MApredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Model Averaging Prediction for lifetime data with time-dependent effects.
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Forecast conditional survival probabilities of new subjects by model averaging approach.
}

\usage{
MApredict(MAtdcCOX.object, newdata, t_star)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAtdcCOX.object}{
object of class MAtdcCOX, see \code{\link{MAtdcCOX}}.
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
data to be predicted, which can be a vector, matrix or data.frame.
}
  \item{t_star}{
%%     ~~Describe \code{t_star} here~~
a vector which includes the time points corresponding to the survival probabilities one wants to predict.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The newdata contains only covariate vectors. If there are missing values in the newdata, the corresponding prompt information will be output.
}

\value{
%%  ~Describe the value returned
\item{n_test}{sample size of newdata}
\item{t_star.order}{the t_star arranged in increasing order}
\item{S_MA_t_star}{survival probabilities of subjects in newdata at time point t_star.order}
\item{time.order}{the follow-up time arranged in increasing order in MAtdcCOX.object}
\item{delta.order}{the status indicator vector in increasing order of follow-up time in MAtdcCOX.object}
\item{S_MA_test}{survival probabilities of subjects in newdata at time.order}
}

\note{
%%  ~~further notes~~
The dimension of the covariate vectors in the newdata should be the same as the one in the MAtdcCOX.object, which is p.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MAtdcCOX}}
}

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

# Estimate the candidate models and their weights based on the train data set

submodels=MAtdcCOX(data,c(5:8),test.plot=F,compute.S=T)

# Predict the conditional survival functions at time points "{1,2,3,4,5}" given the covariate matrix

results_pre=MApredict(submodels,submodels$z.order[1:3,],c(1:5))
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
