\name{w.t.i.G}
\alias{w.t.i.G}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute the Inverse of Probability of Censoring Weighting
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the Inverse of Probability of Censoring Weighting for a subject at some time point.
}
\usage{
w.t.i.G(t, i, time, delta, G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
the specified time point in the Inverse of Probability of Censoring Weighting
}
  \item{i}{
the index of subject
}
  \item{time}{
the follow-up time in increasing order
}
  \item{delta}{
status indicators in increasing order of follow-up time
}
  \item{G}{
a matrix which each column denotes the estimated survival function of a subject for censoring time.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The survival function of censoring time can be estimated by different methods and the Kaplan-Meier estimator is used in our model averaging approach.
}

\value{
%%  ~Describe the value returned
\item{weight.t.i}{Inverse of Probability of Censoring Weighting at time point t}
}

\references{
%% ~put references to the literature/web site here ~
Gerds, T. A. and Schumacher, M. (2006), Consistent estimation of the expected brier score in general survival models with right-censored event times. Biometrical Journal, 48, 1029–1040.
}

\note{
%%  ~~further notes~~
The dimension of "time" need be the same as the order of square matrix "G". The status indicator is 1 for the dead and 0 for the alive. It should be noted that the survival functions in "G" correspond to censoring time rather than survival time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rm(list = ls())
library(survival)

# Create a simple data set

set.seed(1)
n=100
time=rexp(n,0.5)
delta=rbinom(n,1,0.4)
delta=delta[order(time)]
time=time[order(time)]
data=data.frame(time,delta)
model_c=survfit(Surv(time,1-delta)~1,data = data)

# Estimate the survival function of censoring time

G=matrix(rep(model_c$surv,n),byrow = F,ncol = n)

# Compute the Inverse of Probability of Censoring Weighting for the first subject at time point "time[1]"

weight_temp=w.t.i.G(time[1],1,time,delta,G)
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
