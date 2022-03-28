\name{MAtdcCOX}
\alias{MAtdcCOX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Model averaging estimation for time-dependent coefficient COX model.
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fit a final model based on the model averaging for time-dependent coefficient COX model. Each candidate model is selected as a semi-parametric time-varying coefficient Cox model. For the mth candidate model, the coefficient of the mth covariate is considered time-varying and others are time-invariant. The optimal weight is solved by quadratic programming technique.
}

\usage{
MAtdcCOX(data,K.set=c(5:10),test.plot=F,compute.S=T)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
a data.frame which includes covariates, follow up time and status indicators. Here, the data is right-censored, with a status indicator of 1 indicating death, and 0 otherwise.
}
\item{K.set}{
a vector containing the range of the number of basis functions, with a default value of 5 to 10.
}
\item{test.plot}{
a logical variable (defaults to False) to determine whether to plot the coefficient estimators based on scaled Schoenfeld residuals.
}
\item{compute.S}{
a logical variable (defaults to True) to determine whether to compute the estimators of survival functions for subjects in the data.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The last two columns of the data.frame are the follow-up time and status indicators, respectively. When there are tied events, we deal with them by adding random noise. Without loss of generality, the cubic B-splines are employed to approximate the time-dependent coefficients. The data after preprocessing is expanded in a repeated measurement format and the common counting-process style is used, see \code{\link{data.expand}}. We determine the best number of basis functions by minimizing the sum of the AICs of candidate models. If the concern is the prediction of survival probabilities of new subjects, we recommend setting compute.S to False when the computation burden is heavy.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{data_train}{the data after preprocessing}
\item{z.order}{covariate matrix in increasing order of follow-up time}
\item{Bspline}{B-spline basis matrix for the time-dependent coefficients}
\item{time.order}{the follow-up time in increasing order}
\item{delta.order}{the status indicator vector in increasing order of follow-up time}
\item{n}{sample size}
\item{p}{dimension of covariate vector}
\item{K_n}{the number of B-spline basis functions determined by AIC}
\item{test}{results of time-dependent effects tests based on scaled Schoenfeld residuals}
\item{candidate_models}{candidate models in model averaging estimation}
\item{MA_weights}{the estimated weight vector of candidate models}
\item{S_MA_train}{model averaging estimators of conditional survival functions for subjects in data_train}
}

\references{
%% ~put references to the literature/web site here ~
P. Grambsch and T. Therneau (1994), Proportional hazards tests and diagnostics based on weighted residuals. Biometrika, 81, 515-26.
}


\note{
%%  ~~further notes~~
The data containing missing values will be removed and the tied events will be broken by adding random noise. Therefore, the data returned is not necessarily same as the input data. It should be noted that the minimum number of B-spline basis functions can not be less than 5 since cubic B-splines are used by default.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{data.expand},\link{w.t.i.G},\link{Breslow.S.m}}
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

# Estimate the candidate models and their weights

submodels=MAtdcCOX(data,c(5:8),test.plot=F,compute.S=T)
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
