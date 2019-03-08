\name{mortmod.ep}
\alias{mortmod.ep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Age-specific mortality rate model as a function of HIV prevalence and a user-supplied entry parameter (nqx)
}
\description{
This function produces a complete set of mortality rates given region, sex, HIV prevalence and a user-supplied entry parameter}
\usage{
mortmod.ep(recal, entry.par, prev, child.art=NULL, adult.art=NULL, 
  region=1, sex=1, opt=TRUE, entry.nax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{recal}{A character string indicating a different workspace (including file path if not in the working directory) from which to take the fixed parameters. The function \code{hiv.mlt.recal} will re-calibrate the fixed model parameters and save a workspace to the working directory that can be used for this argument.}
  \item{entry.par}{A user-supplied entry parameter (must be some nqx derived from the age groups in the calibration data)}
  \item{prev}{HIV prevalence expressed as a percentage (e.g six percent as \code{prev=6})}
  \item{child.art}{Child ART coverage expressed a percentage (e.g five percent as \code{child.art=5}); optional}
  \item{adult.art}{Adult ART coverage expressed a percentage (e.g. 30 percent as \code{adult.art=30}); optional}
  \item{region}{An integer to indicate which regional model to use. 1 for African and 0 for all other countries}
  \item{sex}{An integer to indicate the sex of the desired output mortality rates. 1 for female, 0 for male}
  \item{opt}{If opt=TRUE, the model varies the model intercept so as to match the user-supplied entry parameter.}
  \item{entry.nax}{A numeric vector of user-defined nax values of the length of \code{ages} defined in \code{hiv.mlt.recal}. If \code{opt=TRUE}, the model calculates a life table in order to match the user-defined entry parameter. The default nax values (\code{entry.nax=NULL}) are half the length of the age interval unless mortality under age 5 is split to 0-1 and 1-4. The user can supply their own nax values for this process with this argument.} 

}
\details{
Further options allow for the inclusion of adult and child ART coverage as predictors of the weights along with prevalence and the user-supplied entry parameter.
}
\value{A set of estimated age-specific mortality rates (nmx) produced for a given set of model inputs}
\references{
Sharrow, David J., Samuel J. Clark, Adrian E. Raftery. 2014 \dQuote{Modeling Age-Specific Mortality for Countries with Generalized HIV Epidemics} PLoS ONE 9(5).
}
\author{David J Sharrow}
\note{
The user-supplied entry parameter must some probability of death over age, e.g. 35q15.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hiv.mlt.recal}}
}
\examples{
\dontrun{
# produce a set of mortality rates for African females 
# with 35q15 of 0.40 and prevalence at 2.5 percent
# and no information on ART coverage
mortmod.ep(recal="SomeRecal.RData", entry.par=0.40, prev=2.5)

# produce a set of mortality rates for African males 
# with 35q15 of 0.15 and prevalence at 6 percent
# and 15 percent adult ART coverage
mortmod.ep(recal="SomeRecal.RData", entry.par=0.15, prev=6, adult.art=15)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ misc }% __ONLY ONE__ keyword per line
