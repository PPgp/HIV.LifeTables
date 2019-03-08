\name{mortmod.e0}
\alias{mortmod.e0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Age-specific mortality rate model as a function of HIV prevalence and life expectancy at birth (e0)}
\description{
This function produces a complete set of mortality rates given a region, sex, HIV prevalence and life expectancy at birth (e0)}
\usage{
mortmod.e0(e0, prev, child.art=NULL, adult.art=NULL, region=1, 
  sex=1, opt=TRUE, recal=NULL, determ=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e0}{Life expectancy at birth}
  \item{prev}{HIV prevalence expressed as a percentage (i.e. six percent as \code{prev=6})}
  \item{child.art}{Child ART coverage expressed a percentage (i.e. five percent as \code{child.art=5}); optional}
  \item{adult.art}{Adult ART coverage expressed a percentage (i.e. 30 percent as \code{adult.art=30}); optional}
  \item{region}{An integer to indicate which regional model to use. 1 for African and 0 for all other countries}
  \item{sex}{An integer to indicate the sex of the desired output mortality rates. 1 for female, 0 for male}
  \item{opt}{\code{opt=TRUE} will vary the intercept in the model such that the resulting estimated mortality rates produce a life expectancy at birth that matches the input life expectancy.}
  \item{determ}{\code{determ=FALSE} will add a random deviation to each age-specific mortality rate. The deviations are drawn from a normal distribution with mean=0. The standard deviation is age-specific and is the standard deviation from the distribution of residuals at each each age after fitting the calibration data with the deterministic model.}
  \item{recal}{\code{recal=NULL} uses the model calibration (fixed parameters) saved in the workspace in the 'data' folder of this package. This argument could also be set to a character string indicating a different workspace (including file path if not in the working directory) from which to take the fixed parameters. The function \code{hiv.mlt.recal} will re-calibrate the fixed model parameters and save a workspace to the working directory that can be used for this argument.}
}
\details{
Further options allow for the inclusion of adult and child ART coverage as predictors of the weights along with prevalence and some demographic indicator.
}
\value{A set of estimated age-specific mortality rates (nmx) produced for a given set of model inputs}
\references{
Sharrow, David J., Samuel J. Clark, Adrian E. Raftery. 2014 \dQuote{Modeling Age-Specific Mortality for Countries with Generalized HIV Epidemics} PLoS ONE 9(5).
}
\author{David J Sharrow}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hiv.mortmod}}, \code{\link{hiv.mlt.recal}}
}
\examples{
# produce a set of mortality rates for African females 
# with life expectancy at birth of 55 and prevalence at 2.5 percent
# and no information on ART coverage
mortmod.e0(e0=55, prev=2.5)

# produce a set of mortality rates for African males 
# with life expectancy at birth of 51 and prevalence at 12 percent
# and 15 percent adult ART coverage
mortmod.e0(e0=51, prev=12, adult.art=15, sex=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ misc }% __ONLY ONE__ keyword per line