\name{HIV.LifeTables-package}
\alias{HIV.LifeTables-package}
\alias{HIV.LifeTables}
\docType{package}
\title{HIV Prevalence-Calibrated Model Life Tables for Countries with Generalized HIV epidemics}
\description{The functions in this package produce a complete set of mortality rates as a function of a combination of HIV prevalence and either life expectancy at birth (e0), child mortality alone (5q0), child mortality with adult mortality (45q15), or a user-defined nqx entry parameter.}
\details{
\tabular{ll}{
Package: \tab HIV.LifeTables\cr
Type: \tab Package\cr
Version: \tab UN\cr
Date: \tab 2015-07-05\cr
License: \tab GPL-2\cr
}
The central function is \code{\link{hiv.mortmod}} which takes prevalence and some mortality indicator (e0, 5q0 alone, or 5q0 with 45q15) to produce set of mortality rates. The user also needs to select a region (Africa or all other countries) and sex. \code{\link{hiv.mortmod}} will also produce a life table based on the estimated mortality rates. \code{\link{mortmod.e0}}, \code{\link{mortmod.5q0}}, and \code{\link{mortmod.45q15}} produce mortality rates for the various input combinations.

Each of the aforementioned functions have additional options to include child and adult ART coverage.

By default the functions listed above access fixed model parameters necessary to produce predicted mortality rates from a workspace in the 'data' folder of this package, but the user can re-calibrate the model with appropriate demographic and epidemiological data (and optionally with another user-defined entry mortality parameter) with the function \code{\link{hiv.mlt.recal}}. That function can produce and save a workspace with re-calibrated model parameters that can be input to the functions listed above. 

}
\author{
David J Sharrow \cr
}
\references{
Sharrow, David J., Samuel J. Clark, Adrian E. Raftery. 2014 \dQuote{Modeling Age-Specific Mortality for Countries with Generalized HIV Epidemics} PLoS ONE 9(5).
}

\keyword{ package }
%\seealso{

%}
\examples{
## a life table at 1.5 percent prevalence, life expectancy of 60, 
## for Africa Females
hiv.mortmod(prev=1.5, e0=60, model=1, region=1, sex=1, lt=TRUE)

## a set of mortality rates at 2.5 percent prevalence, life expectancy of 
## 53, for non-African Males
hiv.mortmod(prev=2.5, e0=53, model=1, region=0, sex=0, lt=FALSE)
}
