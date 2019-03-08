\name{hiv.mortmod}
\alias{hiv.mortmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HIV Prevalence-Calibrated Model Life Tables for Countries with Generalized HIV Epidemics
}
\description{
This function produces a complete set of mortality rates and optional life table given a region, sex, and set of inputs which must include HIV prevalence and one of the following mortality indicators: life expectancy at birth (e0), child mortality alone (5q0), or child mortality with adult mortality (45q15).  
}
\usage{
hiv.mortmod(e0=NULL, child.mort=NULL, adult.mort=NULL, prev=NULL, 
child.art=NULL, adult.art=NULL, model=1, region=1, sex=1, lt=FALSE, 
opt=TRUE, recal=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e0}{Life expectancy at birth (used only if model=1)
}
  \item{child.mort}{child mortality expressed as the probability a new born will die before reaching age 5, 5q0 (used only if model=2 or model=3)
}
  \item{adult.mort}{adult mortality expressed as the probability a 15-year-old will die before reaching age
60, 45q15 (used only if model=3)
}
  \item{prev}{HIV prevalence expressed as a percentage (i.e. six percent as \code{prev=6})
}
  \item{child.art}{Child ART coverage expressed a percentage (i.e. five percent as \code{child.art=5}); optional
}
  \item{adult.art}{Adult ART coverage expressed a percentage (i.e. 30 percent as \code{adult.art=30}); optional
}
  \item{model}{An integer to indicate which mortality indicator(s) are provided. \code{model=1} for life expectancy, \code{model=2} for child mortality alone, \code{model=3} for child mortality with adult mortality
}
  \item{region}{An integer to indicate which regional model to use. 1 for African and 0 for all other countries
}
  \item{sex}{An integer to indicate the sex of the desired output mortality rates. 1 for female, 0 for male}
  \item{lt}{Logical. If TRUE a life table is calculated based on the estimated mortality rates.
}
  \item{opt}{If \code{model=1} and \code{opt=TRUE} the model will vary the intercept such that the resulting estimated mortality rates  produce a life expectancy at birth that closely approximates the input life expectancy. If \code{model=2} and \code{opt=TRUE}, the model adjusts the earliest two probabilities of death (1q0, 4q1) in a proportional way so as to match the input 5q0. If \code{model=3} and \code{opt=TRUE}, the model will vary the intercept such that the resulting estimated mortality rates produce a 45q15 that closely approximates the input 45q15. The model will also adjust the earliest two mortality rates (1q0, 4q1) in a proportional way so as to match the input 5q0.
}
  \item{recal}{\code{recal=NULL} uses the model calibration (fixed parameters) saved in the workspace in the 'data' folder of this package. This argument could also be set to a character string indicating a different workspace (including file path if not in the working directory) from which to take the fixed parameters. The function \code{\link{hiv.mlt.recal}} will re-calibrate the fixed model parameters and save a workspace to the working directory that can be used for this argument.}
}
\value{
If \code{lt=FALSE}, the function produces a set of estimated age-specific mortality rates (nmx) produced for a given set of inputs.

If \code{lt=TRUE}, the function produces a life table with age intervals matching those of the nmx schedule on which the table is built and with columns for age, Lx (person years lived in the interval), lx (cumulative survival), qx (probability of death), ax (average number of years by those who died in the interval; set to half the length of the age interval except when under-5 mortality is split into 0-1 and 1-4), mx (mortality rate), Tx (person years life to live in the life table), and ex (life expectancy).
}
\details{
Further options allow for the inclusion of adult and child ART coverage as predictors of the weights along with prevalence and some demographic indicator.
}
\references{
Sharrow, David J., Samuel J. Clark, Adrian E. Raftery. 2014 \dQuote{Modeling Age-Specific Mortality for Countries with Generalized HIV Epidemics} PLoS ONE 9(5).
}
\author{
David J Sharrow
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mortmod.e0}}, \code{\link{mortmod.5q0}}, \code{\link{mortmod.45q15}}, \code{\link{hiv.mlt.recal}}
}
\examples{
## a life table at 1.5 percent prevalence, life expectancy of 60, 
## for Africa Females
hiv.mortmod(prev=1.5, e0=60, model=1, region=1, sex=1, lt=TRUE)

## a set of mortality rates at 2.5 percent prevalence, life expectancy of 53, 
## for Caribbean Males
hiv.mortmod(prev=2.5, e0=53, model=1, region=0, sex=0, lt=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ misc }% __ONLY ONE__ keyword per line
