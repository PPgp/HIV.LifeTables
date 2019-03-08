\name{mortmod.5q0}
\alias{mortmod.5q0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Age-specific mortality rate model as a function of HIV prevalence and child mortality (5q0)
}
\description{
This function produces a complete set of mortality rates given region, sex, HIV prevalence and child mortality (5q0, the probability a newborn will die before reaching age 5)}
\usage{
mortmod.5q0(child.mort, prev, child.art=NULL, adult.art=NULL, region=1, sex=1, 
  opt=TRUE, recal=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{child.mort}{child mortality expressed as the probability a new born will die before reaching age 5, 5q0}
  \item{prev}{HIV prevalence expressed as a percentage (i.e. six percent as \code{prev=6})}
  \item{child.art}{Child ART coverage expressed a percentage (i.e. five percent as \code{child.art=5}); optional}
  \item{adult.art}{Adult ART coverage expressed a percentage (i.e. 30 percent as \code{adult.art=30}); optional}
  \item{region}{An integer to indicate which regional model to use. 1 for African and 0 for all other countries}
  \item{sex}{An integer to indicate the sex of the desired output mortality rates. 1 for female, 0 for male}
  \item{opt}{If opt=TRUE, the model adjusts the earliest two probabilities of death (1q0, 4q1) in a proportional way so as to match the input 5q0.}
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
# with 5q0 of 0.06 and prevalence at 2.5 percent
# and no information on ART coverage
mortmod.5q0(child.mort=0.06, prev=2.5)

# produce a set of mortality rates for African males 
# with 5q0 of 0.15 and prevalence at 6 percent
# and 15 percent adult ART coverage
mortmod.5q0(child.mort=0.15, prev=6, adult.art=15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ misc }% __ONLY ONE__ keyword per line
