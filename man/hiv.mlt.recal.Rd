\name{hiv.mlt.recal}
\alias{hiv.mlt.recal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HIV Life Tables fixed parameter re-calibration
}
\description{
This function takes sex-specific mortality and epidemiological data of matching periods and re-calibrates (estimates a new set of fixed parameters) for the HIV MLT model
}
\usage{
hiv.mlt.recal(mx.m, mx.f, prev, c.art, a.art, ages=c(0,1,seq(5,100,5)), 
  la.countries=c("Bahamas", "Belize", "Cambodia", "Dominican Republic", "Estonia", 
  "Guyana", "Haiti", "Honduras", "Jamaica", "Panama", "Russian Federation", "Suriname", 
  "Thailand", "Trinidad and Tobago", "Ukraine"), midperiod=seq(1973,2008,5), determ=TRUE,
  user.par=NULL, nax=NULL, save.output="Ws-models.Rdata")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mx.m}{A tab delimited text file of age-specific mortality rates for males readable with \code{read.table()}. The file must have columns titles: "country" (a string for country name); "country code" - a numeric for the UN country code; "age" - a numeric indicating age groups; then column names for each period in the calibration data set. Five year intervals should be "1975-1980", "1980-1985", etc. Single year intervals should be "1980", "1981", etc.
}
  \item{mx.f}{A tab delimited text file of age-specific mortality rates for females readable with \code{read.table()}. The file must have columns titles: "country" (a string for country name); "country code" - a numeric for the UN country code; "age" - a numeric indicating age groups; then column names for each period in the calibration data set. Five year intervals should be "1975-1980", "1980-1985", etc. Single year intervals should be "1980", "1981", etc.
}
  \item{prev}{A tab delimited text file of HIV prevalence rates readable with \code{read.table()}. The file must have column titles: "country" - a string for country name; "country code" - a numeric for the UN country code; then column names for each period in the calibration data set. Even if \code{mx.f} and \code{mx.m} are in 5-year intervals, this data can be 1-year and the argument midperiod will indicate which column should be associated with which 5-year period of mx data. If mx data is 1-year and prevalence data is as well, set \code{midperiod} to the 1-year intervals.
}
  \item{c.art}{A tab delimited text file of child ART coverage rates readable with \code{read.table()}. The file must have column titles: "country" - a string for country name; "country code" - a numeric for the UN country code; then column names for each period in the calibration data set. Even if \code{mx.f} and \code{mx.m} are in 5-year intervals, this data can be 1-year and the argument midperiod will indicate which column should be associated with which 5-year period of mx data. If mx data is 1-year and prevalence data is as well, set \code{midperiod} to the 1-year intervals. 
}
  \item{a.art}{A tab delimited text file of adult ART coverage rates readable with \code{read.table()}. The file must have column titles: "country" - a string for country name; "country code" - a numeric for the UN country code; then column names for each period in the calibration data set. Even if \code{mx.f} and \code{mx.m} are in 5-year intervals, this data can be 1-year and the argument midperiod will indicate which column should be associated with which 5-year period of mx data. If mx data is 1-year and prevalence data is as well, set \code{midperiod} to the 1-year intervals. 
}
  \item{ages}{A numeric vector indicating the beginning age of age intervals of the mx data (i.e. 5-year age intervals with childhood ages broken into two groups would be \code{ages=c(0,1,seq(5,100,5))}). Age intervals must match 'ages' column from mx input data.
}
  \item{la.countries}{A character vector of non-African countries in the calibration data set. Separate models are calibrated for countries in Africa and countries not in Africa. This argument identifies those countries during calibration.
}
  \item{midperiod}{A numeric vector to indicate which columns of the prevalence, adult ART, and child ART inputs should be used for calibration. If mortality data are in greater than 1-year intervals (e.g. 5-year intervals) but epidemiological data are in 1-year intervals, \code{midperiod} tells the function which periods in the epidemiological data should correspond to the multi-year periods in mortality data. If mx data are 1-year and epidemiological data are also 1-year, set \code{midperiod} to the 1-year periods from the calibration data sets.
}
  \item{determ}{\code{determ=FALSE} will calculate the standard deviation of the distribution of residuals at each age after fitting the input data with the calibrated deterministic model. These standard deviations can be used to add random noise to the mortality rate schedule by simulating a random deviation for each age from a normal distribution with mean=0 and the standard deviation for a specific age to be added (or subtracted in the case of a negative deviation) to the predicted mortality rate at a given age. \code{determ=TRUE} will set both the mean and standard deviation to 0, so no deviation will be added--the result being just the deterministic model output. This functionality is only available for the model predicting mortality rates as a function of life expectancy and HIV prevalence, \code{\link{mortmod.e0}}.}
    \item{user.par}{This argument provides for a user defined mortality input parameter that is not one of the three default model inputs. This user defined input must be some nqx that can be derived from the calibration data age groups. Must be entered as a numeric vector of length 2 where the first element indicates the beginning age of the nqx and the second argument is the end of the interval (e.g. if one wants to calibrate a model with 35q15 enter \code{user.par=c(15,50)}).
  }
  \item{nax}{User-defined nax values. When calculating life tables from which the entry mortality parameters (e0, 5q0, 45q15 or a user-defined entry parameter) are derived, the default nax values are half the length of the age interval unless mortality under age 5 is split to 0-1 and 1-5. The user can provide alternative nax values as a matrix with rows for each country-period and columns for each nax. If provided, the matrix must include a set of nax values for each age and each country-period.} 
  \item{save.output}{A character string giving the name of the output workspace file. If \code{NULL}, no output workspace is saved.
}
}
\details{
If \code{user.par} is invoked, the function will calculate a life table for each input mortality rate schedule and derive the user defined entry parameter from those life tables based on the the ages indicated in \code{user.par}. The arguments in \code{user.par} must indicate the age bounds of the nqx interval and must be elements in \code{ages}. 
}
\value{
\item{svd.comps.m}{The first three left-singular vectors (LSV) from a Singular Value Decomposition (SVD) of the matrix of male mortality rates.}
\item{svd.comps.f}{The first three left-singular vectors (LSV) from a Singular Value Decomposition (SVD) of the matrix of female mortality rates.}
\item{weights.m}{A matrix of OLS regression coefficients from regressing the male mortality rates schedules on the LSV from the male SVD.}
\item{weights.f}{A matrix of OLS regression coefficients from regressing the female mortality rates schedules on the LSV from the female SVD.}
\item{params.a.m}{An array containing the fixed model parameters for African males. Rows show the included predictor variables with columns for the various models. Each matrix in the array is for predicting the weight for one of the first three LSV.}
\item{params.a.f}{An array containing the fixed model parameters for African females. Rows show the included predictor variables with columns for the various models. Each matrix in the array is for predicting the weight for one of the first three LSV.}
\item{params.na.m}{An array containing the fixed model parameters for non-African males. Rows show the included predictor variables with columns for the various models. Each matrix in the array is for predicting the weight for one of the first three LSV.}
\item{params.na.f}{An array containing the fixed model parameters for non-African females. Rows show the included predictor variables with columns for the various models. Each matrix in the array is for predicting the weight for one of the first three LSV.}
\item{input.epi}{A data frame containing the input epidemiological information for calibrating the model. This data frame is included as a check to be sure the function is accessing the correct periods for epidemiological data.}
\item{input.mort.m}{A data frame containing the male input mortality information for calibrating the male models. These mortality indicators are derived from life tables calculated from the input mortality rates and this data frame is included as a check to be sure the function is calculating plausible numbers.}
\item{input.mort.f}{A data frame containing the female input mortality information for calibrating the female models. These mortality indicators are derived from life tables calculated from the input mortality rates and this data frame is included as a check to be sure the function is calculating plausible numbers.}
}
\references{
Sharrow, David J., Samuel J. Clark, Adrian E. Raftery. 2014 \dQuote{Modeling Age-Specific Mortality for Countries with Generalized HIV Epidemics} PLoS ONE 9(5).

Sharrow, David J., Yanjun He, Samuel J. Clark, Adrian E. Raftery. 2015. \dQuote{Probabilistic Projections of Mortality in Countries with Generalized HIV Epidemics for Use in Total Population Projection} paper presented at the annual meeting of the Population Association of America, San Diego, CA.\url{http://paa2015.princeton.edu/abstracts/152851}
}
\author{
David J Sharrow
}
\note{
The saved output workspace from this function can be used in the \code{recal} argument in \code{\link{hiv.mortmod}}, \code{\link{mortmod.e0}}, \code{\link{mortmod.5q0}}, \code{\link{mortmod.45q15}}, or \code{\link{mortmod.ep}}.

Setting \code{determ=FALSE} can substantially increase the re-calibration time depending on the number of input mortality rate schedules because all of the input mortality data is fit with the calibrated model before the residuals (and thus the standard deviation of the distribution of residuals at each age) can be calculated. When simulating many future trajectories of mortality rates, adding this random noise will widen subsequent predictive intervals, but if the user is only interested in the deterministic model output, it is advised to set \code{determ=TRUE}, the default.
}

\examples{
\dontrun{
# ## WPP 2012 and UNAIDS 2012
#1. ages
ages.wpp <- c(0,1,seq(5,100,5)) ## default input

#2. prevalence (from Patrick Gerland files sent 2015-03-19)
prev.2012 <- read.table("HIV_PrevalenceRate.txt", header=T)

#3a. ART (adult)
aart.2012 <- read.table("ART_CoverageAdult.txt", header=T)

#3b. ART (child)
cart.2012 <- read.table("ART_CoverageChild.txt", header=T)

#4. nmx
mxF.wpp <- read.table("mxF.txt", header=T)
mxM.wpp <- read.table("mxM.txt", header=T)

#5. non-African countries
# create vectors to indicate which rows are the Latin American (non-African) 
# and African Countries
la.countries=c("Bahamas", "Belize", "Cambodia", "Dominican Republic", 
"Estonia", "Guyana", "Haiti", "Honduras", "Jamaica", "Panama", 
"Russian Federation", "Suriname", "Thailand", "Trinidad and Tobago", 
"Ukraine")

#6. midperiod
midperiod.wpp=seq(1973,2008,5)

hiv.mlt.recal(mx.m=mxM.wpp, mx.f=mxF.wpp, c.art=cart.2012, a.art=aart.2012, 
  la.countries=la.countries, midperiod=midperiod)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{programming}
\keyword{IO}% __ONLY ONE__ keyword per line
