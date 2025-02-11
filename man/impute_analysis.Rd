% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sccwrpwbe_functions.R
\name{impute_analysis}
\alias{impute_analysis}
\title{Impute Analysis}
\usage{
impute_analysis(
  ts_comp,
  rate,
  seeds,
  option_locf = "locf",
  option_interpol = "linear",
  model = "StructTS",
  k = 2,
  option_mean = "mean"
)
}
\arguments{
\item{ts_comp}{Complete Time Series to run algorithm}

\item{rate}{Rate of missing data to simulate}

\item{seeds}{Numnber of iterations to run}

\item{option_locf}{Optional parameters for \code{\link{na_locf}}}

\item{option_interpol}{Optional parameters for \code{\link{na_interpol}}}

\item{model}{model}

\item{k}{K}

\item{option_mean}{Optional Parameter for \code{\link{na_mean}}}
}
\value{
tibble
}
\description{
Function takes in a complete time series (\code{ts_comp}),
a rate of missing data (\code{rate}), and an integer (\code{seeds}) of the number of trials to
run an imputation algorithm which will compare 6 methods of imputing missing time series data.
Both MAPE and RMSE will be returned. This function is used in \code{\link{ggimpute_eval}} and
\code{\link{ggimpute}}
}
\examples{
Hyperion \%>\% date_range(start='2020-04-20' , end = '2020-07-01') \%>\% covid_lag() \%>\% pull(Cases_Offset)\%>\%  impute_analysis(rate = 0.6,seeds = 30)
}
