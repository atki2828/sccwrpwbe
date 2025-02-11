% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sccwrpwbe_functions.R
\name{minmaxDF}
\alias{minmaxDF}
\title{Minmax DF}
\usage{
minmaxDF(data)
}
\arguments{
\item{data}{Input df}
}
\value{
minmax scaled df
}
\description{
Apply minmax function to 0-1 scale all numeric columns of a dataframe
}
\examples{
data <- data.frame('A' = runif(n=25,min=50,max=100),'B'=runif(n=25,min=50,max=100))
minmaxDF(data = data )
}
