#' @title Minmax Scale
#'
#' @description Scales values of numeric vector to be between 0 and 1
#'
#' @param \code{x}     numeric vector
#'
#' @return \code{x} values from input vector scaled to be between 0 and 1
#' @export
#'
#'
#' @examples
#' x <- runif(10, min =0 ,max=100)
#' minmax(x)
#'
#'
minmax <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}




#' @title Minmax DF
#'
#' @description Apply minmax function to 0-1 scale all numeric columns of a dataframe
#'
#' @param data Input df
#'
#' @return minmax scaled df
#' @export
#'
#' @import dplyr
#' @examples
#' data <- data.frame('A' = runif(n=25,min=50,max=100),'B'=runif(n=25,min=50,max=100))
#' minmaxDF(data = data )
minmaxDF <- function(data) {
  if(require(dplyr)) {
    tmp <- data %>% dplyr::mutate(across(where(is.numeric), ~minmax(.x)))
    return(tmp)
  }
}




#' @title LagLead
#'
#' @description Runs \code{\link[dplyr]{lag}} or \code{\link[dplyr]{lead}} functions
#' depending on if offset is positive or negative. Positive values of offset
#' shift data forward. Negative values shift data backwards.
#'
#'
#' @param x Data
#' @param offset Number of intervals index data is to be shifted
#'
#' @return For offset > 0, data with index shifted forward by offset.
#' For Offset < 0, data shifted back by offset
#'
#' @export
#'
#' @examples
#' #Shift data forward by offset = 2
#' data <- data.frame('A' = runif(n=25,min=50,max=100),'B'=runif(n=25,min=50,max=100))
#' laglead(x = data , 2)
#'
#' #Shift data backward by offset = -2
#' laglead(x = data , -2)
#' @import dplyr
laglead <- function(x,offset){
  if(offset < 0 ){return(dplyr::lead(x,abs(offset)))}
  else{return(dplyr::lag(x,offset))}
}




#' @title Date_Range
#'
#' @description Takes df or ts, start, and end date range and returns df filtered to observations
#' in that time period. Additionally, can take in weekday argument and will only return observations
#' occurring on that weekday. Also will work with abbreviated weekday such as 'Tue'. Function is also pipe
#' compatible with Tidy Verse
#'
#' @param data df or time series with Date column
#' @param start Character data type expressing date of form 'YYYY-MM-DD'
#' @param end Character data type of form 'YYYY-MM-DD'
#' @param weekday Character name of weekday e.g. 'Tuesday' or 'Tue'
#'
#'
#' @return Returns df filtered to date range defined by inputs
#' @export
#'
#' @examples
#' data <- data.frame(Date= seq(as.Date("2021/04/20"), by = "day", length.out = 100), Obs = rnorm(100))
#'
#' #Filter for all days between start and end date
#' date_range(data , start = '2021-05-10', end = '2021-06-15')
#'
#' #Filter for every Tuesday between start and end
#' date_range(data , start = '2021-05-10', end = '2021-06-15',weekday = 'Tuesday')
#'
#' #Example with pipe
#' data %>%
#'
#'
#'
date_range <- function(data, start = NULL, end = NULL, weekday=NULL) {

  if(is.null(start)) start <- data$Date[1]
  if(is.null(end)) end <- data$Date[nrow(data)]

  tmp <- data %>% filter(between(Date, ymd(start), ymd(end)))

  if(!is.null(weekday)) {
    tmp <- tmp %>% filter((weekdays(Date) == weekday | weekdays(Date, abbreviate = T) == weekday))
  }
  return(tmp)
}




#' @title Covid Lag
#'
#' @description Covid Lag takes in Main 'data' which has an aggregate total of cases and will
#' calculate lag differences between cases over 'n' days. Covid lag will also take in an offset
#' value and shift the concentration measurements by 'offset' days to line up with the new
#' case count. For example if offset = 7 concentration measurements recorded on '2020-04-20'
#' will be shifted to align with new cases recorded on '2020-04-27'. In addition the CC
#' parameter can be set to TRUE to take the difference in concentrations as well.
#'
#'
#' @param data Main data for Pt loma or Hyperion
#' @param offset Number of days to shift concentrations
#' @param n Number of days of differencing for running cases count to go to new cases
#' @param CC CC = TRUE to also difference concentration n days
#'
#' @return Df with Cases_Offset column recording the new case count for the given day
#' @export
#'
#'
#' @importFrom tsibble difference
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' need to attach data for example
#' }
covid_lag <- function(data, offset = 0, n = 1L, CC = F) {

    tmp <- data %>% arrange(Date) %>%
      mutate(across(starts_with("N"), ~laglead(.x, offset)),
             Cases_Offset = difference(cases, n),
             Lag = as.character(n), Offset = as.character(offset), CC = CC,
             Cases_Offset = floor(na_interpolation(ifelse(Cases_Offset < 1, NA, Cases_Offset)))) %>%
      select(Date, starts_with("N"), Cases_Offset, Lag, Offset, CC)
    if(CC) {
      tmp <- tmp %>% arrange(Date) %>% mutate(across(starts_with("N"), ~difference(.x, n)))
    }
    return(tmp)

}
