#' @title Minmax Scale
#'
#' @description Scales values of numeric vector to be between 0 and 1
#'
#' @param \code{x}     numeric vector
#'
#' @return \code{x} values from input vector scaled to be between 0 and 1
#' @noRd
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
#' @importFrom imputeTS na_interpolation
#' @importFrom tsibble difference
#' @importFrom lubridate ymd
#' @importFrom imputeTS na_interpolation
#'
#' @import dplyr
#'
#'
#' @examples
#' # The following code shows the common way to use date_range and covid_lag to transform the data
#' Hyperion %>% date_range(start ='2020-04-20' , end = '2020-05-25') %>% covid_lag(offset =5 ,n=1, CC=F)
#'
#'
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




#'
#' @title ggcorr_offset
#'
#' @description The function ggcor_offset takes in a vector of \code{offsets} and other
#'  parameters from \code{\link{covid_lag}} as well as ... parameters from \code{\link{date_range}}
#'  and creates a correlation plot matrix. The rows of the matrix are the 4 normalization methods and the
#'  cols are the number of days the concentrations have been offset.
#'
#' @param data Formatted Hyperion or Pt.Loma df
#' @param range Vector of offset values to shift concentrations
#' @param target Character 'N1' or 'N2' "assays"
#' @param n Number of days to difference
#' @param CC Boolean for if concentrations are to be differenced
#' @param ... Parameters from \code{\link{date_range}}
#'
#' @return Prints plot with normalization method across rows and offsets across cols
#' @export
#'
#' @import zoo
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import ggplot2
#'
#' @examples
#' #Plot offset concentrations and cases by 1,3,5,7 days
#' Hyperion %>% ggcorr_offset(c(1,3,5,7) , target = 'N1' , start = '2020-05-25' , end = '2020-10-20')
ggcorr_offset <- function(data, range, target, n = 1L, CC = F,...) {

    var.title = ifelse(CC, " Change", "")
    lag.title <- ifelse(CC, "Lag: ", "")
    norm_methods <- c("Unadjusted", "BoCoV", "BoCoV + PMMV", "PMMV")
    target_labels <- glue::glue("{target} ({norm})", target = target, norm = norm_methods)
    names(target_labels) <- data %>% select(starts_with(target)) %>% colnames() %>% sort()

    data <- data %>% select(Date, cases, starts_with(target))
    data <- Vectorize(function(i) covid_lag(data, n = n, offset = i, CC = CC) %>% date_range(...),
                      vectorize.args = "i", SIMPLIFY = FALSE)(i = range) %>%
      lapply(minmaxDF) %>%
      bind_rows() %>% pivot_longer(cols = !c(Date, Offset, Lag, Cases_Offset, CC)) %>%
      mutate(Offset = as.numeric(Offset)) %>%
      filter(str_detect(name, target))
    g <- data %>% ggpubr::ggscatter(x = "value", y = "Cases_Offset",  add = "reg.line", size = 0.9,
                                    conf.int = TRUE, cor.coef = TRUE, facet.by = c("name", "Offset"),
                                    panel.labs = list(name = target_labels),
                                    add.params = list(fill = "steelblue"),
                                    repel = TRUE, cor.coef.size = 3,
                                    xlab = "Concentration", ylab =  "Case Counts",
                                    title = glue::glue("COVID-19 Case Counts  vs. {target} Concentration{var.title} Offsets",
                                                       target = target, var.title = var.title),
                                    subtitle = glue::glue("{lag.title}{Lag} ({start} to {end})",
                                                          Lag = ifelse(CC, n, ""),
                                                          start = zoo::as.yearmon(min(data$Date)),
                                                          end = zoo::as.yearmon(max(data$Date))),
                                    ylim= c(0,1),
                                    ggtheme = theme_bw())
    print(g)

}


# Input:
# - dataframe with date and time series
# - target: N1, N2
# - k: number of days for window
# - lag: lagged difference
# - offset: number of days to shift
# Note: must provide completed time series or imputed
# Output:
# - ggplot

#' @title ggkeystone
#'
#' @description The ggkeystone function takes in arguments from covid_lag as well as
#' date_range. It also takes in a parameter k which represents the number of days to take a rolling
#' average over. It returns the "keystone" graphic showing the relationship between waste water concentration
#' and new Covid cases.
#'
#' @param data Formatted Hyperion or Pt.Loma df
#' @param target Character 'N1' or 'N2' "assays"
#' @param k Number of days for rolling average.
#' @param lag Number of days for lagged difference from n parameter in \code{\link{covid_lag}}
#' @param offset Number of days to shift concentrations from \code{covid_lag}
#' @param CC CC parameter from \code{covid_lag}
#' @param ... \code{\link{date_range}}
#'
#' @return Returns keystone plot across all 4 normalization methods
#' @export
#'
#' @importFrom glue glue
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @importFrom imputeTS na_kalman
#'
#' @examples
#' #note should be used with some imputation method such as na_kalman()
#' library(imputeTS)
#' Hyperion %>% na_kalman() %>% ggkeystone(target = "N1") + ylab("")
ggkeystone<- function(data, target, k = 7, lag = 1L, offset = 0, CC = F, ...){

  # Add: Raise a warning if too many missing values
  # Do we want it so that we can get a flexible dataset?
  # or is it always going to be comparing the four methods?
  # for now it requires all 4 norm methods to be provided

    data <- data %>% select(Date, cases, starts_with(target)) %>%
      covid_lag(n = lag, offset = offset, CC = CC) %>%
      # left_join(data %>% select(Date, starts_with(target))) %>%
      mutate(Cases = Cases_Offset) %>%
      minmaxDF() %>%
      mutate(across(where(is.numeric) & !Cases, ~rollmean(.x, k = k, fill = NA))) %>%
      date_range(...)

    norm_methods <- c("Unadjusted", "BoCoV", "BoCoV + PMMV", "PMMV")
    target_labels <- glue::glue("{target} ({norm})", target = target, norm = norm_methods)
    names(target_labels) <- data %>% select(starts_with(target)) %>% colnames() %>% sort()

    data %>%
      pivot_longer(cols = !c(Date, Cases_Offset, Lag, Offset, CC, Cases),
                   names_to = "Target", values_to = "Values") %>%
      group_by(Target) %>%
      # mutate(corr = cor.test(Values, Cases)$estimate, p.val = cor.test(Values, Cases)$p.value) %>%
      ggplot(aes(Date, Cases)) +
      geom_col(fill = "lightblue", alpha = 0.9, position = "identity") +
      geom_line(aes(Date, Cases_Offset, color = "steelblue"), size = 0.8) +
      geom_line(aes(Date, Values, color = "indianred3"), size = 0.8, linetype = 'F1') +
      # geom_text(aes(x = ymd('2020-05-30'), y = 0.9,
      # label = glue::glue('R = {r}\np {p}', r = round(corr, 2),
      # p = ifelse(p.val < 0.05, '< 0.05', paste0("=", round(p.val, 4))))),
      # size = 4, check_overlap = TRUE) +
      facet_wrap(~Target, labeller = labeller(Target = target_labels)) +
      labs(subtitle = glue::glue("COVID-19 cases    Lag: {Lag}    Offset: {Offset}    Moving Average: {k}   Date Range: ({start} to {end})",
                                 Lag = lag, Offset = offset, k = k,
                                 start = zoo::as.yearmon(min(data$Date)),
                                 end = zoo::as.yearmon(max(data$Date))),
           title = "Daily new cases (blue) and SARS-CoV-2 concentrations (red)
           with different normalization methods and scaled",
           x = '', y = '') +
      scale_colour_manual(name = '',
                          values = c('steelblue' = 'steelblue','indianred3' = 'indianred3'),
                          labels = c( 'Concentration','COVID-19 Cases')) +
      # scale_x_date(breaks = scales::breaks_pretty(n = 13),
      #              labels = scales::label_date_short()) +
      theme_bw()+
      theme(legend.position = 'bottom', legend.margin = margin(t = 0, unit = 'cm'),
            legend.box.margin = margin(-10, 0, 0, 0))

}


