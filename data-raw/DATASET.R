## code to prepare `DATASET` dataset goes here
library(readr)
library(tibble)
library(magrittr)


Hyperion = read_csv("C:/Users/atki2/OneDrive/Documents/SCCWRP_PROJECT/SCCWRP_PROJECT/Data/main.csv") %>%
            as.data.frame()
Pt.Loma = read_csv("C:/Users/atki2/OneDrive/Documents/SCCWRP_PROJECT/SCCWRP_PROJECT/Data/main_pl.csv") %>%
            as.data.frame()



usethis::use_data(Hyperion, overwrite = TRUE , compress = 'xz')
usethis::use_data(Pt.Loma, overwrite = TRUE, compress='xz')
