library(tidyverse, warn.conflicts = FALSE)

# library to read matlab data formats into R
library(reshape2)
library(lubridate)

# set strings as factors to false
options(stringsAsFactors = FALSE)

# Machine Learning Analysis Packages
# library(dismo)
# library(gbm)

data <- read.csv("../DATA/January/Combined_CARIACO_data_v6.csv")

gmldat <- data %>%
  select('mcc', 'mwp', 'sp', 'cdir', 'sst','tp', 'Satellite_chla') %>%
  drop_na()

gmldat2 <- data %>%
  select('mcc', 'mwp', 'sp', 'cdir', 'sst','tp', 'Chlorophyll_35m') %>%
  drop_na()

library(fitdistrplus)
descdist(gmldat$Satellite_chla, discrete = FALSE)

descdist(log(gmldat$Satellite_chla), discrete = FALSE)

descdist(gmldat2$Chlorophyll_35m, discrete = FALSE)

descdist(log(gmldat2$Chlorophyll_35m), discrete = FALSE)
