library(tidyverse, warn.conflicts = FALSE)

# library to read matlab data formats into R
library(reshape2)
library(lubridate)

# set strings as factors to false
options(stringsAsFactors = FALSE)

# Machine Learning Analysis Packages
library(dismo)
library(gbm)

install.packages( pkgs = "terra"
                  , method = "curl"
                  , configure.args = c(
                    "--with-gdal-config=/Library/Frameworks/GDAL.framework/Programs/gdal-config"
                    , "--with-proj-include=/bpo/home/bin/proj4/include"
                    , "--with-proj-lib=/bpo/home/bin/proj4/lib"
                  ) 
)