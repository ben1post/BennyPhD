# CARICAO
library(ncdf4)
library(chron)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(scales)

# read data from ncdf and plot
filename <- "../../../../Uni/Research/PhD/First Year/TimeSeriesData/BCO CARIACO/Phytoplankton.netcdf"
# large file not on github, link to data: 'https://www.bco-dmo.org/dataset/3095'

phydat<-nc_open(filename)


CruiseID <- ncvar_get(phydat, 'Cruise_ID')
lon<-ncvar_get(phydat, 'lon')
lat<-ncvar_get(phydat, 'lat')
dateRaw <-ncvar_get(phydat, 'Date')
taxon <- ncvar_get(phydat, 'taxon')
group <- ncvar_get(phydat, 'taxon_group')
depth <- ncvar_get(phydat, 'depth')
abundance <- ncvar_get(phydat,'abundance')
month <- ncvar_get(phydat,'Mon')
year <- ncvar_get(phydat, 'Year')


date <- as.POSIXlt(as.character(dateRaw), format = "%Y%m%d", tz="GMT")

# prepare final data set
dat.phy <- data.frame(CruiseID,lon,lat,date,month,year,depth,taxon,group,abundance)
nc_close(phydat)
