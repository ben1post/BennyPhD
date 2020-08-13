library(tidyverse, warn.conflicts = FALSE)
library(cowplot, warn.conflicts = FALSE)
library(scales, warn.conflicts = FALSE)

# library to read matlab data formats into R
library(R.matlab)

library(reshape2)

library(lubridate)

library(tidyverse)
library(worrms)

setwd("~/Documents/GitHub/BennyPhD/CaRiaco/SCRIPPS/CARIACO_data/PrelimDataExploration")
ds <- read.csv("../VeryNEWESTCariacoData/phytoplankton.csv")

ds$Datetime_UTC = parse_date_time(ds$Datetime_UTC, orders = "%Y-%m-%d H:M:S")
ds$Datetime_local = parse_date_time(ds$Datetime_local, orders = "%Y-%m-%d H:M:S")

ds$date = ds$Datetime_UTC

head(ds$Datetime_UTC)

head(ds$AphiaID)

blob <- rbind(data.frame(val=ds$d_1m, depth=1, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_7m, depth=7, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_15m, depth=15, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_25m, depth=25, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_55m, depth=55, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_75m, depth=75, date=ds$date, AphiaID=ds$AphiaID),
      data.frame(val=ds$d_100m, depth=100, date=ds$date, AphiaID=ds$AphiaID))

head(blob)

str(blob)

min(blob$AphiaID)

test = merge(blob, ds[-19:-26])

head(test)

str(test)

test$AphiaID[test$AphiaID <0] <- NA

library(worrms)

wm_attr_data(id = test$AphiaID[200], include_inherited=TRUE)

wm_children(id = test$AphiaID[200])

is_tibble(wm_classification(id = id))

tryCatch({
    wm_data = wm_classification(id = id)
    phylum = wm_data[wm_data$rank=="Infraphylum",]$scientificname
}, error=function(e) phylum=NA)

phylum

ds$AphiaID

id = test$AphiaID[1]

if (is_tibble(wm_classification(id = id))){
    data = wm_classification(id = id)
    phylum = data[data$rank=="Infraphylum",]$scientificname}

data

allID = length(ds$AphiaID)

as.numeric(allID)

for (i in 1:allID){
    print(i)}




###############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

getclassification <- function(x){
    if (x > 0){
      wm_classification = wm_classification(id = x)
      wm_class = paste(wm_classification$scientificname, collapse="_")
    } else {
      wm_class = NA
    }
  
    #print(c('wm_class output',wm_class))
    return(wm_class)
}

getattrs <- function(x){
  
    #print(c('getattrs', x))
  
    wm_attr="NA"
    
    result <- try(wm_attr_data(x), include_inherited = TRUE)#, silent=TRUE)
    #print(class(result))
    
    # Process any error messages
    if (all( class(result) != "try-error")){
      wm_attr = paste(result$measurementType, result$measurementValue, sep=":", collapse=", ") }        
    
    #print(c('wm_attr output',wm_attr, class(wm_attr)))
    return(wm_attr)
}

worms_dat <- ds %>%
  group_by(AphiaID) %>%
  summarise(aphiaid = AphiaID[1], speciesorig = SpeciesNameOriginal[1], speciesclean = SpeciesNameCleaned[1],scientificname = ScientificName_accepted[1],wm_attrs=getattrs(AphiaID[1]), wm_classes=getclassification(AphiaID[1]))

worms_dat


write.csv(worms_dat, "worms_dat_01.csv")




# When this code is done, copy it to a clean jupyter notebook!




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spp_attr <- ds %>%
  #rowwise() %>%
  group_by(AphiaID) %>%
  mutate(wm_class = getclassification(AphiaID[1]))

# wm_class
sum(grepl("Dinofl",spp_attr$phylum, fixed=TRUE))/length(spp_attr$phylum)

sum(grepl("Ochro",spp_attr$phylum, fixed=TRUE))/length(spp_attr$phylum)


spp_attr <- ds %>%
  group_by(AphiaID) %>% 
  group_map(~ get_worms_fgrp(AphiaID=.x$AphiaID)) %>% 
  bind_rows()

testfunc <- function(){
  print("hello")
  print("XXY")
}

testfunc()


source("get_functional_group.r")

library(reshape2)

result <- melt(blob, id.vars=c('depth','val'))

result#[order(result$value),]

data.frame(id=rep(1:4,each=2),
                    z2=sample(letters,8, replace=TRUE))

data.frame(D=as.array(t(ds$d.1m)))

trygetclassification <- function(x){
  phylum = "Uninit"
  print(x)
  #print(x<0)
  if (x > 0){
    #tryCatch({
    wm_data = wm_classification(id = x)
    print(wm_data)
    #phylum = wm_data[wm_data$rank=="Infraphylum",]$scientificname
    phylum = "XXX"
    #}, error= {phylum="NA"}) 
  } else {
    phylum="Aphia number missing"
  }
  return(phylum)
}
