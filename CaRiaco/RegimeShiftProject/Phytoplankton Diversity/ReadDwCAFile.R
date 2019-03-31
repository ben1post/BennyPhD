library('finch')
library('worrms')

zip = 'Phytoplankton Diversity/dwca-cariaco-phytoplankton-v1.2.zip'
(out = dwca_read(zip, read=TRUE))


str(out$data[[1]])
str(out$data[[2]])

species = out$data[[1]]
counts = out$data[[2]]

species$scientificNameID

wm_distribution(id = 235923)
wm_classification(id = 235802)

diversitydata <- merge(counts, species)

diversitydata$scientificName[diversitydata$taxonRank=="Phylum"]

diversitydata$scientificName[diversitydata$taxonRank=="Class"]

levels(as.factor(diversitydata$taxonRank))

which(diversitydata$scientificName == "Nanoflagellate", arr.ind = T)

#PROBLEM: diversitydata does not include nanoflagellates, 
#and does not include functional types! that is why I will for the 1st MS, for now, will use the BCO data
