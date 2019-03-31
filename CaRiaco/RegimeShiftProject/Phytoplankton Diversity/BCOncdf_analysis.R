source('Phytoplankton Diversity/ReadBCOncdf.R')

#' What do I need?
#' Abundances over time per functional group, to split into regime 1 and 2
#' For now I will sum up all depth, but potentially I will have to look at phytoplankton above MLD only
#' but sum up over depth at the last step, so I can plot abundance per depth as well
#' 

# remove ChlMax value from $depth
dp1_0 <- dat.phy[!dat.phy$depth == "ChlMax",]
dp1_0$depth <- as.numeric(dp1_0$depth)
# remove clear outliers (10 abundance values at 9.96921e+36)
dp1_1 <- dp1_0[!dp1_0$abundance > 1e5,]

# remove non-sensical data, that might interfere with sums of Groups
dp1_2 <- dp1_1[!dp1_1$taxon=='SUBTOTAL',]
dp1_3 <- dp1_2[!dp1_2$taxon=='COCCOLITHOPHORIDS_SUBTOTAL',]
dp1_4 <- dp1_3[!dp1_3$taxon=='0',]
dp1_5 <- dp1_4[!dp1_4$taxon=='Number_of_species',]
dp1_6 <- dp1_5[!dp1_5$taxon=='NUMBER_OF_OTHERS',]
dp1_7 <- dp1_6[!dp1_6$taxon=='Others_abundance',]
dp1_8 <- dp1_7[!is.na(dp1_7$abundance),]

dp_filtered <- dp1_8

#sum up abundances per taxons for each date
dp2 <- dp_filtered %>% 
  group_by(date, group, taxon) %>% 
  summarise(Abundance = sum(abundance)) 

# sum up total abundances per group per date
dpGroupPerTime <- dp_filtered %>% 
  group_by(date, group) %>% 
  summarise(Abundance = sum(abundance)) 

dpGroupPerDepthPerTime <- dp_filtered %>% 
  group_by(date,depth, group) %>% 
  summarise(Abundance = sum(abundance)) 

#plot all abundances per group per date
dpGroupPerTime_notot <- dpGroupPerTime[!dpGroupPerTime$group=='TOTAL',]

head(dpGroupPerTime)
ggplot(data=dpGroupPerTime_notot,aes(x=date,y=Abundance,col=group))+geom_point()+scale_y_log10()

#plot all abundances per group per date split by depth
dpGroupPerDepthPerTime_notot <- dpGroupPerDepthPerTime[!dpGroupPerDepthPerTime$group=='TOTAL',]

head(dpGroupPerDepthPerTime_notot)
ggplot(data=dpGroupPerDepthPerTime_notot,aes(x=date,y=Abundance,col=group))+
  geom_point()+scale_y_log10() + facet_grid(depth ~ .)


#split data into two regimes:
dpR1 <- dpGroupPerDepthPerTime_notot[!dpGroupPerDepthPerTime_notot$date>as.Date('2000-10-08'),]
dpR2 <- dpGroupPerDepthPerTime_notot[!dpGroupPerDepthPerTime_notot$date<as.Date('2006-07-08'),]

dpR1$dyear <- as.Date(format(dpR1$date, format="%m-%d"),format="%m-%d")
dpR2$dyear <- as.Date(format(dpR2$date, format="%m-%d"),format="%m-%d")
dpR1$month <- as.numeric(format(dpR1$date, format="%m"))
dpR2$month <- as.numeric(format(dpR2$date, format="%m"))

#plot aggregated yearly abundances per group (all depths)
ggplot(data=dpR1,aes(x=dyear,y=Abundance,col=group))+geom_point()+scale_y_log10() + 
  geom_smooth()+ facet_grid(group ~ .,scales="free_y")
ggplot(data=dpR2,aes(x=dyear,y=Abundance,col=group))+geom_point()+scale_y_log10() + 
  geom_smooth()+ facet_grid(group ~ .)#,scales="free_y")

#plot aggregated abundances in upper MLD
dpR1upper <- dpR1[dpR1$depth>20,]
ggplot(data=dpR1upper,aes(x=dyear,y=Abundance,col=group))+geom_point()+scale_y_log10() + 
  geom_smooth()+ facet_grid(group ~ .,scales="free_y")
dpR2upper <- dpR2[dpR2$depth>20,]
ggplot(data=dpR2upper,aes(x=dyear,y=Abundance,col=group))+geom_point()+#scale_y_log10() + 
  geom_smooth()+ facet_grid(group ~ .,scales="free_y")


# plot boxplot of distribution over the year per group
ggplot(dpR1) + geom_boxplot(aes(x=reorder(month,dyear), y=Abundance, col=group))+scale_y_log10() + 
  facet_grid(group ~ .,scales="free_y")

# plot median of abundance over the year per group
ggplot(dpR1) + geom_boxplot(aes(x=reorder(month,dyear), y=Abundance, col=group))+scale_y_log10() + 
  facet_grid(group ~ .,scales="free_y")

dpR1_monthly <- dpR1 %>%
group_by(month,group) %>%
  summarise(median = median(Abundance),mean = mean(Abundance), sd = sd(Abundance))
dpR2_monthly <- dpR2 %>%
  group_by(month,group) %>%
  summarise(median = median(Abundance),mean = mean(Abundance), sd = sd(Abundance))

ggplot(dpR1_monthly) + geom_point(aes(x=month, y=mean, col=group),shape=1)+#scale_y_log10() + 
  facet_grid(group ~ .,scales="free_y") + scale_y_continuous(limits=c(0,max(dpR1_monthly$mean))) + 
  geom_point(data = dpR2_monthly,aes(x=month, y=mean, col=group),shape=2)

ggplot(dpR1_monthly) + geom_point(aes(x=month, y=median, col=group),shape=1)+#scale_y_log10() + 
  facet_grid(group ~ .,scales="free_y") + scale_y_continuous(limits=c(0,max(dpR1_monthly$median))) + 
  geom_point(data = dpR2_monthly,aes(x=month, y=median, col=group),shape=2)

ggplot(dpR1_monthly) + geom_point(aes(x=month, y=sd, col=group),shape=1)+#scale_y_log10() + 
  facet_grid(group ~ .,scales="free_y") + scale_y_continuous(limits=c(0,max(dpR1_monthly$sd))) + 
  geom_point(data = dpR2_monthly,aes(x=month, y=sd, col=group),shape=2)



#sum up total abundances per groups
groupabunds <- dp2 %>% 
  ungroup() %>% 
  group_by(group) %>% 
  summarise(abundance = sum(Abundance)) %>% 
  arrange(desc(abundance))

#sum up total abundances per taxon 
taxonabunds <- dp2 %>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

#Plot TOP taxons
taxonabunds$taxon <- factor(taxonabunds$taxon, levels = taxonabunds$taxon)

ggplot(taxonabunds[1:26,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) per Species (top 26)")

#plotGROUPS
groupabunds$group <- factor(groupabunds$group, levels = groupabunds$group)

ggplot(groupabunds, aes(x=group,y=abundance)) + geom_bar(stat = "identity")+ xlab('Group') + 
  ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+ 
  coord_flip() + ggtitle("Sum of Abundances (1996-2010) per Group")+
  theme(plot.title = element_text(hjust = 0.5))

#OTHERS
othersabunds <- dp2[dp2$group == "OTHERS",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

othersabunds$taxon <- factor(othersabunds$taxon, levels = othersabunds$taxon)

ggplot(othersabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

#DIATOMS
diatomsabunds <- dp2[dp2$group == "DIATOMS",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

diatomsabunds$taxon <- factor(diatomsabunds$taxon, levels = diatomsabunds$taxon)

ggplot(diatomsabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

#COCCOLITHOPHORIDS
coccosabunds <- dp2[dp2$group == "COCCOLITHOPHORIDS",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

coccosabunds$taxon <- factor(coccosabunds$taxon, levels = coccosabunds$taxon)

ggplot(coccosabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")


#DINOFLAGELLATES
dinosabunds <- dp2[dp2$group == "DINOFLAGELLATES",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

dinosabunds$taxon <- factor(dinosabunds$taxon, levels = dinosabunds$taxon)

ggplot(dinosabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

#NANOFLAGELLATES
nanosabunds <- dp2[dp2$group == "NANOFLAGELLATES",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

nanosabunds$taxon <- factor(nanosabunds$taxon, levels = nanosabunds$taxon)

ggplot(nanosabunds[,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

#CYANOBACTERIA
cyanosabunds <- dp2[dp2$group == "CYANOBACTERIA",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

cyanosabunds$taxon <- factor(cyanosabunds$taxon, levels = cyanosabunds$taxon)

ggplot(cyanosabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

#TOTAL
totalabunds <- dp2[dp2$group == "TOTAL",]%>% 
  ungroup() %>% 
  group_by(taxon) %>% 
  summarise(abundance = sum(Abundance), group = max(group)) %>% 
  arrange(desc(abundance))

totalabunds$taxon <- factor(totalabunds$taxon, levels = totalabunds$taxon)

ggplot(totalabunds[1:36,], aes(x=taxon,y=abundance)) + geom_bar(aes(fill=group),stat = "identity")+ 
  xlab('Species') + ylab(expression(paste("Abundance [cells ",ml^{-1},']')))+
  coord_flip() + labs(fill = "Groups") + ggtitle("Sum of Abundances (1996-2010) in OTHERS group")

