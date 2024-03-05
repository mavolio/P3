library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(codyn)
library(lme4)
library(car)

theme_set(theme_bw(12))

#set wd
setwd("C:/Users/mavolio2/Dropbox/Konza Research")
setwd("E:/Dropbox/Konza Research")

#read in data
spdata<-read.csv("pplots\\Sppcomp\\Species Comp_to Use\\Compiling Datasets in R\\Spp_Data_2002_2021.csv") %>% 
  mutate(present=1)

traits<-read.csv("pplots\\Sppcomp\\Species Comp_to Use\\traits_2023_categoricalonly.csv")

trts<-read.csv('pplots\\treatments.csv') %>% 
  rename(plot_id=plotnum)

pretrtabund<-spdata %>% 
  filter(calendar_year==2002)

##spdata wide

spdatawide<-spdata %>% 
  select(-abundance) %>% 
  pivot_wider(names_from = calendar_year, names_prefix = "y", values_from = present, values_fill=0) %>% 
  select(plot_id, treatment, spnum, genus_species, y2002, y2003, y2004, y2005, y2006, y2007, y2008, y2009, y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019, y2021) %>% 
  mutate(yrspres=rowSums(spdatawide[5:23])) %>% 
  mutate(outcome=ifelse(y2002>0&yrspres==1|y2002>0&y2003>0&yrspres==2|y2002>0&y2003>0&y2004>0&yrspres==3|y2002>0&y2003>0&y2004>0&y2005>0&yrspres==4|y2002>0&y2003>0&y2004>0&y2005>0&y2006>0&yrspres==5|y2002>0&y2003>0&y2004>0&y2005>0&y2006>0&y2007>0&yrspres==6|y2002>0&y2003>0&y2004>0&y2005>0&y2006>0&y2007>0&y2008>0&yrspres==7|y2002>0&y2003>0&y2004>0&y2005>0&y2006>0&y2007>0&y2008>0&y2009>0&yrspres==8|y2002>0&y2003>0&y2004>0&y2005>0&y2006>0&y2021==0&y2019==0&y2018==0&y2017==0&y2016==0, "loss", 
                 ifelse(yrspres>17, 'persist', 
                  ifelse(y2002==0&y2021>0&y2019>0&y2018>0&y2017>0&y2016>0|y2002==0&y2021>0&y2019>0&y2018>0&y2017>0, 'gain', 'blip')))) %>% 
  mutate(loss=ifelse(outcome=='loss', 1, 0),
         gain=ifelse(outcome=='gain', 1, 0),
         persist=ifelse(outcome=="persist", 1, 0)) %>% 
  left_join(pretrtabund) %>% 
  left_join(trts) %>% 
  left_join(traits) %>% 
  filter(outcome!='blip') %>% 
  filter(lifefrom!="S")


spdatawide<-spdata %>% 
  select(-abundance) %>% 
  pivot_wider(names_from = calendar_year, names_prefix = "y", values_from = present, values_fill=0) %>% 
  select(plot_id, treatment, spnum, genus_species, y2002, y2003, y2004, y2005, y2006, y2007, y2008, y2009, y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019, y2021) %>% 
  mutate(yrspres=rowSums(spdatawide[5:23])) %>% 
  mutate(outcome=ifelse(y2002>0&yrspres==1|y2002>0&y2003>0&yrspres==2|y2002>0&y2003>0&y2004>0&yrspres==3|y2002>0&y2003>0&y2004>0&y2005>0&yrspres==4, "loss", 
                        ifelse(yrspres>17, 'persist', 
                               ifelse(y2002==0&y2021>0&y2019>0&y2018>0&y2017>0&y2016>0|y2002==0&y2021>0&y2019>0&y2018>0&y2017>0, 'gain', 'blip')))) %>% 
  mutate(loss=ifelse(outcome=='loss', 1, 0),
         gain=ifelse(outcome=='gain', 1, 0),
         persist=ifelse(outcome=="persist", 1, 0)) %>% 
  left_join(pretrtabund) %>% 
  left_join(trts) %>% 
  left_join(traits) %>% 
  filter(outcome!='blip') %>% 
  filter(lifefrom!="S")
####
##Lifeform
###

lost<-spdatawide %>% 
  select(genus_species, outcome, nitro, phos) %>% 
  group_by(genus_species, outcome, nitro, phos) %>% 
  summarise(n=length(outcome))

#NP
loss_lifeformNP <- glm(loss ~ lifefrom*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(loss_lifeformNP)

plotlifeform<-spdatawide %>% 
  mutate(loss2=ifelse(lifefrom=="F"&loss==0, 0.04, ifelse(lifefrom=='G'&loss==0, 0.02, ifelse(lifefrom=="G"&loss==1, 0.98, loss)))) %>% 
  filter(!is.na(lifefrom))

ggplot(data=subset(plotlifeform, nitro==2&phos>0), aes(x=abundance, y=loss2, color=lifefrom))+
  geom_point()+
  scale_color_manual(name='Lifeform', values=c('deeppink', 'forestgreen', 'brown'))+
  geom_smooth(aes(y=loss), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Loss Probability")

####
##Lifespan
###

#NP
loss_lifespanNP <- glm(loss ~ lifespan*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(loss_lifespanNP)

plotlifespan<-spdatawide %>% 
  mutate(loss2=ifelse(lifespan=="A"&loss==0, 0.02, ifelse(lifespan=='P'&loss==1, 0.98, loss)))

ggplot(data=subset(plotlifespan, nitro==2&phos>0), aes(x=abundance, y=loss2, color=lifespan))+
  geom_point()+
  scale_color_manual(name='Lifespan', values=c('deeppink', 'forestgreen'))+
  geom_smooth(aes(y=loss), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Loss Probability")

####
#N-fix
###

#NP
loss_nifxNP <- glm(loss ~ N_fix*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(loss_nifxNP)

plotnfix<-spdatawide %>% 
  mutate(loss2=ifelse(N_fix=="N"&loss==0, 0.02, ifelse(N_fix=='Y'&loss==1, 0.98, loss)))

ggplot(data=subset(plotnfix, nitro==2&phos>0), aes(x=abundance, y=loss2, color=N_fix))+
  geom_point()+
  scale_color_manual(name='N Fixer', values=c('orange', 'forestgreen'))+
  geom_smooth(aes(y=loss), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Loss Probability")



####
#photo path
###

#NP
loss_ppNP <- glm(loss ~ photo_path*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(loss_ppNP)

plotps<-spdatawide %>% 
  mutate(loss2=ifelse(photo_path=="C3"&loss==0, 0.02, ifelse(photo_path=='C4'&loss==1, 0.98, loss)))

ggplot(data=subset(plotps, nitro==2&phos>0), aes(x=abundance, y=loss2, color=photo_path))+
  geom_point()+
  scale_color_manual(name='Photo Path', values=c('lightgreen', 'forestgreen'))+
  geom_smooth(aes(y=loss), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Loss Probability")
########################################################3
#######################################################33
###PERSISTENCE

####
##Lifeform
###

#NP
persist_lifeformP <- glm(persist ~ lifefrom*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro==1))
#summary(loss_lifeformNP)
Anova(persist_lifeformP)

persist_lifeformNP <- glm(loss ~ lifefrom*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(persist_lifeformNP)

plotlifeform<-spdatawide %>% 
  mutate(persist2=ifelse(lifefrom=="F"&persist==0, 0.04, ifelse(lifefrom=='G'&persist==0, 0.02, ifelse(lifefrom=="G"&persist==1, 0.98, persist)))) %>% 
  filter(!is.na(lifefrom))

ggplot(data=subset(plotlifeform, nitro==1&phos>0), aes(x=abundance, y=persist2, color=lifefrom))+
  geom_point()+
  scale_color_manual(name='Lifeform', values=c('deeppink', 'forestgreen', 'brown'))+
  geom_smooth(aes(y=persist), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Persistence Probability")

####
##Lifespan
###

#NP
persist_lifespanP <- glm(persist ~ lifespan*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro==1))
#summary(loss_lifeformNP)
Anova(persist_lifespanP)

persist_lifespanNP <- glm(persist ~ lifespan*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(persist_lifespanNP)

plotlifespan<-spdatawide %>% 
  mutate(persist2=ifelse(lifespan=="A"&persist==0, 0.02, ifelse(lifespan=='P'&persist==1, 0.98, persist)))

ggplot(data=subset(plotlifespan, nitro==2&phos>0), aes(x=abundance, y=persist2, color=lifespan))+
  geom_point()+
  scale_color_manual(name='Lifespan', values=c('deeppink', 'forestgreen'))+
  geom_smooth(aes(y=persist), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Persistence Probability")

####
#N-fix
###

#NP
persist_nifxP <- glm(persist ~ N_fix*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro==1))
#summary(loss_lifeformNP)
Anova(persist_nifxP)

persist_nifxNP <- glm(persist ~ N_fix*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(persist_nifxNP)

plotnfix<-spdatawide %>% 
  mutate(persist2=ifelse(N_fix=="N"&persist==0, 0.02, ifelse(N_fix=='Y'&persist==1, 0.98, persist)))

ggplot(data=subset(plotnfix, nitro==1&phos>0), aes(x=abundance, y=persist2, color=N_fix))+
  geom_point()+
  scale_color_manual(name='N Fixer', values=c('orange', 'forestgreen'))+
  geom_smooth(aes(y=persist), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Persistence Probability")



####
#photo path
###

#NP
persist_ppP <- glm(persist ~ photo_path*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro==1))
#summary(loss_lifeformNP)
Anova(persist_ppP)

persist_ppNP <- glm(persist ~ photo_path*abundance, family = binomial(), data = subset(spdatawide, phos!=0&nitro!=1))
#summary(loss_lifeformNP)
Anova(persist_ppNP)

plotps<-spdatawide %>% 
  mutate(persist2=ifelse(photo_path=="C3"&persist==0, 0.02, ifelse(photo_path=='C4'&persist==1, 0.98, persist)))

ggplot(data=subset(plotps, nitro==1&phos>0), aes(x=abundance, y=persist2, color=photo_path))+
  geom_point()+
  scale_color_manual(name='Photo Path', values=c('lightgreen', 'forestgreen'))+
  geom_smooth(aes(y=persist), method = 'glm', method.args=list(family='binomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Pre-treatment Abundance")+
  ylab("Persistence Probability")
