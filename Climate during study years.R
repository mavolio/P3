library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(SPEI)

theme_set(theme_bw(12))


ppt<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/APT011.csv")%>%
  filter(ppt!="."&watershed=="HQ")%>%
  separate(RecDate, into=c("month", "day", "year"), sep="/")%>%
  mutate(ppt2=as.numeric(ppt))%>%
  filter(year>2009&year<2016)%>%
  group_by(year)%>%
  summarize(map=sum(ppt2))

MAP<-mean(ppt$map)


climate<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/AWE011.csv")%>%
  select(RECYEAR, RECMONTH, RECDAY, WATERSHED, RECHOUR, TAIR)%>%
  mutate(temp=as.numeric(TAIR))%>%
  filter(RECYEAR>2009&RECYEAR<2016)%>%
  na.omit%>%
  group_by(RECYEAR)%>%
  summarize(mat=mean(temp))

MAT<-mean(climate$mat)


