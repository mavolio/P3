library(tidyverse)


theme_set(theme_bw(12))

ppt<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/APT011.csv")%>%
  filter(ppt!=".")%>%
  separate(RecDate, into=c("month", "day", "year"), sep="/")%>%
  mutate(ppt2=as.numeric(ppt))


#based on 2C growing season precip
annprecip<-ppt%>%
    filter(watershed=="HQ")%>%
    group_by(year, watershed)%>%
    summarize(sppt=sum(ppt2))%>%
    mutate(year=as.integer(year))

mean<-mean(annprecip$sppt)
sd<-sd(annprecip$sppt)

temp<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/AWE011.csv")%>%
  rename(month=RECMONTH,
         year=RECYEAR)%>%
  filter(TAIR!=""&TAIR!="."&TAIR!=-116)%>%
  mutate(TAIR2=as.numeric(TAIR))

temp2<-temp%>%
  group_by(month, year, RECDAY)%>%
  summarize(tave=mean(TAIR2))%>%
  group_by(year, month)%>%
  summarize(Tmed=median(tave))


knz<-ppt%>%
  filter(watershed=="HQ")%>%
  group_by(year, month)%>%
  summarize(ppt=sum(ppt2))%>%
  mutate(year=as.integer(year),
         month=as.integer(month))%>%
  left_join(temp2)%>%
  na.omit%>%
  arrange(year, month)


##what was rainfall during our experiment drought year
knzSub<-knz %>% 
  filter(year %in% c(2009, 2010, 2011, 2012, 2013, 2014, 2015)) %>% 
  mutate(time=ifelse(year==2010&month>3&month<11|year==2011&month>3&month<11|year==2012&month>3&month<11, "on", "off"),
         drop=ifelse(year==2009&month<11, 1, 0)) %>% 
  filter(drop!=1)

##For table 2
#2009-2010
ppt1<-knzSub %>% 
  mutate(keep=ifelse(year==2009, 1, ifelse(year==2010&month<4&month, 1, 0))) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt2<-knzSub %>% 
  mutate(keep=ifelse(year==2010&time=='on', 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt3<-knzSub %>% 
  mutate(keep=ifelse(year==2010&month>10|year==2011&month<4, 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt4<-knzSub %>% 
  mutate(keep=ifelse(year==2011&time=='on', 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt5<-knzSub %>% 
  mutate(keep=ifelse(year==2011&month>10|year==2012&month<4, 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt6<-knzSub %>% 
  mutate(keep=ifelse(year==2012&time=='on', 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt7<-knzSub %>% 
  mutate(keep=ifelse(year==2012&month>10|year==2013&month<11, 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt8<-knzSub %>% 
  mutate(keep=ifelse(year==2013&month>10|year==2014&month<11, 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))

ppt9<-knzSub %>% 
  mutate(keep=ifelse(year==2014&month>10|year==2015&month<11, 1, 0)) %>% 
  filter(keep==1) %>% 
  ungroup() %>% 
  mutate(sum=sum(ppt))




