library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(SPEI)

treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")

p2<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/pplots/Biomass/To use/Compiling Files in R/Biomass_2002_2019.csv")%>%
  select(-X)%>%
  filter(calendar_year>2005,
         treatment=="N1P0"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P3")%>%
  mutate(drought="n")%>%
  rename(plotnum=plot_id)%>%
  left_join(treats)%>%
  select(-treatment, -rep)

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

ggplot(data=subset(annprecip, year>2005&year<2020), aes(x=year, y=sppt))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = mean, linetype="solid")+
  scale_x_continuous(breaks=c(2006, 2008, 2010, 2012, 2014, 2016, 2018))


#trying SPEI
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

meant<-mean(temp2$Tmed)


knz<-ppt%>%
  group_by(year, month)%>%
  summarize(ppt=sum(ppt2))%>%
  mutate(year=as.integer(year),
         month=as.integer(month))%>%
  left_join(temp2)%>%
  na.omit%>%
  arrange(year, month)%>%
  filter(year>2004)

knz$PET <- thornthwaite(knz$Tmed, 39.10216)
knz$BAL <- knz$ppt-knz$PET

knzspei<-spei(knz$BAL, 12)
plot(knzspei)

speivalues<-as.data.frame(knzspei$fitted)%>%
  bind_cols(knz)%>%
  mutate(day=paste(year, month, "01", sep="-"),
         date=as.Date(day))%>%
  filter(year>2005)%>%
  mutate(color=ifelse(PET_tho>-1, 1, 0))
 

ggplot(data=speivalues, aes(x=date, y=PET_tho, color=as.factor(color)))+
  geom_point(show.legend = F)+
  scale_color_manual(values=c("Red", "Black"))+
  geom_hline(yintercept=-1)+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("SPEI")

knzclimate<-knz%>%
  gather(climate, value, ppt:Tmed)%>%
  group_by(year, climate)%>%
  summarize(mean=mean(value))
  

ggplot(data=subset(knzclimate, year>2005), aes(x=year, y=mean))+
  geom_point(stat="identity")+
  #geom_hline(yintercept = mean, linetype="solid")+
  facet_wrap(~climate, scale="free_y", ncol=1)+
  scale_x_continuous(breaks=c(2006, 2008, 2010, 2012, 2014, 2016, 2018))

bioave<-p2%>%
  filter(Trt=="Control")%>%
  group_by(calendar_year)%>%
  summarize(manpp=mean(anpp))

ggplot(data=bioave, aes(x=calendar_year, y=manpp))+
  geom_bar(stat="identity")


###using 2012 and 2018 as drought years. these are non-burned years, flanked by burn years


drt<-p2%>%
  mutate(drt=ifelse(calendar_year==2011|calendar_year==2017, "pre", ifelse(calendar_year==2012|calendar_year==2018, "drt", ifelse(calendar_year==2013|calendar_year==2019,"post", "drop"))),
         drtint=ifelse(calendar_year %in% c(2011, 2012, 2013), 2012, ifelse(calendar_year%in% c(2017, 2018, 2019), 2018, 999)))%>%
  filter(drt!="drop")

hist(log(drt$anpp))

#appraoch 1

m1<-lmer(log(anpp)~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(drt, drtint=="2012"))
anova(m1)

m2<-lmer(log(anpp)~as.factor(calendar_year)*nitro*phos + (1|plotnum), data=subset(drt, drtint=="2018"))
anova(m2)

toplot1<-drt%>% 
  group_by(calendar_year, drt, drtint, nitro, phos, Trt)%>%
  summarize(manpp=mean(anpp), sd=sd(anpp), n=length(anpp))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=toplot1, aes(x=as.factor(calendar_year), y=manpp, color=Trt, group=Trt))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=manpp-se, ymax=manpp+se), width=0.2)+
  facet_wrap(~drtint, scales="free_x")

#approach 2
drt2<-drt%>%
  select(-calendar_year)%>%
  group_by(plotnum, drtint, nitro, phos, Trt)%>%
  spread(drt, anpp)%>%
  mutate(resist=(drt-pre)/pre,
         recov=(post-drt)/drt,
         resil=(post-pre)/pre)

summary(aov(resist~nitro*phos, data=subset(drt2, drtint==2012)))
summary(aov(recov~nitro*phos, data=subset(drt2, drtint==2012)))
summary(aov(resil~nitro*phos, data=subset(drt2, drtint==2012)))

summary(aov(resist~nitro*phos, data=subset(drt2, drtint==2018)))
summary(aov(recov~nitro*phos, data=subset(drt2, drtint==2018)))
summary(aov(resil~nitro*phos, data=subset(drt2, drtint==2018)))

toplot2<-drt2%>%
  gather(yr, value, resist:resil)%>%
  group_by(drtint, yr,Trt)%>%
  summarize(mean=mean(value), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=toplot2, aes(x=yr, y=mean, fill=Trt))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position = position_dodge(0.9))+
  scale_x_discrete(limits=c("resist", "recov", "resil"))+
  facet_wrap(~drtint)


##look at V for all three year intervals - i think the spikes in NP is actually a drought/fire signal and not fire alone

##looking for burn signal
toplot3<-p2%>% 
  group_by(calendar_year, Trt)%>%
  summarize(manpp=mean(anpp), sd=sd(anpp), n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(fire=ifelse(calendar_year %in% c(2007, 2009, 2011, 2013, 2015, 2017, 2019), "*", ""))

ggplot(data=toplot3, aes(x=as.factor(calendar_year), y=manpp, color=Trt, group=Trt, label=fire))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=manpp-se, ymax=manpp+se), width=0.2)+
  facet_wrap(~Trt)+
  geom_text(y=300)

###things to think on: can i get residuals from burn trts and do above stats on that - try this same story, but not sure i'm doing right.

burn<-drt%>% 
   mutate(fire=ifelse(calendar_year %in% c(2007, 2009, 2011, 2013, 2015, 2017, 2019), 1, 0))

burnsig<-lmer(log(anpp)~fire*as.factor(calendar_year) + (1|plotnum), data=burn)
anova(burnsig)

res<-as.data.frame(resid(burnsig))%>%
  bind_cols(burn)%>%
  mutate(residual=resid(burnsig))


m3<-lmer(residual~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(res, drtint=="2012"))
anova(m3)

m4<-lmer(residual~as.factor(calendar_year)*nitro*phos + (1|plotnum), data=subset(res, drtint=="2018"))
anova(m4)

##compare drought years to long term non-burned years only. the 2012 and 2018 drought reduced biomass. 2018 makes more sense. Not sure how to do stats on this. the recovery was not there in 2013 but did recover in 2018
unburnmean<-burn%>%
  filter(fire==0&calendar_year!=2012|calendar_year!=2018)%>%
  group_by(nitro, phos, Trt)%>%
  summarise(mean=mean(anpp))

compareburn<-burn%>%
  filter(drt=="drt")%>%
  group_by(calendar_year, drtint, nitro, phos, Trt)%>%
  summarise(dmean=mean(anpp))%>%
  left_join(unburnmean)%>%
  mutate(resist=(dmean-mean)/mean)

ggplot(data=compareburn, aes(x=Trt, y=resist))+
  geom_bar(stat="identity")+
  facet_wrap(~drtint)+
  ggtitle("Biomass reduction")


burnmean<-burn%>%
  filter(fire==1)%>%
  group_by(nitro, phos, Trt)%>%
  summarise(mean=mean(anpp))

compareburn<-burn%>%
  filter(drt=="post")%>%
  group_by(calendar_year, drtint, nitro, phos, Trt)%>%
  summarise(pmean=mean(anpp))%>%
  left_join(burnmean)%>%
  mutate(recover=(pmean-mean)/mean)

ggplot(data=compareburn, aes(x=Trt, y=recover))+
  geom_bar(stat="identity")+
  facet_wrap(~drtint)
