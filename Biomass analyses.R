library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(nlme)



p3<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/p_cubed_biomass_AllYears.csv")%>%
  mutate(drought="y")%>%
  group_by(drought, year, plot, row, replicate)%>%
  mutate(anpp=sum(grass, forb, woody)*10)%>%
  group_by(year, plot, row, drought)%>%
  summarise(anpp=mean(anpp))%>%
  rename(calendar_year=year)%>%
  left_join(treats)%>%
  select(-row, -plot, -rep)
  
treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")

p2<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/pplots/Biomass/To use/Compiling Files in R/Biomass_2002_2019.csv")%>%
  select(-X)%>%
  filter(calendar_year<2016 & calendar_year>2009,
         treatment=="N1P0"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P3")%>%
  mutate(drought="n")%>%
  rename(plotnum=plot_id)%>%
  left_join(treats)%>%
  select(-treatment, -rep)


##analysis approach one. four way anova
biomass<-p2%>%
  bind_rows(p3)%>%
  mutate(treat=ifelse(calendar_year<2013,"Drought years","Recovery years"),
         nitro=as.factor(nitro),
         phos=as.factor(phos),
         calendar_year=as.factor(calendar_year))

hist(log(biomass$anpp))


##double check nested plots design with a fixed factor here I am nestling drought within plotnum (drought/plotnum).

m.d <- lmer(log(anpp)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
anova(m.d)

m.r <- lmer(log(anpp)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
anova(m.r)

biomave<-biomass%>%
  group_by(calendar_year, drought, treat)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=biomave, aes(x=calendar_year, y=mbio, color=drought, group=drought))+
  geom_point(size=5)+
  geom_line()+
  scale_color_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se),width=.2)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  facet_wrap(~treat, scales="free_x")+
  xlab("Year")



###logRR is not the correct approach we do not think.


###now looking at disc pasture data
dp2011<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/DiscPasture/2011_dispachmeter.csv")%>%
  group_by(Year, plot, row, type)%>%
  summarize(disc=mean(dischgt))
dp2013<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/DiscPasture/p_cubed_disc_2013.csv")%>%
  select(Year, plot, row, type, disc)
dp2014<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/DiscPasture/p_cubed_disc_2014.csv")%>%
  select(Year, plot, row, type, disc)
dp2015<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/P-cubed/DiscPasture/p_cubed_disc_2015.csv")%>%
  select(Year, plot, row, type, disc)

dp<-rbind(dp2011, dp2013, dp2014, dp2015)

dp2<-merge(dp, treats, by=c("plot", "row"))%>%
  select(-rep)%>%
  mutate(nitro=as.factor(nitro),
         phos=as.factor(phos),
         year=as.factor(Year),
         treatment=ifelse(Year<2013,"treatment","recovery"),
         sqrt.hgt=sqrt(disc),
         anpp=(1805*sqrt.hgt-2065)/10)

mdp.t <- lmer(anpp~nitro*phos*type + (1|plotnum), data=subset(dp2, treatment=="treatment"))
Anova(mdp.t)

mdp.r <- lmer(anpp~year*nitro*phos*type + (1|plotnum), data=subset(dp2, treatment=="recovery"))
Anova(mdp.r)

library(doBy)

funcs<-function(x)c(mn=mean(x),se=sd(x)/sqrt(length(x)))
tograph.t<-summaryBy(anpp~Trt+type,data=subset(dp2, treatment=="treatment"), FUN=funcs)
tograph.r<-summaryBy(anpp~Trt+type+year,data=subset(dp2, treatment=="recovery"), FUN=funcs)

ggplot(data=tograph.t, aes(x=type, y=anpp.mn, color=Trt, group=Trt))+
  geom_point(size=5)+
  geom_line()+
  scale_color_discrete(name="Fertilizer\nTreatment")+
  geom_errorbar(aes(ymax=anpp.mn+anpp.se, ymin=anpp.mn-anpp.se), width=.2)+
  scale_x_discrete(name="Water Treatment", label=c("Control","Drought"))+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))

ggplot(data=tograph.r, aes(x=year, y=anpp.mn, color=Trt, shape=type))+
  geom_point(size=5)+
  geom_line()+
  scale_color_discrete(name="Fertilizer\nTreatment")+
  geom_errorbar(aes(ymax=anpp.mn+anpp.se, ymin=anpp.mn-anpp.se), width=.2)+
  #scale_x_discrete(name="Water Treatment", label=c("Control","Drought"))+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))
