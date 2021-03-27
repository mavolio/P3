library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(nlme)


#make figure with all DP and ANPP and say bad blow in 2010/2011. We think it is questionable in those years and had edge effects. Make years continuous on x-axis and note with vertical line drought verus non-drought years.

treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")

p3<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/p_cubed_biomass_AllYears.csv")%>%
  mutate(drought="y")%>%
  group_by(drought, year, plot, row, replicate)%>%
  mutate(anpp=sum(grass, forb, woody)*10)%>%
  group_by(year, plot, row, drought)%>%
  summarise(anpp=mean(anpp), grass=mean(grass*10), forb=mean(forb*10))%>%
  rename(calendar_year=year)%>%
  left_join(treats)%>%
  ungroup()%>%
  select(-row, -plot, -rep)
  

p2<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/pplots/Biomass/To use/Compiling Files in R/Biomass_2002_2019_grassforb.csv")%>%
  select(-X)%>%
  filter(calendar_year<2016 & calendar_year>2009,
         treatment=="N1P0"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P3")%>%
  mutate(drought="n")%>%
  select(-treatment)%>%
  left_join(treats)%>%
  select(-row, -plot, -rep, -plot_id)


##analysis approach one. four way anova
biomass<-p2%>%
  bind_rows(p3)%>%
  mutate(treat=ifelse(calendar_year<2013,"Drought years","Recovery years"),
         nitro=as.factor(nitro),
         phos=as.factor(phos),
         calendar_year=as.factor(calendar_year))

hist(log(biomass$anpp))
hist(log(biomass$grass))
hist(log(biomass$forb))


##double check nested plots design with a fixed factor here I am nestling drought within plotnum (drought/plotnum).

m.d <- lmer(log(anpp)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
summary(m.d)
anova(m.d)
md<-anova(m.d)

m.r <- lmer(log(anpp)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
summary(m.r)
anova(m.r)

biomave<-biomass%>%
  group_by(calendar_year, drought, treat, Trt)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Clipping")

dpave<-dp2%>%
  rename(calendar_year=year, treat=treatment)%>%
  mutate(drought=ifelse(type=="control", "n", "y"))%>%
  group_by(calendar_year, drought, treat, Trt)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Disc Pasture")
  
biomasstoplot<-biomave%>%
  bind_rows(dpave)%>%
  mutate(group=paste(calendar_year, Trt, sep=""))

ggplot(data=biomasstoplot, aes(x=calendar_year, y=mbio, color=Trt, shape=drought))+
  geom_point(aes(group=Trt), size=5, position=position_dodge(0.5))+
  #geom_line(aes(group=group))+
  scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-(1.96*se), ymax=mbio+(1.96*se), group=Trt),position=position_dodge(0.5),width=.1)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  facet_wrap(~type, scales="free_y", ncol = 1)+
  xlab("Year")+
  geom_vline(xintercept = 3.5)


###for grasses
m.dg <- lmer(log(grass)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
summary(m.dg)
anova(m.dg)

mdg<-anova(m.dg)
p.mdg<-mdg$`Pr(>F)`

m.rg <- lmer(log(grass)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
anova(m.rg)

biomaveg<-biomass%>%
  group_by(calendar_year, drought, treat, Trt)%>%
  summarize(mbio=mean(grass),
            sd=sd(anpp),
            n=length(grass))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Grasses")


hist(biomass$forb)
###for forbs
m.df <- lmer(log(forb+1)~nitro*phos*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
anova(m.df)

m.rf <- lmer(log(forb+1)~nitro*phos*drought*calendar_year + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
anova(m.rf)
summary(m.rf)
ls_means(m.rf, pairwise = T)

biomavetype<-biomass%>%
  group_by(calendar_year, drought, treat, Trt)%>%
  summarize(mbio=mean(forb),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Forbs")%>%
  bind_rows(biomaveg)

ggplot(data=biomavetype, aes(x=calendar_year, y=mbio, color=Trt, shape=drought, group=Trt))+
  geom_point(size=5, position = position_dodge(0.5))+
  scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-(1.96*se), ymax=mbio+(1.96*se)), position=position_dodge(0.5), width=.2)+
  ylab(expression(paste("Biomass (g ","m"^"-2",")")))+
  facet_wrap(~type, scales="free", ncol=1)+
  geom_vline(xintercept=3.5)+
  xlab("Year")

biomavef<-biomass%>%
  filter(treat=="Recovery years")%>%
  group_by(drought, treat, Trt)%>%
  summarize(mbio=mean(forb),
            sd=sd(forb),
            n=length(forb))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=biomavef, aes(x=Trt, y=mbio, fill=drought))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se),position=position_dodge(0.9),width=.2)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  xlab("Year")


###resist/recovery framework is not the correct approach we do not think.
biomassdiff<-biomass%>%
  select(-grass, -forb)%>%
  spread(drought, anpp)%>%
  mutate(logrr=log(y/n))%>%
  #filter(calendar_year!=2010)%>%
  na.omit()%>%
  mutate(trtyr=paste(calendar_year, Trt, sep="_"))

hist(biomassdiff$logrr)

m.d <- lmer(logrr~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(biomassdiff, treat=="Drought years"))
summary(m.d)

summary(aov(logrr~nitro*phos, data=subset(biomassdiff, calendar_year==2012)))

m.d <- lmer(logrr~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(biomassdiff, treat=="Recovery years"))
summary(m.d)

ltrtyr<-unique(biomassdiff$trtyr)

ttests<-data.frame()

for (i in 1:length(ltrtyr)){
  
  subset<-biomassdiff%>%
    filter(trtyr==ltrtyr[i])
  
  tt<-t.test(subset$logrr, mu=0)
  
  output<-data.frame(
    Trt=unique(subset$Trt),
    calendar_year=unique(subset$calendar_year),
    treat=unique(subset$treat),
    p.val=tt$p.value
  )

  ttests<-rbind(ttests, output)
  
}

ttests2<-ttests%>%
  mutate(adj=p.adjust(ttests$p.val, method="BH"))%>%
  group_by(treat)%>%
  mutate(sig=ifelse(p.val<0.05, "*",""))

bdave<-biomassdiff%>%
  group_by(calendar_year, treat, Trt)%>%
  summarise(mean=mean(logrr), sd=sd(logrr), n=length(logrr))%>%
  mutate(se=sd/sqrt(n))%>%
  right_join(ttests2)

ggplot(data=bdave, aes(x=calendar_year, y=mean, color=Trt, label=sig, group=Trt))+
  geom_point(size=5)+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~treat, scales = "free")+
    geom_hline(yintercept = 0)




##not going to include this data
##now looking at disc pasture data
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


hist(log(dp2$anpp))

mdp.t <- lmer(log(anpp)~nitro*phos*type + (1|plotnum), data=subset(dp2, treatment=="treatment"))
summary(mdp.t)
anova(mdp.t)
mdp.r <- lmer(anpp~year*nitro*phos*type + (1|plotnum), data=subset(dp2, treatment=="recovery"))
anova(mdp.r)

dp3<-dp2%>%
  select(-disc, -sqrt.hgt)%>%
  spread(type, anpp)%>%
  mutate(resist=log(drought/control))

hist(dp3$resist)

summary(aov(resist~nitro*phos, data=subset(dp3, treatment=="treatment")))

dpave<-dp3%>%
  filter(year==2011)%>%
  group_by(Trt)%>%
  summarise(mean=mean(resist), sd=sd(resist))%>%
  mutate(se=sd/sqrt(6))

ggplot(data=dpave, aes(x=Trt, y=mean))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=.2)+
  ylab("Diff in ANPP")


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
