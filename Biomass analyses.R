library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)
library(emmeans)
library(car)

theme_set(theme_bw(12))

#make figure with all DP and ANPP and say bad blow in 2010/2011. We think it is questionable in those years and had edge effects. Make years continuous on x-axis and note with vertical line drought verus non-drought years.

treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")

p3<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/p_cubed_biomass_AllYears.csv")%>%
  mutate(type="drought")%>%
  group_by(type, year, plot, row, replicate)%>%
  mutate(anpp=sum(grass, forb, woody)*10)%>%
  group_by(year, plot, row, type)%>%
  summarise(anpp=mean(anpp), grass=mean(grass*10), forb=mean(forb*10))%>%
  rename(calendar_year=year)%>%
  left_join(treats)%>%
  ungroup()%>%
  select(-row, -plot, -rep)
  

p2<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/pplots/Biomass/To use/Compiling Files in R/Biomass_2002_2019_grassforb.csv")%>%
  select(-X)%>%
  filter(calendar_year<2016 & calendar_year>2009,
         treatment=="N1P0"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P3")%>%
  mutate(type="control")%>%
  select(-treatment)%>%
  left_join(treats)%>%
  select(-row, -plot, -rep, -plot_id)


##analysis approach one. four way anova
biomass<-p2%>%
  bind_rows(p3)%>%
  mutate(treat=ifelse(calendar_year<2013,"Drought years","Recovery years"),
         nitro=as.factor(nitro),
         phos=as.factor(phos))%>%
  select(calendar_year, anpp, Trt, type, plotnum)%>%
  filter(calendar_year==2011)

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
         treatment=ifelse(Year<2013,"Drought years","Recovery years"),
         sqrt.hgt=sqrt(disc),
         anpp=(1805*sqrt.hgt-2065)/10)%>%
  select(-disc, -sqrt.hgt, -plot, -row, -nitro, -phos, -year, -ForPCube)%>%
  rename(calendar_year=Year, treat=treatment, canpp=anpp)%>%
  filter(calendar_year==2011)


##correlating DP with ANPP
cor<-biomass%>%
  left_join(dp2)%>%
  mutate(Type=ifelse(type=="control", "Un-droughted", "Droughted"))

test<-cor%>%
  group_by(Type)%>%
  summarize(r=round(cor.test(anpp, canpp)$estimate,3),
            p=round(cor.test(anpp, canpp)$p.value,3))%>%
  mutate(rpvalue=paste("r = ", r, ", p = ", p))

ggplot(data=cor, aes(x=canpp, y=anpp))+
  geom_point()+
  geom_abline(intercept=0)+
  scale_y_continuous(limits=c(180,830))+
  scale_x_continuous(limits=c(180,830))+
  ylab("Clipped ANPP")+
  xlab("Disc Pasture ANPP")+
  facet_wrap(~Type)+
  geom_text(data=test, mapping=aes(x=800, y = 250, label = rpvalue), hjust=1.05, vjust=1.5)

#for ESA

ggplot(data=subset(biomasstoplot, type=="Clipping"), aes(x=calendar_year, y=mbio, color=Trt, shape=drought))+
  geom_point(aes(group=Trt), size=5, position=position_dodge(0.5))+
  scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se, group=Trt),position=position_dodge(0.5), width=.1)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  xlab("Year")+
  geom_vline(xintercept = 3.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~Trt, scales="free_y")+
  geom_line(aes(group=drought))


+
  annotate("text", x=3, y=700, label="*", size=8)+
  annotate("text", x=4, y=900, label="*", size=8)+
  annotate("text", x=5, y=900, label="*", size=8)


###for grasses
m.dg <- lmer(log(grass)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
summary(m.dg)
anova(m.dg, ddf="Kenward-Roger")
emmeans(m.dg, pairwise~drought|as.factor(calendar_year), adjust="holm")

mdg<-anova(m.dg)
p.mdg<-mdg$`Pr(>F)`

m.rg <- lmer(log(grass)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
anova(m.rg, ddf="Kenward-Roger")

biomaveg<-biomass%>%
  group_by(calendar_year, drought, treat, Trt)%>%
  summarize(mbio=mean(grass),
            sd=sd(anpp),
            n=length(grass))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Grasses")


hist(biomass$forb)
###for forbs
m.df <- lmer(log(forb+1)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
anova(m.df, ddf="Kenward-Roger")

m.rf <- lmer(log(forb+1)~Trt*drought*calendar_year + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
anova(m.rf, ddf="Kenward-Roger")
summary(m.rf)
emmeans(m.rf, pairwise~drought|Trt, adjust="holm")

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
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), position=position_dodge(0.5), width=.2)+
  ylab(expression(paste("Biomass (g ","m"^"-2",")")))+
  facet_grid(type~Trt, scales="free")+
  geom_vline(xintercept=3.5)+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  geom_line(aes(group=drought))
  

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
  ylab(expression(paste("Forb ANPP (g ","m"^"-2",")")))+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


######



###resist/recovery framework is not the correct approach we do not think.
biomassdiff<-biomass%>%
  select(-grass, -forb)%>%
  spread(drought, anpp)%>%
  mutate(logrr=log(y/n),
         pd=(n-y)/y)%>%
  #filter(calendar_year!=2010)%>%
  na.omit()%>%
  mutate(trtyr=paste(calendar_year, Trt, sep="_"))

hist(biomassdiff$logrr)

m.d <- lmer(logrr~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(biomassdiff, treat=="Drought years"))
summary(m.d)
anova(m.d)

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






###linking biomass to community compositon
##need to run lf code in community_analyses

diff<-biomassdiff%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  left_join(lf)%>%
  filter(trait_cat!="NA", precip=="drought")

ggplot(data=diff, aes(x=relcov, y=logrr))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trait_cat, scales='free')

with(subset(diff, trait_cat=="Annual Forb"), cor.test(logrr, relcov))#notsig
with(subset(diff, trait_cat=="Annual Gram."), cor.test(logrr, relcov))#not sig
with(subset(diff, trait_cat=="C3 Gram."), cor.test(logrr, relcov))#not sig
with(subset(diff, trait_cat=="C4 Gram."), cor.test(logrr, relcov))#sig - but will go away
with(subset(diff, trait_cat=="N-Fixing Forb"), cor.test(logrr, relcov))#not sig
with(subset(diff, trait_cat=="Non-N-Fixing Forb"), cor.test(logrr, relcov))#sig


diff2<-biomassdiff%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  left_join(lf2)%>%
  filter(trait_cat!="NA", precip=="drought", treat=="Drought years")

ggplot(data=diff2, aes(x=relcov, y=pd))+
  geom_point(aes(color=Trt),size=3)+
  geom_smooth(data=subset(diff2, trait_cat=="Forb"), method="lm", se=F, color="black")+
  geom_smooth(data=subset(diff2, trait_cat=="Grass"), method="lm", se=F, color="black")+
  facet_wrap(~trait_cat, scales='free')+
  ylab("% Diff Control-Drought ANPP")+
  xlab("Relative Cover of Functional Type")+
  scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

with(subset(diff2, trait_cat=="Annual"), cor.test(pd, relcov))#notsig
with(subset(diff2, trait_cat=="Forb"), cor.test(pd, relcov))#sig
with(subset(diff2, trait_cat=="Grass"), cor.test(pd, relcov))#sig


diff3<-biomassdiff%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  left_join(lf2)%>%
  filter(trait_cat!="NA", precip=="drought", treat=="Recovery years")

ggplot(data=diff3, aes(x=relcov, y=pd))+
  geom_point(aes(color=Trt),size=3)+
  facet_wrap(~trait_cat, scales='free')+
  ylab("% Diff Control-Drought ANPP")+
  xlab("Relative Cover of Functional Type")+
  scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

with(subset(diff3, trait_cat=="Annual"), cor.test(pd, relcov))#notsig
with(subset(diff3, trait_cat=="Forb"), cor.test(pd, relcov))#notsig
with(subset(diff3, trait_cat=="Grass"), cor.test(pd, relcov))#notsig
