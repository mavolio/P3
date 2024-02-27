library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)
library(emmeans)
library(car)

#test

theme_set(theme_bw(12))

#make figure with all DP and ANPP and say bad blow in 2010/2011. We think it is questionable in those years and had edge effects. Make years continuous on x-axis and note with vertical line drought verus non-drought years.

treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")


#Clipped biomass data
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

biomass<-p2%>%
  bind_rows(p3)%>%
  mutate(treat=ifelse(calendar_year<2013,"Drought years","Recovery years"),
         nitro=as.factor(nitro),
         phos=as.factor(phos))%>%
  select(calendar_year, grass, forb, anpp, Trt, type, treat, plotnum, nitro,phos) %>% 
  rename(drought=type) %>% 
  mutate(year=as.factor(calendar_year)) %>% 
  mutate(plot=as.factor(paste("p", plotnum, sep="-")))

write.csv(biomass, "C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/BiomassforSAS.csv", row.names=F)

##Disc pasture data data
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
  rename(calendar_year=Year, treat=treatment, canpp=anpp, drought=type) %>% 
  mutate(year=as.factor(calendar_year)) %>% 
  mutate(plot=as.factor(paste("p", plotnum, sep="-"))) %>% 
  mutate(ploid=paste(drought, plotnum, sep="_"))

write.csv(dp2, "C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/StandingBiomassforSAS.csv", row.names=F)

##correlating DP with ANPP
cor<-biomass%>%
  left_join(dp2)%>%
  mutate(Type=ifelse(drought=="control", "Un-droughted", "Droughted")) %>%
  filter(calendar_year==2011)

test<-cor%>%
  group_by(Type)%>%
  summarize(r=round(cor.test(anpp, canpp)$estimate,3),
            p=round(cor.test(anpp, canpp)$p.value, 3))%>%
  mutate(rpvalue=paste("r = ", r, ", p = ", p))

ggplot(data=cor, aes(x=canpp, y=anpp))+
  geom_point()+
  geom_abline(intercept=0)+
  scale_y_continuous(limits=c(180,830))+
  scale_x_continuous(limits=c(180,830))+
  ylab(expression(paste("Clipped ANPP (g ","m"^"-2",")")))+
  xlab("Disc Pasture, Standing biomass (g)")+
  facet_wrap(~Type)+
  geom_text(data=test, mapping=aes(x=800, y = 250, label = rpvalue), hjust=1.05, vjust=1.5)

###Running a two-way repated measures anova with Nutrient treatment, drought treatment, and year as fixes effects

#total biomass ANPP
hist(log(biomass$anpp))

m.drt<- lmer(log(anpp)~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(biomass, treat=="Drought years"))
summary(m.drt)
anova(m.drt, ddf="Kenward-Roger")
emmeans(m.drt, pairwise~drought, adjust="holm")


m.rec <- lmer(log(anpp)~Trt*drought*year + (1|plot)+ (1|plot:drought)+(1|plot:year), data=subset(biomass, treat=="Recovery years"))
anova(m.rec, ddf="Kenward-Roger")

emmeans(m.rec, pairwise~drought|Trt, adjust="holm")

biomassave<-biomass %>% 
  group_by(drought, treat) %>% 
  summarize(manpp=mean(anpp))

#disc pasture as measure of standing biomass
hist(log(dp2a$canpp))
#one way anova in drought year
m.drt_dp <- lmer(log(canpp)~Trt*drought + (1|plot), data=subset(dp2, treat=="Drought years"))
anova(m.drt_dp, ddf="Kenward-Roger")
emmeans(m.drt_dp, pairwise~drought|Trt, adjust="holm")

m.rec_dp <- lmer(log(canpp)~Trt*drought*year + (1|plot)+(1|plot:drought) + (1|plot:year), data=subset(dp2, treat=="Recovery years"))
anova(m.rec_dp, ddf="Kenward-Roger")

emmeans(m.rec_dp, pairwise~Trt, adjust="holm")

dpave<-dp2 %>% 
  group_by(drought, treat) %>% 
  summarize(manpp=mean(canpp))


#Figure of standing biomass
dpave2<-dp2%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought2=ifelse(drought=="control", "n", "y"))%>%
  group_by(calendar_year, drought2, treat, Trt)%>%
  summarize(mbio=mean(canpp),
            sd=sd(canpp),
            n=length(canpp))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(group=paste(drought2, Trt, sep = "_"))
  # mutate(label=ifelse(Trt=='Control'&treat=='Recovery years'&drought2=='n', '*', ifelse(Trt=='P&N'&treat=='Recovery years'&drought2=='n','*', ifelse(treat=='Drought years'&drought2=='n', '*', ''))))

ggplot(data=subset(dpave2, treat=='Recovery years'), aes(x=as.factor(calendar_year), y=mbio, color=drought2, shape=Trt, group=group))+
  geom_point(size=4)+
  geom_line()+
  scale_shape_manual(name="Nutrient\nTreatment", breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"), values=c(15, 16, 17, 18))+
  scale_color_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), width=.1)+
  ylab("Standing biomass (g)")+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #scale_x_discrete(limits=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  # geom_text(aes(y=(mbio)+100))+
  facet_wrap(~Trt)

#Figure of biomass
biomassave<-biomass%>%
  mutate(drought2=ifelse(drought=="control", "n", "y"))%>%
  group_by(drought2, treat, Trt)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n)) 

ggplot(data=subset(biomassave, treat=="Recovery years"), aes(x=Trt, y=mbio, fill=drought2))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), position=position_dodge(0.9), width=.2)+
  ylab(expression(paste("Clipped ANPP (g ","m"^"-2",")")))+
  xlab("Nutrient Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))

###need to consider whether or not I can do this.
###linking biomass to community composition
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
