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

dp<-rbind(dp2011, dp2013, dp2015)

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
  mutate(ploid=paste(drought, plotnum, sep="_")) %>% 
  select(-disc, -sqrt.hgt)

# write.csv(dp2, "C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/StandingBiomassforSAS.csv", row.names=F)
write.csv(dp2, "C:\\Users\\mavolio2\\Dropbox\\Konza Research\\P-cubed\\Public data\\StandingBiomassforPublication.csv", row.names=F)


###Running a two-way repeated measures anova with Nutrient treatment, drought treatment, and year as fixes effects

hist(log(dp2$canpp))
#one way anova in drought year
m.drt_dp <- lmer(log(canpp)~Trt*drought + (1|plot), data=subset(dp2, treat=="Drought years"))
anova(m.drt_dp, ddf="Kenward-Roger")
emmeans(m.drt_dp, pairwise~Trt, adjust="tukey")

m.rec_dp <- lmer(log(canpp)~Trt*drought*year + (1|plot)+(1|plot:drought) + (1|plot:year), data=subset(dp2, treat=="Recovery years"))
anova(m.rec_dp, ddf="Kenward-Roger")

emmeans(m.rec_dp, pairwise~Trt, adjust="tukey")

#Figure of standing biomass
dpave2<-dp2%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought2=ifelse(drought=="control", "n", "y"))%>%
  group_by(drought2, treat, Trt)%>%
  summarize(mbio=mean(canpp),
            sd=sd(canpp),
            n=length(canpp))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(label=ifelse(treat=='Drought years'&Trt=="Control"&drought2=='n'|treat=='Recovery years'&Trt=="Control"&drought2=='n', "C",
               ifelse(treat=='Drought years'&Trt=='P'&drought2=='n', 'BC',
               ifelse(treat=='Drought years'&Trt=="N"&drought2=='n', 'AB',
               ifelse(treat=='Drought years'&Trt=='P&N'&drought2=='n'|treat=='Recovery years'&Trt=='P&N'&drought2=='n', 'A', 
               ifelse(treat=='Recovery years'&Trt=="P"&drought2=='n'|treat=='Recovery years'&Trt=="N"&drought2=='n', 'B', ""))))))

collabel<-c("Drought years"="Drought year (2011)", 'Recovery years' = 'Recovery years (2013, 2015)')

ggplot(data=subset(dpave2), aes(x=Trt, y=mbio, fill=drought2, label=label))+
  geom_bar(stat="identity", position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), position = position_dodge(0.9), width=0.1)+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  ylab("Standing biomass (g)")+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  facet_wrap(~treat, labeller = labeller(treat=collabel))+
  geom_text(aes(y=(mbio)+120))
