library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)
library(emmeans)
library(car)
library(gridExtra)

theme_set(theme_bw(12))

#make figure with all DP and ANPP and say bad blow in 2010/2011. We think it is questionable in those years and had edge effects. Make years continuous on x-axis and note with vertical line drought verus non-drought years.

###revised thinking, the ANPP data is just no good and relegate to the appendix. 

treats<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")

# p3<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/Analyses/July 2015 Analyses/p_cubed_biomass_AllYears.csv")%>%
#   mutate(drought="y")%>%
#   group_by(drought, year, plot, row, replicate)%>%
#   mutate(anpp=sum(grass, forb, woody)*10)%>%
#   group_by(year, plot, row, drought)%>%
#   summarise(anpp=mean(anpp), grass=mean(grass*10), forb=mean(forb*10))%>%
#   rename(calendar_year=year)%>%
#   left_join(treats)%>%
#   ungroup()%>%
#   select(-row, -plot, -rep)
#   
# 
# p2<-read.csv("C:/Users/mavolio2/Dropbox/Konza Research/pplots/Biomass/To use/Compiling Files in R/Biomass_2002_2019_grassforb.csv")%>%
#   select(-X)%>%
#   filter(calendar_year<2016 & calendar_year>2009,
#          treatment=="N1P0"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P3")%>%
#   mutate(drought="n")%>%
#   select(-treatment)%>%
#   left_join(treats)%>%
#   select(-row, -plot, -rep, -plot_id)


# ##analysis approach one. four way anova
# # biomass<-p2%>%
# #   bind_rows(p3)%>%
# #   mutate(treat=ifelse(calendar_year<2013,"Drought years","Recovery years"),
# #          nitro=as.factor(nitro),
# #          phos=as.factor(phos),
# #          calendar_year=as.factor(calendar_year))
# # 
# # hist(log(biomass$anpp))
# # hist(log(biomass$grass))
# # hist(log(biomass$forb))
# 
# 
# ##double check nested plots design with a fixed factor here I am nestling drought within plotnum
# 
# # m.d <- lmer(log(anpp)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
# # summary(m.d)
# # anova(m.d, ddf="Kenward-Roger")
# # md<-anova(m.d)
# # 
# # #doing contrasts - not doing these as there are no interactions i am interested in
# # emmeans(m.d, pairwise~drought|as.factor(calendar_year), adjust="holm")
# # 
# # m.r <- lmer(log(anpp)~Trt*drought*as.factor(calendar_year) + (1|plotnum/drought), data=subset(biomass, treat=="Recovery years"))
# # summary(m.r)
# # anova(m.r, ddf="Kenward-Roger")
# # 
# # emmeans(m.r, pairwise~drought|as.factor(calendar_year), adjust="holm")
# # 
# # biomave<-biomass%>%
# #   group_by(calendar_year, drought, treat, Trt)%>%
# #   summarize(mbio=mean(anpp),
# #             sd=sd(anpp),
# #             n=length(anpp))%>%
# #   mutate(se=sd/sqrt(n))%>%
# #   mutate(type="Clipping")
# 
# dpave<-dp2%>%
#   rename(calendar_year=year, treat=treatment)%>%
#   mutate(drought=ifelse(type=="control", "n", "y"))%>%
#   group_by(calendar_year, drought, treat, Trt)%>%
#   summarize(mbio=mean(anpp),
#             sd=sd(anpp),
#             n=length(anpp))%>%
#   mutate(se=sd/sqrt(n))%>%
#   mutate(type="Disc Pasture")
#   
# biomasstoplot<-biomave%>%
#   bind_rows(dpave)%>%
#   mutate(group=paste(calendar_year, Trt, sep=""))
# 
# ggplot(data=biomasstoplot, aes(x=calendar_year, y=mbio, color=Trt, shape=drought))+
#   geom_point(aes(group=Trt), size=5, position=position_dodge(0.5))+
#   #geom_line(aes(group=group))+
#   scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
#   scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
#   geom_errorbar(aes(ymin=mbio-(1.96*se), ymax=mbio+(1.96*se), group=Trt),position=position_dodge(0.5),width=.1)+
#   ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
#   facet_wrap(~type, scales="free_y", ncol = 1)+
#   xlab("Year")+
#   geom_vline(xintercept = 3.5)


#for ESA

# ggplot(data=subset(biomasstoplot, type=="Clipping"), aes(x=calendar_year, y=mbio, color=Trt, shape=drought))+
#   geom_point(aes(group=Trt), size=5, position=position_dodge(0.5))+
#   scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
#   scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
#   geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se, group=Trt),position=position_dodge(0.5), width=.1)+
#   ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
#   xlab("Year")+
#   geom_vline(xintercept = 3.5)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~Trt, scales="free_y")+
#   geom_line(aes(group=drought))
# 
# 
# +
#   annotate("text", x=3, y=700, label="*", size=8)+
#   annotate("text", x=4, y=900, label="*", size=8)+
#   annotate("text", x=5, y=900, label="*", size=8)


# ###for grasses
# m.dg <- lmer(log(grass)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
# summary(m.dg)
# anova(m.dg, ddf="Kenward-Roger")
# emmeans(m.dg, pairwise~drought|as.factor(calendar_year), adjust="holm")
# 
# mdg<-anova(m.dg)
# p.mdg<-mdg$`Pr(>F)`
# 
# m.rg <- lmer(log(grass)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
# anova(m.rg, ddf="Kenward-Roger")
# 
# biomaveg<-biomass%>%
#   group_by(calendar_year, drought, treat, Trt)%>%
#   summarize(mbio=mean(grass),
#             sd=sd(anpp),
#             n=length(grass))%>%
#   mutate(se=sd/sqrt(n))%>%
#   mutate(type="Grasses")
# 
# 
# hist(biomass$forb)
# ###for forbs
# m.df <- lmer(log(forb+1)~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(biomass, treat=="Drought years"))
# anova(m.df, ddf="Kenward-Roger")
# 
# m.rf <- lmer(log(forb+1)~Trt*drought*calendar_year + (1|plotnum), data=subset(biomass, treat=="Recovery years"))
# anova(m.rf, ddf="Kenward-Roger")
# summary(m.rf)
# emmeans(m.rf, pairwise~drought|Trt, adjust="holm")
# 
# biomavetype<-biomass%>%
#   group_by(calendar_year, drought, treat, Trt)%>%
#   summarize(mbio=mean(forb),
#             sd=sd(anpp),
#             n=length(anpp))%>%
#   mutate(se=sd/sqrt(n))%>%
#   mutate(type="Forbs")%>%
#   bind_rows(biomaveg)
# 
# ggplot(data=biomavetype, aes(x=calendar_year, y=mbio, color=Trt, shape=drought, group=Trt))+
#   geom_point(size=5, position = position_dodge(0.5))+
#   scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
#   scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
#   geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), position=position_dodge(0.5), width=.2)+
#   ylab(expression(paste("Biomass (g ","m"^"-2",")")))+
#   facet_grid(type~Trt, scales="free")+
#   geom_vline(xintercept=3.5)+
#   xlab("Year")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
#   geom_line(aes(group=drought))
#   
# 
# biomavef<-biomass%>%
#   filter(treat=="Recovery years")%>%
#   group_by(drought, treat, Trt)%>%
#   summarize(mbio=mean(forb),
#             sd=sd(forb),
#             n=length(forb))%>%
#   mutate(se=sd/sqrt(n))
# 
# ggplot(data=biomavef, aes(x=Trt, y=mbio, fill=drought))+
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
#   geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se),position=position_dodge(0.9),width=.2)+
#   ylab(expression(paste("Forb ANPP (g ","m"^"-2",")")))+
#   xlab("Year")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ######
# 
# 
# 
# ###resist/recovery framework is not the correct approach we do not think.
# biomassdiff<-biomass%>%
#   select(-grass, -forb)%>%
#   spread(drought, anpp)%>%
#   mutate(logrr=log(y/n),
#          pd=(n-y)/y)%>%
#   #filter(calendar_year!=2010)%>%
#   na.omit()%>%
#   mutate(trtyr=paste(calendar_year, Trt, sep="_"))
# 
# hist(biomassdiff$logrr)
# 
# m.d <- lmer(logrr~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(biomassdiff, treat=="Drought years"))
# summary(m.d)
# anova(m.d)
# 
# summary(aov(logrr~nitro*phos, data=subset(biomassdiff, calendar_year==2012)))
# 
# m.d <- lmer(logrr~nitro*phos*as.factor(calendar_year) + (1|plotnum), data=subset(biomassdiff, treat=="Recovery years"))
# summary(m.d)
# 
# 
# 
# 
# 
# ltrtyr<-unique(biomassdiff$trtyr)
# 
# ttests<-data.frame()
# 
# for (i in 1:length(ltrtyr)){
#   
#   subset<-biomassdiff%>%
#     filter(trtyr==ltrtyr[i])
#   
#   tt<-t.test(subset$logrr, mu=0)
#   
#   output<-data.frame(
#     Trt=unique(subset$Trt),
#     calendar_year=unique(subset$calendar_year),
#     treat=unique(subset$treat),
#     p.val=tt$p.value
#   )
# 
#   ttests<-rbind(ttests, output)
#   
# }
# 
# ttests2<-ttests%>%
#   mutate(adj=p.adjust(ttests$p.val, method="BH"))%>%
#   group_by(treat)%>%
#   mutate(sig=ifelse(p.val<0.05, "*",""))
# 
# bdave<-biomassdiff%>%
#   group_by(calendar_year, treat, Trt)%>%
#   summarise(mean=mean(logrr), sd=sd(logrr), n=length(logrr))%>%
#   mutate(se=sd/sqrt(n))%>%
#   right_join(ttests2)
# 
# ggplot(data=bdave, aes(x=calendar_year, y=mean, color=Trt, label=sig, group=Trt))+
#   geom_point(size=5)+
#   geom_line()+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
#   facet_wrap(~treat, scales = "free")+
#     geom_hline(yintercept = 0)

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
         anpp=(1805*sqrt.hgt-2065)/10)


c<-mean(subset(dp2, year==2011&type=="control")$anpp)
d<-mean(subset(dp2, year==2011&type=="drought")$anpp)

(d-c)/c


hist(log(dp2$anpp))

mdp.t <- lmer(log(anpp)~Trt*type + (1|plotnum), data=subset(dp2, treatment=="Drought years"))
summary(mdp.t)
anova(mdp.t, ddf="Kenward-Roger")

mdp.r <- lmer(anpp~as.factor(year)*Trt*type + (1|plotnum), data=subset(dp2, treatment=="Recovery years"))
anova(mdp.r, ddf="Kenward-Roger")
emmeans(mdp.r, pairwise~Trt|type, adjust="holm")
emmeans(mdp.r, pairwise~Trt, adjust="holm")

#means by year and nutrient
dpave<-dp2%>%
  rename(calendar_year=year, treat=treatment)%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought=ifelse(type=="control", "n", "y"))%>%
  group_by(calendar_year, treat, Trt)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(type="Disc Pasture")%>%
  mutate(Trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

#means by drought years/recovery years
dpave2<-dp2%>%
  rename(calendar_year=year, treat=treatment)%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought=ifelse(type=="control", "n", "y"))%>%
  group_by(drought, treat, Trt)%>%
  summarize(mbio=mean(anpp),
            sd=sd(anpp),
            n=length(anpp))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(label=ifelse(Trt=="Control"&treat=="Recovery years"&drought=="y", "*", ""))

#Figure 3
ggplot(data=dpave2, aes(x=Trt, y=mbio, fill=drought, label=label))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se), position=position_dodge(0.9), width=.2)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  xlab("Nutrient Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  facet_wrap(~treat)+
  geom_text(aes(y=(mbio+se)+0.05))

#appendix through time

ggplot(data=subset(dpave, treat=="Recovery years"), aes(x=as.factor(calendar_year), y=mbio, color=Trt, group=Trt))+
  geom_point(aes(group=Trt), size=5)+
  scale_color_manual(name="Treatment", values=c("Black", "Blue", "Red", "Purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  geom_line()+
  scale_shape_manual(name="Droughted", values=c(19, 17),labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mbio-se, ymax=mbio+se, group=Trt), width=.1)+
  ylab(expression(paste("ANPP (g ","m"^"-2",")")))+
  xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##looking at resistance, not sure of this.
# dp3<-dp2%>%
#   select(-disc, -sqrt.hgt)%>%
#   spread(type, anpp)%>%
#   mutate(resist=log(drought/control))
# 
# hist(dp3$resist)
# 
# summary(aov(resist~nitro*phos, data=subset(dp3, treatment=="treatment")))
# 
# dpave<-dp3%>%
#   filter(year==2011)%>%
#   group_by(Trt)%>%
#   summarise(mean=mean(resist), sd=sd(resist))%>%
#   mutate(se=sd/sqrt(6))
# 
# ggplot(data=dpave, aes(x=Trt, y=mean))+
#   geom_bar(stat="identity")+
#   geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=.2)+
#   ylab("Diff in ANPP")

###doing the ANPP as differences

dpdiff<-dp2%>%
    select(-disc, -sqrt.hgt)%>%
    spread(type, anpp)%>%
    mutate(diff=drought-control)%>%
  rename(calendar_year=Year)%>%
  mutate(plotnum=as.character(plotnum))%>%
  select(calendar_year, plotnum, Trt, treatment, diff)

dpdiff.t <- aov(diff~Trt, data=subset(dpdiff, treatment=="Drought years"))
summary(dpdiff.t)


dpdiff.r <- lmer(diff~Trt*as.factor(calendar_year) + (1|plotnum), data=subset(dpdiff, treatment=="Recovery years"))
anova(dpdiff.r, ddf="Kenward-Roger")


dpdiffmeans<-dpdiff%>%
  group_by(Trt, treatment)%>%
  summarize(mean=mean(diff), sd=sd(diff), n=length(diff))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

ggplot(dpdiffmeans, aes(x=trt2, y=mean, fill=Trt))+
  geom_bar(stat="identity", position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width=0.1)+
  facet_wrap(~treatment, ncol=1)+
  scale_fill_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"),  values=c("black", "blue", "red", "purple"))+
  scale_x_discrete(labels=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept = 0)+
  ylab("ANPP differences due to drought")+
  xlab("Nutrient Treatment")+
  theme(legend.position = "none")

##by treatment and time - can put this in an appendix
dpdiffmeans_time<-dpdiff%>%
  group_by(Trt, treatment, calendar_year)%>%
  summarize(mean=mean(diff), sd=sd(diff), n=length(diff))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

ggplot(subset(dpdiffmeans_time, treatment=="Recovery years"), aes(x=as.factor(calendar_year), y=mean, color=trt2))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
  scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept = 0)+
  ylab("ANPP differences due to drought")+
  xlab("Year")+
  geom_line(aes(group=trt2))

###linking biomass to community compositon
##need to run lf code in community_analyses
#can I think of otherways to do this. Maybe chat with Kim.
##also need to get multdiff and multdiff_func

#linking cover composition of plots with control-treat production differences

#looking at all compositonal and fucntional (life forms) differences

##here I am correlating differences in production to differences in composition based on species and composition based on functional types. THere is nothing here. 
diff_comp<-dpdiff%>%
  left_join(multdiff)%>%
  mutate(type="species")

diff_func<-dpdiff%>%
  left_join(multdiff_func)%>%
  mutate(type="functional")%>%
  bind_rows(diff_comp)


with(subset(diff_comp, treatment=="recovery"&type=="species"), cor.test(diff, composition_diff))
with(subset(diff_comp, treatment=="treatment"&type=="species"), cor.test(diff, composition_diff))
with(subset(diff_func, treatment=="recovery"&type=="functional"), cor.test(diff, composition_diff))
with(subset(diff_func, treatment=="treatment"&type=="functional"), cor.test(diff, composition_diff))

ggplot(data=diff_func, aes(x=composition_diff, y=diff))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(type~treatment)+
  geom_hline(yintercept=0)



#tried doing this with the diff in production and the cover of LF in the treated plots. No longer want to do it this way. 

# diff<-dpdiff%>%
#   rename(calendar_year=Year)%>%
#   mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
#   left_join(lf)%>%
#   filter(trait_cat!="NA", precip=="drought")
# 
# ggplot(data=diff, aes(x=relcov, y=diff))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   facet_grid(~trait_cat~treatment, scales='free')
# 
# with(subset(diff, trait_cat=="Annual Forb"&treatment=="recovery"), cor.test(diff, relcov))#sig
# with(subset(diff, trait_cat=="Annual Gram."), cor.test(diff, relcov))#not sig
# with(subset(diff, trait_cat=="C3 Gram."), cor.test(diff, relcov))#sig
# with(subset(diff, trait_cat=="C4 Gram."), cor.test(diff, relcov))#not sig - but will go away
# with(subset(diff, trait_cat=="N-Fixing Forb"), cor.test(diff, relcov))#not sig
# with(subset(diff, trait_cat=="Non-N-Fixing Forb"), cor.test(diff, relcov))#sig
# 
# 
# diff2<-dpdiff%>%
#   rename(calendar_year=Year)%>%
#   mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
#   left_join(lf2)%>%
#   filter(trait_cat!="NA", precip=="drought", treatment=="treatment")
# 
# ggplot(data=diff2, aes(x=relcov, y=diff))+
#   geom_point(aes(color=Trt),size=3)+
#   #geom_smooth(data=subset(diff2, trait_cat=="Forb"), method="lm", se=F, color="black")+
#   #geom_smooth(data=subset(diff2, trait_cat=="Grass"), method="lm", se=F, color="black")+
#   facet_wrap(~trait_cat, scales='free')+
#   ylab("% ANPP Difference\nDroughted-Non-Droughted")+
#   xlab("Relative Cover of Functional Type")+
#   scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# with(subset(diff2, trait_cat=="Annual"), cor.test(diff, relcov))#notsig
# with(subset(diff2, trait_cat=="Forb"), cor.test(diff, relcov))#notsig
# with(subset(diff2, trait_cat=="Grass"), cor.test(diff, relcov))#notsig
# 
# 
# diff3<-dpdiff%>%
#   rename(calendar_year=Year)%>%
#   mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
#   left_join(lf2)%>%
#   filter(trait_cat!="NA", precip=="drought", treatment=="recovery")
# 
# ggplot(data=diff3, aes(x=relcov, y=diff))+
#   geom_point(aes(color=Trt),size=3)+
#   facet_wrap(~trait_cat, scales='free')+
#   ylab("% ANPP Difference\nDroughted-Non-Droughted")+
#   geom_smooth(data=subset(diff3, trait_cat=="Annual"), method="lm", se=F, color="black")+
#   geom_smooth(data=subset(diff3, trait_cat=="Grass"), method="lm", se=F, color="black")+
#   xlab("Relative Cover of Functional Type")+
#   scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# with(subset(diff3, trait_cat=="Annual"), cor.test(diff, relcov))#sig
# with(subset(diff3, trait_cat=="Forb"), cor.test(diff, relcov))#notsig
# with(subset(diff3, trait_cat=="Grass"), cor.test(diff, relcov))#sig
