library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(codyn)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(gridExtra)

theme_set(theme_bw(12))

#set wd
my.wd<-setwd("C:/Users/mavolio2/Dropbox/Konza Research")
my.wd<-setwd("E:/Dropbox/Konza Research")

#read in data
p3plotcomp<-read.csv(paste(my.wd, "/p-cubed/SpComp/pcubed_sp_data2010-2015.csv", sep = ""))%>%
  select(year, plotnum, spnum, genus_species, cover, precip, life, form, C3_C4, n_fixer)%>%
  rename(calendar_year=year, abundance=cover)

treats<-read.csv(paste(my.wd, "/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv", sep=""))

comp<-p3plotcomp%>%
  select(calendar_year, plotnum, spnum, genus_species, abundance, precip)%>%
  left_join(treats)%>%
  mutate(unid=paste(plotnum, precip, sep="_"), 
         unidt=paste(Trt, precip, sep="_"),
         unid3=paste(plotnum, Trt, precip, sep="_"))

##looking at richness and evenness for structure analyses
          
richeven<-community_structure(comp, time.var="calendar_year", abundance.var="abundance", replicate.var = "unid")%>%
  separate(unid, into=c("plotnum", "drought"), sep="_")%>%
  mutate(plotnum=as.integer(plotnum))%>%
  left_join(treats)%>%
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years")) %>% 
  mutate(year=as.factor(calendar_year)) %>% 
  mutate(plot=as.factor(paste("p", plotnum, sep="-")))

hist(richeven$richness)#this is normal
hist(log(richeven$Evar))#log transfrom even

#richness in plots in 2010
rich2010<-richeven%>%
  filter(calendar_year==2010&drought=="control")%>%
  group_by(Trt)%>%
  summarize(mean=mean(richness))


##richness in drought
##this is where i left off this model might not be correct.
rich.d <- lmer(richness~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Drought years"))

#summary(rich.d)
anova(rich.d, ddf="Kenward-Roger")

#doing contrasts
emmeans(rich.d, pairwise~drought|Trt, adjust="holm")

#richness during recovery
rich.r <- lmer(richness~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Recovery years"))

anova(rich.r, ddf="Kenward-Roger")

#doing contrasts
emmeans(rich.r, pairwise~Trt, adjust="holm")

##evenness in drought
even.d <- lmer(log(Evar)~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Drought years"))
anova(even.d, ddf="Kenward-Roger")
#doing contrasts
emmeans(even.d, pairwise~Trt, adjust="holm")

#evenness during recovery
even.r <- lmer(log(Evar)~Trt*drought*year +(1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Recovery years"))
anova(even.r, ddf="Kenward-Roger")

richtoplot<-richeven%>%
  group_by(drought, Trt, treat)%>%
  summarise(mean=mean(richness), sd=sd(richness), n=length(richness))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

eventoplot<-richeven%>%
  group_by(drought, treat)%>%
  summarise(mean=mean(Evar), sd=sd(Evar), n=length(Evar))%>%
  mutate(se=sd/sqrt(n))

trtlab<-c(Control="Control", P="P", N="N", 'P&N'="N+P")

ggplot(richtoplot, aes(x=trt2, y=mean, fill=drought))+
  geom_bar(stat='identity', position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position = position_dodge(0.9))+
  facet_wrap(~trt2, labeller =labeller(trt2=trtlab))+
  geom_line()+
  #scale_color_manual(name="Treatment", values=c("black", "blue", "red", "purple"), breaks = c("Control", "P", "N","P&N"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+
  ylab("Species Richness")+
  xlab("Nutrient Treatment")+
  facet_wrap(~treat)


#Given the experimental design, looking at change does not seem like the best approach because there is not pre-treatment data.

# ###looking at community changes
# deltarac<-RAC_change(comp, time.var = "calendar_year", abundance.var = "abundance", replicate.var = "unid", species.var = "genus_species")%>%
#   separate(unid, into=c("plotnum", "drought"), sep="_")%>%
#   mutate(plotnum=as.integer(plotnum))%>%
#   left_join(treats)%>%
#   mutate(treat=ifelse(calendar_year<2013, "drought", "recovery"))

# ##multivariate change
# deltamult<-multivariate_change(comp, time.var="calendar_year", treatment.var = "unidt", replicate.var = "unid", species.var = "genus_species", abundance.var = "abundance")%>%
#   separate(unidt, into=c("Trt", "drought"), sep="_")
#
# ggplot(data=deltamult, aes(x=calendar_year2, y=composition_change, color=Trt, shape=drought))+
#   geom_point()+
#   facet_wrap(~Trt)
#
# ##multivariate change
# deltamult<-multivariate_change(comp, time.var="calendar_year", treatment.var = "unidt", replicate.var = "unid", species.var = "genus_species", abundance.var = "abundance")%>%
#   separate(unidt, into=c("Trt", "drought"), sep="_")
#
# ggplot(data=deltamult, aes(x=calendar_year2, y=composition_change, color=Trt, shape=drought))+
#   geom_point()+
#   facet_wrap(~Trt)

##mult diff - making this so that it is the distance between drought and control plots for each plot for each year. then I can do stats on this b/c plots are truly paired.

#dropping unknown species
comp2<-comp%>%
  filter(spnum<900)
  
multdiff<-multivariate_difference(comp2, time.var="calendar_year", species.var="genus_species", abundance.var="abundance", replicate.var = "unid", treatment.var = "unid3")%>%
  separate(unid3, into=c("plotnum", "Trt", "drought"), sep="_")%>%
  separate(unid32, into=c("plotnum2", "Trt2", "drought2"), sep="_")%>%
  filter(plotnum==plotnum2)%>%
  select(-Trt2, -plotnum2)%>%
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years"))%>% 
  mutate(year=as.factor(calendar_year)) %>% 
  mutate(plot=as.factor(paste("p", plotnum, sep="-")))


##

mult.d <- lmer(composition_diff~Trt*year + (1|plot), data=subset(multdiff, treat=="Drought years"))
anova(mult.d, ddf="Kenward-Roger")


mult.r <- lmer(composition_diff~Trt*year + (1|plot), data=subset(multdiff, treat=="Recovery years"))
anova(mult.r, ddf="Kenward-Roger")
emmeans(mult.r, pairwise~Trt, adjust="holm")

multdiff_means2<-multdiff%>%
  group_by(Trt, treat)%>%
  summarize(mean=mean(composition_diff), sd=sd(composition_diff))%>%
  mutate(se=sd/sqrt(6))%>%
  mutate(label=ifelse(treat=="Drought years", "A", ifelse(treat=="Recovery years"&Trt=="Control", "B", ifelse(treat=="Recovery years"&Trt=="P&N", "A", "AB"))))
  
dist<-ggplot(data=multdiff_means2, aes(x=Trt, y=mean, label=label))+
  geom_bar(stat='identity', size=3, fill="gray")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Distance between\nCtrl-Drt Centroids")+
  xlab("Nutirient Treatment")+
  facet_wrap(~treat)+
  geom_text(aes(y=mean+se+0.05))+
  ggtitle("B)")+
  theme(plot.title = element_text(hjust=-0.1, vjust = -8))

###doing NMDS. taking the average cover of plot over all three years during drought and control
ave<-comp%>%
  filter(spnum<900) %>% 
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years" )) %>% 
  group_by(calendar_year, plotnum, treat, Trt, precip, genus_species)%>%
  summarize(mabund=mean(abundance))

drought<-ave %>% 
  filter(treat=="Drought years")

#make wide for vegan
dwide<-drought %>% 
  pivot_wider(names_from = genus_species, values_from = mabund, values_fill = 0)

#pull out plot info
dplots<-dwide[,1:5]

#run nmds for drought years
dmds<-metaMDS(dwide[,6:75], trymax = 1000, autotransform = F)
dmds

droughtscores<-dplots%>%
  bind_cols(as.data.frame(dmds$points))

d2010<-dwide %>% 
  filter(calendar_year==2010)
d2011<-dwide %>% 
  filter(calendar_year==2011)
d2012<-dwide %>% 
  filter(calendar_year==2012)

adonis2(d2010[,6:75]~Trt*precip, data=d2010)
adonis2(d2011[,6:75]~Trt*precip, data=d2011)
adonis2(d2012[,6:75]~Trt*precip, data=d2012)
#recovery
recovery<-ave %>% 
  filter(treat=="Recovery years")

rwide<-recovery %>% 
  pivot_wider(names_from = genus_species, values_from = mabund, values_fill = 0)

#pull out plot info
rplots<-rwide[,1:5]

#run nmds for drought years
rmds<-metaMDS(rwide[,6:69], trymax = 3000, autotransform = F)
rmds

r2013<-rwide %>% 
  filter(calendar_year==2013)
r2014<-rwide %>% 
  filter(calendar_year==2014)
r2015<-rwide %>% 
  filter(calendar_year==2015)

adonis2(r2013[,6:69]~Trt*precip, data=r2013)
adonis2(r2014[,6:69]~Trt*precip, data=r2014)
adonis2(r2015[,6:69]~Trt*precip, data=r2015)

recoveryscores<-rplots%>%
  bind_cols(as.data.frame(rmds$points)) %>% 
  bind_rows(droughtscores) %>% 
  group_by(treat, calendar_year, Trt, precip) %>% 
  summarize(NMDS1=mean(MDS1), NMDS2=mean(MDS2), sd1=sd(MDS1), sd2=sd(MDS2)) %>% 
  mutate(se1=sd1/sqrt(6), se2=sd2/sqrt(6),
         Droughted=ifelse(precip=='drought', "Yes", "No"),
         group=paste(Trt, calendar_year, sep="::"),
         yr=substr(calendar_year, 3,4))

nmds<-ggplot(data=recoveryscores, aes(x=NMDS1, y=NMDS2, color=Droughted, shape=Trt, group=group))+
  geom_errorbar(aes(ymin=NMDS2-se2, ymax=NMDS2+se2), width=0.01)+
  geom_errorbarh(aes(xmin=NMDS1-se1, xmax=NMDS1+se1), height=0.01)+
  scale_shape_manual(name="Nutrient\nTreatment", breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"), values=c(21, 22, 23, 24))+
  #scale_fill_manual(values="white")+
  scale_color_manual(name="Droughted", values=c("Blue", "Orange"))+
  geom_point(size=4.5, stroke=1.5, fill="white")+
  geom_path(color="black")+
  geom_text(aes(label=yr), color="black", size=3)+
  scale_x_continuous(limits=c(-0.7, 0.8))+
  scale_y_continuous(limits = c(-0.7, 0.8))+
  facet_wrap(~treat)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("A)")+
  theme(plot.title = element_text(hjust=-0.1, vjust = -8))
  
grid.arrange(nmds, dist, ncol=1)



##life form
totcov<-comp%>%
  group_by(plotnum, calendar_year)%>%
  summarize(tot=sum(abundance))

lf<-p3plotcomp%>%
  filter(spnum<990) %>% #omits unknowns
mutate(trait_cat=ifelse(form=="F"&life=="A", "Annual Forb",
                 ifelse(form=="G"&life=="A", "Annual Grass",
                 ifelse(form=="G"&C3_C4=="C3", "C3 Gram.",
                 ifelse(form=="G"&C3_C4=="C4", "C4 Grass",
                 ifelse(form=='F'|form=="S"&n_fixer=="N", "Non-N-Fixing Forb",
                 ifelse(form=='F'|form=="S"&n_fixer=="Y", "N-Fixing Forb","UNK")))))))%>%
  left_join(treats)%>%
  group_by(plotnum, Trt, precip, calendar_year, trait_cat)%>%
  summarise(cov=sum(abundance))%>%
  #left_join(totcov)%>%
  #mutate(relcov=cov/tot)%>%
  mutate(unid=paste(plotnum, precip, sep="_"), 
         unidt=paste(Trt, precip, sep="_"),
         unid3=paste(plotnum, Trt, precip, sep="_"))

check<-lf %>% 
  group_by(Trt, precip, trait_cat) %>% 
  summarize(mean=mean(cov))

lf_stat<-lf%>% 
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years")) %>% 
  select(plotnum, Trt, precip, calendar_year, treat, trait_cat, cov) %>% 
  pivot_wider(names_from="trait_cat", values_from="cov", values_fill=0) %>% 
  pivot_longer("Annual Forb":"Annual Grass", names_to="trait_cat", values_to = "cov") %>% 
  mutate(year=as.factor(calendar_year)) %>% 
  mutate(plot=as.factor(paste("p", plotnum, sep="-")))


###looking at cover differences of lf through time
lfcov.d <- lmer(cov~Trt*precip*trait_cat*year + (1|plot:precip:year) + (1|plot:precip) + (1|plot:year)+ (1|plot:trait_cat) + (1|plot:year:trait_cat) + (trait_cat), data=subset(lf_stat, treat=="Drought years"))
anova(lfcov.d, ddf="Kenward-Roger")

(1 | plot / drought / year) +
  (1 | plot:year) +
  (1 | plot:drought:PFT) +
  (1 | plot:PFT) +
  (1 | plot:year:PFT),


#doing contrasts
emmeans(lfcov.d, pairwise~precip|Trt|trait_cat, adjust="holm")
emmeans(lfcov.d, pairwise~precip|trait_cat, adjust="holm")

lfcov.r <- lmer(cov~Trt*precip*trait_cat*year + (1|plot) + (1|plot:precip) + (1|plot:year), data=subset(lf_stat, treat=="Recovery years"))
anova(lfcov.r, ddf="Kenward-Roger")
emmeans(lfcov.r, pairwise~precip|Trt|trait_cat, adjust="holm")
emmeans(lfcov.r, pairwise~precip|trait_cat, adjust="holm")

lfave<-lf_stat%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought=ifelse(precip=="control", "n", "y"))%>%
  group_by(drought, treat, Trt, trait_cat)%>%
  summarize(mcov=mean(cov),
            sd=sd(cov),
            n=length(cov))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(label=ifelse(Trt=="N"&treat=="Drought years"&trait_cat=="C3 Gram."&drought=='n','*', ifelse(Trt=='P'&treat=='Drought years'&trait_cat=='C3 Gram.'&drought=='n', '*', ifelse(Trt=='Control'&treat=='Drought years'&trait_cat=='C4 Grass'&drought=='n', '*', ifelse(Trt=='N'&treat=='Drought years'&trait_cat=='C4 Grass'&drought=='n', '*', ifelse(Trt=='P'&treat=='Drought years'&trait_cat=='C4 Grass'&drought=='n', '*', ifelse(Trt=='P&N'&treat=='Drought years'&trait_cat=='Non-N-Fixing Forb'&drought=='n', '*', "")))))))


a<-ggplot(data=subset(lfave, treat=="Drought years"), aes(x=Trt, y=mcov, fill=drought, label=label))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mcov-se, ymax=mcov+se), position=position_dodge(0.9), width=.2)+
  ylab('Total Cover')+
  xlab("Nutrient Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  facet_wrap(~trait_cat, nrow=2)+
  geom_text(aes(y=(mcov+se)+5))+
  ggtitle("A) Drought years")



####doing averages by drought only
lfave_drt<-lf_stat%>%
  mutate(calendar_year=as.integer(as.character(calendar_year)))%>%
  mutate(drought=ifelse(precip=="control", "n", "y"))%>%
  group_by(drought, treat, trait_cat)%>%
  summarize(mcov=mean(cov),
            sd=sd(cov),
            n=length(cov))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(label=ifelse(trait_cat=="C3 Gram."&treat=="Drought years"&drought=="n", '*', 
                      ifelse(trait_cat=="C3 Gram."&treat=="Recovery years"&drought=='n', "*", ifelse(trait_cat=='C4 Grass'&treat=='Drought years'&drought=='n', "*", 
       ifelse(trait_cat=='C4 Grass'&treat=='Recovery years'&drought=='n', '*', 
              ifelse(trait_cat=="Non-N-Fixing Forb"&treat=='Drought years'&drought=='n', '*', ""))))))

b<-ggplot(data=filter(lfave_drt, treat=="Recovery years"), aes(x=trait_cat, y=mcov, fill=drought, label=label))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mcov-se, ymax=mcov+se), position=position_dodge(0.9), width=.2)+
  ylab('Total Cover')+
  xlab("Plant Functional Type")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90))+
  geom_text(aes(y=(mcov+se)+5))+
  ggtitle("B) Recovery years")

grid.arrange(a,b, ncol=1)

# ##making figure through time of life forms
# meanlf<-lf%>%
#   group_by(calendar_year, Trt, precip, trait_cat)%>%
#   summarize(mean=mean(relcov), sd=sd(relcov), n=length(relcov))%>%
#   mutate(se=sd/sqrt(n))%>%
#   filter(trait_cat!="NA")%>%
#   mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))
# 
# trtlab<-c(Control="Control", P="P", N="N", 'P&N'="N+P")
# drtlab<-c(control="Not\nDroughted", drought= "Droughted")
# 
# ggplot(subset(meanlf, trt2=="Control"|trt2=="P&N"), aes(x=calendar_year, y=mean, color=trait_cat))+
#   geom_point(size=4)+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
#   facet_grid(precip~trt2, labeller =labeller(precip=drtlab, trt2=trtlab))+
#   geom_line()+
#   scale_color_manual(name="Functional Type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Gram.", "C3 Gram.",  "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
#   geom_vline(xintercept=2012.5)+
#   theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
#   ylab("Relative Cover")+
#   xlab("Year")


###looking at functional compositions differences through time
multdiff_func<-multivariate_difference(lf, time.var="calendar_year", species.var="trait_cat", abundance.var="cov", replicate.var = "unid", treatment.var = "unid3")%>%
  separate(unid3, into=c("plotnum", "Trt", "drought"), sep="_")%>%
  separate(unid32, into=c("plotnum2", "Trt2", "drought2"), sep="_")%>%
  filter(plotnum==plotnum2)%>%
  select(-Trt2, -plotnum2)%>%
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years"))

##
multf.d <- lmer(composition_diff~Trt*as.factor(calendar_year) + (1|plotnum), data=subset(multdiff_func, treat=="Drought years"))
anova(multf.d, ddf="Kenward-Roger")
emmeans(multf.d, pairwise~as.factor(calendar_year), adjust="holm")

multf.r <- lmer(composition_diff~Trt*as.factor(calendar_year) + (1|plotnum), data=subset(multdiff_func, treat=="Recovery years"))
anova(multf.r, ddf="Kenward-Roger")
emmeans(multf.r, pairwise~Trt|as.factor(calendar_year), adjust="holm")

multdiff_func_means<-multdiff_func%>%
  group_by(Trt, calendar_year)%>%
  summarize(mean=mean(composition_diff), sd=sd(composition_diff))%>%
  mutate(se=sd/sqrt(6))%>%
  mutate(data="Func. type")%>%
  bind_rows(multdiff_means)%>%
  mutate(type=factor(data, levels=c("Species", "Func. type")))

#fig for appendix
ggplot(data=multdiff_func_means, aes(x=calendar_year, y=mean, color=Trt))+
  geom_point(size=3)+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  scale_color_manual(name="Nutrient treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  geom_vline(xintercept=2012.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Compositional Difference")+
  xlab("Year")+
  facet_wrap(~type, ncol=1)

multdiff_func_means2<-multdiff_func%>%
  group_by(Trt, treat)%>%
  summarize(mean=mean(composition_diff), sd=sd(composition_diff))%>%
  mutate(se=sd/sqrt(6))%>%
  mutate(data="Func. type")%>%
  bind_rows(multdiff_means2)%>%
  mutate(type=factor(data, levels=c("Species", "Func. type")))%>%
  mutate(stat=ifelse(type=="Species"&treat=="Recovery years"&Trt=="Control"|type=="Func. type"&treat=="Recovery years"&Trt=="Control"|type=="Func. type"&treat=="Recovery years"&Trt=="P", "B", ifelse(type=="Species"&treat=="Recovery years"&Trt=="P&N"|type=="Func. type"&treat=="Recovery years"&Trt=="N"|type=="Func. type"&treat=="Recovery years"&Trt=="P&N", "A", ifelse(type=="Species"&treat=="Recovery years"&Trt=="P"|type=="Species"&treat=="Recovery years"&Trt=="N", "AB", ""))))

a<-ggplot(data=multdiff_func_means2, aes(x=Trt, y=mean, fill=Trt, label=stat))+
  geom_bar(stat="identity", position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2)+
  scale_fill_manual(name="Nutrient treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"))+
  facet_grid(type~treat)+
  xlab("Nutrient treatment")+
  ylab("Compositional difference")+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text(aes(y=(mean+se)+0.05))+
  ggtitle("A)")


##overall fig of RAC of funcational gropus and N+P and Control treatments - this was for ESA talk
lfoverall<-lf%>%
  filter(precip=="control")%>%
  group_by(Trt, trait_cat)%>%
  summarize(mean=mean(relcov))%>%
  filter(trait_cat!="NA")%>%
  mutate(rank=rank(-mean))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

trtlab<-c("Control"="Control", 'P'="P", 'N'='N', 'P&N'="N+P")

rac_func<-ggplot(lfoverall, aes(x=rank, y=mean))+
  geom_line()+
  geom_point(size=5, aes(color=trait_cat))+
  scale_color_manual(name="Functional Type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Grass", "C3 Gram.",  "Annual Grass","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  xlab("Rank")+
  ylab("Relative Cover")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt2, labeller = labeller(trt2=trtlab))+
  scale_x_continuous(limits=c(1,6), breaks=c(1:6))

lf2<-p3plotcomp%>%
  mutate(trait_cat=ifelse(life=="A", 'Annual', ifelse(form=="F"|form=="S", "Forb",
                          ifelse(form=="G", "Grass",
                                  "UNK"))))%>%
  filter(trait_cat!="UNK")%>%
  left_join(treats)%>%
  group_by(plotnum, Trt, precip, calendar_year, trait_cat)%>%
  summarise(cov=sum(abundance))%>%
  left_join(totcov)%>%
  mutate(relcov=cov/tot)

####exploring trait category changes - based on categoreis.....
traits<-p3plotcomp%>%
  select(genus_species, life, form, C3_C4, n_fixer)%>%
  unique()

diff_abund<-abundance_difference(lf, time.var='calendar_year', abundance.var = "cov", replicate.var = "unid", species.var="trait_cat", treatment.var = "unid3")%>%
  separate(unid3, into=c("plotnum", "Trt", "drought"), sep="_")%>%
  separate(unid32, into=c("plotnum2", "Trt2", "drought2"), sep="_")%>%
  filter(plotnum==plotnum2)%>%
  select(-Trt2, -plotnum2)%>%
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years"))

#need to add zeros for everything
diff_abund_fill<-diff_abund%>%
  select(plotnum, Trt, difference, treat, calendar_year, trait_cat)%>%
  group_by(plotnum, Trt, treat, calendar_year)%>%
  spread(trait_cat, difference, fill=0)%>%
  gather(trait_cat, difference, 'Annual Forb':'Non-N-Fixing Forb')
  
meandiffabund<-diff_abund_fill%>%
  group_by(calendar_year, treat, Trt, trait_cat)%>%
  summarize(mean=mean(difference), sd=sd(difference), n=length(difference))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

# ggplot(meandiffabund, aes(x=calendar_year, y=mean, color=Trt))+
#   geom_point()+
#   geom_line()+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
#   facet_wrap(~trait_cat, ncol=2)+
#   scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
#   theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
#   geom_vline(xintercept = 2012.5)+
#   geom_hline(yintercept = 0)+
#   ylab("Cover differences due to drought")+
#   xlab("Nutrient Treatment")

#other approach to figure
meandiffabund2<-diff_abund_fill%>%
  group_by(treat, Trt, trait_cat)%>%
  summarize(mean=mean(difference), sd=sd(difference), n=length(difference))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))%>%
  mutate(stat=ifelse(trait_cat=="C3 Gram."&treat=="Drought years"&Trt=="Control"|trait_cat=="C4 Grass"&treat=="Drought years"&Trt=="P&N"|trait_cat=="Non-N-Fixing Forb"&treat=="Drought years"&Trt=="Control"|trait_cat=="Non-N-Fixing Forb"&treat=="Drought years"&Trt=="P"|trait_cat=="C4 Grass"&treat=="Recovery years"&Trt=="Control"|trait_cat=="C4 Grass"&treat=="Recovery years"&Trt=="P"|trait_cat=="C4 Grass"&treat=="Recovery years"&Trt=="P&N", "A",                                                         ifelse(trait_cat=="C3 Gram."&treat=="Drought years"&Trt=="N"|trait_cat=="C4 Grass"&treat=="Drought years"&Trt=="Control"|trait_cat=="C4 Grass"&treat=="Drought years"&Trt=="P"|trait_cat=="C4 Grass"&treat=="Drought years"&Trt=="N"|trait_cat=="C4 Grass"&treat=="Recovery years"&Trt=="N"|trait_cat=="Non-N-Fixing Forb"&treat=="Drought years"&Trt=="P&N", "B",  
        ifelse(trait_cat=="C3 Gram."&treat=="Drought years"&Trt=="P&N"|trait_cat=="C3 Gram."&treat=="Drought years"&Trt=="P"|trait_cat=="Non-N-Fixing Forb"&treat=="Drought years"&Trt=="N","AB", ""))))

b<-ggplot(meandiffabund2, aes(x=trt2, y=mean, fill=Trt, label=stat))+
  geom_bar(stat="identity", position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width=0.1)+
  facet_grid(treat~trait_cat)+
  scale_fill_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"))+
  geom_hline(yintercept = 0)+
  ylab("Cover differences")+
  xlab("Nutrient treatment")+
  theme(legend.position = "none")+
  geom_text(aes(y=(mean-se)-4))+
  ggtitle("B)")

grid.arrange(a,b, heights=c(1.5, 2))

abund.d <- lmer(difference~Trt*trait_cat*as.factor(calendar_year) + (1|plotnum), data=subset(diff_abund_fill, treat=="Drought years"))
anova(abund.d, ddf="Kenward-Roger")
emmeans(abund.d, pairwise~Trt|trait_cat, adjust="holm")

abund.r <- lmer(difference~Trt*trait_cat*as.factor(calendar_year) + (1|plotnum), data=subset(diff_abund_fill, treat=="Recovery years"))
anova(abund.r, ddf="Kenward-Roger")
emmeans(abund.r, pairwise~Trt|trait_cat, adjust="holm")

##contingency analsysis

anngram<-lf%>%
  mutate(treat=ifelse(calendar_year<2013, "Drought years", "Recovery years"))%>%
  group_by(Trt, plotnum, treat, precip, trait_cat)%>%
  summarize(mean=mean(cov))%>%
  filter(trait_cat!="NA")%>%
  spread(trait_cat, mean, fill=0)%>%
  gather(trait_cat, cov, 'Annual Forb':'Non-N-Fixing Forb')%>%
  filter(trait_cat=="Annual Gram.",
         treat=="Recovery years", 
         precip=="drought")%>%
  mutate(p=ifelse(cov>0, 1, 0))%>%
  group_by(Trt)%>%
  summarize(present=sum(p))%>%
  mutate(absent=6-present)

prop.test(x=as.matrix(anngram[c('present', 'absent')]), alternative='two.sided')  
