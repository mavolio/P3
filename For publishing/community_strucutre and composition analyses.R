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
setwd("C:/Users/mavolio2/Dropbox/Konza Research/")

#read in data
p3plotcomp<-read.csv("p-cubed/SpComp/pcubed_sp_data2010-2015.csv")%>%
  select(year, plotnum, spnum, genus_species, cover, precip, life, form, C3_C4, n_fixer)%>%
  rename(calendar_year=year, abundance=cover)

treats<-read.csv("p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv") 

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
emmeans(rich.d, pairwise~Trt, adjust="tukey")

#richness during recovery
rich.r <- lmer(richness~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Recovery years"))

anova(rich.r, ddf="Kenward-Roger")

#doing contrasts
emmeans(rich.r, pairwise~Trt, adjust="tukey")

##evenness in drought
even.d <- lmer(log(Evar)~Trt*drought*year + (1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Drought years"))
anova(even.d, ddf="Kenward-Roger")
#doing contrasts
emmeans(even.d, pairwise~Trt, adjust="holm")

#evenness during recovery
even.r <- lmer(log(Evar)~Trt*drought*year +(1|plot) + (1|plot:drought)+(1|plot:year), data=subset(richeven, treat=="Recovery years"))
anova(even.r, ddf="Kenward-Roger")


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
  
ggplot(data=multdiff_means2, aes(x=Trt, y=mean, label=label))+
  geom_bar(stat='identity', size=3, fill="gray")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
  scale_x_discrete(limits=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Euclidean distance between\ndroughted and un-droughted subplots")+
  xlab("Nutrient Treatment")+
  facet_wrap(~treat)+
  geom_text(aes(y=mean+se+0.05))+
  theme(plot.title = element_text(hjust=-0.1, vjust = -8))


##looking at cover of functional types over time
totcov<-comp%>%
  group_by(plotnum, calendar_year)%>%
  summarize(tot=sum(abundance))

lf<-p3plotcomp%>%
  filter(spnum<990) %>% #omits unknowns
mutate(trait_cat=ifelse(form=="F"&life=="A", "Annual Forb",
                 ifelse(form=="G"&life=="A", "Annual Grass",
                 ifelse(form=="G"&C3_C4=="C3", "C3 Grass",
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
lfcov.d <- lmer(cov~Trt*precip*trait_cat*year + (1|plot/precip/year) + (1|plot:year)+ (1|plot:precip:trait_cat) + (1|plot:trait_cat) + (1|plot:year:trait_cat), data=subset(lf_stat, treat=="Drought years"))
anova(lfcov.d, ddf="Kenward-Roger")

#doing contrasts
em<-emmeans(lfcov.d, pairwise~precip|Trt|trait_cat|year, adjust="tukey")
emmeans(lfcov.d, pairwise~year|Trt|trait_cat|precip, adjust="tukey")

lfcov.r <- lmer(cov~Trt*precip*trait_cat*year +  (1|plot/precip/year) + (1|plot:year)+ (1|plot:precip:trait_cat) + (1|plot:trait_cat) + (1|plot:year:trait_cat), data=subset(lf_stat, treat=="Recovery years"))
anova(lfcov.r, ddf="Kenward-Roger")
emmeans(lfcov.r, pairwise~precip|Trt|trait_cat|year, adjust="tukey")
emmeans(lfcov.r, pairwise~year|Trt|trait_cat|precip, adjust="tukey")


lfave<-lf_stat%>%
  mutate(year=as.integer(as.character(year)))%>%
  mutate(drought=ifelse(precip=="control", "n", "y"))%>%
  group_by(drought, treat, Trt, year, trait_cat)%>%
  summarize(mcov=mean(cov),
            sd=sd(cov),
            n=length(cov))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P")),
         cat2=factor(trait_cat, levels=c('C4 Grass', 'C3 Grass', 'Annual Grass', 'Non-N-Fixing Forb', 'N-Fixing Forb', 'Annual Forb'))) %>% 
  group_by(treat, Trt, trait_cat, year) %>% 
  mutate(y=max(mcov)) %>% 
  mutate(label=ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="C3 Grass"&year==2010, "*", 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="C4 Grass"&year==2010, "*", 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P"&trait_cat=="C4 Grass"&year==2010, "*", 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P"&trait_cat=="C3 Grass"&year==2011, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="Non-N-Fixing Forb"&year==2011, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P&N"&trait_cat=="Non-N-Fixing Forb"&year==2011, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="Control"&trait_cat=="C4 Grass"&year==2012, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="C4 Grass"&year==2012, '*', 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P"&trait_cat=="C4 Grass"&year==2012, '*', 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P&N"&trait_cat=="C4 Grass"&year==2012, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="P&N"&trait_cat=="Annual Grass"&year==2013,'*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="N"&trait_cat=="C3 Grass"&year==2013, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="P"&trait_cat=="C3 Grass"&year==2014, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="N"&trait_cat=="C4 Grass"&year==2014, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="N"&trait_cat=="Non-N-Fixing Forb"&year==2014, '*',
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="P&N"&trait_cat=="C3 Grass"&year==2015, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="N"&trait_cat=="C4 Grass"&year==2015, '*', 
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="P&N"&trait_cat=="Non-N-Fixing Forb"&year==2015, '*', ""))))))))))))))))))) 


#figure for the paper
ggplot(data=lfave, aes(x=as.factor(year), y=mcov, color=drought, label=label))+
  geom_rect(data=NULL,aes(xmin=3.5,xmax=6.56,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  geom_point()+
  geom_line(aes(group=drought))+
  scale_color_manual(name="Droughted", values=c("Blue", "Orange"), labels=c("No", "Yes"))+
  geom_errorbar(aes(ymin=mcov-se, ymax=mcov+se), width=.2)+
  ylab('Cover')+
  xlab("Year")+
  scale_x_discrete(labels=c("'10", "'11", "'12","'13", "'14", "'15"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept=3.5)+
  geom_text(aes(y=y+15), color="red", size=5)+
  facet_grid(cat2~trt2, scales='free')
  
#focusing on three cover types
ggplot(data=subset(lfave, cat2=='Non-N-Fixing Forb'&drought=='y'|cat2=='C4 Grass'&drought=='y'|cat2=='C3 Grass'&drought=='y'), aes(x=as.factor(year), y=mcov, color=Trt))+
  geom_rect(data=NULL,aes(xmin=3.5,xmax=6.56,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  geom_point()+
  geom_line(aes(group=Trt))+
  scale_color_manual(name="Nutrient treatment", values = c("black", "blue", "red", "purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  geom_errorbar(aes(ymin=mcov-se, ymax=mcov+se), width=.2)+
  ylab('Cover')+
  xlab("Year")+
  scale_x_discrete(labels=c("'10", "'11", "'12","'13", "'14", "'15"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept=3.5)+
  #geom_text(aes(y=y+15), color="red", size=5)+
  facet_wrap(~cat2, scales='free', ncol=1)



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
