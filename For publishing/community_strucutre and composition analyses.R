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
emmeans(lfcov.d, pairwise~precip|Trt|trait_cat|year, adjust="holm")


lfcov.r <- lmer(cov~Trt*precip*trait_cat*year +  (1|plot/precip/year) + (1|plot:year)+ (1|plot:precip:trait_cat) + (1|plot:trait_cat) + (1|plot:year:trait_cat), data=subset(lf_stat, treat=="Recovery years"))
anova(lfcov.r, ddf="Kenward-Roger")
emmeans(lfcov.r, pairwise~precip|Trt|trait_cat|year, adjust="holm")


lfave<-lf_stat%>%
  mutate(year=as.integer(as.character(year)))%>%
  mutate(drought=ifelse(precip=="control", "n", "y"))%>%
  group_by(drought, treat, Trt, year, trait_cat)%>%
  summarize(mcov=mean(cov),
            sd=sd(cov),
            n=length(cov))%>%
  mutate(se=sd/sqrt(n)) %>% 
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))) %>% 
  group_by(treat, Trt, trait_cat, year) %>% 
  mutate(y=max(mcov)) %>% 
  mutate(label=ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="C3 Grass"&year==2010, "*", 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="C4 Grass"&year==2010, "*", 
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P"&trait_cat=="C3 Grass"&year==2011, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="N"&trait_cat=="Non-N-Fixing Forb"&year==2011, '*',
              ifelse(drought=="n"&treat=="Drought years"&Trt=="P&N"&trait_cat=="Non-N-Fixing Forb"&year==2011, '*',
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
              ifelse(drought=="n"&treat=="Recovery years"&Trt=="P&N"&trait_cat=="Non-N-Fixing Forb"&year==2015, '*', "")))))))))))))))))


#figure for the paper
ggplot(data=subset(lfave), aes(x=as.factor(year), y=mcov, color=drought, label=label))+
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
  geom_text(aes(y=y+14), color="red", size=5)+
  facet_grid(trait_cat~trt2, scales='free')
  


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
