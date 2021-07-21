library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(codyn)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

theme_set(theme_bw(20))

#set wd
my.wd<-setwd("C:/Users/mavolio2/Dropbox/Konza Research")

#read in data
p3plotcomp<-read.csv(paste(my.wd, "/p-cubed/SpComp/pcubed_sp_data2010-2015.csv", sep = ""))%>%
  select(year, plotnum, genus_species, cover, precip, life, form, C3_C4, n_fixer)%>%
  rename(calendar_year=year, abundance=cover)

treats<-read.csv(paste(my.wd, "/p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv", sep=""))

comp<-p3plotcomp%>%
  select(calendar_year, plotnum, genus_species, abundance, precip)%>%
  left_join(treats)%>%
  mutate(unid=paste(plotnum, precip, sep="_"), 
         unidt=paste(Trt, precip, sep="_"))

##looking at richness and evenness - i am just not interested in this
          
richeven<-community_structure(comp, time.var="calendar_year", abundance.var="abundance", replicate.var = "unid")%>%
  separate(unid, into=c("plotnum", "drought"), sep="_")%>%
  mutate(plotnum=as.integer(plotnum))%>%
  left_join(treats)%>%
  mutate(treat=ifelse(calendar_year<2013, "drought", "recovery"))

hist(richeven$richness)#this is normal
hist(log(richeven$Evar))#log transfrom even

#richness in plots in 2010
rich2010<-richeven%>%
  filter(calendar_year==2010&drought=="control")%>%
  group_by(Trt)%>%
  summarize(mean=mean(richness))


##richness in drought
rich.d <- lmer(richness~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(richeven, treat=="drought"))

anova(rich.d, ddf="Kenward-Roger")

#doing contrasts
emmeans(rich.d, pairwise~Trt, adjust="holm")

#richness during recovery
rich.r <- lmer(richness~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(richeven, treat=="recovery"))

anova(rich.r, ddf="Kenward-Roger")

#doing contrasts
emmeans(rich.r, pairwise~Trt|as.factor(calendar_year), adjust="holm")

richtoplot<-richeven%>%
  group_by(calendar_year, drought, Trt, treat)%>%
  summarise(mean=mean(richness), sd=sd(richness), n=length(richness))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

trtlab<-c(Control="Control", P="P", N="N", 'P&N'="N+P")

ggplot(richtoplot, aes(x=calendar_year, y=mean, color=trt2, shape=drought))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~trt2, labeller =labeller(trt2=trtlab))+
  geom_line()+
  scale_color_manual(name="Treatment", values=c("black", "blue", "red", "purple"), breaks = c("Control", "P", "N","P&N"))+
  scale_shape_manual(name="Droughted", labels=c("No", 'Yes'), values=c(19,17))+
  geom_vline(xintercept=2012.5)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  ylab("Species Richness")+
  xlab("Year")





###looking at community changes
deltarac<-RAC_change(comp, time.var = "calendar_year", abundance.var = "abundance", replicate.var = "unid", species.var = "genus_species")%>%
  separate(unid, into=c("plotnum", "drought"), sep="_")%>%
  mutate(plotnum=as.integer(plotnum))%>%
  left_join(treats)%>%
  mutate(treat=ifelse(calendar_year<2013, "drought", "recovery"))

ractoplot<-deltarac%>%
  select(calendar_year2, plotnum, drought, Trt, treat, richness_change:losses)%>%
  gather(measure, value, richness_change:losses)%>%
  group_by(calendar_year2, drought, Trt, treat, measure)%>%
  summarise(mean=mean(value), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))

ggplot(ractoplot, aes(x=calendar_year2, y=mean, color=Trt, shape=drought))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_grid(Trt~measure)



##rank in drought
rank.d <- lmer(rank_change~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(deltarac, treat=="drought"))

anova(rank.d, ddf="Kenward-Roger")
emmeans(rank.d, pairwise~drought|Trt|as.factor(calendar_year), adjust="holm")

##rank in recovery
rank.r <- lmer(rank_change~Trt*drought*as.factor(calendar_year) + (1|plotnum), data=subset(deltarac, treat=="recovery"))

anova(rank.r, ddf="Kenward-Roger")
emmeans(rank.r, pairwise~drought|Trt|as.factor(calendar_year), adjust="holm")

##multivariate change
deltamult<-multivariate_change(comp, time.var="calendar_year", treatment.var = "unidt", replicate.var = "unid", species.var = "genus_species", abundance.var = "abundance")%>%
  separate(unidt, into=c("Trt", "drought"), sep="_")

ggplot(data=deltamult, aes(x=calendar_year2, y=composition_change, color=Trt, shape=drought))+
  geom_point()+
  facet_wrap(~Trt)

##multivariate change
deltamult<-multivariate_change(comp, time.var="calendar_year", treatment.var = "unidt", replicate.var = "unid", species.var = "genus_species", abundance.var = "abundance")%>%
  separate(unidt, into=c("Trt", "drought"), sep="_")

ggplot(data=deltamult, aes(x=calendar_year2, y=composition_change, color=Trt, shape=drought))+
  geom_point()+
  facet_wrap(~Trt)

##mult diff
multdiff<-multivariate_difference(comp, time.var="calendar_year", species.var="genus_species", abundance.var="abundance", replicate.var = "unid", treatment.var = "unidt")%>%
  separate(unidt, into=c("Trt", "drought"), sep="_")%>%
  separate(unidt2, into=c("Trt2", "drought2"), sep="_")%>%
  filter(Trt==Trt2)
  
ggplot(data=multdiff, aes(x=calendar_year, y=composition_diff, color=Trt))+
  geom_point(size=3)+
  geom_line()+
  scale_color_manual(name="Treatment", breaks=c("Control", "P", "N", "P&N"), label=c("Control", "P", "N", "N+P"), values=c("black", "blue", "red", "purple"))+
  geom_vline(xintercept=2012.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Compostitional Diff.\nControl-Drought")+
  xlab("Year")
    
  
#rac diff
rankdiff<-RAC_difference(comp, time.var="calendar_year", species.var="genus_species", abundance.var="abundance", replicate.var = "unid", treatment.var = "unidt")

racdiff2<-rankdiff%>%
  separate(unidt, into=c("Trt", "drought"), sep="_")%>%
  separate(unidt2, into=c("Trt2", "drought2"), sep="_")%>%
  separate(unid, into=c("plotnum", "drought2"), sep="_")%>%
  separate(unid2, into=c("plotnum2", "drought2"), sep="_")%>%
  filter(plotnum==plotnum2)%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(Trt, measure, calendar_year)%>%
  summarise(mean=mean(value), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))
  
ggplot(data=racdiff2, aes(x=calendar_year, y=mean, color=Trt))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))+
  facet_grid(measure~Trt, scales="free")


##life form
totcov<-comp%>%
  group_by(plotnum, calendar_year)%>%
  summarize(tot=sum(abundance))

lf<-p3plotcomp%>%
mutate(trait_cat=ifelse(form=="F"&life=="A", "Annual Forb",
                        ifelse(form=="G"&life=="A", "Annual Gram.",
                               ifelse(form=="G"&C3_C4=="C3", "C3 Gram.",
                                      ifelse(form=="G"&C3_C4=="C4", "C4 Gram.",
                                             ifelse(form=='F'|form=="S"&n_fixer=="N", "Non-N-Fixing Forb",
                                                    ifelse(form=='F'|form=="S"&n_fixer=="Y", "N-Fixing Forb","UNK")))))))%>%
  left_join(treats)%>%
  group_by(plotnum, Trt, precip, calendar_year, trait_cat)%>%
  summarise(cov=sum(abundance))%>%
  left_join(totcov)%>%
  mutate(relcov=cov/tot)

meanlf<-lf%>%
  group_by(calendar_year, Trt, precip, trait_cat)%>%
  summarize(mean=mean(relcov), sd=sd(relcov), n=length(relcov))%>%
  mutate(se=sd/sqrt(n))%>%
  filter(trait_cat!="NA")%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

trtlab<-c(Control="Control", P="P", N="N", 'P&N'="N+P")
drtlab<-c(control="Not\nDroughted", drought= "Droughted")

ggplot(subset(meanlf, trt2=="Control"|trt2=="P&N"), aes(x=calendar_year, y=mean, color=trait_cat))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1)+
  facet_grid(precip~trt2, labeller =labeller(precip=drtlab, trt2=trtlab))+
  geom_line()+
  scale_color_manual(name="Functional Type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Gram.", "C3 Gram.",  "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  geom_vline(xintercept=2012.5)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))+
  ylab("Relative Cover")+
  xlab("Year")

###overall fig
lfoverall<-lf%>%
  filter(Trt=="Control"&precip=="control")%>%
  group_by(trait_cat)%>%
  summarize(mean=mean(relcov), sd=sd(relcov), n=length(relcov))%>%
  mutate(se=sd/sqrt(n))%>%
  filter(trait_cat!="NA")%>%
  mutate(rank=rank(-mean))

ggplot(lfoverall, aes(x=rank, y=mean))+
  geom_line()+
  geom_point(size=5, aes(color=trait_cat))+
  scale_color_manual(name="Functional Type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Gram.", "C3 Gram.",  "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  xlab("Rank")+
  ylab("Relative Cover")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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

####exploring trait category changes
traits<-p3plotcomp%>%
  select(genus_species, life, form, C3_C4, n_fixer)%>%
  unique()

deltabund<-abundance_change(comp, time.var='calendar_year', abundance.var = "abundance", replicate.var = "unid", species.var="genus_species")

deltabund2<-deltabund%>%
  separate(unid, into=c("plotnum", "drought"), sep="_")%>%
  mutate(plotnum=as.integer(as.character(plotnum)))%>%
  left_join(traits)%>%
  mutate(trait_cat=ifelse(form=="F"&life=="A", "Annual Forb",
                          ifelse(form=="G"&life=="A", "Annual Gram.",
                                 ifelse(form=="G"&C3_C4=="C3", "C3 Gram.",
                                        ifelse(form=="G"&C3_C4=="C4", "C4 Gram.",
                                               ifelse(form=='F'|form=="S"&n_fixer=="N", "Non-N-Fixing Forb",
                                                      ifelse(form=='F'|form=="S"&n_fixer=="Y", "N-Fixing Forb","UNK")))))))%>%
  left_join(treats)%>%
  group_by(plotnum, trait_cat, drought, Trt, calendar_year2)%>%
  summarise(abund=mean(change))%>%
  mutate(treat=ifelse(calendar_year2<2013, "drought", "recovery"))

abund.d <- lmer(abund~Trt*drought*as.factor(calendar_year2) + (1|plotnum), data=subset(deltabund2, treat=="drought"&trait_cat=="C4 Gram."))

anova(rank.r, ddf="Kenward-Roger")
emmeans(rank.r, pairwise~drought|Trt|as.factor(calendar_year), adjust="holm")


###doing NMDS
ave<-comp%>%
  group_by(calendar_year, Trt, precip, genus_species)%>%
  summarize(mabund=mean(abundance))

#make the dataset wide for vegan
compwide<-ave%>%
  spread(genus_species, mabund, fill=0)

#pull out plot info
plots<-compwide[,1:3]

#run nmds
mds<-metaMDS(compwide[,4:99], trymax = 100)
mds

#extract NMDS coordinates and bind to plot info for graphs
#the label step adds the label for the first and last year of data
scores<-plots%>%
  bind_cols(as.data.frame(mds$points))%>%
  mutate(label=ifelse(calendar_year==2010|calendar_year==2015, as.character(calendar_year),""),
         Trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

trtlab<-c(Control="Control", P="P", N="N", 'P&N'="N+P")

#make figure with first and last year labeled and path between points connected
ggplot(data=scores, aes(x=MDS1, y=MDS2, color=Trt2, shape=precip, label=label))+
  geom_point(size=4)+
  geom_path()+
  scale_shape_manual(name="Droughted", labels=c("No", "Yes"), values=c(17, 6))+
  geom_text_repel(show.legend = F)+
  scale_color_manual(name="Treatment", values = c("black", "blue", "red", "purple"), labels=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color=F)+
  ylab("NMDS2")+
  xlab("NMDS1")+
  facet_wrap(~Trt2, labeller = labeller(Trt2=trtlab))
