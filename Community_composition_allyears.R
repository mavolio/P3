#community composition through time

library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)

theme_set(theme_bw(14))

#set wd
my.wd<-setwd("C:/Users/mavolio2/Dropbox/Konza Research/")

#read in data
treats<-read.csv("p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv") %>% 
  rename(plot_id=plotnum)


pplotcomp<-read.csv("pplots\\Sppcomp\\Species Comp_to Use\\Compiling Datasets in R\\Spp_Data_2002_2021.csv") %>% 
  left_join(treats) %>% 
  filter(!is.na(Trt))


cattraits<-read.csv(paste(my.wd, "/pplots/traits_2021.csv", sep=""))%>%
  filter(Genus!="")%>%
  rename(spnum=Species.Number) %>% 
  select(spnum, Genus, Species, Annual.Peren.Bi, Forb.grass.shrub, C3.C4, N.fixer..Y.N...)%>%
  mutate(genus_species=paste(tolower(Genus), Species, sep="_"))%>%
  rename(lifespan=Annual.Peren.Bi,
         growthform=Forb.grass.shrub,
         photopath=C3.C4,
         Nfix=N.fixer..Y.N...)%>%
  select(-Genus, -Species)

#get average cover of each species in a treatment for each year of the experiment (average over the plots)
ave<-pplotcomp%>%
  filter(calendar_year<2016) %>% 
  group_by(calendar_year, Trt, spnum, genus_species)%>%
  summarize(mabund=mean(abundance))

#make the dataset wide for vegan
compwide<-ave%>%
  spread(spnum, mabund, fill=0)

#pull out plot info
plots<-compwide[,1:2]

#run nmds
mds<-metaMDS(compwide[,3:76], trymax = 100)
mds

#extract NMDS coordinates and bind to plot info for graphs
#the label step adds the label for the first and last year of data
scores<-plots%>%
  bind_cols(as.data.frame(mds$points))%>%
  mutate(label=ifelse(calendar_year==2002, as.character(calendar_year),""),
         year=ifelse(calendar_year>2009, "n", 'y'))

#make figure with first and last year labeled and path between points connected
nmds<-ggplot(data=scores, aes(x=MDS1, y=MDS2, color=Trt, label=label, shape=year, group=Trt))+
  geom_point(size=3)+
  geom_path()+
  scale_shape_manual(guide='none', values=c(2, 17))+
  geom_label_repel(show.legend = F)+
  scale_color_manual(name="Nutrient treatment", values = c("black", "blue", "red", "purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("NMDS2")+
  xlab("NMDS1")+
  annotate("text", x = 0.5, y = -.5, label="Stress = 0.18")

#permanova
adonis(compwide[3:73]~compwide$Trt)

#perm dispersion
dist<-vegdist(compwide[3:73])
betadisp<-betadisper(dist,compwide$Trt,type="centroid")
permutest(betadisp)

#get average cover of each species in each treatment over all years (average over plots and years). 
#join categorical triats and create trait categories
racave<-ave%>%
  group_by(Trt, spnum)%>%
  summarize(mabund=mean(mabund))%>%
  mutate(rank=rank(-mabund, ties.method = "first"))%>%
  left_join(cattraits)%>%
  mutate(trait_cat=ifelse(growthform=="F"&lifespan=="A", "Annual Forb",
                   ifelse(growthform=="G"&lifespan=="A", "Annual Grass",
                   ifelse(growthform=="G"&photopath=="C3", "C3 Gram.",
                   ifelse(growthform=="G"&photopath=="C4", "C4 Grass",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="N", "Non-N-Fixing Forb",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="Y", "N-Fixing Forb","UNK")))))))%>%
  separate(genus_species, into=c("genus", "species"), sep="_")%>%
  mutate(genera=toupper(substr(genus, 1, 1)),
         sp=paste(genera, species, sep=". "))%>%
  mutate(name=ifelse(rank<6, sp, ""))%>%
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N")))

#make new labels for facet_wrap step  
collabel<-c("Control"="Control", "P"="P","N"="N","P&N"="N+P")

#great rac figure
rac<-
ggplot(data=racave, aes(x=rank, y=mabund, label=name))+
  geom_line()+
  geom_point(aes(color=trait_cat), size=2)+
  scale_color_manual(name="Functional type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Gram.", "C3 Gram.",  "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  geom_text_repel(max.overlaps = 10, size=3)+
  facet_wrap(~trt2, labeller = labeller(trt2=collabel))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Abundance")+
  xlab("Rank")

#bind both figures together.
grid.arrange(nmds, rac)

#run RAC func code in community analyses.

grid.arrange(NMDS, rac_func)



#plot changes over time ESA talk

allpp<-read.csv(paste(my.wd, "/pplots/Sppcomp/Species Comp_to Use/Compiling Datasets in R/Spp_Data_2002_2019.csv", sep = ""))%>%
  filter(treatment=="N1P0"|treatment=="N2P0"|treatment=="N1P3"|treatment=="N2P3")

cattraits<-read.csv(paste(my.wd, "/pplots/traits_2021.csv", sep=""))%>%
  filter(Genus!="")%>%
  select(Genus, Species, Annual.Peren.Bi, Forb.grass.shrub, C3.C4, N.fixer..Y.N...)%>%
  mutate(genus_species=paste(tolower(Genus), Species, sep="_"))%>%
  rename(lifespan=Annual.Peren.Bi,
         growthform=Forb.grass.shrub,
         photopath=C3.C4,
         Nfix=N.fixer..Y.N...)%>%
  select(-Genus, -Species)

#get average cover of each species in a treatment for each year of the experiment (average over the plots)
aveall<-allpp%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarize(mabund=mean(abundance))


                     
#get average cover of each species in each treatment over all years (average over plots and years). 
#join categorical triats and create trait categories
ractime<-aveall%>%
  left_join(cattraits)%>%
  mutate(trait_cat=ifelse(growthform=="F"&lifespan=="A", "Annual Forb",
                          ifelse(growthform=="G"&lifespan=="A", "Annual Grass",
                                 ifelse(growthform=="G"&photopath=="C3", "C3 Gram.",
                                        ifelse(growthform=="G"&photopath=="C4", "C4 Grass",
                                               ifelse(growthform=='F'|growthform=="S"&Nfix=="N", "Non-N-Fixing Forb",
                                                      ifelse(growthform=='F'|growthform=="S"&Nfix=="Y", "N-Fixing Forb","UNK")))))))%>%
  group_by(trait_cat, calendar_year, treatment)%>%
  summarise(abund=sum(mabund))

#make new labels for facet_wrap step  
collabel<-c(
  N1P0="Control", 
  N1P3="P",
  N2P0="N",
  N2P3="N+P")

ggplot(data=ractime, aes(x=calendar_year, y=abund, group=trait_cat))+
  geom_line()+
  geom_point(aes(color=trait_cat), size=4)+
  scale_color_manual(name="Functional Type", values=c("darkgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Gram.", "C3 Gram.",  "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  #geom_text_repel(max.overlaps = 10, size=5)+
  facet_wrap(~treatment, labeller = labeller(treatment=collabel))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Abundance")+
  xlab("Year")
