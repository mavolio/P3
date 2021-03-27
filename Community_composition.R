#community composition through time

library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)

theme_set(theme_bw(12))

#set wd
my.wd<-setwd("C:/Users/mavolio2/Dropbox/Konza Research")

#read in data
pplotcomp<-read.csv(paste(my.wd, "/pplots/Sppcomp/Species Comp_to Use/Compiling Datasets in R/Spp_Data_2002_2019.csv", sep = ""))%>%
  filter(calendar_year>2009&calendar_year<2016)%>%
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
ave<-pplotcomp%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarize(mabund=mean(abundance))

#make the dataset wide for vegan
compwide<-ave%>%
  spread(genus_species, mabund, fill=0)

#pull out plot info
plots<-compwide[,1:2]

#run nmds
mds<-metaMDS(compwide[,3:69])
mds

#extract NMDS coordinates and bind to plot info for graphs
#the label step adds the label for the first and last year of data
scores<-plots%>%
  bind_cols(as.data.frame(mds$points))%>%
  mutate(label=ifelse(calendar_year==2010|calendar_year==2015, as.character(calendar_year),""))

#make figure with first and last year labeled and path between points connected
NMDS<-
ggplot(data=scores, aes(x=MDS1, y=MDS2, color=treatment, label=label))+
  geom_point(size=3)+
  geom_path()+
  geom_text_repel(show.legend = F)+
  scale_color_manual(name="Treatment", values = c("black", "blue", "red", "purple"), labels=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("NMDS2")+
  xlab("NMDS1")

#permanova
adonis(compwide[3:69]~compwide$treatment)

#perm dispersion
dist<-vegdist(compwide[3:69])
betadisp<-betadisper(dist,compwide$treatment,type="centroid")
permutest(betadisp)

#get average cover of each species in each treatment over all years (average over plots and years). 
#join categorical triats and create trait categories
racave<-ave%>%
  group_by(treatment, genus_species)%>%
  summarize(mabund=mean(mabund))%>%
  mutate(rank=rank(-mabund, ties.method = "first"))%>%
  left_join(cattraits)%>%
  mutate(trait_cat=ifelse(growthform=="F"&lifespan=="A", "Annual Forb",
                   ifelse(growthform=="G"&lifespan=="A", "Annual Gram.",
                   ifelse(growthform=="G"&photopath=="C3", "C3 Gram.",
                   ifelse(growthform=="G"&photopath=="C4", "C4 Gram.",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="N", "Non-N-Fixing Forb",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="Y", "N-Fixing Forb","UNK")))))))%>%
  separate(genus_species, into=c("genus", "species"), sep="_")%>%
  mutate(genera=toupper(substr(genus, 1, 1)),
         sp=paste(genera, species, sep=". "))%>%
  mutate(name=ifelse(rank<4, sp, ""))

#make new labels for facet_wrap step  
collabel<-c(
  N1P0="Control", 
  N1P3="P",
  N2P0="N",
  N2P3="N+P")

#great rac figure
rac<-
ggplot(data=racave, aes(x=rank, y=mabund, label=name))+
  geom_line()+
  geom_point(aes(color=trait_cat), size=2)+
  scale_color_manual(name="Functional Type", values=c("forestgreen", "chartreuse3", "green", "darkblue", "lightblue", "deepskyblue"), breaks = c("C3 Gram.", "C4 Gram.", "Annual Gram.","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  geom_text_repel(max.overlaps = 8, size=3)+
  facet_wrap(~treatment, labeller = labeller(treatment=collabel))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Abundance")+
  xlab("Rank")

#bind both figures together.
grid.arrange(NMDS, rac, ncol=2)
                     
