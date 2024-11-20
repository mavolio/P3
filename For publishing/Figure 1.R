#community composition through time

library(vegan)
library(tidyverse)
library(ggrepel)
library(gridExtra)

theme_set(theme_bw(16))

#set wd
setwd("C:/Users/mavolio2/Dropbox/Konza Research/")

#read in data
treats<-read.csv("p-cubed/Analyses/July 2015 Analyses/PPlot_PlotList.csv")%>% 
  rename(plot_id=plotnum)

pplotcomp<-read.csv("pplots\\Sppcomp\\Species Comp_to Use\\Compiling Datasets in R\\Spp_Data_2002_2021.csv") %>% 
  left_join(treats) %>% 
  filter(!is.na(Trt)) %>% 
  mutate(spnum=ifelse(genus_species=='solidago_canadensis', 533, spnum))


cattraits<-read.csv("pplots/traits_2021.csv")%>%
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
  group_by(calendar_year, Trt, genus_species)%>%
  summarize(mabund=mean(abundance))

#make the dataset wide for vegan
compwide<-ave%>%
  spread(genus_species, mabund, fill=0)

#pull out plot info
plots<-compwide[,1:2]

#run nmds
mds<-metaMDS(compwide[,3:90], trymax = 100)
mds

#extract NMDS coordinates and bind to plot info for graphs
#the label step adds the label for the first and last year of data
scores<-plots%>%
  bind_cols(as.data.frame(mds$points))%>%
  mutate(label=ifelse(calendar_year==2002, "start", "nope"),
         Drought=ifelse(calendar_year<2010, "No", ifelse(calendar_year>2009&calendar_year<2013, "Yes, Drought year", 'Yes, Recovery year')))

#make figure with first and last year labeled and path between points connected
NMDS<-
ggplot(data=scores, aes(x=MDS1, y=MDS2, color=Trt, fill=label, shape=Drought, group=Trt))+
  geom_point(size=4)+
  geom_path()+
  scale_shape_manual(name="Year in\nDrought experiment", values=c(21, 15, 17))+
  scale_fill_manual(guide=F, values = c('white', 'gray'))+
  scale_color_manual(name="Nutrient treatment", values = c("black", "blue", "red", "purple"), breaks=c("Control", "P", "N", "P&N"), labels=c("Control", "P", "N", "N+P"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("NMDS2")+
  xlab("NMDS1")+
  annotate("text", x = 0.9, y = -.5, label="Stress = 0.18")
NMDS

#permanova
adonis2(compwide[3:90]~compwide$Trt)

#perm dispersion
dist<-vegdist(compwide[3:90])
betadisp<-betadisper(dist,compwide$Trt,type="centroid")
permutest(betadisp)

#Doing PFTs
#get average cover of each species in each treatment over all years (average over plots and years). 

PFTcover<-pplotcomp%>%
  left_join(cattraits, by="spnum") %>% 
  filter(spnum<900) %>% 
  filter(calendar_year==2010) %>% 
  mutate(trait_cat=ifelse(growthform=="F"&lifespan=="A", "Annual Forb",
                   ifelse(growthform=="G"&lifespan=="A", "Annual Grass",
                   ifelse(growthform=="G"&photopath=="C3", "C3 Grass",
                   ifelse(growthform=="G"&photopath=="C4", "C4 Grass",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="N", "Non-N-Fixing Forb",
                   ifelse(growthform=='F'|growthform=="S"&Nfix=="Y", "N-Fixing Forb","UNK")))))))%>%
  left_join(treats)%>%
  group_by(plot_id, Trt, calendar_year, trait_cat)%>%
  summarise(cov=sum(abundance))

racave<-PFTcover%>%
  group_by(Trt, trait_cat) %>% 
  summarize(meancover=mean(cov)) %>% 
  mutate(rank=rank(-meancover, ties.method = "first")) %>% 
  mutate(trt2=factor(Trt, levels=c("Control", "P", "N", "P&N"))) 

#make new labels for facet_wrap step  
collabel<-c("Control"="Control", "P"="P","N"="N","P&N"="N+P")

#great rac figure
rac<-
ggplot(data=racave, aes(x=rank, y=meancover))+
  geom_line()+
  geom_point(aes(color=trait_cat), size=4)+
  scale_color_manual(name="Functional type", values=c("darkgreen", "chartreuse3", "darkolivegreen1", "darkblue", "lightblue", "deepskyblue"), breaks = c("C4 Grass", "C3 Grass",  "Annual Grass","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"), labels=c(bquote('C'[4]~'Grass'), bquote('C'[3]~'Grass'),"Annual Grass","Non-N-Fixing Forb", "N-Fixing Forb", "Annual Forb"))+
  facet_wrap(~trt2, labeller = labeller(trt2=collabel))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Cover")+
  xlab("Rank")
rac


#bind both figures together.
grid.arrange(NMDS, rac)

