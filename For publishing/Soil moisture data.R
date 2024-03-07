library(tidyverse)
library(gridExtra)
library(lme4)
library(lmerTest)

setwd("C:/Users/mavolio2/Dropbox/Konza Research/p-cubed/soil moisture/SM_Analyses")

theme_set(theme_bw(12))

###import all data
gravi2012<-read.csv("Gravi_2012.csv")
cali2010.2011<-read.csv("probe calibrate_2010_2011.csv")
cali2012<-read.csv("probe_calibrate_2012.csv")
info2010<-read.csv("Info2010.csv")
info2011<-read.csv("Info2011.csv")
info2012<-read.csv("Info2012.csv")
knza2010<-read.csv("KNZA_2010.csv")%>%
  gather(Port, vwc, Port2:Port5)%>%
  filter(vwc!=999)%>%
  mutate(year=2010)
knza2011<-read.csv("KNZA_2011.csv")%>%
  gather(Port, vwc, Port1:Port4)%>%
  filter(vwc!=999)%>%
  mutate(year=2011)
knza2012<-read.csv("KNZA_2012.csv")%>%
  gather(Port, vwc, Port1:Port5)%>%
  filter(vwc!=999)%>%
  mutate(year=2012)
knzc2010<-read.csv("KNZC_2010.csv")%>%
  gather(Port, vwc, Port1:Port4)%>%
  filter(vwc!=999)%>%
  mutate(year=2010)
knzd2010<-read.csv("KNZD_2010.csv")%>%
  gather(Port, vwc, Port1:Port4)%>%
  filter(vwc!=999)%>%
  mutate(year=2010)
knzd2011<-read.csv("KNZD_2011.csv")%>%
  gather(Port, vwc, Port1:Port4)%>%
  filter(vwc!=999)%>%
  mutate(year=2011)
knzd2012<-read.csv("KNZD_2012.csv")%>%
  gather(Port, vwc, Port1:Port5)%>%
  filter(vwc!=999)%>%
  mutate(year=2012)
sgsd2010<-read.csv("SGSD_2010.csv")%>%
  gather(Port, vwc, Port1:Port5)%>%
  filter(vwc!=999)%>%
  mutate(year=2010)

##do calibrations 2010/2011
cali10.11.aveprobe<-cali2010.2011%>%
  gather(probe, vwc, pp1:D8)%>%
  filter(vwc!=999)%>%
  group_by(sm, probe)%>%
  summarize(raw=mean(vwc))

cali10.11.avesm<-cali2010.2011%>%
  gather(probe, vwc, pp1:D8)%>%
  filter(vwc!=999)%>%
  group_by(sm)%>%
  summarize(ave=mean(vwc))

calibrate10.11<-merge(cali10.11.avesm, cali10.11.aveprobe, by="sm")
#y~x
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp3")))
p3<-data.frame(Probe="P3",
               intercept = -0.02105,
               slope = 0.96130)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp11")))
p11<-data.frame(Probe="P11",
                intercept = 0.006207,
                slope = 0.963220)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp2")))
p2<-data.frame(Probe="P2",
               intercept = -0.02235,
               slope = 1.13126)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp17")))
p17<-data.frame(Probe="P17",
                intercept = -0.009771,
                slope = 1.089083)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp14")))
p14<-data.frame(Probe="P14",
                intercept = -0.005979,
                slope = 0.876904)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp5")))
p5<-data.frame(Probe="P5",
               intercept = 0.006034,
               slope = 0.901583)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp9")))
p9<-data.frame(Probe="P9",
               intercept = -0.002676,
               slope = 1.077390)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp12")))
p12<-data.frame(Probe="P12",
                intercept = -0.02544,
                slope = 1.26919)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp10")))
p10<-data.frame(Probe="P10",
                intercept = -0.02166,
                slope = 1.16341)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp1")))
p1<-data.frame(Probe="P1",
               intercept = -0.01567,
               slope = 1.05902)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp20")))
p20<-data.frame(Probe="P20",
                intercept = -0.05282,
                slope = 1.33588)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp4")))
p4<-data.frame(Probe="P4",
               intercept = -0.02746,
               slope = 1.19306)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp16")))
p16<-data.frame(Probe="P16",
                intercept = -0.03131,
                slope = 1.03413)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp13")))
p13<-data.frame(Probe="P13",
                intercept = -0.01667,
                slope = 1.12708)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp15")))
p15<-data.frame(Probe="P15",
                intercept = -0.004587,
                slope = 1.162162)
summary(lm(raw~ave, subset(calibrate10.11, probe=="pp18")))
p18<-data.frame(Probe="P18",
                intercept = -0.03214,
                slope = 1.27366)
summary(lm(raw~ave, subset(calibrate10.11, probe=="D7")))
pd7<-data.frame(Probe="D7",
                intercept = 0.001531,
                slope = 1.108501)

probes.2010.2011<-rbind(p1, p10, p11, p12, p13, p14, p15, p16, p17, p18, p2, p20, p3, p4, p5, p9, pd7)

##Calibrate 2012
cali2012.ave<-cali2012%>%
  gather(potrep, vwc, pot0_1:pot5_3)%>%
  separate(potrep, into=c("Pot","rep"), sep="_")%>%
  group_by(probe, Pot)%>%
  summarize(raw=mean(vwc))

gravi2012.ave<-gravi2012%>%
  group_by(Pot)%>%
  summarize(sm=mean(Gravi))

calibrate12<-merge(cali2012.ave, gravi2012.ave, by="Pot")

#y~x
summary(lm(raw~sm, subset(calibrate12, probe=="202")))
p202<-data.frame(Probe="p202",
                 intercept = 18.314,
                 slope = 4.291)
summary(lm(raw~sm, subset(calibrate12, probe=="201")))
p201<-data.frame(Probe="p201",
                 intercept = 15.583,
                 slope = 4.369)
summary(lm(raw~sm, subset(calibrate12, probe=="208")))
p208<-data.frame(Probe="p208",
                 intercept = 17.722,
                 slope = 4.304)
summary(lm(raw~sm, subset(calibrate12, probe=="203")))
p203<-data.frame(Probe="p203",
                 intercept = 18.516,
                 slope = 4.243)
summary(lm(raw~sm, subset(calibrate12, probe=="210")))
p210<-data.frame(Probe="p210",
                 intercept = 14.656,
                 slope = 4.353)
summary(lm(raw~sm, subset(calibrate12, probe=="207")))
p207<-data.frame(Probe="p207",
                 intercept = 15.23,
                 slope = 4.38)
summary(lm(raw~sm, subset(calibrate12, probe=="211")))
p211<-data.frame(Probe="p211",
                 intercept = 14.568,
                 slope = 4.373)
summary(lm(raw~sm, subset(calibrate12, probe=="206")))
p206<-data.frame(Probe="p206",
                 intercept = 8.642,
                 slope = 4.488)
summary(lm(raw~sm, subset(calibrate12, probe=="204")))
p204<-data.frame(Probe="p204",
                 intercept = 13.02,
                 slope = 4.41)
summary(lm(raw~sm, subset(calibrate12, probe=="205")))
p205<-data.frame(Probe="p205",
                 intercept = 15.223,
                 slope = 4.296)
probes.2012<-rbind(p202,p201,p208,p210,p207,p211, p206, p204, p203, p205)

probes_equations<-rbind(probes.2012, probes.2010.2011)

###merge all soil moisture data
sm2010.raw<-rbind(knza2010, knzc2010, knzd2010, sgsd2010)
sm2010<-merge(sm2010.raw, info2010, by=c("Logger", "Port"))

sm2011.raw<-rbind(knza2011, knzd2011)
sm2011<-merge(sm2011.raw, info2011, by=c("Logger","Port"))

sm2012.raw<-rbind(knza2012, knzd2012)
sm2012<-merge(sm2012.raw, info2012, by=c("Logger","Port"))

sm.raw<-rbind(sm2010, sm2011, sm2012)

sm.correct<-merge(probes_equations, sm.raw, by="Probe")%>%
  mutate(SoilMoist=vwc*slope+intercept)%>%
  mutate(date=as.POSIXct(Time, format="%m/%d/%Y %H:%M", tz="GMT"))%>%
  separate(date, into=c("day","time"), sep=" ")%>%
  mutate(plot.id=paste(Plot, Precip, sep="_"))%>%
  group_by(Probe, Precip, plot.id, day, year)%>%
  summarise(SoilMoist=mean(SoilMoist))

sm.plot<-sm.correct%>%
  mutate(rep=paste(plot.id, year, sep="_"))

# ggplot(data=sm.plot, aes(x=day, y=SoilMoist))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~rep, ncol=5, scales="free")

#drop bad probes and keep only 2010
sm.clean2010<-sm.correct%>%
  filter(Probe!="p211",Probe!="p206", Probe!="p204",Probe!="p203", Probe!="p205", Probe!="D7", Probe!="P13", year==2010)

sm.ave10<-sm.clean2010%>%
  group_by(Precip, day)%>%
  summarise(sm=mean(SoilMoist))%>%
  separate(day, into=c("year","month","day"), sep="-")%>%
  mutate(yr=2010,
         Day=paste(yr, month, day, sep="-"),
         Date=as.Date(Day))

smaverages10<-sm.ave10%>%
  group_by(Precip)%>%
  summarise(sm=mean(sm))

###2011
sm.clean11<-sm.correct%>%
  filter(Probe!="p211",Probe!="p206", Probe!="p204",Probe!="p203", Probe!="p205", Probe!="D7", Probe!="P13", year==2011, SoilMoist>0)

sm.ave11<-sm.clean11%>%
  group_by(Precip, day)%>%
  summarise(sm=mean(SoilMoist))%>%
  separate(day, into=c("year","month","day"), sep="-")%>%
  mutate(yr=2011,
         Day=paste(yr, month, day, sep="-"),
         Date=as.Date(Day))

smaverages11<-sm.ave11%>%
  group_by(Precip)%>%
  summarise(sm=mean(sm))

##2012 raw

sm.12nocor<-sm2012%>%
  filter(Probe!="p211",Probe!="p206", Probe!="p204",Probe!="p203", Probe!="p205", Probe!="D7", Probe!="P13")%>%
  mutate(date=as.POSIXct(Time, format="%m/%d/%Y %H:%M", tz="GMT"))%>%
  separate(date, into=c("day","time"), sep=" ")

sm.ave12nocor<-sm.12nocor%>%
  group_by(Precip, day)%>%
  summarise(sm=mean(vwc))%>%
  separate(day, into=c("year","month","day"), sep="-")%>%
  mutate(yr=2012,
         Day=paste(yr, month, day, sep="-"),
         Date=as.Date(Day))%>%
  filter(Date!="2012-09-27")

smaverages12<-sm.ave12nocor%>%
  group_by(Precip)%>%
  summarise(sm=mean(sm))

#Figures
fig<-sm.ave10%>%
  bind_rows(sm.ave11, sm.ave12nocor)


ggplot(data=fig, aes(x=Date, y=sm, group=Precip, color=Precip))+
  geom_point(aes(color=Precip))+
  scale_color_manual(name="Droughted", values=c("blue", "orange"), labels=c("No", "Yes"))+
  geom_line()+
  scale_y_continuous(name="Soil Moisture (VWC)",limits = c(-0.01,0.3))+
  #scale_x_date(name="", limits=as.Date(c("2010-05-01","2010-10-01")))+
  #annotate("text", x=as.Date("2010-05-15"), y = 0.05, label="2010", size= 6)+
  facet_wrap(~yr, scale="free_x", ncol=1)

ggsave("SoilMoisture.jpeg", plot=smTime, units = "in", width=9, height=6,dpi=300)

#repeated measures anova

sm.12nocor2<-sm.12nocor%>%
  group_by(Precip, Plot, Probe, day)%>%
  summarise(sm=mean(vwc))%>%
  mutate(year=2012)%>%
  filter(day!="2012-09-27")%>%
  rename(SoilMoist=sm)%>%
  mutate(plot.id=paste(Plot, Precip, sep="_"))%>%
  select(Probe, Precip, plot.id, day, year, SoilMoist)
  

small<-sm.clean2010%>%
  bind_rows(sm.clean11, sm.12nocor2)

# m<-lmer(SoilMoist~Precip*as.factor(year) + (1|plot.id), data=small)
# anova(m)

##overall effect


smaveragesall<-sm.ave12nocor%>%
  bind_rows(sm.ave11, sm.ave10)%>%
  group_by(Precip)%>%
  summarise(sm=mean(sm))



