#binning and summarizing SR_Time dataframe so as to calculate mean time to extinction and make the 'wave of extinction style plots'
MeanExtTimeBin <- SR_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/5)) %>%
  group_by(TimeStepRound,Dispersal,Patch_remove, Scale, Rep)%>%
  summarize(Mean_SR = mean(SR, na.rm = T)) %>%
  group_by(Dispersal,Patch_remove, Scale, Rep)%>%
  mutate(NumExt = lag(Mean_SR) - Mean_SR) %>%
  group_by(Dispersal, Patch_remove, Scale, TimeStepRound) %>%
  summarize(Mean_NumExt = mean(NumExt, na.rm = T), SD_NumExt = sd(NumExt, na.rm = T))

ggplot(MeanExtTimeBin,aes(x=(TimeStepRound)*5,y=Mean_NumExt,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale,alpha = 0.1))+
  geom_line()+
  #geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1,alpha = 0.1)+
  geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1, color = NA)+
  facet_grid(Dispersal~Patch_remove)+
  xlab("Time Step")+
  ylab("Mean Number of Extinctions")+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  geom_vline(x=20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#log'd version
#this is figure 3 (4/14/2016)
ggplot(MeanExtTimeBin,aes(x=TimeStepRound*5,y=Mean_NumExt,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale,alpha = 0.1))+
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1,alpha = 0.1)+
  geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1, color = NA)+
  facet_grid(Dispersal~Patch_remove)+
  xlab("Time Step")+
  ylab("Mean Number of Extinctions")+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  geom_vline(x=20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#unbinned version of the above 
MeanExtTime <- SR_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Rep) %>%
  mutate(NumExt = lag(SR) - SR) %>%
  group_by(Dispersal, Patch_remove, Scale, TimeStep) %>%
  summarize(Mean_NumExt = mean(NumExt, na.rm = T), SD_NumExt = sd(NumExt, na.rm = T))

ggplot(MeanExtTime,aes(x=TimeStep,y=Mean_NumExt,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale,alpha = 0.1))+
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1,alpha = 0.1)+
  geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1, color = NA)+
  facet_grid(Dispersal~Patch_remove)+
  xlab("Time Step (unbinned)")+
  ylab("Mean Number of Extinctions")+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  geom_vline(x=20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

##Calculating proportion of species richness remaining, as calculated from the initial equilibrium number of species in the community
PropSR_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2),
                          Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2),
                          Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2),
                          TimeStep = rep(1:length(sampleV)),Scale=rep(c("Local","Regional"), each = length(sampleV)), SR = NA)

#need to figure a better way, that doesn't involve forloops, to add in elements of the above dataframe
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(j in 1:reps){
      PropSR_Time$SR[PropSR_Time$Scale == "Regional" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Rep == j]<- SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j]/SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][20]
      PropSR_Time$SR[PropSR_Time$Scale == "Local" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Rep == j]<- SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j]/SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][20]
      
    }
  }	
}

PropSRTimeSummd <- summarise(group_by(PropSR_Time, Dispersal, Patch_remove, TimeStep, Scale), Mean_SR = mean(SR, na.rm=T), SD_SR = sd(SR, na.rm = T))

#this is figure 2 (4/14/2016) (mean proportional species richness over time, across all scenarios)
require(ggplot2)
#SR over time plots
ggplot(PropSRTimeSummd,aes(x=TimeStep,y=Mean_SR,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Mean Proportion of Species Richness")+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


#for making species richness over time plot
SRTimeSummd <- summarise(group_by(SR_Time, Dispersal, Patch_remove, TimeStep, Scale), Mean_SR = mean(SR, na.rm=T), SD_SR = sd(SR, na.rm = T))

require(ggplot2)
#SR over time plots
ggplot(SRTimeSummd,aes(x=TimeStep,y=Mean_SR,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Mean Species Richness")+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

EDdata_avg <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), 
                        Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T))

#number of species lost vs time until last extinction plot, split into local and regional plots
ggplot(EDdata_avg,aes(x=Mean_LastDebtTime,y=Mean_SRLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+
  scale_color_brewer("Dispersal Level", palette = "BrBG")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  scale_shape_manual(values=c(15,19, 17))+ #25 = upside-down triangle
  xlab("Time Until Last Extinction")+
  ylab("Number of Species Lost")+
  geom_errorbar(aes(ymin=Mean_SRLoss-SD_SRLoss,ymax=Mean_SRLoss+SD_SRLoss),width=0.1)+
  geom_errorbarh(aes(xmin=Mean_LastDebtTime-SD_LastDebtTime,xmax=Mean_LastDebtTime+SD_LastDebtTime),width=0.1)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#adding in percent species loss metric into ED_data dataframe
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(j in 1:reps){
      Numat20 <- SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][20]
      
      ED_data$PercentLoss[ED_data$Scale == "Regional" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Rep == j]<- (Numat20 - SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][length(sampleV)])/Numat20
      
      Numat20 <- SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][20]
      
      ED_data$PercentLoss[ED_data$Scale == "Local" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Rep == j]<- (Numat20 - SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Rep == j][length(sampleV)])/Numat20 
    }
  }	
}


EDdata_avg2 <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss = range(SRLoss, na.rm = T)[1], Highest_SRLoss = range(SRLoss, na.rm = T)[2],
                         Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T),Lowest_LastDebtTime = range(LastDebtTime, na.rm = T)[1], Highest_LastDebtTime = range(LastDebtTime, na.rm = T)[2], Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), 
                         Lowest_PercentLoss = range(PercentLoss, na.rm = T)[1], Highest_PercentLoss = range(PercentLoss, na.rm = T)[2])

#percent of species lost vs time until last extinction plot
ggplot(EDdata_avg2,aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ #alpha = Scale
  scale_color_brewer("Dispersal Level", palette = "BrBG")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Percentage of Species Lost")+
  geom_errorbar(aes(ymin=Mean_PercentLoss-SD_PercentLoss,ymax=Mean_PercentLoss+SD_PercentLoss),width=0.1)+
  geom_errorbarh(aes(xmin=Mean_LastDebtTime-SD_LastDebtTime,xmax=Mean_LastDebtTime+SD_LastDebtTime),width=0.1)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#figure 4 (4/14/2016)
#percent of species lost vs time until last extinction plot
ggplot(EDdata_avg2,aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Percentage of Species Lost")+
  geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Supplementary Material Figure 1 (4.15.2016)
ggplot(EDdata_avg2,aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal), alpha = Scale))+ 
  #alpha = Scale +
  scale_color_brewer("Dispersal Level", palette = "BrBG")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Percentage of Species Lost")+
  geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  #facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines



#looking at the proportion of biomass due to the different metacommunity processes

Metadyn_avg <- summarise(group_by(Meta_dyn_reps,Dispersal,Patch_remove,TimeStep,Dynamic), Mean_Proportion = mean(Proportion, na.rm=T), SD_Proportion = sd(Proportion, na.rm = T))

ggplot(Metadyn_avg,aes(x=TimeStep,y=Mean_Proportion,color=factor(Dynamic),group=interaction(Dynamic, Patch_remove, Dispersal),alpha = 0.1))+
  #scale_color_brewer("Process", palette = "BrBG")+
  xlab("Time Step")+
  ylab("Proportion of Biomass")+
  geom_vline(x=20)+
  geom_ribbon(aes(ymin=Mean_Proportion-SD_Proportion,ymax=Mean_Proportion+SD_Proportion),width=0.1)+
  facet_grid(Patch_remove~Dispersal)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#binning the metadynamics data because otherwise just see massively fluctuating sine waves
MetaDynAvg_Bin <- Meta_dyn_reps %>%
  group_by(Dispersal, Patch_remove, Dynamic, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Dynamic, Rep)%>%
  summarize(Mean_Proportion = mean(Proportion, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Dynamic, TimeStepRound) %>%
  summarize(SD_Proportion = sd(Mean_Proportion, na.rm = T), Mean_Proportion = mean(Mean_Proportion, na.rm = T))

#supplementary materials figure 2 (4/3/2016)
ggplot(MetaDynAvg_Bin,aes(x=TimeStepRound,y=Mean_Proportion,color=Dynamic,fill = Dynamic))+
  xlab("Time Step")+
  ylab("Proportion of Biomass")+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  geom_ribbon(aes(ymin=Mean_Proportion-SD_Proportion,ymax=Mean_Proportion+SD_Proportion), alpha = 0.2, color = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

##Individual Patch Level Stuff
ggplot(IndivPatch,aes(x=Betweenness,y=LastExtTime,color=factor(Rep),group=interaction(Patch_remove, Dispersal, Rep)))+ 
  #scale_color_brewer("Dispersal Level", palette = "BrBG")+
  #geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  geom_point()+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  #scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Betweenness")+
  ylab("Time Until Last Extinction")+
  #geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Dispersal~Patch_remove)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

##Plotting Biomass
require(ggplot2)
#Raw Biomass Plot
ggplot(Biomass_Time,aes(x=TimeStep,y=Biomass,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#binning the biomass data because otherwise just see massively fluctuating sine waves
BiomassTime_Bin <- Biomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/5)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_EffDiv = mean(EffDiv, na.rm = T), Mean_EffDivBetaMult = mean(EffDivBetaMult, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Scale, TimeStepRound) %>%
  summarize(SD_Biomass = sd(Mean_Biomass, na.rm = T), Mean_BiomassFinal = mean(Mean_Biomass, na.rm = T), Mean_EffDivFinal = mean(Mean_EffDiv, na.rm = T), SD_EffDiv = sd(Mean_EffDiv, na.rm = T), Mean_EffDivBetaMultFinal = mean(Mean_EffDivBetaMult, na.rm = T), SD_EffDivBetaMult = sd(Mean_EffDivBetaMult, na.rm = T))

ggplot(BiomassTime_Bin,aes(x=TimeStepRound,y=Mean_BiomassFinal,color=Scale,fill = Scale))+
  xlab("(Binned) Time Step")+
  ylab("Biomass")+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/5)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  geom_ribbon(aes(ymin=Mean_BiomassFinal-SD_Biomass,ymax=Mean_BiomassFinal+SD_Biomass), alpha = 0.2, color = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(BiomassTime_Bin,aes(x=TimeStepRound,y=Mean_EffDivFinal,color=Scale,fill = Scale))+
  xlab("(Binned) Time Step")+
  ylab("Effective Diversity")+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/5)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  geom_ribbon(aes(ymin=Mean_EffDivFinal-SD_EffDiv,ymax=Mean_EffDivFinal+SD_EffDiv), alpha = 0.2, color = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Biomass_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2),
#Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2),
#Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2),
#TimeStep = rep(1:length(sampleV)),Scale=rep(c("Local","Regional"), each = length(sampleV)), 
#EffDiv = NA, EffDivBetaMult = NA, EffDivBetaAdd = NA, Biomass = NA)

#EffectiveDiv_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*3),
#Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*3),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*3),
#TimeStep = rep(1:length(sampleV)),Metric=rep(c("Alpha","Gamma","Beta"), each = length(sampleV)), ExpShannon = NA)

EffectiveDiv_Bin <- EffectiveDiv_Time %>%
  group_by(Dispersal, Patch_remove, Metric, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/5)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Metric, Rep)%>%
  summarize(Mean_ExpShannon = mean(ExpShannon, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Metric, TimeStepRound) %>%
  summarize(SD_ExpShannon = sd(Mean_ExpShannon, na.rm = T), Mean_ExpShannonFinal = mean(Mean_ExpShannon, na.rm = T))

ggplot(EffectiveDiv_Bin,aes(x=TimeStepRound,y=Mean_ExpShannonFinal,color=Metric,fill = Metric))+
  xlab("(Binned) Time Step")+
  ylab("Effective Diversity")+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/5)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  geom_ribbon(aes(ymin=Mean_ExpShannonFinal-SD_ExpShannon,ymax=Mean_ExpShannonFinal+SD_ExpShannon), alpha = 0.2, color = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


#Fourier Transform Plot (plotted for 1 replicate)
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main=paste("Dispersal Level", dispV[i], removeV[j]), 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

par(mfrow=c(length(dispV),length(removeV)))
for(i in 1:length(dispV)){
  for(j in 1:length(removeV)){
    plot.frequency.spectrum(fft(Biomass_Time$Biomass[Biomass_Time$Rep==r & Biomass_Time$Dispersal==dispV[i] & Biomass_Time$Patch_remove==removeV[j] & Biomass_Time$Scale=="Local"]))
  }
}

###old plots

###Tester Plots
#will plot the local biomass of each individual patch for the last scenario run
plot(L_Bmass_sep$X30, type = 'l')
#7, 11, 18, 24, 30
plot(L_Bmass_sep$X30, type = 'l', xlab = "Time Step",ylab = "Biomass", main = paste("Biomass of Patch", sep = " ",30, "over time [Dispersal = ",dispV[i], ", removal sequence = ", removeV[j], "]"))

