#binning and summarizing SR_Time dataframe so as to calculate mean time to extinction and make the 'wave of extinction style plots'
MeanExtTimeBin <- SR_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/5)) %>%
  group_by(TimeStepRound,Dispersal,Patch_remove, Scale, Rep, Species, DelPatches)%>%
  summarize(Mean_SR = mean(SR, na.rm = T)) %>%
  group_by(Dispersal,Patch_remove, Scale, Rep, Species, DelPatches)%>%
  mutate(NumExt = lag(Mean_SR) - Mean_SR) %>%
  group_by(Dispersal, Patch_remove, Scale, TimeStepRound, Species, DelPatches) %>%
  summarize(Mean_NumExt = mean(NumExt, na.rm = T), SD_NumExt = sd(NumExt, na.rm = T))

#log'd version
#this is figure 3 (4/14/2016)
ggplot(MeanExtTimeBin[MeanExtTimeBin$Species == nSpeciesMult[s] & MeanExtTimeBin$DelPatches==nPatchDel[p],],aes(x=TimeStepRound*5,y=Mean_NumExt,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale,alpha = 0.1))+
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1,alpha = 0.1)+
  geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1, color = NA)+
  facet_grid(Dispersal~Patch_remove)+
  xlab("Time Step")+
  ylab("Mean Number of Extinctions")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  geom_vline(x=20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#graphing all of the SR and delpatch variants of figure 3, the forloop isn't working for some reason (did manually)
par(mfrow=c(length(nSpeciesMult),length(nPatchDel)))
for(s in 1:length(nSpeciesMult)){
  for(p in 1:length(nPatchDel)){
    ggplot(MeanExtTimeBin[MeanExtTimeBin$Species==nSpeciesMult[s] & MeanExtTimeBin$DelPatches==nPatchDel[p],],aes(x=TimeStepRound*5,y=Mean_NumExt,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale,alpha = 0.1))+
      geom_line()+
      scale_x_log10()+
      #geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1,alpha = 0.1)+
      geom_ribbon(aes(ymin=Mean_NumExt-SD_NumExt,ymax=Mean_NumExt+SD_NumExt),width=0.1, color = NA)+
      facet_grid(Dispersal~Patch_remove)+
      xlab("Time Step")+
      ylab("Mean Number of Extinctions")+
      ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
      #facet_grid(Dispersal~Patch_remove,scale="free")+
      #facet_grid(Scale~Patch_remove,scale="free")+
      geom_vline(x=20)+
      theme_bw(base_size = 18)+ #gets rid of grey background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  }
}

##Calculating proportion of species richness remaining, as calculated from the initial equilibrium number of species in the community
PropSR_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
    Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), 
    Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),
    Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), SR = NA)
#need to figure a better way, that doesn't involve forloops, to add in elements of the above dataframe
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(s in 1:length(nSpeciesMult)){
      for(p in 1:length(nPatchDel)){
        for(j in 1:reps){
          PropSR_Time$SR[PropSR_Time$Scale == "Regional" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Species == nSpeciesMult[s] & PropSR_Time$DelPatches == nPatchDel[p] & PropSR_Time$Rep == j]<- SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j]/SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][20]
          PropSR_Time$SR[PropSR_Time$Scale == "Local" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Species == nSpeciesMult[s] & PropSR_Time$DelPatches == nPatchDel[p] & PropSR_Time$Rep == j]<- SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j]/SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][20]
          
        } 
      }
    }
  }	
}

PropSRTimeSummd <- summarise(group_by(PropSR_Time, Dispersal, Patch_remove, TimeStep, Species, DelPatches, Scale), Mean_SR = mean(SR, na.rm=T), SD_SR = sd(SR, na.rm = T))

#this is figure 2 (4/14/2016) (mean proportional species richness over time, across all scenarios)
require(ggplot2)
#SR over time plots
ggplot(PropSRTimeSummd[PropSRTimeSummd$Species==nSpeciesMult[s]&PropSRTimeSummd$DelPatches==nPatchDel[p],],aes(x=TimeStep,y=Mean_SR,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Mean Proportion of Species Richness")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#adding in percent species loss metric into ED_data dataframe
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(j in 1:reps){
      for(s in 1:length(nSpeciesMult)){
        for(p in 1:length(nPatchDel)){
          Numat20 <- SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][20]
          
          ED_data$PercentLoss[ED_data$Scale == "Regional" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Numat20 - SR_Time$SR[SR_Time$Scale == "Regional" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][length(sampleV)])/Numat20
          
          Numat20 <- SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][20]
          
          ED_data$PercentLoss[ED_data$Scale == "Local" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Numat20 - SR_Time$SR[SR_Time$Scale == "Local" & SR_Time$Dispersal == dispV[o] & SR_Time$Patch_remove == removeV[w] & SR_Time$Species == nSpeciesMult[s] & SR_Time$DelPatches == nPatchDel[p] & SR_Time$Rep == j][length(sampleV)])/Numat20 
          
        }
      }
    }
  }	
}


EDdata_avg2 <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss = range(SRLoss, na.rm = T)[1], Highest_SRLoss = range(SRLoss, na.rm = T)[2],
                         Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T),Lowest_LastDebtTime = range(LastDebtTime, na.rm = T)[1], Highest_LastDebtTime = range(LastDebtTime, na.rm = T)[2], Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), 
                         Lowest_PercentLoss = range(PercentLoss, na.rm = T)[1], Highest_PercentLoss = range(PercentLoss, na.rm = T)[2])

#figure 4 (4/14/2016)
#percent of species lost vs time until last extinction plot
ggplot(EDdata_avg2[EDdata_avg2$Species == nSpeciesMult[s] & EDdata_avg2$DelPatches == nPatchDel[p],],aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Percentage of Species Lost")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/20,y=BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 2)+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Time Until Re-Equilibrium")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  #geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~., scales = "free_y")+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#this graph is an attempt at a line plot, but it really doesn't work very well
ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/20,y=BiomassChange,color=interaction(factor(Dispersal), Patch_remove),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+ #or "Paired"
  geom_line(aes(group = interaction(factor(Dispersal), Patch_remove)))+
  xlab("Time Until Re-Equilibrium")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  #geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#also doesn't really work :(
ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/20,y=BiomassChange,color=interaction(factor(Dispersal), Patch_remove),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+ #or "Paired"
  geom_point(aes(group = interaction(factor(Dispersal), Patch_remove)))+
  stat_smooth(method = "lm")+
  xlab("Time Until Re-Equilibrium")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  #geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

EDdata_avg <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T),
          Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T), Mean_Mean_Bmass = mean(Mean_Bmass, na.rm=T), SD_Mean_Bmass = sd(Mean_Bmass, na.rm=T),
          Mean_CVBmass = mean(CV_Bmass, na.rm=T), SD_CVBmass = sd(CV_Bmass, na.rm = T), Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=range(BiomassChange,na.rm=T)[1], Highest_BiomassChange=range(BiomassChange,na.rm=T)[2], Mean_LagTime = mean(LagTime, na.rm=T), SD_LagTime = sd(LagTime, na.rm=T), Lowest_LagTime=range(LagTime,na.rm=T)[1], Highest_LagTime=range(LagTime,na.rm=T)[2])

ggplot(EDdata_avg[EDdata_avg$Species == nSpeciesMult[s] & EDdata_avg$DelPatches == nPatchDel[p],],aes(x=Mean_LagTime,y=Mean_BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Time Until Re-Equilibrium")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Mean_BiomassChange - SD_BiomassChange,ymax=Mean_BiomassChange + SD_BiomassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Mean_LagTime - SD_LagTime,xmax=Mean_LagTime + SD_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.,scales = "free_y")+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(EDdata_avg[EDdata_avg$Species == nSpeciesMult[s] & EDdata_avg$DelPatches == nPatchDel[p],],aes(x=Mean_LagTime,y=Mean_BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Time Until Re-Equilibrium")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_BiomassChange,ymax= Highest_BiomassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LagTime,xmax=Highest_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.,scales = "free_y")+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


ggplot(EDdata_avg[EDdata_avg$Species == nSpeciesMult[s] & EDdata_avg$DelPatches == nPatchDel[p],],aes(x=Mean_Mean_Bmass,y=Mean_CVBmass,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Mean Biomass")+
  ylab("CV Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Mean_CVBmass - SD_CVBmass,ymax=Mean_CVBmass + SD_CVBmass),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Mean_Mean_Bmass - SD_Mean_Bmass,xmax=Mean_Mean_Bmass + SD_Mean_Bmass),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Scale~.,scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#looking at the proportion of biomass due to the different metacommunity processes

#binning the metadynamics data because otherwise just see massively fluctuating sine waves
MetaDynAvg_Bin <- Meta_dyn_reps %>%
  group_by(Dispersal, Patch_remove, Dynamic, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Dynamic, Species, DelPatches, Rep)%>%
  summarize(Mean_Proportion = mean(Proportion, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Dynamic, Species, DelPatches, TimeStepRound) %>%
  summarize(SD_Proportion = sd(Mean_Proportion, na.rm = T), Mean_Proportion = mean(Mean_Proportion, na.rm = T))

#supplementary materials figure 2 (4/3/2016)
ggplot(MetaDynAvg_Bin[MetaDynAvg_Bin$Species == nSpeciesMult[s] & MetaDynAvg_Bin$DelPatches == nPatchDel[p],],aes(x=TimeStepRound,y=Mean_Proportion,color=Dynamic,fill = Dynamic))+
  xlab("Time Step")+
  ylab("Proportion of Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/20)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  geom_ribbon(aes(ymin=Mean_Proportion-SD_Proportion,ymax=Mean_Proportion+SD_Proportion), alpha = 0.2, color = NA)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

##Plotting Biomass
require(ggplot2)
#Raw Biomass Plot
ggplot(Biomass_Time[Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p],],aes(x=TimeStep,y=Biomass,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(Biomass_Time[Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Scale == "Regional",],aes(x=TimeStep,y=CVTime,group=interaction(Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("CV Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Biomass, SR, Indiv Biomass on the same plot 
ggplot(Biomass_Time[Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Scale == "Local" & Biomass_Time$Rep == r,],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Biomass/10, colour = "green")) + geom_line(aes(y = SR, colour = "light blue")) + geom_line(aes(y = IndivBiomass/10, colour = "pink"))  +
  #divide biomass and indivbiomass by 100 if using
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

Biomass_TimeSummd <- summarise(group_by(Biomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), SD_SR = sd(SR, na.rm = T), Mean_Biomass = mean(Biomass, na.rm=T), SD_Biomass = sd(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), SD_IndivBiomass = sd(IndivBiomass, na.rm = T), Mean_CVTime = mean(CVTime, na.rm=T), SD_CVTime = sd(CVTime, na.rm = T))

ggplot(Biomass_TimeSummd[Biomass_TimeSummd$Species == nSpeciesMult[s] & Biomass_TimeSummd$DelPatches == nPatchDel[p] & Biomass_TimeSummd$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/100")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Average Biomass per Species /100"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, fill = "blue", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=(Mean_Biomass/10)-(SD_Biomass/10),ymax=(Mean_Biomass/10)+(SD_Biomass/10)),width=0.1, fill = "green", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=(Mean_IndivBiomass/10)-(SD_IndivBiomass/10),ymax=(Mean_IndivBiomass/10)+(SD_IndivBiomass/10)),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(Biomass_TimeSummd[Biomass_TimeSummd$Species == nSpeciesMult[s] & Biomass_TimeSummd$DelPatches == nPatchDel[p],],aes(x=TimeStep,y=Mean_CVTime,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_CVTime-SD_CVTime,ymax=Mean_CVTime+SD_CVTime),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("CV Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20)+ #makes sense because now the first 20 points are just the CV of the last period before patch deletion
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#binning the biomass data because otherwise just see massively fluctuating sine waves
BiomassTime_Bin <- Biomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(SD_Biomass = sd(Mean_Biomass, na.rm = T), Mean_BiomassFinal = mean(Mean_Biomass, na.rm = T), Mean_IndivBiomassFinal = mean(Mean_IndivBiomass, na.rm = T), SD_IndivBiomass = sd(Mean_IndivBiomass, na.rm = T), Mean_SRFinal = mean(Mean_SR, na.rm = T), SD_SR = sd(Mean_SR, na.rm = T), Mean_CVTimeFinal = mean(Mean_CVTime, na.rm=T), SD_CVTime = sd(Mean_CVTime, na.rm=T))

ggplot(BiomassTime_Bin[BiomassTime_Bin$Species == nSpeciesMult[s] & BiomassTime_Bin$DelPatches == nPatchDel[p] & BiomassTime_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_BiomassFinal/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_SRFinal, colour = "Species Richness")) + geom_line(aes(y = Mean_IndivBiomassFinal/10, colour = "avg Biomass per Species /10"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  #geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1, color = NA)+
  xlab("Time Step")+
  #ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(BiomassTime_Bin[BiomassTime_Bin$Species == nSpeciesMult[s] & BiomassTime_Bin$DelPatches == nPatchDel[p] & BiomassTime_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_BiomassFinal/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_SRFinal, colour = "Species Richness")) + geom_line(aes(y = Mean_IndivBiomassFinal/10, colour = "Biomass per Species /10"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_SRFinal-SD_SR,ymax=Mean_SRFinal+SD_SR),width=0.1, color = NA, alpha = 0.1, fill = "blue")+
  geom_ribbon(aes(ymin=(Mean_BiomassFinal/10)-(SD_Biomass/10),ymax=(Mean_BiomassFinal/10)+(SD_Biomass/10)),width=0.1,alpha = 0.1, fill = "green", color = NA)+    
  geom_ribbon(aes(ymin=(Mean_IndivBiomassFinal/10)-(SD_Biomass/10),ymax=(Mean_IndivBiomassFinal/10)+(SD_Biomass/10)),width=0.1,alpha = 0.1, fill = "red", color = NA)+
  xlab("Time Step")+
  #ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(BiomassTime_Bin[BiomassTime_Bin$Species == nSpeciesMult[s] & BiomassTime_Bin$DelPatches == nPatchDel[p],],aes(x=TimeStepRound,y=Mean_CVTimeFinal,color=Scale,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line()+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_CVTimeFinal-SD_CVTime,ymax=Mean_CVTimeFinal+SD_CVTime),width=0.1, color = NA)+
  xlab("Time Step")+
  ylab("CV Biomass (Binned)")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=20/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
 
 #[2016-06-01, 3:33:50 PM] Patrick Thompson: group_by()%>%
#[2016-06-01, 3:34:11 PM] Patrick Thompson:mutate(Stand=decostand(.,method="standardize")
#BiomassTime_Stand <- Biomass_Time %>%
  #group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  #mutate(StandBiomass = decostand(Biomass,method="standardize"),StandIndivBiomass = decostand(IndivBiomass,method="standardize"), StandSR = decostand(SR,method="standardize"), StandCVTime = decostand(CVTime,method="standardize", na.rm=T)) %>%
 #^ doesn't work, would also need to add a summarizing thing after this within the pipeline

BiomassTime_Stand <- Biomass_Time
BiomassTime_Stand$StandSR <- (BiomassTime_Stand$SR - mean(BiomassTime_Stand$SR))/sd(BiomassTime_Stand$SR)
BiomassTime_Stand$StandBiomass <- (BiomassTime_Stand$Biomass - mean(BiomassTime_Stand$Biomass))/sd(BiomassTime_Stand$Biomass)
BiomassTime_Stand$StandIndivBiomass <- (BiomassTime_Stand$IndivBiomass - mean(BiomassTime_Stand$IndivBiomass))/sd(BiomassTime_Stand$IndivBiomass)
BiomassTime_Stand$StandCVTime <- (BiomassTime_Stand$CVTime - mean(BiomassTime_Stand$CVTime, na.rm=T))/sd(BiomassTime_Stand$CVTime, na.rm=T)
#BiomassTime_Stand$StandSR <- decostand(BiomassTime_Stand$SR,method="standardize")
#BiomassTime_Stand$StandBiomass <- decostand(BiomassTime_Stand$Biomass, method = "standardize")
#BiomassTime_Stand$StandIndivBiomass <- decostand(BiomassTime_Stand$IndivBiomass, method = "standardize")
#BiomassTime_Stand$StandCVTime <- decostand(BiomassTime_Stand$CVTime, method = "standardize",na.rm=T)
#^ didn't use decostand because it introduces weird formatting that was making other things fail

BiomassTime_StandSummd <- summarise(group_by(BiomassTime_Stand,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_StandSR = mean(StandSR, na.rm=T), SD_StandSR = sd(StandSR, na.rm = T), Mean_StandBiomass = mean(StandBiomass, na.rm=T), SD_StandBiomass = sd(StandBiomass, na.rm = T), Mean_StandIndivBiomass = mean(StandIndivBiomass, na.rm=T), SD_StandIndivBiomass = sd(StandIndivBiomass, na.rm = T), Mean_StandCVTime = mean(StandCVTime, na.rm=T), SD_StandCVTime = sd(StandCVTime, na.rm = T))

ggplot(BiomassTime_StandSummd[BiomassTime_StandSummd$Species == nSpeciesMult[s] & BiomassTime_StandSummd$DelPatches == nPatchDel[p] & BiomassTime_StandSummd$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_StandBiomass, colour = "Biomass")) + geom_line(aes(y = Mean_StandSR, colour = "Species Richness")) + geom_line(aes(y = Mean_StandIndivBiomass, colour = "Average Biomass per Species"))  +  geom_line(aes(y = Mean_StandCVTime, colour = "CV Biomass"))+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_StandSR-SD_StandSR,ymax=Mean_StandSR+SD_StandSR),width=0.1, fill = "purple", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=(Mean_StandBiomass)-(SD_StandBiomass),ymax=(Mean_StandBiomass)+(SD_StandBiomass)),width=0.1, fill = "green", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=(Mean_StandIndivBiomass)-(SD_StandIndivBiomass),ymax=(Mean_StandIndivBiomass)+(SD_StandIndivBiomass)),width=0.1, fill = "red", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Mean_StandCVTime-SD_StandCVTime,ymax=Mean_StandCVTime+SD_StandCVTime),width=0.1, fill = "cyan", alpha = 0.2, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  
    
  BiomassTimeStand_Bin <- BiomassTime_Stand %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_StandBiomass = mean(StandBiomass, na.rm = T), Mean_StandIndivBiomass = mean(StandIndivBiomass, na.rm = T), Mean_StandSR = mean(StandSR, na.rm = T), Mean_StandCVTime = mean(StandCVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(SD_StandBiomass = sd(Mean_StandBiomass, na.rm = T), Mean_StandBiomassFinal = mean(Mean_StandBiomass, na.rm = T), Mean_StandIndivBiomassFinal = mean(Mean_StandIndivBiomass, na.rm = T), SD_StandIndivBiomass = sd(Mean_StandIndivBiomass, na.rm = T), Mean_StandSRFinal = mean(Mean_StandSR, na.rm = T), SD_StandSR = sd(Mean_StandSR, na.rm = T), Mean_StandCVTimeFinal = mean(Mean_StandCVTime, na.rm=T), SD_StandCVTime = sd(Mean_StandCVTime, na.rm=T))

ggplot(BiomassTimeStand_Bin[BiomassTimeStand_Bin$Species == nSpeciesMult[s] & BiomassTimeStand_Bin$DelPatches == nPatchDel[p] & BiomassTimeStand_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_StandBiomassFinal, colour = "Biomass")) + geom_line(aes(y = Mean_StandSRFinal, colour = "Species Richness")) + geom_line(aes(y = Mean_StandIndivBiomassFinal, colour = "Biomass per Species"))  + geom_line(aes(y = Mean_StandCVTimeFinal, colour = "CV Biomass")) +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Mean_StandSRFinal-SD_StandSR,ymax=Mean_StandSRFinal+SD_StandSR),width=0.1, color = NA, alpha = 0.1, fill = "purple")+
  geom_ribbon(aes(ymin=(Mean_StandBiomassFinal)-(SD_StandBiomass),ymax=(Mean_StandBiomassFinal)+(SD_StandBiomass)),width=0.1,alpha = 0.1, fill = "red", color = NA)+    
  geom_ribbon(aes(ymin=(Mean_StandIndivBiomassFinal)-(SD_StandBiomass),ymax=(Mean_StandIndivBiomassFinal)+(SD_StandBiomass)),width=0.1,alpha = 0.1, fill = "green", color = NA)+
  geom_ribbon(aes(ymin=Mean_StandCVTimeFinal-SD_StandCVTime,ymax=Mean_StandCVTimeFinal+SD_StandCVTime),width=0.1,alpha = 0.1, fill = "cyan", color = NA)+
  xlab("Time Step")+
  #ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local Scale"))+
  geom_vline(x=20/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
  BiomassTime_StandSummd2 <- summarise(group_by(BiomassTime_Stand,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_StandSR = mean(StandSR, na.rm=T), Upper_StandSR = quantile(StandSR, probs=.975, na.rm = T, names = F),Lower_StandSR = quantile(StandSR, probs=.025, na.rm = T, names = F), Mean_StandBiomass = mean(StandBiomass, na.rm=T), Upper_StandBiomass = quantile(StandBiomass, probs=.975, na.rm = T, names = F), Lower_StandBiomass = quantile(StandBiomass, probs=.025, na.rm = T, names = F), Mean_StandIndivBiomass = mean(StandIndivBiomass, na.rm=T), Upper_StandIndivBiomass = quantile(StandIndivBiomass, probs=.975, na.rm = T, names = F),Lower_StandIndivBiomass = quantile(StandIndivBiomass, probs=.025, na.rm = T, names = F), Mean_StandCVTime = mean(StandCVTime, na.rm=T), Upper_StandCVTime = quantile(StandCVTime, probs=.975, na.rm = T, names = F), Lower_StandCVTime = quantile(StandCVTime, probs=.025, na.rm = T, names = F))

ggplot(BiomassTime_StandSummd2[BiomassTime_StandSummd2$Species == nSpeciesMult[s] & BiomassTime_StandSummd2$DelPatches == nPatchDel[p] & BiomassTime_StandSummd2$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_StandBiomass, colour = "Biomass")) + geom_line(aes(y = Mean_StandSR, colour = "Species Richness")) + geom_line(aes(y = Mean_StandIndivBiomass, colour = "Average Biomass per Species"))  +  geom_line(aes(y = Mean_StandCVTime, colour = "CV Biomass"))+
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_StandSR,ymax=Upper_StandSR),width=0.1, fill = "purple", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_StandBiomass,ymax=Upper_StandBiomass),width=0.1, fill = "green", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_StandIndivBiomass,ymax=Upper_StandIndivBiomass),width=0.1, fill = "red", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_StandCVTime,ymax=Upper_StandCVTime),width=0.1, fill = "cyan", alpha = 0.2, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  

#Upper_StandSR = quantile(StandSR, probs=.975, na.rm = T, names = F)
#Lower_StandSR = quantile(StandSR, probs=.025, na.rm = T, names = F)

  BiomassTimeStand_Bin2 <- BiomassTime_Stand %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_StandBiomass = mean(StandBiomass, na.rm = T), Mean_StandIndivBiomass = mean(StandIndivBiomass, na.rm = T), Mean_StandSR = mean(StandSR, na.rm = T), Mean_StandCVTime = mean(StandCVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Upper_StandBiomass = quantile(Mean_StandBiomass, probs=0.975, na.rm = T,names=F),Lower_StandBiomass = quantile(Mean_StandBiomass, probs=0.025, na.rm = T,names=F), Mean_StandBiomassFinal = mean(Mean_StandBiomass, na.rm = T), Mean_StandIndivBiomassFinal = mean(Mean_StandIndivBiomass, na.rm = T), Upper_StandIndivBiomass = quantile(Mean_StandIndivBiomass, probs=0.975, na.rm = T,names=F),Lower_StandIndivBiomass = quantile(Mean_StandIndivBiomass, probs=0.025, na.rm = T,names=F), Mean_StandSRFinal = mean(Mean_StandSR, na.rm = T), Upper_StandSR = quantile(Mean_StandSR, probs=0.975, na.rm = T,names=F),Lower_StandSR = quantile(Mean_StandSR, probs=0.025, na.rm = T,names=F), Mean_StandCVTimeFinal = mean(Mean_StandCVTime, na.rm=T), Upper_StandCVTime = quantile(Mean_StandCVTime, probs=0.975, na.rm = T,names=F),Lower_StandCVTime = quantile(Mean_StandCVTime, probs=0.025, na.rm = T,names=F))

ggplot(BiomassTimeStand_Bin2[BiomassTimeStand_Bin2$Species == nSpeciesMult[s] & BiomassTimeStand_Bin2$DelPatches == nPatchDel[p] & BiomassTimeStand_Bin2$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_StandBiomassFinal, colour = "Biomass")) + geom_line(aes(y = Mean_StandSRFinal, colour = "Species Richness")) + geom_line(aes(y = Mean_StandIndivBiomassFinal, colour = "Biomass per Species"))  + geom_line(aes(y = Mean_StandCVTimeFinal, colour = "CV Biomass")) +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_StandSR,ymax=Upper_StandSR),width=0.1, color = NA, alpha = 0.1, fill = "purple")+
  geom_ribbon(aes(ymin=Lower_StandBiomass,ymax=Upper_StandBiomass),width=0.1,alpha = 0.1, fill = "red", color = NA)+    
geom_ribbon(aes(ymin=Lower_StandIndivBiomass,ymax=Upper_StandIndivBiomass),width=0.1,alpha = 0.1, fill = "green", color = NA)+
  geom_ribbon(aes(ymin=Lower_StandCVTime,ymax=Upper_StandCVTime),width=0.1,alpha = 0.1, fill = "cyan", color = NA)+
  xlab("Time Step")+
  #ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local Scale"))+
  geom_vline(x=20/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
###haven't adjusted the plots below here to deal with the extra levels of species richness and patch deletion, e.g. EffDivFinal moved
ggplot(BiomassTime_Bin,aes(x=TimeStepRound,y=Mean_BiomassFinal,color=Scale,fill = Scale))+
  xlab("(Binned) Time Step")+
  ylab("Biomass")+
  geom_line()+
  scale_x_log10()+
  facet_grid(Dispersal~Patch_remove)+	  
  geom_vline(x=20/20)+
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

