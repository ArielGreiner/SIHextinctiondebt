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

#Biomass and Indiv Biomass on the same plot
ggplot(Biomass_Time[Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p],],aes(x=TimeStep,group=interaction(Scale, Patch_remove, Dispersal),fill=Scale, alpha = 0.1))+
  #geom_point()+ 
  geom_line(aes(y = Biomass, colour = "light blue")) + geom_line(aes(y = IndivBiomass, colour = "pink")) +   
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

#binning the biomass data because otherwise just see massively fluctuating sine waves
BiomassTime_Bin <- Biomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(SD_Biomass = sd(Mean_Biomass, na.rm = T), Mean_BiomassFinal = mean(Mean_Biomass, na.rm = T), Mean_IndivBiomassFinal = mean(Mean_IndivBiomass, na.rm = T), SD_IndivBiomass = sd(Mean_IndivBiomass, na.rm = T), Mean_SRFinal = mean(Mean_SR, na.rm = T), SD_SR = sd(Mean_SR, na.rm = T))

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

