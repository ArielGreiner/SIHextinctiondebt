#binning and summarizing SR_Time dataframe so as to calculate mean time to extinction and make the 'wave of extinction style plots'
MeanExtTimeBin <- Biomass_Time %>%
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
  geom_vline(x=predel_collecttime)+
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
      geom_vline(x=predel_collecttime)+
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
          PropSR_Time$SR[PropSR_Time$Scale == "Regional" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Species == nSpeciesMult[s] & PropSR_Time$DelPatches == nPatchDel[p] & PropSR_Time$Rep == j]<- Biomass_Time$SR[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j]/Biomass_Time$SR[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          PropSR_Time$SR[PropSR_Time$Scale == "Local" & PropSR_Time$Dispersal == dispV[o] & PropSR_Time$Patch_remove == removeV[w] & PropSR_Time$Species == nSpeciesMult[s] & PropSR_Time$DelPatches == nPatchDel[p] & PropSR_Time$Rep == j]<- Biomass_Time$SR[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j]/Biomass_Time$SR[SR_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          
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
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#adding in percent species loss metric into ED_data dataframe (and percent change in biomass, percent change in CV)
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(j in 1:reps){
      for(s in 1:length(nSpeciesMult)){
        for(p in 1:length(nPatchDel)){
          Numpredel <- Biomass_Time$SR[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          Biomasspredel <- Biomass_Time$Biomass[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          CVpredel <- Biomass_Time$CVTime[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          
          ED_data$PercentLoss[ED_data$Scale == "Regional" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Numpredel - Biomass_Time$SR[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV)])/Numpredel
          #^ i feel like this line could be replaced with ED_data$SRLoss[...]/Numpredel now that the SRLoss metric only looks at what's going on after patch deletion though not sure about this
          ED_data$PercentBmasschange[ED_data$Scale == "Regional" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Biomasspredel - Biomass_Time$Biomass[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV)])/Biomasspredel
          ED_data$PercentCVchange[ED_data$Scale == "Regional" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (CVpredel - Biomass_Time$CVTime[Biomass_Time$Scale == "Regional" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV) - ePeriod/samplelength])/CVpredel
          
          
          Numpredel <- Biomass_Time$SR[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          Biomasspredel <- Biomass_Time$Biomass[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          CVpredel <- Biomass_Time$CVTime[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][predel_collecttime]
          
          ED_data$PercentLoss[ED_data$Scale == "Local" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Numpredel - Biomass_Time$SR[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV)])/Numpredel 
          #^ i feel like this line could be replaced with ED_data$SRLoss[...]/Numpredel now that the SRLoss metric only looks at what's going on after patch deletion though not sure about this
          ED_data$PercentBmasschange[ED_data$Scale == "Local" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (Biomasspredel - Biomass_Time$Biomass[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV)])/Biomasspredel
          ED_data$PercentCVchange[ED_data$Scale == "Local" & ED_data$Dispersal == dispV[o] & ED_data$Patch_remove == removeV[w] & ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p] & ED_data$Rep == j]<- (CVpredel - Biomass_Time$CVTime[Biomass_Time$Scale == "Local" & Biomass_Time$Dispersal == dispV[o] & Biomass_Time$Patch_remove == removeV[w] & Biomass_Time$Species == nSpeciesMult[s] & Biomass_Time$DelPatches == nPatchDel[p] & Biomass_Time$Rep == j][length(sampleV) - ePeriod/samplelength])/CVpredel
          
        }
      }
    }
  }	
}


EDdata_avg2 <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss = range(SRLoss, na.rm = T)[1], Highest_SRLoss = range(SRLoss, na.rm = T)[2],
Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T),Lowest_LastDebtTime = range(LastDebtTime, na.rm = T)[1], Highest_LastDebtTime = range(LastDebtTime, na.rm = T)[2], Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), Lowest_PercentLoss = range(PercentLoss, na.rm = T)[1], Highest_PercentLoss = range(PercentLoss, na.rm = T)[2])

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

#summarizing all of the different metrics (biomass, cv, sr) to capture the mean percent change/loss and their 95% quantiles   
EDavg <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), 
Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), Lowest_PercentLoss = quantile(PercentLoss, probs = 0.025, na.rm=T, names = F), Highest_PercentLoss = quantile(PercentLoss, probs = 0.975, na.rm=T, names = F),
Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T), Lowest_LastDebtTime = quantile(LastDebtTime, probs = 0.025, na.rm=T, names = F), Highest_LastDebtTime = quantile(LastDebtTime, probs = 0.975, na.rm=T, names = F),
Mean_PercentBmassChange = mean(PercentBmassChange, na.rm=T), SD_PercentBmassChange = sd(PercentBmassChange, na.rm = T), Lowest_PercentBmassChange = quantile(PercentBmassChange, probs = 0.025, na.rm=T, names = F), Highest_PercentBmassChange = quantile(PercentBmassChange, probs = 0.975, na.rm=T, names = F),
Mean_LagTime = mean(LagTime, na.rm=T), SD_LagTime = sd(LagTime, na.rm=T), Lowest_LagTime=quantile(LagTime, probs = 0.025, na.rm=T, names = F), Highest_LagTime=quantile(LagTime, probs = 0.975, na.rm=T, names = F),
Mean_PercentCVChange = mean(PercentCVChange, na.rm=T), SD_PercentCVChange = sd(PercentCVChange, na.rm = T), Lowest_PercentCVChange = quantile(PercentCVChange, probs = 0.025, na.rm=T, names = F), Highest_PercentCVChange = quantile(PercentCVChange, probs = 0.975, na.rm=T, names = F),
Mean_CVLagTime = mean(CVLagTime, na.rm=T), SD_CVLagTime = sd(CVLagTime, na.rm=T), Lowest_CVLagTime=quantile(CVLagTime, probs = 0.025, na.rm=T, names = F), Highest_CVLagTime=quantile(CVLagTime, probs = 0.975, na.rm=T, names = F))

#SR plot
ggplot(EDavg[EDavg$Species == nSpeciesMult[s] & EDavg$DelPatches == nPatchDel[p],],aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
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

#SR tryptic-style plot
ggplot(EDavg[EDavg$Species == nSpeciesMult[s],],aes(x=Mean_LastDebtTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Percentage of Species Lost")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Biomass plot
ggplot(EDavg[EDavg$Species == nSpeciesMult[s] & EDavg$DelPatches == nPatchDel[p],],aes(x=Mean_LagTime,y=Mean_PercentBmassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Percentage of Biomass Lost")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentBmassChange,ymax=Highest_PercentBmassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LagTime,xmax=Highest_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Biomass plot (tryptic style)
ggplot(EDavg[EDavg$Species == nSpeciesMult[s],],aes(x=Mean_LagTime,y=Mean_PercentBmassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Percentage of Biomass Lost")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentBmassChange,ymax=Highest_PercentBmassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LagTime,xmax=Highest_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+
  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#CV plot
ggplot(EDavg[EDavg$Species == nSpeciesMult[s] & EDavg$DelPatches == nPatchDel[p],],aes(x=Mean_CVLagTime,y=-Mean_PercentCVChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("Percentage of CV Gained")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=-Lowest_PercentCVChange,ymax=-Highest_PercentCVChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#CV plot tryptic-style
ggplot(EDavg[EDavg$Species == nSpeciesMult[s],],aes(x=Mean_CVLagTime,y=-Mean_PercentCVChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("Percentage of CV Gained")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=-Lowest_PercentCVChange,ymax=-Highest_PercentCVChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


#not looking at percent changes...(but otherwise, same as the non-tryptic graphs above)
EDavg2 <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches),
                    Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F),
                    Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T),
                    Lowest_LastDebtTime=quantile(LastDebtTime, probs = 0.025, na.rm=T, names = F), Highest_LastDebtTime=quantile(LastDebtTime, probs = 0.975, na.rm=T, names = F),
                    
                    Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(CVLagTime, probs = 0.975, na.rm=T, names = F), 
                    Mean_LagTime = mean(LagTime, na.rm=T), SD_LagTime = sd(LagTime, na.rm=T), Lowest_LagTime=quantile(LagTime, probs = 0.025, na.rm=T, names = F), Highest_LagTime=quantile(LagTime, probs = 0.975, na.rm=T, names = F),
                    
                    Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), 
                    Mean_CVLagTime = mean(CVLagTime, na.rm=T), SD_CVLagTime = sd(CVLagTime, na.rm=T), Lowest_CVLagTime=quantile(CVLagTime, probs = 0.025, na.rm=T, names = F), Highest_CVLagTime=quantile(CVLagTime, probs = 0.975, na.rm=T, names = F))

#SR plot
ggplot(EDavg2[EDavg2$Species == nSpeciesMult[s] & EDavg2$DelPatches == nPatchDel[p],],aes(x=Mean_LastDebtTime,y=Mean_SRLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Species Lost")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_SRLoss,ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Biomass plot
ggplot(EDavg2[EDavg2$Species == nSpeciesMult[s] & EDavg2$DelPatches == nPatchDel[p],],aes(x=Mean_LagTime,y=Mean_BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Change in Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_BiomassChange,ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LagTime,xmax=Highest_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#CV plot
ggplot(EDavg2[EDavg2$Species == nSpeciesMult[s] & EDavg2$DelPatches == nPatchDel[p],],aes(x=Mean_CVLagTime,y=Mean_CVChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("CV Gained")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_CVChange,ymax=Highest_CVChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/(ePeriod/samplelength),y=BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
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
ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/(ePeriod/samplelength),y=BiomassChange,color=interaction(factor(Dispersal), Patch_remove),group=interaction(Scale, Patch_remove, Dispersal)))+ 
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
ggplot(ED_data[ED_data$Species == nSpeciesMult[s] & ED_data$DelPatches == nPatchDel[p],],aes(x=LagTime/(ePeriod/samplelength),y=BiomassChange,color=interaction(factor(Dispersal), Patch_remove),group=interaction(Scale, Patch_remove, Dispersal)))+ 
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
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
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
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
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

#change in _____ vs # of patches deleted
#CV change shows up underneath SR loss
ggplot(ED_data[ED_data$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove)))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = CVChange, color = "CV Change", size = 4))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = BiomassChange, color = "Biomass Change", size = 4))+
  geom_point(aes(y = SRLoss, color = "SR Loss", size = 4))+
  #scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Change")+
  ggtitle(paste(nSpeciesMult[s], "Species Initially"))+
  #geom_errorbar(aes(ymin=Mean_CVBmass - SD_CVBmass,ymax=Mean_CVBmass + SD_CVBmass),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Mean_Mean_Bmass - SD_Mean_Bmass,xmax=Mean_Mean_Bmass + SD_Mean_Bmass),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#change in _____ (greater than that due to fake patch deletion) vs # of patches deleted
ggplot(ED_data[ED_data$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove)))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = CVChange-CVChangeNoDel, color = "CV Change", size = 4))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = BiomassChange-BmassLossNoDel, color = "Biomass Change", size = 4))+
  geom_point(aes(y = SRLoss-SRLossNoDel, color = "SR Loss", size = 4))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Change")+
  ggtitle(paste(nSpeciesMult[s], "Species Initially"))+
  #geom_errorbar(aes(ymin=Mean_CVBmass - SD_CVBmass,ymax=Mean_CVBmass + SD_CVBmass),width=0.1, linetype = 2)+
  #geom_errorbarh(aes(xmin=Mean_Mean_Bmass - SD_Mean_Bmass,xmax=Mean_Mean_Bmass + SD_Mean_Bmass),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


EDdata_avgchange <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F),
                              Mean_SRLossNoDel = mean(SRLossNoDel, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.975, na.rm=T, names = F),
                              Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), 
                              Mean_BmassLossNoDel = mean(BmassLossNoDel, na.rm=T), SD_BmassLossNoDel = sd(BmassLossNoDel, na.rm = T), Lowest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.975, na.rm=T, names = F), 
                              Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_CVChangeNoDel = mean(CVChangeNoDel, na.rm=T), SD_CVChangeNoDel = sd(CVChangeNoDel, na.rm = T), Lowest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.025, na.rm=T, names = F), Highest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.975, na.rm=T, names = F))

#change in _____ vs # of patches deleted
#CV change shows up underneath SR loss, need to get the error bars to work urgh
ggplot(EDdata_avgchange[EDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_CVChange, color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_BiomassChange, color = "Biomass Change"))+
  geom_point(aes(y = Mean_SRLoss, color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Total Change")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#change in _____ (greater than that due to fake patch deletion) vs # of patches deleted
ggplot(EDdata_avgchange[EDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_CVChange - Mean_CVChangeNoDel, color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_BiomassChange - Mean_BmassLossNoDel, color = "Biomass Change"))+
  geom_point(aes(y = Mean_SRLoss - Mean_SRLossNoDel, color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Change due to Indirect Effects")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#abs(change in _____) (greater than that due to fake patch deletion) vs # of patches deleted
ggplot(EDdata_avgchange[EDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = abs(Mean_CVChange - Mean_CVChangeNoDel), color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = abs(Mean_BiomassChange - Mean_BmassLossNoDel), color = "Biomass Change"))+
  geom_point(aes(y = abs(Mean_SRLoss - Mean_SRLossNoDel), color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("abs(Change)")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#change in ______ (due to fake patch deletion only)
ggplot(EDdata_avgchange[EDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_CVChangeNoDel, color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_BmassLossNoDel, color = "Biomass Change"))+
  geom_point(aes(y = Mean_SRLossNoDel, color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Change Due to Patch Del")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#summary of all measures in the ED_data dataframe other than the Mean_Bmass and CV_Bmass measures
EDdata_avgall <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F),
                           Mean_SRLossNoDel = mean(SRLossNoDel, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.975, na.rm=T, names = F),
                           Mean_LastDebtTime = mean(LastDebtTime, na.rm=T), SD_LastDebtTime = sd(LastDebtTime, na.rm = T), Lowest_LastDebtTime=quantile(LastDebtTime, probs = 0.025, na.rm=T, names = F), Highest_LastDebtTime=quantile(LastDebtTime, probs = 0.975, na.rm=T, names = F),
                           Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), Lowest_PercentLoss=quantile(PercentLoss, probs = 0.025, na.rm=T, names = F), Highest_PercentLoss=quantile(PercentLoss, probs = 0.975, na.rm=T, names = F),
                           
                           Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), 
                           Mean_BmassLossNoDel = mean(BmassLossNoDel, na.rm=T), SD_BmassLossNoDel = sd(BmassLossNoDel, na.rm = T), Lowest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.975, na.rm=T, names = F), 
                           Mean_LagTime = mean(LagTime, na.rm=T), SD_LagTime = sd(LagTime, na.rm = T), Lowest_LagTime=quantile(LagTime, probs = 0.025, na.rm=T, names = F), Highest_LagTime=quantile(LagTime, probs = 0.975, na.rm=T, names = F), 
                           Mean_PercentBmassChange = mean(PercentBmassChange, na.rm=T), SD_PercentBmassChange = sd(PercentBmassChange, na.rm = T), Lowest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.025, na.rm=T, names = F), Highest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.975, na.rm=T, names = F),
                           
                           Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_CVChangeNoDel = mean(CVChangeNoDel, na.rm=T), SD_CVChangeNoDel = sd(CVChangeNoDel, na.rm = T), Lowest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.025, na.rm=T, names = F), Highest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.975, na.rm=T, names = F), Mean_CVLagTime = mean(CVLagTime, na.rm=T), SD_CVLagTime = sd(CVLagTime, na.rm = T), Lowest_CVLagTime=quantile(CVLagTime, probs = 0.025, na.rm=T, names = F), Highest_CVLagTime=quantile(CVLagTime, probs = 0.975, na.rm=T, names = F), Mean_PercentCVChange = mean(PercentCVChange, na.rm=T), SD_PercentCVChange = sd(PercentCVChange, na.rm = T), Lowest_PercentCVChange=quantile(PercentCVChange, probs = 0.025, na.rm=T, names = F), Highest_PercentCVChange=quantile(PercentCVChange, probs = 0.975, na.rm=T, names = F))

#SR tryptic-style plot (change due to effect alone)
ggplot(EDdata_avgall[EDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_LastDebtTime,y=(Mean_SRLoss - Mean_SRLossNoDel),color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Species Lost due to Patch Deletion Effects")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=(Lowest_SRLoss - Lowest_SRLossNoDel),ymax=(Highest_SRLoss - Highest_SRLossNoDel)),width=0.1, linetype = 2)+ 
  geom_errorbarh(aes(xmin=Lowest_LastDebtTime,xmax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#SR tryptic-style plot (change due to deletion alone) <- removed because doesn't make sense to plot this way, no delay in how the change in SR is calculated

#Biomass plot (tryptic style)
ggplot(EDdata_avgall[EDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_LagTime,y=(Mean_BiomassChange - Mean_BmassLossNoDel),color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Biomass Lost due to Patch Deletion Effects")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=(Lowest_BiomassChange - Lowest_BmassLossNoDel),ymax=(Highest_BiomassChange - Highest_BmassLossNoDel)),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_LagTime,xmax=Highest_LagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+
  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#CV plot tryptic-style
ggplot(EDdata_avgall[EDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_CVLagTime,y=(Mean_CVChange - Mean_CVChangeNoDel),color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("Change in CV Due to Effects of Patch Deletion")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=(Lowest_CVChange - Lowest_CVChangeNoDel),ymax=(Highest_CVChange - Highest_CVChangeNoDel)),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#lag time vs # of patches deleted
ggplot(EDdata_avgall[EDdata_avgall$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_CVLagTime, color = "CV"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_LagTime, color = "Biomass"))+
  geom_point(aes(y = Mean_LastDebtTime, color = "SR"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Lag Time")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  geom_errorbar(aes(ymin=Lowest_CVLagTime, ymax=Highest_CVLagTime),width=0.1, linetype = 2)+
  geom_errorbar(aes(ymin=Lowest_LagTime, ymax=Highest_LagTime),width=0.1, linetype = 2)+
  geom_errorbar(aes(ymin=Lowest_LastDebtTime, ymax=Highest_LastDebtTime),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
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
  geom_vline(x=predel_collecttime/20)+
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
  geom_vline(x=predel_collecttime)+
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
  geom_vline(x=predel_collecttime)+
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
  geom_vline(x=predel_collecttime)+
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
  geom_vline(x=predel_collecttime)+
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
  geom_vline(x=predel_collecttime)+ #makes sense because now the first 20 points are just the CV of the last period before patch deletion
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
  geom_vline(x=predel_collecttime/20)+
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
  geom_vline(x=predel_collecttime/20)+
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
  geom_vline(x=predel_collecttime/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
 
Biomass_TimeSummd2 <- summarise(group_by(Biomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

#CV, Biomass, SR on one graph
ggplot(Biomass_TimeSummd2[Biomass_TimeSummd2$Species == nSpeciesMult[s] & Biomass_TimeSummd2$DelPatches == nPatchDel[p] & Biomass_TimeSummd2$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "blue", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime/10),width=0.1, fill = "green", alpha = 0.4, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#CV, Biomass, Indiv Biomass, SR on one graph
ggplot(Biomass_TimeSummd2[Biomass_TimeSummd2$Species == nSpeciesMult[s] & Biomass_TimeSummd2$DelPatches == nPatchDel[p] & Biomass_TimeSummd2$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime/10),width=0.1, fill = "green", alpha = 0.4, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#copied from above
Biomass_TimeSummd2 <- summarise(group_by(Biomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

EDdata_avgchange <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F),
                              Mean_SRLossNoDel = mean(SRLossNoDel, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.975, na.rm=T, names = F),
                              Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), 
                              Mean_BmassLossNoDel = mean(BmassLossNoDel, na.rm=T), SD_BmassLossNoDel = sd(BmassLossNoDel, na.rm = T), Lowest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.975, na.rm=T, names = F),
                              Mean_IndivBmassLossNoDel = mean(IndivBmassLossNoDel, na.rm=T), SD_IndivBmassLossNoDel = sd(IndivBmassLossNoDel, na.rm = T), Lowest_IndivBmassLossNoDel=quantile(IndivBmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_IndivBmassLossNoDel=quantile(IndivBmassLossNoDel, probs = 0.975, na.rm=T, names = F), 
                              Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_CVChangeNoDel = mean(CVChangeNoDel, na.rm=T), SD_CVChangeNoDel = sd(CVChangeNoDel, na.rm = T), Lowest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.025, na.rm=T, names = F), Highest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.975, na.rm=T, names = F))
#^ should probably change this in the for loops below to EDdata_avgall...EDdata_avgchange is sort of a redundant data frame in the context of EDdata_avgall

Biomass_TimeSummd3 <- Biomass_TimeSummd2
for(i in 1:length(dispV)){
  for(j in 1:length(removeV)){
    for(s in 1:length(nSpeciesMult)){
      for(p in 1:length(nPatchDel)){
        Biomass_TimeSummd3$Mean_BmassLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_Biomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        Biomass_TimeSummd3$Mean_IndivBmassLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_IndivBiomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        Biomass_TimeSummd3$Mean_SRLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"] <- Biomass_TimeSummd3$Mean_SR[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"] 
        
        Biomass_TimeSummd3$Mean_CVChangeNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_CVTime[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        Biomass_TimeSummd3$Mean_BmassLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_Biomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
        
        Biomass_TimeSummd3$Mean_IndivBmassLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_IndivBiomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
        
        Biomass_TimeSummd3$Mean_SRLossNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"] <- Biomass_TimeSummd3$Mean_SR[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"] 
        
        Biomass_TimeSummd3$Mean_CVChangeNoDel[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_CVTime[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
      }}}}


#CV, Biomass, Indiv Biomass, SR on one graph <- Goal: plot this with the no-del variants as lines in the same colours
ggplot(Biomass_TimeSummd3[Biomass_TimeSummd3$Species == nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches == nPatchDel[p] & Biomass_TimeSummd3$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  geom_hline(aes(yintercept=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#copied from above
Biomass_TimeSummd2 <- summarise(group_by(Biomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

EDdata_avgchange <- summarise(group_by(ED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F),
                              
                              Mean_SRLossNoDel = mean(SRLossNoDel, na.rm=T), SD_SRLossNoDel = sd(SRLossNoDel, na.rm = T), Lowest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_SRLossNoDel=quantile(SRLossNoDel, probs = 0.975, na.rm=T, names = F),
                              Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), 
                              Mean_BmassLossNoDel = mean(BmassLossNoDel, na.rm=T), SD_BmassLossNoDel = sd(BmassLossNoDel, na.rm = T), Lowest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_BmassLossNoDel=quantile(BmassLossNoDel, probs = 0.975, na.rm=T, names = F),
                              Mean_IndivBmassLossNoDel = mean(IndivBmassLossNoDel, na.rm=T), SD_IndivBmassLossNoDel = sd(IndivBmassLossNoDel, na.rm = T), Lowest_IndivBmassLossNoDel=quantile(IndivBmassLossNoDel, probs = 0.025, na.rm=T, names = F), Highest_IndivBmassLossNoDel=quantile(IndivBmassLossNoDel, probs = 0.975, na.rm=T, names = F), 
                              Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_CVChangeNoDel = mean(CVChangeNoDel, na.rm=T), SD_CVChangeNoDel = sd(CVChangeNoDel, na.rm = T), Lowest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.025, na.rm=T, names = F), Highest_CVChangeNoDel=quantile(CVChangeNoDel, probs = 0.975, na.rm=T, names = F),
                              MeanMean_SRLossNoDel = mean(Mean_SRLossNoDel, na.rm=T),MeanMean_CVChangeNoDel = mean(Mean_CVChangeNoDel, na.rm=T), MeanMean_BmassLossNoDel = mean(Mean_BmassLossNoDel, na.rm=T),MeanMean_IndivBmassLossNoDel = mean(Mean_IndivBmassLossNoDel, na.rm=T) )
#^ should probably change this in the for loops below to EDdata_avgall...EDdata_avgchange is sort of a redundant data frame in the context of EDdata_avgall

Biomass_TimeSummd4 <- Biomass_TimeSummd2
for(i in 1:length(dispV)){
  for(j in 1:length(removeV)){
    for(s in 1:length(nSpeciesMult)){
      for(p in 1:length(nPatchDel)){
        Biomass_TimeSummd4$Mean_BmassLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"] <- 
          c(rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_IndivBmassLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"] <- c(rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_SRLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"] <- c(rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_CVChangeNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"] <- 
          c(rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_BmassLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"] <- 
          c(rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_IndivBmassLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"] <- c(rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_SRLossNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"] <- c(rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)))
        
        Biomass_TimeSummd4$Mean_CVChangeNoDel[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"] <- 
          c(rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime),rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)))
      }}}}

#CV, Biomass, Indiv Biomass, SR on one graph <- Goal: plot this with the no-del variants as lines in the same colours
ggplot(Biomass_TimeSummd4[Biomass_TimeSummd4$Species == nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches == nPatchDel[p] & Biomass_TimeSummd4$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  geom_line(aes(y=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_line(aes(y=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_line(aes(y=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_line(aes(y=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#another way to plot everything, used when only had data for one dispersal level...      
ggplot(Biomass_TimeSummd4[Biomass_TimeSummd4$Species == nSpeciesMult[s] & Biomass_TimeSummd4$Dispersal == 5e-4 & Biomass_TimeSummd4$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, DelPatches)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", 5e-4, "dispersal level", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  geom_hline(aes(yintercept=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(DelPatches~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

BiomassTime_Bin2 <- Biomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Mean_SR_Final = mean(Mean_SR, na.rm=T), Upper_SR = quantile(Mean_SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(Mean_SR, probs=.025, na.rm = T, names = F), Mean_Biomass_Final = mean(Mean_Biomass, na.rm=T), Upper_Biomass = quantile(Mean_Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Mean_Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass_Final = mean(Mean_IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(Mean_IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(Mean_IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime_Final = mean(Mean_CVTime, na.rm=T), Upper_CVTime = quantile(Mean_CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(Mean_CVTime, probs=.025, na.rm = T, names = F))

for(i in 1:length(dispV)){
  for(j in 1:length(removeV)){
    for(s in 1:length(nSpeciesMult)){
      for(p in 1:length(nPatchDel)){
        BiomassTime_Bin2$Mean_BmassLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_Biomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        BiomassTime_Bin2$Mean_IndivBmassLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_IndivBiomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        BiomassTime_Bin2$Mean_SRLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Local"] <- Biomass_TimeSummd3$Mean_SR[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"] 
        
        BiomassTime_Bin2$Mean_CVChangeNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Local"] <- 
          Biomass_TimeSummd3$Mean_CVTime[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Local"][predel_collecttime] - EDdata_avgchange$Mean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"]
        
        BiomassTime_Bin2$Mean_BmassLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_Biomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
        
        BiomassTime_Bin2$Mean_IndivBmassLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_IndivBiomass[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
        
        BiomassTime_Bin2$Mean_SRLossNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Regional"] <- Biomass_TimeSummd3$Mean_SR[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"] 
        
        BiomassTime_Bin2$Mean_CVChangeNoDel[BiomassTime_Bin2$Dispersal==dispV[i] & BiomassTime_Bin2$Patch_remove==removeV[j] & BiomassTime_Bin2$Species==nSpeciesMult[s] & BiomassTime_Bin2$DelPatches==nPatchDel[p] & BiomassTime_Bin2$Scale=="Regional"] <- 
          Biomass_TimeSummd3$Mean_CVTime[Biomass_TimeSummd3$Dispersal==dispV[i] & Biomass_TimeSummd3$Patch_remove==removeV[j] & Biomass_TimeSummd3$Species==nSpeciesMult[s] & Biomass_TimeSummd3$DelPatches==nPatchDel[p] & Biomass_TimeSummd3$Scale=="Regional"][predel_collecttime] - EDdata_avgchange$Mean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"]
      }}}}

#CV, Biomass, Indiv Biomass, SR on one graph <- Goal: plot this with the no-del variants as lines in the same colours
ggplot(BiomassTime_Bin2[BiomassTime_Bin2$Species == nSpeciesMult[s] & BiomassTime_Bin2$DelPatches == nPatchDel[p] & BiomassTime_Bin2$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass_Final/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass_Final/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR_Final, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime_Final, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime/20)+
  geom_hline(aes(yintercept=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_hline(aes(yintercept=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

Biomass_TimeSummd2 <- summarise(group_by(Biomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

Biomass_TimeSummd4 <- Biomass_TimeSummd2

BiomassTime_Bin3 <- Biomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Mean_SR_Final = mean(Mean_SR, na.rm=T), Upper_SR = quantile(Mean_SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(Mean_SR, probs=.025, na.rm = T, names = F), Mean_Biomass_Final = mean(Mean_Biomass, na.rm=T), Upper_Biomass = quantile(Mean_Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Mean_Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass_Final = mean(Mean_IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(Mean_IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(Mean_IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime_Final = mean(Mean_CVTime, na.rm=T), Upper_CVTime = quantile(Mean_CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(Mean_CVTime, probs=.025, na.rm = T, names = F))

for(i in 1:length(dispV)){
  for(j in 1:length(removeV)){
    for(s in 1:length(nSpeciesMult)){
      for(p in 1:length(nPatchDel)){
        BiomassTime_Bin3$Mean_BmassLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Local"] <- 
c(rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)/20))
        
        BiomassTime_Bin3$Mean_IndivBmassLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Local"] <- 
c(rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)/20))
        
        BiomassTime_Bin3$Mean_SRLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Local"] <- c(rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)/20))
        
        BiomassTime_Bin3$Mean_CVChangeNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Local"] <-           c(rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Local"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Local"],(length(sampleV) - predel_collecttime)/20))

        
        BiomassTime_Bin3$Mean_BmassLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Regional"] <- 
c(rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_Biomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_BmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)/20))
        
        BiomassTime_Bin3$Mean_IndivBmassLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Regional"] <- 
c(rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_IndivBiomass[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_IndivBmassLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)/20))
        
        BiomassTime_Bin3$Mean_SRLossNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Regional"] <- c(rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_SR[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_SRLossNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)/20)) 
        
        BiomassTime_Bin3$Mean_CVChangeNoDel[BiomassTime_Bin3$Dispersal==dispV[i] & BiomassTime_Bin3$Patch_remove==removeV[j] & BiomassTime_Bin3$Species==nSpeciesMult[s] & BiomassTime_Bin3$DelPatches==nPatchDel[p] & BiomassTime_Bin3$Scale=="Regional"] <- 
          c(rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]),predel_collecttime/20),rep(mean(Biomass_TimeSummd4$Mean_CVTime[Biomass_TimeSummd4$Dispersal==dispV[i] & Biomass_TimeSummd4$Patch_remove==removeV[j] & Biomass_TimeSummd4$Species==nSpeciesMult[s] & Biomass_TimeSummd4$DelPatches==nPatchDel[p] & Biomass_TimeSummd4$Scale=="Regional"][1:predel_collecttime]) - EDdata_avgchange$MeanMean_CVChangeNoDel[EDdata_avgchange$Dispersal==dispV[i] & EDdata_avgchange$Patch_remove==removeV[j] & EDdata_avgchange$Species==nSpeciesMult[s] & EDdata_avgchange$DelPatches==nPatchDel[p] & EDdata_avgchange$Scale=="Regional"],(length(sampleV) - predel_collecttime)/20))
      }}}}
      
      
#CV, Biomass, Indiv Biomass, SR on one graph <- Goal: plot this with the no-del variants as lines in the same colours
ggplot(BiomassTime_Bin3[BiomassTime_Bin3$Species == nSpeciesMult[s] & BiomassTime_Bin3$DelPatches == nPatchDel[p] & BiomassTime_Bin3$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass_Final/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass_Final/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR_Final, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime_Final, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime/20)+
  geom_line(aes(y=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_line(aes(y=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_line(aes(y=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_line(aes(y=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#binned graph for only one dispersal level at a time
ggplot(BiomassTime_Bin3[BiomassTime_Bin3$Species == nSpeciesMult[s] & BiomassTime_Bin3$Dispersal == 5e-4 & BiomassTime_Bin3$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, DelPatches)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass_Final/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass_Final/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR_Final, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime_Final, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nSpeciesMult[s], "Species and Dispersal level", 5e-4, "Local scale"))+
  geom_vline(x=predel_collecttime/20)+
  geom_line(aes(y=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  geom_line(aes(y=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  geom_line(aes(y=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  geom_line(aes(y=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(DelPatches~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


PercentBiomass_Time_Summd <- summarise(group_by(PercentBiomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F), Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F), Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

ggplot(PercentBiomass_Time_Summd[PercentBiomass_Time_Summd$Species == nSpeciesMult[s] & PercentBiomass_Time_Summd$DelPatches == nPatchDel[p] & PercentBiomass_Time_Summd$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass, colour = "Biomass")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_IndivBiomass, colour = "Average Biomass per Species"))  +  geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))+
  scale_x_log10()+
  scale_y_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, fill = "green", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass,ymax=Upper_IndivBiomass),width=0.1, fill = "red", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "cyan", alpha = 0.2, color = NA)+
  xlab("Time Step")+
  ylab("Percentage Measure")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  

PercentBiomass_Time_Bin <- PercentBiomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Upper_Biomass = quantile(Mean_Biomass, probs=0.975, na.rm = T,names=F),Lower_Biomass = quantile(Mean_Biomass, probs=0.025, na.rm = T,names=F), Mean_BiomassFinal = mean(Mean_Biomass, na.rm = T), Mean_IndivBiomassFinal = mean(Mean_IndivBiomass, na.rm = T), Upper_IndivBiomass = quantile(Mean_IndivBiomass, probs=0.975, na.rm = T,names=F),Lower_IndivBiomass = quantile(Mean_IndivBiomass, probs=0.025, na.rm = T,names=F), Mean_SRFinal = mean(Mean_SR, na.rm = T), Upper_SR = quantile(Mean_SR, probs=0.975, na.rm = T,names=F),Lower_SR = quantile(Mean_SR, probs=0.025, na.rm = T,names=F), Mean_CVTimeFinal = mean(Mean_CVTime, na.rm=T), Upper_CVTime = quantile(Mean_CVTime, probs=0.975, na.rm = T,names=F),Lower_CVTime = quantile(Mean_CVTime, probs=0.025, na.rm = T,names=F))

ggplot(PercentBiomass_Time_Bin[PercentBiomass_Time_Bin$Species == nSpeciesMult[s] & PercentBiomass_Time_Bin$DelPatches == nPatchDel[p] & PercentBiomass_Time_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_BiomassFinal, colour = "Biomass")) + geom_line(aes(y = Mean_SRFinal, colour = "Species Richness")) + geom_line(aes(y = Mean_IndivBiomassFinal, colour = "Average Biomass per Species"))  +  geom_line(aes(y = Mean_CVTimeFinal, colour = "CV Biomass"))+
  scale_x_log10()+
  scale_y_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, fill = "green", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass,ymax=Upper_IndivBiomass),width=0.1, fill = "red", alpha = 0.1, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "cyan", alpha = 0.2, color = NA)+
  xlab("Time Step")+
  ylab("Percentage Measure")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale", "Binned by 20"))+
  geom_vline(x=predel_collecttime/20)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  

EffectiveDiv_TimeSummd <- summarise(group_by(EffectiveDiv_Time,Dispersal,Patch_remove,Metric, Species, DelPatches, TimeStep), Mean_ExpShannon = mean(ExpShannon, na.rm=T), Upper_ExpShannon = quantile(ExpShannon, probs=.975, na.rm = T, names = F),Lower_ExpShannon = quantile(ExpShannon, probs=.025, na.rm = T, names = F), Mean_ExpShannonMult = mean(ExpShannonMult, na.rm=T), Upper_ExpShannonMult = quantile(ExpShannonMult, probs=.975, na.rm = T, names = F),Lower_ExpShannonMult = quantile(ExpShannonMult, probs=.025, na.rm = T, names = F))

ggplot(EffectiveDiv_TimeSummd[EffectiveDiv_TimeSummd$Species == nSpeciesMult[s] & EffectiveDiv_TimeSummd$DelPatches == nPatchDel[p],],aes(y = Mean_ExpShannon, x=TimeStep,group=interaction(Patch_remove, Dispersal, Metric), color = Metric, fill = Metric))+
  #geom_point()+ 
  geom_line() +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_ExpShannon,ymax=Upper_ExpShannon),width=0.1, color = NA, alpha = 0.1)+
  xlab("Time Step")+
  ylab("Effective Diversity Measure")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Option 1: 0's show up as 0's
SpeciesBiomass_TimeSummd <- summarise(group_by(SpeciesBiomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_Biomass = mean(Biomass, na.rm=T), SD_Biomass = sd(Biomass, na.rm = T),Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F))

#Option 2: 0's show up as NA's
SpeciesBiomass_Time2 <- SpeciesBiomass_Time
SpeciesBiomass_Time2$Biomass[which(SpeciesBiomass_Time2$Biomass==0)] <- NA
SpeciesBiomass_TimeSummd2 <- summarise(group_by(SpeciesBiomass_Time2,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_Biomass = mean(Biomass, na.rm=T), SD_Biomass = sd(Biomass, na.rm = T),Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F))

ggplot(SpeciesBiomass_TimeSummd2[SpeciesBiomass_TimeSummd2$DelPatches == nPatchDel[p] & SpeciesBiomass_TimeSummd2$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal, Species), fill = factor(Species), color = factor(Species)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass))+ ####not done fixing
  scale_x_log10()+
  scale_y_log10()+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, alpha = 0.1, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  

#one species at a time
ggplot(SpeciesBiomass_TimeSummd2[SpeciesBiomass_TimeSummd2$DelPatches == nPatchDel[p] & SpeciesBiomass_TimeSummd2$Scale == "Local" & SpeciesBiomass_TimeSummd2$Species == 11,],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass))+
  #scale_x_log10()+
  #scale_y_log10()+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, alpha = 0.1, color = NA)+
  xlab("Time Step")+
  ylab("Biomass")+
  ggtitle(paste(nPatchDel[p], "patches deleted", "Local scale", "Species 11"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines  

#####Indirect Effects
IndirectBiomass_TimeSummd <- summarise(group_by(IndirectBiomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

ggplot(IndirectBiomass_TimeSummd[IndirectBiomass_TimeSummd$Species == nSpeciesMult[s] & IndirectBiomass_TimeSummd$DelPatches == nPatchDel[p] & IndirectBiomass_TimeSummd$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Indirect Effects")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime)+
  #geom_line(aes(y=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  #geom_line(aes(y=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  #geom_line(aes(y=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  #geom_line(aes(y=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

IndirectPropBiomass_TimeSummd <- summarise(group_by(IndirectPropBiomass_Time,Dispersal,Patch_remove,Scale, Species, DelPatches, TimeStep), Mean_SR = mean(SR, na.rm=T), Upper_SR = quantile(SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(SR, probs=.025, na.rm = T, names = F), Mean_Biomass = mean(Biomass, na.rm=T), Upper_Biomass = quantile(Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass = mean(IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime = mean(CVTime, na.rm=T), Upper_CVTime = quantile(CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(CVTime, probs=.025, na.rm = T, names = F))

ggplot(IndirectPropBiomass_TimeSummd[IndirectPropBiomass_TimeSummd$Species == nSpeciesMult[s] & IndirectPropBiomass_TimeSummd$DelPatches == nPatchDel[p],],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal, Scale)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass, colour = "Biomass", linetype = Scale)) + geom_line(aes(y = Mean_IndivBiomass, colour = "Indiv Biomass", linetype = Scale)) + geom_line(aes(y = Mean_SR, colour = "Species Richness", linetype = Scale)) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass", linetype = Scale))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass,ymax=Upper_IndivBiomass),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("% (Indirect Effects)")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_vline(x=predel_collecttime)+
  #facet_grid(Dispersal~Patch_remove)+
  facet_grid(Dispersal~Patch_remove,scale="free_y")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(IndirectPropBiomass_TimeSummd[IndirectPropBiomass_TimeSummd$Species == nSpeciesMult[s] & IndirectPropBiomass_TimeSummd$DelPatches == nPatchDel[p] & IndirectPropBiomass_TimeSummd$Scale == "Local",],aes(x=TimeStep,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass, colour = "Biomass")) + geom_line(aes(y = Mean_IndivBiomass, colour = "Indiv Biomass")) + geom_line(aes(y = Mean_SR, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass,ymax=Upper_Biomass),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass,ymax=Upper_IndivBiomass),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("% (Indirect Effects)")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted, Local scale"))+
  geom_vline(x=predel_collecttime)+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free_y")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

IndirectBiomassTime_Bin <- IndirectBiomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Mean_SR_Final = mean(Mean_SR, na.rm=T), Upper_SR = quantile(Mean_SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(Mean_SR, probs=.025, na.rm = T, names = F), Mean_Biomass_Final = mean(Mean_Biomass, na.rm=T), Upper_Biomass = quantile(Mean_Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Mean_Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass_Final = mean(Mean_IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(Mean_IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(Mean_IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime_Final = mean(Mean_CVTime, na.rm=T), Upper_CVTime = quantile(Mean_CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(Mean_CVTime, probs=.025, na.rm = T, names = F))

#CV, Biomass, Indiv Biomass, SR on one graph 
ggplot(IndirectBiomassTime_Bin[IndirectBiomassTime_Bin$Species == nSpeciesMult[s] & IndirectBiomassTime_Bin$DelPatches == nPatchDel[p] & IndirectBiomassTime_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass_Final/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass_Final/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR_Final, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime_Final, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Indirect Effect")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime/20)+
  #geom_hline(aes(yintercept=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  facet_grid(Dispersal~Patch_remove)+
  #facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

IndirectPropBiomassTime_Bin <- IndirectPropBiomass_Time %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/20)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_IndivBiomass = mean(IndivBiomass, na.rm = T), Mean_SR = mean(SR, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) %>%
  group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  summarize(Mean_SR_Final = mean(Mean_SR, na.rm=T), Upper_SR = quantile(Mean_SR, probs=.975, na.rm = T, names = F),Lower_SR = quantile(Mean_SR, probs=.025, na.rm = T, names = F), Mean_Biomass_Final = mean(Mean_Biomass, na.rm=T), Upper_Biomass = quantile(Mean_Biomass, probs=.975, na.rm = T, names = F),Lower_Biomass = quantile(Mean_Biomass, probs=.025, na.rm = T, names = F), Mean_IndivBiomass_Final = mean(Mean_IndivBiomass, na.rm=T), Upper_IndivBiomass = quantile(Mean_IndivBiomass, probs=.975, na.rm = T, names = F),Lower_IndivBiomass = quantile(Mean_IndivBiomass, probs=.025, na.rm = T, names = F), Mean_CVTime_Final = mean(Mean_CVTime, na.rm=T), Upper_CVTime = quantile(Mean_CVTime, probs=.975, na.rm = T, names = F),Lower_CVTime = quantile(Mean_CVTime, probs=.025, na.rm = T, names = F))

#CV, Biomass, Indiv Biomass, SR on one graph 
ggplot(IndirectPropBiomassTime_Bin[IndirectPropBiomassTime_Bin$Species == nSpeciesMult[s] & IndirectPropBiomassTime_Bin$DelPatches == nPatchDel[p] & IndirectPropBiomassTime_Bin$Scale == "Local",],aes(x=TimeStepRound,group=interaction(Patch_remove, Dispersal)))+
  #geom_point()+ 
  geom_line(aes(y = Mean_Biomass_Final/10, colour = "Biomass/10")) + geom_line(aes(y = Mean_IndivBiomass_Final/10, colour = "Indiv Biomass/10")) + geom_line(aes(y = Mean_SR_Final, colour = "Species Richness")) + geom_line(aes(y = Mean_CVTime_Final, colour = "CV Biomass"))  +
  #divide biomass and indivbiomass by 100 if 'regional'
  scale_x_log10()+
  geom_ribbon(aes(ymin=Lower_SR,ymax=Upper_SR),width=0.1, fill = "purple", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_Biomass/10,ymax=Upper_Biomass/10),width=0.1, fill = "red", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_IndivBiomass/10,ymax=Upper_IndivBiomass/10),width=0.1, fill = "cyan", alpha = 0.4, color = NA)+
  geom_ribbon(aes(ymin=Lower_CVTime,ymax=Upper_CVTime),width=0.1, fill = "green", alpha = 0.4, color = NA)+ #read Upper_CVTime/10 until 6.23.2016
  xlab("Time Step")+
  ylab("Indirect Effect (%)")+
  ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted", "Local scale"))+
  geom_vline(x=predel_collecttime/20)+
  #geom_hline(aes(yintercept=Mean_BmassLossNoDel/10), color = "red", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_IndivBmassLossNoDel/10), color = "cyan", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_CVChangeNoDel), color = "green", linetype = "dashed")+
  #geom_hline(aes(yintercept=Mean_SRLossNoDel), color = "purple", linetype = "dashed")+
  #facet_grid(Dispersal~Patch_remove)+
  facet_grid(Dispersal~Patch_remove,scale="free_y")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


IndirectEDdata_avgchange <- summarise(group_by(IndirectED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F), Mean_SRLoss2 = mean(SRLoss2, na.rm=T), Lowest_SRLoss2=quantile(SRLoss2, probs = 0.025, na.rm=T, names = F), Highest_SRLoss2=quantile(SRLoss2, probs = 0.975, na.rm=T, names = F), Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), Lowest_PercentLoss=quantile(PercentLoss, probs = 0.025, na.rm=T, names = F), Highest_PercentLoss=quantile(PercentLoss, probs = 0.975, na.rm=T, names = F), Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), Mean_PercentBmassChange = mean(PercentBmassChange, na.rm=T), SD_PercentBmassChange = sd(PercentBmassChange, na.rm = T), Lowest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.025, na.rm=T, names = F), Highest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.975, na.rm=T, names = F), Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_PercentCVChange = mean(PercentCVChange, na.rm=T), SD_PercentCVChange = sd(PercentCVChange, na.rm = T), Lowest_PercentCVChange=quantile(PercentCVChange, probs = 0.025, na.rm=T, names = F), Highest_PercentCVChange=quantile(PercentCVChange, probs = 0.975, na.rm=T, names = F))

#change in _____ vs # of patches deleted (only indirect effect)
#need to get the error bars to work...
ggplot(IndirectEDdata_avgchange[IndirectEDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_CVChange, color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_BiomassChange, color = "Biomass Change"))+
  geom_point(aes(y = Mean_SRLoss2, color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("Indirect Effect")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#% change in _____ vs # of patches deleted (only indirect effect)
ggplot(IndirectEDdata_avgchange[IndirectEDdata_avgchange$Species == nSpeciesMult[s],],aes(x=factor(DelPatches),group=interaction(Scale, Patch_remove, Dispersal),shape = factor(Patch_remove), size = 2))+
  #scale_color_brewer("Dispersal Level", palette = "Paired")+
  geom_point(aes(y = Mean_PercentCVChange, color = "CV Change"))+
  #geom_point(aes(size = 4))+
  geom_point(aes(y = Mean_PercentBmassChange, color = "Biomass Change"))+
  geom_point(aes(y = Mean_PercentLoss, color = "SR Loss"))+
  scale_shape_manual(values=c(15,19, 17))+
  xlab("Number of Patches Deleted")+
  ylab("%")+
  scale_y_log10()+
  #ggtitle(paste(nSpeciesMult[s], "Species Initially"))+ <- I think I'm sticking with 11 species for the time being
  #geom_errorbar(aes(ymin=Lowest_CVChange, ymax=Highest_CVChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_BiomassChange, ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  #geom_errorbar(aes(ymin=Lowest_SRLoss, ymax=Highest_SRLoss),width=0.1, linetype = 2)+
  #facet_grid(Scale~.,scales = "free_y")+	
  facet_grid(Dispersal~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#summary of all measures in the IndirectED_data dataframe other than the Mean_Bmass and CV_Bmass measures

IndirectEDdata_avgall <- summarise(group_by(IndirectED_data,Dispersal,Patch_remove,Scale, Species, DelPatches), Mean_SRLoss = mean(SRLoss, na.rm=T), SD_SRLoss = sd(SRLoss, na.rm = T), Lowest_SRLoss=quantile(SRLoss, probs = 0.025, na.rm=T, names = F), Highest_SRLoss=quantile(SRLoss, probs = 0.975, na.rm=T, names = F), Mean_SRLoss2 = mean(SRLoss2, na.rm=T), Lowest_SRLoss2=quantile(SRLoss2, probs = 0.025, na.rm=T, names = F), Highest_SRLoss2=quantile(SRLoss2, probs = 0.975, na.rm=T, names = F), Mean_PercentLoss = mean(PercentLoss, na.rm=T), SD_PercentLoss = sd(PercentLoss, na.rm = T), Lowest_PercentLoss=quantile(PercentLoss, probs = 0.025, na.rm=T, names = F), Highest_PercentLoss=quantile(PercentLoss, probs = 0.975, na.rm=T, names = F), Mean_BiomassChange = mean(BiomassChange, na.rm=T), SD_BiomassChange = sd(BiomassChange, na.rm = T), Lowest_BiomassChange=quantile(BiomassChange, probs = 0.025, na.rm=T, names = F), Highest_BiomassChange=quantile(BiomassChange, probs = 0.975, na.rm=T, names = F), Mean_PercentBmassChange = mean(PercentBmassChange, na.rm=T), SD_PercentBmassChange = sd(PercentBmassChange, na.rm = T), Lowest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.025, na.rm=T, names = F), Highest_PercentBmassChange=quantile(PercentBmassChange, probs = 0.975, na.rm=T, names = F), Mean_CVChange = mean(CVChange, na.rm=T), SD_CVChange = sd(CVChange, na.rm = T), Lowest_CVChange=quantile(CVChange, probs = 0.025, na.rm=T, names = F), Highest_CVChange=quantile(CVChange, probs = 0.975, na.rm=T, names = F), Mean_PercentCVChange = mean(PercentCVChange, na.rm=T), SD_PercentCVChange = sd(PercentCVChange, na.rm = T), Lowest_PercentCVChange=quantile(PercentCVChange, probs = 0.025, na.rm=T, names = F), Highest_PercentCVChange=quantile(PercentCVChange, probs = 0.975, na.rm=T, names = F),
                                   Mean_SRLagTime = mean(SRLagTime, na.rm=T), SD_SRLagTime = sd(SRLagTime, na.rm = T), Lowest_SRLagTime=quantile(SRLagTime, probs = 0.025, na.rm=T, names = F), Highest_SRLagTime=quantile(SRLagTime, probs = 0.975, na.rm=T, names = F),
                                   Mean_SRLagTime2 = mean(SRLagTime2, na.rm=T), SD_SRLagTime2 = sd(SRLagTime2, na.rm = T), Lowest_SRLagTime2=quantile(SRLagTime2, probs = 0.025, na.rm=T, names = F), Highest_SRLagTime2=quantile(SRLagTime2, probs = 0.975, na.rm=T, names = F),
                                   Mean_BmassLagTime = mean(BmassLagTime, na.rm=T), SD_BmassLagTime = sd(BmassLagTime, na.rm = T), Lowest_BmassLagTime=quantile(BmassLagTime, probs = 0.025, na.rm=T, names = F), Highest_BmassLagTime=quantile(BmassLagTime, probs = 0.975, na.rm=T, names = F), 
                                   Mean_CVLagTime = mean(CVLagTime, na.rm=T), SD_CVLagTime = sd(CVLagTime, na.rm = T), Lowest_CVLagTime=quantile(CVLagTime, probs = 0.025, na.rm=T, names = F), Highest_CVLagTime=quantile(CVLagTime, probs = 0.975, na.rm=T, names = F))

#SR tryptic-style plot (change due to indirect effect alone)
ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_SRLagTime,y=Mean_SRLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Species Lost due to Patch Deletion Effects")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_SRLoss,ymax=Highest_SRLoss),width=0.1, linetype = 2)+ 
  geom_errorbarh(aes(xmin=Lowest_SRLagTime,xmax=Highest_SRLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_SRLagTime,y=Mean_PercentLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Last Extinction")+
  ylab("Species Lost due to Patch Deletion Effects (%)")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentLoss,ymax=Highest_PercentLoss),width=0.1, linetype = 2)+ 
  geom_errorbarh(aes(xmin=Lowest_SRLagTime,xmax=Highest_SRLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


#Biomass plot (tryptic style)
ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_BmassLagTime,y=Mean_BiomassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Biomass Lost due to Patch Deletion Effects")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_BiomassChange,ymax=Highest_BiomassChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_BmassLagTime,xmax=Highest_BmassLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+
  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_BmassLagTime,y=Mean_PercentBmassChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until Biomass Stabilizes")+
  ylab("Biomass Lost due to Patch Deletion Effects (%)")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentBmassChange,ymax=Highest_PercentBmassChange),width=0.1, linetype = 2)+ 
  geom_errorbarh(aes(xmin=Lowest_BmassLagTime,xmax=Highest_BmassLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines


#CV plot tryptic-style
#CVChange #CVLagTime
#"Time Until CV Stabilizes"
#"Change in CV Due to Effects of Patch Deletion"
ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_CVLagTime,y=Mean_CVChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("Change in CV Due to Effects of Patch Deletion")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_CVChange,ymax=Highest_CVChange),width=0.1, linetype = 2)+
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+
  #optional: facet_grid(Scale~DelPatches, scales = "free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(IndirectEDdata_avgall[IndirectEDdata_avgall$Species == nSpeciesMult[s],],aes(x=Mean_CVLagTime,y=Mean_PercentCVChange,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal, DelPatches)))+ 
  scale_color_brewer("Dispersal Level", palette = "BrBG")+ #or "Paired"
  geom_point(aes(shape = factor(Patch_remove)), size = 4)+
  #scale_x_log10()+
  #scale_shape_manual(values=c(25,19, 17))+
  scale_shape_manual(values=c(15,19, 17))+
  #scale_alpha_discrete(range = c(0.4,1))+
  xlab("Time Until CV Stabilizes")+
  ylab("Change in CV Due to Effects of Patch Deletion (%)")+
  #ggtitle(paste(nSpeciesMult[s], "Species and", nPatchDel[p], "patches deleted"))+
  geom_errorbar(aes(ymin=Lowest_PercentCVChange,ymax=Highest_PercentCVChange),width=0.1, linetype = 2)+ 
  geom_errorbarh(aes(xmin=Lowest_CVLagTime,xmax=Highest_CVLagTime),width=0.1, linetype = 2)+
  facet_grid(Scale~DelPatches)+	  #optional: facet_grid(Scale~DelPatches, scales = "free_x")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

