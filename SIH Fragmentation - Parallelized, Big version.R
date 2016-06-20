setwd("~/GitHub/SIH Extinction Debt")
source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)
require(data.table)
require(vegan)
require(doParallel)
require(foreach)

#set up parallel
cl <- makeCluster(detectCores())
#cl <- makeCluster(2) 
registerDoParallel(cl)
getDoParWorkers()

reps<- 4 #10
print.plots<-F # set this to true if you want to see the network as the sim runs - it makes it slower
set.seed(2)


nSpeciesMult <- 11 #c(7,11) #c(7,11,15)
nPatchDel <-c(5,10,20) #c(10,20)
#nSpecies <- 11
numCom<-30
randV<- 50 #c(10,50,90)#seq(10,90,by=20) #randV/100 = % random links 
#dispV <- 0.005
#dispV <- c(0.0005, 0.005, 0.015)
dispV<- c(0.0005,0.005,0.015,0.05)
dd<-1 #distance decay
numLinks<-numCom*2
sumby <- 20 #19


rInput<-150 #resource input
rLoss<-10 #resource loss 
eff<-0.2 #conversion efficiency
mort<-0.2 #mortality
Ext<- 0.1 #extinction Threshold

ePeriod<-40000 #period of env sinusoidal fluctuations
eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
initial_time <- 250000
debtcollect_time <- 2000000 # number of time steps after patch deletion
samplelength <- 2000

#Tmax<-250000+40000+debtcollect_time #+40,000 added to ensure that an entire sine wave is taken of the intact network
Tmax<-(initial_time)+(ePeriod*3)+debtcollect_time #+40,000 added to ensure that an entire sine wave is taken of the intact network (+ another 40,000 added so can take moving window CV for an entire sine wave + another so that, when binned, things look sort of normal)
predel_collecttime <- (Tmax - initial_time - debtcollect_time)/(samplelength)
Tdata<- seq(1, Tmax)
DT<- 0.08 # % size of discrete "time steps"
sampleV<-seq(initial_time + samplelength,Tmax,by=samplelength) #controls which time points are sampled from, this ensures that the first 20 samples (1 sine wave worth) taken are of the intact network <- 60 now
removeV<-c("Max betweenness","Min betweenness","Random")

Component_data_reps<-data.frame(Rep=rep(1:reps, each = length(nSpeciesMult)*length(nPatchDel)*length(dispV)*length(removeV)),Dispersal=rep(dispV,each=length(nSpeciesMult)*length(nPatchDel)*length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)), DelPatches=rep(nPatchDel), Component_num=NA,Component_size=NA, Component_range=NA)

#Data frame for recording the proportion of biomass accounted for by each of species sorting, mass effects and base growth at each sampled time point 
Meta_dyn_reps<- data.frame(Rep=rep(1:reps,each=3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Dispersal=rep(dispV,each=reps*3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)*3*length(sampleV)), DelPatches=rep(nPatchDel, each = 3*length(sampleV)),Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")), each = length(sampleV)),TimeStep = rep(1:length(sampleV)),Proportion=NA)

#Data frame recording the time at which the last extinction happens + the number of extinctions that happen, in each scenario
ED_data<-data.frame(Rep=rep(1:reps,each=2*length(dispV)*length(removeV)*length(nSpeciesMult)*length(nPatchDel)),Dispersal=rep(dispV,each=2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = 2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = 2), Scale=rep(c("Local","Regional")), LastDebtTime = NA, SRLoss = NA, Mean_Bmass = NA, CV_Bmass = NA, BiomassChange = NA, LagTime = NA, CVChange = NA, CVLagTime = NA, BmassLossNoDel = NA, SRLossNoDel = NA, CVChangeNoDel = NA, PercentLoss = NA, PercentBmassChange = NA, PercentCVChange = NA)

#SR_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
#Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), SR = NA)

Biomass_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)

PropBiomass_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)

PercentBiomass_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)


#Biomass_Time_Summd <- data.frame(Rep=rep(1:reps, each = (length(sampleV)/sumby)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),Dispersal=rep(dispV, each = length(removeV)*(length(sampleV)/sumby)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = (length(sampleV)/sumby)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = (length(sampleV)/sumby)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = (length(sampleV)/sumby)*2), Scale=rep(c("Local","Regional"), each = (length(sampleV)/sumby)),TimeStep = rep(1:(length(sampleV)/sumby)), Biomass = NA, IndivBiomass = NA, SR = NA)

IndivPatch <- data.frame(Rep=rep(1:reps, each = numCom*length(removeV)*length(dispV)*length(nSpeciesMult)*length(nPatchDel)), Dispersal=rep(dispV, each = length(removeV)*numCom*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = numCom*length(nSpeciesMult)*length(nPatchDel)), Species = rep(nSpeciesMult, each = numCom*length(nPatchDel)), DelPatches = rep(nPatchDel, each = numCom), Patch = rep(1:numCom), Betweenness = NA, LastExtTime = NA, iBiomass = NA)

EffectiveDiv_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*3*length(nSpeciesMult)*length(nPatchDel)), Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*3*length(nSpeciesMult)*length(nPatchDel)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*3*length(nSpeciesMult)*length(nPatchDel)), Species = rep(nSpeciesMult, each = 3*length(nPatchDel)*length(sampleV)), DelPatches = rep(nPatchDel, each = length(sampleV)*3),
Metric=rep(c("Alpha","Gamma","Beta"), each = length(sampleV)), TimeStep = rep(1:length(sampleV)), ExpShannon = NA, ExpShannonMult = NA)

#start of simulations
#initialize community network use rewire for lattice or small world - use random for random
pb <- txtProgressBar(min = 0, max = reps, style = 3)
SIH_frag<- function(){
  
Component_data_noreps<-data.frame(Rep=r,Dispersal=rep(dispV,each=length(nSpeciesMult)*length(nPatchDel)*length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)), DelPatches=rep(nPatchDel), Component_num=NA,Component_size=NA, Component_range=NA)

  #Data frame for recording the proportion of biomass accounted for by each of species sorting, mass effects and base growth at each sampled time point 
  Meta_dyn_noreps<- data.frame(Rep=r,Dispersal=rep(dispV,each=3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)*3*length(sampleV)), DelPatches=rep(nPatchDel, each = 3*length(sampleV)),Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")), each = length(sampleV)),TimeStep = rep(1:length(sampleV)),Proportion=NA) 
  
  #Data frame recording the time at which the last extinction happens + the number of extinctions that happen, in each scenario
ED_data_noreps<-data.frame(Rep=r,Dispersal=rep(dispV,each=2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = 2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = 2), Scale=rep(c("Local","Regional")), LastDebtTime = NA, SRLoss = NA, Mean_Bmass = NA, CV_Bmass = NA, BiomassChange = NA, LagTime = NA, CVChange = NA, CVLagTime = NA, BmassLossNoDel = NA, SRLossNoDel = NA, CVChangeNoDel = NA, PercentLoss = NA, PercentBmassChange = NA, PercentCVChange = NA)
  
#SR_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), SR = NA)

Biomass_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)

PropBiomass_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)

PercentBiomass_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA, CVTime = NA, SR = NA)
  
IndivPatch_noreps <- data.frame(Rep=r, Dispersal=rep(dispV, each = length(removeV)*numCom*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = numCom*length(nSpeciesMult)*length(nPatchDel)), Species = rep(nSpeciesMult, each = numCom*length(nPatchDel)), DelPatches = rep(nPatchDel, each = numCom), Patch = rep(1:numCom), Betweenness = NA, LastExtTime = NA, iBiomass = NA)
  
  EffectiveDiv_Time_noreps <- data.frame(Rep=r, Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*3*length(nSpeciesMult)*length(nPatchDel)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*3*length(nSpeciesMult)*length(nPatchDel)), Species = rep(nSpeciesMult, each = 3*length(nPatchDel)*length(sampleV)), DelPatches = rep(nPatchDel, each = length(sampleV)*3),
Metric=rep(c("Alpha","Gamma","Beta"), each = length(sampleV)), TimeStep = rep(1:length(sampleV)), ExpShannon = NA, ExpShannonMult = NA)

  for(s in 1:length(nSpeciesMult)){
  	for(p in 1:length(nPatchDel)){
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    nSpecies <- nSpeciesMult[s]
    rand<-randV[1]
    #create initial graph
    numEdgesRewired<-rand/100*(numCom*2) 
    success<-FALSE
    while(!success){unweightedgraph<- if(rand==100) create_random_net(numCom, numLinks) else rewire(numCom,numLinks,numEdgesRewired)
    success<-length(V(unweightedgraph))==30}
    for(j in 1:3){ 
      #counter <- 0
      weightedgraph<-addweights(unweightedgraph,numLinks,numCom)
      holdgraph<-weightedgraph
      if(print.plots==T){plot(holdgraph, ylim=c(-1,1),xlim=c(-1,1), main = paste("Original Graph", removeV[j]))} 
      #create dispersal matrix
      d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
      d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
      dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
      
      #define vectors for model
      eOptimum<-1-seq(0,eAMP, by=eAMP/(nSpecies-1)) #species environmental optima
      
      calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=length(R))
      
      Prod<-array(NA,dim=c(numCom,nSpecies,length(sampleV)))
      Abund<-Prod
      
      N0<-N<- matrix(10,ncol=nSpecies,nrow=numCom) # Community x Species abundance matrix
      R0<-R<-rep(10*(nSpecies/10),numCom)
      
      Meta_dyn<-data.frame(Species_sorting=rep(NA,length(sampleV)),Mass_effects=NA,Base_growth=NA)
      #Species_data<-array(NA,dim=c(length(sampleV),nSpecies,2),dimnames = list(sampleV,1:nSpecies,c("Abundance","Occupancy")))
      Effdiv_data<-array(NA,dim=c(length(sampleV),5),dimnames = list(sampleV,c("AddAlpha","MultAlpha","AddBeta","MultBeta","Gamma")))
      Components<-data.frame(Number_components=rep(NA, length(sampleV)),Component_size=NA,Component_envt_range=NA)
      
      
      for(TS in 1:Tmax){ #run through time steps
        #print(TS)
        Immigrants<-calc.immigration(N,disp,dispersal_matrix)
        envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:numCom)*2*pi/numCom)+1)
        #left over bit of code from when all patches were being deleted
        if(is.null(rownames(dispersal_matrix))){ 
          envt.v<-envt.v[as.numeric(names(dispersal_matrix))]
          consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
          Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
          Rt <- DT*rInput+R*(1-DT*(rLoss + sum(consume*N))) #resource step    
          
          Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
          Rt0 <- DT*rInput+R0*(1-DT*(rLoss + sum(consume*N0)))} else { #resource step  
            envt.v<-envt.v[as.numeric(rownames(dispersal_matrix))]
            consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
            Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
            Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
            
            Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
            Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0)))
          }
        
        #when the sampling happens...
        #sampling every time step that corresponds to an element in sampleV
        if(sum(TS==sampleV)==1){ 
          sample_id<-which(TS==sampleV)
          Components$Number_components[sample_id]<-components(weightedgraph)$no
          Components$Component_size[sample_id]<-mean(components(weightedgraph)$csize)
          members<-components(weightedgraph)$membership
          envt.ranges<-sapply(unique(members),function(x){range(envt.v[members==x])})
          Components$Component_envt_range[sample_id]<-mean(envt.ranges[2,]-envt.ranges[1,])
          #metacommunity process calculations
          if(is.null(rownames(N))){
            Prod[as.numeric(names(dispersal_matrix)),,sample_id] <- eff*consume*R*N
            Abund[as.numeric(names(dispersal_matrix)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*R*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
            fitness0<-(N0*(1+DT*(eff*R0*consume - mort))-N0)*(Nt0>Ext)
            home_prod<-mean(rowSums(fitness_w_disp*(fitness>0)))
            disp_prod_ME<-mean(rowSums(fitness_w_disp*(fitness<0 & fitness_w_disp>=0)))
            
            base_prod<-mean(sum(fitness0*(fitness0>0)))
            total_prod<-home_prod+disp_prod_ME
            
            home_prod_prop<-home_prod/total_prod
            SS_prod<-home_prod-base_prod
            SS_prod[SS_prod<0]<-0
            if(mean(sum(N>0))<=1){SS_prod<-0}
            SS<-(SS_prod/home_prod)*home_prod_prop
            SS[is.nan(SS)]<-0
            if(total_prod==0){SS<-NA}
            Meta_dyn$Species_sorting[sample_id]<-SS
            
            ME<-(disp_prod_ME)/total_prod
            ME[is.nan(ME)]<-0
            if(total_prod==0){ME<-NA}
            Meta_dyn$Mass_effects[sample_id]<-ME
            
            BP<-home_prod_prop*(1-(SS_prod/home_prod))
            BP[is.nan(BP)]<-0
            if(total_prod==0){BP<-NA}
            Meta_dyn$Base_growth[sample_id]<-BP
            
            #Meta_dyn$Patches[sample_id]<-1
            
            #Species_data[sample_id,,1]<-N
            #Species_data[sample_id,,2]<-N>0
          } else{
            Prod[as.numeric(rownames(N)),,sample_id] <- eff*consume*R*N
            Abund[as.numeric(rownames(N)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*R*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
            fitness0<-(N0*(1+DT*(eff*R0*consume - mort))-N0)*(Nt0>Ext)
            home_prod<-mean(rowSums(fitness_w_disp*(fitness>0)))
            disp_prod_ME<-mean(rowSums(fitness_w_disp*(fitness<0 & fitness_w_disp>=0)))
            
            base_prod<-mean(rowSums(fitness0*(fitness0>0)))
            total_prod<-home_prod+disp_prod_ME
            
            home_prod_prop<-home_prod/total_prod
            SS_prod<-home_prod-base_prod
            SS_prod[SS_prod<0]<-0
            if(mean(rowSums(N>0))<=1){SS_prod<-0}
            SS<-(SS_prod/home_prod)*home_prod_prop
            SS[is.nan(SS)]<-0
            if(total_prod==0){SS<-NA}
            Meta_dyn$Species_sorting[sample_id]<-SS
            
            ME<-(disp_prod_ME)/total_prod
            ME[is.nan(ME)]<-0
            if(total_prod==0){ME<-NA}
            Meta_dyn$Mass_effects[sample_id]<-ME
            
            BP<-home_prod_prop*(1-(SS_prod/home_prod))
            BP[is.nan(BP)]<-0
            if(total_prod==0){BP<-NA}
            Meta_dyn$Base_growth[sample_id]<-BP
            
            #Meta_dyn$Patches[sample_id]<-nrow(N)
            
            #Species_data[sample_id,,1]<-colSums(N)
            #Species_data[sample_id,,2]<-colSums(N>0)
          }
          
          #Diversity Metrics
          #counter <- counter + 1
          if(sample_id < length(sampleV)+1){
            #if(counter < 1021){
            #com_data <- Abund[,,counter] 
            com_data <- Abund[,,sample_id]
            com_data[is.na(com_data)] = 0
            renyi_avgshannon_a_mult <- prod(renyi(com_data,scales=1, hill=T))^(1/30)
            renyi_avgshannon_a <- mean(renyi(com_data,scales=1, hill=T))
            
            regional_data <- colSums(com_data)
            renyi_shannon_gamma <- renyi(regional_data,scales=1, hill=T)
            
            #Effdiv_data<-array(NA,dim=c(length(sampleV),5),dimnames = list(sampleV,c("AddAlpha","MultAlpha","AddBeta","MultBeta","Gamma"))) 
            Effdiv_data[sample_id,1] <- renyi_avgshannon_a
            Effdiv_data[sample_id,2] <- renyi_avgshannon_a_mult
            Effdiv_data[sample_id,3] <- renyi_shannon_gamma - renyi_avgshannon_a
            Effdiv_data[sample_id,4] <- renyi_shannon_gamma/renyi_avgshannon_a_mult
            Effdiv_data[sample_id,5] <- renyi_shannon_gamma
            
          }
        } 
        
        
        N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
        R <- Rt
        
        N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
        R0 <- Rt0
        
        #delete nPatchDel[p] patches according to scheme of choice (defined by 'j' value)
        if(TS == ((predel_collecttime*samplelength) + initial_time)){
          
          #deletes nPatchDel[p] patches at time step = 250,000+120,000 according to whatever scheme you choose
          if(j==1){btw<-betweenness(weightedgraph)
          if(sum(btw==0)){
            patch.delete<-order(degree(weightedgraph),decreasing = T)[1:nPatchDel[p]]
          } else{patch.delete<-order(btw,decreasing = T)[1:nPatchDel[p]] }
          } else{
            if(j==2){patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1:nPatchDel[p]]} else{
              patch.delete<-sample(nrow(N),nPatchDel[p])}}    
          weightedgraph<-delete.vertices(weightedgraph,patch.delete)
          d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
          d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
          dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
          dispersal_matrix[is.nan(dispersal_matrix)]<-0
          
          
          
          if(print.plots==T){
            if(length(V(weightedgraph))>1){plot(weightedgraph,layout=layout.circle(holdgraph)[as.numeric(colnames(dispersal_matrix)),],ylim=c(-1,1),xlim=c(-1,1), main = paste("Altered Graph", removeV[j]))} else{plot(weightedgraph)}}
          N<-N[-patch.delete,]
          R<-R[-patch.delete]
          N0<-N0[-patch.delete,]
          R0<-R0[-patch.delete]
        }  
      } 
      L_Bmass<-colMeans(apply(Abund,3,rowSums),na.rm=T)
      ##over here: add something to do with collecting information regarding BmassLossNoDel in the ED_data_noreps dataframe
      
      L_Bmass_sep<-data.frame(t(apply(Abund,3,rowSums)))
      R_Bmass<-apply(Abund,3,sum,na.rm=T)
      R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      L_SR<-colMeans(apply((Abund>0),3,rowSums),na.rm=T) #colMeans(apply((Abund>0),3,rowSums, na.rm=T)) <- changed for the same reason as the local SR measure was changed on 6.14.2016
      
      Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"] <- rowMeans(L_Bmass_sep, na.rm = T)
      ED_data_noreps$BmassLossNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- sum(L_Bmass_sep[predel_collecttime,patch.delete])
      #line below: absolute difference in average local biomass as calculated with vs without the deleted patches
      ED_data_noreps$BmassLossNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- abs(rowMeans(L_Bmass_sep, na.rm = T)[predel_collecttime] - rowMeans(L_Bmass_sep[,-patch.delete], na.rm = T)[predel_collecttime])
      Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"] <- R_Bmass #rowMeans(L_Bmass_sep, na.rm = T)
      
      #Average Individual Biomass (regional level)
      Biomass_Time_noreps$IndivBiomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"] <- colMeans(apply(Abund,3,colSums,na.rm = T))
     
      #Average Individual Biomass (local level)
      ##average individual biomass in each patch at each time point
      localibiomass <- matrix(nrow = numCom, ncol = length(sampleV))
      for(w in 1:numCom){
        localibiomass[w,] <- colMeans(Abund[w,,])	
      }
      #then average that value over all 30 patches
      Biomass_Time_noreps$IndivBiomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"] <- colMeans(localibiomass,na.rm=T)
      
      #or average that value over all time steps...
      IndivPatch_noreps$iBiomass[IndivPatch_noreps$Rep==r & IndivPatch_noreps$Dispersal==dispV[i] & IndivPatch_noreps$Patch_remove==removeV[j] & IndivPatch_noreps$Species== nSpeciesMult[s] & IndivPatch_noreps$DelPatches == nPatchDel[p]] <- rowMeans(localibiomass,na.rm=T)
      
      #Effdiv_data<-array(NA,dim=c(length(sampleV),5),dimnames = list(sampleV,c("AddAlpha","MultAlpha","AddBeta","MultBeta","Gamma")))
  
#local CV over time (repeats 1 value over the first 20 data points and then starts calculating the CV after the patch deletion)
x <- rowMeans(L_Bmass_sep, na.rm = T)
CV<-vector(length = length(sampleV)-sumby)
coeff_var<-function(i)
{
  select<-c(i:(i+sumby))
  CV[i]<-sd(x[select])/mean(x[select])
}
coeff_var_right<-function(i)
{
  select<-c(i:(i+sumby))
  CV[i+sumby]<-sd(x[select])/mean(x[select])
}

#firstCV <- sd(x[1:20])/mean(x[1:20])
Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"] <- c(rep(NA,(ePeriod/(samplelength))),tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var)) 
#<-c(tapply(1:predel_collecttime,1:(predel_collecttime - ePeriod/(samplelength)),coeff_var_right), tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
#<- c(rep(firstCV,20), tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
#note on tapply: the first element specifies a vector that will become 'i', the 2nd specifies the elements of that vector you want to use - but ultimately coeff_var is using the 'x' as defined above (as x in the function) -- done this way because tapply technically can only utilize functions that take one input

#calculating the CV change, locally, that would occur due singularly to fake patch deletion
allpatchversion <- c(rep(NA,(ePeriod/(samplelength))),tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
x <- rowMeans(L_Bmass_sep[,-patch.delete], na.rm = T)
missingpatchversion <- c(rep(NA,(ePeriod/(samplelength))),tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
ED_data_noreps$CVChangeNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- abs(allpatchversion[predel_collecttime] - missingpatchversion[predel_collecttime])


#regional CV over time (seems to not quite work, not sure how to fix)
x <- R_Bmass
CV<-vector(length = length(sampleV)-sumby)
coeff_var<-function(i)
{
  select<-c(i:(i+sumby))
  CV[i]<-sd(x[select])/mean(x[select])
}
coeff_var_right<-function(i)
{
  select<-c(i:(i+sumby))
  CV[i+sumby]<-sd(x[select])/mean(x[select])
}

#firstCV <- sd(x[1:20])/mean(x[1:20])

#will need to change the '21' if change the number of samples taken before patch deletion...
Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"] <- c(rep(NA,(ePeriod/(samplelength))), tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
#<-c(tapply(1:predel_collecttime,1:(predel_collecttime - ePeriod/(samplelength)),coeff_var_right), tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var)) <- didn't work because the first 2 arguments can't be the same length
#<- c(rep(firstCV,20), tapply(21:length(sampleV),21:length(sampleV),coeff_var))

#calculating the CV change, regionally, that would occur due singularly to fake patch deletion
allpatchversion <- c(rep(NA,(ePeriod/(samplelength))), tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
x <- apply(Abund[-patch.delete,,],3,sum,na.rm=T) 
missingpatchversion <- c(rep(NA,(ePeriod/(samplelength))), tapply(1:predel_collecttime,1:predel_collecttime,coeff_var_right)[1:(predel_collecttime-(ePeriod/(samplelength)))], tapply((predel_collecttime+1):length(sampleV),(predel_collecttime+1):length(sampleV),coeff_var))
ED_data_noreps$CVChangeNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- abs(allpatchversion[predel_collecttime] - missingpatchversion[predel_collecttime])


  #summarizing biomass, cvtime over time for BiomassChange and LagTime parameters      
  Biomass_Time_Summd <- Biomass_Time_noreps %>%
  group_by(Dispersal, Patch_remove, Scale, Species, DelPatches, Rep) %>%
  mutate(TimeStepRound = ceiling(TimeStep/sumby)) %>%
  group_by(TimeStepRound, Dispersal,Patch_remove, Scale, Species, DelPatches, Rep)%>%
  summarize(Mean_Biomass = mean(Biomass, na.rm = T), Mean_CVTime = mean(CVTime, na.rm = T)) 
  #group_by(Dispersal, Patch_remove, Species, DelPatches, Scale, TimeStepRound) %>%
  #summarize(SD_Biomass = sd(Mean_Biomass, na.rm = T), Mean_BiomassFinal = mean(Mean_Biomass, na.rm = T), Mean_IndivBiomassFinal = mean(Mean_IndivBiomass, na.rm = T), SD_IndivBiomass = sd(Mean_IndivBiomass, na.rm = T))

FinalBiomass <- Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][(predel_collecttime/(ePeriod/samplelength))+1]    
Lagt <- (predel_collecttime/(ePeriod/samplelength))+1
for(w in ((predel_collecttime/(ePeriod/samplelength))+1):(length(sampleV)/sumby)){
	CurrentBiomass <- Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][w]
	if(abs(CurrentBiomass - FinalBiomass) > 0.001){
		FinalBiomass <- CurrentBiomass
		Lagt <- w
	}
}

 ED_data_noreps$BiomassChange[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- abs(Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][(predel_collecttime/(ePeriod/samplelength))] - FinalBiomass)
 
ED_data_noreps$LagTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- (Lagt-(predel_collecttime/(ePeriod/samplelength)))*sumby

#same thing as above, but for CV
FinalCV <- Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][(predel_collecttime/(ePeriod/samplelength))+1]    
Lagt <- (predel_collecttime/(ePeriod/samplelength))+1
#'-1' added in the line below to deal with the fact that the last element of Biomass_Time_Summd$Mean_CVTime will always be NaN due to the finnicky-ness of the average over time nonsense 
for(w in ((predel_collecttime/(ePeriod/samplelength))+1):(length(sampleV)/sumby)-1){
	CurrentCV <- Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][w]
	if(abs(CurrentCV - FinalCV) > 0.000001){
		FinalCV <- CurrentCV
		Lagt <- w
	}
}


 ED_data_noreps$CVChange[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- abs(Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Regional"][(predel_collecttime/(ePeriod/samplelength))] - FinalCV)
 
ED_data_noreps$CVLagTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"] <- (Lagt-(predel_collecttime/(ePeriod/samplelength)))*sumby

 FinalBiomass <- Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][(predel_collecttime/(ePeriod/samplelength))+1]    
Lagt <- (predel_collecttime/(ePeriod/samplelength))+1
for(w in ((predel_collecttime/(ePeriod/samplelength))+1):(length(sampleV)/sumby)){
	CurrentBiomass <- Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][w]
	#'0.000001' fairly arbitrarily chosen, just seemed to match with when things stabilized
	if(abs(CurrentBiomass - FinalBiomass) > 0.001){
		FinalBiomass <- CurrentBiomass
		Lagt <- w
	}	
}
  
 ED_data_noreps$BiomassChange[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- abs(Biomass_Time_Summd$Mean_Biomass[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][(predel_collecttime/(ePeriod/samplelength))] - FinalBiomass)
 
ED_data_noreps$LagTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- (Lagt-(predel_collecttime/(ePeriod/samplelength)))*sumby

#same thing as above, but for CV
 FinalCV <- Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][(predel_collecttime/(ePeriod/samplelength))+1]    
Lagt <- (predel_collecttime/(ePeriod/samplelength))+1
for(w in ((predel_collecttime/(ePeriod/samplelength))+1):(length(sampleV)/sumby)-1){
	CurrentCV <- Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][w]
	#'0.000001' fairly arbitrarily chosen, just seemed to match with when things stabilized
	if(abs(CurrentCV - FinalCV) > 0.000001){ 
		FinalCV <- CurrentCV
		Lagt <- w
	}	
}
  
 ED_data_noreps$CVChange[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- abs(Biomass_Time_Summd$Mean_CVTime[Biomass_Time_Summd$Rep==r & Biomass_Time_Summd$Dispersal==dispV[i] & Biomass_Time_Summd$Patch_remove==removeV[j] & Biomass_Time_Summd$Species == nSpeciesMult[s] & Biomass_Time_Summd$DelPatches == nPatchDel[p] & Biomass_Time_Summd$Scale=="Local"][(predel_collecttime/(ePeriod/samplelength))] - FinalCV)
 
ED_data_noreps$CVLagTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- (Lagt-(predel_collecttime/(ePeriod/samplelength)))*sumby

###Effective Diversity Section           
EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Alpha"] <- Effdiv_data[,1]
      
      EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Gamma"] <- Effdiv_data[,5]
      
      EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Beta"] <- Effdiv_data[,3]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Alpha"] <- Effdiv_data[,2]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Gamma"] <- Effdiv_data[,5]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Beta"] <- Effdiv_data[,4]
      
      cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
      
      L_Bmass_sep_adel <- t(apply(Abund,3,rowSums))[-c(1:predel_collecttime),]
      L_Bmass_adel <- colMeans(apply(Abund,3,rowSums),na.rm=T)[-c(1:predel_collecttime)]
      localcv <- vector(length = numCom)
      #take the cv of each local community, from patch deletion to the end
      for(w in 1:numCom){
        localcv[w] <- cv(L_Bmass_sep_adel[,w])
      }
      #then take the average of all of those (see below)
      
      R_Bmass_adel<-apply(Abund,3,sum,na.rm=T)[-c(1:predel_collecttime)]
      
      
      
      ED_data_noreps$Mean_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-mean(R_Bmass_adel)
      ED_data_noreps$Mean_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"]<-mean(L_Bmass_adel)
      
      ED_data_noreps$CV_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"]<-mean(localcv, na.rm = T)
      ED_data_noreps$CV_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-cv(R_Bmass_adel)
      
      #Random coding note: "Component_data_noreps[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Patch_remove==removeV[j],-c(1:3)]" lets you add to all of the columns not mentioned that aren't the first 3 columns (which are rep, dispersal and patch removal columns)
      Component_data_noreps$Component_num[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:predel_collecttime),1]) 
      Component_data_noreps$Component_range[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:predel_collecttime),3]) 
      Component_data_noreps$Component_size[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:predel_collecttime),2]) 
      
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Species sorting"] <- Meta_dyn$Species_sorting
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Mass effects"] <- Meta_dyn$Mass_effects
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Base growth"] <- Meta_dyn$Base_growth
      
      #for species richness over time plots
      #SR_Time_noreps$SR[SR_Time_noreps$Rep==r & SR_Time_noreps$Dispersal==dispV[i] & SR_Time_noreps$Species == nSpecies & SR_Time_noreps$DelPatches==nPatchDel[p] & SR_Time_noreps$Patch_remove==removeV[j] & SR_Time_noreps$Scale=="Regional"]<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      #SR_Time_noreps$SR[SR_Time_noreps$Rep==r & SR_Time_noreps$Dispersal==dispV[i] & SR_Time_noreps$Species == nSpecies & SR_Time_noreps$DelPatches==nPatchDel[p] & SR_Time_noreps$Patch_remove==removeV[j] & SR_Time_noreps$Scale=="Local"]<-rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T)))
      Biomass_Time_noreps$SR[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Species == nSpecies & Biomass_Time_noreps$DelPatches==nPatchDel[p] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Scale=="Regional"]<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      #line below corrected because realized that the time step right after patch deletion should have the same avg local SR when calculated over all the patches as when calculated over just the non-deleted patches, was just a 0 vs NA problem
      Biomass_Time_noreps$SR[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Species == nSpecies & Biomass_Time_noreps$DelPatches==nPatchDel[p] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Scale=="Local"]<-rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T) #rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T))) changed on 6.14.2016
      #Regional Species Richness loss just due to fake loss of patches (i.e. the regional species richness of the community just before patch deletion - the regional species richness of the community just before patch deletion not counting those 10 patches to be deleted)
      ED_data_noreps$SRLossNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<- colSums(apply(Abund,3,colSums, na.rm=T)>0)[predel_collecttime] - colSums(apply(Abund[-patch.delete,,],3,colSums, na.rm=T)>0)[predel_collecttime]
      #Local Species Richness loss just due to fake loss of patches (i.e. the local species richness of the community just before patch deletion - the local species richness of the community just before patch deletion not counting those 10 patches to be deleted)
      ED_data_noreps$SRLossNoDel[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"]<- rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T)[predel_collecttime] - rowMeans(t(apply((Abund[-patch.delete,,]>0),3,rowSums)),na.rm=T)[predel_collecttime]
      
      ##Dataframe dealing with time to last extinction and SR lost 
      R_SR.df<-data.table(R_SR=colSums(apply(Abund,3,colSums,na.rm=T)>0))
      
      R_lastdebt<-R_SR.df%>%
      	summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR[-c(1:predel_collecttime-1)]!=last(R_SR)),Loss=R_SR[predel_collecttime]-last(R_SR))
        #summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR!=last(R_SR)),Loss=first(R_SR)-last(R_SR))
        #^ changed from what is directly above on 6.3.2016 just in case any changes in SR happen before patch deletion, and to adjust the debt time so that isn't always being artificially inflated by the pre-patch deletion times (because the biomass one isn't, so need to make them comparable) [also changed locally obvs]
      
      ED_data_noreps$LastDebtTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-R_lastdebt$Debt_t
      
      ED_data_noreps$SRLoss[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-R_lastdebt$Loss
      
      
      #L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums, na.rm=T)))
      #^ Modified 5.5.2016 because otherwise for the deleted patches, the local total SR loss and local last extinction time were 'wrong'
      L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums)))
      
      ldebt.f<-function(x){sum(x[-c(1:predel_collecttime-1)]!=last(x))}
      #old version: ldebt.f<-function(x){sum(x!=last(x))}
      
      L_lastdebt<-L_SR.df%>%
        summarise_each(funs(ldebt.f))
      
      loss.f<-function(x){sum(x[predel_collecttime]-last(x))}
      #old version: loss.f<-function(x){sum(first(x)-last(x))}
      
      L_loss<-L_SR.df%>%
        summarise_each(funs(loss.f))
      
      L_SRlast.df<-gather(L_lastdebt,key = Patch,value=Debt_t,L_SR.V1:L_SR.V30) #wide -> long format
      L_loss2<-gather(L_loss,key = Patch, value=Loss, L_SR.V1:L_SR.V30)      
      
      #below: adds the local, mean data to the big data frame.
      
      ED_data_noreps$LastDebtTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- mean(L_SRlast.df$Debt_t, na.rm = T)
      
      ED_data_noreps$SRLoss[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- mean(L_loss2$Loss, na.rm = T)

#Percent Change Metrics Calculation - Regional
Numpredel <- Biomass_Time_noreps$SR[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
Biomasspredel <- Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
CVpredel <- Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
          
          ED_data_noreps$PercentLoss[ED_data_noreps$Scale == "Regional" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (Numpredel - Biomass_Time_noreps$SR[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV)])/Numpredel
          #^ i feel like this line could be replaced with ED_data_noreps$SRLoss[...]/Numpredel now that the SRLoss metric only looks at what's going on after patch deletion though not sure about this
          ED_data_noreps$PercentBmassChange[ED_data_noreps$Scale == "Regional" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (Biomasspredel - Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV)])/Biomasspredel
          ED_data_noreps$PercentCVChange[ED_data_noreps$Scale == "Regional" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (CVpredel - Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Scale == "Regional" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV) - ePeriod/samplelength])/CVpredel
          
#Percent Change Metrics Calculation - Local
          Numpredel <- Biomass_Time_noreps$SR[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
          Biomasspredel <- Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
          CVpredel <- Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][predel_collecttime]
          
          ED_data_noreps$PercentLoss[ED_data_noreps$Scale == "Local" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (Numpredel - Biomass_Time_noreps$SR[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV)])/Numpredel 
          #^ i feel like this line could be replaced with ED_data_noreps$SRLoss[...]/Numpredel now that the SRLoss metric only looks at what's going on after patch deletion though not sure about this
          ED_data_noreps$PercentBmassChange[ED_data_noreps$Scale == "Local" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (Biomasspredel - Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV)])/Biomasspredel
          ED_data_noreps$PercentCVChange[ED_data_noreps$Scale == "Local" & ED_data_noreps$Dispersal == dispV[i] & ED_data_noreps$Patch_remove == removeV[j] & ED_data_noreps$Species == nSpeciesMult[s] & ED_data_noreps$DelPatches == nPatchDel[p] & ED_data_noreps$Rep == r]<- (CVpredel - Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Scale == "Local" & Biomass_Time_noreps$Dispersal == dispV[i] & Biomass_Time_noreps$Patch_remove == removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Rep == r][length(sampleV) - ePeriod/samplelength])/CVpredel
                    
      #Individual Patch Type Metrics
      ##need to fix this metric because 'btw' has length 25, doesn't include the missing patches...but it doesn't crash or anything
      IndivPatch_noreps$Betweenness[IndivPatch_noreps$Rep==r & IndivPatch_noreps$Dispersal==dispV[i] & IndivPatch_noreps$Species == nSpecies & IndivPatch_noreps$DelPatches==nPatchDel[p] & IndivPatch_noreps$Patch_remove==removeV[j]] <- btw
      
      IndivPatch_noreps$LastExtTime[IndivPatch_noreps$Rep==r & IndivPatch_noreps$Dispersal==dispV[i] & IndivPatch_noreps$Patch_remove==removeV[j]] <- L_SRlast.df$Debt_t
      
    #Proportional Versions of Biomass_Time Metrics
    PropBiomass_Time_noreps$Biomass[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Local"] <- rowMeans(L_Bmass_sep, na.rm = T)/Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime]

PropBiomass_Time_noreps$Biomass[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Regional"] <- R_Bmass/Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime] 

      PropBiomass_Time_noreps$IndivBiomass[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Regional"] <- colMeans(apply(Abund,3,colSums,na.rm = T))/Biomass_Time_noreps$IndivBiomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime]
      
            PropBiomass_Time_noreps$IndivBiomass[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Local"] <- colMeans(localibiomass,na.rm=T)/Biomass_Time_noreps$IndivBiomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime]
            
PropBiomass_Time_noreps$CVTime[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Local"] <- Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"]/Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime]


PropBiomass_Time_noreps$CVTime[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Species == nSpeciesMult[s] & PropBiomass_Time_noreps$DelPatches == nPatchDel[p] & PropBiomass_Time_noreps$Scale=="Regional"] <- Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"]/Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime]

PropBiomass_Time_noreps$SR[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Species == nSpecies & PropBiomass_Time_noreps$DelPatches==nPatchDel[p] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Scale=="Regional"]<-colSums(apply(Abund,3,colSums, na.rm=T)>0)/Biomass_Time_noreps$SR[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Species == nSpecies & Biomass_Time_noreps$DelPatches==nPatchDel[p] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime]

PropBiomass_Time_noreps$SR[PropBiomass_Time_noreps$Rep==r & PropBiomass_Time_noreps$Dispersal==dispV[i] & PropBiomass_Time_noreps$Species == nSpecies & PropBiomass_Time_noreps$DelPatches==nPatchDel[p] & PropBiomass_Time_noreps$Patch_remove==removeV[j] & PropBiomass_Time_noreps$Scale=="Local"]<-rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T)/Biomass_Time_noreps$SR[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Species == nSpecies & Biomass_Time_noreps$DelPatches==nPatchDel[p] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime]  

#Percentage Versions of Biomass_Time Metrics
    PercentBiomass_Time_noreps$Biomass[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Local"] <- abs(rowMeans(L_Bmass_sep, na.rm = T)/Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime])

PercentBiomass_Time_noreps$Biomass[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Regional"] <- abs((R_Bmass[predel_collecttime] - R_Bmass)/(R_Bmass[predel_collecttime])) 

      PercentBiomass_Time_noreps$IndivBiomass[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Regional"] <- abs((colMeans(apply(Abund,3,colSums,na.rm = T))[predel_collecttime] - colMeans(apply(Abund,3,colSums,na.rm = T)))/(colMeans(apply(Abund,3,colSums,na.rm = T))[predel_collecttime]))
      
            PercentBiomass_Time_noreps$IndivBiomass[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Local"] <- abs((colMeans(localibiomass,na.rm=T)[predel_collecttime] - colMeans(localibiomass,na.rm=T))/(colMeans(localibiomass,na.rm=T)[predel_collecttime]))
            
PercentBiomass_Time_noreps$CVTime[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Local"] <- abs((Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime] - Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"])
/Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"][predel_collecttime])

PercentBiomass_Time_noreps$CVTime[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Species == nSpeciesMult[s] & PercentBiomass_Time_noreps$DelPatches == nPatchDel[p] & PercentBiomass_Time_noreps$Scale=="Regional"] <- abs((Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime] - Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"])/
Biomass_Time_noreps$CVTime[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Regional"][predel_collecttime])

PercentBiomass_Time_noreps$SR[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Species == nSpecies & PercentBiomass_Time_noreps$DelPatches==nPatchDel[p] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Scale=="Regional"]<-abs((colSums(apply(Abund,3,colSums, na.rm=T)>0)[predel_collecttime] - colSums(apply(Abund,3,colSums, na.rm=T)>0))/colSums(apply(Abund,3,colSums, na.rm=T)>0)[predel_collecttime])

PercentBiomass_Time_noreps$SR[PercentBiomass_Time_noreps$Rep==r & PercentBiomass_Time_noreps$Dispersal==dispV[i] & PercentBiomass_Time_noreps$Species == nSpecies & PercentBiomass_Time_noreps$DelPatches==nPatchDel[p] & PercentBiomass_Time_noreps$Patch_remove==removeV[j] & PercentBiomass_Time_noreps$Scale=="Local"]<-abs((rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T)[predel_collecttime] - rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T))/rowMeans(t(apply((Abund>0),3,rowSums)),na.rm=T)[predel_collecttime])  
      
    }}
    }
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
  #return(list(Component_data_noreps,Meta_dyn_noreps,ED_data_noreps, SR_Time_noreps, Biomass_Time_noreps, IndivPatch_noreps, EffectiveDiv_Time_noreps))
  return(list(Component_data_noreps,Meta_dyn_noreps,ED_data_noreps, Biomass_Time_noreps, IndivPatch_noreps, EffectiveDiv_Time_noreps, PropBiomass_Time_noreps, PercentBiomass_Time_noreps))
}

#run simulation function in parallel
Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr","vegan","data.table")) %dopar% SIH_frag()
for(r in 1:reps){
  Sim_data<-Sim_data_parallel[[r]]
  Component_data_reps[Component_data_reps$Rep==r,]<-Sim_data[[1]]
  Meta_dyn_reps[Meta_dyn_reps$Rep==r,]<-Sim_data[[2]]
  ED_data[ED_data$Rep==r,]<-Sim_data[[3]]
  #SR_Time[SR_Time$Rep==r,]<-Sim_data[[4]]
  Biomass_Time[Biomass_Time$Rep==r,]<-Sim_data[[4]]
  IndivPatch[IndivPatch$Rep==r,]<-Sim_data[[5]]
  EffectiveDiv_Time[EffectiveDiv_Time$Rep==r,]<-Sim_data[[6]]
  PropBiomass_Time[PropBiomass_Time$Rep==r,]<-Sim_data[[7]]
  PercentBiomass_Time[PercentBiomass_Time$Rep==r,]<-Sim_data[[8]]
}  

#save(Component_data_reps, Meta_dyn_reps, ED_data, SR_Time, Biomass_Time, IndivPatch, EffectiveDiv_Time, file = "FullFragmentationDataSet.RData")
save(Component_data_reps, Meta_dyn_reps, ED_data, Biomass_Time, IndivPatch, EffectiveDiv_Time,PropBiomass_Time,PercentBiomass_Time, file = "FullFragmentationDataSet.RData")




