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
#cl <- makeCluster(detectCores())
cl <- makeCluster(2) 
registerDoParallel(cl)
getDoParWorkers()

reps<- 2 #10
print.plots<-F # set this to true if you want to see the network as the sim runs - it makes it slower
set.seed(2)


nSpeciesMult <- c(7,11) #c(7,11) #c(7,11,15)
nPatchDel <- c(10,20) #c(5,10,20)
#nSpecies <- 11
numCom<-30
randV<- 50 #c(10,50,90)#seq(10,90,by=20) #randV/100 = % random links 
dispV <- 0.005
#dispV <- c(0.0005, 0.005, 0.015)
#dispV<- c(0.0005,0.005,0.015,0.05)#c(0.0005,0.005,0.015)
dd<-1 #distance decay
numLinks<-numCom*2


rInput<-150 #resource input
rLoss<-10 #resource loss 
eff<-0.2 #conversion efficiency
mort<-0.2 #mortality
Ext<- 0.1 #extinction Threshold

ePeriod<-40000 #period of env sinusoidal fluctuations
eAMP<-1 #amplitude of envrionment sinusoidal fluctuations

debtcollect_time <- 2000000 # number of time steps after patch deletion

Tmax<-250000+40000+debtcollect_time #+40,000 added to ensure that an entire sine wave is taken of the intact network
Tdata<- seq(1, Tmax)
DT<- 0.08 # % size of discrete "time steps"
sampleV<-seq(252000,Tmax,by=2000) #controls which time points are sampled from, this ensures that the first 20 samples (1 sine wave worth) taken are of the intact network
removeV<-c("Max betweenness","Min betweenness","Random")

Component_data_reps<-data.frame(Rep=rep(1:reps, each = length(nSpeciesMult)*length(nPatchDel)*length(dispV)*length(removeV)),Dispersal=rep(dispV,each=length(nSpeciesMult)*length(nPatchDel)*length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)), DelPatches=rep(nPatchDel), Component_num=NA,Component_size=NA, Component_range=NA)

#Data frame for recording the proportion of biomass accounted for by each of species sorting, mass effects and base growth at each sampled time point 
Meta_dyn_reps<- data.frame(Rep=rep(1:reps,each=3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Dispersal=rep(dispV,each=reps*3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*3*length(sampleV)*length(nSpeciesMult)*length(nPatchDel)),Species=rep(nSpeciesMult, each = length(nPatchDel)*3*length(sampleV)), DelPatches=rep(nPatchDel, each = 3*length(sampleV)),Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")), each = length(sampleV)),TimeStep = rep(1:length(sampleV)),Proportion=NA)

#Data frame recording the time at which the last extinction happens + the number of extinctions that happen, in each scenario
ED_data<-data.frame(Rep=rep(1:reps,each=2*length(dispV)*length(removeV)*length(nSpeciesMult)*length(nPatchDel)),Dispersal=rep(dispV,each=2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = 2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = 2), Scale=rep(c("Local","Regional")), LastDebtTime = NA, SRLoss = NA, Mean_Bmass = NA, CV_Bmass = NA)

SR_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), SR = NA)

Biomass_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),
Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA)

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
ED_data_noreps<-data.frame(Rep=r,Dispersal=rep(dispV,each=2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = 2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = 2), Scale=rep(c("Local","Regional")), LastDebtTime = NA, SRLoss = NA, Mean_Bmass = NA, CV_Bmass = NA)
  
  SR_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), SR = NA)

Biomass_Time_noreps <- data.frame(Rep=r,Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)), Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2*length(nSpeciesMult)*length(nPatchDel)),Species = rep(nSpeciesMult, each = length(sampleV)*2*length(nPatchDel)), DelPatches = rep(nPatchDel, each = length(sampleV)*2), Scale=rep(c("Local","Regional"), each = length(sampleV)),TimeStep = rep(1:length(sampleV)), Biomass = NA, IndivBiomass = NA)
  
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
        if(TS == 290000){
          
          #deletes nPatchDel[p] patches at time step = 290,000 according to whatever scheme you choose
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
      L_Bmass_sep<-data.frame(t(apply(Abund,3,rowSums)))
      R_Bmass<-apply(Abund,3,sum,na.rm=T)
      R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      L_SR<-colMeans(apply((Abund>0),3,rowSums, na.rm=T))
      
      Biomass_Time_noreps$Biomass[Biomass_Time_noreps$Rep==r & Biomass_Time_noreps$Dispersal==dispV[i] & Biomass_Time_noreps$Patch_remove==removeV[j] & Biomass_Time_noreps$Species == nSpeciesMult[s] & Biomass_Time_noreps$DelPatches == nPatchDel[p] & Biomass_Time_noreps$Scale=="Local"] <- rowMeans(L_Bmass_sep, na.rm = T)
      
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
      
      EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Alpha"] <- Effdiv_data[,1]
      
      EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Gamma"] <- Effdiv_data[,5]
      
      EffectiveDiv_Time_noreps$ExpShannon[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Beta"] <- Effdiv_data[,3]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Alpha"] <- Effdiv_data[,2]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Gamma"] <- Effdiv_data[,5]
      
      EffectiveDiv_Time_noreps$ExpShannonMult[EffectiveDiv_Time_noreps$Rep==r & EffectiveDiv_Time_noreps$Dispersal==dispV[i] & EffectiveDiv_Time_noreps$Patch_remove==removeV[j] & EffectiveDiv_Time_noreps$Species == nSpeciesMult[s] & EffectiveDiv_Time_noreps$DelPatches == nPatchDel[p] & EffectiveDiv_Time_noreps$Metric=="Beta"] <- Effdiv_data[,4]
      
      cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
      
      L_Bmass_sep_adel <- t(apply(Abund,3,rowSums))[-c(1:20),]
      L_Bmass_adel <- colMeans(apply(Abund,3,rowSums),na.rm=T)[-c(1:20)]
      localcv <- vector(length = numCom)
      #take the cv of each local community, from patch deletion to the end
      for(w in 1:numCom){
        localcv[w] <- cv(L_Bmass_sep_adel[,w])
      }
      #then take the average of all of those (see below)
      
      R_Bmass_adel<-apply(Abund,3,sum,na.rm=T)[-c(1:20)]
      
      
      
      ED_data_noreps$Mean_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-mean(R_Bmass_adel)
      ED_data_noreps$Mean_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"]<-mean(L_Bmass_adel)
      
      ED_data_noreps$CV_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"]<-mean(localcv, na.rm = T)
      ED_data_noreps$CV_Bmass[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species==nSpeciesMult[s] & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-cv(R_Bmass_adel)
      
      #Random coding note: "Component_data_noreps[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Patch_remove==removeV[j],-c(1:3)]" lets you add to all of the columns not mentioned that aren't the first 3 columns (which are rep, dispersal and patch removal columns)
      Component_data_noreps$Component_num[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:20),1]) 
      Component_data_noreps$Component_range[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:20),3]) 
      Component_data_noreps$Component_size[Component_data_noreps$Rep==r & Component_data_noreps$Dispersal==dispV[i] & Component_data_noreps$Species== nSpecies & Component_data_noreps$DelPatches == nPatchDel[p] & Component_data_noreps$Patch_remove==removeV[j]]<- mean(Components[-c(1:20),2]) 
      
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Species sorting"] <- Meta_dyn$Species_sorting
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Mass effects"] <- Meta_dyn$Mass_effects
      Meta_dyn_noreps$Proportion[Meta_dyn_noreps$Rep==r & Meta_dyn_noreps$Dispersal==dispV[i] & Meta_dyn_noreps$Patch_remove==removeV[j] & Meta_dyn_noreps$Species==nSpecies & Meta_dyn_noreps$DelPatches == nPatchDel[p] & Meta_dyn_noreps$Dynamic=="Base growth"] <- Meta_dyn$Base_growth
      
      #for species richness over time plots
      SR_Time_noreps$SR[SR_Time_noreps$Rep==r & SR_Time_noreps$Dispersal==dispV[i] & SR_Time_noreps$Species == nSpecies & SR_Time_noreps$DelPatches==nPatchDel[p] & SR_Time_noreps$Patch_remove==removeV[j] & SR_Time_noreps$Scale=="Regional"]<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      SR_Time_noreps$SR[SR_Time_noreps$Rep==r & SR_Time_noreps$Dispersal==dispV[i] & SR_Time_noreps$Species == nSpecies & SR_Time_noreps$DelPatches==nPatchDel[p] & SR_Time_noreps$Patch_remove==removeV[j] & SR_Time_noreps$Scale=="Local"]<-rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T)))
      #SR_overtime[j,i,] <- rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T)))
      
      ##Dataframe dealing with time to last extinction and SR lost 
      R_SR.df<-data.table(R_SR=colSums(apply(Abund,3,colSums,na.rm=T)>0))
      
      R_lastdebt<-R_SR.df%>%
        summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR!=last(R_SR)),Loss=first(R_SR)-last(R_SR))
      
      ED_data_noreps$LastDebtTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpecies & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-R_lastdebt$Debt_t
      
      ED_data_noreps$SRLoss[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpecies & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Regional"]<-R_lastdebt$Loss
      
      
      #L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums, na.rm=T)))
      #^ Modified 5.5.2016 because otherwise for the deleted patches, the local total SR loss and local last extinction time were 'wrong'
      L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums)))
      
      ldebt.f<-function(x){sum(x!=last(x))}
      
      L_lastdebt<-L_SR.df%>%
        summarise_each(funs(ldebt.f))
      
      loss.f<-function(x){sum(first(x)-last(x))}
      
      L_loss<-L_SR.df%>%
        summarise_each(funs(loss.f))
      
      L_SRlast.df<-gather(L_lastdebt,key = Patch,value=Debt_t,L_SR.V1:L_SR.V30) #wide -> long format
      L_loss2<-gather(L_loss,key = Patch, value=Loss, L_SR.V1:L_SR.V30)      
      
      #below: adds the local, mean data to the big data frame.
      
      ED_data_noreps$LastDebtTime[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpecies & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- mean(L_SRlast.df$Debt_t, na.rm = T)
      
      ED_data_noreps$SRLoss[ED_data_noreps$Rep==r & ED_data_noreps$Dispersal==dispV[i] & ED_data_noreps$Patch_remove==removeV[j] & ED_data_noreps$Species == nSpecies & ED_data_noreps$DelPatches==nPatchDel[p] & ED_data_noreps$Scale=="Local"] <- mean(L_loss2$Loss, na.rm = T)
      
      #Individual Patch Type Metrics
      ##need to fix this metric because 'btw' has length 25, doesn't include the missing patches...but it doesn't crash or anything
      IndivPatch_noreps$Betweenness[IndivPatch_noreps$Rep==r & IndivPatch_noreps$Dispersal==dispV[i] & IndivPatch_noreps$Species == nSpecies & IndivPatch_noreps$DelPatches==nPatchDel[p] & IndivPatch_noreps$Patch_remove==removeV[j]] <- btw
      
      IndivPatch_noreps$LastExtTime[IndivPatch_noreps$Rep==r & IndivPatch_noreps$Dispersal==dispV[i] & IndivPatch_noreps$Patch_remove==removeV[j]] <- L_SRlast.df$Debt_t
    }}
    }
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
  return(list(Component_data_noreps,Meta_dyn_noreps,ED_data_noreps, SR_Time_noreps, Biomass_Time_noreps, IndivPatch_noreps, EffectiveDiv_Time_noreps))
}

#run simulation function in parallel
Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr","vegan","data.table")) %dopar% SIH_frag()
for(r in 1:reps){
  Sim_data<-Sim_data_parallel[[r]]
  Component_data_reps[Component_data_reps$Rep==r,]<-Sim_data[[1]]
  Meta_dyn_reps[Meta_dyn_reps$Rep==r,]<-Sim_data[[2]]
  ED_data[ED_data$Rep==r,]<-Sim_data[[3]]
  SR_Time[SR_Time$Rep==r,]<-Sim_data[[4]]
  Biomass_Time[Biomass_Time$Rep==r,]<-Sim_data[[5]]
  IndivPatch[IndivPatch$Rep==r,]<-Sim_data[[6]]
  EffectiveDiv_Time[EffectiveDiv_Time$Rep==r,]<-Sim_data[[7]]
}  

save(Component_data_reps, Meta_dyn_reps, ED_data, SR_Time, Biomass_Time, IndivPatch, EffectiveDiv_Time, file = "FullFragmentationDataSet.RData")

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




