setwd("~/GitHub/SIH Extinction Debt")
source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)
require(data.table)

reps<-20
print.plots<-F # set this to true if you want to see the network as the sim runs - it makes it slower
set.seed(2)

nSpecies<-11
numCom<-30
randV<-50#seq(10,90,by=20) #randV/100 = % random links
#dispV <- 0.005
dispV<- c(0.0005,0.005,0.015,0.05)#c(0.0005,0.005,0.015)
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

#Patrick's data frames, left-over from when all patches were being deleted sequentially 
SIH_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Regional_SR=NA,Local_SR=NA,Biomass=NA,Regional_CV=NA,Local_CV=NA)
Component_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Component_num=NA,Component_size=NA, Component_range=NA)

#Data frame for recording the proportion of biomass accounted for by each of species sorting, mass effects and base growth at each sampled time point 
Meta_dyn_reps<- data.frame(Rep=rep(1:reps,each=3*length(sampleV)),Dispersal=rep(dispV,each=reps*3*length(sampleV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*3*length(sampleV)),Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")), each = length(sampleV)),TimeStep = rep(1:length(sampleV)),Proportion=NA) 

#Data frame recording the time at which the last extinction happens + the number of extinctions that happen, in each scenario
ED_data<-data.frame(Rep=rep(1:reps,each=2),Dispersal=rep(dispV,each=reps*2),
                    Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*2),Scale=rep(c("Local","Regional")), LastDebtTime = NA, SRLoss = NA)
                    
#old extinction debt dataframes from when all 30 patches were being deleted sequentially, left in just in case come in use later 
ETime_Localdata<-data.frame(Rep=rep(1:reps,each = length(dispV)*length(removeV)*numCom*nSpecies),
    Dispersal=rep(dispV, each = length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T)),Species = rep(1:nSpecies, each = numCom*length(dispV)*length(removeV)), 
    Patches = rep(1:numCom, each = length(dispV)*length(removeV)),TimeStep = NA)
ETime_Regionaldata<-data.frame(Rep=rep(1:reps, each = length(dispV)*length(removeV)*nSpecies),
   Dispersal=rep(dispV, each = length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T)),Species = rep(1:nSpecies, each = length(dispV)*length(removeV)), TimeStep = NA)
SR_Time <- data.frame(Rep=rep(1:reps, each = length(sampleV)*length(removeV)*length(dispV)*2),
        Dispersal=rep(dispV, each = length(removeV)*length(sampleV)*2),
        Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each = length(sampleV)*2),
        TimeStep = rep(1:length(sampleV)),Scale=rep(c("Local","Regional"), each = length(sampleV)), SR = NA)

#start of simulations
#initialize community network use rewire for lattice or small world - use random for random
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(r in 1:reps){
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    rand<-randV[1]
    #create initial graph
    numEdgesRewired<-rand/100*(numCom*2) 
    success<-FALSE
    while(!success){unweightedgraph<- if(rand==100) create_random_net(numCom, numLinks) else rewire(numCom,numLinks,numEdgesRewired)
    success<-length(V(unweightedgraph))==30}
    for(j in 1:3){ 
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
      Species_data<-array(NA,dim=c(length(sampleV),nSpecies,2),dimnames = list(sampleV,1:nSpecies,c("Abundance","Occupancy")))
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
            
            Species_data[sample_id,,1]<-N
            Species_data[sample_id,,2]<-N>0
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
            
            Species_data[sample_id,,1]<-colSums(N)
            Species_data[sample_id,,2]<-colSums(N>0)
          }
        } 
        
        N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
        R <- Rt
        
        N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
        R0 <- Rt0
        
     #delete 10 patches according to scheme of choice (defined by 'j' value)
        if(TS == 290000){
          
          #deletes 10 patches at time step = 290,000 according to whatever scheme you choose
          if(j==1){btw<-betweenness(weightedgraph)
          if(sum(btw==0)){
            patch.delete<-order(degree(weightedgraph),decreasing = T)[1:10]
          } else{patch.delete<-order(btw,decreasing = T)[1:10] }
          } else{
            if(j==2){patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1:10]} else{
              patch.delete<-sample(nrow(N),10)}}    
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
      #the arrays didn't work, possible that they were too big
      #L_Bmass[j,i,]<-colMeans(apply(Abund,3,rowSums),na.rm=T)
      #L_Bmass_sep[j,i,,]<-data.frame(t(apply(Abund,3,rowSums)))
      L_Bmass<-colMeans(apply(Abund,3,rowSums),na.rm=T)
      L_Bmass_sep<-data.frame(t(apply(Abund,3,rowSums)))
      R_Bmass<-apply(Abund,3,sum,na.rm=T)
      R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      L_SR<-colMeans(apply((Abund>0),3,rowSums, na.rm=T))
      
      cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
      
      #commented out the section below because set drop_length = 0 and so I believe that I would get an error message, fix if necessary later
      #CVdf<-cbind(L_Bmass_sep,data.frame(R_Bmass=R_Bmass,Patches=rep(30:1,each=drop_length/2000))) %>%
        #group_by(Patches) %>%
        #summarise_each(funs(cv))
     # L_CV<-rowMeans(CVdf[,2:31],na.rm=T)
      #R_CV<-CVdf$R_Bmass
      
      # SIH_data_means<-data.frame(R_SR=R_SR,L_SR=L_SR,L_Bmass=L_Bmass,Patches=numCom-colMeans(apply(is.na(Abund),3,colSums))) %>%
      #   group_by(Patches) %>%
      #   summarise_each(funs(mean))
      # SIH_data_means$R_CV <- NA #R_CV
      # SIH_data_means$L_CV<-NA #L_CV
      # 
      Component_data_means<-data.frame(Patches=numCom-colMeans(apply(is.na(Abund),3,colSums)),Component_num=Components$Number_components,Component_size=Components$Component_size,Component_range=Components$Component_envt_range)%>%
        group_by(Patches) %>%
        summarise_each(funs(mean))
      # 
      # SIH_data_reps[SIH_data_reps$Rep==r & SIH_data_reps$Dispersal==dispV[i] & 
      #                 SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-SIH_data_means
      Component_data_reps[SIH_data_reps$Rep==r &
                            SIH_data_reps$Dispersal==dispV[i] & SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-Component_data_means

      Meta_dyn_reps$Proportion[Meta_dyn_reps$Rep==r & Meta_dyn_reps$Dispersal==dispV[i] & Meta_dyn_reps$Patch_remove==removeV[j] & Meta_dyn_reps$Dynamic=="Species sorting"] <- Meta_dyn$Species_sorting
      Meta_dyn_reps$Proportion[Meta_dyn_reps$Rep==r & Meta_dyn_reps$Dispersal==dispV[i] & Meta_dyn_reps$Patch_remove==removeV[j] & Meta_dyn_reps$Dynamic=="Mass effects"] <- Meta_dyn$Mass_effects
      Meta_dyn_reps$Proportion[Meta_dyn_reps$Rep==r & Meta_dyn_reps$Dispersal==dispV[i] & Meta_dyn_reps$Patch_remove==removeV[j] & Meta_dyn_reps$Dynamic=="Base growth"] <- Meta_dyn$Base_growth

#leftover code from when all patches were being deleted...
      for(o in 1:nSpecies){
        ETime_Regionaldata$TimeStep[ETime_Regionaldata$Rep==r & ETime_Regionaldata$Dispersal==dispV[i] & ETime_Regionaldata$Patch_remove==removeV[j] & ETime_Regionaldata$Species==o]<- max(which((apply(Abund,3,colSums, na.rm=T)>0)[o,]))
      }
      
      L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums, na.rm=T)))
      
      LocalSR_timestep <- L_SR.df
      
      #keeping track of when each species goes extinct for the last time, in each patch
      #ETime.df <- data.frame(Species = rep(1:nSpecies, each = numCom), Patches = rep(1:numCom), TimeStep = NA)
      Ext_Time <- function(x){temp <- max(which(x>0))+1
      if(temp==min(which(is.na(x)))){ #if the extinction time is when the patch is deleted, ignore this extinction
        temp <- NA
      } 
      return(temp)
      }
      for(o in 1:nSpecies){
        #ETime.df$TimeStep[ETime.df$Species==o] <- apply(Abund[,o,],1,Ext_Time)
        ETime_Localdata$TimeStep[ETime_Localdata$Rep == r & ETime_Localdata$Dispersal==dispV[i] & ETime_Localdata$Patch_remove==removeV[j] & ETime_Localdata$Species == o] <- apply(Abund[,o,],1,Ext_Time)
      }
      
      
      #for species richness over time plots
      SR_Time$SR[SR_Time$Rep==r & SR_Time$Dispersal==dispV[i] & SR_Time$Patch_remove==removeV[j] & SR_Time$Scale=="Regional"]<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      SR_Time$SR[SR_Time$Rep==r & SR_Time$Dispersal==dispV[i] & SR_Time$Patch_remove==removeV[j] & SR_Time$Scale=="Local"]<-rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T)))
      #SR_overtime[j,i,] <- rowMeans(t(apply((Abund>0),3,rowSums, na.rm=T)))
      
      ##Dataframe dealing with time to last extinction and SR lost 
      R_SR.df<-data.table(R_SR=colSums(apply(Abund,3,colSums,na.rm=T)>0))
      
      R_lastdebt<-R_SR.df%>%
        summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR!=last(R_SR)),Loss=first(R_SR)-last(R_SR))
      
      ED_data$LastDebtTime[ED_data$Rep==r & ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Regional"]<-R_lastdebt$Debt_t
      
      ED_data$SRLoss[ED_data$Rep==r & ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Regional"]<-R_lastdebt$Loss
      
      
      L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums, na.rm=T)))
      
      ldebt.f<-function(x){sum(x!=last(x))}
      
      L_lastdebt<-L_SR.df%>%
        summarise_each(funs(ldebt.f))
      
      loss.f<-function(x){sum(first(x)-last(x))}
      
      L_loss<-L_SR.df%>%
        summarise_each(funs(loss.f))
      
      L_SRlast.df<-gather(L_lastdebt,key = Patch,value=Debt_t,L_SR.V1:L_SR.V30) #wide -> long format
      L_loss2<-gather(L_loss,key = Patch, value=Loss, L_SR.V1:L_SR.V30)
      
      
      #copy of the framework of the big dataframe included below for clarity
      #ED_data<-data.frame(Rep=rep(1:reps,each=(numCom-0)*2),Dispersal=rep(dispV,each=reps*(numCom-0)*2),
      #Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*2),
      #Scale=rep(c("Local","Regional"),each=numCom), Patches=NA, FirstDebtTime = NA, LastDebtTime = NA, SRLoss = NA)
      
      #below: adds the local, mean data to the big data frame.
      
      ED_data$LastDebtTime[ED_data$Rep==r & ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Local"] <- mean(L_SRlast.df$Debt_t)
      
      ED_data$SRLoss[ED_data$Rep==r & ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Local"] <- mean(L_loss2$Loss)
      
    }}
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
}

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
  #facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Supplementary Material Figure 2 (4.15.2016)
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

#supplementary materials figure 1 (4/3/2016)
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
  


###old plots
ggplot(ED_data,aes(x=LastDebtTime,y=SRLoss,color=factor(Dispersal),group=interaction(Scale, Patch_remove, Dispersal)))+
  geom_point()+ 
  geom_point(aes(shape = factor(Patch_remove)))+
  xlab("Time Until Last Extinction")+
  ylab("Number of Species Lost")+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#number of species lost vs time until last extinction plot, split into local and regional plots <- alternate, more colourful version
require(ggplot2)
#number of species lost vs time until last extinction plot
ggplot(EDdata_avg,aes(x=Mean_LastDebtTime,y=Mean_SRLoss,color=interaction(Dispersal, Patch_remove),group=interaction(Scale, Patch_remove, Dispersal)))+
  geom_point()+ 
  xlab("Time Until Last Extinction")+
  ylab("Number of Species Lost")+
  geom_errorbar(aes(ymin=Mean_SRLoss-SD_SRLoss,ymax=Mean_SRLoss+SD_SRLoss),width=0.1)+
  facet_grid(Scale~.)+	  
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

###Tester Plots
#will plot the local biomass of each individual patch for the last scenario run
plot(L_Bmass_sep$X30, type = 'l')
#7, 11, 18, 24, 30
plot(L_Bmass_sep$X30, type = 'l', xlab = "Time Step",ylab = "Biomass", main = paste("Biomass of Patch", sep = " ",30, "over time [Dispersal = ",dispV[i], ", removal sequence = ", removeV[j], "]"))


###Plots no longer in use (some of which were from the 30 patch deletion days)

#copy of the framework of the big dataframe included below for clarity
#ED_data<-data.frame(Rep=rep(1:reps,each=(numCom-0)*2),Dispersal=rep(dispV,each=reps*(numCom-0)*2),
#Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*2),
#Scale=rep(c("Local","Regional"),each=numCom), Patches=NA, FirstDebtTime = NA, LastDebtTime = NA, SRLoss = NA)

#if want to deal with the data not on a rep-by-rep basis
#ED_data_summarized<-summarise(group_by(ED_data, Dispersal, Patch_remove, Patches, Scale), 
                              #Mean_FirstDebtTime=mean(FirstDebtTime,na.rm=T), SD_FirstDebtTime=sd(FirstDebtTime,na.rm=T), Mean_LastDebtTime=mean(LastDebtTime,na.rm=T), 
                              #SD_LastDebtTime=sd(LastDebtTime,na.rm=T), Mean_SRLoss=mean(SRLoss, na.rm=T), SD_SRLoss=sd(SRLoss, na.rm=T))


#regional faunal decay plot
plot(x = c(1:600),y = R_SR.df$R_SR,pch = 20, cex = 0.1, abline(v=seq(0,600, by = 20),col=3,lty=3), xlab = "Time Step", ylab = "Number of Species")

#the regional data file has -Inf if the species never appears in the community (only looking after the first 100k time steps), this is to replace
#that with an NA
#if regional data time to go extinct = 220, or local data time to go extinct = 221 - that means that the species doesn't go extinct during 
#the elapsed time in the simulation
temp1 <- ETime_Regionaldata
temp1[temp1 == -Inf] <- NA
temp1[temp1 == length(sampleV)] <- NA 

temp2 <- ETime_Localdata
temp2[temp2 == -Inf] <- -1
temp2[temp2 == length(sampleV)+1] <- NA

#summarize over all patches
temp2.1 <- summarise(group_by(temp2, Dispersal, Patch_remove, Rep, Species), 
                     Mean_TimeStep=mean(TimeStep,na.rm=T), SD_TimeStep=sd(TimeStep,na.rm=T))

ETime.df <- data.frame(Dispersal=rep(dispV, each = length(removeV)*2*reps),
                       Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T), each=2*reps),
                       Scale=rep(c("Local","Regional"), each = reps), Rep = rep(1:reps),Species = rep(1:nSpecies, each = length(dispV)*length(removeV)*2*reps), 
                       TimeStep = NA)

#sorry Patrick
for(o in 1:length(dispV)){
  for(w in 1:length(removeV)){
    for(i in 1:nSpecies){
      for(j in 1:reps){
        ETime.df$TimeStep[ETime.df$Scale == "Regional" & ETime.df$Dispersal == dispV[o] & ETime.df$Patch_remove == removeV[w] & ETime.df$Species == i & ETime.df$Rep == j] <- temp1$TimeStep[temp1$Dispersal == dispV[o] & temp1$Patch_remove == removeV[w] & temp1$Species == i & temp1$Rep == j]
        
        ETime.df$TimeStep[ETime.df$Scale == "Local" & ETime.df$Dispersal == dispV[o] & ETime.df$Patch_remove == removeV[w] & ETime.df$Species == i & ETime.df$Rep == j] <- temp2.1$Mean_TimeStep[temp2.1$Dispersal == dispV[o] & temp2.1$Patch_remove == removeV[w] & temp2.1$Species == i & temp2.1$Rep == j]
      }
    }	
  }
}

#getting rid of the data points that I called '-1' up above in the histogram
ETime2.df <- ETime.df
ETime2.df$TimeStep[ETime2.df$Scale == "Local" & ETime2.df$TimeStep == -1] <- NA

#histogram - need to reset ylim each time, set so that there is a bar for every time step  
par(mfrow=c(length(removeV),length(dispV)))
for(w in 1:length(removeV)){
  for(o in 1:length(dispV)){
    hist(ETime2.df$TimeStep[ETime2.df$Scale == "Regional" & ETime2.df$Dispersal == dispV[o] & ETime2.df$Patch_remove == removeV[w]], 
         xlab = "Time to go Extinct", main = paste("Dispersal Level", dispV[o], removeV[w]), xlim = c(1,length(sampleV)+1), ylim = c(0,50),
    breaks = seq(1,length(sampleV)+1, by = 20))
    
  }
}
ETimeSum.df <- summarise(group_by(ETime2.df, Dispersal, Patch_remove, Species, Scale), Mean_TimeStep = mean(TimeStep, na.rm=T), SD_TimeStep = sd(TimeStep, na.rm = T))

par(mfrow=c(length(removeV),length(dispV)))
for(w in 1:length(removeV)){
  for(o in 1:length(dispV)){
    plot(ETimeSum.df$Mean_TimeStep[ETimeSum.df$Scale == "Regional" & ETimeSum.df$Dispersal == dispV[o] & ETimeSum.df$Patch_remove == removeV[w]], 
         ylab = "Time to go Extinct", main = paste("Dispersal Level", dispV[o], removeV[w]), xlim = c(1,nSpecies), ylim = c(0,length(sampleV)+5),
         type = 'l')
    
  }
}




###Odd foray into box plots that didn't really work###
#makes the boxes transparent (the last 2 digits of the hex code define the level of transparency)
stripchart(ETime2.df$TimeStep[ETime2.df$Scale == "Local" & ETime2.df$Dispersal == dispV[o] & ETime2.df$Patch_remove == removeV[w] & ETime2.df$Rep == 1],
           ylab="Time to go Extinct",main = paste("Dispersal Level", dispV[o], removeV[w]),
           col=c("light blue"), vertical = TRUE, pch = 19)
#col=c("pink","light blue","pink","light blue","pink","light blue"), vertical = TRUE, pch = 19)
boxplot(ETime2.df$TimeStep[ETime2.df$Scale == "Local" & ETime2.df$Dispersal == dispV[o] & ETime2.df$Patch_remove == removeV[w] & ETime2.df$Rep == 1],
        ylab="Time to go Extinct",main = paste("Dispersal Level", dispV[o], removeV[w]), 
        col=c("light blue"), vertical = TRUE, pch = 19)
#col=c("#FF003322","#9AC0CD22","#FF003322","#9AC0CD22","#FF003322","#9AC0CD22"), add = TRUE)

par(mfrow=c(length(removeV),length(dispV)))
for(w in 1:length(removeV)){
  for(o in 1:length(dispV)){
    stripchart(ETime2.df$TimeStep[ETime2.df$Scale == "Local" & ETime2.df$Dispersal == dispV[o] & ETime2.df$Patch_remove == removeV[w] & ETime2.df$Rep == 1],
               xlab="Time to go Extinct",main = paste("Dispersal Level", dispV[o], removeV[w]),
               col=c("light blue"), vertical = FALSE, pch = 19)
    
  }
}

#not sure if need the line below, because the one time that i tested things i saw that there were 4 extinctions...and when i checked the SR
#over time plot i saw that there was an extinction in that last time step...(0.05 dispersal, min betweenness, regional)
#ED_data$LastDebtTime == length(sampleV) - 1 <- NA

#plot below here is no longer relevant, and if did want to use it again - would need to fix it so that the scale was the same for all of them
require(ggplot2)
#number of species lost vs time until last extinction plot
ggplot(ED_data,aes(x=LastDebtTime,y=SRLoss,color=Scale,group=interaction(Scale, Patch_remove, Dispersal, Rep),fill=Scale),alpha=0.1)+
  geom_point()+ 
  #geom_line()+
  #stat_smooth(method = 'lm', formula = y ~ poly(x,2))+
  stat_smooth(method = 'lm')+
  #geom_ribbon(aes(ymin=Mean_SRLoss-SD_SRLoss,ymax=Mean_SRLoss+SD_SRLoss),width=0.1)+
  facet_grid(Dispersal~Patch_remove,scale="free")+
  #facet_grid(Scale~Patch_remove,scale="free")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines




                         


