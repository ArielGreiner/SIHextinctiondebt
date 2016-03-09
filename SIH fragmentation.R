setwd("~/GitHub/SIH Extinction Debt")
source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)
require(data.table)

reps<-5
print.plots<-F # set this to true if you want to see the network as the sim runs - it makes it slower

nSpecies<-15
numCom<-30
randV<-50#seq(10,90,by=20) #randV/100 = % random links
#dispV <- 0.005
dispV<- c(0.0005,0.005,0.015)
dd<-1 #distance decay
numLinks<-numCom*2


rInput<-150 #resource input
rLoss<-10 #resource loss 
eff<-0.2 #conversion efficiency
mort<-0.2 #mortality
Ext<- 0.1 #extinction Threshold

ePeriod<-40000 #period of env sinusoidal fluctuations
eAMP<-1 #amplitude of envrionment sinusoidal fluctuations

drop_length<-ePeriod*4 

Tmax<-100000+drop_length*(numCom-0) #number of time steps in Sim, drop_length = # of iterations b/w patch deletions
Tdata<- seq(1, Tmax)
DT<- 0.08 # % size of discrete "time steps"
sampleV<-seq(102000,Tmax,by=2000)
removeV<-c("Max betweenness","Min betweenness","Random")

Meta_dyn_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)*3),Dispersal=rep(dispV,each=reps*(numCom-0)*3),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*3),Patches=NA,Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")),each=numCom-0),Proportion=NA)
SIH_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Regional_SR=NA,Local_SR=NA,Biomass=NA,Regional_CV=NA,Local_CV=NA)
Component_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Component_num=NA,Component_size=NA, Component_range=NA)
#Extinction Debt data frames
ED_data<-data.frame(Rep=rep(1:reps,each=(numCom-0)*2),Dispersal=rep(dispV,each=reps*(numCom-0)*2),
  Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*2),
  Scale=rep(c("Local","Regional"),each=numCom), Patches=NA, FirstDebtTime = NA, LastDebtTime = NA, SRLoss = NA)
ETime_Localdata<-data.frame(Rep=rep(1:reps, each = length(dispV)*length(removeV)*numCom*nSpecies),
  Dispersal=rep(dispV, each = length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T)),
  Species = rep(1:nSpecies, each = numCom*length(dispV)*length(removeV)), Patches = rep(1:numCom, each = length(dispV)*length(removeV)),TimeStep = NA)
ETime_Regionaldata<-data.frame(Rep=rep(1:reps, each = length(dispV)*length(removeV)*nSpecies),
  Dispersal=rep(dispV, each = length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T)),
  Species = rep(1:nSpecies, each = length(dispV)*length(removeV)), TimeStep = NA)

#keeping track of local biomass...
L_Bmass_sep <- array(data = NA, dim = c(length(removeV),length(dispV),numCom*(drop_length/2000),numCom))
L_Bmass <- array(data = NA, dim = c(length(removeV),length(dispV),numCom*(drop_length/2000)))

#initialize community network use rewire for lattice or small world - use random for random
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(r in 1:reps){
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    rand<-randV[1]
    numEdgesRewired<-rand/100*(numCom*2) 
    success<-FALSE
    while(!success){unweightedgraph<- if(rand==100) create_random_net(numCom, numLinks) else rewire(numCom,numLinks,numEdgesRewired)
    success<-length(V(unweightedgraph))==30}
    for(j in 1:3){
      weightedgraph<-addweights(unweightedgraph,numLinks,numCom)
      holdgraph<-weightedgraph
      if(print.plots==T){plot(holdgraph, ylim=c(-1,1),xlim=c(-1,1))}
      d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
      d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
      dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
      
      #vectors####
      eOptimum<-1-seq(0,eAMP, by=eAMP/(nSpecies-1)) #species environmental optima
      
      calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=length(R))
      
      Prod<-array(NA,dim=c(numCom,nSpecies,length(sampleV)))
      Abund<-Prod
      
      N0<-N<- matrix(10,ncol=nSpecies,nrow=numCom) # Community x Species abundance matrix
      R0<-R<-rep(10*(nSpecies/10),numCom)
      
      Meta_dyn<-data.frame(Species_sorting=rep(NA,length(sampleV)),Mass_effects=NA,Base_growth=NA,Patches=NA)
      Species_data<-array(NA,dim=c(length(sampleV),nSpecies,2),dimnames = list(sampleV,1:nSpecies,c("Abundance","Occupancy")))
      Components<-data.frame(Number_components=rep(NA, length(sampleV)),Component_size=NA,Component_envt_range=NA)
      
      for(TS in 1:Tmax){
        #print(TS)
        Immigrants<-calc.immigration(N,disp,dispersal_matrix)
        envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:numCom)*2*pi/numCom)+1)
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
        
        if(sum(TS==sampleV)==1){
          sample_id<-which(TS==sampleV)
          Components$Number_components[sample_id]<-components(weightedgraph)$no
          Components$Component_size[sample_id]<-mean(components(weightedgraph)$csize)
          members<-components(weightedgraph)$membership
          envt.ranges<-sapply(unique(members),function(x){range(envt.v[members==x])})
          Components$Component_envt_range[sample_id]<-mean(envt.ranges[2,]-envt.ranges[1,])
          
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
            
            Meta_dyn$Patches[sample_id]<-1
            
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
            
            Meta_dyn$Patches[sample_id]<-nrow(N)
            
            Species_data[sample_id,,1]<-colSums(N)
            Species_data[sample_id,,2]<-colSums(N>0)
          }
        }
        
        N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
        R <- Rt
        
        N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
        R0 <- Rt0
        
        if(max(TS==seq(100000+drop_length,Tmax-1,by=drop_length))){
          if(j==1){btw<-betweenness(weightedgraph)
          if(sum(btw==0)){
            patch.delete<-order(degree(weightedgraph),decreasing = T)[1]
          } else{patch.delete<-order(btw,decreasing = T)[1] }
          } else{
            if(j==2){patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1]} else{
              patch.delete<-sample(nrow(N),1)}}    
          weightedgraph<-delete.vertices(weightedgraph,patch.delete)
          d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
          d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
          dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
          dispersal_matrix[is.nan(dispersal_matrix)]<-0
          if(print.plots==T){
            if(length(V(weightedgraph))>1){plot(weightedgraph,layout=layout.circle(holdgraph)[as.numeric(colnames(dispersal_matrix)),],ylim=c(-1,1),xlim=c(-1,1))} else{plot(weightedgraph)}}
          N<-N[-patch.delete,]
          R<-R[-patch.delete]
          N0<-N0[-patch.delete,]
          R0<-R0[-patch.delete]
        }  
      } 

      L_Bmass[j,i,]<-colMeans(apply(Abund,3,rowSums),na.rm=T)
      L_Bmass_sep[j,i,,]<-data.frame(t(apply(Abund,3,rowSums)))
      R_Bmass<-apply(Abund,3,sum,na.rm=T)
      R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      L_SR<-colMeans(apply((Abund>0),3,rowSums, na.rm=T))
      
      cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
      
      CVdf<-cbind(L_Bmass_sep,data.frame(R_Bmass=R_Bmass,Patches=rep(30:1,each=drop_length/2000))) %>%
        group_by(Patches) %>%
        summarise_each(funs(cv))
      L_CV<-rowMeans(CVdf[,2:31],na.rm=T)
      R_CV<-CVdf$R_Bmass
      
      SIH_data_means<-data.frame(R_SR=R_SR,L_SR=L_SR,L_Bmass=L_Bmass,Patches=numCom-colMeans(apply(is.na(Abund),3,colSums))) %>%
        group_by(Patches) %>%
        summarise_each(funs(mean))
      SIH_data_means$R_CV<-R_CV
      SIH_data_means$L_CV<-L_CV
      
      Component_data_means<-data.frame(Patches=numCom-colMeans(apply(is.na(Abund),3,colSums)),Component_num=Components$Number_components,Component_size=Components$Component_size,Component_range=Components$Component_envt_range)%>%
        group_by(Patches) %>%
        summarise_each(funs(mean))

      SIH_data_reps[SIH_data_reps$Rep==r & SIH_data_reps$Dispersal==dispV[i] & 
        SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-SIH_data_means
      Component_data_reps[SIH_data_reps$Rep==r & 
         SIH_data_reps$Dispersal==dispV[i] & SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-Component_data_means
      
      mean.df<-summarise(group_by(Meta_dyn,Patches),Species_sorting=mean(Species_sorting,na.rm=T),Mass_effects=mean(Mass_effects,na.rm=T),Base_growth=mean(Base_growth,na.rm=T))
      Meta.dyn.long<-gather(mean.df,key = Dynamic,value=Proportion,-Patches)
      
      Meta_dyn_reps[Meta_dyn_reps$Rep==r & Meta_dyn_reps$Dispersal==dispV[i] & Meta_dyn_reps$Patch_remove==removeV[j],-c(1:3,5)]<-Meta.dyn.long[,-2]

##need to replace temp with a dataframe, temp right now gives the timestep at which a particular species goes extinct regionally    

#ETime_Regionaldata<-data.frame(Rep=rep(1:reps, each = length(dispV)*length(removeV)*nSpecies),
  #Dispersal=rep(dispV, each = length(removeV)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T)),
  #Species = rep(1:nSpecies, each = length(dispV)*length(removeV)), TimeStep = NA)      
      
for(o in 1:nSpecies){
  ETime_Regionaldata$TimeStep[ETime_Regionaldata$Rep==r & ETime_Regionaldata$Dispersal==dispV[i] & ETime_Regionaldata$Patch_remove==removeV[j] & ETime_Regionaldata$Species==o]<- max(which((apply(Abund,3,colSums, na.rm=T)>0)[o,]))
}

#Extinction Debt Things      
R_SR.df<-data.table(R_SR=colSums(apply(Abund,3,colSums, na.rm=T)>0),Patches=rep(30:1,each=drop_length/2000))

R_debt<-R_SR.df%>%
  group_by(Patches)%>%
  summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR>=first(R_SR)),Loss=first(R_SR)-last(R_SR))
  
R_debt$Debt_t[R_debt$Debt_t==(drop_length/2000)]<-NA

ED_data$FirstDebtTime[ED_data$Rep==r & ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale == "Regional"]<-R_debt$Debt_t

R_lastdebt<-R_SR.df%>%
  group_by(Patches)%>%
  summarise(Mean_SR=mean(R_SR),Debt_t=sum(R_SR!=last(R_SR)),Loss=first(R_SR)-last(R_SR))
  
R_lastdebt$Debt_t[R_lastdebt$Debt_t==(drop_length/2000)]<-0

ED_data$FirstDebtTime[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Regional"]<-R_debt$Debt_t

ED_data$LastDebtTime[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Regional"]<-R_lastdebt$Debt_t

ED_data$Patches[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j]]<-c(30:1)

ED_data$SRLoss[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Regional"]<-R_debt$Loss

L_SR.df<-data.table(L_SR=t(apply((Abund>0),3,rowSums, na.rm=T)),Patches=rep(30:1,each=drop_length/2000))

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

debt.f<-function(x){sum(x>=first(x))}

L_debt<-L_SR.df%>%
  group_by(Patches)%>%
  summarise_each(funs(debt.f))
  
ldebt.f<-function(x){sum(x!=last(x))}

L_lastdebt<-L_SR.df%>%
  group_by(Patches)%>%
  summarise_each(funs(ldebt.f))

loss.f<-function(x){sum(first(x)-last(x))}

L_loss<-L_SR.df%>%
  group_by(Patches)%>%
  summarise_each(funs(loss.f))

L_SR.df<-gather(L_debt,key = Patch,value=Debt_t,L_SR.V1:L_SR.V30) #wide -> long format
L_SR.df$Debt_t[L_SR.df$Debt_t==(drop_length/2000)]<-NA
L_SRlast.df<-gather(L_lastdebt,key = Patch,value=Debt_t,L_SR.V1:L_SR.V30) #wide -> long format
L_SRlast.df$Debt_t[L_SRlast.df$Debt_t==drop_length/2000]<-0
L_loss2<-gather(L_loss,key = Patch, value=Loss, L_SR.V1:L_SR.V30)
L_SR.df$LastDebt <- L_SRlast.df$Debt_t
L_SR.df$SRLoss <- L_loss2$Loss
#above: created data frame with SR Loss, Last Debt Time and First Debt Time, separated by patch
LocalSum<-summarise(group_by(L_SR.df, Patches), Mean_LocalFirstDebt=mean(Debt_t,na.rm=T), 
  Mean_LocalLastDebt=mean(LastDebt,na.rm=T), Mean_LocalSRLoss=mean(SRLoss,na.rm=T))

#copy of the framework of the big dataframe included below for clarity
#ED_data<-data.frame(Rep=rep(1:reps,each=(numCom-0)*2),Dispersal=rep(dispV,each=reps*(numCom-0)*2),
  #Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*2),
  #Scale=rep(c("Local","Regional"),each=numCom), Patches=NA, FirstDebtTime = NA, LastDebtTime = NA, SRLoss = NA)

#below: adds the local, mean data to the big data frame.
ED_data$FirstDebtTime[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Local"] <- LocalSum$Mean_LocalFirstDebt

ED_data$LastDebtTime[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Local"] <- LocalSum$Mean_LocalLastDebt

ED_data$SRLoss[ED_data$Rep==r & 
   ED_data$Dispersal==dispV[i] & ED_data$Patch_remove==removeV[j] & ED_data$Scale=="Local"] <- LocalSum$Mean_LocalSRLoss

    }}
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
}

#copy of the framework of the big dataframe included below for clarity
#ED_data<-data.frame(Rep=rep(1:reps,each=(numCom-0)*2),Dispersal=rep(dispV,each=reps*(numCom-0)*2),
  #Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*2),
  #Scale=rep(c("Local","Regional"),each=numCom), Patches=NA, FirstDebtTime = NA, LastDebtTime = NA, SRLoss = NA)

#if want to deal with the data not on a rep-by-rep basis
ED_data_summarized<-summarise(group_by(ED_data, Dispersal, Patch_remove, Patches, Scale), 
  Mean_FirstDebtTime=mean(FirstDebtTime,na.rm=T), SD_FirstDebtTime=sd(FirstDebtTime,na.rm=T), Mean_LastDebtTime=mean(LastDebtTime,na.rm=T), 
  SD_LastDebtTime=sd(LastDebtTime,na.rm=T), Mean_SRLoss=mean(SRLoss, na.rm=T), SD_SRLoss=sd(SRLoss, na.rm=T))

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

#regional faunal decay plot
plot(x = c(1:600),y = R_SR.df$R_SR,pch = 20, cex = 0.1, abline(v=seq(0,600, by = 20),col=3,lty=3), xlab = "Time Step", ylab = "Number of Species")

#the regional data file has -Inf if the species never appears in the community (only looking after the first 100k time steps), this is to replace
#that with an NA
  temp1 <- ETime_Regionaldata
  temp1[temp1 == -Inf] <- NA
  
  temp2 <- ETime_Localdata
  temp2[temp2 == -Inf] <- -1
  
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
        xlab = "Time to go Extinct", main = paste("Dispersal Level", dispV[o], removeV[w]), xlim = c(1,1200), breaks = 0:1201,ylim = c(0,6))
   abline(v=seq(0,1200, by = 40),col=3,lty=3)
                                                                                        
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


