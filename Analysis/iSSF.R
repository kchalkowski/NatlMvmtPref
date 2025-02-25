# Pipeline --------------------------------------------------------------------

#MakeFalseSteps > MakeAvailabilityGroups > iSSF > IDsPigsSSF

# Purpose --------------------------------------------------------------------

#Run iSSF on full dataset
#Format steps for iSSF
#Run full model
#Check intxn with step length on full model
#Run leave-one-out loop by study

# Script Setup --------------------------------------------------------------------

#set home directory
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_SSF/"

# load libraries
#remotes::install_github("hrbrmstr/mactheknife") #install if using aliasing for large spatial files
#library(mactheknife) #only load this if using alias as pointer to spatial files
library(raster)
library(amt)
library(sf)
library(dplyr)
library(udpipe)
library(lubridate)
library(dismo)
library(survival)
library(pROC)
library(tidyr)
library(tidyverse)
library(forcats)

#set directory paths
indir=file.path(home,"2_Data","Input") #initial input read from here
objdir=file.path(home,"2_Data","Objects") #intermediate objects, lookup tables, etc. go here
outdir=file.path(home,"4_Output") #intermediate objects, lookup tables, etc. go here

#load data
steps_g=readRDS(file.path(objdir,"steps_grouped.rds",fsep=.Platform$file.sep))

#source functions
source(file.path(home,"1_Scripts","Analysis","Functions","iSSF_Function.R"))

# Loop iSSF Model --------------------------------------------------------------------
groups=unique(steps_g$avail_group)

#Set up loop by group for each conditional logistic model using contr.sum. 
for(i in 1:length(groups)) {
  print(i)
  #Subset to group i
  Dat_T<-subset(steps_g,avail_group==groups[i])
  
  #Get summary of animal id nums in group i
  IDT<-table(Dat_T$animalnum)
  
  #drop factor levels
  Dat_T<-droplevels(Dat_T)
  
  #get summaries of step/lc counts
  NLCDT<-as.data.frame(table(Dat_T$nlcd,Dat_T$case_,dnn=c("nlcd","case")))
  NLCDT<-subset(NLCDT,NLCDT$case=="TRUE"&NLCDT$Freq>0)
  NLCD_ct<-dplyr::count(NLCDT)
  
  #Run iSSF model with contra sums process
  SSF_out=contrSum_iSSF(Dat_T,NLCDT,"animal_p_stp","nlcd_only","parms")
  SSF_params_i=as.data.frame(SSF_out$params)
  SSF_params_i$parms=rownames(SSF_params_i)
  SSF_params_i$aic=SSF_out$aic
  
  SSF_sl_out=contrSum_iSSF(Dat_T,NLCDT,"animal_p_stp","nlcd_sl","parms")
  SSF_sl_params_i=as.data.frame(SSF_sl_out$params)
  SSF_sl_params_i$parms=rownames(SSF_sl_params_i)
  SSF_sl_params_i$aic=SSF_sl_out$aic
  
  SSF_params_i$model="nlcd_only"
  SSF_sl_params_i$model="nlcd_sl"
  
  SSF_params_i$group=groups[i]
  SSF_sl_params_i$group=groups[i]
  
  SSF_i=rbind(SSF_params_i,SSF_sl_params_i)
  
  if(i==1){
    SSF_full=SSF_i
  } else{
    SSF_full=rbind(SSF_full,SSF_i)
  }
  
}


# Tables of study, ID, and group and land use type
Table5<- table(steps_g$study,steps_g$avail_group,dnn=c("study","group"))
Table2<-table(steps_g$nlcd,steps_g$avail_group,steps_g$case_,dnn=c("nlcd","avail_group","case"))
Table3<-table(steps_g$animalnum,steps_g$avail_group,dnn=c("animalnum","avail_group"))

# Leave one out sensitivity analysis ------------------------------------------

#subset to groups with >1 study
grp_singles=table(steps_g$study,steps_g$avail_group,dnn=c("study","avail_group"))
grp_singles[grp_singles>0]=1
l1o_steps=steps_g[(steps_g$avail_group%in%which(colSums(grp_singles)>1)),]

#get lookup with 1/0 if study is in group
study_grp_avail=as.data.frame(table(steps_g$study,steps_g$avail_group,dnn=c("study","avail_group")))
study_grp_avail$Freq[study_grp_avail$Freq>0]=1
study_grp_availw=study_grp_avail %>% pivot_wider(values_from=Freq,names_from=avail_group)
study_grp_availw=study_grp_availw[,c(1,(which(colSums(grp_singles)>1))+1)]

#drop levels
l1o_steps=droplevels(l1o_steps)

#Set up vectors needed for loop
studies <- unique(l1o_steps[, c("study")])
groups <- unique(l1o_steps[, c("avail_group")])

#create empty data frame for output
l1o_out_all=data.frame()

#Loop through studies
for(j in 1:length(groups)) {
  
  #Remove all but group j
  #l1o_steps_rem<-subset(l1o_steps,study!=studies[j])
  print(paste0("group ",groups[j]))
  l1o_steps_rem<-subset(l1o_steps,avail_group==groups[j])
  l1o_steps_rem=droplevels(l1o_steps_rem)
  
  #set group
  group=as.integer(as.character(groups[j]))
  
  #Get vector of studies in group j
  Dat_Grp<-study_grp_avail[study_grp_avail$avail_group==group&study_grp_avail$Freq!=0,]
  studies_j<-Dat_Grp$study
  studies_j=droplevels(studies_j)
  
  #Loop through availability groups that contained study j
  for(i in 1:length(studies_j)) {
   
    print(paste0(studies_j[i]))
    #subset leave1out to group i
    Dat_Trem<-l1o_steps_rem[(l1o_steps_rem$study!=studies_j[i]),]
    
    #drop levels
    Dat_Trem<-droplevels(Dat_Trem)
    
    #get summaries of step/lc counts
    NLCDT<-as.data.frame(table(Dat_Trem$nlcd,Dat_Trem$case_,dnn=c("nlcd","case")))
    NLCDT<-subset(NLCDT,NLCDT$case=="TRUE"&NLCDT$Freq>0)
    
    #this should be same between Dat_T/Dat_Trem
    NLCD_ct<-dplyr::count(NLCDT)
    
    if(NLCD_ct$n<1) {stop("no land classes in group")}
    
    SSF_out=contrSum_iSSF(Dat_Trem,NLCDT,"animal_p_stp","nlcd_only","modelandparms")
    Trem1=SSF_out$model
    Dat_Trem2=SSF_out$dat
    Trem_parms=SSF_out$params
    Trem_parms<-as.data.frame(Trem_parms)
    Trem_parms$params=rownames(Trem_parms)
    
    #Format Trem_parms
    Trem_parms$group=group
    Trem_parms$study=studies_j[i]
    
    #use predict function
    preds.clog_rem=predict(Trem1,Dat_Trem2,type="expected")
    
    #bind results to data frame
    Dat_Trem2=cbind(Dat_Trem2,preds.clog_rem)
    
    #compare MSE and R2
    MSE.clog_rem=sum(na.omit((Dat_Trem2$case_-Dat_Trem2$preds.clog_rem)^2))/nrow(Dat_Trem2) 
    R2.clogit_rem=stats::cor(Dat_Trem2$case_,Dat_Trem2$preds.clog_rem)^2 

    l1o_out=data.frame("group"=group,"study"=studies_j[i],"MSE"=MSE.clog_rem,"R2"=R2.clogit_rem)
    
    ###########
    if(i==1){
    #Run iSSF model with contra sums process
    #Get nlcd frequencies/levels to do contra sums
      Dat_T=l1o_steps_rem
      Dat_T<-droplevels(Dat_T)
      
    #use NLCD_ct count to get contra sum for nlcd
      #get summaries of step/lc counts
      NLCDT<-as.data.frame(table(Dat_T$nlcd,Dat_T$case_,dnn=c("nlcd","case")))
      NLCDT<-subset(NLCDT,NLCDT$case=="TRUE"&NLCDT$Freq>0)
      
    #Run clogit model
    SSF_out=contrSum_iSSF(Dat_T,NLCDT,"animal_p_stp","nlcd_only","modelandparms")
    T1=SSF_out$model
    Dat_T2=SSF_out$dat
    T_parms=SSF_out$params
    T_parms=as.data.frame(T_parms)
    T_parms$params=rownames(T_parms)
    
    #Format T_parms
    T_parms$group=group
    T_parms$study="full"
    
    #use predict function
    preds.clog=predict(T1,Dat_T2,type="expected")

    #bind results to data frame
    Dat_T2=cbind(Dat_T2,preds.clog)

    #compare MSE and R2
    MSE.clog=sum(na.omit((Dat_T2$case_-Dat_T2$preds.clog)^2))/nrow(Dat_T2) 
    R2.clogit=stats::cor(Dat_T2$case_,Dat_T2$preds.clog)^2 
    full_out=data.frame("group"=group,"study"="full","MSE"=MSE.clog,"R2"=R2.clogit)
    }
    
    #info for parm output:
    #group, full_or_study_left_out, rest of parm table
    
    if(i==1&j==1){
    l1o_out_all=rbind(full_out,l1o_out)
    parms_out_all=rbind(Trem_parms,T_parms)
    } else{
      #if full group already got MSE/R2 in output
      if((paste0(group,"full")%in%
          paste0(l1o_out_all$group,l1o_out_all$study))){
        l1o_out_all=rbind(l1o_out_all,l1o_out)
        
        #group, full_or_study_left_out, rest of parm table
        parms_out_all=rbind(parms_out_all,Trem_parms)
        
        #if full group i doesn't already exist in full df, add both
        } else{
        l1o_out_sub=rbind(full_out,l1o_out)
        l1o_out_all=rbind(l1o_out_all,l1o_out_sub)
        
        parms_out_sub=rbind(T_parms,Trem_parms)
        parms_out_all=rbind(parms_out_all,parms_out_sub)
        }
      
    }
    
    
    } #for loop group
} #for loop study


#Write SSF Full output?
#SSF_full

#Write out leave-one-out analysis tables
write.csv(l1o_out_all,file.path(objdir,"l1o_model_summaries.csv",fsep=.Platform$file.sep))
write.csv(parms_out_all,file.path(objdir,"l1o_parm_summaries.csv",fsep=.Platform$file.sep))

# Tables of study, ID, and group and land use type
Table5<- table(steps_g$study,steps_g$avail_group,dnn=c("study","group"))
write.csv(Table5,paste0(objdir,"StudyGroups.csv"))

Table2<-table(steps_g$nlcd,steps_g$avail_group,steps_g$case_,dnn=c("nlcd","avail_group","case"))
write.csv(Table2,paste0(objdir,"Group_NLCD.csv"))

Table3<-table(steps_g$animalnum,steps_g$avail_group,dnn=c("animalnum","avail_group"))
write.csv(Table3,paste0(objdir,"Group_ID.csv"))

