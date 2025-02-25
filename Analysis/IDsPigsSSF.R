# Pipeline --------------------------------------------------------------------

#MakeFalseSteps > MakeAvailabilityGroups > iSSF > IDsPigsSSF > PlottingSSFOutput

# Purpose --------------------------------------------------------------------

#Run iSSF by looping through individual animals

# Script Setup --------------------------------------------------------------------

#set directory paths
indir=file.path(home,"2_Data","Input") #initial input read from here
objdir=file.path(home,"2_Data","Objects") #intermediate objects, lookup tables, etc. go here
outdir=file.path(home,"4_Output") #intermediate objects, lookup tables, etc. go here

#input file
steps_g=readRDS(file.path(objdir,"steps_grouped.rds"))

#load libraries
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
library(fs)
library(stringr)

# Loop through pig/periods -----------------------------------------------------

#Set up loop by pig/periods for each conditional logistic model using contr.sum. 
groups=unique(steps_g$avail_group)

all_parms_g=data.frame()
for(j in 1:length(groups)) {
  print(paste0("avail group ",groups[j]))
  Group_T<-steps_g[steps_g$avail_group==groups[j],]
  Group_T<-droplevels(Group_T)
  
  #get all pig period groupings
  pp_group <- unique(Group_T[, c("pigperID")])

  for(i in 1:length(pp_group)) {
    
  Dat_T=steps_g[steps_g$pigperID==pp_group[i],]
      Dat_T<-droplevels(Dat_T)
      NLCDT<-as.data.frame(table(Dat_T$nlcd,Dat_T$case_,dnn=c("nlcd","case")))
      NLCDT<-subset(NLCDT,NLCDT$case=="TRUE"&NLCDT$Freq>0)
      NLCDT_ct<-dplyr::count(NLCDT)
      
      if(NLCDT_ct$n<1) {stop("no land classes in group")}
      if(NLCDT_ct$n==1) {
        warning(paste0("only 1 land class in group, skipping SSF for pig period ",pp_group[i],", group ",j))
        if(j==1&i==1){
          all_parms_g=data.frame()
        }
        
        } else{
      
      
      #need set strata (animal_p_stp)
      SSF_out_pp=contrSum_iSSF(Dat_T,NLCDT,"animal_p_stp","nlcd_only","modelandparms")
      
      #Format parm table
      pp_parms=SSF_out_pp$params
      pp_parms<-as.data.frame(pp_parms)
      pp_parms$nlcd=rownames(pp_parms)
      pp_parms$group=groups[j]
      pp_parms$study=Dat_T$study[1] #need change how study fed in
      pp_parms$animalnum=Dat_T$animalnum[1]
      pp_parms$pig_period=pp_group[i] #maybe change this too, pig/period separate
      
      if(i==1){
        all_parms_g=pp_parms
      } else{
        all_parms_g=rbind(all_parms_g,pp_parms)
      }
      
        } #close if NLCDT_ct ==1 else bracket
      
} #end for loop through pig periods within group
  
  #combine all_parms across groups, full table
  if(j==1){
    all_parms_total=all_parms_g
  } else{
    all_parms_total=rbind(all_parms_total,all_parms_g)
  }
  
  #summarize for group
  GroupTable<-all_parms_g %>%
    group_by(nlcd) %>%
    dplyr::summarize(Mean = mean(coef),
                     Median=median(coef),
                     min=min(coef),
                     max=max(coef),
                     q05=quantile(coef,0.05,na.rm=TRUE),
                     q95=quantile(coef,0.95,na.rm=TRUE),
                     n=n())
  
  GroupTable$avail_group=groups[j]
  
  #get subset with only significant vals
  Groupsig<-subset(all_parms_g,all_parms_g$`Pr(>|z|)`<0.05)
  
  #Do separate summary for sig values only
  if(nrow(Groupsig)>0){
  GroupTablesig<-Groupsig %>%
    group_by(nlcd) %>%
    dplyr::summarize(Mean = mean(coef),
                     Median=median(coef),
                     min=min(coef),
                     max=max(coef),
                     q05=quantile(coef,0.05,na.rm=TRUE),
                     q95=quantile(coef,0.95,na.rm=TRUE),
                     n=n())
  
  GroupTablesig$avail_group=groups[j]
  } else{
    GroupTablesig=tibble()
  }
  
  #combine as go through loop
  if(j==1){
    AllGroupParms=GroupTable
    AllGroupParmsSig=GroupTablesig
  } else{
    AllGroupParms=rbind(AllGroupParms,GroupTable)
    AllGroupParmsSig=rbind(AllGroupParmsSig,GroupTablesig)
  }
  
}

# Format outputs -----------------------------------------------------

#relabel
AllGroupParms<-AllGroupParms %>%
  mutate(nlcd_str = case_when(nlcd == 11 ~ 'Open_Water',
                              nlcd== 21 ~ 'Developed',
                              nlcd == 22 ~ 'Developed',
                              nlcd == 23 ~ 'Developed',
                              nlcd == 24 ~ 'Developed',
                              nlcd == 31 ~ 'Barren',
                              nlcd == 41 ~ 'Deciduous',
                              nlcd == 42 ~ 'Evergreen',
                              nlcd == 43 ~ 'Mixed',
                              nlcd == 52 ~ 'Shrub',
                              nlcd == 71 ~ 'Grassland',
                              nlcd == 81 ~ 'Pasture',
                              nlcd == 82 ~ 'Crops',
                              nlcd == 90 ~ 'Woody_Wetlands',
                              nlcd == 95 ~ 'Emergent_Herbaceous_Wetlands'
  ))

AllGroupParmsSig<-AllGroupParmsSig %>%
  mutate(nlcd_str = case_when(nlcd == 11 ~ 'Open_Water',
                              nlcd== 21 ~ 'Developed',
                              nlcd == 22 ~ 'Developed',
                              nlcd == 23 ~ 'Developed',
                              nlcd == 24 ~ 'Developed',
                              nlcd == 31 ~ 'Barren',
                              nlcd == 41 ~ 'Deciduous',
                              nlcd == 42 ~ 'Evergreen',
                              nlcd == 43 ~ 'Mixed',
                              nlcd == 52 ~ 'Shrub',
                              nlcd == 71 ~ 'Grassland',
                              nlcd == 81 ~ 'Pasture',
                              nlcd == 82 ~ 'Crops',
                              nlcd == 90 ~ 'Woody_Wetlands',
                              nlcd == 95 ~ 'Emergent_Herbaceous_Wetlands'
  ))

all_parms_total<-all_parms_total %>%
  mutate(nlcd_str = case_when(nlcd == 11 ~ 'Open_Water',
                              nlcd== 21 ~ 'Developed',
                              nlcd == 22 ~ 'Developed',
                              nlcd == 23 ~ 'Developed',
                              nlcd == 24 ~ 'Developed',
                              nlcd == 31 ~ 'Barren',
                              nlcd == 41 ~ 'Deciduous',
                              nlcd == 42 ~ 'Evergreen',
                              nlcd == 43 ~ 'Mixed',
                              nlcd == 52 ~ 'Shrub',
                              nlcd == 71 ~ 'Grassland',
                              nlcd == 81 ~ 'Pasture',
                              nlcd == 82 ~ 'Crops',
                              nlcd == 90 ~ 'Woody_Wetlands',
                              nlcd == 95 ~ 'Emergent_Herbaceous_Wetlands'
  ))


# Write outputs -----------------------------------------------------

#write out averaged parm tables
write.csv(AllGroupParms,paste0(objdir,"AllGroupParms.csv"))
write.csv(AllGroupParmsSig,paste0(objdir,"AllGroupParmsSig.csv"))
write.csv(all_parms_total,paste0(objdir,"all_parms_total.csv"))

#rds versions for next script
saveRDS(AllGroupParms,paste0(objdir,"AllGroupParms.rds"))
saveRDS(AllGroupParmsSig,paste0(objdir,"AllGroupParmsSig.rds"))
saveRDS(all_parms_total,paste0(objdir,"all_parms_total.rds"))





