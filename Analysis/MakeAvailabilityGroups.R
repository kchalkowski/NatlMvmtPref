#Pipeline:
#MakeFalseSteps > MakeAvailabilityGroups > iSSF > IDsPigsSSF

#########################
######## Purpose ######## 
#########################

#Set availability groups based on nlcd land class in true/false steps

##############################
######## Script Setup ######## 
##############################

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
indir=paste0(home,"2_Data/Input/") #initial input read from here
objdir=paste0(home,"2_Data/Objects/") #intermediate objects, lookup tables, etc. go here
outdir=paste0(home,"4_Output/") #intermediate objects, lookup tables, etc. go here

#load data
steps=readRDS(paste0(objdir,"real_false_steps.rds"))

#################################
######## Data formatting ######## 
#################################

#set all developed to same level
steps$nlcd[steps$nlcd==22]<-21
steps$nlcd[steps$nlcd==23]<-21
steps$nlcd[steps$nlcd==24]<-21

#Label NLCD present in data set
steps<-steps %>%
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

#Start sub-setting pigs by access to land types
Ids2b<-as.data.frame(table(steps$animalnum,steps$nlcd_str,steps$case_))
colnames(Ids2b)<-c("animalnum","nlcd","stp_typ","freq")

#separate false/true to isolate frequencies
Fal<-subset(Ids2b,Ids2b$stp_typ=="FALSE")
colnames(Fal)[which(colnames(Fal)=="freq")]<-"freq_false"
Tr<-subset(Ids2b,Ids2b$stp_typ=="TRUE")
colnames(Tr)[which(colnames(Tr)=="freq")]<-"freq_true"

#Remove step type, not needed
Fal=Fal[,c(1,2,4)]
Tr=Tr[,c(1,2,4)]

#merge back to get wide format
Access=merge(Fal,Tr,by=c("animalnum","nlcd"))

#Set up dummy coding for false steps
Access$avail=0
Access$avail[Access$freq_false>0|Access$freq_true>0]<-1

#Set up dummy coding for true/used steps
Access$used=0
Access$used[Access$freq_true>0]<-1

#subset needed cols
Available <-subset(Access,select=c(animalnum,nlcd,avail))

#get wide version
Available_w<-spread(Available,key=nlcd,value=avail)

#Group_Access_Key<-unique(Available_w[,14:17])
#View(Group_Access_Key)

#Make general columns for forest, wetland and urban
# Make new variable for pigs that have access to forest, urban areas and wetlands (Crops already a single land type)
Available_w$forest<-0
Available_w$forest[Available_w$Deciduous=="1"|
                     Available_w$Evergreen=="1"|
                     Available_w$Mixed=="1"]<-1

Available_w$wetland<-0
Available_w$wetland[Available_w$Woody_Wetlands=="1"|
                      Available_w$Emergent_Herbaceous_Wetlands=="1"]<-1

Available_w$urban<-0
Available_w$urban[Available_w$Developed=="1"]<-1

#Assign group number
#Need to verify with Lexi what class IDs supposed to be here
#four classes, considered sig for pigs based on past lit
Available_w$avail_group<-
  unique_identifier(Available_w,
                    c(which(colnames(Available_w)=="Crops"),
                      which(colnames(Available_w)=="forest"),
                      which(colnames(Available_w)=="wetland"),
                      which(colnames(Available_w)=="urban")))

#Make as factor, order levels
Available_w$avail_group<-factor(Available_w$avail_group,levels=as.character(1:length(unique(Available_w$avail_group))))

#Get counts of pigs in each avail group
Avail_Group_Counts<-as.data.frame(table(Available_w$avail_group))

colnames(Avail_Group_Counts)<-c("avail_group","freq")

#Write out availability group df
#saveRDS(Available_w,paste0("Avail_Group_df.rds"))

#Get key for animalnum/avail group matching
Avail_Group_Key <- Available_w[, c("animalnum","avail_group")]

#Get key for breakdown of groups/availabiltiy by land class
Avail_Group_Key_LC<-unique(Available_w[,c(17,3,14:16)])

#merge availability group IDs with steps to loop through groups
steps_g<-merge(steps,Avail_Group_Key, by="animalnum")

#Order steps by group
#convert group to numeric temporarily
steps_g$avail_group<-as.integer(steps_g$avail_group)

#Order by group, then by animalnum, then by t
steps_g=steps_g[with(steps_g,order(avail_group,animalnum,t1_,case_)),]

#change avail group back to factor
steps_g$avail_group<-factor(steps_g$avail_group,levels=as.character(1:length(unique(steps_g$avail_group))))

#Data for which pig is in which group for table
PigInfo_Group <- steps_g[, c("animalnum","state","study","sex","avail_group")]
PigInfo_Group<-unique(PigInfo_Group)

#Get unique groups
groups <- unique(steps_g[, c("avail_group")])

#Get reversed order nlcd
steps_g$nlcd_rev<-fct_rev(steps_g$nlcd)

#Write out objects needed for IDsPigsSSF.R
saveRDS(steps_g,paste0(objdir,"steps_grouped.rds"))

