#Pipeline:
#MakeFalseSteps > MakeAvailabilityGroups > iSSF > IDsPigsSSF

#########################
######## Purpose ######## 
#########################

#Initializes the population matrix for spatial meta-population model
#Also includes option to initialize individual/group in single location on grid
#Option to incorporate heterogeneous landscape preference is in progress

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
library(purrr)

#set directory paths
indir=paste0(home,"2_Data/Input/") #initial input read from here
objdir=paste0(home,"2_Data/Objects/") #intermediate objects, lookup tables, etc. go here

#load data
geo <- read.csv(paste0(indir,"geolocsnatl.csv")) #geo=PigsTotal2
#NLCD <- raster(paste0(indir,"nlcd_2016_land_cover_l48_20210604")) #use this if not using aliasing
LCD <- raster(mactheknife::resolve_alias(paste0(indir,"nlcd_2016_land_cover_l48_20210604"))) #aliasing read method

#################################
######## Data formatting ######## 
#################################

# Fix MI pigs that ended up in MS
geo$state[geo$study=="PhD_MI"]<-"MI"

#convert to sf format
geosf <- st_as_sf(geo, coords = c("longitude", "latitude"), crs="EPSG:4326")

#Transform CRS to match NLCD raster
geosf<-st_transform(geosf, crs=crs(LCD))

#pull x/y coords in aea projection
geosf$x_ <- sapply(geosf$geometry, "[", 1)
geosf$y_ <- sapply(geosf$geometry, "[", 2)

#set times to posixct format
geosf$t_<-as.POSIXct(geosf$datetime, "%Y-%m-%d %H:%M:%S", tz="UTC")

# Set numeric unique id for each pig
geosf$animalnum<-unique_identifier(geosf$animalid)

geosf$hour=hour(geosf$t_)

#Set up data set for steps
steps <- geosf[, c("x_","y_","t_","animalnum","hour")]

#Summary of num steps by animalid
Alln <- steps %>%
  group_by(animalnum) %>%
  dplyr::summarise(
    n=n())

##########################################
######## Subsampling geolocations ######## 
##########################################

#iSSF setting up steps by every 6 hours as a separate track. 
#Loop each pig four times for each set of 6 hours.
#Looking at daily step length.

#get vector of animalids
ids <- unique(steps$animalnum)

#set up empty dataframe for output
full<-data.frame()

#Loop through numeric animalids
for(i in seq_along(ids)) {
  print(i)
  #get subset of geolocs for pig i
  dat_1<-subset(steps, animalnum==ids[[i]][1])
  
  #ensure datetime is in posixct format
  dat_1$t_<-as.POSIXct(dat_1$t_, format="%Y-%m-%d %H:%M:%S")
  
  #stop if any NA datetimes
  if(any(is.na(dat_1$t_))){
    stop(paste0("NAs found in datetimes for animalnum ",ids[i]))
  }
  
  #convert to amt track format
  tr1<-amt::make_track(dat_1,x_,y_,t_, hour=hour, id=animalnum)
  
  #resample to 6 hr intervals to match all tracks to same sampling rate
  tr1 <- tr1 %>% track_resample(rate = hours(6), tolerance = minutes(60))
  
  #get subset for each period
  trp1<-subset(tr1,hour=="23"|hour=="22"|hour=="21"|hour=="20"|hour=="19"|hour=="18")
  trp2<-subset(tr1,hour=="17"|hour=="16"|hour=="15"|hour=="14"|hour=="13"|hour=="12")
  trp3<-subset(tr1,hour=="11"|hour=="10"|hour=="9"|hour=="8"|hour=="7"|hour=="6")
  trp4<-subset(tr1,hour=="0"|hour=="1"|hour=="2"|hour=="3"|hour=="4"|hour=="5")
  
  #make list of trp's, loop through
  trp.list=list(trp1,trp2,trp3,trp4)
  
  for(t in 1:length(trp.list)){
    trp=trp.list[[t]]
  
    if(nrow(trp)>1){  
    #convert to steps format, not bursted
    trp_steps<-trp%>%steps(.,diff_time_units="hours")
    
    #rejoin info
    trp_steps<-trp_steps%>%mutate(animalnum=ids[i]) #pull animalnum back
    trp_steps<-trp_steps%>%dplyr::mutate(period=t) #add period
    
    #bind together
    if(t==1&i==1){
      full=trp_steps
    } else{
      full<-rbind(full,trp_steps)
    }

    } else{ #else, if nrow(trp) not > 1
      if(t==1&i==1){ 
      full<-data.frame()
      }
    }
    
  } #end for loop through list of trp's

} #end for loop through numeric animalids

#Subsampling to daily steps
#Subset hour intervals between 21-27 hours (3 hour window on either side of 24)
realsteps<-subset(full, full$dt_<27.0001 & full$dt_>20.9999)

nrow(geo) #1904392 all geolocs
nrow(full) #234487 subsampled to 6 hrs
nrow(realsteps) #148776 subsampled to daily steps

#Use this to remove those with less than 10 steps
rs.summary <- realsteps %>% #was AllFull
  dplyr::group_by(animalnum,period) %>% 
  dplyr::summarise(
    n=n(),
  ) %>%
  as.data.frame()

#summarise ids with <10 steps
idslt10=rs.summary[rs.summary$n<10,]$animalnum

#remove ids with <10 steps
realsteps_gt10=realsteps[!(realsteps$animalnum%in%idslt10),]

#Summaries
length(ids) #563 ids before resampling
length(unique(realsteps$animalnum)) #543 ids after resampling
length(unique(realsteps_gt10$animalnum)) #464 ids after removing pigs without at least 10 steps

#write out subsampled steps
saveRDS(realsteps_gt10, paste0(objdir,"realsteps_subsampled.rds"))

####################################
######## Create False steps ######## 
####################################
#Need edit this to loop through pig/periods
#might need to undo steps by burst, just do by steps

#get smallest nonzero step length
#smallest_nonzero_sl=min(realsteps_gt10$sl_[which(realsteps_gt10$sl_>0)])

#Add smallest nonzero step length to realsteps_gt10 to try to force random steps to work where sl=0
#realsteps_gt10$sl_[which(realsteps_gt10$sl_==0)]=smallest_nonzero_sl

#ADD 15 random steps
#Loop through pig/period, get random steps and covs
#create unique ID for each pig/period
realsteps_gt10$pigperID<-unique_identifier(realsteps_gt10,fields=c("animalnum","period"))
pigperIDs=unique(realsteps_gt10$pigperID)

for(p in 1:length(unique(pigperIDs))){
  print(p)
  dat=realsteps_gt10[realsteps_gt10$pigperID==pigperIDs[p],]
  #where direction_p is na, take previous step direction_p
  #dat$direction_p[is.na(dat$direction_p)]<-dat$direction_p[(which(is.na(dat$direction_p))-1)]
  dat$direction_p[is.na(dat$direction_p)]<-runif(length((which(is.na(dat$direction_p)))),min=-3.14,max=3.14)
  
    #Tested issue around NA random steps
    #They are produced in the step following a step with 0 length
    #step with 0 length still gets random steps, bc it is fitted with >0 length within the function
    #NAs are propagated due to something happening in the random_steps function,
    #processing of output when direction_p=NA in the previous step due to 0 step length, 
    #even though direction_p is not actually used in any of the calculations for random step generation.
    #developer for amt suggested selecting random value for turn angle in these cases
  dat2<-dat%>%random_steps(n=15)%>%
                          extract_covariates(LCD)
  if(any(dat2$Height==0)){
    stop("nonsensical nlcd value")
  }
  #NAs happen when step has direction_p=NA
  #dat=dat[-which(colnames(dat)=="direction_p")]
  
  if(p==1){
    rf_all=dat2
  } else{
    rf_all=rbind(rf_all,dat2)
  }
}

colnames(rf_all)[ncol(rf_all)]<-"nlcd"

#should no longer be NAs with direction_p=1 fix
CheckNA<-subset(rf_all,is.na(rf_all$nlcd)==TRUE)
nrow(CheckNA) #0

# Create new id for unique animal, step id and time of day (period)
rf_all$animal_p_stp<-unique_identifier(rf_all, fields = c("animalnum","period","step_id_"), start_from = 1)

#set NLCD to factor class
rf_all$nlcd<-as.factor(rf_all$nlcd)

#Pull year from date
rf_all$year<-year(rf_all$t1_)

#Pull back info on state, sex, study and original animalid for each animal
PigsInfo  <-subset(geosf,select=c("animalnum","sex","study","state","animalid"))
PigsInfo <- st_drop_geometry(PigsInfo)
PigsInfo <- unique(PigsInfo)

#save for later
saveRDS(PigsInfo,paste0(objdir,"PigsInfo.rds"))

#Pull animal metadata back to false steps df
rf_all<-merge(rf_all,PigsInfo, by="animalnum")
nrow(rf_all)

#Write false steps df
saveRDS(rf_all, paste0(objdir,"real_false_steps.rds"))


