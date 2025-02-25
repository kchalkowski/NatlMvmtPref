# Pipeline --------------------------------------------------------------------

#MakeFalseSteps > MakeAvailabilityGroups > iSSF > IDsPigsSSF > PlottingSSFOutput.R

# Purpose --------------------------------------------------------------------

#Run iSSF by looping through individual animals

# Script Setup --------------------------------------------------------------------

#set directory paths
indir=paste0(home,"2_Data/Input/") #initial input read from here
objdir=paste0(home,"2_Data/Objects/") #intermediate objects, lookup tables, etc. go here
outdir=paste0(home,"4_Output/") #intermediate objects, lookup tables, etc. go here

#input files
AllGroupParms=readRDS(paste0(objdir,"AllGroupParms.rds"))
AllGroupParmsSig=readRDS(paste0(objdir,"AllGroupParmsSig.rds"))
all_parms_total=readRDS(paste0(objdir,"all_parms_total.rds"))

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
library(hrbrthemes)
library(lme4)
library(sjPlot)


#################################################
######## Make dotplots for each LC/group ########
#################################################

ggplot(data=AllGroupParmsSig)+
  geom_segment(
    mapping=aes(x=q05,xend=q95,y=avail_group,yend=avail_group,color=avail_group),alpha=0.5)+
  #geom_point(mapping=aes(x=min,y=avail_group,color=avail_group),alpha=0.5)+
  #geom_point(mapping=aes(x=max,y=avail_group,color=avail_group),alpha=0.5)+
  facet_wrap(~nlcd_str)+
  geom_point(mapping=aes(x=Mean,y=avail_group,color=avail_group,size=n))+
  geom_vline(xintercept=0,linetype="dotdash",color="red",linewidth=0.6)+
  theme_ipsum()

ggplot(data=AllGroupParms)+
  geom_segment(
    mapping=aes(x=q05,xend=q95,y=avail_group,yend=avail_group,color=avail_group),alpha=0.5)+
  #geom_point(mapping=aes(x=min,y=avail_group,color=avail_group),alpha=0.5)+
  #geom_point(mapping=aes(x=max,y=avail_group,color=avail_group),alpha=0.5)+
  facet_wrap(~nlcd_str)+
  geom_point(mapping=aes(x=Mean,y=avail_group,color=avail_group,size=n))+
  geom_vline(xintercept=0,linetype="dotdash",color="red",linewidth=0.6)+
  theme_ipsum()

colours <- c("#ff5b42","#ff7530","#ffffff","#73e5ff","#4b9ff2")
colour_breaks <- c(-4,-1,0,1,4)

ggplot(AllGroupParmsSig, aes(avail_group, nlcd_str, fill=Mean)) + 
  geom_tile(linewidth=0.8)+theme_ipsum()+
  scale_fill_gradientn(
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = range(AllGroupParmsSig$Mean[!is.na(AllGroupParmsSig$Mean)])), 1))

View(AllGroupParmsSig)


