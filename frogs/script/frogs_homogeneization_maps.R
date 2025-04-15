##Frogs maps - homogeneization with trees datasets 
## Tali Guez - 2025

rm(list = ls())
gc() #Garbage collection

#### Librairies ####
#library(raster)
library(terra)
library(tidyverse)

#### Initialization of the raster ####
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/OR_AUC_10k")
mocklayer<-rast("AF_10k_AF_Scinax_alter_OR_AUC_T10.asc")
mocklayer<-init(mocklayer,"cell") #??
names(mocklayer) <- "Grid"
r<-mocklayer
values(r)<-NA

###### Add indices in maps ########## 
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/results")
all_metrics<-read_csv("REV_obs_results_Frogs_comb_6.csv") #change with the right file 1 to 10

###Create all metric rasters for observed
pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r

pd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$pd #
mpd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$mpd #plot(mpd_ras)
raoq_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$raoq #plot(raoq_ras)
SR_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$n #plot(SR_ras)
fdr_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fdr #plot(fdr_ras)
fd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fd #plot(fd_ras)

plot(SR_ras)
plot(pd_ras)
plot(mpd_ras)
plot(fd_ras)
plot(fdr_ras)
plot(raoq_ras) #need to take off the pixels where there is less then 3 species

### More than 3 species by pixels - TEST ###
combined_ras <- c(SR_ras, raoq_ras, fdr_ras, fd_ras) #, pd_ras, mpd_ras
names(combined_ras) <- c("Species_Richness", "RaoQ", "FDR", "FD") #, "PD", "MPD"
plot(combined_ras)
combined_ras[SR_ras < 4] <- NA
plot(combined_ras) #compare to verify
####

# Remplacer les pixels avec moins de 3 espèces (exclus) par NA
plot(raoq_ras)
raoq_ras[SR_ras < 4] <- NA
plot(raoq_ras) #compare to verify

#### Write rasters ####
####TEST files
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_homogeneization") 
#writeRaster(combined_ras,"combined_test.grd",filetype="RRASTER", overwrite=TRUE) #

writeRaster(pd_ras,filename="comb10/PD_comb10.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(mpd_ras,"comb10/MPD_comb10.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(SR_ras,"comb10/SR_comb10.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(raoq_ras,"comb6/RaoQ_comb6_V2.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fdr_ras,"comb10/fdr_comb10.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fd_ras,"comb10/fd_comb10.grd",filetype="RRASTER", overwrite=TRUE) #
