#####################################################################################
#########MAPPING OBSERVED METRICS AND RESULTS FROM NULL MODELS AND SIMULATIONS #######
##############################June 2023###############################################
#####################################################################################
####Mapping dimensions of diversity based on results from index calculation, simulations and bootstrapping
#Modified by Tali

rm(list = ls())

#library(raster)
library(terra)
library(tidyverse)
#setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Manuscript/NatCommsRev/")
setwd("~/Master/M2/Internship_M2/analyse/trees/datas")
###Create empty rasters based on richness map
mocklayer <- rast("raw/test_base_r.grd") #raster("SR_all_trees_observed.tif")
#mocklayer <- raster::raster("raw/SR_all_trees_observed.tif") #raster("SR_all_trees_observed.tif")
mocklayer<-init(mocklayer,"cell") #??
names(mocklayer) <- "Grid"
r<-mocklayer
values(r)<-NA
#crs(r) <- "EPSG:3857"
#r<-rast(ext(r), resolution = c(0.8333333333333333703, 0.8333333333317306524), crs = "EPSG:3857")
#r<-resample(r, r, method = "bilinear")

#wgs84 so epsg 3857
#Pixel Size	0.8333333333333333703,-0.8333333333317306524
############################
######OBERVED DATA########## 
########################### VS angio only script after 

##Choose the right file: 
#all_metrics<-read.csv("REV_obs_results_200.csv") #All observed results #old
all_metrics<-read_csv("clean/results/REV_obs_results_ANGIO_eco_based_200.csv") #eco based traits
all_metrics<-read_csv("clean/results/REV_obs_results_ANGIO_random_1_200.csv") #random
all_metrics<-read_csv("clean/results/REV_obs_results_ANGIO_random_10_200.csv") #random
all_metrics<-read_csv("clean/results/REV_obs_results_ANGIO_pca_based_200.csv") #pca 
all_metrics<-read_csv("clean/results/REV_obs_results_ANGIO_cluster_random_200.csv") #wihtin clusters random

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

plot(pd_ras)
plot(SR_ras)
plot(fd_ras)
plot(raoq_ras)
plot(fdr_ras)
plot(mpd_ras)

setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps/random_clusters") #eco_based_200 / randomX / random_clusters / PCA_selection
#randomX = randomX at 200 filtering
#writeRaster(x=pd_ras,filename="PD_eco_based.tif",filetype="GTiff", overwrite=TRUE) #PD_all_trees.tif
writeRaster(pd_ras,filename="PD_clusters.grd",filetype="RRASTER", overwrite=TRUE) #PD_all_trees.tif
writeRaster(mpd_ras,"MPD_clusters.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(SR_ras,"SR_clusters.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(raoq_ras,"RaoQ_clusters.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fdr_ras,"fdr_clusters.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fd_ras,"fd_clusters.grd",filetype="RRASTER", overwrite=TRUE) #

### More than 8 species by pixels - TEST ###
combined_ras <- c(SR_ras, raoq_ras, fdr_ras, fd_ras) #, pd_ras, mpd_ras
names(combined_ras) <- c("Species_Richness", "RaoQ", "FDR", "FD") #, "PD", "MPD"
plot(combined_ras)
combined_ras[SR_ras < 9] <- NA
plot(combined_ras) #compare to verify
####

# Remplacer les pixels avec moins de 3 espèces par NA
plot(raoq_ras)
raoq_ras[SR_ras < 9] <- NA
#fdr_ras[SR_ras < 9] <- NA
#fd_ras[SR_ras < 9] <- NA 
plot(raoq_ras)#compare to verify

####TEST
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps/test") #eco_based_200 / randomX / random_clusters / PCA_selection
writeRaster(combined_ras,"combined_test.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(raoq_ras,"raoq_test.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fdr_ras,"fdr_test.grd",filetype="RRASTER", overwrite=TRUE) #

##write again the rasters for RaoQ
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps/PCA_selection") #eco_based_200 / randomX / random_clusters / PCA_selection
writeRaster(raoq_ras,"RaoQ_pca_9sp.grd",filetype="RRASTER", overwrite=TRUE) #
#RaoQ_clusters.grd
#_9sp

##### OLD SCRIPT ######################
###########angiosperms only###############
all_metrics<-read.csv("REV_obs_results_ANGIO_200.csv") #Observed Angio only


###Create all metric rasters for observed
pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r


pd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$pd
mpd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$mpd
raoq_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$raoq
SR_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$n
fdr_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fdr
fd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fd

writeRaster(pd_ras,"PD_angio_trees.tif",format="GTiff")
writeRaster(mpd_ras,"MPD_angio_trees.tif",format="GTiff")
writeRaster(SR_ras,"SR_angio_trees.tif",format="GTiff")
writeRaster(raoq_ras,"RaoQ_angio_trees.tif",format="GTiff")
writeRaster(fdr_ras,"fdr_angio_trees.tif",format="GTiff")
writeRaster(fd_ras,"fd_angio_trees.tif",format="GTiff")


############################
######BOOTSTRAP DATA (mean and CV)##########
###########################
all_metrics<-read.csv("REV_downsample_gatti_200.csv") 

pd_rasB<-r
mpd_rasB<-r
raoq_rasB<-r
SR_rasB<-r  
fdr_rasB<-r
fd_rasB<-r

cv_pd_rasB<-r
cv_mpd_rasB<-r
cv_raoq_rasB<-r
cv_fdr_rasB<-r
cv_fd_rasB<-r

boot<-subset(all_metrics,metric=="pd")
pd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_pd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="mpd")
mpd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_mpd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="raoq")
raoq_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_raoq_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="fdr")
fdr_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_fdr_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="fd")
fd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_fd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean

boot<-subset(all_metrics,metric=="n")
SR_rasB[as.numeric(boot$grid_id)]<-boot$mean

writeRaster(pd_rasB,"Bootstrap_PD_all_trees.tif",format="GTiff")
writeRaster(mpd_rasB,"Bootstrap_MPD_all_trees.tif",format="GTiff")
writeRaster(SR_rasB,"Bootstrap_SR_all_trees.tif",format="GTiff")
writeRaster(raoq_rasB,"Bootstrap_RaoQ_all_trees.tif",format="GTiff")
writeRaster(fdr_rasB,"Bootstrap_fdr_all_trees.tif",format="GTiff")
writeRaster(fd_rasB,"Bootstrap_fd_all_trees.tif",format="GTiff")

writeRaster(cv_pd_rasB,"Bootstrap_CV_PD_all_trees.tif",format="GTiff")
writeRaster(cv_mpd_rasB,"Bootstrap_CV_MPD_all_trees.tif",format="GTiff")
#writeRaster(cv_SR_rasB,"Bootstrap_CV_SR_all_trees.tif",format="GTiff")
writeRaster(cv_raoq_rasB,"Bootstrap_CV_RaoQ_all_trees.tif",format="GTiff")
writeRaster(cv_fdr_rasB,"Bootstrap_CV_fdr_all_trees.tif",format="GTiff")
writeRaster(cv_fd_rasB,"Bootstrap_CV_fd_all_trees.tif",format="GTiff")

plot(cv_pd_rasB)
plot(cv_fdr_rasB)
plot(cv_mpd_rasB)
plot(cv_raoq_rasB)


############################
######NULL MODEL DATA##########
###########################
observed<-read.csv("REV_obs_results_200.csv")
simulations<-read.csv("REV_null_results_200.csv") 
#Z as (obs-mean/sd)

pdSim<-subset(simulations,metric=="pd")
observed <- observed %>% left_join(pdSim, by = c("n"="n")) 
observed$Zpd <- (observed$pd-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd)
mpdSim<-subset(simulations,metric=="mpd")
observed <- observed %>% left_join(mpdSim, by = c("n"="n")) 
observed$Zmpd <- (observed$mpd-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd)

raoSim<-subset(simulations,metric=="raoq")
observed <- observed %>% left_join(raoSim, by = c("n"="n")) 
observed$Zraoq <- (observed$raoq-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd,Zraoq)

fdrSim<-subset(simulations,metric=="fdr")
observed <- observed %>% left_join(fdrSim, by = c("n"="n")) 
observed$Zfdr <- (observed$fdr-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd,Zraoq,Zfdr)

##write file with Z values for Random Forest
write.csv(observed,"REV_obs_results_200_zvals.csv")

Zpd<-r
Zmpd<-r
Zraoq<-r
Zfdr<-r

Zpd[as.numeric(observed$grid_id)]<-observed$Zpd
Zmpd[as.numeric(observed$grid_id)]<-observed$Zmpd
Zraoq[as.numeric(observed$grid_id)]<-observed$Zraoq
Zfdr[as.numeric(observed$grid_id)]<-observed$Zfdr

plot(Zpd)
####WriteRasters to files for maps

writeRaster(Zpd,"Zpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zmpd,"Zmpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zraoq,"Zraoq_all_trees_observed.tif",format="GTiff")
writeRaster(Zfdr,"Zfdr_all_trees_observed.tif",format="GTiff")

