
########################################### LIBRAIRIES ###########################################################################################################
#### For Figure 4 from article
library(tidyverse) #careful some functions can overlap
library(gridExtra) #combine plot

#### For Residuals
#library(FD) #to compute pd
library(terra)
library(mgcv) #gam
library(lme4)


################################################ Residuals GAM ######################################################################
rm(list=ls())

##### Load diversity rasters (10km resolution) ######
setwd("~/Master/M2/Internship_M2/analyse/figures/maps_V3")

pd_raster<-rast("hylids_PD10k_V3.tif")
sr_raster<-rast("richness_hylids_complete_0_083_V3.tif")
#for comb1
fdRic_raster<-rast("hylid_functional_richness10k_cont_V3.tif")
fdDis_raster<-rast("hylid_functional_dispersion10k_cont_V3.tif")

sr_raster_res <- resample(sr_raster, fdRic_raster, method = "near")


stack_SR_PD<-c(pd_raster,sr_raster_res) #merge=combine SpatRasters with different extents (but same origin and resolution)
stack_SR_FD<-c(fdRic_raster,sr_raster_res)
stack_PD_FD<-c(fdRic_raster,pd_raster)
stack_PD_FDdis<-c(fdDis_raster,pd_raster)

#SR and PD
srpd_df<-as.data.frame(stack_SR_PD)
names(srpd_df)[c(1,2)]<- c("PD", "SR")
srpd_df$grid<-1:length(srpd_df$PD)
srpd_df_clean<-na.omit(srpd_df) #ex srpd_df1

#SR and FDric
srfd_df<-as.data.frame(stack_SR_FD)
names(srfd_df)[c(1,2)]<- c("FdRich", "SR")
srfd_df$grid<-1:length(srfd_df$FdRich)
srfd_df_clean<-na.omit(srfd_df) #ex srfd_df1

#PD and FDric
pdfd_df<-as.data.frame(stack_PD_FD)
names(pdfd_df)[c(1,2)]<- c("FdRich", "PD")
pdfd_df$grid<-1:length(pdfd_df$PD)
pdfd_df_clean<-na.omit(pdfd_df)  #ex pdfd_df1

#PD and FDisp
pdfddis_df<-as.data.frame(stack_PD_FDdis)
names(pdfddis_df)[c(1,2)]<- c("FdDisp", "PD")
pdfddis_df$grid<-1:length(pdfddis_df$FdDisp)
pdfddis_df_clean<-na.omit(pdfddis_df)

###LINEAR 
#Generalized additive model GAM 
linear_model_sr_pd<-glm(PD~SR,data=srpd_df_clean)
linear_model_sr_fd<-glm(FdRich~SR,data=srfd_df_clean)
linear_model_pd_fd<-glm(FdRich~PD,data=pdfd_df_clean)
linear_model_pd_fddis<-glm(FdDisp~PD,data=pdfddis_df_clean)

hist(linear_model_sr_pd$residuals) #not normal
plot(srpd_df_clean$SR,linear_model_sr_pd$residuals) #variance not constant
plot(linear_model_sr_pd, which=1) #no independance in data
plot(linear_model_sr_pd, which=2) #not normally distrib

hist(linear_model_sr_fd$residuals) #not normal but better
plot(srfd_df_clean$SR,linear_model_sr_fd$residuals) #variance not constant
plot(linear_model_sr_fd, which=1) #no independance in data
plot(linear_model_sr_fd, which=2) #not normally distrib

hist(linear_model_pd_fd$residuals) #not normal 
plot(pdfd_df_clean$PD,linear_model_pd_fd$residuals) #variance not constant
plot(linear_model_pd_fd, which=1) #no independance in data
plot(linear_model_pd_fd, which=2) #not normally distrib


hist(linear_model_pd_fddis$residuals) #not normal 
plot(pdfddis_df_clean$PD,linear_model_pd_fddis$residuals) #variance not constant but less correlation in data 
plot(linear_model_pd_fddis, which=1) #no homogeneity + no independance in data =  plot(linear_model_pd_fddis$fitted.values, linear_model_pd_fddis$residuals) #not normal 
plot(linear_model_pd_fddis, which=2) #not normally distrib


#### Generate regression models ####
#Generalized additive model GAM 
gam_model_sr_pd<-gam(PD~te(SR),data=srpd_df_clean)
gam_model_sr_fd<-gam(FdRich~te(SR),data=srfd_df_clean)
gam_model_pd_fd<-gam(FdRich~te(PD),data=pdfd_df_clean)
gam_model_pd_fddis<-gam(FdDisp ~te(PD),data=pdfddis_df_clean)



#### GAM RASTERS ####
#Compute and save residuals from gam models
gam_residuals_sr_pd<-resid(gam_model_sr_pd)
gam_residuals_sr_fd<-resid(gam_model_sr_fd)
gam_residuals_pd_fd<-resid(gam_model_pd_fd)
gam_residuals_pd_fddis<-resid(gam_model_pd_fddis)

#New datasets
srpd_df_resid<-srpd_df #with NA
srfd_df_resid<-srfd_df
pdfd_df_resid<-pdfd_df 
pdfddis_df_resid<-pdfddis_df  

###ADD GAM residuals to data frame 
#Add residuals to dataframe OK
srpd_df_resid2<-srpd_df_resid
srpd_df_resid2$GAM_residuals<-rep(NA,length(srpd_df[,1]))
srpd_df_resid2$GAM_residuals[srpd_df_clean$grid]<-gam_residuals_sr_pd

#add residuals of fd to data frame
srfd_df_resid2<-srfd_df_resid
srfd_df_resid2$GAM_residuals<-rep(NA,length(srfd_df[,1]))
srfd_df_resid2$GAM_residuals[srfd_df_clean$grid]<-gam_residuals_sr_fd

#add residuals of pd/fd to data frame
pdfd_df_resid2<-pdfd_df_resid
pdfd_df_resid2$GAM_residuals<-rep(NA,length(pdfd_df[,1]))
pdfd_df_resid2$GAM_residuals[pdfd_df_clean$grid]<-gam_residuals_pd_fd

#add residuals of pd/fddis to data frame ##check
pdfddis_df_resid2<-pdfddis_df_resid
pdfddis_df_resid2$GAM_residuals<-rep(NA,length(pdfddis_df[,1]))
pdfddis_df_resid2$GAM_residuals[pdfddis_df_clean$grid]<-gam_residuals_pd_fddis

####Create raster 
setwd("~/Master/M2/Internship_M2/analyse/datas/cleaned/residuals")

#Generate raster with residual values for PD~TD
gam_residuals_sr_pd_raster<-pd_raster[[1]] #reference raster #gam_residuals_sr_pd_raster<-raster(pd_raster)
values(gam_residuals_sr_pd_raster)<-NA #replace values by NA
srpd_df_resid3<-srpd_df_resid2 #new table
srpd_df_resid3$grid<-rownames(srpd_df_resid2) #change grid column into rownames
srpd_df_resid3$grid<- as.integer(srpd_df_resid3$grid) #change into integer to be read by the raster
gam_residuals_sr_pd_raster[srpd_df_resid3$grid]<-srpd_df_resid3$GAM_residuals #add residuals to raster
names(gam_residuals_sr_pd_raster) <- "GAM_residuals" #change name of the layer for clarity
plot(gam_residuals_sr_pd_raster, main="Residuals PD~TD") #checking

name_group<-"gam_residuals_SR_PD_10km.tif" #gam_residuals_SR_PD_1km_T10
writeRaster(gam_residuals_sr_pd_raster,name_group, filetype = "GTiff", overwrite=TRUE) #raster creation

#Generate raster with residual values for FDRich~TD
gam_residuals_sr_fd_raster<-pd_raster[[1]] #raster(pd_raster)
values(gam_residuals_sr_fd_raster)<-NA
srfd_df_resid3<-srfd_df_resid2
srfd_df_resid3$grid<-rownames(srfd_df_resid2)
srfd_df_resid3$grid<- as.integer(srfd_df_resid3$grid)
gam_residuals_sr_fd_raster[srfd_df_resid3$grid]<-srfd_df_resid3$GAM_residuals 
names(gam_residuals_sr_fd_raster) <- "GAM_residuals" #change name of the layer for clarity
plot(gam_residuals_sr_fd_raster, main="Residuals FDRich~TD") #checking

name_group<-"gam_residuals_SR_FD_10km.tif" #"gam_residuals_SR_FD_1km_T10"
writeRaster(gam_residuals_sr_fd_raster,name_group, filetype = "GTiff")

#Generate raster with residual values for FDRich~PD
gam_residuals_pd_fd_raster<-pd_raster[[1]]
values(gam_residuals_pd_fd_raster)<-NA

pdfd_df_resid3<-pdfd_df_resid2
pdfd_df_resid3$grid<-rownames(pdfd_df_resid2)
pdfd_df_resid3$grid<- as.integer(pdfd_df_resid3$grid)
gam_residuals_pd_fd_raster[pdfd_df_resid3$grid]<-pdfd_df_resid3$GAM_residuals 
names(gam_residuals_pd_fd_raster) <- "GAM_residuals" #change name of the layer for clarity
plot(gam_residuals_pd_fd_raster, main="Residuals FDRich~PD") #checking

name_group<-"gam_residuals_PD_FD_10km.tif" #gam_residuals_PD_FD_1km_T10
writeRaster(gam_residuals_pd_fd_raster,name_group, filetype = "GTiff", overwrite=TRUE)

#Generate raster with residual values for FDisp~PD
gam_residuals_pd_fddis_raster<-pd_raster[[1]]
values(gam_residuals_pd_fddis_raster)<-NA
pdfddis_df_resid3<-pdfddis_df_resid2
pdfddis_df_resid3$grid<-rownames(pdfddis_df_resid2)
pdfddis_df_resid3$grid<- as.integer(pdfddis_df_resid3$grid)
gam_residuals_pd_fddis_raster[pdfddis_df_resid3$grid]<-pdfddis_df_resid3$GAM_residuals 
names(gam_residuals_pd_fddis_raster) <- "GAM_residuals" #change name of the layer for clarity
plot(gam_residuals_pd_fddis_raster, main="Residuals FDRich~PD") #checking

name_group<-"gam_residuals_PD_FDdis_10km.tif" #"gam_residuals_PD_FDdis_1km_T10"
writeRaster(gam_residuals_pd_fddis_raster,name_group, filetype = "GTiff", overwrite=TRUE)
