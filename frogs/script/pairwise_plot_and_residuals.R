################################################################################################################################################## My script
####################################### SCRIPT - Internship M2 ##########################################################################################
#########################################  Tali - Guez ############################################################################################
##################################################################################################################################################My script
#Amphibian 10k
#Pairwise comparison plot
#FIGURE 4 from article # GAM residuals

########################################### LIBRAIRIES ###########################################################################################################
#### For Figure 4 from article
library(tidyverse) #careful some functions can overlap
library(gridExtra) #combine plot

#### For Residuals
#library(FD) #to compute pd
library(terra)
library(mgcv) #gam

################################ Figure 4 ######################################################################
#Import files
PD_TD<-read.table("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/PD_TD_V4.txt", header=TRUE)
df_comb1<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3.txt", header=TRUE)
#df_comb2<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3_comb2.txt", header=TRUE)

setdiff(PD_TD$grilla, df_comb2$grilla)

###refaire les tableaux avec PD et TD!!!
df <- merge(PD_TD, df_comb1, by.x = "grilla", by.y = "grilla")
df<-df %>% 
  select(grilla,PD,SR,fdRic,fdDis)

plot(df$PD,df$fdRic)
#grilla: number of the layer from the layer
#Taxonomic diversity: "SR"
#Phylogenetic diversity:"PD"
#Functional richness: fdRic
#Functional dispersion:fdDis

p1<-ggplot(df, mapping=aes(x=SR, y=fdRic))+
  geom_point()+
  theme_bw() + #black and white 
  labs(x="Taxonomic diversity", y="Functional richness")+ 
  geom_smooth(method = "loess", fill='blue2', level=0.90)

p2<-ggplot(df, mapping=aes(x=SR, y=PD))+
  geom_point()+
  theme_bw() + #black and white 
  labs(x="Taxonomic diversity", y="Phylogentic diversity")+
  geom_smooth(method = "loess", fill='blue2', level=0.90)

p3<-ggplot(df, mapping=aes(x=PD, y=fdRic))+
  geom_point()+
  theme_bw() + #black and white 
  labs(x="Phylogentic diversity", y="Functional richness")+ 
  geom_smooth(method = "loess", fill='blue2', level=0.90)

p4<-ggplot(df, mapping=aes(x=PD, y=fdDis))+
  geom_point()+
  theme_bw() + #black and white 
  labs(x="Phylogentic diversity", y="Functional dispersion")+ 
  geom_smooth(method = "loess", fill='blue2', level=0.90)

comb_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2) #Combine 4 graphs on one image

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/pairwise_plots")
#Save it 
ggsave(
  filename = paste0("Combination_1_plots.png"),
  plot = comb_plot,
  width = 12, height = 8 ) # Ajuste la taille selon tes préférences


##### le reste
pathway <- "~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/"
df_list <- list() #create empty list

#create a list with each table of combinations of traits within (merge everything)
for (i in 2:10) {
  file_name <- paste0(pathway, "communities_hylids_and_fd10k_cont_bis_V3_comb", i, ".txt") #Built name file
  df_comb <- read.table(file_name, header = TRUE) #read table
  df1 <- merge(PD_TD, df_comb, by.x = "grilla", by.y = "grilla") #merge with PD_TD
  
  #Select only 5 column
  df1 <- df1 %>% 
    select(grilla, PD, SR, fdRic, fdDis)
  
  df_list[[paste0("df_final_comb", i)]] <- df1 #stock table into a list
}

#df_list[[1]][, 1]  #Par index
#df_final_list[[1]]$grilla  #Par nom, si c'est "grilla"


#Define the output directory for saving plots
path_plots <- "~/Master/M2/Internship_M2/analyse/frogs/figures/pairwise_plots/"

# Function to create and save the plots
generate_plots <- function(data, table_name, combination_number) {
  # Plot 1: Taxonomic diversity vs Functional richness
  plot1 <- ggplot(data, mapping=aes(x =SR, y=fdRic))+
    geom_point()+
    theme_bw()+
    labs(x="Taxonomic diversity", y="Functional richness")+
    geom_smooth(method = "loess", fill='blue2', level=0.90)
  
  #Save Plot 1
  #ggsave(filename = paste0(path_plots, table_name, "_TD_vs_FR.png"), plot = plot1)
  
  # Plot 2: Taxonomic diversity vs Phylogenetic diversity
  plot2 <- ggplot(data, mapping=aes(x=SR, y=PD))+
    geom_point()+
    theme_bw()+
    labs(x="Taxonomic diversity", y="Phylogenetic diversity")+
    geom_smooth(method = "loess", fill='blue2', level=0.90)
  
  # Save Plot 2
  #ggsave(filename = paste0(path_plots, table_name, "_TD_vs_PD.png"), plot = plot2)
  
  # Plot 3: Phylogenetic diversity VS Functional richness
  plot3 <- ggplot(data, mapping=aes(x=PD, y=fdRic))+
    geom_point()+
    theme_bw()+
    labs(x="Phylogenetic diversity", y="Functional richness")+
    geom_smooth(method = "loess", fill='blue2', level=0.90)
  
  # Save Plot 3
  #ggsave(filename = paste0(path_plots, table_name, "_PD_vs_FR.png"), plot = plot3)
  
  # Plot 4: Phylogenetic diversity VS Functional dispersion
  plot4 <- ggplot(data, mapping=aes(x=PD, y=fdDis))+
    geom_point()+
    theme_bw()+
    labs(x="Phylogenetic diversity", y="Functional dispersion")+
    geom_smooth(method = "loess", fill='blue2', level=0.90)
  
  # Save Plot 4
  #ggsave(filename = paste0(path_plots, table_name, "_PD_vs_Fdisp.png"), plot = plot4)
  
  
  combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2) #Combine 4 graphs on one image
  
  #title = paste("Traits Combination", combination_number), to add for each plot instead of doing combined plots
  
  #Save it 
  ggsave(
    filename = paste0(path_plots, "Combination_", combination_number, "_plots.png"),
    plot = combined_plot,
    width = 12, height = 8  # Ajuste la taille selon tes préférences
  )
}

# Loop through the list and generate plots for each table
for (i in seq_along(df_list)) {
  table_name <- names(df_list)[i] #Get table name (e.g., "df_comb2","df_comb3", ...)
  
  # Calculate the combination number based on the table name
  combination_number <- i + 1 # If the list starts at "df_comb2", shift the index by 1
  
  # Get the data
  data <- df_list[[i]]
  
  # Generate and save the plots (pass the combination number as 'i')
  generate_plots(data, table_name, combination_number = combination_number)
}

################################################ Residuals GAM ######################################################################
rm(list=ls())

##### Load diversity rasters (10km resolution) ######
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

pd_raster<-rast("hylids_PD10k_V3.tif")
sr_raster<-rast("richness_hylids_complete_0_083_V3.tif")
#for comb1
fdRic_raster<-rast("hylid_functional_richness10k_cont_V3.tif")
fdDis_raster<-rast("hylid_functional_dispersion10k_cont_V3.tif")

sr_raster_res <- resample(sr_raster, fdRic_raster, method = "near")

#checking
plot(sr_raster)
plot(sr_raster_res)

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

#### Generate regression models ####
#Generalized additive model GAM #with te 
gam_model_sr_pd_te<-gam(PD~te(SR),data=srpd_df_clean) #R-sq.(adj) =  0.994 
gam_model_sr_fd_te<-gam(FdRich~te(SR),data=srfd_df_clean) #R-sq.(adj) =  0.857
gam_model_pd_fd_te<-gam(FdRich~te(PD),data=pdfd_df_clean) #R-sq.(adj) =  0.854 
gam_model_pd_fddis_te<-gam(FdDisp~te(PD),data=pdfddis_df_clean) #R-sq.(adj) = 0.3  
#summary()

#linear
linear_model_sr_pd<-gam(PD~SR,data=srpd_df_clean)
linear_model_sr_fd<-gam(FdRich~SR,data=srfd_df_clean)
linear_model_pd_fd<-gam(FdRich~PD,data=pdfd_df_clean)
linear_model_pd_fddis<-gam(FdDisp~PD,data=pdfddis_df_clean)

#s()
gam_model_sr_pd_s<-gam(PD~s(SR),data=srpd_df_clean)
gam_model_sr_fd_s<-gam(FdRich~s(SR),data=srfd_df_clean)
gam_model_pd_fd_s<-gam(FdRich~s(PD),data=pdfd_df_clean)
gam_model_pd_fddis_s<-gam(FdDisp~s(PD),data=pdfddis_df_clean)

#ti()
gam_model_sr_pd_ti<-gam(PD~ti(SR),data=srpd_df_clean)
gam_model_sr_fd_ti<-gam(FdRich~ti(SR),data=srfd_df_clean)
gam_model_pd_fd_ti<-gam(FdRich~ti(PD),data=pdfd_df_clean)
gam_model_pd_fddis_ti<-gam(FdDisp~ti(PD),data=pdfddis_df_clean)

#ti()
gam_model_sr_pd_t2<-gam(PD~t2(SR),data=srpd_df_clean)
gam_model_sr_fd_t2<-gam(FdRich~t2(SR),data=srfd_df_clean)
gam_model_pd_fd_t2<-gam(FdRich~t2(PD),data=pdfd_df_clean)
gam_model_pd_fddis_t2<-gam(FdDisp~t2(PD),data=pdfddis_df_clean)

AIC(gam_model_sr_pd_te,linear_model_sr_pd,gam_model_sr_pd_s, gam_model_sr_pd_ti, gam_model_sr_pd_t2)
AIC(gam_model_sr_fd_te,linear_model_sr_fd,gam_model_sr_fd_s, gam_model_sr_fd_ti, gam_model_sr_fd_t2 )
AIC(gam_model_pd_fd_te,linear_model_pd_fd,gam_model_pd_fd_s, gam_model_pd_fd_ti, gam_model_pd_fd_t2 )
AIC(gam_model_pd_fddis_te,linear_model_pd_fddis,gam_model_pd_fddis_s, gam_model_pd_fddis_ti, gam_model_pd_fddis_t2 )

#AIC always best with s() #written again to not change the script
gam_model_sr_pd<-gam(PD~s(SR),data=srpd_df_clean)
gam_model_sr_fd<-gam(FdRich~s(SR),data=srfd_df_clean)
gam_model_pd_fd<-gam(FdRich~s(PD),data=pdfd_df_clean)
gam_model_pd_fddis<-gam(FdDisp~s(PD),data=pdfddis_df_clean)

####PLOTS pairwise comparison ####
preds <- predict(gam_model_sr_pd,se.fit=TRUE)
ggplot()+
  geom_point(data=srpd_df_clean,aes(x=SR,y=PD))+
  geom_line( aes(x=srpd_df_clean$SR, y=preds$fit), linewidth=1, col="blue")+
  #  geom_line(aes(x=linear_model_sr_pd$model$richness_hylids_complete_0_0083,y=linear_model_sr_pd$fitted.values),col="red")+
  theme_bw()+labs(x="Taxonomic diversity",y="Phylogenetic diversity")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

preds <- predict(gam_model_sr_fd,se.fit=TRUE)
ggplot()+
geom_point(data=srfd_df_clean,aes(x=SR,y=FdRich))+
  geom_line(aes(x=srfd_df_clean$SR, y=preds$fit), linewidth=1, col="blue")+
  #  geom_line(aes(x=linear_model_sr_fd$model$richness_hylids_complete_0_0083,y=linear_model_sr_fd$fitted.values),col="red")+
  theme_bw()+labs(x="Taxonomic diversity",y="Functional richness")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

preds <- predict(gam_model_pd_fd,se.fit=TRUE)
ggplot()+
  geom_point(data=pdfd_df_clean,aes(x=PD,y=FdRich))+
  geom_line(aes(x=pdfd_df_clean$PD, y=preds$fit), linewidth=1, col="blue")+
  #  geom_line(aes(x=linear_model_pd_fd$model$hylids_PD,y=linear_model_pd_fd$fitted.values),col="red")+
  theme_bw()+labs(x="Phylogenetic diversity",y="Functional richness")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

preds <- predict(gam_model_pd_fddis,se.fit=TRUE)
ggplot()+
  geom_point(data=pdfddis_df_clean,aes(x=PD,y=FdDisp))+
  geom_line( aes(x=pdfddis_df_clean$PD, y=preds$fit), linewidth=1, col="blue")+
  #  geom_line(aes(x=linear_model_pd_fddis$model$hylids_PD,y=linear_model_pd_fddis$fitted.values),col="red")+
  theme_bw()+labs(x="Phylogenetic diversity",y="Functional dispersion")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

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
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/residuals/with_s")

#Generate raster with residual values for PD~TD
gam_residuals_sr_pd_raster<-pd_raster[[1]] #reference raster #gam_residuals_sr_pd_raster<-raster(pd_raster)
values(gam_residuals_sr_pd_raster)<-NA #replace values by NA
srpd_df_resid3<-srpd_df_resid2 #new table
srpd_df_resid3$grid<-rownames(srpd_df_resid2) #change grid column into rownames
srpd_df_resid3$grid<- as.integer(srpd_df_resid3$grid) #change into integer to be read by the raster
gam_residuals_sr_pd_raster[srpd_df_resid3$grid]<-srpd_df_resid3$GAM_residuals #add residuals to raster
names(gam_residuals_sr_pd_raster) <- "GAM_residuals" #change name of the layer for clarity
plot(gam_residuals_sr_pd_raster, main="Residuals PD~TD") #checking

name_group<-"gam_residuals_SR_PD_10km_s.tif" #gam_residuals_SR_PD_1km_T10
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

name_group<-"gam_residuals_SR_FD_10km_s.tif" #"gam_residuals_SR_FD_1km_T10"
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

name_group<-"gam_residuals_PD_FD_10km_s.tif" #gam_residuals_PD_FD_1km_T10
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

name_group<-"gam_residuals_PD_FDdis_10km_s.tif" #"gam_residuals_PD_FDdis_1km_T10"
writeRaster(gam_residuals_pd_fddis_raster,name_group, filetype = "GTiff", overwrite=TRUE)


###################################### Other combinations #######################################################################
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

# Charger les rasters de base
pd_raster <- rast("hylids_PD10k_V3.tif")
sr_raster <- rast("richness_hylids_complete_0_083_V3.tif")
sr_raster_res <- resample(sr_raster, pd_raster, method = "near") #same dimension
plot(sr_raster)
plot(sr_raster_res)

#Function to create the rasters with the residuals
create_residual_raster <- function(ref_raster, df_resid, var_name, filename, comb_name) {
  res_raster <- ref_raster[[1]]
  values(res_raster) <- NA
  res_raster[df_resid$grid] <- df_resid$GAM_residuals
  names(res_raster) <- "GAM_residuals"
  plot(res_raster, main = paste("Residuals", comb_name, var_name))
  writeRaster(res_raster, filename, filetype = "GTiff", overwrite = TRUE)
}


# Fonction pour traiter les combinaisons
treat_comb <- function(x) {
  cat("Processing combination", x, "\n")
  
  setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")
  #Open rasters
  fdRic_raster <- rast(paste0("hylid_functional_richness10k_cont_V3_comb", x, ".tif"))
  fdDis_raster <- rast(paste0("hylid_functional_dispersion10k_cont_V3_comb", x, ".tif"))
  
  #Stack rasters
  stack_SR_FD <- c(fdRic_raster, sr_raster_res)
  stack_PD_FD <- c(fdRic_raster, pd_raster)
  stack_PD_FDdis <- c(fdDis_raster, pd_raster)
  
  ##Convert into dataframes
  #SR and FDric
  srfd_df <- as.data.frame(stack_SR_FD)
  names(srfd_df)[c(1,2)] <- c("FdRich", "SR")
  srfd_df$grid <- 1:nrow(srfd_df) #1:length ? 
  srfd_df_clean <- na.omit(srfd_df)
  
  #PD and FDric
  pdfd_df <- as.data.frame(stack_PD_FD)
  names(pdfd_df)[c(1,2)] <- c("FdRich", "PD")
  pdfd_df$grid <- 1:nrow(pdfd_df)
  pdfd_df_clean <- na.omit(pdfd_df)
  
  #PD and FDisp
  pdfddis_df <- as.data.frame(stack_PD_FDdis)
  names(pdfddis_df)[c(1,2)] <- c("FdDisp", "PD")
  pdfddis_df$grid <- 1:nrow(pdfddis_df)
  pdfddis_df_clean <- na.omit(pdfddis_df)
  
  #Generate Generalized additive models GAM 
  gam_model_sr_fd <- gam(FdRich ~ s(SR), data = srfd_df_clean) #te
  gam_model_pd_fd <- gam(FdRich ~ s(PD), data = pdfd_df_clean)
  gam_model_pd_fddis <- gam(FdDisp ~ s(PD), data = pdfddis_df_clean)
  
  #Compute and save residuals from gam models
  gam_residuals_sr_fd <- resid(gam_model_sr_fd)
  gam_residuals_pd_fd <- resid(gam_model_pd_fd)
  gam_residuals_pd_fddis <- resid(gam_model_pd_fddis)
  
  #New datasets
  srfd_df_resid<-srfd_df #with NA
  srfd_df_resid2<-srfd_df_resid
  srfd_df_resid2$GAM_residuals<-rep(NA,length(srfd_df[,1]))
  srfd_df_resid2$GAM_residuals[srfd_df_clean$grid]<-gam_residuals_sr_fd
  srfd_df_resid2$grid<-rownames(srfd_df_resid2) #prepare for the maps
  srfd_df_resid2$grid<- as.integer(srfd_df_resid2$grid)
  
  pdfd_df_resid<-pdfd_df 
  pdfd_df_resid2<-pdfd_df_resid
  pdfd_df_resid2$GAM_residuals<-rep(NA,length(pdfd_df[,1]))
  pdfd_df_resid2$GAM_residuals[pdfd_df_clean$grid]<-gam_residuals_pd_fd
  pdfd_df_resid2$grid<-rownames(pdfd_df_resid2)
  pdfd_df_resid2$grid<- as.integer(pdfd_df_resid2$grid)
  
  
  pdfddis_df_resid<-pdfddis_df  
  pdfddis_df_resid2<-pdfddis_df_resid
  pdfddis_df_resid2$GAM_residuals<-rep(NA,length(pdfddis_df[,1]))
  pdfddis_df_resid2$GAM_residuals[pdfddis_df_clean$grid]<-gam_residuals_pd_fddis
  pdfddis_df_resid2$grid<-rownames(pdfddis_df_resid2)
  pdfddis_df_resid2$grid<- as.integer(pdfddis_df_resid2$grid)
  
  # Save residuals as rasters
  setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/residuals/with_s")
  create_residual_raster(pd_raster, srfd_df_resid2, "FDRich~TD", paste0("gam_residuals_SR_FD_10km_s_comb", x, ".tif"), x)
  create_residual_raster(pd_raster, pdfd_df_resid2, "FDRich~PD", paste0("gam_residuals_PD_FD_10km_s_comb", x, ".tif"), x)
  create_residual_raster(pd_raster, pdfddis_df_resid2, "FDisp~PD", paste0("gam_residuals_PD_FDdis_s_10km_comb", x, ".tif"), x)
}

 
# Boucle sur les combinaisons 2 à 10
for (i in 2:10) {
  treat_comb(i)
}



