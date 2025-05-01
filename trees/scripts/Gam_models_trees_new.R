####Gam models for global trees, making the plots
#Modified by Tali Guez 2025
#Trees #new version

##First get coordinates from raster
library(terra)
library(ggplot2)
library(nlme)
library(mgcv)

rm(list = ls())
gc()

##Read in the maps
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps")

sr<-rast("eco_based_200/SR_eco_based_200.grd")
pd<-rast("eco_based_200/PD_eco_based_200.grd")
MPD<-rast("eco_based_200/MPD_eco_based_200.grd")

##### FDR ~ SR #####
###### Eco based - fdr ~ sr ######### 
fdr<-rast("eco_based_200/fdr_eco_based_200.grd")
raoq<-rast("eco_based_200/RaoQ_eco_based_200.grd")

stackRas<-c(raoq,sr,MPD,fdr,pd)
names(stackRas)<-c("raoq","sr","MPD","fdr","pd")

stackRas<-as.data.frame(stackRas)
#stackRas<-subset(stackRas,sr>8)
stackRas<-na.omit(stackRas)

#correcting for the intercept
#stackRas <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas)

#fdr (FDR)
m1_eco_based_gam <- gam(fdr ~ s(sr), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam)

m2_eco_based_gam <- gam(fdr ~ te(sr), family = "gaussian", data = stackRas)
m3_eco_based_gam <- gam(fdr ~ ti(sr), family = "gaussian", data = stackRas)
m4_eco_based_gam <- gam(fdr ~ t2(sr), family = "gaussian", data = stackRas)
m5_eco_based_gam <- gam(fdr ~ sr, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam,m2_eco_based_gam, m3_eco_based_gam,m4_eco_based_gam,m5_eco_based_gam) ##BEST=fdr~s(sr)

model_p<-predict(m1_eco_based_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data <- data.frame(
  sr = stackRas$sr,
  fit = model_p$fit,
  se.fit = model_p$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data
#pred_data <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data)

###### PCA - fdr ~ sr ###### 
fdr_PCA<-rast("PCA_selection/fdr_pca.grd")
raoq_PCA<-rast("PCA_selection/RaoQ_pca.grd")

stackRas_PCA<-c(raoq_PCA,sr,MPD,fdr_PCA,pd)
names(stackRas_PCA)<-c("raoq","sr","MPD","fdr","pd")

stackRas_PCA<-as.data.frame(stackRas_PCA)
#stackRas_PCA<-subset(stackRas_PCA,sr>8)
stackRas_PCA<-na.omit(stackRas_PCA)

#correcting for the intercept
#stackRas_PCA <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_PCA)

#For fdr
#Choose best model
m1_pca_gam <- gam(fdr ~ s(sr),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam)

m2_pca_gam <- gam(fdr ~ te(sr), family = "gaussian", data = stackRas_PCA)
m3_pca_gam <- gam(fdr ~ ti(sr), family = "gaussian", data = stackRas_PCA)
m4_pca_gam <- gam(fdr ~ t2(sr), family = "gaussian", data = stackRas_PCA)
m5_pca_gam <- gam(fdr ~ sr, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam,m2_pca_gam, m3_pca_gam,m4_pca_gam,m5_pca_gam) ##BEST=fdr~s(sr)

model_p_pca<-predict(m1_pca_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA <- data.frame(
  sr = stackRas_PCA$sr,
  fit = model_p_pca$fit,
  se.fit = model_p_pca$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data
#pred_data_PCA <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data_PCA)

###### Random within clusters - fdr ~ sr ###### 
fdr_within_clusters<-rast("random_clusters/fdr_clusters.grd")
raoq_within_clusters<-rast("random_clusters/RaoQ_clusters.grd")

stackRas_clusters<-c(raoq_within_clusters,sr,MPD,fdr_within_clusters,pd)
names(stackRas_clusters)<-c("raoq","sr","MPD","fdr","pd")

stackRas_clusters<-as.data.frame(stackRas_clusters)
#stackRas_clusters<-subset(stackRas_clusters,sr>8)
stackRas_clusters<-na.omit(stackRas_clusters)

#correcting for the intercept
#stackRas_clusters <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_clusters)

#Choose best model
m1_cluster_gam <- gam(fdr ~ s(sr),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam)

m2_cluster_gam <- gam(fdr ~ te(sr), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam <- gam(fdr ~ ti(sr), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam <- gam(fdr ~ t2(sr), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam <- gam(fdr ~ sr, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam,m2_cluster_gam, m3_cluster_gam,m4_cluster_gam,m5_cluster_gam) ##BEST=fdr~s(sr)

model_p_cluster<-predict(m1_cluster_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters <- data.frame(
  sr = stackRas_clusters$sr,
  fit = model_p_cluster$fit,
  se.fit = model_p_cluster$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data_clusters
#pred_data_clusters <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data_clusters)

###### Random X - fdr ~ sr #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  fdr_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  raoq_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random<-c(raoq_random, sr, MPD, fdr_random, pd)
  names(stackRas_random)<- c("raoq","sr","MPD","fdr","pd")
  
  stackRas_random <- as.data.frame(stackRas_random)
  stackRas_random <- subset(stackRas_random, sr>8)
  stackRas_random <- na.omit(stackRas_random)
  
  #Correcting for the intercept
  #stackRas_random <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_random)
  
  #calculate all models
  models <- list(
    s_model = gam(fdr ~ s(sr), family = "gaussian", data = stackRas_random),
    te_model = gam(fdr ~ te(sr), family = "gaussian", data = stackRas_random),
    ti_model = gam(fdr ~ ti(sr), family = "gaussian", data = stackRas_random),
    t2_model = gam(fdr ~ t2(sr), family = "gaussian", data = stackRas_random),
    linear_model = gam(fdr ~ sr, family = "gaussian", data = stackRas_random)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("Random_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values
  )
  
  cat("Random", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}

#Create lists to store models and data
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

#Loop for 10 randoms
for (i in 1:10) {
  #Load data
  fdr_random <- rast(paste0("random", i, "/fdr_random", i, ".grd"))
  raoq_random <- rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random <- c(raoq_random, sr, MPD, fdr_random, pd)
  names(stackRas_random) <- c("raoq", "sr", "MPD", "fdr", "pd")
  
  stackRas_random <-as.data.frame(stackRas_random)
  #stackRas_random<-subset(stackRas_random, sr>8)
  stackRas_random <-na.omit(stackRas_random)
  
  #save 
  stackRas_list[[paste0("stackRas_random", i)]] <- stackRas_random
  
  #Gam models
  model_gam <- gam(fdr ~ s(sr), family = "gaussian", data = stackRas_random)
  model_list[[paste0("model_gam_random", i)]] <- model_gam
  
  model_p_random <- predict(model_gam, se.fit = TRUE)
  
  #Create dataframe with predictions
  pred_data_random <- data.frame(
    sr = stackRas_random$sr,
    fit = model_p_random$fit,
    se.fit = model_p_random$se.fit
  )
  
  #Save predictions
  pred_data_list[[paste0("pred_data_random", i)]] <- pred_data_random
}

#stackRas_list$stackRas_random1 # stocke les datasets pour chaque random
#model_list$model_gam_random1 #stocke les modèles gam
#pred_data_list$pred_data_random1 #stocke les prédictions

#Correcting afterwards - Add manually (0,0) to pred_data_list
#for (i in 1:10) {
#  pred_data_list[[paste0("pred_data_random_raoq_fd", i)]] <- rbind(
#    data.frame(sr = 0, fit = 0, se.fit = 0),
#    pred_data_list[[paste0("pred_data_random_raoq_fd", i)]]
#  )
#}

###### Making the plots - fdr ~ sr #####

#1.all - no subset, (not (0;0) at the origin/not (0;0) at the origin)
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=sr, y=fdr, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random1, aes(x=sr, y=fdr, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random1, aes(x=sr, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random2, aes(x=sr, y=fdr, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random2, aes(x=sr, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random3, aes(x=sr, y=fdr, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random3, aes(x=sr, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random4, aes(x=sr, y=fdr, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random4, aes(x=sr, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random5, aes(x=sr, y=fdr, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random5, aes(x=sr, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random6, aes(x=sr, y=fdr, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random6, aes(x=sr, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random7, aes(x=sr, y=fdr, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random7, aes(x=sr, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random8, aes(x=sr, y=fdr, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random8, aes(x=sr, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random9, aes(x=sr, y=fdr, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random9, aes(x=sr, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random10, aes(x=sr, y=fdr, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random10, aes(x=sr, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=sr, y=fdr, color="Within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters, aes(x=sr, y=fit, color="Within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=sr, y=fdr, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Richness") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR") + 
  ggtitle("Trees - Functional Richness ~ SR (0:0)") #(0:0)

g

#2. subset, (not (0;0) at the origin/not (0;0) at the origin)
g <- ggplot() +
  # eco_based
  geom_point(data=sample_frac(stackRas, 0.05), aes(x=sr, y=fdr, color="Ecologically based"), alpha=0.08) +
  geom_smooth(data=pred_data, aes(x=sr, y=fit, color="Ecologically based")) +
  
  # random 1
  geom_point(data=sample_frac(stackRas_list$stackRas_random1, 0.05), aes(x=sr, y=fdr, color="Random 1"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random1, aes(x=sr, y=fit, color="Random 1")) +
  
  # random 2
  geom_point(data=sample_frac(stackRas_list$stackRas_random2, 0.05), aes(x=sr, y=fdr, color="Random 2"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random2, aes(x=sr, y=fit, color="Random 2")) +
  
  # random 3
  geom_point(data=sample_frac(stackRas_list$stackRas_random3, 0.05), aes(x=sr, y=fdr, color="Random 3"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random3, aes(x=sr, y=fit, color="Random 3")) +
  
  # random 4
  geom_point(data=sample_frac(stackRas_list$stackRas_random4, 0.05), aes(x=sr, y=fdr, color="Random 4"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random4, aes(x=sr, y=fit, color="Random 4")) +
  
  # random 5
  geom_point(data=sample_frac(stackRas_list$stackRas_random5, 0.05), aes(x=sr, y=fdr, color="Random 5"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random5, aes(x=sr, y=fit, color="Random 5")) +
  
  # random 6
  geom_point(data=sample_frac(stackRas_list$stackRas_random6, 0.05), aes(x=sr, y=fdr, color="Random 6"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random6, aes(x=sr, y=fit, color="Random 6")) +
  
  # random 7
  geom_point(data=sample_frac(stackRas_list$stackRas_random7, 0.05), aes(x=sr, y=fdr, color="Random 7"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random7, aes(x=sr, y=fit, color="Random 7")) +
  
  # random 8
  geom_point(data=sample_frac(stackRas_list$stackRas_random8, 0.05), aes(x=sr, y=fdr, color="Random 8"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random8, aes(x=sr, y=fit, color="Random 8")) +
  
  # random 9
  geom_point(data=sample_frac(stackRas_list$stackRas_random9, 0.05), aes(x=sr, y=fdr, color="Random 9"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random9, aes(x=sr, y=fit, color="Random 9")) +
  
  # random 10
  geom_point(data=sample_frac(stackRas_list$stackRas_random10, 0.05), aes(x=sr, y=fdr, color="Random 10"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random10, aes(x=sr, y=fit, color="Random 10")) +
  
  # within_clusters
  geom_point(data=sample_frac(stackRas_clusters, 0.05), aes(x=sr, y=fdr, color="Within clusters"), alpha=0.08) +
  geom_smooth(data=pred_data_clusters, aes(x=sr, y=fit, color="Within clusters")) +
  
  # PCA
  geom_point(data=sample_frac(stackRas_PCA, 0.05), aes(x=sr, y=fdr, color="PCA"), alpha=0.08) +
  geom_smooth(data=pred_data_PCA , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Richness") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Dataset") +
  coord_cartesian(ylim=c(0,4),xlim=c(0,NA),expand=FALSE) + #expand=FALSE, xlim=c(0,NA)), ylim=c(1,NA)
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=15.3),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.15, "cm")) +
  labs(x="Taxonomic diversity SR") #+ 
  ggtitle("Trees - Functional Richness~SR") #(0:0)
g

### RaoQ ~ SR #####
##First get coordinates from raster
library(terra)
library(ggplot2)
library(nlme)
library(mgcv)

rm(list = ls())

##Read in the maps
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps")

sr<-rast("eco_based_200/SR_eco_based_200.grd")
pd<-rast("eco_based_200/PD_eco_based_200.grd")
MPD<-rast("eco_based_200/MPD_eco_based_200.grd")

###### Eco based - raoq ~ sr ######### 
fdr<-rast("eco_based_200/fdr_eco_based_200.grd")
raoq<-rast("eco_based_200/RaoQ_eco_based_200.grd")

stackRas<-c(raoq,sr,MPD,fdr,pd)
names(stackRas)<-c("raoq","sr","MPD","fdr","pd")

stackRas<-as.data.frame(stackRas)
#stackRas<-subset(stackRas,sr>8)
stackRas<-na.omit(stackRas)

#correcting for the intercept
#stackRas <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas)

#RaoQ ~ SR
m1_eco_based_gam_raoq <- gam(raoq ~ s(sr), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam_raoq)

m2_eco_based_gam_raoq  <- gam(raoq ~ te(sr), family = "gaussian", data = stackRas)
m3_eco_based_gam_raoq  <- gam(raoq ~ ti(sr), family = "gaussian", data = stackRas)
m4_eco_based_gam_raoq  <- gam(raoq ~ t2(sr), family = "gaussian", data = stackRas)
m5_eco_based_gam_raoq  <- gam(raoq ~ sr, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam_raoq,m2_eco_based_gam_raoq , m3_eco_based_gam_raoq ,m4_eco_based_gam_raoq ,m5_eco_based_gam_raoq ) ##BEST=raoq~s(sr)
#s best

model_p_raoq<-predict(m1_eco_based_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_raoq <- data.frame(
  sr = stackRas$sr,
  fit = model_p_raoq$fit,
  se.fit = model_p_raoq$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data
#pred_data_raoq <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data_raoq)

###### PCA - raoq ~ sr ###### 
fdr_PCA<-rast("PCA_selection/fdr_pca.grd")
raoq_PCA<-rast("PCA_selection/RaoQ_pca.grd")

stackRas_PCA<-c(raoq_PCA,sr,MPD,fdr_PCA,pd)
names(stackRas_PCA)<-c("raoq","sr","MPD","fdr","pd")

stackRas_PCA<-as.data.frame(stackRas_PCA)
#stackRas_PCA<-subset(stackRas_PCA,sr>8)
stackRas_PCA<-na.omit(stackRas_PCA)

#correcting for the intercept
#stackRas_PCA <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_PCA)

#Choose best model
m1_pca_gam_raoq <- gam(raoq ~ s(sr),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam_raoq )

m2_pca_gam_raoq  <- gam(raoq ~ te(sr), family = "gaussian", data = stackRas_PCA)
m3_pca_gam_raoq  <- gam(raoq ~ ti(sr), family = "gaussian", data = stackRas_PCA)
m4_pca_gam_raoq  <- gam(raoq ~ t2(sr), family = "gaussian", data = stackRas_PCA)
m5_pca_gam_raoq  <- gam(raoq ~ sr, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam_raoq ,m2_pca_gam_raoq , m3_pca_gam_raoq ,m4_pca_gam_raoq ,m5_pca_gam_raoq ) ##BEST=raoq~s(sr)

model_p_pca_raoq<-predict(m1_pca_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA_raoq <- data.frame(
  sr = stackRas_PCA$sr,
  fit = model_p_pca_raoq$fit,
  se.fit = model_p_pca_raoq$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data
#pred_data_PCA_raoq <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data_PCA_raoq)

###### Random within clusters - raoq ~ sr ###### 
fdr_within_clusters<-rast("random_clusters/fdr_clusters.grd")
raoq_within_clusters<-rast("random_clusters/RaoQ_clusters.grd")

stackRas_clusters<-c(raoq_within_clusters,sr,MPD,fdr_within_clusters,pd)
names(stackRas_clusters)<-c("raoq","sr","MPD","fdr","pd")

stackRas_clusters<-as.data.frame(stackRas_clusters)
#stackRas_clusters<-subset(stackRas_clusters,sr>8)
stackRas_clusters<-na.omit(stackRas_clusters)

#correcting for the intercept
#stackRas_clusters <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_clusters)

#Choose best model
m1_cluster_gam_raoq <- gam(raoq ~ s(sr),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam_raoq)

m2_cluster_gam_raoq <- gam(raoq ~ te(sr), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam_raoq <- gam(raoq ~ ti(sr), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam_raoq <- gam(raoq ~ t2(sr), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam_raoq <- gam(raoq ~ sr, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam_raoq,m2_cluster_gam_raoq, m3_cluster_gam_raoq,m4_cluster_gam_raoq,m5_cluster_gam_raoq) ##BEST=raoq~s(sr)

model_p_cluster_raoq<-predict(m1_cluster_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters_raoq <- data.frame(
  sr = stackRas_clusters$sr,
  fit = model_p_cluster_raoq$fit,
  se.fit = model_p_cluster_raoq$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data
#pred_data_clusters_raoq <- rbind(data.frame(sr = 0, fit = 0, se.fit = 0),pred_data_clusters_raoq)

###### Random X - raoq ~ sr #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  fdr_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  raoq_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq<-c(raoq_random, sr, MPD, fdr_random, pd)
  names(stackRas_random_raoq)<- c("raoq","sr","MPD","fdr","pd")
  
  stackRas_random_raoq <- as.data.frame(stackRas_random_raoq)
  #stackRas_random_raoq <- subset(stackRas_random_raoq, sr > 8)
  stackRas_random_raoq <- na.omit(stackRas_random_raoq)
  
  #correcting for the intercept
  #stackRas_random_raoq <- rbind(data.frame(raoq = 0, sr = 0, MPD = 0, fdr = 0, pd = 0), stackRas_random_raoq)
  
  
  #calculate all models
  models <- list(
    s_model = gam(raoq ~ s(sr), family = "gaussian", data = stackRas_random_raoq),
    te_model = gam(raoq ~ te(sr), family = "gaussian", data = stackRas_random_raoq),
    ti_model = gam(raoq ~ ti(sr), family = "gaussian", data = stackRas_random_raoq),
    t2_model = gam(raoq ~ t2(sr), family = "gaussian", data = stackRas_random_raoq),
    linear_model = gam(raoq ~ sr, family = "gaussian", data = stackRas_random_raoq)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("Random_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values,
    data = stackRas_random_raoq # store the data 
  )
  
  cat("Random", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}
#s best  

# Create lists to store results
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

for (i in 1:10) {
  key <- paste0("Random_", i)
  
  # Retrieve data and model from the best_models list
  data_i <- best_models[[key]]$data
  model_i <- best_models[[key]]$best_model
  
  # Save to list (optional, if needed later)
  stackRas_list[[paste0("stackRas_random_raoq", i)]] <- data_i
  model_list[[paste0("model_gam_random_raoq", i)]] <- model_i
  
  # Predict
  pred_i <- predict(model_i, se.fit = TRUE)
  
  # Create prediction data frame
  pred_data_i <- data.frame(
    sr = data_i$sr,
    fit = pred_i$fit,
    se.fit = pred_i$se.fit
  )
  
  pred_data_list[[paste0("pred_data_random_raoq", i)]] <- pred_data_i
}

#stackRas_list$stackRas_random_raoq1 # stocke les datasets pour chaque random
#model_list$model_gam_random_raoq1 #stocke les modèles gam
#pred_data_list$pred_data_randomraoq1 #stocke les prédictions

#Correcting afterwards - Add manually (0,0) to pred_data
#for (i in 1:10) {
#  pred_data_list[[paste0("pred_data_random_raoq_fd", i)]] <- rbind(
#    data.frame(sr = 0, fit = 0, se.fit = 0),
#    pred_data_list[[paste0("pred_data_random_raoq_fd", i)]]
#  )
#}

###### Making the plots - raoq ~ sr #####

#1.all - no subset, (not (0;0) at the origin/(0;0) at the origin)
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=sr, y=raoq, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data_raoq, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random_raoq1, aes(x=sr, y=raoq, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq1, aes(x=sr, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random_raoq2, aes(x=sr, y=raoq, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq2, aes(x=sr, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random_raoq3, aes(x=sr, y=raoq, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq3, aes(x=sr, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random_raoq4, aes(x=sr, y=raoq, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq4, aes(x=sr, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random_raoq5, aes(x=sr, y=raoq, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq5, aes(x=sr, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random_raoq6, aes(x=sr, y=raoq, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq6, aes(x=sr, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random_raoq7, aes(x=sr, y=raoq, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq7, aes(x=sr, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random_raoq8, aes(x=sr, y=raoq, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq8, aes(x=sr, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random_raoq9, aes(x=sr, y=raoq, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq9, aes(x=sr, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random_raoq10, aes(x=sr, y=raoq, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq10, aes(x=sr, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=sr, y=raoq, color="Within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters_raoq, aes(x=sr, y=fit, color="Within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=sr, y=raoq, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA_raoq , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR")+
  ggtitle("Trees -  RaoQ~SR (0:0)") #(0:0)

g

#2. subset, (not (0;0) at the origin/(0;0) at the origin)
g <- ggplot() +
  # eco_based
  geom_point(data=sample_frac(stackRas, 0.05), aes(x=sr, y=raoq, color="Ecologically based"), alpha=0.08) +
  geom_smooth(data=pred_data_raoq, aes(x=sr, y=fit, color="Ecologically based")) +
  
  # random 1
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq1, 0.05), aes(x=sr, y=raoq, color="Random 1"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq1, aes(x=sr, y=fit, color="Random 1")) +
  
  # random 2
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq2, 0.05), aes(x=sr, y=raoq, color="Random 2"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq2, aes(x=sr, y=fit, color="Random 2")) +
  
  # random 3
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq3, 0.05), aes(x=sr, y=raoq, color="Random 3"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq3, aes(x=sr, y=fit, color="Random 3")) +
  
  # random 4
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq4, 0.05), aes(x=sr, y=raoq, color="Random 4"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq4, aes(x=sr, y=fit, color="Random 4")) +
  
  # random 5
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq5, 0.05), aes(x=sr, y=raoq, color="Random 5"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq5, aes(x=sr, y=fit, color="Random 5")) +
  
  # random 6
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq6, 0.05), aes(x=sr, y=raoq, color="Random 6"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq6, aes(x=sr, y=fit, color="Random 6")) +
  
  # random 7
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq7, 0.05), aes(x=sr, y=raoq, color="Random 7"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq7, aes(x=sr, y=fit, color="Random 7")) +
  
  # random 8
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq8, 0.05), aes(x=sr, y=raoq, color="Random 8"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq8, aes(x=sr, y=fit, color="Random 8")) +
  
  # random 9
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq9, 0.05), aes(x=sr, y=raoq, color="Random 9"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq9, aes(x=sr, y=fit, color="Random 9")) +
  
  # random 10
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq10, 0.05), aes(x=sr, y=raoq, color="Random 10"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq10, aes(x=sr, y=fit, color="Random 10")) +
  
  # within_clusters
  geom_point(data=sample_frac(stackRas_clusters, 0.05), aes(x=sr, y=raoq, color="Within clusters"), alpha=0.08) +
  geom_smooth(data=pred_data_clusters_raoq, aes(x=sr, y=fit, color="Within clusters")) +
  
  # PCA
  geom_point(data=sample_frac(stackRas_PCA, 0.05), aes(x=sr, y=raoq, color="PCA"), alpha=0.08) +
  geom_smooth(data=pred_data_PCA_raoq , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position="bottom") +
  labs(x="Taxonomic diversity SR") +
  ggtitle("Trees - RaoQ ~ SR") #(0:0)
g

##### RaoQ ~ FDR #####
##First get coordinates from raster
library(terra)
library(ggplot2)
library(nlme)
library(mgcv)

rm(list = ls())

##Read in the maps
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps")

sr<-rast("eco_based_200/SR_eco_based_200.grd")
pd<-rast("eco_based_200/PD_eco_based_200.grd")
MPD<-rast("eco_based_200/MPD_eco_based_200.grd")

###### Eco based - raoq ~ fdr ######### 
fdr<-rast("eco_based_200/fdr_eco_based_200.grd")
raoq<-rast("eco_based_200/RaoQ_eco_based_200.grd")

stackRas<-c(raoq,sr,MPD,fdr,pd)
names(stackRas)<-c("raoq","sr","MPD","fdr","pd")

stackRas<-as.data.frame(stackRas)
#stackRas<-subset(stackRas,sr>8)
stackRas<-na.omit(stackRas)

#Without correcting
m1_eco_based_gam_raoq_fd <- gam(raoq ~ s(fdr), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam_raoq_fd)

m2_eco_based_gam_raoq_fd<-gam(raoq ~ te(fdr), family = "gaussian", data = stackRas)
m3_eco_based_gam_raoq_fd<-gam(raoq ~ ti(fdr), family = "gaussian", data = stackRas)
m4_eco_based_gam_raoq_fd<-gam(raoq ~ t2(fdr), family = "gaussian", data = stackRas)
m5_eco_based_gam_raoq_fd<-gam(raoq ~ fdr, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam_raoq_fd,m2_eco_based_gam_raoq_fd , m3_eco_based_gam_raoq_fd ,m4_eco_based_gam_raoq_fd ,m5_eco_based_gam_raoq_fd ) ##BEST=raoq~s(fdr)
#best 1
model_p_raoq_fd<-predict(m1_eco_based_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_raoq_fd <- data.frame(
  fdr = stackRas$fdr,
  fit = model_p_raoq_fd$fit,
  se.fit = model_p_raoq_fd$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data_raoq_fd
#pred_data_raoq_fd <- rbind(data.frame(fdr = 0, fit = 0, se.fit = 0),pred_data_raoq_fd)

###### PCA - raoq ~ fdr ###### 
fdr_PCA<-rast("PCA_selection/fdr_pca.grd")
raoq_PCA<-rast("PCA_selection/RaoQ_pca.grd")

stackRas_PCA<-c(raoq_PCA,sr,MPD,fdr_PCA,pd)
names(stackRas_PCA)<-c("raoq","sr","MPD","fdr","pd")

stackRas_PCA<-as.data.frame(stackRas_PCA)
#stackRas_PCA<-subset(stackRas_PCA,sr>8)
stackRas_PCA<-na.omit(stackRas_PCA)

#Choose best model
m1_pca_gam_raoq_fd <- gam(raoq ~ s(fdr),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam_raoq_fd)

m2_pca_gam_raoq_fd <- gam(raoq ~ te(fdr), family = "gaussian", data = stackRas_PCA)
m3_pca_gam_raoq_fd  <- gam(raoq ~ ti(fdr), family = "gaussian", data = stackRas_PCA)
m4_pca_gam_raoq_fd  <- gam(raoq ~ t2(fdr), family = "gaussian", data = stackRas_PCA)
m5_pca_gam_raoq_fd  <- gam(raoq ~ fdr, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam_raoq_fd ,m2_pca_gam_raoq_fd , m3_pca_gam_raoq_fd ,m4_pca_gam_raoq_fd ,m5_pca_gam_raoq_fd) ##BEST=raoq~s(fdr)

model_p_pca_raoq_fd<-predict(m1_pca_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA_raoq_fd <- data.frame(
  fdr = stackRas_PCA$fdr,
  fit = model_p_pca_raoq_fd$fit,
  se.fit = model_p_pca_raoq_fd$se.fit
)

#Correcting afterwards - Add manually (0,0) to pred_data_PCA_raoq_fd
#pred_data_PCA_raoq_fd <- rbind(data.frame(fdr = 0, fit = 0, se.fit = 0), pred_data_PCA_raoq_fd)

###### Random within clusters - raoq ~ fdr ###### 
fdr_within_clusters<-rast("random_clusters/fdr_clusters.grd")
raoq_within_clusters<-rast("random_clusters/RaoQ_clusters.grd")

stackRas_clusters<-c(raoq_within_clusters,sr,MPD,fdr_within_clusters,pd)
names(stackRas_clusters)<-c("raoq","sr","MPD","fdr","pd")

stackRas_clusters<-as.data.frame(stackRas_clusters)
#stackRas_clusters<-subset(stackRas_clusters,sr>8)
stackRas_clusters<-na.omit(stackRas_clusters)

#Choose best model
m1_cluster_gam_raoq_fd <- gam(raoq ~ s(fdr),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam_raoq_fd)

m2_cluster_gam_raoq_fd <- gam(raoq ~ te(fdr), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam_raoq_fd <- gam(raoq ~ ti(fdr), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam_raoq_fd <- gam(raoq ~ t2(fdr), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam_raoq_fd <- gam(raoq ~ fdr, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam_raoq_fd,m2_cluster_gam_raoq_fd, m3_cluster_gam_raoq_fd,m4_cluster_gam_raoq_fd,m5_cluster_gam_raoq_fd) ##BEST=raoq~s(fdr)
#m1 best
model_p_cluster_raoq_fd<-predict(m1_cluster_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters_raoq_fd <- data.frame(
  fdr = stackRas_clusters$fdr,
  fit = model_p_cluster_raoq_fd$fit,
  se.fit = model_p_cluster_raoq_fd$se.fit
)

#Correcting afterwards -  Add manually (0,0) to pred_data_clusters_raoq_fd
#pred_data_clusters_raoq_fd <- rbind(data.frame(fdr = 0, fit = 0, se.fit = 0), pred_data_clusters_raoq_fd)

###### Random X -raoq ~ fdr #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  fdr_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  raoq_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq_fd<-c(raoq_random, sr, MPD, fdr_random, pd)
  names(stackRas_random_raoq_fd)<- c("raoq","sr","MPD","fdr","pd")
  
  stackRas_random_raoq_fd<-as.data.frame(stackRas_random_raoq_fd)
  #stackRas_random_raoq_fd<-subset(stackRas_random_raoq_fd, sr>8)
  stackRas_random_raoq_fd<-na.omit(stackRas_random_raoq_fd)
  
  #calculate all models
  models <- list(
    s_model = gam(raoq ~ s(fdr), family = "gaussian", data = stackRas_random_raoq_fd),
    te_model = gam(raoq ~ te(fdr), family = "gaussian", data = stackRas_random_raoq_fd),
    ti_model = gam(raoq ~ ti(fdr), family = "gaussian", data = stackRas_random_raoq_fd),
    t2_model = gam(raoq ~ t2(fdr), family = "gaussian", data = stackRas_random_raoq_fd),
    linear_model = gam(raoq ~ fdr, family = "gaussian", data = stackRas_random_raoq_fd)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("Random_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values,
    data = stackRas_random_raoq_fd # store the data 
  )
  
  #Print results
  cat("Random", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}
#s model= the best 

# Create lists to store results
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

for (i in 1:10) {
  key <- paste0("Random_", i)
  
  # Retrieve data and model from the best_models list
  data_i <- best_models[[key]]$data
  model_i <- best_models[[key]]$best_model
  
  # Save to list (optional, if needed later)
  stackRas_list[[paste0("stackRas_random_raoq_fd", i)]] <- data_i
  model_list[[paste0("model_gam_random_raoq_fd", i)]] <- model_i
  
  # Predict
  pred_i <- predict(model_i, se.fit = TRUE)
  
  # Create prediction data frame
  pred_data_i <- data.frame(
    fdr = data_i$fdr,
    fit = pred_i$fit,
    se.fit = pred_i$se.fit
  )
  
  pred_data_list[[paste0("pred_data_random_raoq_fd", i)]] <- pred_data_i
}


#Correcting afterwards - Add manually (0,0) to pred_data_list
#for (i in 1:10) {
 # pred_data_list[[paste0("pred_data_random_raoq_fd", i)]] <- rbind(
  #  data.frame(fdr = 0, fit = 0, se.fit = 0),
   # pred_data_list[[paste0("pred_data_random_raoq_fd", i)]]
  #)
#}


###### Making the plots - raoq ~ fdr #####

#1.all - no subset, (not (0;0) at the origin/(0;0) at the origin)
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=fdr, y=raoq, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data_raoq_fd, aes(x=fdr, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random_raoq_fd1, aes(x=fdr, y=raoq, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd1, aes(x=fdr, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random_raoq_fd2, aes(x=fdr, y=raoq, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd2, aes(x=fdr, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random_raoq_fd3, aes(x=fdr, y=raoq, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd3, aes(x=fdr, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random_raoq_fd4, aes(x=fdr, y=raoq, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd4, aes(x=fdr, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random_raoq_fd5, aes(x=fdr, y=raoq, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd5, aes(x=fdr, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random_raoq_fd6, aes(x=fdr, y=raoq, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd6, aes(x=fdr, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random_raoq_fd7, aes(x=fdr, y=raoq, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd7, aes(x=fdr, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random_raoq_fd8, aes(x=fdr, y=raoq, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd8, aes(x=fdr, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random_raoq_fd9, aes(x=fdr, y=raoq, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd9, aes(x=fdr, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random_raoq_fd10, aes(x=fdr, y=raoq, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd10, aes(x=fdr, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=fdr, y=raoq, color="Within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters_raoq_fd, aes(x=fdr, y=fit, color="Within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=fdr, y=raoq, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA_raoq_fd , aes(x=fdr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Functional Richness (FDR)") +
  ggtitle("Trees - RaoQ ~ FDR (0:0)") #(0:0) 

g

#2. subset, (not (0;0) at the origin/(0;0) at the origin)
g <- ggplot() +
  # eco_based
  geom_point(data=sample_frac(stackRas, 0.05), aes(x=fdr, y=raoq, color="Ecologically based"), alpha=0.08) +
  geom_smooth(data=pred_data_raoq_fd, aes(x=fdr, y=fit, color="Ecologically based")) +
  
  # random 1
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd1, 0.05), aes(x=fdr, y=raoq, color="Random 1"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd1, aes(x=fdr, y=fit, color="Random 1")) +
  
  # random 2
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd2, 0.05), aes(x=fdr, y=raoq, color="Random 2"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd2, aes(x=fdr, y=fit, color="Random 2")) +
  
  # random 3
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd3, 0.05), aes(x=fdr, y=raoq, color="Random 3"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd3, aes(x=fdr, y=fit, color="Random 3")) +
  
  # random 4
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd4, 0.05), aes(x=fdr, y=raoq, color="Random 4"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd4, aes(x=fdr, y=fit, color="Random 4")) +
  
  # random 5
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd5, 0.05), aes(x=fdr, y=raoq, color="Random 5"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd5, aes(x=fdr, y=fit, color="Random 5")) +
  
  # random 6
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd6, 0.05), aes(x=fdr, y=raoq, color="Random 6"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd6, aes(x=fdr, y=fit, color="Random 6")) +
  
  # random 7
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd7, 0.05), aes(x=fdr, y=raoq, color="Random 7"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd7, aes(x=fdr, y=fit, color="Random 7")) +
  
  # random 8
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd8, 0.05), aes(x=fdr, y=raoq, color="Random 8"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd8, aes(x=fdr, y=fit, color="Random 8")) +
  
  # random 9
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd9, 0.05), aes(x=fdr, y=raoq, color="Random 9"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd9, aes(x=fdr, y=fit, color="Random 9")) +
  
  # random 10
  geom_point(data=sample_frac(stackRas_list$stackRas_random_raoq_fd10, 0.05), aes(x=fdr, y=raoq, color="Random 10"), alpha=0.08) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd10, aes(x=fdr, y=fit, color="Random 10")) +
  
  # within_clusters
  geom_point(data=sample_frac(stackRas_clusters, 0.05), aes(x=fdr, y=raoq, color="Within clusters"), alpha=0.08) +
  geom_smooth(data=pred_data_clusters_raoq_fd, aes(x=fdr, y=fit, color="Within clusters")) +
  
  # PCA
  geom_point(data=sample_frac(stackRas_PCA, 0.05), aes(x=fdr, y=raoq, color="PCA"), alpha=0.08) +
  geom_smooth(data=pred_data_PCA_raoq_fd, aes(x=fdr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c(
      "Ecologically based" = "darkgreen",
      "Random 1" = "red",
      "Random 2" = "green",
      "Random 3" = "blue",
      "Random 4" = "pink",
      "Random 5" = "brown",
      "Random 6" = "grey",
      "Random 7" = "black",
      "Random 8" = "royalblue1",
      "Random 9" = "orange",
      "Random 10" = "yellow",
      "PCA" = "purple",
      "Within clusters" = "aquamarine"
    ),
    name = "Dataset") +
  coord_cartesian(expand=FALSE, ylim=c(0,NA), xlim=c(0,4)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=15.3),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.15, "cm")) +
  labs(x = "Functional Richness (FDR)")# +
  #ggtitle("Trees - RaoQ~FDR") #(0:0)
g
