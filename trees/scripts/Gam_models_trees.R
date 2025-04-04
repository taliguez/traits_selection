####Gam models for trees
#Modified by Tali Guez 2025
#Trees

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

##### Loading data + create models #####
###### For eco_based ######### 
Frich<-rast("eco_based_200/fdr_eco_based_200.grd")
Fdiv<-rast("eco_based_200/RaoQ_eco_based_200.grd")

stackRas<-c(Fdiv,sr,MPD,Frich,pd)
names(stackRas)<-c("Fdiv","sr","MPD","Frich","pd")

stackRas<-as.data.frame(stackRas)
stackRas<-subset(stackRas,sr>3)
stackRas<-na.omit(stackRas)

#Frich (FDR)
m1_eco_based_gam <- gam(Frich ~ s(sr), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam)

m2_eco_based_gam <- gam(Frich ~ te(sr), family = "gaussian", data = stackRas)
m3_eco_based_gam <- gam(Frich ~ ti(sr), family = "gaussian", data = stackRas)
m4_eco_based_gam <- gam(Frich ~ t2(sr), family = "gaussian", data = stackRas)
m5_eco_based_gam <- gam(Frich ~ sr, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam,m2_eco_based_gam, m3_eco_based_gam,m4_eco_based_gam,m5_eco_based_gam) ##BEST=Frich~s(sr)

model_p<-predict(m1_eco_based_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data <- data.frame(
  sr = stackRas$sr,
  fit = model_p$fit,
  se.fit = model_p$se.fit
)

#RaoQ ~ SR
m1_eco_based_gam_raoq <- gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam_raoq)

m2_eco_based_gam_raoq  <- gam(Fdiv ~ te(sr), family = "gaussian", data = stackRas)
m3_eco_based_gam_raoq  <- gam(Fdiv ~ ti(sr), family = "gaussian", data = stackRas)
m4_eco_based_gam_raoq  <- gam(Fdiv ~ t2(sr), family = "gaussian", data = stackRas)
m5_eco_based_gam_raoq  <- gam(Fdiv ~ sr, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam_raoq,m2_eco_based_gam_raoq , m3_eco_based_gam_raoq ,m4_eco_based_gam_raoq ,m5_eco_based_gam_raoq ) ##BEST=Fdiv~s(sr)

model_p_raoq<-predict(m1_eco_based_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_raoq <- data.frame(
  sr = stackRas$sr,
  fit = model_p_raoq$fit,
  se.fit = model_p_raoq$se.fit
)

#RaoQ ~ FDR
m1_eco_based_gam_raoq_fd <- gam(Fdiv ~ s(Frich), family = "gaussian", data = stackRas)
summary(m1_eco_based_gam_raoq_fd)

m2_eco_based_gam_raoq_fd  <- gam(Fdiv ~ te(Frich), family = "gaussian", data = stackRas)
m3_eco_based_gam_raoq_fd <- gam(Fdiv ~ ti(Frich), family = "gaussian", data = stackRas)
m4_eco_based_gam_raoq_fd  <- gam(Fdiv ~ t2(Frich), family = "gaussian", data = stackRas)
m5_eco_based_gam_raoq_fd  <- gam(Fdiv ~ Frich, family = "gaussian", data = stackRas)
AIC(m1_eco_based_gam_raoq_fd,m2_eco_based_gam_raoq_fd , m3_eco_based_gam_raoq_fd ,m4_eco_based_gam_raoq_fd ,m5_eco_based_gam_raoq_fd ) ##BEST=Fdiv~s(Frich)

model_p_raoq_fd<-predict(m1_eco_based_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_raoq_fd <- data.frame(
  Frich = stackRas$Frich,
  fit = model_p_raoq_fd$fit,
  se.fit = model_p_raoq_fd$se.fit
)

###### PCA ###### 
Frich_PCA<-rast("PCA_selection/fdr_pca.grd")
Fdiv_PCA<-rast("PCA_selection/RaoQ_pca.grd")

stackRas_PCA<-c(Fdiv_PCA,sr,MPD,Frich_PCA,pd)
names(stackRas_PCA)<-c("Fdiv","sr","MPD","Frich","pd")

stackRas_PCA<-as.data.frame(stackRas_PCA)
stackRas_PCA<-subset(stackRas_PCA,sr>3)
stackRas_PCA<-na.omit(stackRas_PCA)

#For Frich
#Choose best model
m1_pca_gam <- gam(Frich ~ s(sr),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam)

m2_pca_gam <- gam(Frich ~ te(sr), family = "gaussian", data = stackRas_PCA)
m3_pca_gam <- gam(Frich ~ ti(sr), family = "gaussian", data = stackRas_PCA)
m4_pca_gam <- gam(Frich ~ t2(sr), family = "gaussian", data = stackRas_PCA)
m5_pca_gam <- gam(Frich ~ sr, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam,m2_pca_gam, m3_pca_gam,m4_pca_gam,m5_pca_gam) ##BEST=Frich~s(sr)

model_p_pca<-predict(m1_pca_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA <- data.frame(
  sr = stackRas_PCA$sr,
  fit = model_p_pca$fit,
  se.fit = model_p_pca$se.fit
)

#For RaoQ ~ SR
#Choose best model
m1_pca_gam_raoq <- gam(Fdiv ~ s(sr),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam_raoq )

m2_pca_gam_raoq  <- gam(Fdiv ~ te(sr), family = "gaussian", data = stackRas_PCA)
m3_pca_gam_raoq  <- gam(Fdiv ~ ti(sr), family = "gaussian", data = stackRas_PCA)
m4_pca_gam_raoq  <- gam(Fdiv ~ t2(sr), family = "gaussian", data = stackRas_PCA)
m5_pca_gam_raoq  <- gam(Fdiv ~ sr, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam_raoq ,m2_pca_gam_raoq , m3_pca_gam_raoq ,m4_pca_gam_raoq ,m5_pca_gam_raoq ) ##BEST=Fdiv~s(sr)

model_p_pca_raoq<-predict(m1_pca_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA_raoq <- data.frame(
  sr = stackRas_PCA$sr,
  fit = model_p_pca_raoq$fit,
  se.fit = model_p_pca_raoq$se.fit
)

#For RaoQ ~ Frich
#Choose best model
m1_pca_gam_raoq_fd <- gam(Fdiv ~ s(Frich),family = "gaussian",data = stackRas_PCA)
summary(m1_pca_gam_raoq_fd )

m2_pca_gam_raoq_fd <- gam(Fdiv ~ te(Frich), family = "gaussian", data = stackRas_PCA)
m3_pca_gam_raoq_fd  <- gam(Fdiv ~ ti(Frich), family = "gaussian", data = stackRas_PCA)
m4_pca_gam_raoq_fd  <- gam(Fdiv ~ t2(Frich), family = "gaussian", data = stackRas_PCA)
m5_pca_gam_raoq_fd  <- gam(Fdiv ~ Frich, family = "gaussian", data = stackRas_PCA)
AIC(m1_pca_gam_raoq_fd ,m2_pca_gam_raoq_fd , m3_pca_gam_raoq_fd ,m4_pca_gam_raoq_fd ,m5_pca_gam_raoq_fd) ##BEST=Fdiv~s(Frich)

model_p_pca_raoq_fd<-predict(m1_pca_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_PCA_raoq_fd <- data.frame(
  Frich = stackRas_PCA$Frich,
  fit = model_p_pca_raoq_fd$fit,
  se.fit = model_p_pca_raoq_fd$se.fit
)

###### Random within clusters ###### 
Frich_within_clusters<-rast("random_clusters/fdr_clusters.grd")
Fdiv_within_clusters<-rast("random_clusters/RaoQ_clusters.grd")

stackRas_clusters<-c(Fdiv_within_clusters,sr,MPD,Frich_within_clusters,pd)
names(stackRas_clusters)<-c("Fdiv","sr","MPD","Frich","pd")

stackRas_clusters<-as.data.frame(stackRas_clusters)
stackRas_clusters<-subset(stackRas_clusters,sr>3)
stackRas_clusters<-na.omit(stackRas_clusters)

#For Frich 
#Choose best model
m1_cluster_gam <- gam(Frich ~ s(sr),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam)

m2_cluster_gam <- gam(Frich ~ te(sr), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam <- gam(Frich ~ ti(sr), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam <- gam(Frich ~ t2(sr), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam <- gam(Frich ~ sr, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam,m2_cluster_gam, m3_cluster_gam,m4_cluster_gam,m5_cluster_gam) ##BEST=Frich~s(sr)

model_p_cluster<-predict(m1_cluster_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters <- data.frame(
  sr = stackRas_clusters$sr,
  fit = model_p_cluster$fit,
  se.fit = model_p_cluster$se.fit
)

#For Fdiv~sr
#Choose best model
m1_cluster_gam_raoq <- gam(Fdiv ~ s(sr),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam_raoq)

m2_cluster_gam_raoq <- gam(Fdiv ~ te(sr), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam_raoq <- gam(Fdiv ~ ti(sr), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam_raoq <- gam(Fdiv ~ t2(sr), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam_raoq <- gam(Fdiv ~ sr, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam_raoq,m2_cluster_gam_raoq, m3_cluster_gam_raoq,m4_cluster_gam_raoq,m5_cluster_gam_raoq) ##BEST=Fdiv~s(sr)

model_p_cluster_raoq<-predict(m1_cluster_gam_raoq, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters_raoq <- data.frame(
  sr = stackRas_clusters$sr,
  fit = model_p_cluster_raoq$fit,
  se.fit = model_p_cluster_raoq$se.fit
)

#For Fdiv~Fric
#Choose best model
m1_cluster_gam_raoq_fd <- gam(Fdiv ~ s(Frich),family = "gaussian",data = stackRas_clusters)
summary(m1_cluster_gam_raoq_fd)

m2_cluster_gam_raoq_fd <- gam(Fdiv ~ te(Frich), family = "gaussian", data = stackRas_clusters)
m3_cluster_gam_raoq_fd <- gam(Fdiv ~ ti(Frich), family = "gaussian", data = stackRas_clusters)
m4_cluster_gam_raoq_fd <- gam(Fdiv ~ t2(Frich), family = "gaussian", data = stackRas_clusters)
m5_cluster_gam_raoq_fd <- gam(Fdiv ~ Frich, family = "gaussian", data = stackRas_clusters)
AIC(m1_cluster_gam_raoq_fd,m2_cluster_gam_raoq_fd, m3_cluster_gam_raoq_fd,m4_cluster_gam_raoq_fd,m5_cluster_gam_raoq_fd) ##BEST=Fdiv~s(Frich)

model_p_cluster_raoq_fd<-predict(m1_cluster_gam_raoq_fd, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_clusters_raoq_fd <- data.frame(
  Frich = stackRas_clusters$Frich,
  fit = model_p_cluster_raoq_fd$fit,
  se.fit = model_p_cluster_raoq_fd$se.fit
)

###### Try for random_1 ###### 
Frich_random_1<-rast("random1/fdr_random1.grd")
Fdiv_random_1<-rast("random1/RaoQ_random1.grd")

stackRas_random1<-c(Fdiv_random_1,sr,MPD,Frich_random_1,pd)
names(stackRas_random1)<-c("Fdiv","sr","MPD","Frich","pd")

stackRas_random1<-as.data.frame(stackRas_random1)
stackRas_random1<-subset(stackRas_random1,sr>3)
stackRas_random1<-na.omit(stackRas_random1)

#Choose best model
m1_random_1_gam <- gam(Frich ~ s(sr),family = "gaussian",data = stackRas_random1)
summary(m1_random_1_gam)

m2_random_1_gam <- gam(Frich ~ te(sr), family = "gaussian", data = stackRas_random1)
m3_random_1_gam <- gam(Frich ~ ti(sr), family = "gaussian", data = stackRas_random1)
m4_random_1_gam <- gam(Frich ~ t2(sr), family = "gaussian", data = stackRas_random1)
m5_random_1_gam <- gam(Frich ~ sr, family = "gaussian", data = stackRas_random1)
AIC(m1_random_1_gam,m2_random_1_gam, m3_random_1_gam,m4_random_1_gam,m5_random_1_gam) ##BEST=Frich~s(sr)

model_p_random1<-predict(m1_random_1_gam, se.fit=TRUE)

#Create dataframe instead of list for ggplot
pred_data_random1 <- data.frame(
  sr = stackRas_random1$sr,
  fit = model_p_random1$fit,
  se.fit = model_p_random1$se.fit
)

##### TRY with Random X - FRich #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random<-c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_random <- as.data.frame(stackRas_random)
  stackRas_random <- subset(stackRas_random, sr > 3)
  stackRas_random <- na.omit(stackRas_random)
  
  #calculate all models
  models <- list(
    s_model = gam(Frich ~ s(sr), family = "gaussian", data = stackRas_random),
    te_model = gam(Frich ~ te(sr), family = "gaussian", data = stackRas_random),
    ti_model = gam(Frich ~ ti(sr), family = "gaussian", data = stackRas_random),
    t2_model = gam(Frich ~ t2(sr), family = "gaussian", data = stackRas_random),
    linear_model = gam(Frich ~ sr, family = "gaussian", data = stackRas_random)
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
  Frich_random <- rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random <- rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random <- c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_random <- as.data.frame(stackRas_random)
  stackRas_random <- subset(stackRas_random, sr > 3)
  stackRas_random <- na.omit(stackRas_random)
  
  #save 
  stackRas_list[[paste0("stackRas_random", i)]] <- stackRas_random
  
  #Gam models
  model_gam <- gam(Frich ~ s(sr), family = "gaussian", data = stackRas_random)
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

##### TRY with Random X - Fdiv~SR #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq<-c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random_raoq)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_random_raoq <- as.data.frame(stackRas_random_raoq)
  stackRas_random_raoq <- subset(stackRas_random_raoq, sr > 3)
  stackRas_random_raoq <- na.omit(stackRas_random_raoq)
  
  #calculate all models
  models <- list(
    s_model = gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_random_raoq),
    te_model = gam(Fdiv ~ te(sr), family = "gaussian", data = stackRas_random_raoq),
    ti_model = gam(Fdiv ~ ti(sr), family = "gaussian", data = stackRas_random_raoq),
    t2_model = gam(Fdiv ~ t2(sr), family = "gaussian", data = stackRas_random_raoq),
    linear_model = gam(Fdiv ~ sr, family = "gaussian", data = stackRas_random_raoq)
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
#s model= the best 

#Create lists to store models and data
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

#Loop for 10 randoms
for (i in 1:10) {
  #Load datas
  Frich_random <- rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random <- rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq <- c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random_raoq) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_random_raoq <- as.data.frame(stackRas_random_raoq)
  stackRas_random_raoq <- subset(stackRas_random_raoq, sr > 3)
  stackRas_random_raoq <- na.omit(stackRas_random_raoq)
  
  #Save data
  stackRas_list[[paste0("stackRas_random_raoq", i)]] <- stackRas_random_raoq
  
  #Gam models
  model_gam_raoq <- gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_random_raoq)
  model_list[[paste0("model_gam_random_raoq", i)]] <- model_gam_raoq
    model_p_random_raoq <- predict(model_gam_raoq, se.fit = TRUE)
  
  #Create dataset
  pred_data_random_raoq <- data.frame(
    sr = stackRas_random_raoq$sr,
    fit = model_p_random_raoq$fit,
    se.fit = model_p_random_raoq$se.fit
  )
  
  #Save dataset
  pred_data_list[[paste0("pred_data_random_raoq", i)]] <- pred_data_random_raoq
}

#stackRas_list$stackRas_random_raoq1 # stocke les datasets pour chaque random
#model_list$model_gam_random_raoq1 #stocke les modèles gam
#pred_data_list$pred_data_randomraoq1 #stocke les prédictions

##### TRY with Random X - Fdiv~Frich #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_random<-rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random<-rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq_fd<-c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random_raoq_fd)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_random_raoq_fd <- as.data.frame(stackRas_random_raoq_fd)
  stackRas_random_raoq_fd <- subset(stackRas_random_raoq_fd, sr > 3)
  stackRas_random_raoq_fd <- na.omit(stackRas_random_raoq_fd)
  
  #calculate all models
  models <- list(
    s_model = gam(Fdiv ~ s(Frich), family = "gaussian", data = stackRas_random_raoq_fd),
    te_model = gam(Fdiv ~ te(Frich), family = "gaussian", data = stackRas_random_raoq_fd),
    ti_model = gam(Fdiv ~ ti(Frich), family = "gaussian", data = stackRas_random_raoq_fd),
    t2_model = gam(Fdiv ~ t2(Frich), family = "gaussian", data = stackRas_random_raoq_fd),
    linear_model = gam(Fdiv ~ Frich, family = "gaussian", data = stackRas_random_raoq_fd)
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
  
  #Print results
  cat("Random", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}
#s model= the best 

##Create lists to store models and data
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

#Loop for 10 randoms
for (i in 1:10) {
  #Load data
  Frich_random <- rast(paste0("random", i, "/fdr_random", i, ".grd"))
  Fdiv_random <- rast(paste0("random", i, "/RaoQ_random", i, ".grd"))
  
  stackRas_random_raoq_fd <- c(Fdiv_random, sr, MPD, Frich_random, pd)
  names(stackRas_random_raoq_fd) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_random_raoq_fd <- as.data.frame(stackRas_random_raoq_fd)
  stackRas_random_raoq_fd <- subset(stackRas_random_raoq_fd, sr > 3)
  stackRas_random_raoq_fd <- na.omit(stackRas_random_raoq_fd)
  
  #Save
  stackRas_list[[paste0("stackRas_random_raoq_fd", i)]] <- stackRas_random_raoq_fd
  
  #Gam models
  model_gam_raoq_fd <- gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_random_raoq_fd)
  model_list[[paste0("model_gam_raoq_fd", i)]] <- model_gam_raoq_fd
  model_p_random_raoq_fd <- predict(model_gam_raoq_fd, se.fit = TRUE)
  
  #Create dataset
  pred_data_random_raoq_fd <- data.frame(
    Frich = stackRas_random_raoq_fd$Frich,
    fit = model_p_random_raoq_fd$fit,
    se.fit = model_p_random_raoq_fd$se.fit
  )
  
  #Save dataset
  pred_data_list[[paste0("pred_data_random_raoq_fd", i)]] <- pred_data_random_raoq_fd
}

#stackRas_list$stackRas_random_raoq_fd1 # stocke les datasets pour chaque random
#model_list$model_gam_random_raoq_fd1 #stocke les modèles gam
#pred_data_list$pred_data_randomraoq_fd1 #stocke les prédictions


##### Making the plots - Frich - with all #####
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=sr, y=Frich, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random1, aes(x=sr, y=Frich, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random1, aes(x=sr, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random2, aes(x=sr, y=Frich, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random2, aes(x=sr, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random3, aes(x=sr, y=Frich, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random3, aes(x=sr, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random4, aes(x=sr, y=Frich, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random4, aes(x=sr, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random5, aes(x=sr, y=Frich, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random5, aes(x=sr, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random6, aes(x=sr, y=Frich, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random6, aes(x=sr, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random7, aes(x=sr, y=Frich, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random7, aes(x=sr, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random8, aes(x=sr, y=Frich, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random8, aes(x=sr, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random9, aes(x=sr, y=Frich, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random9, aes(x=sr, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random10, aes(x=sr, y=Frich, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random10, aes(x=sr, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=sr, y=Frich, color="Random within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters, aes(x=sr, y=fit, color="Random within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=sr, y=Frich, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Richness") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen", 
               "Random 1"="royalblue", 
               "Random 2"="purple", 
               "Random 3"="pink",
               "Random 4"="hotpink",
               "Random 5"="cyan",
               "Random 6"="lightskyblue",
               "Random 7"="red",
               "Random 8"="tomato3",
               "Random 9"="darkgoldenrod2",
               "Random 10"="blue", 
               "PCA"="springgreen2",
               "Random within clusters"="olivedrab2"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR")

g

##### Making the plots - Fdiv~SR - with all #####
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=sr, y=Fdiv, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data_raoq, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random_raoq1, aes(x=sr, y=Fdiv, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq1, aes(x=sr, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random_raoq2, aes(x=sr, y=Fdiv, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq2, aes(x=sr, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random_raoq3, aes(x=sr, y=Fdiv, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq3, aes(x=sr, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random_raoq4, aes(x=sr, y=Fdiv, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq4, aes(x=sr, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random_raoq5, aes(x=sr, y=Fdiv, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq5, aes(x=sr, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random_raoq6, aes(x=sr, y=Fdiv, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq6, aes(x=sr, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random_raoq7, aes(x=sr, y=Fdiv, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq7, aes(x=sr, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random_raoq8, aes(x=sr, y=Fdiv, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq8, aes(x=sr, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random_raoq9, aes(x=sr, y=Fdiv, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq9, aes(x=sr, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random_raoq10, aes(x=sr, y=Fdiv, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq10, aes(x=sr, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=sr, y=Fdiv, color="Random within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters_raoq, aes(x=sr, y=fit, color="Random within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=sr, y=Fdiv, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA_raoq , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen", 
               "Random 1"="royalblue", 
               "Random 2"="purple", 
               "Random 3"="pink",
               "Random 4"="hotpink",
               "Random 5"="cyan",
               "Random 6"="lightskyblue",
               "Random 7"="red",
               "Random 8"="tomato3",
               "Random 9"="darkgoldenrod2",
               "Random 10"="blue", 
               "PCA"="springgreen2",
               "Random within clusters"="olivedrab2"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR")

g

##### Making the plots - Fdiv~Frich - with all #####
g <- ggplot()+
  #eco_based
  geom_point(data=stackRas, aes(x=Frich, y=Fdiv, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data_raoq_fd, aes(x=Frich, y=fit, color="Ecologically based")) +
  
  #random 1
  geom_point(data=stackRas_list$stackRas_random_raoq_fd1, aes(x=Frich, y=Fdiv, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd1, aes(x=Frich, y=fit, color="Random 1")) +
  
  #random 2
  geom_point(data=stackRas_list$stackRas_random_raoq_fd2, aes(x=Frich, y=Fdiv, color="Random 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd2, aes(x=Frich, y=fit, color="Random 2")) +
  
  #random 3
  geom_point(data=stackRas_list$stackRas_random_raoq_fd3, aes(x=Frich, y=Fdiv, color="Random 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd3, aes(x=Frich, y=fit, color="Random 3")) +
  
  #random 4
  geom_point(data=stackRas_list$stackRas_random_raoq_fd4, aes(x=Frich, y=Fdiv, color="Random 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd4, aes(x=Frich, y=fit, color="Random 4")) +
  
  #random 5
  geom_point(data=stackRas_list$stackRas_random_raoq_fd5, aes(x=Frich, y=Fdiv, color="Random 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd5, aes(x=Frich, y=fit, color="Random 5")) +
  
  #random 6
  geom_point(data=stackRas_list$stackRas_random_raoq_fd6, aes(x=Frich, y=Fdiv, color="Random 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd6, aes(x=Frich, y=fit, color="Random 6")) +
  
  #random 7
  geom_point(data=stackRas_list$stackRas_random_raoq_fd7, aes(x=Frich, y=Fdiv, color="Random 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd7, aes(x=Frich, y=fit, color="Random 7")) +
  
  #random 8
  geom_point(data=stackRas_list$stackRas_random_raoq_fd8, aes(x=Frich, y=Fdiv, color="Random 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd8, aes(x=Frich, y=fit, color="Random 8")) +
  
  #random 9
  geom_point(data=stackRas_list$stackRas_random_raoq_fd9, aes(x=Frich, y=Fdiv, color="Random 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd9, aes(x=Frich, y=fit, color="Random 9")) +
  
  #random 10
  geom_point(data=stackRas_list$stackRas_random_raoq_fd10, aes(x=Frich, y=Fdiv, color="Random 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_random_raoq_fd10, aes(x=Frich, y=fit, color="Random 10")) +
  
  #within_clusters
  geom_point(data=stackRas_clusters, aes(x=Frich, y=Fdiv, color="Random within clusters"), alpha=0.01) +
  geom_smooth(data=pred_data_clusters_raoq_fd, aes(x=Frich, y=fit, color="Random within clusters")) +
  
  #PCA
  geom_point(data=stackRas_PCA, aes(x=Frich, y=Fdiv, color="PCA"), alpha=0.01) +
  geom_smooth(data=pred_data_PCA_raoq_fd , aes(x=Frich, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen", 
               "Random 1"="royalblue", 
               "Random 2"="purple", 
               "Random 3"="pink",
               "Random 4"="hotpink",
               "Random 5"="cyan",
               "Random 6"="lightskyblue",
               "Random 7"="red",
               "Random 8"="tomato3",
               "Random 9"="darkgoldenrod2",
               "Random 10"="blue", 
               "PCA"="springgreen2",
               "Random within clusters"="olivedrab2"),
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Functional Richness (Frich - FDR)")

g

