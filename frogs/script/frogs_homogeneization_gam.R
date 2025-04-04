#####Modified by Tali Guez 2025
#Frogs

##First get coordinates from raster
library(terra)
library(ggplot2)
library(nlme)
library(mgcv)
#library(tidymv) #non compatible avec ma version de R
#install.packages("tidygam")
#library(tidygam)

rm(list = ls())

##Read in the maps
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_homogeneization") 
#setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps")

#sr<-rast("eco_based_200/SR_eco_based_200.grd")
#pd<-rast("eco_based_200/PD_eco_based_200.grd")
#MPD<-rast("eco_based_200/MPD_eco_based_200.grd")

sr<-rast("comb1/SR_comb1.grd")
pd<-rast("comb1/PD_comb1.grd")
MPD<-rast("comb1/MPD_comb1.grd")

##### Loading data + create models #####
##### TRY with Random X - FRich #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_comb<-rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb<-rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb<-c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_comb <- as.data.frame(stackRas_comb)
  stackRas_comb <- subset(stackRas_comb, sr > 3)
  stackRas_comb <- na.omit(stackRas_comb)
  
  #calculate all models
  models <- list(
    s_model = gam(Frich ~ s(sr), family = "gaussian", data = stackRas_comb),
    te_model = gam(Frich ~ te(sr), family = "gaussian", data = stackRas_comb),
    ti_model = gam(Frich ~ ti(sr), family = "gaussian", data = stackRas_comb),
    t2_model = gam(Frich ~ t2(sr), family = "gaussian", data = stackRas_comb),
    linear_model = gam(Frich ~ sr, family = "gaussian", data = stackRas_comb)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("Comb_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values
  )
  
  #Print results
  cat("Combination", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}

# Créez des listes pour stocker les données et modèles
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

# Boucle sur les 10 randoms
for (i in 1:10) {
  # Charger les données
  Frich_comb <- rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb <- rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb <- c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_comb <- as.data.frame(stackRas_comb)
  stackRas_comb <- subset(stackRas_comb, sr > 3)
  stackRas_comb <- na.omit(stackRas_comb)
  
  # Sauvegarde des données
  stackRas_list[[paste0("stackRas_comb", i)]] <- stackRas_comb
  
  # Ajuster le modèle GAM
  model_gam <- gam(Frich ~ s(sr), family = "gaussian", data = stackRas_comb)
  model_list[[paste0("model_gam_comb", i)]] <- model_gam
  
  # Prédictions
  model_p_comb <- predict(model_gam, se.fit = TRUE)
  
  # Créer un DataFrame de prédictions
  pred_data_comb <- data.frame(
    sr = stackRas_comb$sr,
    fit = model_p_comb$fit,
    se.fit = model_p_comb$se.fit
  )
  
  # Sauvegarde des prédictions
  pred_data_list[[paste0("pred_data_comb", i)]] <- pred_data_comb
}

# Vous pouvez ensuite accéder aux objets comme ceci :
# stackRas_list$stackRas_comb1 # stocke les datasets pour chaque random
# model_list$model_gam_comb1 #stocke les modèles gam
# pred_data_list$pred_data_comb1 #stocke les prédictions

##### TRY with Random X - Fdiv~SR #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_comb<-rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb<-rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb_raoq<-c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb_raoq)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_comb_raoq <- as.data.frame(stackRas_comb_raoq)
  stackRas_comb_raoq <- subset(stackRas_comb_raoq, sr > 3)
  stackRas_comb_raoq <- na.omit(stackRas_comb_raoq)
  
  #calculate all models
  models <- list(
    s_model = gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_comb_raoq),
    te_model = gam(Fdiv ~ te(sr), family = "gaussian", data = stackRas_comb_raoq),
    ti_model = gam(Fdiv ~ ti(sr), family = "gaussian", data = stackRas_comb_raoq),
    t2_model = gam(Fdiv ~ t2(sr), family = "gaussian", data = stackRas_comb_raoq),
    linear_model = gam(Fdiv ~ sr, family = "gaussian", data = stackRas_comb_raoq)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("comb_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values
  )
  
  #Print results
  cat("comb", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}
#s model= the best 

# Créez des listes pour stocker les données et modèles
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

# Boucle sur les 10 combs
for (i in 1:10) {
  # Charger les données
  Frich_comb <- rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb <- rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb_raoq <- c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb_raoq) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_comb_raoq <- as.data.frame(stackRas_comb_raoq)
  stackRas_comb_raoq <- subset(stackRas_comb_raoq, sr > 3)
  stackRas_comb_raoq <- na.omit(stackRas_comb_raoq)
  
  # Sauvegarde des données
  stackRas_list[[paste0("stackRas_comb_raoq", i)]] <- stackRas_comb_raoq
  
  # Ajuster le modèle GAM
  model_gam_raoq <- gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_comb_raoq)
  model_list[[paste0("model_gam_comb_raoq", i)]] <- model_gam_raoq
  
  # Prédictions
  model_p_comb_raoq <- predict(model_gam_raoq, se.fit = TRUE)
  
  # Créer un DataFrame de prédictions
  pred_data_comb_raoq <- data.frame(
    sr = stackRas_comb_raoq$sr,
    fit = model_p_comb_raoq$fit,
    se.fit = model_p_comb_raoq$se.fit
  )
  
  # Sauvegarde des prédictions
  pred_data_list[[paste0("pred_data_comb_raoq", i)]] <- pred_data_comb_raoq
}

# stackRas_list$stackRas_comb_raoq1 # stocke les datasets pour chaque comb
# model_list$model_gam_comb_raoq1 #stocke les modèles gam
# pred_data_list$pred_data_comb_raoq1 #stocke les prédictions

##### TRY with Random X - Fdiv~Frich #####
#Calcul AIC, choose best model 
best_models<-list() #list to store the results

for (i in 1:10){
  #download data
  Frich_comb<-rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb<-rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb_raoq_fd<-c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb_raoq_fd)<- c("Fdiv","sr","MPD","Frich","pd")
  
  stackRas_comb_raoq_fd <- as.data.frame(stackRas_comb_raoq_fd)
  stackRas_comb_raoq_fd <- subset(stackRas_comb_raoq_fd, sr > 3)
  stackRas_comb_raoq_fd <- na.omit(stackRas_comb_raoq_fd)
  
  #calculate all models
  models <- list(
    s_model = gam(Fdiv ~ s(Frich), family = "gaussian", data = stackRas_comb_raoq_fd),
    te_model = gam(Fdiv ~ te(Frich), family = "gaussian", data = stackRas_comb_raoq_fd),
    ti_model = gam(Fdiv ~ ti(Frich), family = "gaussian", data = stackRas_comb_raoq_fd),
    t2_model = gam(Fdiv ~ t2(Frich), family = "gaussian", data = stackRas_comb_raoq_fd),
    linear_model = gam(Fdiv ~ Frich, family = "gaussian", data = stackRas_comb_raoq_fd)
  )
  
  #AIC
  aic_values <- sapply(models, AIC)
  
  #Lowest AIC
  best_model_name <- names(aic_values)[which.min(aic_values)]
  best_model <- models[[best_model_name]]
  
  #Store the best models
  best_models[[paste0("comb_", i)]] <- list(
    best_model_name = best_model_name,
    best_model = best_model,
    aic_values = aic_values
  )
  
  #Print results
  cat("comb", i, "\n")
  print(aic_values)
  cat("Best model:", best_model_name, "with AIC =", min(aic_values), "\n\n")
}
#s model= the best 

# Créez des listes pour stocker les données et modèles
model_list <- list()
pred_data_list <- list()
stackRas_list <- list()

# Boucle sur les 10 combs
for (i in 1:10) {
  # Charger les données
  Frich_comb <- rast(paste0("comb", i, "/fdr_comb", i, ".grd"))
  Fdiv_comb <- rast(paste0("comb", i, "/RaoQ_comb", i, ".grd"))
  
  stackRas_comb_raoq_fd <- c(Fdiv_comb, sr, MPD, Frich_comb, pd)
  names(stackRas_comb_raoq_fd) <- c("Fdiv", "sr", "MPD", "Frich", "pd")
  
  stackRas_comb_raoq_fd <- as.data.frame(stackRas_comb_raoq_fd)
  stackRas_comb_raoq_fd <- subset(stackRas_comb_raoq_fd, sr > 3)
  stackRas_comb_raoq_fd <- na.omit(stackRas_comb_raoq_fd)
  
  # Sauvegarde des données
  stackRas_list[[paste0("stackRas_comb_raoq_fd", i)]] <- stackRas_comb_raoq_fd
  
  # Ajuster le modèle GAM
  model_gam_raoq_fd <- gam(Fdiv ~ s(sr), family = "gaussian", data = stackRas_comb_raoq_fd)
  model_list[[paste0("model_gam_raoq_fd", i)]] <- model_gam_raoq_fd
  
  # Prédictions
  model_p_comb_raoq_fd <- predict(model_gam_raoq_fd, se.fit = TRUE)
  
  # Créer un DataFrame de prédictions
  pred_data_comb_raoq_fd <- data.frame(
    Frich = stackRas_comb_raoq_fd$Frich,
    fit = model_p_comb_raoq_fd$fit,
    se.fit = model_p_comb_raoq_fd$se.fit
  )
  
  # Sauvegarde des prédictions
  pred_data_list[[paste0("pred_data_comb_raoq_fd", i)]] <- pred_data_comb_raoq_fd
}

# stackRas_list$stackRas_comb_raoq_fd1 # stocke les datasets pour chaque comb
# model_list$model_gam_comb_raoq_fd1 #stocke les modèles gam
# pred_data_list$pred_data_comb_raoq_fd1 #stocke les prédictions


##### Making the plots - Frich - with all #####
g <- ggplot()+
  #eco_based
  #geom_point(data=stackRas, aes(x=sr, y=Frich, color="Ecologically based"), alpha=0.01) +
  #geom_smooth(data=pred_data, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #comb 1
  geom_point(data=stackRas_list$stackRas_comb1, aes(x=sr, y=Frich, color="comb 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb1, aes(x=sr, y=fit, color="comb 1")) +
  
  #comb 2
  geom_point(data=stackRas_list$stackRas_comb2, aes(x=sr, y=Frich, color="comb 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb2, aes(x=sr, y=fit, color="comb 2")) +
  
  #comb 3
  geom_point(data=stackRas_list$stackRas_comb3, aes(x=sr, y=Frich, color="comb 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb3, aes(x=sr, y=fit, color="comb 3")) +
  
  #comb 4
  geom_point(data=stackRas_list$stackRas_comb4, aes(x=sr, y=Frich, color="comb 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb4, aes(x=sr, y=fit, color="comb 4")) +
  
  #comb 5
  geom_point(data=stackRas_list$stackRas_comb5, aes(x=sr, y=Frich, color="comb 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb5, aes(x=sr, y=fit, color="comb 5")) +
  
  #comb 6
  geom_point(data=stackRas_list$stackRas_comb6, aes(x=sr, y=Frich, color="comb 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb6, aes(x=sr, y=fit, color="comb 6")) +
  
  #comb 7
  geom_point(data=stackRas_list$stackRas_comb7, aes(x=sr, y=Frich, color="comb 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb7, aes(x=sr, y=fit, color="comb 7")) +
  
  #comb 8
  geom_point(data=stackRas_list$stackRas_comb8, aes(x=sr, y=Frich, color="comb 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb8, aes(x=sr, y=fit, color="comb 8")) +
  
  #comb 9
  geom_point(data=stackRas_list$stackRas_comb9, aes(x=sr, y=Frich, color="comb 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb9, aes(x=sr, y=fit, color="comb 9")) +
  
  #comb 10
  geom_point(data=stackRas_list$stackRas_comb10, aes(x=sr, y=Frich, color="comb 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb10, aes(x=sr, y=fit, color="comb 10")) +
  
  #within_clusters
  #geom_point(data=stackRas_clusters, aes(x=sr, y=Frich, color="Random within clusters"), alpha=0.01) +
  #geom_smooth(data=pred_data_clusters, aes(x=sr, y=fit, color="Random within clusters")) +
  
  #PCA
  #geom_point(data=stackRas_PCA, aes(x=sr, y=Frich, color="PCA"), alpha=0.01) +
  #geom_smooth(data=pred_data_PCA , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Richness") +
  scale_color_manual(
    values = c(#"Ecologically based" = "darkgreen", 
               "comb 1"="darkgreen", 
               "comb 2"="purple", 
               "comb 3"="pink",
               "comb 4"="hotpink",
               "comb 5"="cyan",
               "comb 6"="lightskyblue",
               "comb 7"="red",
               "comb 8"="tomato3",
               "comb 9"="darkgoldenrod2",
               "comb 10"="royalblue"#, 
               #"PCA"="springgreen2",
               ),#"Random within clusters"="olivedrab2"
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR")

g

##### Making the plots - Fdiv~SR - with all #####
g <- ggplot()+
  #eco_based
  #geom_point(data=stackRas, aes(x=sr, y=Fdiv, color="Ecologically based"), alpha=0.01) +
  #geom_smooth(data=pred_data_raoq, aes(x=sr, y=fit, color="Ecologically based")) +
  
  #comb 1
  geom_point(data=stackRas_list$stackRas_comb_raoq1, aes(x=sr, y=Fdiv, color="comb 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq1, aes(x=sr, y=fit, color="comb 1")) +
  
  #comb 2
  geom_point(data=stackRas_list$stackRas_comb_raoq2, aes(x=sr, y=Fdiv, color="comb 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq2, aes(x=sr, y=fit, color="comb 2")) +
  
  #comb 3
  geom_point(data=stackRas_list$stackRas_comb_raoq3, aes(x=sr, y=Fdiv, color="comb 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq3, aes(x=sr, y=fit, color="comb 3")) +
  
  #comb 4
  geom_point(data=stackRas_list$stackRas_comb_raoq4, aes(x=sr, y=Fdiv, color="comb 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq4, aes(x=sr, y=fit, color="comb 4")) +
  
  #comb 5
  geom_point(data=stackRas_list$stackRas_comb_raoq5, aes(x=sr, y=Fdiv, color="comb 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq5, aes(x=sr, y=fit, color="comb 5")) +
  
  #comb 6
  geom_point(data=stackRas_list$stackRas_comb_raoq6, aes(x=sr, y=Fdiv, color="comb 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq6, aes(x=sr, y=fit, color="comb 6")) +
  
  #comb 7
  geom_point(data=stackRas_list$stackRas_comb_raoq7, aes(x=sr, y=Fdiv, color="comb 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq7, aes(x=sr, y=fit, color="comb 7")) +
  
  #comb 8
  geom_point(data=stackRas_list$stackRas_comb_raoq8, aes(x=sr, y=Fdiv, color="comb 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq8, aes(x=sr, y=fit, color="comb 8")) +
  
  #comb 9
  geom_point(data=stackRas_list$stackRas_comb_raoq9, aes(x=sr, y=Fdiv, color="comb 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq9, aes(x=sr, y=fit, color="comb 9")) +
  
  #comb 10
  geom_point(data=stackRas_list$stackRas_comb_raoq10, aes(x=sr, y=Fdiv, color="comb 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq10, aes(x=sr, y=fit, color="comb 10")) +
  
  #within_clusters
  #geom_point(data=stackRas_clusters, aes(x=sr, y=Fdiv, color="Random within clusters"), alpha=0.01) +
  #geom_smooth(data=pred_data_clusters_raoq, aes(x=sr, y=fit, color="Random within clusters")) +
  
  #PCA
  #geom_point(data=stackRas_PCA, aes(x=sr, y=Fdiv, color="PCA"), alpha=0.01) +
  #geom_smooth(data=pred_data_PCA_raoq , aes(x=sr, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c(#"Ecologically based" = "royalblue", 
               "comb 1"="darkgreen", 
               "comb 2"="purple", 
               "comb 3"="pink",
               "comb 4"="hotpink",
               "comb 5"="cyan",
               "comb 6"="lightskyblue",
               "comb 7"="red",
               "comb 8"="tomato3",
               "comb 9"="darkgoldenrod2",
               "comb 10"="blue"),#, 
               #"PCA"="springgreen2",
              # ), #"Random within clusters"="olivedrab2"
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Taxonomic diversity SR")

g

##### Making the plots - Fdiv~Frich - with all #####
g <- ggplot()+
  #eco_based
  #geom_point(data=stackRas, aes(x=Frich, y=Fdiv, color="Ecologically based"), alpha=0.01) +
  #geom_smooth(data=pred_data_raoq_fd, aes(x=Frich, y=fit, color="Ecologically based")) +
  
  #comb 1
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd1, aes(x=Frich, y=Fdiv, color="comb 1"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd1, aes(x=Frich, y=fit, color="comb 1")) +
  
  #comb 2
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd2, aes(x=Frich, y=Fdiv, color="comb 2"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd2, aes(x=Frich, y=fit, color="comb 2")) +
  
  #comb 3
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd3, aes(x=Frich, y=Fdiv, color="comb 3"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd3, aes(x=Frich, y=fit, color="comb 3")) +
  
  #comb 4
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd4, aes(x=Frich, y=Fdiv, color="comb 4"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd4, aes(x=Frich, y=fit, color="comb 4")) +
  
  #comb 5
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd5, aes(x=Frich, y=Fdiv, color="comb 5"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd5, aes(x=Frich, y=fit, color="comb 5")) +
  
  #comb 6
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd6, aes(x=Frich, y=Fdiv, color="comb 6"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd6, aes(x=Frich, y=fit, color="comb 6")) +
  
  #comb 7
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd7, aes(x=Frich, y=Fdiv, color="comb 7"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd7, aes(x=Frich, y=fit, color="comb 7")) +
  
  #comb 8
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd8, aes(x=Frich, y=Fdiv, color="comb 8"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd8, aes(x=Frich, y=fit, color="comb 8")) +
  
  #comb 9
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd9, aes(x=Frich, y=Fdiv, color="comb 9"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd9, aes(x=Frich, y=fit, color="comb 9")) +
  
  #comb 10
  geom_point(data=stackRas_list$stackRas_comb_raoq_fd10, aes(x=Frich, y=Fdiv, color="comb 10"), alpha=0.01) +
  geom_smooth(data=pred_data_list$pred_data_comb_raoq_fd10, aes(x=Frich, y=fit, color="comb 10")) +
  
  #within_clusters
  #geom_point(data=stackRas_clusters, aes(x=Frich, y=Fdiv, color="Random within clusters"), alpha=0.01) +
  #geom_smooth(data=pred_data_clusters_raoq_fd, aes(x=Frich, y=fit, color="Random within clusters")) +
  
  #PCA
  #geom_point(data=stackRas_PCA, aes(x=Frich, y=Fdiv, color="PCA"), alpha=0.01) +
  #geom_smooth(data=pred_data_PCA_raoq_fd , aes(x=Frich, y=fit, color="PCA")) +
  
  scale_y_continuous(name="Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c(#"Ecologically based" = "", 
               "comb 1"="darkgreen", #royalblue
               "comb 2"="purple", 
               "comb 3"="pink",
               "comb 4"="hotpink",
               "comb 5"="cyan",
               "comb 6"="lightskyblue",
               "comb 7"="red",
               "comb 8"="tomato3",
               "comb 9"="darkgoldenrod2",
               "comb 10"="blue"#, 
               #"PCA"="springgreen2",
              ), # "Random within clusters"="olivedrab2"
    name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") + #, legend.key=element_rect(fill="white", colour="black")
  labs(x="Functional Richness (Frich - FDR)")

g


##### Making the plots - Frich - test #####
g <- ggplot()+
  # eco_based
  geom_point(data=stackRas, aes(x=sr, y=Frich, color="Ecologically based"), alpha=0.01) +
  geom_smooth(data=pred_data, aes(x=sr, y=fit, color="Ecologically based")) +
  
  # random1
  geom_point(data=stackRas_random1, aes(x=sr, y=Frich, color="Random 1"), alpha=0.01) +
  geom_smooth(data=pred_data_random1, aes(x=sr, y=fit, color="Random 1")) +
  
  scale_y_continuous(name="Functional Richness") +
  scale_color_manual(values=c("Ecologically based"="darkgreen", "Random 1"="blue"), 
                     name = "Dataset") +
  coord_cartesian() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="bottom") +
  labs(x="Taxonomic diversity SR")

g

scale=1/60
pred_data %>% 
  ggplot(aes(sr, fit/scale)) + geom_smooth(color="darkgreen")+
  geom_point(data=stackRas,aes(y = Frich/scale),color="darkgreen",alpha=0.01)+
  scale_y_continuous(name="Functional Richness")+
  coord_cartesian()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y.right = element_text(color = "darkgreen"),axis.title.y.right = element_text(color = "darkgreen"),axis.line.y.right=element_line(color = "darkgreen"))+
  labs( x = "Taxonomic diversity SR")

pred_data %>% 
  ggplot(aes(sr, fit)) + geom_smooth(color="darkgreen")+
  geom_point(data=stackRas,aes(y = Frich),color="darkgreen",alpha=0.01)+
  scale_y_continuous(name="Functional Richness")+
  coord_cartesian()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y.right = element_text(color = "darkgreen"),axis.title.y.right = element_text(color = "darkgreen"),axis.line.y.right=element_line(color = "darkgreen"))+
  labs(x = "Taxonomic diversity SR")

#add random1
# Visualisation combinée pour eco_based et random_1
pred_data %>%
  ggplot(aes(sr, fit)) + geom_smooth(color="darkgreen") +
  geom_point(data=stackRas, aes(y = Frich),  color="darkgreen", alpha=0.01) +
  
  # Ajouter l'analyse pour random_1
  geom_smooth(data=pred_data_random1, aes(x = sr, y = fit), color="blue") +
  geom_point(data=stackRas_random1, aes(y = Frich), color="blue", alpha=0.01) +
  
  scale_y_continuous(name="Functional Richness") +
  coord_cartesian() +
  theme_bw() +
  scale_color_discrete(name = "ltitle") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(x = "Taxonomic diversity SR")



