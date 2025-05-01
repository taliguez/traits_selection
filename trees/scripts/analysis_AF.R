#####Tali Guez - 2025
###Crop map of the trees to the extend of the Atlantic forest 

####LIBRARIES ####
library(terra)
library(tidyverse)

library(geometry) #For convex hull calculations (FRic)
library(doParallel) # For parallel computing
library(vegan) #For functional diversity analysis
library(picante) #For phylogenetic diversity analysis
library(feather) #Efficient data storage format

rm(list=ls())

#### Creation of the AF raster for the trees ####
######Initialization######
mocklayer<-rast("~/Master/M2/Internship_M2/analyse/trees/figures/maps/eco_based_200/SR_eco_based_200.grd")
plot(mocklayer)
commu<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/raw/REV_Community_matrix_200.csv")
mocklayer<-init(mocklayer,"cell") #??
names(mocklayer) <- "Grid"
#r<-mocklayer
#values(r)<-NA

trees<-mocklayer
values(trees)<-NA
trees[as.numeric(commu$grid_id)]<-commu$grid_id #
plot(trees)

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")
hylids<-rast("richness_hylids_all.grd")
names(hylids) <- "Grid"
plot(hylids)

#try_projected<-project(trees, hylids, method="near") #, method="near"
#plot(try_projected)

#try_projected<-project(hylids, trees, method="near") #, method="near"
#plot(try_projected)

#extent_hylids <- ext(hylids)
#try_cropped <- crop(try_projected, extent_hylids)
#plot(try_cropped)
#try_masked <- mask(try_cropped, r)
#plot(try_masked)


###### Crop the world map to the Atlantic forest ######
trees_masking<-project(hylids, trees, method="near")
plot(trees_masking)
extent_hylids <- ext(hylids)
trees_extent <- crop(trees_masking, extent_hylids)
plot(trees_extent)

trees_AF<-crop(trees, trees_extent)
plot(trees_AF)
trees_AF<-mask(trees_AF, trees_extent)
plot(trees_AF)

###### Get the grids ######
grids_in_forest <- values(trees_AF, na.rm = TRUE)
grids_in_forest <- as.numeric(grids_in_forest)
grids_in_forest

#Filtering commu dataset that match grids_in_forest
commu_filtered <- commu %>%
  filter(grid_id %in% grids_in_forest)

#Verification of the grids
result_raster <- trees
plot(result_raster)
values(result_raster) <- 0 #ocean
result_raster[as.numeric(commu$grid_id)] <- 1 #rest of the continent
result_raster[as.numeric(commu_filtered$grid_id)] <- 2 #atlantic forest

plot(result_raster, col=c("white", "grey", "green4")) #match OK 

#setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps_AF") 
#writeRaster(trees_AF,"raster_base.grd",filetype="RRASTER", overwrite=TRUE) #save raster base for after
plot(trees_AF)

####Calculations of the indices ####
#read in the data
setwd("~/Master/M2/Internship_M2/analyse/trees/datas")
#suffix_200<-"_200"
#comm0 <- read_csv(paste0("raw/REV_Community_matrix",suffix_200,".csv"))  #Community matrix ,save_suffix,
comm0<-commu_filtered
#setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results_AF")
#write_csv(comm0, "REV_Community_matrix_AF_trees.csv")

pruned_tree <- read.tree("raw/REV_Pruned_tree.csv") #Pruned phylogenetic tree
ag <- read_csv("raw/ANGIO_GYMNO_lookup.csv") %>% 
  mutate(accepted_bin = gsub(" " , "_", accepted_bin)) %>% select(accepted_bin, group) #Species classification (Angiosperms/Gymnosperms)

dt_mat <- read_csv("clean/results/datasets_traits_comm/eco_based_scaled_200.csv")  #Trait matrix -> eco based
dt_mat <- read_csv("clean/results/datasets_traits_comm/random_based_200_10.csv")  #Trait matrix -> eco based
dt_mat <- read_csv("clean/results/datasets_traits_comm/cluster_random_200.csv")
dt_mat <- read_csv("clean/results/datasets_traits_comm/cluster_6_200.csv") #5 & 6
dt_mat <- read_csv("clean/results/datasets_traits_comm/pca_based_200.csv")
dt_mat <- read_csv("clean/results/datasets_traits_comm/pca_based_90_200.csv") #90

save_suffix <- "AF_pca_90" #_200 #change depending on the dataset used to not overwrite
#randomX
#pca_based
#cluster_random
#cluster_3 (5 & 6)
#

# register the parallel cores
registerDoParallel(round(detectCores()/2)) #Uses half of the available cores
getDoParWorkers()

## raoq formula. assumes equal abundances
pairdist <- function(myd){
  return(sum((myd / nrow(myd))^2))
}

#Only analysing Angiosperms
fit_seq <- c(2)  
set.seed(555)

for(gg in fit_seq){
  
  print("Fitting ANGIOSPERMS")
  comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
    filter(group == "Angiosperms")
  outf <- paste0("clean/results_AF/REV_obs_results_ANGIO_", save_suffix, ".csv") #à changer pour adapter aux traits
  
  # Get the plots and shuffle
  plts <- sample(unique(comm$grid_id), length(unique(comm$grid_id)), replace = FALSE)
  
  # Get number of traits
  K <- ncol(dt_mat)-1
  
  # Get the number of species per plot
  nspp <- comm %>% 
    select(grid_id, accepted_bin) %>% 
    distinct() %>% 
    group_by(grid_id) %>% 
    tally() %>% 
    ungroup %>% 
    arrange(desc(n))
  
  # Get the unique species
  spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()
  
  ###################################################.
  # Compute observed FD/PD in parallel
  ###################################################.
  outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE, .packages = c("dplyr", "tidyverse", "geometry", "picante", "ape") ) %dopar%{ ##dopar
    
    set.seed(555)
    
    #Get the focal plot
    my_dt <- comm %>% filter(grid_id == plts[j])
    
    #Initialize variables
    my_dt_sub <- NULL
    pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA ####
    out_res <- c()
    
    if(nrow(my_dt) > 1){
      
      #Subset the data and get the phylometrics
      my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% 
        data.frame() %>% 
        as_tibble() %>% 
        setNames(my_dt$accepted_bin) %>% 
        as.matrix()
      
      pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
      mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))
      
      # Get the trait matrix
      my_tr <- dt_mat %>% 
        filter(accepted_bin %in% my_dt$accepted_bin) %>% 
        select(-accepted_bin) %>% 
        as.matrix()
      
      my_tr <- my_tr + matrix(runif(nrow(my_tr) * ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))
      
      # Compute RaoQ
      my_raoq <- pairdist(as.matrix(dist(my_tr)))
      
      # Compute convex hull (FDRich) if enough unique species
      if(nrow(unique(my_dt)) > ncol(my_tr)){
        my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
        #print("my_fric is calculated")
      }

    } else {
      my_raoq <- pd_cur <- mpd_cur <- 0
    }
    
    # Save results
    out_res <- c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fd" = my_fric, "fdr" = my_fric^(1/K))
    
    return(out_res)
  }
  
  # Save the results
  outmat_check <- outmat_check %>% data.frame() %>% as_tibble() 
  write_csv(outmat_check, outf)
}

stopImplicitCluster()  

#####
rm(list=ls())
#Read file
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results_AF")
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_eco_based.csv") 
#summary(all_metrics)
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_random10.csv") #change names for the rest
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_cluster_random.csv") 
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_cluster_6.csv")
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_pca_based.csv") 
all_metrics<-read_csv("REV_obs_results_ANGIO_AF_pca_40.csv") 

all_metrics$raoq[all_metrics$n<9]<-NA #if n<9, put NA in raoq as fd and fdr

setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps_AF") 
trees_AF<-rast("raster_base.grd")
plot(trees_AF)

#### MAPPING ####
r<-trees_AF
r<-init(r,"cell")
names(r)<-"Grid"
values(r)<-NA

###Create all metric rasters for observed
pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r

cell_df<-as.data.frame(trees_AF, cells=TRUE, na.rm=FALSE) #make a df with the grid_id to the actual grid_number
colnames(cell_df)[2]<-"grid_id" #change name of the column into grid_id
all_metrics_complete<-left_join(cell_df, all_metrics, by = "grid_id") #Create dataset with all_metrics for each pixels

#Create vector NA with the same size of the raster base
pd_vals<-rep(NA, ncell(trees_AF))
mpd_vals<-rep(NA, ncell(trees_AF))
raoq_vals<-rep(NA, ncell(trees_AF))
SR_vals<-rep(NA, ncell(trees_AF))
fdr_vals<-rep(NA, ncell(trees_AF))
fd_vals<-rep(NA, ncell(trees_AF))

#Fill up the values to the corresponding cells
pd_vals[all_metrics_complete$cell]<-all_metrics_complete$pd
mpd_vals[all_metrics_complete$cell]<-all_metrics_complete$mpd
raoq_vals[all_metrics_complete$cell]<-all_metrics_complete$raoq
SR_vals[all_metrics_complete$cell]<-all_metrics_complete$n
fdr_vals[all_metrics_complete$cell]<-all_metrics_complete$fdr
fd_vals[all_metrics_complete$cell]<-all_metrics_complete$fd

#Values in the rasters
values(pd_ras)<-pd_vals
values(mpd_ras)<-mpd_vals
values(raoq_ras)<-raoq_vals
values(SR_ras)<-SR_vals
values(fdr_ras)<-fdr_vals
values(fd_ras)<-fd_vals

###### Plots ######
plot(pd_ras)
plot(SR_ras)
plot(fd_ras)
plot(raoq_ras)
plot(fdr_ras)
plot(mpd_ras)

###### Write rasters ######
setwd("~/Master/M2/Internship_M2/analyse/trees/figures/maps_AF/pca_40") #eco_based_200 / randomX / random_clusters / PCA_selection
#randomX = randomX at 200 filtering
#writeRaster(x=pd_ras,filename="PD_eco_based.tif",filetype="GTiff", overwrite=TRUE) #PD_all_trees.tif
writeRaster(pd_ras,filename="PD_AF.grd",filetype="RRASTER", overwrite=TRUE) #PD_all_trees.tif
writeRaster(mpd_ras,"MPD_AF.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(SR_ras,"SR_AF.grd",filetype="RRASTER", overwrite=TRUE) #

writeRaster(raoq_ras,"RaoQ_AF_pca40.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fdr_ras,"fdr_AF_pca40.grd",filetype="RRASTER", overwrite=TRUE) #
writeRaster(fd_ras,"fd_AF_pca40.grd",filetype="RRASTER", overwrite=TRUE) #



