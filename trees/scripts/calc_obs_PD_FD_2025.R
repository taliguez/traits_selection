######Script calculation diversity indices for trees
#Modified by Tali Guez 2025 (M2)
#from https://github.com/andrepazv/Tree_diversity/blob/main/code/calc_diversity/calc_obs_PD_FD.R

rm(list = ls())
gc() #Garbage collection

library(geometry) #For convex hull calculations (FRic)
library(doParallel) # For parallel computing
library(vegan) #For functional diversity analysis
library(picante) #For phylogenetic diversity analysis
library(feather) #Efficient data storage format
library(tidyverse) #Data manipulation and visualization
library(dplyr)

# register the parallel cores
registerDoParallel(round(detectCores()/2)) #Uses half of the available cores
#registerDoParallel(round(detectCores()-2)) 
getDoParWorkers()
#registerDoParallel(2) #Uses half of the available cores

## raoq formula. assumes equal abundances
pairdist <- function(myd){
  return(sum((myd / nrow(myd))^2))
}

#setwd("~/Git/Tree_diversity/data")
setwd("~/Master/M2/Internship_M2/analyse/trees/datas")

#randomX (X: 1:10)
suffix_200<-"_200"

# read in the data
comm0 <- read_csv(paste0("raw/REV_Community_matrix",suffix_200,".csv"))  #Community matrix ,save_suffix,
#comm0 <- comm0 %>% filter(grid_id %in% sample(unique(grid_id), 5)) #Only 5 sites --> test
#dt_mat <- read_csv("clean/results/datasets_traits_comm/eco_based_scaled.csv")  #Trait matrix 
pruned_tree <- read.tree("raw/REV_Pruned_tree.csv") #Pruned phylogenetic tree
ag <- read_csv("raw/ANGIO_GYMNO_lookup.csv") %>% 
  mutate(accepted_bin = gsub(" " , "_", accepted_bin)) %>% select(accepted_bin, group) #Species classification (Angiosperms/Gymnosperms)
#table(ag$group)

###test 
#test_matrix <- matrix(runif(20, -1, 1), ncol = 2) # 10 points en 2D
#print(convhulln(test_matrix, "FA")$vol) # Vérifier si ça marche

#Choose dt_mat in function of the traits chosen
dt_mat <- read_csv("clean/results/datasets_traits_comm/eco_based_scaled.csv")  #Trait matrix -> eco based
dt_mat <- read_csv("clean/results/datasets_traits_comm/eco_based_scaled_200.csv")  #Trait matrix -> eco based
dt_mat <- read_csv("clean/results/datasets_traits_comm/random_based_200_10.csv")  #Trait matrix -> random X
dt_mat <- read_csv("clean/results/datasets_traits_comm/pca_based_200.csv")  #Trait matrix -> pca based (8 first PC)
dt_mat <- read_csv("clean/results/datasets_traits_comm/cluster_random_200.csv")  #Trait matrix -> pca based (8 first PC)

save_suffix <- "cluster_random" #_200 #change depending on the dataset used to not overwrite
#eco_based
#pca_based
#cluster_random

#Only analysing Angiosperms
fit_seq <- c(2)  
set.seed(555)

for(gg in fit_seq){
  
  print("Fitting ANGIOSPERMS")
  comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
    filter(group == "Angiosperms")
  outf <- paste0("clean/results/REV_obs_results_ANGIO_", save_suffix, suffix_200, ".csv") #à changer pour adapter aux traits
  
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
  
  ###################################################
  # Compute observed FD/PD in parallel
  ###################################################
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
      
      #print(paste("grid_id:", plts[j], 
       #           "- n espèces uniques:", nrow(unique(my_dt)), 
        #          "- n traits:", ncol(my_tr)))
      
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

##### OPEN FILES ##### 
#test<-read_csv("clean/results/REV_obs_results_ANGIO_test.csv") #change names for the rest

eco_based_results<-read_csv("clean/results/REV_obs_results_ANGIO_eco_based_200.csv") #change names for the rest
summary(eco_based_results$fd)
random1<-read_csv("clean/results/REV_obs_results_ANGIO_random_9_200.csv") #change names for the rest
summary(random1$fd)


#########################################
###### Whole original script ############
#########################################

#inv <- read_csv("raw_data/Invasive_glonaf.csv") #List of invasive species
#print(paste0("Using: ", "REV_Community_matrix",save_suffix,".csv")) #Prints the name of the community matrix file being used.

# which null models to do -- 1=all, 2=angio, 3=gymno, 4=no invasives
#fit_seq <- c(4) #Only analyzing the case where invasive species are excluded
fit_seq<-c(2) #only analyzing angiospermes
gg <- 2

#Filtering species based on the group being studied
#This section filters species based on the group being analyzed:
#Angiosperms (gg == 2)
#Gymnosperms (gg == 3)
#Excluding invasive species (gg == 4)
#All species (else)

for(gg in fit_seq){
  
  if(gg == 2){
    print("fitting ANGIOS")
    comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
      filter(group == "Angiosperms")
    outf <- paste0("results/REV_obs_results_ANGIO", save_suffix, ".csv")
  }else if(gg == 3){
    print("fitting GYMNOS")
    comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
      filter(group == "Gymnosperms")
    outf <- paste0("results/REV_obs_results_GYMNO", save_suffix, ".csv")
  }else if(gg == 4){
    print("removing INVASIVES")
    comm <- comm0 %>% filter(!accepted_bin%in%inv$x)
    outf <- paste0("results/REV_obs_results_INV", save_suffix, ".csv")
  }else{
    print("fitting ALL")
    comm <- comm0
    outf <- paste0("results/REV_obs_results", save_suffix, ".csv")
  }
  
  
  # get the plots and shuffle
  plts <- sample(unique(comm$grid_id), length(unique(comm$grid_id)), replace = FALSE)
  
  # get no. of traits #exclude species name
  K <- ncol(dt_mat)-1
  
  # get the number of species per plot
  nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))
  
  # get the unique species
  spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()
  
  ###################################################
  # get the observed fd/pd in parallel
  ###################################################
  #Parallel loop foreach to compute diversity metrics for each grid.
  outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
    
    # get the focal plot #extract data from the current grid 
    my_dt <- comm %>% filter(grid_id == plts[j])
    
    # initialize the variables
    my_dt_sub <- NULL
    pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA
    
    out_res <- c()
    
    if(nrow(my_dt)>1){
      
      # subset the data and get the phylometrics
      my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% data.frame() %>% as_tibble() %>% setNames(my_dt$accepted_bin) %>% as.matrix()
      pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
      mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))
      
      # get the trait matrix
      my_tr <- dt_mat %>% filter(accepted_bin%in%my_dt$accepted_bin) %>% select(-accepted_bin) %>% as.matrix()
      my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))
      
      # get raoq
      my_raoq <- pairdist(as.matrix(dist(my_tr)))
      
      # get convex hull, but only if we have more unique species than traits #FDRich compute
      if(nrow(unique(my_dt)) > ncol(my_tr)){
        my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
      }
      
    }else{
      my_raoq <- pd_cur <- mpd_cur <- 0
    }
    
    # save results
    out_res <- c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fd" = my_fric, "fdr" = my_fric^(1/K))
    
    return(out_res)
  }
  
  
  # save the results
  outmat_check <- outmat_check %>% data.frame() %>% as_tibble() 
  write_csv(outmat_check, outf)
}
