##Frogs maps - homogeneization with trees datasets 
## Tali Guez - 2025

rm(list = ls())
gc() #Garbage collection

#### Librairies ####
library(terra)
library(ape)
library(geiger)
library(picante)
library(ade4)
library(geometry) #For convex hull calculations (FRic)
library(doParallel) #For parallel computing
library(vegan) #For functional diversity analysis
library(picante) #For phylogenetic diversity analysis
library(feather) #Efficient data storage format
library(tidyverse) #Data manipulation and visualization

###### Trait Data ############
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species
summary(trait_filter) #no NA

#species_names <- trait_filter$Species #extract species name 

### By scaled and normalized
#c("size_Celio","HL","HD","TL","Hlimb")#147 entries
comb1<-c("Species","size_Celio","HD","TL") #from the article size_Celio, head_width, tibia_mm
comb2<-c("Species","size_Celio","HL","TL") #size_Celio, head_lenght, tibia_mm 
comb3<-c("Species","size_Celio","HL","Hlimb") #size_Celio, head_lenght,leg_mm
comb4<-c("Species","size_Celio","HD","HL") #size_Celio, head_width, head_lenght
comb5<-c("Species","size_Celio","HD","Hlimb") #size_Celio, head_width, leg_mm
comb6<-c("Species","size_Celio","TL","Hlimb") #size_Celio, tibia_mm, leg_mm
comb7<-c("Species","HD","HL","TL") #head_width, head_lenght, tibia_mm
comb8<-c("Species","HD","HL","Hlimb") #head_width, head_lenght, leg_mm
comb9<-c("Species","HL","TL","Hlimb") #head_lenght, tibia_mm, leg_mm
comb10<-c("Species","HD","TL","Hlimb") #head_width, tibia_mm, leg_mm

#scale them
dt_mat <- trait_filter %>%
  select(all_of(comb10)) %>%  #change combination here (1 to 10)
  #mutate(across(where(is.numeric), scale)) #gives the same values
  gather(PC, value, -Species) %>% 
  group_by(PC) %>%
  mutate(value = (value - mean(value))/sd(value)) %>%
  ungroup %>%
  spread(PC, value) %>% 
  rename(accepted_bin=Species)

###### Community matrix ########
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")
a<-read.table(file="communities_hylids10k_V3.txt", header=TRUE)
a[is.na(a)]<-0

a_comm<-as_tibble(a) %>% 
  rename(grid_id=grilla) %>% 
  pivot_longer(cols=2:148, names_to = "accepted_bin", values_to = "present") %>% 
  relocate(present, .before = accepted_bin) 
comm0<-a_comm%>% 
  filter(present == 1)

#See the mistmatches
sum(!dt_mat$accepted_bin%in%comm0$accepted_bin) #0
unique(comm0$accepted_bin[!comm0$accepted_bin%in%dt_mat$accepted_bin]) #0

###### Pruned matrix #########
#setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned")
#pruned_tree<-read.tree("pruned_tree_frogs.csv") #Pruned phylogenetic tree

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas")
##read phylogeny
user_phylogeny<-ape::read.tree("raw/PD_data/Jetz_Pyron_amph_shl_new_Consensus_7238.tre")
#sp_names<-read.csv("cleaned/species_names_full.csv")
#species_names_full
pruned_tree<-user_phylogeny

####Make sure the names match between maps and phylogeny
setdiff(user_phylogeny$tip.label, sp_names$species_names_full)
setdiff(user_phylogeny$tip.label, comm0$accepted_bin)
setdiff(comm0$accepted_bin, user_phylogeny$tip.label)

##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
#pruned_tree<-ape::drop.tip(user_phylogeny, setdiff(user_phylogeny$tip.label, sp_names$species_names_full))
#pruned_tree<-ape::drop.tip(user_phylogeny, setdiff(user_phylogeny$tip.label, comm0$accepted_bin))

#pruned_tree$tip.label <- gsub("_", " ", pruned_tree$tip.label)

#setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned")
#write.tree(pruned_tree,"pruned_tree_frogs.csv")
##check that all disrtribution data is in tree
#test<-as.data.frame(species_names_full)
#rownames(sp_names)<-sp_names$species_names_full
#test<-as.data.frame(comm0)
#rownames(test)<-test$accepted_bin

#check_names<-geiger::name.check(pruned_tree, sp_names, data.names=NULL) #this must be OK
#check_names<-geiger::name.check(pruned_tree,comm0$accepted_bin, data.names=NULL) #this must be OK

#### Run same script as trees #####
# register the parallel cores
registerDoParallel(round(detectCores()/2)) #Uses half of the available cores
getDoParWorkers()

## raoq formula. assumes equal abundances
pairdist <- function(myd){
  return(sum((myd / nrow(myd))^2))
}


save_suffix <- "comb_10" #1 à 10 #to change

##analyzing frogs
fit_seq <- c(2)  
set.seed(555)

for(gg in fit_seq){
  
  print("Fitting Frogs")
  comm <- comm0 #%>% left_join(ag, by = "accepted_bin") %>% filter(group == "Angiosperms") #no need because frogs
  outf <- paste0("cleaned/results/REV_obs_results_Frogs_", save_suffix, ".csv") #à changer pour adapter aux traits
  
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
    
    # Get the focal plot
    my_dt <- comm %>% filter(grid_id == plts[j])
    
    # Initialize variables
    my_dt_sub <- NULL
    pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA ####
    out_res <- c()
    
    if(nrow(my_dt) > 1){
      
      # Subset the data and get the phylometrics
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
comb1_results<-read_csv("cleaned/results/REV_obs_results_Frogs_comb_1.csv") #change names for the rest
summary(comb1_results)
results<-read_csv("cleaned/results/REV_obs_results_Frogs_comb_3.csv") #change names for the rest
summary(results)

lignes_a_zero <- a %>%
  rowwise() %>%  #sum each line
  mutate(somme = sum(c_across(-grilla))) %>%  #sum each column without "grilla"
  filter(somme == 0) %>%  #keep line where sum=0
  ungroup()  #
#301 pixels with no species
#a : 19293 entries
#comb1_results : 18992 entries (301 entries/pixels with no species missing)
write_csv(lignes_a_zero, "cleaned/pixels_with_no_sp.csv")
