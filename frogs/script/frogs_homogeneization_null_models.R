#### Maps homogeneization - Null models
#Tali Guez

########################################### LIBRAIRIES ###########################################################################################################
library(tidyverse)
#library(FD) #to compute pd
#library(parallel) #this code uses parallel computing
library(terra)
library(vegan)
library(doParallel)
library(dplyr)
library(geometry)
#library(picante)
#library(ape)
library(foreach)

rm(list=ls())

###################################### Functions ##################################################################################################
#### Functional richness ####

## raoq formula. assumes equal abundances
#pairdist <- function(myd){
 # return(sum((myd / nrow(myd))^2))
#}

calculate_indices <- function(trait, community_matrix) {
  require(geometry)
  #results<-data.frame(nbsp = integer(nrow(community_matrix)),fd = numeric(nrow(community_matrix)),fdr = numeric(nrow(community_matrix)), raoq = numeric(nrow(community_matrix)))
  #results<-list()
  
  #Definition of vectors for results, with communities' names as given in 'community_matrix'
  nbsp<-rep(NA, nrow(community_matrix))
  names(nbsp)<-row.names(community_matrix)

  fd<-rep(NA, nrow(community_matrix))
  names(fd)<-row.names(community_matrix)

  fdr<-rep(NA, nrow(community_matrix))
  names(fdr)<-row.names(community_matrix)

  raoq<-rep(NA, nrow(community_matrix))
  names(raoq)<-row.names(community_matrix)
  
  #Number of traits
  K<-ncol(trait) #-1 if matrix has accepted_bin instead of name of species in row_name (here not useful)
  
  for (i in 1:nrow(community_matrix)) {
    sp_present<-colnames(community_matrix)[community_matrix[i, ] == 1]
    n_sp<-length(sp_present)  #total nb of species
    
    #Initialize variables
    my_raoq<-my_fric<-my_fdr<-NA
    
    if (n_sp<2) {
      #Not enough species
      my_raoq<-my_fric<-my_fdr<-0 #
      
    } else {
      trait_subset<-trait[sp_present, , drop = FALSE] #extract a subset of the trait matrices given the sp_present in the given community
      trait_subset<-trait_subset+matrix(runif(length(trait_subset), -1e-8, 1e-8),nrow=nrow(trait_subset),ncol=ncol(trait_subset)) #add noise
      
      #Compute RaoQ
      dist_mat<-as.matrix(dist(trait_subset))
      my_raoq<-sum((dist_mat/nrow(dist_mat))^2)
      
      #Calculation of fd (volume convexe) #Compute convex hull (FDRich) if enough unique species
      if (nrow(trait_subset) > K) {
        my_fric<-tryCatch(convhulln(as.matrix(trait_subset), "FA")$vol, error = function(e) return(NA))
      }
      
      #Calculation of fdr 
      my_fdr<-ifelse(!is.na(my_fric), my_fric^(1/K), NA)
      
      #results[i,]<-list(nbsp=n_sp,fd=my_fric, fdr=my_fdr, raoq=my_raoq)
      #results[i,]<-c(n_sp,my_fric, my_fdr,my_raoq)
    }
    #Results in the corresponding vectors
    nbsp[i]<-n_sp
    fd[i]<-my_fric
    fdr[i]<-my_fdr
    raoq[i]<-my_raoq
  }
  results<-list(nbsp=nbsp, fd=fd, fdr=fdr, raoq=raoq) #save the vectors in a list (as dbfd from FD package)
  return(results)
}


##In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better <- function(m = 1000, n, sn, trait) {
  # Generate community matrix
  v.names <- paste("random_community", 1:m, sep = "") #row names (generate unique names for each row for each comm)
  random_communities <- matrix(0, nrow = m, ncol = length(sn), dimnames = list(v.names, sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  # Calculate functional diversity for the new matrix
  FD.result<-calculate_indices(trait = trait, community_matrix = random_communities) #calcul functional div for this new matrice
  return(FD.result)
}

########################################### Merge datasets ###########################################################################################################
communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3.txt", header=TRUE)
communities<-as_tibble(communities) %>% 
  select(-(149:150)) 

####OLD DATA
#Change here the combination 1 to 10
indices<-read_csv("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/results/REV_obs_results_Frogs_comb_1.csv")
indices<-indices%>% 
  rename(grilla=grid_id)

#from the previous calculated indices from frogs_homogeneization_calculation.R script
new_comm<-communities %>% 
  left_join(indices, by=join_by(grilla)) %>% 
  column_to_rownames(var="grilla") %>% 
  select(-(148:152))

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

dt_mat_st<-dt_mat %>% 
  column_to_rownames(var="accepted_bin") # %>% 
  #relocate(size_Celio, .before="HD") #for comb1 and 4 and 5
  #relocate(size_Celio, .before="HL") #for comb2 and 3
  #relocate(size_Celio, .before="Hlimb") #for comb6
  

#################################### Calculation of the indices for observed commmunities ###########################################################################
obs_comm<-communities %>% 
  column_to_rownames(var="grilla")

obs_indices<-calculate_indices(trait=dt_mat_st, community_matrix=obs_comm) #calculation of the observed values

#transform the list object into a df
df_obs_indices<-obs_indices %>% 
  map(~ enframe(.x, name = "grilla", value="value")) %>%
  bind_rows(.id="indice") %>%
  pivot_wider(names_from="indice", values_from="value")

df_obs_indices<-df_obs_indices %>%  
  mutate(fd=if_else(nbsp<4, NA_real_, fd),
  fdr=if_else(nbsp<4, NA_real_, fdr),
  raoq=if_else(nbsp<4, NA_real_, raoq)) %>%  #NA in the indices when nbsp<4
  mutate(grilla=as.integer(grilla)) #grilla into integer to match previous dataset

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/results/new")
write.table(df_obs_indices, file = "indices_comb10.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#test1<-read.table("indices_comb1.txt", header=TRUE) #df_obs_indices
#test2<-read.table("indices_comb2.txt", header=TRUE)
#test3<-read.table("indices_comb3.txt", header=TRUE)

new_comm<-communities %>% 
  left_join(df_obs_indices, by=join_by(grilla)) %>% 
  column_to_rownames(var="grilla") %>% 
  #select(-nbsp, -fd, -raoq) #select fdr
  #select(-nbsp, -fdr, -raoq) #select fd
  select(-nbsp, -fd, -fdr)  #select raoq


########################################### Creation of random communities  ###########################################################################################################

species_per_pixel <- unique(apply(new_comm[,1:length(new_comm[1,]) -1],1,sum)) #List with number of unique species observed by pixel.
#1 to 60
species_names<-names(new_comm)[1:(length(new_comm[1,])-2)] #147 species ok

#MAC
#b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )

num_cores <- detectCores() - 1 #Define the number of cores to use
cl <- makeCluster(num_cores) #Create a parallele cluster
clusterExport(cl, varlist = c("give_complete_community_better", "calculate_indices", "species_names", "dt_mat_st"), envir = environment())

b <- parLapply(cl, species_per_pixel, function(x) {
  give_complete_community_better(n = x, sn = species_names, trait = dt_mat_st)
})

stopCluster(cl)


setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models_homogeneization")
#Fdr
#save(b,file="random_hylid10k_FDrich_comb10") #don't need to calculate b each time
#load("random_hylid10k_FDrich")
#Raoq
save(b,file="random_hylid10k_FDdisp")
load("random_hylid10k_FDdisp")
#From before
#setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models")
#load("random_hylid10k_FDrich") #to compare


#Create randomization maps
#First compare observed PD with expected PD
communities1<-new_comm
communities1$sp_number<-apply(new_comm[,1:length(new_comm[1,]) -1],1,sum) #add TD per pixel
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel)){
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"raoq"] #change here 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i]) #go to the right community with number of different species (if i==1, go the expected community of 1 species)
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[4]]>=observed_values[j]))/(1000+1))*2 #change here
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[4]]<=observed_values[j]))/(1000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#b[unlist(cond)][[1]][[3]]  #fdr
#1:nbsp
#2:fd
#3:fdr
#4:raoq

hist(b[unlist(40)][[1]][[3]]) #1 à 4.5
#hist(communities1$fdRic) #0.11.. à 8....

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)
comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 
#1-Create an empty raster using a base to ensure correct extent, size and resolution

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/OR_AUC_10k")
r<-rast("AF_10k_AF_Scinax_alter_OR_AUC_T10.asc")
values(r)<-NA #erase all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3

#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$p_values#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$p_values #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$p_values #3

#3-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models_homogeneization/comb10")
#4- Save rasters to file
#writeRaster(p_value_ras,"P_value_randomization_fdrich_hylids.asc")
#writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdrich_hylids.asc")
#writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich_hylids.asc")

#Fdr
#writeRaster(p_value_ras, "P_value_randomization_fdrich_hylids_comb10_homogene.grd", filetype = "RRASTER", overwrite=TRUE)  
#writeRaster(p_value_sup_ras, "P_value_sup_randomization_fdrich_hylids_comb10_homogene.grd", filetype = "RRASTER", overwrite=TRUE)  
#writeRaster(p_value_inf_ras, "P_value_inf_randomization_fdrich_hylids_comb10_homogene.grd", filetype = "RRASTER", overwrite=TRUE)

#Raoq
writeRaster(p_value_ras, "P_value_randomization_hylids_raoq_comb10.grd", filetype = "RRASTER")  
writeRaster(p_value_sup_ras, "P_value_sup_randomization_hylids_raoq_comb10.grd", filetype = "RRASTER")  
writeRaster(p_value_inf_ras, "P_value_inf_randomization_hylids_raoq_comb10.grd", filetype = "RRASTER")  


#5 Reclassify rasters to only show significant pixels
##assuming alpha of 0.05
##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
m<-c(0,0.05,1,0.05,1,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_sup<-classify(p_value_sup_ras, rclmat, include.lowest=T)
m<-c(0,0.05,-1,0.05,1,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_inf<-classify(p_value_inf_ras, rclmat,include.lowest=T)

significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)

par(mfrow=c(1,1))
plot(significant_pixels)

#writeRaster(significant_pixels,"Significant_pixels.asc")
#fdr
#writeRaster(significant_pixels, "Significant_pixels_comb10_homogene.grd", filetype = "RRASTER", overwrite=TRUE)  ###in tif format
#raoq
writeRaster(significant_pixels, "Significant_pixels_raoq_comb10.grd", filetype = "RRASTER", overwrite=TRUE)  ###in tif format


#setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models_homogeneization/comb2")
#write.table(communities1, file = "fdrich_pvalues_comb2.txt", append = FALSE,row.names=F,quote=F,sep="\t")


