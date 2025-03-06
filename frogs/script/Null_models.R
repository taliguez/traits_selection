################################################################################################################################################## My script
####################################### SCRIPT - Internship M2 ##########################################################################################
#########################################  Tali - Guez ############################################################################################
##################################################################################################################################################My script
#Amphibian 10k
#NULL models calculation

########################################### LIBRAIRIES ###########################################################################################################
#### For Null models
library(FD) #to compute pd
library(parallel) #this code uses parallel computing
library(terra)
library(tidyverse)

###################################### NULL models ################################################
#### Functional richness ####

###Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina 

rm(list=ls())
##PD randomizations for each community

####Functions

#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=1000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names (generate unique names for each row for each comm)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F) #calcul functional div for this new matrice
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}

####
##load and process trait database, this must march community data
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
#merged_trait_all<-merge(new_trait, trait, by.x = "X", by.y = "Species", all.x=TRUE, all.y=TRUE) ###keep all lines
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species
#View(trait_filter)

#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
trait1<-trait_filter

rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(2,4,5)]
for (i in 1:3){
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community composition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
#communities<-read.table("communities_hylids_and_fd1k.txt",h=T,row.names=1)#1k
#communities<-read.table("communities_hylids_and_fd10k_cont.txt",h=T,row.names=1) #10k
communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3.txt", header=TRUE)

#communities<-communities[,c(1:158,160)] ##leaving only species and Fdric
#communities<-communities[,1:159] ##for 1 k matrix is different

rownames(communities)<-communities[,1] 

communities<-communities[,c(2:148,150)] ##leaving only species and Fdric (take off fdDisp)

#names(communities)[159]<-"fdRic"##for 1k only
species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum)) #List with number of unique species observed by pixel.
#1 to 60
species_names<-names(communities)[1:(length(communities[1,])-1)] #147 species ok

#MAC
#b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )
#windows doesn t work
#b<-parLapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1)}, mc.cores=16 )

num_cores <- detectCores() - 1 #Define the number of cores to use
cl <- makeCluster(num_cores) #Create a parallele cluster
clusterExport(cl, varlist = c("give_complete_community_better", "species_names", "trait1"), envir = environment()) #Export the necessary object in the cluster

#Execute the function in parallele
b <- parLapply(cl, species_per_pixel, function(x) {
  give_complete_community_better(n = x, sn = species_names, trait = trait1)
})

# Arrêter le cluster
stopCluster(cl)

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models")

save(b,file="random_hylid10k_FDrich")
load("random_hylid10k_FDrich")
#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum) #add TD per pixel
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel)){
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdRic"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i]) #go to the right community with number of different species (if i==1, go the expected community of 1 species)
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[3]]>=observed_values[j]))/(1000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[3]]<=observed_values[j]))/(1000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}


#hist(b[unlist(cond)][[1]][[3]]) #1 à 5
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

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models")
#3- Save rasters to file
#writeRaster(p_value_ras,"P_value_randomization_fdrich_hylids.asc")
#writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdrich_hylids.asc")
#writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich_hylids.asc")

writeRaster(p_value_ras, "P_value_randomization_fdrich_hylidsV2.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_sup_ras, "P_value_sup_randomization_fdrich_hylidsV2.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_inf_ras, "P_value_inf_randomization_fdrich_hylidsV2.tif", filetype = "GTiff")  ###in tif format


#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

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
#writeRaster(significant_pixels,"Significant_pixels.asc")
writeRaster(significant_pixels, "Significant_pixelsV2.tif", filetype = "GTiff")  ###in tif format


plot(significant_pixels)


setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/table_comparaison")
write.table(communities1, file = "fdrich_pvalues", append = FALSE,row.names=F,quote=F,sep="\t")

############################# Frichness - other Comb #######################################################
###### Change only with the right combinations
rm(list=ls())
##PD randomizations for each community

####Functions
#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=1000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names (generate unique names for each row for each comm)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F) #calcul functional div for this new matrice
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}

####
##load and process trait database, this must march community data
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species

#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
#trait_Comb<-trait_filter 
#rownames(trait_Comb)<-trait_Comb[,1] #147 entries
#trait_Comb2<-trait_Comb[,c(2,3,5)] #size_Celio, head_lenght, tibia_mm 
#trait_Comb3<-trait_Comb[,c(2,3,6)] #size_Celio, head_lenght,leg_mm
#trait_Comb4<-trait_Comb[,c(2,3,4)] #size_Celio, head_width, head_lenght
#trait_Comb5<-trait_Comb[,c(2,4,6)] #size_Celio, head_width, leg_mm
#trait_Comb6<-trait_Comb[,c(2,5,6)] #size_Celio, tibia_mm, leg_mm
#trait_Comb7<-trait_Comb[,c(3,4,5)] #head_width, head_lenght, tibia_mm
#trait_Comb8<-trait_Comb[,c(3,4,6)] #head_width, head_lenght, leg_mm
#trait_Comb9<-trait_Comb[,c(3,5,6)] #head_lenght, tibia_mm, leg_mm
#trait_Comb10<-trait_Comb[,c(4,5,6)] #head_width, tibia_mm, leg_mm

trait1<-trait_filter

rownames(trait1)<-trait1[,1]
#Change traits !!!
trait1<-trait1[,c(4,5,6)]
for (i in 1:3){
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community composition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
#Change here the number !!!
communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3_comb10.txt", header=TRUE)
rownames(communities)<-communities[,1] 
communities<-communities[,c(2:148,150)] ##leaving only species and Fdric (take off fdDisp)

species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum)) #List with number of unique species observed by pixel.
#1 to 60
species_names<-names(communities)[1:(length(communities[1,])-1)] #147 species ok

#Parallel
num_cores <- detectCores() - 1 #Define the number of cores to use
cl <- makeCluster(num_cores) #Create a parallele cluster
clusterExport(cl, varlist = c("give_complete_community_better", "species_names", "trait1"), envir = environment()) #Export the necessary object in the cluster

#Execute the function in parallele
b <- parLapply(cl, species_per_pixel, function(x) {
  give_complete_community_better(n = x, sn = species_names, trait = trait1)
})

# Arrêter le cluster
stopCluster(cl)

#change number here !!!
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/Comb10") 

save(b,file="random_hylid10k_FDrich")
#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum) #add TD per pixel
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel)){
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdRic"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[3]]>=observed_values[j]))/(1000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[3]]<=observed_values[j]))/(1000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

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

#Change Number here !!!
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/Comb10")
#3- Save rasters to file ### CHANGE NUMBER HERE!!!
writeRaster(p_value_ras, "P_value_randomization_fdrich_hylidsV2_comb10.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_sup_ras, "P_value_sup_randomization_fdrich_hylidsV2_comb10.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_inf_ras, "P_value_inf_randomization_fdrich_hylidsV2_comb10.tif", filetype = "GTiff")  ###in tif format

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

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
writeRaster(significant_pixels, "Significant_pixelsV2_comb10.tif", filetype = "GTiff")  ###in tif format
#CHANGE NUMBER HERE !!! 

plot(significant_pixels)


############################### Fdispersion ##################################################################

####Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina ##

rm(list=ls())
##PD randomizations for each community
####Functions
give_complete_community_better= function(m=1000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names (generate unique names for each row for each comm)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F) #calcul functional div for this new matrice
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}

####
##load and process trait database, this must march community data
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
#merged_trait_all<-merge(new_trait, trait, by.x = "X", by.y = "Species", all.x=TRUE, all.y=TRUE) ###keep all lines
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species
#View(trait_filter)

#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
trait1<-trait_filter

rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(2,4,5)]
for (i in 1:3){
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
#communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_pd10k_V3.txt", header=TRUE)
communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3.txt", header=TRUE)
rownames(communities)<-communities[,1] 

communities<-communities[,c(2:149)] ##leaving only species and fdDisp (take off Fdric)

species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum)) #List with number of unique species observed by pixel.
species_names<-names(communities)[1:(length(communities[1,])-1)] #147

#MAC
#b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )

#Windows
num_cores <- detectCores() - 1 #Define the number of cores to use
cl <- makeCluster(num_cores) #Create a parallele cluster
clusterExport(cl, varlist = c("give_complete_community_better", "species_names", "trait1"), envir = environment()) #Export the necessary object in the cluster

#Execute the function in parallele
b <- parLapply(cl, species_per_pixel, function(x) {
  give_complete_community_better(n = x, sn = species_names, trait = trait1)
})

# Arrêter le cluster
stopCluster(cl)

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models")

save(b,file="random_hylid10k_FDdisp")
load("random_hylid10k_FDdisp")

#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum) #add TD per pixel
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA
#communities1<-communities1 %>% 
  #filter(sp_number>=3) #filter 


for (i in 1:length(species_per_pixel)){
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdDis"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[6]]>=observed_values[j]))/(1000+1))*2 #1000 instead of 10000 + 6 instead of 3 to have disp
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[6]]<=observed_values[j]))/(1000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}


#hist(b[unlist(cond)][[1]][[6]])
hist(b[unlist(40)][[1]][[6]])
#hist(communities1$fdDis)

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

communities1 <- communities1 %>%
  mutate(
    p_values = if_else(sp_number < 3, NA_real_, p_values),
    p_values_lower = if_else(sp_number < 3, NA_real_, p_values_lower),
    p_values_higher = if_else(sp_number < 3, NA_real_, p_values_higher),
    desviacion = if_else(sp_number < 3, NA_character_, desviacion)
  )


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

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models")
#3- Save rasters to file
#writeRaster(p_value_ras,"P_value_randomization_fdrich_hylids.asc")
#writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdrich_hylids.asc")
#writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich_hylids.asc")

writeRaster(p_value_ras, "P_value_randomization_fdrich_hylids_fDdispV4.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_sup_ras, "P_value_sup_randomization_fdrich_hylids_fDdispV4.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_inf_ras, "P_value_inf_randomization_fdrich_hylids_fDdispV4.tif", filetype = "GTiff")  ###in tif format
#V2, V3 upgraded version
#V4, good version

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#5 Reclassify rasters to only show significant pixels
##assuming alpha of 0.05
##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
m<-c(0,0.05,1,0.05,1,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_sup<-classify(p_value_sup_ras, rclmat, include.lowest=T)
m<-c(0,0.05,-1,0.05,1,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_inf<-classify(p_value_inf_ras, rclmat,include.lowest=T)

#5 Reclassify rasters to only show significant pixels
##assuming alfa of 0.05
##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
#m<-c(0,0.05,1,0.05,2,NA)
#rclmat<-matrix(m,ncol=3,byrow=TRUE)
#p_value_sup<-classify(p_value_sup_ras, rclmat, include.lowest=T)
#m<-c(0,0.05,-1,0.05,2,NA)
#rclmat<-matrix(m,ncol=3,byrow=TRUE)
#p_value_inf<-classify(p_value_inf_ras, rclmat,include.lowest=T)

significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
#writeRaster(significant_pixels,"Significant_pixels.asc")
writeRaster(significant_pixels, "Significant_pixels_fDdispV4.tif", filetype = "GTiff", overwrite=TRUE)  ###in tif format
#V3 wrong

par(mfrow=c(1,1))
plot(significant_pixels)

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/table_comparaison")
write.table(communities1, file = "fdisp_pvalues_V4", append = FALSE,row.names=F,quote=F,sep="\t")



############################# Fdispersion - other Comb #######################################################
###### Change only with the right combinations
rm(list=ls())
##PD randomizations for each community

####Functions
#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=1000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names (generate unique names for each row for each comm)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F) #calcul functional div for this new matrice
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}

give_complete_community_better= function(m=1000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names (generate unique names for each row for each comm)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn)) #generate random matrices
  for (k in 1:m) { #fill the matrices with random species
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F, calc.FRic=F) #calcul functional div for this new matrice
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
} 
#pour comb8 parce que sinon mets trop de temps à converger sur windows et ca marche plus

####
##load and process trait database, this must march community data
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species

#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
#trait_Comb<-trait_filter 
#rownames(trait_Comb)<-trait_Comb[,1] #147 entries
#trait_Comb2<-trait_Comb[,c(2,3,5)] #size_Celio, head_lenght, tibia_mm 
#trait_Comb3<-trait_Comb[,c(2,3,6)] #size_Celio, head_lenght,leg_mm
#trait_Comb4<-trait_Comb[,c(2,3,4)] #size_Celio, head_width, head_lenght
#trait_Comb5<-trait_Comb[,c(2,4,6)] #size_Celio, head_width, leg_mm
#trait_Comb6<-trait_Comb[,c(2,5,6)] #size_Celio, tibia_mm, leg_mm
#trait_Comb7<-trait_Comb[,c(3,4,5)] #head_width, head_lenght, tibia_mm
#trait_Comb8<-trait_Comb[,c(3,4,6)] #head_width, head_lenght, leg_mm
#trait_Comb9<-trait_Comb[,c(3,5,6)] #head_lenght, tibia_mm, leg_mm
#trait_Comb10<-trait_Comb[,c(4,5,6)] #head_width, tibia_mm, leg_mm

trait1<-trait_filter

rownames(trait1)<-trait1[,1]
#Change traits !!!
trait1<-trait1[,c(4,5,6)]
for (i in 1:3){
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community composition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
#Change here the number !!!
communities<-read.table("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3/communities_hylids_and_fd10k_cont_bis_V3_comb10.txt", header=TRUE)
rownames(communities)<-communities[,1] 
communities<-communities[,c(2:149)] ##leaving only species and fdDisp (take off Fdric)

species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum)) #List with number of unique species observed by pixel.
#1 to 60
species_names<-names(communities)[1:(length(communities[1,])-1)] #147 species ok

#Parallel
num_cores <- detectCores() - 1 #Define the number of cores to use
cl <- makeCluster(num_cores) #Create a parallel cluster
clusterExport(cl, varlist = c("give_complete_community_better", "species_names", "trait1"), envir = environment()) #Export the necessary object in the cluster

#Execute the function in parallel
b <- parLapply(cl, species_per_pixel, function(x) {
  give_complete_community_better(n = x, sn = species_names, trait = trait1)
})

# Arrêter le cluster
stopCluster(cl)

#change number here !!!
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/Comb10") 

save(b,file="random_hylid10k_FDdisp") 
#load("random_hylid10k_FDdisp")
#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum) #add TD per pixel
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel)){
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdDis"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[6]]>=observed_values[j]))/(1000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[6]]<=observed_values[j]))/(1000+1))*2 #6 si calcul tout sinon mettre 4 (comme comb8)
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

communities1 <- communities1 %>%
  mutate(
    p_values = if_else(sp_number < 3, NA_real_, p_values),
    p_values_lower = if_else(sp_number < 3, NA_real_, p_values_lower),
    p_values_higher = if_else(sp_number < 3, NA_real_, p_values_higher),
    desviacion = if_else(sp_number < 3, NA_character_, desviacion)
  )

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

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#Change Number here !!!
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned/null_models/Comb10")
#3- Save rasters to file ### CHANGE NUMBER HERE!!!
writeRaster(p_value_ras, "P_value_randomization_fdrich_hylids_fDdisp_comb10.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_sup_ras, "P_value_sup_randomization_fdrich_hylids_fDdisp_comb10.tif", filetype = "GTiff")  ###in tif format
writeRaster(p_value_inf_ras, "P_value_inf_randomization_fdrich_hylids_fDdisp_comb10.tif", filetype = "GTiff")  ###in tif format

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

writeRaster(significant_pixels, "Significant_pixels_fDdisp_comb10.tif", filetype = "GTiff")  ###in tif format
#CHANGE NUMBER HERE !!! 


