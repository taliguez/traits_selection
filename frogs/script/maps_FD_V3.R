################################################################################################################################################## My script
####################################### SCRIPT - Internship M2 ##########################################################################################
#########################################  Tali - Guez ############################################################################################
##################################################################################################################################################My script
#Amphibian maps - 10k
#V3 - Nouveau tableau + merge + cleaned

########################################### LIBRAIRIES ###########################################################################################################

### Shapefile (replace rgdal and raster)
#install.packages("terra")
library(terra)

###Phylogenetic diversity 
library(ape)
library(geiger)
library(picante)

# Function to calculate the four Functional diversity indices (So far we only use the FRic index since no abundance data is available from distribution maps)
library(geometry)
library(ape)
library(ade4)
#install.packages('FD')
library(FD)

########################################### Load map - atlantic forest ###########################################################################################################

###load mask to use (AF) ##change path
#selected_mask <- vect(""/Users/andreapaz/Dropbox/Guano's_files/AF_shapefile_from_Guano", "atlantic_forest")
selected_mask <- vect("C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/frogs/datas/raw/AF_shapefile_from_Guano/AF_shapefile_from_Guano", "atlantic_forest")
selected_mask <- aggregate(selected_mask, dissolve = TRUE) #aggregate and dissolve polygons of a spatial object (selected_mask) into a single polygon
print(selected_mask)
plot(selected_mask)

########################################### Function to rasterize shapefile into ascii format ###########################################################################################################
####First stack rasters with points ####
### function to rasterize points

#Updated function to rasterize points from shapefile
rasterize_species2 <- function(x, mask = selected_mask, output_dir = getwd()) {
  #r<-raster(ncol=300,nrow=400,resolution=resolution,ext=extent(selected_mask))
  #res(r)<-resolution 
  #r<-extend(r,selected_mask)
  #r<-crop(r,extent(selected_mask))
  #r<-mask(r,selected_mask)
  values(r)<-0 #raster initialization
  map<-vect(file.path(maps_folder, paste0(x, ".shp"))) #lire shapefile terra package: shapefile is in the maps_folder
  r<-rasterize(x = map, y = r, field = 1, update = TRUE, background = 0) # Rasterization
  r<-mask(r, mask) #Mask values in a SpatRaster
  valor <- unique(values(r)) #Get the unique cell values of a SpatRaster
  
  if (length(valor) == 1 && is.na(valor)) {
    return(NULL) # Verification
  } else if (length(valor) == 2 && valor[2] == 0) {
    return(NULL)
  } else {
    if (!dir.exists(output_dir)) { #put the file in new directory
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, paste0(x, ".asc"))
    writeRaster(x = r, output_file, filetype = "AAIGrid", overwrite = TRUE)
    return(rast(output_file))
  }
}

#### INITIALIZATION 
#resolution<-0.08333333 ##10km
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/points_14model")
maps_folder <- "C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/frogs/datas/points_14model"

distribution_files<-list.files(path=getwd(), pattern= "*.shp$") #list every file in the actual maps_folder with *.shp$ motifs
species_names_points<-sub(".shp","",distribution_files) ##really file names # extract species name from the files

output_dir <- "C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/frogs/datas/points_14model/points_14model_test" #output directory

##Function to do all 14 files
rasters<-lapply(species_names_points, function(species) {
  rasterize_species2(x = species, mask = selected_mask, output_dir = output_dir)
})
#Example: rasterize_species2(x = "Points_Scinax_uruguayus", mask = selected_mask, output_dir = output_dir)

############################ MAPS - PREPARATION OF RASTERS FOR ANALYSIS ###########################################################################################################

###load SDMs and make correct size
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/") ##to folder with SDMs
#setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap\\OR_AUC_T10/")
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/OR_AUC_10k") ##to folder with SDMs

##Create a new loop to not overwrite the existing file => new directory : output_rasters
#Define the output directory
output_dir <- "output_rasters"

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#Read distribution files and models
#models<-list.files(pattern="AF_")
#models<-list.files() #everything in working directory
models<-list.files(pattern="AF_.*\\.asc$") # List only the .asc files

#The loop is not useful here I guess
for(i in 1:length(models)){ 
  r <- rast(models[i]) #transform into a raster
  #par(mfrow=c(1,4))
  #plot(r)
  
  # Extend, crop, mask, and aggregate the raster i
  r <- extend(r, selected_mask)
  #plot(r)
  r <- crop(r, selected_mask)
  #plot(r)
  r <- mask(r, selected_mask)
  #plot(r)
  r <- aggregate(r, 10, fun=max) #fact=10
  
  # Define the output file path
  output_file <- file.path(output_dir, paste("AF_10k_", models[i], sep=""))
  #output_file <- file.path(output_dir, paste(models[i], sep=""))
  # Write the raster to an ASCII file in the new directory
  writeRaster(r, output_file, filetype="AAIGrid", overwrite=TRUE)
}
#print(models) #104 files 

##when everything is created just:
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/")
#setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap\\OR_AUC_T10/")
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/OR_AUC_10k")

# Load models
#models <- list.files(pattern = "AF_10k_")
models<-list.files(pattern = "AF_10k_.*\\.asc$") #from the working directory

#To Stack
#use c(x1, x2, x3) when providing separate arguments
#use rast(list(x1,x2,x3)) when providing the arguments as a list.

# Load each raster file individually and combine them into a single SpatRaster
raster_list<-lapply(models, rast) #transform models into list

# Combine the rasters into a single SpatRaster
Stack_models<-rast(raster_list) #create a SpatRaster with 3 layers
plot(Stack_models)
#names(Stack_models)

#To clean the files name
species_names_models<-sub("_OR_AUC_T10", "", models) #take _OR_AUC_T10 off
species_names_models<-sub("AF_10k_AF_", "", species_names_models) #take AF_10k_AF_ off
#species_names_models<-sub("AF_10k_AF_10k_AF_", "", species_names_models) #take AF_10k_AF_10k_AF_ off
species_names_models<-sub(".asc", "", species_names_models, fixed = TRUE) #take off .asc --> only have the species name
#104 species name

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/Points_only_10k_and_14")

# Load the points
points<-list.files(pattern = ".asc")
Stack_points<-rast(points) #into 1 spatRaster
plot(Stack_points)
species_names_points<-sub(".asc", "", points, fixed = TRUE) #take off .asc 
species_names_points<-sub("Points_", "", species_names_points) #take off Points_ 
#54 species name
species_names_full<-c(species_names_points, species_names_models) #have every species names 
#158 names of species in total


############################ MAPS - TAXONOMIC DIVERSITY - NEW datasets ###########################################################################################################

#### To get richness and community composition ####
Stack_full <- c(Stack_points, Stack_models) 
plot(Stack_full)
print(global(Stack_full, "sum", na.rm = TRUE))

#### Select species with all information from table
#Trait data: select file containing trait data
trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
new_trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/imputed_data_BM.csv", header=TRUE, sep=";")

###merge both dataset
#merged_trait_all<-merge(new_trait, trait, by.x = "X", by.y = "Species", all.x=TRUE, all.y=TRUE) ###keep all lines
merged_trait <- merge(new_trait, trait, by.x = "X", by.y = "Species")
list_col<-c("X","size_Celio","HL","HD","TL","Hlimb")#147 entries

trait_filter<-merged_trait[list_col] #147 entries
colnames(trait_filter)[1] = "Species" #rename col X into Species
summary(trait_filter) #no NA

species_names <- trait_filter$Species #extract species name 

#Stack with clean names
tabla2 <- as.data.frame(species_names_full)
colnames(tabla2) <- "Grilla" #nom colonne = Grilla

#Change names of models tu species names and add the empty raster in the beggining (pixel number)
Stack_full4<-Stack_full #create new stack : c(Stack_points,Stack_models) 
names(Stack_full4) #To check
names(Stack_full4) <- as.vector(tabla2$Grilla) #clean names from this stack
names(Stack_full4) #To check again 

#Filter raster
Stack_filtre <- Stack_full4[[names(Stack_full4) %in% species_names]]
names(Stack_filtre) #To check (147 species ok)

## To get Taxonomic diversity
plot(Stack_filtre)
print(global(Stack_filtre, "sum", na.rm = TRUE))

TD_filtre<-app(Stack_filtre, sum, na.rm = TRUE) #Sum species of each layers
#app= apply a function to the values of each cell of a SpatRaster
plot(TD_filtre)

TD2_filtre<-crop(TD_filtre, selected_mask)
plot(TD2_filtre)

TD3_filtre<-mask(TD2_filtre, selected_mask)
plot(TD3_filtre, main="TD filtered") #147 species in total 

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")
writeRaster(TD3_filtre, "richness_hylids_complete_0_083_V3.tif", filetype = "GTiff")  ###in tif format
#writeRaster(TD3_filtre, "richness_hylids_all.grd", filetype = "RRASTER")  ###in tif format

par(mfrow=c(1,1)) #to compare
plot(TD3, main="Non-filtered")
plot(TD3_filtre, main="Filtered")

############################ MAPS - PHYLOGENETIC DIVERSITY - With new datasets ###########################################################################################################

####Load Phylogeny, generate list of species
##To read a phylogeny in newick format use read.tree instead of read.nexus
#In working directory always a file with the same name: phylogeny.nex if not then change name here
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/PD_data/")
#setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap")
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/PD_data")

##read phylogeny
user_phylogeny<-ape::read.tree("Jetz_Pyron_amph_shl_new_Consensus_7238.tre")

####MMake sure the names match between maps and phylogeny
setdiff(user_phylogeny$tip.label, species_names_full)
##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
pruned.tree<-ape::drop.tip(user_phylogeny, setdiff(user_phylogeny$tip.label, species_names_full))
##check that all disrtribution data is in tree
test<-as.data.frame(species_names_full)
rownames(test)<-species_names_full
check_names<-geiger::name.check(pruned.tree, test, data.names=NULL) #this must be OK

#########################.
####START ANALYSES#######
#########################.

#Stack with clean names
tabla2 <- as.data.frame(species_names_full)
colnames(tabla2) <- "Grilla" #nom colonne = Grilla

#Change names of models tu species names and add the empty raster in the beggining (pixel number)
Stack_full5<-c(Stack_points,Stack_models) #create new stack ##To get richness and community composition
#names(Stack_full5) #To check
names(Stack_full5) <- as.vector(tabla2$Grilla) #clean names from this stack
#names(Stack_full5) #To check again 

#Create empty raster of the desired area and resolution to assign pixel numbers
r2_filter <- Stack_full5[[1]]
res(r2_filter)<-res(Stack_full5) #resolution

grilla2_filter <- r2_filter
names(grilla2_filter) <- "grilla"
grilla2_filter[1:ncell(grilla2_filter)] <- 1:ncell(grilla2_filter)
lista_grilla2_filter <- list(grilla2_filter)
names(lista_grilla2_filter) <- "Grilla"

#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_full5) <- as.vector(tabla2$Grilla)

#Filter raster
Stack_filtre_1 <- Stack_full5[[names(Stack_full5) %in% species_names]]
names(Stack_filtre_1) #To check 147 species

Stack2_filter <- c(lista_grilla2_filter$Grilla, Stack_filtre_1) ###Full stack of all maps

#Turn maps into dataframe for computation of PD
community_data_filter <- as.data.frame(Stack2_filter)
##remove rows that are NA 
community_data1_filter<-community_data_filter[-which(rowSums(!is.na(community_data_filter[, 2:109])) == 0),]
species_names_filter<-colnames(community_data1_filter)[2:length(community_data1_filter)]  #Store species names

#Store 
setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")
write.table(community_data1_filter, file = "communities_hylids10k_V3.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara

community_data1_filter[is.na(community_data1_filter)]<-0
#In the community data frame NA must be eliminated done before is this if you load community data and not maps?

#III-Phylogenetic diversity computation - computes only Faith's PD others may be added
pd.result_filter <-picante::pd(community_data1_filter[,2:ncol(community_data1_filter)],pruned.tree,include.root = F) 
##if it fails switch to #PhyloMeasures::pd.query(tree, matrix,)

#Add the pixel PD value to data frame
community_data1_filter$pd<-pd.result_filter[[1]]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data1_filter, file = "communities_hylids_and_pd10k_V3.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area
values(r2_filter)<-0
pd_ras_filter<-r2_filter
values(pd_ras_filter)<-NA #Eliminate every value they will be replaced by PD values further down

#2- Assign PD values to raster
pd_ras_filter[community_data1_filter$grilla]<-community_data1_filter$pd

#3- Save raster to file 
#writeRaster(pd_ras,"hylids_PD10k.asc",format="ascii")
writeRaster(pd_ras_filter, "hylids_PD10k_V3.tif", filetype = "GTiff")  ###in tif format

setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/cleaned")
#Write the new matrix to a file to avoid rerunning the script for potential further analyses
df_pd.result_filter<-pd.result_filter
df_pd.result_filter$grilla<-rownames(pd.result_filter)
write.table(df_pd.result_filter, file = "PD_TD_V4.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
plot(pd_ras_filter)

#par(mar = c(4, 4, 2, 2))
#par(mfrow=c(1,2))
plot(pd_ras, main="PD all")
plot(pd_ras_filter, main="PD filtered")

################### MAPS - FD + DISPERSION - new datasets - Comparative models ###########################################################################################################

### Author: Sebastien Villéger, adapted by Claire Fortunel (Please acknowledge as appropriate)  

# Notations corresponds with Villéger et al. (2008) Ecology, 89: 2290-2301 for FRic, FEve, FDiv; and Bellwood et al. (2006) Proc. R. Soc. B., 273: 101-107 for FSpe

##Load maps exactly as for PD if not yet loaded
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/")
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/OR_AUC_10k")
models<-list.files(pattern="AF_") 
#Stack_models<-rast(raster_list) #create a SpatRaster with 3 layers #Fait avant
#plot(Stack_models)
species_names_models<-sub("_OR_AUC_T10","",models) 
species_names_models<-sub("AF_10k_AF_","",species_names_models)
species_names_models<-sub(".asc","",species_names_models,fixed=T) #To have all the species names

#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/Points_only/") 
setwd("~/Master/M2/Internship_M2/analyse/frogs/datas/Points_only_10k_and_14")
points<-list.files(pattern=".asc")
#Stack_points <- rast(points) #fait avant
species_names_points<-sub(".asc","",points,fixed=T) 
species_names_points<-sub("Points_","",species_names_points) #fait avant 
species_names_full<-c(species_names_points,species_names_models) #avoir tous les noms #158 

##To get richness and community composition
Stack_full2<-c(Stack_points,Stack_models) 
plot(Stack_full2)
print(global(Stack_full2, "sum", na.rm = TRUE))

##polygon with mask is called selected_maks from above

####Communities####

#Trait data: select file containing trait data
#trait<-read.csv("~/Master/M2/Internship_M2/analyse/frogs/datas/raw/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
#trait<-read.csv("Species_traits_imputed_lilian.csv")

###Trim trait data to match distribution data #2 don t match the ocean ones
trait_filter2<-trait_filter[trait_filter$Species %in% species_names_full,] 
### Trim the distribution data to match trait data, This is not needed in this case
#species_names_full<-species_names_full[species_names_full %in% trait1$Species]

#tabla<-as.data.frame(species_names_full) #table with names of all species (158)
#colnames(tabla)<-"Grilla" #nom col= Grilla

#Load all distribution maps and make them same extent as will be used
#r1<-Stack_full2[[1]]
#grilla=r1 
#names(grilla)="grilla" #
#grilla[1:ncell(grilla)]<-1:ncell(grilla)
#lista_grilla<-list(grilla)
#names(lista_grilla)<-"Grilla"
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
#names(Stack_full2)<-as.vector(tabla$Grilla)
#Stack<-c(lista_grilla$Grilla,Stack_full2) ###Full stack of all maps

Stack<-Stack2_filter

#Turn maps into dataframe for computation of FD
marco<-as.data.frame(Stack)
marco<-marco[-which(rowSums(!is.na(marco[,2:148]))==0),] ## does not work yet #remove all values in col 2 to 159 when only NA in it
marco[is.na(marco)] = 0 #replace every NA by 0
species_names<-colnames(marco)[2:length(marco)] #Store species names
marco1<-marco[-which(rowSums(marco[,2:148])==0),] #remove rows where the sum of it is equal to 0

ordered_columns <- order(names(marco1)[-1]) #alphabetical order of species (col)
marco1_bis <- marco1[, c(1, ordered_columns + 1)] #reorganize dataframe but keep grilla at first place
#148 col/Species

###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
trait_filter2<-trait_filter
rownames(trait_filter2)<-trait_filter2[,1]
trait_filter2<-trait_filter2[,c(2,4,5)] #from the article size_Celio, head_width, tibia_mm

for (i in 1:3){
  trait_filter2[,i]<-log(trait_filter2[,i]) ## This could be updated to another mathematic transformation if needed
}

rownames(trait_filter2)==colnames(marco1_bis[,2:148]) #if false: not the same
FDindices<-FD::dbFD(trait_filter2,marco1_bis[,2:148],calc.FDiv=F) 

fddis_ras<-r2_filter
fdric_ras<-r2_filter
values(fddis_ras)<-NA # Erase all values from the distribution map
values(fdric_ras)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco4<-marco1_bis
marco4$fdDis<-FDindices$FDis
marco4$fdRic<-FDindices$FRic

fddis_ras[marco4$grilla]<-marco4$fdDis
fdric_ras[marco4$grilla]<-marco4$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_ras, "hylid_functional_dispersion10k_cont_V3.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_ras, "hylid_functional_richness10k_cont_V3.tif", filetype = "GTiff")  ###in tif format
#should be in datas/cleaned
write.table(marco1_bis, file = "communities_hylids_and_fd10k_cont_V3.txt", append = FALSE,row.names=F,quote=F,sep="\t")
write.table(marco4, file = "communities_hylids_and_fd10k_cont_bis_V3.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
plot(fddis_ras, main="dispersion filtered") #dispersion
plot(fdric_ras, main="richness filtered") #FRic

#### ALL PLOTS - null models ####
#par(mfrow=c(2,2))
plot(TD3, main="Taxonomic diversity") #TD
plot(pd_ras, main="Phylogenetic diversity") #PD
plot(fdric_ras, main="Functional richness") #FRic
plot(fddis_ras, main="Functional dispersion") #dispersion

#par(mfrow=c(1,2)) #to compare
plot(TD3, main="Non-filtered")
plot(TD3_filtre, main="Filtered")
plot(pd_ras, main="Phylogenetic diversity")
plot(pd_ras_filter, main="Phylogenetic diversity filtered")

############################ MAPS - Functional diversity - Random traits ###########################################################################################################
#### Need to charge marco1bis_filter_sans_na2 to calculate the FD

# Charger le package utils
library(utils)
# Liste des traits
traits <- names(trait_filter[,c(2:6)])
# Générer toutes les combinaisons de 3 traits parmi les 5
combinations <- combn(traits, 3)

##Comb1/filter: size_Celio, head_width, tibia_mm --> Comparative models
##Comb2: size_Celio, head_lenght, tibia_mm  
##Comb3: size_Celio, head_lenght,leg_mm
##Comb4: size_Celio, head_width, head_lenght
##Comb5: size_Celio, head_width, leg_mm
##Comb6: size_Celio, tibia_mm, leg_mm
##Comb7: head_width, head_lenght, tibia_mm
##Comb8: head_width, head_lenght, leg_mm
##Comb9: head_lenght, tibia_mm, leg_mm
##Comb10: head_width, tibia_mm, leg_mm

###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
#trait_filter (c(2,4,5)) size_Celio, head_width, tibia_mm
trait_Comb<-trait_filter 
rownames(trait_Comb)<-trait_Comb[,1] #147 entries
trait_Comb2<-trait_Comb[,c(2,3,5)] #size_Celio, head_lenght, tibia_mm 
trait_Comb3<-trait_Comb[,c(2,3,6)] #size_Celio, head_lenght,leg_mm
trait_Comb4<-trait_Comb[,c(2,3,4)] #size_Celio, head_width, head_lenght
trait_Comb5<-trait_Comb[,c(2,4,6)] #size_Celio, head_width, leg_mm
trait_Comb6<-trait_Comb[,c(2,5,6)] #size_Celio, tibia_mm, leg_mm
trait_Comb7<-trait_Comb[,c(3,4,5)] #head_width, head_lenght, tibia_mm
trait_Comb8<-trait_Comb[,c(3,4,6)] #head_width, head_lenght, leg_mm
trait_Comb9<-trait_Comb[,c(3,5,6)] #head_lenght, tibia_mm, leg_mm
trait_Comb10<-trait_Comb[,c(4,5,6)] #head_width, tibia_mm, leg_mm


#### COMB 2 ####

for (i in 1:3){
  trait_Comb2[,i]<-log(trait_Comb2[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb2),colnames(marco1_bis[,2:148])) #checking

rownames(trait_Comb2)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb2<-FD::dbFD(trait_Comb2,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

any(rowSums(marco1_bis[, 2:ncol(marco1_bis)]) == 0) #Check
marco1bis_filter_sans_na2 <- marco1bis_filter_sans_na1[-which(rowSums(marco1bis_filter_sans_na1[, 2:ncol(marco1bis_filter_sans_na1)]) == 0), ]
rownames(trait2_Comb2)==colnames(marco1bis_filter_sans_na2[,2:length(marco1bis_filter_sans_na2)]) #if false: not the same
FDindices_comb2<-FD::dbFD(trait2_Comb2,marco1bis_filter_sans_na2[,2:length(marco1bis_filter_sans_na2)],calc.FDiv=F) 

#map2<-rast(Stack_full2[[1]])
fddis_rasc2<-r2_filter
fdric_rasc2<-r2_filter
values(fddis_rasc2)<-NA # Erase all values from the distribution map
values(fdric_rasc2)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c2<-marco1_bis
marco_c2$fdDis<-FDindices_comb2$FDis
marco_c2$fdRic<-FDindices_comb2$FRic

fddis_rasc2[marco_c2$grilla]<-marco_c2$fdDis
fdric_rasc2[marco_c2$grilla]<-marco_c2$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc2, "hylid_functional_dispersion10k_cont_V3_comb2.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc2, "hylid_functional_richness10k_cont_V3_comb2.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c2, file = "communities_hylids_and_fd10k_cont_bis_V3_comb2.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc2, main="FRic filter Comb2") #FRic
plot(fddis_rasc2, main="dispersion filter Comb2") #dispersion


#### COMB 3 ####

for (i in 1:3){
  trait_Comb3[,i]<-log(trait_Comb3[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb3),colnames(marco1_bis[,2:148])) #checking
rownames(trait_Comb3)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb3<-FD::dbFD(trait_Comb3,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc3<-r2_filter
fdric_rasc3<-r2_filter
values(fddis_rasc3)<-NA # Erase all values from the distribution map
values(fdric_rasc3)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c3<-marco1_bis
marco_c3$fdDis<-FDindices_comb3$FDis
marco_c3$fdRic<-FDindices_comb3$FRic

fddis_rasc3[marco_c3$grilla]<-marco_c3$fdDis
fdric_rasc3[marco_c3$grilla]<-marco_c3$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc3, "hylid_functional_dispersion10k_cont_V3_comb3.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc3, "hylid_functional_richness10k_cont_V3_comb3.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c3, file = "communities_hylids_and_fd10k_cont_bis_V3_comb3.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc3, main="FRic filter Comb3") #FRic
plot(fddis_rasc3, main="dispersion filter Comb3") #dispersion

#### COMB 4 ####

for (i in 1:3){
  trait_Comb4[,i]<-log(trait_Comb4[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb4),colnames(marco1_bis[,2:148])) #checking

rownames(trait_Comb4)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb4<-FD::dbFD(trait_Comb4,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

#map2<-rast(Stack_full2[[1]])
fddis_rasc4<-r2_filter
fdric_rasc4<-r2_filter
values(fddis_rasc4)<-NA # Erase all values from the distribution map
values(fdric_rasc4)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c4<-marco1_bis
marco_c4$fdDis<-FDindices_comb4$FDis
marco_c4$fdRic<-FDindices_comb4$FRic

fddis_rasc4[marco_c4$grilla]<-marco_c4$fdDis
fdric_rasc4[marco_c4$grilla]<-marco_c4$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc4, "hylid_functional_dispersion10k_cont_V3_comb4.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc4, "hylid_functional_richness10k_cont_V3_comb4.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c4, file = "communities_hylids_and_fd10k_cont_bis_V3_comb4.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc4, main="FRic filter Comb4") #FRic
plot(fddis_rasc4, main="dispersion filter Comb4") #dispersion

#### COMB 5 ####

for (i in 1:3){
  trait_Comb5[,i]<-log(trait_Comb5[,i]) ## This could be updated to another mathematic transformation if needed
}

# rows in alphabetical order
#ordered_columns_comb5 <- order(rownames(trait1_Comb5))
# reorganize dataframe
#trait2_Comb5 <- trait1_Comb5[ordered_columns_comb5,]
##109 entries by alphabetical order

setdiff(rownames(trait_Comb5),colnames(marco1_bis[,2:148])) #checking
rownames(trait_Comb5)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb5<-FD::dbFD(trait_Comb5,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc5<-r2_filter
fdric_rasc5<-r2_filter
values(fddis_rasc5)<-NA # Erase all values from the distribution map
values(fdric_rasc5)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c5<-marco1_bis
marco_c5$fdDis<-FDindices_comb5$FDis
marco_c5$fdRic<-FDindices_comb5$FRic

fddis_rasc5[marco_c5$grilla]<-marco_c5$fdDis
fdric_rasc5[marco_c5$grilla]<-marco_c5$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc5, "hylid_functional_dispersion10k_cont_V3_comb5.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc5, "hylid_functional_richness10k_cont__V3_comb5.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c5, file = "communities_hylids_and_fd10k_cont_bis_V3_comb5.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc5, main="FRic filter Comb5") #FRic
plot(fddis_rasc5, main="dispersion filter Comb5") #dispersion

#### COMB 6 ####

for (i in 1:3){
  trait_Comb6[,i]<-log(trait_Comb6[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb6),colnames(marco1_bis[,2:148])) #checking
rownames(trait_Comb6)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb6<-FD::dbFD(trait_Comb6,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc6<-r2_filter
fdric_rasc6<-r2_filter
values(fddis_rasc6)<-NA # Erase all values from the distribution map
values(fdric_rasc6)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c6<-marco1_bis
marco_c6$fdDis<-FDindices_comb6$FDis
marco_c6$fdRic<-FDindices_comb6$FRic

fddis_rasc6[marco_c6$grilla]<-marco_c6$fdDis
fdric_rasc6[marco_c6$grilla]<-marco_c6$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc6, "hylid_functional_dispersion10k_cont_V3_comb6.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc6, "hylid_functional_richness10k_cont_V3_comb6.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c6, file = "communities_hylids_and_fd10k_cont_bis_V3_comb6.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc6, main="FRic filter Comb6") #FRic
plot(fddis_rasc6, main="dispersion filter Comb6") #dispersion

#### COMB 7 ####
for (i in 1:3){
  trait_Comb7[,i]<-log(trait_Comb7[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb7),colnames(marco1_bis[,2:148])) #checking
rownames(trait_Comb7)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb7<-FD::dbFD(trait_Comb7,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 
##donne des warnings !!

fddis_rasc7<-r2_filter
fdric_rasc7<-r2_filter
values(fddis_rasc7)<-NA # Erase all values from the distribution map
values(fdric_rasc7)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c7<-marco1_bis
marco_c7$fdDis<-FDindices_comb7$FDis
marco_c7$fdRic<-FDindices_comb7$FRic

fddis_rasc7[marco_c7$grilla]<-marco_c7$fdDis
fdric_rasc7[marco_c7$grilla]<-marco_c7$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc7, "hylid_functional_dispersion10k_cont_V3_comb7.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc7, "hylid_functional_richness10k_cont_V3_comb7.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c7, file = "communities_hylids_and_fd10k_cont_bis_V3_comb7.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc7, main="FRic filter Comb7") #FRic
plot(fddis_rasc7, main="dispersion filter Comb7") #dispersion

#### COMB 8 ####
for (i in 1:3){
  trait_Comb8[,i]<-log(trait_Comb8[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb8),colnames(marco1_bis[,2:148]))
rownames(trait_Comb8)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb8<-FD::dbFD(trait_Comb8,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc8<-r2_filter
fdric_rasc8<-r2_filter
values(fddis_rasc8)<-NA # Erase all values from the distribution map
values(fdric_rasc8)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c8<-marco1_bis
marco_c8$fdDis<-FDindices_comb8$FDis
marco_c8$fdRic<-FDindices_comb8$FRic

fddis_rasc8[marco_c8$grilla]<-marco_c8$fdDis
fdric_rasc8[marco_c8$grilla]<-marco_c8$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc8, "hylid_functional_dispersion10k_cont_V3_comb8.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc8, "hylid_functional_richness10k_cont_V3_comb8.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c8, file = "communities_hylids_and_fd10k_cont_bis_V3_comb8.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc8, main="FRic filter Comb8") #FRic
plot(fddis_rasc8, main="dispersion filter Comb8") #dispersion

#### COMB 9 ####
for (i in 1:3){
  trait_Comb9[,i]<-log(trait_Comb9[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb9),colnames(marco1_bis[,2:148]))
rownames(trait_Comb9)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb9<-FD::dbFD(trait_Comb9,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc9<-r2_filter
fdric_rasc9<-r2_filter
values(fddis_rasc9)<-NA # Erase all values from the distribution map
values(fdric_rasc9)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c9<-marco1_bis
marco_c9$fdDis<-FDindices_comb9$FDis
marco_c9$fdRic<-FDindices_comb9$FRic

fddis_rasc9[marco_c9$grilla]<-marco_c9$fdDis
fdric_rasc9[marco_c9$grilla]<-marco_c9$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc9, "hylid_functional_dispersion10k_cont_V3_comb9.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc9, "hylid_functional_richness10k_cont_V3_comb9.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c9, file = "communities_hylids_and_fd10k_cont_bis_V3_comb9.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc9, main="FRic filter Comb9") #FRic
plot(fddis_rasc9, main="dispersion filter Comb9") #dispersion

#### COMB 10 ####
for (i in 1:3){
  trait_Comb10[,i]<-log(trait_Comb10[,i]) ## This could be updated to another mathematic transformation if needed
}

setdiff(rownames(trait_Comb10),colnames(marco1_bis[,2:148]))
rownames(trait_Comb10)==colnames(marco1_bis[,2:length(marco1_bis)]) #if false: not the same
FDindices_comb10<-FD::dbFD(trait_Comb10,marco1_bis[,2:length(marco1_bis)],calc.FDiv=F) 

fddis_rasc10<-r2_filter
fdric_rasc10<-r2_filter
values(fddis_rasc10)<-NA # Erase all values from the distribution map
values(fdric_rasc10)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco_c10<-marco1_bis
marco_c10$fdDis<-FDindices_comb10$FDis
marco_c10$fdRic<-FDindices_comb10$FRic

fddis_rasc10[marco_c10$grilla]<-marco_c10$fdDis
fdric_rasc10[marco_c10$grilla]<-marco_c10$fdRic

setwd("~/Master/M2/Internship_M2/analyse/frogs/figures/maps_V3")

writeRaster(fddis_rasc10, "hylid_functional_dispersion10k_cont_V3_comb10.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_rasc10, "hylid_functional_richness10k_cont_V3_comb10.tif", filetype = "GTiff")  ###in tif format
write.table(marco_c10, file = "communities_hylids_and_fd10k_cont_bis_V3_comb10.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
#par(mfrow=c(1,2))
plot(fdric_rasc10, main="FRic filter Comb10") #FRic
plot(fddis_rasc10, main="dispersion filter Comb10") #dispersion

###FAIRE UNE FONCTION QUI FAIT TOUT CA !! 


