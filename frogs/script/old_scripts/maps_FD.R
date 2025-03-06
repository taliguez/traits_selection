################################################################################################################################################## My script
####################################### SCRIPT - Internship M2 ##########################################################################################
#########################################  Tali - Guez ############################################################################################
##################################################################################################################################################My script


########################################### LIBRAIRIES ###########################################################################################################

##Calculate FD
#install.packages('FD')
library('FD')

### Shapefile (replace rgdal)
#install.packages("terra")
library(terra)

########################################### Load map - atlantic forest ###########################################################################################################

setwd("~/Master/M2/Internship_M2/analyse")

###load mask to use (AF) ##change path
#selected_mask <- vect(""/Users/andreapaz/Dropbox/Guano's_files/AF_shapefile_from_Guano", "atlantic_forest")
selected_mask <- vect("C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/AF_shapefile_from_Guano/AF_shapefile_from_Guano", "atlantic_forest")
selected_mask <- aggregate(selected_mask, dissolve = TRUE) #aggregate and dissolve polygons of a spatial object (selected_mask) into a single polygon
print(selected_mask)
plot(selected_mask)


########################################### Function to rasterize shapefile into ascii format ###########################################################################################################
####First stack rasters with points ####
### function to rasterize points
#rasterize_species <- function(x, mask=selected_mask, r=raster_base) {
  ##r<-raster(ncol=300,nrow=400,resolution=resolution,ext=extent(selected_mask))
  ##res(r)<-resolution 
  ##r<-extend(r,selected_mask)
  ##r<-crop(r,extent(selected_mask))
  ##r<-mask(r,selected_mask)
  #values(r)<-0 #raster initialization
  #map<-vect(file.path(maps_folder, paste0(x, ".shp"))) #lire shapefile terra package: shapefile is in the maps_folder
  #r <- rasterizeGeom(x=map, y= r, field = 1, update = TRUE, background = 0) # Rasterization of geometric properties of vector data. 
  #r <- mask(r, selected_mask) #mask values in a SpatRaster
  #valor <- unique(values(r)) #get the unique cell values of a SpatRaster of the attributes of a SpatVector 
  
  #if (length(valor) == 1 && is.na(valor)) {
    #return(NULL) #verification
  #} else if (length(valor) == 2 && valor[2] == 0) {
    #return(NULL)
  #} else {
    #writeRaster(x=r, paste0(x, ".asc"), filetype = "AAIGrid") 
    #return(rast(paste0(x, ".asc")))
  #}
#}
#### INITIALIZATION 
#resolution<-0.08333333 ##10km
#maps_folder<-getwd()
#distribution_files<-list.files(path=getwd(), pattern= "*.shp$") #list every file in the actual maps_folder with *.shp$ motifs
#species_names_points<-sub(".shp","",distribution_files) ##really file names # extract species name from the files
#raster_base<-Stack_models[[1]] ### this needs to be created with a different path section must be moved after having the stack
#raster_base<-raster(models[139]) ### this needs to be created with a different path


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
setwd("~/Master/M2/Internship_M2/analyse/points_14model")
maps_folder <- "C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/points_14model"
#maps_folder<-getwd()

distribution_files<-list.files(path=getwd(), pattern= "*.shp$") #list every file in the actual maps_folder with *.shp$ motifs
species_names_points<-sub(".shp","",distribution_files) ##really file names # extract species name from the files

output_dir <- "C:/Users/talig/OneDrive/Documents/Master/M2/Internship_M2/analyse/points_14model/points_14model_test" #output directory
#rasterize_species2(x = "Points_Scinax_uruguayus", mask = selected_mask, output_dir = output_dir)

##Function to do all 14 files
rasters<-lapply(species_names_points, function(species) {
  rasterize_species2(x = species, mask = selected_mask, output_dir = output_dir)
})


#Rasterize distribution files and keep only those of interest
#layers<-lapply(species_names_points,rasterize_species2)
##get species names from with grep removing the Points part. 
#species_names_points<-sub("Points_","",species_names_points)
#assign species names to layers with species distributions
#layers[sapply(layers,    is.null)]<-NULL
#Stack_points<-rast(layers)

############################ MAPS - TAXONOMIC DIVERSITY ###########################################################################################################

###load SDMs and make correct size
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/") ##to folder with SDMs
#setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap\\OR_AUC_T10/")
setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_10k") ##to folder with SDMs
#setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_T10")

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

#setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_10k/output_rasters")
setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_10k")

# Load models
#models <- list.files(pattern = "AF_10k_")
models<-list.files(pattern = "AF_10k_.*\\.asc$") #from the working directory
#models2<-as.list(models) #transform models into a list

#To Stack
#use c(x1, x2, x3) when providing separate arguments
#use rast(list(x1,x2,x3)) when providing the arguments as a list.

# Load each raster file individually and combine them into a single SpatRaster
raster_list<-lapply(models, rast) #transform models into list

# Combine the rasters into a single SpatRaster
Stack_models<-rast(raster_list) #create a SpatRaster with 3 layers
plot(Stack_models)
#names(Stack_models)
#Stack_models2 <- c(raster_list) #create list of 3 separate SpatRasters #not a SpatRaster
#plot(Stack_models2[[1]])
#Stack_models3 <- c(models2)
#Stack_models4 <- rast(models2) #ne marche pas

##Verification for each layers
#For Stack_models: one layer at a time #it will be 104 maps
for (i in 1:nlyr(Stack_models)) {
  print(paste("Layer", i, ":", names(Stack_models)[i]))
  print(global(Stack_models[[i]], "sum", na.rm = TRUE))  # Somme des valeurs
  plot(Stack_models[[i]])  # Visualisez la couche
}

#See if there are values in the raster (different than NA)
for (i in 1:nlyr(Stack_models)) {
  print(values(raster_list[[i]])[!is.na(values(raster_list[[i]]))]) #there are some values with 1 or 0
}

#values(raster_list[[1]])[!is.na(values(raster_list[[1]]))] #there are some values with 1 or 0

#To clean the files name
species_names_models<-sub("_OR_AUC_T10", "", models) #take _OR_AUC_T10 off
species_names_models<-sub("AF_10k_AF_", "", species_names_models) #take AF_10k_AF_ off
#species_names_models<-sub("AF_10k_AF_10k_AF_", "", species_names_models) #take AF_10k_AF_10k_AF_ off
species_names_models<-sub(".asc", "", species_names_models, fixed = TRUE) #take off .asc --> only have the species name

setwd("~/Master/M2/Internship_M2/analyse/Points_only_10k_and_14")
#setwd("~/Master/M2/Internship_M2/analyse/Points_only_10k")

# Load the points
points<-list.files(pattern = ".asc")
Stack_points<-rast(points) #into 1 spatRaster
plot(Stack_points)
species_names_points<-sub(".asc", "", points, fixed = TRUE) #take off .asc 
species_names_points<-sub("Points_", "", species_names_points) #take off Points_ 
species_names_full<-c(species_names_points, species_names_models) #have every species =37 species 
#158 names of species in total

#### To get richness and community composition ####
#Stack_full<-stack(Stack_points,Stack_models) #stack the rasters together (combine 2 stacks of rasters into a single one)
#TD<-calc(Stack_full,sum,na.rm=T) #calculate species richness for each cell in the stacked raster
#TD<-crop(TD,selected_mask) #crop the raster to the extent of selected_mask
#TD<-mask(TD,selected_mask) #masks the raster using selected_mask, setting values outside the mask to NA
#plot(TD)
#writeRaster(TD,"richness_hylids_complete_0_083.asc") #write the raster to an ASCII file 

Stack_full <- c(Stack_points, Stack_models) 
plot(Stack_full)
#Stack_full2 <- c(Stack_points, Stack_models2) 
#plot(Stack_full)
print(global(Stack_full, "sum", na.rm = TRUE))

TD<-app(Stack_full, sum, na.rm = TRUE) #Sum species of each layers
#app= apply a function to the values of each cell of a SpatRaster
plot(TD)

# Values in TD before cropping
print(values(TD))
# Check if there are any non-NaN values in TD before cropping
non_nan_values <- values(TD)[!is.na(values(TD))]
print(global(TD, "sum", na.rm = TRUE))

TD2<-crop(TD, selected_mask)
plot(TD2)

#Values in TD after cropping
print(values(TD2))
# Check the extent of the raster and the mask
print(ext(TD2)) #Extent of TD after cropping
print(ext(selected_mask)) #Extent of selected_mask
print(global(TD2, "sum", na.rm = TRUE))


TD3<-mask(TD2, selected_mask)

#Values in TD after masking
print(values(TD3))
print(global(TD3, "sum", na.rm = TRUE))

plot(TD3)
setwd("~/Master/M2/Internship_M2/analyse")
writeRaster(TD3, "richness_hylids_complete_0_083.tif", filetype = "GTiff")  ###in tif format

TD_test <- app(Stack_models, sum, na.rm = TRUE)
plot(TD_test)
print(global(TD_test, "sum", na.rm = TRUE))

TD_test2 <- app(Stack_points, sum, na.rm = TRUE)
plot(TD_test2)
print(global(TD_test2, "sum", na.rm = TRUE))


############################ MAPS - PHYLOGENETIC DIVERSITY ###########################################################################################################

library(ape)
library(geiger)
library(picante)

####Load Phylogeny, generate list of species
##To read a phylogeny in newick format use read.tree instead of read.nexus
#In working directory always a file with the same name: phylogeny.nex if not then change name here
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/PD_data/")
#setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap")
setwd("~/Master/M2/Internship_M2/analyse/PD_data")

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
### Trim the distribution data to match phylogeny
#species_names_full<-species_names_full[species_names_full %in% user_phylogeny$tip,label]
#distribution_files<-sub("*$","_OR_AUC_T10.asc",species_names)


#########################
####START ANALYSES#######
#########################

##To get richness and community composition
Stack_full3<-c(Stack_points,Stack_models) 
#plot(Stack_full3)
#print(global(Stack_full3, "sum", na.rm = TRUE))

##polygon with mask is called selected_maks from above


# Calcul de la diversité phylogénétique
tabla2 <- as.data.frame(species_names_full)
colnames(tabla2) <- "Grilla" #nom colonne = Grilla

#Create empty raster of the desired area and resolution to assign pixel numbers
r2 <- Stack_full3[[1]]
res(r2)<-res(Stack_full3) #resolution
#r[is.na(r[])]<-0
#r<-crop(r,selected_mask)
#r<-mask(r,selected_mask)
#r<-aggregate(r,fact=10,fun=max)

grilla2 <- r2
names(grilla2) <- "grilla"
grilla2[1:ncell(grilla2)] <- 1:ncell(grilla2)
lista_grilla2 <- list(grilla2)
names(lista_grilla2) <- "Grilla"

#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_full3) <- as.vector(tabla2$Grilla)
#Stack <- rast(lista_grilla$Grilla, Stack_full) ###Full stack of all maps
Stack2 <- c(lista_grilla2$Grilla, Stack_full3) ###Full stack of all maps


#Turn maps into dataframe for computation of PD
community_data <- as.data.frame(Stack2)
##remove rows that are NA 
community_data1<-community_data[-which(rowSums(!is.na(community_data[, 2:159])) == 0),]
species_names<-colnames(community_data1)[2:length(community_data1)]  #Store species names

#Store 
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
setwd("~/Master/M2/Internship_M2/analyse")
write.table(community_data1, file = "communities_hylids10k.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara

community_data1[is.na(community_data1)]<-0
#In the community data frame NA must be eliminated done before is this if you load community data and not maps?
#community_data=na.omit(community_data)
#head(community_data)

#III-Phylogenetic diversity computation 
#computes only Faith's PD others may be added

pd.result <-picante::pd(community_data1[,2:ncol(community_data1)],pruned.tree,include.root = F) 
##if it fails switch to
##PhyloMeasures::pd.query(tree, matrix,)

#Add the pixel PD value to data frame
community_data1$pd<-pd.result[[1]]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data1, file = "communities_hylids_and_pd10k.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area
values(r2)<-0
pd_ras<-r2
values(pd_ras)<-NA #Eliminate every value they will be replaced by PD values further down

#2- Assign PD values to raster
pd_ras[community_data1$grilla]<-community_data1$pd

#3- Save raster to file 
#writeRaster(pd_ras,"hylids_PD10k.asc",format="ascii")
writeRaster(pd_ras, "hylids_PD10k.tif", filetype = "GTiff")  ###in tif format

#4-Optional plotting map in R 
plot(pd_ras)

#par(mar = c(4, 4, 2, 2))
#par(mfrow=c(1,1))


############################ MAPS - Functional diversity - Null models ###########################################################################################################

### Author: Sebastien Vill�ger, adapted by Claire Fortunel (Please acknowledge as appropriate)  

# Notations corresponds with Vill�ger et al. (2008) Ecology, 89: 2290-2301 for FRic, FEve, FDiv; and Bellwood et al. (2006) Proc. R. Soc. B., 273: 101-107 for FSpe

# Function to calculate the four Functional diversity indices (So far we only use the FRic index since no abundance data is available from distribution maps)
library(geometry)
library(ape)
library(ade4)
library(FD)

##Load maps exactly as for PD if not yet loaded
#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/")
#setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_T10/output_rasters")
setwd("~/Master/M2/Internship_M2/analyse/OR_AUC_10k")
models<-list.files(pattern="AF_") 
#Stack_models<-stack(models)
#Stack_models<-rast(raster_list) #create a SpatRaster with 3 layers #Fait avant
#plot(Stack_models)
species_names_models<-sub("_OR_AUC_T10","",models) 
#species_names_models<-sub("AF_","",species_names_models)
species_names_models<-sub("AF_10k_AF_","",species_names_models)
species_names_models<-sub(".asc","",species_names_models,fixed=T) #To have all the species names

#setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/Points_only/") 
setwd("~/Master/M2/Internship_M2/analyse/Points_only_10k_and_14")
points<-list.files(pattern=".asc")
#Stack_points<-stack(points)
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
trait<-read.csv("~/Master/M2/Internship_M2/analyse/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
#trait<-read.csv("Species_traits_imputed_lilian.csv")
###Trim trait data to match distribution data #2 don t match the ocean ones
trait<-trait[trait$Species %in% species_names_full,] 
#trait<-na.omit(trait)
### Trim the distribution data to match trait data, This is not needed in this case
#species_names_full<-species_names_full[species_names_full %in% trait1$Species]

tabla<-as.data.frame(species_names_full) #table with names of all species (158)
colnames(tabla)<-"Grilla" #nom col= Grilla

#Load all distribution maps and make them same extent as will be used
r1<-Stack_full2[[1]]
#print(global(Stack_full2, "sum", na.rm = TRUE))
#r<-extend(r,AF,value=0)
#r[is.na(r[])]<-0
#r<-crop(r,AF)
#r<-mask(r,AF)
#r<-aggregate(r,fact=10,fun=max)
grilla=r1 
names(grilla)="grilla" #
grilla[1:ncell(grilla)]<-1:ncell(grilla)
lista_grilla<-list(grilla)
names(lista_grilla)<-"Grilla"
#layers<-lapply(distribution_files,load_species)
#names(layers)<-as.vector(tabla$Grilla)
#layers[sapply(layers,is.null)]<-NULL
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_full2)<-as.vector(tabla$Grilla)
Stack<-c(lista_grilla$Grilla,Stack_full2) ###Full stack of all maps

#Turn maps into dataframe for computation of FD
marco<-as.data.frame(Stack)
marco<-marco[-which(rowSums(!is.na(marco[,2:159]))==0),] ## does not work yet #remove all values in col 2 to 159 when only NA in it
marco[is.na(marco)] = 0 #replace every NA by 0
species_names<-colnames(marco)[2:length(marco)] #Store species names
marco1<-marco
marco1<-marco[-which(rowSums(marco[,2:159])==0),] #remove rows where the sum of it is equal to 0

#marco2<-marco1
#rownames(marco2)<-marco1$grilla
#marco2<-marco2[,-1]

# Obtenir l'ordre alphabétique des colonnes
#ordered_columns <- order(names(marco2))
# Réorganiser le dataframe
#marco3 <- marco2[, ordered_columns]

# Obtenir l'ordre alphabétique des noms de colonnes, en excluant la première colonne
ordered_columns <- order(names(marco1)[-1])
# Réorganiser le dataframe, en conservant la première colonne à sa place
marco1_bis <- marco1[, c(1, ordered_columns + 1)]

###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
#trait1<-read.csv("Species_traits_imputed_lilian.csv")
trait1<-read.csv("~/Master/M2/Internship_M2/analyse/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")
#FD_dis<-FD::gowdis(trait1[,c(3,4,7,8)])
#FDDis<-FD::fdisp(FD_dis,marco)
##Trim trait data to match distribution data #2 don t match the ocean ones
#trait1<-trait1[trait1$Species %in% names(marco)[2:159],] #don't need that maybe
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(4,7,8)] #null model (from the article size_Celio, head_width, tibia_mm)

for (i in 1:3){
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

trait2<-trait1
# Obtenir l'ordre alphabétique des colonnes
ordered_columns2 <- order(rownames(trait1))
# Réorganiser le dataframe
trait2 <- trait1[ordered_columns2,]

rownames(trait2)==colnames(marco1_bis[,2:159]) #if false: not the same
#FDindices<-FD::dbFD(trait1,marco[,2:159],calc.FDiv=F) 
#FDindices<-FD::dbFD(trait1,marco1[,2:159],calc.FDiv=F) 
#FDindices<-FD::dbFD(trait2,marco3,calc.FDiv=F) 
FDindices<-FD::dbFD(trait2,marco1_bis[,2:159],calc.FDiv=F) 


map<-rast(Stack_full2[[1]])
fddis_ras<-r1
fdric_ras<-r1
values(fddis_ras)<-NA # Erase all values from the distribution map
values(fdric_ras)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
#marco$fdDis<-FDindices$FDis
#marco$fdRic<-FDindices$FRic
marco4<-marco1_bis
marco4$fdDis<-FDindices$FDis
marco4$fdRic<-FDindices$FRic

fddis_ras[marco4$grilla]<-marco4$fdDis
fdric_ras[marco4$grilla]<-marco4$fdRic
#writeRaster(fddis_ras,"hylid_functional_dispersion10k_cont.asc", format="ascii")
#writeRaster(fdric_ras,"hylid_functional_richness10k_cont.asc", format="ascii")
#write.table(marco, file = "communities_hylids_and_fd10k_cont.txt", append = FALSE,row.names=F,quote=F,sep="\t")

setwd("~/Master/M2/Internship_M2/analyse")

writeRaster(fddis_ras, "hylid_functional_dispersion10k_cont.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_ras, "hylid_functional_richness10k_cont.tif", filetype = "GTiff")  ###in tif format
write.table(marco1_bis, file = "communities_hylids_and_fd10k_cont.txt", append = FALSE,row.names=F,quote=F,sep="\t")
write.table(marco4, file = "communities_hylids_and_fd10k_cont_bis.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
plot(fddis_ras) #dispersion
plot(fdric_ras) #FRic


#### ALL PLOTS - null models ####
#par(mfrow=c(2,2))
plot(TD3, main="Taxonomic diversity") #TD
plot(pd_ras, main="Phylogenetic diversity") #PD
plot(fdric_ras, main="Functional richness") #FRic
plot(fddis_ras, main="Functional dispersion") #dispersion

####### NOT USED #####
###try with categorical
trait2<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
trait2<-read.csv("Species_traits_imputed_lilian.csv")
rownames(trait2)<-trait2[,1]
trait2<-trait2[,c(3,4,7,8)] #rep mode, size, head w tibiaL
for (i in 2:4)
{
  trait2[,i]<-log(trait2[,i]) ## This could be updated to another mathematic transformation if needed
}
##check if cat is cat
FDindices_cat<-FD::dbFD(trait2,marco[,2:159],calc.FDiv=F,corr="cailliez")

map<-raster(Stack_full[[1]])
fddisC_ras<-r
fdricC_ras<-r
values(fddisC_ras)<-NA # Erase alll values from the distribution map
values(fdricC_ras)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco$fdDisC<-FDindices_cat$FDis
marco$fdRicC<-FDindices_cat$FRic
fddisC_ras[marco$grilla]<-marco$fdDisC
fdricC_ras[marco$grilla]<-marco$fdRicC
writeRaster(fddisC_ras,"hylid_functional_dispersion10k_cat.asc", format="ascii")
writeRaster(fdricC_ras,"hylid_functional_richness10k_cat.asc", format="ascii")


FDind=function(trait,abund) {
  # T = number of traits
  T=dim(trait)[2]
  # c = number of communities
  C=dim(abund)[1]
  # check coherence of number of species in 'traits' and 'abundances'
  if (dim(abund)[2]!=dim(trait)[1]) stop(" Error : different number of species in 'trait' and 'abund' matrices ")
  # check format of traits values
  if (ncol(trait)<2) stop ("'Trait' must have at least 2 columns")
  if (is.numeric(trait)==F) stop ("Traits values must be numeric")
  # check absence of NA in 'traits'
  if (length(which(is.na(trait)==T))!=0) stop(" Error : NA in 'trait' matrix ")
  # replacement of NA in 'abund' by '0'
  abund[which(is.na(abund))]=0
  # definition of vector for results, with communities'names as given in 'abund'
  Nbsp=rep(NA,C) ; names(Nbsp)=row.names(abund)
  FRic=rep(NA,C) ; names(FRic)=row.names(abund)
  FEve=rep(NA,C) ; names(FEve)=row.names(abund)
  FDiv=rep(NA,C) ; names(FDiv)=row.names(abund)
  FSpe=rep(NA,C) ; names(FSpe)=row.names(abund)
  # scaling and centering of each trait according to all species values
  traitCS=scale(trait, center=TRUE, scale=TRUE)
  # functional specialization of each species (distance to point 0,0 in the standardized functional space)
  FSpeS=(apply(traitCS, 1, function(x) {x%*%x}))^0.5
  # loop to compute FRic, FEve, FDiv and FSpe on each community
  for (i in 1:C){
    # selection of species present in the community
    esppres=which(abund[i,]>0)
    #  number of species in the community
    S=length(esppres) ; Nbsp[i]=S
    # check if more species than traits
    # if (S<=T) stop(paste("Number of species must be higher than number of traits in community:",row.names(abund)[i]))
    # filter on 'trait' and 'abund' to keep only values of species present in the community
    tr=traitCS[esppres,] ; ab=as.matrix(abund[i,esppres])
    # scaling of abundances
    abondrel=ab/sum(ab)
    # Functional Diversity Indices
    # FRic
    # Using convhulln function
    # volume
    FRic[i]=round(convhulln(tr,"FA")$vol,6)
    FD_dis[i]<-FD::gowdis(tr)
    FDDis<-fdisp(FD_dis,marco)
    # identity of vertices
    vert0=convhulln(tr,"Fx TO 'vert.txt'")
    vert1=scan("vert.txt",quiet=T)
    vert2=vert1+1
    vertices=vert2[-1]
    # FEve
    # computation of inter-species euclidian distances
    distT=dist(tr, method="euclidian")
    # computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
    linkmst=mst(distT) ; mstvect=as.dist(linkmst)
    # computation of the pairwise cumulative relative abundances and conversion into 'dist' class
    abond2=matrix(0,nrow=S,ncol=S)
    for (q in 1:S)
      for (r in 1:S)
        abond2[q,r]=abondrel[q]+abondrel[r]
    abond2vect=as.dist(abond2)  # end of q,r
    # computation of weighted evenness (EW) for the (S-1) branches to link S species
    EW=rep(0,S-1)
    flag=1
    for (m in 1:((S-1)*S/2)){if (mstvect[m]!=0) {EW[flag]=distT[m]/(abond2vect[m]) ; flag=flag+1}}  # end of m
    # computation of the partial weighted evenness (PEW) and comparison with 1/S-1, and computation of FEve
    minPEW=-rep(0,S-1) ; OdSmO=1/(S-1)
    for (l in 1:(S-1))
      minPEW[l]=min((EW[l]/sum(EW)), OdSmO)  # end of l
    FEve[i]=round(((sum(minPEW))- OdSmO)/(1-OdSmO),6)
    # FDiv
    # traits values of vertices of the convex hull
    trvertices=tr[vertices,]
    # coordinates of the center of gravity of the vertices (Gv) of the convex hull
    baryv=apply(trvertices,2,mean)
    # euclidian distances to Gv (dB) of each of S species (centro de gravedad)
    distbaryv=rep(0,S)
    for (j in 1:S)
      distbaryv[j]=(sum((tr[j,]-baryv)^2) )^0.5  # end of j
    # mean euclidian distance to the center of gravity of the S species (i.e. mean of dB values)  (mean centro gravedad)
    # bigger effect if you are very abundant (andrea)
    meandB=mean(distbaryv)
    # deviation of each species dB from mean dB
    devdB=distbaryv-meandB
    # relative abundances-weighted mean deviation
    abdev=abondrel*devdB
    # relative abundances-weighted mean of absolute deviations
    ababsdev=abondrel*abs(devdB)
  } # end of i
  # result storage
  res=data.frame(Nbsp=Nbsp, FRic=FRic, FDis) ; row.names(res)=row.names(abund)
  invisible(res)
}# end of function


############################ MAPS - Functional diversity - Random traits ###########################################################################################################
#### Need to charge marco1_bis to calculate the FD

# Charger le package utils
library(utils)
# Liste des traits
traits <- names(trait[,c(4,6,7,8,9)])
# Générer toutes les combinaisons de 3 traits parmi les 5
combinations <- combn(traits, 3)

##Comb1: size_Celio, head_width, head_lenght--> NULL models
##Comb2: size_Celio, head_lenght, tibia_mm    
##Comb3: size_Celio, head_lenght,leg_mm
##Comb4: size_Celio, head_width, tibia_mm   
##Comb5: size_Celio, head_width, leg_mm
##Comb6: size_Celio, tibia_mm, leg_mm
##Comb7: head_width, head_lenght, tibia_mm
##Comb8: head_width, head_lenght, leg_mm
##Comb9: head_lenght, tibia_mm, leg_mm
##Comb10: head_width, tibia_mm, leg_mm

###compute Functional dispersion (including categorical data)
#Trait data: select file containing trait data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
#trait1<-read.csv("Species_traits_imputed_lilian.csv")
trait1<-read.csv("~/Master/M2/Internship_M2/analyse/Species_traits_imputed_lilian.csv",header=TRUE, sep = ";")

##Trim trait data to match distribution data #2 don t match the ocean ones
#trait1<-trait1[trait1$Species %in% names(marco)[2:159],] #don't need that maybe
rownames(trait1)<-trait1[,1] #row names= species

#### COMB 2 - size_Celio, head_lenght, tibia_mm   ####
trait1_Comb2<-trait1[,c(4,6,8)]
names(trait1_Comb2) #verify good columns

for (i in 1:3){
  trait1_Comb2[,i]<-log(trait1_Comb2[,i]) ## This could be updated to another mathematic transformation if needed
}

ordered_columns2_Comb2 <- order(rownames(trait1_Comb2)) #Column in alphabetical order
trait2_Comb2 <- trait1_Comb2[ordered_columns2_Comb2,] #Reorganize df
 
rownames(trait2_Comb2)==colnames(marco1_bis[,2:159]) #Verification: if false: not the same


##Calculation indices
FDindices_Comb2<-FD::dbFD(trait2_Comb2,marco1_bis[,2:159],calc.FDiv=F) 


map<-rast(Stack_full2[[1]])
fddis_ras<-r1
fdric_ras<-r1
values(fddis_ras)<-NA # Erase all values from the distribution map
values(fdric_ras)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
#marco$fdDis<-FDindices$FDis
#marco$fdRic<-FDindices$FRic
marco4<-marco1_bis
marco4$fdDis<-FDindices$FDis
marco4$fdRic<-FDindices$FRic

fddis_ras[marco4$grilla]<-marco4$fdDis
fdric_ras[marco4$grilla]<-marco4$fdRic
#writeRaster(fddis_ras,"hylid_functional_dispersion10k_cont.asc", format="ascii")
#writeRaster(fdric_ras,"hylid_functional_richness10k_cont.asc", format="ascii")
#write.table(marco, file = "communities_hylids_and_fd10k_cont.txt", append = FALSE,row.names=F,quote=F,sep="\t")

setwd("~/Master/M2/Internship_M2/analyse")

writeRaster(fddis_ras, "hylid_functional_dispersion10k_cont.tif", filetype = "GTiff")  ###in tif format
writeRaster(fdric_ras, "hylid_functional_richness10k_cont.tif", filetype = "GTiff")  ###in tif format
write.table(marco1_bis, file = "communities_hylids_and_fd10k_cont.txt", append = FALSE,row.names=F,quote=F,sep="\t")
write.table(marco4, file = "communities_hylids_and_fd10k_cont_bis.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#4-Optional plotting map in R 
plot(fddis_ras) #dispersion
plot(fdric_ras) #FRic


#### ALL PLOTS - null models ####
#par(mfrow=c(2,2))
plot(TD3, main="Taxonomic diversity") #TD
plot(pd_ras, main="Phylogenetic diversity") #PD
plot(fdric_ras, main="Functional richness") #FRic
plot(fddis_ras, main="Functional dispersion") #dispersion

