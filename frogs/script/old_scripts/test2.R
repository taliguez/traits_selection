setwd("~/Master/M2/Internship_M2/analyse/figures/maps_V3")
x<-read.table("communities_hylids_and_fd10k_cont_V3.txt", header=TRUE)
unique(x)
summary(x)

x2<-read.table("communities_hylids_and_fd10k_cont_bis_V3_comb2.txt", header=TRUE)
hist(x2$fdRic)

x3<-read.table("communities_hylids_and_fd10k_cont_bis_V3_comb3.txt", header=TRUE)
hist(x3$fdRic)

x6<-read.table("communities_hylids_and_fd10k_cont_bis_V3_comb6.txt", header=TRUE)
summary(x6)
hist(x6$fdRic)

x7<-read.table("communities_hylids_and_fd10k_cont_bis_V3_comb7.txt", header=TRUE)
summary(x7)
hist(x7$fdRic)

x8<-read.table("communities_hylids_and_fd10k_cont_bis_V3_comb8.txt", header=TRUE)
summary(x8)
hist(x8$fdRic)



# Calculate the distance matrix
dist_matrix <- dist(trait_Comb7, method = "euclidean")
dist_matrix

# Check for zero distances
zero_distances <- which(dist_matrix == 0, arr.ind = TRUE)
print("Paires d'espèces avec des distances nulles :")
print(zero_distances)

library(vegan)  # For PCoA
# Visualize the distances using PCoA
pcoa_result <- pcoa(dist_matrix)
plot(pcoa_result, main = "PCoA of Trait Data")


# Calculate the distance matrix
dist_matrix <- dist(trait_Comb7, method = "euclidean")

# Perform PCoA
pcoa_result <- pcoa(dist_matrix)

# Print the PCoA results
print(pcoa_result)

# Plot the PCoA results using the vegan package
plot(pcoa_result, main = "PCoA of Trait Data")

# Alternatively, extract the coordinates and plot manually
scores <- scores(pcoa_result)
plot(scores[, 1], scores[, 2], main = "PCoA of Trait Data", xlab = "PCoA1", ylab = "PCoA2", pch = 19, col = "blue")



# Load necessary libraries
library(vegan)  # For PCoA

# Example data
trait_Comb7 <- trait_Comb[, c(3, 4, 5)] # Combinaison 7

# Calculate the distance matrix
dist_matrix <- dist(trait_Comb7, method = "euclidean")

# Perform PCoA
pcoa_result <- pcoa(dist_matrix)
str(pcoa_result)

# Print the PCoA results
print(pcoa_result)

# Extract the scores from the PCoA result
scores <- pcoa_result$vectors

# Plot the scores manually
plot(scores[, 1], scores[, 2], main = "PCoA of Trait Data", xlab = "PCoA1", ylab = "PCoA2", pch = 19, col = "blue")




#####
# Load necessary libraries
library(FD)

# Example data
trait_Comb6 <- trait_Comb[, c(2, 5, 6)] # Combinaison 6
trait_Comb7 <- trait_Comb[, c(3, 4, 5)] # Combinaison 7


for (i in 1:3){
  trait_Comb5[,i]<-log(trait_Comb5[,i]) ## This could be updated to another mathematic transformation if needed
}


for (i in 1:3){
  trait_Comb7[,i]<-log(trait_Comb7[,i]) ## This could be updated to another mathematic transformation if needed
}


# Check the traits in combinaison 6
print("Traits in combinaison 6:")
print(trait_Comb6)

# Check the traits in combinaison 7
print("Traits in combinaison 7:")
print(trait_Comb7)

# Check for NA values in combinaison 7
na_values_comb7 <- is.na(trait_Comb7)
print("NA values in combinaison 7:")
print(na_values_comb7)

# Check for unique values in combinaison 7
unique_values_comb7 <- apply(trait_Comb7, 2, unique)
print("Unique values in combinaison 7:")
print(unique_values_comb7)

# Check the community data
community <- x[, 2:length(x)]
print("Community data:")
print(community)


setdiff(rownames(trait_Comb7),colnames(x[,2:148])) #checking
rownames(trait_Comb7)==colnames(x[,2:length(x)]) #if false: not the same
# Calculate functional diversity indices for combinaison 7 with debug information
FDindices_comb7<-FD::dbFD(trait_Comb7,x[,2:length(x)],calc.FDiv=F) 

# Print the results for combinaison 7
print("FD indices for combinaison 7:")
print(FDindices_comb7)

# Check the warnings
warnings_comb7 <- warnings()
print("Warnings for combinaison 7:")
print(warnings_comb7)


# Calculate the distance matrix
dist_matrix <- dist(trait_Comb7, method = "euclidean")

# Check for zero distances in the distance matrix
zero_distances <- which(dist_matrix == 0, arr.ind = TRUE)
print("Zero distances in the distance matrix:")
print(zero_distances)

# Identify the species involved in zero distances
species_names <- rownames(trait_Comb7)
species_index1 <- zero_distances[1, 1]
species_index2 <- zero_distances[1, 2]
species_name1 <- species_names[species_index1]
species_name2 <- species_names[species_index2]
print(paste("Zero distance between:", species_name1, "and", species_name2))

# Check the trait data for the species involved
trait_data1 <- trait_Comb7[species_index1, ]
trait_data2 <- trait_Comb7[species_index2, ]
print(paste("Trait data for", species_name1, ":"))
print(trait_data1)
print(paste("Trait data for", species_name2, ":"))
print(trait_data2)




# Example data
trait_Comb7 <- trait_Comb[, c(3, 4, 5)] # Combinaison 7

for (i in 1:3){
  trait_Comb7[,i]<-log(trait_Comb7[,i]) ## This could be updated to another mathematic transformation if needed
}


# Calculate the distance matrix
dist_matrix <- dist(trait_Comb7, method = "euclidean")

# Check for zero distances in the distance matrix
zero_distances <- which(dist_matrix == 0, arr.ind = TRUE)
print("Zero distances in the distance matrix:")
print(zero_distances)

# Identify the species involved in zero distances
species_names <- rownames(trait_Comb7)
species_index1 <- zero_distances[1, 1]
species_index2 <- zero_distances[1, 2]
species_name1 <- species_names[species_index1]
species_name2 <- species_names[species_index2]
print(paste("Zero distance between:", species_name1, "and", species_name2))

# Check the trait data for the species involved
trait_data1 <- trait_Comb7[species_index1, ]
trait_data2 <- trait_Comb7[species_index2, ]
print(paste("Trait data for", species_name1, ":"))
print(trait_data1)
print(paste("Trait data for", species_name2, ":"))
print(trait_data2)

# Check for functionally singular species
community <- x[, 2:length(x)]
functionally_singular_species <- apply(community, 2, function(x) sum(x > 0) == 1)
print("Functionally singular species:")
print(names(which(functionally_singular_species)))

# Calculate functional diversity indices for combinaison 7 with debug information
FD_indices_comb7 <- dbFD(trait_Comb = trait_Comb7, marcol_bis = community, calc_FDiv = FALSE, debug = TRUE)

# Print the results for combinaison 7
print("FD indices for combinaison 7:")
print(FD_indices_comb7)

# Check the warnings
warnings_comb7 <- warnings()
print("Warnings for combinaison 7:")
print(warnings_comb7)




#############
setwd("~/Master/M2/Internship_M2/analyse/datas/cleaned/null_models/table_comparaison")
fdrich<-read.table("fdrich_pvalues", header=TRUE)
fddisp<-read.table("fdisp_pvalues", header=TRUE)
