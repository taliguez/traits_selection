########################################## Script Tali ###########################################################
#merging datasets traits and community
#2025

rm(list = ls())

########################################## LIBRAIRIES ################################################################
library("tidyverse")

#install.packages("devtools")
#library("devtools")
#devtools::install_github("speckerf/treemendous")
library("treemendous")


########################################## INPUT DATASETS ################################################################
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/raw")
df<- read_csv("Maynard_trait_table_no_monocots.csv") #traits database
df_comm<- read_csv("REV_Community_matrix.csv") #community matrix 

df_comm<-read_csv("REV_Community_matrix_200.csv") #200

###Filter dataset for the fit --> to have only phylogenetic imputation measurement (phy)
df_filter<-df %>% 
  filter(fit=="phy")

#Transforming the traits dataset into columns and name of species with _ 
df_col_pred_traits<- df_filter %>%
  mutate(accepted_bin = str_replace_all(accepted_bin, " ", "_")) %>% #Replace spaces by underscores
  mutate(trait = str_replace_all(trait, " ", "_"))%>%  
  group_by(accepted_bin, fit, LAT, LON) %>%  #Group species lines
  summarise(pred_value = list(pred_value), trait = list(trait), .groups = "drop") %>%  #put values into list to avoid duplication of species
  unnest(c(trait, pred_value)) %>%  #Break the list into columns
  pivot_wider(names_from = trait, values_from = pred_value)  #transform into "widen data"
#Select by imputed values but could also be by observation

#create the input dataset for translate trees = input
df_col_pred_traits_treemendous<- df_col_pred_traits %>% #input or df
  distinct(accepted_bin) %>%
  separate(accepted_bin, into = c("Genus", "Species"), sep = "_", remove = FALSE) %>%  #create two columns
  filter(!(accepted_bin %in% c("Hovea_elliptica", "Pelea_elliptica"))) #remove those ambiguous species

#target
df_comm_treemendous<-df_comm %>%#target ###df_comm_200
  distinct(accepted_bin) %>%
  separate(accepted_bin, into = c("Genus", "Species"), sep = "_", remove = FALSE) 


#### Treemendous package - translate_trees ####
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean")
df_treemendous<-translate_trees(df = df_col_pred_traits_treemendous, target = df_comm_treemendous)
#write_csv(df_treemendous, "treemendous_species_match.csv")
df_treemendous<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/clean/species_names_treemendous//treemendous_species_match.csv")

#species_match<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/clean/treemendous_species_match.csv)
#Hovea elliptica and Pelea elliptica only in trait dataset so take them off to not have the treemendous_ambiguous_genera.csv dataset
#get list of correct species names: treemendous_ambiguous_genera.csv

df_treemendous2 <- df_treemendous %>% 
  select(Orig.Genus, Orig.Species, Matched.Genus,Matched.Species,accepted_bin, direct_match)%>%
  mutate(accepted_bin_old = if_else( 
    is.na(Orig.Genus) | is.na(Orig.Species), 
    NA_character_, 
    str_c(Orig.Genus, Orig.Species, sep = "_"))) %>% 
  mutate(accepted_bin_new = if_else(
    is.na(Matched.Genus) | is.na(Matched.Species), 
    NA_character_, 
    str_c(Matched.Genus, Matched.Species, sep = "_")
  ))

#write_csv(x=df_treemendous2, file="from_treemendous.csv", append=FALSE) #_200

#Add measurements
df_col_pred_traits_names<- df_col_pred_traits %>%
  left_join(df_treemendous2 %>% select(accepted_bin_old, accepted_bin_new, direct_match), 
            by = c("accepted_bin" = "accepted_bin_old")) %>% 
  relocate(accepted_bin_new, .after=accepted_bin) %>% 
  relocate(direct_match, .after = accepted_bin_new) 
#write_csv(x=df_col_pred_traits_names, file="traits_with_new_names.csv", append=FALSE)

setdiff(df_comm_treemendous$accepted_bin,df_col_pred_traits_names$accepted_bin_new) #rows that appear in x but not y

#Species that do not find a match from community to traits database
x <- df_comm_treemendous %>% 
  anti_join(df_col_pred_traits_names, by = c("accepted_bin" = "accepted_bin_new")) %>%
  mutate(present_in_treemendous = accepted_bin %in% df_treemendous2$accepted_bin_old)

setdiff(x$accepted_bin,df_comm$accepted_bin) #rows that appear in x but not y
#write_csv(x, "presence_in_comm_matrix_not_in_traits.csv")

#Species with duplicates
df_col_pred_traits_names_filtered <- df_col_pred_traits_names %>%
  filter(accepted_bin_new %in% df_comm$accepted_bin) #39,412 entries VS 38,064 entries in df_comm_treemendous
#Extract duplicates names
duplicates<-df_col_pred_traits_names_filtered %>%
  group_by(accepted_bin_new) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  distinct(accepted_bin_new) %>%
  count(accepted_bin_new)  #Count number of occurences of same sp

####### Add measurements to the matching species (df_new)  ###########
##I take the traits for the exact match and don't do the average, 
#if we don't have an exact match and it's not a duplicate I take those traits 
#and we'll see next week for the one that are not find by the translation from the community matrix ?

df_new <- df_comm_treemendous %>%
  left_join(df_col_pred_traits_names, by = c("accepted_bin" = "accepted_bin_new")) %>%
  mutate(direct_match = ifelse(is.na(direct_match), FALSE, direct_match)) %>%
  group_by(accepted_bin) %>%
  arrange(accepted_bin, desc(direct_match)) %>%
  slice(1) %>% #keep only the first line of the group to avoid duplicates
  ungroup() %>%
  select(accepted_bin, direct_match, everything())

#Identifier les lignes où accepted_bin.y est NA #same as dataset as x 
na_rows <- df_new %>%
  filter(is.na(accepted_bin.y))

#Verification
setdiff(x$accepted_bin,na_rows$accepted_bin) #rows that appear in x but not y

#### Dealing with non-matching species #####
#### Not matching species: try another method with matching(), resolve_synonyms(), and translate_trees() 
#either choose the backbone (WFO, BGCI...) or do them all at the same time
input <- df_col_pred_traits_treemendous %>% 
  matching() %>% #backbone = 'BGCI'
  resolve_synonyms()  #'BGCI'

flags_input<-input %>% highlight_flags('WFO')

input2 <- input %>%
  select(Accepted.Genus, Accepted.Species) %>%
  rename(Genus = Accepted.Genus, Species = Accepted.Species) %>% 
  drop_na() %>% #enlève 103 lignes! #52,092 entries
  distinct(Genus, Species) #51,627 entries
which(is.na(input2$Genus))
which(is.na(input2$Species))

#setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/species_names_treemendous")
#write_csv(x=input, file="species_names_traits.csv", append=FALSE)
#write_csv(x=flags_input, file="species_names_traits_complicated.csv", append=FALSE)

target<- x %>% 
  matching() %>% #backbone = 'BGCI'
  resolve_synonyms() #'BGCI'

flags_target<-target %>% highlight_flags('WFO')

target2 <- target %>% #from community
  select(Accepted.Genus, Accepted.Species) %>%
  rename(Genus = Accepted.Genus, Species = Accepted.Species)
which(is.na(target2$Genus))
which(is.na(target2$Species))

df_treemendous_tropical<-translate_trees(df = input2, target = target2)
df_treemendous_tropical2<-df_treemendous_tropical %>% 
  inner_join(target2, by = c("Orig.Genus" = "Genus", "Orig.Species" = "Species"))

df_treemendous_tropical_bis<-translate_trees(df =target2 , target =input2)
df_treemendous_tropical_bis2<-df_treemendous_tropical_bis %>% 
  inner_join(target2, by = c("Orig.Genus" = "Genus", "Orig.Species" = "Species")) %>% 
  mutate(accepted_bin_old = if_else( 
    is.na(Orig.Genus) | is.na(Orig.Species), 
    NA_character_, 
    str_c(Orig.Genus, Orig.Species, sep = "_"))) %>% 
  mutate(accepted_bin_new = if_else(
    is.na(Matched.Genus) | is.na(Matched.Species), 
    NA_character_, 
    str_c(Matched.Genus, Matched.Species, sep = "_")
  )) #have the new matching (old=comm, new=trait)

target_try<-target %>% 
  select(1,2,5,6,15) %>% 
  rename(accepted_bin_comm=accepted_bin) %>% 
  mutate(accepted_bin_trait = if_else( 
    is.na(Accepted.Genus) | is.na(Accepted.Species), 
    NA_character_, 
    str_c(Accepted.Genus, Accepted.Species, sep = "_")))

setdiff(target_try$accepted_bin_comm, df_new$accepted_bin)
setdiff(target_try$accepted_bin_comm, df_col_pred_traits$accepted_bin)
setdiff(target_try$accepted_bin_trait, df_col_pred_traits$accepted_bin)

na_rows_test <- na_rows %>% 
  left_join(df_treemendous_tropical_bis2 %>% select(accepted_bin_old, accepted_bin_new),
            by = c("accepted_bin" = "accepted_bin_old")) %>% #create a temporary column to stock the corresponding value from accepted_bin_new
  relocate(accepted_bin_new, .after=accepted_bin.y) %>% 
  mutate(accepted_bin.y = if_else(is.na(accepted_bin_new), accepted_bin.y, accepted_bin_new)) %>% #Replace accepted_bin.y by the values from temporary column
  select(-accepted_bin_new)   #Supresss temporary column

setdiff(na_rows_test$accepted_bin.y, df_col_pred_traits$accepted_bin)
setdiff(na_rows_test$accepted_bin, df_col_pred_traits$accepted_bin)


#Verify if values from accepted_bin.y exist in df_col_pred_traits_treemendous$accepted_bin
#and replace by NA if it is not the case
na_rows_test2 <- na_rows_test %>%
  mutate(accepted_bin.y = if_else(accepted_bin.y %in% df_col_pred_traits_treemendous$accepted_bin,
                                  accepted_bin.y, NA_character_)) %>% 
  select(-(6:26)) %>% 
  select(-direct_match) #suppress unecessary columns

na_rows_test3 <- na_rows_test2 %>%
  left_join(target_try, by = c("accepted_bin" = "accepted_bin_comm"))%>%
  left_join(df_col_pred_traits, by = c("accepted_bin_trait" = "accepted_bin"))%>% 
  select(-(5:8), -(10:30))%>%
  mutate(accepted_bin_trait = if_else(accepted_bin_trait %in% df_col_pred_traits_treemendous$accepted_bin,
                                  accepted_bin_trait, NA_character_)) %>% 
  mutate(accepted_bin_traits = coalesce(accepted_bin.y, accepted_bin_trait)) %>% #combiner les deux colonnes en une
  select(-(4:5))

# Étape 2 : Ajouter les colonnes correspondantes de df_col_pred_traits_names aux lignes où les espèces correspondent
#na_rows_test4<- na_rows_test3 %>%
 # left_join(df_col_pred_traits, by = c("accepted_bin_traits" = "accepted_bin"))
#69 species

#write_csv(x=df_col_pred_traits_names, file="traits_with_new_names.csv", append=FALSE)

####finish manually the rest ####
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/raw")
matched_trait_names<-read_csv("matched_trait_names.csv") #3354 sp, manually done datasets

matched_trait_names_x<- matched_trait_names %>% #missing species from traits database
  mutate(raw_name=str_replace_all(raw_name, " ", "_"),trait_name = str_replace_all(trait_name, " ", "_")) %>% 
  inner_join(x, by = c("raw_name" = "accepted_bin"))

manually<-na_rows_test3 %>% #4
  left_join(matched_trait_names_x %>% select(raw_name, trait_name),
            by = c("accepted_bin" = "raw_name")) %>% 
  distinct() #Pinus_henryi have 3 names which are false so #Pinus_tabuliformis & Pinus_henryi same measurement 

manually2 <- manually %>%
  mutate(accepted_bin_traits = coalesce(accepted_bin_traits, trait_name))%>% #add corresponding values of trait_name in the NA in accepted_bin_traits
  select(-trait_name) %>% #Supress temporary column trait_name
  distinct(accepted_bin, .keep_all = TRUE) %>% #Deal Pinus_henryi (Supress duplicates for Pinus_henryi to keep only one row)
  mutate(accepted_bin_traits = if_else(accepted_bin == "Pinus_henryi", "Pinus_tabuliformis", accepted_bin_traits)) #remplacer par Pinus_tabuliformis

#remplir les données des traits
manually3<-manually2 %>%
  left_join(df_col_pred_traits, by = c("accepted_bin_traits" = "accepted_bin"))
#final dataset with traits


####na_rows - The final ####
na_rows_manually<-na_rows %>% #4
  select(-(2:4), -(6:26)) %>% 
  rename(accepted_bin_traits=accepted_bin.y) %>% 
  left_join(matched_trait_names_x %>% select(raw_name, trait_name),
            by = c("accepted_bin" = "raw_name")) %>% 
  mutate(accepted_bin_traits = coalesce(accepted_bin_traits, trait_name))%>% #add corresponding values of trait_name in the NA in accepted_bin_traits
  select(-trait_name) %>% #Supress temporary column trait_name
  distinct() %>%   #Pinus_henryi have 3 names which are false so #Pinus_tabuliformis & Pinus_henryi same measurement 
  distinct(accepted_bin, .keep_all = TRUE) %>% #Deal Pinus_henryi (Supress duplicates for Pinus_henryi to keep only one row)
  mutate(accepted_bin_traits = if_else(accepted_bin == "Pinus_henryi", "Pinus_tabuliformis", accepted_bin_traits)) #remplacer par Pinus_tabuliformis

#remplir les données des traits
na_rows_manually2<-na_rows_manually %>%
  left_join(df_col_pred_traits, by = c("accepted_bin_traits" = "accepted_bin"))
#final dataset with traits
setdiff(na_rows_manually2$accepted_bin, df_new$accepted_bin )
setdiff(na_rows_manually2$accepted_bin_traits, df_col_pred_traits$accepted_bin )
setdiff(na_rows_manually2$accepted_bin, df_col_pred_traits$accepted_bin)
###ok

####merged dataset (na_rows_manually2) df_new and ####
# Étape 1 : Effectuer une jointure entre df_new et na_rows_manually2 pour obtenir les valeurs correspondantes
df_new_updated <- df_new %>%
  left_join(na_rows_manually2 %>% select(accepted_bin, accepted_bin_traits),
            by = c("accepted_bin" = "accepted_bin"))%>%
  mutate(accepted_bin.y = coalesce(accepted_bin.y, accepted_bin_traits)) %>%
  select(-accepted_bin_traits) %>% 
  rename(accepted_bin_traits=accepted_bin.y)
  
df_new_updated2<-df_new_updated %>% 
  select(-(6:length(df_new_updated))) %>% 
  left_join(df_col_pred_traits, by = c("accepted_bin_traits"="accepted_bin")) %>% 
  select(-Genus, -Species)

setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm")
#write_csv(x=df_new_updated2, file="all_traits_comm.csv", append=FALSE) #_200

###choisir les traits du dt_mat dataset: 
#But maybe scale them at each run so that its scaled for each analysis
#So that would be more like what people would do if only measuring those traits
#Col df_new
#accepted_bin #direct_match #Genus #Species #accepted_bin.y #fit #LAT #LON
#Wood_density #Root_depth #Leaf_N_per_mass #Leaf_P_per_mass #Stem_diameter 
#Bark_thickness #Seed_dry_mass #Leaf_K_per_mass #Stomatal_conductance 
#Leaf_thickness #Leaf_density #Leaf_Vcmax_per_dry_mass #Stem_conduit_diameter 
#Crown_diameter #Crown_height #Tree_height #Leaf_area #Specific_leaf_area

##### BUILDING DATASETS #####
rm(list = ls())
trait<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/raw/REV_Trait_matrix.csv")
#38,064 entries

all_traits<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/all_traits_comm.csv")
#38,064 rows and 24 columns
all_traits_200<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/all_traits_comm_200.csv")


##1. Selected traits, ecologically based ####
eco_based<-c("Bark_thickness", "Leaf_P_per_mass", "Root_depth","Seed_dry_mass", "Specific_leaf_area",
             "Stem_conduit_diameter", "Tree_height", "Wood_density")

#scale them
df_eco_based_scaled_200 <- all_traits_200 %>%
  select((1:6),all_of(eco_based)) %>% 
  select(-(2:6)) %>% 
  #mutate(across(where(is.numeric), scale)) #gives the same values
  gather(PC, value, -accepted_bin) %>% 
  group_by(PC) %>%
  mutate(value = (value - mean(value))/sd(value)) %>%
  ungroup %>%
  spread(PC, value)

setdiff(df_eco_based_scaled_200$accepted_bin, trait$accepted_bin) #same species

#write_csv(df_eco_based_scaled_200, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/eco_based_scaled_200.csv")

##2. Completely RANDOM x10####
traits_names <- c("Bark_thickness", "Crown_diameter", "Crown_height", "Leaf_area", 
            "Leaf_density", "Leaf_K_per_mass", "Leaf_N_per_mass", "Leaf_P_per_mass", 
            "Leaf_thickness", "Leaf_Vcmax_per_dry_mass", "Root_depth", "Seed_dry_mass", 
            "Specific_leaf_area", "Stem_conduit_diameter", "Stem_diameter", 
            "Stomatal_conductance", "Tree_height", "Wood_density")

set.seed(123)
samples <- map(1:10, ~ sample(traits_names, 8, replace = FALSE)) #Generate 10 randoms selection of 8 traits each

random_datasets <- map(samples, function(selected_traits) { #Create and scale 10 datasets in a list 
  all_traits_200 %>%
    select(accepted_bin, all_of(selected_traits)) %>%  # Garde l'ID (ex: accepted_bin)
    mutate(across(where(is.numeric), ~ (. - mean(.)) / sd(.)))  # Standardisation
}) #random_datasets <- map(samples, ~ select(all_traits_200, all_of(.x))) 

walk2(random_datasets, 1:10, ~ write_csv(.x, paste0("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/random_based_200_", .y, ".csv")))
#save datasets (purr package)

#test<-read_csv("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/REV_Traits_ANGIO_200_random_3.csv")
for (i in 1:10) {
  cat("Dataset", i, ":", paste(samples[[i]], collapse = ", "), "\n\n")
}

## 3. PCA ####
#Select 8 first component 
##### a. Scale and centered them ####
df_pca_scaled <- all_traits_200 %>%
  select(-(2:6)) %>% 
  #mutate(across(where(is.numeric), scale)) #gives the same values
  gather(PC, value, -accepted_bin) %>% 
  group_by(PC) %>%
  mutate(value = (value - mean(value))/sd(value)) %>%
  ungroup %>%
  spread(PC, value)

df_pca_scaled_rownames <- df_pca_scaled %>%
  column_to_rownames(var = "accepted_bin") #put accepted_bin in rownames

##### b. Specific libraries to add ####
library(stats)
library(vegan)
library(factoextra)
library("ggplot2")
library("ggfortify")
library("gridExtra")
library("corrplot")

##### c. PCA + get table ####
pca_traits<-prcomp(df_pca_scaled_rownames, center=FALSE, scale.=FALSE) #not centered and scaled because already done before
pca_traits
summary(pca_traits)
eig.val<-get_eigenvalue(pca_traits)
eig.val #dim 1 - dim 8 = 69.9% captured
#dim1-2:29.8%
#dim1-3:39.8%
#dim1-15:93.8%
#dim1-16:96.2%

#Visualization by a screeplot of which of the 8 PCAs contribute the most to the variation of the community variables (which PCA explains the best the data)
screeplot(pca_traits, bstick = TRUE, main="PCA on the traits variables") 
fviz_eig(pca_traits, col.var="blue", addlabels = TRUE)

pca_dataset<- as_tibble(pca_traits$x)%>%
  rownames_to_column("accepted_bin")

pca_dataset<- as.data.frame(pca_traits$x)
pca_dataset$accepted_bin <- rownames(pca_traits$x)

#8 first PC
pca_dataset<- as_tibble(pca_dataset)%>% #transform into tibble
  relocate(accepted_bin)%>%
  select(1:9) #select species names and 8 first principal component
#write_csv(pca_dataset, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/pca_based_200.csv")

#2 first PC
pca_dataset<- as_tibble(pca_dataset)%>% #transform into tibble
  relocate(accepted_bin)%>%
  select(1:3) #select species names and 8 first principal component
#write_csv(pca_dataset, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/pca_based_30_200.csv")

#3 first PC
pca_dataset<- as_tibble(pca_dataset)%>% #transform into tibble
  relocate(accepted_bin)%>%
  select(1:4) #select species names and 8 first principal component
#write_csv(pca_dataset, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/pca_based_40_200.csv")

#15 first PC
pca_dataset<- as_tibble(pca_dataset)%>% #transform into tibble
  relocate(accepted_bin)%>%
  select(1:16) #select species names and 8 first principal component
#write_csv(pca_dataset, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/pca_based_90_200.csv")

#16 first PC
pca_dataset<- as_tibble(pca_dataset)%>% #transform into tibble
  relocate(accepted_bin)%>%
  select(1:17) #select species names and 8 first principal component
#write_csv(pca_dataset, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/pca_based_96_200.csv")

##### d. Contrib of traits - Optional + corrplot####
var <- get_pca_var(pca_traits)
var
head(var$cos2, n=8)
corrplot(var$cos2, is.corr=FALSE)
head(var$contrib, n=8)
corrplot(var$contrib, is.corr=FALSE)    

fviz_cos2(pca_traits, choice = "var", axes = 1:2)
fviz_contrib(pca_traits, choice = "var", axes = 1:2)

fviz_pca_var(pca_traits,
             col.var = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "red"), #darkorange
             repel = TRUE
)
#fviz_pca_var(pca_traits, col.var = "black")
fviz_pca_var(pca_traits, repel = TRUE, select.var = list(cos2 = 8))

fviz_pca_var(pca_traits, repel = TRUE, select.var = list(contrib = 8))

fviz_pca_var(pca_traits,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "red"), #darkorange
             repel = TRUE
)

# Coordinates
head(var$coord, n=8)
# Cos2: quality on the factore map
head(var$cos2, n=8)
# Contributions to the principal components
head(var$contrib, n=8)

# Contributions of variables to PC1
a<-fviz_contrib(pca_traits, choice = "var", axes = 1)
# Contributions of variables to PC2
b<-fviz_contrib(pca_traits, choice = "var", axes = 2)
grid.arrange(a,b, ncol=2, top='Contribution of the variables to the first two PCs')

ind <- get_pca_ind(pca_traits)
ind
fviz_pca_ind(pca_traits,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE
)
# Total contribution on PC1 and PC2
fviz_contrib(pca_traits, choice = "ind", axes = 1:2)
autoplot(pca_traits, loadings=TRUE, loadings.colour='darkorchid4', loadings.label=TRUE, loadings.label.size=3)

#Representation of the 2 PCAs (PC1 and PC2) which explains the better the dataset
plot(pca_traits$x, main="PCA results on traits variables") #Plot the PCA results

biplot_traits <-envfit(pca_traits$x, df_pca_scaled_rownames, w=NULL, permutations = 1000)
plot(biplot_traits, add = TRUE, at = c (0,0)) 
plot(biplot_traits, p.max=0.05, at = c (0,0), col='red') 
biplot_traits #Get p-values

##4. random in CLUSTERS ####
trait_clusters <- tibble(
  trait = c("Leaf_density", "Wood_density", "Root_depth", "Specific_leaf_area", 
            "Leaf_thickness", "Leaf_N_per_mass", "Leaf_K_per_mass", "Leaf_P_per_mass", 
            "Stem_conduit_diameter", "Leaf_Vcmax_per_dry_mass", "Stomatal_conductance", 
            "Leaf_area", "Crown_height", "Crown_diameter", "Tree_height", 
            "Seed_dry_mass", "Bark_thickness", "Stem_diameter"),
  cluster = c(1, 1, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 8, 8)
)

#Random selection by cluster
set.seed(999)  # Pour la reproductibilité
selected_traits <- trait_clusters %>%
  group_by(cluster) %>%
  slice_sample(n = 1) %>%  #select randomly 1 trait by cluster 
  pull(trait)  

#Create dataset
cluster_random <- all_traits_200 %>%
  select(accepted_bin, all_of(selected_traits))%>%
  gather(PC, value, -accepted_bin) %>% 
  group_by(PC) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%  
  ungroup() %>%
  spread(PC, value)  
#write_csv(cluster_random, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/cluster_random_200.csv")

##5. Datasets by clusters : 3, 5 and 6 ####
trait = c("Leaf_density", "Wood_density", "Root_depth", "Specific_leaf_area", 
          "Leaf_thickness", "Leaf_N_per_mass", "Leaf_K_per_mass", "Leaf_P_per_mass", 
          "Stem_conduit_diameter", "Leaf_Vcmax_per_dry_mass", "Stomatal_conductance", 
          "Leaf_area", "Crown_height", "Crown_diameter", "Tree_height", 
          "Seed_dry_mass", "Bark_thickness", "Stem_diameter")

cluster3<-c("Specific_leaf_area","Leaf_thickness","Leaf_N_per_mass")
cluster5<-c("Stem_conduit_diameter","Leaf_Vcmax_per_dry_mass","Stomatal_conductance","Leaf_area")
cluster6<-c("Crown_height","Crown_diameter","Tree_height")

#Create dataset
clusters<- all_traits_200 %>%
  select(accepted_bin, all_of(cluster6))%>% #choisir cluster
  gather(PC, value, -accepted_bin) %>% 
  group_by(PC) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%  
  ungroup() %>%
  spread(PC, value)  
#write_csv(clusters, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/results/datasets_traits_comm/cluster_6_200.csv")


##### TABLE WITH ALL TRAITS #####
noms_datasets <- c("ecologically_based", "within_clusters", paste0("random_", 1:10))

ecologically_based_sorted <- sort(eco_based)
random_samples_sorted <- map(samples, sort)
selected_traits_sorted<-sort(selected_traits)

trait_combinations <- tibble(
  Dataset = noms_datasets,
  Trait_1 = c(ecologically_based_sorted[1], selected_traits_sorted[1], map_chr(random_samples_sorted, 1)),
  Trait_2 = c(ecologically_based_sorted[2], selected_traits_sorted[2], map_chr(random_samples_sorted, 2)),
  Trait_3 = c(ecologically_based_sorted[3], selected_traits_sorted[3], map_chr(random_samples_sorted, 3)),
  Trait_4 = c(ecologically_based_sorted[4], selected_traits_sorted[4], map_chr(random_samples_sorted, 4)),
  Trait_5 = c(ecologically_based_sorted[5], selected_traits_sorted[5], map_chr(random_samples_sorted, 5)),
  Trait_6 = c(ecologically_based_sorted[6], selected_traits_sorted[6], map_chr(random_samples_sorted, 6)),
  Trait_7 = c(ecologically_based_sorted[7], selected_traits_sorted[7], map_chr(random_samples_sorted, 7)),
  Trait_8 = c(ecologically_based_sorted[8], selected_traits_sorted[8], map_chr(random_samples_sorted, 8))
)

write_csv(trait_combinations, "~/Master/M2/Internship_M2/analyse/trees/datas/clean/trait_combinations.csv")
