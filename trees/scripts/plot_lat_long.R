#Script making plots trees functional indices depending on the latitude
#Tali Guez
#Trees Global and AF

library(tidyverse)
rm(list = ls())

#### GLOBAL TREES - INITIALIZATION ####
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/raw")
coord<-read_csv("grid_coordinates.csv")

setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results")
eco_based<-read_csv("REV_obs_results_ANGIO_eco_based_200.csv")
random1<-read_csv("REV_obs_results_ANGIO_random_1_200.csv")
random2<-read_csv("REV_obs_results_ANGIO_random_2_200.csv")
random3<-read_csv("REV_obs_results_ANGIO_random_3_200.csv")
random4<-read_csv("REV_obs_results_ANGIO_random_4_200.csv")
random5<-read_csv("REV_obs_results_ANGIO_random_5_200.csv")
random6<-read_csv("REV_obs_results_ANGIO_random_6_200.csv")
random7<-read_csv("REV_obs_results_ANGIO_random_7_200.csv")
random8<-read_csv("REV_obs_results_ANGIO_random_8_200.csv")
random9<-read_csv("REV_obs_results_ANGIO_random_9_200.csv")
random10<-read_csv("REV_obs_results_ANGIO_random_10_200.csv")
pca<-read_csv("REV_obs_results_ANGIO_pca_based_200.csv")
clusters<-read_csv("REV_obs_results_ANGIO_cluster_random_200.csv")

#Combine all datasets into one and add a 'source' column
combined_data <- bind_rows(
  eco_based %>% mutate(source = "Ecologically based"),
  random1 %>% mutate(source = "Random 1"),
  random2 %>% mutate(source = "Random 2"),
  random3 %>% mutate(source = "Random 3"),
  random4 %>% mutate(source = "Random 4"),
  random5 %>% mutate(source = "Random 5"),
  random6 %>% mutate(source = "Random 6"),
  random7 %>% mutate(source = "Random 7"),
  random8 %>% mutate(source = "Random 8"),
  random9 %>% mutate(source = "Random 9"),
  random10 %>% mutate(source = "Random 10"),
  pca %>% mutate(source = "PCA"),
  clusters %>% mutate(source = "Within clusters")
) %>%
  # Add NA to raoq, fdr, and fd if n < 9
  mutate(
    raoq = ifelse(n < 9, NA, raoq),
    fdr = ifelse(n < 9, NA, fdr),
    fd = ifelse(n < 9, NA, fd)
  )

# Merge with coordinates
combined_data <- left_join(combined_data, coord, by = "grid_id")

#Sampling for the plots points
set.seed(123)
sampled_points <- combined_data %>% 
  group_by(source) %>%
  sample_frac(0.1) #10% of the points represented 

##### Latitude ######
###### fdr~lat ######
#1. Trees - FDR ~ Latitude (no subset)
fdr_lat <- ggplot(combined_data, aes(x = Latitude, y = fdr, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - Functional Richness ~ Latitude") +
  theme(legend.position = "bottom")

fdr_lat

#2. Trees - FDR ~ Absolute Latitude (subset)
fdr_lat <- ggplot() +
  geom_point(data=sampled_points, aes(x = abs(Latitude), y = fdr, color = source),alpha = 0.03) +
  geom_smooth(data=combined_data, aes(x = abs(Latitude), y = fdr, color = source), method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets" )  +
  theme_bw() +
  labs(x = "Absolute latitude") +
  coord_cartesian(xlim=c(0,70), ylim=c(0,3.5), expand=FALSE)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=15.3),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.15, "cm")) #+
 # ggtitle("Trees - Functional Richness ~ Absolute Latitude") #

fdr_lat

###### raoq~lat ####
#1.Trees - RaoQ ~ Latitude (no subset) 
raoq_lat<-ggplot(combined_data, aes(x = Latitude, y = raoq, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "RaoQ") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - RaoQ ~ Latitude") +
  theme(legend.position = "bottom")

raoq_lat

#2.  Trees - Raoq ~ Absolute Latitude (subset)
raoq_lat <- ggplot() +
  geom_point(data=sampled_points, aes(x = abs(Latitude), y = raoq, color = source),alpha = 0.03) +
  geom_smooth(data=combined_data, aes(x = abs(Latitude), y = raoq, color = source), method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )  +
  theme_bw() +
  labs(x = "Absolute latitude")+
  #coord_cartesian(xlim=c(0,70), ylim=c(0,30), expand=FALSE)+ #first zoom
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
   #     axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=15.3),
  #      axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.15, "cm"))
  coord_cartesian(xlim=c(10,40), ylim=c(10,30), expand=FALSE) + #choose to make the zoom better
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position ="none" ,
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=20),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.3, "cm"), plot.margin=margin(r=15, t=10)) #+
  #ggtitle("Trees - RaoQ ~ Latitude") 

raoq_lat

#plot.margin = margin(t = 20,  # Top margin
 #                    r = 50,  # Right margin
  #                   b = 40,  # Bottom margin
  #                   l = 10)) #
#theme_get()$plot.margin 
#5.5points 5.5points 5.5points 5.5points

##### Longitude - NOT USED ######
###### fdr~long #####
fdr_long <- ggplot(combined_data, aes(x = Longitude, y = fdr, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - Functional Richness ~ Longitude") +
  theme(legend.position = "bottom")

fdr_long

###### raoq~long ####
raoq_long <- ggplot(combined_data, aes(x = Longitude, y = raoq, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "RaoQ") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - RaoQ ~ Longitude") +
  theme(legend.position = "bottom")

raoq_long


#### ATLANTIC FOREST TREES - INITIALIZATION ####
library(tidyverse)
rm(list = ls())
setwd("~/Master/M2/Internship_M2/analyse/trees/datas/raw")
coord<-read_csv("grid_coordinates.csv")

setwd("~/Master/M2/Internship_M2/analyse/trees/datas/clean/results_AF")
eco_based<-read_csv("REV_obs_results_ANGIO_AF_eco_based.csv")
random1<-read_csv("REV_obs_results_ANGIO_AF_random1.csv")
random2<-read_csv("REV_obs_results_ANGIO_AF_random2.csv")
random3<-read_csv("REV_obs_results_ANGIO_AF_random3.csv")
random4<-read_csv("REV_obs_results_ANGIO_AF_random4.csv")
random5<-read_csv("REV_obs_results_ANGIO_AF_random5.csv")
random6<-read_csv("REV_obs_results_ANGIO_AF_random6.csv")
random7<-read_csv("REV_obs_results_ANGIO_AF_random7.csv")
random8<-read_csv("REV_obs_results_ANGIO_AF_random8.csv")
random9<-read_csv("REV_obs_results_ANGIO_AF_random9.csv")
random10<-read_csv("REV_obs_results_ANGIO_AF_random10.csv")
pca<-read_csv("REV_obs_results_ANGIO_AF_pca_based.csv")
clusters<-read_csv("REV_obs_results_ANGIO_AF_cluster_random.csv")

#Combine all datasets into one and add a 'source' column
combined_data <- bind_rows(
  eco_based %>% mutate(source = "Ecologically based"),
  random1 %>% mutate(source = "Random 1"),
  random2 %>% mutate(source = "Random 2"),
  random3 %>% mutate(source = "Random 3"),
  random4 %>% mutate(source = "Random 4"),
  random5 %>% mutate(source = "Random 5"),
  random6 %>% mutate(source = "Random 6"),
  random7 %>% mutate(source = "Random 7"),
  random8 %>% mutate(source = "Random 8"),
  random9 %>% mutate(source = "Random 9"),
  random10 %>% mutate(source = "Random 10"),
  pca %>% mutate(source = "PCA"),
  clusters %>% mutate(source = "Within clusters")
) %>%
  # Add NA to raoq, fdr, and fd if n < 9
  mutate(
    raoq = ifelse(n < 9, NA, raoq),
    fdr = ifelse(n < 9, NA, fdr),
    fd = ifelse(n < 9, NA, fd)
  )

# Merge with coordinates
combined_data <- left_join(combined_data, coord, by = "grid_id")

#Sampling for the plots points
set.seed(123)
sampled_points <- combined_data %>% 
  group_by(source) %>%
  sample_frac(0.1) #10% of the points represented 

##### Latitude ######
###### fdr~lat ######
#1. AF - FDR ~ Latitude (no subset)
fdr_lat <- ggplot(combined_data, aes(x = Latitude, y = fdr, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_minimal() +
  ggtitle("Trees - Functional Richness ~ Latitude") +
  theme(legend.position = "bottom")

fdr_lat

#2. AF - FDR ~ Absolute latitude (no subset) #NON SENS
fdr_lat <- ggplot(data=combined_data, aes(x = abs(Latitude), y = fdr, color = source)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees AF - Functional Richness ~ Latitude") +
  theme(legend.position = "bottom")+
  labs(x = "Absolute latitude")

fdr_lat

#3. AF - FDR ~ Latitude (no subset) less transparent
fdr_lat <- ggplot(data=combined_data, aes(x = Latitude, y = fdr, color = source)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(x = "Latitude")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=20),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.3, "cm"), plot.margin=margin(r=15, t=10))+
  coord_cartesian(expand=FALSE) #choose to make the zoom better
  #coord_cartesian(xlim=c(10,40), ylim=c(10,30), expand=FALSE) #choose to make the zoom better
  #coord_cartesian(xlim=c(0,70), ylim=c(0,40), expand=FALSE) #first zoom 
  #ggtitle("Trees AF - Functional Richness ~ Latitude")
  

fdr_lat

###### raoq~lat ####
#1. AF - RaoQ ~ Latitude (no subset)
raoq_lat <- ggplot(data=combined_data, aes(x = Latitude, y = raoq, color = source)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=20),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.3, "cm"), plot.margin=margin(r=15, t=10))+
  labs(x = "Latitude")+
  coord_cartesian(expand=FALSE)#+
  #ggtitle("Trees AF - RaoQ ~ Latitude")

raoq_lat

##### Longitude ######
###### fdr~long #####
#1. AF - FDR ~ Longitude (no subset)
fdr_long <- ggplot(combined_data, aes(x = Longitude, y = fdr, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - Functional Richness ~ Longitude") +
  theme(legend.position = "bottom")

fdr_long

#2. AF - FDR ~ Longitude (no subset) less transparent
fdr_long <- ggplot(data=combined_data, aes(x = Longitude, y = fdr, color = source)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Richness (FDR)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=20),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.3, "cm"), plot.margin=margin(r=15, t=10))+
  labs(x = "Longitude")+
  coord_cartesian(ylim=c(2,3.5),expand=FALSE) #
  #ggtitle("Trees AF - Functional Richness ~ Longitude") 

fdr_long

###### raoq~long ####
#1. AF - RaoQ ~ Longitude (no subset)
raoq_long <- ggplot(combined_data, aes(x = Longitude, y = raoq, color = source)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "RaoQ") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("Trees - RaoQ ~ Longitude") +
  theme(legend.position = "bottom")

raoq_long

#2. AF - RaoQ ~ Longitude (no subset) less transparent
raoq_long <- ggplot(data=combined_data, aes(x = Longitude, y = raoq, color = source)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_continuous(name = "Functional Dispersion (RaoQ)") +
  scale_color_manual(
    values = c("Ecologically based" = "darkgreen",
               "Random 1" = "red",
               "Random 2" = "green",
               "Random 3" = "blue",
               "Random 4" = "pink",
               "Random 5" = "brown",
               "Random 6" = "grey",
               "Random 7" = "black",
               "Random 8" = "royalblue1",
               "Random 9" = "orange",
               "Random 10" = "yellow",
               "PCA" = "purple",
               "Within clusters" = "aquamarine"),
    name = "Datasets"
  )+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",
        axis.line=element_line(size=1.5, linewidth ="solid"), axis.text=element_text(color="black",size=20),
        axis.ticks=element_line(size=1.5), axis.ticks.length = unit(.3, "cm"), plot.margin=margin(r=15, t=10))+
  labs(x = "Longitude")+
  coord_cartesian(ylim=c(10,26),expand=FALSE) #
  #ggtitle("Trees AF - RaoQ ~ Longitude") 

raoq_long
