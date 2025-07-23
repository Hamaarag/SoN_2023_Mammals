# this code is for exploring the threshold for minimum number of camera nights
# there is a minimum of 7 properly operating cameras per transect, under 7 the transect must be re-deployed.
# This code is meant to explore the number of properly operating cameras over a range of thresholds for the
# minimal number of nights with proper operation per camera (under which the camera is disqualified)
library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
# this is the metadata file exported from fulcrum
A <- read_excel(path = "c:/Users/Ron Chen 2/Desktop/tmp/20211122 mammals_2021 - fulcrum output.xlsx")
A$point_name <- gsub("Hashofat","Hashofet",A$point_name)
A$`10_days`[which(A$n_nights=="yes")] <- "yes"
A$n_nights[which(A$n_nights=="yes")] <- NA
A$n_nights <- as.numeric(A$n_nights)
A$include_camera = rep(TRUE,nrow(A),1)
A$transect_id <- paste(A$unit,str_match(A$point_name,"^([^0-9]+\\s(\\d_)?)\\d{1}$")[,2],sep = "_")
A <- A[-grep(pattern = "NA",x = A$transect_id),]
B_all <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
A$include_camera[A$n_nights<4] <- FALSE
B4 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p4 <- ggplot(B4,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p4
print(paste("number of transects with 7 cameras or more for a minimum of 4 nights of operation:",sum(B4$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 4 nights of operation:",sum(B4$tot_camera<7)))
A$include_camera[A$n_nights<5] <- FALSE
B5 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p5 <- ggplot(B5,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p5
print(paste("number of transects with 7 cameras or more for a minimum of 5 nights of operation:",sum(B5$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 5 nights of operation:",sum(B5$tot_camera<7)))
A$include_camera[A$n_nights<6] <- FALSE
B6 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p6 <- ggplot(B6,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p6
print(paste("number of transects with 7 cameras or more for a minimum of 6 nights of operation:",sum(B6$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 6 nights of operation:",sum(B6$tot_camera<7)))
A$include_camera[A$n_nights<7] <- FALSE
B7 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p7 <- ggplot(B7,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p7
print(paste("number of transects with 7 cameras or more for a minimum of 7 nights of operation:",sum(B7$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 7 nights of operation:",sum(B7$tot_camera<7)))
A$include_camera[A$n_nights<8] <- FALSE
B8 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p8 <- ggplot(B8,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p8
print(paste("number of transects with 7 cameras or more for a minimum of 8 nights of operation:",sum(B8$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 8 nights of operation:",sum(B8$tot_camera<7)))
A$include_camera[A$n_nights<9] <- FALSE
B9 <- A %>% group_by(transect_id) %>% summarise(tot_camera = sum(include_camera))
p9 <- ggplot(B9,aes(tot_camera)) + geom_histogram(binwidth = 1) + scale_x_continuous(0:9)
p9
print(paste("number of transects with 7 cameras or more for a minimum of 9 nights of operation:",sum(B9$tot_camera>=7)))
print(paste("number of transects with less than 7 cameras for a minimum of 9 nights of operation:",sum(B9$tot_camera<7)))
