require(data.table)
require(dplyr)
require(lubridate)
A <- fread(file = "Data/2023-03-28_Mammals_abundance_matrix.csv")
A[,good_date:=as_date(Start.date)]
A[,td:=good_date-min(good_date)]
A[,td_sc:=scale(td)]

# explore raw data----
abund <- A[,`Canis aureus`:`Vulpes vulpes`] %>% unlist()
plot(sort(abund))
abund[abund>150]
A[`Canis aureus`>150,]
A[Site=="Kisalon" & Proximity=="Near",.(Project.ID,`Canis aureus`)] # this needs to checked further
A[Site=="Ofer" & Proximity=="Near",.(Project.ID,`Canis aureus`)]  # seems reasonable

A[`Sus scrofa`>200,]
A[Site=="Goren" & Proximity=="Near",.(Project.ID,Proximity,`Sus scrofa`)][order(Project.ID,Proximity)] # seems reasonable
A[Site=="Ein Yaakov" & Proximity=="Near",.(Project.ID,Proximity,`Sus scrofa`)][order(Project.ID,Proximity)] # need to check
A[Site=="Kfar Shamai" & Proximity=="Near",.(Project.ID,Proximity,`Sus scrofa`)][order(Project.ID,Proximity)] # reasonable

# explore dates
