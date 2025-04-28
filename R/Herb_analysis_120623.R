
save(all.data.mammals.final.final.without.rare,Herb.abu.data,Herb.mammals.abundance.matrix,Herb_abundmat,Herb_spp,Herb_abu,env_data.herb,Rare.sp_by.two.conditions,Herb_abu_no.wolf,Herb_spp2,mva_herb_nb_no.wolf, anov_mva_herb_nb_no.wolf,coefp_Herb_no.wolf,Herb,Herb.data,fox.herb,glm.fox.herb,file="Herb_analysis_2023.RData")


#Load packages
library(devtools)
library(rlang)
library(ggplot2)
library(data.table)
library(dplyr)
library(mvabund)
library(kableExtra)
library(lme4)
library(MASS)
library(vegan)
library(jtools)
library(interactions)
library(car)
library(ecoCopula)
library(chron)
library(data.table)
library(readxl)
library(performance)
library(glmmTMB) 
library(emmeans)
library(Cairo)
library(extrafont)
library(knitr)
library(chron)

#Load function
source("functions/plot_model_effect_mammals.R")
source("functions/plot_interaction_mammals.R")

#data exploration:
#find non- regular values with abundance:
D <- Herb_abu %>% as.data.table()
D$TotAb <- rowSums(D[,1:7])
D$TotRi <- rowSums(D[,1:7]>0)
View(D)

my_cols <- c("blue", "red")
dotchart(D$TotAb, ylab = "Site",
         groups = env_data.herb$Site, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")


my_cols <- c("blue", "red","darkgreen","brown","purple")
dotchart(D$TotAb, ylab = "Campaign",gcolor = my_cols,
         groups = env_data.herb$Project.ID, 
         color = my_cols[env_data.herb$Project.ID],
         cex = 0.9,  pch = 1, xlab = "abundance")

#There are exception in T3: karei deshe, because of high number of boars, which was changed,
#but still very high. 
#need to think if we need to remove them from the analysis.


Da <- Herb.mammals.abundance.matrix %>% as.data.table()
Da$TotAb <- rowSums(Da[,18:24])


kable(summary(Da[,.(TotAb,Project.ID, Settlements, Site,rescaled_Time.Diff_new, cosinus_Monitoring.Time.Diff, sinus_Monitoring.Time.Diff,Transect, Distance, Distance.Agri)]),"pipe")

pairs(Da[,lapply(X = .SD,FUN = as.numeric),.SDcols=c("Settlements","Project.ID","Site","rescaled_Time.Diff_new", "sinus_Monitoring.Time.Diff", "cosinus_Monitoring.Time.Diff", "Transect", "Distance","Distance.Agri")])

#seems that Project.ID and rescaled.time.diff has strong correaltion, obviously.
#Site and transect has strong correlation

kable(cor(Da[,lapply(X = .SD,FUN = as.numeric),.SDcols=c("Settlements","Project.ID","Site","rescaled_Time.Diff_new", "sinus_Monitoring.Time.Diff", "cosinus_Monitoring.Time.Diff", "Transect", "Distance","Distance.Agri")]),"pipe")

#the distance has some correlation with site (-0.7), and transect (-0.61/-0.68), and with distance.Agri (-0.82).
#distance.agri has correlation with site (-0.82) and transect (-0.79/-0.82) and with distance (0.75).
#it might be problematic to put site and distance together in the model, because of the high correlation.

# Many glm

#Make Site as factor
env_data.herb <- env_data.herb[,Site:=factor(Site)]

#Factor levels
unique(env_data.herb$Site)

#Compare manyglm model distributions
mva_herb_poiss <- manyglm(Herb_spp~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                      sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff , family = "poisson", data=env_data.herb)

plot.manyglm(mva_herb_poiss, which = 1:3)

mva_herb_nb <- manyglm(Herb_spp~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                            sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.herb)

plot.manyglm(mva_herb_nb, which = 1:3)

#mva_herb_nb looks much better.

#model selection:

drop1(mva_herb_nb)


#nothing has to be dropped. 

#the final model is mva_herb_nb

anov_mva_herb_nb <- anova(mva_herb_nb, p.uni = "adjusted")
print(anov_mva_herb_nb)

#graph:
coefp_Herb <- merge(data.table(t(coef(mva_herb_nb)),keep.rownames=TRUE), data.table(t(anov_mva_herb_nb$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_Herb <- merge(coefp_Herb,as.data.table(colSums(Herb_spp),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_Herb)

colnames(coefp_Herb) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_Karei.Deshe.coef","Site_Natur.coef","Site_Shaal.coef","Site_Yavneel.coef", "sinus_monitoring.time.coef","cosinus_monitoring.time.coef","Distance_Time.coef", "Intercept.p","Distance.p", "Time_diff.p","Site.p", "sin_time.p","Cosin_time.p","Distance_Time.p" ,"species_abundance")


#add hebrew name
A <- as.data.table(all.data.mammals.final.final.without.rare)
sci_heb <- A[,.(Species_Heb=unique(Species_Heb)), keyby=Species]
sci_heb[,Species:=make.names(Species)]
coefp_Herb <- merge(coefp_Herb, sci_heb, by.x = "SciName", by.y = "Species")

write.csv(coefp_Herb, "coefficients_Herb.csv")

#only golden jackal, grey wolf and red fox are affected by the parameters in the model.



golden.jackal.herb <- Herb_abundmat[,c(1:18,25:26)]

names(golden.jackal.herb)[names(golden.jackal.herb) == 'Canis aureus'] <- 'Canis.aureus'

#only distance is significant, so I run model without intercation:


glm.golden.jackal.herb <- glm.nb(Canis.aureus ~ Distance_rescaled+rescaled_Time.Diff_new+Site + 
                                   sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff, data = golden.jackal.herb, maxit = 999)


summary(glm.golden.jackal.herb)


#no dots:
effect_plot_g.Jackal.herb <- effect_plot(model = glm.golden.jackal.herb, data=golden.jackal.herb, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
effect_plot_g.Jackal.herb+ylim(0,200)

#Change over distance:

#Mean distance rescaled herb
# mean_distance_rescaled <- mean(Herb.mammals.abundance.matrix$Distance)
# sd_distance_rescaled <- sd(Herb.mammals.abundance.matrix$Distance)
# 
# # Calculate for 100 and 2000 meters
# (100-mean_distance_rescaled)/sd_distance_rescaled
# (2000-mean_distance_rescaled)/sd_distance_rescaled
# 
# # Calculate for 500, 1000, 1500
# (500-mean_distance_rescaled)/sd_distance_rescaled
# (1000-mean_distance_rescaled)/sd_distance_rescaled
# (1500-mean_distance_rescaled)/sd_distance_rescaled
# 
# jackal.distance <- data.frame(Site = "Gamla",
#                            sinus_Monitoring.Time.Diff = 0.1007873, 
#                            cosinus_Monitoring.Time.Diff = -0.05638926, 
#                            rescaled_Time.Diff_new = 0.1421484, 
#                            Distance_rescaled = c(-1.455879, -1.132998, 
#                                                  -0.7293972, -0.3257959,
#                                                  0.07780538
#                            ))
# 
# predict(glm.golden.jackal.herb, newdata = jackal.distance,
#         type = "response")
# #Decrease from 100 to 2000 m of 543.06% or 55.3 fold
# (3983.50352/72.02639) - 1  
# 
# 3983.50352/72.02639
# 
# 72.02639/207.06964
# 207.06964/595.30730

my_cols <- c("blue", "red")
dotchart(golden.jackal.herb$Canis.aureus, ylab = "Site",
         groups = golden.jackal.herb$Site, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")

#gcolor = my_cols,


my_cols <- c("blue", "red","darkgreen","brown","purple")
dotchart(golden.jackal.herb$Canis.aureus, ylab = "Campaign",gcolor = my_cols,
         groups = golden.jackal.herb$Project.ID, 
         color = my_cols[golden.jackal.herb$Project.ID],
         cex = 0.9,  pch = 1, xlab = "abundance")

#find how many jackals are in each site and each campaign:
jackal <- as.data.table(golden.jackal.herb)
jackal.counts <- jackal[,.(total.count = sum(Canis.aureus)), keyby=c("Project.ID","Site")]


p <- ggplot(golden.jackal.herb, aes(x=Proje, y=Canis.aureus))+ geom_point(aes(color=Site),size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

write.csv(jackal.counts, "jackal.counts_Herb.csv")

#there are some sites with less jackal obs, but not completely zeros. So, I do not know if there is
#justification to remove sites. 

p7 <- ggplot(golden.jackal.herb, aes(x=Project.ID, y=Canis.aureus))+ geom_point(aes(color=Site),size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# scale_y_continuous(limits = c(0,100))

p7

p5 <- ggplot(golden.jackal.herb, aes(fill=Site,y=Canis.aureus,x=Project.ID))+ geom_bar(position="dodge", stat="identity")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# scale_y_continuous(limits = c(0,100))

p5

rare.species.wolf <- Rare.sp_by.two.conditions[Rare.sp_by.two.conditions$Species=="Canis lupus",]


write.csv(rare.species.wolf, "rare.species.wolf.csv")

#grey wolf:
wolf.herb <- Herb_abundmat[,c(1:17,19,25:26)]

unique(wolf.herb$Site)

names(wolf.herb)[names(wolf.herb) == 'Canis lupus'] <- 'Canis.lupus'

#only the interaction between distance and time is significant:


glm.wolf.herb <- glm.nb(Canis.lupus ~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                                   sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff, data = wolf.herb, maxit = 999)


summary(glm.wolf.herb)

#the estimate numbers are really low.

my_cols <- c("blue", "red")
dotchart(wolf.herb$Canis.lupus, ylab = "Site",
         groups = wolf.herb$Site, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")

#gcolor = my_cols,


my_cols <- c("blue", "red","darkgreen","brown","purple")
dotchart(wolf.herb$Canis.lupus, ylab = "Campaign",gcolor = my_cols,
         groups = wolf.herb$Project.ID, 
         color = my_cols[wolf.herb$Project.ID],
         cex = 0.9,  pch = 1, xlab = "abundance")

#remove Karei Deshe:
wolf.herb_no.KD <- wolf.herb[!Site %in% "Karei Deshe",]

wolf.herb_no.KD <- wolf.herb_no.KD[,Site:=factor(Site)]

unique(wolf.herb_no.KD$Site)

glm.wolf.herb2 <- glm.nb(Canis.lupus ~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                          sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff, data = wolf.herb_no.KD, maxit = 999)


summary(glm.wolf.herb2)

#the estimates are better.

wolf_interaction <- interact_plot(model = glm.wolf.herb2, data=wolf.herb_no.KD, pred = rescaled_Time.Diff_new, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus")
wolf_interaction 
#no interval:
interact_plot(model = glm.wolf.herb2, data=wolf.herb_no.KD, pred = rescaled_Time.Diff_new, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = F,modx.values = "plus-minus")

#for wolf- try to remove Karei Deshe and Natur too:

wolf.herb_thin <- subset(wolf.herb, Site!="Karei Deshe" & Site!="Natur")

wolf.herb_thin <- as.data.table(wolf.herb_thin)

wolf.herb_thin <- wolf.herb_thin[,Site:=factor(Site)]

unique(wolf.herb_thin$Site)

glm.wolf.herb3 <- glm.nb(Canis.lupus ~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                           sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff, data = wolf.herb_thin, maxit = 999)

summary(glm.wolf.herb3)

wolf.herb_thin$Site <- relevel(wolf.herb_thin$Site,ref="Yavneel")


interact_plot(model = glm.wolf.herb3, data=wolf.herb_thin, pred = rescaled_Time.Diff_new, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus")

#no interval:
wolf_inetraction_plot <- interact_plot(model = glm.wolf.herb3, data=wolf.herb_thin, pred = rescaled_Time.Diff_new, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = F,modx.values = "plus-minus")

unique(wolf.herb$rescaled_Time.Diff_new)
#red fox:
fox.herb <- Herb_abundmat[,c(1:17,24:26)]

unique(fox.herb$Site)

names(fox.herb)[names(fox.herb) == 'Vulpes vulpes'] <- 'Vulpes.vulpes'

#the time and distance are significant separately, so model without interaction:

glm.fox.herb <- glm.nb(Vulpes.vulpes ~ Distance_rescaled+rescaled_Time.Diff_new+Site + 
                          sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff, data = fox.herb, maxit = 999)

summary(glm.fox.herb)


#Calculate increase from July 1st 2012 to July 1st 2022

# fox.time <- data.frame(Site = "Gamla",
#                        sinus_Monitoring.Time.Diff = 0.1007873, 
#                        cosinus_Monitoring.Time.Diff = -0.05638926, 
#                        Distance_rescaled = 2.522537e-17, 
#                        rescaled_Time.Diff_new = c(-1.129633336,
#                                                   1.218933532
#                        ))
#                        
# predict(glm.fox.herb, newdata = fox.time,
#         type = "response")
# 
# 
# # Increase of 260.5% between 2014 - 2022
# (2.2065974/0.6710217) - 1
# 
# #or 3.28 fold increase
# 
# (2.2065974/0.6710217) 



#Change over distance:

#Mean distance rescaled herb
# mean_distance_rescaled <- mean(Herb.mammals.abundance.matrix$Distance)
# sd_distance_rescaled <- sd(Herb.mammals.abundance.matrix$Distance)
# 
# # Calculate for 100 and 2000 meters
# (100-mean_distance_rescaled)/sd_distance_rescaled
# (500-mean_distance_rescaled)/sd_distance_rescaled
# (1000-mean_distance_rescaled)/sd_distance_rescaled
# (1500-mean_distance_rescaled)/sd_distance_rescaled
# (2000-mean_distance_rescaled)/sd_distance_rescaled
# 
# fox.distance <- data.frame(Site = "Gamla",
#                            sinus_Monitoring.Time.Diff = 0.1007873, 
#                            cosinus_Monitoring.Time.Diff = -0.05638926, 
#                            rescaled_Time.Diff_new = 0.1421484, 
#                            Distance_rescaled = c(-1.455879, -1.132998, -0.7293972,
#                                                  -0.3257959, 0.07780538))
# 
# predict(glm.fox.herb, newdata = fox.distance,
#         type = "response")
# 
# #Decrease from 100 to 2000 m of 133.45% or 2.33 fold
# (2.859044/1.224651) - 1
# (2.859044/1.224651)
# 
# # decrease of 80% every 500 m
# 
# 1.224651 / 1.530768
# 1.530768 / 1.913403

#for distance:
effect_plot_fox.herb <- effect_plot(model = glm.fox.herb, data=fox.herb, 
    pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), 
    point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
effect_plot_fox.herb+ylim(0,40)


#time:
effect_plot_fox.herb2 <- effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = T, plot.points = T, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
effect_plot_fox.herb+ylim(0,40)

#no dots:
effect_plot_fox.herb2 <- effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

effect_plot_fox.herb2 <- effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = F, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

save(all.data.mammals,all.data.mammals_below.10th.day,all.data.mammals.final,all.data.mammals.final_duplicates,all.data.mammals.final_Inland.sands_new,all.data.mammals.final.final,corrected.sites,corrected.sites_below.10th.day,counting.duplicates,date.per.Unit.Site_corrected,deployments_Negev_Traznisition_T4,Duplicate_obs_new,Duplicate_obs,duplicate_obs_all_by.deployment,list_of_duplicates,list_of_duplicates_prop,list_of_duplicates_prop_sum,list_of_duplicates_prop_freq1, list_of_duplicates_proportion,list_of_duplicates_prop_sum_all, list_of_duplicates_sum,list_of_duplicates_proportion_by.count,list_of_transects,species_Latin_Heb_names,T4.deployments.data,T4.mammals_below.10th.day,T4.mammals_below.10th.day_Negev_Transition,T4.mammals.data,T4.mammals.data_Negev_Transition,T4.mammals.data_seq.correction,T4.mammals.data_seq.correction_Negev_Transition,T4.mammals.data.combined, T4.mammals.data.combined_Negev_Transition,T4.mammals.data.new,T4.mammals.data.new_Negev_Transition,Unit_per_year_corrected,Unit.site_per_project,Unit.site_per_year_corrected,Unit.site_per_year_Negev_Transition_T4,Time,Time.diff,Time_diff_below_30.min,Time_diff_proportion,Time_diff_proportion_new,all.data.mammals_cast,mammals.abundance_matrix,all.data.mammals_cast_2,mammals.indcidence.matrix,mammals_indcidence_matrix,Rare,rare.species,Maquis.data,Maquis_total_counts_site_campaign,Maquis_total_counts,Maquis_total_counts_campaign_all_sp,Maquis_total_counts_campaign_proximity_all_sp,mammals.abundance.matrix,sp.richness.site,all.data.mammals.final.final.without.rare,Rare.campaign,rare.species.per.campaign,Abund_total_counts_campaign,sp.richness.Maquis,Maquis.data,Rare_by_count,Rare_by_count_deployment,Rare.sp_by.two.conditions,all.data.mammals.final.final.for.richness,sp.richness.Maquis.per.cam,Maquis.richness.matrix,Maquis.incidence.long,list_of_deployments,list_of_existing_deployments,Distance_data, distance_data_T4, distance_data_T0_T3,Agri_Distance_data,Maquis.abu.data,env_data.maquis,Maquis_abu,Maquis_abu2,site.scrs,spp.scrs,nMDS,nMDS2, deployments_Yatir_Mirsham, T4.mammals.data.new_yatir_mirsham,T4.mammals.data_yatir_mirsham, T4.mammals.data.combined_yatir_mirsham, T4.mammals_below.10th.day_yatir_mirsham,Maquis_spp,Maquis.mammals.abundance.matrix,nMDS4,Maquis_abu_without.T2.GorenFar7,Maquis_abu_without.T2.GorenNear4,env_data.maquis_without.T2.GorenFar7,env_data.maquis_without.T2.GorenNear4, Herb.data,Herb.mammals.abundance.matrix,Herb_abu,Herb_abundmat,Herb_spp,sp.richness.Herb.per.cam,Herb.abu.data,env_data.herb,Distance_diff,mvabund4,mvabund43,anov_mvabund43,mvabund53,mvabund5, mvabund6, mvabund63,anov_mvabund63, coefp,Hyena_Maquis,Hyena_Maquis_per.unit,Maquis.sample.size.per.species,Maquis.sample.size.per.species.no.campaign, model.coef, Maquis.cam.number,Maquis.cam.number.per.campaign,Maquis.cam.number.per.campaign.per.species,Maquis.cam.number.per.species,Jackal_Maquis,Jackal_Maquis_per.unit, coefp_Herb,coefp_Herb2,envd.galil,envd.jud,galil_spp,subunit_spp,galilee.abundance.matrix,galil_abu,jud_abu,judea.abundance.matrix,mva.galil2,mva.jud2,mva0.jud,mva1.jud,mva1.jud2,mva0, anov_mva.jud2, mva.carmel2,jud_abundmat,coefp_jud,model.coef_jud,carmel.abundance.matrix,carmel_abundmat,carmel_spp,envd.carmel,carmel_abu,anov_mva.carmel2, coefp_carmel, mva.carmel2, Maquis_abundmat,maqi_mvabund_nb2,maqi_mvabund_nb21,anov_maqi_mvabund21,anov_mva.jud.nb1,mva.jud.nb1,coefp_maqi,coefp_maqi_judea,anov_mva.galil2,coefp_maqi_galil,coefp_maqi_carmel,mva.carmel2,anov_mva.carmel2, anov_mva_herb_nb,mva_herb_nb,coefp_Herb, file="Mammals_analysis_2023.RData")

effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = T, plot.points = T, jitter = c(0.1,0), 
            partial.residuals = F, colors = "Qual1") + scale_x_continuous(breaks=c(-1.233466888,-0.601840099,0.028922631,0.66054942,1.291312151), 
                                                                                       labels=c("June 2014", "June 2016", "June 2018", "June 2020", "June 2022"), name = "Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylab("Red fox abundance")

#run the mvabund without wolf (37 indvs and 14 cameras):

Herb_abu_no.wolf <- Herb_abundmat[,c(18,20:24)]
Herb_spp2 <- mvabund(Herb_abu_no.wolf)

mva_herb_nb_no.wolf <- manyglm(Herb_spp2~ Distance_rescaled*rescaled_Time.Diff_new+Site + 
                         sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.herb)

plot.manyglm(mva_herb_nb_no.wolf, which = 1:3)

#model selection:

drop1(mva_herb_nb_no.wolf)

 
#nothing has to be dropped. 

#the final model is mva_herb_nb_no.wolf

anov_mva_herb_nb_no.wolf <- anova(mva_herb_nb_no.wolf, p.uni = "adjusted")
print(anov_mva_herb_nb_no.wolf)

coefp_Herb_no.wolf <- merge(data.table(t(coef(mva_herb_nb_no.wolf)),keep.rownames=TRUE), data.table(t(anov_mva_herb_nb_no.wolf$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_Herb_no.wolf <- merge(coefp_Herb_no.wolf,as.data.table(colSums(Herb_spp2),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_Herb_no.wolf)

colnames(coefp_Herb_no.wolf) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_Karei.Deshe.coef","Site_Natur.coef","Site_Shaal.coef","Site_Yavneel.coef", "sinus_monitoring.time.coef","cosinus_monitoring.time.coef","Distance_Time.coef", "Intercept.p","Distance.p", "Time_diff.p","Site.p", "sin_time.p","Cosin_time.p","Distance_Time.p" ,"species_abundance")

write.csv(coefp_Herb_no.wolf, "coefficients_Herb_no.wolf.csv")

Herb <- all.data.mammals.final.final.for.richness[all.data.mammals.final.final.for.richness$Unit=="Herbaceous and Dwarf Shrub Vegetation",]

unique(Herb$Species)

#to get the numbers of obs and inds:

Herb.data <- all.data.mammals.final.final.without.rare[all.data.mammals.final.final.without.rare$Unit=="Herbaceous and Dwarf Shrub Vegetation",]

unique(Herb.data$Species)

#remove canis lupus and hyena hyena:

library(dplyr)
Herb.data.no.rare <- Herb.data %>% filter(!Species %in% c('Canis lupus','Hyaena hyaena'))
unique(Herb.data.no.rare$Species)

sum(Herb.data.no.rare$Sum)

## Exporting Hebrew vectorized PDFs----

# The plot:

effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = T, plot.points = T, jitter = c(0.1,0), 
            partial.residuals = F, colors = "Qual1",point.size = 5,line.thickness = 2,point.alpha = 0.25,line.colors = "black") + scale_x_continuous(breaks=c(-1.233466888,-0.601840099,0.028922631,0.66054942,1.291312151), 
     labels=c("2014", "2016", "2018", "2020", "2022"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"))+
  ylab("עפש")+xlab("ןמז")



#Red fox plot using function


redfox_shrub <- effect_plot(model = glm.fox.herb, data=fox.herb, pred = rescaled_Time.Diff_new, interval = T, plot.points = T, jitter = c(0.1,0), 
                            partial.residuals = F, colors = "Qual1",point.size = 5,line.thickness = 2,point.alpha = 0.25,line.colors = "black") 
                                                                                                                                                                   
effplot_vulpes_time <- plot_model_effect(P.anal = fox.herb, 
                               m = glm.fox.herb, 
                               eff2plot = "rescaled_Time.Diff_new", 
                               plot_points=FALSE, plot_residuals=FALSE, 
                               export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                               fontsize=22, pdf_width=160, outpath = "output/batha/")

mean(all.data.mammals.final.final.without.rare$Time.Diff)
sd(all.data.mammals.final.final.without.rare$Time.Diff)
