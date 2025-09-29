
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
library(reshape2)

save(all.data.mammals.final.final.without.rare,Maquis.abu.data,Maquis.mammals.abundance.matrix,Maquis_abundmat,Maquis_abu,Maquis_spp,env_data.maquis,maqi_mvabund_nb21,anov_maqi_mvabund21,coefp_maqi,golden_jackal.maquis,glm.golden.jackal_maquis,gazelle.maquis,glm.gazelle_maquis2,hyena.maquis,glm.hyena_maquis,porcupine.maquis,glm.porcupine_maquis,badger.maquis,glm.badger_maquis,w.boar.maqi_thin,glm.w.boar_maquis2,fox.maquis,glm.fox_maquis,w.boar.jud,mva.jud.nb1,jud_abundmat,subunit_spp,envd.jud,w.boar.jud_thin2,glm.w.boar_jud3,envd.galil,galil_spp,mva.galil2,anov_mva.galil2,coefp_maqi_galil,galil_abundmat,golden_jackal.galil,glm.golden.jackal_galil, file="Maqi_analysis_2023.RData")

Maquis.abu.data <- all.data.mammals.final.final.without.rare[all.data.mammals.final.final.without.rare$Unit=="Mediterranean Maquis",]

#remove rare species of Maquis: According to "rare.species":
#Lepus capensis need to be removed:

Maquis.abu.data <- Maquis.abu.data %>% filter(!Species %in% c('Lepus capensis'))
unique(Maquis.abu.data$Species_Heb)

#need to remove Dama mesopotamica from the data:
#unique(Maquis.abu.data$Species)

Maquis.abu.data <- Maquis.abu.data %>% filter(!Species %in% c('Dama mesopotamica'))
unique(Maquis.abu.data$Species_Heb)


Maquis.mammals.abundance.matrix <- dcast(Maquis.abu.data,
                                         Project.ID+Deployment.id_new+Unit+Subunit+Site+Settlements+Lon+Lat+Start.date+Time.Diff+rescaled_Time.Diff+Monitoring.Time.Diff_numeric+sinus_Monitoring.Time.Diff+cosinus_Monitoring.Time.Diff+Distance+Transect+Distance.Agri ~ Species, value.var = 'Sum', fun = sum)


#remove duplicate rows:
#There a duplicate lines of "T1_kerem maharal_carmel_kerem maharal Far CA17" and CA18- need to check!
#line 332
#line 334
Maquis.mammals.abundance.matrix <- Maquis.mammals.abundance.matrix%>%  filter(!row_number() %in% c(332,334))

unique(Maquis.abu.data$Species)


Maquis_abundmat <- as.data.table(Maquis.mammals.abundance.matrix)
#Maquis_abu <- Maquis_abundmat[,`Canis aureus`:`Vulpes vulpes`]
#abu.maquis <- abundmat[grepl("Maquis",Unit),`Canis aureus`:`Vulpes vulpes`] %>% mvabund()
env_data.maquis <- Maquis_abundmat[,1:17]
#Maquis_spp <- mvabund(Maquis_abu)

write.csv(env_data.maquis, "env_data.maquis.csv")

mean(env_data.maquis$Time.Diff)
sd(env_data.maquis$Time.Diff)

#data exploration:
#find non- regular values with abundance:
D <- Maquis_abu %>% as.data.table()
D$TotAb <- rowSums(D[,1:8])
D$TotRi <- rowSums(D[,1:8]>0)


dotchart(D$TotAb, ylab = "Site",
         groups = env_data.maquis$Site, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")


my_cols <- c("blue", "red","darkgreen")
dotchart(D$TotAb, ylab = "Subunit",gcolor = my_cols,
         groups = env_data.maquis$Subunit, 
         color = my_cols[env_data.maquis$Subunit],
         cex = 0.9,  pch = 1, xlab = "abundance")

my_cols <- c("blue", "red","darkgreen","brown","purple")
dotchart(D$TotAb, ylab = "Campaign",gcolor = my_cols,
         groups = env_data.maquis$Project.ID, 
         color = my_cols[env_data.maquis$Project.ID],
         cex = 0.9,  pch = 1, xlab = "abundance")

#There are some exceptions, mostly in the Galilee, on T1,T4.
#need to think if we need to remove them from the analysis.

Da <- Maquis.mammals.abundance.matrix %>% as.data.table()
Da$TotAb <- rowSums(Da[,19:26])


kable(summary(Da[,.(TotAb,Project.ID, Settlements,Subunit, Site,rescaled_Time.Diff, cosinus_Monitoring.Time.Diff, sinus_Monitoring.Time.Diff,Transect, Distance, Distance.Agri)]),"pipe")

pairs(Da[,lapply(X = .SD,FUN = as.numeric),.SDcols=c("Settlements","Project.ID","Subunit","Site","rescaled_Time.Diff", "sinus_Monitoring.Time.Diff", "cosinus_Monitoring.Time.Diff", "Transect", "Distance", "Distance.Agri")])

#seems that Project.ID and rescaled.time.diff has strong correaltion, obviously.
#Site and transect has strong correlation
#Distance and distance.agri does not have correlation.
#Distance seems strange, mainly because Kerem Maharal Far, which is 2 km far from settle, and some far transects in other sites 
#are 200-300m from settle. 

kable(cor(Da[,lapply(X = .SD,FUN = as.numeric),.SDcols=c("Settlements","Project.ID","Subunit","Site","rescaled_Time.Diff", "sinus_Monitoring.Time.Diff", "cosinus_Monitoring.Time.Diff", "Transect", "Distance","Distance.Agri")]),"pipe")
#Distance and Settlements has pretty strong correlation (-0.67), but we use only the distance parameter.
#distance and distance.agri does not have strong correlation (0.34), so can be both included in the model.

#with rescaling the distance:
env_data.maquis$Distance_rescaled = (env_data.maquis$Distance - mean(env_data.maquis$Distance))/sd(env_data.maquis$Distance)
env_data.maquis$Distance.Agri_rescaled = (env_data.maquis$Distance.Agri - mean(env_data.maquis$Distance.Agri))/sd(env_data.maquis$Distance.Agri)
env_data.maquis[,Subunit_rescaled:=scale(as.numeric(Subunit))]

env_data.maquis <- env_data.maquis[,Site:=factor(Site)]
unique(env_data.maquis$Site)


maqi_mvabund_poiss <- manyglm(Maquis_spp~Distance_rescaled*rescaled_Time.Diff+ Subunit+ Site+
                      sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff , family = "poisson", data=env_data.maquis)


maqi_mvabund_nb <- manyglm(Maquis_spp~Distance_rescaled*rescaled_Time.Diff+ Subunit+ Site+
                                sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.maquis)

plot.manyglm(maqi_mvabund_poiss, which = 1:3)
plot.manyglm(maqi_mvabund_nb, which = 1:3)

#the nb looks better. 

maqi_mvabund_nb1 <- manyglm(Maquis_spp~Distance_rescaled*rescaled_Time.Diff+ Subunit+ 
                             sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.maquis)

maqi_mvabund_nb2 <- manyglm(Maquis_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                             sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.maquis)

aic_dt <- data.table(nb1 = maqi_mvabund_nb1$aic, nb2=maqi_mvabund_nb2$aic)
colMeans(aic_dt)

#The AIC is better with site, not subunit.

#model selection:
drop1(maqi_mvabund_nb2)


#drop sinus:
maqi_mvabund_nb21 <- manyglm(Maquis_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                              cosinus_Monitoring.Time.Diff , family = "negative.binomial", data=env_data.maquis)

drop1(maqi_mvabund_nb21)

#no need to drop more parameters. 

#final model is maqi_mvabund_nb21

plot.manyglm(maqi_mvabund_nb21, which = 1:3)

save(all.data.mammals,all.data.mammals_below.10th.day,all.data.mammals.final,all.data.mammals.final_duplicates,all.data.mammals.final_Inland.sands_new,all.data.mammals.final.final,corrected.sites,corrected.sites_below.10th.day,counting.duplicates,date.per.Unit.Site_corrected,deployments_Negev_Traznisition_T4,Duplicate_obs_new,Duplicate_obs,duplicate_obs_all_by.deployment,list_of_duplicates,list_of_duplicates_prop,list_of_duplicates_prop_sum,list_of_duplicates_prop_freq1, list_of_duplicates_proportion,list_of_duplicates_prop_sum_all, list_of_duplicates_sum,list_of_duplicates_proportion_by.count,list_of_transects,species_Latin_Heb_names,T4.deployments.data,T4.mammals_below.10th.day,T4.mammals_below.10th.day_Negev_Transition,T4.mammals.data,T4.mammals.data_Negev_Transition,T4.mammals.data_seq.correction,T4.mammals.data_seq.correction_Negev_Transition,T4.mammals.data.combined, T4.mammals.data.combined_Negev_Transition,T4.mammals.data.new,T4.mammals.data.new_Negev_Transition,Unit_per_year_corrected,Unit.site_per_project,Unit.site_per_year_corrected,Unit.site_per_year_Negev_Transition_T4,Time,Time.diff,Time_diff_below_30.min,Time_diff_proportion,Time_diff_proportion_new,all.data.mammals_cast,mammals.abundance_matrix,all.data.mammals_cast_2,mammals.indcidence.matrix,mammals_indcidence_matrix,Rare,rare.species,Maquis.data,Maquis_total_counts_site_campaign,Maquis_total_counts,Maquis_total_counts_campaign_all_sp,Maquis_total_counts_campaign_proximity_all_sp,mammals.abundance.matrix,sp.richness.site,all.data.mammals.final.final.without.rare,Rare.campaign,rare.species.per.campaign,Abund_total_counts_campaign,sp.richness.Maquis,Maquis.data,Rare_by_count,Rare_by_count_deployment,Rare.sp_by.two.conditions,all.data.mammals.final.final.for.richness,sp.richness.Maquis.per.cam,Maquis.richness.matrix,Maquis.incidence.long,list_of_deployments,list_of_existing_deployments,Distance_data, distance_data_T4, distance_data_T0_T3,Agri_Distance_data,Maquis.abu.data,env_data.maquis,Maquis_abu,Maquis_abu2,site.scrs,spp.scrs,nMDS,nMDS2, deployments_Yatir_Mirsham, T4.mammals.data.new_yatir_mirsham,T4.mammals.data_yatir_mirsham, T4.mammals.data.combined_yatir_mirsham, T4.mammals_below.10th.day_yatir_mirsham,Maquis_spp,Maquis.mammals.abundance.matrix,nMDS4,Maquis_abu_without.T2.GorenFar7,Maquis_abu_without.T2.GorenNear4,env_data.maquis_without.T2.GorenFar7,env_data.maquis_without.T2.GorenNear4, Herb.data,Herb.mammals.abundance.matrix,Herb_abu,Herb_abundmat,Herb_spp,sp.richness.Herb.per.cam,Herb.abu.data,env_data.herb,Distance_diff,mvabund4,mvabund43,anov_mvabund43,mvabund53,mvabund5, mvabund6, mvabund63,anov_mvabund63, coefp,Hyena_Maquis,Hyena_Maquis_per.unit,Maquis.sample.size.per.species,Maquis.sample.size.per.species.no.campaign, model.coef, Maquis.cam.number,Maquis.cam.number.per.campaign,Maquis.cam.number.per.campaign.per.species,Maquis.cam.number.per.species,Jackal_Maquis,Jackal_Maquis_per.unit, coefp_Herb,coefp_Herb2,envd.galil,envd.jud,galil_spp,subunit_spp,galilee.abundance.matrix,galil_abu,jud_abu,judea.abundance.matrix,mva.galil2,mva.jud2,mva0.jud,mva1.jud,mva1.jud2,mva0, anov_mva.jud2, mva.carmel2,jud_abundmat,coefp_jud,model.coef_jud,carmel.abundance.matrix,carmel_abundmat,carmel_spp,envd.carmel,carmel_abu,anov_mva.carmel2, coefp_carmel, mva.carmel2, Maquis_abundmat,maqi_mvabund_nb2,maqi_mvabund_nb21,anov_maqi_mvabund21, file="Mammals_analysis_2023.RData")

anov_maqi_mvabund21 <- anova(maqi_mvabund_nb21, p.uni = "adjusted")
print(anov_maqi_mvabund21)

#Coefficients
coefp_maqi <- merge(data.table(t(coef(maqi_mvabund_nb21)),keep.rownames=TRUE), data.table(t(anov_maqi_mvabund21$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_maqi <- merge(coefp_maqi,as.data.table(colSums(Maquis_spp),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_maqi)
colnames(coefp_maqi) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_BeitOren.coef","Site_Ein Yaakov.coef","Site_GivatYearim.coef", "Site_GivatYeshayahu.coef","Site_Goren.coef","Site_Iftach.coef", "Site_KeremMaharl.coef","Site_KfarShamai.coef","Site_Kisalon.coef","Site_Margaliot.coef","Site_Nehusha.coef","Site_NirEtzion.coef","Site_Ofer.coef","Site_Yagur.coef","cosinus_monitoring.time.diff.coef","Distance_Time.coef", "Intercept.p","Distance.p", "Time_diff.p", "Site.p","Cosin_time.p","Distance_Time.p" ,"species_abundance")


#add hebrew name
A <- as.data.table(all.data.mammals.final.final.without.rare)
sci_heb <- A[,.(Species_Heb=unique(Species_Heb)), keyby=Species]
sci_heb[,Species:=make.names(Species)]
coefp_maqi <- merge(coefp_maqi, sci_heb, by.x = "SciName", by.y = "Species")

write.csv(coefp_maqi, "coefficients_Med.Maquis.csv")

#golden jackal:
Maquis_abundmat <- Maquis_abundmat[,Site:=factor(Site)]

unique(Maquis_abundmat$Site)

golden_jackal.maquis <- Maquis_abundmat[,c(1:19,27:28)]

unique(golden_jackal.maquis$Site)

glm.golden.jackal_maquis <- glm.nb(Canis.aureus ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                                     cosinus_Monitoring.Time.Diff, data = golden_jackal.maquis, maxit = 999)

summary(glm.golden.jackal_maquis)

#according to the model distance from settlements and the interaction between distance and time are significant. 

# interact_plot(model = glm.golden.jackal_maquis, data=golden_jackal.maquis, pred = rescaled_Time.Diff, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus",colors = "Dark2")+ 
# scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
# theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
#   ylab("Golden Jackal abundance")

source("functions/plot_interaction_mammals.R")
plot_model_interaction(P.anal = golden_jackal.maquis, m = glm.golden.jackal_maquis, eff2plot = "rescaled_Time.Diff", 
                       modvar2plot = "Distance_rescaled", plot_points=FALSE, 
                       plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                       fontsize=22, pdf_width=160, outpath = "output/maquis/")


#no dots:
effect_plot_g.Jackal.maquis <- effect_plot(model = glm.golden.jackal_maquis, data=golden_jackal.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#Testing whether trends and means are different from zero
#Test whether temporal trend is different from zero
#Calculate mean and sd of distance
sd_distance <- sd(env_data.maquis$Distance_rescaled)
print(sd_distance)

#1

mean_distance <- mean(env_data.maquis$Distance_rescaled)
print(mean_distance)

#-2.324816e-17

#calculate +/1 1 SD
mean_distance - sd_distance 

#-1

mean_distance + sd_distance

#1



mm_rescaled_time <- emtrends(glm.golden.jackal_maquis, specs = "Distance_rescaled", var = "rescaled_Time.Diff", type = "response", at = list(Distance_rescaled = c(-1,1)))
print(mm_rescaled_time)


#Results are averaged over the levels of: Site 
#Confidence level used: 0.95

test_results_mm_rescaled_time <- test(mm_rescaled_time, null = 0, adjust = "fdr")
print(test_results_mm_rescaled_time)

#Distance_rescaled rescaled_Time.Diff.trend     SE  df z.ratio p.value
#-1                    0.124 0.0728 Inf   1.700  0.0891
#1                   -0.169 0.0809 Inf  -2.089  0.0735

#Results are averaged over the levels of: Site 
#P value adjustment: fdr method for 2 tests 

#meaning: temporal trend is almost significantly different from zero in small (p=0.089) and large (p=0.073) distances from settlements. 


#Testing for pairwise differences between the distances levels:
mm_distance <- emmeans(object = glm.golden.jackal_maquis, ~Distance_rescaled*rescaled_Time.Diff)
test_results_distance <- test(pairs(mm_distance, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)


mm_distance <- emmeans(object = glm.golden.jackal_maquis, specs = "Distance_rescaled", by = "rescaled_Time.Diff", at = list(Distance_rescaled = c(-1,1)))
print(mm_distance)


#Results are averaged over the levels of: Site 
#Results are given on the log (not the response) scale. 
#Confidence level used: 0.95 

test_results_distance <- test(pairs(mm_distance, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)


#meaning:#in average, the difference of jackal abundance between near and far is significant (p<0.0001).


#m. gazelle:
gazelle.maquis <- Maquis_abundmat[,c(1:18,21,27:28)]

unique(gazelle.maquis$Site)

glm.gazelle_maquis <- glm.nb(Gazella.gazella ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                                     cosinus_Monitoring.Time.Diff, data = gazelle.maquis, maxit = 999)

summary(glm.gazelle_maquis)

#model without the interaction: according to the model distance and time are significant, separately. 
glm.gazelle_maquis2 <- glm.nb(Gazella.gazella ~ Distance_rescaled+rescaled_Time.Diff+ Site+
                               cosinus_Monitoring.Time.Diff, data = gazelle.maquis, maxit = 999)

summary(glm.gazelle_maquis2)


# Distance:


source("functions/plot_model_effect_mammals.R")
gazella.distance <- plot_model_effect(P.anal = gazelle.maquis, m = glm.gazelle_maquis2, 
                                      eff2plot = "Distance_rescaled", plot_points=FALSE, 
                                      plot_residuals=FALSE, export_plot=TRUE, 
                  ylabel = NULL, fontname = "Almoni ML v5 AAA",
                  fontsize=22, pdf_width=160, outpath = "output/maquis/")

gazella.distance <- effect_plot(model = glm.gazelle_maquis2, data=gazelle.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1") + 
   scale_x_continuous(breaks=c(-0.70378,1.134963,3.178015,5.221068), labels=c("100", "1000", "2000", "3000"), name = "Distance (m)")

gazella.distance <- data.frame(Site = "Aderet",
                              cosinus_Monitoring.Time.Diff = -0.08646667, 
                              rescaled_Time.Diff = -0.226265, 
                              Distance_rescaled = c(-1.455879, -1.132998, 
                                                    -0.7293972, -0.3257959,
                                                    0.07780538))

predict(glm.gazelle_maquis2, newdata = gazella.distance,
        type = "response")

#Increase from 100 to 2000 m of 224.08% or 3.24 fold

(1.6346021/0.5043732 ) - 1  

1.6346021/0.5043732 

#Increase of 136.26% every 500 m

1.1995761/0.8803260
0.8803260/0.6460400 


#time:
gazelle.time <- effect_plot_gazelle.maquis3 <- effect_plot(model = glm.gazelle_maquis2, data=gazelle.maquis, pred = rescaled_Time.Diff, interval = T, plot.points = T, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
 effect_plot_gazelle.maquis3+ scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
   theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
   ylab("Mountain Gazelle abundance")

 gazelle.time <- data.frame(Site = "Aderet",
                              sinus_Monitoring.Time.Diff = 0.1007873, 
                              cosinus_Monitoring.Time.Diff = -0.08646667, 
                              rescaled_Time.Diff = c(-1.691456827,0.950479961),
                              Distance_rescaled = -2.317273e-17)

predict(glm.gazelle_maquis2, newdata = gazelle.time,
        type = "response")

#Increase of 60.95% during monitoring time, or 1.6 fold

(1.903606 /1.182676 ) - 1  

1.903606/1.182676 





source("functions/plot_model_effect_mammals.R")
gazella.time <- plot_model_effect(P.anal = gazelle.maquis, m = glm.gazelle_maquis2, 
                                      eff2plot = "rescaled_Time.Diff", plot_points=FALSE, 
                                      plot_residuals=FALSE, export_plot=TRUE, 
                                      ylabel = NULL, fontname = "Almoni ML v5 AAA",
                                      fontsize=22, pdf_width=160, outpath = "output/maquis/")


#hyena:
hyena.maquis <- Maquis_abundmat[,c(1:18,22,27:28)]

#unique(gazelle.maquis$Site)

glm.hyena_maquis <- glm.nb(Hyaena.hyaena ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                               cosinus_Monitoring.Time.Diff, data = hyena.maquis, maxit = 999)

summary(glm.hyena_maquis)

#according to the model, distance and the interaction between distance and time are significant. 

hyaena.interact <- interact_plot(model = glm.hyena_maquis, data=hyena.maquis, pred = rescaled_Time.Diff, 
              modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus",colors="Dark2")+
  



#distance near
hyaena.interact <- data.frame(Site = "Aderet",
                           sinus_Monitoring.Time.Diff = 0.1007873, 
                           cosinus_Monitoring.Time.Diff = -0.08646667, 
                           rescaled_Time.Diff = c(-1.691456827,0.950479961),
                           Distance_rescaled = -2.317273e-17)

predict(glm.hyena_maquis, newdata = hyaena.interact,
        type = "response")

#Increase of 60.95% during monitoring time, or 1.6 fold

(1.903606 /1.182676 ) - 1  

1.903606/1.182676 



source("functions/plot_interaction_mammals.R")
plot_model_interaction(P.anal = hyena.maquis, m = glm.hyena_maquis, eff2plot = "rescaled_Time.Diff", 
                       modvar2plot = "Distance_rescaled", plot_points=FALSE, 
                       plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                       fontsize=22, pdf_width=160, outpath = "output/maquis/")

#to understand what happens today with hyena abundance:
interact_plot(model = glm.hyena_maquis, data=hyena.maquis, pred =Distance_rescaled , modx =rescaled_Time.Diff , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus")

#distance:
effect_plot_hyena.maquis <- effect_plot(model = glm.hyena_maquis, data=hyena.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = T, colors = "Qual1",line.colors = "black") 
effect_plot_hyena.maquis+ylim(0,160)

#no dots:
effect_plot_hyena.maquis <- effect_plot(model = glm.hyena_maquis, data=hyena.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#Test whether temporal trend is different from zero
install.packages("emmeans")
require(emmeans)

mm_rescaled_time5 <- emtrends(glm.hyena_maquis, specs = "Distance_rescaled", var = "rescaled_Time.Diff", type = "response", at = list(Distance_rescaled = c(-1,1)))
print(mm_rescaled_time5)

 
test_results_mm_rescaled_time <- test(mm_rescaled_time5, null = 0, adjust = "fdr")
print(test_results_mm_rescaled_time)


#meaning: temporal trend is not significantly different from zero in small distances from settlements (p=0.43) 
#and significantly different in large distances from settlements (p=0.01). 

#Testing whether the abundance is different between near and far: 
mm_distance5 <- emmeans(object = glm.hyena_maquis, specs = "Distance_rescaled", by = "rescaled_Time.Diff", at = list(Distance_rescaled = c(-1,1)))
print(mm_distance5)

 
test_results_distance <- test(pairs(mm_distance5, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)

small_distance <- mean_distance - sd_distance 
large_distance <- mean_distance + sd_distance

pairwise_distance <- emmeans(object = glm.hyena_maquis, ~Distance_rescaled*rescaled_Time.Diff, at = list(Distance_rescaled = c(small_distance,large_distance)))
test_results_distance <- test(pairs(pairwise_distance, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)



#porcupine:
porcupine.maquis <- Maquis_abundmat[,c(1:18,23,27:28)]

#unique(gazelle.maquis$Site)

glm.porcupine_maquis <- glm.nb(Hystrix.indica ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                             cosinus_Monitoring.Time.Diff, data = porcupine.maquis, maxit = 999)

summary(glm.porcupine_maquis)

interact_plot(model = glm.porcupine_maquis, data=porcupine.maquis, pred = rescaled_Time.Diff, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus", colors="Dark2")+
scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylab("Indian porcupine abundance")

source("functions/plot_interaction_mammals.R")
plot_model_interaction(P.anal = porcupine.maquis, m = glm.porcupine_maquis, eff2plot = "rescaled_Time.Diff", 
                       modvar2plot = "Distance_rescaled", plot_points=FALSE, 
                       plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                       fontsize=22, pdf_width=160, outpath = "output/maquis/")
#distance:
effect_plot_porcupine.maquis <- effect_plot(model = glm.porcupine_maquis, data=porcupine.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = T, colors = "Qual1",line.colors = "black") 
effect_plot_porcupine.maquis+ylim(0,210)

#no dots:
effect_plot_porcupine.maquis <- effect_plot(model = glm.porcupine_maquis, data=porcupine.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#time:
effect_plot_porcupine.maquis2 <- effect_plot(model = glm.porcupine_maquis, data=porcupine.maquis, pred = rescaled_Time.Diff, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = T, colors = "Qual1",line.colors = "black") 
effect_plot_porcupine.maquis2+ylim(0,200)

#no dots:
effect_plot_porcupine.maquis2 <- effect_plot(model = glm.porcupine_maquis, data=porcupine.maquis, pred = rescaled_Time.Diff, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#Test whether temporal trend is different from zero
install.packages("emmeans")
require(emmeans)


mm_rescaled_time2 <- emtrends(glm.porcupine_maquis, specs = "Distance_rescaled", var = "rescaled_Time.Diff", type = "response", at = list(Distance_rescaled = c(-1,1)))
print(mm_rescaled_time)

test_results_mm_rescaled_time <- test(mm_rescaled_time2, null = 0, adjust = "fdr")
print(test_results_mm_rescaled_time)

#meaning: temporal trend is significantly different from zero in small distances from settlements (p=0.0001) 
#and not significantly different in large distances from settlements (p=0.5). 

#Testing whether the abundance is differenct between near and far; 
mm_distance2 <- emmeans(object = glm.porcupine_maquis, specs = "Distance_rescaled", by = "rescaled_Time.Diff", at = list(Distance_rescaled = c(-1,1)))
print(mm_distance2)

test_results_distance <- test(pairs(mm_distance2, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)


#Results are averaged over the levels of: Site 
#Results are given on the log (not the response) scale. 
 
pairwise_distance <- emmeans(object = glm.porcupine_maquis, ~Distance_rescaled*rescaled_Time.Diff, at = list(Distance_rescaled = c(small_distance,large_distance)))
test_results_distance <- test(pairs(pairwise_distance, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)

#badger:
badger.maquis <- Maquis_abundmat[,c(1:18,24,27:28)]


#model without interaction:
glm.badger_maquis <- glm.nb(Meles.meles ~ Distance_rescaled+rescaled_Time.Diff+ Site+
                                 cosinus_Monitoring.Time.Diff, data = badger.maquis, maxit = 999)

summary(glm.badger_maquis)

#distance:
effect_plot_badger.maquis <- effect_plot(model = glm.badger_maquis, data=badger.maquis, pred = Distance_rescaled, interval = T, plot.points = T, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
effect_plot_badger.maquis+ylim(0,20)

#no dots:
effect_plot_badger.maquis <- effect_plot(model = glm.badger_maquis, data=badger.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#with actual values of distance:
# effect_plot(model = glm.badger_maquis, data=badger.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1") + 
#   scale_x_continuous(breaks=c(-0.70378,1.134963,3.178015,5.221068), labels=c("100", "1000", "2000", "3000"), name = "Distance (m)")

source("functions/plot_model_effect_mammals.R")
plot_model_effect(P.anal = badger.maquis, m = glm.badger_maquis, eff2plot = "Distance_rescaled", 
                  plot_points=FALSE, plot_residuals=FALSE, 
                  export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                  fontsize=22, pdf_width=160, outpath = "output/maquis/")
#wild boar:
w.boar.maquis <- Maquis_abundmat[,c(1:18,25,27:28)]


glm.w.boar_maquis <- glm.nb(Sus.scrofa ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                              cosinus_Monitoring.Time.Diff, data = w.boar.maquis, maxit = 999)

summary(glm.w.boar_maquis)

interact_plot(model = glm.w.boar_maquis, data=w.boar.maquis, pred = rescaled_Time.Diff, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus")


#Error in qr.solve(qr.R(qr.lm(object))[p1, p1]) : 
#singular matrix 'a' in solve
#In addition: Warning message:
#  -1 is outside the observed range of Distance_rescaled 

#maybe remove sites that there are no boar`s obs:
dotchart(w.boar.maquis$Sus.scrofa, ylab = "Site",
         groups = w.boar.maquis$Site, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")

#Aderet and Givat Yashahu lack boar obs. 
dotchart(w.boar.maquis$Sus.scrofa, ylab = "Project.ID",
         groups = w.boar.maquis$Project.ID, 
         # color = my_cols[sp.richness.Maquis.per.cam$Settlements],
         cex = 0.9,  pch = 1, xlab = "abundance")

#remove Aderet and Givat yeshayahu:
w.boar.maqi_thin <- subset(w.boar.maquis, Site!="Givat Yeshayahu" & Site!="Aderet")

unique(w.boar.maqi_thin$Site)      

w.boar.maqi_thin <- w.boar.maqi_thin[,Site:=factor(Site)]
                       
glm.w.boar_maquis2 <- glm.nb(Sus.scrofa ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                              cosinus_Monitoring.Time.Diff, data = w.boar.maqi_thin, maxit = 999)

summary(glm.w.boar_maquis2)

# interact_plot(model = glm.w.boar_maquis2, data=w.boar.maqi_thin, pred = rescaled_Time.Diff, modx =Distance_rescaled , partial.residuals = F, jitter = c(0.2), point.size = 3, rug = F, point.shape = F, interval = T,modx.values = "plus-minus", colors="Dark2")+
# scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
#   theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
#   ylab("Wild boar abundance")

source("functions/plot_interaction_mammals.R")
plot_model_interaction(P.anal = w.boar.maqi_thin, m = glm.w.boar_maquis2, eff2plot = "rescaled_Time.Diff", 
                       modvar2plot = "Distance_rescaled", plot_points=FALSE, 
                       plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                       fontsize=22, pdf_width=160, outpath = "output/maquis/")

#no dots:
effect_plot_boar.maquis <- effect_plot(model = glm.w.boar_maquis2, data=w.boar.maqi_thin, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#Test whether temporal trend is different from zero

mm_rescaled_time3 <- emtrends(glm.w.boar_maquis2, specs = "Distance_rescaled", var = "rescaled_Time.Diff", type = "response", at = list(Distance_rescaled = c(-1,1)))
print(mm_rescaled_time)

 
test_results_mm_rescaled_time <- test(mm_rescaled_time3, null = 0, adjust = "fdr")
print(test_results_mm_rescaled_time)
 

#meaning: temporal trend is significantly different from zero in small distances from settlements (p=0.003) 
#and not significantly different in large distances from settlements (p=0.91). 

#Testing whether the abundance of boars is different between near and far:
mm_distance3 <- emmeans(object = glm.w.boar_maquis2, specs = "Distance_rescaled", by = "rescaled_Time.Diff", at = list(Distance_rescaled = c(-1,1)))
print(mm_distance3)

test_results_distance <- test(pairs(mm_distance3, by="rescaled_Time.Diff"), by=NULL, adjust="fdr")
print(test_results_distance)


#red fox: model with no interaction
fox.maquis <- Maquis_abundmat[,c(1:18,26:28)]

#unique(gazelle.maquis$Site)

#model without interaction:
glm.fox_maquis <- glm.nb(Vulpes.vulpes~ Distance_rescaled+rescaled_Time.Diff+ Site+
                              cosinus_Monitoring.Time.Diff, data = fox.maquis, maxit = 999)

summary(glm.fox_maquis)

source("functions/plot_model_effect_mammals.R")
plot_model_effect(P.anal = fox.maquis, m = glm.fox_maquis, eff2plot = "Distance_rescaled", 
                  plot_points=FALSE, plot_residuals=FALSE, 
                  export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                  fontsize=22, pdf_width=160, outpath = "output/maquis/")

#with actual values of distance:
# effect_plot(model = glm.fox_maquis, data=fox.maquis, pred = Distance_rescaled, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1") + 
#   scale_x_continuous(breaks=c(-0.70378,1.134963,3.178015,5.221068), labels=c("100", "1000", "2000", "3000"), name = "Distance (m)")




 ############
#special story- boars in Judea:
w.boar.jud <- jud_abundmat[,c(1:18,25,27:29)]

p5 <- ggplot(w.boar.jud, aes(fill=Project.ID,y=Sus.scrofa,x=Project.ID))+ geom_bar(position="dodge", stat="identity")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# scale_y_continuous(limits = c(0,100))

p5

p6 <- ggplot(w.boar.jud, aes(x=Project.ID, y=Sus.scrofa))+ geom_point(aes(),size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# scale_y_continuous(limits = c(0,100))

p6

#the model fits the separated model if judea:
mva.jud.nb1 <- manyglm(subunit_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                         cosinus_Monitoring.Time.Diff,
                       data = envd.jud, family="negative.binomial")



glm.w.boar_jud <- glm.nb(Sus.scrofa ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                           cosinus_Monitoring.Time.Diff, data = w.boar.jud, maxit = 999)

summary(glm.w.boar_jud)

effect_plot_boar.jud <- effect_plot(model = glm.w.boar_jud, data=w.boar.jud, pred = rescaled_Time.Diff, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = T, colors = "Qual1",line.colors = "black") 
effect_plot_boar.jud+ylim(0,650)

#Error in qr.solve(qr.R(qr.lm(object))[p1, p1]) : 
#singular matrix 'a' in solve

#remove Aderet and Givat Yashayahu, and t0-t1:

unique(w.boar.jud_thin$Site)
unique(w.boar.jud_thin$Project.ID)

glm.w.boar_jud2 <- glm.nb(Sus.scrofa ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                           cosinus_Monitoring.Time.Diff, data = w.boar.jud_thin, maxit = 999)

summary(glm.w.boar_jud2)

effect_plot_boar.jud <- effect_plot(model = glm.w.boar_jud2, data=w.boar.jud_thin, pred = rescaled_Time.Diff, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = T, colors = "Qual1",line.colors = "black") 
effect_plot_boar.jud+ylim(0,15)

effect_plot_boar.jud <- effect_plot(model = glm.w.boar_jud2, data=w.boar.jud_thin, pred = rescaled_Time.Diff, interval = T, plot.points = F, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 

#
w.boar.jud_thin2 <- subset(w.boar.jud, Site!="Givat Yeshayahu" & Site!="Aderet")

unique(w.boar.jud_thin2$Site)      

w.boar.jud_thin2 <- w.boar.jud_thin2[,Site:=factor(Site)]

unique(w.boar.jud_thin2$Project.ID)

glm.w.boar_jud3 <- glm.nb(Sus.scrofa ~ Distance_rescaled*rescaled_Time.Diff+ Site+
                            cosinus_Monitoring.Time.Diff, data = w.boar.jud_thin2, maxit = 999)

summary(glm.w.boar_jud3)

#Wild boar judea - final plot

source("functions/plot_model_effect_mammals.R")
plot_model_effect(P.anal = w.boar.jud_thin2, m = glm.w.boar_jud3, eff2plot = "rescaled_Time.Diff", 
                  plot_points=FALSE, plot_residuals=FALSE, 
                  export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                  fontsize=22, pdf_width=160, outpath = "output/maquis/")

# effect_plot_boar.jud <- effect_plot(model = glm.w.boar_jud3, data=w.boar.jud_thin2, pred = rescaled_Time.Diff, interval = T, plot.points = T, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
# effect_plot_boar.jud+scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
#   theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
#   ylab("Wild boar abundance")


#golden jackal in Galil:
library(data.table)
galil_abundmat <- as.data.table(galil_abundmat)

galil_abundmat <- Maquis.mammals.abundance.matrix[Maquis.mammals.abundance.matrix$Subunit=="Galilee",]

unique(galil_abundmat$Site)

galil_abundmat <- galil_abundmat[,Site:=factor(Site)]

unique(galil_abundmat$Site)

golden_jackal.galil <- galil_abundmat[,c(1:19,27:28)]

unique(golden_jackal.galil$Site)

glm.golden.jackal_galil <- glm.nb(Canis.aureus ~ Distance_rescaled+rescaled_Time.Diff+ Site+
                                    sinus_Monitoring.Time.Diff, data = golden_jackal.galil, maxit = 999)

summary(glm.golden.jackal_galil)

#according to mva_galil2 distance from settlements and time have significant effect, separately. So, effect graph for time:

#time:
effect_plot_g.Jackal2 <- effect_plot(model = glm.golden.jackal_galil, data=golden_jackal.galil, pred = rescaled_Time.Diff, interval = T, plot.points = T, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = F, colors = "Qual1",line.colors = "black") 
effect_plot_g.Jackal2+scale_x_continuous(breaks=c(-1.86423,-0.91809,-0.28646,0.344304,0.975931), labels=c("June 2012", "June 2015", "June 2017", "June 2019", "June 2021"), name = "Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylab("Golden jackal abundance")


source("functions/plot_model_effect_mammals.R")
plot_model_effect(P.anal = golden_jackal.galil, m = glm.golden.jackal_galil, eff2plot = "rescaled_Time.Diff", 
                  plot_points=FALSE, plot_residuals=FALSE, 
                  export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                  fontsize=22, pdf_width=160, outpath = "output/maquis/")

#analysis separately for each subunit:
#First, Judea:
#envd.jud <- env_data.maquis[Subunit=="Judea"][,Site:=factor(Site)]

unique(envd.jud$Site)

#subunit_ind <- env_data.maquis$Subunit=="Judea"

#subunit_spp <- Maquis_spp[subunit_ind,colSums(Maquis_spp[subunit_ind,])>0]

#data exploration: was done in Maquis_analysis_300523

mva.jud.poiss <- manyglm(subunit_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                           sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                         data = envd.jud, family="poisson")


mva.jud.nb <- manyglm(subunit_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                           sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                         data = envd.jud, family="negative.binomial")

plot.manyglm(mva.jud.poiss, which = 1:3)
plot.manyglm(mva.jud.nb, which = 1:3)

#mva.jud.nb looks much better.

drop1(mva.jud.nb)

#
#Model:
 # subunit_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + 
  #sinus_Monitoring.Time.Diff + cosinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  5378.7
#Site                                 28 5552.1
#sinus_Monitoring.Time.Diff            7 5371.1
#cosinus_Monitoring.Time.Diff          7 5389.2
#Distance_rescaled:rescaled_Time.Diff  7 5401.9

#drop sinus:
mva.jud.nb1 <- manyglm(subunit_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                        cosinus_Monitoring.Time.Diff,
                      data = envd.jud, family="negative.binomial")

drop1(mva.jud.nb1)

#
#Model:
 # subunit_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + 
  #cosinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  5371.1
#Site                                 28 5545.7
#cosinus_Monitoring.Time.Diff          7 5382.9
#Distance_rescaled:rescaled_Time.Diff  7 5395.1

#mva.jud.nb1 is the final model.

anov_mva.jud.nb1 <- anova(mva.jud.nb1, p.uni = "adjusted")
print(anov_mva.jud.nb1)

save(all.data.mammals,all.data.mammals_below.10th.day,all.data.mammals.final,all.data.mammals.final_duplicates,all.data.mammals.final_Inland.sands_new,all.data.mammals.final.final,corrected.sites,corrected.sites_below.10th.day,counting.duplicates,date.per.Unit.Site_corrected,deployments_Negev_Traznisition_T4,Duplicate_obs_new,Duplicate_obs,duplicate_obs_all_by.deployment,list_of_duplicates,list_of_duplicates_prop,list_of_duplicates_prop_sum,list_of_duplicates_prop_freq1, list_of_duplicates_proportion,list_of_duplicates_prop_sum_all, list_of_duplicates_sum,list_of_duplicates_proportion_by.count,list_of_transects,species_Latin_Heb_names,T4.deployments.data,T4.mammals_below.10th.day,T4.mammals_below.10th.day_Negev_Transition,T4.mammals.data,T4.mammals.data_Negev_Transition,T4.mammals.data_seq.correction,T4.mammals.data_seq.correction_Negev_Transition,T4.mammals.data.combined, T4.mammals.data.combined_Negev_Transition,T4.mammals.data.new,T4.mammals.data.new_Negev_Transition,Unit_per_year_corrected,Unit.site_per_project,Unit.site_per_year_corrected,Unit.site_per_year_Negev_Transition_T4,Time,Time.diff,Time_diff_below_30.min,Time_diff_proportion,Time_diff_proportion_new,all.data.mammals_cast,mammals.abundance_matrix,all.data.mammals_cast_2,mammals.indcidence.matrix,mammals_indcidence_matrix,Rare,rare.species,Maquis.data,Maquis_total_counts_site_campaign,Maquis_total_counts,Maquis_total_counts_campaign_all_sp,Maquis_total_counts_campaign_proximity_all_sp,mammals.abundance.matrix,sp.richness.site,all.data.mammals.final.final.without.rare,Rare.campaign,rare.species.per.campaign,Abund_total_counts_campaign,sp.richness.Maquis,Maquis.data,Rare_by_count,Rare_by_count_deployment,Rare.sp_by.two.conditions,all.data.mammals.final.final.for.richness,sp.richness.Maquis.per.cam,Maquis.richness.matrix,Maquis.incidence.long,list_of_deployments,list_of_existing_deployments,Distance_data, distance_data_T4, distance_data_T0_T3,Agri_Distance_data,Maquis.abu.data,env_data.maquis,Maquis_abu,Maquis_abu2,site.scrs,spp.scrs,nMDS,nMDS2, deployments_Yatir_Mirsham, T4.mammals.data.new_yatir_mirsham,T4.mammals.data_yatir_mirsham, T4.mammals.data.combined_yatir_mirsham, T4.mammals_below.10th.day_yatir_mirsham,Maquis_spp,Maquis.mammals.abundance.matrix,nMDS4,Maquis_abu_without.T2.GorenFar7,Maquis_abu_without.T2.GorenNear4,env_data.maquis_without.T2.GorenFar7,env_data.maquis_without.T2.GorenNear4, Herb.data,Herb.mammals.abundance.matrix,Herb_abu,Herb_abundmat,Herb_spp,sp.richness.Herb.per.cam,Herb.abu.data,env_data.herb,Distance_diff,mvabund4,mvabund43,anov_mvabund43,mvabund53,mvabund5, mvabund6, mvabund63,anov_mvabund63, coefp,Hyena_Maquis,Hyena_Maquis_per.unit,Maquis.sample.size.per.species,Maquis.sample.size.per.species.no.campaign, model.coef, Maquis.cam.number,Maquis.cam.number.per.campaign,Maquis.cam.number.per.campaign.per.species,Maquis.cam.number.per.species,Jackal_Maquis,Jackal_Maquis_per.unit, coefp_Herb,coefp_Herb2,envd.galil,envd.jud,galil_spp,subunit_spp,galilee.abundance.matrix,galil_abu,jud_abu,judea.abundance.matrix,mva.galil2,mva.jud2,mva0.jud,mva1.jud,mva1.jud2,mva0, anov_mva.jud2, mva.carmel2,jud_abundmat,coefp_jud,model.coef_jud,carmel.abundance.matrix,carmel_abundmat,carmel_spp,envd.carmel,carmel_abu,anov_mva.carmel2, coefp_carmel, mva.carmel2, Maquis_abundmat,maqi_mvabund_nb2,maqi_mvabund_nb21,anov_maqi_mvabund21,anov_mva.jud.nb1,mva.jud.nb1,coefp_maqi,coefp_maqi_judea, file="Mammals_analysis_2023.RData")

#coef
coefp_maqi_judea <- merge(data.table(t(coef(mva.jud.nb1)),keep.rownames=TRUE), data.table(t(anov_mva.jud.nb1$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_maqi_judea <- merge(coefp_maqi_judea,as.data.table(colSums(subunit_spp),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_maqi_judea)
colnames(coefp_maqi_judea) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_GivatYearim.coef", "Site_GivatYeshayahu.coef","Site_Kisalon.coef","Site_Nehusha.coef","cosinus_monitoring.time.diff.coef","Distance_Time.coef", "Intercept.p","Distance.p", "Time_diff.p", "Site.p","Cosin_time.p","Distance_Time.p" ,"species_abundance")


#add hebrew name
#A <- as.data.table(all.data.mammals.final.final.without.rare)
#sci_heb <- A[,.(Species_Heb=unique(Species_Heb)), keyby=Species]
#sci_heb[,Species:=make.names(Species)]
coefp_maqi_judea <- merge(coefp_maqi_judea, sci_heb, by.x = "SciName", by.y = "Species")

write.csv(coefp_maqi_judea, "coefficients_Med.Maquis_judea.csv")


##Galil:
envd.galil <- env_data.maquis[Subunit=="Galilee"][,Site:=factor(Site)]

unique(envd.galil$Site)

galil_ind <- env_data.maquis$Subunit=="Galilee"

galil_spp <- Maquis_spp[galil_ind,colSums(Maquis_spp[galil_ind,])>0]

#data exploration:in Maqis_analysis_300523

mva.galil.poiss <- manyglm(galil_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                             sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                           data = envd.galil, family="poisson")


mva.galil <- manyglm(galil_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                       sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                     data = envd.galil)

plot.manyglm(mva.galil.poiss, which = 1:3)
plot.manyglm(mva.galil, which = 1:3)

#mva.galil looks better.

#model selection:
drop1(mva.galil)

#
#Model:
 # galil_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + sinus_Monitoring.Time.Diff + 
  #cosinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  9046.0
#Site                                 32 9320.1
#sinus_Monitoring.Time.Diff            8 9049.6
#cosinus_Monitoring.Time.Diff          8 9037.4
#Distance_rescaled:rescaled_Time.Diff  8 9037.9

#drop cosinus:
mva.galil1 <- manyglm(galil_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                       sinus_Monitoring.Time.Diff,
                     data = envd.galil)
drop1(mva.galil1)

#
#Model:
 # galil_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + sinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  9037.4
#Site                                 32 9345.7
#sinus_Monitoring.Time.Diff            8 9040.7
#Distance_rescaled:rescaled_Time.Diff  8 9029.1

#drop distance:time:
mva.galil2 <- manyglm(galil_spp~Distance_rescaled+rescaled_Time.Diff+ Site+
                        sinus_Monitoring.Time.Diff,
                      data = envd.galil)
drop1(mva.galil2)

#
#Model:
 # galil_spp ~ Distance_rescaled + rescaled_Time.Diff + Site + sinus_Monitoring.Time.Diff
#Df    AIC
#<none>                        9029.1
#Distance_rescaled           8 9132.8
#rescaled_Time.Diff          8 9041.5
#Site                       32 9336.1
#sinus_Monitoring.Time.Diff  8 9035.0

#mva.galil2 is the final model.

anov_mva.galil2 <- anova(mva.galil2, p.uni = "adjusted")
print(anov_mva.galil2)

#coef
coefp_maqi_galil <- merge(data.table(t(coef(mva.galil2)),keep.rownames=TRUE), data.table(t(anov_mva.galil2$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_maqi_galil <- merge(coefp_maqi_galil,as.data.table(colSums(galil_spp),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_maqi_galil)
colnames(coefp_maqi_galil) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_Goren.coef", "Site_Iftach.coef","Site_Kfar Shamai.coef","Site_Margaliot.coef","sinus_monitoring.time.diff.coef", "Intercept.p","Distance.p", "Time_diff.p", "Site.p","sin_time.p" ,"species_abundance")


#add hebrew name
#A <- as.data.table(all.data.mammals.final.final.without.rare)
#sci_heb <- A[,.(Species_Heb=unique(Species_Heb)), keyby=Species]
#sci_heb[,Species:=make.names(Species)]
coefp_maqi_galil <- merge(coefp_maqi_galil, sci_heb, by.x = "SciName", by.y = "Species")

write.csv(coefp_maqi_galil, "coefficients_Med.Maquis_galil.csv")



#carmel:
envd.carmel <- env_data.maquis[Subunit=="Carmel"][,Site:=factor(Site)]

unique(envd.carmel$Site)

carmel_ind <- env_data.maquis$Subunit=="Carmel"

carmel_spp <- Maquis_spp[carmel_ind,colSums(Maquis_spp[carmel_ind,])>0]

#data exploration:in Maquis_analysis_300523

mva.carmel.poiss <- manyglm(carmel_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                              sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                            data = envd.carmel, family="poisson")

plot.manyglm(mva.carmel.poiss, which = 1:3)


mva.carmel <- manyglm(carmel_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                        sinus_Monitoring.Time.Diff+ cosinus_Monitoring.Time.Diff,
                      data = envd.carmel)

plot.manyglm(mva.carmel, which = 1:3)

unique(envd.carmel$Site)

drop1(mva.carmel)

#
#Model:
 # carmel_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + 
#  sinus_Monitoring.Time.Diff + cosinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  7650.7
#Site                                 28 7727.1
#sinus_Monitoring.Time.Diff            7 7641.6
#cosinus_Monitoring.Time.Diff          7 7647.9
#Distance_rescaled:rescaled_Time.Diff  7 7663.3

#drop sinus:
mva.carmel1 <- manyglm(carmel_spp~Distance_rescaled*rescaled_Time.Diff+ Site+
                         cosinus_Monitoring.Time.Diff,
                      data = envd.carmel)

drop1(mva.carmel1)

#
#Model:
 # carmel_spp ~ Distance_rescaled * rescaled_Time.Diff + Site + 
  #cosinus_Monitoring.Time.Diff
#Df    AIC
#<none>                                  7641.6
#Site                                 28 7717.4
#cosinus_Monitoring.Time.Diff          7 7639.1
#Distance_rescaled:rescaled_Time.Diff  7 7656.6

#drop cosinus:
mva.carmel2 <- manyglm(carmel_spp~Distance_rescaled*rescaled_Time.Diff+ Site,
                       data = envd.carmel)

drop1(mva.carmel2)

#
#Model:
 # carmel_spp ~ Distance_rescaled * rescaled_Time.Diff + Site
#Df    AIC
#<none>                                  7639.1
#Site                                 28 7727.9
#Distance_rescaled:rescaled_Time.Diff  7 7652.4

#mva.carmel2 is the final model.

anov_mva.carmel2 <- anova(mva.carmel2, p.uni = "adjusted")
print(anov_mva.carmel2)

#coef
coefp_maqi_carmel <- merge(data.table(t(coef(mva.carmel2)),keep.rownames=TRUE), data.table(t(anov_mva.carmel2$uni.p), keep.rownames = TRUE), by = "rn") # coefficients are exponentiated because the link function used for poisson is log -> above 1 is positive and below 1 is negative

# add total species abundance
coefp_maqi_carmel <- merge(coefp_maqi_carmel,as.data.table(colSums(carmel_spp),keep.rownames = TRUE), by.x = "rn", by.y = "V1")

colnames(coefp_maqi_carmel)
colnames(coefp_maqi_carmel) <- c("SciName","Intercept.coef","Distance.coef","Time_diff.coef","Site_Kerem Maharal.coef", "Site_Nir Eztion.coef","Site_Ofer.coef","Site_Yagur.coef","Distance_time.coef", "Intercept.p","Distance.p", "Time_diff.p", "Site.p","Distance_time.p" ,"species_abundance")


#add hebrew name
#A <- as.data.table(all.data.mammals.final.final.without.rare)
#sci_heb <- A[,.(Species_Heb=unique(Species_Heb)), keyby=Species]
#sci_heb[,Species:=make.names(Species)]
coefp_maqi_carmel <- merge(coefp_maqi_carmel, sci_heb, by.x = "SciName", by.y = "Species")

write.csv(coefp_maqi_carmel, "coefficients_Med.Maquis_carmel.csv")

sum(Maquis_abundmat$`Canis aureus`)
sum(Maquis_abundmat$`Hystrix indica`)
sum(Maquis_abundmat$`Gazella gazella`)

nrow(Maquis_abundmat[Maquis_abundmat$Subunit == "Judea",])

sum(Maquis_abundmat$Subunit == "Judea")

sum(jud_abundmat$Sus.scrofa)
sum(galil_abundmat$Canis.aureus)

write.csv(env_data.maquis, "env_data.maquis.csv")

sd(env_data.maquis$Distance)
mean(env_data.maquis$Distance)

#Kerem Maharal: try to differenciate between distance from settlements and distance from agriculture:
Kerem.Maharal <- Maquis_abundmat[Maquis_abundmat$Site=="Kerem Maharal",]

write.csv(Kerem.Maharal, "Kerem.maharal.csv")

#put both distances in the same column

p6 <- ggplot(fox.maquis, aes(x=Distance_rescaled, y=Vulpes.vulpes))+ geom_point(aes(color=Subunit),size=5)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# scale_y_continuous(limits = c(0,100))

p6