
save(all.data.mammals.final.final.without.rare,Herb.abu.data,Herb.mammals.abundance.matrix,Herb_abundmat,Herb_spp,Herb_abu,env_data.herb,Rare.sp_by.two.conditions,Herb_abu_no.wolf,Herb_spp2,mva_herb_nb_no.wolf, anov_mva_herb_nb_no.wolf,coefp_Herb_no.wolf,Herb,Herb.data,fox.herb,glm.fox.herb,file="Herb_analysis_2023.RData")










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
