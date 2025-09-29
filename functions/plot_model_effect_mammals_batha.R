plot_model_effect <- function (P.anal, m, eff2plot, plot_points=TRUE, plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                               fontsize=22, pdf_width=160, outpath = "../output/") {
  # this function is for plotting and exporting to pdf of categorical effects or continuous year effect, no interactions. using the function effect_plot from jtools.
  # P.anal is the data table
  # m is the model object, outcome of glm or glmer
  # eff2plot is the factor name for plotting (independent variable)
  # plot_points, plot_residuals: booleans. if plot_residuals is TRUE overrides then residuals will be plotted, regardless of the value of plot_points.
  # export_plot is a boolean indicating whether to export the file to PDF. The file name is created automatically
  # ylabel is an optional character vector specifying the label of the y-axis. If not given then the label will be the response variable in Hebrew
  #   (for richness, gma and abundance) or English (for other variables)
  # fontname is the name of the font for all text in the plot, one of the fonts included in extrafont::fonts()
  # fontsize is the size of all text in the plot
  # pdf_width is the export plot size in mm
  # outpath is the path for saving exported pdf files
  
  # Define a function to reverse Hebrew text in case markdown is knit to pdf
  strReverse <- function(x){sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")}
  
  pdf_aspect_ratio <- 2/3
  
  require(jtools)
  require(Cairo)
  require(ggplot2)
  require(extrafont)
  
  P.anal <- copy(P.anal)
  
  # determine what is the response and set the y label
  response_var <- as.character(formula(m))[2]
  if (is.null(ylabel)) {
    ylabel <- switch(response_var,
                     "Canis.aureus" = "שפע",
                     "Canis.lupus" = "שפע",
                     "Vulpes.vulpes" = "שפע",
                     "Lepus.capensis" = "שפע",
                     "Hystrix.indica" = "שפע",
                     "Gazella.gazella" = "שפע",
                     "Gazella.dorcas" = "שפע",
                     "Meles.meles" = "שפע",
                     "Hyaena.hyaena" = "שפע",
                     "Sus.scrofa" = "שפע",
                     "`Sus scrofa`" = "שפע",
                     "`Canis aureus`" = "שפע",
                     "gma" = "שפע ממוצע למין",
                     "abundance" = "שפע כולל",
                     
                     response_var)
  }

  # # if near vs. far factor, then relevel to get Near on the left hand side
  # if (eff2plot %in% c("agriculture","settlements")) {
  #   P.anal[,c(eff2plot):=factor(get(eff2plot),levels = c("Far","Near"))]
  # 
  #   # if year covariate (continuous) then get the data labels to reverse them if needed
  # } else if (eff2plot=="year_ct") {
  #   tmp_effplot <- effect_plot(model = m, pred = year_ct, data = P.anal)
  #   xticklabs <- ggplot_build(tmp_effplot)$layout$panel_params[[1]]$x$get_labels()
  # }
  #  
  
  # determine what is the plotted effect and set the x label
  xlabel <- switch(eff2plot,
                   "settlements" = "ישוב",
                   "agriculture" = "חקלאות",
                   "dunes" = "סיווג הדיונה",
                   "Subunit" = "תת-יחידה",
                   "rescaled_Time.Diff_new" = "זמן",
                   "Distance_rescaled" = "מרחק מיישוב )מטר(",
                   "Distance_agri_rescaled" = "מרחק מחקלאות )מטר(")
  
  # if (eff2plot %in% c("settlements","agriculture","dunes","habitat","land_use","Subunit")) {
  # 
  #   #  Set the x tick mark labels (labels of levels of the factor plotted)
  #   # create a table containing English level names and Hebrew translations
  #   fac_lev <- data.table(fac = c("settlements","settlements","agriculture","agriculture","dunes","dunes","habitat","habitat",
  #                                 "land_use","land_use","land_use","Subunit","Subunit","Subunit"),
  #                         eng_lev = c("Near","Far","Near","Far","Shifting","Semi-shifting","Slope","Wadi",
  #                                     "Beduoin Agriculture","KKL Plantings","Loess","Carmel","Galilee","Judea"),
  #                         heb_lev = c("קרוב","רחוק","קרוב","רחוק","נודדת","מיוצבת למחצה","מדרון","ערוץ",
  #                                     "חקלאות מסורתית","נטיעות קקל","לס","כרמל","גליל","הרי-יהודה"))
  #   #get the levels of the current factor, in the order in which they appear
  #   curr_lev <- data.table(eng_lev = levels(P.anal[,get(eff2plot)]))
  #   #add ordering variable
  #   curr_lev[,ord:=seq(1,.N)]
  #   #get the Hebrew translations, in the original order of levels
  #   hebrew_levels <- curr_lev[fac_lev[fac==eff2plot], on = .(eng_lev)][,.(ord,heb_lev)][order(ord)][,heb_lev]


  # predlabels <- hebrew_levels }
   if (eff2plot=="rescaled_Time.Diff_new") {
    predlabels <- as.character(as.numeric(c("2014", "2016", "2018", "2020", "2022")))
      #as.character(as.numeric(xticklabs) + 2012)
  } else if (eff2plot %in% c("Distance_rescaled","Distance_agri_rescaled")) {
    predlabels <- as.character(c("100", "1,000", "2,000", "3,000"))
  } else {
    predlabels <- as.character(c("כרמל","גליל","הרי-יהודה"))
  }

  # set the main title and filename 
  main_title <- paste("Effect of",eff2plot,sep = " ")
  filename <- paste("mammals",unique(P.anal$unit),response_var,eff2plot,sep = "_")
  if (plot_points==TRUE) {
    main_title <- paste(main_title,"(points are observations)",sep = " ")
    filename <- paste(filename,"w_observations.pdf",sep = "_")
  } else if (plot_residuals==TRUE) {
    main_title <- paste(main_title,"(points are partial residuals)",sep = " ")
    filename <- paste(filename,"w_residuals.pdf",sep = "_")
  } else {
    filename <- paste(filename,"pdf",sep = ".")
  }
  
  # if markdown is knit to pdf, reverse hebrew text
  curr_device <- knitr::opts_chunk$get("dev")
  if (export_plot || (!is.null(curr_device) && grepl("pdf",curr_device))) {
    xlabel <- strReverse(xlabel)
    ylabel <- strReverse(ylabel)
  #  if (eff2plot!="year_ct") {
  #    predlabels <- strReverse(predlabels)
  # }
   }
  
  eval(parse(text=eval(expr = paste0("effplot <- effect_plot(model = m, data=P.anal, pred = ",
                                     eff2plot,
                                     ", interval = T,plot.points = plot_points, jitter = c(0.1,0), point.size = 3, cat.pred.point.size = 4, partial.residuals = plot_residuals, colors = 'Qual1', main.title = main_title, line.colors = 'black', point.alpha = 0.25, x.label = xlabel, y.label = ylabel)"))))
  if (eff2plot=="rescaled_Time.Diff_new") {
    #All other units
    effplot <- effplot + scale_x_continuous(breaks=c(-1.233466888,-0.601840099,0.028922631,0.66054942,1.291312151), labels = predlabels)
    #y axis for negev highlands
    #effplot <- effplot + scale_y_continuous(breaks=c(0.5,1,1.5,2))
  } else if (eff2plot %in% c("Distance_rescaled","Distance_agri_rescaled")) {
    effplot <- effplot + scale_x_continuous(breaks=c(-0.820445737,0.018246686,0.950127155,1.882007625), labels = predlabels)
  } else {
    effplot <- effplot + scale_x_discrete(labels = predlabels)
  }
  
  # if font name specified then apply
  if (!is.null(fontname)) {
    effplot <- effplot + theme(text=element_text(family = fontname))
    effplot$theme$plot.title$family <- fontname
    effplot$theme$axis.title$family <- fontname
    effplot$theme$axis.text.x$family <- fontname
    effplot$theme$axis.text.y$family <- fontname
  }
  
  #apply default font size
  effplot$theme$text$size <- fontsize
  effplot$theme$plot.title$size <- fontsize
  effplot$theme$axis.title$size <- fontsize-1
  
  # justify text to right (Hebrew)
  effplot$theme$plot.title$hjust <- 1
  
  print(effplot)
  
  if (export_plot) {
    Cairo(file = paste0(outpath,filename), width = pdf_width, height = pdf_width*pdf_aspect_ratio, type = "PDF", units = "mm")
    print(effplot)
    dev.off()
  }
  
  return(effplot)

}
