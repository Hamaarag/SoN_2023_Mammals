plot_model_interaction <- function (P.anal, m, eff2plot, modvar2plot, plot_points=TRUE, plot_residuals=FALSE, export_plot=TRUE, ylabel = NULL, fontname = "Almoni ML v5 AAA",
                                    fontsize=22, pdf_width=160, outpath = "../output/", colors = "Dark2", legend_position = "right") {
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
  # colors: if length of this argument is 1 then it is interpreted as color palette name. If longer thank 1 the it is interpreted as a vector of color names / codes, as received by interactions::interact_plot
  # legend_position is a string indicating position of the legend as accepted by theme(legend.position = legend_position), e.g. "right", "bottom", "inside", etc.
  
  # Define a function to reverse Hebrew text in case markdown is knit to pdf
  strReverse <- function(x){sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")}
  
  pdf_aspect_ratio <- 2/3
  
  require(jtools)
  require(Cairo)
  require(ggplot2)
  require(extrafont)
  
  # set colors string if colors were specified by user
  if (length(colors)==1) {
    colors <- paste0("'",colors,"'")
  } else {
    colors <- paste0('c("',paste(colors, collapse = '","'),'")')
  }
  
  # determine what is the response and set the y label
  response_var <- as.character(formula(m))[2]
  if (is.null(ylabel)) {
    ylabel <- switch(response_var,
                     "richness" = "עושר מינים",
                     "Canis.aureus" = "שפע",
                     "`Capra nubiana`" = "שפע",
                     "Gazella.dorcas" = "שפע",
                     "Hystrix.indica" = "שפע",
                     "Meles.meles" = "שפע", 
                     "Sus.scrofa" = "שפע",
                     "Vulpes.vulpes" = "שפע",
                     "Hyaena.hyaena" = "שפע",
                     "gma" = "שפע ממוצע למין",
                     "abundance" = "שפע כולל",
                     response_var)
  }
  
  # if near vs. far factor, then relevel to get Near on the left hand side
  # if (eff2plot %in% c("agriculture","settlements")) {
  #     P.anal[,c(eff2plot):=factor(get(eff2plot),levels = c("Far","Near"))]
  #   
  #      # if year covariate (continuous) then get the data labels to reverse them if needed
  #    } else if (eff2plot=="year_ct") {
  #      tmp_effplot <- interact_plot(model = m, pred = year_ct, data = P.anal)
  #      xticklabs <- ggplot_build(tmp_effplot)$layout$panel_params[[1]]$x$get_labels()
  #xticklabs <- as.character(as.numeric(c("0", "2", "4", "6", "8")))
  #  }
  
  # determine what is the plotted effect and set the x label
  xlabel <- switch(eff2plot,
                   "settlements" = "קרבה ליישובים",
                   "agriculture" = "חקלאות",
                   "dunes" = "סיווג הדיונה",
                   "habitat" = "בית הגידול",
                   "land_use" = "שימוש הקרקע",
                   "Subunit" = "תת-יחידה",
                   "year_ct" = "",                                      # do not print x label for temporal trend
                   "rescaled_Time.Diff" = "",                           # do not print x label for temporal trend
                   "Distance_rescaled" = "מרחק מיישוב",
                   "Distance_agri_rescaled" = "מרחק מחקלאות")
  
  
  if (eff2plot %in% c("settlements","agriculture","Dunes","habitat","land_use","Subunit")) {

    #  Set the x tick mark labels (labels of levels of the factor plotted)
    # create a table containing English level names and Hebrew translations
    fac_lev <- data.table(fac = c("settlements","settlements","agriculture","agriculture","dunes","dunes","habitat","habitat",
                                  "land_use","land_use","land_use","Subunit","Subunit","Subunit"),
                          eng_lev = c("Near","Far","Near","Far","Shifting","Semi-shifting","Slope","Wadi",
                                      "Beduoin Agriculture","KKL Plantings","Loess","Judean Highlands","Carmel","Galilee"),
                          heb_lev = c("קרוב","רחוק","קרוב","רחוק","נודדת","מיוצבת למחצה","מדרון","ערוץ",
                                      "חקלאות מסורתית","נטיעות קקל","לס","הרי יהודה","הכרמל","הגליל"))
    #get the levels of the current factor, in the order in which they appear
    curr_lev <- data.table(eng_lev = levels(P.anal[,get(eff2plot)]))
    #add ordering variable
    curr_lev[,ord:=seq(1,.N)]
    #get the Hebrew translations, in the original order of levels
    hebrew_levels <- curr_lev[fac_lev[fac==eff2plot], on = .(eng_lev)][,.(ord,heb_lev)][order(ord)][,heb_lev]
    predlabels <- hebrew_levels
    
  } else if (eff2plot=="rescaled_Time.Diff") {
    #Maquis and conifer forest
    predlabels <- as.character(as.numeric(c("2013", "2015", "2017", "2019", "2021")))
  } else if (eff2plot %in% c("Distance_rescaled","Distance_agri_rescaled")) {
    predlabels <- as.character(c("100", "1,000", "2,000", "3,000"))
  }
  
 
  if (modvar2plot %in% c("Distance_rescaled","Distance_agri_rescaled")){
    modxvals = 'plus-minus'
  } else if (modvar2plot=="Subunit") {
    modxvals = NULL
  }
  
  modx_title <- switch(modvar2plot,
                       "settlements" = "קרבה ליישובים",
                       "agriculture" = "חקלאות",
                       "dunes" = "סיווג הדיונה",
                       "habitat" = "בית הגידול",
                       "land_use" = "שימוש הקרקע",
                       "year_ct" = "שנה",
                       "Subunit" = "תת-יחידה",
                       "rescaled_Time.Diff" = "זמן",
                       "Distance_rescaled" = "מרחק מיישוב",
                       "Distance_agri_rescaled" = "מרחק מחקלאות")
  
  # set xmod lables
  
  if (modvar2plot %in% c("settlements","agriculture","Distance_rescaled","Distance_agri_rescaled")) { 
    modlabs = as.character(c("קרוב", "רחוק")) 
  } else if (modvar2plot=="habitat") {
    modlabs = as.character(c("ערוץ", "מדרון"))
  } else if (modvar2plot=="dunes") {
    modlabs = as.character(c("נודדות", "מיוצבת למחצה"))
  } else if (modvar2plot=="Subunit") {
    modlabs = as.character(c("הכרמל","הגליל","הרי יהודה"))
  }
  
  # set the main title and filename 
  main_title <- paste("Effect of",eff2plot,"by",modvar2plot,sep = " ")
  filename <- paste("mammals",unique(P.anal$unit),response_var,eff2plot,modvar2plot,"interaction",sep = "_")
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
    modx_title <- strReverse(modx_title)
    modlabs <- strReverse(modlabs)
    # if (eff2plot!="Distance_rescaled") {
    #   predlabels <- strReverse(predlabels)
    # }
  }
  
  eval(parse(text=eval(expr = paste0("effplot <- interact_plot(model = m, data=P.anal, pred = ",
                                     eff2plot, 
                                     ", legend.main = modx_title, modx = ", modvar2plot,
                                     ", interval = T,plot.points = plot_points, jitter = c(0.1,0), 
                                     point.size = 3, cat.pred.point.size = 4, partial.residuals = plot_residuals, 
                                     main.title = main_title, colors = ",colors,
                                     ", modxvals = modxvals, modx.labels = modlabs,
                                     line.colors = 'black', point.alpha = 0.25, x.label = xlabel, 
                                     y.label = ylabel)"))))
  
  if (eff2plot=="year_ct") {
    effplot <- effplot + scale_x_continuous(labels = predlabels)
  } else if (eff2plot=="rescaled_Time.Diff") {
    #Maquis and conifer forest
    effplot <- effplot + scale_x_continuous(breaks=c(-1.496140683,-0.865377952,-0.233751163,0.397011567,1.028638356), labels = predlabels)
  } else if (eff2plot=="Distance_rescaled") {
    effplot <- effplot + scale_x_continuous(breaks=c(-0.820445737,0.018246686,0.950127155,1.882007625), labels = predlabels)
    #effplot$layers[[2]]$geom_params$width <- 0.4
  } else { effplot <- effplot + scale_x_discrete(labels = predlabels) 
  effplot$layers[[2]]$geom_params$width <- 0.4 
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
  
  # set legend position
  effplot <- effplot + theme(legend.position = legend_position)
  if (legend_position=="bottom") {
    effplot <- effplot + guides(color = guide_legend(title.position = "right"))
  }
  
  print(effplot)
  
  if (export_plot) {
    Cairo(file = paste0(outpath,filename), width = pdf_width, height = pdf_width*pdf_aspect_ratio, type = "PDF", units = "mm")
    print(effplot)
    dev.off()
  }
  
  return(effplot)
}

