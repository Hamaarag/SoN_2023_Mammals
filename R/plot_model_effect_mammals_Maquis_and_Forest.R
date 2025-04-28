plot_model_effect <- function (P.anal, m, eff2plot, plot_points=FALSE, plot_residuals=FALSE, export_plot=FALSE, ylabel = NULL, fontname = NULL,
                               fontsize=22, pdf_width=160, outpath = "../output/", point_size=3, point_alpha=0.25,
                               horz_jitter=0.1, at_list = NULL, y_cutoff = NULL) {
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
  # point_size is the size of the observation points (in case plot_points=TRUE)
  # point_alpha is the opacity of the observation points - 0 is transparent, 1 is opaque (in case plot_points=TRUE)
  # horz_jitter is the degree to which observation points should be jittered (in case plot_points=TRUE)
  # at_list: If you want to manually set the values of other variables in the model, do so by providing a named list where the names are the variables and the list values are vectors of the values. This can be useful especially when you are exploring interactions or other conditional predictions. 
  # y_cutoff is the max limit of the y-axis. Points excluded by this cutoff will be enumerated and indicated in the plot title.
  
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
  
  
  # determine what is the plotted effect and set the x label
  xlabel <- switch(eff2plot,
                   "settlements" = "יישוב",
                   "agriculture" = "חקלאות",
                   "dunes" = "סיווג הדיונה",
                   "habitat" = "בית הגידול",
                   "land_use" = "שימוש הקרקע",
                   "is_loess" = "שימוש הקרקע",
                   "Subunit" = "תת-יחידה",
                   "rescaled_Time.Diff" = "זמן",
                   "Distance_rescaled" = "מרחק מיישוב )מטר(",
                   "Distance_agri_rescaled" = "מרחק מחקלאות )מטר(")
  
  # data table with subunit names
  dt_subunits <- data.table(subunit_eng = c("Carmel","Galilee","Judea"),
                            subunit_heb = c("הכרמל","הגליל","הרי יהודה"))
  
  if (eff2plot=="rescaled_Time.Diff") {
    #Maquis and conifer forest
    predlabels <- as.character(as.numeric(c("2013", "2015", "2017", "2019", "2021")))
    #as.character(as.numeric(xticklabs) + 2012)
  } else if (eff2plot %in% c("Distance_rescaled","Distance_agri_rescaled")) {
    predlabels <- as.character(c("100", "1,000", "2,000", "3,000"))
  } else {
    predlabels <- dt_subunits[subunit_eng %in% unique(P.anal$Subunit),subunit_heb]
  }
  
  # set the main title and filename 
  main_title <- paste("Effect of",eff2plot,sep = " ")
  filename <- paste("mammals",unique(P.anal$unit),response_var,eff2plot,sep = "_")
  if (plot_points==TRUE) {
    main_title <- paste(main_title,"(points are observations)",sep = " ")
    filename <- paste(filename,"w_observations",sep = "_")
  } else if (plot_residuals==TRUE) {
    main_title <- paste(main_title,"(points are partial residuals)",sep = " ")
    filename <- paste(filename,"w_residuals",sep = "_")
  }
  
  # add y_cutoff to figure title----
  if (!is.null(y_cutoff)) {
    if (plot_points) {
      main_title <- paste0(main_title,"\n",paste(P.anal[,sum(.SD>y_cutoff, na.rm = T),.SDcols=response_var],"points ommited above",y_cutoff))
    } else if (plot_residuals) {
      main_title <- paste0(main_title,"\n",paste(as.data.table(effplot$layers[[1]]$data)[,..response_var][,sum(.SD>y_cutoff, na.rm = T),.SDcols=1],"points ommited above",y_cutoff))
    }
    filename <- paste(filename,"y_cutoff",y_cutoff,sep = "_")
  }
  
  # add extension to file name----
  filename <- paste(filename,"pdf",sep = ".")
  
  # if markdown is knit to pdf, reverse hebrew text
  curr_device <- knitr::opts_chunk$get("dev")
  if (export_plot || (!is.null(curr_device) && grepl("pdf",curr_device))) {
    xlabel <- strReverse(xlabel)
    ylabel <- strReverse(ylabel)
    if (eff2plot!="rescaled_Time.Diff" & eff2plot!="Distance_rescaled") {
      predlabels <- strReverse(predlabels)
    }
  }
  
  eval(parse(text=eval(expr = paste0("effplot <- effect_plot(model = m, data=P.anal, pred = ",
                                     eff2plot,
                                     ", interval = T,plot.points = plot_points, jitter = c(horz_jitter,0), point.size = point_size, cat.pred.point.size = 4, partial.residuals = plot_residuals, colors = 'Qual1', main.title = main_title, line.colors = 'black', point.alpha = point_alpha, x.label = xlabel, y.label = ylabel, at = at_list)"))))
  if (eff2plot=="rescaled_Time.Diff") {
    #maquis and conifer forest
    effplot <- effplot + scale_x_continuous(breaks=c(-1.496140683,-0.865377952,-0.233751163,0.397011567,1.028638356), labels = predlabels)
  } else if (eff2plot %in% c("Distance_rescaled","Distance_agri_rescaled")) {
    effplot <- effplot + scale_x_continuous(breaks=c(-0.820445737,0.018246686,0.950127155,1.882007625), labels = predlabels)
  } else {
    effplot <- effplot + scale_x_discrete(labels = predlabels)
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
  
  # apply y cutoff
  if (!is.null(y_cutoff)) {
    effplot <- effplot + ylim(0,y_cutoff)
  }
  
  print(effplot)
  
  if (export_plot) {
    Cairo(file = paste0(outpath,filename), width = pdf_width, height = pdf_width*pdf_aspect_ratio, type = "PDF", units = "mm")
    print(effplot)
    dev.off()
  }
  
  return(effplot)
  
}
