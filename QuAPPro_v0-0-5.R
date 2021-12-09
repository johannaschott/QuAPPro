####################################
###     Script Shiny R 2021      ###
###  Polysome profile analysis   ###
####################################


library(shiny)
library(shinythemes)
library(colourpicker)
library(stringr)
library(shinyFiles)
library(Cairo)
library(grDevices)
library(plyr)
library(colorspace)

###################################################

### FUNCTIONS FOR APP

find_closest_minmax <- function(y, a, range)
{
  d <- diff(y)
  d[ d > 0 ] <- 1
  d[ d < 0 ] <- -1
  seg <- rle(d)
  seg_ends <- cumsum(seg$lengths) # last nt of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)] + 1) # first nt of a segment
  ## get plateaus and their middle positions
  flat <- which( seg$lengths > range & seg$values == 0 )
  flat_center <- ( seg_ends[flat] + seg_starts[flat] ) / 2
  up_down <- which( seg$lengths > range & seg$values != 0 )
  up_down_change <- abs(diff( seg$values[up_down]) ) == 2
  change_seg <- up_down[which(up_down_change) + 1]
  before_change_seg <-  up_down[which(up_down_change)]
  minmax <- (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  all_minmax <- c(minmax, flat_center) + 0.5
  closest <- order( abs(all_minmax -  a), decreasing = F)[1]
  return(all_minmax[closest])
}

closest_point <- function(x, y, X, Y){
  all_dist <- sqrt( (x-X)^2 + (y-Y)^2)
  closest <- order(all_dist, decreasing = F)[1]
  return(X[closest])
}
find_closest_inflections <- function(y, a, range)
{
  d <- diff(y, difference = 2)
  d[ d > 0 ] <- 1
  d[ d < 0 ] <- -1
  seg <- rle(d)
  seg_ends <- cumsum(seg$lengths) # last nt of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)] + 1) # first nt of a segment
  up_down <- which( seg$lengths > range & seg$values != 0 )
  up_down_change <- abs(diff( seg$values[up_down]) ) == 2
  change_seg <- up_down[which(up_down_change) + 1]
  before_change_seg <-  up_down[which(up_down_change)]
  minmax <- (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  all_minmax <- minmax + 1.5
  closest <- order( abs(all_minmax -  a), decreasing = F)[1]
  return(all_minmax[closest])
}

smooth_profile <- function(y, half_window)
{
  window <- 2*half_window + 1
  extension <- rep(NA, half_window)
  y_extended <- c(extension, y, extension)
  y_smooth <- vector()
  for(i in 1:(length(y_extended) - (window - 1) ) )
  {
    y_smooth[i] <- mean( y_extended[i:(i+window)], na.rm = T )
  }
  return(y_smooth)
}


############################
### SHINY APP:


### USER INTERFACE ###

ui <- fluidPage(
  # add title for webpage tab
  tags$head(tags$title("QuAPPro at MI3")),
  # add shiny theme 
  theme = shinythemes::shinytheme(("simplex")),
  #create title and subtitles visible in every tab
  titlePanel(h2(tags$strong("QuAPPro"), align = "center")),
  tags$h4(" - Quantification and Alignment of Polysome Profiles -", align = "center"),
  #tags$hr(tags$h5(tags$strong("Load"), "your files.",tags$strong("Select"), "your profiles.",
  #                tags$strong("Analyze"), "using the buttons.")),
  tags$head(
    tags$style(
      HTML(".shiny-notification {
           position:fixed;
           top: calc(30%);
           left: calc(10%);
           }
           ")
      )
  ),
  # create fluid layout with several tabs for displaying outputs
  tabsetPanel(
    # FIRST TAB
    tabPanel(tags$strong("Analysis profiles"), icon = icon("area-chart"),
             fluidRow(column(6, tags$h4(tags$strong("Single profiles for quantification"))),
                      column(6, tags$h4(tags$strong("Multiple aligned profiles")))),
             
             fluidRow(
               column(2,
                      # create area for uploading .pks file
                      fileInput("input_data", "Upload pks or csv file", multiple = T, accept = c(".pks",".csv")),
                      
                      #Let user select their loaded files and set x-anchor and baseline
                      selectInput("select", "Select files", choices = c(), width = '100%'),
                      
      
                      
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow(tags$h5(tags$strong("Fluorescence")),
                               
                        # show fluorescence signal or not?
                        checkboxInput("show_fl", "Show fluorescence signal", value = FALSE, width = NULL),
                        
                        column(6,
                               numericInput("axis1_fl", "Set y min", value = NULL)),
                        column(6,
                               numericInput("axis2_fl", "Set y max", value = NULL))),
                      
                      # introduce a slider for smoothing of the fluorescence signal
                      fluidRow(sliderInput("slider1", label = "Smooth profile", min = 0, 
                                           max = 100, value = 0)),
                      
                      # set a separate baseline for the fluorescence signal
                      fluidRow(
                              column(6,actionButton("baseline_fl", "Set baseline", width = '100%'))),
                      
                      fluidRow(tags$h5(tags$strong("Polysome profile")),
                               column(6,
                                      numericInput("axis1", "Set y min", value = NULL)),
                               column(6,
                                      numericInput("axis2", "Set y max", value = NULL))),
                      fluidRow(),
                      fluidRow(
                        column(6,
                               numericInput("axis3", "Set x min", value = NULL)),
                        column(6,
                               numericInput("axis4", "Set x max", value = NULL))),
                      
                      fluidRow(
                               column(6,actionButton("baseline", "Set baseline", width = '100%')),
                               column(6,actionButton("x_anchor", "Set x-anchor", width = '100%')) ),
                      # checkbox determining if lines are displayed in the plot
                      checkboxInput("red_lines", "Show baseline and x-anchor", value = TRUE, width = NULL)),
               # display plot area with click option
               column(4, plotOutput("plot_single", click = "click")),
               column(4, plotOutput("plot_align")),
               column(2,
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow(tags$h5(tags$strong("Fluorescence")),
                               
                               # show fluorescence signal or not?
                                checkboxInput("show_fl_al", "Show fluorescence signal", value = FALSE, width = NULL),
                               
                               column(6,
                                      numericInput("axis1_a_fl", "Set y min", value = NULL)),
                               column(6,
                                      numericInput("axis2_a_fl", "Set y max", value = NULL))),
                      fluidRow(tags$h5(tags$strong("Polysome profile")),
                               column(6,
                                      numericInput("axis1_a", "Set y min", value = NULL)),
                               column(6,
                                      numericInput("axis2_a", "Set y max", value = NULL))),
                      fluidRow(),
                      fluidRow(
                        column(6,
                               numericInput("axis3_a", "Set x min", value = NULL)),
                        column(6,
                               numericInput("axis4_a", "Set x max", value = NULL))))),
             
             
             # Create new area directly below the plot
             fluidRow(
               # Quantification options on the left
               column(3,
                      tags$h4(tags$strong("Quantification")),
                      tags$h6("(Baseline needs to be set)"),
                      fluidRow(column(6, actionButton("file_start", "Set area start", icon = icon("caret-square-right"), width = '100%')),
                               column(6,actionButton("file_end", "Set area end", icon = icon("pause-circle"), width = '100%'))),
                      selectInput("select_area", "Select area to quantify", choices = c("Total", "Monosomes", "Polysomes", "40S", "60S"), width = '100%'),
                      actionButton("quantify_area", "Quantify selected area", icon = icon("calculator"), width = '100%'),
                      checkboxInput("green_lines", "Show area and lines", value = TRUE, width = NULL),
                      textInput("name_area", "Optional: Name new area", value = "", width = NULL, placeholder = "Unknown peak"),
                      actionButton("take_over_name", "Add area name to list", icon = icon("check-circle"))
               ),
               column(3,
                      radioButtons("helper_functions", "Use helper functions for setting vertical lines:",
                                   c("Help find min/max" = "peak_help",
                                     "Help find inflections" = "inflection_help",
                                     "No usage of helper functions" = "no_help"
                                   ), selected = "peak_help")),
               # Alignment options in the middle
               column(3,
                      # let user select files ready to align 
                      tags$h4(tags$strong("Alignment")),
                      tags$h6("(Baseline and x-anchor need to be set)"),
                      selectInput("select_alignment", "Files ready to align", choices = c(), width = '100%'),
                      # actionButton("align", "Align profiles", icon = icon("align-center"), width = '100%'),
                      fluidRow(
                        column(6, actionButton("up", "up", icon = icon("angle-double-up"), width = '100%')),
                        column(6, actionButton("down", "down", icon = icon("angle-double-down"), width = '100%'))),
                      # aligne profiles by pressing the button
                      #changing their color and linetype seing the changes directly in the plot
                      
                      
                      colourInput("color", "Change colors for aligned profiles", palette = "square"),
                      selectInput("linetype", "Change linetypes for aligned profiles", choices = c("solid", "dashed", "dotted",
                                                                                                   "dotdash", "longdash", "twodash"), width = '100%'),
                      selectInput("linewidth", "Change linewidth for aligned profiles", choices = c("thin", "medium", "bold"), width = '100%'),
                      fluidRow(column(6, actionButton("remove_file", "Remove file", icon = icon("trash"), width = '100%')),
                               column(6, actionButton("add_file", "Restore file", icon = icon("trash-restore"), width = '100%')))),
               # Plot display options on the right with plot download button
               column(3,
                      #checkboxInput("peak_help", "Help find max/min for setting lines", value = TRUE, width = NULL),
                      checkboxInput("anchor_line", "Display x-anchor in alignment", value = TRUE, width = NULL),
                      checkboxInput("normalize_height", HTML("Normalize <b>height</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL),
                      checkboxInput("normalize_length", HTML("Normalize <b>length</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL),
                      tags$h4(tags$strong("Download plot")),
                      #shinyDirButton("dir", "Input directory", "Upload"),
                      #verbatimTextOutput("dir", placeholder = TRUE),
                      #textInput("filename_user", "Filename", value = "", width = NULL, placeholder = NULL),
                      downloadButton("downloadPlot", "Download Alignment plot", icon = icon("file-download"), width = '100%'),
                      downloadButton("downloadSinglePlot", "Download Single plot", icon = icon("file-download"), width = '100%')
               ))
    ),
    # SECOND TAB
    # shows updating table of quantified areas for respective plots, download possible 
    tabPanel(tags$strong("Quantification summary"), icon = icon("list-alt"), 
             tags$h4("Table of quantified polysome profile areas"),
             tableOutput("quantification"),
             tableOutput("ad_quantification"),
             #textInput("filen_quant_user", "Insert quantification filename", value = "", width = NULL, placeholder = NULL),
             downloadButton("downloadQuant", "Download .csv file")),
    # THIRD TAB
    # shows updating table of aligned (and normalized) profiles, download possible 
    tabPanel(tags$strong("Alignment table"), icon = icon("table"), 
             tags$h4("Table of all aligned (and normalized) profiles"),
             downloadButton("downloadData", "Download .csv file"),
             tableOutput("csv_file")),
    # FOURTH TAB 
    # general information and impressum
    tabPanel(tags$strong("Contact", style = "color:blue"), icon = icon("exclamation-circle", style = "color:blue"), 
             tags$h4("About Us"),
             tags$div("This shiny app was developed by Chiara Schiller and",
                      tags$a(href="https://www.umm.uni-heidelberg.de/biochemie/research/research-johanna-schott/", "Johanna Schott", style = "color:blue"),
                      "at the",
                      tags$a(href="https://www.umm.uni-heidelberg.de/medical-faculty-mannheim/home/", "Medical Faculty Mannheim", style = "color:blue"),
                      "(", 
                      tags$a(href="https://www.umm.uni-heidelberg.de/mi3/", "MI3", style = "color:blue"), 
                      ",",
                      tags$a(href="https://www.umm.uni-heidelberg.de/biochemie/", "Stoecklin lab", style = "color:blue"),
                      ")."
             ),
             tags$h4("Acknowledgements"),
             tags$div("We would like to thank Dr. Andreas Bohne-Lang at the",
                      tags$a(href="https://www.umm.uni-heidelberg.de/fakultaet/edv/", "IT department", style = "color:blue"),
                      "of the Medical Faculty Mannheim and all members of the Stoecklin lab for support and input."),
             tags$h4(tags$a(href="https://www.umm.uni-heidelberg.de/impressum/", "Impressum", style = "color:blue") )
    )
  ))


#### SERVER ####

server <- function(input, output, session) {
  
  ### REACTIVE VALUES FOR LISTS BUILT ALONG THE SESSION
  # create reactive values for collecting values in a list/vector/data.frame
  
  val <- reactiveValues(
    paths_collected = vector()
  )
  val <- reactiveValues(
    colors_collected = list()
  )
  val <- reactiveValues(
    linetype_collected = list()
  )
  val <- reactiveValues(
    linewidth_collected = list()
  )
  val <- reactiveValues(
    file_starts = list()
  )
  val <- reactiveValues(
    file_ends = list()
  )
  val <- reactiveValues(
    sum_areas = list()
  )
  val <- reactiveValues(
    sum_additive_areas = list()
  )
  val <- reactiveValues(
    file_starts_add = list()
  )
  val <- reactiveValues(
    file_ends_add = list()
  )
  val <- reactiveValues(
    files_to_align = vector()
  )
  val <- reactiveValues(
    remove_files_list = vector()
  )
  val <- reactiveValues(
    csv_file_df = data.frame()
  )
  val <- reactiveValues(
    factors_list = list()
  )
  
  ### INITIAL LOADED FILES
  # update select output field with the name of every newly loaded file
  observe({
    session$userData$files_list <- c(input$input_data$name, session$userData$files_list)
    updateSelectInput(session, "select",
                      choices = session$userData$files_list)
  })
  observe({
    req(input$select)
    
    if((grepl(".pks",as.character(input$select)) == T)){
      val$factors_list[[input$select]] <- (0.1/60)}else if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
        val$factors_list[[input$select]] <- (0.32/60)}else if((grepl(".csv",as.character(input$select)) == T) & !("SampleFluor" %in% colnames(file_plot()))){
          val$factors_list[[input$select]] <- (0.2/60)
          print(val$factors_list[[input$select]])
        }
    })
  # when an input file is selected store the path to the selected file name as current_path
  current_path <- reactive({
    val$paths_collected[input$input_data$name] <- input$input_data$datapath
    val$paths_collected[input$select]
  })
  
  # read selected datapath from the current selected file (input$select) for input data storage
  file_plot <- reactive({
    #require that the input is available
    req(input$select) 
    if(grepl(".pks",as.character(input$select)) == T){
      read.table(current_path(), dec = ",", header = F)}else if(grepl(".csv",as.character(input$select)) == T){
        start <- grep("Data Columns:", readLines(current_path()))
        read.csv(current_path(), skip = start)}
    
  })
  
  # Setting y and x column values for plotting from the loaded dataset
  yvalue <- reactive({
    req(input$select) 
    
    if((grepl(".pks",as.character(input$select)) == T)){
      file_plot()[ ,3]}else if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
        file_plot()[ ,5]}else if((grepl(".csv",as.character(input$select)) == T) & (ncol(file_plot()) == 4)){
          head(file_plot()[,2], -1)
        }
  })
  
  fluorescence <- reactive({
    req(input$select) 
    req(input$show_fl)
    # one condition missing!!
    if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
      smooth_profile( file_plot()[ ,3], input$slider1 )}
    
  })
  
  xvalue <- reactive({
    req(input$select) 
    if(grepl(".pks",as.character(input$select)) == T){
      (file_plot()[ ,1]+1)*val$factors_list[[input$select]]}else if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
        (1:length((file_plot()[ ,4])))*val$factors_list[[input$select]]
      }else if((grepl(".csv",as.character(input$select)) == T) & !("SampleFluor" %in% colnames(file_plot()))){
        (head(1:length(file_plot()[,1]),-1))*val$factors_list[[input$select]]
      }
    
    
  })
  
  # Create axis limits for plot area, starting with 0 and max limit values of the plotted dataset
  
  # axis labels are modified if there are more than three numerals
  lost_num_fl <- reactive({
    max_numeral <- max( floor(log10(abs( fluorescence() ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_pol <- reactive({
    max_numeral <- max( floor(log10(abs( yvalue() ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When user selects other axis limits, limit values change (works for both single and aligned (+ normalized) plots)
  ymin_single <- reactive({
    if(is.na(input$axis1)){
      0
    }else{
      input$axis1*lost_num_pol()
    }})
  
  ymin_single_fl <- reactive({
    if(is.na(input$axis1_fl)){
      0
    }else{
      input$axis1_fl*lost_num_fl()
    }})
  
  ymax_single <- reactive({
    if(is.na(input$axis2)){
      max(yvalue())
    }else{
      input$axis2*lost_num_pol()
    }})
  
  ymax_single_fl <- reactive({
    if(is.na(input$axis2_fl)){
      max(fluorescence())
    }else{
      input$axis2_fl*lost_num_fl()
    }})
  
  xmin_single <- reactive({
    if(is.na(input$axis3)){
      0
    }else{
      input$axis3
    }})

  xmax_single <- reactive({
    if(is.na(input$axis4)){
      max(xvalue())
    }else{
      input$axis4
    }})
  
  #### ALIGNMENT:
  
  # value of val$buttons changes by clicking on different action buttons triggering the respective action(s).
  # initial value of val$buttons and value after selecting a new file is reset to 10 (value 10 not used and also not NA)
  val <- reactiveValues(
    buttons = 10)
  observeEvent(input$select,{
    val$buttons = 10
  })
  
  # specific values set for val$buttons after clicking on specific buttons
  observeEvent(input$x_anchor, {val$buttons = 1})
  observeEvent(input$baseline, {val$buttons = 2})
  observeEvent(input$align, {val$buttons = 3})
  observeEvent(input$file_start, {val$buttons = 4})
  observeEvent(input$file_end, {val$buttons = 5})
  observeEvent(input$baseline_fl, {val$buttons = 6})
  
  # setting profile start and end, x-anchor and baseline triggered by clicking into the plot.
  # Values get stored as list or vector
  # Either function for helping to find max and min is used or not (tick box input$peak_help)
  observeEvent(input$click$x, {
    # only runs when button was pressed
    req(val$buttons)
    # setting profile start
    if(val$buttons == 4)
      if(input$helper_functions == "peak_help"){
        val$file_starts[[input$select]] <- (find_closest_minmax(smooth_profile(yvalue(), 2),
                                                                (closest_point(input$click$x, 
                                                                               input$click$y, 
                                                                               xvalue(),
                                                                               yvalue()))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$file_starts[[input$select]] <- (find_closest_inflections(smooth_profile(yvalue(), 2),
                                                                     (closest_point(input$click$x, 
                                                                                    input$click$y, 
                                                                                    xvalue(),
                                                                                    yvalue()))/val$factors_list[[input$select]], 5))
      }else{
        val$file_starts[[input$select]] <-(closest_point(input$click$x,
                                                         input$click$y,
                                                         xvalue(),
                                                         yvalue())/val$factors_list[[input$select]])
      }
    # setting profile ends
    if(val$buttons == 5)
      if(input$helper_functions == "peak_help"){
        val$file_ends[[input$select]] <- (find_closest_minmax(smooth_profile(yvalue(), 2),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             xvalue(),
                                                                             yvalue()))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$file_ends[[input$select]] <- (find_closest_inflections(smooth_profile(yvalue(), 2),
                                                                   (closest_point(input$click$x, 
                                                                                  input$click$y, 
                                                                                  xvalue(),
                                                                                  yvalue()))/val$factors_list[[input$select]], 5))
      }else{
        val$file_ends[[input$select]] <-(closest_point(input$click$x,
                                                       input$click$y,
                                                       xvalue(),
                                                       yvalue())/val$factors_list[[input$select]])
      }
    # setting baseline for polysome profile
    if(val$buttons == 2){
      val$baseline[input$select] <- yvalue()[(input$click$x)/val$factors_list[[input$select]]]
    }
    
    # setting baseline for fluorescence profile
    if(val$buttons == 6){
      par(xpd = T)
      val$baseline_fl[input$select] <- fluorescence()[(input$click$x)/val$factors_list[[input$select]]]
    }
    
    # setting x-anchor
    if(val$buttons == 1)
      if(input$helper_functions == "peak_help"){
        val$anchor[input$select] <- round(find_closest_minmax(smooth_profile(yvalue(), 2),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             xvalue(),
                                                                             yvalue()))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$anchor[input$select] <- (find_closest_inflections(smooth_profile(yvalue(), 2),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             xvalue(),
                                                                             yvalue()))/val$factors_list[[input$select]], 5))
      }else{
        val$anchor[input$select] <-(closest_point(input$click$x,
                                                  input$click$y,
                                                  xvalue(),
                                                  yvalue())/val$factors_list[[input$select]])
      }
    
    # if both, anchor and baseline, exist for one plot, select_alignment button gets updated with this file
    session$userData$align_files_list <- c((intersect(names(val$baseline),
                                                      names(val$anchor))),
                                           session$userData$align_files_list)
    updateSelectInput(session, "select_alignment",
                      choices = session$userData$align_files_list)
    # reactive vector with respective file names gets updated
    val$files_to_align <- intersect(names(val$baseline), names(val$anchor))
  })
  
  # update choices of areas to be quantified with names given by the user
  observeEvent(input$take_over_name,{
    session$userData$area_list <- c(input$name_area, session$userData$area_list)
    updateSelectInput(session, "select_area",
                      choices = c(session$userData$area_list, "Total", "Monosomes", "Polysomes", "40S", "60S"))
  })
  
  # show notification if one of the 3 required values is missing for area quantification, but quantify button is pressed
  observeEvent(input$quantify_area,{
    # only execute when all required values are set
    if(is.null(val$file_starts[[input$select]])|is.null(val$file_ends[[input$select]])|is.null(val$baseline[input$select])){
      showNotification(paste("Please be sure to have set start, end and baseline of the area."),
                       duration = NULL, type = "message")
    }
  })
  
  # Calculate sum of yvalues to quantify areas when button is pressed and required values were set before
  observeEvent(input$quantify_area,{
    req(val$file_starts[[input$select]], val$file_ends[[input$select]], val$baseline[input$select])
    x_first <- round(val$file_starts[[input$select]])*val$factors_list[[input$select]]
    x_last <- round(val$file_ends[[input$select]])*val$factors_list[[input$select]]
    # sum up all values from start to end area and subtract baseline, assign area name selected by user
    # store in list of areas including quantification with respective file name
    val$sum_areas[[input$select_area]][input$select] <- (sum(yvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))]-val$baseline[input$select]))
    
    if(input$show_fl){
    val$sum_areas[[paste(input$select_area, "_fluo", sep = "")]][input$select] <- (sum(fluorescence()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))]-val$baseline_fl[input$select]))
    }
    # create data frame with quantified areas containing NAs for files without respective area
    non.null.list <- lapply(val$sum_areas, lapply, function(x)ifelse(is.null(x), NA, x))
    df_quant <- rbind.fill(lapply(non.null.list, as.data.frame))
    df_quant <- as.data.frame(df_quant, row.names = names(val$sum_areas))
    val$df_quant <- df_quant
    
    
    
    # store not only area but also start and stop values in the same ways as areas!
    val$area_starts[[input$select_area]][input$select] <- val$file_starts[[input$select]]
    val$area_ends[[input$select_area]][input$select] <- val$file_ends[[input$select]]
  })
  #select_area
  
  # create "remove vector" containing files selected by remove button. Only take files into the vector that are not already there
  observeEvent(input$remove_file,{
    req(val$files_to_align)
    if(sum(str_detect(val$remove_files_list,
                      intersect(names(val$baseline[input$select_alignment]),
                                names(val$anchor[input$select_alignment]))))== 0)
      val$remove_files_list <- c((intersect(names(val$baseline[input$select_alignment]),
                                            names(val$anchor[input$select_alignment]))), val$remove_files_list)
  })
  # Create "add vector" for removing files again from remove list (so addition to alignment again)
  observeEvent(input$add_file,{
    req(val$files_to_align)
    add_vector <- c((intersect(names(val$baseline[input$select_alignment]),
                               names(val$anchor[input$select_alignment]))), val$add_files_list)
    val$remove_files_list <- val$remove_files_list[!val$remove_files_list %in% add_vector] 
  })
  
  # Let user change colors of the different plots (selected out of aligned files)
  observeEvent(input$color, {
    if(input$color != "#FFFFFF"){
      val$colors_collected[[input$select_alignment]] <- input$color}
  })
  
  # Let user change linetypes of the different plots (selected out of aligned files)
  observeEvent(input$linetype, {
    if(input$linetype == "solid")line_selected <- 1
    if(input$linetype == "dashed")line_selected <- 2
    if(input$linetype == "dotted")line_selected <- 3
    if(input$linetype == "dotdash")line_selected <- 4
    if(input$linetype == "longdash")line_selected <- 5
    if(input$linetype == "twodash")line_selected <- 6
    val$linetype_collected[[input$select_alignment]] <- line_selected
  })
  
  # Let user change linewidths of the different plots (selected out of aligned files)
  observeEvent(input$linewidth, {
    if(input$linewidth == "thin")linewidth_selected <- 1
    if(input$linewidth == "medium")linewidth_selected <- 2
    if(input$linewidth == "bold")linewidth_selected <- 3
    val$linewidth_collected[[input$select_alignment]] <- linewidth_selected
  })
  
  ### OUTPUT
  
  # show notification if align is pressed but no files to align were set (by x- anchor and baseline)
  observeEvent(input$align,{
    if(is.null(val$files_to_align)){
      showNotification(paste("Please be sure to have set x-anchor and baseline of the plots you want to align"),
                       duration = Null, type = "error")
    }
  })
  
  # show notification if show fluorescence is selected but there is no fluorescence signal
  observeEvent(input$show_fl,{
    if(is.null(fluorescence())){
      showNotification(paste("Your file does not contain the column SampleFluor."),
                       duration = NULL, type = "error")
    }
  })
  
  ## plot individual profiles
  plot_singleFl <-function(){
    if(lost_num_fl() == 1 ){
      ylab <- "Fluorescence"
    }else{
      ylab <- paste("Fluoresc. (x ", lost_num_fl(), ")", sep = "")
    }
    
    par(mar = c(0, 4, 0.5, 2)) 
    plot(xvalue(), fluorescence(), type = "l", xaxt = "n", col = "darkgreen", 
         xlim =c(xmin_single(),xmax_single()), ylim = c(ymin_single_fl(), ymax_single_fl()),
         las = 1, ylab = ylab, mgp = c(2, 0.6, 0), xlab = "", yaxt = "n")
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_fl(), las = 1, mgp = c(2, 0.6, 0))
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      polygon(c(x_first, xvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], x_last),
              c(val$baseline_fl[input$select], fluorescence()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], val$baseline_fl[input$select]),
              col = "#c7e9c0", border = "darkgreen")
    }
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2)
    }
    if(input$red_lines){
      abline(h = val$baseline_fl[input$select],col = "red", lty=2)
    }
  }
  
  plot_singlePol <-function(){
    if(lost_num_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_pol(), ")", sep = "")
    }
    
    par(mar = c(4, 4, 0, 2))
    plot(xvalue(), yvalue(), type = "l", las = 1,
         ylab = ylab, xlab = "time (min)",
         ylim = c(ymin_single(),ymax_single()), 
         xlim =c(xmin_single(),xmax_single()), mgp = c(2, 0.6, 0),
         yaxt = "n"
    )
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = c(2, 0.6, 0))
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      polygon(c(x_first, xvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], x_last),
              c(val$baseline[input$select], yvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], val$baseline[input$select]),
              col = "#c7e9c0", border = "black")
    }
    # selected x-anchor and baseline are displayed if box is ticked
    if(input$red_lines){
      abline(v = val$anchor[input$select]*val$factors_list[[input$select]],col = "red", lty=2)
      abline(h = val$baseline[input$select],col = "red", lty=2)
    }
    # selected area starts and ends are displayed if box is ticked
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2)
    }
  }


  
  plot_singleInput <- function(){
    if(!is.null(xvalue())){
      if(input$show_fl){
        if(val$buttons == 6){
          layout(matrix(2:1, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
          plot_singlePol()
          plot_singleFl()
        }else{
          layout(matrix(1:2, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
          plot_singleFl()
          plot_singlePol()
        }
      }else{
        plot_singlePol()
      }
    }
  }
  
  output$plot_single <- renderPlot({
    plot_singleInput()
  })
  
  # Enable download of current plot as pdf
  
  output$downloadSinglePlot <- downloadHandler(
    filename = function(){ sub(pattern = "(.*?)\\..*$", replacement = "\\1", as.character(input$select) ) },
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_singleInput() )
      dev.off()
    })  
  
  # create table output of aligned files with option to download the table with given name
  output$csv_file <- renderTable(
    val$csv_file_df
  )
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$filename_user,".csv", sep = "")
    },
    content = function(file) {
      write.csv2(val$csv_file_df, file, row.names = FALSE)
    }
  )
  
  # create output table showing quantification data with option to download table with given name
  output$quantification <- renderTable(
    t(val$df_quant), rownames = T
  )
  
  output$downloadQuant <- downloadHandler(
    filename = function() {
      paste(input$filen_quant_user,".csv", sep = "")
    },
    content = function(file) {
      write.csv2( t(val$df_quant), file, row.names = T)
    }
  )




  # Create alignment plot
  
  # Opportunity to remove and add files again to files_to_plot vector
  files_to_plot <- reactive({
    val$files_to_align[!val$files_to_align %in% val$remove_files_list]
  })
    
  # move files up or down in the val$files_to_align vector
  observeEvent(input$up, {
    to_be_shifted <- which( val$files_to_align == input$select_alignment )
    if(to_be_shifted > 1)
    {
      files_order <- 1:length(val$files_to_align)
      files_order[to_be_shifted] <- files_order[to_be_shifted] - 1
      files_order[to_be_shifted - 1] <- files_order[to_be_shifted - 1] + 1
      val$files_to_align <- val$files_to_align[files_order]
    }
  })
  
  # move files up or down in the val$files_to_align vector
  observeEvent(input$down, {
    to_be_shifted <- which( val$files_to_align == input$select_alignment )
    if(to_be_shifted < length(val$files_to_align))
    {
      files_order <- 1:length(val$files_to_align)
      files_order[to_be_shifted] <- files_order[to_be_shifted] + 1
      files_order[to_be_shifted + 1] <- files_order[to_be_shifted + 1] - 1
      val$files_to_align <- val$files_to_align[files_order]
    }
  })
   
  # create color values, rainbow palette as initial colors, further replaced by selected color if selected
  colors_vector <- reactive({
    dummy <- qualitative_hcl(length(files_to_plot()), palette = "Dark 3")
    names(dummy) <- sort(files_to_plot())
    if (length(val$colors_collected) > 0){
      dummy[names(val$colors_collected) ] <- unlist(val$colors_collected)
    }
    dummy
  })
  
  # create linetype values, solid lines as initial linetype, further replaced by selected linetype
  lines_vector <- reactive({
    dummy <- rep(1, length(files_to_plot()))
    names(dummy) <- files_to_plot()
    if (length(val$linetype_collected) > 0){
      dummy[names(val$linetype_collected)] <- unlist(val$linetype_collected)
    }
    dummy
  })
  
  # create linewidth vector, linewidth 1 as initial linewidth, further replaced by selected linetype
  linewidth_vector <- reactive({
    dummy <- rep(1, length(files_to_plot()))
    names(dummy) <- files_to_plot()
    if (length(val$linewidth_collected) > 0){
      dummy[names(val$linewidth_collected)] <- unlist(val$linewidth_collected)
    }
    dummy
  })
  
  # for normalization of surfaces:
  norm_factor <- reactive({
    dummy <- vector()
    for(f in files_to_plot()){
      if((!is.null(val$sum_areas[["Total"]][f])) & input$normalize_height){
        # normalization y-values
        dummy[f] <- val$sum_areas[["Total"]][f]/val$sum_areas[["Total"]][files_to_plot()[1]]
      }else{
        dummy[f] <- 1
      }
    }
    dummy
  })
  
  # for normalization of length:
  norm_factor_x <- reactive({
    dummy <- vector()
    for(f in files_to_plot()){
      if(!is.null(val$sum_areas[["Total"]][f]) & input$normalize_length){
        dummy[f] <- (round(val$area_ends[["Total"]][f] - val$area_starts[["Total"]][f])/(round(val$area_ends[["Total"]][files_to_plot()[1]] - val$area_starts[["Total"]][files_to_plot()[1]])))
      }else{
        dummy[f] <- 1
      }
    }
    dummy
  })
  
  # determine normalized y-values of fluorescence
  values_fluorescence <- reactive({
    dummy <- list()
    for(f in files_to_plot())
    {
      start <- grep("Data Columns:", readLines(val$paths_collected[f]))
      fluo <- read.csv(val$paths_collected[f], skip = start)$SampleFluor
      dummy[[f]] <- ( ( fluo - val$baseline_fl[[f]] )/norm_factor()[f])*norm_factor_x()[f]
    }
    dummy
  })
  
  
  values_list <- reactive({
    dummy <- list()
    for(f in files_to_plot())
    {
      if(val$factors_list[[f]] == (0.1/60)){
        next_y <- read.table(val$paths_collected[f], dec = ",", header = F)$V3
      }else if(val$factors_list[[f]] == (4.8/15)/60){
        start <- grep("Data Columns:", readLines(val$paths_collected[f]))
        next_y <- read.csv(val$paths_collected[f], skip = start)$AbsA
      }else if(val$factors_list[f] == (0.2)/60){
        start <- grep("Data Columns:", readLines(val$paths_collected[f]))
        next_y <- head(read.csv(val$paths_collected[f], skip = start)$Absorbance, -1)
      }
      
      dummy[[f]] <- ( ( next_y - val$baseline[f])/norm_factor()[f])*norm_factor_x()[f]
    }
    dummy
  })
  
  
  # determine shift along x-axis for all:
  shifts <- reactive({
    # find max anchor of all files to be aligned
    all_anchors <- (val$anchor[files_to_plot()])/norm_factor_x()[files_to_plot()]
    super_anchor <- max(all_anchors)
    dummy <- super_anchor - all_anchors
    names(dummy) <- files_to_plot()
    dummy
  })
  
  # find starts and ends of aligned profiles
  aligned_starts <- reactive({
    dummy <- ((1/norm_factor_x()[files_to_plot()] + (shifts())))
    names(dummy) <- files_to_plot()
    dummy
  })
  
  aligned_ends <- reactive({
    # store lengths of the profiles 
    profile_lengths <- sapply(values_list(), FUN = length)
    names(profile_lengths) <- files_to_plot()
    dummy <- ((profile_lengths/norm_factor_x()[files_to_plot()] + (shifts())))
    names(dummy) <- files_to_plot()
    dummy
  })
  
  xmax_all <- reactive({
    max(aligned_ends())
  })
  
  # find maximum of y-axis values for polysome profiles
  ymax_all <- reactive({
    profile_heights <- sapply(values_list()[files_to_plot()], FUN = max)
    max(profile_heights )
  })
    
  # find maximum of y-axis values for fluorescence profiles
  ymax_all_fl <- reactive({
    profile_heights <- sapply(values_fluorescence()[files_to_plot()], FUN = max)
    max(profile_heights )
  })
  
  # axis labels are modified if there are more than three numerals
  lost_num_al_fl <- reactive({
    max_numeral <- max( floor(log10(abs( unlist(values_fluorescence()[files_to_plot()] ) ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_al_pol <- reactive({
    max_numeral <- max( floor(log10(abs( unlist(values_list()[files_to_plot()]) ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When user selects other axis limits, limit values change 
  # (works for both single and aligned (+ normalized) plots)
  ymin <- reactive({
    if(is.na(input$axis1_a)){
      0
    }else{
      input$axis1_a
    }})
  
  ymax <- reactive({
    if(is.na(input$axis2_a)){
      ymax_all()  
    }else{
      input$axis2_a*lost_num_al_pol()
    }})
  
  ymin_fl <- reactive({
    if(is.na(input$axis1_a_fl)){
      0
    }else{
      input$axis1_a_fl
    }})
  
  ymax_fl <- reactive({
    if(is.na(input$axis2_a_fl)){
      ymax_all_fl()  
    }else{
      input$axis2_a_fl*lost_num_al_fl()
    }})
  
  xmin <- reactive({
    if(is.na(input$axis3_a)){
      0
    }else{
      input$axis3_a
    }})
  
  xmax <- reactive({
    if(is.na(input$axis4_a)){
      xmax_all()   
    }else{
      input$axis4_a
    }
  })    
  
  
  #### plotting function for aligned profiles
  
  plot_alignedPol <- function(){
    
    if(lost_num_al_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_al_pol(), ")", sep = "")
    }
    
    f <- files_to_plot()[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    par(mar = c(4, 4, 0, 2)) 
    plot(x, values_list()[[f]], type = "l", lty = lines_vector()[f],
         lwd = linewidth_vector()[f],
         ylab = ylab, xlab = "Index", las = 1,
         col = colors_vector()[f], mgp = c(2, 0.6, 0), 
         ylim = c(ymin(),ymax()), xlim = c(xmin(),xmax()),
         yaxt = "n"
    )
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_pol(), las = 1, mgp = c(2, 0.6, 0))
    
    for(f in files_to_plot()[-1])
    {
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      points(x, values_list()[[f]], type = "l", lty = lines_vector()[f], 
           lwd = linewidth_vector()[f], col = colors_vector()[f])
    }
    legend("topright", legend = files_to_plot(), lty = lines_vector()[files_to_plot()], 
           lwd = linewidth_vector()[files_to_plot()], col = colors_vector()[files_to_plot()],
           bty = "n"
           )
  }
  
  plot_alignedFluo <- function(){
    
    if(lost_num_al_fl() == 1 ){
      ylab <- "Fluorescence"
    }else{
      ylab <- paste("Fluo. (x ", lost_num_al_fl(), ")", sep = "")
    }
    
    f <- files_to_plot()[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    par(mar = c(0, 4, 0.5, 2)) 
    plot(x, values_fluorescence()[[f]], type = "l", lty = lines_vector()[f],
         lwd = linewidth_vector()[f],
         ylab = ylab, xlab = "", las = 1,
         col = colors_vector()[f], mgp = c(2, 0.6, 0), 
         ylim = c(ymin_fl(),ymax_fl()), xlim = c(xmin(),xmax()),
         yaxt = "n", yaxt = "n"
    )
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_fl(), las = 1, mgp = c(2, 0.6, 0))
    
    for(f in files_to_plot()[-1])
    {
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      points(x, values_fluorescence()[[f]], type = "l", lty = lines_vector()[f], 
             lwd = linewidth_vector()[f], col = colors_vector()[f])
    }
  }

  plot_alignment <- function(){
      if(input$show_fl_al){
        layout(matrix(1:2, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
        plot_alignedFluo()
        plot_alignedPol()
      }else{
        plot_alignedPol()
      }
  }
  
  output$plot_align <- renderPlot({
    req(files_to_plot())
    if(length(files_to_plot()) >= 1)
    {
      plot_alignment()
    }
  })

  # Enable download of current plot as pdf
  
  output$downloadPlot <- downloadHandler(
    filename = "alignment",
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_alignedPol() )
      dev.off()
    })  
  
  
  
}

shinyApp(ui = ui, server = server)

###############
