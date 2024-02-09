#####################################################
#     ____                 _____  _____             #
#    / __ \          /\   |  __ \|  __ \            #
#   | |  | |_   _   /  \  | |__) | |__) | __ ___    #
#   | |  | | | | | / /\ \ |  ___/|  ___/ '__/ _ \   #
#   | |__| | |_| |/ ____ \| |    | |   | | | (_) |  #
#    \___\_\\__,_/_/    \_\_|    |_|   |_|  \___/   #
#                                                   #                                                
#####################################################                                                 

# Note: The tabs for de-convolution and statistics can be found
# in their most recent version in QuAPPro_v0-0-6_240205.R.

# Problems and Ideas:
# Quantification when total area does not exist (line 2423) -> prevent 
# When is spar NULL? Should smoothing be performed in this case, or not?
# Determine time factor from time column
# Set minimum of profiles not to zero, but to minimum of values
# Introduce logos and message for optimal display

#### REQUIRED PACKAGES ####

library(shiny)
library(shinythemes)
library(colourpicker)
library(stringr)
library(shinyFiles)
library(plyr)
library(colorspace)
library(markdown)
library(emg)
library(DT)


#### FUNCTIONS ####

# Given coordinates of a selected point (x, y),
# this function returns the index of the closest point
# among the coordinates X and Y according to the Euclidean distance:

closest_point <- function(x, y, X, Y){
  all_dist <- sqrt( (x-X)^2 + (y-Y)^2)
  closest <- order(all_dist, decreasing = F)[1]
  return(closest)
}

# Given a set of values (y), this function identifies 
# the indices of all local minima
# that lie between at least r consecutive decreasing values
# and at least r consecutive increasing values.
all_min <- function(y, r){
  d <- diff(y )
  d0 <- d
  d0[ d > 0 ] <- 1
  d0[ d < 0 ] <- -1
  seg <- rle(d0)
  seg_ends <- cumsum(seg$lengths) + 1 # last position of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)]) # first position of a segment
  up_down <- which( seg$lengths >= r & seg$values != 0  ) # segments with increasing or decreasing values
  up_down_change <- diff( seg$values[up_down]) == 2 
  change_seg <- up_down[which(up_down_change) + 1]  # segments with a positive direction change compared to previous non-zero segment?
  before_change_seg <-  up_down[which(up_down_change)]
  all_min <- (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  return(all_min)
}

# Given a set of values (y), this function identifies 
# the indices of all local maxima
# that lie between at least r consecutive increasing values
# and at least r consecutive decreasing values.
all_max <- function(y, r){
  d <- diff(y )
  d0 <- d
  d0[ d > 0 ] <- 1
  d0[ d < 0 ] <- -1
  seg <- rle(d0)
  seg_ends <- cumsum(seg$lengths) + 1 # last position of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)]) # first position of a segment
  up_down <- which( seg$lengths >= r & seg$values != 0  ) # segments with increasing or decreasing values
  up_down_change <- diff( seg$values[up_down]) == -2 
  change_seg <- up_down[which(up_down_change) + 1]  # segments with a negative direction change compared to previous non-zero segment?
  before_change_seg <-  up_down[which(up_down_change)]
  all_max <- (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  return(all_max)
}

# This function identifies the closest local minimum or maximum 
# in a vector of values (y) to a selected index a of the vector.
find_closest_minmax <- function(y, a, r){
  minmax <- c(all_min(y, r), all_max(y, r) )
  closest <- order( abs(minmax -  a), decreasing = F )[1]
  return(minmax[closest])
}

# This function returns the "second derivative" of a vector of values (y).
second_deriv <- function(y){
  return( diff(y, differences = 2) )
}

# This function applies smooth.spline() to a given vector of values y
# with the smoothing parameter spar. 
# If spar is zero, the values are returned without smoothing.
smooth_profile <- function(y, spar){
  if(is.null(spar) )  {
    return( smooth.spline(1:length(y), y, spar = spar)$y )
  }else{
    if(spar == 0){
      return(y)
    }else{
      return( smooth.spline(1:length(y), y, spar = spar)$y )
    }
  }
}

#### USER INTERFACE ####

ui <- fluidPage(
  # Browser tab title:
  title = tags$head(tags$title("QuAPPro at MI3")),
  
  # Header with logo and title:
  fluidRow(
    column(2),
    column(8, 
           titlePanel(
             div( h2(tags$strong("QuAPPro"), align = "center" ),
                  h4("- Quantification and Alignment of Polysome Profiles -", align = "center" )
             )
           )
    ),
    column(2,
           style = "padding-right:10px;padding-top:5px",
           tags$img(src='UMM_MEDMA.png', align = "right", width = '150%', height = '150%')
    )
  ),
  
  # Further settings:
  # Add shiny theme: 
  theme = shinythemes::shinytheme(("simplex")),
  # Settings of notifications (error messages and warnings):
  tags$head(
    tags$style(
      HTML(".shiny-notification {
           position:fixed;
           top: calc(50%);
           left: calc(50%);
           }
           ")
    )
  ),
  
  # Create fluid layout with several tabs:
  tabsetPanel(
    # FIRST TAB: ALIGNMENT AND QUANTIFICATION OF AREAS
    tabPanel(tags$strong("Alignments/Areas"), icon = icon("chart-area"),
             fluidRow(column(2),
                      column(4, tags$h4(tags$strong("Individual profiles for quantification") ) ),
                      column(6, tags$h4(tags$strong("Aligned profiles") ) )
             ),
             fluidRow(
               # Section for uploading and selecting files:
               column(2, style = "padding-left:20px",
                      tags$h4(tags$strong("File import and export")),
                      
                      # fileInput for uploading files
                      fileInput("input_data", "Load profiles or QuAPPro analysis", multiple = T, accept = c(".pks",".csv", ".txt", ".RData")),
                      #textOutput("test"),
    
                      # Download entire analysis
                      downloadButton("export", "Export analysis", width = "100%"),
                      
                      # Menu for selecting which of the uploaded files should be displayed
                      fluidRow(style = "padding-left:15px;padding-top:20px;padding-right:15px",
                      selectInput("select", "Select a profile", choices = c(), width = '100%')
                      ),
                      
                      # Section for setting baseline and x-anchor of individual polysome profiles
                      tags$h4(tags$strong("Polysome profile")),
                      fluidRow(
                        column(6, style = "padding-top:10px", 
                               actionButton("baseline", "Baseline", width = '100%')
                        ),
                        column(6,style = "padding-top:10px", 
                               actionButton("x_anchor", "X-anchor", width = '100%')
                        )
                      ),                     
                      
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               numericInput("axis1", "Set y min", value = NULL, step = 1) ),
                        column(6, style = "padding-top:10px",
                               numericInput("axis2", "Set y max", value = NULL, step = 1) )
                      ),
                      
                      fluidRow(
                        column(6, style = "padding-top:0px",
                               numericInput("axis3", "Set x min", value = NULL, step = 1) ),
                        column(6, style = "padding-top:0px",
                               numericInput("axis4", "Set x max", value = NULL, step = 1) )
                      )
               ),
       
               # Section for plots of individual profiles: 
               column(4, style = "padding-top:22px",
                      plotOutput("plot_single", click = "click"),
                      column(8,
                             checkboxInput("helper_functions", "Automatic min/max recognition", value = TRUE, width = NULL)
                      ),
                      column(4,
                             downloadButton("downloadPlot_single", "Download plot", icon = icon("file-download"), width = '100%')
                      )
               ),
                 
               # Section for plots of aligned profiles:
               column(4, style = "padding-top:22px",
                      plotOutput("plot_align"),
                      column(8),
                      column(4,
                             downloadButton("downloadPlot_aligned", "Download plot", icon = icon("file-download"), width = '100%')
                      )
               ),
               
               column(2, style = "padding-right:22px",
                      # Section for axis settings of aligned polysome profiles:
                      tags$h4(tags$strong("Polysome profile")),
                               
                      # Set y-axis limits for polysome profiles:
                      fluidRow(
                        column(6,
                               numericInput("axis1_a", "Set y min", value = NULL, step = 1)),
                        column(6,
                               numericInput("axis2_a", "Set y max", value = NULL, step = 1))
                      ),
                      
                     
                      fluidRow(
                        # Set x-axis limits for polysome profiles:
                        column(6, 
                               numericInput("axis3_a", "Set x min", value = NULL, step = 1) ),
                        column(6,
                               numericInput("axis4_a", "Set x max", value = NULL, step = 1) )
                      ),
                      
                      # Section for axis settings of aligned fluorescence profiles:
                      tags$h4(tags$strong("Fluorescence profile")),
                      
                      # Show fluorescence signal or not?
                      fluidRow(
                        column(12,  
                               checkboxInput("show_fl_al", "Show fluorescence signal", value = FALSE, width = NULL) )
                      ),
                      
                      # Set y-axis limits for fluorescence profiles:
                      fluidRow(
                        column(6, 
                               numericInput("axis1_a_fl", "Set y min", value = NULL, step = 1) ),
                        column(6, 
                               numericInput("axis2_a_fl", "Set y max", value = NULL, step = 1) )
                      ),
                     
                      # Section for normalization options:
                      tags$h4(tags$strong("Normalization")),
                      tags$h5("Only available when area \"Total\" was quantified for all profiles"),
                      checkboxInput("normalize_height", HTML("Normalize <b>area</b> (y-values)"), value = FALSE, width = NULL),
                      checkboxInput("normalize_length", HTML("Normalize <b>length</b> (x-values)"), value = FALSE, width = NULL)
                )
             ),
             
             fluidRow(
               # Section for individual fluorescence profiles:
               column(2, style = "padding-left:20px",
                        column(8, 
                               fluidRow( tags$h4(tags$strong("Fluorescence profile")) ) 
                        ),
                      # Show fluorescence signal or not?
                        column(4, checkboxInput("show_fl", "Show", value = FALSE, width = '100%')
                        ),
                      
                        fluidRow( 
                          # Set a separate baseline for the fluorescence signal
                          column(6, style = "padding-top:23px", actionButton("baseline_fl", "Baseline", width = '100%') 
                          ),
                          # For smoothing of the fluorescence signal
                          column(6, style = "padding-top:0px",
                                 numericInput("slider1", "Smooth profile", value = 0, min = 0, max = 1, step = 0.01) 
                          )
                        ),
                        fluidRow(
                          column(6, style = "padding-top:0px",
                                 numericInput("axis1_fl", "Set y min", value = NULL, step = 1) ),
                          column(6, style = "padding-top:0px",
                                 numericInput("axis2_fl", "Set y max", value = NULL, step = 1) )
                        )
               ),
                     
               # Section for quantification of areas
               column(4, 
                      column(8, tags$h4(tags$strong("Quantification") ),
                             
                             fluidRow(style = "padding-top:30px",
                               column(6,
                                      actionButton("file_start", "Start", width = '100%')
                               ),
                               column(6,
                                      actionButton("file_end", "End", width = '100%')
                               )
                             ),
                             fluidRow(
                               column(6, style = "padding-top:15px",
                                      selectizeInput("select_area", "Select area", options = list(create=TRUE, plugins = list('restore_on_backspace')),
                                              choices = c("Total", "80S", "Polysomes", "40S", "60S"), width = '100%')
                               ),
                               column(6, style = "padding-top:35px",
                                      actionButton("quantify_area", "Quantify", width = '100%')
                               )
                             )
                      ),
                      column(4,
                             checkboxInput("red_lines", "Show x-anchor", value = TRUE, width = NULL),
                             checkboxInput("red_lines", "Show baseline", value = TRUE, width = NULL),
                             checkboxInput("green_lines", "Show quantified area", value = TRUE, width = NULL)
                      )
               ),
            
               # Section for alignment
               column(6,
                      fluidRow(
                        column(8,
                               fluidRow(
                                 column(4, tags$h4(tags$strong("Alignment") ) ),
                                 column(4),
                                 column(4, tags$h4(tags$strong("Colors and Lines") ) )
                               )
                        ),
                        column(4)
                      ),
                      
                      fluidRow(
                        column(8,
                          fluidRow(
                            column(4,
                                   selectInput("select_alignment", "Available profiles", choices = c(), width = '100%')
                            ),
                            column(1, style = "padding-top:22px",
                                   actionButton("up", label = NULL, icon = icon("angle-double-up"), width = '100%',
                                                style = "padding-left:5px") 
                            ),
                            column(1, style = "padding-top:22px",
                                   actionButton("down", label = NULL, icon = icon("angle-double-down"), width = '100%',
                                                style = "padding-left:5px") 
                            ),
                            column(2, style = "padding-top:22px",
                                   checkboxInput("show_in_al", "Display", value = TRUE)
                            ),
                            column(4,
                                   colourInput("color", "Color", palette = "square", value = "#FFFFFF")
                            )
                          )
                        ),
                      
                        column(4,
                               fluidRow(
                                   column(8,
                                          selectInput("linetype", "Line Type", choices = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                                                      selected = "solid")
                                   ),
                                   column(4,
                                          numericInput("linewidth", "Width", value = 1, min = 0) 
                                   )
                               )
                        )
                      ),
            
                      fluidRow(
                        column(8,
                          fluidRow(
                            column(4, style = "padding-top:10px",
                                   textInput("new_name", "Name in Legend", value = "", width = '100%', placeholder = "Name") 
                            ),
                            column(4, style = "padding-top:32px",
                                   actionButton("rename", "Rename", width = '100%')
                            ),
                            column(4,
                                     radioButtons("color_palette", "Select default color palette:",
                                                  c("Dark" = "dark_palette",
                                                    "Rainbow" = "rainbow_palette",
                                                    "Color blind friendly" = "color_blind"), 
                                                  selected = "dark_palette"),
                                   style = "padding-left:30px"
                            )
                          )
                        ),
                        
                        column(4,
                               checkboxInput("anchor_line", "Show x-anchor", value = TRUE, width = NULL)
                        )
                      ),
                   
                      fluidRow(
                        column(12,
                               style = "padding-left:0px",
                               tags$img(src='mi3_logo.png', align = "right", width = '10%', height = '10%')
                        )
                      )
               )
             )
    ),

    
    # SECOND TAB
    # Shows updating table of aligned (and normalized) profiles, 
    # which can be download 
    tabPanel(tags$strong("Alignment table"), icon = icon("table"), 
             tags$h4("Table of all aligned (and normalized) profiles"),
             downloadButton("downloadData", "Download .csv file"),
             tableOutput("csv_file")
    ),
    
    # THIRD TAB
    # Shows updating table of quantified areas for respective plots, 
    # which can be downloaded 
    tabPanel(tags$strong("Quantification summary"), icon = icon("list-alt"), 
             tags$h4("Table of quantified areas"),
             tableOutput("quantification"),
             downloadButton("downloadQuant", "Download .csv file")
    ),
    
    # FOURTH TAB 
    # Shows manual  
    tabPanel(tags$strong("QuAPPro manual"), icon = icon("question-circle"),
             fluidRow(
               column(2),
               column(8, htmltools::includeMarkdown("QuAPPro_v0-0-7_manual.Rmd")),
               column(2)
             ) 
    ),
    
    # FIFTH TAB
    # Release notes with link to GitHub repository:
    tabPanel(tags$strong("Release notes", style = "color:blue"), icon = icon("info", style = "color:blue"), 
             tags$h4("Release notes"),
             tags$div("This is version QuAPPro v0.0.7. 
                      To run it locally, you can download the code from our GitHub repository:",
                      tags$a(href="https://github.com/johannaschott/QuAPPro", "https://github.com/johannaschott/QuAPPro", style = "color:blue")
                      
             )
    ),
    
    # SIXTH TAB
    # Contact information, acknowledgments and impressum with links
    tabPanel(tags$strong("Contact", style = "color:blue"), icon = icon("address-card", style = "color:blue"), 
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
    )
)


#### SERVER ####

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 50 * 1024^2) # file size limit: 50 MB
  
  # Collection of reactive values: 
  val <- reactiveValues(
    colors_collected = list(), # profile-specific colors set by the user
    color_list = list(), # colors defined initially from a color palette
    linetype_collected = list(), # profile-specific line types set by the user
    linewidth_collected = list(), # profile-specific line widths set by the user
    file_starts = list(), # selected start points for area quantification
    file_ends = list(), # selected end points for area quantification
    sum_areas = list(), # sum of data points from start to end point of a quantification
    files_to_align = vector(), # profiles with baseline and x-anchor (which are automatically part of the alignment)
    files_to_plot = vector(), # profiles that should be shown in the alignment plot (individual profiles can be omitted by the user)
    show_in_alignment = list(), # logical value whether a profile in files_to_align should be plotted
    csv_file_df = data.frame(), # data frame with x- and y-values of all loaded profiles
    df_quant = data.frame(),# data frame with quantified areas
    factors_list = list(), # time-factor for each profile (how far are the data points apart?)
    control_baseline = list(), # baseline that was used for a quantification
    fl_control_baseline = list(), # baseline that was used for a quantification in a fluorescence profile
    file_types = vector(), # list of source file types
    xmin_collected = list(), # x-axis minima as defined by the user
    xmax_collected = list(), # x-axis maxima as defined by the user
    ymin_collected = list(), # y-axis minima as defined by the user
    ymax_collected = list(), # y-axis maxima as defined by the user
    ymin_collected_fl = list(), # y-axis minima of fluorescence profiles as defined by the user 
    ymax_collected_fl = list(), # y-axis maxima of fluorescence profiles as defined by the user 
    polysome_data = list(), # 
    xvalues = list(),
    files_list = list(),
    align_files_list = vector(),
    legend_names = vector()
  )
  
  
  
  ### INITIAL LOADED FILES
  
  # function for loading previous analysis from .RData file:
  import <- function(){
    load(input$input_data$datapath)
    
    for(i in names(forExport))
    {
      val[[ i ]] <- forExport[[i]]
    }
    
    updateSelectInput(session, "select",
                      choices = val$files_list, selected = val$files_list[[1]])
    
    updateSelectInput(session, "select2",
                      choices = val$files_list, selected = val$files_list[[1]])
    
    updateSelectInput(session, "select_alignment",
                      choices = val$align_files_list)
    
    updateColourInput(session, "color", value = val$color_list[[input$select_alignment]])
    updateNumericInput(session, "linewidth", value = val$linewidth_collected[[input$select_alignment]]  )
    updateSelectInput(session, "linetype", 
                      selected = c("solid", "dashed", "dotted",
                                   "dotdash", "longdash", "twodash")[ val$linetype_collected[[input$select_alignment]] ])
  }
      
 # function for reading polysome profile data:
  read_data <- function(){
    
    new_names <- input$input_data$name
    new_paths <- input$input_data$datapath
    
    for(i in 1:length(new_names))
    {
      this_name <- new_names[i]
      this_path <- new_paths[i]
      
      if( grepl(".pks",as.character(this_name) ) ){
        
        data <- read.table(this_path, dec = ",", header = F) 
        log(data[,3]) # returns an error when data is not numeric
        val$polysome_data[[this_name]] <- data[,3]
        
        val$factors_list[[this_name]] <- (0.1/60)
        val$file_types[this_name] <- "pks"
        val$xvalues[[this_name]] <- (data[ ,1]+1)*val$factors_list[[this_name]]
      }
      
      if( grepl(".csv",as.character(this_name) )  ){
        start <- grep("Data Columns:", readLines(this_path))
        data <- read.csv(this_path, skip = start)
        
        if( "SampleFluor" %in% colnames(data) )
        {
          log(data[,5]) # returns an error when data is not numeric
          log(data$SampleFluor) # returns an error when data is not numeric
          val$polysome_data[[this_name]] <- data[ ,5]
          val$fluo_data[[this_name]] <- data$SampleFluor
          
          val$factors_list[[this_name]] <- (0.32/60)
          val$file_types[this_name] <- "csv_fluo" 
          val$xvalues[[this_name]] <- (1:length(data[,5]))*val$factors_list[[this_name]]
        }
        
        if(!"SampleFluor" %in% colnames(data) )
        {
          log(data[,5]) # returns an error when the data is not numeric
          val$polysome_data[[this_name]] <- data[ ,5]
          
          val$factors_list[[this_name]] <- (0.2/60)
          val$file_types[this_name] <- "csv" 
          val$xvalues[[this_name]] <- (1:length(data[ ,5]))*val$factors_list[[this_name]]
        }
      }
      
      if( grepl(".txt",as.character(this_name) ) ){
        data <- read.delim(this_path, dec = ",")
        
        val$factors_list[[this_name]] <- (0.1/60)
        val$file_types[this_name] <- "txt"
        
        log(data[,2])
        log(data[,1])
        val$polysome_data[[this_name]] <- data[,2]
        val$xvalues[[this_name]] <-  val$factors_list[[this_name]]*(1:length(val$polysome_data[[this_name]]))
          # (data[ ,1]+1)*val$factors_list[[this_name]]
      }
    }
    
    val$paths_collected[input$input_data$name] <- input$input_data$datapath
    val$files_list <- c(input$input_data$name[input$input_data$size != 0], val$files_list)
    
    updateSelectInput(session, "select",
                      choices = val$files_list)
    
    updateSelectInput(session, "select2",
                      choices = val$files_list)
    
  }
  
  # Load data and show error message when loading failed:
  observeEvent(input$input_data, {
    #output$test <- renderText(paste("These files were loaded:", input$input_data$name) )
    
    if( length(input$input_data$name) > 1 & any(grepl("RData", input$input_data$name ) ) )
    {
      showNotification("Please select only one .RData file.",
                       duration = NULL, type = "error")
    }else{
    
      if( length(input$input_data$name) == 1 & any( grepl("RData", input$input_data$name ) ) )
      {
        loading_attempt <- try( import(), silent = T )
      }else{
        loading_attempt <- try( read_data(), silent = T )
      }
      if( any( class(loading_attempt) == "try-error") )
      {
        showNotification("File format not accepted.",
                         duration = NULL, type = "error")
      }
    }
  })
  
  
  fluorescence <- reactive({
    req(input$select) 
    req(input$show_fl)
    if( val$file_types[input$select] == "csv_fluo" ){
      smooth_profile( val$fluo_data[[input$select]], input$slider1 )}
    
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
    max_numeral <- max( floor(log10(abs( val$polysome_data[[input$select]] ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When user selects other axis limits, limit values change 
  ymin_single <- reactive({
    if(isTruthy(val$ymin_collected[[input$select]])){
      val$ymin_collected[[input$select]]*lost_num_pol()
    }else{
      0
    }})
  
  ymin_single_fl <- reactive({
    if(isTruthy(val$ymin_collected_fl[[input$select]])){
      val$ymin_collected_fl[[input$select]]*lost_num_fl()
    }else{
      0
    }})
  
  ymax_single <- reactive({
    if(isTruthy(val$ymax_collected[[input$select]])){
      val$ymax_collected[[input$select]]*lost_num_pol()
    }else{
      max(val$polysome_data[[input$select]])
    }})
  
  ymax_single_fl <- reactive({
    if(isTruthy(val$ymax_collected_fl[[input$select]])){
      val$ymax_collected_fl[[input$select]]*lost_num_fl()
    }else{
      max(fluorescence())
    }})
  
  xmin_single <- reactive({
    req(input$select)
    if(isTruthy(val$xmin_collected[[input$select]])){
      val$xmin_collected[[input$select]]
    }else{
      0
    }})
  
  xmax_single <- reactive({
    if(isTruthy(val$xmax_collected[[input$select]])){
      val$xmax_collected[[input$select]]
    }else{
      max(val$xvalues[[input$select]])
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
      if(input$helper_functions){
        val$file_starts[[input$select]] <- (find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                                (closest_point(input$click$x, 
                                                                               input$click$y, 
                                                                               val$xvalues[[input$select]],
                                                                               val$polysome_data[[input$select]])), 5))
      }else{
        val$file_starts[[input$select]] <-(closest_point(input$click$x,
                                                         input$click$y,
                                                         val$xvalues[[input$select]],
                                                         val$polysome_data[[input$select]]) )
      }
    # setting profile ends
    if(val$buttons == 5)
      if(input$helper_functions){
        val$file_ends[[input$select]] <- (find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             val$xvalues[[input$select]],
                                                                             val$polysome_data[[input$select]])), 5))
      }else{
        val$file_ends[[input$select]] <-(closest_point(input$click$x,
                                                       input$click$y,
                                                       val$xvalues[[input$select]],
                                                       val$polysome_data[[input$select]]))
      }
    # setting baseline for polysome profile
    if(val$buttons == 2){
      val$baseline[input$select] <- val$polysome_data[[input$select]][(input$click$x)/val$factors_list[[input$select]]]
    }
    
    # setting baseline for fluorescence profile
    if(val$buttons == 6){
      par(xpd = T)
      val$baseline_fl[input$select] <- fluorescence()[(input$click$x)/val$factors_list[[input$select]]]
    }
    
    # setting x-anchor
    if(val$buttons == 1)
      if(input$helper_functions){
        val$anchor[input$select] <- round(find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             val$xvalues[[input$select]],
                                                                             val$polysome_data[[input$select]])), 5))
      }else{
        val$anchor[input$select] <-(closest_point(input$click$x,
                                                  input$click$y,
                                                  val$xvalues[[input$select]],
                                                  val$polysome_data[[input$select]]))
      }
    
    # if both, anchor and baseline, exist for one plot, select_alignment button gets updated with this file
    val$align_files_list <- c((intersect(names(val$baseline),
                                                      names(val$anchor))),
                                           val$align_files_list)
    updateSelectInput(session, "select_alignment",
                      choices = val$align_files_list)
    # reactive vector with respective file names gets updated
    val$files_to_align <- intersect(names(val$baseline), names(val$anchor))
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
  # Calculate sum of yvalues to quantify areas when button is pressed and required values were set before
  observeEvent(input$quantify_area,{
    req(val$file_starts[[input$select]], val$file_ends[[input$select]], val$baseline[input$select])
    
    # store file baselines for quantified areas to be displayed in the plotting area even when baseline is changed again by the user
    if(isTruthy(val$baseline_fl[input$select])){
      val$fl_control_baseline[[input$select_area]][input$select] <- val$baseline_fl[input$select]
    }
    
    val$control_baseline[[input$select_area]][input$select] <- val$baseline[input$select]
    
    x_first <- round(val$file_starts[[input$select]])*val$factors_list[[input$select]]
    x_last <- round(val$file_ends[[input$select]])*val$factors_list[[input$select]]
    # sum up all values from start to end area and subtract baseline, assign area name selected by user
    # store in list of areas including quantification with respective file name
    val$sum_areas[[input$select_area]][input$select] <- (sum(val$polysome_data[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))]-val$baseline[input$select]))
    
    if(input$show_fl){
      val$sum_areas[[paste(input$select_area, "_fluo", sep = "")]][input$select] <- (sum(fluorescence()[(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))]-val$baseline_fl[input$select]))
    }
    
    # store not only area but also start and stop values in the same ways as areas!
    val$area_starts[[input$select_area]][input$select] <- val$file_starts[[input$select]]
    val$area_ends[[input$select_area]][input$select] <- val$file_ends[[input$select]]
  })
  #select_area
  
  
  # create data frame with quantified areas containing NAs for files without respective area
  observeEvent(val$sum_areas, {
    areas_with_quant <- names(val$sum_areas)
    
    a <- areas_with_quant[1]
    df_quant <- data.frame(File = names(val$sum_areas[[ a ]]), val$sum_areas[[ a ]])
    colnames(df_quant)[2] <- a
    for(a in areas_with_quant[-1])
    {
      df_quant_new <- data.frame(File = names(val$sum_areas[[ a ]]), val$sum_areas[[ a ]])
      colnames(df_quant_new)[2] <- a
      df_quant <- merge(df_quant, df_quant_new, by = "File", all = T)
    }
    val$df_quant <- df_quant
  })
  
  # save axis limits set by the user for single profiles
  observeEvent(input$axis3, {
    val$xmin_collected[[input$select]] <- input$axis3
  })
  observeEvent(input$axis4, {
    val$xmax_collected[[input$select]] <- input$axis4
  })
  observeEvent(input$axis1, {
    val$ymin_collected[[input$select]] <- input$axis1
  })
  observeEvent(input$axis2, {
    val$ymax_collected[[input$select]] <- input$axis2
  })
  observeEvent(input$axis1_fl, {
    if(is.numeric(input$axis1_fl)){
      val$ymin_collected_fl[[input$select]] <- input$axis1_fl 
    }
  })
  observeEvent(input$axis2_fl, {
    if(is.numeric(input$axis2_fl)){
      val$ymax_collected_fl[[input$select]] <- input$axis2_fl
    }
  })
  
  # When a new profile is available for alignment, automatically show in alignment
  observeEvent(val$files_to_align, {
    val$show_in_alignment[[ rev( val$files_to_align)[1] ]] <- TRUE
    val$files_to_plot <- val$files_to_align[unlist(val$show_in_alignment[val$files_to_align])] # only the names of those that are true
  })
  
  # When a new profile is available for alignment, set the linetype to 1
  observeEvent(val$files_to_align, {
    val$linetype_collected[[rev( val$files_to_align)[1] ]] <- 1
  })
  
  # When a new profile is available for alignment, set the line width to 1
  observeEvent(val$files_to_align, {
    val$linewidth_collected[[rev( val$files_to_align)[1] ]] <- 2
  })
  
  # Let user change linetypes of the different plots (selected out of aligned files)
  observeEvent(input$linetype, {
    req(input$select_alignment)
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
    req(input$select_alignment)
    val$linewidth_collected[[input$select_alignment]] <- input$linewidth
  })
  
  # Let user change colors of the different plots (selected out of aligned files)
  observeEvent(input$color, {
    req(input$select_alignment)
    val$colors_collected[[input$select_alignment]] <- input$color
  })
  
  
  # update color, line type and width to the values of the selected profiles
  # Caution: updateColourInput seems to require double square brackets [[]]!
  observe({
    req(input$select_alignment)
    updateColourInput(session, "color", value = val$color_list[[input$select_alignment]])
    updateNumericInput(session, "linewidth", value = val$linewidth_collected[[input$select_alignment]]  )
    updateSelectInput(session, "linetype", 
                      selected = c("solid", "dashed", "dotted",
                                   "dotdash", "longdash", "twodash")[ val$linetype_collected[[input$select_alignment]] ])
  })
  
  
  ### OUTPUT
  
  # show_fl is deselected when a new file is selected that does not contain fluorescence
  # no warning or error is shown
  observeEvent(input$select,{
    if(input$show_fl && val$file_types[input$select] != "csv_fluo"){
      updateCheckboxInput(session, "show_fl", value = FALSE)
    }
  })
  
  #observe({
  #  updateNumericInput(session, "axis3", value = round(xmin_single(), digits = 2 ) )
  #  updateNumericInput(session, "axis4", value = round(xmax_single(), digits = 2 ) )
  #  updateNumericInput(session, "axis1", value = round(ymin_single()/lost_num_pol(), digits = 2) )
  #  updateNumericInput(session, "axis2", value = round(ymax_single()/lost_num_pol(), digits = 2) )
  #})
  
  # show notification if show fluorescence is selected but there is no fluorescence signal in the currently selected file
  observeEvent(input$show_fl,{
    req(input$select)
    if( val$file_types[input$select] != "csv_fluo" & input$show_fl){
      showNotification("Your file does not contain a fluorescence signal (column 'SampleFluor').",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl", value = FALSE)
    }
  })
  
  # let axis limits update for fluorescence axis
  #observe({
  #  req(input$select)
  #  if( val$file_types[input$select] == "csv_fluo" & input$show_fl){
  #    updateNumericInput(session, "axis1_fl", value = round( ymin_single_fl()/lost_num_fl(), digits = 2 ) )
  #    updateNumericInput(session, "axis2_fl", value = round( ymax_single_fl()/lost_num_fl(), digits = 2 ) )
  #  }else{
  #    updateNumericInput(session, "axis1_fl", value = "")
  #    updateNumericInput(session, "axis2_fl", value = "")
  #  } 
  #})
  
  # show notification if "Show fluorescence" was selected for alignment 
  # but there is not fluorescence signal
  # or there was not baseline set for at least on fluorescence signal
  observeEvent(input$show_fl_al,{
    if( !any( val$file_types[val$files_to_plot] == "csv_fluo") & input$show_fl_al){
      showNotification("None of your aligned profiles contains a fluorescence signal (column 'SampleFluor').",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( !all( val$file_types[val$files_to_plot] == "csv_fluo") & input$show_fl_al  ){
        showNotification("Some of your aligned profiles do not contain a fluorescence signal (column 'SampleFluor').",
                         duration = NULL, type = "warning")
      }
    }
    
    if(any( val$file_types[val$files_to_plot] == "csv_fluo") & is.null(val$baseline_fl) & input$show_fl_al){
      showNotification("You did not set a baseline for your fluorescence profiles.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( any( val$file_types[val$files_to_plot] == "csv_fluo") & !all( names(val$file_types[val$file_types[val$files_to_plot] == "csv_fluo"]) %in% names(val$baseline_fl) ) & input$show_fl_al  ){
        showNotification("Some of your fluorescence profiles do not have a baseline.",
                         duration = NULL, type = "warning")
      }
    }
  })
  
  # show notification if how fluorescence was selected already but a new file in val$files_to_plot
  # does not contain a fluorescence signal
  observeEvent(val$files_to_plot, {
    if( !all( names(val$file_types[val$files_to_plot] == "csv_fluo") %in% names(val$baseline_fl) ) & input$show_fl_al  ){
      showNotification("Some of your fluorescence profiles do not have a baseline.",
                       duration = NULL, type = "warning")
    }
  })
  
  
  
  # show notification when normalization to length or height is set but total area is missing
  observe({
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( val$files_to_plot %in% files_with_total ) & (input$normalize_length | input$normalize_height) ){
      showNotification("Please select a total area for all profiles in the alignment.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "normalize_length", value = FALSE)
      updateCheckboxInput(session, "normalize_height", value = FALSE)
    }
  })
  
  ## plot individual profiles
  plot_singleFl <-function(cex_lab, cex_axis, lwd, mgp2){
    if(lost_num_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_fl(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select]], fluorescence(), type = "l", xaxt = "n", col = "black", lwd = lwd,
         xlim =c(xmin_single(),xmax_single()), ylim = c(ymin_single_fl(), ymax_single_fl()),
         las = 1, ylab = ylab, xlab = "", yaxt = "n", cex.lab = cex_lab, mgp = mgp2)
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_fl(), las = 1, mgp = mgp2, cex.axis = cex_axis)
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      # quantified area is colored in green. If baseline, area_end or area_satrt lines are changed, the colored area stays the same for the selected quantified area 
      # until the button "quantify area" is pressed again.
      if(isTruthy(val$fl_control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, val$xvalues[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], x_last),
                c(val$fl_control_baseline[[input$select_area]][input$select], fluorescence()[(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], val$fl_control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "darkgreen", lwd = 2)
      }
    }
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = lwd)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = lwd)
    }
    if(input$red_lines){
      abline(h = val$baseline_fl[input$select],col = "red", lty=2, lwd = lwd)
    }
  }
  
  plot_singlePol <-function(cex_lab, cex_axis, lwd, mgp1, mgp2){
    req(val$xvalues[[input$select]])
    req(val$polysome_data[[input$select]])
    if(lost_num_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_pol(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select]], val$polysome_data[[input$select]], type = "l", las = 1, lwd = lwd,
         ylab = ylab, xlab = "", mgp = mgp2,
         ylim = c(ymin_single(),ymax_single()), 
         xlim =c(xmin_single(),xmax_single()),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab
    )
    
    
    axis(1, las = 1, mgp = mgp1, cex.axis = cex_axis)
    mtext("Time (min)", side = 1, line = mgp1[1], cex = cex_lab)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = mgp2, cex.axis = cex_axis)
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      if(isTruthy(val$control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, val$xvalues[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], x_last),
                c(val$control_baseline[[input$select_area]][input$select], val$polysome_data[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], val$control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "black", lwd = lwd)
      }
    }
    # selected x-anchor and baseline are displayed if box is ticked
    if(input$red_lines){
      abline(v = val$anchor[input$select]*val$factors_list[[input$select]],col = "red", lty=2, lwd = lwd)
      abline(h = val$baseline[input$select],col = "red", lty=2, lwd = lwd)
    }
    # selected area starts and ends are displayed if box is ticked
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = lwd)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = lwd)
    }
  }
  
  
  
  plot_singleInput <- function(cex_lab, cex_axis, lwd, mgp1, mgp2, mar_factor){
    if(!is.null(val$xvalues[[input$select]])){
      if(input$show_fl & val$file_types[input$select] == "csv_fluo" ){
        if(val$buttons == 6){
          layout(matrix(2:1, 2, 1), height = c(0.6, 0.9) ) # divides the plotting area into 2 rows
          par(mar = c(5, 5, 0, 2)*mar_factor)
          plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2)
          par(mar = c(0, 5, 0.5, 2)*mar_factor)
          plot_singleFl(cex_lab, cex_axis, lwd, mgp2)
        }else{
          layout(matrix(1:2, 2, 1), height = c(0.6, 0.9) ) # divides the plotting area into 2 rows
          par(mar = c(0, 5, 0.5, 2)*mar_factor)
          plot_singleFl(cex_lab, cex_axis, lwd, mgp2)
          par(mar = c(5, 5, 0, 2)*mar_factor)
          plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2)
        }
      }else{
        par(mar = c(5, 5, 0.5, 2)*mar_factor)
        plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2)
      }
    }
  }
  
  output$plot_single <- renderPlot({
    req(val$xvalues[[input$select]])
    plot_singleInput(cex_lab = 1.6, cex_axis = 1.6, lwd = 2, 
                     mgp1 = c(3.5, 1.2, 0), mgp2 = c(3.5, 0.8, 0), 
                     mar_factor = 1)
  })
  
  # Enable download of current plot as pdf
  
  output$downloadPlot_single <- downloadHandler(
    filename = function(){ 
      paste( 
        sub(pattern = "(.*?)\\..*$", replacement = "\\1", as.character(input$select) ),
        ".pdf", sep = ""
      )
    },
    content = function(file) {
      pdf(file, width = 2, height = 1.2, pointsize = 6 )
      print( plot_singleInput(cex_lab = 1, cex_axis = 0.8, lwd = 0.8, 
                              mgp1 = c(1.7, 0.5, 0), mgp2 = c(2, 0.8, 0), 
                              mar_factor = 0.6) )
      dev.off()
    })  
  
  # create output table showing quantification data with option to download table with given name
  output$quantification <- renderTable(
    val$df_quant
  )
  
  output$downloadQuant <- downloadHandler(
    filename = "AreaQuant.csv",
    content = function(file) {
      write.csv2( val$df_quant, file, row.names = F)
    }
  )

  # Create alignment plot
  observeEvent(input$show_in_al, {
    if(input$show_in_al)
    {
      val$show_in_alignment[[input$select_alignment]] <- TRUE
    }else{
      val$show_in_alignment[[input$select_alignment]] <- FALSE
    }
    val$files_to_plot <- val$files_to_align[unlist( val$show_in_alignment[val$files_to_align]) ] # only the names of those that are true
  })
  
  observeEvent(input$select_alignment, {
    updateCheckboxInput(session, "show_in_al", value = val$show_in_alignment[[input$select_alignment]])
  })
  
  
  # move files up or down in the val$files_to_align vector
  observeEvent(input$up, {
    if(input$select_alignment %in% val$files_to_plot)
    {
      to_be_shifted <- which( val$files_to_plot == input$select_alignment )
      if(to_be_shifted > 1)
      {
        files_order <- 1:length(val$files_to_plot)
        files_order[to_be_shifted] <- files_order[to_be_shifted] - 1
        files_order[to_be_shifted - 1] <- files_order[to_be_shifted - 1] + 1
        val$files_to_plot <- val$files_to_plot[files_order]
      }
    }
  })
  
  # move files up or down in the val$files_to_align vector
  observeEvent(input$down, {
    if(input$select_alignment %in% val$files_to_plot)
    {
      to_be_shifted <- which( val$files_to_plot == input$select_alignment )
      if(to_be_shifted < length(val$files_to_plot) )
      {
        files_order <- 1:length(val$files_to_plot)
        files_order[to_be_shifted] <- files_order[to_be_shifted] + 1
        files_order[to_be_shifted + 1] <- files_order[to_be_shifted + 1] - 1
        val$files_to_plot <- val$files_to_plot[files_order]
      }
    }
  })
  
  # if color palette choice gets changed by the user, the originally selected colors are reset to an empty list:
  observeEvent(
    input$color_palette,              
    {val$colors_collected = list()}, 
    ignoreInit = TRUE
  )
  
  # create color values, rainbow palette as initial colors, further replaced by selected color if selected
  observe({
    if(input$color_palette == "dark_palette"){
      dummy <- qualitative_hcl(length(val$files_to_plot), palette = "Dark 3") 
    }
    if(input$color_palette == "rainbow_palette"){
      dummy <- rainbow(length(val$files_to_plot)) 
    }
    if(input$color_palette == "color_blind"){
      pal <- c("#000000", "#ff6db6", "#006ddb", "#920000",
               "#b66dff", "#004949", "#ffb6db", "#6db6ff", "#924900",
               "#24ff24", "#929200", "#490092", "#b6dbff", "#db6d00", "#ffff6d")
      dummy <- rep(pal, ceiling( length( val$files_to_plot )/length(pal ) ) )
    }
    names(dummy) <- sort(val$files_to_plot)
    if (length(val$colors_collected) > 0){
      dummy[names(val$colors_collected) ] <- unlist(val$colors_collected)
    }
    val$color_list <- as.list(dummy)
  })
  
  # When a new file is added to the alignment, the file name is used for the figure legend:
  observeEvent(val$files_to_align, {
    newest <- tail(val$files_to_align, 1)
    basename <- gsub(".csv|.txt|.pks", "", newest)
    val$legend_names[newest] <- basename
  })
  
  # Change names of profiles for legend:
  observeEvent(input$rename, {
    val$legend_names[input$select_alignment] <- input$new_name 
  })
  
  # for normalization of surfaces:
  norm_factor <- reactive({
    dummy <- vector()
    files_with_total <- names( val$sum_areas[["Total"]])
    for(f in val$files_to_plot){
      if((f %in% files_with_total) & input$normalize_height){
        # normalization y-values
        dummy[f] <- val$sum_areas[["Total"]][f]/val$sum_areas[["Total"]][files_with_total[1]]
      }else{
        dummy[f] <- 1
      }
    }
    dummy
  })
  
  # for normalization of length:
  norm_factor_x <- reactive({
    dummy <- vector()
    files_with_total <- names( val$sum_areas[["Total"]])
    for(f in val$files_to_plot){
      if( (f %in% files_with_total) & input$normalize_length){
        dummy[f] <- (round(val$area_ends[["Total"]][f] - val$area_starts[["Total"]][f])/(round(val$area_ends[["Total"]][files_with_total[1]] - val$area_starts[["Total"]][files_with_total[1]])))
      }else{
        dummy[f] <- 1
      }
    }
    dummy
  })
  
  # determine normalized y-values of fluorescence, only for files with baseline
  values_fluorescence <- reactive({
    if(any(val$file_types[val$files_to_plot] == "csv_fluo") )
    {
      dummy <- list()
      for(f in names(val$baseline_fl)) # only go through fluorescence signals with baseline
      {
        fluo <- val$fluo_data[[f]]
        dummy[[f]] <- smooth_profile( ( ( fluo - val$baseline_fl[[f]] )/norm_factor()[f])*norm_factor_x()[f], input$slider1)
      }
      dummy
    }else{  
      NULL
    }
  })
  
  
  values_list <- reactive({
    dummy <- list()
    for(f in val$files_to_plot)
    {
      next_y <- val$polysome_data[[f]]
      dummy[[f]] <- ( ( next_y - val$baseline[f])/norm_factor()[f])*norm_factor_x()[f]
    }
    dummy
  })
  
  
  # determine shift along x-axis for all:
  shifts <- reactive({
    # find max anchor of all files to be aligned
    all_anchors <- (val$anchor[val$files_to_plot])/norm_factor_x()[val$files_to_plot]
    super_anchor <- max(all_anchors)
    dummy <- super_anchor - all_anchors
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # find starts and ends of aligned profiles
  aligned_starts <- reactive({
    dummy <- ((1/norm_factor_x()[val$files_to_plot] + (shifts())))
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  aligned_ends <- reactive({
    # store lengths of the profiles 
    profile_lengths <- sapply(values_list(), FUN = length)
    names(profile_lengths) <- val$files_to_plot
    dummy <- ((profile_lengths/norm_factor_x()[val$files_to_plot] + (shifts())))
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  xmax_all <- reactive({
    max(aligned_ends())
  })
  
  # find maximum of y-axis values for polysome profiles
  ymax_all <- reactive({
    profile_heights <- sapply(values_list()[val$files_to_plot], FUN = max)
    max(profile_heights )
  })
  
  # find maximum of y-axis values for fluorescence profiles
  ymax_all_fl <- reactive({
    profile_heights <- sapply(values_fluorescence()[val$files_to_plot], FUN = max)
    max(profile_heights )
  })
  
  xmin_all <- reactive({
    min(aligned_starts())
  })
  
  # find maximum of y-axis values for polysome profiles
  ymin_all <- reactive({
    profile_mins <- sapply(values_list()[val$files_to_plot], FUN = min)
    min(profile_mins )
  })
  
  # find maximum of y-axis values for fluorescence profiles
  ymin_all_fl <- reactive({
    profile_mins <- sapply(values_fluorescence()[val$files_to_plot], FUN = min)
    min(profile_mins )
  })
  
  # axis labels are modified if there are more than three numerals
  lost_num_al_fl <- reactive({
    max_numeral <- max( floor(log10(abs( unlist(values_fluorescence()[intersect( val$files_to_plot, names(val$baseline_fl) ) ] ) ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_al_pol <- reactive({
    max_numeral <- max( floor(log10(abs( unlist(values_list()[val$files_to_plot]) ))) + 1 ) # how many numerals has the maximum UV absorption?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_al_Index <- reactive({
    max_numeral <- max( floor(log10(abs( aligned_ends() ))) + 1 ) # how many numerals has the maximum Index?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When user selects other axis limits, limit values change 
  observeEvent(input$axis1_a, {
    val$ymin <- input$axis1_a*lost_num_al_pol()
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis2_a, {
    val$ymax <- input$axis2_a*lost_num_al_pol()
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis1_a_fl, {
    val$ymin_fl <- input$axis1_a_fl*lost_num_al_fl()
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis2_a_fl, {
    val$ymax_fl <- input$axis2_a_fl*lost_num_al_fl()
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis3_a, {
    val$xmin <- input$axis3_a*lost_num_al_Index()
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis4_a, {
    val$xmax <- input$axis4_a*lost_num_al_Index()
  }, 
  ignoreInit = T)
  
  ymin_aligned <- reactive({
    if(isTruthy(val$ymin) ){
      val$ymin
    }else{
      ymin_all()
    }})
  
  ymin_aligned_fl <- reactive({
    if(isTruthy(val$ymin_fl)){
      val$ymin_fl
    }else{
      ymin_all_fl()
    }})
  
  ymax_aligned <- reactive({
    if(isTruthy(val$ymax)){
      val$ymax
    }else{
      ymax_all()
    }})
  
  ymax_aligned_fl <- reactive({
    if(isTruthy(val$ymax_fl)){
      val$ymax_fl
    }else{
      ymax_all_fl()
    }})
  
  xmin_aligned <- reactive({
    req(input$select)
    if(isTruthy(val$xmin)){
      val$xmin
    }else{
      xmin_all()
    }})
  
  xmax_aligned <- reactive({
    if(isTruthy(val$xmax)){
      val$xmax
    }else{
      xmax_all()
    }})
  
  
  #observe({
  #  if(length(val$files_to_plot >= 1))
  #  {
  #    updateNumericInput(session, "axis3_a", value = round( xmin_aligned() /lost_num_al_Index(), digits = 2 ) )
  #    updateNumericInput(session, "axis4_a", value = round( xmax_aligned() /lost_num_al_Index(), digits = 2 ) )
  #    updateNumericInput(session, "axis1_a", value = round( ymin_aligned()/lost_num_al_pol(), digits = 2 ) )
  #    updateNumericInput(session, "axis2_a", value = round( ymax_aligned()/lost_num_al_pol(), digits = 2 ) )
  #  }
  #})
  
  # let axis limits update for fluorescence axis
  observe({
    if(length(val$files_to_plot >= 1))
    {
      req(val$baseline_fl)
      req(input$show_fl_al)
      updateNumericInput(session, "axis1_a_fl", value = round( ymin_aligned_fl() /lost_num_al_fl(), digits = 2 ) )
      updateNumericInput(session, "axis2_a_fl", value = round( ymax_aligned_fl() /lost_num_al_fl(), digits = 2 ) )
    }
  })
  
  
  #### plotting function for aligned profiles
  
  plot_alignedPol <- function(cex_lab, cex_axis, lwd, 
                              mgp1, mgp2, 
                              cex_legend, lwd_factor){
    if(lost_num_al_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_al_pol(), ")", sep = "")
    }
    
    f <- val$files_to_plot[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    plot(x, values_list()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]]*lwd_factor, 
         ylab = ylab, xlab = "", las = 1,
         col = val$color_list[[f]], mgp = mgp2, 
         ylim = c(ymin_aligned(),ymax_aligned()), xlim = c(xmin_aligned(),xmax_aligned()),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab
    )
    
    
    # store values in df for creating alignment table
    y_aligned <- values_list()[[f]]
    x_aligned <- x
    df = list(x_aligned = x_aligned, y_aligned = y_aligned)
    attributes(df) = list(names = names(df),
                          row.names=1:max(length(x_aligned), length(y_aligned)), class='data.frame')
    colnames(df) <- c("Index", as.character(str_remove(f, ".pks|.csv|.txt")))
    csv_file_df <- df
    
    ### Shows anchor (when "Display x-anchor in alignment" is selected)
    if(input$anchor_line == TRUE){
      anchor_line <- val$anchor[val$files_to_plot[1]] / norm_factor_x()[val$files_to_plot[1]] + shifts()[val$files_to_plot[1]]
      abline(v = anchor_line, col = "red", lty=2, lwd = lwd)
    }
    
    a <- axTicks(1)
    axis(1, at = a, labels = a/lost_num_al_Index(), las = 1, mgp = mgp1, cex.axis = cex_axis)
    mtext("Relative position", side = 1, line = mgp1[1], cex = cex_lab)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_pol(), las = 1, mgp = mgp2, cex.axis = cex_axis)
    
    for(f in val$files_to_plot[-1])
    {
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      
      # create dataframe with x and y values to merge together in each for loop round. Col names are data names.
      y_aligned <- values_list()[[f]]
      x_aligned <- x
      df2 = list(x_aligned = x_aligned, y_aligned=y_aligned)
      attributes(df2) = list(names = names(df2),
                             row.names=1:max(length(x_aligned), length(y_aligned)), class='data.frame')
      colnames(df2) <- c("Index", as.character(str_remove(f, ".pks|.csv|.txt"))) 
      csv_file_df <- merge(csv_file_df, df2, by="Index", all = T)
      
      # plot files in alignment
      points(x, values_list()[[f]], type = "l", lty = val$linetype_collected[[f]], 
             lwd = val$linewidth_collected[[f]]*lwd_factor, col = val$color_list[[f]])
    }
    
    #create reactive dataframe of all plots to have access outside of renderPlot function
    val$csv_file_df <- csv_file_df
  }
  
  plot_alignedFluo <- function(cex_lab, cex_axis, lwd, 
                               mgp2, lwd_factor){
    
    if(lost_num_al_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_al_fl(), ")", sep = "")
    }
    
    files_in_al_fluo <- val$files_to_plot[val$files_to_plot %in% intersect(names(val$baseline_fl),  names(val$baseline))] # this is necessary to maintain the order as given by val$files_to_plot
    f <-  files_in_al_fluo[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    plot(x, values_fluorescence()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]]*lwd_factor,
         ylab = ylab, xlab = "", las = 1,
         col = val$color_list[[f]], mgp = mgp2, 
         ylim = c(val$ymin_fl,val$ymax_fl), xlim = c(val$xmin,val$xmax),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab
    )
    # store values in df for creating fluo alignment table
    y_aligned <- values_fluorescence()[[f]]
    x_aligned <- x
    df = list(x_aligned = x_aligned, y_aligned = y_aligned)
    attributes(df) = list(names = names(df),
                          row.names=1:max(length(x_aligned), length(y_aligned)), class='data.frame')
    colnames(df) <- c("Index", paste0(as.character(str_remove(f, ".pks|.csv|.txt")),"_fluo"))
    csv_file_df_fluo <- df
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_fl(), las = 1, mgp = mgp2, cex.axis = cex_axis)
    
    for(f in  files_in_al_fluo[-1])
    {
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      
      # create dataframe with x and y values to merge together in each for loop round. Col names are data names plus "_fluo".
      y_aligned <- values_fluorescence()[[f]]
      x_aligned <- x
      df2 = list(x_aligned = x_aligned, y_aligned=y_aligned)
      attributes(df2) = list(names = names(df2),
                             row.names=1:max(length(x_aligned), length(y_aligned)), class='data.frame')
      
      colnames(df2) <- c("Index", paste0(as.character(str_remove(f, ".pks|.csv|.txt")),"_fluo"))
      csv_file_df_fluo <- merge(csv_file_df_fluo, df2, by="Index", all = T)
      
      points(x, values_fluorescence()[[f]], type = "l", lty = val$linetype_collected[[f]], 
             lwd = val$linewidth_collected[[f]]*lwd_factor, col = val$color_list[[f]])
    }
    
    #create reactive dataframe of all plots to have access outside of renderPlot function
    val$csv_file_df_all <- merge(val$csv_file_df, csv_file_df_fluo, by="Index", all = T)
  }
  
  legend_alignment <- function(cex_legend, lwd_factor){
    par(mar = c(0, 0, 0, 0) )
    plot(0, 0, yaxt = "n", xaxt = "n", xlab = "", ylab = "", type = "n", bty = "n")
    legend_names <- val$legend_names[val$files_to_plot]
    legend_adjustment <- sqrt( min( c(1, 8/max( nchar(legend_names) ) ) ) )
    legend("topright", legend = val$legend_names[val$files_to_plot], lty = unlist(val$linetype_collected[val$files_to_plot]), 
           lwd = unlist(val$linewidth_collected[val$files_to_plot])*lwd_factor, col = unlist(val$color_list[val$files_to_plot]),
           bty = "n", cex = cex_legend*legend_adjustment)
  }
  
  plot_alignment <- function(cex_lab, cex_axis, lwd, 
                             mgp1, mgp2, 
                             cex_legend, lwd_factor,
                             mar_factor){
    if(input$show_fl_al){
      layout(matrix(1:4, 2, 2), height = c(0.6, 0.9), width = c(0.8, 0.2) ) # divides the plotting area into 2 rows
      par(mar = c(0, 5, 0.5, 2)*mar_factor)
      plot_alignedFluo(cex_lab, cex_axis, lwd, 
                       mgp2, lwd_factor)
      par(mar = c(5, 5, 0, 2)*mar_factor) 
      plot_alignedPol(cex_lab, cex_axis, lwd, 
                      mgp1, mgp2, 
                      cex_legend, lwd_factor)
      legend_alignment(cex_legend, lwd_factor)
    }else{
      layout(matrix(1:2, 1, 2), height =  c(0.6, 0.9), width = c(0.8, 0.2) ) # divides the plotting area into 2 columns
      par(mar = c(5, 5, 0.5, 2)*mar_factor)
      plot_alignedPol(cex_lab, cex_axis, lwd, 
                      mgp1, mgp2, 
                      cex_legend, lwd_factor)
      legend_alignment(cex_legend, lwd_factor)
    }
  }
  
  output$plot_align <- renderPlot({
    req(val$files_to_plot)
    if(length(val$files_to_plot) >= 1)
    {
      plot_alignment(cex_lab = 1.6, cex_axis = 1.6, lwd = 2, 
                     mgp1 = c(3.5, 1.2, 0), mgp2 = c(3.5, 0.8, 0), 
                     cex_legend = 1.5, lwd_factor = 1,
                     mar_factor = 1)
    }
  })
  
  # Enable download of current plot as pdf
  output$downloadPlot_aligned <- downloadHandler(
    filename = "alignment.pdf",
    content = function(file) {
      pdf(file,width = 2, height = 1.2, pointsize = 6 )
      print( plot_alignment(cex_lab = 1, cex_axis = 0.8, lwd = 0.8, 
                            mgp1 = c(1.7, 0.5, 0), mgp2 = c(2, 0.8, 0),  
                            cex_legend = 0.5, lwd_factor = 0.5,
                            mar_factor = 0.6) )
      dev.off()
    })  
  
  # create table output of aligned files with option to download table with given name
  output$csv_file <- renderTable(
    if(input$show_fl_al){
      val$csv_file_df_all
    }else{
      val$csv_file_df
    }
  )
  output$downloadData <- downloadHandler(
    filename = "alignment.csv",
    content = function(file) {
      if(input$show_fl_al){
        write.csv2(val$csv_file_df_all, file, row.names = FALSE)
      }else{
        write.csv2(val$csv_file_df, file, row.names = FALSE)
      }
    }
  )
  
  # Export entire analysis:
  output$export <- downloadHandler(
    filename <- function(){
      paste("QuAPPro.RData")
    },
    content = function(file) {
      forExport <- list()
      for(i in names(val) )
      {
        forExport[[ i ]] <- val[[i]]
      }
      save(forExport, file = file)
    }
  )
}

###############
shinyApp(ui = ui, server = server)