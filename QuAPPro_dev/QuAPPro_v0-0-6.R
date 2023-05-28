#####################################################
#     ____                 _____  _____             #
#    / __ \          /\   |  __ \|  __ \            #
#   | |  | |_   _   /  \  | |__) | |__) | __ ___    #
#   | |  | | | | | / /\ \ |  ___/|  ___/ '__/ _ \   #
#   | |__| | |_| |/ ____ \| |    | |   | | | (_) |  #
#    \___\_\\__,_/_/    \_\_|    |_|   |_|  \___/   #
#                                                   #                                                
#####################################################                                                 

# LOAD REQUIRED PACKAGES:

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


# FUNCTIONS:

render_dt = function(data, editable = 'cell', server = TRUE, ...) {
  renderDT(data, selection = 'none', server = server, editable = editable, ...)
}

# Find closest point in a two-dimensional coordinate system:

# This function 

closest_point <- function(x, y, X, Y){
  all_dist <- sqrt( (x-X)^2 + (y-Y)^2)
  closest <- order(all_dist, decreasing = F)[1]
  return(X[closest])
}

# Find closest minimum or maximum:

# This function identifies the closest local minimum or maximum 
# in a vector of values to a selected position in the vector.

# y: vector of values
# a: index of selected point
# range: Minimal number of consecutive values that show increase or decrease

find_closest_minmax <- function(y, a, range)
{
  d <- diff(y)
  d[ d > 0 ] <- 1
  d[ d < 0 ] <- -1
  seg <- rle(d)
  seg_ends <- cumsum(seg$lengths) # last position of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)] + 1) # first position of a segment
  up_down <- which( seg$lengths > range & seg$values != 0 )
  up_down_change <- abs(diff( seg$values[up_down]) ) == 2
  change_seg <- up_down[which(up_down_change) + 1]
  before_change_seg <-  up_down[which(up_down_change)]
  minmax <- (seg_ends[before_change_seg] + 1 + seg_starts[change_seg]) / 2
  closest <- order( abs(minmax -  a), decreasing = F)[1]
  return(minmax[closest])
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

smooth_profile <- function(y, spar)
{
  if(is.null(spar) )
  {
    return( smooth.spline(1:length(y), y, spar = spar)$y )
  }else{
    if(spar == 0)
    {
      return(y)
    }else{
      return( smooth.spline(1:length(y), y, spar = spar)$y )
    }
    
  }
}

second_deriv <- function(y)
{
  return( diff(y, differences = 2) )
}

all_min <- function(y, range)
{
  d <- diff(y )
  d0 <- d
  d0[ d > 0 ] <- 1
  d0[ d < 0 ] <- -1
  seg <- rle(d0)
  
  seg_names <- 1:length(seg$values)
  seg_names <- unlist( mapply(seg_names, seg$lengths, FUN = rep) )
  seg_sum <- by(d, seg_names, sum)
  
  seg_ends <- cumsum(seg$lengths) # last nt of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)] + 1) # first nt of a segment
  up_down <- which( seg$lengths >= range & seg$values != 0  )
  up_down_change <- diff( seg$values[up_down]) == 2
  change_seg <- up_down[which(up_down_change) + 1]
  before_change_seg <-  up_down[which(up_down_change)]
  all_min <- 0.5 + (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  return(all_min)
}

all_max <- function(y, range)
{
  d <- diff(y )
  d0 <- d
  d0[ d > 0 ] <- 1
  d0[ d < 0 ] <- -1
  seg <- rle(d0)
  
  seg_names <- 1:length(seg$values)
  seg_names <- unlist( mapply(seg_names, seg$lengths, FUN = rep) )
  seg_sum <- by(d, seg_names, sum)
  
  seg_ends <- cumsum(seg$lengths) # last nt of a segment
  seg_starts <- c(1, seg_ends[-length(seg_ends)] + 1) # first nt of a segment
  up_down <- which( seg$lengths >= range & seg$values != 0  )
  up_down_change <- diff( seg$values[up_down]) == -2
  change_seg <- up_down[which(up_down_change) + 1]
  before_change_seg <-  up_down[which(up_down_change)]
  all_max <- 0.5 + (seg_ends[before_change_seg] + seg_starts[change_seg]) / 2
  return(all_max)
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
           top: calc(50%);
           left: calc(50%);
           }
           ")
    )
  ),
  # create fluid layout with several tabs for displaying outputs
  tabsetPanel(
    # FIRST TAB
    tabPanel(tags$strong("Alignments/Areas"), icon = icon("chart-area"), # updated to newer version
             fluidRow(column(2),
                      column(4, tags$h4(tags$strong("Single profiles for quantification")
                      )),
                      column(6, tags$h4(tags$strong("Multiple aligned profiles")))),
             
             fluidRow(
               column(2, style = "padding-left:20px",
                      # create area for uploading .pks file
                      fileInput("input_data", "Upload polysome profile or previous analysis", multiple = T, accept = c(".pks",".csv", ".txt", ".RData")),
                      
                      #Let user select their loaded files and set x-anchor and baseline
                      selectInput("select", "Select files", choices = c(), width = '100%'),
                      
                      
                      
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow( 
                        column(12, tags$h4(tags$strong("Fluorescence profile")))),
                      fluidRow( 
                        column(12, tags$h6("Set fluorescence baseline"))),
                      fluidRow(         
                        # show fluorescence signal or not?
                        # set a separate baseline for the fluorescence signal
                        
                        column(6, style = "padding-top:10px",actionButton("baseline_fl", "Baseline", width = '100%')
                        ),
                        
                        column(6,style = "padding-top:8px",
                               checkboxInput("show_fl", "Show fluorescence signal", value = FALSE, width = NULL))
                        
                        
                      ),
                      
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               numericInput("axis1_fl", "Set y min", value = NULL, step = 1)
                        ),
                        column(6, style = "padding-top:10px",
                               numericInput("axis2_fl", "Set y max", value = NULL, step = 1))),
                      
                      # introduce a slider for smoothing of the fluorescence signal
                      column(12,
                             fluidRow(sliderInput("slider1", label = "Smooth fluorescence profile", min = 0, 
                                                  max = 1, value = 0, step = 0.01)))
                      
                      
               ),
               column(4, style = "padding-top:22px",
                      textOutput("test"),
                      plotOutput("plot_single", click = "click")
               ),
               column(4, style = "padding-top:22px",
                      plotOutput("plot_align")
               ),
               column(2, style = "padding-right:22px",
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow(tags$h5(tags$strong("Fluorescence profile")),
                               
                               # show fluorescence signal or not?
                               column(12,  
                                      checkboxInput("show_fl_al", "Show fluorescence signal", value = FALSE, width = NULL)
                               ),
                               
                               column(6, 
                                      numericInput("axis1_a_fl", "Set y min", value = NULL, step = 1)),
                               column(6, 
                                      numericInput("axis2_a_fl", "Set y max", value = NULL, step = 1))
                      ),
                      fluidRow(tags$h5(tags$strong("Polysome profile")),
                               column(6,
                                      numericInput("axis1_a", "Set y min", value = NULL, step = 1)),
                               column(6,
                                      numericInput("axis2_a", "Set y max", value = NULL, step = 1))
                      ),
                      fluidRow(
                        column(6, 
                               numericInput("axis3_a", "Set x min", value = NULL, step = 1)
                        ),
                        column(6,
                               numericInput("axis4_a", "Set x max", value = NULL, step = 1)
                        )
                      )
               )
             ),
             
             fluidRow(
               column(4, style = "padding-left:20px",
                      fluidRow(
                        column(6,
                               tags$h4(tags$strong("Polysome profile"))
                        ),
                        column(6,
                               tags$h4(tags$strong("Quantification"))
                        )
                      ),
                      
                      fluidRow(
                        column(6,
                               tags$h6("Set baseline and alignment anchor")
                        ),
                        column(6,
                               tags$h6("Set area start and end")
                        )
                      ),
                      
                      fluidRow(
                        
                        column(3,style = "padding-top:10px", actionButton("baseline", "Baseline", width = '100%'
                        )
                        ),
                        column(3,style = "padding-top:10px", actionButton("x_anchor", "X-anchor", width = '100%')
                        ),
                        
                        column(3, style = "padding-top:10px",
                               actionButton("file_start", "Start", width = '100%')
                        ),
                        column(3, style = "padding-top:10px",
                               actionButton("file_end", "End", width = '100%')
                        )
                      ),
                      fluidRow(
                        column(3, style = "padding-top:10px",
                               numericInput("axis1", "Set y min", value = NULL, step = 1)
                        ),
                        column(3, style = "padding-top:10px",
                               numericInput("axis2", "Set y max", value = NULL, step = 1)
                        ),
                        
                        column(3, style = "padding-top:10px",
                               selectInput("select_area", "Select area to quantify", choices = c("Total", "Monosomes", "Polysomes", "40S", "60S"), width = '100%')
                        ),
                        column(3, style = "padding-top:10px",
                               textInput("name_area", "Optional: Name area", value = "", width = NULL, placeholder = "Unknown peak")
                        )
                      ),
                      
                      fluidRow(
                        column(3, style = "padding-top:0px",
                               numericInput("axis3", "Set x min", value = NULL, step = 1)
                        ),
                        column(3, style = "padding-top:0px",
                               numericInput("axis4", "Set x max", value = NULL, step = 1)
                        ),
                        column(3, style = "padding-top:22px",
                               actionButton("quantify_area", "Quantify", width = '100%')
                        ),
                        column(3, style = "padding-top:22px",
                               actionButton("take_over_name", "Add name", width = '100%')
                        )
                      ),
                      fluidRow(
                        column(6,
                               checkboxInput("red_lines", "Show baseline and x-anchor", value = TRUE, width = NULL)
                        ),
                        column(6,
                               checkboxInput("green_lines", "Show area and lines", value = TRUE, width = NULL)
                        )
                      )
               ),
               
               column(2,
                      fluidRow(
                        downloadButton("downloadSinglePlot", "Download plot", icon = icon("file-download"), width = '100%')
                      ),
                      fluidRow(
                        radioButtons("helper_functions", "Use helper functions for setting vertical lines:",
                                     c("Help find min/max" = "peak_help",
                                       "Help find inflections" = "inflection_help",
                                       "No usage of helper functions" = "no_help"), selected = "peak_help"),
                        style = "padding-top:20px"
                      ),
                      style = "padding-left:35px"
               ),
               
               column(2,
                      # let user select files ready to align 
                      fluidRow(
                        column(12,
                               tags$h4(tags$strong("Alignment")) 
                        )
                      ),
                      fluidRow(
                        column(12,
                               tags$h6("(Baseline and x-anchor need to be set)")
                        )
                      ),
                      
                      
                      fluidRow(
                        column(8,
                               selectInput("select_alignment", "Available profiles", choices = c(), width = '100%')
                        ),
                        column(1, style = "padding-top:20px",
                               actionButton("up", label = NULL, icon = icon("angle-double-up"), width = '100%',
                                            style = "padding-left:5px"
                               )
                        ),
                        column(1, style = "padding-top:20px",
                               actionButton("down", label = NULL, icon = icon("angle-double-down"), width = '100%',
                                            style = "padding-left:5px"
                               )
                        ),
                        column(2)
                        
                      ),
                      
                      fluidRow(
                        column(12,
                               tags$h6("Remove and restore profiles from alignment"))
                      ),
                      
                      fluidRow(
                        column(6, 
                               actionButton("remove_file", "Hide", icon = icon("trash"), width = '100%')
                        ),
                        column(6, 
                               actionButton("add_file", "Show", icon = icon("trash-restore"), width = '100%')
                        ),
                        style = "padding-right:35px"
                      )
                      
                      
                      
                      
               ),
               
               column(2, #changing their color and linetype seing the changes directly in the plot
                      
                      
                      fluidRow(
                        radioButtons("color_palette", "Select default color palette:",
                                     c("Dark palette" = "dark_palette",
                                       "Rainbow palette" = "rainbow_palette",
                                       "Color blind friendly palette" = "color_blind"
                                     ), selected = "dark_palette"),
                        style = "padding-top:20px"
                      ),
                      
                      fluidRow(
                        tags$h6("Adjust selected profile yourself:")
                      ),
                      
                      fluidRow(
                        column(5, 
                               colourInput("color", "Colors", palette = "square", value = "#FFFFFF")
                        ),
                        column(4, style = "padding-left:0px", 
                               selectInput("linetype", "Type", choices = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                                           selected = "solid")
                        ),
                        column(3, style = "padding-left:0px",
                               numericInput("linewidth", "Width", value = 1, min = 0)
                        )
                      )
                      
               ),
               
               column(2,
                      downloadButton("downloadPlot", "Download plot", icon = icon("file-download"), width = '100%'),
                      checkboxInput("anchor_line", "Display x-anchor in alignment", value = TRUE, width = NULL),
                      checkboxInput("normalize_height", HTML("Normalize <b>height</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL),
                      checkboxInput("normalize_length", HTML("Normalize <b>length</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL),
                      downloadButton("export", "Export analysis", width = "100%")
                  
               )
             )
             
    ),
    # SECOND TAB
    # for user-guided peak-deconvolution
    tabPanel(tags$strong("Deconvolution"), icon = icon("cut"),
             fluidRow(column(2),
                      column(4, tags$h4(tags$strong("Single profile for deconvolution"))
                      ),
                      column(6)
             ),
             
             fluidRow(
               column(2, style = "padding-left:20px",
                      
                      #Let user select their loaded files and set x-anchor and baseline
                      selectInput("select2", "Select files", choices = c(), width = '100%'),
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow( 
                        column(12, tags$h4(tags$strong("Identify peaks")))),
                      
                      fluidRow(
                        column(12, style = "padding-left:20px",
                               tags$h5(tags$strong("Peaks can be identified from local maxima with different resolution:")))
                      ),
                      column(12,
                             fluidRow(sliderInput("slider3", 
                                                  label = "Resolution", 
                                                  min = 50, 
                                                  max = 100, value = 95)
                             )
                      ),
                      
                      fluidRow(
                        column(12, style = "padding-left:20px",
                               tags$h5(tags$strong("Hidden peaks can be identified from 2nd deriv. minima after smoothing.")))
                      ),
                      # introduce a slider for smoothing of the fluorescence signal
                      column(12,
                             fluidRow(sliderInput("slider2", 
                                                  label = "Smoothing", 
                                                  min = 0, 
                                                  max = 1, value = 0.3, step = 0.01))),
                      
                      
                      fluidRow(
                        column(6,
                               checkboxInput("show_deriv", "Show 2nd deriv. minima", value = FALSE, width = NULL)
                        ),
                        column(6,
                               checkboxInput("show_local_max", "Show local maxima", value = TRUE, width = NULL)
                        )
                      ),
                      fluidRow(
                        column(12, style = "padding-left:20px",
                               tags$h4(tags$strong("Deconvolution"))
                        )
                      ),
                      
                      fluidRow(
                        column(12, style = "padding-left:20px",
                               tags$h5(tags$strong("Select or add peaks by clicking into the plot."))
                        )
                      ),
                      
                      
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               numericInput("height", "Height", value = NULL, step = 0.5)
                        ),
                        column(6, style = "padding-top:10px",
                               numericInput("SD", "Width", value = 0.1, min = 0.01, step = 0.01)
                        )
                      ),
                      
                      fluidRow(
                        column(12, sliderInput(
                          "slider4", label = "Degree of asymmetry", value = 0, min = -0.2, max = 0.2,
                          step = 0.01)
                        )
                      ),
                      fluidRow(
                        column(12,style = "padding-top:32px", actionButton("delete_peak", "Delete peak", width = '100%',
                                                                           icon = icon("trash"))
                        )
                      )
               ),
               
               column(8, style = "padding-top:22px", 
                      plotOutput("plot_deconv", click = "click_deconv", height = 700)
               ),
               
               column(2, style = "padding-right:22px",
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      
                      fluidRow(tags$h4(tags$strong("Axis limits")),
                               column(6,
                                      numericInput("axis1_d", "Set y min", value = NULL)),
                               column(6,
                                      numericInput("axis2_d", "Set y max", value = NULL))
                      ),
                      fluidRow(
                        column(6, 
                               numericInput("axis3_d", "Set x min", value = NULL)),
                        column(6,
                               numericInput("axis4_d", "Set x max", value = NULL))
                      ),
                      fluidRow(
                        column(12,
                               tags$h4(tags$strong("Peak quantification"))
                        )
                      ),
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               selectInput("select_peak", "Peak to quantify", choices = c("40S", "60S", "Monosomes", "Halfmers"), width = '100%')
                        ),
                        column(6, style = "padding-top:10px",
                               textInput("name_peak", "Optional: Name peak", value = "", width = NULL, placeholder = "Unknown peak")
                        )
                      ),
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               actionButton("quantify_peaks", "Quantify", width = '100%')
                        ),
                        column(6, style = "padding-top:10px",
                               actionButton("take_over_peak", "Add name", width = '100%')
                        )
                      ),
                      fluidRow(
                        column(6,
                               downloadButton("downloadDeconv", "Download plot", icon = icon("file-download"), width = '100%')
                        ),
                        column(6),
                        style = "padding-top:22px"
                      )
               )
             )
    ),
    # THIRD TAB
    # shows updating table of aligned (and normalized) profiles, download possible 
    tabPanel(tags$strong("Alignment table"), icon = icon("table"), 
             tags$h4("Table of all aligned (and normalized) profiles"),
             downloadButton("downloadData", "Download .csv file"),
             tableOutput("csv_file")),
    # FOURTH TAB
    # shows updating table of quantified areas for respective plots, download possible 
    tabPanel(tags$strong("Quantification summary"), icon = icon("list-alt"), 
             fluidRow(
               column(10,
                      tags$h4("Table of quantified areas or peaks"),
                      tableOutput("quantification"),
                      downloadButton("downloadQuant", "Download .csv file")
               ),
               column(2, style = "padding-top:20px", 
                      actionButton("to_stats", "Transfer to Stats", width = '80%', icon = icon("chart-line") )
               )
             )
    ),
    # FIFTH TAB
    # for performing beta regression on replicates
    tabPanel(tags$strong("Statistics"), icon = icon("chart-line"), 
             fluidRow(
               column(4, style = "padding-top:10px",
                      selectInput("select_region", "Area or peak", choices = c(), width = '80%')
               ),
               column(8, style = "padding-top:10px",
                      tableOutput("stats_tab") 
               )
             ),
             fluidRow(
               column(4,
                      fluidRow( column(12, DTOutput("files_conditions") ) ),
                      fluidRow(
                        column(12, actionButton("add_variable", "Add column", width = "20%", icon = icon("plus") ) )
                      )
                      
               )
             )
    ),
    # SIXTH TAB 
    # shows Rmd manual  
    tabPanel(tags$strong("QuAPPro manual"), icon = icon("question-circle"),
             fluidRow(
               column(2),
               column(8, htmltools::includeMarkdown("RMD_menu_template_2022.Rmd")),
               column(2))
    ),
    
    # SIXTH TAB
    # general information and impressum
    tabPanel(tags$strong("Release notes", style = "color:blue"), icon = icon("info", style = "color:blue"), 
             tags$h4("Release notes"),
             tags$div("This is version QuAPPro v0.0.5. 
                      To run it locally, you can download the code from our GitHub repository:",
                      tags$a(href="https://github.com/johannaschott/QuAPPro", "https://github.com/johannaschott/QuAPPro", style = "color:blue")
                      
             )
    ),
    
    # SEVENTH TAB
    # general information and impressum
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
  ))



#### SERVER ####

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 50 * 1024^2)
  
  ### REACTIVE VALUES FOR LISTS BUILT ALONG THE SESSION
  # create reactive values for collecting values in a list/vector/data.frame
  val <- reactiveValues(
    paths_collected = vector(),
    colors_collected = list(),
    color_vector = vector(),
    linetype_collected = list(),
    linewidth_collected = list(),
    file_starts = list(),
    file_ends = list(),
    sum_areas = list(),
    sum_additive_areas = list(),
    file_starts_add = list(),
    file_ends_add = list(),
    files_to_align = vector(),
    files_to_plot = vector(),
    remove_files_list = vector(),
    csv_file_df = data.frame(),
    df_quant = data.frame(),
    factors_list = list(),
    control_baseline = list(),
    fl_control_baseline = list(),
    file_types = vector(),
    xmin_collected = list(),
    xmax_collected = list(),
    ymin_collected = list(),
    ymax_collected = list(),
    ymin_collected_fl = list(),
    ymax_collected_fl = list(),
    polysome_data = list(),
    files_list = list(),
    peak_list = vector(),
    align_files_list = vector(),
    area_list = vector()
    
  )
  
  
  
  ### INITIAL LOADED FILES
  # update select output field with the name of every newly loaded file
  observeEvent(input$input_data, {
    
    if(grepl("RData", input$input_data$name ) )
      
    {
      load(input$input_data$datapath)
      for(i in names(forExport))
      {
        val[[ i ]] <- forExport[[i]]
      }
      updateSelectInput(session, "select",
                        choices = val$files_list)
      updateSelectInput(session, "select2",
                        choices = val$files_list)
    }else{
    
      val$paths_collected[input$input_data$name] <- input$input_data$datapath
      val$files_list <- c(input$input_data$name[input$input_data$size != 0], val$files_list)
      
      updateSelectInput(session, "select",
                        choices = val$files_list)
      
      updateSelectInput(session, "select2",
                        choices = val$files_list)
      
      new_names <- input$input_data$name
      new_paths <- input$input_data$datapath
      
      for(i in 1:length(new_names))
      {
        this_name <- new_names[i]
        this_path <- new_paths[i]
        
        if( grepl(".pks",as.character(this_name) ) ){
          val$factors_list[[this_name]] <- (0.1/60)
          val$file_types[this_name] <- "pks"
          data <- read.table(this_path, dec = ",", header = F) 
          val$polysome_data[[this_name]] <- data[,3]
          val$xvalues[[this_name]] <- (data[ ,1]+1)*val$factors_list[[this_name]]
        }
        
        if( grepl(".csv",as.character(this_name) )  ){
          start <- grep("Data Columns:", readLines(this_path))
          data <- read.csv(this_path, skip = start)
          
          if( "SampleFluor" %in% colnames(data) )
          {
            val$factors_list[[this_name]] <- (0.32/60)
            val$file_types[this_name] <- "csv_fluo" 
            val$polysome_data[[this_name]] <- data[ ,5]
            val$fluo_data[[this_name]] <- data$SampleFluor
            val$xvalues[[this_name]] <- (1:length(data[ ,4]))*val$factors_list[[this_name]]
          }
          
          if(!"SampleFluor" %in% colnames(data) )
          {
            val$factors_list[[this_name]] <- (0.2/60)
            val$file_types[this_name] <- "csv" 
            val$polysome_data[[this_name]] <- data[ ,5]
            val$xvalues[[this_name]] <- (1:length(data[ ,5]))*val$factors_list[[this_name]]
          }
        }
        
        if( grepl(".txt",as.character(this_name) ) ){
          val$factors_list[[this_name]] <- (0.1/60)
          val$file_types[this_name] <- "pks"
          data <- read.delim(this_path)[,2]
          val$polysome_data[[this_name]] <- data[,2]
          val$xvalues[[this_name]] <-  (data[ ,1]+1)*val$factors_list[[this_name]]
        }
        
      }
    }
  })
  
  output$test <- renderText(input$input_data$name)
  
  # notification if empty file is loaded
  observeEvent(input$input_data,{
    if(length(input$input_data$name[input$input_data$size != 0]) != length(input$input_data$name)){
      if((length(input$input_data$name) - length(input$input_data$name[input$input_data$size != 0])) == 1)
      {showNotification(paste0("File ",input$input_data$name[input$input_data$size == 0]," is empty and will be ignored."),
                        duration = NULL, type = "message")
      }else if((length(input$input_data$name) - length(input$input_data$name[input$input_data$size != 0])) > 1){
        showNotification(paste0("Files ",paste(input$input_data$name[input$input_data$size == 0], collapse = ", ")," are empty and will be ignored."),
                         duration = NULL, type = "message")}
    }
  })
  
  
  
  fluorescence <- reactive({
    req(input$select) 
    req(input$show_fl)
    if( val$file_types[[input$select]] == "csv_fluo" ){
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
      if(input$helper_functions == "peak_help"){
        val$file_starts[[input$select]] <- (find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                                (closest_point(input$click$x, 
                                                                               input$click$y, 
                                                                               val$xvalues[[input$select]],
                                                                               val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$file_starts[[input$select]] <- (find_closest_inflections(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                                     (closest_point(input$click$x, 
                                                                                    input$click$y, 
                                                                                    val$xvalues[[input$select]],
                                                                                    val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else{
        val$file_starts[[input$select]] <-(closest_point(input$click$x,
                                                         input$click$y,
                                                         val$xvalues[[input$select]],
                                                         val$polysome_data[[input$select]])/val$factors_list[[input$select]])
      }
    # setting profile ends
    if(val$buttons == 5)
      if(input$helper_functions == "peak_help"){
        val$file_ends[[input$select]] <- (find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             val$xvalues[[input$select]],
                                                                             val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$file_ends[[input$select]] <- (find_closest_inflections(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                                   (closest_point(input$click$x, 
                                                                                  input$click$y, 
                                                                                  val$xvalues[[input$select]],
                                                                                  val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else{
        val$file_ends[[input$select]] <-(closest_point(input$click$x,
                                                       input$click$y,
                                                       val$xvalues[[input$select]],
                                                       val$polysome_data[[input$select]])/val$factors_list[[input$select]])
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
      if(input$helper_functions == "peak_help"){
        val$anchor[input$select] <- round(find_closest_minmax(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             val$xvalues[[input$select]],
                                                                             val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else if(input$helper_functions == "inflection_help"){
        
        val$anchor[input$select] <- (find_closest_inflections(smooth_profile(val$polysome_data[[input$select]], NULL),
                                                              (closest_point(input$click$x, 
                                                                             input$click$y, 
                                                                             val$xvalues[[input$select]],
                                                                             val$polysome_data[[input$select]]))/val$factors_list[[input$select]], 5))
      }else{
        val$anchor[input$select] <-(closest_point(input$click$x,
                                                  input$click$y,
                                                  val$xvalues[[input$select]],
                                                  val$polysome_data[[input$select]])/val$factors_list[[input$select]])
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
  
  # update choices of areas to be quantified with names given by the user
  observeEvent(input$take_over_name,{
    val$area_list <- c(input$name_area, val$area_list)
    updateSelectInput(session, "select_area",
                      choices = c(val$area_list, "Total", "Monosomes", "Polysomes", "40S", "60S"))
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
  
  
  
  # create "remove vector" containing files selected by remove button. Only take files into the vector that are not already there
  observeEvent(input$remove_file,{
    req(val$files_to_align) 
    val$remove_files_list <- unique( c(val$remove_files_list, input$select_alignment) )
    
  })
  # Create "add vector" for removing files again from remove list (so addition to alignment again)
  observeEvent(input$add_file,{
    req(val$files_to_align)
    add_vector <- c((intersect(names(val$baseline[input$select_alignment]),
                               names(val$anchor[input$select_alignment]))), val$add_files_list)
    val$remove_files_list <- val$remove_files_list[!val$remove_files_list %in% add_vector] 
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
    updateColourInput(session, "color", value = val$color_vector[[input$select_alignment]])
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
  
  observe({
    updateNumericInput(session, "axis3", value = round(xmin_single(), digits = 2 ) )
    updateNumericInput(session, "axis4", value = round(xmax_single(), digits = 2 ) )
    updateNumericInput(session, "axis1", value = round(ymin_single()/lost_num_pol(), digits = 2) )
    updateNumericInput(session, "axis2", value = round(ymax_single()/lost_num_pol(), digits = 2) )
  })
  
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
  observe({
    req(input$select)
    if( val$file_types[input$select] != "csv_fluo" & input$show_fl){
      updateNumericInput(session, "axis1_fl", value = round( ymin_single_fl()/lost_num_fl(), digits = 2 ) )
      updateNumericInput(session, "axis2_fl", value = round( ymax_single_fl()/lost_num_fl(), digits = 2 ) )
    }else{
      updateNumericInput(session, "axis1_fl", value = "")
      updateNumericInput(session, "axis2_fl", value = "")
    } 
  })
  
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
  plot_singleFl <-function(){
    if(lost_num_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_fl(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select]], fluorescence(), type = "l", xaxt = "n", col = "darkgreen", lwd = 2,
         xlim =c(xmin_single(),xmax_single()), ylim = c(ymin_single_fl(), ymax_single_fl()),
         las = 1, ylab = ylab, mgp = c(3.5, 0.8, 0), xlab = "", yaxt = "n", cex.lab = 1.6)
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_fl(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 1.6)
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
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = 2)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = 2)
    }
    if(input$red_lines){
      abline(h = val$baseline_fl[input$select],col = "red", lty=2, lwd = 2)
    }
  }
  
  plot_singlePol <-function(){
    if(lost_num_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_pol(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select]], val$polysome_data[[input$select]], type = "l", las = 1, lwd = 2,
         ylab = ylab, xlab = "Time (min)",
         ylim = c(ymin_single(),ymax_single()), 
         xlim =c(xmin_single(),xmax_single()), mgp = c(3.5, 0.8, 0),
         yaxt = "n", xaxt = "n", cex.lab = 1.6
    )
    
    
    axis(1, las = 1, mgp = c(3.5, 1.2, 0), cex.axis = 1.6)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 1.6)
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      if(isTruthy(val$control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, val$xvalues[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], x_last),
                c(val$control_baseline[[input$select_area]][input$select], val$polysome_data[[input$select]][(which(val$xvalues[[input$select]] == (x_first))):(which(val$xvalues[[input$select]] == (x_last)))], val$control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "black", lwd = 2)
      }
    }
    # selected x-anchor and baseline are displayed if box is ticked
    if(input$red_lines){
      abline(v = val$anchor[input$select]*val$factors_list[[input$select]],col = "red", lty=2, lwd = 2)
      abline(h = val$baseline[input$select],col = "red", lty=2, lwd = 2)
    }
    # selected area starts and ends are displayed if box is ticked
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = 2)
      abline(v = val$file_ends[[input$select]]*val$factors_list[[input$select]],col = "#238b45", lty=2, lwd = 2)
    }
  }
  
  
  
  plot_singleInput <- function(){
    if(!is.null(val$xvalues[[input$select]])){
      if(input$show_fl && val$file_types[input$select] == "csv_fluo" ){
        if(val$buttons == 6){
          layout(matrix(2:1, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
          par(mar = c(5, 5, 0, 2))
          plot_singlePol()
          par(mar = c(0, 5, 0.5, 2))
          plot_singleFl()
        }else{
          layout(matrix(1:2, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
          par(mar = c(0, 5, 0.5, 2))
          plot_singleFl()
          par(mar = c(5, 5, 0, 2))
          plot_singlePol()
        }
      }else{
        par(mar = c(5, 5, 0.5, 2))
        plot_singlePol()
      }
    }
  }
  
  output$plot_single <- renderPlot({
    plot_singleInput()
  })
  
  # Enable download of current plot as pdf
  
  output$downloadSinglePlot <- downloadHandler(
    filename = function(){ 
      paste( 
        sub(pattern = "(.*?)\\..*$", replacement = "\\1", as.character(input$select) ),
        ".pdf", sep = ""
      )
    },
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_singleInput() )
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
  
  # Opportunity to remove and add files again to files_to_plot vector
  observe({
    val$files_to_plot <- val$files_to_align[!val$files_to_align %in% val$remove_files_list]
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
  
  # if color palette choice gets changed by the user, the original color_vector is reset to an empty list again
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
    val$color_vector <- as.list(dummy)
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
  
  
  observe({
    if(length(val$files_to_plot >= 1))
    {
      updateNumericInput(session, "axis3_a", value = round( xmin_aligned() /lost_num_al_Index(), digits = 2 ) )
      updateNumericInput(session, "axis4_a", value = round( xmax_aligned() /lost_num_al_Index(), digits = 2 ) )
      updateNumericInput(session, "axis1_a", value = round( ymin_aligned()/lost_num_al_pol(), digits = 2 ) )
      updateNumericInput(session, "axis2_a", value = round( ymax_aligned()/lost_num_al_pol(), digits = 2 ) )
    }
  })
  
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
  
  plot_alignedPol <- function(){
    if(lost_num_al_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_al_pol(), ")", sep = "")
    }
    
    f <- val$files_to_plot[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    plot(x, values_list()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]], 
         ylab = ylab, xlab = "Relative position", las = 1,
         col = val$color_vector[[f]], mgp = c(3.5, 1, 0), 
         ylim = c(val$ymin,val$ymax), xlim = c(val$xmin,val$xmax),
         yaxt = "n", xaxt = "n", cex.lab = 1.6
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
      abline(v = anchor_line, col = "red", lty=2, lwd = 2)
    }
    
    a <- axTicks(1)
    axis(1, at = a, labels = a/lost_num_al_Index(), las = 1, mgp = c(3.5, 1.2, 0), cex.axis = 1.6)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_pol(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 1.6)
    
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
             lwd = val$linewidth_collected[[f]], col = val$color_vector[[f]])
    }
    legend("topright", legend = val$files_to_plot, lty = unlist(val$linetype_collected[val$files_to_plot]), 
           lwd = unlist(val$linewidth_collected[val$files_to_plot]), col = unlist(val$color_vector[val$files_to_plot]),
           bty = "n", cex = 1.5
    )
    
    #create reactive dataframe of all plots to have access outside of renderPlot function
    val$csv_file_df <- csv_file_df
  }
  
  plot_alignedFluo <- function(){
    
    if(lost_num_al_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_al_fl(), ")", sep = "")
    }
    
    files_in_al_fluo <- val$files_to_plot[val$files_to_plot %in% intersect(names(val$baseline_fl),  names(val$baseline))] # this is necessary to maintain the order as given by val$files_to_plot
    f <-  files_in_al_fluo[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    plot(x, values_fluorescence()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]],
         ylab = ylab, xlab = "", las = 1,
         col = val$color_vector[[f]], mgp = c(3.5, 0.8, 0), 
         ylim = c(val$ymin_fl,val$ymax_fl), xlim = c(val$xmin,val$xmax),
         yaxt = "n", xaxt = "n", cex.lab = 1.6
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
    axis(2, at = a, labels = a/lost_num_al_fl(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 1.6)
    
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
             lwd = val$linewidth_collected[[f]], col = val$color_vector[[f]])
    }
    
    #create reactive dataframe of all plots to have access outside of renderPlot function
    val$csv_file_df_all <- merge(val$csv_file_df, csv_file_df_fluo, by="Index", all = T)
  }
  
  plot_alignment <- function(){
    if(input$show_fl_al){
      layout(matrix(1:2, 2, 1), height = c(0.5, 1) ) # divides the plotting area into 2 rows
      par(mar = c(0, 5, 0.5, 2))
      plot_alignedFluo()
      par(mar = c(5, 5, 0, 2)) 
      plot_alignedPol()
    }else{
      par(mar = c(5, 5, 0.5, 2))
      plot_alignedPol()
    }
  }
  
  output$plot_align <- renderPlot({
    req(val$files_to_plot)
    if(length(val$files_to_plot) >= 1)
    {
      plot_alignment()
    }
  })
  
  # Enable download of current plot as pdf
  
  output$downloadPlot <- downloadHandler(
    filename = "alignment.pdf",
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_alignment() )
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
  
  
  ### Peak de-convolution:
  yvalue_deconv <- function()
  {
    val$polysome_data[[input$select2]] - val$baseline[input$select2]
  }
  
  observeEvent(input$select2, {
    updateSelectInput(session, "select", selected = input$select2)
  })
  
  observeEvent(input$select, {
    updateSelectInput(session, "select2", selected = input$select)
  })
  
  # update choices of peaks to be quantified with names given by the user
  observeEvent(input$take_over_peak,{
    val$peak_list <- c(input$name_peak, val$peak_list)
    updateSelectInput(session, "select_peak",
                      choices = c(val$peak_list, "40S", "60S", "80S", "Halfmers"))
  })
  
  second_deriv_smoothed <- function()
  {
    smooth_profile( second_deriv( smooth_profile(yvalue_deconv(), input$slider2 )), input$slider2)
  }
  
  second_deriv_min <- function()
  {
    resolution <- (1 - input$slider3/100 )*100
    val$xvalues[[input$select2]][-c(1:2)] [ all_min( second_deriv_smoothed(), resolution ) ]
  }
  
  profile_smoothed <- function()
  {
    smooth_profile(yvalue_deconv(), input$slider2 )
  }
  
  local_max <- function()
  {
    resolution <- (1 - input$slider3/100 )*100
    val$xvalues[[input$select2]][ all_max( yvalue_deconv(), resolution ) ]
  }
  
  
  plot_2nd_deriv <-function(){
    plot(val$xvalues[[input$select2]][-c(1:2)], second_deriv_smoothed(), type = "l", lwd = 2, xaxt = "n", col = "black", 
         las = 1, ylab = "2nd derivative", mgp = c(3.5, 0.8, 0), xlab = "", yaxt = "n", cex.lab = 1.6,
         xlim = c(xmin_d(), xmax_d()) )
    a <- axTicks(2)
  }
  
  # When user selects other axis limits, limit values change 
  val <- reactiveValues(
    xmin_d_collected = list(),
    xmax_d_collected = list(),
    ymin_d_collected = list(),
    ymax_d_collected = list()
  )
  
  # save axis limits set by the user for single profiles
  observeEvent(input$axis3_d, {
    val$xmin_d_collected[[input$select2]] <- input$axis3_d
  })
  observeEvent(input$axis4_d, {
    val$xmax_d_collected[[input$select2]] <- input$axis4_d
  })
  observeEvent(input$axis1_d, {
    val$ymin_d_collected[[input$select2]] <- input$axis1_d
  })
  observeEvent(input$axis2_d, {
    val$ymax_d_collected[[input$select2]] <- input$axis2_d
  })
  
  ymin_d <- reactive({
    if(isTruthy(val$ymin_d_collected[[input$select2]])){
      val$ymin_d_collected[[input$select2]]*lost_num_pol()
    }else{
      0
    }})
  
  
  ymax_d <- reactive({
    if(isTruthy(val$ymax_d_collected[[input$select2]])){
      val$ymax_d_collected[[input$select2]]*lost_num_pol()
    }else{
      max(yvalue_deconv())
    }})
  
  
  xmin_d <- reactive({
    if(isTruthy(val$xmin_d_collected[[input$select2]])){
      val$xmin_d_collected[[input$select2]]
    }else{
      0
    }})
  
  xmax_d <- reactive({
    if(isTruthy(val$xmax_d_collected[[input$select2]])){
      val$xmax_d_collected[[input$select2]]
    }else{
      max(val$xvalues[[input$select2]])
    }})
  
  observe({
    req(yvalue_deconv())
    updateNumericInput(session, "axis3_d", value = round(xmin_d(), digits = 2 ) )
    updateNumericInput(session, "axis4_d", value = round(xmax_d(), digits = 2 ) )
    updateNumericInput(session, "axis1_d", value = round(ymin_d()/lost_num_pol(), digits = 2) )
    updateNumericInput(session, "axis2_d", value = round(ymax_d()/lost_num_pol(), digits = 2) )
  })
  
  # add peak or change active peak by clicking into the plot:
  # store peak information:
  val <- reactiveValues(
    peak_pos = list(), 
    peak_height = list(),
    peak_sd = list(),
    active_peak = list(),
    peak_asym = list(),
    peak_type = list(),
    peak_values = list()
  )
  
  
  
  generate_peak <- function(x, peak_height, peak_pos, peak_sd, peak_asym)
  {
    # Normal distribution for symmetric peaks
    if( peak_asym == 0)
    {
      Gauss <- dnorm(x, mean = peak_pos, sd = peak_sd )
      cor_factor <- peak_height/max(Gauss)
      return(Gauss*cor_factor)
    }
    
    
    # Exponentially modified Gauss, unless asymmetry is zero
    if(peak_asym > 0)
    {
      lambda <- 1/peak_asym
      y1 <- demg(x, mu = peak_pos, sigma = peak_sd, lambda = lambda)
      pos_of_max <- x[ order(y1, decreasing = TRUE)[1] ]
      diff_to_mean <- pos_of_max - peak_pos
      Gauss <- demg(x, mu = peak_pos - diff_to_mean, sigma = peak_sd, lambda = lambda)
    }
    
    if(peak_asym < 0)
    {
      lambda <- -1/peak_asym
      # mirror position of peak along middle of profile
      peak_pos_mirror <- median(x) + ( median(x) - peak_pos )
      y1 <- demg(x, mu = peak_pos_mirror, sigma = peak_sd, lambda = lambda)
      pos_of_max <- x[ order(y1, decreasing = TRUE)[1] ]
      diff_to_mean <- pos_of_max - peak_pos_mirror
      Gauss <- rev( demg(x, mu = peak_pos_mirror - diff_to_mean, sigma = peak_sd, lambda = lambda) )
    }
    cor_factor <- peak_height/max(Gauss)
    return(Gauss*cor_factor)
  }
  
  
  sum_of_Peaks <- function(x, peak_height, peak_pos, peak_sd, peak_asym, peak_type){
    y <- matrix(NA, length(peak_pos), length(x))
    for(i in 1:length(peak_pos))
    {
      y[i,] <- val$peak_values[[input$select2]][[ i ]]
    }
    return(apply(y, 2, sum))
  }
  
  current_model <- function(){
    sum_of_Peaks(val$xvalues[[input$select2]], val$peak_height[[input$select2]], 
                 val$peak_pos[[input$select2]], val$peak_sd[[input$select2]],
                 val$peak_asym[[input$select2]],
                 val$peak_type[[input$select2]])
  }
  
  
  # Which peak positions are available?
  # Choose from local maxima or second derivative minima, depending on 
  # the selection
  pos_to_choose <- function()
  {
    if(input$show_local_max)
    {
      return( local_max() )
    }
    if(input$show_deriv)
    {
      return( second_deriv_min() )
    }
    if( !input$show_local_max & !input$show_deriv )
    {
      return(val$peak_pos[[input$select2]])
    }
  }
  
  observeEvent(input$click_deconv$x, {
    req(pos_to_choose())
    if(isTruthy(val$baseline[input$select2]))
    {
      # determine peak position that is closest to the click:
      closest <- order( abs(pos_to_choose() - input$click_deconv$x) )[1]
      closest_peak <- pos_to_choose()[closest]
      # When closest peak has not been generated before:
      if( !(closest_peak %in% val$peak_pos[[input$select2]]) )
      {
        position_on_list <- length( val$peak_values[[input$select2]] ) + 1
        val$active_peak[[input$select2]] <- position_on_list
        
        if(length(val$peak_pos[[input$select2]]) > 0 )
        {
          val$peak_height[[input$select2]] <- c(val$peak_height[[input$select2]] ,( yvalue_deconv() - current_model() )[val$xvalues[[input$select2]] == closest_peak ])
        }else{
          val$peak_height[[input$select2]] <- c(val$peak_height[[input$select2]], yvalue_deconv()[val$xvalues[[input$select2]] == closest_peak ])
        }
        val$peak_sd[[input$select2]] <- c(val$peak_sd[[input$select2]], input$SD)
        val$peak_pos[[input$select2]] <- c(val$peak_pos[[input$select2]], closest_peak)
        val$peak_asym[[input$select2]] <- c(val$peak_asym[[input$select2]], input$slider4)
        val$peak_values[[input$select2]][[position_on_list]] <- generate_peak(val$xvalues[[input$select2]], 
                                                                              val$peak_height[[input$select2]][val$active_peak[[input$select2]] ], 
                                                                              val$peak_pos[[input$select2]][val$active_peak[[input$select2]] ],
                                                                              val$peak_sd[[input$select2]][val$active_peak[[input$select2]] ], 
                                                                              val$peak_asym[[input$select2]][val$active_peak[[input$select2]] ])
        
      }else{
        val$active_peak[[input$select2]] <- which(val$peak_pos[[input$select2]] == closest_peak)
        updateNumericInput(session, "SD", value = val$peak_sd[[input$select2]][val$active_peak[[input$select2]]])
        updateNumericInput(session, "height", value = round( val$peak_height[[input$select2]][val$active_peak[[input$select2]]], digits = 2) )
        updateSliderInput(session, "slider4", value = val$peak_asym[[input$select2]][val$active_peak[[input$select2]]])
      }
    }
  })
  
  
  
  
  # update parameters of active peak, when user input changes
  observeEvent(input$SD, {
    req(val$active_peak[[input$select2]])
    val$peak_sd[[input$select2]][val$active_peak[[input$select2]]] <- input$SD
  })
  
  observeEvent(input$height, {
    req(val$active_peak[[input$select2]])
    val$peak_height[[input$select2]][val$active_peak[[input$select2]]] <- input$height
  })
  
  observeEvent(input$slider4, {
    req(val$active_peak[[input$select2]])
    val$peak_asym[[input$select2]][ val$active_peak[[input$select2]] ] <- input$slider4
  })
  
  # update values of peak when parameters change
  observe({
    req(val$active_peak[[input$select2]])
    val$peak_values[[input$select2]][[ val$active_peak[[input$select2]] ]] <- generate_peak(val$xvalues[[input$select2]], 
                                                                                            val$peak_height[[input$select2]][val$active_peak[[input$select2]] ], 
                                                                                            val$peak_pos[[input$select2]][val$active_peak[[input$select2]] ],
                                                                                            val$peak_sd[[input$select2]][val$active_peak[[input$select2]] ], 
                                                                                            val$peak_asym[[input$select2]][val$active_peak[[input$select2]] ]
    )
  })
  
  observeEvent(input$delete_peak, {
    val$peak_values[[input$select2]] <- val$peak_values[[input$select2]][-val$active_peak[[input$select2]]]
    val$peak_pos[[input$select2]] <- val$peak_pos[[input$select2]][-val$active_peak[[input$select2]] ]
    val$peak_height[[input$select2]] <- val$peak_height[[input$select2]][-val$active_peak[[input$select2]] ]
    val$peak_sd[[input$select2]] <- val$peak_sd[[input$select2]][-val$active_peak[[input$select2]] ]
    val$peak_asym[[input$select2]] <- val$peak_asym[[input$select2]][-val$active_peak[[input$select2]] ]
    val$peak_type[[input$select2]] <- val$peak_type[[input$select2]][-val$active_peak[[input$select2]] ]
    val$active_peak[[input$select2]] <- NULL
  })
  
  # only 2nd derivative minima or local maxima can be shown
  observeEvent(input$show_deriv, {
    if(input$show_deriv)
    {
      updateCheckboxInput(session, "show_local_max", value = FALSE)
    }
  })
  
  observeEvent(input$show_local_max, {
    if(input$show_local_max)
    {
      updateCheckboxInput(session, "show_deriv", value = FALSE)
    }
  })
  
  
  # Calculate sum of peak values to quantify peaks
  observeEvent(input$quantify_peaks,{
    req(val$active_peak[[input$select2]])
    curve <- val$peak_values[[input$select2]][[ val$active_peak[[input$select2]]  ]]
    val$sum_areas[[input$select_peak]][input$select2] <- sum(curve)
    # store peak positions of quantified peaks to display them again when the user selects them 
    val$quantified_peaks_pos[[input$select_peak]][input$select2] <- val$peak_pos[[input$select]]
  })
  
  
  # Change active peak when selected from the menu by the user:
  observeEvent(input$select_peak, {
    req(val$quantified_peaks_pos[[input$select_peak]][input$select2])
    val$active_peak[[input$select2]] <- which(val$peak_pos[[input$select2]] == val$quantified_peaks_pos[[input$select_peak]][input$select2]  )
  })
  
  plot_singlePol_deconv <-function(){
    if(lost_num_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_pol(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select2]], yvalue_deconv(), type = "l", lwd = 2, las = 1,
         ylab = ylab, xlab = "Time (min)",
         ylim = c(ymin_d(),ymax_d()), 
         xlim =c(xmin_d(),xmax_d()), mgp = c(3.5, 0.8, 0),
         yaxt = "n", xaxt = "n", cex.lab = 1.6
    )
    
    
    axis(1, las = 1, mgp = c(3.5, 1.2, 0), cex.axis = 1.6)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 1.6)
  }
  
  # add curves of peaks and model:
  # Plot active peak last
  plot_peaks <- function(){
    req(val$peak_pos[[input$select2]])
    if(length(val$peak_pos[[input$select2]]) > length( val$active_peak[[input$select2]] ) )
    {
      # plot inactive peaks
      inactive_peaks <- which( !( ( 1:length(val$peak_pos[[input$select2]] ) ) %in% val$active_peak[[input$select2]] ) )
      for(peak in inactive_peaks )
      {
        curve <- val$peak_values[[input$select2]][[ peak ]]
        points(val$xvalues[[input$select2]], curve,
               type = "l")
        polygon( c(val$xvalues[[input$select2]][1], val$xvalues[[input$select2]], tail(val$xvalues[[input$select2]], 1) ), 
                 c(0, curve, 0), col = rgb(0, 0, 0.9, alpha = 0.5),
                 border = NA)
      }
    }
    
    if(!is.null(val$active_peak[[input$select2]]))
    {
      # plot active peak last
      peak <- val$active_peak[[input$select2]]
      curve <- val$peak_values[[input$select2]][[ peak ]]
      points(val$xvalues[[input$select2]], curve,
             type = "l")
      polygon(c(val$xvalues[[input$select2]][1], val$xvalues[[input$select2]], tail(val$xvalues[[input$select2]], 1) ), 
              c(0, curve, 0), col = rgb(0, 0.9, 0, alpha = 0.5),
              border = NA)
    }
    
    points(val$xvalues[[input$select2]], current_model(), type = "l", 
           col = "grey3", lwd = 2, lty = 2)
    
  }
  
  plot_singleInput_deconv <- function(){
    if(!is.null(val$xvalues[[input$select2]])){
      if(input$show_deriv & !input$show_local_max)
      {
        layout(matrix(1:2, 2, 1), height = c(0.5, 1) ) 
        par(mar = c(0, 5, 0.5, 2))
        plot_2nd_deriv()
        abline( v = second_deriv_min(), col = "red", lty = 2, lwd = 2 )
        par(mar = c(5, 5, 0, 2))
        plot_singlePol_deconv()
        abline( v = second_deriv_min(), col = "red", lty = 2, lwd = 2 )
        if(length(val$peak_pos[[input$select2]]) > 0){
          plot_peaks()
        }
        
      }
      if(input$show_local_max & !input$show_deriv)
      {
        par(mar = c(5, 5, 0.5, 2))
        plot_singlePol_deconv()
        abline( v = local_max(), col = "blue", lty = 2, lwd = 2 )
        if(length(val$peak_pos[[input$select2]]) > 0){
          plot_peaks()
        }
      }
      if(!input$show_local_max & !input$show_deriv)
      {
        par(mar = c(5, 5, 0.5, 2))
        plot_singlePol_deconv()
        if(length(val$peak_pos[[input$select2]]) > 0){
          plot_peaks()
        }
      }
    }
  }
  
  show_note <- function(){
    plot(0,type='n',axes=FALSE,ann=FALSE)
    legend("topleft", legend = "Please select a baseline for your profile.",
           bty = "n", cex = 2)
  }
  
  
  output$plot_deconv <- renderPlot({
    req(val$xvalues[[input$select2]])
    if( isTruthy( val$baseline[input$select2] ) )
    {
      plot_singleInput_deconv()
    }else{
      show_note()
    }
  })
  
  # Enable download of current plot as pdf
  
  output$downloadDeconv <- downloadHandler(
    filename = function(){ 
      paste( 
        sub(pattern = "(.*?)\\..*$", replacement = "\\1", as.character(input$select) ),
        ".pdf", sep = ""
      )
    },
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_singleInput_deconv() )
      dev.off()
    })  
  
  # Transfer quantification to statistics panel:
  observeEvent(input$to_stats,{
    files <- unique( val$df_quant[1] )
    val$conditions_tab <- cbind( files, 
                                 rep("", length(files) ) ,
                                 rep("", length(files) )
    )
    colnames(val$conditions_tab) <- c("File", "Variable1", "Variable2")
    
    available_regions <- setdiff( colnames(val$df_quant), c("File", "Total") )
    updateSelectInput(session, "select_region",
                      choices = available_regions)
  })
  
  observeEvent(input$add_variable, {
    new_col <- rep("", dim(val$conditions_tab)[1] )
    val$conditions_tab <- cbind(val$conditions_tab, new_col)
    colnames(val$conditions_tab)[dim(val$conditions_tab)[2]] <- paste("Variable", dim(val$conditions_tab)[2] - 1, sep = "")
  })
  
  observeEvent(input$select_region, {
    Proportions <- val$df_quant[, input$select_region]/val$df_quant[,"Total"]
    File <- val$df_quant[,1]
    Conditions <- val$conditions_tab[,-1]
    val$stats_tab <- cbind(File, Proportions, Conditions)
  })
  
  
  #### Replace data 
  observeEvent(input$files_conditions_cell_edit, {
    val$conditions_tab <<- editData(val$conditions_tab, 
                                    input$files_conditions_cell_edit,
                                    "files_conditions")
    
    Proportions <- val$df_quant[, input$select_region]/val$df_quant[,"Total"]
    File <- val$df_quant[,1]
    Conditions <- val$conditions_tab[,-1]
    val$stats_tab <- cbind(File, Proportions, Conditions)
  })
  
  output$files_conditions <- renderDT(
    val$conditions_tab, editable = "cell",
    options = list(searching=FALSE)
  )
  
  # create output table showing data for statistics
  output$stats_tab <- renderTable(
    val$stats_tab
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
