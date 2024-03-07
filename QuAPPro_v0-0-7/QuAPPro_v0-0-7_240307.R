#####################################################
#     ____                 _____  _____             #
#    / __ \          /\   |  __ \|  __ \            #
#   | |  | |_   _   /  \  | |__) | |__) | __ ___    #
#   | |  | | | | | / /\ \ |  ___/|  ___/ '__/ _ \   #
#   | |__| | |_| |/ ____ \| |    | |   | | | (_) |  #
#    \___\_\\__,_/_/    \_\_|    |_|   |_|  \___/   #
#                                                   #                                                
#####################################################                                                 

# This is version 0.0.7 of QuAPPro, a shiny app for 
# interactive alignment and quantification
# of polysome profiles. It is available under the 
# GNU General Public License v3.0.

# Authors: Chiara Schiller, Johanna Schott

# Structure:

#### REQUIRED PACKAGES ####
#### FUNCTIONS ####
#### USER INTERFACE ####
#### SERVER ####
## COLLECTION OF REACTIVE VALUES
## LOADING FILES
## USER INPUT FOR PLOTTING INDIVIDUAL PROFILES
## PLOT INDIVIDUAL PROFILES
## INTERACTION OF THE USER WITH THE PROFILE PLOT
## GENERATE ALIGNMENT
## DISPLAY AND MODIFY ALIGNMENT
## PLOT ALIGNMENTS
## DOWNLOADS
## NOTIFICATIONS

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
           tags$img(src='logo_mfm_ukm.svg', align = "right", width = '80%', height = '80%')
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
             fluidRow(column(2, tags$h4(tags$strong("File import")), 
                             style = "padding-left:20px"),
                      column(4, tags$h4(tags$strong("Individual profiles for quantification") ) ),
                      column(4, tags$h4(tags$strong("Aligned profiles") ) ),
                      column(2, tags$h4(tags$strong("File export")) )
             ),
             fluidRow(
               # Section for uploading and selecting files:
               column(2, style = "padding-left:20px",
                      # fileInput for uploading files
                      fileInput("input_data", label = NULL, multiple = T, accept = c(".pks",".csv", ".txt", ".RData")),
                      
                      # Menu for selecting file type
                      fluidRow(
                        column(5, 
                               selectInput("type", "Profile type", 
                                           choices = c("PeakTrak", "Biocomp", "Äkta", "Custom"),
                                           selected = "PeakTrak",
                                           width = '100%')
                        ),
                        column(4,
                               selectInput("sep", "Delim.", 
                                           choices = c("space", "tab", ",", ";"), 
                                           selected = "space",
                                           width = '100%')
                        ),
                        column(3,
                               selectInput("dec", "Dec.", 
                                           choices = c(".", ","), selected = ",",
                                           width = '100%')
                        )
                      ),
                      
                      fluidRow( 
                        column(6, 
                               numericInput("skip", "Skip lines", value = 0, step = 1)
                        ),
                        column(3,
                               numericInput("UV_col", "UV", 
                                            value = 3, step = 1)
                        ),
                        column(3,
                               numericInput("fluo_col", "Fluor.", 
                                            value = NA, step = 1)
                        )
                      ),
                      
                      # Menu for selecting which of the uploaded files should be displayed
                      fluidRow(style = "padding-left:15px;padding-top:0px;padding-right:15px",
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
                      
                      # Download entire analysis
                      fluidRow(
                        style = "padding-bottom:10px;padding-left:15px",
                        tags$h6("Export your entire analysis as .RData. You can re-import it via File import."),
                        downloadButton("export", "Export analysis", width = "100%")
                      ), 
                      
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
                      column(8, 
                             fluidRow( tags$h4(tags$strong("Fluorescence profile")) ) 
                      ),
                      # Show fluorescence signal or not?
                      column(4, checkboxInput("show_fl_al", "Show", value = FALSE, width = '100%')
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
                               numericInput("smooth_spar", "Smooth profile", value = 0, min = 0, max = 1, step = 0.01) 
                        )
                      ),
                      fluidRow(
                        column(6, style = "padding-top:0px",
                               numericInput("axis1_fl", "Set y min", value = NULL, step = 1) ),
                        column(6, style = "padding-top:0px",
                               numericInput("axis2_fl", "Set y max", value = NULL, step = 1) )
                      ),
                      fluidRow(
                        column(12,
                               style = "padding-left:0px;padding-top:20px",
                               tags$img(src='logo_umm.svg', align = "left", width = '75%', height = '75%')
                        )
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
                             checkboxInput("red_lines2", "Show baseline", value = TRUE, width = NULL),
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
                               tags$img(src='mi3_logo.png', align = "right", width = '12%', height = '12%')
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

  options(shiny.maxRequestSize = 50 * 1024^2) # file size limit: 50 MB; previous analyses can be loaded as .RData and might be rather large
  
  ## COLLECTION OF REACTIVE VALUES
  val <- reactiveValues(
    baseline = c(), # baseline of polysome profiles as vector
    baseline_fl = list(), # baseline of fluorescence profiles
    colors_collected = list(), # profile-specific colors set by the user
    color_list = list(), # colors defined initially from a color palette
    linetype_collected = list(), # profile-specific line types set by the user
    linewidth_collected = list(), # profile-specific line widths set by the user
    file_starts = list(), # selected start points for area quantification
    file_ends = list(), # selected end points for area quantification
    sum_areas = list(), # sum of data points from start to end point of a quantification
    files_to_align = vector(), # profiles with baseline and x-anchor (which are automatically part of the alignment)
    files_to_plot = vector(), # profiles that should be shown in the alignment plot (individual profiles can be omitted by the user)
    csv_file_df = data.frame(), # data frame with x- and y-values of all aligned polysome profiles
    csv_file_df_fluo = data.frame(), # data frame with x- and y-values of all aligned fluorescence profiles
    csv_file_df_all = data.frame(), # merged data frame with x- and y-values of all aligned polysome and fluorescence profiles
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
    polysome_data = list(), # y-values of the polysome profiles 
    xvalues = list(), # x-values of the polysome profiles in minutes 
    files_list = list(), # all files that have been loaded so far
    files_to_align = vector(), # profiles with baseline and x-anchor (which are automatically part of the alignment)
    # redundant to files_to_align? Is a vector, not a list!
    legend_names = vector(), # profile names as they should appear in the legend
    fluo_data = list(), # y-values of fluorescence profile data
    smooth_fluo = 0, # spar parameter for smoothing fluorescence profiles with smooth.spline
    checkboxes = list(), # list of all checkbox inputs (for restoring them when a previous analysis is loaded as .RData)
    color_palette = "dark_palette", # selected color palette
    type = "PeakTrak", # selected file type for import
    buttons = NA
    )
  
  ## LOADING FILES
  
  # Adjust import parameters (column separator, decimal separator 
  # and lines to skip) to selected file type:
  observeEvent( input$type, {
    if(input$type == "PeakTrak"){
      updateSelectInput(session, "sep", selected = "space")
      updateSelectInput(session, "dec", selected = ",")
      updateNumericInput(session, "skip", value = 0)
      updateNumericInput(session, "UV_col", value = 3)
      updateNumericInput(session, "fluo_col", value = NA)
    }
    if(input$type == "Biocomp"){
      updateSelectInput(session, "sep", selected = ",")
      updateSelectInput(session, "dec", selected = ".")
      updateNumericInput(session, "skip", value = 51)
      updateNumericInput(session, "UV_col", value = 5)
      updateNumericInput(session, "fluo_col", value = 3)
    }
    if(input$type == "Äkta"){
      updateSelectInput(session, "sep", selected = "tab")
      updateSelectInput(session, "dec", selected = ".")
      updateNumericInput(session, "skip", value = 2)
      updateNumericInput(session, "UV_col", value = 2)
      updateNumericInput(session, "fluo_col", value = NA)
    }
    if(input$type == "Custom"){
      updateSelectInput(session, "sep", selected = "tab")
      updateSelectInput(session, "dec", selected = ".")
      updateTextInput(session, "skip", value = 0)
      updateNumericInput(session, "UV_col", value = 2)
      updateNumericInput(session, "fluo_col", value = NA)
    }
    val$type <- input$type
  })
  
  # Function for importing previous analysis from .RData file:
  import <- function(){
    load(input$input_data$datapath)
    
    for(i in names(forExport)){
      val[[i]] <- forExport[[i]]
    }
    
    updateSelectInput(session, "select",
                      choices = val$files_list, selected = val$files_list[[1]])
    updateSelectInput(session, "select_alignment",
                      choices = val$files_to_align)
    updateColourInput(session, "color", value = val$color_list[[input$select_alignment]])
    updateNumericInput(session, "linewidth", value = val$linewidth_collected[[input$select_alignment]]  )
    updateNumericInput(session, "smooth_spar", value = val$smooth_fluo )
    updateSelectInput(session, "linetype", 
                      selected = c("solid", "dashed", "dotted",
                                   "dotdash", "longdash", "twodash")[ val$linetype_collected[[input$select_alignment]] ])
    updateCheckboxInput(session, "red_lines", val = val$checkboxes[["red_lines"]])
    updateCheckboxInput(session, "red_lines2", val = val$checkboxes[["red_lines2"]])
    updateCheckboxInput(session, "green_lines", val = val$checkboxes[["green_lines"]])
    updateCheckboxInput(session, "show_fl", val = val$checkboxes[["show_fl"]])
    updateCheckboxInput(session, "show_fl_al", val = val$checkboxes[["show_fl_al"]])
    updateCheckboxInput(session, "normalize_length", val = val$checkboxes[["normalize_length"]])
    updateCheckboxInput(session, "normalize_height", val = val$checkboxes[["normalize_height"]])
    updateCheckboxInput(session, "anchor_line", val = val$checkboxes[["anchor_line"]])
    updateCheckboxInput(session, "helper_functions", val = val$checkboxes[["helper_functions"]])
    updateRadioButtons(session, "color_palette", selected = val$color_palette)
    updateSelectInput(session, "type", selected = val$type)
  }
      
  # Function for reading polysome profile data from individual files
  read_data <- function(){
    
    new_names <- input$input_data$name
    new_paths <- input$input_data$datapath
    
    for(i in 1:length(new_names)){
      this_name <- new_names[i]
      this_path <- new_paths[i]
      
      # How many lines need to be skipped from the beginning of the file?
      if(isTruthy(input$skip)){
        skip <- input$skip
      }else{
        skip <- 0
      }
      
      # Determine column separator
      sep <- c("", "\t", ",", ";")[ match(input$sep, c("space", "tab", ",", ";") ) ]
      
      # Import data 
      data <- read.delim(this_path, sep = sep, dec = input$dec, skip = skip)
      UV <- data[ , input$UV_col] 
      log(UV) # returns an error when the data is not numeric and triggers notification
      val$polysome_data[[this_name]] <- UV
      val$xvalues[[this_name]] <- 1:length(UV)
      
      if(isTruthy(input$fluo_col)){
        fluo <- data[ , input$fluo_col] 
        log(fluo) # returns an error when the data is not numeric and triggers notification
        val$fluo_data[[this_name]] <- fluo 
      }
    }
    
    val$paths_collected[input$input_data$name] <- input$input_data$datapath
    val$files_list <- c(input$input_data$name[input$input_data$size != 0], val$files_list)
    
    updateSelectInput(session, "select",
                      choices = val$files_list)
  }
  
  # Load data and show error message when loading failed:
  observeEvent(input$input_data, {
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
      if( any( class(loading_attempt) == "try-error") ){
        showNotification("File format not accepted.",
                         duration = NULL, type = "error")
      }
    }
  })
  
  ## USER INPUT FOR PLOTTING INDIVIDUAL PROFILES
  
  # Change smoothing of the fluorescence profile to spar value entered by user:
  observeEvent(input$smooth_spar, {
    val$smooth_fluo <- input$smooth_spar
  })
  
  # Smoothing of fluorescence profiles
  fluorescence <- reactive({
    req(input$select) 
    req(input$show_fl)
    if( isTruthy(val$fluo_data[[input$select]]) ){
      smooth_profile( val$fluo_data[[input$select]], val$smooth_fluo )
    }
  })
  
  # Store all checkBox inputs to restore them when a previous analysis is loaded as .RData:
  observeEvent(input$red_lines, {
    val$checkboxes[["red_lines"]] <- input$red_lines
  })
  observeEvent(input$red_lines2, {
    val$checkboxes[["red_lines2"]] <- input$red_lines2
  })
  observeEvent(input$green_lines, {
    val$checkboxes[["green_lines"]] <- input$green_lines
  })
  observeEvent(input$show_fl, {
    val$checkboxes[["show_fl"]] <- input$show_fl
  })
  observeEvent(input$show_fl_al, {
    val$checkboxes[["show_fl_al"]] <- input$show_fl_al
  })
  observeEvent(input$normalize_length, {
    val$checkboxes[["normalize_length"]] <- input$normalize_length
  })
  observeEvent(input$normalize_height, {
    val$checkboxes[["normalize_height"]] <- input$normalize_height
  })
  observeEvent(input$anchor_line, {
    val$checkboxes[["anchor_line"]] <- input$anchor_line
  })
  observeEvent(input$helper_functions, {
    val$checkboxes[["helper_functions"]] <- input$helper_functions
  })
  observeEvent(input$color_palette, {
    val$color_palette <- input$color_palette
  })
  
  # Create axis limits for plot of individual profiles:
  # By default, the minimum and maximum of the profile data is used as limits.
  # If the user has selected other limits for a profile, they are used instead: 
  ymin_single <- reactive({
    req(input$select)
    if(isTruthy(val$ymin_collected[[input$select]])){
      val$ymin_collected[[input$select]]
    }else{
      floor( min( val$polysome_data[[input$select]] ) )
    }})
  
  ymin_single_fl <- reactive({
    req(input$select)
    if(isTruthy(val$ymin_collected_fl[[input$select]])){
      val$ymin_collected_fl[[input$select]]
    }else{
      floor( min( fluorescence() ) )
    }})
  
  ymax_single <- reactive({
    req(input$select)
    if(isTruthy(val$ymax_collected[[input$select]])){
      val$ymax_collected[[input$select]]
    }else{
      ceiling( max( val$polysome_data[[input$select]] ) )
    }})
  
  ymax_single_fl <- reactive({
    req(input$select)
    if(isTruthy(val$ymax_collected_fl[[input$select]])){
      val$ymax_collected_fl[[input$select]]
    }else{
      ceiling( max( fluorescence() ) )
    }})
  
  xmin_single <- reactive({
    req(input$select)
    if(isTruthy(val$xmin_collected[[input$select]])){
      val$xmin_collected[[input$select]]
    }else{
      floor( min( val$xvalues[[input$select]] ) )
    }})
  
  xmax_single <- reactive({
    req(input$select)
    if(isTruthy(val$xmax_collected[[input$select]])){
      val$xmax_collected[[input$select]]
    }else{
      ceiling( max( val$xvalues[[input$select]] ) )
    }})
  
  # The axis limits are then shown in the respective input fields.
  observeEvent( input$select, {
    updateNumericInput(session, "axis1", value = ymin_single(), step = lost_num_pol() )
  })
  
  observeEvent( input$select, {
    updateNumericInput(session, "axis2", value = ymax_single(), step = lost_num_pol() )
  })
  
  observeEvent( input$select, {
    updateNumericInput(session, "axis3", value = xmin_single(), step = lost_num_x() )
  })
  
  observeEvent( input$select, {
    updateNumericInput(session, "axis4", value = xmax_single(), step = lost_num_x() )
  })
  
  observeEvent( input$select, {
    if(input$show_fl & isTruthy(val$fluo_data[[input$select]]) ){
      updateNumericInput(session, "axis1_fl", value = ymin_single_fl(), step = lost_num_fl() )
    }else{
      updateNumericInput(session, "axis1_fl", value = NA )
    }
  })
  
  observeEvent( input$select, {
    if(input$show_fl & isTruthy(val$fluo_data[[input$select]]) ){
      updateNumericInput(session, "axis2_fl", value = ymax_single_fl(), step = lost_num_fl() )
    }else{
      updateNumericInput(session, "axis2_fl", value = NA )
    }
  })
  
  observeEvent( input$show_fl, {
    if(input$show_fl & isTruthy(val$fluo_data[[input$select]]) ){
      updateNumericInput(session, "axis1_fl", value = ymin_single_fl(), step = lost_num_fl() )
    }else{
      updateNumericInput(session, "axis1_fl", value = NA )
    }
  })
  
  observeEvent( input$show_fl, {
    if(input$show_fl & isTruthy(val$fluo_data[[input$select]]) ){
      updateNumericInput(session, "axis2_fl", value = ymax_single_fl(), step = lost_num_fl() )
    }else{
      updateNumericInput(session, "axis2_fl", value = NA )
    }
  })
  
  # The functions lost_num_fl(), lost_num_pol() and lost_num_x() determine how many numerals 
  # should be omitted from the tick mark labels for readability
  # and return a factor by which the original values have to be divided to obtain the labels. 
  lost_num_x <- reactive({
    # how many numerals has the maximum x-value?
    max_numeral <- max( floor(log10(abs( xmax_single() ))) + 1 ) 
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_fl <- reactive({
    # how many numerals has the maximum fluorescence signal?
    max_numeral <- max( floor( log10( abs( ymax_single_fl() ) ) ) + 1 ) 
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_pol <- reactive({
    # how many numerals has the maximum UV absorption?
    max_numeral <- max( floor(log10(abs( ymax_single() ))) + 1 ) 
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When a profile is selected that does not contain fluorescence data,
  # but fluorescence was shown for the previously selected profile,
  # the check box "Show" for fluorescence is automatically de-selected.
  observeEvent(input$select,{
    if( input$show_fl & !isTruthy(val$fluo_data[[input$select]]) ){
      updateCheckboxInput(session, "show_fl", value = FALSE)
    }
  })
  
  ## PLOT INDIVIDUAL PROFILES
  
  # The function plot_singleFl() is used to plot individual fluorescence profiles.
  plot_singleFl <-function(cex_lab, cex_axis, lwd, mgp2, lwd_box){
    if(lost_num_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_fl(), ")", sep = "")
    }
    plot(val$xvalues[[input$select]], fluorescence(), type = "l", xaxt = "n", col = "black", lwd = lwd,
         xlim = c( xmin_single(), xmax_single() ), ylim = c( ymin_single_fl(), ymax_single_fl() ),
         las = 1, ylab = ylab, xlab = "", yaxt = "n", cex.lab = cex_lab, mgp = mgp2, bty = "n")
    box(lwd = lwd_box)
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_fl(), las = 1, mgp = mgp2, cex.axis = cex_axis, lwd = lwd_box)
    
    # If an area has been quantified by the user in the polysome profile, 
    # also show the area as green polygon also in the fluorescence profile.
    # If baseline, area_end or area_start lines are changed, 
    # the colored area stays the same for the selected quantified area 
    # until the button "quantify area" is pressed again.
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) & isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      # Rounding of the start and end point might be necessary,
      # in case they are located between two data points
      # (because a local minimum or maximum might be in the middle between two data points).
      x_first <- round(val$area_starts[[input$select_area]][input$select])
      x_last <- round(val$area_ends[[input$select_area]][input$select])
      
      if(isTruthy(val$fl_control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, val$xvalues[[input$select]][ x_first : x_last ], x_last),
                c(val$fl_control_baseline[[input$select_area]][input$select], fluorescence()[ x_first : x_last ], val$fl_control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "black", lwd = 2)
      }
    }
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]], col = "#238b45", lty = 2, lwd = lwd)
      abline(v = val$file_ends[[input$select]], col = "#238b45", lty = 2, lwd = lwd)
    }
    if(input$red_lines2){
      abline(h = val$baseline_fl[input$select], col = "red", lty = 2, lwd = lwd)
    }
  }
  
  # The function plot_singlePol() plots individual polysome profiles.
  plot_singlePol <-function(cex_lab, cex_axis, lwd, mgp1, mgp2, lwd_box){
    req(val$xvalues[[input$select]])
    req(val$polysome_data[[input$select]])
    
    if(lost_num_x() == 1 ){
      xlab <- "Data Points"
    }else{
      xlab <- paste("Data points (x ", lost_num_x(), ")", sep = "")
    }
    
    if(lost_num_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_pol(), ")", sep = "")
    }
    
    plot(val$xvalues[[input$select]], val$polysome_data[[input$select]], type = "l", las = 1, lwd = lwd,
         ylab = ylab, xlab = "", mgp = mgp2,
         ylim = c(ymin_single(),ymax_single()), 
         xlim =c(xmin_single(),xmax_single()),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab, bty = "n"
    )
    box(lwd = lwd_box)
    # The axis tick marks are corrected by the factor determined with the functions lost_num_pol()
    # and lost_num_x() to prevent large tick mark labels:
    a <- axTicks(1)
    axis(1, at = a, labels = a/lost_num_x(), las = 1, mgp = mgp1, cex.axis = cex_axis, lwd = lwd_box)
    mtext(xlab, side = 1, line = mgp1[1], cex = cex_lab)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = mgp2, cex.axis = cex_axis, lwd = lwd_box)
    
    # If an area has been quantified by the user in the polysome profile, 
    # it is shown as a green polygon.
    # If baseline, area_end or area_start lines are changed, 
    # the colored area stays the same for the selected quantified area 
    # until the button "quantify area" is pressed again.
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) & isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      
      # Rounding of the start and end point might be necessary,
      # in case they are located between two data points
      # (because a local minimum or maximum might be in the middle between two data points).
      x_first <- round(val$area_starts[[input$select_area]][input$select])
      x_last <- round(val$area_ends[[input$select_area]][input$select])
      
      if(isTruthy(val$control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, val$xvalues[[input$select]][ x_first : x_last ], x_last),
                c(val$control_baseline[[input$select_area]][input$select], val$polysome_data[[input$select]][x_first : x_last], val$control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "black", lwd = lwd)
      }
    }
    # selected x-anchor and baseline are displayed if box is ticked
    if(input$red_lines){
      abline(v = val$anchor[input$select], col = "red", lty = 2, lwd = lwd)
    }
    
    if(input$red_lines2){
      abline(h = val$baseline[input$select], col = "red", lty = 2, lwd = lwd)
    }
    
    if(input$green_lines){
      abline(v = val$file_starts[[input$select]], col = "#238b45", lty=2, lwd = lwd)
      abline(v = val$file_ends[[input$select]], col = "#238b45", lty=2, lwd = lwd)
    }
  }
  
  # The function plot_singleInput() divides the plotting area into two parts 
  # when there is a fluorescence profile to show and adjusts the margins accordingly.
  plot_singleInput <- function(cex_lab, cex_axis, lwd, mgp1, mgp2, mar_factor, lwd_box){
    if(!is.null(val$xvalues[[input$select]])){
      if(input$show_fl & isTruthy(val$fluo_data[[input$select]])  ){
        if(val$buttons == 5){
          layout(matrix(2:1, 2, 1), height = c(0.6, 0.9) ) # divides the plotting area into 2 rows
          par(mar = c(5, 5, 0, 0.1)*mar_factor)
          plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2, lwd_box)
          par(mar = c(0, 5, 0.5, 0.1)*mar_factor)
          plot_singleFl(cex_lab, cex_axis, lwd, mgp2, lwd_box)
        }else{
          layout(matrix(1:2, 2, 1), height = c(0.6, 0.9) ) # divides the plotting area into 2 rows
          par(mar = c(0, 5, 0.5, 0.1)*mar_factor)
          plot_singleFl(cex_lab, cex_axis, lwd, mgp2, lwd_box)
          par(mar = c(5, 5, 0, 0.1)*mar_factor)
          plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2, lwd_box)
        }
      }else{
        par(mar = c(5, 5, 0.5, 0.1)*mar_factor)
        plot_singlePol(cex_lab, cex_axis, lwd, mgp1, mgp2, lwd_box)
      }
    }
  }
  
  output$plot_single <- renderPlot({
    req(val$xvalues[[input$select]])
    plot_singleInput(cex_lab = 1.6, cex_axis = 1.6, lwd = 2, 
                     mgp1 = c(3.5, 1.2, 0), mgp2 = c(3.5, 0.8, 0), 
                     mar_factor = 1, lwd_box = 1)
  })

  ## INTERACTION OF THE USER WITH THE PROFILE PLOT
  
  # The value of val$buttons changes by clicking on different action buttons that trigger the respective action(s).
  # The initial value of val$buttons and value after selecting a new file is 10 (which stands fo "no action selected")
  val <- reactiveValues(
    buttons = 10)
  observeEvent(input$select,{
    val$buttons = 10
  })
  
  # Specific values are set for val$buttons encoding for the current action initiated by the user
  observeEvent(input$x_anchor, {val$buttons = 1})
  observeEvent(input$baseline, {val$buttons = 2})
  observeEvent(input$file_start, {val$buttons = 3})
  observeEvent(input$file_end, {val$buttons = 4})
  observeEvent(input$baseline_fl, {val$buttons = 5})
  
  # By clicking into the plot, the user can set a baseline (val$buttons = 2) and an x-anchor (val$buttons = 1) for alignment,
  # or start (val$buttons = 3) and end (val$buttons = 4) of the area to quantify.
  # Also for fluorescence profiles, a baseline needs to be set for alignment and quantification (val$buttons = 5)
  # The selected values are stored for each profile.
  # A function to identify the closest local minimum or maximum is used by default,
  # but can be inactivated (tick box input$peak_help).
  observeEvent(input$click$x, {
    req(val$buttons)
    
    # Setting start of area to quantify
    if(val$buttons == 3)
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
    
    # Setting ends of area to quantify
    if(val$buttons == 4)
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
    
    # Setting baseline for polysome profile
    if(val$buttons == 2){
      val$baseline[input$select] <- val$polysome_data[[input$select]][input$click$x]
    }
    
    # Setting baseline for fluorescence profile
    if(val$buttons == 5){
      par(xpd = T)
      val$baseline_fl[input$select] <- fluorescence()[input$click$x]
    }
    
    # Setting x-anchor
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
    
    # If both x-anchor and baseline exist for a polysome profile, 
    # the name of the file is added to the select_alignment input.
    val$files_to_align <- intersect( names(val$baseline), names(val$anchor) )
    updateSelectInput(session, "select_alignment",
                      choices = val$files_to_align)
  })
  
  # Calculate sum of values to quantify an area when the quantify button is pressed 
  # and all required values (start, end and baseline) were set before
  observeEvent(input$quantify_area,{
    req(val$file_starts[[input$select]], val$file_ends[[input$select]], val$baseline[input$select])
    # Store profile baselines for quantified areas to be correctly displayed in the plotting area 
    # even when the baseline is changed after the quantification by the user
    if(isTruthy(val$baseline_fl[input$select])){
      val$fl_control_baseline[[input$select_area]][input$select] <- val$baseline_fl[input$select]
    }
    val$control_baseline[[input$select_area]][input$select] <- val$baseline[input$select]
    
    # Rounding of the start and end values is necessary, because they might lie between two points
    x_first <- round(val$file_starts[[input$select]])
    x_last <- round(val$file_ends[[input$select]])
    
    # Estimate the surface of the selected area by summing up
    # all values from start to end area after subtracting the baseline.
    val$sum_areas[[input$select_area]][input$select] <-
      sum( val$polysome_data[[input$select]][x_first:x_last] - val$baseline[input$select] )
    
    if(input$show_fl){
      val$sum_areas[[paste(input$select_area, "_fluo", sep = "")]][input$select] <- 
        sum( fluorescence()[x_first:x_last] - val$baseline_fl[input$select] )
    }
    # store not only area but also start and stop values in the same ways as areas!
    val$area_starts[[input$select_area]][input$select] <- val$file_starts[[input$select]]
    val$area_ends[[input$select_area]][input$select] <- val$file_ends[[input$select]]
  })
  
  # Create a data frame with quantified areas,
  # filling in NAs when the respective area is missing for a file.
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
  
  # Save axis limits set by the user for individual polysome profiles
  observeEvent(input$axis1, {
    val$ymin_collected[[input$select]] <- input$axis1
  })
  observeEvent(input$axis2, {
    val$ymax_collected[[input$select]] <- input$axis2
  })
  observeEvent(input$axis3, {
    val$xmin_collected[[input$select]] <- input$axis3
  })
  observeEvent(input$axis4, {
    val$xmax_collected[[input$select]] <- input$axis4
  })
  
  # Save axis limits set by the user for individual fluorescence profiles
  observeEvent(input$axis1_fl, {
    if(!is.na(input$axis1_fl) ){
      val$ymin_collected_fl[[input$select]] <- input$axis1_fl 
    }
  })
  
  observeEvent(input$axis2_fl, {
    if(!is.na(input$axis2_fl) ){
      val$ymax_collected_fl[[input$select]] <- input$axis2_fl
    }
  })
 
  ## GENERATE ALIGNMENT
  
  # For alignment of polysome profiles along the y-axis, the baseline is subtracted.
  # If the user selects "Normalize area", the result is divided by a normalization 
  # factor to obtain equal total areas for all profiles. 
  # If the user selects "Normalize length", it is additionally necessary to multiply
  # by the normalization factor of the x-values to ensure that the total area is not affected
  # by the length normalization along the the x-axis. 
  values_list <- reactive({
    dummy <- list()
    for(f in val$files_to_plot){
      next_y <- val$polysome_data[[f]]
      dummy[[f]] <- ( ( next_y - val$baseline[f])/norm_factor()[f])*norm_factor_x()[f]
    }
    dummy
  })
  
  # For alignment of fluorescence profiles along the y-axis, the baseline is subtracted.
  # If the user selects "Normalize area" or "Normalize length", 
  # the fluorescence profiles are adjusted by the same normalization factor as the polysome profiles.
  # In addition, the profiles are smoothed using the function smooth_profile.
  values_fluorescence <- reactive({
      dummy <- list()
      for(f in intersect(names(val$baseline_fl), val$files_to_plot ) ){ # Only go through fluorescence signals with baseline
        fluo <- val$fluo_data[[f]]
        dummy[[f]] <- smooth_profile( ( ( fluo - val$baseline_fl[f] )/norm_factor()[f])*norm_factor_x()[f], val$smooth_fluo)
      }
      dummy
  })
  
  # The x-values of each profile are generated as sequence from the first to the last
  # x-value in steps of one, if there is no length normalization, and in steps of 
  # one divided by the normalization factor, if the user selected "Normalize length".
  x_values <- reactive({
    dummy <- list()
    for(f in val$files_to_plot){
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      dummy[[f]] <- x
    }
    dummy
  })
  
  # For alignment of profiles along the x-axis, the shift has to be determined for all
  # profiles from the position of the x-anchors that have been set by the user.
  # The x-anchor with the highest value is used as a reference.
  # If length normalization is selected by the user, 
  # the position of the x-anchors is adjusted accordingly.
  shifts <- reactive({
    all_anchors <- val$anchor[val$files_to_plot] / norm_factor_x()[val$files_to_plot]
    super_anchor <- max(all_anchors)
    dummy <- super_anchor - all_anchors
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # The function aligned_starts() determines the first x-value of each profile,
  # taking the shift for alignment and, if selected by the user, 
  # the length normalization into account.
  aligned_starts <- reactive({
    dummy <- 1/norm_factor_x()[val$files_to_plot] + shifts()  
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # The function aligned_ends() determines the last x-value of each profile,
  # taking the shift for alignment and, if selected by the user, 
  # the length normalization into account.
  aligned_ends <- reactive({
    profile_lengths <- sapply(values_list(), FUN = length)
    names(profile_lengths) <- val$files_to_plot
    dummy <- profile_lengths/norm_factor_x()[val$files_to_plot] + shifts()
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # The function norm_factor() determines a normalization factor 
  # to ensure that all polysome profiles have the same total area.
  # As long as "Normalize area" is not selected by the user, the factor is set to 1
  # for all profiles and no normalization is performed.
  # Otherwise, the profile with the highest total area is used as reference.
  norm_factor <- reactive({
    files_with_total <- names( val$sum_areas[["Total"]] )
    if( all(val$files_to_plot %in% files_with_total) & input$normalize_height ){
      ref_area <- max( val$sum_areas[["Total"]][val$files_to_plot] )
      dummy <- val$sum_areas[["Total"]][val$files_to_plot]/ref_area
    }else{
      dummy <- rep( 1, length(val$files_to_plot) )
    }
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # The function norm_factor_x() determines a normalization factor 
  # to ensure that all polysome profiles have the same length of their total area.
  # As long as "Normalize area" is not selected by the user, the factor is set to 1
  # for all profiles and no normalization is performed.
  # Otherwise, the profile with the longest total area is used as reference.
  norm_factor_x <- reactive({
    files_with_total <- names( val$sum_areas[["Total"]] )
    if( all(val$files_to_plot %in% files_with_total) & input$normalize_length ){
      all_lengths <- round(val$area_ends[["Total"]][val$files_to_plot] - val$area_starts[["Total"]][val$files_to_plot] )
      ref_length <- max(all_lengths)
      dummy <- round(val$area_ends[["Total"]][val$files_to_plot] - val$area_starts[["Total"]][val$files_to_plot] ) / ref_length
    }else{
      dummy <- rep( 1, length(val$files_to_plot) )
    }
    names(dummy) <- val$files_to_plot
    dummy
  })
  
  # For the tab "Alignment table", a data frame is generated 
  # that shows the aligned profile values with a common set of x-values. 
  # Empty cells are filled with NAs.
  create_df_pol_alignment <- reactive({
    y_aligned <- values_list()
    x_aligned <- x_values()
    
    y <- y_aligned[[ val$files_to_plot[1] ]]
    x <- x_aligned[[ val$files_to_plot[1] ]]
    df <- data.frame(x, y)
    name <- str_remove(val$files_to_plot[1], ".pks|.csv|.txt|.tsv")
    colnames(df) <- c("Position", name)
    
    if( length(val$files_to_plot) > 1){
      for( f in val$files_to_plot[-1] ){
        y <- y_aligned[[f]]
        x <- x_aligned[[f]]
        df_new <- data.frame(x, y)
        name <- str_remove(f, ".pks|.csv|.txt|.tsv")
        colnames(df_new) <- c("Position", name)
        df <- merge(df, df_new, by = "Position", all = T)
      }
    }
    val$csv_file_df <- df
  })
  
  # For the tab "Alignment table", a data frame is generated 
  # that shows the aligned fluorescence values with a common set of x-values. 
  # Empty cells are filled with NAs.
  create_df_fluo_alignment <- reactive({
    y_aligned <- values_fluorescence()
    x_aligned <- x_values()
    
    fluo_files <- intersect( names( values_fluorescence() ), val$files_to_plot)
    y <- y_aligned[[ fluo_files[1] ]]
    x <- x_aligned[[ fluo_files[1] ]]
    df <- data.frame(x, y)
    name <- str_remove(fluo_files[1], ".pks|.csv|.txt|.tsv")
    colnames(df) <- c("Position", paste(name, "_fluo") )
    
    if(length(fluo_files) > 1){
      for( f in fluo_files[-1] ){
        y <- y_aligned[[f]]
        x <- x_aligned[[f]]
        df_new <- data.frame(x, y)
        name <- str_remove(f, ".pks|.csv|.txt|.tsv")
        colnames(df_new) <- c("Position", paste(name, "_fluo") )
        df <- merge(df, df_new, by = "Position", all = T)
      }
    }
    val$csv_file_df_fluo <- df
  })
  
  # The data frames with aligned polysome profile data 
  # and aligned fluorescence profile data have to be merged:
  merge_df_alignment <- reactive({
    req(val$csv_file_df)
    req(val$csv_file_df_fluo)
    val$csv_file_df_all <- merge(val$csv_file_df, val$csv_file_df_fluo, by = "Position", all = T)
  })
  
  ## DISPLAY AND MODIFY ALIGNMENT
  
  # Determine default x-axis limits for the alignment plot 
  # from the minima and maxima of the aligned data.
  
  xmax_all <- reactive({
    ceiling( max(aligned_ends()) )
  })
  
  xmin_all <- reactive({
    floor( min(aligned_starts()) )
  })
  
  # Determine default y-axis limits for the alignment plot 
  # from the minima and maxima of the aligned polysome profile data.
  ymax_all <- reactive({
    profile_heights <- sapply(values_list()[val$files_to_plot], FUN = max)
    ceiling( max( profile_heights ) )
  })
  
  ymin_all <- reactive({
    profile_mins <- sapply(values_list()[val$files_to_plot], FUN = min)
    floor( min( profile_mins ) )
  })
  
  # Determine default y-axis limits for the alignment plot 
  # from the minima and maxima of the aligned fluorescence profile data.
  ymax_all_fl <- reactive({
    profile_heights <- sapply(values_fluorescence()[val$files_to_plot], FUN = max)
    ceiling( max( unlist( profile_heights ) ) )
  })
  
  ymin_all_fl <- reactive({
    profile_mins <- sapply(values_fluorescence()[val$files_to_plot], FUN = min)
    floor( min( unlist( profile_mins ) ) )
  })
  
  # The functions lost_num_al_fl(), lost_num_al_pol() and lost_num_al_Index() determine how many numerals 
  # should be omitted from the tick mark labels for readability
  # and return a factor by which the original values have to be divided to obtain the labels. 
  lost_num_al_fl <- reactive({
    max_numeral <- floor( log10( abs( ymax_aligned_fl() ) ) + 1 )
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_al_pol <- reactive({
    max_numeral <- floor( log10( abs( ymax_aligned() ) ) + 1 )
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  lost_num_al_Index <- reactive({
    max_numeral <- floor(log10(abs( xmax_aligned() ) ) + 1 )
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When the user selects other axis limits, the limits are adapted.
  observeEvent(input$axis1_a, {
    val$ymin <- input$axis1_a
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis2_a, {
    val$ymax <- input$axis2_a
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis1_a_fl, {
    val$ymin_fl <- input$axis1_a_fl
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis2_a_fl, {
    val$ymax_fl <- input$axis2_a_fl
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis3_a, {
    val$xmin <- input$axis3_a
  }, 
  ignoreInit = T)
  
  observeEvent(input$axis4_a, {
    val$xmax <- input$axis4_a
  }, 
  ignoreInit = T)
  
  # Axis limits are either determined automatically by the functions 
  # xmin_all() and xmax_all() for the x-values,
  # ymin_all() and ymax_all() for the polysome profile y-values and
  # ymin_all_fl() and ymax_all_fl() for the fluorescence profile y-values.
  ymin_aligned <- reactive({
    if(isTruthy(val$ymin) ){
      val$ymin
    }else{
      ymin_all()
    }
  })
  
  ymin_aligned_fl <- reactive({
    if(isTruthy(val$ymin_fl)){
      val$ymin_fl
    }else{
      ymin_all_fl()
    }
  })
  
  ymax_aligned <- reactive({
    if(isTruthy(val$ymax)){
      val$ymax
    }else{
      ymax_all()
    }
  })
  
  ymax_aligned_fl <- reactive({
    if(isTruthy(val$ymax_fl)){
      val$ymax_fl
    }else{
      ymax_all_fl()
    }
  })
  
  xmin_aligned <- reactive({
    if(isTruthy(val$xmin)){
      val$xmin
    }else{
      xmin_all()
    }
  })
  
  xmax_aligned <- reactive({
    if(isTruthy(val$xmax)){
      val$xmax
    }else{
      xmax_all()
    }
  })
  
  # When the set of aligned profiles is changed,
  # the minima or maxima are automatically set as axis limits
  # and are displayed in the respective input field.
  observeEvent( val$files_to_plot, {
    updateNumericInput(session, "axis1_a", value = ymin_all(), step = lost_num_al_pol() )
  })
  
  observeEvent( val$files_to_plot, {
    updateNumericInput(session, "axis2_a", value = ymax_all(), step = lost_num_al_pol() )
  })
  
  observeEvent( val$files_to_plot, {
    if(input$show_fl_al & isTruthy( values_fluorescence()[val$files_to_plot] ) ){
      updateNumericInput(session, "axis1_a_fl", value = ymin_all_fl(), step = lost_num_al_fl() )
    }
  })
  
  observeEvent( val$files_to_plot, {
    if(input$show_fl_al & isTruthy( values_fluorescence()[val$files_to_plot] ) ){
    updateNumericInput(session, "axis2_a_fl", value = ymax_all_fl(), step = lost_num_al_fl() )
    }
  })
  
  observeEvent( val$files_to_plot, {
    req( xmin_aligned() )
    updateNumericInput(session, "axis3_a", value = xmin_all(), step = lost_num_al_Index() )
  })
  
  observeEvent( val$files_to_plot, {
    req( xmax_aligned() )
    updateNumericInput(session, "axis4_a", value = xmax_all(), step = lost_num_al_Index() )
  })
  
  # When "Show" is selected for aligned fluorescence profiles,
  # the minima or maxima are automatically set as axis limits
  # and are displayed in the respective input field.
  observeEvent( input$show_fl_al, {
    if(input$show_fl_al & isTruthy( values_fluorescence()[val$files_to_plot] ) ){
      updateNumericInput(session, "axis1_a_fl", value = ymin_all_fl(), step = lost_num_al_fl() )
    }
  })
  
  observeEvent( input$show_fl_al, {
    if(input$show_fl_al & isTruthy( values_fluorescence()[val$files_to_plot] ) ){
      updateNumericInput(session, "axis2_a_fl", value = ymax_all_fl(), step = lost_num_al_fl() )
    }
  })
 
  # When a new profile is available for alignment, 
  # because both baseline and x-anchor have been set by the user,
  # it is automatically shown in the alignment.
  observeEvent(val$files_to_align, {
    val$files_to_plot <- setdiff( val$files_to_align, val$files_to_exclude) 
  })
  
  # When a new profile is available for alignment, 
  # the line type is set to 1.
  observeEvent(val$files_to_align, {
    newest <- tail(val$files_to_align, 1)
    # To prevent overwriting of the last profile linetype
    # when a previous analysis is imported as .RData:
    if( !isTruthy( val$linetype_collected[[newest]] ) ){
      val$linetype_collected[[ newest ]] <- 1
    }
  })
  
  # When a new profile is available for alignment, 
  # the line width is set to 1.
  observeEvent(val$files_to_align, {
    newest <- tail(val$files_to_align, 1)
    # To prevent overwriting of the last profile linetype
    # when a previous analysis is imported as .RData
    if( !isTruthy( val$linewidth_collected[[newest]] ) ){
      val$linewidth_collected[[newest]] <- 2
    }
  })
  
  # The user can change the line type of a selected profile in the alignment.
  observeEvent(input$linetype, {
    req(input$select_alignment)
    if(input$linetype == "solid") line_selected <- 1
    if(input$linetype == "dashed") line_selected <- 2
    if(input$linetype == "dotted") line_selected <- 3
    if(input$linetype == "dotdash") line_selected <- 4
    if(input$linetype == "longdash") line_selected <- 5
    if(input$linetype == "twodash") line_selected <- 6
    val$linetype_collected[[input$select_alignment]] <- line_selected
  })
  
  # The user can change the line width of a selected profile in the alignment.
  observeEvent(input$linewidth, {
    req(input$select_alignment)
    val$linewidth_collected[[input$select_alignment]] <- input$linewidth
  })
  
  # Initially, the rainbow palette is applied to the aligned profiles.
  # Alternatively, the palette "Dark 3" or a color blind friendly palette can be selected.
  # If colors of individual profiles are changed via input$color,
  # they replace the respective color in the palette.
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
    
    # The names of the color vector correspond to the SORTED file names,
    # because the order of the plotted files can change, while the assigned
    # colors should remain the same irrespective of the order.
    names(dummy) <- sort(val$files_to_plot) 
    # If individual colors have been set already by the user via input$color:
    if (length(val$colors_collected) > 0){
      dummy[names(val$colors_collected) ] <- unlist(val$colors_collected)
    }
    val$color_list <- as.list(dummy)
  })
  
  # If the selected color palette is changed by the user, 
  # the individually selected colors are reset to an empty list.
  observeEvent(input$color_palette, {
    val$colors_collected = list()
  }, 
  ignoreInit = TRUE
  )
  
  # The user can change the color of a selected profile in the alignment.
  observeEvent(input$color, {
    req(input$select_alignment)
    val$colors_collected[[input$select_alignment]] <- input$color
  })
  
  # When a profile is selected among the aligned profiles,
  # the color, line type and line width are updated to the values of the selected profile.
  # Caution: updateColourInput seems to require a list and double square brackets [[]]!
  observe({
    req(input$select_alignment)
    updateColourInput(session, "color", value = val$color_list[[input$select_alignment]])
    updateNumericInput(session, "linewidth", value = val$linewidth_collected[[input$select_alignment]]  )
    updateSelectInput(session, "linetype", 
                      selected = c("solid", "dashed", "dotted",
                                   "dotdash", "longdash", "twodash")[ val$linetype_collected[[input$select_alignment]] ])
  })
  
  # Individual profiles can be excluded from the alignment.
  observeEvent(input$show_in_al, {
    if(input$show_in_al){
      val$files_to_exclude <- setdiff(val$files_to_exclude, input$select_alignment)
    }else{
      val$files_to_exclude <- c(val$files_to_exclude, input$select_alignment)
    }
    val$files_to_plot <- setdiff( val$files_to_align, val$files_to_exclude) 
  })
  
  # The tick mark for showing a profile in the alignment is updated
  # when a profile is selected among the aligned profiles.
  observeEvent(input$select_alignment, {
    value <- input$select_alignment %in% val$files_to_plot
    updateCheckboxInput(session, "show_in_al", value = value)
  })
  
  # The order of aligned profiles can be changed by the user.
  # If a profile is moved further up, it will swap positions
  # with the profile on the position above.
  observeEvent(input$up, {
    if(input$select_alignment %in% val$files_to_plot){
      to_be_shifted <- which( val$files_to_plot == input$select_alignment )
      if(to_be_shifted > 1){
        files_order <- 1:length(val$files_to_plot)
        files_order[to_be_shifted] <- files_order[to_be_shifted] - 1
        files_order[to_be_shifted - 1] <- files_order[to_be_shifted - 1] + 1
        val$files_to_plot <- val$files_to_plot[files_order]
      }
    }
  })
  
  # If a profile is moved further down, it will swap positions
  # with the profile on the position below.
  observeEvent(input$down, {
    if(input$select_alignment %in% val$files_to_plot){
      to_be_shifted <- which( val$files_to_plot == input$select_alignment )
      if(to_be_shifted < length(val$files_to_plot) ){
        files_order <- 1:length(val$files_to_plot)
        files_order[to_be_shifted] <- files_order[to_be_shifted] + 1
        files_order[to_be_shifted + 1] <- files_order[to_be_shifted + 1] - 1
        val$files_to_plot <- val$files_to_plot[files_order]
      }
    }
  })
  
  # When a new file is added to the alignment, 
  # the file name is used for the figure legend.
  observeEvent(val$files_to_align, {
    newest <- tail(val$files_to_align, 1)
    # To prevent overwriting of the last profile name
    # when a previous analysis is imported as .RData:
    if( !isTruthy( val$legend_names[newest] ) ){
      basename <- gsub(".csv|.txt|.pks|.tsv", "", newest)
      val$legend_names[newest] <- basename
    }
  })
  
  # The user can change the profile names in the legend of the alignment.
  observeEvent(input$rename, {
    val$legend_names[input$select_alignment] <- input$new_name 
  })
  
  ## PLOT ALIGNMENTS
  
  # This function plots the aligned polysome profiles,
  # and takes a set of graphical parameters as input.
  
  plot_alignedPol <- function(cex_lab, cex_axis, lwd, 
                              mgp1, mgp2, 
                              cex_legend, lwd_factor, lwd_box){
    if(lost_num_al_pol() == 1 ){
      ylab <- "UV abs."
    }else{
      ylab <- paste("UV abs. (x ", lost_num_al_pol(), ")", sep = "")
    }
    
    if(lost_num_al_Index() == 1 ){
      xlab <- "Relative position"
    }else{
      xlab <- paste("Relative position (x ", lost_num_al_Index(), ")", sep = "")
    }
    
    f <- val$files_to_plot[1]
    
    plot(x_values()[[f]], values_list()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]]*lwd_factor, 
         ylab = ylab, xlab = "", las = 1,
         col = val$color_list[[f]], mgp = mgp2, 
         ylim = c( ymin_aligned(), ymax_aligned() ), xlim = c(xmin_aligned(), xmax_aligned() ),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab, bty = "n"
    )
    box(lwd = lwd_box)
    # Shows the common anchor (when "Show x-anchor" is selected)
    if(input$anchor_line == TRUE){
      anchor_line <- val$anchor[val$files_to_plot[1]] / norm_factor_x()[val$files_to_plot[1]] + shifts()[val$files_to_plot[1]]
      abline(v = anchor_line, col = "grey", lty = 2, lwd = lwd)
    }
    
    a <- axTicks(1)
    axis(1, at = a, labels = a/lost_num_al_Index(), las = 1, mgp = mgp1, cex.axis = cex_axis, lwd = lwd_box)
    mtext(xlab, side = 1, line = mgp1[1], cex = cex_lab)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_pol(), las = 1, mgp = mgp2, cex.axis = cex_axis, lwd = lwd_box)
    
    for(f in val$files_to_plot[-1]){
      points(x_values()[[f]], values_list()[[f]], type = "l", lty = val$linetype_collected[[f]], 
             lwd = val$linewidth_collected[[f]]*lwd_factor, col = val$color_list[[f]])
    }
  }
  
  # This function plots the aligned fluorescence profiles,
  # and takes a set of graphical parameters as input.
  plot_alignedFluo <- function(cex_lab, cex_axis, lwd, 
                               mgp2, lwd_factor, lwd_box){
    if(lost_num_al_fl() == 1 ){
      ylab <- "Fluo."
    }else{
      ylab <- paste("Fluo. (x ", lost_num_al_fl(), ")", sep = "")
    }
    
    files_in_al_fluo <- val$files_to_plot[val$files_to_plot %in% intersect(names(val$baseline_fl),  names(val$baseline))] # This is necessary to maintain the order as given by val$files_to_plot.
    f <-  files_in_al_fluo[1]
    
    plot(x_values()[[f]], values_fluorescence()[[f]], type = "l", lty = val$linetype_collected[[f]],
         lwd = val$linewidth_collected[[f]]*lwd_factor,
         ylab = ylab, xlab = "", las = 1,
         col = val$color_list[[f]], mgp = mgp2, 
         ylim = c(ymin_aligned_fl(), ymax_aligned_fl() ), xlim = c(xmin_aligned(), xmax_aligned() ),
         yaxt = "n", xaxt = "n", cex.lab = cex_lab, bty = "n"
    )
    box(lwd = lwd_box)
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_fl(), las = 1, mgp = mgp2, cex.axis = cex_axis, lwd = lwd_box)
    
    for(f in  files_in_al_fluo[-1]){
      points(x_values()[[f]], values_fluorescence()[[f]], type = "l", lty = val$linetype_collected[[f]], 
             lwd = val$linewidth_collected[[f]]*lwd_factor, col = val$color_list[[f]])
    }
  }
  
  # This function creates a legend next to the alignment plot.
  # The size of the legend is adjusted to the length of the profile names.
  legend_alignment <- function(cex_legend, lwd_factor){
    req(val$files_to_plot)
    par(mar = c(0, 0, 0, 0) )
    plot(0, 0, yaxt = "n", xaxt = "n", xlab = "", ylab = "", type = "n", bty = "n")
    legend_names <- val$legend_names[val$files_to_plot]
    legend_adjustment <- sqrt( min( c(1, 8/max( nchar(legend_names) ) ) ) )
    legend("topright", legend = legend_names, lty = unlist(val$linetype_collected[val$files_to_plot]), 
           lwd = unlist(val$linewidth_collected[val$files_to_plot])*lwd_factor, col = unlist(val$color_list[val$files_to_plot]),
           bty = "n", cex = cex_legend*legend_adjustment)
  }
  
  # This function divides the plotting area into two columns for alignment plot and legend,
  # and two rows when polysome and fluorescence alignment plots have to be combined.
  # At the same time, it generates a table of aligned values which is shown in the tab "Alignment table".
  plot_alignment <- function(cex_lab, cex_axis, lwd, 
                             mgp1, mgp2, 
                             cex_legend, lwd_factor,
                             mar_factor, lwd_box){
    if(input$show_fl_al & any( names(val$baseline_fl) %in% val$files_to_plot  ) ){
      layout(matrix(1:4, 2, 2), height = c(0.6, 0.9), width = c(0.8, 0.2) ) 
      par(mar = c(0, 5, 0.5, 0.1)*mar_factor)
      plot_alignedFluo(cex_lab, cex_axis, lwd, 
                       mgp2, lwd_factor, lwd_box)
      par(mar = c(5, 5, 0, 0.1)*mar_factor) 
      plot_alignedPol(cex_lab, cex_axis, lwd, 
                      mgp1, mgp2, 
                      cex_legend, lwd_factor, lwd_box)
      legend_alignment(cex_legend, lwd_factor)
      create_df_pol_alignment()
      create_df_fluo_alignment()
      merge_df_alignment()
    }else{
      layout(matrix(1:2, 1, 2), height = 1, width = c(0.8, 0.2) ) 
      par(mar = c(5, 5, 0.5, 0.1)*mar_factor)
      plot_alignedPol(cex_lab, cex_axis, lwd, 
                      mgp1, mgp2, 
                      cex_legend, lwd_factor, lwd_box)
      legend_alignment(cex_legend, lwd_factor)
      create_df_pol_alignment()
    }
  }
  
  # If there is at least one file in the alignment,
  # the alignment plot is generated as output.
  output$plot_align <- renderPlot({
    req(val$files_to_plot)
    if(length(val$files_to_plot) >= 1)
    {
      plot_alignment(cex_lab = 1.6, cex_axis = 1.6, lwd = 2, 
                     mgp1 = c(3.5, 1.2, 0), mgp2 = c(3.5, 0.8, 0), 
                     cex_legend = 1.5, lwd_factor = 1,
                     mar_factor = 1, lwd_box = 1)
    }
  })
  
  ## DOWNLOADS
  
  # The current plot of an individual polysome profile can be downloaded.
  # The width is reduced by 20% compared to the file of the aligned profiles,
  # to ensure the same size of the plot, although the aligned profiles also
  # contain a legend.
  output$downloadPlot_single <- downloadHandler(
    filename = function(){ 
      paste( 
        sub(pattern = "(.*?)\\..*$", replacement = "\\1", as.character(input$select) ),
        ".pdf", sep = ""
      )
    },
    content = function(file) {
      height_factor <- c(1, 1.75)[input$show_fl + 1]
      pdf(file, width = 1.92, height = 1.44*height_factor, pointsize = 8 )
      print( plot_singleInput(cex_lab = 1, cex_axis = 0.8, lwd = 0.6, 
                              mgp1 = c(1.7, 0.5, 0), mgp2 = c(2, 0.8, 0), 
                              mar_factor = 0.6, lwd_box = 0.5) )
      dev.off()
    })  
  
  # The current alignment plot can be downloaded.
  output$downloadPlot_aligned <- downloadHandler(
    filename = "alignment.pdf",
    content = function(file) {
      height_factor <- c(1, 1.75)[input$show_fl_al + 1]
      pdf(file, width = 2.4, height = 1.44*height_factor, pointsize = 8 )
      print( plot_alignment(cex_lab = 1, cex_axis = 0.8, lwd = 0.6, 
                            mgp1 = c(1.7, 0.5, 0), mgp2 = c(2, 0.8, 0),  
                            cex_legend = 0.6, lwd_factor = 0.4,
                            mar_factor = 0.6, lwd_box = 0.5) )
      dev.off()
    })  
  
  # A table with quantified areas is shown as output in the tab "Quantification".
  output$quantification <- renderTable(
    val$df_quant
  )
  
  # The table with quantified areas can be downloaded from the "Quantification" tab.
  output$downloadQuant <- downloadHandler(
    filename = "AreaQuant.csv",
    content = function(file) {
      write.csv2( val$df_quant, file, row.names = F)
    }
  )
  
  # A table with aligned profile values is shown as output in the tab "Alignment table".
  output$csv_file <- renderTable(
    if(input$show_fl_al){
      val$csv_file_df_all
    }else{
      val$csv_file_df
    }
  )
  
  # The table with aligned profile values can be downloaded from the "Alignment table" tab.
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
  
  # The entire analysis can be exported as .RData file, 
  # which allows the user to re-import it
  # and continue the analysis. 
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
  
  ## NOTIFICATIONS
  
  # Two additional notifications are part of the functions to import data.
  
  # A notification is shown if one of the three required values for area quantification (start, end and baseline) is missing, 
  # but the user activates the "Quantify" button.
  observeEvent(input$quantify_area, {
    if( !isTruthy(val$file_starts[[input$select]]) | !isTruthy(val$file_ends[[input$select]]) | !isTruthy(val$baseline[input$select]) ){
      showNotification(paste("Please set start, end and baseline of the area that you would like to quantify."),
                       duration = NULL, type = "message")
    }
  })
  
  # A notification is shown if "Show" is selected for fluorescence profiles,
  # but there is no fluorescence signal in the currently selected file.
  observeEvent(input$show_fl,{
    req(input$select)
    if( input$show_fl & !isTruthy(val$fluo_data[[input$select]]) ){
      showNotification("Your file does not contain a fluorescence signal.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl", value = FALSE)
    }
  })
  
  # A notification is shown if "Show" is selected for fluorescence profiles 
  # in the alignment but there is not fluorescence signal,
  # or there was no baseline set for at least one of the fluorescence profiles.
  observeEvent(input$show_fl_al,{
    if( !any( val$files_to_plot %in% names(val$fluo_data) ) & input$show_fl_al){
      showNotification("None of your aligned profiles contains a fluorescence signal.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( !all( val$files_to_plot %in% names(val$fluo_data)  ) & input$show_fl_al  ){
        showNotification("Some of your aligned profiles do not contain a fluorescence signal.",
                         duration = NULL, type = "warning")
      }
    }
    
    fluo_to_plot <- intersect(names(val$fluo_data), val$files_to_plot ) 
    if(length(fluo_to_plot) >= 1){
      if( !any( fluo_to_plot %in% names(val$baseline_fl[fluo_to_plot]) ) & input$show_fl_al ){
        showNotification("You did not set a baseline for your fluorescence profiles.",
                         duration = NULL, type = "error")
        updateCheckboxInput(session, "show_fl_al", value = FALSE)
      }else{
        if( !all( fluo_to_plot %in% names(val$baseline_fl[fluo_to_plot])  ) & input$show_fl_al ){
          showNotification("Some of your fluorescence profiles do not have a baseline.",
                           duration = NULL, type = "warning")
        }
      }
    }
  })
  
  # A notification is shown if "Show" was selected already but a new file in val$files_to_plot
  # does not contain a baseline for the fluorescence signal.
  observeEvent(val$files_to_plot, {
    fluo_to_plot <- intersect(names(val$fluo_data), val$files_to_plot ) 
    if( length(fluo_to_plot) == 0 & input$show_fl_al){
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( !all( fluo_to_plot %in% names(val$baseline_fl[fluo_to_plot])  ) & input$show_fl_al ){
        showNotification("Some of your fluorescence profiles do not have a baseline.",
                         duration = NULL, type = "warning")
      }
    }
  })
 
  # A notification is shown when normalization to length or height is selected but the total area 
  # was not quantified for ate least one of the profiles that are shown in the alignment.
  observe({
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( val$files_to_plot %in% files_with_total ) & (input$normalize_length | input$normalize_height) ){
      showNotification("Please select a total area for all profiles in the alignment.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "normalize_length", value = FALSE)
      updateCheckboxInput(session, "normalize_height", value = FALSE)
    }
  })
  
}

###############
shinyApp(ui = ui, server = server)