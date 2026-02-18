#####################################################
#     ____                 _____  _____             #
#    / __ \          /\   |  __ \|  __ \            #
#   | |  | |_   _   /  \  | |__) | |__) | __ ___    #
#   | |  | | | | | / /\ \ |  ___/|  ___/ '__/ _ \   #
#   | |__| | |_| |/ ____ \| |    | |   | | | (_) |  #
#    \___\_\\__,_/_/    \_\_|    |_|   |_|  \___/   #
#                                                   #                                                
#####################################################                                                 

# This is version 1.0.0 of QuAPPro, a shiny app for 
# interactive Quantification and Alignment
# of Polysome Profiles. It is available under the 
# GNU General Public License v3.0.

# Authors: Chiara Schiller, Matthias Lemmer, Johanna Schott

# Structure:

#### REQUIRED PACKAGES ####
#### FUNCTIONS ####
#### USER INTERFACE ####
#### SERVER ####
## LOADING FILES
## USER INPUT FOR PLOTTING INDIVIDUAL PROFILES
## PLOT INDIVIDUAL PROFILES
## INTERACTION OF THE USER WITH THE PROFILE PLOT
## GENERATE ALIGNMENT
## DISPLAY AND MODIFY ALIGNMENT
## PLOT ALIGNMENTS
## DECONVOLUTION
## BAR PLOT
## DOWNLOADS
## NOTIFICATIONS
## GUIDED TOURS

#### REQUIRED PACKAGES ####

library(shiny)
library(shinythemes)
library(colourpicker)
library(colorspace)
library(stringr)
library(shinyalert)
library(DT)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(cicerone)
library(shinybusy)

#### FUNCTIONS ####

# Given coordinates of a selected point (x, y),
# this function returns the index of the closest point
# among the coordinates X and Y according to the Euclidean distance.
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

# This function identifies the closest local minimum  
# in a vector of values (y) to a selected index a of the vector.
find_closest_min <- function(y, a, r){
  minima <- c( all_min( y, r) )
  closest <- order( abs(minima -  a), decreasing = F )[1]
  return(minima[closest])
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

# Given a set of values (y), this function identifies 
# the indices of all inflection points
# defined by local minima or maxima of the 
# first derivative of y.
all_inflections <- function(y, r)
{
  d <- diff(y)
  minmax <- c( all_min(d, r), all_max(d, r) )
  all_minmax <- minmax + 1
  return(all_minmax)
}

# This function identifies the closest inflection point 
# in a vector of values (y) to a selected index a of the vector.
find_closest_inflections <- function(y, a, r){
  infl <- all_inflections(y, r)
  closest <- order( abs(infl -  a), decreasing = F )[1]
  return(infl[closest])
}

# This function applies smooth.spline() to a given vector of values y
# with the smoothing parameter spar. 
# If spar is zero, the values are returned without smoothing.
smooth_profile <- function(y, spar){
  if(!isTruthy(spar) )  {
    return(y)
  }else{
    if(spar == 0){
      return(y)
    }else{
      return( smooth.spline(1:length(y), y, spar = spar)$y )
    }
  }
}

# This function determines the sum of several normal density distributions
# at the positions pos with the given heights and standard deviations,
# provided as vectors.
gauss_sum <- function(x, pos, height, sd){
  gauss_overlay <- rep(0, length(x))
  for(i in 1:length(pos)){
    gauss <- dnorm(x, mean = pos[i], sd = sd[i])
    norm_factor <- height[i]/max(gauss)
    gauss_overlay <- gauss_overlay + gauss*norm_factor
  }
  gauss_overlay
}

#### USER INTERFACE ####

ui <- fluidPage(
  
  # Use the cicerone guide for the guided tour:
  use_cicerone(),
  
  # Show a spinner when the server is busy:
  add_busy_spinner(spin = "fading-circle", position = "full-page", timeout = "1000"),
  
  # Browser tab title:
  title = tags$head(tags$title("QuAPPro at MI3")),
  # Use version 5.7 of fontawesome
  tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
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
           tags$img(src='logo_mfm_ukm.svg', align = "right", width = '100%', height = '100%')
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
    tabPanel(tags$strong("Alignment/Areas"), icon = icon("chart-area"),
             fluidRow(
               column(2, 
                      wellPanel( style = "padding-top:0px",
                                 fluidRow( style = "padding-bottom:10px;padding-left:10px",
                                           tags$h3(tags$strong("Import and Export")),
                                 ),
                                 # Tutorial with example data and guided tour:
                                 tags$h4(tags$strong("Take a Tour")),
                                 fluidRow( style = "padding-top:0px;padding-bottom:10px",
                                           column(6, 
                                                  actionButton("example1", label = "Example", icon = icon("chart-area"),
                                                               width = "100%")
                                           ),
                                           column(6,
                                                  actionButton("tour1", label = "Tour", icon = icon("circle-info"),
                                                               width = "100%")
                                           )
                                 ),
                                 tags$h4(tags$strong("File Import")), 
                                 # fileInput for uploading files
                                 div(
                                   id = "input_data_wrapper",
                                   fileInput("input_data", label = NULL, multiple = T, accept = c(".pks",".csv", ".txt", ".asc", ".tsv", ".RData"))
                                 ),
                                 # Menu for selecting file type
                                 fluidRow(
                                   column(6, 
                                          div(
                                            id = "type_wrapper",
                                            selectInput("type", "Profile type", 
                                                        choices = c("PeakTrak", "TRIAX (UV)", "TRIAX (+ Fl.)", "PrimeView", "CustomTSV"),
                                                        selected = "PeakTrak",
                                                        width = '110%')                               )
                                          
                                   ),
                                   column(6,
                                          numericInput("skip", "Skip lines", value = 0, step = 1, min = 0)
                                   )
                                 ),
                                 
                                 fluidRow( 
                                   column(6, 
                                          selectInput("sep", "Column delim.", 
                                                      choices = c("space", "tab", ",", ";"), 
                                                      selected = "space",
                                                      width = '110%')
                                   ),
                                   
                                   column(6,
                                          selectInput("dec", "Dec. separator", 
                                                      choices = c(".", ","), selected = ",",
                                                      width = '110%')
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          numericInput("UV_col", "UV column", 
                                                       value = 3, step = 1, min = 1, width = '110%')
                                   ),
                                   column(6,
                                          numericInput("fluo_col", "Fluor. column", 
                                                       value = NA, step = 1, min = 1, width = '110%')
                                   )
                                 ),
                                 
                                 # Menu for selecting which of the uploaded files should be displayed
                                 fluidRow(
                                   column(12, style = "padding-left:15px;padding-top:0px;padding-right:15px",
                                          selectInput("select", "Select a profile", choices = c(), width = '100%')
                                   )
                                 ), 
                                 
                                 # Download entire analysis
                                 tags$h4(tags$strong("Export Analysis")),
                                 fluidRow(
                                   style = "padding-bottom:0px;padding-left:15px",
                                   tags$h6("Export your entire analysis as .RData. You can re-import it via File Import."),
                                   downloadButton("export", "Export Analysis", width = "100%")
                                 ),
                                 fluidRow(  style = "padding-left:0px;padding-top:65px",
                                            tags$img(src='Logo-MI3-RGB.png', align = "left", width = '90%', height = '90%')
                                 )
                      )
               ),
               
               # Section for plots of individual profiles
               column(5, 
                      wellPanel( style = "padding-top:0px;",
                                 fluidRow( style = "padding-bottom:10px;padding-left:10px",
                                           tags$h3(tags$strong("Process and Quantify Profiles")),
                                 ),
                                 plotOutput("plot_single", click = "click"),
                                 fluidRow( style = "padding-top:10px",
                                           column(3,
                                                  checkboxInput("red_lines", "Show x-anchor", value = TRUE, width = NULL)
                                           ),
                                           column(3,
                                                  checkboxInput("red_lines2", "Show baseline", value = TRUE, width = NULL)
                                           ),
                                           column(3,
                                                  checkboxInput("green_lines", "Show quantified area", value = TRUE, width = NULL)
                                           ),
                                           column(3, 
                                                  checkboxInput("start_end", "Show start/end", value = TRUE, width = NULL)
                                           )
                                 ),
                                 fluidRow( style = "border-bottom: 1px dotted #d3d3d3;padding-top:14px",
                                           column(2, style = "padding-top:0px",
                                                  numericInput("axis3", "Set x min", value = NULL, step = 1) ),
                                           column(2, style = "padding-top:0px",
                                                  numericInput("axis4", "Set x max", value = NULL, step = 1) ),
                                           column(2,
                                                  checkboxInput("helper_functions", "Min/max detection", value = TRUE, width = NULL)
                                           ),
                                           column(2,
                                                  checkboxInput("inflection_functions", "Inflection detection", value = FALSE, width = NULL)
                                           ),
                                           column(4, align = "right", style = "padding-top:15px", 
                                                  downloadButton("downloadPlot_single", "Download Plot", 
                                                                 icon = icon("file-download"), width = '100%')
                                           ),
                                 ),
                                 fluidRow(
                                   # Section for UV profile
                                   column(4, style = "border-right: 1px dotted #d3d3d3",
                                          fluidRow(
                                            column(6, style = "padding-top:10px",
                                                   tags$h4(tags$strong("UV Profile") )   
                                            ),
                                            column(6, style = "padding-top:20px",
                                                   numericInput("smooth_pol", "Smoothing", value = 0, min = 0, max = 1, step = 0.01) 
                                            )  
                                          ),
                                          fluidRow( style = "padding-bottom:10px",
                                                    column(6, style = "padding-top:10px", 
                                                           actionButton("baseline", "Baseline", width = '100%')
                                                    ),
                                                    column(6,style = "padding-top:10px", 
                                                           actionButton("x_anchor", "X-anchor", width = '100%')
                                                    )
                                          ),
                                          fluidRow(
                                            column(6, style = "padding-top:10px",
                                                   numericInput("axis1", "Set y min", value = NULL, step = 1) 
                                            ),
                                            column(6, style = "padding-top:10px",
                                                   numericInput("axis2", "Set y max", value = NULL, step = 1) 
                                            )
                                          )
                                   ),
                                   # Section for fluorescence profile
                                   column(4, style = "border-right: 1px dotted #d3d3d3",
                                          fluidRow(
                                            column(6, style = "padding-top:10px",
                                                   tags$h4(tags$strong("Fluorescence") ) 
                                            ),
                                            column(6, style = "padding-top:20px",
                                                   numericInput("smooth_fluo", "Smoothing", value = 0, min = 0, max = 1, step = 0.01) 
                                            )
                                          ),
                                          fluidRow( style = "padding-bottom:10px",
                                                    # Set a separate baseline for the fluorescence signal
                                                    column(6, style = "padding-top:10px", 
                                                           actionButton("baseline_fl", "Baseline", width = '100%') 
                                                    ),
                                                    # Show fluorescence signal or not?
                                                    column(6, style = "padding-top:0px;padding-bottom:0px",
                                                           checkboxInput("show_fl", "Show", value = FALSE, width = '100%')
                                                    )
                                          ),
                                          fluidRow(
                                            column(6, style = "padding-top:10px",
                                                   numericInput("axis1_fl", "Set y min", value = NULL, step = 1) ),
                                            column(6, style = "padding-top:10px",
                                                   numericInput("axis2_fl", "Set y max", value = NULL, step = 1) )
                                          )
                                   ),
                                   # Section for quantification
                                   column(4,
                                          fluidRow(
                                            column(6, style = "padding-top:10px",
                                                   tags$h4(tags$strong("Quantify") )
                                            ),
                                            column(6, style = "padding-top:20px",
                                                   div( id = "select_area_wrapper",
                                                        selectizeInput("select_area", "Set name", options = list(create=TRUE, plugins = list('restore_on_backspace')),
                                                                       choices = c("Total", "80S", "Polysomes", "40S", "60S"), width = '100%')
                                                   )
                                            )
                                          ),
                                          fluidRow( style = "padding-bottom:10px",
                                                    column(6, style = "padding-top:10px",
                                                           actionButton("file_start", "Start", width = '100%')
                                                    ),
                                                    column(6, style = "padding-top:10px",
                                                           actionButton("file_end", "End", width = '100%')
                                                    )
                                          ),
                                          fluidRow(
                                            column(12, style = "padding-top:32px",
                                                   actionButton("quantify_area", "Quantify", width = '100%')
                                            )
                                          )
                                   )
                                 )
                      )),
               # Section for showing and modifying the alignment
               column(5,
                      wellPanel( 
                        style = "padding-top:0px;",
                        fluidRow( style = "padding-bottom:10px;padding-left:10px",
                                  tags$h3(tags$strong("Visualize Alignments")),
                        ),
                        plotOutput("plot_align"),
                        fluidRow( style = "padding-top:10px",
                                  column(3,
                                         checkboxInput("anchor_line", "Show x-anchor", value = TRUE, width = NULL)    
                                  ),
                                  column(6,
                                         radioButtons("color_palette", "Default color palette:",
                                                      c("Dark" = "dark_palette",
                                                        "Rainbow" = "rainbow_palette",
                                                        "Color blind friendly" = "color_blind"), 
                                                      selected = "dark_palette", inline = TRUE),
                                         style = "padding-left:30px"
                                  )
                        ),
                        fluidRow( style = "border-bottom: 1px dotted #d3d3d3;",
                                  column(2, style = "padding-top:0px",
                                         numericInput("axis3_a", "Set x min", value = NULL, step = 1) ),
                                  column(2, style = "padding-top:0px",
                                         numericInput("axis4_a", "Set x max", value = NULL, step = 1) ),
                                  column(2,
                                         checkboxInput("normalize_height", HTML("Normalize <b>area</b>"), 
                                                       value = FALSE, width = NULL)
                                  ),
                                  column(2,
                                         checkboxInput("normalize_length", HTML("Normalize <b>length</b>"), 
                                                       value = FALSE, width = NULL)
                                  ),
                                  column(4, align = "right", style = "padding-top:15px", 
                                         downloadButton("downloadPlot_aligned", "Download Plot", 
                                                        icon = icon("file-download"), width = '100%')
                                  ),
                        ),
                        fluidRow(
                          # Section for UV profile
                          column(4, style = "border-right: 1px dotted #d3d3d3;border-bottom: 1px dotted #d3d3d3",
                                 fluidRow(
                                   column(12, style = "padding-top:10px",
                                          tags$h4(tags$strong("UV Profile") )   
                                   )
                                 ),
                                 fluidRow(
                                   column(6, style = "padding-top:10px;padding-bottom:15px",
                                          numericInput("axis1_a", "Set y min", value = NULL, step = 1) 
                                   ),
                                   column(6, style = "padding-top:10px;padding-bottom:15px",
                                          numericInput("axis2_a", "Set y max", value = NULL, step = 1) 
                                   )
                                 )
                          ),
                          # Section for fluorescence profile
                          column(4, style = "border-right: 1px dotted #d3d3d3;border-bottom: 1px dotted #d3d3d3",
                                 fluidRow(
                                   column(6, style = "padding-top:10px",
                                          tags$h4(tags$strong("Fluorescence") ) 
                                   ),
                                   column(6, style = "padding-top:0px",
                                          checkboxInput("show_fl_al", "Show", value = FALSE, width = '100%')
                                   )
                                 ),
                                 fluidRow(
                                   column(6, style = "padding-top:10px;padding-bottom:15px",
                                          numericInput("axis1_a_fl", "Set y min", value = NULL, step = 1) ),
                                   column(6, style = "padding-top:10px;padding-bottom:15px",
                                          numericInput("axis2_a_fl", "Set y max", value = NULL, step = 1) )
                                 )
                          ),
                          column(4,
                                 fluidRow(
                                   column(6, style = "padding-top:10px",
                                          tags$h4(tags$strong("Colors/Lines") ) 
                                   ),
                                   column(6, style = "padding-top:10px",
                                          selectInput("linetype", "Line type", choices = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                                                      selected = "solid")
                                   )
                                 ),
                                 fluidRow(
                                   column(6, style = "padding-top:0px",
                                          colourInput("color", "Color", palette = "square", value = "#E16A86")
                                   ),
                                   column(6, style = "padding-top:0px",
                                          numericInput("linewidth", "Width", value = 1, min = 0)  
                                   )
                                 )
                          )
                        ),
                        # Section for order of profiles 
                        fluidRow(
                          column(4,
                                 div(
                                   id = "select_alignment_wrapper",
                                   selectInput("select_alignment", "Select a profile", choices = c(), width = '100%')
                                 )
                          ),
                          column(2,
                                 column(6, style = "padding-top:22px",
                                        actionButton("up", label = "\u21D1", width = '140%',
                                        ) 
                                 ),
                                 column(6, style = "padding-top:22px",
                                        actionButton("down", label = "\u21D3", width = '140%'
                                        ) 
                                 ),
                          ),
                          column(2, style = "padding-top:22px",
                                 checkboxInput("show_in_al", "Display", value = TRUE)
                          ),
                          column(2, style = "padding-top:0px",
                                 textInput("new_name", "Label", value = "", width = '100%', placeholder = "Name") 
                          ),
                          column(2, style = "padding-top:22px",
                                 actionButton("rename", "Rename", width = '100%')
                          )
                        )
                        
                      ))
             )
    ),
    
    # SECOND TAB: DECONVOLUTION
    # Shows the selected polysome profile with its second derivative 
    # for deconvolution and quantification of individual peaks.
    tabPanel(tags$strong("Deconvolution"), icon = icon("cut"),
             fluidRow(
               # Section for uploading and selecting files:
               column(2,
                      wellPanel( style = "padding-top:0px",
                                 fluidRow( style = "padding-bottom:10px;padding-left:10px",
                                           tags$h3(tags$strong("Deconvolution")),
                                 ),
                                 # Tutorial with example data and guided tour:
                                 tags$h4(tags$strong("Take a Tour")),
                                 fluidRow( style = "padding-top:0px;padding-bottom:10px",
                                           column(6, 
                                                  actionButton("example2", label = "Example", icon = icon("chart-area"),
                                                               width = "100%")
                                           ),
                                           column(6,
                                                  actionButton("tour2", label = "Tour", icon = icon("circle-info"),
                                                               width = "100%")
                                           )
                                 ),
                                 # Section settings of UV profile
                                 fluidRow(
                                   column(12, tags$h4(tags$strong("UV Profile") )
                                   )
                                 ),
                                 # Menu for selecting which of the profiles should be displayed
                                 # (among all profiles with a baseline)
                                 fluidRow(
                                   column(12, style = "padding-left:15px;padding-top:0px;padding-right:15px",
                                          div(id = "select2_wrapper",
                                              selectInput("select2", "Select a profile", choices = c(), width = '100%')
                                          )
                                   )
                                 ), 
                                 
                                 # Set axis limits for UV profile
                                 fluidRow(
                                   column(6, style = "padding-top:10px",
                                          numericInput("axis1_d", "Set y min", value = NULL, step = 1) ),
                                   column(6, style = "padding-top:10px",
                                          numericInput("axis2_d", "Set y max", value = NULL, step = 1) )
                                 ),
                                 
                                 fluidRow(
                                   column(6, style = "padding-top:0px",
                                          numericInput("axis3_d", "Set x min", value = NULL, step = 1) ),
                                   column(6, style = "padding-top:0px",
                                          numericInput("axis4_d", "Set x max", value = NULL, step = 1) )
                                 ),
                                 
                                 # Section for deconvolution settings
                                 fluidRow(
                                   column(12, tags$h4(tags$strong("Peak Detection") )
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          numericInput("smooth_pol_decon", "Smoothing", value = 0.2, min = 0, max = 1, step = 0.01) 
                                   ),
                                   column(6,
                                          numericInput("resolution", "Resolution", value = 5, min = 1, max = 100, step = 1) 
                                   )
                                 ),
                                 
                                 fluidRow(
                                   column(12, tags$h4(tags$strong("Deconvolution") )
                                   )
                                 ),
                                 
                                 fluidRow(style = "padding-top:10px",
                                          column(4,
                                                 actionButton("decon_start", "Start", width = '100%')
                                          ),
                                          column(4,
                                                 actionButton("decon_end", "End", width = '100%')
                                          ),
                                          column(4,
                                                 actionButton("decon_go", "Go", width = '100%')
                                          )
                                 ),
                                 fluidRow(
                                   column(12,
                                          checkboxInput("helper_functions_decon", "Min/max detection", value = TRUE, width = NULL),
                                          checkboxInput("inflection_functions_decon", "Inflection detection", value = TRUE, width = NULL),
                                          checkboxInput("derivMin_functions_decon", "2nd deriv. minima detection", value = FALSE, width = NULL)
                                   )
                                 ),
                                 
                                 # Section for quantification of individual peaks
                                 fluidRow(
                                   column(12, tags$h4(tags$strong("Quantification") )
                                   )
                                 ),
                                 fluidRow(
                                   column(6, style = "padding-top:15px",
                                          div(
                                            id = "select_peak_wrapper",
                                            selectizeInput("select_peak", "Select a name", options = list(create=TRUE, plugins = list('restore_on_backspace')),
                                                           choices = c("40S", "60S", "80S"), width = '100%')
                                          )
                                   ),
                                   column(6, style = "padding-top:35px",
                                          actionButton("quantify_peak", "Quantify", width = '100%')
                                   )
                                 )
                      )
               ),
               
               # Section for plots of individual profiles: 
               column(8, style = "padding-top:22px",
                      plotOutput("plot_deconv", click = "click_deconv", height = 600),
                      column(8, 
                             checkboxInput("show_deriv", "Show 2nd deriv.", value = TRUE, width = NULL)
                      ),
                      column(4, align = "right", style = "padding-right:0px;padding-top:15px",
                             downloadButton("downloadPlot_decon", "Download Plot", 
                                            icon = icon("file-download"), width = NULL)
                      )
               ),
               column(2, style = "padding-right:22px",
                      fluidRow( style = "padding-top:20px",
                                tags$h4(tags$strong("Model Parameters") ),
                      ),
                      textOutput("MSE"),
                      tableOutput("decon_param"),
                      fluidRow(
                        column(12,
                               downloadButton("downloadParam", "Download .csv file")
                        )
                      )
                      
               )
             )
    ),
    
    # THIRD TAB: ALIGNMENT TABLE
    # Shows updating table of aligned (and normalized) profiles, 
    # which can be download.
    tabPanel(tags$strong("Alignment Table"), icon = icon("table"), 
             tags$h4("Table of all aligned (and normalized) profiles"),
             downloadButton("downloadData", "Download .csv file"),
             tableOutput("csv_file")
    ),
    
    # FOURTH TAB: QUANTIFICATION TABLE
    # Shows updating table of quantified areas for respective plots, 
    # which can be downloaded. 
    tabPanel(tags$strong("Quantification Summary"), icon = icon("list-alt"), 
             tags$h4("Table of quantified areas"),
             tableOutput("quantification"),
             downloadButton("downloadQuant", "Download .csv file")
    ),
    
    # FIFTH TAB: BAR PLOTS
    tabPanel(tags$strong("Bar Plot"), icon = icon("chart-bar"),
             fluidRow(
               column(4,
                      wellPanel( style = "padding-top:0px",
                                 fluidRow( style = "padding-bottom:10px;padding-left:10px",
                                           tags$h3(tags$strong("Bar Plot Summary")),
                                 ),
                                 # Tutorial with example data and guided tour:
                                 tags$h4(tags$strong("Take a Tour")),
                                 # Tutorial with example data and guided tour:
                                 fluidRow( style = "padding-top:10px;padding-bottom:10px",
                                           column(3, 
                                                  actionButton("example3", label = "Example Data", icon = icon("chart-area"),
                                                               width = "100%")
                                           ),
                                           column(3,
                                                  actionButton("tour3", label = "Guided Tour", icon = icon("circle-info"),
                                                               width = "100%")
                                           )
                                 ),
                                 tags$h4(tags$strong("Select Values")),
                                 fluidRow(
                                   column(6, 
                                          div( id = "input_table_barplot_wrapper",
                                               selectInput("input_table_barplot", label = "Area or peak", choices = NULL)
                                          )
                                   ),
                                   column(6,
                                          div( id = "normalize_barplot_wrapper",
                                               selectInput("normalize_barplot", label = "Normalization", choices = "None")
                                          )
                                   )
                                 ),
                                 tags$h4(tags$strong("Condition Assignment") ),
                                 DTOutput("barplot_table"),
                                 downloadButton("download_barplot_table", "Download Table")
                      )
               ),
               column(6, style = "padding-top:22px",
                      plotOutput("barplot"),
                      column(8),
                      column(4, align = "right", style = "padding-right:0px;padding-top:15px",
                             downloadButton("download_barplot", "Download Plot")
                      ),
                      fluidRow(
                        column(6,
                               wellPanel(
                                 tags$h4(tags$strong("Colors")),
                                 fluidRow(
                                   column(4, 
                                          div( id = "selected_cond1_wrapper",
                                               selectInput("selected_cond1", label = "Select 'Condition1'", choices = NULL) 
                                          )
                                   ),
                                   column(4, colourInput("color_cond1", "Bars", palette = "square", value = "#858181")
                                   ),
                                   column(4, colourInput("color_dots_cond1", "Data points", palette = "square", value = "#000000")
                                   )
                                 ),
                                 fluidRow(
                                   column(4, 
                                          div( id = "selected_cond2_wrapper",
                                               selectInput("selected_cond2", label = "Select 'Condition2'", choices = NULL) 
                                          )
                                   ),
                                   column(4, colourInput("color_cond2", "Bars", palette = "square", value = "#858181")
                                   ),
                                   column(4, colourInput("color_dots_cond2", "Data points", palette = "square", value = "#000000")
                                   )
                                 )
                               )
                        ),
                        column(6)
                      )
               ),
               column(2, style = "padding-right:22px",
                      
               )
             )
    ),
    
    # SIXTH TAB: MANUAL 
    # Shows manual  
    tabPanel(tags$strong("QuAPPro Manual"), icon = icon("question-circle"),
             fluidRow(
               column(1),
               column(10,
                      tags$div(
                        style = "font-size: 18px;",
                        htmltools::includeMarkdown("QuAPPro_v1-0-0_manual.Rmd")
                      )
               ),
               column(1)
             ) 
    ),
    
    # SEVENTH TAB: RELEASE NOTES
    # Release notes with link to GitHub repository:
    tabPanel(tags$strong("Release Notes", style = "color:blue"), icon = icon("info", style = "color:blue"), 
             fluidRow(
               column(1),
               column(10, 
                      tags$div(
                        style = "font-size: 18px;",
                        htmltools::includeMarkdown("QuAPPro_v1-0-0_releaseNotes.Rmd"))
               ), 
               column(1)
             ) 
    ),
    
    # EIGHTH TAB: CONTACT
    # Contact information and acknowledgments
    tabPanel(tags$strong("Info", style = "color:blue"), icon = icon("address-card", style = "color:blue"), 
             tags$h4("About Us"),
             tags$div("This shiny app was developed by Chiara Schiller, Matthias Lemmer and",
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
             tags$div("QuAPPro is running on cloud infrastructure of Heidelberg University",
                      tags$a(href="https://www.urz.uni-heidelberg.de/de/service-katalog/cloud/heicloud", "(HeiCLOUD)", style = "color:blue"),
                      "."),
             tags$h4(tags$a(href="https://www.umm.uni-heidelberg.de/impressum/", "Impressum", style = "color:blue") )
    )
  )
)
