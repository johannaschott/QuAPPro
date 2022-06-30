####################################
###     Script Shiny R 2021      ###
###  Polysome profile analysis   ###
####################################


library(shiny)
library(shinythemes)
library(colourpicker)
library(stringr)
library(shinyFiles)
library(plyr)
library(colorspace)
library(markdown)

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
           top: calc(50%);
           left: calc(50%);
           }
           ")
      )
      ),
  # create fluid layout with several tabs for displaying outputs
  tabsetPanel(
    # FIRST TAB
    tabPanel(tags$strong("Analysis profiles"), icon = icon("chart-area"), # updated to newer version
             fluidRow(column(2),
                      column(4, tags$h4(tags$strong("Single profiles for quantification"))),
                      column(6, tags$h4(tags$strong("Multiple aligned profiles")))),
             
             fluidRow(
               column(2, style = "padding-left:20px",
                      # create area for uploading .pks file
                      fileInput("input_data", "Upload pks, txt or csv file", multiple = T, accept = c(".pks",".csv", ".txt")),
                      
                      #Let user select their loaded files and set x-anchor and baseline
                      selectInput("select", "Select files", choices = c(), width = '100%'),
                      
                      
                      
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow( 
                        column(12, tags$h4(tags$strong("Fluorescence profile")))),
                      fluidRow(         
                        # show fluorescence signal or not?
                        # set a separate baseline for the fluorescence signal
                        
                        column(6, style = "padding-top:10px",actionButton("baseline_fl", "Set Fluo. baseline", width = '90%')
                        ),
                        
                        column(6,style = "padding-top:8px",
                               checkboxInput("show_fl", "Show fluorescence signal", value = FALSE, width = NULL))
                        
                        
                      ),
                      
                      fluidRow(
                        column(6, style = "padding-top:10px",
                               numericInput("axis1_fl", "Set y min", value = NULL)
                        ),
                        column(6, style = "padding-top:10px",
                               numericInput("axis2_fl", "Set y max", value = NULL))),
                      
                      # introduce a slider for smoothing of the fluorescence signal
                      column(12,
                             fluidRow(sliderInput("slider1", label = "Smooth fluorescence profile", min = 0, 
                                                  max = 100, value = 0)))
                      
                      
               ),
               column(4,
                      plotOutput("plot_single", click = "click")
               ),
               column(4,
                      plotOutput("plot_align")
               ),
               column(2, style = "padding-right:22px",
                      # create inputs for user to set x and y axis limits, should start with initially set values before user input 
                      fluidRow(tags$h5(tags$strong("Fluorescence profile")),
                               
                               # show fluorescence signal or not?
                               column(12,  
                                      checkboxInput("show_fl_al", "Show fluorescence signal", value = FALSE, width = NULL),
                                             ),
                               
                               column(6, 
                                      numericInput("axis1_a_fl", "Set y min", value = NULL)),
                               column(6, 
                                      numericInput("axis2_a_fl", "Set y max", value = NULL))
                      ),
                      fluidRow(tags$h5(tags$strong("Polysome profile")),
                               column(6,
                                      numericInput("axis1_a", "Set y min", value = NULL)),
                               column(6,
                                      numericInput("axis2_a", "Set y max", value = NULL))
                      ),
                      fluidRow(
                        column(6, 
                               numericInput("axis3_a", "Set x min", value = NULL)
                        ),
                        column(6,
                               numericInput("axis4_a", "Set x max", value = NULL)
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
                        
                        column(3,style = "padding-top:10px", actionButton("baseline", "Set baseline", width = '100%')
                        ),
                        column(3,style = "padding-top:10px", actionButton("x_anchor", "Set x-anchor", width = '100%')
                        ),
                        
                        column(3, style = "padding-top:10px",
                               actionButton("file_start", "Set area start", icon = icon("caret-square-right"), width = '100%')
                        ),
                        column(3, style = "padding-top:10px",
                               actionButton("file_end", "Set area end", icon = icon("pause-circle"), width = '100%')
                        )
                      ),
                      fluidRow(
                        column(3, style = "padding-top:10px",
                               numericInput("axis1", "Set y min", value = NULL)
                        ),
                        column(3, style = "padding-top:10px",
                               numericInput("axis2", "Set y max", value = NULL)
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
                               numericInput("axis3", "Set x min", value = NULL)
                        ),
                        column(3, style = "padding-top:0px",
                               numericInput("axis4", "Set x max", value = NULL)
                        ),
                        column(3, style = "padding-top:22px",
                               actionButton("quantify_area", "Quantify area", icon = icon("calculator"), width = '100%')
                        ),
                        column(3, style = "padding-top:22px",
                               actionButton("take_over_name", "Add name", icon = icon("check-circle"), width = '100%')
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
                        column(6, 
                               actionButton("remove_file", "Hide profile", icon = icon("trash"), width = '100%')
                        ),
                        column(6, 
                               actionButton("add_file", "Show profile", icon = icon("trash-restore"), width = '100%')
                        ),
                        style = "padding-right:35px",
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
                               colourInput("color", "Colors", palette = "square")
                        ),
                        column(4, style = "padding-left:0px", 
                               selectInput("linetype", "Type", choices = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                                           selected = "solid")
                        ),
                        column(3, style = "padding-left:0px",
                               numericInput("linewidth", "Width", value = 1)
                        )
                      )
                      
               ),
               
               column(2,
                      downloadButton("downloadPlot", "Download plot", icon = icon("file-download"), width = '100%'),
                      checkboxInput("anchor_line", "Display x-anchor in alignment", value = TRUE, width = NULL),
                      checkboxInput("normalize_height", HTML("Normalize <b>height</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL),
                      checkboxInput("normalize_length", HTML("Normalize <b>length</b> in alignment <br/> (Requirement: ALL total areas)"), value = FALSE, width = NULL))
             )
             
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
    # shows Rmd manual  
    tabPanel(tags$strong("QuAPPro manual"), icon = icon("question-circle"),
             fluidRow(
               column(2),
               column(8, htmltools::includeMarkdown("RMD_menu_template_2022.Rmd")),
               column(2))
    ),
    
    # FIFTH TAB
    # general information and impressum
    tabPanel(tags$strong("Release notes", style = "color:blue"), icon = icon("info", style = "color:blue"), 
             tags$h4("Release notes"),
             tags$div("The current version QuAPPro v0.1.0 can be downlaoded",
                      downloadLink("DownloadQuAPPro", label = "here"),
                      "to run it locally."
             )
    ),
    
    # SIXTH TAB
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
  
  ### REACTIVE VALUES FOR LISTS BUILT ALONG THE SESSION
  # create reactive values for collecting values in a list/vector/data.frame
  
  # 
  output$DownloadQuAPPro <- downloadHandler(
    filename = function() {
      paste("QuAPPro_v0-1-0", ".zip")
    },
    content = function(file) {
      file.copy("QuAPPro_v0-1-0.zip", file)
    }
  )
  
  
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
    df_quant = data.frame()
  )
  val <- reactiveValues(
    factors_list = list()
  )
  val <- reactiveValues(
    control_baseline = list()
  )
  val <- reactiveValues(
    fl_control_baseline = list()
  )
  
  val <- reactiveValues(
    file_types = vector()
  )
  
  val <- reactiveValues(
    xmin_collected = list()
  )
  val <- reactiveValues(
    xmax_collected = list()
  )
  val <- reactiveValues(
    ymin_collected = list()
  )
  val <- reactiveValues(
    ymax_collected = list()
  )
  val <- reactiveValues(
    ymin_collected_fl = list()
  )
  val <- reactiveValues(
    ymax_collected_fl = list()
  )
  
  
  ### INITIAL LOADED FILES
  # update select output field with the name of every newly loaded file
  observe({
    session$userData$files_list <- c(input$input_data$name[input$input_data$size != 0], session$userData$files_list)
    updateSelectInput(session, "select",
                      choices = session$userData$files_list)
  })
  
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
  
  observe({
    req(input$select)
    
    if((grepl(".pks",as.character(input$select)) == T)){
      val$factors_list[[input$select]] <- (0.1/60)
      val$file_types[input$select] <- "pks"
    }else if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
      val$factors_list[[input$select]] <- (0.32/60)
      val$file_types[input$select] <- "csv_fluo" 
    }else if((grepl(".csv",as.character(input$select)) == T) & !("SampleFluor" %in% colnames(file_plot()))){
      val$factors_list[[input$select]] <- (0.2/60)
      val$file_types[input$select] <- "csv" 
      print(val$factors_list[[input$select]])
    }else if((grepl(".txt",as.character(input$select)) == T)){
      val$factors_list[[input$select]] <- (0.1/60)
      val$file_types[input$select] <- "pks"
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
      read.table(current_path(), dec = ",", header = F)
    }else if(grepl(".csv",as.character(input$select)) == T){
      start <- grep("Data Columns:", readLines(current_path()))
      read.csv(current_path(), skip = start)
    }else if(grepl(".txt",as.character(input$select)) == T){
      read.delim(current_path())}
  })
  
  # Setting y and x column values for plotting from the loaded dataset
  yvalue <- reactive({
    req(input$select) 
    
    if((grepl(".pks",as.character(input$select)) == T)){
      file_plot()[ ,3]}else if((grepl(".csv",as.character(input$select)) == T) & ("SampleFluor" %in% colnames(file_plot()))){
        file_plot()[ ,5]}else if((grepl(".csv",as.character(input$select)) == T) & (ncol(file_plot()) == 4)){
          head(file_plot()[,2], -1)
        }else if((grepl(".txt",as.character(input$select)) == T)){
          file_plot()[ ,2]
        }
  })
  
  fluorescence <- reactive({
    req(input$select) 
    req(input$show_fl)
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
      }else if((grepl(".txt",as.character(input$select)) == T)){
        (file_plot()[ ,1]+1)*val$factors_list[[input$select]]
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
      max(yvalue())
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
      max(xvalue())
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
    
    # store file baselines for quantified areas to be displayed in the plotting area even when baseline is changed again by the user
    if(isTruthy(val$baseline_fl[input$select])){
      val$fl_control_baseline[[input$select_area]][input$select] <- val$baseline_fl[input$select]
    }
    
    val$control_baseline[[input$select_area]][input$select] <- val$baseline[input$select]
    
    x_first <- round(val$file_starts[[input$select]])*val$factors_list[[input$select]]
    x_last <- round(val$file_ends[[input$select]])*val$factors_list[[input$select]]
    # sum up all values from start to end area and subtract baseline, assign area name selected by user
    # store in list of areas including quantification with respective file name
    val$sum_areas[[input$select_area]][input$select] <- (sum(yvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))]-val$baseline[input$select]))
    
    if(input$show_fl){
      val$sum_areas[[paste(input$select_area, "_fluo", sep = "")]][input$select] <- (sum(fluorescence()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))]-val$baseline_fl[input$select]))
    }
    # create data frame with quantified areas containing NAs for files without respective area
    areas_with_quant <- names(val$sum_areas)
    
    a <- areas_with_quant[1]
    df_quant <- data.frame(files = names(val$sum_areas[[ a ]]), val$sum_areas[[ a ]])
    colnames(df_quant)[2] <- a
    for(a in areas_with_quant[-1])
    {
      df_quant_new <- data.frame(files = names(val$sum_areas[[ a ]]), val$sum_areas[[ a ]])
      colnames(df_quant_new)[2] <- a
      df_quant <- merge(df_quant, df_quant_new, by = "files", all = T)
    }
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
    val$linewidth_collected[[input$select_alignment]] <- input$linewidth
  })
  
  # update color, line type and width to the values of the selected profiles
  observeEvent({input$select_alignment
    input$color_palette},  {
      req(input$select_alignment)
      updateColourInput(session, "color", value = colors_vector()[[input$select_alignment]])
      updateNumericInput(session, "linewidth", value = linewidth_vector()[[input$select_alignment]]  )
      updateSelectInput(session, "linetype", 
                        selected = c("solid", "dashed", "dotted",
                                     "dotdash", "longdash", "twodash")[ lines_vector()[[input$select_alignment]] ])
    })
  
  
  ### OUTPUT
  
  # show_fl is deselected when a new file is selected that does not contain fluorescence
  # no warning or error is shown
  observeEvent(input$select,{
    if(input$show_fl && !("SampleFluor" %in% colnames(file_plot()))){
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
    if(!("SampleFluor" %in% colnames(file_plot())) & input$show_fl){
      showNotification("Your file does not contain a fluorescence signal (column 'SampleFluor').",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl", value = FALSE)
    }
  })
  
  # let axis limits update for fluorescence axis
  observe({
    if(("SampleFluor" %in% colnames(file_plot())) & input$show_fl){
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
    if( !any( val$file_types[files_to_plot()] == "csv_fluo") & input$show_fl_al){
      showNotification("None of your aligned profiles contains a fluorescence signal (column 'SampleFluor').",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( !all( val$file_types[files_to_plot()] == "csv_fluo") & input$show_fl_al  ){
        showNotification("Some of your aligned profiles do not contain a fluorescence signal (column 'SampleFluor').",
                         duration = NULL, type = "warning")
      }
    }
    
    if(any( val$file_types[files_to_plot()] == "csv_fluo") & is.null(val$baseline_fl) & input$show_fl_al){
      showNotification("You did not set a baseline for your fluorescence profiles.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "show_fl_al", value = FALSE)
    }else{
      if( any( val$file_types[files_to_plot()] == "csv_fluo") & !all( names(val$file_types[files_to_plot()] == "csv_fluo") %in% names(val$baseline_fl) ) & input$show_fl_al  ){
        showNotification("Some of your fluorescence profiles do not have a baseline.",
                         duration = NULL, type = "warning")
      }
    }
  })
  
  # show notification if how fluorescence was selected already but a new file in files_to_plot()
  # does not contain a fluorescence signal
  observeEvent(files_to_plot(), {
    if( !all( names(val$file_types[files_to_plot()] == "csv_fluo") %in% names(val$baseline_fl) ) & input$show_fl_al  ){
      showNotification("Some of your fluorescence profiles do not have a baseline.",
                       duration = NULL, type = "warning")
    }
  })
  
  
  
  # show notification when normalization to length or height is set but total area is missing
  observeEvent(input$normalize_length,{
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( files_to_plot() %in% files_with_total ) & input$normalize_length ){
      showNotification("Please select a total area for all profiles in the alignment.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "normalize_length", value = FALSE)
    }
  })
  
  # show notification when normalization to length or height is set but total area is missing
  observeEvent(input$normalize_height,{
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( files_to_plot() %in% files_with_total ) & input$normalize_height ){
      showNotification("Please select a total area for all profiles in the alignment.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "normalize_height", value = FALSE)
    }
  })
  
  # show notification when normalization to length or height is set but total area is missing
  observeEvent(files_to_plot(),{
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( files_to_plot() %in% files_with_total ) & input$normalize_length ){
      showNotification("Some profiles in the alignment do not have a total area.",
                       duration = NULL, type = "error")
      updateCheckboxInput(session, "normalize_length", value = FALSE)
    }
  })
  
  observeEvent(files_to_plot(),{
    files_with_total <- names(val$sum_areas[["Total"]])
    if(!all( files_to_plot() %in% files_with_total ) & input$normalize_height ){
      showNotification("Some profiles in the alignment do not have a total area.",
                       duration = NULL, type = "error")
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
    
    par(mar = c(0, 5, 0.5, 2)) 
    plot(xvalue(), fluorescence(), type = "l", xaxt = "n", col = "darkgreen", 
         xlim =c(xmin_single(),xmax_single()), ylim = c(ymin_single_fl(), ymax_single_fl()),
         las = 1, ylab = ylab, mgp = c(3.5, 0.8, 0), xlab = "", yaxt = "n", cex.lab = 2)
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_fl(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 2)
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      # quantified area is colored in green. If baseline, area_end or area_satrt lines are changed, the colored area stays the same for the selected quantified area 
      # until the button "quantify area" is pressed again.
      if(isTruthy(val$fl_control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, xvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], x_last),
                c(val$fl_control_baseline[[input$select_area]][input$select], fluorescence()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], val$fl_control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "darkgreen")
      }
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
    
    par(mar = c(5, 5, 0, 2))
    plot(xvalue(), yvalue(), type = "l", las = 1,
         ylab = ylab, xlab = "Time (min)",
         ylim = c(ymin_single(),ymax_single()), 
         xlim =c(xmin_single(),xmax_single()), mgp = c(3.5, 0.8, 0),
         yaxt = "n", xaxt = "n", cex.lab = 2
    )
    
    
    axis(1, las = 1, mgp = c(3.5, 1.2, 0), cex.axis = 2)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_pol(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 2)
    # show polygon automatically from start to stop value the user quantified
    # if the user selects a quantified area again, the respective polygon gets displayed in the plot again
    if(isTruthy(val$area_starts[[input$select_area]][input$select]) && isTruthy(val$area_ends[[input$select_area]][input$select]) & input$green_lines & isTruthy(val$baseline[input$select])){
      x_first <- round(val$area_starts[[input$select_area]][input$select])*val$factors_list[[input$select]]
      x_last <- round(val$area_ends[[input$select_area]][input$select])*val$factors_list[[input$select]]
      
      if(isTruthy(val$control_baseline[[input$select_area]][input$select])){
        polygon(c(x_first, xvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], x_last),
                c(val$control_baseline[[input$select_area]][input$select], yvalue()[(which(xvalue() == (x_first))):(which(xvalue() == (x_last)))], val$control_baseline[[input$select_area]][input$select]),
                col = "#c7e9c0", border = "black")
      }
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
      if(input$show_fl && ("SampleFluor" %in% colnames(file_plot()))){
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
    filename = "quantification.csv",
    content = function(file) {
      write.csv2( val$df_quant, file, row.names = F)
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
  
  # if color palette choice gets changed by the user, the original color_vector is reset to an empty list again
  observeEvent(
    input$color_palette,              
    {val$colors_collected = list()}, 
    ignoreInit = TRUE
  )
  
  # create color values, rainbow palette as initial colors, further replaced by selected color if selected
  colors_vector <- reactive({
    if(input$color_palette == "dark_palette"){
      dummy <- qualitative_hcl(length(files_to_plot()), palette = "Dark 3") 
    }
    if(input$color_palette == "rainbow_palette"){
      dummy <- rainbow(length(files_to_plot())) 
    }
    if(input$color_palette == "color_blind"){
      pal <- c("#000000", "#ff6db6", "#006ddb", "#920000",
               "#b66dff", "#004949", "#ffb6db", "#6db6ff", "#924900",
               "#24ff24", "#929200", "#490092", "#b6dbff", "#db6d00", "#ffff6d")
      dummy <- rep(pal, ceiling( length( files_to_plot() )/length(pal ) ) )
    }
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
    files_with_total <- names( val$sum_areas[["Total"]])
    for(f in files_to_plot()){
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
    for(f in files_to_plot()){
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
    if(any(val$file_types[files_to_plot()] == "csv_fluo") )
    {
      dummy <- list()
      for(f in names(val$baseline_fl)) # only go through fluorescence signals with baseline
      {
        start <- grep("Data Columns:", readLines(val$paths_collected[f]))
        fluo <- read.csv(val$paths_collected[f], skip = start)$SampleFluor
        dummy[[f]] <- smooth_profile( ( ( fluo - val$baseline_fl[[f]] )/norm_factor()[f])*norm_factor_x()[f], input$slider1)
      }
      dummy
    }else{  
      NULL
    }
  })
  
  
  values_list <- reactive({
    dummy <- list()
    for(f in files_to_plot())
    {
      if(val$factors_list[[f]] == (0.1/60) & grepl(".pks",as.character(f))){
        next_y <- read.table(val$paths_collected[f], dec = ",", header = F)$V3
      }else if(val$factors_list[[f]] == (4.8/15)/60){
        start <- grep("Data Columns:", readLines(val$paths_collected[f]))
        next_y <- read.csv(val$paths_collected[f], skip = start)$AbsA
      }else if(val$factors_list[[f]] == (0.2)/60){
        start <- grep("Data Columns:", readLines(val$paths_collected[f]))
        next_y <- head(read.csv(val$paths_collected[f], skip = start)$Absorbance, -1)
      }else if(val$factors_list[[f]] == (0.1/60) & grepl(".txt",as.character(f))){
        next_y <- read.delim(val$paths_collected[f])[,2]
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
  
  xmin_all <- reactive({
    min(aligned_starts())
  })
  
  # find maximum of y-axis values for polysome profiles
  ymin_all <- reactive({
    profile_mins <- sapply(values_list()[files_to_plot()], FUN = min)
    min(profile_mins )
  })
  
  # find maximum of y-axis values for fluorescence profiles
  ymin_all_fl <- reactive({
    profile_mins <- sapply(values_fluorescence()[files_to_plot()], FUN = min)
    min(profile_mins )
  })
  
  # axis labels are modified if there are more than three numerals
  lost_num_al_fl <- reactive({
    max_numeral <- max( floor(log10(abs( unlist(values_fluorescence()[intersect( files_to_plot(), names(val$baseline_fl) ) ] ) ))) + 1 ) # how many numerals has the maximum UV absorption?
    
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
  
  lost_num_al_Index <- reactive({
    max_numeral <- max( floor(log10(abs( aligned_ends() ))) + 1 ) # how many numerals has the maximum Index?
    
    if(max_numeral > 2 ){
      10^(max_numeral - 2)
    }else{
      1
    }
  })
  
  # When user selects other axis limits, limit values change 
  # (works for both single and aligned (+ normalized) plots)
  
  ymin <- reactive({
      input$axis1_a*lost_num_al_pol()
  })
  
  ymax <- reactive({
      input$axis2_a*lost_num_al_pol()
  })
  
  ymin_fl <- reactive({
      input$axis1_a_fl*lost_num_al_fl()
  })
  
  ymax_fl <- reactive({
      input$axis2_a_fl*lost_num_al_fl()
  })
  
  xmin <- reactive({
      input$axis3_a*lost_num_al_Index()
  })
  
  xmax <- reactive({
      input$axis4_a*lost_num_al_Index()
  }) 
  
  
  
  
  observe({
    req(val$files_to_align)
    updateNumericInput(session, "axis3_a", value = round( xmin_all() /lost_num_al_Index(), digits = 2 ) )
    updateNumericInput(session, "axis4_a", value = round( xmax_all() /lost_num_al_Index(), digits = 2 ) )
    updateNumericInput(session, "axis1_a", value = round( ymin_all()/lost_num_al_pol(), digits = 2 ) )
    updateNumericInput(session, "axis2_a", value = round( ymax_all()/lost_num_al_pol(), digits = 2 ) )
  })
  
  # let axis limits update for fluorescence axis
  observe({
    req(val$files_to_align)
    req(val$baseline_fl)
    req(input$show_fl_al)
    updateNumericInput(session, "axis1_a_fl", value = round( ymin_all_fl() /lost_num_al_fl(), digits = 2 ) )
    updateNumericInput(session, "axis2_a_fl", value = round( ymax_all_fl() /lost_num_al_fl(), digits = 2 ) )
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
    par(mar = c(5, 5, 0, 2)) 
    plot(x, values_list()[[f]], type = "l", lty = lines_vector()[f],
         lwd = linewidth_vector()[f],
         ylab = ylab, xlab = "Relative position", las = 1,
         col = colors_vector()[f], mgp = c(3.5, 1, 0), 
         ylim = c(ymin(),ymax()), xlim = c(xmin(),xmax()),
         yaxt = "n", xaxt = "n", cex.lab = 2
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
      anchor_line <- val$anchor[files_to_plot()[1]] / norm_factor_x()[files_to_plot()[1]] + shifts()[files_to_plot()[1]]
      abline(v = anchor_line, col = "red", lty=2)
    }
    
    a <- axTicks(1)
    axis(1, at = a, labels = a/lost_num_al_Index(), las = 1, mgp = c(3.5, 1.2, 0), cex.axis = 2)
    
    a <- axTicks(2)
    axis(2, at = a, labels = a/lost_num_al_pol(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 2)
    
    for(f in files_to_plot()[-1])
    {
      x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
      
      # create dataframe with x and y values to merge together in each for loop round. Col names are data names.
      y_aligned <- values_list()[[f]]
      x_aligned <- x
      df2 = list(x_aligned = x_aligned, y_aligned=y_aligned)
      attributes(df2) = list(names = names(df2),
                             row.names=1:max(length(x_aligned), length(y_aligned)), class='data.frame')
      colnames(df2) <- c("Index", as.character(str_remove(f, ".pks|.csv|.txt"))) # muss auch noch ".csv" removen kÃÂÃÂ¶nnen
      csv_file_df <- merge(csv_file_df, df2, by="Index", all = T)
      
      # plot files in alignment
      points(x, values_list()[[f]], type = "l", lty = lines_vector()[f], 
             lwd = linewidth_vector()[f], col = colors_vector()[f])
    }
    legend("topright", legend = files_to_plot(), lty = lines_vector()[files_to_plot()], 
           lwd = linewidth_vector()[files_to_plot()], col = colors_vector()[files_to_plot()],
           bty = "n", cex = 2
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
    
    files_in_al_fluo <- files_to_plot()[files_to_plot() %in% intersect(names(val$baseline_fl),  names(val$baseline))] # this is necessary to maintain the order as given by files_to_plot()
    f <-  files_in_al_fluo[1]
    x <- seq(aligned_starts()[f], aligned_ends()[f], by = 1/norm_factor_x()[f])
    par(mar = c(0, 5, 0.5, 2)) 
    plot(x, values_fluorescence()[[f]], type = "l", lty = lines_vector()[f],
         lwd = linewidth_vector()[f],
         ylab = ylab, xlab = "", las = 1,
         col = colors_vector()[f], mgp = c(3.5, 0.8, 0), 
         ylim = c(ymin_fl(),ymax_fl()), xlim = c(xmin(),xmax()),
         yaxt = "n", xaxt = "n", cex.lab = 2
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
    axis(2, at = a, labels = a/lost_num_al_fl(), las = 1, mgp = c(3.5, 0.8, 0), cex.axis = 2)
    
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
      
      points(x, values_fluorescence()[[f]], type = "l", lty = lines_vector()[f], 
             lwd = linewidth_vector()[f], col = colors_vector()[f])
    }
    
    #create reactive dataframe of all plots to have access outside of renderPlot function
    val$csv_file_df_all <- merge(val$csv_file_df, csv_file_df_fluo, by="Index", all = T)
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
    filename = "alignment.pdf",
    content = function(file) {
      pdf(file, width = 10, height = 6 )
      print( plot_alignedPol() )
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
  
}

###############
shinyApp(ui = ui, server = server)
