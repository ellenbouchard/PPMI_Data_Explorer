#
# PPMI Data Explorer
# A Shiny app
#
# Author: Ellen Bouchard
#
# Date Created: February 11 2025
#
# This is a Shiny app that allows for exploration of data sourced from PPMI 
# for the cohort of individuals included in the Ajami Lab PBMC sequencing study.
#
# NOTE: This is specifically a "Demo" version that obfuscates unnecessary information,
# Such as dates, patient IDs, and the patient metrics being reported.


require(shiny)
require(plotly)
require(rlang)
require(dplyr)

################  SETUP  ######################
# Upload datasets
setwd("~/Documents/OHSU/shiny/PPMI_Data_Explorer")
data_dir <- "./data_files_unnamed"
file_paths <- list.files(data_dir, full.names = TRUE)
names <- tools::file_path_sans_ext(basename(file_paths))
data_dfs <- lapply(file_paths, read.csv)
names(data_dfs) <- names

# Define valid list of datasets for user to choose from
datasets <- names

# Quick cleanup:
# Make sure Date (INFODT) is parsed correctly for all datasets
# Add INFODT_NUM column, which includes dates in numeric format, required for smoothing via LOESS
# Change SEX from binary (0, 1) to Female or Male for human readability
# Add GENOTYPE column in which WT is reported as "WT" (instead of 0) for human readability
# Add CONDITION column that includes "Control" and "PD" conditions
data_dfs <- lapply(data_dfs, function(df) {
  df$INFODT <- as.Date(df$INFODT, format = "%m/%d/%y")
  df$INFODT_NUM <- as.numeric(df$INFODT)
  df <- df %>%
    mutate(SEX = ifelse(SEX == 0, "Female", "Male"),
           GENOTYPE = ifelse(VAR_GENE == "0", "WT", VAR_GENE),
           CONDITION = ifelse(GROUP %in% c("HC", "LHC"), "Control", "PD"))
  return(df)
})

# Define list of valid types of plots for user to choose from
plottypes <- c("Line Graph", "Box Plot", "Bar Chart")

# Dictionary of data column names: 
# Includes human readable versions of column names for display purposes
colnames_dict <- c("Years Since Enrollment" = "DELTA_ENROLL", 
                   "Date" = "INFODT", 
                   "Age" = "AGE",
                   "Years Since Diagnosis" = "DELTA_DX", 
                   "Years Since Symptom Onset" = "DELTA_SX",
                   "Group" = "GROUP", 
                   "Condition" = "CONDITION", 
                   "LRRK2 Genotype" = "GENOTYPE", 
                   "Sex" = "SEX", 
                   "APOE Genotype" = "APOE",
                   "Normalized Score" = "SCORE_NORM")
# Make a reversed one for plot labeling purposes
colnames_lookup <- setNames(names(colnames_dict), colnames_dict)
# Break down into which columns users can use for different types of plots
xaxis_options_linegraph <- colnames_dict[1:5]
group_options_linegraph <- colnames_dict[6:10]
xaxis_options_boxplot <- colnames_dict[6:10]
group_options_boxplot <- colnames_dict[6:10]
xaxis_options_barchart <- colnames_dict[6:10]
group_options_barchart <- colnames_dict[6:10]

# Define color codes for groups
# NOTE: This is manually hard-coded for now
# I intend to add in the option for the user to define custom colors
group_colors <- c(
  "HC" = "#2b3a6b",
  "LHC" = "#5167ad",
  "iPD" = "#b56c1f",
  "LPD" = "#f29735",
  "Control" = "#5167ad",
  "PD" = "#f29735",
  "WT" = "#5167ad",
  "LRRK2" = "#f29735",
  "PRKN" = "#4f4f4f",
  "Female" = "#f29735",
  "Male" = "#5167ad",
  "E2/E2" = "#2b3a6b",
  "E2/E3" = "#5167ad",
  "E2/E4" = "#6e2121",
  "E3/E3" = "#f29735",
  "E3/E4" = "#b32e2e",
  "E4/E4" = "#e84d4d"
)

# Define shapes for representation of SEX 
sex_shapes <- c("Female" = "circle", "Male" = "square")

############### PLOT FUNCTION ###################
# Function to generate Plotly plot
# (No defaults since defaults are specified in the app)
# ARGUMENTS:
#   data = R dataframe chosen from "datasets" list
#   type = string indicating type of chart; options are "Line Graph", "Box Plot", or "Bar Chart"
#   xaxis = string indicating column of dataframe to use as X axis variable 
#   groupinput = string indicating column of dataframe to use as color grouping variable
#   filter_closest = boolean used to determine whether data will be filtered to only include patient visits closest in time to the date of sample collection
#     (only used for boxplots)
#   yaxis_bar = string indicating column of dataframe to use as Y axis variable (bar charts only)
#   plot_title = string to be used as plot title
# OUTPUT: 
#   plot = a plotly object
make_plot_plotly <- function(data,
                             type,
                             xaxis,
                             groupinput,
                             filter_closest,
                             yaxis_bar,
                             plot_title) {
  
  # If user chooses Date (column name "INFODT") as x axis variable, must be changed to numeric date (column name INFODT_NUM) to allow for smoothing calculations
  if(xaxis == "INFODT") {
    xaxis_internal <- "INFODT_NUM"
  } else {
    xaxis_internal <- xaxis
  }
  
  # Turn input variables from strings into symbols so that they can be recognized as column names
  xaxissym <- sym(xaxis_internal) # In the case of plotting against DATE, is equal to INFODT_NUM (numeric format)
  xaxisdisp <- sym(xaxis) # In the case of plotting against DATE, is equal to INFODT (date format)
  groupsym <- sym(groupinput)
  
  # Filter to remove any rows that have blank or NA values for xaxis
  # (This is specifically to account for control patients who don't have a diagnosis or symptom onset date)
  data <- data %>%
    filter(!is.na(!!xaxissym), !!xaxissym != "")
  
  # If user chooses DELTA_SX, remove patient who is outlier for symptom onset date 
  # (>30 years between symptom onset and diagnosis)
  if(xaxis == "DELTA_SX") {
    data <- data %>%
      group_by(PATNO) %>%
      filter(!any(DELTA_SX > 30, na.rm = TRUE)) %>%
      ungroup()
  }
  
  # If user is making a boxplot AND "filter for visits closest to collection date" box is checked, 
  #filter the data to only include rows for which "CLOSEST_TO_COLLECTION" is equal to 1
  if(type == "Box Plot" && filter_closest) {
    data <- data %>%
      filter(CLOSEST_TO_COLLECTION == 1)
  }
  
  # BAR CHART
  if(type == "Bar Chart") {
    yaxis_title = "Number of Individuals"
    plot_title = NULL
    
    # If the xaxis groupings and color grouping variables are the same,
    # do not include a color variable 
    if(xaxis == groupinput) {
      # Prep data by counting number of individuals per group
      bar_data <- data %>%
        group_by(!!xaxissym) %>%
        summarise(N_PATNOS = n_distinct(PATNO), .groups = "drop")
      
      # Make bar chart with no color option
      plot <- plot_ly(data = bar_data) %>%
        add_bars(x = as.formula(paste0("~",xaxis)),
                 y = ~N_PATNOS,
                 type = "bar",
                 marker = list(color = "#2b3a6b"))
      
      # If x axis and grouping variables are different, include color variable
      } else {
      # First count number of individuals per GROUP category (for normalization purposes)
      totals <- data %>%
        group_by(!!groupsym) %>%
        summarise(TOTAL_PATNOS_PERGROUP = n_distinct(PATNO), .groups = "drop")
        
      # Prep data by counting number of individuals per group
      # Calculate proportion and normalized proportion
      bar_data <- data %>%
        group_by(!!xaxissym, !!groupsym) %>%
        summarise(N_PATNOS = n_distinct(PATNO), .groups = "drop") %>%
        group_by(!!xaxissym) %>%
        mutate(PROP = N_PATNOS / sum(N_PATNOS)) %>%
        left_join(totals, by = groupinput) %>%
        mutate(NORM_FACTOR = 1 / TOTAL_PATNOS_PERGROUP,
               NORM_N = N_PATNOS * NORM_FACTOR,
               NORM_PROP = NORM_N / sum(NORM_N))
      
      # If Y axis is set to "number":
      if(yaxis_bar == "Number") {
        plot <- plot_ly(data = bar_data) %>%
          add_bars(x = as.formula(paste0("~",xaxis)),
                   y = ~N_PATNOS,
                   color = as.formula(paste0("~",groupinput)),
                   colors = group_colors,
                   type = "bar") %>%
          layout(barmode = "stack")
        
      # If Y axis is set to "proportion":
      } else if(yaxis_bar == "Proportion") {
        yaxis_title = "Proportion of Individuals"
        plot <- plot_ly(data = bar_data) %>%
          add_bars(x = as.formula(paste0("~",xaxis)),
                   y = ~PROP,
                   color = as.formula(paste0("~",groupinput)),
                   colors = group_colors,
                   type = "bar") %>%
          layout(barmode = "stack")
        
        # If Y axis is set to "normalized proportion" :
      } else if(yaxis_bar == "Normalized Proportion") {
        yaxis_title = "Normalized Proportion of Individuals"
        plot <- plot_ly(data = bar_data) %>%
          add_bars(x = as.formula(paste0("~",xaxis)),
                   y = ~NORM_PROP,
                   color = as.formula(paste0("~",groupinput)),
                   colors = group_colors,
                   type = "bar") %>%
          layout(barmode = "stack")
      }
    }
  }
  
  # BOX PLOT
  if(type == "Box Plot") {
    yaxis_title = "Normalized Score"
    
    # Generate boxplot
    plot <- plot_ly(data = data,
                    y = ~SCORE_NORM) %>%
      add_boxplot(x = as.formula(paste0("~",xaxis)),
                  color = as.formula(paste0("~",groupinput)),
                  colors = group_colors,
                  boxpoints = "all",
                  jitter = 0.5,           
                  pointpos = 0,   
                  marker = list(opacity = 0.6, size = 6),
                  text = ~paste(
                    "Patient:", " ",
                    "<br>Condition:", CONDITION,
                    "<br>Genotype:", GENOTYPE,
                    "<br>APOE:", APOE,
                    "<br>Sex:", SEX,
                    "<br>Age:", " ",
                    "<br>nScore:", round(SCORE_NORM, 1)
                  ),
                  hoverinfo = "text") 
    
    # If there are different inputs for x axis and group, alter the box arrangement to be grouped
    # (as oppossed to overlapping)
    if(xaxis != groupinput) {
      plot <- plot %>%
        layout(boxmode = "group")
    }
  }
  
  # LINE GRAPH
  if(type == "Line Graph") {
    yaxis_title = "Normalized Score"
    
    
  # Plotly doesn't have a "geom_smooth" equivalent that automatically generates loess curves,
  # So I need to manually calculate those as a new dataframe.
  # First, calculate average score per x axis category (for smoothing calculation purposes)
  summary_data <- data %>%
    rename(x = !!xaxissym) %>%
    rename(group = !!groupsym) %>%
    group_by(group, x) %>%
    summarise(
      AVGSCORE = mean(SCORE_NORM, na.rm = TRUE),
      .groups = "drop"
    )
  # Then calculate loess smoothing values based on averages 
  loess_data <- summary_data %>%
    group_by(group) %>%
    do({
      loess_model <- loess(AVGSCORE ~ x, data = ., span = 0.3)
      prediction <- predict(loess_model)
      data.frame(x = .$x, group = .$group, SMOOTHED = prediction)
    }) %>%
    ungroup() 
  
  # When plotting dates, need to change dates from numeric back to date format
  # for display purposes
  if(xaxis == "INFODT") {
    loess_data <- loess_data %>%
      mutate(x = as.Date(x, origin = "1970-01-01"))
  }
    
    # Generate plot, adding point then line traces
  plot <- plot_ly(data = data,
                  y = ~SCORE_NORM) %>%
      add_trace(data = data,
                x = as.formula(paste0("~", xaxisdisp)),
                y = ~SCORE_NORM,
                color = as.formula(paste0("~", groupsym)),
                colors = group_colors,
                symbols = sex_shapes,
                symbol = ~SEX,
                type = "scatter",
                mode = "markers",
                marker = list(size = 7),
                showlegend = TRUE,
                alpha = 0.5,
                text = ~paste("Patient:", " ",
                              "<br>Condition:", CONDITION,
                              "<br>Genotype:", GENOTYPE,
                              "<br>APOE:", APOE,
                              "<br>Sex:", SEX,
                              "<br>Age:", " ",
                              "<br>nScore:", round(SCORE_NORM, 1)),
                hoverinfo = "text") %>%
      add_trace(data = loess_data,
                x = ~x, 
                y = ~SMOOTHED, 
                mode = "lines", 
                type = "scatter",
                color = ~group,
                colors = group_colors,
                name = ~group,
                showlegend = FALSE
                )
  }
  # Finalize layout
  plot <- plot %>% layout(title = plot_title, 
                          legend = list(
                            title = list(text = colnames_lookup[[groupinput]])
                          ),
                          xaxis = list(
                            title = colnames_lookup[[xaxis]],
                            showgrid = FALSE,
                            zeroline = FALSE 
                          ),
                          yaxis = list(
                            title = list(text = yaxis_title),
                            showgrid = FALSE,
                            zeroline = FALSE  
                          ),
                          dragmode = FALSE,
                          margin = list(
                            l = 70,
                            r = 70,
                            t = 70,
                            b = 70
                          ))
  return(plot)
}

###################  UI  #######################
ui <- fluidPage(

    # App Title
    titlePanel("PPMI Data Explorer"),

    # Sidebar with input
    sidebarLayout(
        sidebarPanel(
          selectInput("plottype", "Select Plot Type", plottypes, selected = "Line Graph"),
          tabsetPanel(
            id = "parameters",
            type = "hidden",
            tabPanel("Line Graph",
              selectInput("dataset", "Select Patient Outcome", datasets),
              selectInput("xaxis_line", "X axis", choices = xaxis_options_linegraph),
              selectInput("groupings_line", "Group by", choices = group_options_linegraph)
            ),
            tabPanel("Box Plot",
              selectInput("dataset", "Select Patient Outcome", datasets),
              selectInput("xaxis_box", "X axis", choices = xaxis_options_boxplot),
              selectInput("groupings_box", "Group by", choices = group_options_boxplot),
              checkboxInput("filter_closest", "Only Show Patient Visits Closest to Date of Sample Collection", value = TRUE)
            ),
            tabPanel("Bar Chart",
                     selectInput("xaxis_bar", "X axis", choices = xaxis_options_barchart),
                     selectInput("groupings_bar", "Group by", choices = group_options_barchart),
                     selectInput("yaxis_bar", "Y axis", choices = c("Number", "Proportion", "Normalized Proportion"))
            ),
          ),
        ),
        # Show plot and data table
        mainPanel(
          plotlyOutput("plot"),
          DT::DTOutput("dataTable")
        )
    )
)

##################  SERVER LOGIC  ######################
server <- function(input, output, session) {
  
  # Define input variables as reactive values
  dataset <- reactive(input$dataset)
  plottype <- reactive(input$plottype)
  filter_closest <- reactive(input$filter_closest)
  yaxis_bar <- reactive(input$yaxis_bar)
  
  # Define xaxis and grouping input variables as reactive values based on plot type
  xaxis <- reactive({
    if (input$plottype == "Line Graph") input$xaxis_line
    else if (input$plottype == "Box Plot") input$xaxis_box
    else if(input$plottype == "Bar Chart") input$xaxis_bar
  })
  groupings <- reactive({
    if (input$plottype == "Line Graph") input$groupings_line
    else if (input$plottype == "Box Plot") input$groupings_box
    else if (input$plottype == "Bar Chart") input$groupings_bar
  })
  
  # Define data based on user input
  data <- reactive({
    req(input$dataset)
    df <- data_dfs[[dataset()]]
    return(df)
  })
  
  # JUST FOR DEMO: make "display" dataframe that does not show unnecessary information 
  data_demo <- reactive({
    req(input$dataset)
    df <- data_dfs[[dataset()]] %>%
      mutate(PATNO = "PatientID",
             PATNO_VISIT = "Patient Visit",
             INFODT = "Visit Date",
             COLLECTION_DATE = "Collection Date",
             ENROLL_DATE = "Enrollment Date",
             SXDT = "Symptom Onset Date",
             PDDXDT = "Diagnosis Date") %>%
      select(-c(LABEL, SCORE, BARCODE, ST_NUMBER, GROUP_CODE, BIRTHDT, INFODT_NUM, COHORT_DEFINITION,
                AGE_SX, AGE_DX, HYPEREXPANDED,DELTA_COLLECT, DELTA_STATUS, ENROLL_AGE, ENROLL_STATUS,
                STATUS_DATE, AGE_COLLECTION))
    return(df)
  })
  
  # Observe input plot type
  observe({
    req(input$plottype)
    updateTabsetPanel(session, "parameters", selected = input$plottype)
  })
  
  # Generate plot using make_plot function 
  make_plot <- reactive({
    p <- NULL
    if(!is.null(data())) {
      p <- make_plot_plotly(
        data(),
        type = plottype(),
        plot_title = dataset(),
        xaxis = xaxis(),
        groupinput = groupings(),
        filter_closest = filter_closest(),
        yaxis_bar = yaxis_bar())
    }
    return(p)
  })
  
  
  # Render the plot
  output$plot <- renderPlotly({
      make_plot()
    })
  
  # Render data table 
  output$dataTable <- DT::renderDT(data_demo(), options = list(pageLength = 5))
}  

##################  Run the application  ##################
shinyApp(ui = ui, server = server)
