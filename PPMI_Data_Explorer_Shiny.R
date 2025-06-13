#
# PPMI Data Explorer:
# A Shiny app
#
# Author: Ellen Bouchard
#
# Date Created: February 11 2025
#
# This is a Shiny app that allows for exploration of data sourced from PPMI 
# for the cohort of individuals included in the Ajami Lab PBMC sequencing study.


require(shiny)
require(plotly)
require(rlang)
require(dplyr)
require(broom)
require(purrr)
require(emmeans)

################  SETUP  ######################
# Set working directory
setwd("~/Documents/OHSU/shiny/PPMI_Data_Explorer")

# Upload molecular (myeloid/lymphoid percents) data for correlation graphs
molecular_data <- read.csv("./myeloid_lymphoid_pcts.csv")
# Get column names from molecular data to use for correlation plots
molecular_colnames <- colnames(molecular_data)
molecular_names <- molecular_colnames[4:6]
#celltypes <- unique(gene_module_data$CellType)

# Upload datasets of patient outcomes - point in time 
# For correlation graphs
# NOTE: These have the same names as the longitudinal outcome datasets
corr_data_dir <- "./correlation_files"
corr_file_paths <- list.files(corr_data_dir, full.names = TRUE)
corr_names <- tools::file_path_sans_ext(basename(corr_file_paths))
corr_data_dfs <- lapply(corr_file_paths, read.csv)
names(corr_data_dfs) <- corr_names

# Upload datasets of patient outcomes - Longitudinal
data_dir <- "./data_files"
file_paths <- list.files(data_dir, full.names = TRUE)
names <- tools::file_path_sans_ext(basename(file_paths))
data_dfs <- lapply(file_paths, read.csv)
names(data_dfs) <- names

# Define valid list of datasets for user to choose from
datasets <- names

# Quick cleanup/preprocessing:
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

# Preprocessing of correlation data:
corr_data_dfs <- lapply(corr_data_dfs, function(df) {
  df <- df %>%
    mutate(SEX = ifelse(SEX == 0, "Female", "Male"),
           GENOTYPE = ifelse(VAR_GENE == "0", "WT", VAR_GENE),
           CONDITION = ifelse(GROUP %in% c("HC", "LHC"), "Control", "PD"))
  return(df)
})

# Define list of valid types of plots for user to choose from
plottypes <- c("Line Graph", "Box Plot", "Bar Chart", "Correlation")

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
                   "Normalized Score" = "SCORE_NORM",
                   "Estimated Score at Time of Collection" = "EST_SCORE_AT_COLLECTION",
                   "Rate of Change of Score" = "SLOPE",
                   "Age at Symptom Onset" = "AGE_SX",
                   "Age at Diagnosis" = "AGE_DX",
                   "Percent Lymphoid" = "PCT_LYMPHOID",
                   "Percent Myeloid" = "PCT_MYELOID")
# Make a reversed one for plot labeling purposes
colnames_lookup <- setNames(names(colnames_dict), colnames_dict)
# Break down into which columns users can use for different types of plots
xaxis_options_linegraph <- colnames_dict[1:5]
group_options_linegraph <- c(colnames_dict[6:10], "None")
xaxis_options_boxplot <- colnames_dict[6:10]
group_options_boxplot <- c(colnames_dict[6:10])
yaxis_options_boxplot <- colnames_dict[12:15]
xaxis_options_barchart <- colnames_dict[6:10]
group_options_barchart <- c(colnames_dict[6:10])
yaxis_options_corr <- colnames_dict[12:15]

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
  "E4/E4" = "#e84d4d",
  "All" = "black"
)

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
#   yaxis = string indicating column of dataframe to use as Y axis variable
#   plot_title = string to be used as plot title
# OUTPUT: 
#   plot = a plotly object
make_plot_plotly <- function(data,
                             molecular_data,
                             type,
                             xaxis,
                             groupinput,
                             filter_closest,
                             yaxis,
                             celltype,
                             plot_title) {
  # If user chooses Date (column name "INFODT") as x axis variable, 
  # must be changed to numeric date (column name INFODT_NUM) to allow for smoothing calculations
  if(xaxis == "INFODT") {
    xaxis_internal <- "INFODT_NUM"
  } else {
    xaxis_internal <- xaxis
  }
  
  # Turn input variables from strings into symbols so that they can be recognized as column names
  xaxissym <- sym(xaxis_internal) # In the case of plotting against DATE, is equal to INFODT_NUM (numeric format)
  xaxisdisp <- sym(xaxis) # In the case of plotting against DATE, is equal to INFODT (date format)
  if(!is.null(groupinput) && groupinput != "None") groupsym <- sym(groupinput)
  
  # Filter to remove any rows that have blank or NA values for xaxis
  # (This is specifically to account for control patients who don't have a diagnosis or symptom onset date)
  if(type != "Correlation") {
    data <- data %>%
      filter(!is.na(!!xaxissym), !!xaxissym != "")
  }
  
  # If user chooses DELTA_SX, remove patient who is outlier for symptom onset date 
  # (>30 years between symptom onset and diagnosis)
  if(xaxis == "DELTA_SX") {
    data <- data %>%
      group_by(PATNO) %>%
      filter(!any(DELTA_SX > 30, na.rm = TRUE)) %>%
      ungroup()
  }
  
  # BAR CHART
  if(type == "Bar Chart") {
    yaxis_title = "Number of Individuals"
    xaxis_title = colnames_lookup[[xaxis]]
    legend_title = colnames_lookup[[groupinput]]
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
      if(yaxis == "Number") {
        plot <- plot_ly(data = bar_data) %>%
          add_bars(x = as.formula(paste0("~",xaxis)),
                   y = ~N_PATNOS,
                   color = as.formula(paste0("~",groupinput)),
                   colors = group_colors,
                   type = "bar") %>%
          layout(barmode = "stack")
        
      # If Y axis is set to "proportion":
      } else if(yaxis == "Proportion") {
        yaxis_title = "Proportion of Individuals"
        plot <- plot_ly(data = bar_data) %>%
          add_bars(x = as.formula(paste0("~",xaxis)),
                   y = ~PROP,
                   color = as.formula(paste0("~",groupinput)),
                   colors = group_colors,
                   type = "bar") %>%
          layout(barmode = "stack")
        
        # If Y axis is set to "normalized proportion" :
      } else if(yaxis == "Normalized Proportion") {
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
    yaxis_title = colnames_lookup[[yaxis]]
    xaxis_title = colnames_lookup[[xaxis]]
    legend_title = colnames_lookup[[groupinput]]
    
    # Generate boxplot
    plot <- plot_ly(data = data) %>%
      add_boxplot(x = as.formula(paste0("~",xaxis)),
                  y = as.formula(paste0("~",yaxis)),
                  color = as.formula(paste0("~",groupinput)),
                  colors = group_colors,
                  boxpoints = "all",
                  jitter = 0.5,           
                  pointpos = 0,   
                  marker = list(opacity = 0.6, size = 6),
                  text = ~paste(
                    "Patient:", PATNO,
                    "<br>Condition:", CONDITION,
                    "<br>Genotype:", GENOTYPE,
                    "<br>APOE:", APOE,
                    "<br>Sex:", SEX,
                    "<br>Age:", round(AGE_COLLECTION, 0)
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
    yaxis_title = "Score"
    xaxis_title = colnames_lookup[[xaxis]]
    legend_title = colnames_lookup[[groupinput]]
    
  # Plotly doesn't have a "geom_smooth" equivalent that automatically generates loess curves,
  # So I need to manually calculate those as a new dataframe.
  # First, calculate average score per x axis category (for smoothing calculation purposes)
  summary_data <- data %>%
    rename(x = !!xaxissym) %>%
    rename(group = !!groupsym) %>%
    group_by(group, x) %>%
    summarise(
      AVGSCORE = mean(SCORE, na.rm = TRUE),
      .groups = "drop"
    )
  # Then calculate loess smoothing values based on averages 
  loess_data <- summary_data %>%
    group_by(group) %>%
    do({
      loess_model <- loess(AVGSCORE ~ x, data = ., span = 0.6)
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
                  y = ~SCORE) %>%
      add_trace(data = data,
                x = as.formula(paste0("~", xaxisdisp)),
                y = ~SCORE,
                color = as.formula(paste0("~", groupsym)),
                colors = group_colors,
                type = "scatter",
                mode = "markers",
                marker = list(size = 7),
                showlegend = TRUE,
                alpha = 0.5,
                text = ~paste("Patient:", PATNO,
                              "<br>Condition:", CONDITION,
                              "<br>Genotype:", GENOTYPE,
                              "<br>APOE:", APOE,
                              "<br>Sex:", SEX,
                              "<br>Age:", round(AGE, 0),
                              "<br>Score:", round(SCORE, 1)),
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
  
  # CORRELATION PLOT
  if(type == "Correlation") {
    if (groupinput != "None") {
      groupsym <- sym(groupinput)
      legend_title <- colnames_lookup[[groupinput]]
    } else {
      legend_title <- NULL
    }
    
    yaxis_title = colnames_lookup[[yaxis]]
    xaxis_title = colnames_lookup[[xaxis]]
    plot_title = paste0(plot_title, ": ", xaxis)
    
    # For this one, convert y axis to symbol (only graph type for which Y axis is a dataframe variable)
    yaxissym <- sym(yaxis)
    
    # Filter to only include cell type of interest in expression data
    #expression_data <- expression_data %>%
    #  filter(CellType == celltype)
    
    # Add molecular data
    data <- data %>%
      left_join(molecular_data, by = "LABEL")
    
    # Create dataset that fits line to data for plotting
    # First, calculate average score per x axis category (for line calculation purposes)
    summary_data <- data %>%
      rename(x = !!xaxissym,
             y = !!yaxissym) %>%
      {
        if (groupinput != "None") {
          rename(., group = !!groupsym) %>%
            group_by(group, x)
        } else {
          mutate(., group = "All") %>%
            group_by(group, x)
        }
      } %>%
      summarise(AVGSCORE = mean(y, na.rm = TRUE), .groups = "drop")

    # Then calculate linear regression values based on averages 
    # NOTE: need to filter to remove groups that only have NA values
    lm_data <- summary_data %>%
      group_by(group) %>%
      filter(n_distinct(x[!is.na(x)]) >= 2,
             sum(!is.na(AVGSCORE)) >= 2) %>%
      do({
        lm_model <- lm(AVGSCORE ~ x, data = .)
        prediction <- predict(lm_model,  newdata = data.frame(x = .$x))
        data.frame(x = .$x, group = .$group, FITTED = prediction)
      }) %>%
      ungroup()
    
    # Calculate spearmann rank correlation coefficient and P value
    #cor_result <- cor.test(pull(data, !!xaxissym), data$EST_SCORE_AT_COLLECTION, method = "spearman")
    #spearman_rho <- round(cor_result$estimate, 2)
    #p_value <- signif(cor_result$p.value, 2)

    plot <- plot_ly(data = data,
                    y = ~EST_SCORE_AT_COLLECTION) %>%
      # Display points
      add_trace(data = data,
                x = as.formula(paste0("~", xaxisdisp)),
                y = as.formula(paste0("~", yaxissym)),
                color = if (groupinput != "None") as.formula(paste0("~", groupsym)) else NULL,
                colors = group_colors,
                type = "scatter",
                mode = "markers",
                marker = list(size = 7),
                showlegend = TRUE,
                alpha = 0.5,
                text = ~paste("Patient:", PATNO,
                              "<br>Condition:", CONDITION,
                              "<br>Genotype:", GENOTYPE,
                              "<br>APOE:", APOE,
                              "<br>Sex:", SEX,
                              "<br>Age at Collection:", round(AGE_COLLECTION, 0)),
                hoverinfo = "text") %>%
      # Display loess smoothed line on plot
      add_trace(data = lm_data,
                x = ~x, 
                y = ~FITTED, 
                color = ~group,
               colors = group_colors,
                mode = "lines", 
                type = "scatter",
                showlegend = FALSE)# %>%
      # Display Spearman's score on plot
     # add_annotations(
      #  text = paste0("Spearman's r = ", spearman_rho, "<br>p = ", p_value),
      # xref = "paper", yref = "paper",
       # x = 0.95, y = 1.1,
       # showarrow = FALSE,
       # font = list(size = 12),
       # bgcolor = 'rgba(255, 255, 255, 0.5)')
  }
  
  
  # Finalize layout
  plot <- plot %>% layout(title = plot_title, 
                          legend = list(
                            title = list(text = legend_title)
                          ),
                          xaxis = list(
                            title = list(text = xaxis_title),
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

# Function to generate results dataframe
# (No defaults since defaults are specified in the app)
# ARGUMENTS:
#   data = R dataframe chosen from "datasets" list
#   type = string indicating type of chart; options are "Line Graph", "Box Plot", or "Bar Chart"
#   xaxis = string indicating column of dataframe to use as X axis variable 
#   groupinput = string indicating column of dataframe to use as color grouping variable
#   yaxis = string indicating column of dataframe to use as Y axis variable
# OUTPUT: 
#   dataframe = a dataframe of results and statistics
make_dataframe <- function(data,
                             molecular_data,
                             type,
                             xaxis,
                             groupinput,
                             yaxis,
                             celltype) {
  # If user chooses Date (column name "INFODT") as x axis variable, 
  # must be changed to numeric date (column name INFODT_NUM) to allow for smoothing calculations
  if(xaxis == "INFODT") {
    xaxis_internal <- "INFODT_NUM"
  } else {
    xaxis_internal <- xaxis
  }
  
  # Turn input variables from strings into symbols so that they can be recognized as column names
  xaxissym <- sym(xaxis_internal) # In the case of plotting against DATE, is equal to INFODT_NUM (numeric format)
  xaxisdisp <- sym(xaxis) # In the case of plotting against DATE, is equal to INFODT (date format)
  if(!is.null(groupinput) && groupinput != "None") groupsym <- sym(groupinput)
  if(!is.null(yaxis)) yaxissym <- sym(yaxis)
  
  # Filter to remove any rows that have blank or NA values for xaxis
  # (This is specifically to account for control patients who don't have a diagnosis or symptom onset date)
  if(type != "Correlation") {
    data <- data %>%
      filter(!is.na(!!xaxissym), !!xaxissym != "")
  }
  
  # If user chooses DELTA_SX, remove patient who is outlier for symptom onset date 
  # (>30 years between symptom onset and diagnosis)
  if(xaxis == "DELTA_SX") {
    data <- data %>%
      group_by(PATNO) %>%
      filter(!any(DELTA_SX > 30, na.rm = TRUE)) %>%
      ungroup()
  }
  
  # BAR CHART
  if(type == "Bar Chart") {
    
    # If the xaxis groupings and color grouping variables are the same,
    # do not include a color variable 
    if(xaxis == groupinput) {
      # Prep data by counting number of individuals per group
      results_df <- data %>%
        group_by(!!xaxissym) %>%
        summarise(N_PATNOS = n_distinct(PATNO), .groups = "drop") %>%
        mutate(Proportion_of_total = round(N_PATNOS / sum(N_PATNOS), 2)) %>%
        rename(Number = N_PATNOS)
      
      # If x axis and grouping variables are different, include color variable
    } else {
      # First count number of individuals per GROUP category (for normalization purposes)
      totals <- data %>%
        group_by(!!groupsym) %>%
        summarise(TOTAL_PATNOS_PERGROUP = n_distinct(PATNO), .groups = "drop")
      
      # Prep data by counting number of individuals per group
      # Calculate proportion and normalized proportion
      results_df <- data %>%
        group_by(!!xaxissym, !!groupsym) %>%
        summarise(N_PATNOS = n_distinct(PATNO), .groups = "drop") %>%
        group_by(!!xaxissym) %>%
        mutate(PROP = round(N_PATNOS / sum(N_PATNOS), 2)) %>%
        left_join(totals, by = groupinput) %>%
        mutate(NORM_FACTOR = 1 / TOTAL_PATNOS_PERGROUP,
               NORM_N = round(N_PATNOS * NORM_FACTOR, 2),
               NORM_PROP = round(NORM_N / sum(NORM_N),2)) %>%
        select(!!xaxissym,
               !!groupsym,
               N_PATNOS,
               PROP,
               NORM_PROP) %>%
        rename(Number = N_PATNOS,
               Proportion_Of_Group = PROP,
               Normalized_Proportion_Of_Group = NORM_PROP
               )
    }
  return(results_df)
  }
  
  # BOX PLOT
  if(type == "Box Plot") {
  
    # If the "group" variable is different from the "xaxis" variable, 
    # then run pairwise tests between each group for each xaxis category
    if (groupinput != xaxis) {
      # Split data by x axis category
      results_df <- data %>%
        group_by(!!xaxissym) %>%
        group_split() %>%
        purrr::map_dfr(function(subdata) {
          # Remove groups with fewer than 2 values
          groups_use <- subdata %>%
            group_by(!!groupsym) %>%
            summarise(n = sum(!is.na(!!yaxissym)), .groups = "drop") %>%
            filter(n >= 2) %>%
            pull(!!groupsym)
          
          subdata_filtered <- subdata %>%
            filter((!!groupsym) %in% groups_use)
          
          # If there are not enough groups in the category to compare, return NULL
          if (length(unique(subdata_filtered[[groupinput]])) < 2) {
            return(NULL)
          }
          
          # Run pairwise t-tests
          ttest <- pairwise.t.test(
            x = pull(subdata_filtered, !!yaxissym),
            g = pull(subdata_filtered, !!groupsym),
            p.adjust.method = "fdr"
          )
          
          # Create results df
          out_df <- as.data.frame(as.table(ttest$p.value)) %>%
            filter(!is.na(Freq)) %>%
            rename(Group1 = Var1,
                   Group2 = Var2,
                   FDR_pValue_num = Freq) %>%
            mutate(
              FDR_pValue = formatC(FDR_pValue_num, format = "e", digits = 2),
              Significant = FDR_pValue_num < 0.05,
              X_Axis = unique(subdata[[xaxis]])
            ) %>%
            select(X_Axis, Group1, Group2, FDR_pValue, Significant)
          return(out_df)
        })
    } else {
      # If xaxis is the same as group, run once on full dataset
      # Remove groups with fewer than 2 values 
      groups_use <- data %>%
        group_by(!!groupsym) %>%
        summarise(n_present = sum(!is.na(!!yaxissym)), .groups = "drop") %>%
        filter(n_present >= 2) %>%
        pull(!!groupsym)
      
      data <- data %>%
        filter((!!groupsym) %in% groups_use)
      
      # Run pairwise t tests between all groups
      ttest_results <- pairwise.t.test(x = data[[yaxis]], 
                                       g = data[[groupinput]], 
                                       p.adjust.method = "fdr")
      
      # Format for result dataframe
      results_df <- as.data.frame(as.table(ttest_results$p.value)) %>%
        filter(!is.na(Freq)) %>%
        rename(Group1 = Var1, 
               Group2 = Var2,
               FDR_pValue_num = Freq) %>%
        mutate(
          FDR_pValue = formatC(FDR_pValue_num, format = "e", digits = 2),
          Significant = FDR_pValue_num < 0.05
        ) %>%
        select(Group1, Group2, FDR_pValue, Significant)
    }
  }
  
  # LINE GRAPH
  if(type == "Line Graph") {
    print("DEBUG: DATAFRAME FOR LINE GRAPH")
    
    # Fit linear model with interaction
    formula <- paste("SCORE ~", xaxis, "*", groupinput)
    model <- lm(as.formula(formula), data = data)
    contrast_formula <- as.formula(paste("pairwise ~", groupinput))
    results <- emtrends(model, specs = contrast_formula, var = xaxis)
    
    results_df <- as.data.frame(summary(results$contrasts)) %>%
      mutate(
        pValue = formatC(p.value, format = "e", digits = 2),
        Significant = p.value < 0.05,
        Slope_difference = round(estimate, 2), 
        Standard_Error = round(SE, 2)
      ) %>%
      select(contrast, Slope_difference, Standard_Error, pValue, Significant)
    
    
    # Plotly doesn't have a "geom_smooth" equivalent that automatically generates loess curves,
    # So I need to manually calculate those as a new dataframe.
    # First, calculate average score per x axis category (for smoothing calculation purposes)
    summary_data <- data %>%
      rename(x = !!xaxissym) %>%
      rename(group = !!groupsym) %>%
      group_by(group, x) %>%
      summarise(
        AVGSCORE = mean(SCORE, na.rm = TRUE),
        .groups = "drop"
      )
    # Then calculate loess smoothing values based on averages 
    loess_data <- summary_data %>%
      group_by(group) %>%
      do({
        loess_model <- loess(AVGSCORE ~ x, data = ., span = 0.6)
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
  }
  
  # CORRELATION PLOT
  if(type == "Correlation") {
    if (groupinput != "None") {
      groupsym <- sym(groupinput)
      legend_title <- colnames_lookup[[groupinput]]
    } else {
      legend_title <- NULL
    }
    
    # For this one, convert y axis to symbol (only graph type for which Y axis is a dataframe variable)
    yaxissym <- sym(yaxis)
    
    # Filter to only include cell type of interest in expression data
    #expression_data <- expression_data %>%
    #  filter(CellType == celltype)
    
    # Add molecular data
    data <- data %>%
      left_join(molecular_data, by = "LABEL")
    
    # If y axis is AGE_SX or AGE_DX, filter on groups
    if(yaxis %in% c("AGE_SX", "AGE_DX")) {
      data <- data %>% filter(GROUP %in% c("iPD", "LPD"))
    }
    
    # Calculate spearmann rank correlation coefficient and P value for each group
    if (!is.null(groupinput) && groupinput != "None") {
      results_df <- data %>%
        group_by(group = !!groupsym) %>%
        summarise(
          test = list(cor.test(!!xaxissym, !!yaxissym, method = "spearman")),
          .groups = "drop"
        ) %>%
        mutate(tidy_result = lapply(test, broom::tidy)) %>%
        tidyr::unnest(tidy_result) %>%
        select(Group = group, `Spearmann Coefficient` = estimate, pValue = p.value)
      
    } else {
      
      test <- cor.test(data[[as_string(xaxissym)]], data[[as_string(yaxissym)]], method = "spearman")
      
      results_df <- broom::tidy(test) %>%
        mutate(Group = "All") %>%
        select(Group, `Spearmann Coefficient` = estimate, pValue = p.value)
    }
    results_df <- results_df %>%
      mutate(across(where(is.numeric), ~formatC(.x, format = "e", digits = 2))) %>%
      mutate(pValue_num = as.numeric(pValue),
             Significant = pValue_num < 0.05) %>%
      select(-pValue_num)
    
  }
  
return(results_df)
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
              selectInput("yaxis_box", "Y axis", choices = yaxis_options_boxplot)
            ),
            tabPanel("Bar Chart",
                     selectInput("dataset", "Select Patient Outcome", datasets),
                     selectInput("xaxis_bar", "X axis", choices = xaxis_options_barchart),
                     selectInput("groupings_bar", "Group by", choices = group_options_barchart),
                     selectInput("yaxis_bar", "Y axis", choices = c("Number", "Proportion", "Normalized Proportion"))
            ),
            tabPanel("Correlation",
                     selectInput("dataset", "Select Patient Outcome", datasets),
                     selectInput("yaxis_corr", "Y axis", choices = yaxis_options_corr),
                    # selectInput("celltype", "Cell Type", celltypes),
                     selectInput("xaxis_corr", "Cell Type Percent", molecular_names),
                     selectInput("groupings_corr", "Group by", choices = group_options_linegraph)
            )
          ),
        ),
        # Show plot and data table
        mainPanel(
          fluidRow(
            div(
              style = "padding: 10px; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 10px;",
              HTML("
            <ul>
              <br>Explore data for a set of Parkinson's patients (and controls) enrolled in the Parkinson's Progression Markers Initiative (PPMI).
              <br> <b>Choose from the following data representations:</b>
              <li> Line Graph: shows longitudinal data of patient outcome metrics over time. </li>
              <li> Box Plot: shows patient outcome metrics per group. </li>
              <li> Bar Chart: shows number or proportion of individuals per group. </li>
              <li> Correlation: shows correlation between patient outcome and gene module score. </li>
              <br> <b>Hover over data to view details!</b>
            </ul>
          ")
            )
          ),
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
  #celltype <- reactive(input$celltype)
  
  # Define xaxis, yaxis and grouping input variables as reactive values based on plot type
  xaxis <- reactive({
    if (input$plottype == "Line Graph") input$xaxis_line
    else if (input$plottype == "Box Plot") input$xaxis_box
    else if(input$plottype == "Bar Chart") input$xaxis_bar
    else if(input$plottype == "Correlation") input$xaxis_corr
  })
  yaxis <- reactive({
    if(input$plottype == "Bar Chart") input$yaxis_bar
    else if(input$plottype == "Box Plot") input$yaxis_box
    else if(input$plottype == "Correlation") input$yaxis_corr
  })
  groupings <- reactive({
    if (input$plottype == "Line Graph") input$groupings_line
    else if (input$plottype == "Box Plot") input$groupings_box
    else if (input$plottype == "Bar Chart") input$groupings_bar
    else if (input$plottype == "Correlation") input$groupings_corr
  })
  
  # Define data based on user input
  data <- reactive({
    req(input$dataset)
    if(input$plottype %in% c("Correlation", "Box Plot")) {
      df <- corr_data_dfs[[dataset()]] } 
    else {df <- data_dfs[[dataset()]]}
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
        molecular_data = molecular_data,
        type = plottype(),
        plot_title = dataset(),
        xaxis = xaxis(),
        groupinput = groupings(),
        yaxis = yaxis(),
        celltype = celltype())
    }
    return(p)
  })
  
  # Generate results dataframe
  results_dataframe <- reactive({
    df <- NULL
    if(!is.null(data())) {
      df <- make_dataframe(
        data(),
        molecular_data = molecular_data,
        type = plottype(),
        xaxis = xaxis(),
        groupinput = groupings(),
        yaxis = yaxis(),
        celltype = celltype()) 
    }
    return(df)
  })
  
  
  # Render the plot
  output$plot <- renderPlotly({
      make_plot()
    })
  
  # Render data table 
  output$dataTable <- DT::renderDT(results_dataframe(), options = list(pageLength = 5))
}  

##################  Run the application  ##################
shinyApp(ui = ui, server = server)
