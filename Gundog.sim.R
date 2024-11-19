# This script simulates the movement of agents through a grid system based on dead-reckoning (DR) calculations, 
# using empirical distributions for heading, step length, and turn angle. 
# The model is designed to approximate realistic animal movement by incorporating journey phases, 
# grid-cell based distributions, and autocorrelation adjustments. Users can specify parameters for various distributions,
# journey phase handling, autocorrelation, and other options to customize movement behaviors.


# Input data
#ID                <- df$ID                    # Unique identifier for each animal contributing empirical data.
#datetime          <- as.POSIXct(df$DateTime, format="%Y-%m-%d %H:%M:%S")  # Timestamp of each recorded observation.
#heading           <- df$heading.corr          # Heading (in degrees 0-360) for each observation.
#turn.angle        <- df$abs.step.ang          # Absolute turn angle (in degrees) for each observation.
#step.length       <- df$step.dist             # Distance traveled between consecutive observations.
#journey_phase     <- df$trip_status           # Optional - Journey phase label for each observation (Required to be labelled as "outbound" and "inbound").

# Initial conditions
#lo                <- -63.8653                 # Initial longitude coordinate of the agent's starting position.
#la                <- -42.08                   # Initial latitude coordinate of the agent's starting position.
#grid.x            <- df$gridded.X.d           # X-coordinate of the grid cell for each observation. Must be in longlat decimal format.
#grid.y            <- df$gridded.Y.d           # Y-coordinate of the grid cell for each observation. Must be in longlat decimal format.

# Distribution calculation parameters
#min.sector.data   <- 75                       # Minimum data points required per grid sector to compute distributions; if insufficient, fallback 'neighbors' or global options occurs.
#bandwidth         <- "ucv"                    # Smoothing bandwidth for frequency distribution (can also specify "nrd0", "SJ").

# Quantile options
#quantile.turn.angle <- 0.999                  # Default quantile threshold for filtering turn angle data.
#quantile.step.length <- 0.999                 # Default quantile threshold for filtering step length data.
#quantile.scope      <- "per_journey_phase"    # Scope for quantile calculations: "global", "per_journey_phase", or "per_grid" (the latter calculates and applies quantile filtering per journey phase and grid cell).

# Neighborhood search parameters
#radius            <- 1                        # Search radius (in grid cells) for neighboring data in fallback scenarios (1 = directly adjacent cells, incl. diagonals cells). 2 = 2 cells away and so on.
#step_size         <- 1                        # Step size for x-axis in frequency distribution (adjusts granularity).

# Autocorrelation settings
#autocorr.turn.angle <- FALSE                  # Apply 1st order autocorrelation to turn angle if set to TRUE.
#autocorr.step.length <- TRUE                  # Apply 1st order autocorrelation to step length if set to TRUE.
#autocorr.heading    <- FALSE                  # Apply  1st order autocorrelation to heading if set to TRUE.
#auto.corr.heading.method <- "trig"            # Method for calculating autocorrelation for heading: "trig" (trigonometric), "circ.corr" (cor.circular).
#autocorr.scope    <- "per_grid"               # Scope for autocorrelation calculation: "global", "per_journey_phase", or "per_grid" (the latter calculates and applies quantile filtering per journey phase and grid cell). Note, this is calcuated in arranged time sequence per bird and a grand mean computed per journey phase/sector/just a global values, depending on scope used. 

# Randomization parameters
#seed              <- 10                       # Seed for reproducibility in random sampling.

# Movement parameters
#agents            <- 5                        # Number of agents to simulate.
#max.distance.moved <- 270000                  # Maximum cumulative distance (meters) an agent can travel.
#switch_proportion <- 0.45                     # Proportion of maximum distance at which agents switch journey phase (if journey phase is not supplied, then this is irrelevant).
#switch_based_on <- "cumulative"               # Choose switching method ("cumulative" or "straight_line").
#switch_distance <- 50000                      # Straight-line distance in meters for switching (used if switch_based_on == "straight_line")

# Heading offset (for different journey phases)
#head.offset.out   <- 50                       # Heading offset to apply when in "outbound" journey phase. If journey_phase is not specified, then this will be applied to all values.
#head.offset.in    <- -50                      # Heading offset to apply when in "inbound" journey phase.

# Heading distribution and Journey phase settings
#use_global_phase_only <- TRUE                 # TRUE to use journey phase only when computing heading frequency distributions (no grid-based distributions), FALSE for grid + journey phase.
#num_ID_heading <- NULL                  # Control which IDs are used to compute the heading frequency distributions for the simulation. When NULL (default), all available IDs in the are used to compute the heading frequency distributions.
                                               # When a positive integer is provided (e.g., num_ID_heading = 5): The function randomly selects the specified number of IDs from the available data. The same randomly selected IDs are used across all journey phases.
                                               # When a character vector or factor of IDs is provided (e.g., num_ID_heading = c("P2E", "P5D", "P11D")): The function uses only the specified IDs to compute the heading frequency distributions.
                                               # If there is insufficient data for a journey phase, the function falls back to using all available data for that phase and informs the user.

# User-defined behavior for heading and turn angle adjustments.
#user_choice <- "target_heading"                     # user_choice defines turning movement adjustment strategy:
# - "target_heading": Adjust heading to minimize deviation by selecting turn angle direction that when applied to previous heading, results in smallest deviation to the next sampled 'target' heading.
# - "probability_turn": Similar to target_heading, but instead of always choosing the direction that minimizes the deviation, this introduces a probabilistic element where the likelihood of turning left or right depends on how much each direction would minimize the deviation (probabilities inversely proportional to the deviations).
# - "independent_turn" : Sampled turn angle randomly added to or subtracted from selected target heading, regardless of deviation away from previous or next heading
# - "just_heading" : Turn angle not considered, just the sampled heading.

#plot <- TRUE # If TRUE, summary plots are produced

Gundog.sim <- function(ID = 1,
                       datetime,
                       heading,
                       turn.angle,
                       step.length,
                       journey_phase,
                       lo,
                       la,
                       grid.x,
                       grid.y,
                       min.sector.data = 50,
                       bandwidth = "ucv",
                       quantile.turn.angle = 1,
                       quantile.step.length = 1,
                       quantile.scope = "per_journey_phase", 
                       radius = 1,
                       step_size = 1,
                       autocorr.turn.angle = FALSE,
                       autocorr.step.length = FALSE,
                       autocorr.heading  = FALSE,                  
                       auto.corr.heading.method = "circ.corr",            
                       autocorr.scope = "per_grid", 
                       seed = 1,
                       agents = 1,
                       max.distance.moved = 200000,
                       switch_proportion = 0.5,
                       switch_based_on = "cumulative",
                       switch_distance = 50000, 
                       head.offset.out = 0,
                       head.offset.in = 0,
                       num_ID_heading = NULL,
                       use_global_phase_only = TRUE,
                       user_choice = "target_heading",
                       plot = TRUE){
  
  suppressWarnings({
  # Suppress dplyr summarise info
  options(dplyr.summarise.inform = FALSE)
  
  #### User input checks ####
  # Check if 'datetime' is a proper timestamp object
  if (!inherits(datetime, c("POSIXct", "POSIXlt"))) {
    datetime_converted <- as.POSIXct(datetime, tz = "UTC")
    if (any(is.na(datetime_converted))) {
      stop("'datetime' must be a POSIXct/POSIXlt datetime object or convertible to one.", call. = FALSE)
    } else {
      datetime <- datetime_converted
    }
  }
  # List of variables to check
  variables_to_check <- list(
    heading = heading,
    turn.angle = turn.angle,
    step.length = step.length,
    grid.x = grid.x,
    grid.y = grid.y
  )
  # Check if all variables are numeric and of the same length
  variable_lengths <- sapply(variables_to_check, length)
  if (length(unique(variable_lengths)) > 1) {
    stop("All input vectors ('heading', 'turn.angle', 'step.length', 'grid.x', 'grid.y') must be of the same length.", call. = FALSE)
  }
  # Check if all variables are numeric
  non_numeric_vars <- names(variables_to_check)[!sapply(variables_to_check, is.numeric)]
  if (length(non_numeric_vars) > 0) {
    stop(paste("The following variables must be numeric:", paste(non_numeric_vars, collapse = ", ")), call. = FALSE)
  }
  # Check 'grid.x' (longitude)
  if (any(grid.x < -180 | grid.x > 180, na.rm = TRUE)) {
    stop("'grid.x' (longitude) values must be between -180 and 180 degrees.", call. = FALSE)
  }
  
  # Check 'grid.y' (latitude)
  if (any(grid.y < -90 | grid.y > 90, na.rm = TRUE)) {
    stop("'grid.y' (latitude) values must be between -90 and 90 degrees.", call. = FALSE)
  }
  
   ### Automatically Load or Install Required Packages ###
    required_packages <- c("dplyr", "scales", "ggplot2", "gridExtra", "circular", "devtools", "sf", "purrr", "tidyr", "cowplot", "zoo")
    
    install_and_load <- function(packages) {
      for (pkg in packages) {
        if (!require(pkg, character.only = TRUE)) {
          message(paste("Installing missing package:", pkg))
          install.packages(pkg, dependencies = TRUE, type = "source")
          suppressMessages(library(pkg, character.only = TRUE))
        } else {
          suppressMessages(library(pkg, character.only = TRUE))
        }
      }
    }
    missing_pkgs <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
    if (length(missing_pkgs) > 0) {
      message(paste("Installing missing packages:", paste(missing_pkgs, collapse = ", ")))
    }
    # Execute the function to load all required packages
    install_and_load(required_packages)
    
    # Required functions
    # 1)
    #### Create frequency ECDF for a given input variable
    freq.distr <- function(x, bandwidth = "ucv", step_size = 1) {
      suppressWarnings({
        #Remove NA values
        x = na.omit(x)
        if (length(x) == 0) return(NA)  # Return NA if there's no data
        # Determine bandwidth based on user input
        bw <- switch(bandwidth,
                     "ucv" = bw.ucv(x),
                     "nrd0" = bw.nrd0(x),
                     "nrd" = bw.nrd(x),
                     "SJ" = bw.SJ(x),
                     stop("Invalid bandwidth option. Choose 'ucv', 'nrd0', 'nrd', or 'SJ'.")
        )
        # Define the range based on observed min and max values in x, with specified step size
        x_range <- seq(min(x), max(x), by = step_size)
        # Compute density and cumulative density function within the specified range
        x.density <- density(x, bw = bw, from = min(x), to = max(x), n = length(x_range))
        # Compute cumulative distribution function (CDF) based on the density values
        x.cdf <- cumsum(x.density$y) / sum(x.density$y)
        # Return data frame with specified x range and corresponding ECDF values
        df.freq <- data.frame(x = x_range, cdf = x.cdf, density = x.density$y)
        return(df.freq)
      })
    }
    
    # 2)
    # Function to determine the correct UTM zone based on input coordinates
    get_utm_epsg <- function(longitude, latitude) {
      zone_number <- floor((longitude + 180) / 6) + 1 # Calculate UTM zone number based on longitude
      epsg_code <- ifelse(latitude >= 0, 32600 + zone_number, 32700 + zone_number) # Determine if the location is in the Northern [latitude >= 0] or Southern [latitude < 0] Hemisphere 
      return(epsg_code)
    }
    
    # 3)
    # Define a helper function to find neighboring cells with data specific to a given journey_phase
    find_neighbors_by_phase <- function(data, target_x, target_y, journey_phase_check, radius = 1) {
      data %>%
        filter(
          journey_phase == journey_phase_check,
          abs(grid.x.m - target_x) <= radius * cell_size_x, 
          abs(grid.y.m - target_y) <= radius * cell_size_y
        )
    }
    
    # 4)
    # Define 1st order autocorrelation function (for lag-1) 
    compute_acf <- function(x) {
      if(length(na.omit(x)) < 2) return(NA)  # Return NA if not enough data
      return(acf(na.omit(x), plot = FALSE)$acf[2])  # Extract lag-1 autocorrelation
    }
    
    # 5) 
    # Helper function for quantile filtering
    apply_quantile_filter <- function(data, column, quantile_val, scope) {
      if (scope == "global") {
        # Global quantile filtering
        data %>%
          filter(!is.na(.data[[column]])) %>%
          filter(.data[[column]] <= quantile(.data[[column]], quantile_val, na.rm = TRUE))
      } else if (scope == "per_journey_phase") {
        # per_journey_phase quantile filtering
        data %>%
          filter(!is.na(.data[[column]])) %>%
          group_by(journey_phase) %>%
          mutate(threshold = quantile(.data[[column]], quantile_val, na.rm = TRUE)) %>%
          ungroup() %>%
          filter(.data[[column]] <= threshold) %>%
          dplyr::select(-threshold)  # Remove temporary threshold column
      } else if (scope == "per_grid") {
        # Per-grid quantile filtering
        data %>%
          filter(!is.na(.data[[column]])) %>%
          group_by(grid.x.m, grid.y.m, journey_phase) %>%
          mutate(threshold = quantile(.data[[column]], quantile_val, na.rm = TRUE)) %>%
          ungroup() %>%
          filter(.data[[column]] <= threshold) %>%
          dplyr::select(-threshold)  # Remove temporary threshold column
      } else {
        stop("Invalid 'quantile.scope' value. Use 'global', 'per_journey_phase', or 'per_grid'.")
      }
    }
    
    # 6) 
    # Helper function for autocorrelation calculation with grand mean
    calculate_autocorrelation <- function(data, column, scope, min.sector.data = min.sector.data) {
      # Arrange data by ID and datetime for correct temporal sequence
      data <- data %>%
        filter(!is.na(.data[[column]])) %>% # Remove NAs
        arrange(ID, datetime)
      # Define the scope of autocorrelation calculation
      if (scope == "global") {
        # Global autocorrelation by ID, followed by a grand mean across all individuals
        data %>%
          group_by(ID) %>%
          filter(n() >= min.sector.data) %>%  # Ensure minimum data points per ID
          summarise(
            acf_value = compute_acf(.data[[column]]),
            .groups = 'drop'
          ) %>%
          filter(!is.na(acf_value)) %>%
          summarise(grand_mean_acf = mean(acf_value, na.rm = TRUE))
      } else if (scope == "per_journey_phase") {
        # Autocorrelation by journey phase and ID, followed by a grand mean for each journey phase
        data %>%
          group_by(journey_phase, ID) %>%
          filter(n() >= min.sector.data) %>%  # Ensure minimum data points per journey phase and ID
          summarise(
            acf_value = compute_acf(.data[[column]]),
            .groups = 'drop'
          ) %>%
          filter(!is.na(acf_value)) %>%
          group_by(journey_phase) %>%
          summarise(grand_mean_acf = mean(acf_value, na.rm = TRUE))
      } else if (scope == "per_grid") {
        # Autocorrelation by grid, journey phase, and ID, followed by a grand mean for each grid and journey phase
        data %>%
          group_by(grid.x.m, grid.y.m, journey_phase, ID) %>%
          filter(n() >= min.sector.data) %>%  # Ensure minimum data points per grid, journey phase, and ID
          summarise(
            acf_value = compute_acf(.data[[column]]),
            .groups = 'drop'
          ) %>%
          filter(!is.na(acf_value)) %>%
          group_by(grid.x.m, grid.y.m, journey_phase) %>%
          summarise(grand_mean_acf = mean(acf_value, na.rm = TRUE))
        
      } else {
        stop("Invalid 'autocorr.scope' value. Use 'global', 'per_journey_phase', or 'per_grid'.")
      }
    }
    
    #7)
    # Helper function for circular autocorrelation calculation with grand mean
    calculate_circular_autocorrelation <- function(data, column, scope, method = "trig", min.sector.data = min.sector.data, use_global_phase_only = FALSE, selected_ids = NULL) {
      # Define the trigonometric circular lagged autocorrelation function
      circular_lag_acf_trig <- function(angles, lag = 1) {
        angles <- circular(angles, units = "degrees", template = "geographics")
        n <- length(angles)
        if (n <= lag) return(NA)  # Return NA if insufficient data for lag calculation
        # Compute lagged product of cosines and sines
        cos_lag <- cos(angles[-length(angles)]) * cos(angles[-1])
        sin_lag <- sin(angles[-length(angles)]) * sin(angles[-1])
        rho <- mean(cos_lag + sin_lag, na.rm = TRUE)
        return(rho)
      }
      # Define the circular correlation approach for lagged autocorrelation
      circular_lag_acf_circ_corr <- function(angles, lag = 1) {
        angles <- circular(angles, units = "degrees", template = "geographics")
        n <- length(angles)
        if (n <= lag) return(NA)  # Return NA if insufficient data for lag calculation
        # Use `cor.circular` with lagged values
        cor_value <- cor.circular(angles[-length(angles)], angles[-1])
        return(cor_value)
      }
      # Select the desired method
      compute_acf <- if (method == "trig") circular_lag_acf_trig else circular_lag_acf_circ_corr
      # Arrange data by ID and datetime for correct temporal sequence
      data <- data %>%
        filter(!is.na(.data[[column]])) %>% # Remove NAs
        arrange(ID, datetime)
      
      # Apply filters based on 'use_global_phase_only' and 'selected_ids'
      if (!is.null(selected_ids)) {
        data <- data %>% filter(as.character(ID) %in% selected_ids)
      }
      
      if (use_global_phase_only) {
        # Do not use grid-based grouping
        scope <- ifelse(scope == "per_grid", "per_journey_phase", scope)
      }
      
      # Define the scope of circular autocorrelation calculation
      if (scope == "global") {
        # Global circular autocorrelation by ID, followed by a grand mean across all individuals
        data %>%
          group_by(ID) %>%
          filter(n() >= min.sector.data) %>%  # Ensure minimum data points per ID
          summarise(
            circular_acf_value = compute_acf(.data[[column]], lag = 1),
            .groups = 'drop'
          ) %>%
          filter(!is.na(circular_acf_value)) %>%
          summarise(grand_mean_circular_acf = mean.circular(circular(circular_acf_value, units = "degrees")))
      } else if (scope == "per_journey_phase") {
        # Circular autocorrelation by journey phase and ID, followed by a grand mean for each journey phase
        data %>%
          group_by(journey_phase, ID) %>%
          filter(n() >= min.sector.data) %>%  # Ensure minimum data points per journey phase and ID
          summarise(
            circular_acf_value = compute_acf(.data[[column]], lag = 1),
            .groups = 'drop'
          ) %>%
          filter(!is.na(circular_acf_value)) %>%
          group_by(journey_phase) %>%
          summarise(grand_mean_circular_acf = mean.circular(circular(circular_acf_value, units = "degrees")))
      } else if (scope == "per_grid") {
        if (use_global_phase_only) {
          # Change scope to per_journey_phase
          data %>%
            group_by(journey_phase, ID) %>%
            filter(n() >= min.sector.data) %>%
            summarise(
              circular_acf_value = compute_acf(.data[[column]], lag = 1),
              .groups = 'drop'
            ) %>%
            filter(!is.na(circular_acf_value)) %>%
            group_by(journey_phase) %>%
            summarise(grand_mean_circular_acf = mean.circular(circular(circular_acf_value, units = "degrees")))
        } else {
          # Circular autocorrelation by grid, journey phase, and ID, followed by a grand mean for each grid and journey phase
          data %>%
            group_by(grid.x.m, grid.y.m, journey_phase, ID) %>%
            filter(n() >= min.sector.data) %>%  # Ensure minimum data points per grid, journey phase, and ID
            summarise(
              circular_acf_value = compute_acf(.data[[column]], lag = 1),
              .groups = 'drop'
            ) %>%
            filter(!is.na(circular_acf_value)) %>%
            group_by(grid.x.m, grid.y.m, journey_phase) %>%
            summarise(grand_mean_circular_acf = mean.circular(circular(circular_acf_value, units = "degrees")))
        }
      } else {
        stop("Invalid 'scope' value. Use 'global', 'per_journey_phase', or 'per_grid'.")
      }
    }
    
    # 8: Fill missing ACF values using neighbors within specified radius or fallback to global values
    fill_missing_acf <- function(expanded_data, column, radius, global.acf, global_column, data_original, variable_type = "linear", selected_ids = NULL) {
      expanded_data %>%
        rowwise() %>%
        mutate(
          journey_phase_check = journey_phase,
          current_grid_x = grid.x.m,
          current_grid_y = grid.y.m,
          !!sym(column) := if (is.na(!!sym(column))) {
            # For heading data with use_global_phase_only = TRUE, use global ACF
            if (variable_type == "circular" && use_global_phase_only) {
              global.acf %>%
                filter(journey_phase == journey_phase_check) %>%
                pull(global_column)
            } else {
              # Find neighbors with existing ACF values
              neighbors_acf <- find_neighbors_by_phase(
                expanded_data %>% filter(journey_phase == journey_phase_check),
                current_grid_x, current_grid_y, journey_phase = journey_phase_check, radius = radius
              ) %>%
                filter(!is.na(!!sym(column)))
              
              if (nrow(neighbors_acf) >= 1) {
                mean(neighbors_acf[[column]], na.rm = TRUE)
              } else {
                # Fallback to global mean
                global.acf %>%
                  filter(journey_phase == journey_phase_check) %>%
                  pull(global_column)
              }
            }
          } else {
            !!sym(column)
          }
        ) %>%
        dplyr::select(-journey_phase_check, -current_grid_x, -current_grid_y) %>%
        ungroup()
    }
    
    # 9)
    # Circular mean
    circ.mean <-function(x){
      V_east= mean(sin(x * pi/180))
      V_north = mean(cos(x * pi/180))
      mean_MH = (atan2(V_east, V_north))* 180/pi
      mean_MH = (360 + mean_MH) %% 360
      return(mean_MH)
    }
    
    # 10)
    # Haversine distance (m)
    disty = function(long1, lat1, long2, lat2) { #longitude and latitude supplied in degrees
      long1 = long1 * pi/180 ; long2 = long2 * pi/180 ; lat1 = lat1 * pi/180 ; lat2 = lat2 * pi/180 # convert to radians
      a = sin((lat2 - lat1) / 2) * sin((lat2 - lat1) / 2) + cos(lat1) * cos(lat2) * sin((long2 - long1) / 2) * sin((long2 - long1) / 2)
      c = 2 * atan2(sqrt(a), sqrt(1 - a))
      d1 = 6378137 * c
      return(d1)
    }
    
    # 11)
    # Great circular bearing between 2D positions
    beary = function(long1, lat1, long2, lat2) { #Longitude and Latitude supplied in degrees
      long1 = long1 * pi/180 ; long2 = long2 * pi/180 ; lat1 = lat1 * pi/180 ; lat2 = lat2 * pi/180 #Function converts to radians
      a = sin(long2 - long1) * cos(lat2)
      b = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(long2 - long1)
      c = ((atan2(a, b) / pi) * 180)  #Units returned in degrees (-180 to +180 degree scale)
      return(c)
    }
    
    ########################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################
    
    # Create main data frame with all variables
    data <- data.frame(ID = ID,
                       datetime = datetime,
                       heading = heading, 
                       turn.angle = turn.angle,
                       step.length = step.length,
                       grid.x = grid.x,
                       grid.y = grid.y)
    
    # Check if journey_phase is provided and meets the requirements
    if (exists("journey_phase") && !is.null(journey_phase)) {
      # Convert to lowercase and factor
      journey_phase <- factor(tolower(journey_phase))
      # Ensure journey_phase has exactly two levels: 'outbound' and 'inbound'
      unique_levels <- levels(journey_phase)
      if (length(unique_levels) != 2 || !all(sort(unique_levels) == c("inbound", "outbound"))) {
        stop("journey_phase must contain exactly two unique values: 'outbound' and 'inbound'.")
      }
      # Reorder levels to ensure 'outbound' is first
      journey_phase <- factor(journey_phase, levels = c("outbound", "inbound"))
      # Add journey_phase to the data frame
      data$journey_phase <- journey_phase
      print("Journey phase provided")
    } else {
      # Create a default journey_phase column filled with "outbound" if journey_phase is not provided
      data$journey_phase <- "outbound"
      print("Journey phase not provided")
    }
    
    # Apply separate heading offsets (head.offset.in for inbound and head.offset.out for outbound).
    # If this condition isn't met, it applies only head.offset.out to all headings. 
    # Identify unique levels in the 'journey_phase' column
    unique_levels <- na.omit(unique(data$journey_phase))
    # Check if the levels contain exactly 'inbound' and 'outbound'
    if (length(unique_levels) == 2 && all(sort(unique_levels) == c("outbound", "inbound"))) {
      # Apply different offsets based on journey phase
      data <- data %>%
        mutate(
          heading = case_when(
            journey_phase == "outbound" ~ (heading + head.offset.out) %% 360,
            journey_phase == "inbound"  ~ (heading + head.offset.in) %% 360
            #TRUE ~ heading  # For any other case, retain the original heading
          )
        )
    } else {
      # Apply 'head.offset.out' to all headings if the condition isn't met
      data <- data %>%
        mutate(
          heading = (heading + head.offset.out) %% 360
        )
    }
    # Correct circular nature for headings so that they are in range [0, 360)
    data <- data %>%
      mutate(
        heading = ifelse(heading < 0, heading + 360, heading)  # Ensures headings stay positive within 0-360
      )
    
    # Convert data frame to sf object in WGS84 and transform to UTM
    crs = get_utm_epsg(lo, la) # Automatically obtain local UTM CRS based on a reference coordinate
    grid_sf <- st_as_sf(data, coords = c("grid.x", "grid.y"), crs = 4326)
    grid_utm <- st_transform(grid_sf, crs = crs)
    
    # Extract UTM coordinates and add them back to the main data frame
    data$grid.x.m <- st_coordinates(grid_utm)[,1]
    data$grid.y.m <- st_coordinates(grid_utm)[,2]
    data$grid.x.m = round(data$grid.x.m)
    data$grid.y.m = round(data$grid.y.m)
    # Define cell width and height based on the maximum difference between unique x and y values
    cell_size_x <- max(diff(sort(unique(data$grid.x.m))))
    cell_size_y <- max(diff(sort(unique(data$grid.y.m))))
    
    # Print the estimated cell width and height
    print(paste("Estimated cell width (meters):", cell_size_x))
    print(paste("Estimated cell height (meters):", cell_size_y))
    
    # Expand the grid with combinations of grid coordinates and journey phase
    # Identify if journey_phase has only 'outbound'
    if (all(data$journey_phase == "outbound")) {
      # Only 'outbound' phase is present, so no need to duplicate
      grid_template <- expand.grid(grid.x.m = unique(data$grid.x.m), 
                                   grid.y.m = unique(data$grid.y.m),
                                   journey_phase = "outbound")
    } else {
      # Both 'inbound' and 'outbound' are present, apply the full expansion logic
      grid_template <- expand.grid(grid.x.m = unique(data$grid.x.m), 
                                   grid.y.m = unique(data$grid.y.m),
                                   journey_phase = unique(data$journey_phase)) %>%
        group_by(grid.x.m, grid.y.m) %>%
        summarize(
          has_inbound = any(journey_phase == "inbound", na.rm = TRUE),
          has_outbound = any(journey_phase == "outbound", na.rm = TRUE)
        ) %>%
        mutate(
          journey_phase_expanded = case_when(
            !has_inbound & !has_outbound ~ list(c("inbound", "outbound")),
            has_inbound & !has_outbound ~ list(c("inbound", "outbound")),
            !has_inbound & has_outbound ~ list(c("outbound", "inbound")),
            TRUE ~ list(c("inbound", "outbound"))
          )
        ) %>%
        unnest(journey_phase_expanded) %>%
        rename(journey_phase = journey_phase_expanded) %>%
        ungroup()%>%
        # Remove the temporary columns
        dplyr::select(-has_inbound, -has_outbound)
    }
    
    # Merge the template with data to ensure all combinations are included
    data <- grid_template %>%
      left_join(data, by = c("grid.x.m", "grid.y.m", "journey_phase"))
    
    ###########
    
    # 1) Apply quantile filtering for turn angle
    turn_angle_filtered <- apply_quantile_filter(data, "turn.angle", quantile.turn.angle, quantile.scope)
    
    # 2) Apply quantile filtering for step length
    step_length_filtered <- apply_quantile_filter(data, "step.length", quantile.step.length, quantile.scope)
    
    ###########
    
    # Calculate global frequency distributions by journey_phase as a fallback option

    # Compute global heading distributions
    # Initialize a variable to store the selected IDs
    selected_ids <- NULL
    
    if (!is.null(num_ID_heading)) {
      if (is.numeric(num_ID_heading)) {
        # User has provided a number; we will randomly select that many IDs
        num_ids_to_select <- as.integer(num_ID_heading)
        if (num_ids_to_select <= 0) {
          stop("'num_ID_heading' must be a positive integer.", call. = FALSE)
        }
        available_IDs <- unique(as.character(na.omit(data$ID)))
        num_ids_to_select <- min(num_ids_to_select, length(available_IDs))
        selected_ids <- sample(available_IDs, num_ids_to_select)
        # Print message about selected IDs
        message(paste("Randomly selected IDs for heading distributions:", 
                      paste(selected_ids, collapse = ", ")))
      } else if (is.character(num_ID_heading) || is.factor(num_ID_heading)) {
        # User has provided specific IDs (character or factor)
        selected_ids <- unique(as.character(num_ID_heading))
        available_IDs <- unique(as.character(data$ID))
        missing_ids <- setdiff(selected_ids, available_IDs)
        if (length(missing_ids) > 0) {
          stop(paste("The following IDs are not present in the data:", paste(missing_ids, collapse = ", ")), call. = FALSE)
        }
      } else {
        stop("'num_ID_heading' must be either a positive integer or a vector of IDs (character or factor).", call. = FALSE)
      }
      
      # Now compute the global heading distributions using the selected IDs
      global_heading_distributions <- data %>%
        filter(!is.na(heading)) %>%
        group_by(journey_phase) %>%
        do({
          phase_data <- .
          # Filter phase data to include only the selected IDs
          phase_data_filtered <- phase_data %>% filter(as.character(ID) %in% selected_ids)
          # Check if there is sufficient data
          if (nrow(phase_data_filtered) >= min.sector.data) {
            heading_dist <- freq.distr(phase_data_filtered$heading, bandwidth = bandwidth, step_size = step_size)
            data_used <- "Filtered IDs"
          } else {
            # Fall back to using all data for this journey phase
            heading_dist <- freq.distr(phase_data$heading, bandwidth = bandwidth, step_size = step_size)
            data_used <- "All IDs (fallback)"
            message(paste("Insufficient data for journey phase", unique(phase_data$journey_phase), "using filtered IDs. Falling back to all available data for this phase."))
          }
          data.frame(journey_phase = unique(phase_data$journey_phase), heading_dist = I(list(heading_dist)), data_used = data_used)
        }) %>%
        ungroup()
    } else {
      # Original code
      global_heading_distributions <- data %>%
        filter(!is.na(heading)) %>%
        group_by(journey_phase) %>%
        summarise(
          heading_dist = list(freq.distr(heading, bandwidth = bandwidth, step_size = step_size)),
          .groups = "drop"
        ) %>%
        ungroup()
    }
    
    # Global Turn angle frequency
    global_turn_angle_distributions <- turn_angle_filtered %>%
      group_by(journey_phase) %>%
      summarise(
        turn_angle_dist = list(freq.distr(turn.angle, bandwidth = bandwidth, step_size = step_size)), .groups = "drop"
      )
    # Global Step distance frequency
    global_step_length_distributions <- step_length_filtered %>%
      group_by(journey_phase) %>%
      summarise(
        step_length_dist = list(freq.distr(step.length, bandwidth = bandwidth, step_size = step_size)), .groups = "drop"
      )
    
    #### Plot global distributions (according to journey phase) ###
    if(plot == TRUE){
    plot_list <- list()
    counter <- 1
    variables <- c('heading', 'turn.angle', 'step.length')
    
    for (var in variables) {
      # Select the appropriate dataframe
      if (var == 'heading') {
        data_var <- data
         bins <- seq(0, 360, by = step_size)
         # If num_ID_heading is specified, filter data_var to include only selected IDs
         if (!is.null(num_ID_heading)) {
           data_var <- data_var %>% filter(as.character(ID) %in% selected_ids)
         }
       } else if (var == 'turn.angle') {
         data_var <- turn_angle_filtered
         bins <- seq(0, 180, by = step_size)
      } else if (var == 'step.length') {
         data_var <- step_length_filtered
        bins <- seq(0, max(data_var$step.length, na.rm = TRUE) + step_size, by = step_size)
       }
      
      journey_phases <- unique(data_var$journey_phase)
      for (journey_phase_i in journey_phases) {
        x <- subset(data_var, journey_phase == journey_phase_i)
        x <- subset(x, !is.na(x[[var]]))
    
        # Ensure there is data to process
        if (nrow(x) == 0) {
          message(paste("No data available for", var, "in journey phase", journey_phase_i, ". Skipping plot."))
          next
        }
        
        # Add 'x_var' as a column to 'x' data frame
        x <- x %>% mutate(x_var = x[[var]])
        if(step_size > 1){
          x$x_var = as.numeric(as.character(cut(x$x_var, breaks = bins, include.lowest = TRUE, labels = bins[-length(bins)])))
        } else {
          x <- x %>% mutate(x_var = round(x_var))
        }
      
        # Now 'x_var' is available as a column in 'x'
        x_var <- x$x_var
        # Ensure there is data to process
        if (length(x_var) == 0) next
        
        # Use freq.distr function to compute density and cdf
        freq_data <- freq.distr(x_var, bandwidth = bandwidth, step_size = step_size)

        # Interpolate density and cdf values at x_var points using approx
        density_values <- approx(freq_data$x, freq_data$density, xout = x_var, rule = 2)$y
        cdf_values <- approx(freq_data$x, freq_data$cdf, xout = x_var, rule = 2)$y
        
        # Construct x_plot_data
        x_plot_data <- data.frame(x_var = x_var, density = density_values, cdf = cdf_values)
        
        # Compute frequency counts using 'x_var' as a column in 'x'
        freq_counts <- x %>%
          group_by(x_var) %>%
          summarise(freq_count = n())
        x_plot_data <- merge(x_plot_data, freq_counts, by = "x_var", all.x = TRUE)
        
        # Normalize density and freq_count to fit on the cumulative distribution's scale
        max_cdf <- max(x_plot_data$cdf, na.rm = TRUE)
        x_plot_data$density <- x_plot_data$density / max(x_plot_data$density, na.rm = TRUE) * max_cdf
        x_plot_data$freq_count <- x_plot_data$freq_count / max(x_plot_data$freq_count, na.rm = TRUE) * max_cdf
        
        # Remove duplicates
        x_plot_data <- x_plot_data[!duplicated(x_plot_data$x_var), ]
        
        # Plot
        fig <- ggplot(data = x_plot_data, aes(x = x_var, weight = freq_count)) +
          geom_histogram(stat = "count", alpha = 0.4, aes(fill = after_stat(count)), colour = "black") +
          geom_line(aes(y = density), colour = "purple", size = 1) +
          geom_line(aes(y = cdf), colour = "darkgreen", size = 1.4, linetype = "dashed") +
          scale_y_continuous(sec.axis = sec_axis(~ ., name = "ECDF")) +
          xlab(paste0(var, ' (°)')) + ylab('Density') +
          theme_minimal() +
          theme(
            axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "purple", size = 12),
            axis.title.x = element_text(color = "black", size = 14),
            axis.title.y = element_text(color = "purple", size = 14),
            axis.title.y.right = element_text(color = "darkgreen"),
            axis.text.y.right = element_text(color = "darkgreen", size = 14),
            legend.title = element_blank()
          ) +
          scale_fill_gradient('Density', low = 'blue', high = 'red') +
          guides(fill = "none") +
          ggtitle(paste0(journey_phase_i, ' - ', var))
        
        plot_list[[counter]] <- fig
        counter <- counter + 1
      }
    }
    
    # Arrange and display all plots
    print(plot_grid(plotlist = plot_list, ncol = length(unique(data$journey_phase)), align = 'v'))
    }
    
    # Generate frequency distributions for each grid cell, considering neighboring cells if empty
    # This function calculates distributions for each variable based on available data, with fallbacks as needed
    # Generate frequency distributions for each grid cell, with heading flexibility
    extended_distributions <- grid_template %>%
      rowwise() %>%
      mutate(
        journey_phase_check = journey_phase,
        grid.x.m_check = grid.x.m,
        grid.y.m_check = grid.y.m,
        # For 'heading'
        heading_dist = list({
          if (use_global_phase_only) {
            # Use global heading distribution computed earlier
            global_dist <- global_heading_distributions$heading_dist[global_heading_distributions$journey_phase == journey_phase_check][[1]]
            global_dist
          } else {
            # Use grid-based heading distribution
            grid_data <- data %>%
              filter(grid.x.m == grid.x.m_check & grid.y.m == grid.y.m_check & journey_phase == journey_phase_check)
            
            # If num_ID_heading is specified, filter by selected IDs
            if (!is.null(num_ID_heading)) {
              # Use the same selected_ids as before
              ids <- selected_ids
              # Filter grid_data by ids
              grid_data <- grid_data %>% filter(as.character(ID) %in% ids)
            }
            
            # Check if we have enough data in the grid cell
            phase_headings <- grid_data$heading
            if (length(na.omit(phase_headings)) < min.sector.data) {
              # Not enough data, check neighbors
              neighbors <- find_neighbors_by_phase(data %>% filter(journey_phase == journey_phase_check), grid.x.m_check, grid.y.m_check, journey_phase_check, radius = radius)
              
              if (!is.null(num_ID_heading)) {
                # Filter neighbors by ids
                neighbors <- neighbors %>% filter(as.character(ID) %in% ids)
              }
              
              phase_headings <- neighbors$heading
              # Check if we have enough data after including neighbors
              if (length(na.omit(phase_headings)) < min.sector.data) {
                # If still insufficient data, use global distribution
                global_dist <- global_heading_distributions$heading_dist[global_heading_distributions$journey_phase == journey_phase_check][[1]]
                global_dist
              } else {
                # Compute the distribution based on the available data from neighbors
                freq.distr(phase_headings, bandwidth = bandwidth, step_size = step_size)
              }
            } else {
              # Compute the distribution based on the available data in the grid cell
              freq.distr(phase_headings, bandwidth = bandwidth, step_size = step_size)
            }
          }
        }),
        
        # For 'turn.angle' (no change)
        turn_angle_dist = list({
          grid_data <- turn_angle_filtered$turn.angle[turn_angle_filtered$grid.x.m == grid.x.m_check & turn_angle_filtered$grid.y.m == grid.y.m_check & turn_angle_filtered$journey_phase == journey_phase_check]
          if (length(na.omit(grid_data)) < min.sector.data) {
            neighbors = find_neighbors_by_phase(turn_angle_filtered %>% filter(journey_phase == journey_phase_check), grid.x.m_check, grid.y.m_check, journey_phase_check, radius = radius)
            if (nrow(neighbors) > min.sector.data) {
              freq.distr(neighbors$turn.angle, bandwidth = bandwidth, step_size = step_size)
            } else {
              global_turn_angle_distributions$turn_angle_dist[global_turn_angle_distributions$journey_phase == journey_phase_check][[1]]
            }
          } else {
            freq.distr(grid_data, bandwidth = bandwidth, step_size = step_size)
          }
        }),
        
        # For 'step.length' (no change)
        step_length_dist = list({
          grid_data <- step_length_filtered$step.length[step_length_filtered$grid.x.m == grid.x.m_check & step_length_filtered$grid.y.m == grid.y.m_check & step_length_filtered$journey_phase == journey_phase_check]
          if (length(na.omit(grid_data)) < min.sector.data) {
            neighbors = find_neighbors_by_phase(step_length_filtered %>% filter(journey_phase == journey_phase_check), grid.x.m_check, grid.y.m_check, journey_phase_check, radius = radius)
            if (nrow(neighbors) > min.sector.data) {
              freq.distr(neighbors$step.length, bandwidth = bandwidth, step_size = step_size)
            } else {
              global_step_length_distributions$step_length_dist[global_step_length_distributions$journey_phase == journey_phase_check][[1]]
            }
          } else {
            freq.distr(grid_data, bandwidth = bandwidth, step_size = step_size)
          }
        })
      ) %>%
      ungroup()
    
    # Generate metadata table to track data presence and source type
    # This table will provide insights into which data source was used for each grid cell distribution
    # Generate metadata table to track data presence and source type for heading and step length/turn angle
    metadata_table <- grid_template %>%
      rowwise() %>%
      mutate(
        journey_phase_check = journey_phase,
        current_grid_x = grid.x.m,
        current_grid_y = grid.y.m,
        
        ### Heading Data Status and Source Type ###
        heading_data_status = {
          if (use_global_phase_only) {
            "Global Data"
          } else {
            # Attempt to use grid data
            grid_data <- data %>%
              filter(grid.x.m == current_grid_x & grid.y.m == current_grid_y & journey_phase == journey_phase_check)
            
            if (!is.null(num_ID_heading)) {
              # Use selected IDs
              grid_data <- grid_data %>% filter(as.character(ID) %in% selected_ids)
            }
            
            if (length(na.omit(grid_data$heading)) >= min.sector.data) {
              "Grid Data"
            } else {
              # Check neighbors
              neighbors <- find_neighbors_by_phase(
                data %>% filter(journey_phase == journey_phase_check),
                current_grid_x, current_grid_y, journey_phase_check, radius = radius
              )
              
              if (!is.null(num_ID_heading)) {
                # Filter neighbors by selected IDs
                neighbors <- neighbors %>% filter(as.character(ID) %in% selected_ids)
              }
              
              if (length(na.omit(neighbors$heading)) >= min.sector.data) {
                "Neighbor Data"
              } else {
                "Global Data"
              }
            }
          }
        },
        
        ### Step Length and Turn Angle Data Status and Source Type ###
        step_turn_data_status = {
          grid_data <- step_length_filtered %>%
            filter(grid.x.m == current_grid_x & grid.y.m == current_grid_y & journey_phase == journey_phase_check)
          
          if (length(na.omit(grid_data$step.length)) >= min.sector.data) {
            "Grid Data"
          } else {
            # Check neighbors
            neighbors <- find_neighbors_by_phase(
              step_length_filtered %>% filter(journey_phase == journey_phase_check),
              current_grid_x, current_grid_y, journey_phase_check = journey_phase_check, radius = radius
            )
            
            if (nrow(neighbors) >= min.sector.data) {
              "Neighbor Data"
            } else {
              "Global Data"
            }
          }
        }
      ) %>%
      ungroup()
    
    # Bind metadata table to extended_distributions
    # This combines the calculated distributions with metadata information for easy tracking and analysis
    extended_distributions <- cbind(extended_distributions, metadata_table[, c(7:8)])
    
    ######################################################################################################################################################################################################################################
    
    # Summary plot of grid data count and type of data frequency distributions used
    
    # Prepare the count data per grid cell for plotting, including the metadata source type
    # Filter heading data based on selected IDs if specified
    if(plot == TRUE){
    if (!is.null(num_ID_heading)) {
      # Use selected_ids from earlier
      heading_filtered <- data %>%
        filter(!is.na(heading), as.character(ID) %in% selected_ids)
    } else {
      # Use all data
      heading_filtered <- data %>%
        filter(!is.na(heading))
    }
    
    # Prepare the count data per grid cell for heading data
    grid_counts_heading <- heading_filtered %>%
      group_by(grid.x.m, grid.y.m, journey_phase) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        xmin = grid.x.m - cell_size_x / 2,
        xmax = grid.x.m + cell_size_x / 2,
        ymin = grid.y.m - cell_size_y / 2,
        ymax = grid.y.m + cell_size_y / 2
      )
    
    # Prepare the count data per grid cell for step/turn data
    grid_counts_step_turn <- step_length_filtered %>%
      group_by(grid.x.m, grid.y.m, journey_phase) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        xmin = grid.x.m - cell_size_x / 2,
        xmax = grid.x.m + cell_size_x / 2,
        ymin = grid.y.m - cell_size_y / 2,
        ymax = grid.y.m + cell_size_y / 2
      )
    
    ### Step 2: Adjust the Plots
    
    # Define color scale with 10 colors
    colfunc <- c("white", "lightblue", "cyan", "green", "yellowgreen", "yellow", "orange", "darkorange", "red", "darkred")
    
    # Combine counts to find overall max and quantiles
    combined_counts <- bind_rows(
      grid_counts_heading %>% mutate(data_type = "Heading"),
      grid_counts_step_turn %>% mutate(data_type = "Step/Turn")
    )
    
    val_max_combined <- max(combined_counts$count, na.rm = TRUE)
    val_99_combined <- ceiling(as.numeric(quantile(combined_counts$count, na.rm = TRUE, 0.99)))
    val_scaled_combined <- seq(0, val_99_combined, length.out = length(colfunc))
    
    # Plot for Heading Data
    p_heading <- ggplot() +
      theme_bw() +
      theme(
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 18),
        axis.title.y = element_text(color = "black", size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      geom_rect(
        data = grid_counts_heading,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = count),
        color = "grey70"
      ) +
      scale_fill_gradientn(
        colours = colfunc,
        values = rescale(val_scaled_combined),
        limits = c(0, val_99_combined),
        oob = scales::squish,
        name = expression("Data Count [" * P[99] * "]")
      ) +
      geom_point(
        data = extended_distributions,
        aes(x = grid.x.m, y = grid.y.m, shape = heading_data_status),
        size = 2, color = "black"
      ) +
      scale_shape_manual(
        name = "Heading Data Source Type",
        values = c("Grid Data" = 16, "Neighbor Data" = 17, "Global Data" = 4)
      ) +
      coord_equal() +
      geom_hline(yintercept = unique(grid_counts_heading$ymin), color = "grey70", linetype = "dotted") +
      geom_vline(xintercept = unique(grid_counts_heading$xmin), color = "grey70", linetype = "dotted") +
      labs(
        x = "Distance E-W (m)",
        y = "Distance N-S (m)",
        title = "Heading Data"
      ) + facet_grid(. ~ journey_phase)
    
    # Plot for Step/Turn Data
    p_step_turn <- ggplot() +
      theme_bw() +
      theme(
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 18),
        axis.title.y = element_text(color = "black", size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      geom_rect(
        data = grid_counts_step_turn,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = count),
        color = "grey70"
      ) +
      scale_fill_gradientn(
        colours = colfunc,
        values = rescale(val_scaled_combined),
        limits = c(0, val_99_combined),
        oob = scales::squish,
        name = expression("Data Count [" * P[99] * "]")
      ) +
      geom_point(
        data = extended_distributions,
        aes(x = grid.x.m, y = grid.y.m, shape = step_turn_data_status),
        size = 2, color = "black"
      ) +
      scale_shape_manual(
        name = "Step/Turn Data Source Type",
        values = c("Grid Data" = 16, "Neighbor Data" = 17, "Global Data" = 4)
      ) +
      coord_equal() +
      geom_hline(yintercept = unique(grid_counts_step_turn$ymin), color = "grey70", linetype = "dotted") +
      geom_vline(xintercept = unique(grid_counts_step_turn$xmin), color = "grey70", linetype = "dotted") +
      labs(
        x = "Distance E-W (m)",
        y = "Distance N-S (m)",
        title = "Step Length and Turn Angle Data"
      ) + facet_grid(. ~ journey_phase)
    
    ### Step 3: Combine Plots into One Panel
    
    # Suppress legends in both plots
    p_heading <- p_heading + theme(legend.position = "none")
    p_step_turn <- p_step_turn + theme(legend.position = "none")
    
    ### Step 2: Extract Both Legends
    
    # Increase legend sizes and adjust margins before extracting the legends
    p_heading_legend <- p_heading +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5)
      ) +
      guides(fill = FALSE)  # Hide the 'Data Count' legend
     
    p_step_turn_legend <- p_step_turn +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5)
      ) +
      guides(shape = FALSE)  # Hide the 'Source Type' legend
    
    # Extract the legends
    legend_data_count <- get_legend(p_heading_legend)
    legend_source_type <- get_legend(p_step_turn_legend)
    
    ### Step 2: Align the x-axis and y-axis Breaks
    
    # Get the x-axis and y-axis range from both plots
    x_range <- range(
      c(grid_counts_heading$xmin, grid_counts_heading$xmax,
        grid_counts_step_turn$xmin, grid_counts_step_turn$xmax),
      na.rm = TRUE
    )
    y_range <- range(
      c(grid_counts_heading$ymin, grid_counts_heading$ymax,
        grid_counts_step_turn$ymin, grid_counts_step_turn$ymax),
      na.rm = TRUE
    )
    # Define common breaks
    x_breaks <- pretty(x_range, n = 3)
    y_breaks <- pretty(y_range, n = 3)
    
    # Apply the same x-axis breaks to both plots
    p_heading <- p_heading +
      scale_x_continuous(breaks = x_breaks) +
      scale_y_continuous(breaks = y_breaks)
    p_step_turn <- p_step_turn +
      scale_x_continuous(breaks = x_breaks) +
      scale_y_continuous(breaks = y_breaks)
    
    # Add an overall title
    title_text <- if (!is.null(num_ID_heading)) {
      paste0(
        "Data Counts and Source Types \n(counts based on quantile-filtered data and selected IDs for heading)"
      )
    } else {
      "Data Counts and Source Types \n(counts based on quantile-filtered data)"
    }
    
    ### Step 4: Assemble the Final Plot
    
    # Create the title grob
    title <- ggdraw() +
      draw_label(
        title_text,
        fontface = 'bold',
        x = 0.5,
        hjust = 0.5,
        size = 16
      )
    # Combine the plots vertically
    plots_combined <- plot_grid(
      p_heading,
      p_step_turn,
      labels = c("A", "B"),
      ncol = 1,
      align = "v"
    )
    # Combine legends vertically
    # Combine legends vertically
    legend_combined <- plot_grid(
      legend_data_count,
      legend_source_type,
      ncol = 1,
      align = "v"
    )
    
    # Adjust rel_widths to reduce the gap between plots and legends
    plot_with_legend <- plot_grid(
      plots_combined,
      legend_combined,
      ncol = 2,
      rel_widths = c(3, 0.9)  # Adjust the second value as needed
    )
    # Combine the title and plot
    final_plot_with_title <- plot_grid(
      title,
      plot_with_legend,
      ncol = 1,
      rel_heights = c(0.1, 0.9)
    )
    
    # Print the final plot
    print(final_plot_with_title)
    }
    
    ######################################################################################################################################################################################################################################
    
    # Autocorrelation of metrics?
    
    # Step 1: Compute initial ACF values for each grid, as specified by autocorr.scope
    # Only proceed with neighbor filling if autocorr.scope is "per_grid"
    if(autocorr.turn.angle == TRUE || autocorr.step.length == TRUE || autocorr.heading == TRUE){
      if (autocorr.scope == "per_grid") {
        
        # Step 1: Compute per-grid ACF values for turn angle, step length, and heading
        # For turn angle
        if (autocorr.turn.angle == TRUE) {
          autocorr_turn_angle <- calculate_autocorrelation(turn_angle_filtered, "turn.angle", autocorr.scope, min.sector.data)
          autocorr_turn_angle <- autocorr_turn_angle %>%
            rename(acf_turn_angle = grand_mean_acf)
        }
        # For step length
        if (autocorr.step.length == TRUE) {
          autocorr_step_length <- calculate_autocorrelation(step_length_filtered, "step.length", autocorr.scope, min.sector.data)
          autocorr_step_length <- autocorr_step_length %>%
            rename(acf_step_length = grand_mean_acf)
        }
        # For heading
        if (autocorr.heading == TRUE) {
          autocorr_heading <- calculate_circular_autocorrelation(
            data, "heading", autocorr.scope, method = auto.corr.heading.method,
            min.sector.data = min.sector.data, use_global_phase_only = use_global_phase_only, selected_ids = selected_ids
          )
        }

        # Step 2: Merge computed ACF values with expanded grid template to ensure all grids are represented
        expanded_grid_acf <- grid_template
        if (autocorr.turn.angle == TRUE) {
          expanded_grid_acf <- expanded_grid_acf %>%
            left_join(autocorr_turn_angle, by = c("grid.x.m", "grid.y.m", "journey_phase"))
        }
        if (autocorr.step.length == TRUE) {
          expanded_grid_acf <- expanded_grid_acf %>%
            left_join(autocorr_step_length, by = c("grid.x.m", "grid.y.m", "journey_phase"))
        }
        if (autocorr.heading == TRUE) {
          if (use_global_phase_only && autocorr.scope == "per_grid") {
            # Merge only by journey_phase since grid columns are not present
            expanded_grid_acf <- expanded_grid_acf %>%
              left_join(autocorr_heading, by = "journey_phase")
          } else {
            expanded_grid_acf <- expanded_grid_acf %>%
              left_join(autocorr_heading, by = c("grid.x.m", "grid.y.m", "journey_phase"))
          }
        }
        
        # Step 3: Compute global (per_journey_phase) autocorrelation for fallback values
        if (autocorr.turn.angle == TRUE) {
          autocorr_turn_angle.global <- calculate_autocorrelation(turn_angle_filtered, "turn.angle", "per_journey_phase", min.sector.data)
          autocorr_turn_angle.global <- autocorr_turn_angle.global %>%
            rename(acf_turn_angle = grand_mean_acf)
        }
        if (autocorr.step.length == TRUE) {
          autocorr_step_length.global <- calculate_autocorrelation(step_length_filtered, "step.length", "per_journey_phase", min.sector.data)
          autocorr_step_length.global <- autocorr_step_length.global %>%
            rename(acf_step_length = grand_mean_acf)
        }
        if (autocorr.heading == TRUE) {
          autocorr_heading.global <- calculate_circular_autocorrelation(
            data, "heading", "per_journey_phase", method = auto.corr.heading.method,
            min.sector.data = min.sector.data, use_global_phase_only = use_global_phase_only, selected_ids = selected_ids
          )
        }
        
        # Step 4: Fill missing ACF values
        if (autocorr.turn.angle == TRUE) {
          expanded_grid_acf <- fill_missing_acf(
            expanded_grid_acf, "acf_turn_angle", radius,
            global.acf = autocorr_turn_angle.global, global_column = "acf_turn_angle",
            data_original = turn_angle_filtered, variable_type = "linear"
          )
        }
        if (autocorr.step.length == TRUE) {
          expanded_grid_acf <- fill_missing_acf(
            expanded_grid_acf, "acf_step_length", radius,
            global.acf = autocorr_step_length.global, global_column = "acf_step_length",
            data_original = step_length_filtered, variable_type = "linear"
          )
        }
        if (autocorr.heading == TRUE) {
          if (!(use_global_phase_only && autocorr.scope == "per_grid")) {
            expanded_grid_acf <- fill_missing_acf(
              expanded_grid_acf, "grand_mean_circular_acf", radius,
              global.acf = autocorr_heading.global, global_column = "grand_mean_circular_acf",
              data_original = data, variable_type = "circular", selected_ids = selected_ids
            )
          } else {
            # Since we have journey phase level data, set acf_status_heading to "Global Data"
            expanded_grid_acf$acf_status_heading <- "Global Data"
          }
        }
        
        # Create a metadata table for tracking the source of each ACF value
        metadata_acf_table <- grid_template %>%
          rowwise() %>%
          mutate(
            journey_phase_check = journey_phase,
            current_grid_x = grid.x.m,
            current_grid_y = grid.y.m,
            
            # Metadata for turn angle ACF source
            acf_status_turn_angle = if (autocorr.turn.angle == TRUE) {
              # Attempt to use grid data
              grid_data_ta <- turn_angle_filtered %>%
                filter(grid.x.m == current_grid_x & grid.y.m == current_grid_y & journey_phase == journey_phase_check)
              
              if (nrow(grid_data_ta) >= min.sector.data) {
                "Grid Data"
              } else {
                # Check neighbors
                neighbors_ta <- find_neighbors_by_phase(
                  turn_angle_filtered %>% filter(journey_phase == journey_phase_check),
                  current_grid_x, current_grid_y, journey_phase_check, radius = radius
                )
                
                if (nrow(neighbors_ta) >= min.sector.data) {
                  "Neighbor Data"
                } else {
                  "Global Data"
                }
              }
            } else {
              NA_character_
            },
            
            # Metadata for step length ACF source
            acf_status_step_length = if (autocorr.step.length == TRUE) {
              # Attempt to use grid data
              grid_data_sl <- step_length_filtered %>%
                filter(grid.x.m == current_grid_x & grid.y.m == current_grid_y & journey_phase == journey_phase_check)
              
              if (nrow(grid_data_sl) >= min.sector.data) {
                "Grid Data"
              } else {
                # Check neighbors
                neighbors_sl <- find_neighbors_by_phase(
                  step_length_filtered %>% filter(journey_phase == journey_phase_check),
                  current_grid_x, current_grid_y, journey_phase_check, radius = radius
                )
                
                if (nrow(neighbors_sl) >= min.sector.data) {
                  "Neighbor Data"
                } else {
                  "Global Data"
                }
              }
            } else {
              NA_character_
            },
            
            # Metadata for heading ACF source
            acf_status_heading = if (autocorr.heading == TRUE) {
              if (use_global_phase_only && autocorr.scope == "per_grid") {
                "Global Data"
              } else if (use_global_phase_only) {
                "Global Data"
              } else {
                # Attempt to use grid data
                grid_data_h <- data %>%
                  filter(grid.x.m == current_grid_x & grid.y.m == current_grid_y & journey_phase == journey_phase_check)
                
                if (!is.null(selected_ids)) {
                  grid_data_h <- grid_data_h %>% filter(as.character(ID) %in% selected_ids)
                }
                
                if (nrow(grid_data_h) >= min.sector.data) {
                  "Grid Data"
                } else {
                  # Check neighbors
                  neighbors_h <- find_neighbors_by_phase(
                    data %>% filter(journey_phase == journey_phase_check),
                    current_grid_x, current_grid_y, journey_phase_check, radius = radius
                  )
                  
                  if (!is.null(selected_ids)) {
                    neighbors_h <- neighbors_h %>% filter(as.character(ID) %in% selected_ids)
                  }
                  
                  if (nrow(neighbors_h) >= min.sector.data) {
                    "Neighbor Data"
                  } else {
                    "Global Data"
                  }
                }
              }
            } else {
              NA_character_
            }
          ) %>%
          ungroup() %>%
          dplyr::select(-journey_phase_check, -current_grid_x, -current_grid_y)
        
        # Bind the metadata table to the ACF values
        expanded_grid_acf <- expanded_grid_acf %>%
          left_join(metadata_acf_table, by = c("grid.x.m", "grid.y.m", "journey_phase"))
        
      } else if (autocorr.scope != "per_grid") {

        # If not per_grid, simply calculate global or journey-phase level ACFs without expanding
        # For turn angle
        if (autocorr.turn.angle == TRUE) {
          autocorr_turn_angle <- calculate_autocorrelation(turn_angle_filtered, "turn.angle", autocorr.scope, min.sector.data)
        }
        
        # For step length
        if (autocorr.step.length == TRUE) {
          autocorr_step_length <- calculate_autocorrelation(step_length_filtered, "step.length", autocorr.scope, min.sector.data)
        }
        
        # For heading
        if (autocorr.heading == TRUE) {
          # Adjust scope for heading data if 'use_global_phase_only' is TRUE
          heading_autocorr_scope <- if (use_global_phase_only && autocorr.scope == "per_grid") "per_journey_phase" else autocorr.scope
          autocorr_heading <- calculate_circular_autocorrelation(
            data, "heading", heading_autocorr_scope, method = auto.corr.heading.method,
            min.sector.data = min.sector.data, use_global_phase_only = use_global_phase_only, selected_ids = selected_ids
          )
        }
        
        # Since we're not dealing with grids, we can directly use the computed ACF values
        # For plotting or further analysis, we might need to prepare a data frame that includes journey phases
        
        # Prepare a data frame with the ACF values
        if (autocorr.scope == "global") {
          expanded_acf <- data.frame(
            journey_phase = unique(data$journey_phase)
          )
          
          if (autocorr.turn.angle == TRUE) {
            expanded_acf$acf_turn_angle <- autocorr_turn_angle$grand_mean_acf
          }
          
          if (autocorr.step.length == TRUE) {
            expanded_acf$acf_step_length <- autocorr_step_length$grand_mean_acf
          }
          
          if (autocorr.heading == TRUE) {
            expanded_acf$acf_heading <- autocorr_heading$grand_mean_circular_acf
          }
          
        } else if (autocorr.scope == "per_journey_phase") {
          expanded_acf <- data.frame(
            journey_phase = unique(data$journey_phase)
          )
          
          if (autocorr.turn.angle == TRUE) {
            expanded_acf <- expanded_acf %>%
              left_join(autocorr_turn_angle, by = "journey_phase")
          }
          
          if (autocorr.step.length == TRUE) {
            expanded_acf <- expanded_acf %>%
              left_join(autocorr_step_length, by = "journey_phase")
          }
          
          if (autocorr.heading == TRUE) {
            expanded_acf <- expanded_acf %>%
              left_join(autocorr_heading, by = "journey_phase")
          }
        }
        
        # Metadata: Since we are not using grids or neighbors, the source is straightforward
        # For each variable, the source is either "Journey Phase Data" or "Global Data"
        
        expanded_acf <- expanded_acf %>%
          mutate(
            acf_status_turn_angle = if (autocorr.turn.angle == TRUE) {
              if (autocorr.scope == "per_journey_phase") {
                "Journey Phase Data"
              } else if (autocorr.scope == "global") {
                "Global Data"
              } else {
                NA_character_
              }
            } else {
              NA_character_
            },
            acf_status_step_length = if (autocorr.step.length == TRUE) {
              if (autocorr.scope == "per_journey_phase") {
                "Journey Phase Data"
              } else if (autocorr.scope == "global") {
                "Global Data"
              } else {
                NA_character_
              }
            } else {
              NA_character_
            },
            acf_status_heading = if (autocorr.heading == TRUE) {
              if (use_global_phase_only) {
                "Global Data"
              } else if (autocorr.scope == "per_journey_phase") {
                "Journey Phase Data"
              } else if (autocorr.scope == "global") {
                "Global Data"
              } else {
                NA_character_
              }
            } else {
              NA_character_
            }
          )
      }
    }
    
    #####################################################################################################################################################################
    #####################################################################################################################################################################
    ########################### Start agent based model ###############################
    # Assuming DR.lat and DR.lon are storing latitude and longitude
    # Convert expanded grids in metres projection back to longlat format
    grid_sf <- st_as_sf(data, coords = c("grid.x.m", "grid.y.m"), crs = crs)
    grid_wgs84 <- st_transform(grid_sf, crs = 4326)
    # Extract longlat coordinates and add them back to the main data frame
    data$long<- st_coordinates(grid_wgs84)[,1]
    data$lat <- st_coordinates(grid_wgs84)[,2]
    
    # Define the boundary limits in lat/lon for your area
    lat_min <- min(data$lat)
    lat_max <- max(data$lat)
    lon_min <- min(data$long)
    lon_max <- max(data$long)
    lat_min.m <- min(data$grid.y.m)
    lat_max.m <- max(data$grid.y.m)
    lon_min.m <- min(data$grid.x.m)
    lon_max.m <- max(data$grid.x.m)
    
    #Starting grid
    # Find closest grid cell using Haversine distance
    closest_grid <- data %>%
      mutate(
        distance = disty(lo, la, long, lat)
      ) %>%
      slice_min(distance, n = 1)  # Extract the closest grid - First sector 
    
    # Initialize data frame to store all agents' movements
    df_agents <- data.frame(long = numeric(), lat = numeric(), id = numeric(), cum.dist = numeric(), grid_x = numeric(), grid_y = numeric(), grid_x.utm = numeric(), grid_y.utm = numeric())
    switch_proportion <- switch_proportion # Proportion of cumulative distance travelled before journey phase switches (if specified)
    
    # Set initial sector type as "outbound"
    sector_type <- "outbound"
    
    #Set seed?
    if(is.null(seed) != TRUE){
      set.seed(seed) 
    }
    
    #Interactive Plot
    # Define color palette for agents
    agent_colors <- rainbow(agents)  # Generates a different color for each agent
    # Plot the initial position for this agent with assigned color
    
    
    # How many agents?
    DR.lat <- vector(mode = "numeric", length = 0) # Empty vector to be filled with dead-reckoned longitude coordinates
    DR.lon <- vector(mode = "numeric", length = 0) # Empty vector to be filled with dead-reckoned latitude coordinates
    # Set initial position
    DR.lat[1] <- la * pi/180 # Initial latitude in radians
    DR.lon[1] <- lo * pi/180 # Initial longitude in radians
    
    for(a in 1:agents){
     
      #Circular mean of first sector - Set initial heading value to mean circular heading of first sector
      h = circ.mean(na.omit(closest_grid$heading[closest_grid$journey_phase == "outbound"])) 
      # Set the color for this agent
      if(plot == TRUE){
      agent_color <- agent_colors[a]
      if(a == 1){
      plot(DR.lon * 180 / pi, DR.lat * 180 / pi, xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), 
           xlab = "DR Longitude", ylab = "DR Latitude", col = agent_color[a], pch = 16, main =  "Agent Movement Trajectory")
      }
      }
        # Initialize time and distance trackers for plotting
      last_update_time <- Sys.time()   # For time-based updates
      last_dist_update <- 0            # For distance-based updates
      
      # Initialize variables
      cum.dist <- 0 # Initial cumulative distance
      DR.lat <- vector(mode = "numeric", length = 0) # Empty vector to be filled with dead-reckoned longitude coordinates
      DR.lon <- vector(mode = "numeric", length = 0) # Empty vector to be filled with dead-reckoned latitude coordinates
      dist <- vector(mode = "numeric", length = 0) ; dist[1] = 0 # Empty vector to be filled with distance travelled (m)
      straight_line_dist <- vector(mode = "numeric", length = 0); straight_line_dist[1] = 0 # Empty vector to be filled with straight-line distance travelled (m)
      
      # Set initial position
      DR.lat[1] <- la * pi/180 # Initial latitude in radians
      DR.lon[1] <- lo * pi/180 # Initial longitude in radians
      
      # Set initial heading
      simulated.head <- vector(mode = "numeric", length = 0) ; simulated.head[1] = h # Empty vector to be filled with selected heading values from freq. distr. Initial value is h
      simulated.dist <- vector(mode = "numeric", length = 0) ; simulated.dist[1] = 0 # Empty vector to be filled with selected step length (distance) values from freq. distr. Initial value is 0
      simulated.turns <- vector(mode = "numeric", length = 0) ; simulated.turns[1] = 0 # If turn angles supplied; empty vector to be filled with selected turn angle values from freq. distr. Initial value is 0
      
      # Gridded sectors
      grid_x <- vector(mode = "numeric", length = 0)
      grid_y <- vector(mode = "numeric", length = 0)
      grid_x.utm <- vector(mode = "numeric", length = 0)
      grid_y.utm <- vector(mode = "numeric", length = 0)
      
      i <- 2 # Step index
      
      while(cum.dist <= max.distance.moved){
        
        # Convert current latitude and longitude to degrees
        current_lon_deg <- DR.lon[i-1] * 180 / pi
        current_lat_deg <- DR.lat[i-1] * 180 / pi
        
        # Calculate the Haversine distance to all grid points in latitude and longitude
        closest_row <- data[which.min(disty(current_lon_deg, current_lat_deg, data$long, data$lat)),]
        
        # Extract the closest grid coordinates in both long/lat and meters
        grid_x_current <- closest_row$long
        grid_y_current <- closest_row$lat
        grid_x_current_m <- closest_row$grid.x.m
        grid_y_current_m <- closest_row$grid.y.m
        
        df.sub.distr <- extended_distributions %>%
          filter(
            journey_phase == sector_type,
            round(grid.x.m) == round(grid_x_current_m),  # Rounding to avoid floating-point issues
            round(grid.y.m) == round(grid_y_current_m)
          )
        if(autocorr.step.length || autocorr.turn.angle || autocorr.heading){
        df.sub.acf <- expanded_grid_acf %>%
          filter(
            journey_phase == sector_type,
            round(grid.x.m) == round(grid_x_current_m),  # Rounding to avoid floating-point issues
            round(grid.y.m) == round(grid_y_current_m)
          )
        }
        
        # Simulate a heading value from ECDF
        prob_head <- runif(1, 0, 1)
        if(!autocorr.heading){
          simulated.head[i] <- approx(
            x = df.sub.distr$heading_dist[[1]][[2]], 
            y = df.sub.distr$heading_dist[[1]][[1]], 
            xout= prob_head, 
            rule=2)$y
        } else {
          # Incorporate autocorrelation for circular data
          sampled_head <- approx(x = df.sub.distr$heading_dist[[1]][[2]], 
                                 y = df.sub.distr$heading_dist[[1]][[1]], 
                                 xout= prob_head, 
                                 rule=2)$y
          # Ensure autocorrelation coefficient is between 0 and 1
          acf_heading_weight <- pmax(0, pmin(df.sub.acf$acf_heading, 1))
          # Calculate weighted circular mean using the autocorrelation
          simulated.head[i] <- as.numeric(
            circular::mean.circular(
              circular(c(simulated.head[i-1],  sampled_head), type = "angles", units = "degrees", template = "geographics"),  
              w = c(acf_heading_weight, 1 - acf_heading_weight)))
          simulated.head[i] <- (simulated.head[i] + 360) %% 360
        } 
        
        # Simulate a step length value from ECDF
        prob_step <- runif(1, 0, 1)
        if(!autocorr.step.length){
          simulated.dist[i] <- approx(x = df.sub.distr$step_length_dist[[1]][[2]], 
                                      y = df.sub.distr$step_length_dist[[1]][[1]], 
                                      xout= prob_step, 
                                      rule=2)$y
        } else {
          sampled_dist <- approx(x = df.sub.distr$step_length_dist[[1]][[2]],
                                 y = df.sub.distr$step_length_dist[[1]][[1]], 
                                 xout= prob_step, 
                                 rule=2)$y
          # Update step length with autocorrelation
          simulated.dist[i] <- (df.sub.acf$acf_step_length * simulated.dist[i-1] + 
                                  (1 - df.sub.acf$acf_step_length) * sampled_dist)
          simulated.dist[i] <- max(simulated.dist[i], 0) # Ensure non-negative
        }
        q = simulated.dist[i] / 6378137 # Incorporate Earth's radius (m)
        
        # Simulate a turn angle value from ECDF
        prob_step <- runif(1, 0, 1)
        if(!autocorr.turn.angle){
          simulated.turns[i] <- approx(x = df.sub.distr$turn_angle_dist[[1]][[2]], 
                                       y = df.sub.distr$turn_angle_dist[[1]][[1]], 
                                       xout= prob_step, 
                                       rule=2)$y
        } else {
          sampled_turn_angle <- approx(x = df.sub.distr$turn_angle_dist[[1]][[2]], 
                                       y = df.sub.distr$turn_angle_dist[[1]][[1]], 
                                       xout= prob_step, 
                                       rule=2)$y
          
          # Update turn angle with autocorrelation
          simulated.turns[i] <- (df.sub.acf$acf_turn_angle * simulated.turns[i-1] + 
                                   (1 - df.sub.acf$acf_turn_angle) * sampled_turn_angle)
          simulated.turns[i] <- max(simulated.turns[i], 0) # Ensure non-negative
        }
        
        # Sample a new target heading and turn angle from their respective distributions
        target_heading <- simulated.head[i]
        turn_angle <- simulated.turns[i]
        
        if (user_choice == "target_heading") {
        
          # New Approach: Use heading from the distribution, apply turn angle as a natural deviation
          # Adjust turn angle direction to minimize angular deviation from previous heading
          left_heading <- (h - turn_angle) %% 360
          right_heading <- (h + turn_angle) %% 360
          
          # Calculate deviations, correcting for circular nature
          left_diff <- abs(target_heading - left_heading)
          left_diff <- ifelse(left_diff > 180, 360 - left_diff, left_diff)
          
          right_diff <- abs(target_heading - right_heading)
          right_diff <- ifelse(right_diff > 180, 360 - right_diff, right_diff)
          
          # Apply the turn angle in the direction that brings the agent closer to the target
          # Apply the turn angle in the direction that brings the agent closer to the target
          if (left_diff < right_diff) {
            h <- left_heading
          } else if (left_diff > right_diff) {
            h <- right_heading
          } else {
            # Deviations are equal; randomize turn direction
            h <- sample(c(left_heading, right_heading), 1)
          }
          h <- (h + 360) %% 360 # Ensure h is always in the range 0-360
          
          # Calculate the probabilities inversely proportional to the deviations
        } else if (user_choice == "probability_turn") {
          left_heading <- (h - turn_angle) %% 360
          right_heading <- (h + turn_angle) %% 360
          left_diff <- abs(target_heading - left_heading)
          left_diff <- ifelse(left_diff > 180, 360 - left_diff, left_diff)
          right_diff <- abs(target_heading - right_heading)
          right_diff <- ifelse(right_diff > 180, 360 - right_diff, right_diff)
          prob_left <- 1 / (left_diff + 1e-6)  # Adding a small number to prevent division by zero
          prob_right <- 1 / (right_diff + 1e-6)
          # Normalize probabilities
          total_prob <- prob_left + prob_right
          prob_left <- prob_left / total_prob
          prob_right <- prob_right / total_prob
          # Sample turn direction based on probabilities
          turn_direction <- sample(c("left", "right"), size = 1, prob = c(prob_left, prob_right))
          # Apply the turn
          if (turn_direction == "left") {
            h <- left_heading
          } else {
            h <- right_heading
          }
          
          } else if (user_choice == "independent_turn") {           # Random turn direction
          h <- (target_heading + turn_angle * sample(c(-1, 1), 1)) %% 360
        } else{
         h <- target_heading # No turn
        }
          
        # Dead Reckoning (DR) within boundary check
        next_lat <- asin(sin(DR.lat[i-1]) * cos(q) + cos(DR.lat[i-1]) * sin(q) * cos(h * pi/180)) 
        next_lon <- DR.lon[i-1] + atan2(sin(h * pi/180) * sin(q) * cos(DR.lat[i-1]), cos(q) - sin(DR.lat[i-1]) * sin(next_lat))
        # Convert to degrees for boundary checking
        next_lat_deg <- next_lat * 180 / pi
        next_lon_deg <- next_lon * 180 / pi
        
        # Check if the next position is within boundaries
        if (next_lat_deg >= lat_min && next_lat_deg <= lat_max && next_lon_deg >= lon_min && next_lon_deg <= lon_max) {
          # Update DR position if within boundaries
          DR.lat[i] <- next_lat
          DR.lon[i] <- next_lon 
        } else {
          print("Extent of grid boundary hit")
           if (length(unique_levels) == 2 && all(sort(unique_levels) == c("outbound", "inbound"))) {
             if (next_lat_deg >= lat_max) { sector_type <- "inbound" } 
           }
          # If out of bounds, apply a loop until a new direction keeps it within bounds
          while (next_lat_deg < lat_min || next_lat_deg > lat_max || next_lon_deg < lon_min || next_lon_deg > lon_max) {
            # Generate a deviation based on the inverted heading using von Mises distribution
            inverted_heading <- (h + 180) %% 360  # Calculate the inverted heading
            deviation <- as.numeric(rvonmises(1, circular(inverted_heading * pi / 180), 4)) * 180 / pi  # Kappa = 4, Convert to degrees
            h <- (deviation) %% 360  # Apply the deviation
            
            # Recalculate next position with adjusted heading
            next_lat <- asin(sin(DR.lat[i-1]) * cos(q) + cos(DR.lat[i-1]) * sin(q) * cos(h * pi/180)) 
            next_lon <- DR.lon[i-1] + atan2(sin(h * pi/180) * sin(q) * cos(DR.lat[i-1]), cos(q) - sin(DR.lat[i-1]) * sin(next_lat))
            next_lat_deg <- next_lat * 180 / pi
            next_lon_deg <- next_lon * 180 / pi
          }
          # After finding a suitable heading, update the DR position
          DR.lat[i] <- next_lat
          DR.lon[i] <- next_lon
        }
        
        # Calculate distance moved in this step
        dist[i] = disty(DR.lon[i]*180/pi, DR.lat[i]*180/pi, DR.lon[i-1]*180/pi, DR.lat[i-1]*180/pi) # Step-wise distance moved (m)
        straight_line_dist[i] <- disty(DR.lon[i]*180/pi, DR.lat[i]*180/pi, DR.lon[1]*180/pi, DR.lat[1]*180/pi) #Straight-line distance moved (m) 
        cum.dist = cum.dist + dist[i] # Update cumulative distance
        grid_x[i] <- grid_x_current
        grid_y[i] <- grid_y_current
        grid_x.utm[i] <- grid_x_current_m
        grid_y.utm[i] <- grid_y_current_m
        
        # Check if sector switch is needed
        if (length(unique_levels) == 2 && all(sort(unique_levels) == c("outbound", "inbound"))) {
          if (switch_based_on == "cumulative") {
            if (cum.dist >= max.distance.moved * switch_proportion) {
              sector_type <- "inbound"
            }
          } else if (switch_based_on == "straight_line") {
            # Switch based on straight-line distance from starting point
            if (straight_line_dist[i] >= switch_distance) {
              sector_type <- "inbound"
            }
          }
        }
        
        if(plot == TRUE){
        # Condition to check if 60 seconds have passed since the last plot update
        if ((cum.dist - last_dist_update) >= 5000) {
          
          # Plot the current position on the map
          # Plot the current position and path for this agent
          points(DR.lon[i] * 180 / pi, DR.lat[i] * 180 / pi, col = agent_color, pch = 16)  # Point for the current position
          lines(DR.lon[1:i] * 180 / pi, DR.lat[1:i] * 180 / pi, col = agent_color)         # Line showing path
          
          # Update time and distance trackers
          last_update_time <- Sys.time()
          last_dist_update <- cum.dist
        }
        }
        
        # Increment step index
        i = i + 1

        # Print status
        print(paste0("Agent:", a, ". ",
                     "Cumulative distance travelled: ", round(cum.dist), " m. ",
                     "Straight-line distance from start: ", round(straight_line_dist[i-1]), " m. ",
                     "Current sector: ", round(grid_x_current_m), ", ", round(grid_y_current_m), ". ",
                     sector_type))
        flush.console()
       
        
      } # End of while loop
      
      # Update seed if multiple agents
      if(!is.null(seed)){
        seed = seed + 1
        set.seed(seed) 
      }
      # Reset to outbound
      sector_type <- "outbound"
      
      # Add results into 'df_agents'
      df_agents = rbind(df_agents, data.frame(
        long = DR.lon*180/pi, 
        lat = DR.lat*180/pi, 
        id = rep(a, length(DR.lon)), 
        cum.dist = cumsum(dist), 
        straight_line_dist = straight_line_dist,
        grid_x, grid_y, 
        grid_x.utm, grid_y.utm))  
      
    }
  
    #Descriptive parameters
    # Compute the values
    df_agents <- df_agents %>% group_by(id) %>% mutate(step.length = disty(lag(long), lag(lat), long, lat)) %>% ungroup()
    df_agents <- df_agents %>% group_by(id) %>% mutate(heading = beary(lag(long), lag(lat), long, lat)) %>% ungroup()
    df_agents$heading <- ifelse(df_agents$heading < 0, df_agents$heading + 360, df_agents$heading)
    df_agents <- df_agents %>% group_by(id) %>%  mutate(turn.angle = lag(heading) - heading) %>% ungroup()
    df_agents$turn.angle <- ifelse(df_agents$turn.angle < -180, df_agents$turn.angle + 360, df_agents$turn.angle)
    df_agents$turn.angle <- ifelse(df_agents$turn.angle > 180, df_agents$turn.angle - 360, df_agents$turn.angle)   
    ### Plot of the main result ###
    # Convert DR coords to UTM
    grid_sf <- st_as_sf(df_agents, coords = c("long", "lat"), crs = 4326)
    grid_utm <- st_transform(grid_sf, crs = crs)
    # Extract UTM coordinates and add them back to the agents data frame
    df_agents$long.utm <- st_coordinates(grid_utm)[,1]
    df_agents$lat.utm<- st_coordinates(grid_utm)[,2]
    x_gridlines <- seq(lon_min.m, lon_max.m, length.out = (lon_max.m - lon_min.m) / cell_size_x) 
    y_gridlines <- seq(lat_min.m, lat_max.m, length.out = (lat_max.m - lat_min.m) / cell_size_y)   
    
    print(p2 <- ggplot(data = df_agents, aes(x = long.utm, y = lat.utm))+
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#4C00FF") +
      theme_bw() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", size = 14),
            axis.title.y = element_text(color = "black", size = 14)) +  
      stat_density2d(aes(fill = ..level..), geom = "polygon", n = 100, contour = TRUE) +
      scale_fill_gradientn(colors = topo.colors(a)) + guides(fill = "none") + 
      #Annotate horizontal gridlines
      annotate("segment", x = -Inf, xend = Inf, y = y_gridlines, yend = y_gridlines, color = "grey90", linetype = "dashed") +
      # Annotate vertical gridlines
      annotate("segment", x = x_gridlines, xend = x_gridlines, y = -Inf, yend = Inf, color = "grey90", linetype = "dashed") +
      geom_path(data = df_agents, aes(x = long.utm, y = lat.utm, group = factor(id), colour = factor(id))) +
      xlab("Relative distance E-W") + ylab("Relative distance N-S") + guides(colour = "none")  + coord_cartesian(xlim = c(min(x_gridlines), max(x_gridlines)), 
                                                                                                                 ylim = c(min(y_gridlines), max(y_gridlines))))

    # Return data frame
    return(df_agents)
    })
  
}

########## End of function ###########
#############################################################################################################################################################################################################

# e.g.,
sim = Gundog.sim(ID = df$ID,
                 datetime = as.POSIXct(df$DateTime, format="%Y-%m-%d %H:%M:%S"),
                 heading = df$Mag.heading.smoothed,
                 turn.angle = df$abs.step.ang,
                 step.length = df$step.dist,
                 journey_phase = df$trip_status,
                 lo = -63.8653 ,
                 la = -42.08 ,
                 grid.x = df$gridded.X.d,
                 grid.y = df$gridded.Y.d,
                 min.sector.data = 5,
                 bandwidth = "ucv",
                 quantile.turn.angle = 1,
                 quantile.step.length = 0.9995,
                 quantile.scope = "global", 
                 radius = 2,
                 step_size = 1,
                 autocorr.turn.angle = FALSE,
                 autocorr.step.length = FALSE,
                 autocorr.heading  = TRUE,                  
                 auto.corr.heading.method = "circ.corr",            
                 autocorr.scope = "per_grid", 
                 seed = sample(150, 1),
                 agents = 1,
                 max.distance.moved = 370000,
                 switch_proportion = 0.54,
                 switch_based_on = "cumulative",
                 switch_distance <- 50000, 
                 head.offset.out = 0,
                 head.offset.in = 0,
                 num_ID_heading = c("P10C", "P12E", "P14A", "P5A", "P9E", "P2A"),
                 use_global_phase_only = FALSE,
                 user_choice = "target_heading",
                 plot = TRUE)

