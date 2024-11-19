
# Gundog.sim

## Description

**Gundog.sim** is an individual-based model that predicts the movement paths of animals by simulating their movement using dead-reckoning (DR) calculations through a spatial grid system. The model uses high-resolution dead-reckoning data to inform movement decisions, incorporating empirical distributions of fundamental step lengths (Fstepdistance), turn angles (Fturnangle), and compass headings (H). By dividing the movement landscape into grid cells and segmenting movement data by journey phase (outbound or inbound), the model dynamically adjusts agent movements to reflect realistic behavioral patterns observed in empirical studies.

### Key Features

- **Empirical Distributions:** Incorporates empirical frequency distributions (Empirical Cumulative Distribution Functions, ECDFs) for compass heading, step lengths and turn angles specific to each grid cell and journey phase.
- **Journey Phases:** Differentiates between outbound and inbound phases, allowing for phase-specific movement behaviors.
- **Autocorrelation Adjustments:** Applies first-order autocorrelation to movement metrics to capture temporal dependencies in movement patterns.
- **Customizable Parameters:** Offers extensive customization options for distributions, journey phase handling, autocorrelation settings, and more.

## Installation

### Prerequisites

Ensure that you have **R** installed on your system. You can download R from [CRAN](https://cran.r-project.org/).

### Required R Packages

Gundog.sim relies on several R packages. You can install them using the following commands:

```R
install.packages(c("dplyr", "scales", "ggplot2", "gridExtra", "circular", "devtools", "sf", "purrr", "tidyr", "cowplot", "zoo"))
```

### Prepare Your Input Data

Ensure your data frame (df) contains the following columns (though names of columns can vary):
- **ID:** Unique identifier for each animal.
- **DateTime:** Timestamp of each observation.
- **heading:** Animal heading in degrees (0-360).
- **abs.step.ang:** Absolute turn angle in degrees.
- **step.dist:** Distance traveled between observations.
- **trip_status:** Journey phase label ("outbound" or "inbound").
- **gridded.X.d:** X-coordinate of the grid cell.
- **gridded.Y.d:** Y-coordinate of the grid cell.



### Run the Simulation

Call the Gundog.sim function with your defined parameters.

```R
# Example parameters (modify as needed)
result <- Gundog.sim(
  ID = df$ID,
  datetime = as.POSIXct(df$DateTime, format="%Y-%m-%d %H:%M:%S"),
  heading = df$heading.corr,
  turn.angle = df$abs.step.ang,
  step.length = df$step.dist,
  journey_phase = df$trip_status,
  lo = -63.8653,
  la = -42.08,
  grid.x = df$gridded.X.d,
  grid.y = df$gridded.Y.d,
  min.sector.data = 75,
  bandwidth = "ucv",
  quantile.turn.angle = 1,
  quantile.step.length = 0.9995,
  quantile.scope = "per_journey_phase",
  radius = 1,
  step_size = 1,
  autocorr.turn.angle = FALSE,
  autocorr.step.length = TRUE,
  autocorr.heading = FALSE,
  auto.corr.heading.method = "trig",
  autocorr.scope = "per_grid",
  seed = 10,
  agents = 5,
  max.distance.moved = 270000,
  switch_proportion = 0.45,
  switch_based_on = "cumulative",
  switch_distance = 50000,
  head.offset.out = 50,
  head.offset.in = -50,
  num_ID_heading = NULL,
  use_global_phase_only = TRUE,
  user_choice = "target_heading",
  plot = TRUE
)

```

## Customization Options

Users can customize the movement behaviors of agents by adjusting various parameters within the Gundog.sim function. Some of the key customization options include:
- **Quantile Filtering**
- **Autocorrelation Settings**
- **Journey Phase Handling**
- **Heading Adjustment Strategies**
- **Grid Cell Handling**

## License

This project is licensed under the MIT License.

## Contact

For questions, suggestions, or contributions, please contact:
- Richard Gunner
- Email: rgunner@ab.mpg.de
- GitHub: [Richard6195](https://github.com/Richard6195)
