# Assignment_1 Bartek Puzio 1146247 ####

# Load all packages for this script
# I took the time to add proper headers for each section of this assignment, please open and collapse each section as needed for a structured experience.
library(tidyverse)
library(tidyr)
library(vegan)
library(readr)
library(dplyr)
library(ggplot2)
library(maps)
library(styler)

## Housekeeping ####

# Load data file into Rstudio, remove the # from the next line and run styler if you want because I already did:)
# setwd(Assignment_1.zip/Assignment_1/R)
# hopefully the working directory is appropriately set. The pathways used to call the data file are sound but needs to be set in the correct directory to work. I use the environment to find the file with the R script and using the blue gear "more" setting, set working directory.
# styler::style_file("Assignment_1.R")
data.b <- read_tsv(file = "../data/water_bear.tsv")




# Data manipulation + cleaning ####
# The Coordinates in the dataset columns had square brackets, a comma, and were not separated into 2 columns by latitude and longitude. This code removes the special characters and separates the coordinates into 2 columns and recognizes the values as numeric. The original cood column remains in the dataset for reproducabiltiy

data.b <- data.b %>%
  separate(coord, into = c("lat", "long"), sep = ",", remove = FALSE)

data.b$lat <- gsub("\\[", "", data.b$lat)
data.b$long <- gsub("\\]", "", data.b$long)

data.b$lat <- as.numeric(data.b$lat)
data.b$long <- as.numeric(data.b$long)

## Graph - Specimen collection site on Globe ####
# Raw data inputted into the map to identify collection hotspots and investigate distribution patterns by region. Collection site does not equal distribution patterns. This graph was inspired by the maps package using stackoverflow for help.  https://stackoverflow.com/questions/23130604/plot-coordinates-on-map
world_map <- map_data("world")

ggplot() +
  geom_polygon(
    data = world_map, aes(x = long, y = lat, group = group),
    fill = "lightblue", color = "gray70", alpha = 0.5
  ) +
  geom_point(
    data = data.b, aes(x = long, y = lat),
    color = "red", alpha = 0.5
  ) +
  labs(title = "Specimen Locations on Globe", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

## Data objects for Graph 2 ####
# Statistics

data.b.ratio <- data.b %>% filter(!is.na(lat), lat > 0)

n_bins <- n_distinct(data.b.ratio$bin_uri, na.rm = TRUE)
n_species <- n_distinct(data.b.ratio$species, na.rm = TRUE)
ratio <- n_bins / n_species

c(n_bins = n_bins, n_species = n_species, ratio = ratio)

# Result: 1.81. This suggests that there are more distinctive clusters, or mitochondrial lineages, than named species in this dataset. This can be an indicator of cryptic species diversity or populations that have been geographically isolated for some time.



data.b.clean <- data.b %>%
  filter(
    !is.na(bin_uri),
    !is.na(`country/ocean`),
    !is.na(species),
    `country/ocean` != "Unrecoverable",
    !is.na(lat),
    lat >= 0
  )

bin_country <- data.b.clean %>%
  group_by(`country/ocean`) %>%
  summarise(
    n_specimens = n(),
    n_bins = n_distinct(bin_uri),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_bins))


## Graph - BIN richness vs Named species from Countries North of the Equator####
bin_country %>%
  slice_max(n_bins, n = 15) %>%
  select(`country/ocean`, n_bins, n_species) %>%
  pivot_longer(
    cols = c(n_bins, n_species),
    names_to = "metric", values_to = "value"
  ) %>%
  ggplot(aes(x = reorder(`country/ocean`, value), y = value, fill = metric)) + # https://r4ds.hadley.nz/data-visualize for help and Chatgpt Ai for debugging
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Country",
    y = "Count",
    fill = "",
    title = "BIN Richness vs. Named Species (Top 15 Countries, Equator & North)"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


## Graph - Distribution of BIN's by Latitude North of the Equator ####
# https://ggplot2.tidyverse.org/reference/geom_histogram.html

data.b %>%
  filter(!is.na(lat), lat > 0) %>%
  ggplot(aes(x = lat)) +
  geom_histogram(
    bins = 100,
    fill = "blue",
    color = "white"
  ) +
  labs(
    title = "Distribution of BINs by Latitude (North of Equator)",
    x = "Latitude (°N)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

## Graph - distribution of unique bins north of the equator by latitude  ####
# https://ggplot2.tidyverse.org/reference/geom_histogram.html
data.b %>%
  distinct(bin_uri, .keep_all = TRUE) %>%
  filter(lat > 0) %>%
  ggplot(aes(x = lat)) +
  geom_histogram(
    bins = 100,
    fill = "blue",
    color = "white"
  ) +
  labs(
    title = "Distribution of Unique BINs North of the Equator Latitude",
    x = "Latitude (°)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

## Graph - BIN richness across latitude bands north of the equator ####

data.b.bands <- data.b %>%
  filter(!is.na(lat), lat >= 0) %>%
  mutate(
    lat_band = cut(
      lat,
      breaks = seq(0, 90, by = 10),
      include.lowest = TRUE, right = FALSE,
      labels = paste0(seq(0, 80, 10), "-", seq(10, 90, 10))
    )
  )

richness.bands <- data.b.bands %>%
  filter(!is.na(bin_uri)) %>%
  group_by(lat_band) %>%
  summarise(BIN_richness = n_distinct(bin_uri), .groups = "drop")

ggplot(richness.bands, aes(x = lat_band, y = BIN_richness)) +
  geom_col(fill = "blue") +
  theme_minimal() +
  labs(
    title = "Tardigrade BIN richness across latitude bands \n (Northern Hemisphere)",
    x = "Latitude band (°N)",
    y = "Number of unique BINs"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 200))

