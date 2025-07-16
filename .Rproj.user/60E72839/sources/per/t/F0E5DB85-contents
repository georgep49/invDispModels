# https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

library(tidyverse)
library(rgbif)

dispersalTx <- read_csv("gbif/dispersalDistanceSumm_060125.csv")

gbif_taxon_keys <- 
  pull(dispersalTx, species) |>
  unique() |>
  name_backbone_checklist(strict = TRUE) |> # match to backbone
  filter(matchType == "EXACT") |>
  distinct()

gbif_key_lu <- select(gbif_taxon_keys, 
        searchSpecies = verbatim_name, 
        taxonKey = usageKey, 
        speciesKey = speciesKey)

max_dd <- dispersalTx |>
  group_by(species) |>
  summarise(dd = max(log10mdd_best)) |>
  ungroup()

occ_download(
    pred_in("taxonKey", gbif_taxon_keys$usageKey),
    pred("country", "NZ"),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user = "georgep",
    pwd = "zhM6ZuhhXRmzh6K",
    email = "george.perry@auckland.ac.nz")

d <- occ_download_get("0057162-241126133413365", overwrite = TRUE) %>%
  occ_download_import()

write_csv(d, file = "rawGBIFDownload_060125.csv")
# GBIF Occurrence Download https://doi.org/10.15468/dl.d6esbz Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2025-01-05

######################
# remove records without coordinates
library(tidyverse)
library(sf)
library(southernMaps)

load(".gbif/rangeAnalysis.RData")

# d <- read_csv("./gbif/rawGBIFDownload_060125.csv")

d <- d %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

####

dat_df <- data.frame(d)
dat_sf <- dat_df |> 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# This is the points in the main NZ archipelago
dat_nz_sf <- st_join(dat_sf, nzHigh84, join = st_within) |>
  filter(!is.na(name))

# Getting the range of each species
dat_range <- dat_nz_sf |>
  st_coordinates() |>
  bind_cols(species = dat_nz_sf$scientificName, 
            year = dat_nz_sf$year,
            taxonKey = dat_nz_sf$taxonKey,
            speciesKey = dat_nz_sf$speciesKey) |>
  left_join(gbif_key_lu, by = "speciesKey")


# and summarising it
dat_range_summ <- dat_range |>
  group_by(searchSpecies) |>
  summarise(lat_max = max(Y), lat_min = min(Y), lat_range = max(Y) - min(Y), n = n()) |>
  ungroup() |>
  left_join(max_dd, by = c("searchSpecies" = "species"))

write_csv(dat_range_summ, file = "rangeSummaryInvasiveDispersal_060125.csv")

  
#### Visualisation code
load("./gbif/rangeAnalysis.RData")
dat_range_summ$dd_raw <- 10 ^ dat_range_summ$dd

lbl <- dat_range_summ |>
  filter(n < 100 & lat_range < 5 & dd_raw > 100)

# ggplot() +
#   geom_sf(data = nzHigh84) +
#   geom_sf(data = dat_nz_sf)
set.seed(3141)
rangePlot <- ggplot(dat_range_summ) +
  geom_point(aes(x = n, y = lat_range, col = dd), size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(data = lbl, 
                          aes(x = n, y = lat_range, label = searchSpecies),
    min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, size = 3, fontface = "italic") +
  labs(x = "No. of records", y = "Latitudinal range (degrees)") +
  scale_colour_distiller(name = "Max dispersal distance (m)", palette = "YlGnBu", direction = 1, breaks = 0:4, labels = 10^(0:4)) +
  scale_x_continuous(breaks = c(0, 10^(0:5)), transform = "pseudo_log", 
 guide = "axis_logticks") +
  theme_bw()

pdf(file = "gbif/rangeAnalysisAllRestricted_v2.pdf")
rangePlot
dev.off()

save.image("D:/Research/wottonDispersalModels/gbif/rangeAnalysis.RData")

#####
wm <- borders("world", regions = "New Zealand", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = dat %>% filter(speciesKey == 3974372),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "blue",
             size = 0.5) +
  # xlim(165, 180) +
   ylim(-50, -25) +
  theme_bw()

