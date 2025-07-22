
library(tidyverse)
library(readxl)
library(taxize)

load("data/dispeRsal.rda")

# load("data/debTaxonomy_040924.Rdata")

dispersal_db <- read_xlsx("data/weed dispersal data 2025-07-15.xlsx")

# Tally mechanisms, syndromes, non-std mechanisms, mechanism prevalence
# nonNA in mechanism*, nonNA in syndrome*, diff is no. nonstd, tabulate across mechanism* cols

# I just checked the dispersalTable data again, and the species with unspecialised syndrome and 
# no known mechanisms (i.e. mechanism1 = unknown) have a tally of one for the number of 
# mechanisms. This should be zero.

# Species  with an unspecialised syndrome need to have all the mechanisms counted as 
# non-standard (i.e. don’t subtract unspecialised from the tally of mechanisms). 

# Except for J. vulgaris, which is the only species with two syndromes where one of them (syndrome2) 
# is unspecialised (in this case, only subtract the wind syndrome from the number of mechanisms!).

# Nearly there – the species with two syndromes need to have both syndromes subtracted from the 
# number of mechanisms to calculate the number of non-standard mechanisms.

dispersal_db <- dispersal_db |>
  rowwise() |>
  mutate(
    n_mech = sum(!is.na(c_across(starts_with("mechanism")))) - if_else(mechanism1 == "unknown", 1, 0, 0),
    n_synd = sum(!is.na(c_across(starts_with("syndrome")))),
    n_nonstd_mech = if_else(syndrome1 != "unspecialised", n_mech - n_synd, n_mech),
    n_nonstd_mech = if_else(!is.na(syndrome2) & species == "Jacobaea vulgaris", n_mech - (n_synd - 1), 
                    if_else(!is.na(syndrome2) & species != "Jacobaea vulgaris", n_mech - n_synd, n_nonstd_mech)))
# catch for JacVul    

# prevalence of the mechanisms in the flora
mech_table <- lapply(dispersal_db[,paste0("mechanism", 1:5)], function(x) {as.vector(x)}) |> 
  unlist() |>
  table() |>
  data.frame() |>
  mutate(prop = Freq / length(unique(dispersal_db$species)))

names(mech_table) <- c("mechanism", "freq", "prop")

#######
# Sort the taxonomy
# First get standardised taxa lists via <<taxize>>

# rows = 1 here to stop ia sel
cl_gbif <- classification(dispersal_db$species, db = "gbif", rows = 1)
cl_gbif_rb <- rbind(cl_gbif)

# extract family and order for PGLS models
cl_all <- cl_gbif_rb |> 
  filter(rank %in% c("family", "order"))

cl_all_pw <- cl_all |>
  pivot_wider(id_cols = c(query), names_from = rank, values_from = name)

# join back in dispersal database
dispersal_db_tax <- dispersal_db |>
  left_join(cl_all_pw, by = join_by(species == query))

save.image("data/debTaxonomy_250716.RData")


###########
# Get data into format for dispeRsal
# need DS (dispersal syndrome), growth-form (GF), release height (RH, log10), seed mass (SM), terminal vel (TV)
# Seed.mass = SM, DS = mechanism1, height = RH, TV = NA

load("data/debTaxonomy_250716.RData")

source("src/tidyMechanisms.r")
source("src/dispersal_gp.r")

save.image("data/debTaxonomy_250716.RData")

######
#  Now the models themselves (using Tamme code modified by GP)
model2_pred <- dispeRsal_gp(dispersal_db_tax_allMech |> filter(model == 2), model = 2, 
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model3_pred <- dispeRsal_gp(dispersal_db_tax_allMech |> filter(model == 3), model = 3, 
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model4_pred <- dispeRsal_gp(dispersal_db_tax_allMech |> filter(model == 4), model = 4,
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model5_pred <- dispeRsal_gp(dispersal_db_tax_allMech |> filter(model == 5), model = 5, 
    CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

# bind the model output together
all_predictions <- bind_rows(model2_pred$predictions,
                model3_pred$predictions,
                model4_pred$predictions,
                model5_pred$predictions,
                .id = "model_number") |>
      mutate(model_number = as.numeric(model_number) + 1)
      


# tidy up and identify the best prediction for each on nested random family model
all_predictions <- all_predictions |>
  mutate(Species = str_replace(Species, "_mech.", "")) |>
  distinct() |>
  as_tibble()  |>
  mutate(log10MDD_best =
    case_when(
      !is.na(log10MDD_Family) ~ log10MDD_Family,
      !is.na(log10MDD_Order) & is.na(log10MDD_Family) ~ log10MDD_Order,
      is.na(log10MDD_Order)  & is.na(log10MDD_Family) ~ log10MDD,
      .default = NA))

all_predictions <- all_predictions |>
  janitor::clean_names() |>
  left_join(dispersal_db |> select(-c(1, 3, 4, 5)))

## Write to file...

write_csv(dispersal_db, "output/dispersal_db.csv")

write_csv(dispersal_db_tax, "output/dispersalTable_2307125.csv")
write_csv(all_predictions, "output/dispersalDistance_2307125.csv")
write_csv(mech_table, "output/mechanismsTable_2307125.csv")



######

# Also, we were going to try modelling wind.none dispersal distances as wind.special instead, and see what the results looked like. We don’t want to include the two types as different mechanisms in our counts of mechanisms etc., but we may want to include both distances in the SI. 
# So can we just run the model first and then decide how we’re going to proceed?

dispersal_allMech_wind_spec <- dispersal_db_tax_allMech |>
  mutate(DS = str_replace(DS, "wind.none", "wind.special"))

#  Now the models themselves (using Tamme code modified by GP)
model2_ws_pred <- dispeRsal_gp(dispersal_allMech_wind_spec |> filter(model == 2), model = 2, 
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model3_ws_pred <- dispeRsal_gp(dispersal_allMech_wind_spec |> filter(model == 3), model = 3, 
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model4_ws_pred <- dispeRsal_gp(dispersal_allMech_wind_spec |> filter(model == 4), model = 4,
      CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

model5_ws_pred <- dispeRsal_gp(dispersal_allMech_wind_spec |> filter(model == 5), model = 5, 
    CI = FALSE, random = TRUE, tax = "family", write.result = FALSE)

# bind the model output together
all_ws_predictions <- bind_rows(model2_ws_pred$predictions,
                model3_ws_pred$predictions,
                model4_ws_pred$predictions,
                model5_ws_pred$predictions,
                .id = "model_number") |>
      mutate(model_number = as.numeric(model_number) + 1)
      


# tidy up and identify the best prediction for each on nested random family model
all_ws_predictions <- all_ws_predictions |>
  mutate(Species = str_replace(Species, "_mech.", "")) |>
  distinct() |>
  as_tibble()  |>
  mutate(log10MDD_best =
    case_when(
      !is.na(log10MDD_Family) ~ log10MDD_Family,
      !is.na(log10MDD_Order) & is.na(log10MDD_Family) ~ log10MDD_Order,
      is.na(log10MDD_Order)  & is.na(log10MDD_Family) ~ log10MDD,
      .default = NA))

all_ws_predictions <- all_ws_predictions |>
  janitor::clean_names() |>
  left_join(dispersal_db |> select(-c(1, 3, 4, 5)))

## Write to file...
write_csv(all_ws_predictions, "output/dispersalDistanceWindSpec_2307125.csv")
