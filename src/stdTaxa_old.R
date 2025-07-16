library(tidyverse)
library(readxl)
library(taxize)

load("data/dispeRsal.rda")
load("debTaxonomy.Rdata")

deb_data <- read_xlsx("data/weed dispersal data 2023-11-10.xlsx")

# First get standardised taxa lists via <<taxize>>



# Resolve names from the GNR server

# "Erythrina xsykesii" [145]
# "Passiflora xrosea" [255]
# "Reynoutria xbohemica"[306] 
# "Salix xfragilis" [332]
# "Sporobolus xtownsendii" [359]


# standardised list of names via Kew's checklist
# drop 145, 255, 306, 332, 359
names_res <- gnr_resolve(deb_data$species [c(1:144, 146:254, 256:305, 307:331, 333:358, 360:376)]) %>%
  group_by(submitted_name) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(matched_name = word(matched_name, 1,2, sep = " "))

deb_data_j <- deb_data %>%
     left_join(names_res %>% select(submitted_name, matched_name), by = join_by(species == submitted_name))

question_names <- deb_data_j %>%
  filter(is.na(matched_name)) %>%
  pull(species)

question_names <- c(question_names, deb_data$species[c(145, 255, 306, 332, 359)])
  

# Classification via taxonize
cl_gbif <- classification(deb_data$species, db = 'gbif', rows = 1) # rwos here to stop ia sel


deb_data_filt <- deb_data_j %>% filter(!is.na(matched_name))
cl_ncbi <- classification(deb_data_filt$matched_name, db = 'ncbi')

cl_ncbi_rb <- rbind(cl_ncbi)
spp_ncbi <- unique(cl_ncbi_rb$query)
not_ncbi <- setdiff(deb_data_filt$species, spp_ncbi)


cl_gbif <- classification(question_names, db = 'gbif')

cl_gbif_rb <- rbind(cl_gbif)
cl_gbif_rb$source <- "gbif"


cl_all <- rbind(cl_ncbi_rb, cl_gbif_rb) 

cl_all <- cl_all %>% 
  filter(rank %in% c("family", "order")) %>%
  group_by(query) %>%
  
  
  



save.image("debTaxonomy.RData")
               