library(tidyverse)
library(taxize)

pinops_id <- get_gbifid("Pinopsida")
cList <- downstream(pinops_id, "gbif", downto = "species")
