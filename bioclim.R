setwd("bioclim/")

library(rgbif)
library(taxize)
library(purrr)
library(tibble)
library(dplyr)
library(magrittr)
library(roperators)
library (raster)
library (sp)
library (phytools)
library (ape)
library (ggplot2)
library (qdapTools)
library (ggrepel)
library (plot3D)
library (plot3Drgl)
library (gridExtra)
library (gridBase)
library (lattice)
library (geomorph)
library (CoordinateCleaner)
library (maptools)
library (rgdal)
library (rasterVis)
library (terra)
library (sf)
library (geodata)

#Log in to GBIF

user <- "" # your gbif.org username 
pwd <- "" # your gbif.org password
email <- "" # your email

#Download data from GBIF

file_url <- "taxa"
# match the names 
gbif_taxon_keys <- 
  readr::read_csv(file_url) %>% 
  pull("Taxon name") %>% # use fewer names if you want to just test 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(file = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED" | status == "SYNONYM") %>% # get only accepted and matched names
  filter(kingdom == "Plantae" & order == "Apiales") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys

occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred("hasCoordinate", TRUE),
  pred_in("basisOfRecord", c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN")),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_get('0012475-240506114902167') %>%
  occ_download_import()

#Clean and bind GBIF dataframe
data <- read.delim ("download.csv", header = TRUE, sep = "\t", stringsAsFactors=FALSE, quote="", fill=FALSE)
data.cleaned <- clean_coordinates(data, lon = "decimalLongitude", lat = "decimalLatitude", value = "clean")

#Correct synonyms
data.cleaned$species[data.cleaned$species == "Mutellina adonidifolia"] <- "Mutellina purpurea"
data.cleaned$species[data.cleaned$species == "Hippomarathrum vulgare"] <- "Seseli hippomarathrum"
data.cleaned$species[data.cleaned$species == "Laser archangelica"] <- "Laserpitium archangelica"
data.cleaned$species[data.cleaned$species == "Scaligeria allioides"] <- "Elaeosticta allioides"
data.cleaned$species[data.cleaned$species == "Hansenia phaea"] <- "Haplosphaera phaea"
data.cleaned$species[data.cleaned$species == "Psammogeton capillifolium"] <- "Aphanopleura capillifolia"
data.cleaned$species[data.cleaned$species == "Pachypleurum mutellinoides"] <- "Pachypleurum alpinum"
data.cleaned$species[data.cleaned$species == "Zosima absinthiifolia"] <- "Zosima orientalis"
data.cleaned$species[data.cleaned$species == "Pimpinella aromatica"] <- "Pimpinella anisum"

saveRDS (data.cleaned, "bioclim/data.cleaned.RDS")
#data.cleaned <- readRDS("bioclim/data.cleaned.RDS")

#Download worldclim data
wclim <- getData("worldclim",var = "bio", res = 2.5)
#wclim <- worldclim_global(var = "bio", path = "bioclim/", res = 2.5)
gain (wclim) = 0.1

#Extract worldclim data
wc <- data[c(22,23)]
lats <- wc$decimalLatitude
lons <- wc$decimalLongitude
coords <- data.frame(x = lons, y = lats)
points.wclim <- SpatialPoints (coords, proj4string = wclim@crs)
values.clim <- raster::extract(wclim, points.wclim, method = "simple")

#Bind data frame and aggregate by species
data.new <- cbind.data.frame(data.cleaned$species, coordinates(points.wclim), values.clim)
names (data.new) <- c ("species", 
                       "Lats", 
                       "Lons", 
                       names(data.new[,4:ncol(data.new)]))

bioclim <- aggregate(data.new[,4:ncol(data.new)], list (data.new$species), mean, na.rm = TRUE)
colnames(bioclim)[1] <- "species"
bioclim <- na.omit (bioclim)
bioclim$species <- gsub(" ", "_", bioclim$species)

morpho <- read.table("data/morpho", 
                     header = TRUE, sep = "\t")

bioclim <- merge(morpho[,1:5], bioclim, by.x = "species")
write.table (bioclim, "data/bioclim", 
             sep = "\t", quote = FALSE, row.names = FALSE)