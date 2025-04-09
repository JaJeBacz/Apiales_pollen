require(knitr)
knitr::opts_chunk$set(cache=TRUE, dev="pdf")
require(Momocs) 
require(phytools)
require(car)
require(boxplotdbl)
require(reshape2)
require(multcompView)
require(ggplot2)
require(dplyr)
require(coin)
require(kableExtra)
require(rgl)
require(pracma)
require(reticulate)
require(plotly)
require(htmlwidgets)

source("functions.R") #Various mildly-helpful functions

####################
#### POLAR VIEW ####
####################

#List the jpg files

jpg.list_all <- list.files("outlines/polar",pattern = ".jpg", full.names = T)

#Import them

returns_pollen_polar_view <- import_jpg(jpg.paths = jpg.list_all, auto.notcentered = T) 

#Create Out object

pollen <- Out(returns_pollen_polar_view)

#Import grouping variable which is in the folder with the outlines

groups_pollen <- read.csv("outlines/polar/groups.csv", sep = ";")
groups_pollen<-data.frame(groups_pollen[,1:5])

# Assign grouping variable

pollen$fac <- groups_pollen
names(pollen) <- groups_pollen$ID

#Center and template pollen grains.

coo_pollen_centered <- coo_center(pollen) 
coo_pollen_template <- coo_template(coo_pollen_centered) 

#Panel of outlines in polar view.

panel(pollen, col=0, names="ID", cex=0.3)

#Elliptical Fourier Analysis and PCA 

efou_pollen <- efourier(coo_pollen_template, norm=F, nb.h=32, smooth.it=1)
efou_pollen <- rm_sym(efou_pollen) # remove noise

# PCA of the eFa harmonic coefficients

pca_efou_pollen <- PCA(efou_pollen)

#Here we can see the contribution of each component

pc_scores_po <- summary(pca_efou_pollen)
kable_styling(kable(pc_scores_po$importance[,1:5],"latex"),latex_options="hold_position")

#Organize PC scores of polar view

all_data_pcs_efou_pollen_polar_view <- cbind(groups_pollen[1:5], 
                                             data.frame(pca_efou_pollen$x[,1:10])) 

#Contribution to shape in polar view of each PC.

PCcontrib(pca_efou_pollen, nax=1:2, sd.r=c(-3, -1.5, 0, 1.5, 3)) 

#Plot morphospace resulting from the scores of the first two PCs.

plot_PCA (pca_efou_pollen, ~family, axes = c(1, 2), chullfilled = TRUE, legend = TRUE, chull = FALSE, labelpoints = TRUE, zoom = 0.9)

#########################
#### EQUATORIAL VIEW ####
#########################

#First we list the jpg files

jpg.list_all <- list.files("outlines/equatorial/",pattern = ".jpg", full.names = T)

#Import them

returns_pollen_equatorial_view <- import_jpg (jpg.path = jpg.list_all, auto.notcentered = T) 

#Create Out object

pollen_eq <- Out(returns_pollen_equatorial_view)

#Import grouping variable which is in the folder with the outlines

groups_pollen_equatorial_view <- read.csv("outlines/equatorial/groups.csv", sep = ';')
groups_pollen_equatorial_view <- data.frame(groups_pollen_equatorial_view[,1:5])

#Assign grouping variable

pollen_eq$fac <- groups_pollen_equatorial_view
names(pollen_eq) <- groups_pollen_equatorial_view$ID # One can use the species column as well

#Here is the panel of outlines in equatorial view.

panel(pollen_eq, col=0, names="ID", cex=0.3)

#Importance  of each PC 

#Elliptical Fourier Analysis and PCA 

efou_pollen_eq_view <- efourier(pollen_eq, norm=T, nb.h=32, smooth.it=1)
efou_pollen_eq_view <- rm_asym(efou_pollen_eq_view) # Remove noise 

#PCA of the eFa harmonic coefficients

pca_efou_pollen_eq_view <- PCA(efou_pollen_eq_view,  center=T)

#Contribution of each component

pc_scores_eq <- summary(pca_efou_pollen_eq_view)
kable_styling(kable(pc_scores_eq$importance[,1:5],"latex"),latex_options="hold_position")

#Organize PC scores.

all_data_pcs_efou_pollen_eq_view <- cbind(groups_pollen_equatorial_view[1:5], 
                                          data.frame(pca_efou_pollen_eq_view$x[,1:10])) 

#Morphological variation explained in each PC. 

PCcontrib(pca_efou_pollen_eq_view, nax=1:2, sd.r=c(-3, -1.5,0,1.5, 3)) 

#Plot PCA

plot_PCA (pca_efou_pollen_eq_view, ~family, axes = c(1, 2), chullfilled = TRUE, legend = TRUE, chull = FALSE, labelpoints = TRUE, zoom = 0.9)

##########################
#### Raw measurements ####
##########################

size_data_exp <- read.table("data/measurements", header = TRUE, sep = ";")

#######################################################
#### Prepare final table with continuous variables ####
#######################################################

size_data_fin <- cbind.data.frame(
  size_data_exp[,1:6],
  "P/E" = size_data_exp$P/size_data_exp$E,
  "WP/P" = (size_data_exp$SP+size_data_exp$NP)/size_data_exp$P,
  "WE/P" = (size_data_exp$SE+size_data_exp$NE)/size_data_exp$P,
  "NP/WP" = size_data_exp$NP/(size_data_exp$SP+size_data_exp$NP),
  "NE/WE" = size_data_exp$NE/(size_data_exp$SE+size_data_exp$NE),
  "COSW/NP" = size_data_exp$COSW/size_data_exp$NP,
  "ECT/P" = size_data_exp$ECT/size_data_exp$P,
  "COSL/NP" = size_data_exp$COSL/size_data_exp$NP,
  "ENDL/P" = size_data_exp$ENDL/size_data_exp$P,
  "ENDL/ENDW" = size_data_exp$ENDL/size_data_exp$ENDW
)

#Plot examples of solids of revolution
Turgenia_latifolia <- import_jpg("outlines/equatorial/100.jpg")
Turgenia_latifolia_centered <- coo_center(Out(Turgenia_latifolia))
Turgenia_latifolia_plot <- plot_solid_of_revolution(Turgenia_latifolia_centered)
saveWidget (Turgenia_latifolia_plot, file = "plots/Turgenia_sor.html")

Ammoides_pusilla <- import_jpg("outlines/equatorial/112.jpg")
Ammoides_pusilla_centered <- coo_center(Out(Ammoides_pusilla))
Ammoides_pusilla_plot <- plot_solid_of_revolution(Ammoides_pusilla_centered)
saveWidget (Ammoides_pusilla_plot, file = "plots/Ammoides_sor.html")

Bifora_radians <- import_jpg("outlines/equatorial/145.jpg")
Bifora_radians_centered <- coo_center(Out(Bifora_radians))
Bifora_radians_plot <- plot_solid_of_revolution(Bifora_radians_centered)
saveWidget (Bifora_radians_plot, file = "plots/Bifora_sor.html")

Griselinia_lucida <- import_jpg("outlines/equatorial/3.jpg")
Griselinia_lucida_centered <- coo_center(Out(Griselinia_lucida))
Griselinia_lucida_plot <- plot_solid_of_revolution(Griselinia_lucida_centered)
saveWidget (Griselinia_lucida_plot, file = "plots/Griselinia_sor.html")

#Calulate pollen volume by transforming half-image into a solid of revolution
scaling <- size_data_exp[,c(1,6)]
scaling <- scaling [match(pollen_eq$species, scaling$species),]

#Rescale
volumes_rescaled <- mapply(calculate_volume_of_revolution_x_axis_rescaled, 
                           pollen_eq$coo, 
                           scaling$P)

volume_df_rescaled <- data.frame(
  species = groups_pollen_equatorial_view$species, 
  volume_rescaled = volumes_rescaled
)

volume_df_rescaled <- volume_df_rescaled [match(size_data_fin$species, volume_df_rescaled$species),]

#Get shapes
shapes <- cbind.data.frame(species = all_data_pcs_efou_pollen_eq_view$species,
                           P.SHAPE = all_data_pcs_efou_pollen_polar_view[,6],
                           E.SHAPE = all_data_pcs_efou_pollen_eq_view[,6])
shapes <- shapes [(match(size_data_fin$species, shapes$species)),]

#Bind final table
size_data_fin$VOL <- volume_df_rescaled$volume_rescaled

size_data_fin <- cbind.data.frame(
  size_data_fin,
  shapes [,2:3]
)

write.table(size_data_fin, file = "data/morpho", 
            sep = "\t", quote = F, row.names = FALSE)

#############
#### PCA ####
#############

#Create dataset for plotting
pca_data <- cbind.data.frame(PC1eq = pca_efou_pollen_eq_view$x[,1],
                             PC1pol = pca_efou_pollen$x[,1],
                             VOL = volume_df_rescaled$volume_rescaled,
                             clade = groups_pollen_equatorial_view$family,
                             species = groups_pollen$species)

pca_data <- pca_data [match(size_data_exp$species, pca_data$species),]
pca_data$family <- size_data_exp$clade

#Read consensus tree
contree <- read.nexus("tree/MrBayes/concatenated.nex.con.tre")
contree.rooted <- root(contree, "Helianthus_annuus", resolve.root = TRUE)
contree.rooted <- keep.tip(contree.rooted, pca_data$species)
contree.rooted$tip.label <- gsub("_", " ", contree.rooted$tip.label)

pca_data$species <- gsub("_", " ", pca_data$species)

#Make plots
phylomorpho_eq_pol <- ggmorphoJB_color(contree.rooted, pca_data, xvar = PC1eq, yvar = PC1pol,
                                       factorvar = family, labelvar = species,
                                       repel = TRUE, edge.width = 0.5,
                                       point.size = 3,
                                       fontface = "italic",
                                       tree.alpha = 0.3, title = NULL, 
                                       show.labels = F) +
  theme_minimal() +
  xlab("PC1 equatorial") + 
  ylab("PC1 polar")


ggsave("plots/phylomorpho_eq_pol.pdf", 
       phylomorpho_eq_pol, 
       device = "pdf", 
       width = 170, 
       height = 150, 
       units = "mm")

phylomorpho_eq_vol <- ggmorphoJB_color(contree.rooted, pca_data, xvar = PC1eq, yvar = VOL,
                                       factorvar = family, labelvar = species,
                                       repel = TRUE, edge.width = 0.5,
                                       point.size = 3,
                                       fontface = "italic",
                                       tree.alpha = 0.3, title=NULL, 
                                       show.labels = F) +
  theme_minimal() +
  xlab("PC1 equatorial") + 
  ylab("volume")


ggsave("plots/phylomorpho_eq_vol.pdf", 
       phylomorpho_eq_vol, 
       device = "pdf", 
       width = 170, 
       height = 150, 
       units = "mm")

phylomorpho_pol_vol <- ggmorphoJB_color(contree.rooted, pca_data, xvar = PC1pol, yvar = VOL,
                                        factorvar = family, labelvar = species,
                                        repel = TRUE, edge.width = 0.5,
                                        point.size = 3,
                                        fontface = "italic",
                                        tree.alpha = 0.3, title=NULL, 
                                        show.labels = F) +
  theme_minimal() +
  xlab("PC1 polar") + 
  ylab("volume")


ggsave("plots/phylomorpho_pol_vol.pdf", 
       phylomorpho_pol_vol, 
       device = "pdf", 
       width = 170, 
       height = 150, 
       units = "mm")