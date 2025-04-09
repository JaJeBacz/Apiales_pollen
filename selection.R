require(mvMORPH)
require(phytools)
require(ape)
require(factoextra)
require(cluster)
require(webshot2)
require(rgl)
require(magick)
require(viridisLite)
require(dplyr)
require(purrr)
require(ggplot2)
require(tidyr)
require(ggpubr)
require(parallel)
require(reshape2)
library (PhylogeneticEM)
require (l1ou)
require (mvMORPH)
require (phangorn)
require (viridisLite)

source("functions.R") #Various mildly-helpful functions

#Read RDS just in case

folder_path <- "RDS/"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for(file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  data <- readRDS(file_path)
  assign(file_name, data, envir = .GlobalEnv)
}

#############################
#### Bioclimatic regimes ####
#############################

#Read and prepare data
bioclim <- read.table("data/bioclim", header = TRUE, row.names = 1, sep = "\t")

#Read trees

#Read consensus
contree.rooted <- read.tree("tree/TreePL/tree_mb_dated.tre")
contree.rooted <- keep.tip(contree.rooted, row.names(bioclim))

#Correct rounding (the tree is ulrramteric!)
contree.rooted <- force.ultrametric (contree.rooted, method = "extend")

#PCA
PCA <- prcomp(bioclim[,5:ncol(bioclim)], center = TRUE, scale. = TRUE)

#Determine optimal number of components
PCs <- PCA$x[,1:2]

#Variance by component
round(PCA$sdev^2 / sum(PCA$sdev^2)*100, 4)

#Make table with important components
threshold <- 0.3
PCA.vecs <- PCA$rotation[,1:2]
ifelse(abs(PCA.vecs) > threshold, PCA.vecs, NA)

########################
#### PhylogeneticEM ####
########################

t_PCs_con <- t(PCs)

phyloem.con <- PhyloEM(contree.rooted, t_PCs_con,
                       process = "scOU",
                       parallel_alpha = TRUE,
                       Ncores = 16,
                       random.root = TRUE,
                       stationary.root = TRUE)

params.con <- params_process(phyloem.con)

#Equivalent solutions

phyloem.con.eq <- equivalent_shifts(contree.rooted, params.con)

saveRDS(phyloem.con.eq, "RDS/phyloem.con.eq.RDS")

#Plotting

pdf("phyloEM.pdf", width = 20, height = 10)
plot (phyloem.con)
dev.off()

#Plotting

pdf ("plots/phyloEM_eq.pdf", width = 20, height = 10)
plot (phyloem.con.eq,
      show_shifts_values = TRUE, 
      shifts_cex = 0.5)
dev.off()

########################
######### l1OU #########
########################
l1ou.data.con <- adjust_data(contree.rooted, PCs)

l1ou.pBIC.con <- estimate_shift_configuration(
  tree = l1ou.data.con$tree,
  Y = l1ou.data.con$Y,
  criterion = c("pBIC"),
  max.nShifts = 20,
  nCores = 16,
  root.model = c("OUrandomRoot"),
  rescale = FALSE
)

saveRDS(l1ou.pBIC.con, "RDS/l1ou.pBIC.con.RDS")

l1ou.pBIC.con.boot <- l1ou_bootstrap_support (l1ou.pBIC.con, 
                                              nItrs = 200,
                                              multicore = TRUE,
                                              nCores = 10)

saveRDS(l1ou.pBIC.con.boot, "RDS/l1ou.pBIC.con.boot.RDS")

l1ou.BIC.con <- estimate_shift_configuration(
  tree = l1ou.data.con$tree,
  Y = l1ou.data.con$Y,
  criterion = c("BIC"),
  max.nShifts = 20,
  nCores = 16,
  root.model = c("OUrandomRoot"),
  rescale = FALSE
)

saveRDS(l1ou.BIC.con, "RDS/l1ou.BIC.con.RDS")

l1ou.BIC.con.boot <- l1ou_bootstrap_support (l1ou.BIC.con, 
                                             nItrs = 200,
                                             multicore = TRUE,
                                             nCores = 10)

saveRDS(l1ou.BIC.con.boot, "RDS/l1ou.BIC.con.boot.RDS")

#Plotting
bootstrap <- round(l1ou.pBIC.con.boot$detection.rate*100, digits = 1)
bootstrap [bootstrap <70] = NA

#pBIC
pdf("plots/l1ou_pBIC.pdf", height = 6.6929*2, width = 8.8976*2)
palette.esc <- palette(c(viridis(l1ou.pBIC.con$nShifts),"lightgrey"))
plot(l1ou.pBIC.con, 
     edge.label.ann=TRUE, 
     edge.shift.ann = FALSE, 
     edge.label = bootstrap, 
     edge.ann.cex = 0.5, cex=0.5, 
     label.offset=0.02, 
     edge.width=2, 
     palette = palette.esc)
dev.off()

#BIC
bootstrap <- round(l1ou.BIC.con.boot$detection.rate*100, digits = 1)
bootstrap [bootstrap <70] = NA

pdf("plots/l1ou_BIC.pdf", height = 6.6929*2, width = 8.8976*2)
palette.esc <- palette(c(viridis(l1ou.BIC.con$nShifts),"lightgrey"))
plot(l1ou.BIC.con, 
     edge.label.ann=TRUE, 
     edge.shift.ann = FALSE, 
     edge.label = bootstrap, 
     edge.ann.cex = 0.5, cex=0.5, 
     label.offset=0.02, 
     edge.width=2, 
     palette = palette.esc)
dev.off()

###########################################################
###########################################################
################## Hypotheses testing #####################
###########################################################
###########################################################

#Samplers and indices
set.seed(69420)
sampler <- sample(1:100, 1)
indices <- 1:100
treesampler <- 100

#Read and scale data
morpho.raw <- read.table("data/morpho", header = TRUE, row.names = 1, sep = "\t")
morpho <- scale(morpho.raw[,5:ncol(morpho.raw)])

#Read trees
tree1 <- read.nexus("tree/MrBayes/concatenated.nex.run1.t")
tree2 <- read.nexus("tree/MrBayes/concatenated.nex.run2.t")
tree.merge <- c (tree1 [round(length(tree1)/4):length(tree1)],
                 tree2 [round(length(tree2)/4):length(tree2)])
tree.rooted <- root.multiPhylo(tree.merge, 
                               "Helianthus_annuus", 
                               resolve.root = TRUE)
tree.sample <- sample(tree.rooted, treesampler)
tree.sample <- keep.tip.multiPhylo(tree.sample, row.names(morpho))
tree.sample <- lapply(tree.sample, ladderize)

saveRDS(tree.sample, "RDS/tree.sample.RDS")

#Define regimes with addtional shifts for core apioids
tree.sample.regimes.bioclim <- tree.sample.regimes.expanded <- lapply(tree.sample, function(tree) {
  tree <- paintSubTree(
    tree,
    getMRCA(tree, c("Daucus_carota", "Actinotus_helianthi", "Platysace_lanceolata")),
    state = 1,
    anc.state = 0
  )
  tree <- paintSubTree(
    tree,
    getMRCA(tree, c("Daucus_carota", "Ammi_majus", "Chamaesium_paradoxum")),
    state = 2,
    anc.state = 0
  )
  return(tree)
})

tree.sample.regimes.expanded <- lapply(tree.sample, function(tree) {
  tree <- paintSubTree(
    tree,
    getMRCA(tree, c("Daucus_carota", "Actinotus_helianthi", "Platysace_lanceolata", "Hymenosporum_flavum", "Myodocarpus_fraxinifolius")),
    state = 1,
    anc.state = 0
  )
  tree <- paintSubTree(
    tree,
    getMRCA(tree, c("Daucus_carota", "Ammi_majus", "Chamaesium_paradoxum")),
    state = 2,
    anc.state = 0
  )
  tree <- paintSubTree(
    tree,
    getMRCA(tree, c("Harmsiopanax_ingens", "Dendropanax_arboreus", "Fatsia_japonica", "Heptapleurum_heptaphyllum", "Pseudopanax_crassifolius")),
    state = 3,
    anc.state = 0
  )
  return(tree)
})

#Hypothesis 1: Climate with increased desiccation intensity favors larger pollen grains
morpho.H1 <- as.data.frame(morpho[row.names(morpho)%in%row.names(PCs),])[12]
bioclim.H1 <- as.data.frame(PCs)
H1.data <- cbind.data.frame(morpho.H1, bioclim.H1)

tree.sample.H1 <- keep.tip.multiPhylo(tree.sample, row.names(PCs))

H1.fit.BM <- mclapply(indices, function(i) {
  phylolm(VOL~PC1+PC2, 
          data = H1.data,
          model = "BM",
          tree.sample.H1[[i]])
}, mc.cores = 20)

saveRDS(H1.fit.BM, "RDS/H1.fit.BM.RDS")

H1.fit.OU <- mclapply(indices, function(i) {
  phylolm(VOL~PC1+PC2, 
        data = H1.data,
        model = "OUrandomRoot",
        tree.sample.H1[[i]])
}, mc.cores = 20)

saveRDS(H1.fit.OU, "RDS/H1.fit.OU.RDS")

#Hypothesis 2: Transition from more tropical to more temperate climate favors reduction of apertures and thickening 
#of pollen wall to counteract increased desiccation intensity

#Fit models
H2.BM1 <- mclapply(indices, function(i) {
  mvBM(tree.sample.regimes.bioclim[[i]], 
       morpho [,c(3,4,8,10,11)], 
       model = "BM1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 16)

saveRDS(H2.BM1, "RDS/H2.BM1.RDS")

H2.BMM <- mclapply(indices, function(i) {
  mvBM(tree.sample.regimes.bioclim[[i]], 
       morpho [,c(3,4,8,10,11)], 
       model = "BMM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 16)

saveRDS(H2.BMM, "RDS/H2.BMM.RDS")

#H2.BMM.expanded <- mclapply(indices, function(i) {
#  mvBM(tree.sample.regimes.expanded[[5]], 
#       morpho [,c(3,4,8,10,11)], 
#       model = "BMM", 
#       method = "sparse",
#       scale.height = TRUE,
#       control = list(maxit = 300000),
#       param = list(root = FALSE))
#}, mc.cores = 16)

saveRDS(H2.BMM.expanded, "RDS/H2.BMM.expanded.RDS")

H2.OU1 <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.bioclim[[i]], 
       morpho [,c(3,4,8,10,11)], 
       model = "OU1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 16)

saveRDS(H2.OU1, "RDS/H2.OU1.RDS")

H2.OUM <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.bioclim[[i]], 
       morpho [,c(3,4,8,10,11)], 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 16)

saveRDS(H2.OUM, "RDS/H2.OUM.RDS")

H2.OUM.expanded <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.expanded[[5]], 
       morpho [,c(3,4,8,10,11)], 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 16)

saveRDS(H2.OUM.expanded, "RDS/H2.OUM.expanded.RDS")

AICc.compare(H2.BM1)
AICc.compare(H2.BMM)
#AICc.compare(H2.BMM.expanded)
AICc.compare(H2.OU1)
AICc.compare(H2.OUM)
AICc.compare(H2.OUM.expanded)

#Hypothesis 3: Pollen shape is a good predictor of aperture and nexine morphology
#as well as wall thickness

PGLS.data <- list(
  shape = as.matrix (morpho[,13:14]), 
  nexine = as.matrix (morpho [,c(5,6,7,9)]),
  aperture = as.matrix (morpho [,c(8,10,11)]),
  wall = as.matrix (morpho[,3:4]))

H3.fit.BM <- mclapply(indices, function(i) {
  mvgls(shape ~ aperture + wall + nexine, 
        data = PGLS.data, 
        tree.sample[[i]], 
        model="BM",
        method = "LL", 
        penalty="RidgeArch")
}, mc.cores = 16)

saveRDS(H3.fit.BM, "RDS/H3.fit.BM.RDS")

H3.fit.OU <- mclapply(indices, function(i) {
  mvgls(shape ~ aperture + wall + nexine, 
        data = PGLS.data, 
        tree.sample[[i]], 
        model="OU",
        method = "LL", 
        penalty="RidgeArch")
}, mc.cores = 16)

saveRDS(H3.fit.OU, "RDS/H3.fit.OU.RDS")

H3.MANCOVA.BM <- mclapply(indices, function(i) {
  manova.gls(H3.fit.BM[[i]], test = "Wilks", nperm = 9999, type = "II")
}, mc.cores = 16)

saveRDS(H3.MANCOVA.BM, "RDS/H3.MANCOVA.BM.RDS")

H3.MANCOVA.OU <- mclapply(indices, function(i) {
  manova.gls(H3.fit.OU[[i]], test = "Wilks", nperm = 9999, type = "II")
}, mc.cores = 16)

saveRDS(H3.MANCOVA.OU, "RDS/H3.MANCOVA.OU.RDS")

#Hypothesis 4: Pollen shape lags behind aperture length/Aperture length lags behind pollen shape
morpho.lag <- morpho [,c(8,10,11,3,4)]
morpho.lag2 <- morpho [,c(3,4,8,10,11)]

alpha_matrix_lag <- matrix(0, nrow=5, ncol=5)
alpha_matrix_lag[1,] <- c (1,6,7,NA,NA) 
alpha_matrix_lag[2,] <- c (6,2,8,NA,NA)
alpha_matrix_lag[3,] <- c (7,8,3,NA,NA)
alpha_matrix_lag[4,] <- c (9,10,11,4,15) 
alpha_matrix_lag[5,] <- c (12,13,14,15,5) 

alpha_matrix_lag2 <- matrix(0, nrow=5, ncol=5)
alpha_matrix_lag2[1,] <- c (1,6,NA,NA,NA) 
alpha_matrix_lag2[2,] <- c (6,2,NA,NA,NA)
alpha_matrix_lag2[3,] <- c (7,8,3,11,14)
alpha_matrix_lag2[4,] <- c (9,10,11,4,15) 
alpha_matrix_lag2[5,] <- c (12,13,14,15,5) 

H4.OUM.lag <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.bioclim[[i]], 
       morpho.lag, 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE, decomp = alpha_matrix_lag))
}, mc.cores = 50)

saveRDS(H4.OUM.lag, "RDS/H4.OUM.lag.RDS")

H4.OUM.lag2 <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.bioclim[[i]], 
       morpho.lag2, 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE, decomp = alpha_matrix_lag2))
}, mc.cores = 50)

saveRDS(H4.OUM.lag2, "RDS/H4.OUM.lag2.RDS")

H4.OUM.B.lag <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.expanded[[i]], 
       morpho.lag, 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE, decomp = alpha_matrix_lag))
}, mc.cores = 50)

saveRDS(H4.OUM.lag, "RDS/H4.OUM.B.lag.RDS")

H4.OUM.B.lag2 <- mclapply(indices, function(i) {
  mvOU(tree.sample.regimes.expanded[[i]], 
       morpho.lag2, 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE, decomp = alpha_matrix_lag2))
}, mc.cores = 50)

saveRDS(H4.OUM.B.lag2, "RDS/H4.OUM.B.lag2.RDS")

##############################################
##############################################
################## Plots #####################
##############################################
##############################################

###########################################################################
#### Hypothesis 1 - you will see white cells, since CI are around 0 :( ####
###########################################################################

# Step 1: Extract AICc values from both BM and OU models
H1.AICc.summary <- data.frame(
  Tree = indices,
  BM = sapply(H1.fit.BM, function(x) x$aic),
  OU = sapply(H1.fit.OU, function(x) x$aic)
)

# Step 2: Melt to long format
H1.AICc.summary.melt <- H1.AICc.summary %>%
  pivot_longer(cols = -Tree, names_to = "Model", values_to = "AICc")

# Step 3: Compute model order and delta (difference)
H1.AICc.summary.melt <- H1.AICc.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc) %>%
  mutate(
    ModelOrder = max(row_number()) + 1 - row_number(),
    Difference = case_when(
      row_number() == 1 ~ AICc,
      row_number() == 2 ~ AICc - lag(AICc),
      row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1)
    )
  )

# Step 4: Plot raw AICc per tree with model color and group-based bar stacking
H1.AICc.plot <- ggplot(H1.AICc.summary.melt, aes(x = factor(Tree), y = AICc, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"  # consistent with your H3.EIC plot
  ) +
  scale_fill_manual(values = c("firebrick4", "springgreen4")) +
  scale_x_discrete(breaks = seq(1, 100, 10)) +
  coord_cartesian(ylim = c(min(H1.AICc.summary.melt$AICc), NA)) +
  labs(
    title = element_blank(),
    x = "Tree",
    y = "AICc",
    fill = "Model"
  )


# Helper function to generate a smooth CI gradient
expand_ci <- function(lower, upper, n = 100) {
  seq(lower, upper, length.out = n)
}

# Step 1: Extract coefficients and compute 95% CI across 100 trees
results_list <- lapply(H1.fit.OU, function(fit) {
  summary_table <- summary(fit)[[2]]
  data.frame(
    Predictor = rownames(summary_table),
    Estimate = summary_table[, 1],
    PValue = summary_table[, 4]
  )
})

combined_results <- bind_rows(results_list, .id = "Tree") %>%
  filter(Predictor != "(Intercept)")

# Step 2: Compute CI for each predictor
ci_results <- combined_results %>%
  group_by(Predictor) %>%
  summarize(
    LowerCI = quantile(Estimate, probs = 0.025, na.rm = TRUE),
    UpperCI = quantile(Estimate, probs = 0.975, na.rm = TRUE)
  )

# Step 3: Expand CI to gradient rows for each predictor
expanded_df <- ci_results %>%
  rowwise() %>%
  mutate(
    CI_Values = list(expand_ci(LowerCI, UpperCI, n = 100)),
    Position = list(seq(-0.5, 0.5, length.out = 100))
  ) %>%
  unnest(c(CI_Values, Position)) %>%
  ungroup()

# Step 4: Heatmap with CI gradient
H1.heatmap <- ggplot(expanded_df, aes(x = Predictor, y = "VOL")) +
  geom_tile(aes(fill = CI_Values, height = 1, width = 1), color = "white") +
  geom_raster(aes(x = as.numeric(as.factor(Predictor)) + Position, fill = CI_Values)) +
  geom_text(data = ci_results, aes(y = "VOL", label = paste0("[", round(LowerCI, 2), ", ", round(UpperCI, 2), "]")),
            color = "black", size = 4, vjust = 1.5) +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick4",
                       midpoint = 0, limits = c(-1, 1), name = "95% CI Range") +
  theme_minimal() +
  labs(
    x = "Predictor",
    y = "Response"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(),
        strip.text.y = element_text(size = 10))

# Display the heatmap
print(H1.heatmap)

######################
#### Hypothesis 2 ####
######################

H2.summary <- data.frame(Tree = indices,
                         BM1 = AICc.vec(H2.BM1),
                         BMM = AICc.vec(H2.BMM),
                         OU1 = AICc.vec(H2.OU1),
                         OUM.e = AICc.vec(H2.OUM.expanded),
                         OUM = AICc.vec(H2.OUM),
                         OUM.lag = AICc.vec(H4.OUM.lag),
                         OUM.lag.e = AICc.vec(H4.OUM.B.lag),
                         OUM.lag2 = AICc.vec(H4.OUM.lag2),
                         OUM.lag2.e = AICc.vec(H4.OUM.B.lag2))

H2.summary.melt <- H2.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

H2.summary.melt <- H2.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1),
           row_number() == 4 ~ AICc - lag(AICc, order_by = AICc, n = 2),
           row_number() == 5 ~ AICc - lag(AICc, order_by = AICc, n = 3),
           row_number() == 6 ~ AICc - lag(AICc, order_by = AICc, n = 3),
           row_number() == 7 ~ AICc - lag(AICc, order_by = AICc, n = 3),
           row_number() == 8 ~ AICc - lag(AICc, order_by = AICc, n = 3),
           row_number() == 9 ~ AICc - lag(AICc, order_by = AICc, n = 3)
           
         )
  )

H2.summary.plot <- ggplot(H2.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  scale_fill_manual(values = c("azure3",
                              "cornflowerblue", 
                              "darkseagreen3", 
                              "goldenrod2", 
                              "indianred", 
                              "lightslateblue", 
                              "orchid4", 
                              "purple4",
                              "lightsteelblue4")) +
  theme(legend.position="bottom") +
  # scale_y_reverse() +
  labs(title = "",
       x = "Tree", y = "AICc", fill = "Model") +
  scale_x_discrete(breaks = seq(from = 1, to = 100, by = 10)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  coord_cartesian(ylim = c(min(H2.summary.melt$AICc), NA))

ggsave("plots/H2_AICc.pdf", H2.summary.plot, device = "pdf", width = 160, height = 100, units = "mm")

######################
#### Hypothesis 3 ####
######################

#EIC

H3.fit.BM.EIC <- mclapply(H3.fit.BM, EIC, mc.cores = 16, mc.preschedule = FALSE)

H3.fit.OU.EIC <- mclapply(H3.fit.OU, EIC, mc.cores = 16, mc.preschedule = FALSE)

H3.EIC.summary <- data.frame(Tree = indices,
                             BM = EICvec(H3.fit.BM.EIC),
                             OU = EICvec(H3.fit.OU.EIC))

H3.EIC.summary.melt <- H3.EIC.summary %>%
  gather(key = "Model", value = "H3.EIC", -Tree) 

H3.EIC.summary.melt  <- H3.EIC.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, H3.EIC)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ H3.EIC,
           row_number() == 2 ~ H3.EIC - lag(H3.EIC),
           row_number() == 3 ~ H3.EIC - lag(H3.EIC, order_by = H3.EIC, n = 1)
         )
  )

H3.EIC.summary.plot <- ggplot(H3.EIC.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_discrete(breaks = seq(from = 1, to = 100, by = 10)) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick4", "springgreen4")) +
  coord_cartesian(ylim = c(min(H3.EIC.summary.melt$H3.EIC), NA)) +
  # scale_y_reverse() +
  labs(title = element_blank(),
       x = "Tree", y = "EIC", fill = "Model")

ggsave("/scratch/jjb67768/HELIOPHILA_FIN_OF_FINS_KAMIL_NEW/plots/H3.EIC.summary.plot.pdf", 
       H3.EIC.summary.plot, 
       device = "pdf", 
       width = 8, 
       height = 10, 
       units = "in")

#MANOVA

expand_ci <- function(lower, upper, n = 100) {
  seq(lower, upper, length.out = n)
}

# Step 1: Extract coefficients and calculate CIs across 100 trees
results_list <- map2(H3.fit.OU, H3.MANCOVA.OU, function(fit, manova) {
  coefficients <- as.data.frame(coef(fit))
  coefficients <- coefficients[rownames(coefficients) != "(Intercept)", ]
  
  coefficients$Trait <- rownames(coefficients)
  
  coefficients <- coefficients %>%
    pivot_longer(-Trait, names_to = "Predictor", values_to = "Coefficient")
  
  coefficients <- coefficients %>%
    mutate(
      Matrix = case_when(
        grepl("aperture", Trait) ~ "aperture",
        grepl("nexine", Trait) ~ "nexine",
        grepl("wall", Trait) ~ "wall"
      ),
      PValue = case_when(
        Matrix == "aperture" ~ manova$pvalue[1],
        Matrix == "nexine" ~ manova$pvalue[2],
        Matrix == "wall" ~ manova$pvalue[3]
      ),
      Trait = sub("^[a-z]+", "", Trait)
    )
  
  return(coefficients)
})

# Combine results and calculate CI across all trees
results_df <- bind_rows(results_list) %>%
  group_by(Trait, Predictor, Matrix) %>%
  summarize(
    LowerCI = quantile(Coefficient, probs = 0.025, na.rm = TRUE),
    UpperCI = quantile(Coefficient, probs = 0.975, na.rm = TRUE)
  )

# Step 2: Expand each CI into a sequence of values for gradient effect
expanded_df <- results_df %>%
  rowwise() %>%
  mutate(CI_Values = list(expand_ci(LowerCI, UpperCI, n = 100))) %>%
  unnest(CI_Values) %>%
  mutate(Position = seq(-0.5, 0.5, length.out = 100)) %>%
  ungroup()

# Step 3: Create the heatmap with CI gradient and add CI range text
H3.heatmap <- ggplot(expanded_df, aes(x = Predictor, y = Trait)) +
  geom_tile(aes(fill = CI_Values, height = 1, width = 1), color = "white") +
  geom_raster(aes(x = as.numeric(as.factor(Predictor)) + Position, fill = CI_Values)) +
  geom_text(data = results_df, aes(label = paste0("[", round(LowerCI, 2), ", ", round(UpperCI, 2), "]")), 
            color = "black", size = 5, vjust = 1.5) +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick4",
                      limits = c(-1,1), midpoint = 0, name ="95% CI Range") +
  facet_grid(Matrix ~ ., scales = "free", space = "free", switch = "y") +
  theme_minimal() +
  labs(
    x = "Response",
    y = "Predictor"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(size = 10))


ggsave("plots/H3.MANOVA.pdf", 
       H3.heatmap, 
       device = "pdf", 
       width = 8, 
       height = 10, 
       units = "in")

#P-value density plots

# Step 1: Extract p-values from each tree's MANOVA result

p_values_list <- map(H3.MANCOVA.OU, function(manova) {
  data.frame(
    aperture = manova$pvalue[1],
    wall = manova$pvalue[2],
    nexine = manova$pvalue[3]
  )
})

# Combine p-values from all trees into a single data frame

p_values_df <- bind_rows(p_values_list) %>%
  pivot_longer(cols = everything(), names_to = "Matrix", values_to = "PValue")

# Step 2: Plot the distribution of p-values

H3.pvalues <- ggplot(p_values_df, aes(x = PValue)) +
  geom_density() +
  facet_wrap(~ Matrix, scales = "free", ncol = 1) +
  theme_minimal() +
  labs(
    title = "",
    x = "P-Value",
    y = "Density"
  ) +
  theme(strip.text = element_text(size = 10))

ggsave("plots/H3.prob.pdf", 
       H3.pvalues, 
       device = "pdf", 
       width = 2.5, 
       height = 4, 
       units = "in")

######################
#### Hypothesis 4 ####
######################

H4.summary <- data.frame(Tree = indices,
                         H4.nolag = AICc.vec(H2.OUM),
                         H4.lag = AICc.vec(H4.OUM.lag),
                         H4.lag2 = AICc.vec(H4.OUM.lag2))

H4.summary.melt <- H4.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

H4.summary.melt <- H4.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1)
         )
  )

H4.summary.plot <- ggplot(H4.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  scale_fill_manual(values = c("dodgerblue1", "indianred", "springgreen2")) +
  theme(legend.position="bottom") +
  # scale_y_reverse() +
  labs(title = "",
       x = "Tree", y = "AICc", fill = "Model") +
  scale_x_discrete(breaks = seq(from = 1, to = 100, by = 10)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  coord_cartesian(ylim = c(min(H4.summary.melt$AICc), NA))

ggsave("plots//H4_AICc.pdf", H4.summary.plot, device = "pdf", width = 160, height = 50, units = "mm")

###################################################
###################################################
################## Halflife H2 #####################
###################################################
###################################################

H2.hl <- cbind.data.frame (meand = rowMeans(sapply(H4.OUM.lag2, halflife)),
                           sd = apply(sapply(H2.OUM, halflife), 1, sd))

rownames(H2.hl) <- colnames(H2.OUM[[4]]$theta)

###########################################################
###########################################################
################## Evolutionary optima ####################
###########################################################
###########################################################

#Revert scaling

polmean <- colMeans(morpho.raw[,c(7,8,12,14,15)])
polscale <- sqrt(diag(var(morpho.raw[,c(7,8,12,14,15)])))

H2.theta <- theta.compare(H4.OUM.lag2, polscale, polmean)

H2.theta.bind <-  theta_converter(H2.theta, "theta", "climate")

H2.theta.bind <- H2.theta.bind %>%
  mutate(column = factor(climate), 
         panel_id = factor(panel_id),
         list_id = factor(list_id))

H2.theta.bind$row_name <- gsub("_", " ", H2.theta.bind$row_name)
H2.theta.bind$row_name <- paste(H2.theta.bind$row_name)
H2.theta.bind$climate <- gsub("0", "climatic regime - 0", H2.theta.bind$climate)
H2.theta.bind$climate <- gsub("1", "climatic regime - 1", H2.theta.bind$climate)
H2.theta.bind$climate <- gsub("2", "climatic regime - 2", H2.theta.bind$climate)

H2.theta.plot <- ggplot(H2.theta.bind, aes(x = value, fill = interaction(list_id, climate))) + 
  geom_density(alpha = 0.6) + # Add some transparency
  facet_wrap(~row_name, scales = "free", ncol = 1) + # Create a panel for each row, adjust ncol as needed
  labs(title = element_blank(), 
       x = "", 
       y = "Density") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(values = c("grey40" ,"goldenrod1", "purple3")) +
  theme(legend.position="none") +
  labs(fill = "")

H2.theta.meansd <- H2.theta.bind %>%
  group_by(row_name, climate) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

ggsave("plots/H2_theta.pdf", 
       H2.theta.plot, 
       device = "pdf", 
       width = 160, 
       height = 140, 
       units = "mm")