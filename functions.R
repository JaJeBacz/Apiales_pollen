#Function returning mean and sd from AICc scores

AICc.compare <- function (data) {
  AICc <- list()
  for (i in seq_along(data)){
    AICc[[i]] <- data[[i]]$AICc
  }
  return (c (mean(unlist(AICc)), sd(unlist(AICc))))
} 

AICc.vec <- function (data) {
  AICcvec <- list()
  for (i in seq_along(data)){
    AICcvec[[i]] <- data[[i]]$AICc
  }
  return (unlist(AICcvec))
} 

GICvec <- function (data) {
  GIC <- list()
  for (i in seq_along(data)){
    GIC[[i]] <- data[[i]]$GIC
  }
  return (unlist(GIC))
} 

EICvec <- function (data) {
  EIC <- list()
  for (i in seq_along(data)){
    EIC[[i]] <- data[[i]][[1]]
  }
  return (unlist(EIC))
} 

pvalvec <- function (data) {
  pval <- list()
  for (i in seq_along(data)){
    pval[[i]] <- data[[i]]$pvalue[[1]]
  }
  return (unlist(pval))
} 

theta.compare <- function (data, scale, center) {
  theta <- list()
  for (i in seq_along(data)){
    theta[[i]] <- apply (data[[i]]$theta, 1, function(x) x * scale + center)
  }
  return (theta)
} 

theta.compare.BIO <- function (data) {
  theta <- list()
  for (i in seq_along(data)){
    theta[[i]] <- t(data[[i]]$theta)
  }
  return (theta)
} 

theta_converter <- function(df_list, list_id, names) {
  bind_rows(lapply(seq_along(df_list), function(i) {
    as.data.frame(df_list[[i]]) %>%
      mutate(row_name = rownames(.),
             df_id = i,
             list_id = list_id) %>%
      pivot_longer(cols = c('0', '1', '2'), names_to = names, values_to = "value")
  }), .id = "panel_id")
}

#Function calculating mean Q matrix

Q.avg <- function (corHMM) {
  Q <- list()
  for (i in seq_along(corHMM)){
    Q[[i]] <- corHMM[[i]]$solution
  }
  Q.dcb <- do.call (cbind, Q)
  Q.array <- array (Q.dcb, 
                    dim = c(dim(Q[[1]]), length(Q)))
  Q.mean <- apply (Q.array, 
                   c(1,2), 
                   mean, 
                   na.rm = TRUE)
  return(Q)
  return(Q.mean)
} 

#Modifiet corHMM plot function

plotRECONJB <- function (phy, likelihoods, tipcolors, piecolors = NULL, cex = 0.5, pie.cex = 0.25, 
                         file = NULL, height = 11, width = 8.5, show.tip.label = TRUE, 
                         title = NULL, ...) 
{
  if (is.null(piecolors)) {
    piecolors = c("white", "black", "red", "yellow", "forestgreen", 
                  "blue", "coral", "aquamarine", "darkorchid", "gold", 
                  "grey", "yellow", "#3288BD", "#E31A1C")
  }
  if (!is.null(file)) {
    pdf(file, height = height, width = width, useDingbats = FALSE)
  }
  plot(ladderize(phy), cex = cex, show.tip.label = show.tip.label, tip.color = tipcolors, ...)
  if (!is.null(title)) {
    title(main = title)
  }
  nodelabels(pie = likelihoods, piecol = piecolors, cex = pie.cex)
  states <- colnames(likelihoods)
  if (!is.null(file)) {
    dev.off()
  }
}

#Modified phylomorphospace3d plotting function

phylomorphospace3dJB <- function (tree, X, A = NULL, label = TRUE, control = list(), 
                                  method = c("dynamic", "static"), point.colors = NULL, 
                                  sphere.size = 0.02, label.size = 1, label.colors = NULL, ...) 
{
  method <- method[1]
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  con = list(spin = TRUE, axes = TRUE, box = TRUE, simple.axes = FALSE, 
             lwd = 1, ftype = "reg", col.edge = rep("black", nrow(tree$edge)))
  con[(namc <- names(control))] <- control
  if (con$simple.axes) 
    con$box <- con$axes <- FALSE
  con$ftype <- which(c("off", "reg", "b", "i", "bi") == con$ftype) - 1
  if (is.null(A)) 
    A <- apply(X, 2, function(x, tree) fastAnc(tree, x), tree = tree)
  else A <- A[as.character(1:tree$Nnode + length(tree$tip)), ]
  x <- y <- z <- matrix(NA, nrow(tree$edge), 2)
  X <- X[tree$tip.label, ]
  for (i in 1:length(tree$tip)) {
    x[tree$edge[, 2] == i, 2] <- X[i, 1]
    y[tree$edge[, 2] == i, 2] <- X[i, 2]
    z[tree$edge[, 2] == i, 2] <- X[i, 3]
  }
  for (i in length(tree$tip) + 1:tree$Nnode) {
    x[tree$edge[, 1] == i, 1] <- x[tree$edge[, 2] == i, 2] <- A[as.character(i), 1]
    y[tree$edge[, 1] == i, 1] <- y[tree$edge[, 2] == i, 2] <- A[as.character(i), 2]
    z[tree$edge[, 1] == i, 1] <- z[tree$edge[, 2] == i, 2] <- A[as.character(i), 3]
  }
  if (is.null(colnames(X))) 
    colnames(X) <- c("x", "y", "z")
  if (method == "dynamic") {
    chk <- require("rgl")
    if (!chk) {
      cat("  method = \"dynamic\" requires the package \"rgl\"\n  Defaulting to method = \"static\"\n\n")
      method <- "static"
      lines3d <- play3d <- plot3d <- segments3d <- points3d <- spin3d <- text3d <- function(...) NULL
    }
  }
  if (method == "dynamic") {
    params <- get("r3dDefaults")
    plot3d(rbind(X, A), xlab = colnames(X)[1], ylab = colnames(X)[2], 
           zlab = colnames(X)[3], axes = con$axes, box = con$box, 
           params = params)
    for (i in 1:nrow(tree$edge)) segments3d(x[i, ], y[i, ], z[i, ], lwd = con$lwd, col = con$col.edge[i])
    ms <- colMeans(X)
    rs <- apply(rbind(X, A), 2, range)
    if (con$simple.axes) {
      lines3d(x = rs[, 1], y = c(rs[1, 2], rs[1, 2]), z = c(rs[1, 3], rs[1, 3]))
      lines3d(x = c(rs[1, 1], rs[1, 1]), y = rs[, 2], z = c(rs[1, 3], rs[1, 3]))
      lines3d(x = c(rs[1, 1], rs[1, 1]), y = c(rs[1, 2], rs[1, 2]), z = rs[, 3])
    }
    rs <- rs[2, ] - rs[1, ]
    if (is.null(point.colors)) point.colors <- rep("black", nrow(X))
    points3d(X[, 1], X[, 2], X[, 3], color = point.colors, size = sphere.size * 100)
    for (i in 1:length(tree$tip)) {
      adj <- 0.03 * rs * (2 * (X[i, ] > ms) - 1)
      if (con$ftype && !is.null(label.colors)) 
        text3d(X[i, ] + adj, texts = tree$tip.label[i], font = con$ftype, cex = label.size, color = label.colors[i])
      else if (con$ftype)
        text3d(X[i, ] + adj, texts = tree$tip.label[i], font = con$ftype, cex = label.size)
    }
    if (con$spin) {
      xx <- spin3d(axis = c(0, 0, 1), rpm = 10)
      play3d(xx, duration = 5)
      invisible(xx)
    }
    else invisible(NULL)
  }
  else if (method == "static") {
    if (hasArg(angle)) 
      angle <- list(...)$angle
    else angle <- 30
    if (hasArg(xlim)) 
      xlim <- list(...)$xlim
    else xlim <- NULL
    if (hasArg(ylim)) 
      ylim <- list(...)$ylim
    else ylim <- NULL
    if (hasArg(zlim)) 
      zlim = list(...)$zlim
    else zlim = NULL
    point_colors <- if (is.null(point.colors)) rep("black", nrow(X)) else point.colors
    xx <- scatterplot3d(X, xlab = colnames(X)[1], zlab = colnames(X)[3], 
                        pch = 19, angle = angle, ylab = colnames(X)[2], cex.symbols = 1.3, 
                        xlim = xlim, ylim = ylim, zlim = zlim, color = point_colors)
    aa <- xx$xyz.convert(A)
    points(aa$x, aa$y, pch = 19, cex = 0.8)
    for (i in 1:nrow(tree$edge)) {
      aa <- xx$xyz.convert(x[i, ], y[i, ], z[i, ])
      lines(aa$x, aa$y, lwd = 2, col = con$col.edge[i])
    }
    for (i in 1:length(tree$tip.label)) {
      aa <- xx$xyz.convert(x[which(tree$edge[, 2] == i), 2], y[which(tree$edge[, 2] == i), 2], z[which(tree$edge[, 2] == i), 2])
      if (con$ftype && !is.null(label.colors)) 
        text(tree$tip.label[i], x = aa$x, y = aa$y, pos = 2, font = con$ftype, col = label.colors[i])
      else if (con$ftype)
        text(tree$tip.label[i], x = aa$x, y = aa$y, pos = 2, font = con$ftype)
    }
    invisible(xx)
  }
}

#Summarizing model fitting

summarize_model_fit <- function(data) {
  best_models <- data %>% filter(ModelOrder == 9)
  second_best_models <- data %>% filter(ModelOrder == 8)
  
  combined <- best_models %>% 
    inner_join(second_best_models, by = "Tree", suffix = c("_best", "_second_best"))
  
  combined <- combined %>%
    mutate(Difference = AICc_second_best - AICc_best)
  
  best_fit_counts <- combined %>%
    mutate(FitType = ifelse(Difference > 2, "BestFit", "Inconclusive")) %>%
    group_by(Model_best, FitType) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = FitType, values_from = Count, values_fill = list(Count = 0))
  
  summary_table <- best_fit_counts %>%
    mutate(Total = BestFit + Inconclusive)

}

#Wrap splitter

save_facet_plots <- function(plot, facet_var, data, path_prefix) {
  facets <- unique(data[[facet_var]])
  
  for (facet in facets) {
    p_facet <- plot +
      ggforce::facet_wrap_paginate(~row_name, scales = "free", ncol = 1, page = which(facets == facet)) +
      labs(title = paste("Hypothesis 4:", facet))
    
    ggsave(filename = paste0(path_prefix, "_", facet, ".tiff"), plot = p_facet, device = "tiff")
  }
}

#Phylomorphospace 2D

ggmorphoJB <- function (tree, tipinfo, xvar = PC1, yvar = PC2, factorvar = group, 
                        labelvar = taxon, title = "Phylomorphospace", xlab = "PC1", 
                        ylab = "PC2", repel = TRUE, edge.width = 1, fontface = "italic", 
                        tree.alpha = 0.7, point.size = 5, point.shapes = NULL, show.labels = TRUE) 
{
  require(ggplot2)
  require(phytools)
  require(ggrepel)
  
  # Default shapes if none are provided
  if (is.null(point.shapes)) {
    point.shapes <- c(16, 17, 18, 10, 15, 8, 1, 2, 0, 3, 4, 5, 6, 7)
  }
  
  mat <- cbind(eval(substitute(xvar), tipinfo), eval(substitute(yvar), 
                                                     tipinfo))
  rownames(mat) <- eval(substitute(labelvar), tipinfo)
  stopifnot(length(setdiff(tree$tip.label, rownames(mat))) == 
              0)
  
  xAnc <- fastAnc(tree, mat[, 1])
  yAnc <- fastAnc(tree, mat[, 2])
  
  all_node_coords <- data.frame(x = c(mat[tree$tip.label, 1], 
                                      xAnc), y = c(mat[tree$tip.label, 2], yAnc), nodeid = 1:(tree$Nnode + 
                                                                                                length(tree$tip.label)))
  
  edges <- data.frame(tree$edge)
  names(edges) <- c("node1", "node2")
  
  edgecoords <- merge(merge(edges, all_node_coords, by.x = "node1", 
                            by.y = "nodeid"), all_node_coords, by.x = "node2", by.y = "nodeid")
  
  pointsForPlot <- data.frame(x = eval(substitute(xvar), tipinfo), 
                              y = eval(substitute(yvar), tipinfo), shape = eval(substitute(factorvar), 
                                                                                tipinfo), label = eval(substitute(labelvar), tipinfo))
  
  theplot <- ggplot() + geom_segment(data = edgecoords, aes(x = x.x, 
                                                            xend = x.y, y = y.x, yend = y.y), linewidth = edge.width, 
                                     alpha = tree.alpha) + geom_point(data = pointsForPlot, 
                                                                      aes(x = x, y = y, shape = shape), size = point.size) + 
    labs(title = title, x = xlab, y = ylab) + theme_bw(20) + 
    theme(legend.position = "bottom")
  
  theplot <- theplot + scale_shape_manual(values = point.shapes)
  
  if (show.labels) {
    if (repel) {
      theplot <- theplot + geom_text_repel(data = pointsForPlot, 
                                           aes(x = x, y = y, label = label), segment.alpha = 0.5, 
                                           fontface = fontface)
    } else {
      theplot <- theplot + geom_text(data = pointsForPlot, 
                                     aes(x = x, y = y, label = label), fontface = fontface)
    }
  }
  
  return(theplot)
}

#Plotting solids of revolution

plot_solid_of_revolution <- function(outline, n = 100) {
  
  # Extract and center the coordinates
  coords <- coo_center(outline)$coo[[1]]
  x <- coords[, 1]
  y <- coords[, 2]
  
  # Find the midpoint of the x-range
  x_mid <- (max(x) + min(x)) / 2
  
  # Filter to keep only the lower half of the outline (y <= y_mid)
  keep_indices <- y <= 0  # Rotate around y = 0, so we keep the lower half
  x_half <- x[keep_indices]
  y_half <- y[keep_indices]
  
  # Sort by x for consistency
  sorted_indices <- order(x_half)
  x_half <- x_half[sorted_indices]
  y_half <- y_half[sorted_indices]
  
  # Generate points for revolution around the x-axis
  theta <- seq(0, 2*pi, length.out = n)
  Y <- outer(y_half, cos(theta))  # Rotate y coordinates around x-axis
  Z <- outer(y_half, sin(theta))  # Rotate y coordinates around x-axis
  X <- outer(x_half, rep(1, n))   # X stays the same, this is the axis of rotation
  
  # Plot using plotly without legend and axis numbers
  plot_ly(x = X, y = Y, z = Z, type = 'surface', colorscale = list(c(0, 1), c("black", "#ffcc00ff")), showscale = FALSE) %>%
    layout(scene = list(
      xaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE),
      yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE),
      zaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE),
      showlegend = FALSE
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0)
    )
}

#Save plots

save_plotly_as_png_webshot <- function(plot, filename) {
  temp_html <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(as_widget(plot), temp_html, selfcontained = TRUE)
  webshot(temp_html, file = filename)
}

calculate_volume_of_revolution_x_axis_rescaled <- function(outline, length, n = 100) {
  # Extract and center the coordinates
  coords <- coo_center(outline)
  x <- coords[, 1]
  y <- coords[, 2]
  
  # Rescale the coordinates based on the length of the pollen grain
  original_length <- max(x) - min(x)
  scale_factor <- length / original_length
  x_rescaled <- x * scale_factor
  y_rescaled <- y * scale_factor
  
  # Filter to keep only the lower half of the outline (y <= 0)
  keep_indices <- y_rescaled <= 0  # Rotate around y = 0, so we keep the lower half
  x_half <- x_rescaled[keep_indices]
  y_half <- y_rescaled[keep_indices]
  
  # Sort by x for consistency
  sorted_indices <- order(x_half)
  x_half <- x_half[sorted_indices]
  y_half <- y_half[sorted_indices]
  
  # Calculate the volume using the disk method
  volume <- 2 * pi * trapz(x_half, y_half^2)
  
  # Convert volume into cubic micrometers (this is now in real-world units)
  return(volume)
}


phylomorphospace3dJB_plot3D <- function(tree, X, A = NULL, control = list(), 
                                        shapes = NULL, clade_colors = NULL, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  
  # Set the edge color to black for the phylogenetic tree
  con = list(col.edge = rep("grey70", nrow(tree$edge)))
  con[(namc <- names(control))] <- control
  
  # Check and clean the input data X
  X_subset <- X[, c("PC1eq", "VOL", "PC1pol")]
  if (anyNA(X_subset) || any(!is.finite(as.matrix(X_subset)))) {
    stop("The input data X contains NA, NaN, or non-numeric values. Please clean your data before proceeding.")
  }
  
  if (is.null(A)) 
    A <- apply(X_subset, 2, function(x, tree) fastAnc(tree, x), tree = tree)
  else 
    A <- A[as.character(1:tree$Nnode + length(tree$tip)), ]
  
  x <- y <- z <- matrix(NA, nrow(tree$edge), 2)
  X_subset <- X_subset[tree$tip.label, ]
  for (i in 1:length(tree$tip)) {
    x[tree$edge[, 2] == i, 2] <- X_subset[i, 1]
    y[tree$edge[, 2] == i, 2] <- X_subset[i, 2]
    z[tree$edge[, 2] == i, 2] <- X_subset[i, 3]
  }
  for (i in length(tree$tip) + 1:tree$Nnode) {
    x[tree$edge[, 1] == i, 1] <- x[tree$edge[, 2] == i, 2] <- A[as.character(i), 1]
    y[tree$edge[, 1] == i, 1] <- y[tree$edge[, 2] == i, 2] <- A[as.character(i), 2]
    z[tree$edge[, 1] == i, 1] <- z[tree$edge[, 2] == i, 2] <- A[as.character(i), 3]
  }
  
  # Generate colors based on clade
  if (is.null(clade_colors)) {
    stop("Please provide a vector of colors corresponding to each clade.")
  }
  point_colors <- clade_colors[match(X$clade, names(clade_colors))]
  
  # Directly assign shapes based on family
  shapes_vector <- shapes[match(X$family, names(shapes))]
  
  # Get the minimum z value for the vertical lines
  z_min <- min(X_subset[, 3])
  
  # Plot the points with custom shapes and colors
  scatter3D(x = X_subset[, 1], y = X_subset[, 2], z = X_subset[, 3], 
            colvar = NULL, col = point_colors, pch = shapes_vector, 
            cex = 0.5, bty = "f", ticktype = "detailed", main = "",  
            colkey = FALSE, box = TRUE, axes = TRUE, 
            d = 2, theta = 30, phi = 30)
  
  # Add vertical lines to the bottom plane (type = "h")
  for (i in 1:nrow(X_subset)) {
    segments3D(x0 = X_subset[i, 1], y0 = X_subset[i, 2], z0 = z_min, 
               x1 = X_subset[i, 1], y1 = X_subset[i, 2], z1 = X_subset[i, 3], 
               col = "gray90", lwd = 0.5, add = TRUE)  # Gray vertical lines
  }
  
  # Plot the phylogenetic tree edges in black
  for (i in 1:nrow(tree$edge)) {
    segments3D(x0 = x[i, 1], y0 = y[i, 1], z0 = z[i, 1], 
               x1 = x[i, 2], y1 = y[i, 2], z1 = z[i, 2], 
               col = con$col.edge[i], lwd = 2, add = TRUE)
  }
}





ggmorphoJB_color <- function (tree, tipinfo, xvar = PC1, yvar = PC2, factorvar = group, 
                              labelvar = taxon, title = "Phylomorphospace", xlab = "PC1", 
                              ylab = "PC2", repel = TRUE, edge.width = 1, fontface = "italic", 
                              tree.alpha = 0.7, point.size = 5, show.labels = TRUE) 
{
  require(ggplot2)
  require(phytools)
  require(ggrepel)
  
  # Define color palette for clades
  clade_colors <- c(
    "Azorelloideae" = "#1f77b4",
    "Apioideae" = "#ff7f0e",
    "Araliaceae" = "#2ca02c",
    "Griseliniaceae" = "#d62728",
    "Myodocarpaceae" = "#9467bd",
    "Pennantiaceae" = "#8c564b",
    "Pittosporaceae" = "#e377c2",
    "Torricelliaceae" = "#7f7f7f",
    "Hermas" = "#bcbd22",
    "Klotzschia" = "#17becf",
    "Mackinlayoideae" = "#ff9896",
    "Platysace" = "#aec7e8",
    "Saniculoideae" = "#c5b0d5"
  )
  
  mat <- cbind(eval(substitute(xvar), tipinfo), eval(substitute(yvar), tipinfo))
  rownames(mat) <- eval(substitute(labelvar), tipinfo)
  stopifnot(length(setdiff(tree$tip.label, rownames(mat))) == 0)
  
  xAnc <- fastAnc(tree, mat[, 1])
  yAnc <- fastAnc(tree, mat[, 2])
  
  all_node_coords <- data.frame(x = c(mat[tree$tip.label, 1], xAnc),
                                y = c(mat[tree$tip.label, 2], yAnc),
                                nodeid = 1:(tree$Nnode + length(tree$tip.label)))
  
  edges <- data.frame(tree$edge)
  names(edges) <- c("node1", "node2")
  
  edgecoords <- merge(merge(edges, all_node_coords, by.x = "node1", 
                            by.y = "nodeid"), all_node_coords, by.x = "node2", 
                      by.y = "nodeid")
  
  pointsForPlot <- data.frame(
    x = eval(substitute(xvar), tipinfo),
    y = eval(substitute(yvar), tipinfo),
    clade = eval(substitute(factorvar), tipinfo),
    label = eval(substitute(labelvar), tipinfo)
  )
  
  theplot <- ggplot() +
    geom_segment(data = edgecoords, aes(x = x.x, xend = x.y, y = y.x, yend = y.y), 
                 linewidth = edge.width, alpha = tree.alpha) +
    geom_point(data = pointsForPlot, aes(x = x, y = y, color = clade), size = point.size) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw(20) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = clade_colors)
  
  if (show.labels) {
    if (repel) {
      theplot <- theplot + geom_text_repel(data = pointsForPlot, 
                                           aes(x = x, y = y, label = label),
                                           segment.alpha = 0.5, fontface = fontface)
    } else {
      theplot <- theplot + geom_text(data = pointsForPlot, 
                                     aes(x = x, y = y, label = label), fontface = fontface)
    }
  }
  
  return(theplot)
}
