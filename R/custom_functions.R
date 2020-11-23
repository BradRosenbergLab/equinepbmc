## This script contains custom functions used for analysis and visualizations for the equine PBMC manuscript

## Load required functions to run this function. Install prior; packages are a medley of CRAN, github and Bioconductor packages
library(rlang)
library(Seurat)
library(ggnetwork)
library(RColorBrewer)
library(ggrepel)
library(readxl)
library(dplyr)
library(ggpubr)
library(ggforce)
library(lemon)
library(patchwork)
library(ape)
library(jntools)
library(BiocManager)
library(viridis)
library(ggtree)
library(ggforce)
library(legocolors)
library(clustree)
library(readxl)
library(openxlsx)
library(cluster)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(ggtree)

#### Hack Seurat functions to plot geom_point (shape) and geom_text for label clusters, underlayed with a circles. This is key, as it identifies where the labels are plotted on the UMAP
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}

## Custom function to add circles labels
custom.LabelClusters <- function(
  plot, # Use DimPlot to generate base ggplot to apply function
  id,   # Th
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = F,
  shape = 21,
  colors = colors,
  circle.size = circle.size,
  text.size = text.size,
  text.color = text.color,
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape=shape, stroke = 0.66, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  ) + geom_text( 
    size = text.size,
    color = text.color,
    data = labels.loc, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  )
  return(plot)
}

## Custom function to add circles labels
custom.LabelClusters <- function(
  plot, # Use DimPlot to generate base ggplot to apply function
  id,   # Th
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = F,
  colors = colors,
  circle.size = circle.size,
  text.size = text.size,
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape=21, stroke = 0.66, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  ) + geom_text( 
    size = text.size,
    data = labels.loc, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  )
  return(plot)
}

cluster.umap.plot <- function(object, 
                              group.ident = "integrated_snn_res.1",
                              colors = legocolors$hex[-1],
                              pt.size = 0.01,
                              legend.col = 4,
                              ciclelabels = TRUE,
                              byrow = TRUE,
                              label = FALSE,
                              inset = TRUE,
                              shape = 21,
                              legend.title = "Cluster",
                              legend.position = "top left",
                              legend.text.align = 0){
  cellnumber <- dim(object)[2]
  p2 <- 
    DimPlot(object = object, 
            pt.size = pt.size,
            label = label,
            group.by = group.ident) +
    scale_color_manual(values = colors)  + 
    scale_fill_manual(values = colors)  + 
    theme.c +
    xlab(expression('UMAP'[1])) +
    ylab(expression('UMAP'[2])) + 
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(hjust = 0, vjust = 0),
          axis.title.y = element_text(hjust = 0, vjust = 0),
          legend.title = element_text(face = "bold"),
          legend.text.align = legend.text.align,
          legend.text = element_text(color = "black",  inherit.blank = T,
                                     margin = margin(t = 0, r = 0, b = 0, l = -5)),
          legend.key.width = unit(0, "mm"),
          legend.key = element_rect(fill = NA, inherit.blank = T),
          legend.background = element_rect(fill = NA, color = "black", size = 0.25),        
          legend.margin = margin(0,0,0,0)) + 
    labs(color = legend.title) +
    ggtitle(paste("Unsupervised clustering", " (", scales::comma_format(digits = 12)(cellnumber), " cells)", sep = "")) +
    guides(colour = guide_legend(override.aes = list(shape = 21, 
                                                     size = 6,
                                                     color = "black",
                                                     fill = colors[1:length(levels(object@meta.data[,group.ident]))]), 
                                 label = T, 
                                 byrow = byrow,
                                 ncol = legend.col)) 
  if (ciclelabels) {
  p2 <- custom.LabelClusters(plot = p2, 
                             id = group.ident, 
                             clusters = levels(object@meta.data[,group.ident]),
                             circle.size = 7, 
                             text.color = "black",
                             text.size = 4, 
                             shape = shape,
                             fill = colors[1:length(levels(object@meta.data[,group.ident]))],
                             repel = T)
  }
  if (inset) {
    p2.gtable <- reposition_legend(aplot = p2, 
                                   offset = c(0.001,0.001),
                                   position = legend.position) 
    
    p2 <- wrap_ggplot_grob(p2.gtable)
  }
  print(p2)
}

FeaturePlot.c <- function(object, features){FeaturePlot(object = object, 
                                                        features = features, 
                                                        min.cutoff = "q2", 
                                                        max.cutoff = "q98",
                                                        cols = c("grey90",
                                                                 "darkgreen"), 
                                                        order = T) + NoAxes() }

theme.c <- theme(axis.line = element_blank(),
                 aspect.ratio = 1,
                 panel.background = element_blank(),
                 panel.border = element_rect(color = "black", fill = NA, size = 0.25))

Workbook <- function(x, file.name){
  wb <- createWorkbook()
  ## Loop through the list of split tables as well as their names
  ##   and add each one as a sheet to the workbook
  Map(function(data, name){
    
    addWorksheet(wb, name)
    writeData(wb, name, data, keepNA = TRUE)
    
  }, x, names(x))
  ## Save workbook to working directory
  saveWorkbook(wb, file = file.name, overwrite = TRUE)
}

library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

GeneDetector <- function(SeuratObj, threshold, listname) { 
  counts <- SeuratObj
  genes <- c()
  for (i in 1:(length(unique(SeuratObj@active.ident))))
  {
    cluster <- subset(SeuratObj, idents = levels(SeuratObj@active.ident)[i])
    cluster_df <- cluster@assays$RNA@counts
    cluster_df <- as.data.frame(cluster_df)
    print(sum(rowMeans(cluster_df > 0) >= threshold))
    genes_df <- as.data.frame(rowMeans(cluster_df > 0))
    filtered_df <- subset(x = genes_df, subset = rowMeans(cluster_df > 0) > threshold)
    genes <- c(genes, rownames(filtered_df))
  } 
  genes_list <- unique(genes)
  assign(x = listname, value = genes_list, envir = globalenv())
}

contrastConstructor <- 
  function(clusterIDs,
           clusterIDs.2 = TRUE,
           design,
           makeContrast = FALSE,
           name = NULL){
    clusterIDs
    clusterIDs.char <- as.character(clusterIDs)
    clusterIDs.char
    clusterIDs.char.celltype <- paste("celltype", clusterIDs.char, sep = "")
    clusterIDs.char.celltype
    clusterIDs.char.celltype.combined <- paste(clusterIDs.char.celltype, collapse = "+")
    clusterIDs.char.celltype.combined
    clusterIDs.length <- as.character(length(clusterIDs.char.celltype))
    clusterIDs.length
    clusterIDs.char.celltype.combined.parantheses <- paste0("(", clusterIDs.char.celltype.combined, ")")
    clusterIDs.ready <- paste(clusterIDs.char.celltype.combined.parantheses, 
                              clusterIDs.length, 
                              sep = "/", collapse = "")
    clusterIDs.ready
    
    if(missing(clusterIDs.2))
    {
      otherClusterIDs <- setdiff(x = colnames(design), 
                                 y = clusterIDs.char.celltype)
      otherClusterIDs.num <- grep("[[:digit:]]", otherClusterIDs, value = T)
      otherClusterIDs.length <- as.character(length(otherClusterIDs.num))
      otherClusterIDs.length
      otherClusterIDs.combined <- paste(otherClusterIDs.num, collapse = "+")
      otherClusterIDs.combined
      otherClusterIDs.combined.parentheses <- paste0("(", otherClusterIDs.combined, ")")
      otherClusterIDs.combined.parentheses
      otherClusterIDs.ready <- paste(otherClusterIDs.combined.parentheses, otherClusterIDs.length, sep = "/", collapse = "")
      otherClusterIDs.ready
      contrastString <- paste(clusterIDs.ready, otherClusterIDs.ready, sep = "-")
      contrastString
    }else{
      clusterIDs.2
      clusterIDs.2.char <- as.character(clusterIDs.2)
      clusterIDs.2.char
      clusterIDs.2.char.celltype <- paste("celltype", clusterIDs.2.char, sep = "")
      clusterIDs.2.char.celltype
      clusterIDs.2.char.celltype.combined <- paste(clusterIDs.2.char.celltype, collapse = "+")
      clusterIDs.2.char.celltype.combined
      clusterIDs.2.length <- as.character(length(clusterIDs.2.char.celltype))
      clusterIDs.2.length
      clusterIDs.2.char.celltype.combined.parantheses <- paste0("(", clusterIDs.2.char.celltype.combined, ")")
      clusterIDs.2.ready <- paste(clusterIDs.2.char.celltype.combined.parantheses, 
                                  clusterIDs.2.length, 
                                  sep = "/", collapse = "")
      clusterIDs.2.ready
      contrastString <- paste(clusterIDs.ready, clusterIDs.2.ready, sep = "-")
      
    }
    if(makeContrast == T){
      if(missing(name)){
        contrast <- makeContrasts(contrasts = contrastString, levels = design)
      }
      else{
        contrast <- makeContrasts(contrasts = contrastString, levels = design)
        assign(x = name, value = contrast, envir = globalenv())  
      }
    }else{
      return(contrastString)
    }
  }

generateContrasts <- function(cluster.numbers, prefix){
  for(i in 1:length(cluster.numbers)){
    contrastConstructor(
      clusterIDs = cluster.numbers[i], 
      clusterIDs.2 = setdiff(cluster.numbers, cluster.numbers[i]),
      design = design.,
      makeContrast = T, 
      name = paste(prefix, cluster.numbers[i], sep = ""))
  }
}


plottree <- function (object, ...) 
{
  library(ape)
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  plot.phylo(x = data.tree, direction = "downwards", ...)
}

prettyheatmap <- function(object, 
                          gene.list,
                          text.size = 8,
                          flagged_genes, 
                          downsample = 30, 
                          colors, 
                          ident, 
                          order){
  {
    ## Set/order of identities to plot (right to left)
    levels(object) <- order
    ## To create custom heatmaps, extract the assay slot you wish to visualize. Generally, this slot is the natural log normalized counts slot in the Seurat object.
    object <- ScaleData(object, scale.max = 50)
    mat <- object@assays$RNA@scale.data
    ## Filter the extracted matrix against the genes you are interested in visualizing. The rownames are gene names.
    cellnames <- WhichCells(object, downsample = downsample)
    mat <- mat[gene.list, cellnames]
    object.sub <- subset(object, cells = cellnames)
    
    ## Next, extract relevant meta.data of interest
    # It is important that meta.data stays 'connected' to cell barcodes/names
    cluster.data <- object.sub@meta.data[, c("orig.ident", ident)]
    # Give meta.data a name; will be the legend in the pheatmap
    colnames(cluster.data) <- c("source", "cluster")
    # I did this because I kept losing the rownames when I applied the sort. If you know a better way to sort AND retain the rownames let me know.
    # I think its a problem when its a dataframe with only 1 column
    cluster.data$name <- rownames(cluster.data)
    # Sort the cell barcodes, according to a specified order. I am sorting the clusters(descending) by numerical order. You can also explicitly specify this.
    cluster.data.ordered <- cluster.data[order(cluster.data$cluster), ]
    cluster.data.ordered$name <- NULL
    cluster.data.ordered$source <- NULL
    cluster.data.ordered$cluster2 <-
      factor(cluster.data.ordered$cluster, levels = order %>% as.character())
    cluster.data.ordered <-
      cluster.data.ordered[order(cluster.data.ordered$cluster2), ]
    cluster.data.ordered$cluster <- cluster.data.ordered$cluster2
    cluster.data.ordered$cluster2 <- NULL
    ## Order the extracted matrix according to the meta.data
    mat.ordered <-
      mat[, match(x = rownames(cluster.data.ordered), table = cellnames)]
    complex.matrix.rev <- as.matrix(mat.ordered)
    
    ## Extract top n number of non-ENSECAG genes to highlight/flag on the heatmap - unsupervised
    annoGenes.highlight <- flagged_genes
    anno_coords.highlight <-
      match(x = annoGenes.highlight, table = rownames(complex.matrix.rev))
    ## Create annotation object for gene flags using anno_mark
    ha = rowAnnotation(
      foo = anno_mark(
        at = anno_coords.highlight,
        labels = annoGenes.highlight,
        padding = 0.3,
        labels_gp = gpar(fontsize = text.size, fontface = "italic"),
        side = "left",
        link_gp = gpar(lwd = 0.2)
      )
    )
    lha = rowAnnotation(gene_expression_sum = anno_lines(complex.matrix.rev))
    
    names(colors) <- order
    ta = columnAnnotation(
      cluster = cluster.data.ordered$cluster,
      col = list(cluster = colors[na.omit(names(colors))]),
      show_legend = T,
      show_annotation_name = FALSE
    )
    
    ## Extract RColorBrewer 'RdYlBu' color pallete and ramp accordingly to use for expression scale
    pal_RdYlBu <- rev(brewer.pal(n = 11, name = "RdYlBu"))
    hm_breaks <- seq(from =  -2,
                     to = 2,
                     length.out = 11)
    col_fun = colorRamp2(breaks = hm_breaks,
                         space = "LAB",
                         colors = pal_RdYlBu)
    ## Set hieght and width of plots within complex heatmap call, dependent of figure dimensions specified
    map_height <- nrow(complex.matrix.rev) * 0.1
    map_width <- length(unique(cluster.data.ordered$cluster)) * 0.4
    
    heatmap <- Heatmap(
      complex.matrix.rev,
      cluster_rows = T,
      cluster_columns = FALSE, show_row_dend = F,
      #  width = unit(3.1033, "in"), height = unit(3.8867, "in"),
      col = col_fun,
      color_space = 'LAB',
      left_annotation = ha,
      top_annotation = ta,
      show_row_names = T,
      show_column_names = F,
      show_heatmap_legend = T,
      row_names_gp = gpar(fontsize = 6),
      #        column_names_rot = T, column_names_centered = T,
      name = "Scaled\nexpression",
      border = T
    )
    draw(heatmap)
  }
}



prettyDots <-
  function(dotplot,
           colors,
           col.min = -2.5,
           col.max = 2.5) {
    library(cowplot)
    dotplot$data$Scaled.Expression <- dotplot$data$avg.exp.scaled
    dotplot$data$Percent.Expressed <- dotplot$data$pct.exp
    
    plot <- ggplot(dotplot$data,
                   aes(x = factor(id), y = features.plot)) +
      geom_point(aes(color = Scaled.Expression, size = Percent.Expressed)) +
      scale_size(range = c(0.5, 8), limits = c(0, 100)) +
      geom_tile(aes(fill = id, y = 0), size = 1) +
      scale_fill_manual(values = colors) +
      scale_color_gradientn(
        limits = c(col.min, col.max),
        breaks = scales::pretty_breaks()(col.min:col.max),
        colours = (brewer.pal(name = "YlOrRd",  n = 9))
      ) +
      theme_bw() +
      theme(
        legend.title = element_text(size = 12),
        legend.box = "horizontal",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size = 1)
      ) +
      geom_text(aes(y = 0, label = id)) +
      geom_point(aes(y = 0, size = 75), shape = 21, stroke = 0.5) + coord_fixed() 
    
    plot.noleg <- ggplot(dotplot$data,
                         aes(x = id, y = features.plot)) +
      geom_point(aes(color = Scaled.Expression, size = Percent.Expressed)) +
      scale_size(range = c(0.5, 8), limits = c(0, 100)) +
      geom_tile(aes(fill = id, y = 0), size = 1) +
      scale_fill_manual(values = colors) +
      scale_color_gradientn(
        limits = c(col.min, col.max),
        breaks = scales::pretty_breaks()(col.min:col.max),
        colours = (brewer.pal(name = "YlOrRd",  n = 9))
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size = 1)
      ) +
      geom_text(aes(y = 0, label = id)) +
      geom_point(aes(y = 0, size = 75), shape = 21, stroke = 0.5) + coord_fixed() 
    print(plot)
    print(plot.noleg)
  }
