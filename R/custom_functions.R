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
library(edgeR)
library(tidyr)
library(readr)
library(inflection)

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

GeneDetector <- function(SeuratObj, threshold) { 
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
  return(genes_list)
}

contrastConstructor <- 
  function(clusterIDs,
           clusterIDs.2 = TRUE,
           design,
           makeContrast = FALSE,
           name = NULL){
    clusterIDs.char <- as.character(clusterIDs)
    clusterIDs.char.celltype <- paste("celltype", clusterIDs.char, sep = "")
    clusterIDs.char.celltype.combined <- paste(clusterIDs.char.celltype, collapse = "+")
    clusterIDs.length <- as.character(length(clusterIDs.char.celltype))
    clusterIDs.char.celltype.combined.parantheses <- paste0("(", clusterIDs.char.celltype.combined, ")")
    clusterIDs.ready <- paste(clusterIDs.char.celltype.combined.parantheses, 
                              clusterIDs.length, 
                              sep = "/", collapse = "")
    if(missing(clusterIDs.2))
    {
      otherClusterIDs <- setdiff(x = colnames(design), 
                                 y = clusterIDs.char.celltype)
      otherClusterIDs.num <- grep("celltype[[:digit:]]", otherClusterIDs, value = T)
      otherClusterIDs.length <- as.character(length(otherClusterIDs.num))
      otherClusterIDs.combined <- paste(otherClusterIDs.num, collapse = "+")
      otherClusterIDs.combined.parentheses <- paste0("(", otherClusterIDs.combined, ")")
      otherClusterIDs.ready <- paste(otherClusterIDs.combined.parentheses, otherClusterIDs.length, sep = "/", collapse = "")
      contrastString <- paste(clusterIDs.ready, otherClusterIDs.ready, sep = "-")
    }else{
      clusterIDs.2.char <- as.character(clusterIDs.2)
      clusterIDs.2.char.celltype <- paste("celltype", clusterIDs.2.char, sep = "")
      clusterIDs.2.char.celltype.combined <- paste(clusterIDs.2.char.celltype, collapse = "+")
      clusterIDs.2.length <- as.character(length(clusterIDs.2.char.celltype))
      clusterIDs.2.char.celltype.combined.parantheses <- paste0("(", clusterIDs.2.char.celltype.combined, ")")
      clusterIDs.2.ready <- paste(clusterIDs.2.char.celltype.combined.parantheses, 
                                  clusterIDs.2.length, 
                                  sep = "/", collapse = "")
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

edgeR_create_fit <- function(object, gene.detection.thres = 0.10, celltype.ident = "order")
  {
## Set the object identity 
object <- SetIdent(object, value = celltype.ident)
  
## Detect genes that are expressed in at least 10% of any given cluster in the Seurat object and write them out into a vector called "marker.names".
marker.names <- GeneDetector(SeuratObj = object, threshold = gene.detection.thres)

## Filter the raw gene-cell counts matrix(stored in the Seurat object) with the genes identified above
tmp.matrix <- object@assays$RNA@counts[marker.names, ]

## Create the Seurat object to recalculate the nFeature_RNA (unique genes detected/cell)
tmp <- CreateSeuratObject(tmp.matrix)
colnames(tmp@meta.data)[3] <- "Recalculated_nGene"
tmp <- AddMetaData(object = tmp, metadata = object@meta.data)

## Calculate the gdr from the nGene column in the Seurat metadat
tmp@meta.data$gdr <- as.numeric(scale(tmp@meta.data$Recalculated_nGene))

## Create DGE counts object from the edgeR package
counts <- DGEList(counts = tmp@assays$RNA@counts, genes = rownames(tmp))
## Calculate the library sizes of each single cell to be analyzed
counts$samples$lib.PBMSCize <- colSums(counts$counts)
## Calculate the normalization factors to be applied to make single cell comparable
counts <- edgeR::calcNormFactors(counts)

## Reformat metadata from the Seurat object to make compatible with edgeR/design matrix
edgeRmetadata <- tmp@meta.data[, c("Subject", celltype.ident, "gdr")]
colnames(edgeRmetadata)[2] <- "celltype"
edgeRmetadata$celltype <- as.factor(edgeRmetadata$celltype) %>% droplevels()
design <- model.matrix(~ 0 + gdr + celltype + Subject, data = edgeRmetadata)
colnames(design) <- make.names(colnames(design))

## Estimate dispersion
counts <- estimateDisp(counts, design, robust = T)
## Perform the edgeR fit
fit <- glmQLFit(counts, design)

return_list <- list(fit, counts, design)
names(return_list) <- c("fit", "counts", "design")
return(return_list)
}

add.percent.expression <- 
  function(dge_list, object, ident, subset = NULL, group.ident = NULL){
    if(!is.null(subset)){
      object <- SetIdent(object, value = group.ident)
      object <- subset(x = object, ident = subset)
    }
    dge_list <- dge_list[gtools::mixedorder(names(dge_list))]
    my_genes <- lapply(dge_list, tibble::remove_rownames) %>% 
      lapply(tibble::column_to_rownames, var = "genes") %>% 
      lapply(rownames)
    object <- SetIdent(object, value = ident)
    exp <- SplitObject(object, split.by = ident)
    exp <- exp[gtools::mixedorder(names(exp))]
    mat <- 
      mapply(function(x,y) {
        FetchData(object = x, vars = y)
      }, 
      exp, my_genes)
    matrix <- lapply(mat, function(x) as.matrix(colMeans(x  > 0))*100)
    matrix <- lapply(matrix, as.data.frame)
    matrix <- lapply(matrix, setNames, "percent expressed")
    dge_list.add <- mapply(cbind, dge_list, matrix, SIMPLIFY=FALSE)
    return(dge_list.add)
  }

left <- function (string,char) {
  base::substr(string,1,char)
}

right <- function (string, char) {
  base::substr(string,nchar(string)-(char-1),nchar(string))
}

cross_species_dendrogram <- function(species.1.obj, species.2.obj, 
                                     species1.ident, species2.ident,
                                     species1.res, species2.res,
                                     species1.name, species2.name)
{
  ## Set objects to compare here (seuratObject)
  obj_1 <- subset(species.1.obj, idents = species1.ident)
  obj_2 <- subset(species.2.obj, idents = species2.ident)
  
  ### Extract there metadata for cluster ids (seuratObject@meta.data)
  obj_1.meta <- obj_1@meta.data
  obj_2.meta <- obj_2@meta.data
  ## Set resolution to perform comparisons (STRING)
  obj_1.res <- species1.res
  obj_2.res <- species2.res
  ### Set string for species parameter for id on dendrogram (STRING)
  obj_1.species <- species1.name
  obj_2.species <- species2.name
  
  # Extract counts matrix
  obj_1.cMat <- obj_1@assays$RNA@counts
  obj_2.cMat <- obj_2@assays$RNA@counts
  
  obj1.Sobj <- CreateSeuratObject(obj_1.cMat)
  obj2.Sobj <- CreateSeuratObject(obj_2.cMat)
  
  # Find highly variable genes under the SCTransform framework
  obj1.Sobj <- SCTransform(obj1.Sobj)
  obj2.Sobj <- SCTransform(obj2.Sobj)
  
  # Extact highly variable genes
  obj_1.genes <- VariableFeatures(obj1.Sobj)
  obj_2.genes <- VariableFeatures(obj2.Sobj)
  
  obj1.varplot <- VariableFeaturePlot(obj1.Sobj, log = T)$data %>% arrange(residual_variance)
  obj2.varplot <- VariableFeaturePlot(obj2.Sobj, log = T)$data %>% arrange(residual_variance)
  
  library(inflection)
  obj1.inf <- ese(x = 1:length(obj1.varplot[,2]), 
                  y = obj1.varplot[,2], 
                  index = 1)
  obj1.inf.pt <- round(length(obj1.varplot[,2]) - (length(obj1.varplot[,2]) - obj1.inf[,3])/2)
  obj_1.genes <- rownames(obj1.varplot)[obj1.inf.pt:length(rownames(obj1.varplot))]
  
  obj2.inf <- ese(x = 1:length(obj2.varplot[,2]), 
                  y = obj2.varplot[,2] %>% sort(), 
                  index = 1)
  obj2.inf.pt <- round(length(obj2.varplot[,2]) - (length(obj2.varplot[,2]) - obj2.inf[,3])/2)
  obj_2.genes <- rownames(obj2.varplot)[obj2.inf.pt:length(rownames(obj2.varplot))]
  
  
  # Intersect species hvg's to identify shared highly variable orthologs
  orthologs <- intersect(obj_1.genes,obj_2.genes)
  
  # Use RelativeCounts function within Seurat to calculate TPM from counts matrix
  obj_1.TPM <- RelativeCounts(obj_1.cMat, scale.factor = 10^6)
  obj_2.TPM <- RelativeCounts(obj_2.cMat, scale.factor = 10^6)
  
  # Add some cluster meta data to cell names in the TPM matrices; will be used to identify whose 'who'
  obj_1 <- SetIdent(obj_1, value = obj_1.res)
  obj_2 <- SetIdent(obj_2, value = obj_2.res)
  
  colnames(obj_1.TPM) <- paste(colnames(obj_1.TPM), Idents(obj_1), sep = "_")
  colnames(obj_2.TPM) <- paste(colnames(obj_2.TPM), Idents(obj_2), sep = "_")
  
  # Filter matrices by highly var orthologs to make file sizes more maneagable
  obj_1.TPM <- obj_1.TPM[orthologs, ]
  obj_2.TPM <- obj_2.TPM[orthologs, ]
  
  # Convert from sparseMatrix to matrix format to work well with base R functions
  obj_1.TPM <- as.matrix(obj_1.TPM)
  obj_2.TPM <- as.matrix(obj_2.TPM)
  
  # Tidy data and calculate gene-cluster averages across all groups
  library(reshape2)
  obj_1.melt <- melt(obj_1.TPM)
  colnames(obj_1.melt) <- c("gene", "cell", "count")
  obj_1.melt$cluster <- sub(".*[0-9]_", "", obj_1.melt$cell)
  
  obj_2.melt <- melt(obj_2.TPM)
  colnames(obj_2.melt) <- c("gene", "cell", "count")
  obj_2.melt$cluster <- sub(".*_", "", obj_2.melt$cell)
  
  library(dplyr)
  obj_1.groupedbyClusterandGene_means <- obj_1.melt %>%
    group_by(cluster, gene) %>%
    summarise(mean = mean(count))
  obj_2.groupedbyClusterandGene_means <- obj_2.melt %>%
    group_by(cluster, gene) %>%
    summarise(mean = mean(count))
  
  obj_1.groupedbyClusterandGene_means.df <- dcast(data = obj_1.groupedbyClusterandGene_means, formula = gene ~ cluster, fun.aggregate = sum, value.var = "mean")
  rownames(obj_1.groupedbyClusterandGene_means.df) <- obj_1.groupedbyClusterandGene_means.df[, 1]
  obj_1.groupedbyClusterandGene_means.df[, 1] <- NULL
  
  obj_2.groupedbyClusterandGene_means.df <- dcast(data = obj_2.groupedbyClusterandGene_means, formula = gene ~ cluster, fun.aggregate = sum, value.var = "mean")
  rownames(obj_2.groupedbyClusterandGene_means.df) <- obj_2.groupedbyClusterandGene_means.df[, 1]
  obj_2.groupedbyClusterandGene_means.df[, 1] <- NULL
  
  # Add psuedocount TPM of 50 to each matrix of mean expression
  obj_1.groupedbyClusterandGene_means.df.p <- obj_1.groupedbyClusterandGene_means.df + 50
  obj_2.groupedbyClusterandGene_means.df.p <- obj_2.groupedbyClusterandGene_means.df + 50
  
  # Calculate medians of average expression per gene
  obj_1.median.mean <- apply(obj_1.groupedbyClusterandGene_means.df.p, MARGIN = 1, median)
  obj_2.median.mean <- apply(obj_2.groupedbyClusterandGene_means.df.p, MARGIN = 1, median)
  
  # Add psuedocount TPM of 50 to each median
  obj_1.median.mean.p <- obj_1.median.mean + 50
  obj_2.median.mean.p <- obj_2.median.mean + 50
  
  # Perform klein transformation
  obj_1.kleinTrans <- sweep(obj_1.groupedbyClusterandGene_means.df.p, 1, obj_1.median.mean, "/")
  obj_2.kleinTrans <- sweep(obj_2.groupedbyClusterandGene_means.df.p, 1, obj_2.median.mean, "/")
  
  ## Append species identifier to column names
  colnames(obj_1.kleinTrans) <- paste(obj_1.species, colnames(obj_1.kleinTrans), sep = "_")
  colnames(obj_2.kleinTrans) <- paste(obj_2.species, colnames(obj_2.kleinTrans), sep = "_")
  
  # Combine matrices via bind columns(species clusters)
  merged.obj <- cbind(obj_1.kleinTrans, obj_2.kleinTrans)
  
  # Perform log transformation(natural)
  merged.obj.log1p <- log1p(merged.obj)
  
  # Calculate pearson distance metrices
  pearson.dist.df <- as.dist((1 - cor(merged.obj.log1p)) / 2)
  
  # Run hclust and visualize dendogram
  dendo <- hclust(pearson.dist.df, method = "ward.D2")
  return(dendo)
}

