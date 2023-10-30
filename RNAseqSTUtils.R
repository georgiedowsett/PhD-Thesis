#functions 

calc_percent_cluster = function(seurat_object, genes, metadata_name) {
  gene_expr <- FetchData(seurat_object, vars = genes)
  clusters <- seurat_object@meta.data[,metadata_name]
  avex_percluster <- AverageExpression(seurat_object, features = genes, group.by = metadata_name)
  avex_percluster<- data.frame(t(avex_percluster$RNA))
  if (any(grepl("^ENSMUSG", colnames(avex_percluster)))) {
    avex_percluster <- convert_colnames(avex_percluster)
  }
  colnames(avex_percluster) <- paste0("AvEx_", colnames(avex_percluster))
  avex_percluster <- round(avex_percluster, 2)
  gene_expr[gene_expr > 0] <- 1
  gene_expr[gene_expr < 1] <- 0
  total_cells_per_cluster <- table(clusters)
  total_cells_per_cluster <- subset(total_cells_per_cluster, total_cells_per_cluster > 0)
  cells_expressing_gene_per_cluster <- rowsum(gene_expr, clusters)
  percent_expressing_gene_per_cluster <- cells_expressing_gene_per_cluster / total_cells_per_cluster * 100
  percent_expressing_gene_per_cluster <- round(percent_expressing_gene_per_cluster, 2)
  if (any(grepl("^ENSMUSG", colnames(percent_expressing_gene_per_cluster)))) {
    percent_expressing_gene_per_cluster <- convert_colnames(percent_expressing_gene_per_cluster)
  }
  result<-cbind(percent_expressing_gene_per_cluster, avex_percluster)
  return(as.data.frame(result))
}

convert_colnames <- function(df) {
  new_colnames <- sapply(colnames(df), function(ensid) {
    gene_symbol <- gene_names1[gene_names1$Gene.stable.ID == ensid, 'Gene.name']
    if (length(gene_symbol) == 0) {
      return(ensid)  # Return original ENSID if no matching gene symbol
    } else {
      return(as.character(gene_symbol))
    }
  })
  colnames(df) <- new_colnames
  return(df)
}

top_gene_clusters <- function(data, subset, freq_cut, percent_cut, metadata_name = metadata_name) {
  num_cells_per_cluster <- table(data@meta.data[,metadata_name])
  gene_cells_per_cluster <- table(subset@meta.data[,metadata_name])
  gene_clusters <- merge(num_cells_per_cluster, gene_cells_per_cluster, by = "Var1")
  gene_clusters$percent <- gene_clusters$Freq.y / gene_clusters$Freq.x * 100
  gene_clusters_subset <- gene_clusters[gene_clusters$Freq.y > freq_cut & gene_clusters$percent > percent_cut,]
  return(gene_clusters_subset)
}

reformat_table<-function(cluster_table, num_to_name, data, gene){
  x<-merge(cluster_table, num_to_name, by.x = 'Var1', by.y = 'cluster',  all.x = T, all.y = F)
  y<-as.table(t(as.data.frame(AverageExpression(subset(data, idents = as.character(cluster_table$Var1)), features = gene)$RNA)))
  x<-merge(x, y, by= 'Var1')
  names(x)[names(x) == "Freq"] <- "AverageExpression"
  x<-x[,-which(names(x) == 'Var2')]
  names(x)[names(x) == "Freq.x"] <- "Number_Cells"
  names(x)[names(x) == "Freq.y"] <- paste0("Num_",gene,"+Cells")
  names(x)[names(x) == "percent"] <- paste0("Percent_",gene)
  names(x)[names(x) == "Var1"] <- "Cluster"
  names(x)[names(x) == "genes_conc"] <- "ClusterName"
  return(x)
}
add_marker_genes <- function(gene_clusters, clusters, markers){
  for (i in clusters) {
    x <- subset(markers, subset = cluster == i)
    x <- x[!grepl("^LINC|^AC\\d+|^AL\\d+|^MIR\\d+|RIK$", x$gene), ]
    x <- x[order(-x$specificity), ]
    x <- x$gene[1:4]
    row_index <- which(clusters == i)
    x <- paste(x, collapse = "|")
    gene_clusters[row_index, "gene_markers"] <- x
  }
  return(gene_clusters)
}

my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)
huhy_plots<-function(data, gene, percentile = NA, image.alpha = 0, assay.use = data@active.assay, nrow = 1){
  if (is.na(percentile)){
    gene_expression<-data@assays[[assay.use]]@data[gene,]
    a<-SpatialFeaturePlot(data, gene, image = 'slice5A', image.alpha = image.alpha) 
    b<-SpatialFeaturePlot(data, gene, image = 'slice6B', image.alpha = image.alpha) 
    c<-SpatialFeaturePlot(data, gene, image = 'slice2B', image.alpha = image.alpha) 
    d<-SpatialFeaturePlot(data, gene, image = 'slice4B', image.alpha = image.alpha) 
    e<-SpatialFeaturePlot(data, gene, image = 'slice7A', image.alpha = image.alpha) 
    f<-SpatialFeaturePlot(data, gene, image = 'slice3A', image.alpha = image.alpha) 
    g<-SpatialFeaturePlot(data, gene, image = 'slice8B', image.alpha = image.alpha) 
    a+b+c+d+e+f+g & plot_layout(nrow = nrow) & scale_fill_gradientn(colors = my_spatial_palette, oob=squish)
  }
  else {
    gene_expression <- data@assays[[assay.use]]@data[gene,]
    scale <- quantile(gene_expression, percentile)
    scale <- as.numeric(round(scale, 1))
    a<-SpatialFeaturePlot(data, gene, image = 'slice5A', image.alpha = image.alpha)  
    b<-SpatialFeaturePlot(data, gene, image = 'slice6B', image.alpha = image.alpha) 
    c<-SpatialFeaturePlot(data, gene, image = 'slice2B', image.alpha = image.alpha) 
    d<-SpatialFeaturePlot(data, gene, image = 'slice4B', image.alpha = image.alpha) 
    e<-SpatialFeaturePlot(data, gene, image = 'slice7A', image.alpha = image.alpha) 
    f<-SpatialFeaturePlot(data, gene, image = 'slice3A', image.alpha = image.alpha) 
    g<-SpatialFeaturePlot(data, gene, image = 'slice8B', image.alpha = image.alpha) 
    a+b+c+d+e+f+g & plot_layout(nrow = nrow) & scale_fill_gradientn(limits = c(0, scale), breaks = c(0,scale), colors = my_spatial_palette, oob=squish)
  }
}

my_cell2loc_palette<-colorRampPalette(c(magma(7, begin = 0.2)), bias = 1)(n = 256)
cell2loc_plots<-function(data, cluster, scale_cutoff = 0.992, assay = data@active.assay){
  active_assay <- DefaultAssay(data)
  gene_expression <- t(data@assays[[active_assay]]@data[cluster, , drop = FALSE])
  slide_info<-as.matrix(data$slide)
  gene_slide <- as.data.frame(cbind(gene_expression, slide_info))
  quantiles <- aggregate(gene_expression ~ V2, data = gene_slide, FUN = function(x) quantile(x, probs = scale_cutoff))
  x <- max(quantiles[,2])
  if (x < 0.05) {
    scale <- signif(x, digits = 2)
  } else {
    scale <- round(x, digits = 1)
  }
  a<-SpatialFeaturePlot(data, cluster, image = 'slice5A', image.alpha = 0)  
  b<-SpatialFeaturePlot(data, cluster, image = 'slice6B', image.alpha = 0) 
  c<-SpatialFeaturePlot(data, cluster, image = 'slice2B', image.alpha = 0) 
  d<-SpatialFeaturePlot(data, cluster, image = 'slice4B', image.alpha = 0) 
  e<-SpatialFeaturePlot(data, cluster, image = 'slice7A', image.alpha = 0) 
  f<-SpatialFeaturePlot(data, cluster, image = 'slice3A', image.alpha = 0) 
  g<-SpatialFeaturePlot(data, cluster, image = 'slice8B', image.alpha = 0)
  a+b+c+d+e+f+g & plot_layout(nrow = 1) & scale_fill_gradientn(limits = c(0, scale), breaks = c(0,scale), colors = my_cell2loc_palette, oob=squish)
}







##############

DotPlot_cell2loc <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA,
    threshold = 0.05
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  
  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = threshold)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}




SpatialDimPlot_ordered<-function(data, image.alpha, metadata, nrow, crop = FALSE){
  col <- scCustomize_Palette(num_groups = 40, ggplot_default_colors = FALSE)
  names(col) <- unique(st@meta.data[[metadata]])
  a<-SpatialDimPlot(data, image = 'slice5A', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE)  
  b<-SpatialDimPlot(data, image = 'slice6B', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE) 
  c<-SpatialDimPlot(data, image = 'slice2B', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE) 
  d<-SpatialDimPlot(data, image = 'slice4B', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE) 
  e<-SpatialDimPlot(data, image = 'slice7A', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE) 
  f<-SpatialDimPlot(data, image = 'slice3A', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE) 
  g<-SpatialDimPlot(data, image = 'slice8B', image.alpha = image.alpha, group.by = metadata, cols = col, crop = FALSE)
  a+b+c+d+e+f+g + plot_layout(nrow = nrow, guides = 'collect')
}




#Function to calculate the percentage of cells with greater than 0.01 mapping in a set of clusters. 

percent_map_per_cluster<-function(data = data, 
                                  nucseq_clusters = NULL, 
                                  assay = 'C6', 
                                  percentage_threshold = 0.01, 
                                  group_by = "seurat_clusters",
                                  split_by_slide = FALSE){
  #make sure clusters are in correct format
  if (!is.null(nucseq_clusters)) {
    nucseq_clusters <- gsub("-", ".", nucseq_clusters)
  }  
  x<-data.frame(t(data@assays[[assay]]@data))
  y<-data@meta.data[[group_by]]
  if (!is.null(nucseq_clusters)) {
    x <- x[, colnames(x) %in% nucseq_clusters]
  }
  if (split_by_slide == TRUE){
    print("Splitting clusters by tissue sections..")
    z<-data$slide
    x[x > percentage_threshold] <- 1
    x[x <= percentage_threshold] <- 0
    x$cluster <- y
    x$slide <- z
    x$cluster_slide<-paste(x$slide, x$cluster, sep = '_')
    total_cells_per_cluster_slide <- data.frame(table(x$cluster_slide))
    colnames(total_cells_per_cluster_slide)<- c('cluster_slide', 'total')
    x_sum<-aggregate(.~cluster_slide -cluster -slide, data = x, FUN = sum)
    x_sum<-merge(x_sum, total_cells_per_cluster_slide, by = 'cluster_slide', all = TRUE)
    row.names(x_sum)<-x_sum[,1]
    x_sum<-x_sum[-1]
    x_percent <- x_sum
    for (col in tail(colnames(x_sum), -2)) {
      x_percent[[col]] <- round((x_sum[[col]] / x_sum$total * 100), 2)    
    }
    x_percent <- subset(x_percent, select = -c(total))
    return(x_percent)
  }
  else if (split_by_slide == FALSE){
    print("Calculating percentage per cluster..")
    # Convert expression values to binary based on the threshold for numeric columns only
    x[x > percentage_threshold] <- 1
    x[x <= percentage_threshold] <- 0
    x$cluster <- y
    total_cells_per_cluster <- data.frame(table(x$cluster))
    colnames(total_cells_per_cluster)<- c('cluster', 'total')
    x_sum <- aggregate(. ~ cluster, data = x, FUN = sum)
    x_sum<-merge(x_sum, total_cells_per_cluster, by = 'cluster', all = TRUE)
    row.names(x_sum)<-x_sum[,1]
    x_sum<-x_sum[-1]
    x_percent <- x_sum
    for (col in colnames(x_sum)) {
      x_percent[[col]] <- x_sum[[col]]/ x_sum$total * 100
    }
    x_percent<-round(x_percent, 2)
    x_percent <- subset(x_percent, select = -c(total))
    return(x_percent)
  }
}

# top_regions_mapping<-function(data = st,
#                               )
# 
# 
#   find_top_regions<-function(table, cluster_name, number_of_regions = 3){
#     cluster_data<-table[, cluster_name, drop = FALSE]
#     sorted_data <- cluster_data[order(cluster_data[[1]], decreasing = TRUE), , drop = FALSE]
#     top_regions<-sorted_data[1:number_of_regions, , drop = FALSE]
#     return(top_regions)
#   }



cell2loc_avg_abundance<-function(data = data,
                                 assay = assay,
                                 regional_clusters = regional_clusters,
                                 average_function = mean,
                                 split_by_slide = FALSE){
  x<-data.frame(t(data@assays[[assay]]@data))
  y<-data@meta.data[[regional_clusters]]
  if (split_by_slide == TRUE){
    print('Calculating average cell abundance per region, per slide..')
    z<-data$slide
    x$cluster <- y
    x$slide <-z
    x$cluster_slide<-paste(x$slide, x$cluster, sep = '_')
    x <- subset(x, select = -c(slide, cluster))
    avex<-aggregate(.~cluster_slide, data = x, FUN = mean)
    rownames(avex) <- avex$cluster
    avex<-avex[,-c(1)]
    return(avex)
  }
  else if (split_by_slide == FALSE){
    print('Calculating average cell abundance per region..')
    x$cluster <- y
    avex<-aggregate(.~cluster, x, FUN = mean)
    rownames(avex) <- avex$cluster
    avex<-avex[,-c(1)]
    return(avex)
  }
}

#Flipping images in seurat object to allow for easier plotting. these functions were foudn from github... https://github.com/satijalab/seurat/issues/2702

# flip_angle %in% c(180, "R90", "L90", "Hf", "Vf")

rotimat=function(foo,rotation){
  if(!is.matrix(foo)){
    cat("Input is not a matrix")
    return(foo)
  }
  if(!(rotation %in% c("180","Hf","Vf", "R90", "L90"))){
    cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
    return(foo)
  }
  if(rotation == "180"){
    foo <- foo %>% 
      .[, dim(.)[2]:1] %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "Hf"){
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  
  if(rotation == "Vf"){
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "L90"){
    foo = t(foo)
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "R90"){
    foo = t(foo)
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  return(foo)
}

rotateSeuratImage = function(seuratVisumObject, slide = "slice1", rotation="Vf"){
  if(!(rotation %in% c("180","Hf","Vf", "L90", "R90"))){
    cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
    return(NULL)
  }else{
    seurat.visium = seuratVisumObject
    ori.array = (seurat.visium@images)[[slide]]@image
    img.dim = dim(ori.array)[1:2]/(seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- c()  
    # transform the image array
    for (rgb_idx in 1:3){
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # construct new rgb image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]],
                         new.mx[[2]],
                         new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    #swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    ## step4: change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    #swap index
    if(rotation == "Hf"){
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "Vf"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
    }
    
    if(rotation == "180"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "L90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2]-img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if(rotation == "R90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1]-img.index$imagerow
    }
    
    return(seurat.visium)
  }  
}


calc_NTs <- function(data, mouse = TRUE) {
  labeled_data <- data
  
  if (!'NTs' %in% colnames(labeled_data)) {
    labeled_data$NTs <- character(nrow(labeled_data))
  }
  
  for (i in 1:nrow(labeled_data)) {
    row <- labeled_data[i, ]
    labels <- c()
    
    if (mouse) {
      if (row['Slc17a6'] == 1) labels <- c(labels, 'GLUT2')
      if (row['Slc17a7'] == 1) labels <- c(labels, 'GLUT1')
      if (row['Slc17a8'] == 1) labels <- c(labels, 'GLUT3')
      if ((row['Gad1'] == 1 || row['Gad2'] == 1) && row['Slc32a1'] == 1) labels <- c(labels, 'GABA')
      if ((row['Gad1'] == 1 || row['Gad2'] == 1) && (row['Slc6a5'] == 1 || row['Slc6a9'] == 1)) labels <- c(labels, 'GLY')
      if (row['Chat'] == 1 && row['Slc5a7'] == 1) labels <- c(labels, 'Chol')
      if ((row['Slc6a3'] == 1 || row['Ddc'] == 1 || row['Th'] == 1) && row['Dbh'] == 0 && row['Pnmt'] == 0 && row['Tph2'] == 0 && row['Slc6a4'] == 0) labels <- c(labels, 'Dopa')
      if (row['Pnmt'] == 1 || row['Dbh'] == 1) labels <- c(labels, 'Noradr')
      if (row['Tph2'] == 1 && row['Slc6a4'] == 1) labels <- c(labels, '5-HT')
      if (row['Hdc'] == 1) labels <- c(labels, 'His');
    } else {
      if (row['SLC17A6'] == 1) labels <- c(labels, 'GLUT2')
      if (row['SLC17A7'] == 1) labels <- c(labels, 'GLUT1')
      if (row['SLC17A8'] == 1) labels <- c(labels, 'GLUT3')
      if ((row['GAD1'] == 1 || row['GAD2'] == 1) && row['SLC32A1'] == 1) labels <- c(labels, 'GABA')
      if ((row['GAD1'] == 1 || row['GAD2'] == 1) && (row['SLC6A5'] == 1 || row['SLC6A9'] == 1)) labels <- c(labels, 'GLY')
      if (row['CHAT'] == 1 && row['SLC5A7'] == 1) labels <- c(labels, 'Chol')
      if ((row['SLC6A3'] == 1 || row['DDC'] == 1 || row['TH'] == 1) && row['DBH'] == 0 && row['PNMT'] == 0 && row['TPH2'] == 0 && row['SLC6A4'] == 0) labels <- c(labels, 'Dopa')
      if (row['PNMT'] == 1 || row['DBH'] == 1) labels <- c(labels, 'Noradr')
      if (row['TPH2'] == 1 && row['SLC6A4'] == 1) labels <- c(labels, '5-HT')
      if (row['HDC'] == 1) labels <- c(labels, 'His');
    }
    
    if (length(labels) > 0) {
      if (!is.na(row['NTs']) && row['NTs'] != '') {
        existing_labels <- unlist(strsplit(row['NTs'], '_'))
        labels <- c(existing_labels, labels)
      }
      labeled_data[i, 'NTs'] <- paste(labels, collapse = '_')
    } else {
      labeled_data[i, 'NTs'] <- NA
    }
  }
  
  return(labeled_data)
}

