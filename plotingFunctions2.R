# do the enrichment using clusterprofiler for those unique genes in IMCN and MICN 
# ORA
perform_enricher <- function(gene_symbols, TERM2GENE) {
  result <- 
    clusterProfiler::enricher(
      gene_symbols,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = NULL,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      gson = NULL,
      TERM2GENE = TERM2GENE,
      TERM2NAME = NA) #%>% as.data.frame() 
  return(result)
}

# GSEA
perform_GSEA <- function(gene_symbols, TERM2GENE) {
  result <- 
    clusterProfiler::GSEA(
      gene_symbols,
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = 0.05,
      pAdjustMethod = "none",
      gson = NULL,
      TERM2GENE = TERM2GENE,
      TERM2NAME = NA,
      verbose = TRUE,
      seed = FALSE,
      by = "fgsea") #%>% as.data.frame() 
  return(result)
}

# Function to convert gene symbols to gene IDs
symbol_to_geneID <- function(gene_symbols) {
  # Load the human gene annotation database
  hs_db <- org.Hs.eg.db
  
  # Convert gene symbols to gene IDs
  gene_ids <- mapIds(hs_db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL")
  
  return(gene_ids)
}

# Function for disease ontology enrichment(ORA)
disease_enrichment <- function(gene_symbols, ont, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                               minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.05,
                               readable = TRUE) {
  # Convert gene symbols to gene IDs
  # Load the human gene annotation database
  hs_db <- org.Hs.eg.db
  
  # Convert gene symbols to gene IDs
  gene_ids <- mapIds(hs_db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL")
  
  # Perform disease ontology enrichment
  enriched_terms <- DOSE::enrichDO(
    gene = as.character(gene_ids),
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    qvalueCutoff = qvalueCutoff,
    readable = readable
  )
  
  return(enriched_terms)
}

# Function for disease ontology enrichment(GSEA)
disease_gseDO <- function(gene_symbols, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.05,
                          readable = TRUE) {
  # Convert gene symbols to gene IDs
  hs_db <- org.Hs.eg.db
  # Convert gene symbols to gene IDs
  gene_ids <- mapIds(hs_db, keys = names(gene_symbols), column = "ENTREZID", keytype = "SYMBOL")
  gene_ids2 <- setNames(gene_symbols, gene_ids)
  
  # Perform disease ontology enrichment
  enriched_terms <- DOSE::gseDO(
    gene = gene_ids2,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    exponent = 1,
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea"
  )
  
  return(enriched_terms)
}

# enrichment bar plots

# Custom function to split the pathway names into two lines if needed
split_description <- function(description, threshold = 20) {
  description <- as.character(description) # Convert to character (in case of factors)
  words <- strsplit(description, " ")[[1]]
  num_words <- length(words)
  if (nchar(description) > threshold && num_words > 2) {
    first_line <- paste(words[1:2], collapse = " ")
    second_line <- paste(words[3:num_words], collapse = " ")
    return(paste(first_line, second_line, sep = "\n"))
  } else {
    return(description)
  }
}

# enrichment bar plots where enriched pathways are highlighted tissue wise
enrichment_barplot <- function(data,
                               labs_title = NULL, 
                               labs_subtitle = NULL) {
  
  # Split pathway names into two lines if needed using the custom function
  # data$Description <- sapply(data$Description, function(x) paste(strwrap(x, width = 25)[1:2], collapse = "\n"))
  data$Description <- sapply(data$Description, split_description)
  
  # reorder the data by pvalue and assign colors to tissues
  data_reorder <- data %>% 
    mutate(Description = fct_reorder(Description, -log10(p.adjust),.desc = F))
  
  # create the ggplot object with the data and aesthetics
  p <- ggplot(data_reorder, aes(x = Description, y = -log10(p.adjust), fill = Tissue)) +
    # add the bar chart layer
    geom_bar(stat = "identity", position = position_dodge()) + 
    scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                                 'muscle' = 'darkmagenta')) +
    # add a horizontal line layer
    geom_hline(yintercept = 1.3, linetype = 2, colour='black', size = 1) +
    # set the theme
    theme_minimal() +
    # flip the coordinates to create a horizontal bar chart
    coord_flip() +
    # set the fill colors manually
    #scale_fill_manual(values = tissue_colors) +
    #set the axis and legend elements
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position = "right",
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold")) +
    # set the axis labels
    labs(title = labs_title, subtitle = labs_subtitle) +
    labs(x = "",
         y = bquote("" ~ -Log[10] ~ "p.adjust" ~ ""))
  
  return(p)
}


# general enrichment bar plots 
enrichment_1_barplot <- function(data,
                                 labs_title = NULL, 
                                 labs_subtitle = NULL) {
  data$Description <- sapply(data$Description, split_description)
  
  # reorder the data by pvalue and assign colors to tissues
  data_reorder <- data %>% 
    mutate(Description = fct_reorder(Description, -log10(p.adjust),.desc = F))
  
  # create the ggplot object with the data and aesthetics
  p <- ggplot(data_reorder, aes(x = Description, y = -log10(p.adjust))) +
    # add the bar chart layer
    geom_bar(stat = "identity", position = position_dodge()) + 
    # add a horizontal line layer
    geom_hline(yintercept = 1.3, linetype = 2, colour='black', size = 1) +
    # set the theme
    theme_minimal() +
    # flip the coordinates to create a horizontal bar chart
    coord_flip() +
    # set the fill colors manually
    #scale_fill_manual(values = tissue_colors) +
    #set the axis and legend elements
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position = "right",
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold")) +
    # set the axis labels
    labs(title = labs_title, subtitle = labs_subtitle) +
    labs(x = "",
         y = bquote("" ~ -Log[10] ~ "p.adjust" ~ ""))
  
  return(p)
}


# enrichment bar plots enriched pathways are highlighted network wise
enrichment_2_barplot <- function(data,
                                 labs_title = NULL, 
                                 labs_subtitle = NULL) {
  # Split pathway names into two lines if needed using the custom function
  # data$Description <- sapply(data$Description, function(x) paste(strwrap(x, width = 25)[1:2], collapse = "\n"))
  data$Description <- sapply(data$Description, split_description)
  
  # reorder the data by pvalue and assign colors to tissues
  data_reorder <- data %>% 
    mutate(Description = fct_reorder(Description, -log10(p.adjust),.desc = F))
  
  # create the ggplot object with the data and aesthetics
  p <- ggplot(data_reorder, aes(x = Description, y = -log10(p.adjust), fill = network)) +
    # add the bar chart layer
    geom_bar(stat = "identity", position = position_dodge()) + 
    scale_fill_manual(values = c('IMCN'  = 'black' ,
                                 'MICN' = 'gray')) +
    # add a horizontal line layer
    geom_hline(yintercept = 1.3, linetype = 2, colour='black', size = 1) +
    # set the theme
    theme_minimal() +
    # flip the coordinates to create a horizontal bar chart
    coord_flip() +
    # set the fill colors manually
    #scale_fill_manual(values = tissue_colors) +
    #set the axis and legend elements
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position = "right",
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold")) +
    # set the axis labels
    labs(title = labs_title, subtitle = labs_subtitle) +
    labs(x = "",
         y = bquote("" ~ -Log[10] ~ "p.adjust" ~ ""))
  
  return(p)
}

# enrichment bar plots where enriched pathways are highlighted as positive or negative 
enrichment_3_barplot <- function(data,
                                 labs_title = NULL, 
                                 labs_subtitle = NULL) {
  data$Description <- sapply(data$Description, split_description)
  
  # reorder the data by pvalue and assign colors to tissues
  data_reorder <- data %>% 
    mutate(Description = fct_reorder(Description, -log10(p.adjust),.desc = F))
  
  # create the ggplot object with the data and aesthetics
  p <- ggplot(data_reorder, aes(x = Description, y = -log10(p.adjust), fill = genes)) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    scale_fill_manual(values = c('pos'  = 'dodgerblue' ,
                                 'neg' = 'gray')) +
    # add a horizontal line layer
    geom_hline(yintercept = 1.3, linetype = 2, colour='black', size = 1) +
    # set the theme
    theme_minimal() +
    # flip the coordinates to create a horizontal bar chart
    coord_flip() +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position = "right",
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold")) +
    # set the axis labels
    labs(title = labs_title, subtitle = labs_subtitle) +
    labs(x = "",
         y = bquote("" ~ -Log[10] ~ "p.adjust" ~ ""))
  
  return(p)
}

# enrichment bar plots for GSEA
enrichment_barplot_GSEA <- function(data,
                                    labs_title = NULL, 
                                    labs_subtitle = NULL) {
  ggplot(data, aes(x = reorder(Description, 
                               NES,decreasing = TRUE),
                   y = round(NES,digits = 2), 
                   fill = Tissue)) + 
    coord_flip()+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = 
                        c("IMAT" ="darkorange",  
                          "muscle" = "darkmagenta"))+
    coord_flip() +
    labs(x = "", 
         y = "enrichment Score",
         tag = "") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold"),
          legend.position = "left",
          strip.text  = element_text(size = 15))
}


# Heatplot 
heatplot_plot = function(data) {
  enrichplot::heatplot(data) +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold"),
          legend.position = "right",
          strip.text  = element_text(size = 15))
}

ridgeplot_plot = function(data) {
  enrichplot::ridgeplot(data) +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold"),
          legend.position = "right",
          strip.text  = element_text(size = 15))
}

# show it in dot plot
pathwaysDotplot <- function(data, 
                            n, 
                            caption = NULL,
                            title = NULL,
                            hjust = NULL,
                            y_axis_title = NULL, 
                            legend_position = "bottom", 
                            legend_direction = "horizontal", 
                            color_scale_low = "green",
                            color_scale_mid = 'white', 
                            color_scale_high = "red",
                            fill_scale_low = 'lightblue', 
                            fill_scale_mid = "white",
                            fill_scale_high = "dodgerblue",
                            shape_values = c('IMAT' = 21, 'muscle' = 23), 
                            title_size = 20, axis_title_size = 12, 
                            axis_text_size = 8, strip_text_size = 12) {
  
  # Subset data based on n value
  if (length(unique(data$geneID)) > 40) {
    df1 <- data[1:nrow(data)/2,]
    df2 <- data[(nrow(data)/2+1):(nrow(data)),]
    p <- ggplot(data = df1, aes(x = Description, 
                                y = geneID, 
                                shape = Tissue,
                                color = Rewiring,
                                fill = corGIR)) +
      geom_point(alpha = 1, size = 4) +
      scale_colour_gradient2(low = color_scale_low, 
                             mid = color_scale_mid,
                             high = color_scale_high) +
      scale_fill_gradient2(low = fill_scale_low, 
                           mid = fill_scale_mid,
                           high = fill_scale_high) +
      scale_shape_manual(values = shape_values,
                         guide = guide_legend()) +
      coord_flip() +
      theme_bw() +
      theme(axis.title.y = element_text(size = axis_title_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.text.x = element_text(size = 12, angle = 45),
            axis.text.y = element_text(size = axis_title_size),
            strip.text = element_text(size = strip_text_size),
            legend.position = legend_position,
            legend.direction = legend_direction,
            legend.spacing.x = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.box.margin = margin(0, 0, 4, 0)) +
      labs(x = "", y = "") +
      theme(legend.spacing.x = unit(0, "cm")) +
      # guides(
      #   shape = guide_legend(order = 1),
      #   color = guide_legend(order = 2),
      #   fill = guide_legend(order = 3)) +
      patchwork::plot_annotation(title = title,
                                 caption = caption,
                                 theme = theme(plot.title = element_text(size = title_size, hjust = hjust, face = "bold")))
    
    if (!is.null(y_axis_title)) {
      p <- p + ylab(y_axis_title)
    } else {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle("")
    }
    
    # Plot using ggplot2 with pvalue as fill color
    p <-  p + ggplot(data = df2, aes(x = Description, 
                                     y = geneID, 
                                     shape = Tissue,
                                     color = Rewiring, 
                                     fill = corGIR)) +
      geom_point(alpha = 1, size = 4) +
      scale_colour_gradient2(low = color_scale_low, 
                             mid = color_scale_mid,
                             high = color_scale_high) +
      scale_fill_gradient2(low = fill_scale_low, 
                           mid = fill_scale_mid,
                           high = fill_scale_high) +
      scale_shape_manual(values = shape_values,
                         guide = guide_legend()) +
      coord_flip() +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = axis_title_size),
            axis.text.x = element_text(size = 12, angle = 45),
            axis.text.y = element_blank(),
            strip.text = element_text(size = strip_text_size),
            legend.position = legend_position,
            legend.direction = legend_direction,
            legend.spacing.x = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.box.margin = margin(0, 0, 4, 0)) +
      labs(x = "", y = "") +
      theme(legend.spacing.x = unit(0, "cm")) +
      # guides(
      #   shape = guide_legend(order = 1),
      #   color = guide_legend(order = 2),
      #   fill = guide_legend(order = 3)) +
      patchwork::plot_annotation(title = title,
                                 caption = caption,
                                 theme = theme(plot.title = element_text(size = title_size, hjust = hjust, face = "bold")))
    
    if (!is.null(y_axis_title)) {
      p <- p + ylab(y_axis_title)
    } else {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    print(p)
    
  } else {
    p <- ggplot(data, aes(x = Description, 
                          y = geneID, 
                          shape = Tissue,
                          color = Rewiring, 
                          fill = corGIR
    )) +
      geom_point(alpha = 1, size = 4) +
      scale_colour_gradient2(low = color_scale_low, 
                             mid = color_scale_mid,
                             high = color_scale_high) +
      scale_fill_gradient2(low = fill_scale_low, 
                           mid = fill_scale_mid,
                           high = fill_scale_high) +
      scale_shape_manual(values = shape_values,
                         guide = guide_legend()) +
      coord_flip() +
      theme_bw() +
      theme(axis.title.y = element_text(size = axis_title_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.text.x = element_text(size = 12, angle = 45),
            axis.text.y = element_text(size = axis_title_size),
            strip.text = element_text(size = strip_text_size),
            legend.position = legend_position,
            legend.direction = legend_direction,
            legend.spacing.x = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.box.margin = margin(0, 0, 4, 0)) +
      labs(x = "", y = "") +
      theme(legend.spacing.x = unit(0, "cm")) +
      patchwork::plot_annotation(title = title,
                                 caption = caption,
                                 theme = theme(plot.title = element_text(size = title_size, hjust = hjust, face = "bold")))
    print(p)
  }
  
}


# circus plot function 
plotcircosinR_from_to_cor <- function(data,
                                      gene_col = NULL,
                                      transparency = 0.5,
                                      link.arr.lwd = 1,
                                      link.arr.lty = NULL,
                                      link.arr.col = NULL,
                                      link.arr.width = NULL,
                                      link.arr.type = NULL,
                                      facing = "clockwise",
                                      cell_col = NULL,
                                      print.cell = NULL,
                                      print.cell_layer = TRUE,
                                      track.height_1 = uh(4, "mm"),
                                      track.height_2 = uh(14, "mm"),
                                      annotation.height_1 = 0.02,
                                      annotation.height_2 = 0.02,
                                      text.vjust = "0.6cm",
                                      ...)
{
  ## Highlighting this IS genes in the circus plot
  cell_group <- unique(c(data$cell_Sender, data$cell_Receiver))
  genes <- c(
    structure(data$Sender, names = data$cell_Sender),
    structure(data$Receiver, names = data$cell_Receiver)
  )
  genes <- genes[!duplicated(paste(names(genes), genes))]
  genes <- genes[order(names(genes))]
  
  ## IMAT/senders/ correlated to GIR
  genes_2 <-
    c(structure(data$Sender_Cor_to_GIR_estimate, names = data$cell_Sender))
  genes_2 <- genes_2[!is.na(genes_2)]
  genes_2 <- genes_2[!duplicated(paste(names(genes_2), genes_2))]
  genes_2 <- genes_2[order(names(genes_2))]
  
  
  ## muscle/receiver/ correlated to GIR
  genes_4 <-
    c(structure(data$Receiver_Cor_to_GIR_estimate, names = data$cell_Receiver))
  genes_4 <- genes_4[!is.na(genes_4)]
  genes_4 <- genes_4[!duplicated(paste(names(genes_4), genes_4))]
  genes_4 <- genes_4[order(names(genes_4))]
  
  
  if (is.null(link.arr.lty)) {
    link.arr.lty = "solid"
  }
  if (is.null(link.arr.col)) {
    data <-
      data %>% mutate(link_col ="#2166ac")
  }
  else {
    data$link_col = link.arr.col
  }
  if (is.null(link.arr.type)) {
    link.arr.type = "triangle"
  }
  if (is.null(gene_col)) {
    geneSendr_col <- structure(
      c(RColorBrewer::brewer.pal(7, "Paired")),
      names = c(
        "Cytokines",
        "Growth-Factors",
        "Junctional",
        "Neuropeptides",
        "Predicted-Ligands",
        "Secreted_Glycoprotein",
        "Secreted-Proteins"
      )
    )
    
    geneReciver_col <-
      structure(
        c(RColorBrewer::brewer.pal(4, "Dark2")),
        names = c(
          "Adhesion",
          "Enzyme-Receptors",
          "GPCRs",
          "Others"
        )
      )
    gene_col <-
      structure(c(geneSendr_col[data$Sender_gene_Class], 
                  geneReciver_col[data$Receiver_gene_Class]),
                names = c(data$Sender, data$Receiver))
    
  }
  
  if (is.null(cell_col)) {
    cell_col <-
      structure(randomColor(count = length(unique(names(
        genes
      ))), luminosity = "dark"),
      names = unique(names(genes)))
  }
  
  if (is.null(link.arr.lwd)) {
    data <- data %>% mutate(arr_width = 1)
  }
  else if (max(abs(link.arr.lwd)) - min(abs(link.arr.lwd)) == 0 &&
           all(link.arr.lwd != 1e-04)) {
    data <-
      data %>% mutate(arr_width = ifelse(abs(link.arr.lwd < 5), abs(link.arr.lwd), 5))
  }
  else {
    data <-
      data %>% mutate(arr_width = ifelse(link.arr.lwd == 1e-04, 2, 1 + 5 / (max(
        abs(link.arr.lwd)
      ) - min(
        abs(link.arr.lwd)
      )) *
        (abs(link.arr.lwd) - min(
          abs(link.arr.lwd)
        ))))
  }
  if (length(cell_group) != 1) {
    gap.degree <-
      do.call("c", lapply(table(names(genes)), function(i)
        c(rep(1, i - 1), 8)))
  }
  else {
    gap.degree <-
      do.call("c", lapply(table(names(genes)), function(i)
        c(rep(1, i))))
  }
  circos.par(gap.degree = gap.degree)
  if (length(gene_col) == 1) {
    grid.col = gene_col
  }
  else {
    grid.col = gene_col[genes]
    names(grid.col) <- paste(names(genes), genes)
  }
  
  if (is.null(link.arr.width)) {
    data <- data %>% mutate(link.arr.width = data$arr_width / 10)
  }
  else if (max(abs(link.arr.width)) - min(abs(link.arr.width)) == 0 &&
           all(link.arr.width != 1e-04)) {
    data <-
      data %>% mutate(link.arr.width = ifelse(abs(link.arr.width) <
                                                0.5, abs(link.arr.width), 0.5))
  }
  else {
    data <-
      data %>% mutate(link.arr.width = ifelse(link.arr.width == 1e-04, 0.2, (1 + 5 /
                                                                               (
                                                                                 max(abs(link.arr.width)) - min(abs(link.arr.width))
                                                                               ) *
                                                                               (
                                                                                 abs(link.arr.width) - min(abs(link.arr.width))
                                                                               )) / 10))
  }
  chordDiagram(
    as.data.frame(cbind(
      paste(data$cell_Sender, data$Sender),
      paste(data$cell_Receiver, data$Receiver)
    )),
    order = paste(names(genes),  genes),
    grid.col = grid.col,
    transparency = transparency,
    directional = 1,
    direction.type = "arrows",
    link.arr.lwd = data$arr_width,
    link.arr.lty = link.arr.lty,
    link.arr.type = link.arr.type,
    link.arr.width = data$link.arr.width,
    link.arr.col = data$link_col,
    col = "#00000000",
    annotationTrack = c("grid"),
    preAllocateTracks = list(
      list(track.height = track.height_1),
      list(track.height = track.height_2)
    ),
    annotationTrackHeight = c(annotation.height_1, annotation.height_2),
    ...
  )
  
  circos.trackPlotRegion(
    track.index = 2,
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.index = genes[get.cell.meta.data("sector.numeric.index")]
      circos.text(
        mean(xlim),
        mean(ylim),
        labels =  sector.index,
        col =  ifelse(
          sector.index %in% unlist(genes_2),
          'red',
          ifelse(
            sector.index %in% unlist(genes_4),
            'red', 'black')
        ),
        
        cex = 0.9,
        facing = facing,
        niceFacing = TRUE
      )
      
    },
    bg.border = 0
  )
  
  
  if (print.cell) {
    for (c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(
        sector.index = paste(c, gene),
        track.index = 1,
        col = ifelse(length(cell_col) == 1, cell_col,
                     cell_col[c]),
        #if you wan't print the tissue name on the highlighted track, uncommint the code below
        #text = c,
        
        #if you won't print the name tissue on the highlighted track
        text = "",
        #to make the title horizontal for(circos.par(start.degree = 270) # vertically symmetric)
        facing = "downward",
        cex = 2,
        text.vjust = text.vjust,
        niceFacing = TRUE,
        lwd = 1
      )
    }
  }
  
  circos.clear()
}


# hallmark pathway list
hallmark <- msigdbr(species = "Homo sapiens",category = "H") %>% data.frame()%>%
  dplyr::select(gs_name,human_gene_symbol) %>%
  dplyr::filter(gs_name %in% c("HALLMARK_ADIPOGENESIS",
                               "HALLMARK_MYOGENESIS",
                               "HALLMARK_ANGIOGENESIS",
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                               "HALLMARK_APOPTOSIS",
                               "HALLMARK_COAGULATION",
                               "HALLMARK_HYPOXIA",
                               "HALLMARK_INFLAMMATORY_RESPONSE",
                               "HALLMARK_COMPLEMENT",
                               "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                               "HALLMARK_APICAL_JUNCTION",
                               "HALLMARK_APICAL_SURFACE",
                               "HALLMARK_PEROXISOME",
                               "HALLMARK_GLYCOLYSIS",
                               "HALLMARK_FATTY_ACID_METABOLISM",
                               "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                               "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                               "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                               "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                               "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                               "HALLMARK_TGF_BETA_SIGNALING",
                               "HALLMARK_PI3K_AKT_MTOR_SIGNALING"))
hallmark$gs_name <- gsub("HALLMARK_", "", hallmark$gs_name) # Remove "HALLMARK_"
hallmark$gs_name <- gsub("_", " ", hallmark$gs_name)  # Remove underscore

# KEGG pathway list
KEGGpa = readr::read_tsv("/Users/amare.wolide/Desktop/PhD/finalIMATmuscle/percentile10/selected.KEGG.txt") %>%
  dplyr::pull(pathway)

KEGG <- msigdbr(species = "Homo sapiens",category = "C2", 
                subcategory = "CP:KEGG") %>% data.frame()%>%
  dplyr::select(gs_description,human_gene_symbol) %>%
  dplyr::filter(gs_description %in% KEGGpa)

# Reactome pathway list 
list_ReactomePA = readr::read_tsv("/Users/amare.wolide/Desktop/PhD/finalIMATmuscle/percentile10/selected.Reactome.txt") %>%
  dplyr::pull(pathway)
ReactomePA = msigdbr(species = "Homo sapiens", category = "C2",
                     subcategory = "CP:REACTOME") %>% data.frame()%>%
  dplyr::select(gs_description,human_gene_symbol) %>%
  dplyr::filter(gs_description %in% list_ReactomePA)


# disease ontology list 
list_DO <- c(
  "hyperinsulinemic hypoglycemia",
  "insulin resistance",
  "hyperinsulinemia",
  "abnormal insulin level",
  "hypoinsulinemia",
  "insulin insensitivity",
  "diabetes mellitus",
  "type II diabetes mellitus",
  "maturity-onset diabetes of the young",
  "obesity",
  "abdominal obesity",
  "truncal obesity",
  "arteriosclerosis",
  "muscular disease",
  "collagen disease",
  "myopathy",
  "myotonic disease",
  "skeletal muscle atrophy",
  "fatty replacement of skeletal muscle",
  "abnormal skeletal muscle morphology",
  "fatigable weakness of skeletal muscles",
  "skeletal muscle fibrosis",
  "skeletal muscle hypertrophy",
  "muscular dystrophy",
  "atherosclerosis",
  "muscle tissue disease",
  "type 2 diabetes mellitus",
  "lipid storage disease",
  "increased intramyocellular lipid droplets",
  "increased muscle lipid content",
  "abnormality of glycolipid metabolism",
  "abnormality of metabolism/homeostasis",
  "abnormality of vitamin D metabolism",
  "abnormality of vitamin metabolism",
  "abnormality of mucopolysaccharide metabolism",
  "abnormal circulating lipid concentration",
  "hyperlipidemia",
  "hyperinsulinism",
  "muscular atrophy",
  "muscle fiber necrosis",
  "quadriceps muscle weakness",
  "muscle abnormality related to mitochondrial dysfunction",
  "abnormal muscle physiology",
  "insulin sensitivity",
  "Insulin resistance syndrome",
  "hyperglycemia",
  "postprandial hyperglycemia",
  "insulin-resistant diabetes mellitus",
  "insulin resistance in diabetes",
  "increased sarcoplasmic glycoge",
  "muscle damage",
  "muscle degeneration",
  "abnormal muscle glycogen content",
  "skeletal dysplasia",
  "increased muscle glycogen content"
)

# Function to combine data frames with missing pathways
combine_data_frames <- function(df1, df2) {
  # Find the pathways that are present in df1 but not in df2
  missing_pathways_df2 <- setdiff(df1$Description, df2$Description)
  
  # Create a new data frame with missing pathways for df2
  missing_rows_df2 <- data.frame(
    Description = missing_pathways_df2,
    GeneRatio = "0/0",
    BgRatio = "0/0",
    pvalue = 0,
    p.adjust = 0,
    qvalue = 0,
    network = "MICN"
  )
  
  # Combine data frames using rbind, including the missing rows for df2
  combined_df <- rbind(df1, df2, missing_rows_df2)
  
  # Find the pathways that are present in df2 but not in df1
  missing_pathways_df1 <- setdiff(df2$Description, df1$Description)
  
  # Create a new data frame with missing pathways for df1
  missing_rows_df1 <- data.frame(
    Description = missing_pathways_df1,
    GeneRatio = "0/0",
    BgRatio = "0/0",
    pvalue = 0,
    p.adjust = 0,
    qvalue = 0,
    network = "IMCN"
  )
  
  # Combine data frames using rbind, including the missing rows for df1
  combined_df <- rbind(combined_df, missing_rows_df1)
  
  # Return the final combined data frame
  return(combined_df)
}