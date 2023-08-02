## PCA plots 
pca_plot <- function(expression, 
                     Sample) {
  # package 
  require(ggplot2)
  require(ggplotify)
  require(ggpubr)
  require(ggrepel)
  # Check input data
  if (!is.matrix(expression)) {
    stop("gene expression data must be a matrix")
  }
  if (!is.data.frame(Sample)) {
    stop("sample information must be a data frame")
  }
  if (ncol(expression) != nrow(Sample)) {
    stop("gene expression data and sample information must have the same number of rows")
  }
  pca <- prcomp(expression)
  pca_rotation <- as.data.frame(pca[2]$rotation)
  pca_rotation$Tissue <-  Sample$Tissue
  pca_rotation$sample <- rownames(Sample)
  pca_rotation$group <- Sample$category
  percentage_pca <- round(pca$sdev^2 / sum(pca$sdev^2)*100,2)
  percentage_pca <- paste0(colnames(pca)[grep("^PC",colnames(pca))],
                           "(", paste0(as.character(percentage_pca), "%", ")")
  )
  
  plot <- ggplot(data=pca_rotation, aes(x=PC1, y=PC2, color=Tissue)) +
    geom_point(size=3) + theme_minimal() +
    scale_color_manual(values = c("darkorange","darkmagenta"))+
    theme(axis.title=element_text(size= 20),
          axis.text.x  = element_text(size = 18),
          axis.text.y  = element_text(size = 18),
          plot.title = element_text(hjust = 0.5,size=25),
          legend.title = element_text(color = "black", size = 10),
          legend.text=element_text(color = "black", size = 10),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.line = element_line(size = 1)) +
    theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text  = element_text(size = 15)) +
    guides(colour = guide_legend(override.aes = list(size=6)))+
    xlab(paste0("PC 1\n", percentage_pca[1])) +
    ylab(paste0("PC 1\n", percentage_pca[2]))
}

#heatmap
expression_heatmap = function(expression,
                              Sample, 
                              show_row_dend = TRUE, 
                              show_column_dend = FALSE,
                              show_column_names = FALSE,
                              show_row_names = FALSE,
                              cluster_columns = FALSE, 
                              cluster_rows = TRUE,
                              column_title = "Heatmap",
                              row_colors = NULL,
                              column_colors = NULL,
                              ...){
  
  # load package if it require
  require(ComplexHeatmap)
  require(circlize)
  require(grid)
  require(gridExtra)
  require(gridBase)
  
  # Set default row colors
  # if (is.null(row_colors)) {
  #   row_colors = list("Aver.exp" = colorRamp2(c(min(rowMeans(expression)), 
  #                                                  mean(rowMeans(expression)), 
  #                                                  max(rowMeans(expression))),
  #                                                c("white","#FFC0CB","#FF69B4")))
  # }
  
  # Set default column colors
  if (is.null(column_colors)) {
    column_colors = list(Tissue = c(IMAT = "darkorange", Muscle = "darkmagenta"),
                         Group = c("LC" = "lightblue", "ATH" = "#0073C2FF",
                                   "OB" = "#EFC000FF", "T2D" = "#CD534CFF"),
                         GIR = colorRamp2(c(min(Sample$GIR), 
                                            mean(Sample$GIR), 
                                            max(Sample$GIR)),
                                          c("green", "white", "red")))
  }
  
  # make sure sample row and expression column size are identical 
  stopifnot(ncol(expression) == nrow(Sample))
  
  # Check for missing values in expression data
  if (any(is.na(expression))) {
    stop("Expression data contains missing values.")
  }
  
  # Check that row colors are in correct format
  # if (!is.list(row_colors) || !all(names(row_colors) %in% rownames(expression))) {
  #   stop("Row colors are not in the correct format.")
  # }
  # Check that column colors are in correct format
  if (!is.list(column_colors) || !all(names(column_colors) %in% colnames(Sample))) {
    stop("Column colors are not in the correct format.")
  }
  
  # Create heatmap
  heatmapplot = grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        as.matrix(t(scale(t(expression)))),
        row_names_gp = gpar(fontsize = 8, fontface = "bold"),
        show_row_names = show_row_names,
        cluster_columns = cluster_columns,
        show_column_names = FALSE,
        use_raster = TRUE,
        show_column_dend = FALSE,
        show_heatmap_legend = TRUE,
        column_title = column_title,
        column_title_gp = gpar(fontsize = 25),
        raster_by_magick = TRUE,
        heatmap_legend_param = list(
          legend_direction = c("horizontal"),
          title = "Expr.level",
          grid_height = unit(3, "mm"),
          grid_width = unit(1.5, "mm"),
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10)),
        # left_annotation = ComplexHeatmap::rowAnnotation(gene= anno_block(gp = gpar(fill = 4:5),
        #                                                                  labels = c("high",
        #                                                                             "low"),
        #                                                                  labels_gp = gpar(col = c("white","black"),
        #fontsize = 15))),
        # right_annotation = HeatmapAnnotation(
        #   which = "row",
        #   "Aver.exp" = rowMeans(expression), 
        #   bar_plot = anno_barplot(rowMeans(expression)),
        #   show_legend = T,
        #   show_annotation_name = F,
        #   annotation_name_gp = gpar(fontsize = 20, fontface = "bold"),
        #   annotation_legend_param = list("Aver.exp" = list(direction = "horizontal")),
        #   border = TRUE,
        #   col = row_colors),
        top_annotation = HeatmapAnnotation(
          which = "column",
          #GIR = anno_barplot(as.numeric(Sample$GIR),bar_width = 1),
          #GIR = anno_boxplot(as.numeric(Sample$GIR), height = unit(4, "cm"))
          #GIR = as.numeric(Sample$GIR),
          #GIR = as.numeric(Sample$GIR),
          #Group = Sample$Group,
          Tissue = Sample$Tissue,
          col = column_colors,
          simple_anno_size = unit(0.5, "cm"),
          annotation_name_side = "right",
          show_legend = T,
          show_annotation_name = F,
          annotation_name_gp = gpar(fontsize = 15),
          annotation_legend_param = list(Tissue = list(nrow = 2)
                                         #Group = list(nrow = 2),
                                         #GIR  = list(direction = "horizontal")
          ),
          border = TRUE
        ),
        bottom_annotation = HeatmapAnnotation(
          which = "column",
          GIR = as.numeric(Sample$GIR),
          col = column_colors,
          simple_anno_size = unit(0.5, "cm"),
          annotation_name_side = "right",
          show_legend = T,
          show_annotation_name = F,
          annotation_name_gp = gpar(fontsize = 15),
          annotation_legend_param = list(GIR  = list(direction = "horizontal")
          ),
          border = TRUE
        )
      ),
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom"
    )
  )
}

# Heatmap and PCA for the processed data 
# PCA
PCAplot = pca_plot(expression = as.matrix(expression_df), 
                   Sample = Sample) +
  ggtitle("PCA")

# heatmap
annotation_col <- data.frame(Tissue = Sample$Tissue,
                             GIR = Sample$GIR,
                             Group = Sample$Group)
rownames(annotation_col) <- colnames(expression_df)
expressionheatmap <- expression_heatmap(expression = expression_df, 
                                        Sample = annotation_col,
                                        cluster_columns = F,
                                        show_row_names = FALSE,
                                        column_title = "Heatmap")

plot.new()
grid1 = 
  cowplot::plot_grid(
    expressionheatmap,
    PCAplot,
    ncol = 2,
    rel_heights = c(1,1),
    rel_widths = c(1,1),
    labels = c("A","B"))  + 
  patchwork::plot_annotation(title = "Data preprocessing", subtitle = "") &  
  theme(plot.title = element_text(hjust = 0.5,size=25,face="bold")) &
  theme(plot.subtitle = element_text(hjust = 0.5,size=15,face="bold"))

outline1 <- ggdraw() +
  draw_plot(grid1) +
  draw_grob(rectGrob(gp = gpar(col = "black",
                               lwd = 5, fontsize = 25,
                               fill = NA)))
cowplot::plot_grid(outline1)

# figure 1
ggsave("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/Suppfigure1.pdf",
       width = 35, height = 15, units = "cm",
       cowplot::plot_grid(outline1))

ggsave("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/Suppfigure1.svg",
       width = 35, height = 15, units = "cm",
       cowplot::plot_grid(outline1))

ggsave("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/Suppfigure1.png",
       width = 35, height = 15, units = "cm",
       cowplot::plot_grid(outline1))


##-----correlation to GIR-------------------------------------------------------
# histogram for the distribution of IMAT or muscle to GIR
histGIR = 
  ggplot(bindBoth, 
         aes(y = estimate, fill = Tissue)) + 
  geom_histogram(color = "gray", position = "stack") +
  scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                               'muscle' = 'darkmagenta'))+
  theme_bw()+
  ggtitle("Histogram") +
  xlab("Count")+
  ylab("Cor.coff")+
  labs(tag = "") +
  theme(axis.title=element_text(size= 18),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        plot.title = element_text(hjust = 0.5,size=25),
        legend.title = element_text(color = "black", size = 10),
        legend.text=element_text(color = "black", size = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.line = element_line(size = 1)) +
  theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
  theme(strip.background = element_blank(), strip.placement = "outside") +
  theme(strip.text  = element_text(size = 15)) + coord_flip() 

##-------------------for sender and receivers 
SR = read.csv("~/Desktop/PhD/Ligand_Receptor/cellcell_data.csv") %>%
  dplyr::select(Pair.Name,Sender_gene.ApprovedSymbol,Receiver_gene.ApprovedSymbol)
SRlist = SR %>% dplyr::select(Sender_gene.ApprovedSymbol,Receiver_gene.ApprovedSymbol) %>% unlist () %>% unique()

# now filter the genes based on the SRlist and get the sample data 
Sample2 <- data.table::fread("~/Desktop/PhD/IMATmuscle_network/All_SampleAnnotation_Tissue.csv") %>%
  as.data.frame() %>% 
  dplyr::mutate(Tissue = stringr::str_remove(Tissue, "\\s")) %>%
  dplyr::group_by(Pb_ID) %>%
  dplyr::rename(Group = category) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::arrange(desc(GIR)) %>%
  dplyr::mutate_if(.,is.character, factor)  %>%
  tibble::column_to_rownames(var = "V1")  %>%
  dplyr::arrange(desc(Tissue))

expression_df = expression[rownames(expression) %in% SRlist, rownames(Sample2)]
colnames(expression_df) %in% rownames(Sample2)
all(rownames(Sample2) %in% colnames(expression_df))
identical(colnames(expression_df), rownames(Sample2))

# all(rownames(Sample2) == colnames(expression_df))
# Sample2 = Sample2[match(rownames(Sample2), colnames(expression_df)), ]
# expression_df <- expression_df[, rownames(Sample2)]
# Sample2 <- Sample2[colnames(expression_df), ]
# all(rownames(Sample2) == colnames(expression_df))
# identical(colnames(expression_df), rownames(Sample2))

plot.new()
expressionheatmap2 <- expression_heatmap(expression = expression_df, 
                                         Sample = Sample2,
                                         show_row_names = FALSE,
                                         column_title = "Expression level")

cowplot::plot_grid(expressionheatmap2)




# Histogram for the correlation coefficient between gene expression and GIR and for edge weight,
# and the bar plot for number of nodes and edges in the IMAT to muscle communication network
cowplot::plot_grid(
  ggplot(nodeStatimatMuscle, 
         aes(y = corGIR, fill = Tissue)) + 
    geom_histogram(color = "gray", position = "stack") +
    scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                                 'muscle' = 'darkmagenta'))+
    theme_bw()+
    ggtitle("Cor.coefficient") +
    xlab("Count")+
    ylab("Cor.coff")+
    labs(tag = "A") +
    theme(plot.title = element_text(size = 25,hjust = 0.5))+
    theme(axis.title.y = element_text(size = 20))+
    theme(axis.title.x = element_text(size = 20))+
    theme(axis.text = element_text(size = 20))+
    theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(legend.position="left",legend.direction="vertical")+
    theme(strip.text  = element_text(size = 15)) + ggplot2::coord_flip(),
  
  ggplot(nodeStatimatMuscle, aes(x=Tissue, fill = Tissue)) +
    geom_bar(position = "dodge",
             #fill = "white", 
             size = 0.7) +
    scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                                 'muscle' = 'darkmagenta'))+
    theme_bw()+
    ggtitle("Number of Genes") +
    ylab("")+
    xlab("Tissue")+
    labs(tag = "B") +
    theme(axis.title.y = element_text(size = 20))+
    theme(axis.title.x = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 20))+
    theme(axis.text.y = element_text(size = 20))+
    theme(strip.text  = element_text(size = 20)) +
    theme(plot.title = element_text(size = 25,hjust = 0.5)) +
    theme(legend.position="none",legend.direction="vertical"),
  
  ggplot(EdgeStatimatMuscle, aes(x=EdgeWeight_addingCorcoffs, fill = Edges)) +
    geom_histogram(color = "gray", position = "stack") +
    scale_fill_manual(values = c('positive' = 'red',
                                 'negative' = 'blue'))+
    theme_bw()+
    ggtitle("Edge weight") +
    xlab("Sum of Cor.coffs")+
    ylab("Count")+
    labs(tag = "C") +
    theme(plot.title = element_text(size = 25,hjust = 0.5))+
    theme(axis.title.y = element_text(size = 20))+
    theme(axis.title.x = element_text(size = 20))+
    theme(axis.text = element_text(size = 20))+
    theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text  = element_text(size = 20)) +         
    theme(legend.position="left",legend.direction="vertical"),
  
  ggplot(EdgeStatimatMuscle, aes(x=Edges, fill = Edges)) +
    geom_bar(position = "stack",
             #fill = "gray",
             size=0.7) +
    ggtitle("Number of Edges") +
    theme_bw()+
    ylab("")+
    xlab("Edges")+
    labs(tag = "D") +
    scale_fill_manual(values = c('positive' = 'red',
                                 'negative' = 'blue',alpha = 4))+
    scale_y_continuous(breaks = c(0,100,200,300,400,500,600))+
    theme(plot.title = element_text(size = 25,hjust = 0.5))+
    theme(axis.title.y = element_text(size = 20))+
    theme(axis.title.x = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 20))+
    theme(axis.text.y = element_text(size = 20))+
    theme(strip.text  = element_text(size = 20)) +
    theme(legend.position="none",legend.direction="vertical"),
  ncol = 2,
  rel_widths = c(1,0.5,1,0.5))+
  patchwork::plot_annotation(title = "Histograms for gene expression-GIR correlation, edge weight,\n and bar plot for IMAT-muscle network size",
                             subtitle = "") &  
  theme(plot.title = element_text(hjust = 0.5,size=25,face="bold")) &
  theme(plot.subtitle = element_text(hjust = 0.5,size=15,face="bold"))


# number of genes and edges in the muscle to IMAT  
ggsave("~/Desktop/PhD/finalIMATmuscle/percentile10/MuscletoIMATNodeEdgebar.pdf",
       width = 30, height = 25, units = "cm",
       plot = 
         cowplot::plot_grid(
           ggplot(nodeStatMuscleimat, 
                  aes(y = corGIR, fill = Tissue)) + 
             geom_histogram(color = "gray", position = "stack") +
             scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                                          'muscle' = 'darkmagenta'))+
             scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
             theme_bw()+
             ggtitle("Cor.coefficient") +
             xlab("Count")+
             ylab("Cor.coff")+
             labs(tag = "A") +
             theme(plot.title = element_text(size = 25,hjust = 0.5))+
             theme(axis.title.y = element_text(size = 20))+
             theme(axis.title.x = element_text(size = 20))+
             theme(axis.text = element_text(size = 20))+
             theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
             theme(strip.background = element_blank(), strip.placement = "outside") +
             theme(legend.position="left",legend.direction="vertical")+
             theme(strip.text  = element_text(size = 15)) + ggplot2::coord_flip(),
           
           
           ggplot(nodeStatMuscleimat, aes(x=Tissue, fill = Tissue)) +
             geom_bar(position = "dodge",
                      #fill = "white", 
                      size = 0.7) +
             scale_fill_manual(values = c('IMAT'  = 'darkorange' ,
                                          'muscle' = 'darkmagenta'))+
             theme_bw()+
             ggtitle("Number of Genes") +
             ylab("")+
             xlab("Tissue")+
             labs(tag = "B") +
             theme(axis.title.y = element_text(size = 20))+
             theme(axis.title.x = element_text(size = 20))+
             theme(axis.text.x = element_text(size = 20))+
             theme(axis.text.y = element_text(size = 20))+
             theme(strip.text  = element_text(size = 20)) +
             theme(plot.title = element_text(size = 25,hjust = 0.5)) +
             theme(legend.position="none",legend.direction="vertical"),
           
           ggplot(EdgeStatMuscleimat, aes(x=EdgeWeight_addingCorcoffs, fill = Edges)) +
             geom_histogram(color = "gray", position = "stack") +
             scale_fill_manual(values = c('positive' = 'red',
                                          'negative' = 'blue'))+
             scale_y_continuous(breaks = c(0,20,40,60,80))+
             theme_bw()+
             ggtitle("Edge weight") +
             xlab("Sum of Cor.coffs")+
             ylab("Count")+
             labs(tag = "C") +
             theme(plot.title = element_text(size = 25,hjust = 0.5))+
             theme(axis.title.y = element_text(size = 20))+
             theme(axis.title.x = element_text(size = 20))+
             theme(axis.text = element_text(size = 20))+
             theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
             theme(strip.background = element_blank(), strip.placement = "outside") +
             theme(strip.text  = element_text(size = 20)) +         
             theme(legend.position="left",legend.direction="vertical"),
           
           ggplot(EdgeStatMuscleimat, aes(x=Edges, fill = Edges)) +
             geom_bar(position = "stack",
                      #fill = "gray",
                      size=0.7) +
             ggtitle("Number of Edges") +
             theme_bw()+
             ylab("")+
             xlab("Edges")+
             labs(tag = "D") +
             scale_fill_manual(values = c('positive' = 'red',
                                          'negative' = 'blue',alpha = 4))+
             scale_y_continuous(breaks = c(0,100,200,300,400,500,600,800))+
             theme(plot.title = element_text(size = 25,hjust = 0.5))+
             theme(axis.title.y = element_text(size = 20))+
             theme(axis.title.x = element_text(size = 20))+
             theme(axis.text.x = element_text(size = 20))+
             theme(axis.text.y = element_text(size = 20))+
             theme(strip.text  = element_text(size = 20)) +
             theme(legend.position="none",legend.direction="vertical"),
           ncol = 2,
           rel_widths = c(1,0.5,1,0.5))+
         patchwork::plot_annotation(title = "Histograms for gene expression-GIR correlation, edge weight,\n and bar plot for muscle-IMAT network size",
                                    subtitle = "") &  
         theme(plot.title = element_text(hjust = 0.5,size=25,face="bold")) &
         theme(plot.subtitle = element_text(hjust = 0.5,size=15,face="bold"))
)


# How many of the genes are shared in both network 
library("ggvenn")
ggsave("~/Desktop/PhD/finalIMATmuscle/percentile10/venndigram.pdf",
       width = 15, height = 15, units = "cm",
       ggvenn(list(
         "IMAT to muscle" = nodeStatimatMuscle$gene,
         "muscle to IMAT" = nodeStatMuscleimat$gene),
         fill_color = c(#"blue", "yellow", 
           "red", "green"),
         text_size = 6,
         set_name_size = 6) +
         theme(plot.title = element_text(hjust = 0.5,size=25)) +
         labs(title = "Number genes in both networks", tag = "A") 
)



#
## use upset plot to see number genes unique and shared 
upset = list(
  "IMAT_sender" = nodeStatimatMuscle %>% dplyr::filter(Tissue == "IMAT") %>% dplyr::pull(gene),
  "muscle_receiver" = nodeStatimatMuscle %>% dplyr::filter(Tissue == "muscle") %>% dplyr::pull(gene),
  "IMAT_receiver" = nodeStatMuscleimat %>% dplyr::filter(Tissue == "IMAT") %>% dplyr::pull(gene),
  "muscle_sender" = nodeStatMuscleimat %>% dplyr::filter(Tissue == "muscle") %>% dplyr::pull(gene))
UpSet(make_comb_mat(upset))

##  annotations and row/column orders
dev.off()
upset = make_comb_mat(upset)
ss = set_size(upset)
cs = comb_size(upset)
ht = UpSet(upset,
           column_title = "Number of unique and shared\ngenes in both networks",
           column_title_gp = gpar(fontsize = 25),
           set_order = order(ss),
           comb_order = order(comb_degree(upset), -cs),
           top_annotation = HeatmapAnnotation(
             "Unique and shared" = anno_barplot(cs,  
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(6, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Number of genes" = anno_barplot(-ss, 
                                              baseline = 0,
                                              axis_param = list(
                                                at = c(0, -100, -200, -300,400),
                                                labels = c(0, 100, 200, 300,400),
                                                labels_rot = 0),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              width = unit(6, "cm")
             ),
             set_name = anno_text(set_name(upset), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(upset)) + unit(6, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)


ht = draw(ht)
od = column_order(ht)
decorate_annotation("Unique and shared", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 15, col = "#404040"), rot = 45)
})


# ggsave(
#   "~/Desktop/PhD/finalIMATmuscle/Secpercentile10/main_fig1_venndigram.pdf",
#   cowplot::plot_grid(
#     ggvenn(list(
#       "IMC" = nodeStatimatMuscle$gene,
#       "MIC" = nodeStatMuscleimat$gene),
#       fill_color = c(
#         "red", "green"),
#       text_size = 6,
#       set_name_size = 6)+ 
#       theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))+
#       labs(title = "# nodes", tag = ""),
#   ggvenn(list(
#     "IMC" = IMCN$Pair.Name,
#     "MIC" = MIC$Pair.Name),
#     fill_color = c(
#       "red", "green"),
#     text_size = 6,
#     set_name_size = 6)+ 
#     theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))+
#     labs(title = "# edges", tag = ""), labels = "G", ncol = 1) +
#     patchwork::plot_annotation(title = "Number of nodes and edges\n in both networks", subtitle = "") &  
#     theme(plot.title = element_text(hjust = 0.5,size=25,face="bold")) &
#     theme(plot.subtitle = element_text(hjust = 0.5,size=15,face="bold")) 
# )

## use upset plot to see number genes unique and shared 
upset = list(
  "IMAT_sender" = nodeStatimatMuscle %>% dplyr::filter(Tissue == "IMAT") %>% dplyr::pull(gene),
  "muscle_receiver" = nodeStatimatMuscle %>% dplyr::filter(Tissue == "muscle") %>% dplyr::pull(gene),
  "IMAT_receiver" = nodeStatMuscleimat %>% dplyr::filter(Tissue == "IMAT") %>% dplyr::pull(gene),
  "muscle_sender" = nodeStatMuscleimat %>% dplyr::filter(Tissue == "muscle") %>% dplyr::pull(gene))
upsetPlot =
  UpSetR::upset(
    UpSetR::fromList(upset),
    order.by = "freq",
    main.bar.color = "black",
    matrix.color = "black",
    point.size = 6,
    line.size = 1,
    text.scale = 2.5,
    shade.color = c("gray"),
    shade.alpha = 0.25,
    matrix.dot.alpha = 1,
    group.by = "degree") %>% ggplotify::as.ggplot() + 
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))+
  labs(title = "Number of genes in both networks", tag = "B")

ggsave(
  "~/Desktop/PhD/finalIMATmuscle/Secpercentile10/upsetplot.pdf",
  cowplot::plot_grid(
    upsetPlot
  )
)

## Venn and upset plots
upsetPlot =
  UpSetR::upset(
    UpSetR::fromList(upset),
    order.by = "freq",
    main.bar.color = "black",
    matrix.color = "black",
    point.size = 6,
    line.size = 1,
    text.scale = 2.5,
    shade.color = c("gray"),
    shade.alpha = 0.25,
    number.angles = 30,
    matrix.dot.alpha = 1,
    group.by = "degree") %>% ggplotify::as.ggplot() 


##---------Convert i graph to to KNN graph for clustering----------------------- 
imatTOmuscleNetwork_igraph2 <- readRDS("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/IMATmuscleSR_communication.rds") %>%  
  dplyr::filter(!is.na(Sender_Cor_to_GIR_estimate),
                !is.na(Receiver_Cor_to_GIR_estimate),
                !is.na(EdgeWeight_addingCorcoffs))  %>% 
  dplyr::select(Sender,Receiver,EdgeWeight_addingCorcoffs) 
imatTOmuscleNetwork_KNN <- imatTOmuscleNetwork_igraph2 %>% igraph::graph_from_data_frame(directed = F)
imatTOmuscleNetwork_KNN <- igraph::set_edge_attr(imatTOmuscleNetwork_KNN, "edgeWeight", 
                                                 value = c(imatTOmuscleNetwork_igraph2$EdgeWeight_addingCorcoffs))

# look the edge weights of the graph
E(imatTOmuscleNetwork_KNN)$edgeWeight 

# Convert the weighted graph to a weighted adjacency matrix
imatTOmuscleNetwork_KNN <- as_adjacency_matrix(imatTOmuscleNetwork_KNN, attr = "edgeWeight", sparse = FALSE)

# Select Principal components
IMATmuscle_seurat <- CreateSeuratObject(as.matrix(imatTOmuscleNetwork_KNN))
all_genes <-  colnames(imatTOmuscleNetwork_KNN)

## Scale correlation values
IMATmuscle_seurat <- ScaleData(IMATmuscle_seurat, features = all_genes)

## Run pca
n_resolution = 0.5
IMATmuscle_seurat <- RunPCA(IMATmuscle_seurat, features = colnames(imatTOmuscleNetwork_KNN))
IMATmuscle_seurat <- FindNeighbors(IMATmuscle_seurat,reduction = "pca", dims = 1:10)
IMATmuscle_seurat <- FindClusters(IMATmuscle_seurat,resolution = 0.5)
IMATmuscle_seurat <- RunUMAP(IMATmuscle_seurat,reduction = "pca", 
                             dims = 1:10, n.neighbors = 10)

# Number of proteins in each modules 
modules <- IMATmuscle_seurat[[paste0('RNA_snn_res.', n_resolution)]] %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    RNA_snn_res.0.5 = dplyr::case_when(
      group =
        RNA_snn_res.0.5 == 0 ~ 1,
      RNA_snn_res.0.5 == 1 ~ 2,
      RNA_snn_res.0.5 == 2 ~ 3,
      RNA_snn_res.0.5 == 3 ~ 4,
      RNA_snn_res.0.5 == 4 ~ 5,
      RNA_snn_res.0.5 == 5 ~ 6,
      .default = as.numeric(RNA_snn_res.0.5)
    )
  )

# count how many proteins in each modules  
write.csv(modules, "~/Desktop/PhD/finalIMATmuscle/Secpercentile10/KNN_imatTOmuscle.csv")
p_res <- as.data.frame(table(modules))
ggarrange(
  (ElbowPlot(IMATmuscle_seurat) +
     theme_minimal()+
     theme(axis.title=element_text(size= 15),
           axis.text = element_text(size = 13),
           plot.title = element_text(hjust = 0.5,size=20),
           legend.title = element_text(color = "Black", size = 16),
           legend.text=element_text(color = "black", size = 16),
           legend.position = "right",
           legend.direction = "vertical",
           axis.line = element_line(size = 1)) +
     labs(title = "Screen plot") +
     labs(tag = "A")+
     xlab("PCAs") +
     ylab("% of variance"))+
    (DimPlot(IMATmuscle_seurat, reduction = "umap") + 
       theme_dark()+
       ggtitle("KNN-clustering") +
       labs(tag = "B")+
       theme(axis.title=element_text(size= 15),
             axis.text = element_text(size = 13),
             plot.title = element_text(hjust = 0.5,size=20),
             legend.title = element_text(color = "Black", size = 16),
             legend.text=element_text(color = "black", size = 16),
             legend.position = "right",
             legend.direction = "vertical"))+
    (ggplot(p_res ,aes(x = RNA_snn_res.0.5, y = Freq)) +
       geom_bar(stat = 'identity', fill = "steelblue") +
       scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
       labs(x = 'modules', y ='# genes', tag = "C") +
       theme_minimal() +
       theme(axis.title=element_text(size= 15),
             axis.text = element_text(size = 13),
             plot.title = element_text(hjust = 0.5,size=20),
             legend.title = element_text(color = "Black", size = 16),
             legend.text=element_text(color = "black", size = 16),
             legend.position = "bottom",
             legend.direction = "horizontal",
             axis.line = element_line(size = 1)) +
       ggtitle("Number of genes\n in the modules"))) +
  patchwork::plot_annotation(title = "IMAT to muscle network") &  
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))


## module to biology 
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
hallmark_list <- read.gmt('/Users/amare.wolide/ProtomicsNSCLC/h.all.v2023.1.Hs.symbols.gmt')
modules = read.csv("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/KNN_imatTOmuscle.csv")
idx_sorted <- sort(unique(modules[,2]))

# do the enrichment using clusterprofiler 
module_enrichment <- do.call(rbind, lapply(idx_sorted, function(i) {
  print(i)
  res <- clusterProfiler::enricher(gene = modules[,1][modules[,2] %in%  i], 
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   qvalueCutoff = 0.05,
                                   gson = NULL,
                                   TERM2GENE = hallmark_list,
                                   TERM2NAME = NA
  )
  res <- as.data.frame(res)
  if(nrow(res) == 0) {res[1, ] = NA} else {
    res$Module <- i
    return(res)}
}))

module_enrichment <- module_enrichment[order(module_enrichment$Count, decreasing = TRUE), ]
module_enrichment <- module_enrichment[complete.cases(module_enrichment), ]

# Remove prefix "hallmark_", convert to lowercase, and replace underscores with spaces while keeping first letter capitalized
module_enrichment$Description <- str_replace(module_enrichment$Description, "^HALLMARK_", "") # Remove prefix
module_enrichment$Description <- str_replace_all(module_enrichment$Description, "_", " ") # Replace underscores with spaces
module_enrichment <- module_enrichment %>% dplyr::select(-ID)
write.csv(module_enrichment,"~/Desktop/PhD/finalIMATmuscle/Secpercentile10/module_enrichmentimatTOmuscle.csv")

# Function to create bar plot for each cluster
create_bar_plot <- function(df, cluster_num) {
  # Subset data for the current cluster
  df_cluster <- subset(df, Module == cluster_num)
  
  # Create -log10(P-value)
  df_cluster$neg_log10_Pvalue <- -log10(df_cluster$p.adjust)
  
  # Sort data frame by p-value in descending order
  df_cluster <- df_cluster %>% arrange(desc(p.adjust))
  
  # Create bar plot using ggplot2
  p <- ggplot(df_cluster, aes(x = reorder(Description, neg_log10_Pvalue), y = neg_log10_Pvalue)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Module ", cluster_num),
         x = "",
         y = "-log10(p.adjust)") +
    theme_minimal() + coord_flip() +
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 13))+
    theme(plot.title = element_text(size = 15))+
    theme(axis.text.y  = element_text(size = 10))+
    theme(axis.text.x  = element_text(size = 12))+
    theme(panel.grid.major = element_line(colour = "gray",linetype = 2))+
    theme(strip.background = element_blank(), strip.placement = "outside") 
  #theme(strip.text  = element_text(size = 15))
  return(p)
}

plot_list <- list()

# Loop through each cluster and create bar plot
for (cluster_num in unique(module_enrichment$Module)) {
  plot_list[[as.character(cluster_num)]] <- create_bar_plot(module_enrichment, cluster_num)
}

# Reorder the plot_list based on module number
plot_list <- plot_list[order(as.numeric(names(plot_list)))]

# Print the plots
for (cluster_num in unique(module_enrichment$Module)) {
  print(plot_list[[as.character(cluster_num)]])
}

# Combine all the plots into one figure
combined_plot <- wrap_plots(plotlist = plot_list, nrow = 3) 

# Print the combined plot
print(combined_plot)

##------------------------------------------------------------------------------

muscleTOimatNetwork_igraph2 <- readRDS("~/Desktop/PhD/finalIMATmuscle/Secpercentile10/muscleIMATSR_communication.rds") %>%  
  dplyr::filter(!is.na(Sender_Cor_to_GIR_estimate),
                !is.na(Receiver_Cor_to_GIR_estimate),
                !is.na(EdgeWeight_addingCorcoffs))  %>% 
  dplyr::select(Sender,Receiver,EdgeWeight_addingCorcoffs) 
muscleTOimatNetwork_KNN <- muscleTOimatNetwork_igraph2 %>% igraph::graph_from_data_frame(directed = F)
muscleTOimatNetwork_KNN <- igraph::set_edge_attr(muscleTOimatNetwork_KNN, "edgeWeight", 
                                                 value = c(muscleTOimatNetwork_igraph2$EdgeWeight_addingCorcoffs))

# look the edge weights of the graph
E(muscleTOimatNetwork_KNN)$edgeWeight 

# Convert the weighted graph to a weighted adjacency matrix
muscleTOimatNetwork_KNN <- as_adjacency_matrix(muscleTOimatNetwork_KNN, attr = "edgeWeight", sparse = FALSE)

# Select Principal components
muscleIMAT_seurat <- CreateSeuratObject(as.matrix(muscleTOimatNetwork_KNN))
all_genes <-  colnames(muscleTOimatNetwork_KNN)

## Scale correlation values
muscleIMAT_seurat <- ScaleData(muscleIMAT_seurat, features = all_genes)

## Run pca
muscleIMAT_seurat <- RunPCA(muscleIMAT_seurat, features = colnames(muscleTOimatNetwork_KNN))
muscleIMAT_seurat <- FindNeighbors(muscleIMAT_seurat,reduction = "pca", dims = 1:10)
muscleIMAT_seurat <- FindClusters(muscleIMAT_seurat,resolution = 0.5)
muscleIMAT_seurat <- RunUMAP(muscleIMAT_seurat,reduction = "pca", 
                             dims = 1:10, n.neighbors = 10)

# Number of proteins in each modules 
modules2 <- muscleIMAT_seurat[[paste0('RNA_snn_res.', n_resolution)]] %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    RNA_snn_res.0.5 = dplyr::case_when(
      group =
        RNA_snn_res.0.5 == 0 ~ 1,
      RNA_snn_res.0.5 == 1 ~ 2,
      RNA_snn_res.0.5 == 2 ~ 3,
      RNA_snn_res.0.5 == 3 ~ 4,
      RNA_snn_res.0.5 == 4 ~ 5,
      RNA_snn_res.0.5 == 5 ~ 6,
      .default = as.numeric(RNA_snn_res.0.5)
    )
  )

# count how many proteins in each modules  
write.csv(modules2, "~/Desktop/PhD/finalIMATmuscle/Secpercentile10/KNN_muscleTOimat.csv")
p_res2 <- as.data.frame(table(modules2))
ggarrange(
  (ElbowPlot(muscleIMAT_seurat) +
     theme_minimal()+
     theme(axis.title=element_text(size= 15),
           axis.text = element_text(size = 13),
           plot.title = element_text(hjust = 0.5,size=20),
           legend.title = element_text(color = "Black", size = 16),
           legend.text=element_text(color = "black", size = 16),
           legend.position = "right",
           legend.direction = "vertical",
           axis.line = element_line(size = 1)) +
     labs(title = "Screen plot") +
     labs(tag = "A")+
     xlab("PCAs") +
     ylab("% of variance"))+
    (DimPlot(muscleIMAT_seurat, reduction = "umap") + 
       theme_dark()+
       ggtitle("KNN-clustering") +
       labs(tag = "B")+
       theme(axis.title=element_text(size= 15),
             axis.text = element_text(size = 13),
             plot.title = element_text(hjust = 0.5,size=20),
             legend.title = element_text(color = "Black", size = 16),
             legend.text=element_text(color = "black", size = 16),
             legend.position = "right",
             legend.direction = "vertical"))+
    (ggplot(p_res2 ,aes(x = RNA_snn_res.0.5, y = Freq)) +
       geom_bar(stat = 'identity', fill = "steelblue") +
       scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
       labs(x = 'modules', y ='# genes', tag = "C") +
       theme_minimal() +
       theme(axis.title=element_text(size= 15),
             axis.text = element_text(size = 13),
             plot.title = element_text(hjust = 0.5,size=20),
             legend.title = element_text(color = "Black", size = 16),
             legend.text=element_text(color = "black", size = 16),
             legend.position = "bottom",
             legend.direction = "horizontal",
             axis.line = element_line(size = 1)) +
       ggtitle("Number of genes\n in the modules"))) +
  patchwork::plot_annotation(title = "muscle to IMAT network") &  
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))

hallmark_list <- read.gmt('/Users/amare.wolide/ProtomicsNSCLC/h.all.v2023.1.Hs.symbols.gmt')
idx_sorted2 <- sort(unique(modules2[,2]))

# do the enrichment using clusterprofiler 
module_enrichment2 <- do.call(rbind, lapply(idx_sorted2, function(i) {
  print(i)
  res <- clusterProfiler::enricher(gene = modules2[,1][modules2[,2] %in%  i], 
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   qvalueCutoff = 0.05,
                                   gson = NULL,
                                   TERM2GENE = hallmark_list,
                                   TERM2NAME = NA
  )
  res <- as.data.frame(res)
  if(nrow(res) == 0) {res[1, ] = NA} else {
    res$Module <- i
    return(res)}
}))

module_enrichment2 <- module_enrichment2[order(module_enrichment2$Count, decreasing = TRUE), ]
module_enrichment2 <- module_enrichment2[complete.cases(module_enrichment2), ]

# Remove prefix "hallmark_", convert to lowercase, and replace underscores with spaces while keeping first letter capitalized
module_enrichment2$Description <- str_replace(module_enrichment2$Description, "^HALLMARK_", "") # Remove prefix
module_enrichment2$Description <- str_replace_all(module_enrichment2$Description, "_", " ") # Replace underscores with spaces
module_enrichment2 <- module_enrichment2 %>% dplyr::select(-ID)
write.csv(module_enrichment2,"~/Desktop/PhD/finalIMATmuscle/Secpercentile10/module_enrichmentmuscleTOimat.csv")

plot_list2 <- list()

# Loop through each cluster and create bar plot
for (cluster_num in unique(module_enrichment2$Module)) {
  plot_list2[[as.character(cluster_num)]] <- create_bar_plot(module_enrichment2, cluster_num)
}

# Reorder the plot_list based on module number
plot_list2 <- plot_list2[order(as.numeric(names(plot_list2)))]

# Print the plots
for (cluster_num in unique(module_enrichment2$Module)) {
  print(plot_list2[[as.character(cluster_num)]])
}

# Combine all the plots into one figure
combined_plot2 <- wrap_plots(plotlist = plot_list2, nrow = 3) 

# Print the combined plot
print(combined_plot2)
##----------

# enrichment bar plots
library(enrichR)
library(clusterProfiler)
enrichment_barplot <- function(data,
                               labs_title = NULL, 
                               labs_subtitle = NULL,
                               legend_position = NULL) {
  data_reorder <- data %>% 
    mutate(Description = fct_reorder(Description, -log10(p.adjust)))
  
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
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text = element_text(color = "Black", size = 10, face = "bold")) +
    theme(legend.position = legend_position)+
    labs(title = labs_title, subtitle = labs_subtitle) +
    labs(x = "",
         y = bquote("" ~ -Log[10] ~ "p.adjust" ~ ""))
  
  return(p)
}
