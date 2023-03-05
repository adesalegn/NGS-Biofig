##---------------Sender-receiver databases--------------------------------------
rm(list = ls())
gc()
plot.new()
all_sender_receiver_human <-
  read.csv("~/revision_newanalysis/new_sender_receiver_ppi_DB.csv") %>%
  dplyr::select(Sender_gene_Class, Receiver_gene_Class)

table(all_sender_receiver_human$Sender_gene_Class)
table(all_sender_receiver_human$Receiver_gene_Class)
plot.new()
grid.col = structure(c(
  RColorBrewer::brewer.pal(8, "Paired"),
  RColorBrewer::brewer.pal(4, "Dark2")
),
names = unique(
  c(
    all_sender_receiver_human$Sender_gene_Class,
    all_sender_receiver_human$Receiver_gene_Class
  )
))

# Base plot
circos.par(start.degree = 270) # vertically symmetric
chordDiagram(
  all_sender_receiver_human,
  grid.col = grid.col,
  transparency = 0.5,
  annotationTrack = c("grid"),
  directional = 1,
  #direction.type = "arrows",
  diffHeight  = -0.04,
  #link.arr.type = "triangle",
  link.arr.lty = "solid",
  link.arr.lwd = 1,
  link.arr.width =  0.2,
  link.sort = TRUE,
  #link.arr.col = "black",
  #col = "#00000000",
  link.largest.ontop = TRUE,
  #link.largest.ontop = FALSE,
  preAllocateTracks = list(
    list(track.height = uh(5, "mm")),
    list(track.height = uh(5, "mm")),
    list(track.height = uh(5, "mm"))
  )
)


## Tracking the sectors
circos.track(
  track.index = 3,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      mean(xlim),
      mean(ylim),
      col = "black",
      sector.index,
      cex = 2,
      adj = c(1, 0),
      facing = "clockwise",
      niceFacing = TRUE
    )
  },
  track.height = 0.25,
  bg.border = NA
)

##highlighting sectors of sender
highlight.sector(
  ##growth factors, activin and inhibin,Bone morphogenetic proteins,hormones, Neuropeptides, cytokines,chemokines,and interlukines as ligand
  c("Neuropeptides", "Secreted_Glycoprotein"),
  track.index = 2,
  col = "black",
  #text = "Endocrine", cex = 2,  facing = "clockwise", text.col = "green",
  niceFacing = TRUE
)

highlight.sector(
  c(
    "Surface_Kinases",
    "Secreted_Proteins",
    "Cytokines",
    "Growth_Factors"
  ),
  track.index = 2,
  col = "violet",
  #text = "Paracrine",cex = 2, facing = "clockwise",text.col = "black",
  niceFacing = TRUE
)

highlight.sector(
  c("Junctional"),
  track.index = 2,
  col = "blue",
  #text = "Juxtacrine",cex = 2, facing = "clockwise",text.col = "black",
  niceFacing = TRUE
)

highlight.sector(
  c("Predicted_Ligand"),
  track.index = 2,
  col = "#1EFF00",
  #text = "Predicted_Ligand",cex = 2, facing = "clockwise",text.col = "black",
  niceFacing = TRUE
)

## Tracking all sender information
highlight.sector(
  c(
    "Cytokines",
    "Growth_Factors",
    "Junctional",
    "Neuropeptides",
    "Predicted_Ligand",
    "Secreted_Glycoprotein",
    "Secreted_Proteins",
    "Surface_Kinases"
  ),
  track.index = 1,
  col = "darkorange",
  text.vjust = "0.8cm",
  text = "Sender",
  cex = 2.5,
  facing = "downward",
  text.col = "black",
  niceFacing = TRUE
)

## highlighting all sectors of receiver
highlight.sector(
  c(
    "Adhesion_Protein",
    "Catalytic_Receptors",
    "GPCRs",
    "Transport_Protein"
  ),
  track.index = 1,
  col = "darkmagenta",
  text.vjust = "0.8cm",
  text = "Receiver",
  cex = 2.5,
  facing = "downward",
  text.col = "black",
  niceFacing = TRUE
)
abline(v = 0,
       lty = 2,
       col = "#00000080",
       lwd = 5)
title("Sender-Receiver DB", cex.main=2)
my_base <- recordPlot()
cowplot::plot_grid(my_base,
                   rel_widths = 1,
                   rel_heights = 1)

text_tbl <- data.frame(
  Gene = c("Sender", "Receiver"),
  Features = c(
    "Senders are protein coding genes known in signaling function.
     Summarized in into three superclass such as endocrine(Neuropeptides and Secreted Glycoprotein),
     paracrine(Cytokines,Growth Factors,Secreted Proteins, and Surface Kinases)
     Juxtcrine(Junctional Proteins).",
    "Receivers are protein coding genes involved in responding for the incoming signal.
    Summarized in four class such as Adhesion Protein,Catalytic Receptors,GPCRs, and Transporter Proteins."
  )
)

library(kableExtra)
kbl(text_tbl) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "40em", background = "yellow")


##--------------------------------------------------------------------------------------------------
plot.new()
grid.col = structure(c(
  RColorBrewer::brewer.pal(8, "Paired"),
  RColorBrewer::brewer.pal(4, "Dark2")
),
names = unique(
  c(
    all_sender_receiver_human$Sender_gene_Class,
    all_sender_receiver_human$Receiver_gene_Class
  )
))

# Base plot
circos.par(start.degree = 270) # vertically symmetric
chordDiagram(
  all_sender_receiver_human,
  grid.col = grid.col,
  transparency = 0.5,
  annotationTrack = c("grid"),
  directional = 1,
  diffHeight  = -0.04,
  link.arr.lty = "solid",
  link.arr.lwd = 1,
  link.arr.width =  0.2,
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  preAllocateTracks = list(
    list(track.height = uh(5, "mm")),
    list(track.height = uh(5, "mm")),
    list(track.height = uh(5, "mm"))
  )
)

##highlighting sectors of sender
highlight.sector(
  c("Neuropeptides", "Secreted_Glycoprotein"),
  track.index = 2,
  col = "black",
  niceFacing = TRUE
)

highlight.sector(
  c(
    "Surface_Kinases",
    "Secreted_Proteins",
    "Cytokines",
    "Growth_Factors"
  ),
  track.index = 2,
  col = "violet",
  niceFacing = TRUE
)

highlight.sector(
  c("Junctional"),
  track.index = 2,
  col = "blue",
  niceFacing = TRUE
)

highlight.sector(
  c("Predicted_Ligand"),
  track.index = 2,
  col = "#1EFF00",
  niceFacing = TRUE
)

## Tracking all sender information
highlight.sector(
  c(
    "Cytokines",
    "Growth_Factors",
    "Junctional",
    "Neuropeptides",
    "Predicted_Ligand",
    "Secreted_Glycoprotein",
    "Secreted_Proteins",
    "Surface_Kinases"
  ),
  track.index = 1,
  col = "darkorange",
  text.vjust = "0.8cm",
  text = "",
  cex = 2.5,
  facing = "downward",
  text.col = "black",
  niceFacing = TRUE
)

## highlighting all sectors of receiver
highlight.sector(
  c(
    "Adhesion_Protein",
    "Catalytic_Receptors",
    "GPCRs",
    "Transport_Protein"
  ),
  track.index = 1,
  col = "darkmagenta",
  text.vjust = "0.8cm",
  text = "",
  cex = 2.5,
  facing = "downward",
  text.col = "black",
  niceFacing = TRUE
)
abline(v = 0,
       lty = 2,
       col = "#00000080",
       lwd = 5)
plot.new()
table(all_sender_receiver_human$Sender_gene_Class)
#Paracrine(840 + 1663 + 3217 + 557) = 6277
#Juxtacrine = 117
#Endocrine(648 + 893) = 1541
lgd_points_superclass <- Legend(
  at = c(
    "Endocrine(n=1541)",
    "Paracrine(n=6277)",
    "Juxtacrine(n=117)",
    "Predicted_Ligand(n=1942)"
  ),
  type = "points",
  nrow = 4,
  ncol = 1,
  legend_gp = gpar(
    col = c("black", "violet", "blue", "#1EFF00"),
    lwd = 5
  ),
  size =  unit(10, "mm"),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 24, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 30),
  title = "Sender gene\n Super-Class"
)

draw(lgd_points_superclass, just = "centre")

##--------------------------------------------------------------------------------------------------
plot.new()
lgd_points_sr = Legend(
  at = c("Sender", "Receiver"),
  type = "points",
  legend_gp = gpar(col = c("darkorange", "darkmagenta"), lwd = 5),
  size =  unit(10, "mm"),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 24, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 30),
  title = "Cell"
)

lgd_points_superclass <- Legend(
  at = c(
    "Endocrine(n=1541)",
    "Paracrine(n=6277)",
    "Juxtacrine(n=117)",
    "Predicted_Ligand(n=1942)"
  ),
  type = "points",
  nrow = 4,
  ncol = 1,
  legend_gp = gpar(
    col = c("black", "violet", "blue", "#1EFF00"),
    lwd = 5
  ),
  size =  unit(10, "mm"),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 24, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 30),
  title = "Sender gene\n Super-Class"
)

lgd_points_Sender = Legend(
  at = c(
    "Cytokines",
    "Growth_Factors",
    "Junctional",
    "Neuropeptides",
    "Predicted_Ligand",
    "Secreted_Glycoprotein",
    "Secreted_Proteins",
    "Surface_Kinases"
  ),
  type = "points",
  legend_gp = gpar(col = RColorBrewer::brewer.pal(8, "Paired"), lwd = 5),
  size =  unit(10, "mm"),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 24, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 30),
  title = "Sender gene\n Gene Class"
)


lgd_points_Receiver = Legend(
  at =  c(
    "Adhesion_Protein",
    "Catalytic_Receptors",
    "GPCRs",
    "Transport_Protein"
  ),
  type = "points",
  legend_gp = gpar(col = RColorBrewer::brewer.pal(4, "Dark2"), lwd = 5),
  size =  unit(10, "mm"),
  title_position = "topcenter",
  title_gp = gpar(fontsize = 24, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 30),
  title = "Receiver gene\n Gene Class"
)

plot.new()
lgd_list_vertical = packLegend(
  direction = "vertical",
  lgd_points_sr,
  lgd_points_superclass,
  lgd_points_Sender,
  lgd_points_Receiver
)

draw(lgd_list_vertical, just = "centre")
plot.new()
draw(lgd_points_superclass, just = "centre")





##---------data visualization for the correlation
## Visualization---using cytoscape
## Visualization---Circus plot
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(magick)
library(GetoptLong)
library(reshape2)
library(ComplexHeatmap)
library(magrittr)
library(grid)
library(gridExtra)
library(gridBase)
library(randomcoloR)

# plot.new()
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
    c(structure(data$IMAT_corGIR, names = data$cell_Sender))
  genes_2 <- genes_2[!is.na(genes_2)]
  genes_2 <- genes_2[!duplicated(paste(names(genes_2), genes_2))]
  genes_2 <- genes_2[order(names(genes_2))]
  
  ## IMAT/senders/ correlated to FG
  genes_3 <- c(structure(data$IMAT_corFG, names = data$cell_Sender))
  genes_3 <- genes_3[!is.na(genes_3)]
  genes_3 <- genes_3[!duplicated(paste(names(genes_3), genes_3))]
  genes_3 <- genes_3[order(names(genes_3))]
  
  ## muscle/receiver/ correlated to GIR
  genes_4 <-
    c(structure(data$muscle_corGIR, names = data$cell_Receiver))
  genes_4 <- genes_4[!is.na(genes_4)]
  genes_4 <- genes_4[!duplicated(paste(names(genes_4), genes_4))]
  genes_4 <- genes_4[order(names(genes_4))]
  
  ## muscle/receiver/ correlated to FG
  genes_5 <-
    c(structure(data$muscle_corFG, names = data$cell_Receiver))
  genes_5 <- genes_5[!is.na(genes_5)]
  genes_5 <- genes_5[!duplicated(paste(names(genes_5), genes_5))]
  genes_5 <- genes_5[order(names(genes_5))]
  
  if (is.null(link.arr.lty)) {
    link.arr.lty = "solid"
  }
  if (is.null(link.arr.col)) {
    data <-
      data %>% mutate(link_col = ifelse(data$pcc > 0, "#B2182B", "#2166ac"))
  }
  else {
    data$link_col = link.arr.col
  }
  if (is.null(link.arr.type)) {
    link.arr.type = "triangle"
  }
  if (is.null(gene_col)) {
    geneSendr_col <- structure(
      c(RColorBrewer::brewer.pal(8, "Paired")),
      names = c(
        "Cytokines",
        "Growth_Factors",
        "Junctional",
        "Neuropeptides",
        "Predicted_Ligand",
        "Secreted_Glycoprotein",
        "Secreted_Proteins",
        "Surface_Kinases"
      )
    )
    
    geneReciver_col <-
      structure(
        c(RColorBrewer::brewer.pal(4, "Dark2")),
        names = c(
          "Adhesion_Protein",
          "Catalytic_Receptors",
          "GPCRs",
          "Transport_Protein"
        )
      )
    gene_col <-
      structure(c(geneSendr_col[data$Sender_gene_Class], geneReciver_col[data$Receiver_gene_Class]),
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
            sector.index %in% unlist(genes_3),
            'blue',
            ifelse(
              sector.index %in% unlist(genes_4),
              'red',
              ifelse(sector.index %in% unlist(genes_5), 'blue', 'black')
            )
          )
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


## LC
plot.new()
circos.par(start.degree = 270) # vertically symmetric
iLC_IMAT_muscle <-
  read.csv(
    "~/revision_newanalysis/LC/LC005_IMAT_to_muscle_communication_posative.csv",
    row.names = 1
  ) %>%
  dplyr::arrange(Sender_gene_Class) %>%
  dplyr::mutate(
    cell_Receiver = "muscle",
    cell_Sender = "IMAT",
    Sender = IMAT,
    Receiver = muscle
  )
tissue_imat_muscle <-
  structure(c('darkorange', 'darkmagenta'), names = unique(
    c(iLC_IMAT_muscle$cell_Sender, iLC_IMAT_muscle$cell_Receiver)
  ))
cirlc = plotcircosinR_from_to_cor(iLC_IMAT_muscle, cell_col = tissue_imat_muscle, print.cell = T)
text(
  x = 0,
  y = 0,
  "IMAT to Muscle(LC)",
  col = "blue",
  cex = 2
)

## OW
plot.new()
circos.par(start.degree = 270)
iOW_IMAT_muscle <-
  read.csv(
    "~/revision_newanalysis/OW/OW005_IMAT_to_muscle_communication_posative.csv",
    row.names = 1
  ) %>%
  dplyr::arrange(Sender_gene_Class) %>%
  dplyr::mutate(
    cell_Receiver = "muscle",
    cell_Sender = "IMAT",
    Sender = IMAT,
    Receiver = muscle
  )
tissue_imat_muscle <-
  structure(c('darkorange', 'darkmagenta'), names = unique(
    c(iOW_IMAT_muscle$cell_Sender, iOW_IMAT_muscle$cell_Receiver)
  ))

cirow = plotcircosinR_from_to_cor(iOW_IMAT_muscle, cell_col = tissue_imat_muscle, print.cell = T)
text(
  x = 0,
  y = 0,
  "IMAT to Muscle(OW)",
  col = "blue",
  cex = 2
)
## T2D
plot.new()
circos.par(start.degree = 270)
iT2D_IMAT_muscle <-
  read.csv(
    "~/revision_newanalysis/T2D/T2D005_IMAT_to_muscle_communication_posative.csv",
    row.names = 1
  ) %>%
  dplyr::arrange(Sender_gene_Class) %>%
  dplyr::mutate(
    cell_Receiver = "muscle",
    cell_Sender = "IMAT",
    Sender = IMAT,
    Receiver = muscle
  )
tissue_imat_muscle <-
  structure(c('darkorange', 'darkmagenta'), names = unique(
    c(
      iT2D_IMAT_muscle$cell_Sender,
      iT2D_IMAT_muscle$cell_Receiver
    )
  ))
cirt2d = plotcircosinR_from_to_cor(iT2D_IMAT_muscle, cell_col = tissue_imat_muscle, print.cell = T)
text(
  x = 0,
  y = 0,
  "IMAT to Muscle(T2D)",
  col = "blue",
  cex = 2
)

##Legend
#x = circle_size,
plot.new()
lgd_points_tissue = Legend(
  at = c("IMAT", "muscle"),
  type = "points",
  legend_gp = gpar(col = c("darkorange", "darkmagenta"), lwd = 5),
  size =  unit(7, "mm"),
  title_position = "topleft",
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 20),
  title = "Tissue"
)

lgd_points_Sender = Legend(
  at = c(
    "Cytokines",
    "Growth_Factors",
    "Junctional",
    "Neuropeptides",
    "Predicted_Ligand",
    "Secreted_Glycoprotein",
    "Secreted_Proteins",
    "Surface_Kinases"
  ),
  type = "points",
  legend_gp = gpar(col = RColorBrewer::brewer.pal(8, "Paired")),
  size =  unit(7, "mm"),
  title_position = "topleft",
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 20),
  title = "Sender Protein \ncoding Gene Class"
)


lgd_points_Receiver = Legend(
  at =  c(
    "Adhesion_Protein",
    "Catalytic_Receptors",
    "GPCRs",
    "Transport_Protein"
  ),
  type = "points",
  legend_gp = gpar(col = RColorBrewer::brewer.pal(4, "Dark2")),
  size =  unit(7, "mm"),
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 20),
  title_position = "topleft",
  title = "Receiver Protein \ncoding Gene Class"
)


lgd_lines_correlations = Legend(
  at = c("+pcc", "-pcc"),
  type = "lines",
  legend_gp = gpar(col = c("#B2182B", "#2166ac"), lwd = 2),
  size =  unit(7, "mm"),
  title_position = "topleft",
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 20),
  title = "Correlation"
)


lgd_lines_IS = Legend(
  at = c("GIR", "FG", "None"),
  type = "boxplot",
  legend_gp = gpar(col = c("red", "blue", "black"), lwd = 2),
  size =  unit(7, "mm"),
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  title_gap = unit(4, "mm"),
  grid_height = unit(7, "mm"),
  grid_width = unit(5, "mm"),
  column_gap = unit(4, "mm"),
  row_gap = unit(2, "mm"),
  labels_gp = gpar(fontsize = 20),
  title_position = "topleft",
  title = "Genes associated with\nInsulin Sensitivity (IS)"
)

plot.new()
lgd_list_vertical = packLegend(
  direction = "horizontal",
  lgd_points_tissue,
  lgd_points_Sender,
  lgd_points_Receiver,
  lgd_lines_correlations,
  lgd_lines_IS
)
draw(lgd_list_vertical, just = "centre")
