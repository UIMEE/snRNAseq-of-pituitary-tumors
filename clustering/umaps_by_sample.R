library("ggplot2")
library("vroom")
library("readr")
#Creando los UMAPS + cluster individuales + feature plots
scRNA <- lapply(
  scRNA, 
  function(x){
    x <- RunUMAP(x, dims = 1:50)
    #find cluster
    x <- FindNeighbors(x, dims = 1:30, verbose = FALSE)
    x <- FindClusters(x, verbose = FALSE)
    return(x)
  }
)

#Finally some visualization
DimPlot(scRNA$AC_001, group.by = "orig.ident")


DimPlot(scRNA$AC_001, label = TRUE) + NoLegend()

lapply(
  names(scRNA), 
  function(x){
    p <- DimPlot(scRNA[[x]], label = TRUE) + NoLegend()
    ggsave(
      p,
      file = paste0("results/", "UMAP_", x, ".pdf"),
      height = 6,
      width  = 7
    )
  }
)
#TO DO panel with experimental group by samples

#Identificando por Seurat
markers <- lapply(
  scRNA,
  FindAllMarkers
)

#guardando a archivo
lapply(
  names(markers), 
  function(x){
    y <- markers[[x]]
    write_csv(
      y,
      file = paste0("results/", "markersSimple_", x, ".csv")
    )
  }
)
lapply(markers, dim)

#CP_001
p <- FeaturePlot(
  scRNA$CP_001,
  feature = unique(c(
    "POMC", "CRHR1", "AVPR1B", "TBX19",
    "LGR4", "RBPMS", "WWTR1", "CTNNB1", "JUND", "JUN", "PITX1",
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74", "FCGR2A",
    "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206", "CD14", 
    "PTPRC", "CD45", "CD53", "CD2", "IL7R", "GZMA", "NKG7", "LTB",
    "CD96", "TFRC", "CD71", "GYPA", "AHSP", "CA1", "HBB", "HBA1",
    "HBA2"
  ))
)
ggsave(
  p,
  file = "results/feature_CP_001.pdf",
  height = 20,
  width  = 20 
)

#AC_002
p <- FeaturePlot(
  scRNA$AC_002,
  feature = unique(c(
    "POMC", "TBX19", "CRHR1", "GAL", "CDH7", "KCNH7", "AVPR1B",
    "POMC", "TBX19", "RAB3B", "POU1F1", "PRL", "ALK", "PRLR", "ANGPT1",
    "POMC", "TBX19", "POU1F1", "GH1", "GHRHR",
    "PDGFRB", "CSPG4", "NES", "MCAM", "COL4A1", "ADGRD1", 
    "VWF", "CD34", "ENG", "PECAM1", 
    "LGR4", "RBPMS", "WWTR1", "CTNNB1", "JUND", "JUN",
    "PITX1",
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74",
    "FCGR2A", "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206",
    "CD14", "PTPRC", "CD45", "CD53", "CD2", "IL7R", 
    "GZMA", "NKG7", "LTB", "CD96"
  ))
)
ggsave(
  p,
  file = "results/feature_AC_002.pdf",
  height = 20,
  width  = 20 
)

#AC_001
p <- FeaturePlot(
  scRNA$AC_001,
  feature = unique(c(
    "POU1F1", "GH1", "GHRHR", 
    "PBK", "SPC25", "CDCA3", "BIRC5", "CENPF", "UBE2C",
    "HMMR", "TOP2A", "CKAP2L", "CCNB1", "CDK1", "CCNA2",
    "CCNB2", "PTPRC", "CD45", "CD53", "CD163", 
    "CD86", "CD74", "FCGR2A", "TLR2", "ITGAX", "ITGA4", "MRC1",
    "CD206", "CD14", 
    "PTPRC", "CD45", "CD53", "CD2", "IL7R", "GZMA", "NKG7", "LTB",
    "CD96"
  ))
)
ggsave(
  p,
  file = "results/feature_AC_001.pdf",
  height = 20,
  width  = 20 
)

#CU_002
p <- FeaturePlot(
  scRNA$CU_002,
  feature = unique(c(
    "POU1F1", "GH1", "GHRHR", 
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74", 
    "FCGR2A", "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206",
    "CD14", "PTPRC", "CD45", "CD53", "CD2", "IL7R", "GZMA",
    "NKG7", "LTB", "CD96",
    "VWF", "CD34", "ENG", "PECAM1"
  ))
)
ggsave(
  p,
  file = "results/feature_CU_002.pdf",
  height = 20,
  width  = 20 
)

#NF_003
p <- FeaturePlot(
  scRNA$NF_003,
  feature = unique(c(
    "NR5A1", "LHB", "FSHB", "CGA", "GNRHR", "GATA2", "TGFBR3L",
    "LGR4", "RBPMS", "WWTR1", "CTNNB1", "JUND",
    "JUN", "PITX1", 
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74", "FCGR2A",
    "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206", "CD14"
  ))
)
ggsave(
  p,
  file = "results/feature_NF_003.pdf",
  height = 20,
  width  = 20 
)

#NF_002
p <- FeaturePlot(
  scRNA$NF_002,
  feature = unique(c(
    "NR5A1", "LHB", "FSHB", "CGA", "GNRHR", "GATA2", "TGFBR3L", 
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74", "FCGR2A",
    "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206", "CD14", 
    "PDGFRB", "CSPG4", "NES", "MCAM", "COL4A1", "ADGRD1" 
  ))
)
ggsave(
  p,
  file = "results/feature_NF_002.pdf",
  height = 20,
  width  = 20 
)

#NF_004
p <- FeaturePlot(
  scRNA$NF_004,
  feature = unique(c(
    "NR5A1", "LHB", "FSHB", "CGA", "GNRHR", "GATA2", "TGFBR3L" 
  ))
)
ggsave(
  p,
  file = "results/feature_NF_004.pdf",
  height = 20,
  width  = 20 
)

#SC_001
p <- FeaturePlot(
  scRNA$SC_001,
  feature = unique(c(
    "NR5A1", "LHB", "FSHB", "CGA", "GNRHR", "GATA2", "TGFBR3L",
    "PTPRC", "CD45", "CD53", "CD163", "CD86", "CD74", "FCGR2A",
    "TLR2", "ITGAX", "ITGA4", "MRC1", "CD206", "CD14", 
    "PTPRC", "CD45", "CD53", "CD2", "IL7R", "GZMA", "NKG7",
    "LTB", "CD96"
  ))
)
ggsave(
  p,
  file = "results/feature_SC_001.pdf",
  height = 20,
  width  = 20 
)

save(
  scRNA, 
  file="results/scRNAindividual.RData", 
  compress="xz"
)

#features

#Etiquetando clusters
clusterid <- list(
  # Muestra AC-001-2021
  # Clusters 0, 1, 2, 3, 4, 5, 7= POU1F1 Tumoral cells 
  # Cluster 6 = Macrophages 
  # Cluster 8 = NK cells
  "AC_001" =  c(
    rep("POU1F1 Tumoral cells", 6), #0 al 5
    "Macrophages", #6
    "POU1F1 Tumoral cells", #7
    "NK cells"
  ),
  # Muestra AC-002-2021
  # Cluster 7, 16 = Macrophages 
  # Cluster 11 = NK cells 
  # Cluster 8 = Endothelium  
  # Cluster 10 = Pericytes
  # Cluster 15, 4 = Progenitor cells 
  # Cluster 6, 9, 17, = POMC/GH tumor cells
  # Cluster 0, 1, 2, 3, 5, 12, 13, 14, = POMC/PRL tumor cells
  "AC_002" =  c(
    rep("POMC/PRL tumor cells", 4), #0 al 3
    "Progenitor cells", #4
    "POMC/PRL tumor cells", #5
    "POMC/GH tumor cells",  #6
    "Macrophages",          #7
    "Endothelium",          #8
    "POMC/GH tumor cells",  #9
    "Pericytes",            #10
    "NK cells ",            #11
    rep("POMC/PRL tumor cells", 3), #12 al 14
    "Progenitor cells",     #15
    "Macrophages",          #16 
    "POMC/GH tumor cells"   #17
  ),
  # Muestra CU-002-2021
  # Clusters 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 13, = POU1F1 Tumoral cells 
  # Cluster 11 = Pericytes
  # Cluster 14 = Endothelium  
  # Cluster 7 = Macrophages
  # Cluster 9 = NK cells
  "CU_002" = c(
    rep("POU1F1 Tumoral cells", 7),   #0 al 6
    "Macrophages",                    #7
    "POU1F1 Tumoral cells",           #8
    "NK cells",                       #9
    "POU1F1 Tumoral cells",           #10
    "Pericytes",                      #11
    rep("POU1F1 Tumoral cells", 2),   #12 y 13
    "Endothelium"                     #14
  ),
  # Muestra CP-001-2021
  # Clusters 0, 1, 2, 3, 4, 5, 7, 8, 10 = TBX19 Tumoral cells 
  # Cluster 6 = Macrophages 
  # Cluster 9 = NK cells
  # Cluster 13 = Erythroid
  # Cluster 12 = Pericytes
  # Cluster 11 = unidentified
  "CP_001" = c(
    rep("TBX19 Tumoral cells", 6), #0 al 5
    "Macrophages",                 #6
    rep("TBX19 Tumoral cells", 2), #7 y 8
    "NK cells",                    #9
    "TBX19 Tumoral cells",         #10
    "unidentified",                #11
    "Pericytes",                   #12
    "Erythroid"                    #13
  ),
  # Muestra SC-001-2021
  # Clusters 0, 1, 2, 3, 4, 5, 6, 8, 10, 11= NR5A1 Tumoral cells 
  # Cluster 7 = Macrophages 
  # Cluster 9 = NK cells
  "SC_001" = c(
    rep("NR5A1 Tumoral cells", 7),  #0 al 6
    "Macrophages",                  #7
    "NR5A1 Tumoral cells",          #8
    "NK cells",                     #9
    "NR5A1 Tumoral cells"           #10 no hay grupo 11
  ),
  # Muestra NF-002-2021
  # Cluster 9 = Macrophages 
  # Clusters 0, 1, 2, 3, 4, 5, 7, 8, 10 = NR5A1 Tumoral cells
  "NF_002" = c(
    rep("NR5A1 Tumoral cells", 9), #0 al 8 ojo faltaba el 6
    "Macrophages",                 #9
    "NR5A1 Tumoral cells"          #10
  ),
  # Muestra NF-003-2021
  # Cluster 4 = Stem Cells 
  # Clusters 0, 1, 2, 3, 5, = NR5A1 Tumoral cells 
  "NF_003" = c(
    rep("NR5A1 Tumoral cells", 4), #0 al 3
    "Stem Cells",                  #4
    "NR5A1 Tumoral cells"          #5
  ),
  # Muestra NF-004-2021
  # Clusters 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11= NR5A1 Tumoral cells 
  "NF_004" = c(
    rep("NR5A1 Tumoral cells", 13) #todos
  )
)

scRNAEtiquetados <- lapply(
  names(scRNA), 
  function(x){
    names(clusterid[[x]]) <- levels(scRNA[[x]])
    scRNA[[x]] <- RenameIdents(scRNA[[x]], clusterid[[x]])
    p <- DimPlot(scRNA[[x]], reduction = "umap", label = TRUE, pt.size = 0.5)
    ggsave(
      p,
      file = paste0("results/", "Etiqueta_UMAP_", x, ".pdf"),
      height = 6,
      width  = 7
    )
    return(scRNA[[x]])
  }
)
names(scRNAEtiquetados) <- names(scRNA)

#trajectorias modificar objetos
#modificar
subset()

sce <- slingshot(
  data = scRNA$AC_001$umap@cell.embeddings,
  clusterLabels = scRNA$AC_001@active.ident, 
  approx_points = 150,
)

plot(
  scRNA.integrated@reductions$umap@cell.embeddings, 
  col = scRNA.integrated@active.ident,
  pch=16, 
  asp = 1
)
lines(SlingshotDataSet(sce), lwd=2, col='black')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot(
  scRNA.integrated@reductions$umap@cell.embeddings, 
  col = gg_color_hue(16)[scRNA.integrated@active.ident], 
  pch = 16, 
  asp = 1
)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

####Uniendo por grupo tumoral
# 1).- NF-002-2021 + NF-003-2021 + NF-004-2021 + SC-001-2021
# 2.- AC-001-2021 + CU-002-2021 
# 3).- CP-001-2021 + AC-002-2021

#features for integration by Experimental group
features <- list()
features$NFs <- SelectIntegrationFeatures(
  object.list = scRNAEtiquetados[c("NF_002", "NF_003", "NF_004", "SC_001")]
)
features$ACs  <- SelectIntegrationFeatures(
  object.list = scRNAEtiquetados[c("AC_001", "CU_002")]
)
features$CPs  <- SelectIntegrationFeatures(
  object.list = scRNAEtiquetados[c("CP_001", "AC_002")]
)
mascara <- c(
  "NF_002" = "NFs", 
  "NF_003" = "NFs",
  "NF_004" = "NFs",
  "SC_001" = "NFs",
  "AC_001" = "ACs",
  "CU_002" = "ACs",
  "CP_001" = "CPs",
  "AC_002" = "CPs"
)
  

#Rerun the PCAs with the features
scRNAEtiquetados <- lapply(
  names(scRNAEtiquetados), 
  function(x) {
    RunPCA(
      scRNAEtiquetados[[x]], 
      features = features[[mascara[x]]], 
      verbose = FALSE
    )
  }
)
names(scRNAEtiquetados) <- names(scRNA)
#Do not run in parallel, disk dead for sure!!!! Sequential process required
# Go take a coffee
scRNAEtiquetadosGrupos <- list(
  NFs = scRNAEtiquetados[c("NF_002", "NF_003", "NF_004", "SC_001")],
  ACs = scRNAEtiquetados[c("AC_001", "CU_002")],
  CPs = scRNAEtiquetados[c("CP_001", "AC_002")]
)


#finding anchors
anchors <- lapply(
  names(features), 
  function(x){
    FindIntegrationAnchors(
      object.list = scRNAEtiquetadosGrupos[[x]], 
      anchor.features = features[[x]],
      scale = FALSE,
      reduction = "cca", 
      dims = 1:50
    )  
  }
)
#aca
names(anchors)  <-  names(features)

# Warning message:
#   In CheckDuplicateCellNames(object.list = object.list) :
#   Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.


#Create and integrate object
scRNAEtiquetados.integrated <- lapply(
  anchors, 
  function(x){
    out <- IntegrateData(anchorset = x, dims = 1:50)
    out <- SCTransform(
      out,
      method = "glmGamPoi", 
      vars.to.regress = "percent.mt", 
      verbose = FALSE
    )
    out <- RunPCA(out, verbose = FALSE)
    out <- RunUMAP(out, dims = 1:50)
    return(out)
  } 
)

#crear los clusters



# pasar a Sergio para los featureplots
# de los integrados por grupo experimental
# trajectoria

sce <- slingshot(
  data = scRNA.integrated@reductions$umap@cell.embeddings,
  clusterLabels = scRNA.integrated@active.ident, 
  approx_points = 150
)

plot(
  scRNA.integrated@reductions$umap@cell.embeddings, 
  col = scRNA.integrated@active.ident,
  pch=16, 
  asp = 1
)
lines(SlingshotDataSet(sce), lwd=2, col='black')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot(
  scRNA.integrated@reductions$umap@cell.embeddings, 
  col = gg_color_hue(16)[scRNA.integrated@active.ident], 
  pch = 16, 
  asp = 1
)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')




# Reintegro todo en 1
# pasar a sergio para featureplots
# trajectoria


#nuevos etiquetados
# expresion diferencial
# heatmaps
# proliferacion 
# SEAs + Venn



#etiquetar
#Revisar como si o no se pierden las etiquetas
scRNAEtiquetados.integrated


#Finally some visualization
DimPlot(
  scRNAEtiquetados.integrated$NFs, 
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.5
)

DimPlot(scRNAEtiquetados.integrated$NFs, group.by = "orig.ident")
DimPlot(scRNAEtiquetados.integrated$NFs, group.by = "orig.ident")





