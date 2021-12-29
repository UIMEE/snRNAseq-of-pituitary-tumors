library("ggplot2")
library("vroom")
library("readr")
library(Seurat)
library(dplyr)
library(patchwork)
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
################################################################################
################################################################################
#Etiquetando clusters por muestra
clusterid_sample <- list(
#  Muestra CU-002-2021
#  Clusters 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 13, = GH Tumoral cells 
#  Cluster 11 = Pericytes
#  Cluster 14 = Endothelium  
#  Cluster 7 = Macrophages
#  Cluster 9 grande= NK cells
#  Cluster 9 pequeño= B cells
"CU_002" = c(
    rep("GH Tumoral cells", 7),
    "Macrophages",
    "GH Tumoral cells",
    "NK cells & B cells",
    "GH Tumoral cells",
    "Pericytes",
    rep("GH Tumoral cells", 2),
    "Endothelium"
  ),
#  Muestra AC-001-2021
#  Clusters 0, 1, 2, 3, 4, 5, 7= GH Tumoral cells 
#  Cluster 6 = Macrophages 
#  Cluster 8 = NK cells
"AC_001" = c(
  rep("GH Tumoral cells", 6),
  "Macrophages",
  "GH Tumoral cells",
  "NK cells"
),
#  Muestra CP-001-2021
#  Clusters 0, 1, 2, 3, 4, 5, 7, 8, 10 = POMC Carcinoma cells 
#  Cluster 6 = Macrophages 
#  Cluster 9 = NK cells
#  Cluster 13 = Erythroid
#  Cluster 12 = Pericytes
#  Cluster 11 = unidentified 
"CP_001" = c(
  rep("POMC Carcinoma cells", 6),
  "Macrophages",
  "POMC Carcinoma cells",
  "POMC Carcinoma cells",
  "NK cells",
  "POMC Carcinoma cells",
  "unidentified",
  "Pericytes",
  "Erythroid"
),
#  Muestra AC-002-2021
#  Cluster 7, 16 = Macrophages 
#  Cluster 11 = NK cells 
#  Cluster 8 = Endothelium  
#  Cluster 10 = Pericytes
#  Cluster 15, 4 = Progenitor cells 
#  Cluster 6, 9, 17, = POMC/GH tumor cells
# Cluster 0, 1, 2, 3, 5, 12, 13, 14, = POMC/PRL tumor cells
"AC_002" = c(
  rep("POMC/PRL tumor cells", 4),
  "Progenitor cells",
  "POMC/PRL tumor cells",
  "POMC/GH tumor cells",
  "Macrophages",
  "Endothelium",
  "POMC/GH tumor cells",
  "Pericytes",
  "NK cells",
  rep("POMC/PRL tumor cells", 3),
  "Progenitor cells",
  "Macrophages",
  "POMC/GH tumor cells"
),
#  Muestra NF-002-2021
#  Cluster 9 = Macrophages 
#  Clusters 0, 1, 2, 3, 4, 5, 7, 8, 10 = NR5A1 Tumoral cells 
"NF_002" = c(
  rep("NR5A1 Tumoral cells", 9),
  "Macrophages",
  "NR5A1 Tumoral cells"
),
#  Muestra NF-003-2021
#  Cluster 4 = Stem Cells 
#  Clusters 0, 1, 2, 3, = NR5A1 Tumoral cells 
#  Cluster 5 = Macrophages 
"NF_003" = c(
  rep("NR5A1 Tumoral cells", 4),
  "Stem Cells",
  "Macrophages"
),
#  Muestra NF-004-2021
#  Clusters 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11= NR5A1 Tumoral cells 
"NF_004" = c(
  rep("NR5A1 Tumoral cells", 13)
), 
#  Muestra SC-001-2021
#  Clusters 0, 1, 2, 3, 4, 5, 6, 8, 10,= NR5A1 Tumoral cells 
#  Cluster 7 = Macrophages 
#  Cluster 9 = NK cells
"SC_001" = c(
  rep("NR5A1 Tumoral cells", 7),
  "Macrophages",
  "NR5A1 Tumoral cells",
  "NK cells",
  "NR5A1 Tumoral cells"
)  
)  

### Construyendo objeto con etiquetas para cada muestra y cada tipo celular
### Construyendo umaps
scRNAEtiquetados_by_sample <- lapply(
  names(scRNA), 
  function(x){
    names(clusterid_sample[[x]]) <- levels(scRNA[[x]])
    scRNA[[x]] <- RenameIdents(scRNA[[x]], clusterid_sample[[x]])
    dir.create(paste0(names(scRNA[x]), "_feature_plots"))
    pdf(file = paste0(names(scRNA[x]), "_feature_plots", "/", "umap.pdf"),
        paper = "letter")
    print(DimPlot(scRNA[[x]], reduction = "umap", label = TRUE))
    dev.off()
    return(scRNA[[x]])
  }
)

names(scRNAEtiquetados_by_sample) <- names(scRNA)

################################################################################
#Etiquetando clusters por grupo

clusterid_group <- list(  
#   UMAP_ACs
# Cluster 15 = B Cells
# Cluster 11 = Pericytes
# Cluster 13 = Endothelium  
# Cluster 7, 14 = Macrophages
# Cluster 10 = NK cells
# Clusters 0, 1, 3, 4, 5, 6, 8, 9, 12, = GH Tumoral cells 1
# Clusters 0, 1, 3, 4, 5, 6, 8, 9, 12, = GH Tumoral cells 2
ACs = c(
  rep("GH Tumoral cells", 7),
  "Macrophages",
  rep("GH Tumoral cells", 2),
  "NK cells",
  "Pericytes",
  "GH Tumoral cells",
  "Endothelium",
  "Macrophages",
  "B Cells"
),
  
# UMAP_CPs
# Cluster 21, 9 = Macrophages 1
# Cluster 17 = Macrophages 2
# Cluster 13 = NK cells
# Cluster 22 = Erythroid
# Cluster 14 = Pericytes
# Cluster 15 = Endothelium  
# Cluster 16 = POMC/GH tumor cells
# Cluster 19 = unidentified 
# Cluster 11 = Progenitor cells 
# Cluster 12, 4, 5, 20 = POMC/PRL tumor cells 
# Clusters 0, 1, 2, 3, 6, 7, 8, 10, 18 = POMC Carcinoma cells 
CPs = c(
  rep("POMC Carcinoma cells", 4),
  "POMC/PRL tumor cells",
  "POMC/PRL tumor cells",
  rep("POMC Carcinoma cells", 3),
  "Macrophages 1",
  "POMC Carcinoma cells",
  "Progenitor cells",
  "POMC/PRL tumor cells",
  "NK cells",
  "Pericytes",
  "Endothelium",
  "POMC/GH tumor cells",
  "Macrophages 2",
  "POMC Carcinoma cells",
  "unidentified",
  "POMC/PRL tumor cells",
  "Macrophages 1",
  "Erythroid"
),
# UMAP_NFs
# Cluster 11 = Macrophages 
# Clusters 2, 3, = NR5A1 Tumoral cells 1
# Clusters 13 = NR5A1 Tumoral cells 2
# Clusters 5, 6, 7, 8, 12 = NR5A1 Tumoral cells 3
# Clusters 0, 1, 4, 9, 10, 14, 15 = NR5A1 Tumoral cells 4
# Cluster pequeño sin numero en la parte superior = Pericytes
# Cluster 17 = NK cells
# Cluster 16 = Stem cells 
NFs = c(
  "NR5A1 Tumoral cells 4",
  "NR5A1 Tumoral cells 4",
  "NR5A1 Tumoral cells 1",
  "NR5A1 Tumoral cells 1",
  "NR5A1 Tumoral cells 4",
  rep("NR5A1 Tumoral cells 3", 4),
  "NR5A1 Tumoral cells 4",
  "NR5A1 Tumoral cells 4",
  "Macrophages",
  "NR5A1 Tumoral cells 3",
  "NR5A1 Tumoral cells 2",
  "NR5A1 Tumoral cells 4",
  "NR5A1 Tumoral cells 4",
  "Stem cells",
  "NK cells"
  )
)

scRNAEtiquetados_by_group <- lapply(
  names(scRNA), # hay que cambiar por el objeto que contenga los seurat por grupo
  function(x){
    names(clusterid_group[[x]]) <- levels(scRNA[[x]])
    scRNA[[x]] <- RenameIdents(scRNA[[x]], clusterid_group[[x]])
    dir.create(paste0(names(scRNA[x]), "_feature_plots"))
    pdf(file = paste0(names(scRNA[x]), "_feature_plots", "/", "umap.pdf"),
        paper = "letter")
    print(DimPlot(scRNA[[x]], reduction = "umap", label = TRUE))
    dev.off()
    return(scRNA[[x]])
  }
)

names(scRNAEtiquetados_by_group) <- names(scRNA)

################################################################################ 
# UMAP_Final
# Cluster 14, 17 = Macrophages 1
# Cluster 23 = Macrophages 1
# Cluster 18 = NK cells
# Cluster E = Erythroid
# Cluster 22 = Endothelium  
# Cluster 20 = Pericytes
# Cluster 24 = Stem cells 
# Clusters 5, 7, 8, 9, 10 = POMC Carcinoma cells 
# Clusters 2, 13 = GH Tumoral cells 1
# Clusters 15 = GH Tumoral cells 2
# Cluster 19 = Progenitor cells 
# Cluster 4 = POMC/PRL tumor cells 
# Cluster dentro del cluster 4 hay un circulo que delimita = POMC/PRL tumor cells 
# Clusters 12, 3 = NR5A1 Tumoral cells 1
# Cluster 21 = NR5A1 Tumoral cells 2
# Clusters 11, 1 = NR5A1 Tumoral cells 3
# Clusters 16, 6, 0 = NR5A1 Tumoral cells 4
# Cluster U = Unidentified

clusterid_final <- c(
  "NR5A1 Tumoral cells 4",
  "NR5A1 Tumoral cells 3",
  "GH Tumoral cells 1",
  "NR5A1 Tumoral cells 1",
  "POMC/PRL tumor cells",
  "POMC Carcinoma cells",
  "NR5A1 Tumoral cells 4",
  rep("POMC Carcinoma cells", 4),
  "NR5A1 Tumoral cells 3",
  "NR5A1 Tumoral cells 1",
  "GH Tumoral cells 1",
  "Macrophages 1",
  "GH Tumoral cells 2",
  "NR5A1 Tumoral cells 4",
  "Macrophages 1",
  "NK cells",
  "Progenitor cells",
  "Pericytes",
  "NR5A1 Tumoral cells 2",
  "Endothelium",
  "Macrophages 1",
  "Stem cells"
  )

#Revisar estos comando para cambiar el nombre de cluster y construcción del umap
#cambiar scRNA por el objeto seurat con el objeto scRNA_final (todas las muestras)
scRNA <- RenameIdents(scRNA, clusterid_final)
dir.create("final_feature_plots")
pdf(file = paste0("final_feature_plots", "/", "umap.pdf"),
    paper = "letter")
print(DimPlot(scRNA, reduction = "umap", label = TRUE))
dev.off()

################################################################################
################################################################################

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


#################


