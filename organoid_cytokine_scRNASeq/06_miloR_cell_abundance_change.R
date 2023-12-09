library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(miloR)
library(SingleCellExperiment)
library(patchwork)
library(dplyr)
library(statmod)
library(ggplot2)
library(RColorBrewer)

Convert("HVG_1000_latent_10_after_qc.h5ad", dest = "h5seurat", overwrite = TRUE)
sj <- LoadH5Seurat('HVG_1000_latent_10_after_qc.h5seurat')
sj$CQ <- paste0(sj$ID, '_', sj$timepoint, '_', sj$activation_status)
sj@reductions[['pca']] = sj@reductions$scVI

ml <- as.SingleCellExperiment(sj)
ml <- Milo(ml)
ml <- buildGraph(ml, k = 15, d = 10)
ml <- makeNhoods(ml, prop = 0.1, k = 15, d = 10, refined = TRUE)
ml <- countCells(ml, meta.data = as.data.frame(colData(ml)), sample="CQ")
ml_design <- data.frame(colData(ml))[, c("CQ", "activation_status")]
ml_design <- distinct(ml_design)
rownames(ml_design) <- ml_design[, 1]
ml_design[, 'activation_status'] <- factor(ml_design$activation_status, levels=c('no_cytokines', 'cytokines'))
ml <- calcNhoodDistance(ml, d=10)
ml_results <- testNhoods(ml, design = ~activation_status, design.df = ml_design)
ml <- buildNhoodGraph(ml)
ml_results <- annotateNhoods(ml, ml_results, coldata_col = "celltype")

pdf('miloR_res.pdf', height=8, width=8.5)
ctys <- c('TOM_VCT_proliferating', 'TOM_VCT', 'TOM_SCT', 'VCT', 'SCT', 'EVT_proliferating', 'EVT_early_1', 'EVT_early_2', 'EVT_early_3',
          'EVT_intermediate_1', 'EVT_intermediate_2', 'EVT_late_1', 'EVT_late_2', 'EVT_late_3')
ml_sres <- ml_results[ml_results$celltype_fraction > 0.7, ]
ml_sres$celltype <- factor(ml_sres$celltype, levels = rev(ctys))
print(plotDAbeeswarm(ml_sres, group.by = "celltype") + scale_colour_gradientn(colours = rev(brewer.pal(9, 'PuOr'))) + scale_x_discrete(position = "top"))
dev.off()
