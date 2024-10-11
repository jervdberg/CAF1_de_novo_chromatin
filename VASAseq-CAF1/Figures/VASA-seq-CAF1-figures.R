library(Seurat)
#library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(Matrix)
library(rlang)
library(leiden)
library(reticulate)
library(tidyverse)
library(Rfast)
library(rstatix)

###########rep1$###########

use_python("/Users/jeroen/miniconda3/bin/python", required=T)

JD <- readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/unfiltcounts.rds')
#JD <- JD$unspliced
ENSG_to_name <- readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/ENSGtoName.rds')
rownames(JD) <- ENSG_to_name[rownames(JD), "gene_out"]

JD_VASA <- CreateSeuratObject(counts = JD, project = "JD")

meta <-  readRDS("~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/unfiltmeta.rds")


JD_VASA <-AddMetaData(JD_VASA, metadata = meta)
JD_VASA@meta.data$line <- "RPE1 CAF1-dTAG"
JD_VASA@meta.data$rep <- "rep1"
JD_VASA@meta.data$cond <- JD_VASA@meta.data$plate

print("Before filtering")
JD_VASA
JD_VASA = subset(JD_VASA, subset = nFeature_RNA >6500)
print("After filtering")
JD_VASA
JD_VASA = subset(JD_VASA, subset = nCount_RNA >800)


#########rep2############


JD2 <- readRDS('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/unfiltcounts.rds')
#JD <- JD$unspliced
ENSG_to_name <- readRDS('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/ENSGtoName.rds')
#rownames(JD2) <- ENSG_to_name[rownames(JD2), "gene_out"]

JD2_VASA <- CreateSeuratObject(counts = JD2, project = "JD2")

meta <-  readRDS("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/unfiltmeta.rds")


JD2_VASA <-AddMetaData(JD2_VASA, metadata = meta)
JD2_VASA@meta.data$line <- "RPE1 CAF1-dTAG"
JD2_VASA@meta.data$rep <- "rep2"
JD2_VASA@meta.data$cond <- JD2_VASA@meta.data$plate

print("Before filtering")
JD2_VASA
JD2_VASA = subset(JD2_VASA, subset = nFeature_RNA >6500)
JD2_VASA = subset(JD2_VASA, subset = nCount_RNA >800)
print("After filtering")
JD2_VASA



DMSO_VASA@meta.data


DMSO <- readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/silv_dmso/counts_dmso_symbols.rds')
#ENSG_to_name1 <- readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/silv_dmso/JD-VASA/compiled/ENSGtoName.rds')
#rownames(DMSO)  <-ENSG_to_name1[rownames(DMSO), "gene_out"]
DMSO_VASA <- CreateSeuratObject(counts = DMSO, project = "DMSO" )
meta_DMSO <-  readRDS("~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/silv_dmso/meta_dmso.rds")
DMSO_VASA <-AddMetaData(DMSO_VASA, metadata = meta_DMSO)
DMSO_VASA@meta.data$line <- "RPE1 WT"
DMSO_VASA@meta.data$cond <- "WT DMSO"
DMSO_VASA@meta.data$rep <- "rep0"


DMSO_VASA <- SplitObject(DMSO_VASA, split.by = "sort_population")
DMSO_VASA <- DMSO_VASA$Interphase

JD_VASA <- merge(JD_VASA, JD2_VASA)
JD_VASA <- merge(JD_VASA, DMSO_VASA)
JD_VASA <- subset(JD_VASA, subset = nFeature_RNA >4000)
JD_VASA <- subset(JD_VASA, subset = nCount_RNA >1000)


unique(JD_VASA@meta.data$sort_population)





JD.list <- SplitObject(JD_VASA, split.by = "rep")

JD.list <- lapply(X = JD.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = JD.list)

JD.anchors <- FindIntegrationAnchors(object.list = JD.list, anchor.features = features)

JD.combined <- IntegrateData(anchorset = JD.anchors)

DefaultAssay(JD.combined) <- "integrated"




p1 = VlnPlot(JD.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 1, pt.size = 0.2) 


p1.2 <- FeatureScatter(JD.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1.2

JD.combined@meta.data %>% group_by(cond)%>% summarize(n=n())

all.genes = rownames(JD.combined)
JD.combined = ScaleData(JD.combined, features = all.genes)

JD.combined = RunPCA(JD.combined, features = VariableFeatures(JD.combined))

p1.4 = VizDimLoadings(JD.combined, dims = 1:2, reduction = "pca")

p1.4

p1.5 = DimPlot(JD.combined, reduction = "pca") + ggtitle("JD.combined")

p1.5


p1.6 = ElbowPlot(JD.combined) + ggtitle("JD.combined")
p1.6


JD.combined = FindNeighbors(JD.combined, dims = 1:15, annoy.metric = 'manhattan')
JD.combined = FindClusters(JD.combined, resolution = 1, algorithm = 4, method = "igraph")

JD.combined <- RunUMAP(JD.combined, dims = 1:15, n.neighbors = 10, min.dist = 0.3)
JD.combined <- RunTSNE(JD.combined, dims = 1:15)
JD.combined@meta.data$plate

p1.7 = DimPlot(JD.combined, reduction = "umap", pt.size = 1,  group.by = 'line', shuffle=T) + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'right')+coord_equal()

p1.7a = DimPlot(JD.combined, reduction = "umap", pt.size = 0.8, group.by = 'cond', seed = 42,shuffle = T) + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'none')+coord_equal()+
  scale_color_manual(values = c("#1F97FF",
                               "#1F97FF",
                               "#1F97FF",
                               "#1F97FF",
                               "#FF7B15", 
                               "#FF7B15" ,
                               "#FF7B15" ,
                               "#FF7B15",
                               "grey80"))

p1.7a[[1]]$layers[[1]]$aes_params$alpha =  .9
p1.7a
p1.8 = DimPlot(JD.combined, reduction = "umap", pt.size = 1) + ggtitle("")+theme(legend.position = 'right')+coord_equal()

DimPlot(JD.combined, reduction = "umap", pt.size = 1,  split.by = 'rep') + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'bottom')+coord_equal()

DimPlot(JD.combined, reduction = "umap", pt.size = 1,split.by = 'cond', ncol = 4) + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'bottom')+coord_equal()

DimPlot(JD.combined, reduction = "umap", pt.size = 1,split.by = 'rep', ncol = 4) + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'bottom')+coord_equal()

DimPlot(JD.combined, reduction = "umap", pt.size = 1,split.by = 'plate', ncol = 4) + 
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'bottom')+coord_equal()

(p1.7/p1.7a/p1.8)+plot_layout(heights = c(4,4,4))

JD_2h = JD.list <- SplitObject(JD.combined, split.by = "cond")
JD_2h1 <- merge(JD_2h$`DMSO-2h`, JD_2h$`V1-2h`, merge.data = TRUE, merge.dr = 'umap')
JD_2h2 <- merge(JD_2h1, JD_2h$`DMSO-2h-pl2`, merge.data = TRUE, merge.dr = 'umap')
JD_2h3 <- merge(JD_2h2, JD_2h$`V1-2h-pl2`, merge.data = TRUE, merge.dr = 'umap')  
  
JD_30h = JD.list <- SplitObject(JD.combined, split.by = "cond")
JD_30h1 <- merge(JD_30h$`DMSO-30h`, JD_30h$`V1-30h`, merge.data = TRUE, merge.dr = 'umap')
JD_30h2 <- merge(JD_30h1, JD_30h$`DMSO-30h-pl2`, merge.data = TRUE, merge.dr = 'umap')
JD_30h3 <- merge(JD_30h2, JD_30h$`V1-30h-pl2`, merge.data = TRUE, merge.dr = 'umap')  

p1.7b <- DimPlot(JD_2h3, reduction = "umap", pt.size = 0.8, group.by = 'cond', seed = 42,shuffle = T) + 
  ggtitle("2h")+theme(legend.position = 'none')+coord_equal()+
  scale_color_manual(values = c("#1F97FF",
                                "#1F97FF",
                                "#FF7B15", 
                                "#FF7B15" 
                             ))

p1.7c <- DimPlot(JD_30h3, reduction = "umap", pt.size = 0.8, group.by = 'cond', seed = 42,shuffle = T) + 
  ggtitle("30h")+theme(legend.position = 'none')+coord_equal()+
  scale_color_manual(values = c("#1F97FF",
                                "#1F97FF",
                                "#FF7B15", 
                                "#FF7B15" 
  ))

p1.7b/p1.7c

p1.7d <- DimPlot(JD_CAF1, reduction = "umap", pt.size = 1,split.by = 'rep', ncol = 1) + 
  ggtitle("")+theme(legend.position = 'none')+coord_equal()



FeaturePlot(JD.combined, pt.size = 1, features = c("nCount_RNA"), order = T)+coord_equal()
FeaturePlot(JD.combined, pt.size = 1, features = c("nFeature_RNA"), order = T)+coord_equal()

p1.8



JD.combined@meta.data$green

DefaultAssay(JD.combined) <- "RNA"



JD.combined_markers = FindAllMarkers(JD.combined, min.pct = 0.9, logfc.threshold = 0.25, only.pos = T)
JD.combined_top10_DEG = JD.combined_markers %>%  dplyr::filter(cluster ==2) %>% 
  dplyr::filter(p_val_adj < 1e-3)%>%
    arrange(avg_log2FC) 

overlap6_9<- JD.combined_markers %>% dplyr::filter(cluster == 6 | cluster == 9)%>% group_by(cluster)%>% nest()%>%
  ungroup() %>% summarise(intersect_feat = map(1:n_distinct(.$cluster), function(x).$data[[x]]$gene) %>% reduce(intersect))%>%
  mutate(new_gene = str_remove(intersect_feat, '.[0-9]{1,2}-.*'))%>% rownames_to_column()%>%dplyr::select(new_gene) %>% as.data.table()



overlap6_9%>% write_tsv('overlap_GO1.txt')


JD.combined_top10_DEG = JD.combined_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 


JD.combined_markers %>% dplyr::filter(cluster == 6 & avg_log2FC >0.25)%>% arrange(desc(avg_log2FC), p_val_adj)
  mutate(new_gene = str_remove(gene, '.[0-9]{1,2}-.*'))%>% rownames_to_column()%>%select(new_gene)
  

JD.combined@meta.data$green

p1.11<- FeaturePlot(JD.combined, pt.size = 1, features = c("red"), order = T)+coord_equal()
p1.11<-p1.11+scale_colour_gradientn(colours = brewer.pal(3,"Reds"), na.value = "grey92", trans = "log10")+ggtitle('hCdt1-mKO')


p1.10<- FeaturePlot(JD.combined, pt.size = 1,  features = c("green"), order = T)+coord_equal()
p1.10 <- p1.10+scale_colour_gradientn(colours = brewer.pal(3,"Greens"), na.value = "grey92", trans = "log10")+ggtitle('hGem-mAG')

(p1.7|p1.8)/(p1.10/p1.11)



library(data.table)
rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., 'MCM2'))

FeaturePlot(JD.combined, features = c("ENSG00000094804.12-CDC6"),pt.size=1)

FeaturePlot(JD.combined, features = c("ENSG00000163950.13-SLBP"),pt.size=1)

FeaturePlot(JD.combined, features = c("ENSG00000148773.14-MKI67","ENSG00000094804.12-CDC6",   "ENSG00000287080.1-H3C3","ENSG00000134057.15-CCNB1"),pt.size=1, ncol=2,
            cols = rev(brewer.pal(n = 10, name = "Spectral")))


DefaultAssay(JD.combined) <- "RNA"




FeaturePlot(JD.combined, pt.size = 1.5,  features = c("ENSG00000158406.4-H4C8", "ENSG00000276368.1-H2AC14",
                                  "ENSG00000184357.4-H1-5","ENSG00000203814.6-H2BC18"),order = T, split.by = 'rep', ncol = 2,
            cols = rev(brewer.pal(n = 8, name = "Spectral"  )))


FeaturePlot(JD.combined, pt.size = 1.5, features = c("ENSG00000158373.8-H2BC5",
                                  "ENSG00000168298.6-H1-4", "ENSG00000286522.1-H3C2", "ENSG00000287080.1-H3C3"),order = T,  split.by = 'cond', ncol = 2,
            cols = rev(brewer.pal(n = 5, name = "Spectral")))


VlnPlot(JD.combined, pt.size = 1.5, features = c("ENSG00000158373.8-H2BC5",
                                                     "ENSG00000168298.6-H1-4", "ENSG00000286522.1-H3C2", "ENSG00000287080.1-H3C3"),  split.by = 'cond', ncol = 2,
            cols = rev(brewer.pal(n = 5, name = "Spectral")))

VlnPlot(JD.combined, features = c( "ENSG00000067369.13-TP53BP1",
                                       "ENSG00000078804.13-TP53INP2",
                                       "ENSG00000115129.14-TP53I3",
                                   "ENSG00000105327.17-BBC3",
                                       "ENSG00000141510.17-TP53"
                                
                                  ),pt.size=1, ncol = 2, group.by='seurat_clusters',
            cols = rev(brewer.pal(n = 10, name = "Spectral")))

DefaultAssay(JD.combined) <- "RNA"
FeaturePlot(JD.combined, features =c(
  "ENSG00000124762.13-CDKN1A",
  "ENSG00000001617.12-SEMA3F",
  "ENSG00000108691.9-CCL2",
  "ENSG00000163453.11-IGFBP7",
  "ENSG00000113328.19-CCNG1",
  "ENSG00000135679.25-MDM2",
  "ENSG00000198625.13-MDM4"
),pt.size=0.75, ncol = 2, order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))

FeaturePlot(JD.combined, features =c(
  "ENSG00000124762.13-CDKN1A",
  "ENSG00000163453.11-IGFBP7",
  "ENSG00000113328.19-CCNG1",
  "ENSG00000135679.25-MDM2",
  "ENSG00000141510.17-TP53"
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))




FeaturePlot(JD.combined, features =c(
  "ENSG00000100867.14-DHRS2" 
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))



FeaturePlot(JD.combined, features =c(
  "ENSG00000163453.11-IGFBP7"
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))



FeaturePlot(JD.combined, features =c(
  "ENSG00000124762.13-CDKN1A"
),pt.size=0.75, order=T, split.by = 'cond',
cols = rev(brewer.pal(n = 10, name = "Spectral")))

VlnPlot(JD.combined, features =c(
  "ENSG00000124762.13-CDKN1A"
),group.by = 'plates')

JD.combined@meta.data$seurat_clusters
RidgePlot(JD.combined, features =c(
  "ENSG00000148773.14-MKI67"
),group.by = 'cond')

"ENSG00000148773.14-MKI67"

ENSG00000167670.16-CHAF1A

VlnPlot(JD.combined, features =c(
  "ENSG00000124762.13-CDKN1A"
),group.by = 'cond', pt.size=0.75)

VlnPlot(JD.combined, features =c(
  "ENSG00000167670.16-CHAF1A"
),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))


VlnPlot(JD.combined, features =c(
  "ENSG00000148773.14-MKI67"
  ),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))

VlnPlot(JD.combined, features =c(
  "ENSG00000132475.10-H3-3B"
),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))


VlnPlot(JD.combined, features =c(
  "ENSG00000132475.10-H3-3B"
),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))



VlnPlot(JD.combined, features =c(
  "ENSG00000275714.1-H3C1"
),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))



VlnPlot(JD.combined, features =c(
  "ENSG00000231617.8-DAXX"
),group.by = 'cond', pt.size=0.75, cols = c("#FFDB6D",  "#D16103", "#52854C", "#4E84C4", "lightgrey"))



FeaturePlot(JD.combined, features =c(
    "ENSG00000124762.13-CDKN1A",
  "ENSG00000001617.12-SEMA3F",
  "ENSG00000108691.9-CCL2",
  "ENSG00000163453.11-IGFBP7",
  "ENSG00000113328.19-CCNG1"
),pt.size=0.75, ncol = 2, order=T, split.by = 'cond',
cols = rev(brewer.pal(n = 10, name = "Spectral")))

DimPlot(JD.combined,pt.size=1.5, split.by = 'cond')

ENSG00000141458.13-NPC1

FeaturePlot(JD.combined, features =c(
  "ENSG00000100867.14-DHRS2" 
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))

FeaturePlot(JD.combined, features =c(
  "ENSG00000187764.11-SEMA4D",
 "ENSG00000197632.9-SERPINB2",
   "ENSG00000103740.10-ACSBG1",
       "ENSG00000134363.12-FST"
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 5, name = "Spectral")))



DHRS2a <- FeaturePlot(JD_CAF1, features =c(
  "ENSG00000100867.14-DHRS2" 
),pt.size=1,order=T,
cols = rev(brewer.pal(n = 10, name = "Spectral")))

DHRS2aa <- VlnPlot(JD_CAF1, features =c(
  "ENSG00000100867.14-DHRS2" 
),pt.size=1, group.by = 'cond2',
cols = rev(brewer.pal(n = 10, name = "Spectral")))



 DHRS2a/DHRS2aa+plot_layout(heights = c(3,3))
 
 
 
 DNAH3a <- FeaturePlot(JD.combined, features =c(
   "ENSG00000158486.13-DNAH3"
 ),pt.size=1,order=T,
 cols = rev(brewer.pal(n = 10, name = "Spectral")))
 
 DNAH3b <- VlnPlot(JD.combined, features =c(
   "ENSG00000158486.13-DNAH3"
 ),pt.size=1, group.by = 'cond',
 cols = rev(brewer.pal(n = 10, name = "Spectral")))


 


 DNAH3a/DNAH3b

 
 
 KRT18a <- FeaturePlot(JD.combined, features =c(
   "ENSG00000205426.10-KRT81"
 ),pt.size=1,order=T,
 cols = rev(brewer.pal(n = 10, name = "Spectral")))
 
KRT81b <- VlnPlot(JD_CAF1, features =c(
  "ENSG00000205426.10-KRT81" 
 ),pt.size=1, group.by = 'cond2',
 cols = rev(brewer.pal(n = 10, name = "Spectral")))

 
 
 saveRDS(JD.combined, "~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/RPE1-dTAG-CAF1-2-30hr.rds")
 
 
 JD_sphase = subset(x=JD.combined, seurat_clusters == c("3", "7"))

 JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")

 table(JD_CAF1@meta.data$cond,JD_CAF1@meta.data$seurat_clusters )%>% as.data.frame()%>%dplyr::filter(Freq>0)%>%
   setNames(c("cond", "cluster", "n"))%>% arrange(n)%>%
   ggplot()+geom_col(aes(x=cond, y=n, fill=cluster), position = 'fill')+

  theme_bw()
 
 JD_sphase@assays$RNA %>% 
   as.data.frame()%>%
   tibble::rownames_to_column(var = 'gene')
 
 
 JD_CAF1@meta.data <- JD_CAF1@meta.data %>% mutate(cond2 = gsub("-pl2", "", cond))
 JD_2h_DMSO <-subset(x=JD_CAF1, cond2 == c("DMSO-2h"))
 JD_2h_dTAG <- subset(x=JD_CAF1, cond2 == c("V1-2h"))
 JD_30h_DMSO <-subset(x=JD_CAF1, cond2 == c("DMSO-30h"))
 JD_30h_dTAG <- subset(x=JD_CAF1, cond2 == c("V1-30h"))

 
 
 volcano2 <-  bind_rows(
     (JD_2h_DMSO[["RNA"]]@counts %>%   as.data.frame()%>%
        tibble::rownames_to_column(var = 'gene')%>%
        pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "DMSO-2h")),
     
     (JD_2h_dTAG[["RNA"]]@counts %>%   as.data.frame()%>%
        tibble::rownames_to_column(var = 'gene')%>%
        pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "V1-2h")))%>% as.data.table()


 
 volcano30 <-  bind_rows(
   (JD_30h_DMSO[["RNA"]]@counts %>%   as.data.frame()%>%
      tibble::rownames_to_column(var = 'gene')%>%
      pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "DMSO-30h")),
   
   (JD_30h_dTAG[["RNA"]]@counts %>%   as.data.frame()%>%
      tibble::rownames_to_column(var = 'gene')%>%
      pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "V1-30h")))%>% as.data.table()
 
 
 cells2 <- volcano2[, .(count_tot = sum(count)), .(cell, cluster)
 ][, clusterN := uniqueN(cell), .(cluster)
 ]
 
 cells30 <- volcano30[, .(count_tot = sum(count)), .(cell, cluster)
 ][, clusterN := uniqueN(cell), .(cluster)
 ]
 
 
 
 volcano2_p <- volcano2[count > 0
 ][, count_norm := log1p(count / cells2[, setNames(count_tot, cell)][cell] * median(cells2$count_tot))
 ][order(gene, cluster), .(count_norm = list(count_norm)), 
   .(gene, cluster, clusterN = as.factor(cells2[,setNames(clusterN, cell)][cell]))
   ][, data.table(map2_dfr(count_norm, as.numeric(as.character(clusterN)), function(x, cN){
     test <- nafill(x[1:cN], fill = 0)  
     rest <- nafill(c(0, unlist(count_norm[clusterN != cN], F, F))[2:(sum(as.numeric(levels(clusterN))) - cN + 1)], fill = 0)
     data.table(pvalue = t.test(test, rest)$p.value,
                logfold = log2(mean(test) / mean(rest)))
   }), cluster = cluster), .(gene)
   ][, pvalue_cor := p.adjust(pvalue, "BH"), .(cluster)]
 
 volcano30_p <- volcano30[count > 0
 ][, count_norm := log1p(count / cells30[, setNames(count_tot, cell)][cell] * median(cells30$count_tot))
 ][order(gene, cluster), .(count_norm = list(count_norm)), 
   .(gene, cluster, clusterN = as.factor(cells30[,setNames(clusterN, cell)][cell]))
 ][, data.table(map2_dfr(count_norm, as.numeric(as.character(clusterN)), function(x, cN){
   test <- nafill(x[1:cN], fill = 0)  
   rest <- nafill(c(0, unlist(count_norm[clusterN != cN], F, F))[2:(sum(as.numeric(levels(clusterN))) - cN + 1)], fill = 0)
   data.table(pvalue = t.test(test, rest)$p.value,
              logfold = log2(mean(test) / mean(rest)))
 }), cluster = cluster), .(gene)
 ][, pvalue_cor := p.adjust(pvalue, "BH"), .(cluster)]
 
 volcano2_p %>% dplyr::filter(str_detect(cluster, "DMSO"))%>%
  select(gene, logfold, pvalue_cor) %>% mutate(comparison = "2hr log2(DMSO/dTAG)")%>%
   write_tsv(.,"~/Desktop/post-doc/Collaberations/Francesca/2hr_LFC_adjP_BH.tsv")
 
 volcano30_p %>% dplyr::filter(str_detect(cluster, "DMSO")) %>%
   select(gene, logfold, pvalue_cor) %>% mutate(comparison = "30hr log2(DMSO/dTAG)")%>% 
   write_tsv(.,"~/Desktop/post-doc/Collaberations/Francesca/30hr_LFC_adjP_BH.tsv")
 
 
 hits <- volcano2_p %>%  
   dplyr::filter(pvalue_cor<0.01 & logfold > 0.25 | pvalue_cor<0.01 & logfold < -0.25)%>%
   dplyr::filter(str_detect(cluster, "DMSO"))%>%
   mutate(Transcript = if_else( str_detect(gene, "H4C|H2AC|H2BC|H3C"), true = "Histone transcripts", false = "Other Transcripts"))

   
   
   
 ns <- volcano2_p%>% 
   dplyr::filter(str_detect(cluster, "DMSO"))
 
   ggplot()+
   geom_point(data = ns, aes(y=-log10(pvalue_cor), x=logfold), col= "grey80", alpha=0.5, size= 0.25)+
   geom_point(data = hits, aes(y=-log10(pvalue_cor), x=logfold, color=Transcript), size= 0.75)+
   theme_bw()+
   geom_vline(xintercept = c(-0.25, 0.25), linetype =2)+
  xlab('Log2 Fold Change (2h DMSO/dTAG)')+ #xlim(-2.5,2.5)+
   ggtitle('CHAF1A-dTAG 2h dTAG')+
   ylab('-log10(adjusted p-value')+coord_fixed(ratio = 0.4)
 
 
volcano30_p[logfold > 0.5, prank := frank(pvalue, ties.method = "first"), .(cluster)]

volcano30_p %>% arrange(prank)

library(Rfast)


volcano[gene %in% volcano_p[(prank < 20)]$gene & count > 0, !"cluster"
 ][, .SD[cells[, .(cell, cluster)], on = .(cell)], .(gene)
 ][is.na(count), count := 0
 ][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
 ][, countz := (count_norm - mean(count_norm)) / sd(count_norm), .(gene)
 ][, cellorder := as.numeric(as.factor(cell)), .(cluster)
 ][, gene := factor(gene,
                    levels = dcast(.SD[,.(gene, cell, countz)], 
                                   gene ~ cell, value.var = "countz", fill = 0) %>%
                      {as.dist(as.matrix(data.frame(Dist((.[,!"gene"]), method = "manhattan"), 
                                                    row.names = .$gene)))} %>%
                      hclust("ward.D") %>% 
                      {.$labels[.$order]})
 ][abs(countz) > 1.5, countz := sign(countz) * 1.5
   #[leidenn == 3
 ][] %>% ggplot()+
   geom_raster(aes(x = cellorder, y = gene, fill = countz)) +
   facet_grid(cols = vars(cluster),  scales = "free", space = "free") +
   scale_fill_gradient2(low = "magenta", mid = "black", high = "green") +
   theme(panel.grid = element_blank(), 
         axis.text.x = element_blank(), 
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(), text = element_text(size = 9)) +
   coord_cartesian(expand = F)+ggtitle('S-phase cells')
 
 
Up_dTAG <-volcano[gene %in% volcano_p[(prank < 10)]$gene & count > 0, !"cluster"
][, .SD[cells[, .(cell, cluster)], on = .(cell)], .(gene)
][is.na(count), count := 0
][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
][, countz := (count_norm - mean(count_norm)) / sd(count_norm), .(gene)
][, cellorder := as.numeric(as.factor(cell)), .(cluster)
][, gene := factor(gene,
                   levels = dcast(.SD[,.(gene, cell, countz)], 
                                  gene ~ cell, value.var = "countz", fill = 0) %>%
                     {as.dist(as.matrix(data.frame(Dist((.[,!"gene"]), method = "manhattan"), 
                                                   row.names = .$gene)))} %>%
                     hclust("ward.D") %>% 
                     {.$labels[.$order]})
][abs(countz) > 1.5, countz := sign(countz) * 1.5
  #[leidenn == 3
][]%>% dplyr::filter(cluster=="DMSO") %>% group_by(gene)%>% summarize(countz = mean(countz))%>% dplyr::filter(countz<0)#%>%
  mutate(new_gene = str_remove(gene, '.[0-9]{1,2}-.*'))%>% rownames_to_column()%>%select(new_gene)
 
  
##############cluster 4vs 9 #########
  
  
  
  JD_eG1phase_DMSO <-subset(x=JD.combined, seurat_clusters == c("5"))
  
  JD_eG1phase_dTAG <- subset(x=JD.combined, seurat_clusters == c( "6"))
  
  JD_eG1phase_DMSO <- subset(x= JD_eG1phase_DMSO, line == "RPE1 CAF1-dTAG")
  JD_eG1phase_dTAG <- subset(x= JD_eG1phase_dTAG, line == "RPE1 CAF1-dTAG")


  cells <- volcano_eG1[, .(count_tot = sum(count)), .(cell, cluster)
  ][, clusterN := uniqueN(cell), .(cluster)
  ]
  
  volcano_p_eG1 <- volcano_eG1[count > 0
  ][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
  ][order(gene, cluster), .(count_norm = list(count_norm)), 
    .(gene, cluster, clusterN = as.factor(cells[,setNames(clusterN, cell)][cell]))
  ][, data.table(map2_dfr(count_norm, as.numeric(as.character(clusterN)), function(x, cN){
    test <- nafill(x[1:cN], fill = 0)  
    rest <- nafill(c(0, unlist(count_norm[clusterN != cN], F, F))[2:(sum(as.numeric(levels(clusterN))) - cN + 1)], fill = 0)
    data.table(pvalue = t.test(test, rest)$p.value,
               logfold = log2(mean(test) / mean(rest)))
  }), cluster = cluster), .(gene)
  ][, pvalue_cor := p.adjust(pvalue, "bonferroni"), .(cluster)][]
  
  
  volcano_p_eG1 %>%  mutate(sig = pvalue_cor<0.05,  updown = logfold >0.5 | logfold < -0.5)%>% dplyr::filter(cluster == 'DMSO')%>%
    ggplot+
    geom_point(aes(y=-log10(pvalue_cor), x=logfold, fill= sig & updown), size= 1.5, alpha=0.8,shape =21, stroke=0.15)+
    theme_bw()+
    scale_fill_manual(values = c('grey',  'blue2', 'yellow2' ))+xlab('log2(DMSO/dTAG)')+ggtitle('30hrs')+
    ylab('Bonferoni corrected T-test\n-log10(p-value)')
  
  
  volcano_p_eG1[logfold > 0.5, prank := frank(pvalue, ties.method = "first"), .(cluster)]
  
  library(Rfast)
  
  
  volcano_eG1[gene %in% volcano_p_eG1[(prank < 50)]$gene & count > 0, !"cluster"
  ][, .SD[cells[, .(cell, cluster)], on = .(cell)], .(gene)
  ][is.na(count), count := 0
  ][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
  ][, countz := (count_norm - mean(count_norm)) / sd(count_norm), .(gene)
  ][, cellorder := as.numeric(as.factor(cell)), .(cluster)
  ][, gene := factor(gene,
                     levels = dcast(.SD[,.(gene, cell, countz)], 
                                    gene ~ cell, value.var = "countz", fill = 0) %>%
                       {as.dist(as.matrix(data.frame(Dist((.[,!"gene"]), method = "manhattan"), 
                                                     row.names = .$gene)))} %>%
                       hclust("ward.D") %>% 
                       {.$labels[.$order]})
  ][abs(countz) > 1.5, countz := sign(countz) * 1.5
    #[leidenn == 3
  ] %>% ggplot()+
    geom_raster(aes(x = cellorder, y = gene, fill = countz)) +
    facet_grid(cols = vars(cluster),  scales = "free", space = "free") +
    scale_fill_gradient2(low = "magenta", mid = "black", high = "green") +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), text = element_text(size = 9)) +
    coord_cartesian(expand = F)+ggtitle('30hrs')
  
  library(MASS)
  library(tidyverse)


  table(JD_CAF1@meta.data$cond,JD_CAF1@meta.data$seurat_clusters) %>% as.data.frame()%>%
    setNames(c("cond", "cluster", "n"))%>% arrange(n)%>% mutate(cond2 = gsub('-pl2', '', cond))%>%
    ggplot()+geom_col(aes(x=cond, y=n, fill=cluster))+
    theme_bw()+facet_grid(cols=vars)
  
  tab <- JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>%
    mutate(cond = gsub('-pl2', '', cond))%>% group_by(rep, cond, seurat_clusters)%>% summarize(n=n())
  
  
  summary(tab)
  
  tab%>% ggplot()+geom_col(aes(x=cond, y=n, fill=seurat_clusters))+
       theme_bw()+facet_grid(cols=vars(seurat_clusters))
  
  tab%>%   dplyr::filter(str_detect(cond, "2h")) %>% ggplot()+
    geom_col(aes(x=cond, y=n, fill=cond), position= position_dodge(w=1))+
    theme_bw()+facet_grid( rows=vars(rep), cols=vars(seurat_clusters))+ ggtitle('2hr')+
    scale_fill_manual(values = c("#1F97FF", "#FF7B15" ))+
   tab%>%   dplyr::filter(str_detect(cond, "30h")) %>% ggplot()+
    geom_col(aes(x=cond, y=n, fill=cond), position= position_dodge(w=1))+ ggtitle('30hr')+
    theme_bw()+facet_grid( rows=vars(rep), cols=vars(seurat_clusters), space = "free")+
    scale_fill_manual(values = c("#1F97FF", "#FF7B15" ))
  
  
two <-JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>% summarize(n=n())%>% ungroup()%>%
  dplyr::filter(str_detect(cond, "2h"))%>% #mutate(n=n())%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%

  mutate( time = 2,
          LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= as.numeric(binom.test(x=DMSO, 
                                  n=total, 
                                  p=prob)[3])
  )%>% as.data.table()

two1<- JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>% dplyr::filter(str_detect(cond, "2h"), rep == "rep1")%>%
  summarize(n=n())%>% ungroup()%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%
  mutate( time = 2,
          rep = "rep1",
          LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= as.numeric(binom.test(x=DMSO, 
                                  n=total, 
                                  p=prob)[3])
  )%>% as.data.table()


two2<- JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>%dplyr::filter(str_detect(cond, "2h"), rep == 'rep2')%>%
  summarize(n=n())%>% ungroup()%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%

  mutate( time = 2,
          rep = "rep2",
          LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= as.numeric(binom.test(x=DMSO, 
                                  n=total, 
                                  p=prob)[3])
  )%>% as.data.table()

thirty1<- JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>% dplyr::filter(str_detect(cond, "30h"), rep == "rep1")%>%
  summarize(n=n())%>% ungroup()%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%
  mutate( time = 30,
          rep = "rep1",
          LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= as.numeric(binom.test(x=DMSO, 
                                  n=total, 
                                  p=prob)[3])
  )%>% as.data.table()


thirty2<- JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>%dplyr::filter(str_detect(cond, "30h"), rep == 'rep2')%>%
  summarize(n=n())%>% ungroup()%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%
  mutate( time = 30,
          rep = "rep2",
          LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= as.numeric(binom.test(x=DMSO, 
                                  n=total, 
                                  p=prob)[3])
  )%>% as.data.table()



bind_rows(two1, two2)%>%  arrange(LFC) %>% 
  ggplot()+
  geom_point(aes(x=LFC, y=rep, fill=-log10(p), group= rep), position=position_dodge(w=0.3), size= 3.5, shape = 21, stroke=0.5)+
  geom_vline(xintercept = 0, color= "blue", size = 0.5, linetype=2 )+
  scale_fill_viridis_c(option='B', 'Binomial Test\n -log10(p)')+ facet_grid(rows=vars(cluster, time))+
  scale_fill_distiller(palette = 'Spectral', direction = -1, type = "div" )+
  theme_bw()+ylab("Clusters")+xlab('Differential Abundance\n log2 Fold Change (DMSO/V-1)')+theme(legend.position = 'bottom')
  
  bind_rows(two1, two2, 
            thirty1, thirty2
  )%>%  arrange(LFC) %>% 
  ggplot()+
  geom_point(aes(x=LFC, y=fct_rev(rep), fill=-log10(p), group= fct_rev(rep)), position=position_dodge(w=0.3), size= 1, shape = 21, stroke=0.5)+
  geom_vline(xintercept = 0, color= "blue", size = 0.5, linetype=2 )+
  scale_fill_viridis_c(option='B', 'Binomial Test\n -log10(p)')+ facet_grid(rows=vars(cluster), cols=vars(time))+
  scale_fill_distiller(palette = 'Spectral', direction = -1, type = "div" )+
  theme_bw()+ylab("Clusters")+xlab('Differential Abundance\n log2 Fold Change (DMSO/V-1)')+theme(legend.position = 'bottom') 
           



class(thirty$p)





JD_CAF1@meta.data%>% dplyr::select(seurat_clusters, cond, rep)%>% remove_rownames()%>% 
  mutate(cond = gsub('-pl2', '', cond))%>% group_by(cond, seurat_clusters)%>% summarize(n=n())%>% ungroup()%>%
  dplyr::filter(str_detect(cond, "2h"))%>% #mutate(n=n())%>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0)%>% 
  setNames(c("cluster", "DMSO", 'V1'))%>% ungroup()%>%
  # mutate(DMSO = DMSO+1,
  #         V1 = V1+1)%>%
  mutate( LFC = log2((DMSO/sum(DMSO))/ (V1/sum(V1))),
          total = DMSO+V1,
          sumDMSO=sum(DMSO),
          sumV1=sum(V1),
          prob = (sumDMSO)/(sumDMSO+sumV1))%>%
  rowwise()%>%
  mutate(p= binom.test(x=DMSO, 
                       n=total, 
                       p=prob)[3]
  )%>% as.data.table()


JD.CAF1_markers = FindAllMarkers(JD_CAF1, min.pct = 0.1, logfc.threshold = 0.25, only.pos = T)
JD.CAF1_top10_DEG = JD.CAF1_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) 

JD_CAF1 <- ScaleData(JD_CAF1, verbose = FALSE)
DoHeatmap(JD_CAF1, features = JD.CAF1_top10_DEG$gene,  size = 2) + NoLegend() + 
  theme(text = element_text(size = 5), legend.title = element_text(size=10), legend.text=element_text(size=5))

#######2hr######

JD_CAF1@meta.data <- JD_CAF1@meta.data %>% mutate(cond2 = gsub("-pl2", "", cond),
                                                  time = gsub("DMSO-|V1-", "", cond2))

JD_CAF1_2h <- subset(JD_CAF1, time == "2h")
JD_CAF1_2h <- Seurat::SplitObject(JD_CAF1_2h, split.by = "cond2")  


volcano_CAF1_2h <-  bind_rows(
  (JD_CAF1_2h$`DMSO-2h`[["RNA"]]@counts %>%   as.data.frame()%>%
     tibble::rownames_to_column(var = 'gene')%>%
     pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "DMSO")),
  
  (JD_CAF1_2h$`V1-2h`[["RNA"]]@counts %>%   as.data.frame()%>%
     tibble::rownames_to_column(var = 'gene')%>%
     pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "dTAG")))%>% as.data.table()


cells <- volcano_CAF1_2h[, .(count_tot = sum(count)), .(cell, cluster)
][, clusterN := uniqueN(cell), .(cluster)
]

volcano_CAF1_2h_p <- volcano_CAF1_2h[count > 0
][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
][order(gene, cluster), .(count_norm = list(count_norm)), 
  .(gene, cluster, clusterN = as.factor(cells[,setNames(clusterN, cell)][cell]))
][, data.table(map2_dfr(count_norm, as.numeric(as.character(clusterN)), function(x, cN){
  test <- nafill(x[1:cN], fill = 0)  
  rest <- nafill(c(0, unlist(count_norm[clusterN != cN], F, F))[2:(sum(as.numeric(levels(clusterN))) - cN + 1)], fill = 0)
  data.table(pvalue = t.test(test, rest)$p.value,
             logfold = log2(mean(test) / mean(rest)))
}), cluster = cluster), .(gene)
][, pvalue_cor := p.adjust(pvalue, "bonferroni"), .(cluster)]

volcano_early <- volcano_CAF1_2h_p %>%  mutate(sig = pvalue_cor<0.05,  updown = logfold >0.5 | logfold < -0.5)%>%dplyr::filter(cluster == "DMSO")%>%
  ggplot+
  geom_point(aes(y=-log10(pvalue_cor), x=logfold, fill= sig & updown), size= 1.5, alpha=0.8,shape =21, stroke=0.15)+
  theme_bw()+
  scale_fill_manual(values = c('grey',  'blue2', 'yellow2' ))+xlab('log2(DMSO/dTAG)')+ ggtitle('')+
  ylab('Bonferoni corrected T-test\n-log10(p-value)')+theme(legend.position = 'none')


volcano_CAF1_2h_p[logfold > 0.5, prank := frank(pvalue, ties.method = "first"), .(cluster)]




heatmap_early <- volcano_CAF1_2h[gene %in% volcano_CAF1_2h_p[(prank < 8 & pvalue_cor<0.05  & logfold>0.5)]$gene & count > 0, !"cluster"
][, .SD[cells[, .(cell, cluster)], on = .(cell)], .(gene)
][is.na(count), count := 0
][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
][, countz := (count_norm - mean(count_norm)) / sd(count_norm), .(gene)
][, cellorder := as.numeric(as.factor(cell)), .(cluster)
][, gene := factor(gene,
                   levels = dcast(.SD[,.(gene, cell, countz)], 
                                  gene ~ cell, value.var = "countz", fill = 0) %>%
                     {as.dist(as.matrix(data.frame(Dist((.[,!"gene"]), method = "manhattan"), 
                                                   row.names = .$gene)))} %>%
                     hclust("ward.D") %>% 
                     {.$labels[.$order]})
][abs(countz) > 1.5, countz := sign(countz) * 1.5

][] %>% ggplot()+
  geom_raster(aes(x = cellorder, y = gene, fill = countz)) +
  facet_grid(cols = vars(cluster),  scales = "free", space = "free") +
  scale_fill_gradient2(low = "magenta", mid = "black", high = "green") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), text = element_text(size = 9)) +
  coord_cartesian(expand = F)+ggtitle('')+theme(legend.position = 'none')

(volcano_early|heatmap_early)+plot_layout(widths = c(3,6))



JD_CAF1_30h <- subset(JD_CAF1, time == "30h")
JD_CAF1_30h <-SplitObject(JD_CAF1_30h, split.by = "cond2")

volcano_CAF1_30h <-  bind_rows(
  (JD_CAF1_30h$`DMSO-30h`[["RNA"]]@counts %>%   as.data.frame()%>%
     tibble::rownames_to_column(var = 'gene')%>%
     pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "DMSO")),
  
  (JD_CAF1_30h$`V1-30h`[["RNA"]]@counts %>%   as.data.frame()%>%
     tibble::rownames_to_column(var = 'gene')%>%
     pivot_longer(!gene,  names_to = "cell", values_to = "count")%>% mutate(cluster = "dTAG")))%>% as.data.table()


cells <- volcano_CAF1_30h[, .(count_tot = sum(count)), .(cell, cluster)
][, clusterN := uniqueN(cell), .(cluster)
]

volcano_CAF1_30h_p <- volcano_CAF1_30h[count > 0
][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
][order(gene, cluster), .(count_norm = list(count_norm)), 
  .(gene, cluster, clusterN = as.factor(cells[,setNames(clusterN, cell)][cell]))
][, data.table(map2_dfr(count_norm, as.numeric(as.character(clusterN)), function(x, cN){
  test <- nafill(x[1:cN], fill = 0)  
  rest <- nafill(c(0, unlist(count_norm[clusterN != cN], F, F))[2:(sum(as.numeric(levels(clusterN))) - cN + 1)], fill = 0)
  data.table(pvalue = t.test(test, rest)$p.value,
             logfold = log2(mean(test) / mean(rest)))
}), cluster = cluster), .(gene)
][, pvalue_cor := p.adjust(pvalue, "bonferroni"), .(cluster)]

late_volcano <- volcano_CAF1_30h_p %>%  mutate(sig = pvalue_cor<0.05,  updown = logfold >0.5 | logfold < -0.5)%>%dplyr::filter(cluster == "DMSO")%>%
  ggplot+
  geom_point(aes(y=-log10(pvalue_cor), x=logfold, fill= sig & updown), size= 1.5, alpha=0.8,shape =21, stroke=0.15)+
  theme_bw()+
  scale_fill_manual(values = c('grey',  'blue2', 'yellow2' ))+xlab('log2(DMSO/dTAG)')+ ggtitle('')+
  ylab('Bonferoni corrected T-test\n-log10(p-value)')+theme(legend.position = 'none')


volcano_CAF1_30h_p[logfold > 0.5, prank := frank(pvalue, ties.method = "first"), .(cluster)]

library(Rfast)

volcano_CAF1_30h_p %>% mutate(sig = pvalue_cor<0.05,  updown = logfold > 0.5) %>% dplyr::filter(sig, updown & cluster =='dTAG')%>%
  dplyr::select(gene)%>% mutate(new_gene = str_remove(gene, '.[0-9]{1,2}-.*'))%>% dplyr::select(new_gene)%>%  
  write_tsv(., '~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/30hrs_up.txt')
  

volcano_CAF1_2h_p %>% mutate(sig = pvalue_cor<0.05,  updown = logfold > 0.25) %>% dplyr::filter(sig, updown & cluster =='dTAG')%>%
  dplyr::select(gene)%>% mutate(new_gene = str_remove(gene, '.[0-9]{1,2}-.*'))%>% dplyr::select(new_gene)%>%  
  write_tsv(., '~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/2hrs_up.txt')

  
lateheatmap <- volcano_CAF1_30h[gene %in% volcano_CAF1_30h_p[(prank < 20 & pvalue_cor<0.05)]$gene & count > 0, !"cluster"
][, .SD[cells[, .(cell, cluster)], on = .(cell)], .(gene)
][is.na(count), count := 0
][, count_norm := log1p(count / cells[, setNames(count_tot, cell)][cell] * median(cells$count_tot))
][, countz := (count_norm - mean(count_norm)) / sd(count_norm), .(gene)
][, cellorder := as.numeric(as.factor(cell)), .(cluster)
][, gene := factor(gene,
                   levels = dcast(.SD[,.(gene, cell, countz)], 
                                  gene ~ cell, value.var = "countz", fill = 0) %>%
                     {as.dist(as.matrix(data.frame(Dist((.[,!"gene"]), method = "manhattan"), 
                                                   row.names = .$gene)))} %>%
                     hclust("ward.D") %>% 
                     {.$labels[.$order]})
][abs(countz) > 1.5, countz := sign(countz) * 1.5
  
] %>% ggplot()+
  geom_raster(aes(x = cellorder, y = gene, fill = countz)) +
  facet_grid(cols = vars(cluster),  scales = "free", space = "free") +
  scale_fill_gradient2(low = "magenta", mid = "black", high = "green") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), text = element_text(size = 9)) +
  coord_cartesian(expand = F)+ggtitle('')+theme(legend.position = 'none')

(late_volcano|lateheatmap)+plot_layout(widths=c(3,6))

library(data.table)
library(tidyverse)
library(Matrix)

rep2_unspliced <- {readRDS('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/unfiltvelo.rds')$unspliced}%>%
  data.table::as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(unspliced = sum(count))
  

rep2_spliced <- {readRDS('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/unfiltvelo.rds')$spliced}%>%
  as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(spliced = sum(count))

rep2_amb <- {readRDS('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/compiled/unfiltvelo.rds')$ambiguous}%>%
  data.table::as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(amb = sum(count))



rep1_unspliced <- {readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/unfiltvelo.rds')$unspliced}%>%
  as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(unspliced = sum(count))


rep1_spliced <- {readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/unfiltvelo.rds')$spliced}%>%
  as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(spliced = sum(count))

rep1_amb <- {readRDS('~/Desktop/post-doc/Experiments/exp137-VASA-CAF1-dTAG-V1-2hr-30hr/JD-VASA/compiled/unfiltvelo.rds')$ambiguous}%>%
  data.table::as.data.table()%>% remove_rownames() %>% pivot_longer(everything(), values_to = "count", names_to= "cell") %>% group_by(cell)%>%
  summarise(amb = sum(count))





rep2 <- full_join(rep2_spliced, rep2_unspliced, c('cell'='cell')) %>% mutate(sum_total = (spliced+unspliced))%>%
  mutate(Txn_activity = unspliced/(unspliced+spliced))%>% mutate(cond = str_split(cell, "\\_", simplify = T)[,1])%>% mutate(cond = str_replace(cond, "-pl2", ""), rep= "rep2")


rep2 <- full_join(rep2, rep2_amb, c('cell'='cell')) 


rep1 <- full_join(rep1_spliced, rep1_unspliced, c('cell'='cell')) %>% mutate(sum_total = (spliced+unspliced))%>%
  mutate(Txn_activity = unspliced/(unspliced+spliced))%>% mutate(cond = str_split(cell, "\\_", simplify = T)[,1], rep = "rep1")

rep1 <- full_join(rep1, rep1_amb, c('cell'='cell'))

bind_rows(rep2, rep1) %>%
  dplyr::filter(spliced >2500) %>%
  pivot_longer(cols = c("spliced", "unspliced", "amb")) %>%
  ggplot()+
  geom_boxplot(aes(x =  factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= log10(value), fill= cond), alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= log10(value), col= cond), alpha=0.5, size =0.75)+
  theme_bw()+
  xlab('')+
  facet_grid(cols=vars(factor(name, c("spliced", "unspliced", "amb"))))+
  ylab('log10 reads per cell')+
  scale_fill_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")

  


bind_rows(rep2, rep1)%>%  
  ggplot()+
  geom_violin(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")) , y= Txn_activity, fill= cond), alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Txn_activity, col= cond), size =0.75)+
  theme_bw()+
  xlab('')+
  ylab('Normalized Nascent Transription')+
  scale_fill_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")



stderror <- function(x) sd(x)/sqrt(length(x))
  
bind_rows(rep2, rep1) %>%
  dplyr::filter(spliced >2500) %>%
  mutate(sum_total = spliced + unspliced + amb) %>%
  pivot_longer(cols = c(sum_total)) %>%
  ggplot()+
  geom_boxplot(aes(x =  factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= (value), fill= cond), alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= (value), col= cond), alpha=0.5, size =0.75)+
  theme_bw()+
  xlab('')+
  ylab('log10 reads per cell')+
  ylim(0,1e5)+
  scale_fill_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")




library(rstatix)

bind_rows(rep2, rep1)%>%
  dplyr::filter(str_detect(cond, "2h"))%>% group_by(rep)%>%
  rstatix::wilcox_test(Txn_activity ~ cond)%>%  adjust_pvalue(method = "bonferroni")%>%
  add_significance("p.adj")

bind_rows(rep2, rep1)%>%
  dplyr::filter(str_detect(cond, "30h"))%>% group_by(rep)%>%
   rstatix::wilcox_test(Txn_activity ~ cond)%>%  adjust_pvalue(method = "bonferroni")

JD.combined <- readRDS("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/RPE1-dTAG-CAF1-2-30hr.rds")
meta <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname) %>% select(cell, seurat_clusters) 
all <- bind_rows(rep2, rep1)

left_join(all, meta, c('cell'='cell'))%>% drop_na()%>%
  ggplot()+
  geom_violin(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")) , y= Txn_activity, fill= cond), scale = 'width',  alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Txn_activity, col= cond), size =1.3)+
  theme_bw()+
  xlab('')+
  ylab('Normalized Nascent Transription')+
  scale_fill_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")+
  facet_grid(rows=vars(seurat_clusters))


S_phase_nascent <- full_join(meta, all, c('cell'='cell'))%>% drop_na()%>% 
  dplyr::filter(str_detect(cond, "2h"))%>%
  dplyr::filter(seurat_clusters %in% c(3,7,9))%>%
  ggplot()+
  geom_violin(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")) , y= Txn_activity, fill= cond),   trim=F, alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Txn_activity, col= cond), size =1.3)+
  theme_bw()+
  xlab('')+
  ylab('Normalized Nascent Transription')+
  scale_fill_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")+ggtitle('Sphase cells')+ylim(0,0.5)





full_join(meta, all, c('cell'='cell'))%>% drop_na()%>% 
  dplyr::filter(str_detect(cond, "2h"))%>%
  dplyr::filter(seurat_clusters %in% c(3,7,9))%>% 
    rstatix::t_test(Txn_activity ~ cond)


rest_nascent <- full_join(meta, all, c('cell'='cell'))%>% drop_na()%>% 
  dplyr::filter(str_detect(cond, "2h"))%>%
  dplyr::filter(seurat_clusters %in% c( 2, 4,5,8))%>%
  ggplot()+
  geom_violin(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Txn_activity, fill= cond), trim=F, alpha=0.1, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Txn_activity, col= cond), size =1.3)+
  theme_bw()+
  xlab('')+
  ylab('Normalized Nascent Transription')+
  scale_fill_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  theme(legend.position = "none")+ggtitle('G1, G2, Mitotic cells')+ylim(0,0.5)



full_join(meta, all, c('cell'='cell'))%>% drop_na()%>% 
  dplyr::filter(str_detect(cond, "2h"))%>%
  dplyr::filter(seurat_clusters %in% c(2,4 ,5,8))%>%
  rstatix::t_test(Txn_activity ~ cond)%>%  adjust_pvalue(method = "bonferroni")


S_phase_nascent | rest_nascent



JD.combined <- readRDS("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/RPE1-dTAG-CAF1-2-30hr.rds")

hist_genes <- rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., '-H3A|-H3B|-H3C|-H4A|-H4B|-H4C|-H2A|-H2C|-H2B|-H1'))
hist_genes <- pull(hist_genes, .)

JD.combined <- AddModuleScore(
  object = JD.combined,
  features = list(hist_genes),
  ctrl = 5,
  name = 'Histone'
)

JD.combined@meta.data$Histone1
library(RColorBrewer)
hist_a<- FeaturePlot(JD.combined, pt.size = 1.1, features = c("Histone1"), order = T)+coord_equal()
hist_a<-hist_a+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdYlBu")), na.value = "grey92")+ggtitle('Histone Gene Expression Score')


meta1 <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname) 

hist_B <- meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= Histone1, fill= cond2),width=1,   trim=F, alpha=0.33)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= Histone1, col= cond2),width=0.35, size =1, alpha=0.5)+

  theme_bw()+
  xlab('')+
  ylab('Histone Gene \n Expression Score')+

  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

(hist_a|hist_B)+plot_layout(widths = c(5,2.5))

library(rstatix)
meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         cond3 = paste0(Sphase, cond2))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>%
  rstatix::pairwise_t_test(Histone1 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")



rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., '-MCM|-CDC6|-SLBP|-H3A|-H3B|-H3C|-H4A|-H4B|-H4C|-H2A|-H2C|-H2B|-H1'))

#######S-phase_split#####


S_genes <- rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., '-MCM|-CDC6|-SLBP|-H3A|-H3B|-H3C|-H4A|-H4B|-H4C|-H2A|-H2C|-H2B|-H1'))
S_genes <- pull(S_genes, .)

JD.combined <- AddModuleScore(
  object = JD.combined,
  features = list(S_genes),
  ctrl = 5,
  name = 'S_genes'
)


hist_c <- FeaturePlot(JD.combined, pt.size = 1, features = c("S_genes1"), order = T)+coord_equal()
hist_c <-hist_c+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdYlBu")), na.value = "grey92")+ggtitle('Histone Gene Expression Score')

hist_c

meta1 <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=S_genes1), bins=100)

 meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% mutate(cond2 = gsub('-pl2', '', cond))%>%
  ggplot()+
  geom_violin(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Histone1, fill= cond2), scale = "width", trim=F, alpha=0.2, draw_quantiles = 0.5)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= Histone1, col= cond2), size =1.3)+
  theme_bw()+
  xlab('')+
  ylab('Histone Gene \n Expression Score')+
  facet_grid(cols = vars(S_genes1 >0.2))+
  scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% mutate(cond2 = gsub('-pl2', '', cond))%>%
  ggplot()+
  geom_violin(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= S_genes1, fill= cond2), scale = "width", trim=F, alpha=0.2, draw_quantiles = T)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= S_genes1, col= cond2), size =0.8)+
  theme_bw()+
  xlab('')+
  ylab('Histone Gene \n Expression Score')+
  facet_grid(cols = vars(S_genes1 >-0.25))+
  scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")





p53_hist_genes <- rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., 'CDKN1A|TP53I3|CD82|PMAIP|MDM2|SFN|DDB2|PPM1D|SESN1|-RRM2|ENSG00000026103.22-FAS'))
p53_hist_genes <- pull(p53_hist_genes, .)

JD.combined <- AddModuleScore(
  object = JD.combined,
  features = list(p53_hist_genes),
  ctrl = 5,
  name = 'p53_hist_genes'
)


hist_e <- FeaturePlot(JD.combined, pt.size = 1, features = c("p53_hist_genes1"), order = T)+coord_equal()
hist_e <-hist_e+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdYlBu")), na.value = "grey92")+ggtitle('TP53 histone dpeletion Score')

hist_e

meta1 <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=p53_hist_genes1), bins=100)





meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% mutate(Histone2  = Histone1 > 0 )%>%  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  ggplot()+
  geom_violin(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= p53_hist_genes1, fill= cond2), scale = "width", trim=F, alpha=0.2, draw_quantiles = 0.5)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= p53_hist_genes1, col= cond2), size =1.3)+
  theme_bw()+
  xlab('')+
  ylab('TP53 response to Histone depletion')+
  scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")


########





all_hist_genes <- rownames(JD.combined@assays$RNA[]) %>% as.data.table()%>%dplyr::filter(str_detect(., 'ENSG00000103197.18-TSC2|TP73|SESN2|CCNG2|BBC3|SERPINE1|CDKN1A|TP53I3|CD82|PMAIP|MDM2|SFN|DDB2|PPM1D|SESN1|-RRM2|ENSG00000026103.22-FAS'))
all_hist_genes <- pull(all_hist_genes, .)

JD.combined <- AddModuleScore(
  object = JD.combined,
  features = list(all_hist_genes),
  ctrl = 5,
  name = 'all_hist_genes'
)

hist_g <- FeaturePlot(JD.combined, pt.size = 1, features = c("all_hist_genes1"), order = T)+coord_equal()
hist_g <-hist_g+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('Histone depletion response Score')

hist_g

meta1 <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=non_p53_hist_genes1), bins=100)


install.packages("ggpol")
library(ggpol)

hist_responseA <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
   dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= p53_hist_genes1, fill= cond2), trim=F, alpha=0.2, draw_quantiles = T)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= p53_hist_genes1, col= cond2), size =1)+
  theme_bw()+
  xlab('')+
  ylab('helleday p53 histone')+

  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))

hist_responseB <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= all_hist_genes1, fill= cond2),width=1,   trim=F, alpha=0.33)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= all_hist_genes1, col= cond2),width=0.35, size =1, alpha=0.5)+
  theme_bw()+
  xlab('')+
  ylab('Perturbed Histone synthesis reponse')+

  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+
  theme(legend.position = 'none')


hist_B <- meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= Histone1, fill= cond2),width=1,   trim=F, alpha=0.33)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= Histone1, col= cond2),width=0.35, size =1, alpha=0.5)+
  theme_bw()+
  xlab('')+
  ylab('Histone Gene set Score')+

  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")




meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         cond3 = paste0(Sphase, cond2))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>%
  rstatix::pairwise_t_test(Histone1 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")


meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         cond3 = paste0(Sphase, cond2))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>%
  rstatix::pairwise_t_test(all_hist_genes1 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")

hist_B|hist_responseB 

meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% mutate(Histone2  = Histone1 > 0 )%>%  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%#, Histone2)%>%
  ggplot()+  
  geom_boxjitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= all_hist_genes1, fill= cond2),
               jitter.shape = 21, jitter.color = NA, 
               jitter.height = 0, jitter.width = 0.04,
               outlier.color = NA, errorbar.draw = TRUE) +
  scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
  theme_bw()+  ggtitle('All cells')+
  xlab('')+
  ylab('Histone depletion response score')+theme(legend.position = "none")
  
  hist_responseB <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% mutate(Histone2  = Histone1 > 0 )%>%  dplyr::filter(line =="RPE1 CAF1-dTAG", Histone2)%>%
    ggplot()+
    geom_violin(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= all_hist_genes1, fill= cond2),
                width=0.35, size =0.8, alpha=0.5)+
    geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= all_hist_genes1, col= cond2),
                width=1,   trim=T, alpha=0.33,size =0.8)+
    theme_bw()+
    ggtitle('S-phase cells')+
  xlab('')+
    ylab('Histone depletion response score')+
    #facet_grid(cols = vars(Histone1 >0))+
    scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
    scale_color_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")
  
  

hist_g/(hist_responseA|hist_responseB)


ATAC_list <- fread('~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/join_overlap_intersect (1).csv')%>% 
  dplyr::filter(gene_type == "protein_coding")%>%
  distinct(gene_name, .keep_all = T)%>%
  arrange(-log2_M)%>% select(gene_id)%>% slice_head(n=500)%>% pull()

ATAC_list

ATAC_list <-  rownames(JD.combined@assays$RNA[]) %>% 
              as.data.table()%>%
              dplyr::filter(str_detect(.,paste(ATAC_list, collapse = "|")))%>% pull()

JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")
JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(ATAC_list),
  ctrl = 5,assay = "RNA",
  bins=20,
  name = 'ATAC'
)


ATAC_a <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("ATAC1"), order = T)+coord_equal()
ATAC_a <-hist_g+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('ATACa')

ATAC_a

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=ATAC1), bins=100)

ATAC_response <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% mutate(Histone2  = Histone1 > 0 )%>%  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  ggplot()+
  geom_violin(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= ATAC1, fill= cond2), scale = "width", trim=F, alpha=0.2, draw_quantiles = 0.5)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= ATAC1, col= cond2), size =0.8)+
  theme_bw()+
  ggtitle('Top 500 log2FOLD ATAC genes in H3k27me3')+
  xlab('')+
  ylab('Gene expression module score')+ ylim(-0.02,0.02)+
  #facet_grid(cols = vars(Histone1 >0))+
  scale_fill_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF","#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

ATAC_response


meta1 %>% mutate(cond2 = gsub('-pl2', '', cond))%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%  dplyr::filter(str_detect(cond2, '30h'))%>%
  rstatix::anova_test(ATAC1 ~ cond2)%>%  adjust_pvalue(method = "BH")%>%
  add_significance("p.adj")
  

meta1 %>% mutate(cond2 = gsub('-pl2', '', cond))%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%  dplyr::filter(str_detect(cond2, '2h'))%>% 
  rstatix::anova_test(ATAC1 ~ cond2)%>%  adjust_pvalue(method = "BH")%>%
  add_significance("p.adj")


nonrep_h3 <- rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(., '-H3-3')) %>% pull()

rep_h3 <- rownames(JD.combined@assays$RNA[]) %>%
  as.data.table()%>%
  dplyr::filter(str_detect(., '-H3C|-H3-2|H3Y|H3-2')) %>% 
  pull()

JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(nonrep_h3),
  ctrl = 5,assay = "RNA",
  bins=20,
  name = 'H3.3'
)


JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(rep_h3),
  ctrl = 5,assay = "RNA",
  bins=20,
  name = 'H3.1-2'
)


non_rep_a <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("H3.31"), order = T)+coord_equal()
non_rep_a <- non_rep_a+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('H3.3 expression score')

non_rep_a

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=H3.31), bins=100)

nonrepA_response <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= H3.31, fill= cond2),
               alpha=0.5)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= H3.31, col= cond2),
         trim=T, alpha=0.33,size =0.8)+
  theme_bw()+
  ylab('H3.3 expression score')+
  xlab('')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

nonrepA_response 


rep_a <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("H3.1.21"), order = T)+coord_equal()
rep_a <- rep_a+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('H3.3 expression score')

rep_a

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=H3.31), bins=100)

repA_response <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  mutate(cond2 = gsub('-pl2', '', cond),
        Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-2h", "V1-2h"))%>% 
  ggplot()+
  geom_violin(aes( x = paste0(Sphase, cond2), y= H3.1.21, fill= cond2),
           alpha=0.5)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= H3.1.21, col= cond2),
               trim=T, alpha=0.33,size =0.8)+
  theme_bw()+
  ylab('H3.1-H3.2 expression score')+
  xlab('')+
  scale_fill_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

repA_response|nonrepA_response 
 



meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         cond3 = paste0(Sphase, cond2))%>%
  dplyr::filter(cond2 %in% c("DMSO-30h", "V1-30h"))%>%
  rstatix::pairwise_t_test(H3.31 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")



meta1 %>% dplyr::filter(line =="RPE1 CAF1-dTAG")%>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         cond3 = paste0(Sphase, cond2))%>%
  dplyr::filter(cond2 %in% c("DMSO-30h", "V1-30h"))%>%
  rstatix::pairwise_t_test(H3.1.21 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")


nonrepA_response30 <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  dplyr::filter(cond2 %in% c("DMSO-30h", "V1-30h"))%>% 
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, cond2), y= H3.31, fill= cond2),
              alpha=0.5)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= H3.31, col= cond2),
              trim=T, alpha=0.33,size =0.8)+
   theme_bw()+
  ylab('H3.3 expression score')+
  xlab('')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

nonrepA_response 


rep_a <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("H3.1.21"), order = T)+coord_equal()
rep_a <- rep_a+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('H3.3 expression score')

rep_a

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)

meta1%>% ggplot()+geom_histogram(aes(x=H3.31), bins=100)

repA_response30 <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(cond2 %in% c("DMSO-30h", "V1-30h"))%>% 
  ggplot()+
  geom_violin(aes( x = paste0(Sphase, cond2), y= H3.1.21, fill= cond2),
              alpha=0.5)+
  geom_jitter(aes(x = paste0(Sphase, cond2), y= H3.1.21, col= cond2),
              trim=T, alpha=0.33,size =0.8)+
   theme_bw()+
  ylab('H3.1-H3.2 expression score')+
  xlab('')+
  scale_fill_manual(values = c("#1F97FF",  "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")

repA_response30|nonrepA_response30

G0 <-rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(., '4.14-LRP|5.3-GPR37|MXD4|GABARAP1|IGFBP5|PDE5A|CDC42EP2|NEAT1')) %>% pull()




JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(G0),
  ctrl = 11,assay = "RNA",
  bins=20,
  name = 'G0_2'
)


G0umap <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("G0_21"), order = T)+coord_equal()
G0umap <- G0umap+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('Quiescence Score')

G0umap

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)



G0_response <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  ggplot()+
  geom_violin(aes(x = paste0(cond2), y= G0_21, fill= cond2),
              alpha=0.5)+
  geom_jitter(aes(x = paste0(cond2), y= G0_21, col= cond2),
              trim=T, alpha=0.7,size =0.8)+
   theme_bw()+
  ylab('Quiescence score')+
  xlab('')+
   ylim(-0.2,0.8)+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF",  "#FF7B15","#1F97FF","#FF7B15"  ))+theme(legend.position = "none")
(G0umap|G0_response)+plot_layout(widths = c(6,3))

KI67 <-rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(., 'KI67')) %>% pull()
JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(KI67),
  ctrl = 11,assay = "RNA",
  bins=20,
  name = 'KI67'
)


meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)


meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")
         )%>%
  dplyr::filter(line =="RPE1 CAF1-dTAG")%>%
  ggplot()+
  geom_point(aes(x=KI671,y=G0_21,color=cond2))+
  theme_bw()+
  ylab('Quiescence score')+
  xlab('KI67')+
  facet_grid(cols = vars(time, treatment))+
  scale_fill_manual(values = c("#1F97FF","#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")



meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>% 
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(seurat_clusters %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG"),
         cond3 = paste0(time, treatment))%>%
  rstatix::pairwise_t_test(G0_21 ~ cond3)%>%  adjust_pvalue(method = "bonferroni")


fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv') %>%
  mutate(log2_M = log2(Mature_ATAC_dTAG/Mature_ATAC_DMSO))%>%
  ggplot()+geom_point(aes(x=log2_N, y=log2_M))

fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv')%>%
mutate(Z_K27 = scale(H3K27me3_score),
       Z_K9 = scale(H3K9me3_score),
       Z_div = Z_K9-Z_K27,
       mod =  Z_K9 < 1 & Z_K27 > 0.5)%>%
  dplyr::filter(gene_type == "protein_coding") %>%
  ggplot()+geom_bin_2d(aes(x=Z_K27, y= Z_K9))+theme_bw()+
  facet_grid(cols=vars(mod))+
  scale_fill_viridis_c(option="B", trans="log2")

k27_list <- fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv')%>%
mutate(Z_K27 = scale(H3K27me3_score),
       Z_K9 = scale(H3K9me3_score))%>%
  dplyr::filter(gene_type == "protein_coding" & Z_K9 < 0 & Z_K27 > 0 & Z_K27< 4)%>%

  distinct(gene_name, .keep_all = T)%>%
  arrange(-log2_N)%>% select(gene_id)%>% slice_head(n=200)%>% 
  pull()

k27_list

k27_list <-  rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(k27_list, collapse = "|")))%>% pull()

JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")
JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(k27_list),
  ctrl = 100,assay = "RNA",
  bins=10,
  name = 'k27_list'
)

ATAC_a <- FeaturePlot(JD_CAF1, pt.size = 1, features = c("k27_list1"), order = T)+coord_equal()
ATAC_a <-ATAC_a+scale_colour_gradientn(colours = rev(brewer.pal(10,"RdBu")), na.value = "grey92")+ggtitle('ATACa')

ATAC_a

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)




k27_list <- fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv')%>%
  mutate(Z_K27 = scale(H3K27me3_score),
         Z_K9 = scale(H3K9me3_score))%>%
  dplyr::filter(gene_type == "protein_coding" & Z_K9 < 0 & Z_K27 > 0 & Z_K27< 4)%>%
  distinct(gene_name, .keep_all = T)%>%
  arrange(-log2_N)%>% select(gene_id)%>% slice_head(n=200)%>% 
  pull()

k27_list <-  rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(k27_list, collapse = "|")))%>% pull()




k9_list <- fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv')%>%
  mutate(Z_K27 = scale(H3K27me3_score),
         Z_K9 = scale(H3K9me3_score))%>%
  dplyr::filter(gene_type == "protein_coding" & Z_K27 < 0 & Z_K9 > 0)%>%
  distinct(gene_name, .keep_all = T)%>%
  arrange(-log2_N)%>% select(gene_id)%>% slice_head(n=200)%>% 
  pull()



k9_list <-  rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(k9_list, collapse = "|")))%>% pull()


set.seed(100)
random_list <- fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv')%>%
  mutate(Z_K27 = scale(H3K27me3_score),
         Z_K9 = scale(H3K9me3_score))%>%
  dplyr::filter(gene_type == "protein_coding")%>%
  distinct(gene_name, .keep_all = T)%>%
  select(gene_id)%>% slice_sample(n=200)%>% 
  pull()



random_list <-  rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(random_list, collapse = "|")))%>% pull()




JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")


JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(k9_list),
  ctrl = 50,assay = "RNA",
  bins=5,
  name = 'k9_list'
)

JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(k27_list),
  ctrl = 50,assay = "RNA",
  bins=5,
  name = 'k27_list'
)

JD_CAF1 <- AddModuleScore(
  object = JD_CAF1,
  features = list(random_list),
  ctrl = 50,assay = "RNA",
  bins=5,
  name = 'random'
)

meta1 <- JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname)


ATAC_response_k9 <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  ggplot()+
  geom_boxplot(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= k9_list1, fill= cond2), alpha=0.2, outlier.shape = NA)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= k9_list1, col= cond2), size =0.8)+
  theme_bw()+
  ggtitle('Top 200 log2FOLD ATAC genes in \nH3k9me3')+
  xlab('')+
  ylab('Gene expression module score')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")



ATAC_response_k27 <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  ggplot()+
  geom_boxplot(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= k27_list1, fill= cond2), alpha=0.2, outlier.shape = NA)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= k27_list1, col= cond2), size =0.8)+
  theme_bw()+
  ggtitle('Top 200 log2FOLD ATAC genes in \nH3k27me3')+
  xlab('')+
  ylab('Gene expression module score')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")



ATAC_response_random <- meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  ggplot()+
  geom_boxplot(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= random1, fill= cond2), alpha=0.2, outlier.shape = NA)+
  geom_jitter(aes(x = factor(cond2, c("DMSO-2h", "V1-2h", "DMSO-30h", "V1-30h")), y= random1, col= cond2), size =0.8)+
  theme_bw()+
  ggtitle('Random 200 genes')+
  xlab('')+
  ylab('Gene expression module score')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15","#FF7B15"  ))+theme(legend.position = "none")


ATAC_response_random |ATAC_response_k9| ATAC_response_k27   

meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  rstatix::pairwise_t_test(random1 ~ cond2)%>%  adjust_pvalue(method = "bonferroni")

meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  rstatix::pairwise_t_test(k27_list1 ~ cond2)%>%  adjust_pvalue(method = "bonferroni")

meta1 %>% mutate(cond2 = gsub('-pl2', '', cond)) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG",
                str_detect(cond2, "30h"))%>%
  rstatix::pairwise_t_test(k9_list1 ~ cond2)%>%  adjust_pvalue(method = "bonferroni")




HuHist <- fread("~/Desktop/post-doc/CAF1-DNA-replication/HuHist_txt") %>% pull()

library(data.table)
HuHist_list <- rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(HuHist, collapse = "|")))%>% pull()
var_list1 <- c("ident", "cond", HuHist_list) 

JD.combined <- readRDS("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/RPE1-dTAG-CAF1-2-30hr.rds")
JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")
JD_WT <- subset(x= JD.combined, line == "RPE1 WT")

MS_hits <- fread("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/Changing genes full proteome_3outof3_noG1_summedintensities.csv") %>%
  select(`Up in dTAG`) %>% pull()


MS_hits_list <- rownames(JD.combined@assays$RNA[]) %>% 
  as.data.table()%>%
  dplyr::filter(str_detect(.,paste(MS_hits, collapse = "|")))%>% pull()


var_list = c("ident", HuHist_list)
FetchData(JD_WT, var_list ) %>% 
rownames_to_column(var = "cell")%>%
  mutate(Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
   pivot_longer(starts_with("ENSG"))%>% 
  group_by(Sphase, name) %>%
    summarize(mean = median(value)) %>%
     pivot_wider(names_from = Sphase, values_from = mean) %>%
     mutate(enrichment = Sphase / `non-Sphase`,
            name  = gsub(".*\\.","",name))%>%
  dplyr::filter(`non-Sphase` > 0 | Sphase > 0) %>%
  ggplot()+
  geom_point(aes(y=Sphase, x = `non-Sphase`, fill= enrichment), shape =21, size =4)+
  ggrepel::geom_text_repel(aes(y=Sphase, x = `non-Sphase`, label = name), size =2.5)+
  geom_abline(intercept = 0, slope = 2, col="blue", linetype=2, size =1)+
  scale_fill_viridis_c(option ="B", trans= "log2")+theme_bw()+coord_fixed(ratio=0.8)+
  ggtitle('Histone expression')






FetchData(JD_WT, var_list ) %>% 
  rownames_to_column(var = "cell")%>%
  mutate(Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  pivot_longer(starts_with("ENSG"))%>% 
  group_by(Sphase, name) %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = Sphase, values_from = mean) %>%
  mutate(enrichment = Sphase / `non-Sphase`) %>%
  arrange(-enrichment)  %>%
  setNames(c("Gene", "MeanSphase", "MeanNonSphase", "enrichment")) %>%
  write.tsv(., "~/Desktop/post-doc/CAF1-DNA-replication/Sphase_histones.txt", quote = F, row.names = F, sep = "\t") 


var_list1 <- c("ident", "cond", HuHist_list) 

FetchData(JD_CAF1, var_list1 ) %>% 
  rownames_to_column(var = "cell")%>%
  pivot_longer(starts_with("ENSG"))%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(str_detect(name, "MACRO|H2AZ|H3-3|H2AX|6-H1-10")) %>%
  dplyr::filter(time == "2hr") %>%
  mutate(name  = gsub(".*\\.","",name)) %>%
  ggplot()+
  geom_violin(aes(x = paste0(Sphase, treatment), y= value, fill= treatment),
              alpha=0.5)+
  geom_jitter(aes(x = paste0(Sphase, treatment), y= value, col= treatment),
              trim=T, alpha=0.7,size =0.8)+
  theme_bw()+
  ylab('Expression')+
  xlab('')+
  facet_wrap(vars(name), nrow = 2, scales = )+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
  scale_color_manual(values = c("#1F97FF",  "#FF7B15","#1F97FF","#FF7B15"  ))+
  theme(legend.position = "none")


FetchData(JD_CAF1, var_list1 ) %>% 
  rownames_to_column(var = "cell")%>%
  pivot_longer(starts_with("ENSG"))%>%
  mutate(Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  distinct(cell, Sphase, time, treatment)%>%
  group_by( Sphase, time, treatment)%>%
  summarize(n=n())



diff <- FetchData(JD_CAF1, var_list1 ) %>% 
  rownames_to_column(var = "cell")%>%
  pivot_longer(starts_with("ENSG"))%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(Sphase == "Sphase")%>%
  group_by(name, treatment) %>%
  summarize(mean = mean(value)) %>%
  ungroup()%>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  dplyr::filter(DMSO > 0.1 | dTAG > 0.1) %>%
  mutate(diff = dTAG/DMSO)%>%
  arrange(diff)
  

absolute <- FetchData(JD_WT, var_list1 ) %>% 
  rownames_to_column(var = "cell")%>%
  mutate(Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"))%>%
  pivot_longer(starts_with("ENSG"))%>% 
  group_by(Sphase, name) %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = Sphase, values_from = mean) %>%
  mutate(enrichment = Sphase / `non-Sphase`)#%>%
  dplyr::filter(`non-Sphase` > 0.1 | Sphase > 0.1)
  
  
  
  absodiff <- full_join(diff, absolute, by = "name") %>%
    mutate(name  = gsub(".*\\.","",name)) %>%
    mutate(name  = gsub(".*\\-","",name))
  
  
  absodiff %>%
  ggplot()+
  geom_abline(intercept = 0, slope = 1.98, col="grey75", linetype=2, size =1)+
  ggrepel::geom_text_repel(aes(y=Sphase, x = `non-Sphase`, label = name), 
                           size =3,
                           #nudge_x = 0.1,
                           #nudge_y = 0.1,
                           segment.curvature = 0.2,
                           segment.ncp = 2,
                           segment.angle = 15,
                           max.time = 5, 
                           max.iter = 10e5,
                           force=30,
                           direction = c("both", "y", "x"),
                           min.segment.length = 1)+
  geom_point(aes(y=Sphase, x = `non-Sphase`, fill= diff, size = enrichment),  shape =21)+
  scale_fill_gradientn( colors = c( "#1F97FF", "#5293FF","#A9C0FA", "white", "#FF7B15"
                                    ),
                        limits = c(0.25, 1.25), oob = scales::squish,
                        name = "dTAG/DMSO \n fold change")+
  theme_bw()+coord_fixed(ratio=0.5)+
  ggtitle('Histone expression')+
  xlab("Mean expression outside of S-phase")+
  ylab("Mean S-phase exression")

  
  absodiff %>%
  
    dplyr::filter(`non-Sphase` > 0.1 | Sphase > 0.1)%>%
    ggplot()+
    geom_point(aes(x=DMSO, y = dTAG, fill= log2(diff)), size =3, shape =21)+
    ggrepel::geom_text_repel(aes(x=DMSO, y = dTAG, label = name), size =2.5)+
    scale_fill_distiller(palette = "RdBu", direction = -1)+theme_bw()+coord_fixed(ratio=0.8)+
    ggtitle('Histone expression')

absodiff %>% 
    arrange(-diff) %>%slice_head(n = 30)%>%select(name) %>% pull()
  
  
    ggplot()+
    geom_point(aes(x=diff, y=enrichment))
    
library(tidyverse)
library(Seurat) 
    

FetchData(JD_CAF1, var_list1 ) %>% rownames_to_column(var= "id") %>% select(id, ident, cond) %>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG"))%>%
  dplyr::filter(time == "2hr")%>%
  write_tsv("~/Desktop/JD_VASA.txt")
 
    
    
    FetchData(JD_CAF1, var_list1 ) %>% 
      rownames_to_column(var = "cell")%>%
      pivot_longer(starts_with("ENSG"))%>%
      mutate(cond2 = gsub('-pl2', '', cond),
             Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
             time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
             treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
      dplyr::filter(str_detect(name, ".2-H2BC8|14.1-H3C1|H4C2|NSG00000183598.3-H3C13|H2BC14|2AC15|1-H2AC15|1-H2AC4|H2AC7")) %>%
      dplyr::filter(time == "2hr") %>%
      mutate(name  = gsub(".*\\.","",name),
             cond3 = paste0(Sphase, "_", treatment)) %>%
            ggplot()+
      geom_violin(aes(x = paste0(Sphase, treatment), y= value, fill= treatment),
                  alpha=0.5, scale = "width", trim = T)+
      geom_jitter(aes(x = paste0(Sphase, treatment), y= value, col= treatment),
                   alpha=0.7,size =0.8)+
      theme_bw()+
      ylab('Expression')+
      xlab('')+
      facet_wrap(vars(name), nrow = 2, scales = "free")+
      scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
      scale_color_manual(values = c("#1F97FF",  "#FF7B15","#1F97FF","#FF7B15"  ))+
      theme(legend.position = "none")
    
    
    
  
    FetchData(JD_CAF1, var_list1 ) %>% 
      rownames_to_column(var = "cell")%>%
      pivot_longer(starts_with("ENSG"))%>%
      mutate(cond2 = gsub('-pl2', '', cond),
             Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
             time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
             treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
      dplyr::filter(str_detect(name, ".2-H2BC8|14.1-H3C1|H4C2|NSG00000183598.3-H3C13|H2BC14|2AC15|1-H2AC15|1-H2AC4|H2AC7")) %>%
      dplyr::filter(time == "2hr") %>%
      mutate(name  = gsub(".*\\.","",name),
             cond3 = paste0(Sphase, "_", treatment)) %>%
      ggplot()+
      geom_violin(aes(x = paste0(Sphase, treatment), y= value, fill= treatment),
                  alpha=0.5, scale = "width", trim = T)+
      geom_jitter(aes(x = paste0(Sphase, treatment), y= value, col= treatment),
                  alpha=0.7,size =0.8)+
      theme_bw()+
      ylab('Expression')+
      xlab('')+
      facet_wrap(vars(name), nrow = 2, scales = "free")+
      scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
      scale_color_manual(values = c("#1F97FF",  "#FF7B15","#1F97FF","#FF7B15"  ))+
      theme(legend.position = "none")
    
    FetchData(JD_CAF1, var_list1 ) %>% 
      rownames_to_column(var = "cell")%>%
      pivot_longer(starts_with("ENSG"))%>%
      mutate(cond2 = gsub('-pl2', '', cond),
             Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
             time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
             treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
      dplyr::filter(str_detect(name, "MACRO|H2AX|H2AZ|1-H2AC4|9-H2AC6|2-H2AC7|NSG00000164508.4-H2AC1|5-H2AC21")) %>%
      dplyr::filter(time == "2hr") %>%
      mutate(name  = gsub(".*\\.","",name),
             cond3 = paste0(Sphase, "_", treatment)) %>%
      ggplot()+
      geom_violin(aes(x = paste0(Sphase, treatment), y= value, fill= treatment),
                  alpha=0.5, scale = "width", trim = T)+
      geom_jitter(aes(x = paste0(Sphase, treatment), y= value, col= treatment),
                  alpha=0.7,size =0.8)+
      theme_bw()+
      ylab('Expression')+
      xlab('')+
      facet_wrap(vars(name), ncol = 3, scales = "free")+
      scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
      scale_color_manual(values = c("#1F97FF",  "#FF7B15","#1F97FF","#FF7B15"  ))+
      theme(legend.position = "none")
    
    
  
    
    
MS_hits <- fread("~/Desktop/post-doc/Experiments/exp140_VASA_CAF1_dTAG2_30hr_rep2/Changing genes full proteome_3outof3_noG1_summedintensities.csv") %>%
      select(`Up in dTAG`) %>% pull()
    
MS_hits_list <- rownames(JD.combined@assays$RNA[]) %>% 
      as.data.table()%>%
      dplyr::filter(str_detect(.,paste(MS_hits, collapse = "|")))%>% pull()
    
    
var_list2 = c("ident","cond", MS_hits_list)
 
FetchData(JD_CAF1, var_list2 ) %>% 
  rownames_to_column(var = "cell")%>%
  pivot_longer(starts_with("ENSG"))%>%
  mutate(cond2 = gsub('-pl2', '', cond),
         Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  mutate(name  = gsub(".*\\.","",name)) %>%
  mutate(name  = gsub(".*\\.","",name)) %>%
  mutate(name  = gsub(".*\\.","",name)) %>%
  dplyr::filter(str_detect("DHRS2", name))%>%
  
  group_by(name) %>%
  mutate(Zvalue = scale(value)) %>%  group_by(name, time, treatment) %>%
  summarise(expression = mean(Zvalue)) %>%
  dplyr::filter(time == "30hr") %>%
  ggplot()+geom_raster(aes(y=treatment, x=name, fill= expression))+
  scale_fill_viridis_c(option="H")



     
FetchData(JD_CAF1, var_list1 ) %>% 
      rownames_to_column(var = "cell")%>%
      pivot_longer(starts_with("ENSG"))%>%
      mutate(cond2 = gsub('-pl2', '', cond),
             Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
             time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
             treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
      dplyr::filter(str_detect(name, "MACRO|H2AZ|H3-3|H2AX|6-H1-10")) %>%
      dplyr::filter(time == "2hr") %>%
      mutate(name  = gsub(".*\\.","",name),
             cond3 = paste0(Sphase, "_", treatment)) %>%
      select(name, value, cond3) %>%
      group_by(name) %>%
      rstatix::pairwise_t_test(value ~ cond3) %>% rstatix::adjust_pvalue(method = "bonferroni") %>% write_tsv("~/Desktop/nonrepli_hist_stats.txt")    
      
      
    
    FetchData(JD_CAF1, var_list1 ) %>% 
      rownames_to_column(var = "cell")%>%
      pivot_longer(starts_with("ENSG"))%>%
      mutate(cond2 = gsub('-pl2', '', cond),
             Sphase = if_else(ident %in% c(3,7,9), true = "Sphase", false = "non-Sphase"),
             time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
             treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
      dplyr::filter(str_detect(name, ".2-H2BC8|14.1-H3C1|H4C2|NSG00000183598.3-H3C13|H2BC14|2AC15|1-H2AC15|1-H2AC4|H2AC7")) %>%
      dplyr::filter(time == "2hr") %>%
      mutate(name  = gsub(".*\\.","",name),
             cond3 = paste0(Sphase, "_", treatment)) %>%
      select(name, value, cond3) %>%
      group_by(name) %>%
      rstatix::pairwise_t_test(value ~ cond3) %>% rstatix::adjust_pvalue(method = "bonferroni") %>% write_tsv("~/Desktop/repli_hist_stats.txt")
    
    
    
meta1 %>%
  mutate(cond2 = gsub('-pl2', '', cond),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG"
                ) %>%
  group_by(seurat_clusters, time) %>%
  mutate(total =  n()) %>%
  group_by(cond2, seurat_clusters, time, treatment, total) %>%
  summarise(n= (n()/total))%>% distinct(seurat_clusters, time, treatment, n) %>%
  ggplot()+
  geom_col(aes(x=seurat_clusters, y= n, fill= treatment ) )+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
  facet_grid(rows= vars(time))+theme_bw()+
  xlab('Clusters')+
  ylab('Fraction of cells')

meta1 %>%
  mutate(cond2 = gsub('-pl2', '', cond),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG"
  ) %>%
  group_by(time) %>%
  mutate(total =  n()) %>%
  group_by(cond2, seurat_clusters, time, treatment, total) %>%
  summarise(n= (n()/total))%>% distinct(seurat_clusters, time, treatment, n) %>%
  ggplot()+
  geom_col(aes(x=seurat_clusters, y= n*100, fill= treatment ) )+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
  facet_grid(rows= vars(time))+theme_bw()+
  xlab('Clusters')+
  ylab('Percentage of cells')






library(tidyverse)
library(data.table)
absodiff
HuHist <- fread("~/Desktop/post-doc/CAF1-DNA-replication/HuHist_txt") %>% pull()

hist_genes <- fread("~/Desktop/post-doc/ForkDOT/NatureMethods/Rebuttal/speed/Homo_sapiens.GRCh38.90.gtf.tsv.gz") %>%
  group_by(type) %>% 
  dplyr::filter(type == 'transcript',
                transcript_biotype == "protein_coding",
                str_detect(gene_id ,paste(HuHist, collapse = "|"))) 



hist_genes %>% mutate(start = floor(start/1.6)*1e6,
                      end = ceiling(end/1.5e6)*1e6) %>%
  group_by(chr, start) %>%
  summarize(n=n()) %>%
  arrange(-n)

absodiff %>% ggplot()+geom_histogram(aes(x=diff), bins =50)


absodiff <- full_join(diff, absolute, by = "name") %>%
  mutate(name  = gsub(".*\\.","",name)) %>%
  mutate(name  = gsub(".*\\-","",name)) %>%
 setNames(c("gene_name", "DMSO",  "dTAG",  "diff", "Sphase", "nonSphase", "enrichment"))

atac_mature <- fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv') %>%
  mutate(histgene = if_else( str_detect(gene_id ,paste(HuHist, collapse = "|")), true = "histone", false = "non-histone")) %>%
  dplyr::filter(gene_type == "protein_coding",
                Mature_ATAC_DMSO > 20)

fread('~/Desktop/post-doc/CAF1-DNA-replication/merged_mature_clone8_h3k27me3_genes.csv') %>%
  dplyr::fi


mCAF1 <- left_join(atac_mature, absodiff, by= "gene_name") %>%
  mutate(CAF1hist = diff <0.8) %>%
  mutate_at(vars(CAF1hist), ~replace_na(., FALSE)) %>%
  group_by(histgene, CAF1hist) %>%
  slice_sample(n=200) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, fill =fct_rev(paste(histgene, CAF1hist))), alpha=0.1)+
  geom_jitter(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, col= fct_rev(paste(histgene, CAF1hist))), alpha=0.8)+
  theme_bw()+
  ylab('log2(dTAG/DMSO) Mature')+
  scale_fill_manual(values = c("grey20", "grey40","grey80" ))+
  scale_color_manual(values = c("grey20", "grey40","grey80" ))+
  ggtitle("CAF1 dependent histone expression")+
  ylim(-1,3)+
  theme(legend.position = "bottom")


mSphase <- left_join(atac_mature, absodiff, by= "gene_name") %>%
  mutate(CAF1hist = enrichment > 2.5) %>%
  mutate_at(vars(CAF1hist), ~replace_na(., FALSE)) %>%
  group_by(histgene, CAF1hist) %>%
  slice_sample(n=200) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, fill =fct_rev(paste(histgene, CAF1hist))), alpha=0.1)+
  geom_jitter(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, col= fct_rev(paste(histgene, CAF1hist))), alpha=0.8)+
  theme_bw()+
  ylab('log2(dTAG/DMSO) Mature')+
  scale_fill_manual(values =c("grey20", "grey40","grey80" ))+
  scale_color_manual(values = c("grey20", "grey40","grey80" ))+
  ggtitle("Sphase-dependent histone expression")+
  ylim(-1,3)+
  theme(legend.position = "none")


mCAF1|mSphase



atac_nascent <- fread('~/Desktop/post-doc/CAF1-DNA-replication/N_clone8_log2_final_genes.csv') %>%
  mutate(histgene = if_else( str_detect(gene_id ,paste(HuHist, collapse = "|")), true = "histone", false = "non-histone")) %>%
  dplyr::filter(gene_type == "protein_coding")

nCAF1 <- left_join(atac_nascent, absodiff, by= "gene_name") %>%
  mutate(CAF1hist = diff <0.8) %>%
  mutate_at(vars(CAF1hist), ~replace_na(., FALSE)) %>%
  group_by(histgene, CAF1hist) %>%
  slice_sample(n=200) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, fill = fct_rev(paste(histgene, CAF1hist))), alpha=0.1)+
  geom_jitter(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, col= fct_rev(paste(histgene, CAF1hist))), alpha=0.8)+
  theme_bw()+
  ylab('log2(dTAG/DMSO) Nascent')+
  scale_fill_manual(values = c("grey20", "grey40","grey80" ))+
  scale_color_manual(values = c("grey20", "grey40","grey80" ))+
  ggtitle("CAF1 dependent histone expression")+
  ylim(-1,3)+
  theme(legend.position = "bottom")


nSphase <- left_join(atac_nascent, absodiff, by= "gene_name") %>%
  mutate(CAF1hist = enrichment > 2.5) %>%
  mutate_at(vars(CAF1hist), ~replace_na(., FALSE)) %>%
  group_by(histgene, CAF1hist) %>%
  slice_sample(n=200) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, fill =fct_rev(paste(histgene, CAF1hist))), alpha=0.1)+
  geom_jitter(aes(x=fct_rev(paste(histgene, CAF1hist)), y= log2_N, col= fct_rev(paste(histgene, CAF1hist))), alpha=0.8)+
  theme_bw()+
  ylab('log2(dTAG/DMSO) Nascent')+
  scale_fill_manual(values =c("grey20", "grey40","grey80" ))+
  scale_color_manual(values = c("grey20", "grey40","grey80" ))+
  ggtitle("Sphase-dependent histone expression")+
  ylim(-1,3)+
  theme(legend.position = "none")


nCAF1|nSphase

nCAF1|mCAF1

nSphase|mSphase




left_join(atac_nascent1, HLB, by= "loc") %>%
  mutate_at(vars(HLB), ~replace_na(., FALSE)) %>%
  group_by(HLB) %>%
  slice_sample(n=200) %>%
  ggplot()+
  geom_boxplot(aes(x=HLB, y= log2_N, fill = HLB), alpha=0.1)+
  geom_jitter(aes(x=HLB, y= log2_N, col = HLB), alpha=0.8)+
  theme_bw()+
  ylab('log2(dTAG/DMSO) Nascent')+
  scale_fill_manual(values =c("grey20", "grey40","grey80" ))+
  scale_color_manual(values = c("grey20", "grey40","grey80" ))+
  ggtitle("Histone Locus Body Nascent")+
  ylim(-1,3)+
  theme(legend.position = "none")
  


absodiff %>% ggplot()+geom_point(aes(x=enrichment, y= diff))


absodiff %>% 
  mutate(category = if_else(diff < 0.8 & enrichment > 2, true = "replicative", false = "nonreplicative")) %>%
  select(name
         ) %>% write_tsv("~/Desktop/histone_category2.txt")
  

Seurat::DimPlot(JD.combined, reduction = "pca", pt.size = 1.2, group.by = 'cond', seed = 42,shuffle = F)+
  ggtitle("RPE-1 dTAG-CAF1 scVASAseq")+theme(legend.position = 'left')+coord_equal()+
  scale_color_manual(values = c("#1F97FF",
                                "#1F97FF",
                                "blue3",
                                "blue3",
                                "#FF7B15", 
                                "#FF7B15" ,
                                "orange4" ,
                                "orange4",
                                "grey80"))

JD.combined
JD_CAF1 <- subset(x= JD.combined, line == "RPE1 CAF1-dTAG")

JD_CAF1@meta.data <- JD_CAF1@meta.data %>% mutate(cond2 = gsub("-pl2", "", cond),
                                                  time = gsub("DMSO-|V1-", "", cond2))

JD_CAF1_2h <- subset(JD_CAF1, time == "2h")



Seurat::Idents(JD_CAF1_2h) <- JD_CAF1_2h@meta.data$cond2


JD_CAF1_2h.markers <- Seurat::FindAllMarkers(JD_CAF1_2h)
JD_CAF1_2h.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25)

JD_CAF1_2h.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  ungroup() -> top100_2hr

JD_CAF1_2h <- Seurat::ScaleData(JD_CAF1_2h)
Seurat::DoHeatmap(JD_CAF1_2h, features = top100_2hr$gene)

top100_2hr %>%
  group_by(cluster) %>%
  summarize(n=n())

JD_CAF1_30h <- subset(JD_CAF1, time == "30h")



Seurat::Idents(JD_CAF1_30h) <- JD_CAF1_30h@meta.data$cond2

JD_CAF1


JD_CAF1_30h.markers <- Seurat::FindAllMarkers(JD_CAF1_30h)
JD_CAF1_30h.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25)

JD_CAF1_30h.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  ungroup() -> top100

JD_CAF1_30h <- Seurat::ScaleData(JD_CAF1_30h)
Seurat::DoHeatmap(JD_CAF1_30h, features = top100$gene)


top100 %>%
  group_by(cluster) %>%
  summarize(n=n())



meta1 <- JD.combined@meta.data %>% rownames_to_column()%>% rename(cell = rowname) 



JD_CAF1@meta.data %>% rownames_to_column()%>% rename(cell = rowname) %>%
  mutate(cond2 = gsub('-pl2', '', cond),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG"
  ) %>%
  group_by(cond2) %>%
  summarize(n=n())




meta1 %>%
  mutate(cond2 = gsub('-pl2', '', cond),
         time = if_else(cond2 %in% c("DMSO-30h", "V1-30h"), true = "30hr", false = "2hr"),
         treatment = if_else(cond2 %in% c("DMSO-30h", "DMSO-2h"), true = "DMSO", false = "dTAG")) %>%
  dplyr::filter(line =="RPE1 CAF1-dTAG"
  ) %>%
  group_by(seurat_clusters, time) %>%
  mutate(total =  n()) %>%
  group_by(cond2, seurat_clusters, time, treatment, total) %>%
  summarise(n= n())%>% distinct(seurat_clusters, time, treatment, n) %>%
  ggplot()+
  geom_col(aes(x=seurat_clusters, y= n, fill= treatment ) )+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15","#1F97FF","#FF7B15" ))+
  facet_grid(rows= vars(time))+theme_bw()+
  xlab('Clusters')+
  ylab('Number of cells')



JD_eG1phase <-subset(x=JD_CAF1, seurat_clusters == c("7", "3"))

Seurat::Idents(JD_eG1phase)

JD_eG1phase.markers <- Seurat::FindAllMarkers(JD_eG1phase)
JD_eG1phase.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.33)

JD_eG1phase.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.33) %>%
  ungroup() -> top100_JD_eG1phase

JD_eG1phase <- Seurat::ScaleData(JD_eG1phase)
Seurat::DoHeatmap(JD_eG1phase, features = top100_JD_eG1phase$gene)

top100_JD_eG1phase %>%
  group_by(cluster) %>%
  summarize(n=n())

