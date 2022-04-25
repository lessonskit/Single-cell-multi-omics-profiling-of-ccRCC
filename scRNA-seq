###scRNA-seq

library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(dplyr)

K81.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC81/RCC81_result/outs/filtered_feature_bc_matrix")
kid81 <- CreateSeuratObject(counts = K81.data, project = "mRCC81", min.cells = 8, min.features = 500)
kid81[["percent.mt"]] <- PercentageFeatureSet(kid81, pattern = "^MT-")
VlnPlot(kid81, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid81 <- subset(kid81, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

K84.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC84/RCC84_result/outs/filtered_feature_bc_matrix")
kid84 <- CreateSeuratObject(counts = K84.data, project = "mRCC84", min.cells = 7, min.features = 500)
kid84[["percent.mt"]] <- PercentageFeatureSet(kid84, pattern = "^MT-")
VlnPlot(kid84, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid84 <- subset(kid84, subset = nFeature_RNA > 500 & nFeature_RNA < 4700 & percent.mt < 10)

K86.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC86/RCC86_result/outs/filtered_feature_bc_matrix")
kid86 <- CreateSeuratObject(counts = K86.data, project = "mRCC86", min.cells = 9, min.features = 500)
kid86[["percent.mt"]] <- PercentageFeatureSet(kid86, pattern = "^MT-")
VlnPlot(kid86, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid86 <- subset(kid86, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K87.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC87/RCC87_result/outs/filtered_feature_bc_matrix")
kid87 <- CreateSeuratObject(counts = K87.data, project = "mRCC87", min.cells = 9, min.features = 500)
kid87[["percent.mt"]] <- PercentageFeatureSet(kid87, pattern = "^MT-")
VlnPlot(kid87, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid87 <- subset(kid87, subset = nFeature_RNA > 500 & nFeature_RNA < 2900 & percent.mt < 10)

K94.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC94/RCC94_result/outs/filtered_feature_bc_matrix")
kid94 <- CreateSeuratObject(counts = K94.data, project = "mRCC94", min.cells = 8, min.features = 500)
kid94[["percent.mt"]] <- PercentageFeatureSet(kid94, pattern = "^MT-")
VlnPlot(kid94, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid94 <- subset(kid94, subset = nFeature_RNA > 500 & nFeature_RNA < 3700 & percent.mt < 10)

K96.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC96/RCC96_result/outs/filtered_feature_bc_matrix")
kid96 <- CreateSeuratObject(counts = K96.data, project = "mRCC96", min.cells = 8, min.features = 500)
kid96[["percent.mt"]] <- PercentageFeatureSet(kid96, pattern = "^MT-")
VlnPlot(kid96, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid96 <- subset(kid96, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K99.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC99/RCC99_result/outs/filtered_feature_bc_matrix")
kid99 <- CreateSeuratObject(counts = K99.data, project = "mRCC99", min.cells = 9, min.features = 500)
kid99[["percent.mt"]] <- PercentageFeatureSet(kid99, pattern = "^MT-")
VlnPlot(kid99, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid99 <- subset(kid99, subset = nFeature_RNA > 500 & nFeature_RNA < 2800 & percent.mt < 10)

K100.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC100/RCC100_result/outs/filtered_feature_bc_matrix")
kid100 <- CreateSeuratObject(counts = K100.data, project = "mRCC100", min.cells = 10, min.features = 500)
kid100[["percent.mt"]] <- PercentageFeatureSet(kid100, pattern = "^MT-")
VlnPlot(kid100, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid100 <- subset(kid100, subset = nFeature_RNA > 500 & nFeature_RNA < 3800 & percent.mt < 10)

K101.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC101/RCC101_result/outs/filtered_feature_bc_matrix")
kid101 <- CreateSeuratObject(counts = K101.data, project = "mRCC101", min.cells = 10, min.features = 500)
kid101[["percent.mt"]] <- PercentageFeatureSet(kid101, pattern = "^MT-")
VlnPlot(kid101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid101 <- subset(kid101, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K103.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC103/RCC103_result/outs/filtered_feature_bc_matrix")
kid103 <- CreateSeuratObject(counts = K103.data, project = "mRCC103", min.cells = 10, min.features = 500)
kid103[["percent.mt"]] <- PercentageFeatureSet(kid103, pattern = "^MT-")
VlnPlot(kid103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid103 <- subset(kid103, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K104.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC104/RCC104_result/outs/filtered_feature_bc_matrix")
kid104 <- CreateSeuratObject(counts = K104.data, project = "mRCC104", min.cells = 10, min.features = 500)
kid104[["percent.mt"]] <- PercentageFeatureSet(kid104, pattern = "^MT-")
VlnPlot(kid104, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid104 <- subset(kid104, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K106.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC106/RCC106_result/outs/filtered_feature_bc_matrix")
kid106 <- CreateSeuratObject(counts = K106.data, project = "mRCC106", min.cells = 10, min.features = 500)
kid106[["percent.mt"]] <- PercentageFeatureSet(kid106, pattern = "^MT-")
VlnPlot(kid106, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid106 <- subset(kid106, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K112.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC112/RCC112_result/outs/filtered_feature_bc_matrix")
kid112 <- CreateSeuratObject(counts = K112.data, project = "mRCC112", min.cells = 10, min.features = 500)
kid112[["percent.mt"]] <- PercentageFeatureSet(kid112, pattern = "^MT-")
VlnPlot(kid112, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid112 <- subset(kid112, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K113.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC113/RCC113_result/outs/filtered_feature_bc_matrix")
kid113 <- CreateSeuratObject(counts = K113.data, project = "mRCC113", min.cells = 9, min.features = 500)
kid113[["percent.mt"]] <- PercentageFeatureSet(kid113, pattern = "^MT-")
VlnPlot(kid113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid113 <- subset(kid113, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K114.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC114/RCC114_result/outs/filtered_feature_bc_matrix")
kid114 <- CreateSeuratObject(counts = K114.data, project = "mRCC114", min.cells = 7, min.features = 500)
kid114[["percent.mt"]] <- PercentageFeatureSet(kid114, pattern = "^MT-")
VlnPlot(kid114, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid114 <- subset(kid114, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 10)

K115.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC115/RCC115_result/outs/filtered_feature_bc_matrix")
kid115 <- CreateSeuratObject(counts = K115.data, project = "mRCC115", min.cells = 9, min.features = 500)
kid115[["percent.mt"]] <- PercentageFeatureSet(kid115, pattern = "^MT-")
VlnPlot(kid115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid115 <- subset(kid115, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K116.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC116/RCC116_result/outs/filtered_feature_bc_matrix")
kid116 <- CreateSeuratObject(counts = K116.data, project = "mRCC116", min.cells = 6, min.features = 500)
kid116[["percent.mt"]] <- PercentageFeatureSet(kid116, pattern = "^MT-")
VlnPlot(kid116, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid116 <- subset(kid116, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K119.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC119/RCC119_result/outs/filtered_feature_bc_matrix")
kid119 <- CreateSeuratObject(counts = K119.data, project = "mRCC119", min.cells = 10, min.features = 500)
kid119[["percent.mt"]] <- PercentageFeatureSet(kid119, pattern = "^MT-")
VlnPlot(kid119, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid119 <- subset(kid119, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K120.data <- Read10X(data.dir = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/mRNA_RCC120/RCC120_result/outs/filtered_feature_bc_matrix")
kid120 <- CreateSeuratObject(counts = K120.data, project = "mRCC120", min.cells = 10, min.features = 500)
kid120[["percent.mt"]] <- PercentageFeatureSet(kid120, pattern = "^MT-")
VlnPlot(kid120, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid120 <- subset(kid120, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

mRCC <- merge(x = kid81, y = list(kid84, kid86, kid87, kid94, kid96, kid99, kid100, kid101, kid103, kid104, kid106, kid112, kid113, kid114, kid115, kid116, kid119, kid120))
mRCC[["percent.mt"]] <- PercentageFeatureSet(mRCC, pattern = "^MT-")
VlnPlot(mRCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size= 0)
plot1 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
mRCC <- NormalizeData(mRCC, normalization.method = "LogNormalize", scale.factor = 10000)
mRCC <- NormalizeData(mRCC)
mRCC <- FindVariableFeatures(mRCC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mRCC), 10)
plot1 <- VariableFeaturePlot(mRCC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
mRCC <- CellCycleScoring(mRCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mRCC <- ScaleData(mRCC, vars.to.regress = c("S.Score", "G2M.Score"))
mRCC <- RunPCA(mRCC, pc.genes = mRCC@var.genes, npcs = 30, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
mRCC <- mRCC %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(mRCC, 'harmony')
harmony_embeddings[1:5, 1:5]
mRCC <- mRCC %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.6) %>% 
    identity()
DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
mRCC.markers <- FindAllMarkers(mRCC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(mRCC.markers,sep="\t",file="/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.xls")
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)
DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
DimPlot(mRCC, reduction = "umap", split.by = "orig.ident",label = TRUE, pt.size = .5, ncol = 5)
##save object
save(mRCC, file= "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.RData")

#Marker genes
features1 <- c("NDUFA4L2","CA9","KRT18","KRT8","KLRD1","KLRB1","GNLY","NKG7","CD3D","CD3E","CD8A","CD8B","IL7R","CD68","CD163","GPNMB","SLC40A1","MSR1","ACTA2","PDGFRB","COL1A2","PECAM1","KDR","CDH5","S100A8","S100A9","LYZ","KCNQ1OT1","CP","CTLA4","FOXP3","CD1C","CD1E","ACKR1","VWF","MKI67","TOP2A","CD79A","CD79B","IGKC","IGLC2","MS4A1","TPSB2","TPSAB1","KIT","KRT19","WFDC2")
DotPlot(mRCC, features = features1, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

###The proportion of each cell group
library(ggplot2)
require(reshape2)
b=c(30339,19819,16444,15115,8720,7742,1728,1979,837)
column=c("T cells","Tumor cells","Monocyte_macrophages","NK cells","Endothelial cells","CAF","Dendritic cells","B cells","Mast cells")
row=c("cell")
A=rbind(b)
dimnames(A)=list(row,column)
lals_pct <- paste(column, "", "(",round(A/sum(A)*100, 2), "%",")", sep="")
pie(A,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 9,name = "Set2"))
title("Proportion of cells in each category",cex.main=2)
#T cells
d=c(14944,12007,1980,1321,87)
column=c("CD8+ T cells","CD4+ T cells","Exhausted T cells","Proliferative T cells 1","Proliferative T cells 2")
row=c("cell")
A=rbind(d)
dimnames(A)=list(row,column)
lals_pct <- paste(column, "", "(",round(A/sum(A)*100, 2), "%",")", sep="")
pie(A,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 9,name = "Set2"))
title("Proportion of T cells in each category",cex.main=2)
#Tumor cells
e=c(15845,2115,1206,653)
column=c("ccRCC 1","ccRCC 2","ccRCC 3","ccRCC 4")
row=c("cell")
A=rbind(e)
dimnames(A)=list(row,column)
lals_pct <- paste(column, "", "(",round(A/sum(A)*100, 2), "%",")", sep="")
pie(A,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 9,name = "Set2"))
title("Proportion of ccRCC cells in each category",cex.main=2)
#EC cells
f=c(7247,1326,121,26)
column=c("Endothelial cells 1","Endothelial cells 2","Endothelial cells 3","Endothelial cells 4")
row=c("cell")
A=rbind(f)
dimnames(A)=list(row,column)
lals_pct <- paste(column, "", "(",round(A/sum(A)*100, 2), "%",")", sep="")
pie(A,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 9,name = "Set2"))
title("Proportion of ccRCC cells in each category",cex.main=2)

###A map of the distribution of various cells
a1=subset(mRCC, idents = 1)
a2=subset(mRCC, idents = 2)
a3=subset(mRCC, idents = 3)
a4=subset(mRCC, idents = 4)
a5=subset(mRCC, idents = 5)
a6=subset(mRCC, idents = 6)
a7=subset(mRCC, idents = 7)
a8=subset(mRCC, idents = 8)
a9=subset(mRCC, idents = 9)
a10=subset(mRCC, idents = 10)
a11=subset(mRCC, idents = 11)
a12=subset(mRCC, idents = 12)
a13=subset(mRCC, idents = 13)
a14=subset(mRCC, idents = 14)
a15=subset(mRCC, idents = 15)
a16=subset(mRCC, idents = 16)
a17=subset(mRCC, idents = 17)
a18=subset(mRCC, idents = 18)
a19=subset(mRCC, idents = 19)
a20=subset(mRCC, idents = 20)
a21=subset(mRCC, idents = 21)
a22=subset(mRCC, idents = 22)

c1=table(a1@meta.data$orig.ident)
c2=table(a2@meta.data$orig.ident)
c3=table(a3@meta.data$orig.ident)
c4=table(a4@meta.data$orig.ident)
c5=table(a5@meta.data$orig.ident)
c6=table(a6@meta.data$orig.ident)
c7=table(a7@meta.data$orig.ident)
c8=table(a8@meta.data$orig.ident)
c9=table(a9@meta.data$orig.ident)
c10=table(a10@meta.data$orig.ident)
c11=table(a11@meta.data$orig.ident)
c12=table(a12@meta.data$orig.ident)
c13=table(a13@meta.data$orig.ident)
c14=table(a14@meta.data$orig.ident)
c15=table(a15@meta.data$orig.ident)
c16=table(a16@meta.data$orig.ident)
c17=table(a17@meta.data$orig.ident)
c18=table(a18@meta.data$orig.ident)
c19=table(a19@meta.data$orig.ident)
c20=table(a20@meta.data$orig.ident)
c21=table(a21@meta.data$orig.ident)
c22=table(a22@meta.data$orig.ident)
C=rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22)

#You take the C, you modify it, and then you plug it in
E=read.csv(file = "/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2/C.csv", header = TRUE)
data<-data.frame(E,celltypes = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
library(reshape2)
mydata <- melt(data,id.vars="celltypes",variable.name="samples",value.name="Proportion", na.rm = TRUE)
mydata$celltypes <- factor(mydata$celltypes, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"), ordered=TRUE)
mydata$samples <- factor(mydata$samples, levels=c("mRCC81", "mRCC84", "mRCC86", "mRCC87", "mRCC94", "mRCC96", "mRCC99","mRCC100","mRCC101","mRCC103","mRCC104","mRCC106","mRCC112","mRCC113","mRCC114","mRCC115","mRCC116","mRCC119","mRCC120"), ordered=TRUE)

library(ggplot2)
library(ggthemes)
p1=ggplot(mydata,aes(celltypes,Proportion,fill=samples))+
geom_bar(stat="identity",position="fill", colour = "white")+
ggtitle("The proportion of samples")+
theme_wsj()+ 
theme(axis.ticks.length=unit(0.5,'cm'))+
guides(fill=guide_legend(title=NULL))

mytheme <- theme_bw()+
theme(legend.position = 'top', panel.border = element_blank(),
panel.grid.major = element_line(linetype = 'dashed'), panel.grid.minor =
element_blank(), legend.text = element_text(size=9,color='#003087',family = "CA"),
plot.title = element_text(size=15,color="#003087",family = "CA"), legend.key =
element_blank(), axis.text = element_text(size=10,color='#003087',family = "CA"),
strip.text = element_text(size=12,color="#EF0808",family = "CA"),
strip.background = element_blank())

p1+mytheme

##### Molecular typing of tumor cells
Tumor=subset(mRCC, idents = c(1,9,15,19))
DimPlot(Tumor, reduction = "umap", label = TRUE, pt.size = .5)

top10gene=c("CA9","NDUFA4L2","CRYAB","TMEM176A","TMEM176B","HILPDA","RARRES2","KRT18","ATP1B1","ANGPTL4","KRT8","KCNQ1OT1","CP","CA12","VMP1","LINC01320","VEGFA","IGFBP3","PAX8","MME","CADM1",
"CYP1B1","TGM2","ELF3","PLEKHA1","CYB5A","CD24","BNIP3","LDHA","ANXA4","PDZK1IP1","FXYD2","TPI1","EGLN3","CMBL","FABP6","NNMT","SPP1","SLPI","C15orf48","IGFBP1","LGALS3","NMB","FHL2","KRT19","LINC00626","DEFB1","TMEM91","WFDC2")
DotPlot(Tumor, features = top10gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
#tumor cells CNV
library(infercnv)
library(tidyverse)
Tumor_cnv=subset(mRCC, idents = c(1,8,9,15,19))
table(Tumor_cnv@active.ident)
Tumor_cnv[["celltype"]]<-Tumor_cnv@active.ident
Tumor_cellAnnota <- subset(Tumor_cnv@meta.data, select='celltype')
Tumor_exprMatrix <- as.matrix(GetAssayData(Tumor_cnv, slot='counts'))

write.table(Tumor_exprMatrix, '/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2/Tumor_exprMatrix.txt', col.names=NA, sep='\t')
write.table(Tumor_cellAnnota, '/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2/Tumor_cellAnnota.txt', col.names=F, sep='\t')

infercnv_Tumor = CreateInfercnvObject(delim = '\t',
                  raw_counts_matrix = Tumor_exprMatrix,
                  annotations_file = Tumor_cellAnnota,
                  gene_order_file = '/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2/geneLocate.txt',
                  ref_group_names = NULL)
dir.create("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2")
#10x data cutoff use 0.1
infercnv_Tumor = infercnv::run(infercnv_Tumor,
                             cutoff=0.1, 
                             out_dir='/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure2', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)


###Trajectory Analysis of tumorcells by monocle3
library(monocle3)
library(SeuratWrappers)
						 
new.cluster.ids <- c("ccRCC 1", "NK cells", "CD8+ T cells", "CD4+ T cells", "TAM 1", "CAF", "Endothelial cells 1", "Monocytes","ccRCC 2","Exhausted T cells","Dendritic cells","TAM 2","Endothelial cells 2","Pro CD8+ T",
"ccRCC 3","Plasma cells","B cells","Mast cells","ccRCC 4","Endothelial cells 3","Exhau_Pro CD8+ T","Endothelial cells 4")
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)

Tumor=subset(mRCC, idents = c("ccRCC 1","ccRCC 2","ccRCC 3","ccRCC 4"))
Tumor.cds <- as.cell_data_set(Tumor)
Tumor.cds <- cluster_cells(cds = Tumor.cds, reduction_method = "UMAP")
Tumor.cds <- learn_graph(Tumor.cds, use_partition = TRUE)

Tumor.cds <- order_cells(Tumor.cds, reduction_method = "UMAP")
plot_cells(
  cds = Tumor.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

DimPlot(Tumor, reduction = "umap", split.by = "orig.ident", pt.size = .1, ncol = 5)


###Analysis between outcomes by TCGA database

ccRCC1=c(44,6)
ccRCC2=c(28,9)
ccRCC3=c(38,5)
ccRCC4=c(19,21)

data<-data.frame(Celltypes = c("ccRCC1","ccRCC2","ccRCC3","ccRCC4"), positive = c(44,28,38,19),negctive = c(6,9,5,21))
library(reshape2)
mydata <- melt(data,id.vars="Celltypes",variable.name="Survival",value.name="Gene_Num")

library(ggplot2)
library(ggthemes)
ggplot(mydata,aes(Celltypes,Gene_Num,fill=Survival))+
geom_bar(stat="identity",position="dodge",colour = "white")+
ggtitle("The relationship between Tumor cell types and prognosis")+
theme_bw()

#chi-square
ccRCC1=c(44,6)
ccRCC2=c(28,9)
ccRCC3=c(38,5)
ccRCC4=c(19,21)

x <- matrix(c(44, 6, 19, 21), ncol = 2)
chisq.test(x)


######Endothelial cell subtype
End=subset(mRCC, idents = c("Endothelial cells 1","Endothelial cells 2","Endothelial cells 3","Endothelial cells 4"))
top10gene=c("PECAM1","KDR","VWF","CDH5","FLT1","PLVAP","ESM1","TIMP3","PLPP1","ACKR1","RNASE1","RAMP3","HYAL2","ENG","NOTCH3","PLAC9","TBX2","COL1A2","PDGFRB","CTSG","TPSB2","TPSAB1","HDC","SSUH2","ENPP7","TMEM174")
DotPlot(End, features = top10gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

attach(mRCC.markers)
kid5cgene=mRCC.markers[which(gene%in% c("PECAM1","KDR","VWF","CDH5","FLT1","PLVAP","ESM1","TIMP3","PLPP1","ACKR1","RNASE1","RAMP3","HYAL2","ENG","NOTCH3","PLAC9","ACTA2","COL1A2","PDGFRB","TAGLN","TPSB2","TPSAB1","HDC","SSUH2","ENPP7","TMEM174")&cluster%in% c(7,13,20,22)),]
kid5cgene$gene <- factor(kid5cgene$gene, levels=c("PECAM1","KDR","VWF","CDH5","FLT1","PLVAP","ESM1","TIMP3","PLPP1","ACKR1","RNASE1","RAMP3","HYAL2","ENG","NOTCH3","PLAC9","ACTA2","COL1A2","PDGFRB","TAGLN","TPSB2","TPSAB1","HDC","SSUH2","ENPP7","TMEM174"), ordered=TRUE)
detach(mRCC.markers)

library(ggplot2)
p1 <- ggplot(kid5cgene, aes(x=cluster,y=gene)) +  xlab("cluster") + theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank())  + theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(angle=45,hjust=1, vjust=1,size=18,color="black"),axis.text.y=element_text(size=15,color="black")) + theme(legend.position="top") +  geom_tile(aes(fill=avg_log2FC)) + scale_fill_gradient2(low = "navy", mid = "white",high = "red")
#Relationship with ligand receptor of tumor cells
EC=subset(mRCC, idents = c("ccRCC 1", "Endothelial cells 1","ccRCC 2","Endothelial cells 2", "ccRCC 3","ccRCC 4", "Endothelial cells 3","Endothelial cells 4"))

setwd("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure3/cellphoneDB")
write.table(as.matrix(EC@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
EC@meta.data$cell_type=EC@active.ident
meta_data <- cbind(rownames(EC@meta.data), EC@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#After generating the file, transfer it to the Linux system to obtain the OUT file
pbmc='/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure3/cellphoneDB/out/' 
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(scales)
library(ggsci)

#Cell interactions
netf<- "count_network.txt"
mynet <- read.delim(paste0(pbmc,"count_network.txt"), check.names = FALSE)
table(mynet$count)
mynet %>% filter(count>0) -> mynet 
head(mynet)
net<- graph_from_data_frame(mynet) 
plot(net)

allcolour=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
"#E64B3533", "#4DBBD533", "#00A08733", "#3C548833", "#F39B7F33", "#8491B433", "#91D1C233", "#DC000033", "#7E614833", "#B09C8533",
"#E64B3566", "#4DBBD566", "#00A08766", "#3C548866", "#F39B7F66", "#8491B466", "#91D1C266", "#DC000066", "#7E614866", "#B09C8566")

							 
E(net)$width  <- E(net)$count/10

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     vertex.label.cex=.7) 

net2 <- net

for (i in 1: length(unique(mynet$SOURCE)) ){ 
        E(net)[map(unique(mynet$SOURCE),function(x) {
                get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
        })%>% unlist()]$color <- allcolour[i]
}

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     vertex.label.cex=.7) 
	 
#Ligand receptor pairs
mypvals <- read.delim(paste0(pbmc,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(pbmc,"means.txt"), check.names = FALSE)

# 

END <- grep("NRP1|ITGB1|KDR|FLT1|CD63|INSR|NRP2|CALCRL|SPP1", 
            mymeans$interacting_pair,value = T)

mymeans %>% dplyr::filter(interacting_pair %in% END)%>%
  dplyr::select("interacting_pair",starts_with("Endothelial cells 1"),ends_with("Endothelial cells 1"))  %>%  
  reshape2::melt() -> meansdf

colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals %>% dplyr::filter(interacting_pair %in% END)%>%
  dplyr::select("interacting_pair",starts_with("Endothelial cells 1"),ends_with("Endothelial cells 1"))%>%  
  reshape2::melt()-> pvalsdf

colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >1))$means)

pldf%>% filter(means >1) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="#DC0000FF",mid = "#4DBBD5FF",low ="#3C5488FF",midpoint = 1)+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 

 




#Tumor cells and Immune cells
.libPaths(c("/usr/local/lib/R/site-library","/usr/local/lib/R/library", .libPaths()))
library(Seurat)
library(reshape2)
library(DT)
library(ggthemes)
library(ggplot2)
require(tidyr)
require(ComplexHeatmap)
require(circlize)

setwd("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure3")
load("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.RData")

immune=subset(mRCC,idents = c('2','3','4','5','8','10','11','12','14','16','17','18','21'))
nonimmune=subset(mRCC,idents = c('1','6','7','9','13','15','19','20','22'))
DimPlot(immune, reduction = "umap", label = TRUE, pt.size = .5)

B1=c("immune","nonimmune")
B2=c(66442,36281)
B=data.frame(B1, B2)
lals_pct <- paste(B1, "", "(",round(B2/sum(B2)*100, 2), "%",")", sep="")
library(RColorBrewer)
pie(B2,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 2,name = "Set2"))

A=table(immune@meta.data$RNA_snn_res.0.6)
barplot(A, col=RColorBrewer::brewer.pal(n = 12,name = "Set2"))

###T cells
TC=subset(mRCC,idents = c('3','4','10','14','21'))
DimPlot(TC, reduction = "umap", label = TRUE, pt.size = .5)

VlnPlot(object = TC, pt.size= 0,features = c("CD3D","CD8A", "FOXP3","PDCD1","MKI67","CD3E","IL7R","IL2RA","HAVCR2","TOP2A"),ncol=5)

##The proportion of T cells expressing PD1 and those expressing PD1 greater than 1 were screened out
BB=subset(TC, PDCD1>1)
FeaturePlot(BB, features = "PDCD1",cols = c("gray", "red"))
B3=c("PD1+ T cell","PD1- T cell")
B4=c(6520,23819)
B=data.frame(B3, B4)
lals_pct <- paste(B3, "", "(",round(B4/sum(B4)*100, 2), "%",")", sep="")
library(RColorBrewer)
pie(B4,labels = lals_pct,col=RColorBrewer::brewer.pal(n = 2,name = "Set2"))

####TAM
TAM=subset(mRCC,idents = c('5','12'))
DimPlot(TAM, reduction = "umap", label = FALSE, pt.size = .5)
VlnPlot(object = TAM, pt.size= 0,features = c("CD68","GPNMB","SLC40A1","MSR1","CD163", "CSF1R","MRC1","CD86"),ncol=4)

##Cellphone DB
new.cluster.ids <- c("ccRCC 1", "NK cells", "Exhau CD8+ T", "CD4+ T cells", "TAM 1", "CAF", "Endothelial cells 1", "Monocytes","ccRCC 2","Treg cells","Dendritic cells","TAM 2","Endothelial cells 2","Pro CD8+ T",
"ccRCC 3","Plasma cells","B cells","Mast cells","ccRCC 4","Endothelial cells 3","Exhau_Pro CD8+ T","Endothelial cells 4")
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)

immune=subset(mRCC,idents = c("NK cells", "Exhau CD8+ T", "CD4+ T cells", "TAM 1", "Monocytes","Treg cells","Dendritic cells","TAM 2","Pro CD8+ T",
"Plasma cells","B cells","Mast cells","Exhau_Pro CD8+ T"))

setwd("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure4/cellphoneDB")
write.table(as.matrix(immune@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
immune@meta.data$cell_type=immune@active.ident
meta_data <- cbind(rownames(immune@meta.data), immune@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#After generating the file, transfer it to the Linux system to obtain the OUT file
pbmc='/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA_result/figure4/cellphoneDB/out/' 
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(scales)
library(ggsci)

#Cell interaction
netf<- "count_network.txt"
mynet <- read.delim(paste0(pbmc,"count_network.txt"), check.names = FALSE) 
table(mynet$count)
mynet %>% filter(count>0) -> mynet  
head(mynet)
net<- graph_from_data_frame(mynet) 
plot(net)

allcolour=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
"#E64B3533", "#4DBBD533", "#00A08733", "#3C548833", "#F39B7F33", "#8491B433", "#91D1C233", "#DC000033", "#7E614833", "#B09C8533",
"#E64B3566", "#4DBBD566", "#00A08766", "#3C548866", "#F39B7F66", "#8491B466", "#91D1C266", "#DC000066", "#7E614866", "#B09C8566")

							 
E(net)$width  <- E(net)$count/10

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     vertex.label.cex=.7) 

net2 <- net

for (i in 1: length(unique(mynet$SOURCE)) ){ 
        E(net)[map(unique(mynet$SOURCE),function(x) {
                get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
        })%>% unlist()]$color <- allcolour[i]
}

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     vertex.label.cex=.7) 
	 
#Cell interactions
mypvals <- read.delim(paste0(pbmc,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(pbmc,"means.txt"), check.names = FALSE)

chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", 
            mymeans$interacting_pair,value = T)
th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", 
            mymeans$interacting_pair,value = T)
th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", 
             mymeans$interacting_pair,value = T)
treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", mymeans$interacting_pair,value = T)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
                      mymeans$interacting_pair,value = T)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", 
                     mymeans$interacting_pair,value = T)
niche <- grep("CSF", mymeans$interacting_pair,value = T)

mymeans %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("Exhau_Pro CD8+ T"),ends_with("Exhau_Pro CD8+ T"))  %>%  
  reshape2::melt() -> meansdf

colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("Exhau_Pro CD8+ T"),ends_with("Exhau_Pro CD8+ T"))%>%  
  reshape2::melt()-> pvalsdf

colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >1))$means)

pldf%>% filter(means >1) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="#DC0000FF",mid = "#4DBBD5FF",low ="#3C5488FF",midpoint = 1)+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 


####Sample mutation genes
library(ComplexHeatmap)
library(ggplot2)
require(reshape2)
require(scales)
mu=read.csv("mutation.csv",header = T)
attach(mu)
mu$Sample <- factor(mu$Sample, levels=c("RCC81", "RCC84", "RCC86", "RCC87", "RCC94", "RCC96", "RCC99", "RCC100","RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113", "RCC114", "RCC115","RCC116", "RCC119", "RCC120"), ordered=TRUE) #调X轴指标顺序
mu$Gene <- factor(mu$Gene, levels=c("VHL","PBRM1","SETD2","NKX2-4","SNTG1","ZUFSP","KCNU1","MTOR","TP53","CD24","KRT14","CA9"), ordered=TRUE)
p1 <- ggplot(mu, aes(x=Sample,y=Gene)) +  xlab("Sample") + theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank())  + theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
axis.text.x=element_text(angle=45,hjust=1, vjust=1,size=18,color="black"),axis.text.y=element_text(size=15,color="black")) + theme(legend.position="top") +  geom_tile(aes(fill=Mutation)) + scale_fill_manual(breaks = c("Nonsense","Missense","Nonframe shift deletion","Frame shift deletion","Nonframe shift insertion","Frame shift insertion","NA"),values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF","white"))

M=c(15.8,15.8,47.3,5.3,10.5,5.3,5.3,5.3,5.3,21.1,15.8,84.2)
barplot(M, col= "#B09C85FF")


####Mass spectrometry and transcriptome correlation analysis
library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(dplyr)

load("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.RData")

new.cluster.ids <- c("ccRCC 1", "NK cells", "CD8+ T cells", "CD4+ T cells", "TAM 1", "CAF", "Endothelial cells 1", "Monocytes","ccRCC 2","Exhausted T cells","Dendritic cells","TAM 2","Endothelial cells 2","Pro CD8+ T",
"ccRCC 3","Plasma cells","B cells","Mast cells","ccRCC 4","Endothelial cells 3","Exhau_Pro CD8+ T","Endothelial cells 4")
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)

features1 <- c("TPM1","TPM4","MYH9","HSPB1","CNBP","MYL6","PUF60","KRT8")
DotPlot(mRCC, features = features1, cols = c("blue4", "yellow"), dot.scale = 8) + RotatedAxis()

library(ggplot2)
library(Cairo)

setwd("/home/yuzhenyuan/single_cell_data/fig7")
CTB <-read.csv("CTB.csv",header = T)

CTB$Gene <- factor(CTB$Gene, levels=c("TPM1","TPM4","MYH9","HSPB1","CNBP","MYL6","PUF60","KRT8"), ordered=TRUE)


p1=ggplot(CTB,aes(x=Gene,y=lncRNA)) + 
  geom_point(aes(size=Coverage,color=1*(Score)))+
  scale_colour_gradient(low="yellow4",high="yellow")+
  labs(
       color=expression(Score),
       size="Coverage",
       x="Gene"
       # y="Coverage")
      )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )



###Tumor cells were extracted and grouped according to different sample sources to explore gene expression characteristics
.libPaths(c("/data_8t/file/teacher/yuzhenyuan/R_packages/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(dplyr)
load("/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.RData")

new.cluster.ids <- c("ccRCC 1", "NK cells", "CD8+ T cells", "CD4+ T cells", "TAM 1", "CAF", "Endothelial cells 1", "Monocytes","ccRCC 2","Exhausted T cells","Dendritic cells","TAM 2","Endothelial cells 2","Proliferative T cells 1",
"ccRCC 3","Plasma cells","B cells","Mast cells","ccRCC 4","Endothelial cells 3","Proliferative T cells 2","Endothelial cells 4")
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)
Tumor=subset(mRCC, idents = c("ccRCC 1","ccRCC 2","ccRCC 3","ccRCC 4"))

Tumor@active.ident <- as.factor(Tumor$orig.ident)
VHL <- subset(Tumor, ident=c("mRCC81","mRCC84","mRCC86","mRCC94","mRCC96","mRCC99","mRCC101","mRCC103","mRCC104","mRCC106","mRCC112","mRCC113","mRCC114","mRCC115","mRCC119","mRCC120"))
NonVHL <- subset(Tumor, ident=c("mRCC87","mRCC100","mRCC116"))

DimPlot(VHL, reduction = "umap", pt.size = .5)
DimPlot(NonVHL, reduction = "umap", pt.size = .5)
avg1=AverageExpression(VHL)
avg1 <- log1p(avg1$RNA)
avg1=as.data.frame(avg1)
avg1$Mean1= rowMeans(avg1)
avg1$gene <- rownames(avg1)

avg2=AverageExpression(NonVHL)
avg2 <- log1p(avg2$RNA)
avg2=as.data.frame(avg2)
avg2$Mean2= rowMeans(avg2)
avg2$gene <- rownames(avg2)

A=merge(avg1,avg2,by = 'gene')

r1 <- cor(A$Mean1,A$Mean2)
r1
attach(A)
library(ggplot2)
p1 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("R=0.984") + stat_smooth(method=lm)
p2 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("0.984") + stat_smooth(method=lm) + geom_text(aes(label=gene), size=3)
detach(A)

CA9 <- subset(Tumor, ident=c("mRCC86","mRCC100","mRCC119"))
NonCA9 <- subset(Tumor, ident=c("mRCC81","mRCC84","mRCC87","mRCC94","mRCC96","mRCC99","mRCC101","mRCC103","mRCC104","mRCC106","mRCC112","mRCC113","mRCC114","mRCC115","mRCC116","mRCC120"))

DimPlot(CA9, reduction = "umap", pt.size = .5)
DimPlot(NonCA9, reduction = "umap", pt.size = .5)
avg1=AverageExpression(CA9)
avg1 <- log1p(avg1$RNA)
avg1=as.data.frame(avg1)
avg1$Mean1= rowMeans(avg1)
avg1$gene <- rownames(avg1)

avg2=AverageExpression(NonCA9)
avg2 <- log1p(avg2$RNA)
avg2=as.data.frame(avg2)
avg2$Mean2= rowMeans(avg2)
avg2$gene <- rownames(avg2)

A=merge(avg1,avg2,by = 'gene')

r1 <- cor(A$Mean1,A$Mean2)
r1
attach(A)
library(ggplot2)
p1 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("R=0.987") + stat_smooth(method=lm)
p2 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("0.987") + stat_smooth(method=lm) + geom_text(aes(label=gene), size=3)
detach(A)


KRT14 <- subset(Tumor, ident=c("mRCC81","mRCC86","mRCC116"))
NonKRT14 <- subset(Tumor, ident=c("mRCC84","mRCC87","mRCC94","mRCC96","mRCC99","mRCC100","mRCC101","mRCC103","mRCC104","mRCC106","mRCC112","mRCC113","mRCC114","mRCC115","mRCC119","mRCC120"))

DimPlot(KRT14, reduction = "umap", pt.size = .5)
DimPlot(NonKRT14, reduction = "umap", pt.size = .5)
avg1=AverageExpression(KRT14)
avg1 <- log1p(avg1$RNA)
avg1=as.data.frame(avg1)
avg1$Mean1= rowMeans(avg1)
avg1$gene <- rownames(avg1)

avg2=AverageExpression(NonKRT14)
avg2 <- log1p(avg2$RNA)
avg2=as.data.frame(avg2)
avg2$Mean2= rowMeans(avg2)
avg2$gene <- rownames(avg2)

A=merge(avg1,avg2,by = 'gene')

r1 <- cor(A$Mean1,A$Mean2)
r1
attach(A)
library(ggplot2)
p1 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("R=0.988") + stat_smooth(method=lm)
p2 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("0.988") + stat_smooth(method=lm) + geom_text(aes(label=gene), size=3)
detach(A)

CD24 <- subset(Tumor, ident=c("mRCC84","mRCC86","mRCC96","mRCC99","mRCC100","mRCC113","mRCC115","mRCC116","mRCC120"))
NonCD24 <- subset(Tumor, ident=c("mRCC81","mRCC87","mRCC94","mRCC101","mRCC103","mRCC104","mRCC106","mRCC112","mRCC114","mRCC119"))

DimPlot(CD24, reduction = "umap", pt.size = .5)
DimPlot(NonCD24, reduction = "umap", pt.size = .5)
avg1=AverageExpression(CD24)
avg1 <- log1p(avg1$RNA)
avg1=as.data.frame(avg1)
avg1$Mean1= rowMeans(avg1)
avg1$gene <- rownames(avg1)

avg2=AverageExpression(NonCD24)
avg2 <- log1p(avg2$RNA)
avg2=as.data.frame(avg2)
avg2$Mean2= rowMeans(avg2)
avg2$gene <- rownames(avg2)

A=merge(avg1,avg2,by = 'gene')

r1 <- cor(A$Mean1,A$Mean2)
r1
attach(A)
library(ggplot2)
p1 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("R=0.988") + stat_smooth(method=lm)
p2 <- ggplot(A, aes(x=Mean1,y=Mean2)) + geom_point() + ggtitle("0.988") + stat_smooth(method=lm) + geom_text(aes(label=gene), size=3)
detach(A)

###end
