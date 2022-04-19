#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(infercnv)
library(AnnotationHub)
library(ensembldb)


# :: Glimmunology :: ----

# A :: snRNA Sample_Y GBM ----


rm(sid, object_1)
gc()


sid <- 'van_Hijfte_Sample_Y'
object_1 <- Read10X(data.dir = "data/Glimmunology_GBM_1/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix")
object_1 <- CreateSeuratObject(counts = object_1,
                               min.cells = 3,
                               min.features = 200,
                               project = "glioma_glim")
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 

    
ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1400,col="red") +
  geom_hline(yintercept = 4500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 2200,col="red") +
  geom_hline(yintercept = 14000,col="red")
  #scale_y_log10()

object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 1400 &
                     nFeature_RNA < 4500 & 
                     nCount_RNA > 2200 &
                     nCount_RNA < 14000 &
                     percent.mito < 0.025)

object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2

## scaling of data ----
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

## cluster the cells ----

object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
head(Idents(object_1), 20)


## UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")

# levels(object_1$seurat_clusters) <- gsub("^21$","\\1.Neuron",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(5|10|13|14)$","Immune + T-cells.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(4|6|12|19)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^16$","Endothelial",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^15$","Pericytes",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|9|8|11)$","Tumor.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^7$","Tumor/Dividing",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^17$","Tumor outlier",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^20$","Astrocyte",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^18$","Tumor/Apoptosis?",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^22$","Immune cell / OD hybrid?",levels(object_1$seurat_clusters))

object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(21), "NE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(5,10,13,14), "TAM/MG", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(4,6,12,19), "OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16), "EN", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15), "PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0,1,2,3,9,8,11), "T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(7), "T", object_1$class) # dividing
object_1$class <- ifelse(object_1$seurat_clusters %in% c(17), "T ?", object_1$class) #  Outlier?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18), "T ?", object_1$class) # Apoptotic?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(22), "TAM/MG|OD", object_1$class)
object_1$class <- ifelse(object_1@reductions$umap@cell.embeddings[,1] >= 10 &
                           object_1@reductions$umap@cell.embeddings[,1] <= 11 &
                           object_1@reductions$umap@cell.embeddings[,2] >= 1.5 &
                           object_1@reductions$umap@cell.embeddings[,2] <= 3,
                         "TC", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","7. T","8. T","9. T","11. T",
  "17. T ?","18. T ?",
  "20. AC",
  "21. NE",
  "4. OD","6. OD","12. OD","19. OD",
  "16. EN",
  "15. PE",
  "22. TAM/MG|OD" ,
  "5. TAM/MG","10. TAM/MG","13. TAM/MG","14. TAM/MG",
  "14. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")  +
  labs(subtitle=sid) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3)))



#ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.pdf"),width=10,height=8)
#ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.png"),width=12,height=10)


od.markers <- FindMarkers(object_1, ident.1 = c(4,6,12,19,22))
View(od.markers)



# 
# 
# tmp.17 <- FindMarkers(object_1, ident.1 = 17)
# head(tmp.17,20)
# 
# tmp.22 <- FindMarkers(object_1, ident.1 = 22)
# head(tmp.22,20)

tmp.15 <- FindMarkers(object_1, ident.1 = 15) # PE
View(tmp.15)




## 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = "PDGFRB") # Tumor/MES


FeaturePlot(object = object_1, features = "EREG") # EREG
FeaturePlot(object = object_1, features = "EGF") # BTC
FeaturePlot(object = object_1, features = "BTC") # BTC

DotPlot(object = object_1, features =list( 'ligands'=c('EGF','EREG','AREG','BTC','EPGN','HBEGF') ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("EGFR and ligands in: ",sid))


FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


# AC
FeaturePlot(object = object_1, features = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1"))
FeaturePlot(object = object_1, features = c("GFAP","FABP7","BCAN","PON2","METTL7B","SPARC","GATM","RAMP1")) # Tumor
FeaturePlot(object = object_1, features = c("PMP2","AQP4","DBI","EDNRB","PTPRZ1","CLU","PMP22","ATP1A2")) # Tumor
FeaturePlot(object = object_1, features = c("S100A16","HEY1","PCDHGC3","TTYH1","NDRG2","PRCP","ATP1B2","AGT","PLTP","GPM6B"))
FeaturePlot(object = object_1, features = c("F3","RAB31","PPAP2B","ANXA5","TSPAN7")) # Tumor

# MES1
FeaturePlot(object = object_1, features = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2")) # Tumor

# MES2
FeaturePlot(object = object_1, features = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT")) # Tumor

# OPC
FeaturePlot(object = object_1, features = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4")) # Tumor

# NPC1
FeaturePlot(object = object_1, features = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST")) # Tumor

# NPC2
FeaturePlot(object = object_1, features = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4")) # Tumor

# NPC1, NPC2, Neuron
FeaturePlot(object = object_1, features = c("BCAN", "NREP", "RBFOX3")) # Tumor
# BCAN


# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")

FeaturePlot(object = object_1, features = "NODAL")


## 2. Astrocyte (+) ----

FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



## 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



## 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


## 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor
FeaturePlot(object = object_1, features = "HBG2") # Tumor


## 3D. ? Mono/Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))



## 4. Neurons (+) ----


tmp <- list('C1'=neuron.genes[neuron.genes %in% NPC2 == F] ,
           'NPC1'=NPC1[NPC1 %in% NPC2 == F] ,
           'NPC1+2' = intersect(NPC1, NPC2),
           'NPC2'=NPC2[NPC2 %in% c(NPC1, neuron.genes) == F],
           'NPC2 + C1' = intersect(neuron.genes, NPC2))

DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C1/NPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.pdf"),width=7.5*3, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.png"),width=7.5*3, height=4,scale=1.2)




DotPlot(object = object_1, features = c("EGFR", "GFAP","MOG", "PLP1", "TMEM144", 
                                        "RBFOX1", "RBFOX2", "RBFOX3", "CD2",
                                        "CD3D", "P2RY12", "CD163", "ABCB1", "RGS5"
                                        ))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


FeaturePlot(object = object_1, features = c("SOX4", "RBFOX3"))


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")


#FeaturePlot(object = object_1, features = "DCN") # DCN
#FeaturePlot(object = object_1, features = "COL1A2") # DCN
FeaturePlot(object = object_1, features = "ANPEP") # DCN


## 5. Oligodendrocytes (+) ----


# + "PEAR1", "HEYL" , "CFH"
DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C2.pdf"),width=7.5*2, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C2.png"),width=7.5*2, height=4,scale=1.2)



FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


FeaturePlot(object = object_1, features = "ST18") # OD?
FeaturePlot(object = object_1, features = "MBP") # OD?
FeaturePlot(object = object_1, features = "CTNNA3") # OD?
FeaturePlot(object = object_1, features = "SLC24A2") # OD?
FeaturePlot(object = object_1, features = "KIRREL3") # OD?
FeaturePlot(object = object_1, features = "NKAIN2") # OD?
FeaturePlot(object = object_1, features = "MAP7") # OD?
FeaturePlot(object = object_1, features = "RNF220") # OD?
FeaturePlot(object = object_1, features = "PEX5L") # OD?
FeaturePlot(object = object_1, features = "TMEM144") # OD?
FeaturePlot(object = object_1, features = "EDIL3") # OD?
FeaturePlot(object = object_1, features = "DOCK5") # OD?
FeaturePlot(object = object_1, features = "MOBP") # OD?
FeaturePlot(object = object_1, features = "UNC5C") # OD?
FeaturePlot(object = object_1, features = "CLDN11") # OD?
FeaturePlot(object = object_1, features = "SPOCK3") # OD?
FeaturePlot(object = object_1, features = "CNTNAP4") # OD?
FeaturePlot(object = object_1, features = "MAN2A1") # OD?
FeaturePlot(object = object_1, features = "PCSK6") # OD?
FeaturePlot(object = object_1, features = "TTLL7") # OD?

FeaturePlot(object = object_1, features = "OLIG2") # OD?


## 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



## 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))




## 7. Cycling cells (?) ----


FeaturePlot(object = object_1, features = "TOP2A" )
FeaturePlot(object = object_1, features = "AURKA" )
FeaturePlot(object = object_1, features = "AURKB" )
FeaturePlot(object = object_1, features = "BUB1" )
FeaturePlot(object = object_1, features = "BUB1B" )
FeaturePlot(object = object_1, features = "CDC20" )
FeaturePlot(object = object_1, features = "CENPF" )
FeaturePlot(object = object_1, features = "FAM64A" )
FeaturePlot(object = object_1, features = "FOXM1" )
FeaturePlot(object = object_1, features = "TACC3" )
FeaturePlot(object = object_1, features = "TMPO" )
FeaturePlot(object = object_1, features = "TPX2" )
FeaturePlot(object = object_1, features = "TUBA1C" )


# FeaturePlot(object = object_1, features = "RRM2" )
FeaturePlot(object = object_1, features = "PCNA" )
FeaturePlot(object = object_1, features = "KIAA0101" )
# FeaturePlot(object = object_1, features = "HIST1H4C" )
# FeaturePlot(object = object_1, features = "MLF1IP" )
# FeaturePlot(object = object_1, features = "GMNN" )
FeaturePlot(object = object_1, features = "RNASEH2A" )
# FeaturePlot(object = object_1, features = "MELK" )
# FeaturePlot(object = object_1, features = "CENPK" )
# FeaturePlot(object = object_1, features = "TK1" )
FeaturePlot(object = object_1, features = "TMEM106C" )
# FeaturePlot(object = object_1, features = "CDCA5" )
FeaturePlot(object = object_1, features = "CKS1B" )
FeaturePlot(object = object_1, features = "CDC45" )
FeaturePlot(object = object_1, features = "MCM3" )
FeaturePlot(object = object_1, features = "CENPM" )
FeaturePlot(object = object_1, features = "AURKB" )
FeaturePlot(object = object_1, features = "PKMYT1" )
FeaturePlot(object = object_1, features = "MCM4" )
FeaturePlot(object = object_1, features = "ASF1B" )
FeaturePlot(object = object_1, features = "GINS2" )
FeaturePlot(object = object_1, features = "MCM2" )
FeaturePlot(object = object_1, features = "FEN1" )
FeaturePlot(object = object_1, features = "RRM1" )
FeaturePlot(object = object_1, features = "DUT" )
FeaturePlot(object = object_1, features = "RAD51AP1" )
FeaturePlot(object = object_1, features = "MCM7" )
FeaturePlot(object = object_1, features = "CCNE2" )
FeaturePlot(object = object_1, features = "ZWINT" )


FeaturePlot(object = object_1, features = "CCNB1" )
FeaturePlot(object = object_1, features = "CDC20" )
FeaturePlot(object = object_1, features = "CCNB2" )
FeaturePlot(object = object_1, features = "PLK1" )
FeaturePlot(object = object_1, features = "CCNA2" )
FeaturePlot(object = object_1, features = "CKAP2" )
FeaturePlot(object = object_1, features = "KNSTRN" )
FeaturePlot(object = object_1, features = "RACGAP1" )
FeaturePlot(object = object_1, features = "CDCA3" )
FeaturePlot(object = object_1, features = "TROAP" )
FeaturePlot(object = object_1, features = "KIF2C" )
FeaturePlot(object = object_1, features = "KPNA2" )
FeaturePlot(object = object_1, features = "KIF20A" )
FeaturePlot(object = object_1, features = "ECT2" )
FeaturePlot(object = object_1, features = "CDCA8" )
FeaturePlot(object = object_1, features = "TTK" )
FeaturePlot(object = object_1, features = "NCAPD2" )
FeaturePlot(object = object_1, features = "ARL6IP1" )
FeaturePlot(object = object_1, features = "KIF4A" )
FeaturePlot(object = object_1, features = "CKAP2L" )
FeaturePlot(object = object_1, features = "MZT1" )
FeaturePlot(object = object_1, features = "KIFC1" )
FeaturePlot(object = object_1, features = "SPAG5" )
FeaturePlot(object = object_1, features = "ANP32E" )
FeaturePlot(object = object_1, features = "KIF11" )
FeaturePlot(object = object_1, features = "PSRC1" )
FeaturePlot(object = object_1, features = "TUBB4B" )
FeaturePlot(object = object_1, features = "SMC4" )
FeaturePlot(object = object_1, features = "MXD3" )
FeaturePlot(object = object_1, features = "CDC25B" )
FeaturePlot(object = object_1, features = "OIP5" )
FeaturePlot(object = object_1, features = "REEP4" )
FeaturePlot(object = object_1, features = "GPSM2" )
FeaturePlot(object = object_1, features = "HMGB3" )
FeaturePlot(object = object_1, features = "ARHGAP11A" )
FeaturePlot(object = object_1, features = "RANGAP1" )
FeaturePlot(object = object_1, features = "H2AFZ" )

## 8. longitudinal sig. chr6 ----


sig <- c("H4C1", "H3C2", "H2AC4", "H3C3",            "HIST1H4A","HIST1H3B", "HIST1H2AB", "HIST1H3C",
         "H1-6", "H3C7", "H2BC9", "H2BC11",          "HIST1H1T","HIST1H3F","HIST1H2BH","HIST1H2BH",
         "H2AC11", "H2BC12", "H2AC12", "H2BC13",     "HIST1H2AG","HIST1H2BK","HIST1H2AH","HIST1H2BL",
         "H2AC13", "H3C10", "H2AC14", "H2BC14",      "HIST1H2AI","HIST1H3H","HIST1H2AJ","H2BC14",
         "H2AC15", "H2AC16", "H1-5", "H3C11",        "HIST1H2AK","HIST1H2AL","HIST1H1B","HIST1H3I",
         "H3C12", "H2BC17",                          "HIST1H3J","HIST1H2BO"
)
sig <- unique(sig)

sig %in% all.genes
sig <- sig[ sig %in% all.genes]
sig


DotPlot(object_1, features=sig, group.by = "seurat_clusters")



FeaturePlot(object_1, features=sig)


# up.3 lgg ----

set <- c("OPRD1","RGS8","RAB7B","SULT1C2","TNFAIP6","DLX1","DLX2","SP9","LINC01792","KCNE4","MLPH","COL8A1","SPATA18","BANK1","TENM3-AS1","PDLIM4","SPRY4-AS1","LINC01411","PTCHD4","SP8","AQP1","MET","SLC13A4","SCARA5","NRG1","SFRP1","FABP5","CCN3","DMRT2","ASS1","TDRD1","VAX1","P2RY2","AP002761.4","PTGDR","ASB2","GREM1","AC087473.1","ISLR","ALDH1A3","HS3ST3A1","SLC47A1","CCL18","AC015909.3","TYMSOS","ADCYAP1","CDH19","PODNL1","ZNF98","GPR143")
set <- c("LGR6","NTSR2","AC012593.1","CD8B2","CHL1-AS2","SIDT1","GRAMD1C","ALDH1L1-AS2","GABRG1","AC107398.3","ETNPPL","FSTL5","PLCXD3","INSYN2B","ENPP5","MPP6","ITPRID1","CALN1","AC006148.1","AC019257.1","AC027117.1","IQANK1","TRPM3","BRINP1","SFTPD-AS1","AL157396.1","ZNF214","NLRP14","SPX","KSR2","AL137139.2","AL157769.1","RGS6","VAC14-AS1","DNAAF1","KCNG4","SHISA6","DNAH9","CYP4F3","DMKN","SDC4","AL021578.1","LANCL3","XACT")

# down.1
set <- c("AL391845.2","AL121992.1","GRHL3-AS1","AL603840.1","AL627316.1","AL355306.2","AL445426.1","IGSF3","ATP1A4","CFAP126","LMX1A","F5","SELL","AL392172.2","AL591686.2","AL591686.1","AC106875.1","STON1-GTF2A1L","MAP3K19","FIGN","AC016766.1","LINC01117","TTN","DNAH7","AC069148.1","SPHKAP","RTP5","AC092053.6","AC098969.1","CNTN3","ST3GAL6-AS1","MAATS1","SLC15A2","AC063919.1","ANKUB1","VEPH1","SERPINI2","WDR49","AC026316.5","AC026316.2","DTHD1","AC097709.1","AC093677.3","LINC01088","NDST4","MYOZ2","JADRR","AC107223.1","AC097450.1","IQCM","SFRP2","LINC02427","RANBP3L","AC091435.2","C5orf64","LINC01338","AC113167.1","AC124854.1","LIX1","SLC25A48","AC002428.2","AL021918.3","AL033519.5","AL096709.2","AL589935.3","OGFRL1","AL357522.1","GJA1","AL513524.1","COL28A1","AC005162.2","AC005999.2","AC002451.1","AC105052.5","FOXP2","AC024730.1","SHH","AF131216.1","AF131216.3","FAM167A-AS1","AC037441.2","AC113133.1","TTPA","CALB1","AP000424.1","AC103726.1","AC104257.1","AL160270.1","CNTFR-AS1","FAM189A2","TMC1","AL807761.4","PAPPA","TLR4","AL160272.1","LCN9","AL732437.3","ARMC3","LINC00836","MKX","ARMC4","AL022344.1","RET","AC244230.2","LINC00842","GDF10","STOX1","AL359844.1","CFAP70","AL096706.1","HOGA1","HPSE2","DNMBP-AS1","CFAP43","TECTB","AL133482.1","AC005383.1","PRLHR","AL157388.1","LINC01165","AC136475.4","AC136475.5","AC136475.9","OR51B5","OR51I1","USH1C","MYOD1","LUZP2","DCDC1","AC131571.1","AP001266.1","SLC29A2","CFAP300","AP001527.1","COLCA1","COLCA2","A2ML1","AC092112.1","RERGL","PIK3C2G","LMNTD1","AC022367.1","SSPN","AC012150.2","MGAT4C","DAO","HSPB8","SGCG","LINC02343","AL354833.1","OBI1-AS1","LINC00446","AL445209.1","AL355916.2","PAPLN","SLC24A4","AC012405.1","AC073941.1","AC009754.2","AC103740.1","LINC00927","ADAMTSL3","PLIN1","C1QTNF8","AC007220.1","ACSM5","VWA3A","SDR42E2","HYDIN","TAT","MARVELD3","HPR","AC099508.1","DYNLRB2","TEKT1","CFAP52","AC100793.4","RAMP2-AS1","AC005746.2","AC002546.1","KCNJ16","LDLRAD4-AS1","LINC01894","AC011825.2","SLC14A1","ACSBG2","CRLF1","LINC01532","AC016590.3","HIF3A","MIR663AHG","SYCP2","AP001172.1","AP001065.3","AP000547.3","AP000365.1","Z80897.1","Z98885.2","ARSD-AS1","LINC01546","DGKK","LINC01496","SYTL4","IGSF1","CCDC160")

DotPlot(object = object_1, features = set, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features in: ",sid))



# B :: snRNA LGG 101P ----

rm(sid, object_1)
gc()

sid <- 'van_Hijfte_LGG_P101'
object_1 <- Read10X(data.dir = "data/Glimmunology_LGG_101P/outs/filtered_feature_bc_matrix")
object_1 <- CreateSeuratObject(counts = object_1,
                               min.cells = 3,
                               min.features = 200,
                               project = "glioma_glim")
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1400,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 2000,col="red") +
  geom_hline(yintercept = 14000,col="red")
#scale_y_log10()



object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 1400 &
                     nFeature_RNA < 4500 & 
                     nCount_RNA > 2000 &
                     nCount_RNA < 14000 &
                     percent.mito < 0.025)

object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
#object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot2


## scaling of data ----
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)



## cluster the cells ----

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca",label = TRUE, pt.size = .8, group.by = "seurat_clusters")


ElbowPlot(object_1, ndims = 45)## estimation of the number of principle components in your dataset

d <- 14 # ? 27 #
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


## UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters",label.size=6)



object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
levels(object_1$seurat_clusters) <- gsub("^(7|8|2|1|5|9|13|4|16)$","\\1. T",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|6|11|12)$","\\1. T",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(3|14|15)$","\\1. TAM",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10)$","\\1. OD",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(18)$","\\1. EN|PE",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(17)$","\\1. T/hybrid?",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T", "1. T", "2. T", "4. T", "5. T", "6. T", "7. T", "8. T", "9. T", "11. T","12. T","13. T","16. T", 
  "17. T/hybrid?",
  "10. OD",
  "3. TAM", "14. TAM", "15. TAM",
  "18. EN|PE"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters",label.size=5)  +
  labs(subtitle=sid) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3)))



## FindMarkers ----


tmp.t1.t2 <- FindMarkers(object_1, ident.1 = c(7,8,2,1,5,9,13,4,16),ident.2 = c(0,11,12,6)) # 17, specifically
FeaturePlot(object = object_1, features = c("NFIA","KCNN3","LINC01088","SPARCL1","ADCY2","SLC1A3","LIFR","SSBP2"))
FeaturePlot(object = object_1, features = c("KAZN","DNM3","TNR","GALNT13","RAPGEF4","CADM2","FGF12","APBB2"))




tmp.17a <- FindMarkers(object_1, ident.1 = '17') # 17, specifically
tmp.17b <- FindMarkers(object_1, ident.1 = '17', ident.2 = '18') # 17, compared to EN/PE


# what is cluster 18? seems partially tumor but close to PE/EN? x-check in pca?
FeaturePlot(object = object_1, features = c("GABRB1","ALDH1A1","CADM1","HPSE2","AQP4","RANBP3L","MGST1","RNF219-AS1"))
FeaturePlot(object = object_1, features = c("DTNA","LSAMP","NRCAM","MGAT4C","FRMD5","CDH20","FMN2","SRPX2"))
FeaturePlot(object = object_1, features = c("CTNND2","GPM6A","APC","ADGRB3","PCDH9","TNIK","PPP2R2B","RP11-384F7.2"))
FeaturePlot(object = object_1, features = c("BMPR1B","MIR99AHG","CADM2","MSI2","MAPK10","ANK2","NCAM2")) # Tumor


## 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

## 2. Astrocyte (+) ----

FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



## 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



## 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


## 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor
FeaturePlot(object = object_1, features = "HBG2") # Tumor


## 3D. ? Mono/Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))



## 4. Neurons (+) ----


tmp <- list('C1'=neuron.genes[neuron.genes %in% NPC2 == F] ,
            'NPC1'=NPC1[NPC1 %in% NPC2 == F] ,
            'NPC1+2' = intersect(NPC1, NPC2),
            'NPC2'=NPC2[NPC2 %in% c(NPC1, neuron.genes) == F],
            'NPC2 + C1' = intersect(neuron.genes, NPC2))

DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C1/NPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.pdf"),width=7.5*3, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C1.png"),width=7.5*3, height=4,scale=1.2)




DotPlot(object = object_1, features = c("EGFR", "GFAP","MOG", "PLP1", "TMEM144", 
                                        "RBFOX1", "RBFOX2", "RBFOX3", "CD2",
                                        "CD3D", "P2RY12", "CD163", "ABCB1", "RGS5"
))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


FeaturePlot(object = object_1, features = c("SOX4", "RBFOX3"))


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")


#FeaturePlot(object = object_1, features = "DCN") # DCN
#FeaturePlot(object = object_1, features = "COL1A2") # DCN
FeaturePlot(object = object_1, features = "ANPEP") # DCN


## 5. Oligodendrocytes (+) ----


FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


FeaturePlot(object = object_1, features = "ST18") # OD?
FeaturePlot(object = object_1, features = "MBP") # OD?
FeaturePlot(object = object_1, features = "CTNNA3") # OD?
FeaturePlot(object = object_1, features = "SLC24A2") # OD?
FeaturePlot(object = object_1, features = "KIRREL3") # OD?
FeaturePlot(object = object_1, features = "NKAIN2") # OD?
FeaturePlot(object = object_1, features = "MAP7") # OD?
FeaturePlot(object = object_1, features = "RNF220") # OD?
FeaturePlot(object = object_1, features = "PEX5L") # OD?
FeaturePlot(object = object_1, features = "TMEM144") # OD?
FeaturePlot(object = object_1, features = "EDIL3") # OD?
FeaturePlot(object = object_1, features = "DOCK5") # OD?
FeaturePlot(object = object_1, features = "MOBP") # OD?
FeaturePlot(object = object_1, features = "UNC5C") # OD?
FeaturePlot(object = object_1, features = "CLDN11") # OD?
FeaturePlot(object = object_1, features = "SPOCK3") # OD?
FeaturePlot(object = object_1, features = "CNTNAP4") # OD?
FeaturePlot(object = object_1, features = "MAN2A1") # OD?
FeaturePlot(object = object_1, features = "PCSK6") # OD?
FeaturePlot(object = object_1, features = "TTLL7") # OD?

FeaturePlot(object = object_1, features = "OLIG2") # OD?


## 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



## 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))




## 7. Cycling cells (?) ----


FeaturePlot(object = object_1, features = "TOP2A" )
FeaturePlot(object = object_1, features = "AURKA" )
FeaturePlot(object = object_1, features = "AURKB" )
FeaturePlot(object = object_1, features = "BUB1" )
FeaturePlot(object = object_1, features = "BUB1B" )
FeaturePlot(object = object_1, features = "CDC20" )
FeaturePlot(object = object_1, features = "CENPF" )
FeaturePlot(object = object_1, features = "FAM64A" )
FeaturePlot(object = object_1, features = "FOXM1" )
FeaturePlot(object = object_1, features = "TACC3" )
FeaturePlot(object = object_1, features = "TMPO" )
FeaturePlot(object = object_1, features = "TPX2" )
FeaturePlot(object = object_1, features = "TUBA1C" )



sig <- c("H4C1", "H3C2", "H2AC4", "H3C3",            "HIST1H4A","HIST1H3B", "HIST1H2AB", "HIST1H3C",
         "H1-6", "H3C7", "H2BC9", "H2BC11",          "HIST1H1T","HIST1H3F","HIST1H2BH","HIST1H2BH",
         "H2AC11", "H2BC12", "H2AC12", "H2BC13",     "HIST1H2AG","HIST1H2BK","HIST1H2AH","HIST1H2BL",
         "H2AC13", "H3C10", "H2AC14", "H2BC14",      "HIST1H2AI","HIST1H3H","HIST1H2AJ","H2BC14",
         "H2AC15", "H2AC16", "H1-5", "H3C11",        "HIST1H2AK","HIST1H2AL","HIST1H1B","HIST1H3I",
         "H3C12", "H2BC17",                          "HIST1H3J","HIST1H2BO")
sig <- unique(sig)
FeaturePlot(object = object_1, features = sig )





