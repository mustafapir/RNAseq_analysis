
library(Seurat)
library(data.table)


raw_counts<-ReadMtx(mtx = "files/GSM6067321_H1N1_D7_matrix.mtx.gz", features = "files/GSM6067321_H1N1_D7_features.tsv.gz",
  cells = "files/GSM6067321_H1N1_D7_barcodes.tsv.gz")

sc167<-ReadMtx(mtx = "files/COPD/GSM5100998_SC167matrix.mtx.gz", features = "files/COPD/GSM5100998_SC167genes.tsv.gz",
              cells = "files/COPD/GSM5100998_SC167barcodes.tsv.gz")

sc170<-ReadMtx(mtx = "files/COPD/GSM5100999_SC170matrix.mtx.gz", features = "files/COPD/GSM5100999_SC170genes.tsv.gz",
               cells = "files/COPD/GSM5100999_SC170barcodes.tsv.gz")

sc200<-ReadMtx(mtx = "files/COPD/GSM5101000_SC200matrix.mtx.gz", features = "files/COPD/GSM5101000_SC200genes.tsv.gz",
               cells = "files/COPD/GSM5101000_SC200barcodes.tsv.gz")

sc167<-CreateSeuratObject(counts = sc167, project = "sc167")
sc170<-CreateSeuratObject(counts = sc170, project = "sc170")
sc200<-CreateSeuratObject(counts = sc200, project = "sc200")

sc167@meta.data$type<-"sc167"
sc170@meta.data$type<-"sc170"
sc200@meta.data$type<-"sc200"

object<-list(sc167, sc170, sc200)

object<- lapply(X = object, FUN = function(x) {
  x[["percent.mt"]]<-PercentageFeatureSet(x, pattern = "^MT-")
  x@meta.data<-x@meta.data[!is.nan(x@meta.data$percent.mt),]
  x<-subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features<-SelectIntegrationFeatures(object.list = object)
anchors<-FindIntegrationAnchors(object.list = object, anchor.features = features)
combined<-IntegrateData(anchorset = anchors)

all.genes<-rownames(combined)
combined<-ScaleData(combined, features = all.genes)
combined<-RunPCA(combined, npcs = 30)

ElbowPlot(combined, ndims = 30)

combined<-RunUMAP(combined, reduction = "pca", dims = 1:24)
combined<-FindNeighbors(combined, reduction = "pca", dims = 1:24)
combined<-FindClusters(combined, resolution = 0.5)

DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)

DefaultAssay(combined)<-"RNA"
markers1<-FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

ciliary1<-markers1[markers1$cluster == 6,]

ciliary2<-FindConservedMarkers(combined, ident.1 = 6)

object[["percent.mt"]]<-PercentageFeatureSet(object, pattern = "^MT-")

object@meta.data<-object@meta.data[!is.nan(object@meta.data$percent.mt),]
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#f<-fread("files/GSM6067321_H1N1_D7_features.tsv.gz")

object<-subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

object<-NormalizeData(object)

object<-FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10<-head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1<-VariableFeaturePlot(object)
plot2<-LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes<-rownames(object)
object<-ScaleData(object, features = all.genes)

object<-RunPCA(object, features = VariableFeatures(object = object))
object<-RunPCA(object)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
DimPlot(object, reduction = "pca")
DimHeatmap(object, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(object, ndims = 50)

object<-FindNeighbors(object, dims = 1:30)
object<-FindClusters(object, resolution = 0.5)

object<-RunUMAP(object, dims = 1:30)
DimPlot(object, reduction = "umap", label = TRUE)

object<-RunTSNE(object, dims = 1:22)
DimPlot(object, reduction = "tsne")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
object.markers<-FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
object.markers<-object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

