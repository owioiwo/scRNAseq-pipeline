## function create_seurat_obj，quantile_QC
# for 2_data_initial.R

#---alt function: quantile_QC
# for quantile gene & count
quantile_QC = function(object){
  cat("-------QC quantiles start-------:\n")
  cat("nCount(cells):\n")
  print(quantile(colSums(object@assays[["RNA"]]@counts)))
  cat("nFeature(cells):\n")
  print(quantile(colSums(object@assays[["RNA"]]@counts >0)))
  cat("nCount(genes):\n")
  print(quantile(rowSums(object@assays[["RNA"]]@counts)))
  cat("nFeature(genes) == min_cells:\n")
  print(quantile(rowSums(object@assays[["RNA"]]@counts >0)))
  cat("-------QC quantiles end-------\n")
  print(VlnPlot(object,features = c("nFeature_RNA","nCount_RNA","percent.mt")))
  print(FeatureScatter(object,feature1 = "nCount_RNA",feature2 = "nFeature_RNA"))
  cat("VlnPlot(nCount,nFeature,mt) & plot(nCount,nFeature) are drawn.\n")
}
#---alt function: end
#----create seurat obj and initial
create_seurat_obj= function(object,project_name,min_cells,min_features,ncount,nfeature){
  cat("dims of inputs:",dim(object),"\n")
  object = CreateSeuratObject(project =project_name,min.cells=min_cells,min.features=min_features,object)
  cat("dims of cr_obj:",dim(object),"\n")
  object[["percent.mt"]]<- PercentageFeatureSet(object = object,pattern = "^mt-|^MT-|^Mt-")
  object<- subset(object,subset = nFeature_RNA >nfeature & nCount_RNA >ncount)
  cat("dims of subset:",dim(object),"\n")
  cat("object created.\n")
  return(object)
}
