## functions seurat_pipeline
# for 3_deepclustering.R : 

seurat_pipeline = function(object,normalize_method,pipeline_method){
  var_nfeatures = nrow(object)
  if(pipeline_method == 1){
    var_nfeatures = ifelse(normalize_method == 1,2000,3000)
  }
  if(normalize_method == 1){
    dims= 1:20
    object <- NormalizeData(object) %>%  FindVariableFeatures(nfeature = var_nfeatures) %>% ScaleData()
  } else {
    dims= 1:30
    object <- SCTransform(object, method = "glmGamPoi", vars.to.regress = "percent.mt",variable.features.n = var_nfeatures)
  }
  object <- RunPCA(object) %>% RunUMAP(dims = dims) %>% FindNeighbors(dims = dims) %>% FindClusters()
  method= ifelse(normalize_method == 1,"LogNormal","SCTrans" )
  print(UMAPPlot(object,label = T) + labs(subtitle=paste0(" method= ",method,"\n var_nfeatures= ",length(VariableFeatures(object)))))
  cat("++++++ seurat pipeline complete,UMAPPlot drawn with dims:",dim(object),"++++++\n")
  return(object)
}
## end