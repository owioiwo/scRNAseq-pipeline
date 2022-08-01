##script 3 : normalize & deepclustering
## input = 2output_seurat_obj.rds ( = processed seurat_obj)
## output = project_name.rds (= prcessed seurat_clustering)

#---inputdata
seurat_obj = readRDS("2output_seurat_obj.rds")

#--select normalize method
cat("choose normalize method:\n 1 = LogNormalize (general method)\n 2 = SCTranform (better for UMIs methods,such as 10X scRNA-seq)\n")
repeat{
  normalize_method = readline("1 or 2:")
  if(normalize_method %in% c("1","2")){break}
}
#--select normalize method end
#--select seurat pipeline mode
cat("choose seurat pipeline mode :\n 1 = standard (this mode use 2000~3000 features as high variable genes)\n 2 = filter (this mode use all features as high variable genes,used in some papers) \n")
repeat{
  pipeline_method = readline("1 or 2:")
  if(pipeline_method %in% c("1","2")){break}
}
#--select seurat pipeline mode end
#==do seurat pipeline & deep clustering to find/delete outliers
seurat_clustering <- seurat_pipeline(seurat_obj,normalize_method,pipeline_method) %>% find_cluster()
repeat{
  cat("Next we can do :\n 1 = QC by deleting bad clusters \n 2 = QC by deleting subcluster's outliers (median ~ 3*sd)\n 3 = refind clusters and choose again \n 4 = result ok,finish at once \n")
  repeat{
    qc_choice = readline("1 or 2 or 3 or 4 : ")
    if(qc_choice %in% c("1","2","3","4")){break}
  }
  if(qc_choice == 1){
    seurat_clustering = delete_cluster(seurat_clustering)
  } else if(qc_choice == 2){
    seurat_clustering = subcluster_sd(seurat_clustering)
  } else if(qc_choice == 4){
    break
  } else {
    seurat_clustering = CreateSeuratObject(seurat_clustering@assays[["RNA"]]@counts) %>% PercentageFeatureSet(pattern = "^mt-|^MT-|^Mt-", col.name = "percent.mt") %>% seurat_pipeline(normalize_method,pipeline_method) %>% find_cluster()
  } 
}

## save result
seurat_clustering$sample = seurat_obj@project.name
seurat_clustering@project.name = seurat_obj@project.name
deep_cluster_result = paste0(seurat_obj@project.name,".rds")
saveRDS(seurat_clustering,deep_cluster_result)
cat("===== Deep clustering complete: \n===== result save as project name set before.")
## == end ==