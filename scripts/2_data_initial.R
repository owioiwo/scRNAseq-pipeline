##script 2 :
## input = 1output_check_gene_names.rds ( = processed inputdata)
## output = 2output_seurat_obj.rds ( = processed seurat_obj)

#---inputdata
inputdata = readRDS("1output_check_gene_names.rds")
#---start here if inputdata from other script
#-----delete duplicated & empty gene names
cat("dim of inputdata:",dim(inputdata),"\n")
inputdata = inputdata[!duplicated(rownames(inputdata)),]
cat("dim without duplicate:",dim(inputdata),"\n")
inputdata = inputdata[rownames(inputdata)!="",]
cat("dim without blank gene:",dim(inputdata),"\n")
#-----delete end
#===set QC vars
project_name = set_project_name()
min_cells= 0.01*ncol(inputdata) # initial filter: genes in >1% cells kept
min_features =200 #initial filter,cell > nfeature to kept
ncount= set_qc_vars("ncount")
nfeature = set_qc_vars("nfeature")
#===set var end
colnames(inputdata) = paste(project_name,colnames(inputdata),sep = "_")
seurat_obj = create_seurat_obj(inputdata,project_name,min_cells,min_features,ncount,nfeature)
quantile_QC(seurat_obj)
## cycles inital
init_loop = 1
repeat{
  seurat_obj = median_3sd_QC(seurat_obj) %>% lmfit_3sd_QC()
  seurat_obj = create_seurat_obj(seurat_obj@assays[["RNA"]]@counts,project_name,min_cells,min_features,ncount,nfeature)
  if(min(seurat_obj$nCount_RNA) >ncount & min(seurat_obj$nFeature_RNA) >nfeature & min(rowSums(seurat_obj@assays[["RNA"]]@counts >0)) >= min_cells ){
    break
    }
  cat("==============3sd QC again==============\n")
  init_loop =init_loop +1
}
quantile_QC(seurat_obj)
cat("#####data initial done by",init_loop,"cycles.#####\n")
cat("dims of initial data:",dim(seurat_obj),"\n")
## cycles inital end
#output ===      = seurat_obj
seurat_obj$sample = project_name
saveRDS(seurat_obj,"2output_seurat_obj.rds")

## pinpeline script end here.
##===there is an extra script to explore abnormal genes,see "extra" folder