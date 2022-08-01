## functions print_plots
# for 3_deepclustering.R

##==function print_plots
print_plots = function(object){
  color = rainbow(length(levels(object)))
  print(VlnPlot(object, c("nFeature_RNA", "nCount_RNA","percent.mt"),cols = color))
  print(UMAPPlot(object,cols = color) + labs(subtitle=paste0("resolution=",object@commands$FindClusters$resolution)))
  print("summary => number of cells in eath cluster:",quote =F)
  print(summary(object$seurat_clusters))
}
##==function print_plots end