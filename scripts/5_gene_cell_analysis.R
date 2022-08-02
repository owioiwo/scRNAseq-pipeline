##script 5 : do gene de analysis and cell annotation
## input = any seurat object list in "5_gene_cell_analysis.txt" at "results" directory
## output 1 = data_name.rds (prcessed seurat.object)
## output 2 = markers_data_name.rds.txt (seurat.markers =all markers)

# 1 import data
options(warn = -1)
inputdata_list = read.table("5_gene_cell_analysis.txt")
options(warn = 0)
cat("data list: \n")
print(inputdata_list)
repeat{
  cat("choose a data to analysis,\n")
  reso= readline(prompt = "input the number of that data:")
  if(reso ==""){
  } else if(!sum(!unlist(strsplit(reso,split = "")) %in% c(0:9))){
    reso= as.numeric(reso)
    if(reso > 0 & reso <= nrow(inputdata_list)){
      break
    } else{
      cat("input number is wrong... \n")
    }
  } 
}
seurat.object = readRDS(inputdata_list[reso,])
print_plots(seurat.object)
cat("data loaded,check summary and two plots.. \n")
# 2 find markers
seurat.markers = FindAllMarkers(seurat.object)
repeat{
  cat("show how many top differential expressed genes in eatch cluster,default is 5: \n")
  top.number = readline(prompt = "input a top exp number:")
  if(top.number ==""){
    top.number = 5
    break
  } else if(!sum(!unlist(strsplit(top.number,split = "")) %in% c(0:9))){
    top.number= as.numeric(top.number)
    if(top.number > 0){
      break
    } else{
      cat("input number is wrong... \n")
    }
  }
}
top.markers = seurat.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(n = top.number, order_by = avg_log2FC)
top.markers = top.markers[,"gene"] %>% unique()
cat("plot top",top.number,"makers in each cluster,check : \n",top.markers)
print(DotPlot(seurat.object,features = top.markers) + RotatedAxis())
write.table(seurat.markers,paste0("markers_",inputdata_list[reso,],".txt"),sep = "\t",quote = F,row.names = F)
cat("find markers done, save all markers as file: ", paste0("markers_",inputdata_list[reso,],".txt"),“ \n")
# 3 plot marker genes and annotate cells
repeat{
  cat("Next we can do: \n 1 = plot gene expression by multi functions \n 2 = annotate cells by mannual \n 3 = annotate cells by R package: SingleR \n 4 = end anaylsis \n")
  repeat{
    qc_choice = readline("1 or 2 or 3 or 4 : ")
    if(qc_choice %in% c("1","2","3","4")){break}
  }
  if(qc_choice == 1){
    gene.plots(seurat.object,top.markers)
  } else if(qc_choice == 2 | qc_choice == 3){
    seurat.object = cell.annotation(seurat.object,qc_choice)
  } else {
    UMAPPlot(seurat.object) -> plot1
    UMAPPlot(seurat.object, group.by = "sample") -> plot2
    print(plot1 + plot2)
    cat("gene de analysis and cell annotation done,check plots..\n")
    break
  } 
}
##
cat("save result as file: ",inputdata_list[reso,],“ \n")
saveRDS(object = seurat.object,inputdata_list[reso,])
# end
