## function cell.annotation
# for 5_gene_cell_analysis.R : 

####------function cell.annotation
cell.annotation = function(seurat.object,qc_choice){
  if(qc_choice == 2){
    repeat{
      new.cluster.ids = as.character()
      for (i in 1:length(levels(seurat.object))) {
        cat("annotate cells manually,input cell type of cluster",i-1,"..\n")
        repeat{
          cell.type = readline(paste0("cell type of cluster ",i-1,": "))
          if(cell.type != ""){
            new.cluster.ids[i] = cell.type
            break
          }
        }
      }
      names(new.cluster.ids) <- levels(seurat.object)
      seurat.object = RenameIdents(seurat.object, new.cluster.ids)
      print_plots(seurat.object)
      cat("cell annotation done,need to do this again?\n y = yes\n n = no \n")
      repeat{
        re.anno = readline("y or n : ")
        if(re.anno %in% c("y","n")){
          break
        }
      }
      if(re.anno == "n"){
        cat("cell annotate manually complete..\n")
        break
      }
    } #manual end here
  } else {
    repeat{
      cat("which species to annotate? \n 1 = Homo sapiens (human) \n 2 = Mus musculus (house mouse)\n")
      repeat{
        anno.species = readline("1 or 2 : ")
        if(anno.species %in% c("1","2")){
          break
        }
      }
      if(anno.species == 1){
        anno.ref = celldex::HumanPrimaryCellAtlasData()
      } else {
        cat("choose a mouse reference database.. \n 1 = MouseRNAseqData \n 2 = ImmGenData \n")
        repeat{
          anno.mouse = readline("1 or 2 : ")
          if(anno.mouse %in% c("1","2")){
            break
          }
        }
        if(anno.mouse == 1){
          anno.ref = celldex::MouseRNAseqData()
        } else {
          anno.ref = celldex::ImmGenData()
        }
      } # mouse ref
      anno.data = GetAssayData(seurat.object,slot = "data")
      anno.anno = SingleR(test = anno.data,ref = anno.ref,labels = anno.ref$label.main)
      cat("check SingleR annotation results & plots : \n")
      print(table(anno.anno$labels))
      print(plotScoreHeatmap(anno.anno))
      # plotScoreDistribution(anno.anno)
      # DoHeatmap(seurat.object,features = unique(unlist(anno.anno@metadata$de.genes$Fibroblasts)))
      Idents(seurat.object) = anno.anno$labels
      print_plots(seurat.object)
      cat("SingleR annotation done,want to change a reference database? \n y = yes \n n = no \n")
      repeat{
        re.anno = readline("y or n : ")
        if(re.anno %in% c("y","n")){
          break
        }
      }
      if(re.anno == "n"){
        cat("SingleR annotation complete..\n")
        break
      }
    }
  } # singleR end here
  return(seurat.object)
}
#end