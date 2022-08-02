## function gene.plots & input.plot.genes
# for 5_gene_cell_analysis.R : 

## alt function: input.plot.genes
input.plot.genes = function(plot.choose,default.genes){
  cat("input gene names to plot one by one,end by press enter with blank..\n Note:press 'enter' without any input will set plot genes to top markers randomly. \n")
  plot.feature = scan(what = character())
  if(length(plot.feature) == 0){
    plot.feature = switch(plot.choose,
           sample(1:length(default.genes),min(10,length(default.genes))),
           sample(1:length(default.genes),min(1,length(default.genes))),
           sample(1:length(default.genes),min(1,length(default.genes))),
           sample(1:length(default.genes),min(20,length(default.genes))),
           sample(1:length(default.genes),min(1,length(default.genes)))
    ) %>% sort()
    plot.feature = default.genes[plot.feature]
  } else {
    plot.feature = unique(plot.feature)
    }
  return(plot.feature)
}
## alt end

####------function gene.plots
gene.plots = function(seurat.object,default.genes){
  repeat{
    cat("choose a plot to draw: \n 1 = DotPlot \n 2 = FeaturePlot \n 3 = VlnPlot \n 4 = Heatmap \n 5 = RidgePlot \n 6 = plots ok,let's quit")
    repeat{
      plot.choose = readline("1 or 2 or 3 or 4 or 5 or 6 : ")
      if(plot.choose %in% c("1","2","3","4","5","6")){
        plot.choose = as.numeric(plot.choose)
        break
        }
    }
    if(plot.choose == 6){
      cat("gene expression plot complete..\n")
      break
    } else {
      plot.feature = input.plot.genes(plot.choose,default.genes)
      }
    if(plot.choose == 1){
      print(DotPlot(seurat.object,features = plot.feature) + RotatedAxis())
    } else if(plot.choose == 2){
      print(FeaturePlot(seurat.object,features = plot.feature))
    } else if(plot.choose == 3){
      print(VlnPlot(seurat.object,features = plot.feature))
    } else if(plot.choose == 4){
      print(DoHeatmap(seurat.object,features = plot.feature) + NoLegend())
    } else {
      print(RidgePlot(seurat.object,features = plot.feature))
    }
    cat("plot drwan.\n")
  }
}
# gene.plot end
