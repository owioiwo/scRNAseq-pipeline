## function integrate.pipeline: integrate sct & lognormal + rpca & cca
# for 4_integrate_data.R : 

####------function integrate.pipeline
integrate.pipeline = function(data.list,normalize.method,integrate.method){
  # 1 normalize
  if(normalize.method == 1){
    data.normalized <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x) %>% FindVariableFeatures()
    })
  } else {
    data.normalized <- lapply(X = data.list, FUN = SCTransform, method = "glmGamPoi")
  }
  # 2 select features & prepare
  var.nfeature = ifelse(normalize.method == 1,2000,3000)
  integrate.feature <- SelectIntegrationFeatures(object.list = data.normalized, nfeatures = var.nfeature)
  if(normalize.method == 2){
    data.normalized <- PrepSCTIntegration(object.list = data.normalized, anchor.features = integrate.feature)
  }
  # 3 rpca integrate mode only: RunPCA
  reduction.method = "cca"
  if(integrate.method == 1){
    reduction.method = "rpca"
    if(normalize.method == 1){
      data.normalized <- lapply(X = data.normalized,FUN = function(x){
        x <- ScaleData(x, features = integrate.feature) %>% RunPCA(features = integrate.feature,verbose = F)
      })
    } else{
      data.normalized <-lapply(X = data.normalized, FUN = RunPCA, features = integrate.feature,verbose = F)
    }
  }
  # 4 perform integration
  Log.SCT = ifelse(normalize.method == 1,"LogNormalize","SCT")
  k.anchor = 5
  repeat{
  data.anchors <- FindIntegrationAnchors(object.list = data.normalized, anchor.features = integrate.feature, reduction = reduction.method,normalization.method = Log.SCT,k.anchor = k.anchor)
  data.integrate <- IntegrateData(anchorset = data.anchors,normalization.method = Log.SCT)
  # 5 Lognormalize only: Scaledata
  if(normalize.method == 1){
    data.integrate <- ScaleData(data.integrate, verbose = FALSE)
  }
  # 6 downstream anlysis 
  data.integrate <- RunPCA(data.integrate,verbose= F) %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters()
  # 7 plots and return result
  sub.title = paste0(" methods = ",ifelse(normalize.method == 1,"log","sct")," + ",reduction.method,"\n k.anchor = ",k.anchor)
  print(UMAPPlot(data.integrate) + labs(subtitle= sub.title))
  print(UMAPPlot(data.integrate, group.by = "sample") + labs(subtitle= sub.title))
  # 8  reset integrate strength or complete 
    cat("=====================================\nintegrate pipeline complete,TWO UMAP-plots drawn,check and decide whether the result is ok or not,if not,we can reset the k.anchor parameter and modify the strength of integration..\n so next we can do: \n 1 = reset k.anchor to try again \n 2 = result is ok,end this pipeline \n")
    repeat{
      reset_par = readline("1 or 2: ")
      if(reset_par %in% c("1","2")){break}
    }
    if(reset_par == 1){
      repeat{
        cat("input a number as new k.anchor: \npress 'enter' directly will set to default '5' \n")
        k.anchor = readline(prompt = "input a number as new k.anchor:")
        if(k.anchor == ""){
          k.anchor = 5
          break
        }
        if(!sum(!unlist(strsplit(k.anchor,split = "")) %in% c(0:9,"."))){
          k.anchor = as.numeric(k.anchor)
          break
        } 
      }
    } else {
      cat("integrate pipeline complete...\n")
      break
    }
  }
  return(data.integrate)
}

# end

