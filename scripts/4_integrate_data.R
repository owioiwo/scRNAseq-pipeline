##script 4 :  Integration
## input = integrate_list of [project_name.rds]s ( = processed [seurat_clustering]s)
## output =  project_name.rds (= prcessed data.integrate)
options(warn = -1)
integrate.list = read.table("4_integrate_list.txt") %>% unique()
options(warn = 0)
cat("check data list to be integrated: \n")
print(integrate.list)
if(nrow(integrate.list) <=1){
  cat("Not enough dataset to integrate.")
} else{
  data.list = list()
  for (i in 1:nrow(integrate.list)) {
    data.list[[i]] = readRDS(integrate.list[i,])
  }
}
project.name = set_project_name()
#--select normalize method
cat("choose normalize method:\n 1 = LogNormalize (general method)\n 2 = SCTranform (better for UMIs methods,such as 10X scRNA-seq)\n")
repeat{
  normalize.method = readline("1 or 2:")
  if(normalize.method %in% c("1","2")){break}
}
#--choose integrate method
cat("choose integrate method:\n 1 = reciprocal PCA (fast,better for different disease states between samples)\n 2 = canonical CCA (slow,maybe overcorrection,better for conserved celltypes)\n")
repeat{
  integrate.method = readline("1 or 2:")
  if(integrate.method %in% c("1","2")){break}
}

#====integrate pipeline
data.integrate = integrate.pipeline(data.list,normalize.method,integrate.method)

#===reset cluster numbers
cat("Need to reset rosulution and refind clusters:? \n y = yes,reset resolution \n n = no,result is ok \n")
repeat{
  reset = readline("y or n:")
  if(reset %in% c("y","n")){break}
}
if(reset== "n"){
  cat("find clusters done...\n")
} else{
  data.integrate = find_cluster(data.integrate)
}
## output
sub.title = paste0(" methods = ", Command(data.integrate,"FindIntegrationAnchors","normalization.method")," + ",Command(data.integrate,"FindIntegrationAnchors","reduction"),"\n k.anchor = ",Command(data.integrate,"FindIntegrationAnchors","k.anchor")," ; ","resolution = ",Command(data.integrate,"FindClusters","resolution"))
print(UMAPPlot(data.integrate, group.by = "sample") + labs(subtitle= sub.title))
cat("======= data integration complete =======\n")
data.integrate@project.name = project.name
saveRDS(data.integrate,paste0(project.name,".rds"))
# end