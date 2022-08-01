## functions find_cluster
# for 3_deepclustering.R

####--function find_cluster
find_cluster = function(object){
  cat("cycles begin to find a best resolution for clustering...\n")
  repeat{
    repeat{
      cat("Input a resolution number to try clustering,\npress 'enter' directly will set resolution to 0.8\n")
      reso= readline(prompt = "input a resolution number:")
      if(reso == ""){
        reso=0.8
        break
      }
      if(!sum(!unlist(strsplit(reso,split = "")) %in% c(0:9,"."))){
        reso= as.numeric(reso)
        break
      } 
    }
    object <- FindClusters(object, resolution = reso)
    print_plots(object)
    cat("Found",length(levels(object)),"clusters by resolution =",reso,", see Vlnplot & UMAPPlot & summary\n")
    cat("Want to refind clusters by reset resolution? \n y = yes \n n = no\n")
    repeat{
      reset = readline("y or n:")
      if(reset %in% c("y","n")){break}
    }
    if(reset== "n"){
      cat("find clusters done.\n")
      break
    }
  }
  return(object)
}
####--function find_cluster end
