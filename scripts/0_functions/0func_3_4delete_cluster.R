## function delete_cluster
# for 3_deepclustering.R : 

####------function delete_cluster
delete_cluster= function(object){
  repeat{
    cat("Need delete some clusters? if yes,input the NO. of that cluster(s),or if not press 'enter' directly.\n")
    repeat{
      badcluster = scan("",what = "")
      if(length(grep("[^0-9]",badcluster)) != 0){
        cat("inputs are wrong, integers only,re input bad cluster No.:\n")
      } else{
        break
      }
    }
    if(length(badcluster)==0){
      cat("no bad cluster No. detected,deleting step complete.\n")
      break
    } 
    object <- subset(object, cells= colnames(object)[!(object$seurat_clusters %in% badcluster)])
    print_plots(object)                   
    cat("new plots drawn,need repeat to delete bad clusters? \n y = yes \n n = no\n")
    repeat{
      loop = readline("y or n:")
      if(loop %in% c("y","n")){break}
    }
    if(loop == "n"){
      cat("cluster deleting step complete.\n")
      break
    }
  }
  return(object)
}
####------function delete_cluster end