## function subcluster_sd
# for 3_deepclustering.R : 

####------function subcluster_sd
subcluster_sd = function(object){
  cat("subcluser QC,median with 3*sd start, dims:",dim(object),"\n")
  repeat{
    cat("which subcluster'outliers in the plots would you like to delete this time: \n")
    repeat{
      cat("subcluster No. : ",levels(object),"\n")
      cluster_NO = readline("choose a subcluster NO. : ")
      if(cluster_NO %in% levels(object)){break} 
    }
    sub_cluster = subset(object,subset = seurat_clusters == cluster_NO)
    sub_embeddings = Embeddings(sub_cluster,reduction = "umap")
    repeat{
      medembed_min1 = median(sub_embeddings[,1]) -3*sd(sub_embeddings[,1])
      medembed_max1 = median(sub_embeddings[,1]) +3*sd(sub_embeddings[,1])
      medembed_min2 = median(sub_embeddings[,2]) -3*sd(sub_embeddings[,2])
      medembed_max2 = median(sub_embeddings[,2]) +3*sd(sub_embeddings[,2])
      #------
      if(max(sub_embeddings[,1]) < medembed_max1 & max(sub_embeddings[,2]) < medembed_max2 & min(sub_embeddings[,1]) > medembed_min1 & min(sub_embeddings[,2]) > medembed_min2){
        break
      } else{
        sub_embeddings = sub_embeddings[sub_embeddings[,1] < medembed_max1 & sub_embeddings[,1] > medembed_min1 & sub_embeddings[,2] < medembed_max2 & sub_embeddings[,2] > medembed_min2,]
        print("subcluser QC,median with 3*sd +1",quote =F)
      }
    }
    sub_outlier = setdiff(colnames(sub_cluster),rownames(sub_embeddings))
    object_inlier = setdiff(colnames(object),sub_outlier)
    object = subset(object,cells = object_inlier)
    print_plots(object)
    cat("subcluster",cluster_NO,"QC done,new plots drawn,dims:",dim(object),"\n")
    cat("any other subclusters' outliers need to be delete? y = yes, n = no \n")
    repeat{
      re_qc = readline("y or n : ")
      if(re_qc %in% c("y","n")){break} 
    }
    if(re_qc == "n" ){break}
  }
  cat("subcluser QC,median with 3*sd done. final dims:",dim(object),"\n")
  return(object)
}
####------function delete_cluster end