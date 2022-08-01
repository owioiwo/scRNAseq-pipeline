##function median_3sd_QC
# for 2_data_initial.R

#---cal & QC median~3*sd
median_3sd_QC = function(object){
  repeat{
    medcount_min = median(object$nCount_RNA) -3*sd(object$nCount_RNA)
    medcount_max = median(object$nCount_RNA) +3*sd(object$nCount_RNA)
    medfeature_min = median(object$nFeature_RNA)- 3*sd(object$nFeature_RNA)
    medfeature_max = median(object$nFeature_RNA)+ 3*sd(object$nFeature_RNA)
    medmt_min = median(object$percent.mt) -3*sd(object$percent.mt)
    medmt_max = median(object$percent.mt) +3*sd(object$percent.mt)
    #------
    if(max(object$nCount_RNA) < medcount_max & max(object$nFeature_RNA) < medfeature_max & min(object$nCount_RNA) > medcount_min & min(object$nFeature_RNA) > medfeature_min & max(object$percent.mt) < medmt_max & min(object$percent.mt) > medmt_min ){
      break
    } else{ print("median with 3*sd QC +1",quote =F)}
    #------
    object<- subset(object,subset = nFeature_RNA >medfeature_min & nFeature_RNA < medfeature_max & nCount_RNA > medcount_min & nCount_RNA < medcount_max & percent.mt < medmt_max & percent.mt > medmt_min)
  }
  cat("nCount:",c(min(object$nCount_RNA),max(object$nCount_RNA)),"\n")
  cat("nFeature:",c(min(object$nFeature_RNA),max(object$nFeature_RNA)),"\n")
  cat("mt%:",c(min(object$percent.mt),max(object$percent.mt)),"\n")
  cat("dims:",dim(object),"\n")
  #---cal & QC end
  cat("median with 3sd QC done.\n")
  return(object)
}
