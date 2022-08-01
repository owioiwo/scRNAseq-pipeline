##function lmfit_3sd_QC
# for 2_data_initial.R
# lmfit with residual's 3sd QC

lmfit_3sd_QC =function(object){
  repeat{
    fit = lm(object$nFeature_RNA~object$nCount_RNA)
    res = residuals(fit)
    fit_res_sd_max = median(res) + 3*sd(res)
    fit_res_sd_min = median(res) - 3*sd(res)
    object$res = res
    #------
    if(max(object$res) < fit_res_sd_max & min(object$res) > fit_res_sd_min ){
      break
    } else{ print("lm res median with 3*sd QC +1..",quote =F)}
    #------
    object<- subset(object,subset = res >fit_res_sd_min & res < fit_res_sd_max)
  }
  cat("dims:",dim(object),"\n")
  cat("lmfit residual:",fit_res_sd_min,fit_res_sd_max,"\n")
  cat("lmfit with residual 3sd QC done.\n")
  return(object)
}
