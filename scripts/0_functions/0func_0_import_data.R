# function: import_data
# for script: main.R
## data should be in "rawdata" directory
# data name in 0_inputdata_list.txt, one line = one data name
#only support gene matrix or seurat obj files.

import_data = function(){
  options(warn = -1)
  inputdata_list = read.table("../rawdata/0_inputdata_list.txt")
  options(warn = 0)
  cat("data list: \n")
  print(inputdata_list)
  repeat{
    cat("choose a data to process,\n")
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
  postfix = inputdata_list[reso,] %>% substr(nchar(inputdata_list[reso,])-2,nchar(inputdata_list[reso,])) %>% toupper()
  if(postfix == "RDS"){
    inputdata= readRDS(paste0("../rawdata/",inputdata_list[reso,]))
    # inputdata= inputdata@assays[["RNA"]]@counts
    inputdata = GetAssayData(inputdata,slot = "counts")
  } else if(postfix == "TXT"){
    inputdata= read.table(paste0("../rawdata/",inputdata_list[reso,]),sep= "\t",header =T)
  } else {
    cat("only support 'rds' or 'txt' files,others can import by hand.. " )
  }
  return(inputdata)
}
# end