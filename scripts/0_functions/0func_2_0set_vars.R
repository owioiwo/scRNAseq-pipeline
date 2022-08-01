## function project_name,set_qc_vars
# for 3_deepclustering.R

# set_project_name
set_project_name = function(){
  cat("input a project name,which will be finally save as the output's name:\n")
  repeat{
    project_name = readline("input a project name,for example 'sepsis1':") 
    if(project_name != ""){break}
  }
  return(project_name)
}
# set_project_name end

#== set_qc_vars
set_qc_vars = function(vars){
  repeat{
    cat("set",vars,"threshold for filter,default is 200:\n blank = default(200) \n y ='ok,200'\n n ='not ok'\n")
    notreset = readline("y or n:")
    if(notreset %in% c("y","n","")){break}
  }
  if(notreset != "n"){
    vars = 200
  } else{
    repeat{
      cat("reset",vars,"QC threshold,press enter with a blank line will set to default (200):\n")
      notreset= readline(prompt = "input a number:")
      if(notreset == ""){
        vars = 200
        break
      }
      if(!sum(!unlist(strsplit(notreset,split = "")) %in% c(0:9))){
        vars= as.numeric(notreset)
        break
      }
    }
  }
  return(vars)
}
#== set_ncount_nfeature end