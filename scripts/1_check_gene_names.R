##script1: 
##input = inputdata; 
##output = 1output_check_gene_names.rds ( = processed inputdata )
gene_names = rownames(inputdata)

#== find special characters
special_char= strsplit(gene_names,"") %>% unlist() %>% unique()
special_char= special_char[!special_char %in% c(0:9,letters,LETTERS)]
cat("special characters in gene names:",special_char,"\n")
special_char_backup = special_char
special_char = special_char[!special_char %in% c("-","(",")",".")]#,"_")]
cat("not good special characters:",special_char,"\n")

#== find special characters end
#-----find out gene names with special_char 
if(!length(special_char)){
  print("no special characters,gene names seems good",quote = F)
} else { 
  special_gene= as.numeric()
  for (i in special_char) {
     temp= grep(i,gene_names,fixed = T)
     special_gene = c(special_gene,temp)
                           }
  print("not good special character gene names detected:",quote = F)
  print(c("special gene names:",gene_names[special_gene]),quote = F)
  print(c("in positons:",special_gene),quote = F)
        }
#-----find out special char end
#====================find gene names transform to date
# b= c("1-Mar","111-Mar","11-Mar","Mar-1","1-sep","Oct-3","Mar1","Oct-ps1","29-Dec","29sep","Sep-29") #test var
# paste(1:31,collapse = "|")
# substr(month.name,1,3) %>%unique()  %>% paste(collapse = "|")
#####sep-1,1-sep ....septin1
septin_maybe1= grep("^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30|31)-SEP$",ignore.case = T,x = gene_names) 
cat("septin_gene_maybe1:",gene_names[septin_maybe1],"\n")
septin_maybe2=grep("^SEP-(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30|31)$",ignore.case = T,x = gene_names) 
cat("septin_gene_maybe2:",gene_names[septin_maybe2],"\n")
if(length(septin_maybe1) | length(septin_maybe2)){
  print("septin gene names need to be transformed by hand maybe",quote=F)
  print(c(septin_maybe1,gene_names[septin_maybe1]))
  print(c(septin_maybe2,gene_names[septin_maybe2]))
} else{
  print("septin gene names detected well",quote=F)
}
#####sep-1,1-sep ....septin1 end
if("+" %in% special_char){
  print("E+ gene names need to be transformed by hand",quote = F)  
  temp = grep("+",gene_names,fixed = T)
  print(c(temp,gene_names[temp]))
} else{print("no gene names with E+,that's good",quote = F)}
##find and transform back
month_D_M = grep("^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30|31)-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$",ignore.case = T,x = gene_names)

month_M_D =grep("^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30|31)$",ignore.case = T,x = gene_names)

if(length(month_D_M) | length(month_M_D)){
  print("date format gene names detected, need to be transformed.",quote=F)
  print(c("positions:",month_D_M,month_M_D),quote = F)
  print(c("date format gene names:",gene_names[c(month_D_M,month_M_D)]),quote = F)
  print("gene names transforming...",quote = F)
  for (i in month_D_M) {
    temp= strsplit(gene_names[i],split = "-") %>% unlist()
    gene_names[i]= paste(temp[2],temp[1],sep = "")
  }
  for (i in month_M_D) {
    temp= strsplit(gene_names[i],split = "-") %>% unlist()
    gene_names[i]= paste(temp[1],temp[2],sep = "")
  }
  print(c("transformed gene names:",gene_names[c(month_D_M,month_M_D)]),quote = F)
  print("gene names transformed complete",quote=F)
} else{
  print("no date format gene names detected",quote=F)
}
#---
print("gene names detecting and transforming script complete.",quote= F)
#===============transform date gene names back end
#---end & push result to inputdata for other scripts
gene_names -> rownames(inputdata)
saveRDS(inputdata,"1output_check_gene_names.rds")