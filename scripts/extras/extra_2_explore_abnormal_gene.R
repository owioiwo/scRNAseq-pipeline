## run this script behind 2_data_initial.R
##===extra_of_step2 : explore abnormal genes
gene_ncount = rowSums(seurat_obj@assays[["RNA"]]@counts)
gene_nfeature = rowSums(seurat_obj@assays[["RNA"]]@counts >0)
boxplot(gene_ncount)
boxplot(gene_nfeature)
ggplot(mapping = aes(gene_ncount,gene_nfeature)) +  geom_point()
cat("each gene's nCount & nFeature plotsss are drawn.\n")
quantile(gene_ncount)
quantile(gene_nfeature)
gene_medcount_min = median(gene_ncount) -3*sd(gene_ncount)
gene_medcount_max = median(gene_ncount) +3*sd(gene_ncount)
cat("gene_ncount med~3sd:",gene_medcount_min,gene_medcount_max,"\n")
gene_medfeature_min = median(gene_nfeature)- 3*sd(gene_nfeature)
gene_medfeature_max = median(gene_nfeature)+ 3*sd(gene_nfeature)
cat("gene_nfeature med~3sd:",gene_medfeature_min,gene_medfeature_max,"\n")
# find top(n) abnormal gene , see 5 genes' counts special in plot  
n = 5
abnormal_gene_bycount = head(sort(gene_ncount,decreasing = T),n)
abnormal_gene_bycount 
#####   find 4 genes' counts > 10000 :mt-Rnr2 Gm42418 Gm26917 mt-Rnr1
###   mt-Rnr1,mt-Rnr2:12S + 16S rRNA, mitochondrial,rRNA 
##  Gm42418  = Rn18s-rs5,pseudo gene
#  Gm26917 ,ncRNA
VlnPlot(seurat_obj,features = names(abnormal_gene_bycount))
##===extra : explore abnormal genes end