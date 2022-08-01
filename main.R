setwd("~/scRNAseq-pipeline/results/")
source("../scripts/0_load_functions.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(sctransform)
library(glmGamPoi)
library(SingleR)
library(celldex)
##++++++++++++++++++++++++++++++++++++++++++++
##### step A : load data, support matrix.txt and seurat.rds
# data should be in "rawdata" directory
# data names should be in "0_inputdata_list.txt",one line = one data name
# only support gene matrix or seurat obj files
# also users can inputdata through RStudio console manually.

inputdata = import_data()

##### step B : single gene exp matrix processing (check and filter)
### warning: belows can only run one line at one time =======
source("../scripts/1_check_gene_names.R")
source("../scripts/2_data_initial.R")
source("../scripts/3_deepclustering.R")
#### repeat step A & B to process multi datasets one by one

##### step optional : integrate datas, using seurat integrate method
# data should be in "results" directory
# data names should be in "4_integrate_list.txt",one line = one data name
source("../scripts/4_integrate_data.R")

##### step C : do gene de analysis and cell annotation
# data should be in "results" directory
# data name should be in "5_gene_cell_analysis.txt",one line = one data name
source("../scripts/5_gene_cell_analysis.R")

### end


