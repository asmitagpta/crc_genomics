# The code was used to download and process gene expression and other datasets from TCGA portal 

library(TCGAbiolinks)
library(TCGAbiolinksGUI.data)
library(ExperimentHub)
library(SummarizedExperiment)
library(stringr)
library(pheatmap)
library(CRCAssigner)
library(CMSclassifier)

coad <- read.csv('coad_tcga_barcodes.txt', header = F)
#read <- read.csv()
eh=ExperimentHub()
mm450 = query(eh, 'HM450.address')
mm450[['EH5966']]

query <- GDCquery(project = 'TCGA-COAD',
                  access = 'open', 
                  legacy = F,
                  barcode = coad$V1,
                  sample.type = "Primary Tumor",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")

query.meth <- GDCquery(project = 'TCGA-COAD',
                  access = 'open', 
                  legacy = F,
                  barcode = coad$V1,
                  sample.type = "Primary Tumor",
                  data.category = "DNA Methylation",
                  data.type = "Masked Intensities",
                  platform = "Illumina Human Methylation 450")                


summ <- getResults(query.rna)
GDCdownload(query = query, method = "api")
coad_meth <- GDCprepare(query.meth)

TCGAanalyze_Preprocessing(coad_rna)
###### analyze rna-seq data #####

ens_ids <- rownames(coad_rna_matrix) 
ens_ids_nover <- str_split(ens_ids,'\\.')
genes5p <- read.csv('coad_5p_gene_ens_ids.txt', header = F, stringsAsFactors = F)
set.seed(21)
gene_id <- character(length(ens_ids_nover))
for(x in 1:length(ens_ids_nover)){
  gene_id[x] <- ens_ids_nover[[x]][1]
 } 

rownames(coad_rna_matrix) <- gene_id
coad_rna_matrix_5p_fusion <-  coad_rna_matrix[rownames(coad_rna_matrix) %in% genes5p$V1,]
coad_rna_matrix_no_fusion_mat <- coad_rna_matrix[!rownames(coad_rna_matrix) %in% genes5p$V1,]
samp_rows <- sample(rownames(coad_rna_matrix_no_fusion_mat), dim(coad_rna_matrix_5p_fusion)[1])
coad_rna_matrix_no_fusion_test <- coad_rna_matrix_no_fusion_mat[samp_rows,]

# remove genes with no expression across samples
tmp=coad_rna_matrix_5p_fusion[apply(coad_rna_matrix_5p_fusion[,-1],1,function(x) length(x[x>=5])>5),]
tmp2=coad_rna_matrix_no_fusion_test[apply(coad_rna_matrix_no_fusion_test[,-1],1,function(x) length(x[x>=5])>5),]
pheatmap(tmp)

#CRC subtype
write.table(coad_rna_matrix_filter, file = 'coad_rna_matrix.txt', sep = '\t', row.names = TRUE, col.names = TRUE)
CRCA(direc = '.', file = 'coad_rna_matrix.txt', PAM = 'PAM786')

