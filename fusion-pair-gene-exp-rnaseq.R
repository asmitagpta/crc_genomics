library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(AnnotationHub)

cdfd_id <- read.csv('CDFD_transcriptome_data_IDs.csv')
cdfd_id <- cdfd_id %>% dplyr::filter(Datatype == 'Transcriptome')
cdfd_id$CDFD.File.IDs <- gsub('/ Tumor','', cdfd_id$CDFD.File.IDs)
cdfd_id$CDFD.File.IDs <- gsub('A[0-9]*/ |E[0-9]*/','', cdfd_id$CDFD.File.IDs)
colnames(cdfd_id) <- c('sample_ID','datatype','alt_id','real_id')

cdfd_id$sample_ID <- gsub('MDB_R[I]*_','', cdfd_id$sample_ID)

hub <- AnnotationHub()
query(hub, "EnsDb.Hsapiens.v104")
edb <- hub[["AH95744"]]
gene_ids <- AnnotationDbi::select(edb, keys=keys(edb), columns=c("GENENAME"))
colnames(gene_ids) <- c('gene_ids','gene_name')

exp_counts <- read.csv('EOSRC_expression_data/eosrc_rlog_normalization_on_fcounts.txt', sep = '\t')
#exp_counts <- read.csv('../../eosrc_project/rsem_exp_counts.txt', sep = '\t')
colnames(exp_counts)[1] <- 'gene_ids'
exp_counts <- left_join(exp_counts, gene_ids, first, by='gene_ids')

fnames <- readxl::read_xlsx('fusion-pair-list-rt.xlsx', sheet = 2)
fnames <- fnames %>% mutate(pair = paste(up_gene,'--',dw_gene)) %>% 
  dplyr::select(up_gene, dw_gene, pair) %>% 
  gather(key = 'key', value = 'value', -pair)

fnames.s <- fnames[order(fnames$value), ]
colnames(fnames.s) <- c('pair','key','gene_ids')

exp_counts_sub <- exp_counts[exp_counts$gene_ids %in% fnames.s$gene_ids, ]
#values <- exp_counts_sub$gene_id # was used when rsem_exp_counts was used

exp_counts_sub <- left_join(exp_counts_sub, fnames.s, by='gene_ids')
#exp_counts_sub_data <- gather(exp_counts_sub, key = 'sample_ID', value = 'counts', -gene_ids,
#                             -genename, -pair, -key)
exp_counts_sub_data <- gather(exp_counts_sub, key = 'sample_ID', value = 'counts',
                                                           -gene_ids, -pair, -key)
colnames(exp_counts_sub_data) <- c('gene_name','pair','gene_loc','sample_ID','counts')
exp_counts_sub_data <- left_join(exp_counts_sub_data, cdfd_id, by='sample_ID')

## sub data for rapsac meeting
#exp_counts.s <- exp_counts %>% filter(gene_id %in% c('ENSG00000104408','ENSG00000147655','ENSG00000144959','ENSG00000144962'))
#exp_counts.s$gene <- c('EIF3E','NCEH1','SPATA16','RSPO2')
#exp_counts.s.long <- exp_counts.s %>% gather(key='samp_id',value='values',-gene_id,-gene)

p1 <- exp_counts_sub_data %>% 
  #dplyr::filter(!gene_name %in% c('PTPRK','RSPO3')) %>% 
  ggplot(aes(x = real_id, y = counts, fill=gene_loc)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  #facet_wrap_paginate(~ pair + gene_name, scales = 'free', labeller = label_wrap_gen(multi_line = FALSE),
   #                  nrow = 6, ncol = 2, page = 1)
  facet_wrap_paginate(~ pair + factor(gene_loc, levels = c('up_gene','dw_gene')),
                      scales = 'free', labeller = label_wrap_gen(multi_line = FALSE),
                                      ncol = 2, page = 4) +
                     
  #facet_wrap( ~ pair + factor(gene_loc, levels = c('up_gene','dw_gene')), 
   #           ncol = 2, scales = 'free_y', labeller = label_wrap_gen(multi_line = F)) +
 # facet_grid(gene_name ~ pair,  scales = 'free_y') +
  
  theme(panel.background = element_blank(), axis.line = element_line(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #scale_x_discrete(labels = seq(1,37)) +
  scale_fill_manual(values = c('steelblue3','red3'), name='')
p1

ggsave(p1, filename = 'fusion-pair-gene-exp-rnaseq-ab-fusions.eps', width = 10, height = 55)
ggsave(p1, filename = 'fusion-pair-gene-exp-rnaseq.eps', device = cairo_ps, width = 10, height = 16, dpi = 300)
ggsave(p1, filename = 'fusion-pair-gene-exp-rnaseq-mod.eps', device = cairo_ps, width = 10, height = 15, dpi = 300)



