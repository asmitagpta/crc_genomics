library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
library(tibble)
library(export)
library(AnnotationHub)

f5 <- 'analysis_data_files/compartment_files/crc_5p_genes_overlap_ab_comp_2.txt'
f3 <- 'analysis_data_files/compartment_files/crc_3p_genes_overlap_ab_comp_2.txt'

dname <- str_split((str_split(f5,'/')[[1]][3]), '_')[[1]][1]

process_ab_overlap_file <- function(f5, f3){
  
  f5 <- read.csv(f5, sep = '\t', header = F)
  f3 <- read.csv(f3, sep = '\t', header = F)
  
  colnames(f5) <- c('pair','sign')
  colnames(f3) <- c('pair','sign')
  
  f5[, c('up','dw')] <- str_split_fixed(f5$pair, ':', 2)
  f3[, c('up','dw')] <- str_split_fixed(f3$pair, ':', 2)
  
  f5.s <- f5[!duplicated(f5[c('sign', 'up')]),]
  f3.s <- f3[!duplicated(f3[c('sign', 'dw')]),]
  
  f5.s <- f5.s %>% dplyr::select(sign, up) %>% dplyr::filter(up != '')
  f3.s <- f3.s %>% dplyr::select(sign, dw) %>% dplyr::filter(dw != '')
  
  f5.s$dir <- 'up'
  f3.s$dir <- 'dw'
  
  colnames(f5.s) <- c('sign','gene','dir')
  colnames(f3.s) <- c('sign','gene','dir')
  
  all_set <- rbind(f5.s, f3.s)
  all_set <- all_set[!duplicated(all_set[c('sign', 'gene')]),]
  
  all_set$dataset <- dname
  return(all_set)  
  
}
eorc <- process_ab_overlap_file(f5, f3)
eorc.s <- eorc %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% filter(n < 2)
eorc.m <- eorc %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% filter(n > 1)
crc <- process_ab_overlap_file(f5, f3)
crc.s <- crc %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% filter(n < 2)
crc.m <- crc %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% filter(n > 1)
all.set.s <- rbind(eorc.s, crc.s) # only for cell line comparison

#paad <- process_ab_overlap_file(f5, f3)
#paad.s <- paad %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% dplyr::filter(n < 2)
#paad.m <- paad %>% group_by(gene) %>% mutate(n=n()) %>% arrange(desc(n)) %>% dplyr::filter(n > 1)
#crc <- process_ab_overlap_file(f5, f3)

#all.set.s <- paad.s

# Perform sample wide analysis for RNA-Seq expression values
exp.mat <- read.csv('EOSRC_expression_data/coad_rlog_normalization_on_starcounts.txt', sep = '\t')

# pancancer RNA matrix specific part; not required for EOSRC and CRC
rowclean <- function(x){rownames(x) <- gsub("\\.[0-9]*", "", rownames(x), fixed = FALSE, perl = FALSE); x }
exp.mat <- rowclean(exp.mat)
exp.mat$gene <- rownames(exp.mat)

hub <- AnnotationHub()
query(hub, "EnsDb.Hsapiens.v104")
edb <- hub[["AH95744"]]
gene_ids <- AnnotationDbi::select(edb, keys=keys(edb), columns=c("GENENAME"))
colnames(gene_ids) <- c('gene','gene_name')
exp.mat <- left_join(exp.mat, gene_ids, first, by='gene')
colnames(exp.mat)[180] <- 'gene'
colnames(exp.mat)[179] <- 'gene_id'

# general section, can be used for any dataset; tcaa or non tcga
exp.mat.s <- left_join(exp.mat, paad.s, first, by='gene')
#tmp <- exp.mat.s[,c(-1107:-1112)]
exp.mat.s$med <- apply(exp.mat.s[,c(-179:-184)], 1, median, na.rm=T)
exp.mat.s <- exp.mat.s %>% mutate(col=ifelse(sign < -0.010, 'red3', 
                                               ifelse(sign >= -0.010 & sign < 0.010, 'yellow3', 'royalblue1')))

# perform cell line based expression analysis and A/B correlation
exp.cl <- read.csv('Expression_Public_22Q4_subsetted.csv')
exp.cl <- as.data.frame(t(exp.cl[, c(-1,-3:-6)]))
colnames(exp.cl) <- exp.cl[1,]
exp.cl <- rownames_to_column(exp.cl, 'gene')
exp.cl <- exp.cl[-1,]
exp.cl[,2:5] <- sapply(exp.cl[,2:5], as.double)
exp.cl <- exp.cl[rowSums(sapply(exp.cl[,2:5], as.double)) > 0,]

exp.mat.s2 <- left_join(exp.cl, all.set.s, by='gene')
exp.mat.s2 <- exp.mat.s2 %>% mutate(col=ifelse(sign < -0.010, 'red3', 
                                                 ifelse(sign >= -0.010 & sign < 0.010, 'yellow3', 'royalblue1')))
#eq <- lm(exp.mat.s2$sign ~ exp.mat.s2$HCT116)
#labl <- paste0('r2=', format(summary(eq)$r.squared, digits = 3))

cairo_ps(file = "analysis_data_files/pca_vs_gene_exp_eosrc_crc.eps", onefile = FALSE, fallback_resolution = 600, width = 10, height = 6)
p <- exp.mat.s2 %>% 
  dplyr::filter(dataset != '') %>% 
  ggplot(aes(x=sign, y=PANC1)) + #cell line based analysis
  #geom_point(aes(fill=col, color=col, shape=factor(dataset, levels = c('eosrc','crc')),
  #               size=factor(dataset, levels = c('eosrc','crc'))),  alpha=0.4) +
  geom_point(aes(fill=col, color=col), shape=21, alpha=0.4) +
  facet_wrap(. ~ dataset) +
  theme(panel.border =  element_rect(fill=NA, colour = 'black'),
        panel.background = element_blank(), panel.grid = element_line(color = 'grey94')) +
  scale_fill_manual(values = c(`red3`='red3',`yellow3`='yellow3', `royalblue1`='royalblue1')) + 
  scale_color_manual(values = c(`red3`='red3',`yellow3`='yellow3', `royalblue1`='royalblue1')) + 
  
  #scale_fill_manual(values = 'tomato2') +
  scale_x_continuous(expand = c(0,0.00)) +
  scale_y_continuous(expand = c(0.0,0.0)) +
  scale_shape_manual(values = c(21, 3), name='') +
  scale_size_manual(values = c(2.5,3.5), name='') +
  labs(x='Principal Eigen vector sign', y='TPM(log2)') +
  stat_smooth(method = 'lm', se = F, linetype='dashed', aes(group=dataset)) 
  
p
#p +   annotate(geom = "rect", xmin = min(exp.mat.s2$sign, na.rm = T), xmax = -0.010, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='red', alpha=0.2) + 
 #     annotate(geom = "rect", xmin = -0.010, xmax = 0.010, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='yellow', alpha=0.2) + 
  #    annotate(geom = "rect", xmin = 0.010, xmax = 0.056, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='blue', alpha=0.2) +
   #   stat_smooth(method = 'lm', color='black',linetype='dashed') +
    #  geom_text(x=-0.030, y=12, label=labl)

dev.off()

cairo_ps(file = "analysis_data_files/panc1-gene-exp-vs-pca-sample-paad.eps", onefile = FALSE, fallback_resolution = 600, width = 7, height = 6)
p <- exp.mat.s %>% 
  dplyr::filter(dataset != '') %>% 
  ggplot(aes(x=sign, y=med)) + #cell line based analysis
  geom_point(aes(fill=col, color=col), shape=21,  alpha=0.4) +
  #geom_point(aes(fill=dataset),shape=21, size=1.5, alpha=0.4) +
  theme(panel.border =  element_rect(fill=NA, colour = 'black'),
        panel.background = element_blank(), panel.grid = element_line(color = 'grey94')) +
  scale_fill_manual(values = c(`red3`='red3',`yellow3`='yellow3', `royalblue1`='royalblue1')) + 
  scale_color_manual(values = c(`red3`='red3',`yellow3`='yellow3', `royalblue1`='royalblue1')) + 
  labs(x='Principal Eigen vector sign', y='TPM(log2)') +
  stat_smooth(method = 'lm', se = F, linetype='dashed', aes(group=dataset)) 

p
#p +   annotate(geom = "rect", xmin = min(exp.mat.s2$sign, na.rm = T), xmax = -0.010, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='red', alpha=0.2) + 
#     annotate(geom = "rect", xmin = -0.010, xmax = 0.010, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='yellow', alpha=0.2) + 
#    annotate(geom = "rect", xmin = 0.010, xmax = 0.056, ymin = 0 , ymax = max(as.double(exp.mat.s2$MCF7), na.rm = T), fill='blue', alpha=0.2) +
#   stat_smooth(method = 'lm', color='black',linetype='dashed') +
#  geom_text(x=-0.030, y=12, label=labl)

dev.off()

