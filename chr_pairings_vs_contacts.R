library(tidyverse)
library(ggplot2)
library(pheatmap)
library(extrafont)
library(hrbrthemes)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(stringr)

options(hrbrthemes.loadfonts = TRUE)
#read files and clean the format
eorc <- read.csv('analysis_data_files/eosrc_chromosome_pairings.csv', sep = '\t')
eorc <- eorc %>%  group_by(cohort) %>% mutate(cohort_size=length(unique(samp_id))) 

crc <- read.csv('analysis_data_files/crc_chromosome_pairings.csv', sep = '\t')
crc <- crc %>%  group_by(cohort) %>% mutate(cohort_size=length(unique(samp_id))) 

tcga <- read.csv('analysis_data_files/pancan_chromosome_pairings.csv', sep = '\t')
tcga <- tcga %>% group_by(cohort, prog) %>% mutate(cohort_size=length(unique(samp_id)))
tcga <-  tcga %>% ungroup() %>% select(-prog)
tcga <- tcga[, colnames(eorc)]

# generate pooled data
all_data <- tcga %>% filter(cohort %in% c('TCGA-BRCA','TCGA-PAAD','TCGA-PRAD','TCGA-COAD'))
#all_data <- tcga %>% filter(cohort %in% c('TCGA-COAD','TCGA-READ','TCGA-BRCA'))
#all_data <- rbind(eorc, crc)
all_data$chr_pair <- str_replace(all_data$chr_pair,'chr','')
all_data$chr_pair <- str_replace(all_data$chr_pair,':chr',':')

length(unique(all_data$chr_pair))  # count all possible pairings observed
tbl <- all_data %>% group_by(cohort, chr_pair) %>% mutate(chr_pair_counts=n()) # get counts of all pairings
tbl <- tbl %>% ungroup()
tbl <- tbl %>% group_by(chr_pair) %>% mutate(chr_pair_mean = mean(chr_pair_counts), chr_pair_std = sd(chr_pair_counts))
#tbl <- tbl %>% mutate(chr_pair_mean = mean(chr_pair_counts), chr_pair_std = sd(chr_pair_counts))
tbl <- tbl %>% ungroup()
tbl <- tbl %>% group_by(chr_pair) %>% mutate(chr_pair_zscore =  (chr_pair_counts - chr_pair_mean)/chr_pair_std)
#tbl <- tbl  %>% mutate(chr_pair_zscore =  (chr_pair_counts - chr_pair_mean)/chr_pair_std)

tbl.short <- tbl %>% select(cohort, chr_pair, chr_pair_counts, chr_pair_mean, chr_pair_zscore) 
tbl.short <- tbl.short[!duplicated(tbl.short),] #temporary
#tbl.short.filt <- tbl.short %>% filter(chr_pair_zscore > 0.5 | chr_pair_zscore < -0.5)

#mat <- data.frame(tbl.short.filt %>% tidyr::pivot_wider(names_from = 'cohort', values_from = 'chr_pair_mean'))
mat <- data.frame(tbl.short.filt %>% tidyr::pivot_wider(names_from = 'cohort', values_from = 'chr_pair_zscore')) #tcga specific

mat.na <- mat %>% filter(if_any(everything(), is.na))
mat.nona <- mat %>% filter(!if_any(everything(), is.na))

rownames(mat) <- mat[,1] #tcga specific
mat[,1] <- NULL #tcga specific
mat <- na.omit(mat) #tcga specific

mat.nona[,1] <- NULL
rownames(mat.nona) <- mat.nona[,1]
mat.nona[,1] <- NULL
rownames(mat.na) <- mat.na[,1]
mat.na[,1] <- NULL

f <- read.csv('analysis_data_files/all_chr_pair_contactcounts_final.txt', sep = '\t', header = F)

#select pairs based on top contributing chr pairs shown in figure in EOSRC and CRC
colnames(f) <- c('pair','hct116','mcf7','pc3','panc1')
f[c('chr1','chr2')] <- str_split_fixed(f$pair,':',2)
f.melt <- f %>% gather(key = 'cohort', value = 'contacts', -pair, -chr1, -chr2)
colnames(f.melt) <- c('chr_pair','chr1','chr2','cohort','contacts')

brca <- tbl.short %>% filter(cohort == 'TCGA-BRCA')
mcf7 <- f.melt %>% filter(cohort == 'mcf7')
brca.set <- left_join(brca, mcf7, by='chr_pair')

prad <- tbl.short %>% filter(cohort == 'TCGA-PRAD')
pc3 <- f.melt %>% filter(cohort == 'pc3')
prad.set <- left_join(prad, pc3, by='chr_pair')

paad <- tbl.short %>% filter(cohort == 'TCGA-PAAD')
panc1 <- f.melt %>% filter(cohort == 'panc1')
paad.set <- left_join(paad, panc1, by='chr_pair')

coad <- tbl.short %>% filter(cohort == 'TCGA-COAD')
hct116 <- f.melt %>% filter(cohort == 'hct116')
coad.set <- left_join(coad, hct116, by='chr_pair')

all.set <- rbind(coad.set, brca.set, paad.set, prad.set)


p <- all.set %>%  ggplot(aes(x=contacts, y=chr_pair_counts)) +
  geom_point() +
  facet_wrap(. ~ cohort.x, scales = 'free') +
  stat_smooth(method = 'lm', se = F, linetype='dashed', aes(group=cohort.x)) +
  theme(panel.background = element_blank(), axis.line = element_line(),
      axis.text = element_text(colour = 'black')) +
  geom_text(label=ifelse((all.set$chr_pair_counts > 200 & all.set$contacts > 25000), all.set$chr_pair, ''))

save_pheatmap_png <- function(x, filename, width=8, height=8, units='in', res = 300) {
  setEPS()
  postscript(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p, "analysis_data_files/chr_contacts_vs_pair_counts_all_sets.eps")
ggsave(p, filename = 'analysis_data_files/chr_contacts_vs_pair_counts_all_sets.eps', device = cairo_ps, width = 8, height = 8)


  
