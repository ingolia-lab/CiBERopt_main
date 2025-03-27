library(ggplot2)
library(dplyr)
library(ggrastr)
library(fastmatch)
library(cowplot)

#load GOterms for plotting
GOpath = setwd('~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/previousCiBERcomparisons/gcn4goTerms/')
file_names = list.files(patter='*.txt')
list_tables <- lapply(file_names, read.table)

ids = c('Amino acid biosynthesis', 'Nuclear Pore' ,'Polymerase III subunit','Sumoylation',
        'Translational control','tRNA amino acyl charging','tRNA processing')


for(i in 1:length(list_tables)){
  list_tables[[i]]$id = ids[i]
}

goMaster = do.call(rbind,list_tables) 

goUTR = goMaster %>% filter(id != 'Nuclear Pore' & id != 'Sumoylation')
goCDS = goMaster %>% filter(id == 'Nuclear Pore' | id == 'Sumoylation')


utr = read.table("~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/previousCiBERcomparisons/UTRscreen_results.txt")
cds = read.table("~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/previousCiBERcomparisons/CDSscreen_results.txt")

#load the general transcription machinery that was aggregated from multiple GOterms
genTF = read.csv("generalTxn_handcur.csv")


#clip pvals and log2fc for adjustment, get yorf name, and match to GO term. 
utr = utr %>% mutate(clipped_pval = ifelse(adj.P.Val < 1e-10,1e-10,adj.P.Val),
                           clipped_log2fc = -1*ifelse(abs(logFC) >= 6, sign(logFC)*6, logFC), #-1 accounts for comparison differences in original paper (pre/post vs post/pre)
                           yorf = sapply(strsplit(Guide,'_'),'[[',1)) %>%
  mutate(goTerm = goUTR$id[fmatch(yorf,goUTR$name)]) %>% mutate(goTerm = ifelse(yorf %in% genTF$yorf, 'genTF', goTerm))
  

cds = cds %>% mutate(clipped_pval = ifelse(adj.P.Val < 1e-10,1e-10,adj.P.Val),
                     clipped_log2fc = -1*ifelse(abs(logFC) >= 4, sign(logFC)*4, logFC), #-1 accounts for comparison differences in original paper (pre/post vs post/pre)
                     yorf = sapply(strsplit(Guide,'_'),'[[',1)) %>%
  mutate(goTerm = goCDS$id[fmatch(yorf,goCDS$name)]) %>% mutate(goTerm = ifelse(yorf %in% genTF$yorf, 'genTF', goTerm))




utrPlot = ggplot() +
  geom_point_rast(data = subset(utr,is.na(goTerm)) , aes(x = clipped_log2fc, y = -log10(clipped_pval)),color = '#bdbdbd', raster.dpi = 600) +
  geom_point_rast(data = subset(utr, goTerm == 'genTF'), aes(x = clipped_log2fc, y = -log10(clipped_pval)), color = '#d95f02', raster.dpi = 600) + 
  geom_point_rast(data = subset(utr, !is.na(goTerm) & goTerm != 'genTF'), aes(x = clipped_log2fc, y = -log10(clipped_pval)),color = '#1b9e77',raster.dpi = 600) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01)) + xlim(-6,6) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_blank())


cdsPlot = ggplot() +
  geom_point_rast(data = subset(cds,is.na(goTerm)) , aes(x = clipped_log2fc, y = -log10(clipped_pval)),color = '#bdbdbd', raster.dpi = 600) +
  geom_point_rast(data = subset(cds, goTerm == 'genTF'), aes(x = clipped_log2fc, y = -log10(clipped_pval)), color = '#d95f02',raster.dpi = 600) +  #scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#E78AC3','#A6D854'))+
  geom_point_rast(data = subset(cds, !is.na(goTerm) & goTerm != 'genTF'), aes(x = clipped_log2fc, y = -log10(clipped_pval)), color = '#1b9e77',raster.dpi = 600) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01)) + xlim(-4,4) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_blank())

figplot = plot_grid(utrPlot, cdsPlot,ncol=2)

if(!(file.exists('~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/previousCiBERcomparisons/previousScreens.pdf'))){
  ggsave('~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/previousCiBERcomparisons/previousScreens.pdf',plot = figplot, height = 5, width = 12, units = 'in')
}






 
