library(mpra)
library(ggplot2)
library(ggrastr)

setwd("~/CiBERopt_main/CiBERopt_counts/cl1_counts/")
f = list.files(pattern="*.txt") 

tables=  lapply(f, function(x){
  f = read.table(x, header=F,sep = '')
  colnames(f) = c(x,'bcs')
  #need to load from position 1
  f = f[1:length(f[,1]),]
  return(f)
})

setwd("~/CiBERopt_main/CiBERopt_counts/")
if (!file.exists("SGD_features.tab")) {
  sgdfile <- 'https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab'
  sgd <- download.file(sgdfile, destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))


setwd("~/CiBERopt_main/CiberOpt_fig3/CL1_degronAnalysis/")


counts= Reduce(function(...) merge(..., all = TRUE, by = "bcs"), tables)
counts[is.na(counts)] <- 0
colnames(counts)[2:ncol(counts)] = sapply(strsplit(colnames(counts)[2:ncol(counts)],'_'),'[[',2)

#remove unassign, add guide/barcode
counts = counts[counts$bcs != '*' &  counts$bcs != 'No_barcode',]
counts$guide = paste0(sapply(strsplit(counts$bcs,"_"),'[[',1),'_',sapply(strsplit(counts$bcs,"_"),'[[',2))
counts$barcode = sapply(strsplit(counts$bcs,"_"),'[[',3)

#apply filter of 5 reads in pre
filtered = counts[counts$`CL1-Pre-1-CFP` > 5 & 
                    counts$`CL1-Pre-1-YFP` > 5 &
                    counts$`CL1-Pre-2-CFP` > 5 &
                    counts$`CL1-Pre-2-YFP` > 5,]


mpra_cl1 = MPRASet(DNA=filtered[,grepl('CFP',colnames(filtered))], 
                   RNA=filtered[,grepl('YFP',colnames(filtered))], 
                   barcode=filtered$barcode, eid=filtered$guide,eseq=NULL)

design_cl1 <- data.frame(intcpt = 1,
                         post = grepl("Post", colnames(mpra_cl1)),
                         turb = grepl('-1', colnames(mpra_cl1))) 

fit_cl1 <- mpralm(object = mpra_cl1, design = design_cl1,
                  aggregate = "sum", normalize = TRUE,
                  model_type = "indep_groups", plot = TRUE)

cl1_tab <- topTable(fit_cl1, coef = 2, number = Inf)
cl1_tab$sgid = sapply(strsplit(rownames(cl1_tab),'_'),'[[',1)
cl1_tab$gene = sgd[match(cl1_tab$sgid, sgd$name),'gene']
cl1_tab$gene = ifelse(!(nzchar(cl1_tab$gene)),cl1_tab$sgid,cl1_tab$gene)

#load functional annotations. Notes, guides assigned to multiple genes were called for the pathway
#they were next to (i.e ALG13 is next to RPT6). #completed for manual
ids = read.csv('annotationsCL1.csv')

cl1_tab$label = ids$annot[match(cl1_tab$gene,ids$id)]

#check number of guides that are invovled in doa10, proteosome, or cdc48
up = cl1_tab[cl1_tab$logFC > 2 & cl1_tab$adj.P.Val < 0.01,]
perc_guides = sum(!(is.na(up$label)))/nrow(up)



cl1_plot = ggplot() +
  geom_point_rast(data = subset(cl1_tab, is.na(label)), aes(x = logFC, y = -log10(adj.P.Val)),color = '#bdbdbd', raster.dpi = 600) +
  geom_point_rast(data = subset(cl1_tab, !is.na(label)), aes(x = logFC, y = -log10(adj.P.Val), color = label),raster.dpi = 600) +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  geom_vline(xintercept = c(-2, 2), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_blank())


if(!(file.exists('cl1_volcano.pdf'))){
  ggsave('cl1_volcano.pdf',plot = cl1_plot, height = 5, width = 5, units = 'in')
}


#for GO for supp
sigup = unique(cl1_tab[cl1_tab$logFC > 2 & cl1_tab$adj.P.Val < 0.01, 'gene'])
sigdown = unique(cl1_tab[cl1_tab$logFC < -2 & cl1_tab$adj.P.Val < 0.01, 'gene'])
all_genes = unique(cl1_tab$gene)

setwd('~/CiBERopt_main/CiBERsupp/ciberOpt_supp3/CL1_goTerms/')

#for analysis, used fisher exact test and bonferroni correction in pantherDB. download results. 
write.table(x = sigup,'cl1_up.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = sigdown,'cl1_down.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = all_genes,'cl1_all.txt',sep = '',quote = F,row.names = F, col.names = F)



