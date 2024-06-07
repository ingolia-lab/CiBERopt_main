library(mpra)
library(ggplot2)
library(ggrastr)
library(ggrepel)

setwd("~/CiBERopt_main/CiBERopt_counts/nmd_counts/")

files = list.files(pattern="*.txt") 

nmd_tables=  lapply(files, function(x){
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

setwd("~/CiBERopt_main/CiberOpt_fig4/nmd_volcanoPlots/")
nmd_counts= Reduce(function(...) merge(..., all = TRUE, by = "bcs"), nmd_tables)
nmd_counts[is.na(nmd_counts)] <- 0
colnames(nmd_counts)[2:ncol(nmd_counts)] = sapply(strsplit(colnames(nmd_counts)[2:ncol(nmd_counts)],'_'),'[[',2)

#combine reads from multiple sequencing runs, didn't reseq all libs
counts = data.frame(nmd_counts$bcs)
for(col in colnames(nmd_counts)[2:ncol(nmd_counts)]){
  if(sum(grepl(col,colnames(nmd_counts))) > 1){
    counts[,col] = rowSums(nmd_counts[,grepl(col,colnames(nmd_counts))])
  }else{
    counts[,col] = nmd_counts[,col]
  }
  colnames(counts)[1] = 'bcs'
}

counts = counts[counts$bcs != '*' &  counts$bcs != 'No_barcode',]
counts$guide = paste0(sapply(strsplit(counts$bcs,"_"),'[[',1),'_',sapply(strsplit(counts$bcs,"_"),'[[',2))
counts$barcode = sapply(strsplit(counts$bcs,"_"),'[[',3)

#apply filter of 5 reads in pre,
filtered = counts[counts$`NMD-Pre-1-Scarlet` > 5 & 
                    counts$`NMD-Pre-2-Scarlet` > 5 &
                    counts$`NMD-Pre-1-YFP` > 5 &
                    counts$`NMD-Pre-2-YFP` >5,]


####perform pre vs post analysis
mpra_nmd = MPRASet(DNA=filtered[,grepl('Scarlet',colnames(filtered)) & !(grepl('CHX',colnames(filtered)))], 
                   RNA=filtered[,grepl('YFP',colnames(filtered)) & !(grepl('CHX',colnames(filtered)))], 
                   barcode=filtered$barcode, eid=filtered$guide,eseq=NULL)

design_nmd <- data.frame(intcpt = 1,
                         post = grepl("Post", colnames(mpra_nmd)),
                         turb = grepl('-1', colnames(mpra_nmd))) 

fit_nmd <- mpralm(object = mpra_nmd, design = design_nmd,
                  aggregate = "sum", normalize = TRUE,
                  model_type = "indep_groups", plot = TRUE)

nmd_tab <- topTable(fit_nmd, coef = 2, number = Inf)
nmd_tab$sgid = sapply(strsplit(rownames(nmd_tab),'_'),'[[',1)
nmd_tab$gene = sgd[match(nmd_tab$sgid, sgd$name),'gene']
nmd_tab$gene = ifelse(!(nzchar(nmd_tab$gene)),nmd_tab$sgid,nmd_tab$gene)

#some pvals are very high, set max pval< 1e-50 to 1e-50 
nmd_tab$adj.P.Val =ifelse(nmd_tab$adj.P.Val<1e-50,1e-50,nmd_tab$adj.P.Val)

#get upregulated genes
nmd_tab$up = ifelse(nmd_tab$logFC > 1 & nmd_tab$adj.P.Val < 0.01, T,F)
#env11 annotation certianly targets adjacent Upf3, reassign name. Set Namy7 to Upf1, Nmd2 to upf2 for ease
nmd_tab$gene = ifelse(nmd_tab$gene == 'ENV11' & nmd_tab$up,'UPF3',nmd_tab$gene)
nmd_tab$gene = ifelse(nmd_tab$gene == 'NAM7' & nmd_tab$up,'UPF1',nmd_tab$gene)
nmd_tab$gene = ifelse(nmd_tab$gene == 'NMD2' & nmd_tab$up,'UPF2',nmd_tab$gene)

nmd_plot = ggplot(data=nmd_tab, aes(x=logFC,y=-log10(adj.P.Val),color=up)) +
  geom_point_rast(raster.dpi=100) +
  scale_color_manual(values=c("#bdbdbd","#f03b20")) + 
  geom_text_repel(data = subset(nmd_tab, nmd_tab$adj.P.Val < 0.01 & nmd_tab$up),aes(label=gene), max.overlaps=Inf) + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,50)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = 'none',
        panel.border = element_blank())


if(!(file.exists('nmd_prePost_volcano.pdf'))){
  ggsave('nmd_prePost_volcano.pdf',plot = nmd_plot, height = 5, width = 5, units = 'in')
}

#write GO terms to supp
setwd("~/CiBERopt_main/CiBERsupp/ciberOpt_supp4/nmd_goTerms/")
sigup = unique(nmd_tab[nmd_tab$logFC > 1 & nmd_tab$adj.P.Val < 0.01, 'gene'])
sigdown = unique(nmd_tab[nmd_tab$logFC < -1 & nmd_tab$adj.P.Val < 0.01, 'gene'])
all_genes = unique(nmd_tab$gene)

#used fisher exact test and bonferroni correction in pantherDB. download results to supp
write.table(x = sigup,'nmd_up.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = sigdown,'nmd_down.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = all_genes,'nmd_all.txt',sep = '',quote = F,row.names = F, col.names = F)



###compare CHX treatment to 'post', using post as baseline for supp
#expectation is NMD guides that  lose activity will have NEGATIVE activity score
#since it won't be further stabilized upon chx treatment 
#normal uninteresting gene will have stabilized reporter and manifest as an positive score'

setwd('~/CiBERopt_main/CiBERsupp/ciberOpt_supp4/nmd_chx_screen/')


mpra_chx = MPRASet(DNA=filtered[,grepl('Scarlet',colnames(filtered)) & !(grepl('Pre',colnames(filtered)))], 
                   RNA=filtered[,grepl('YFP',colnames(filtered)) & !(grepl('Pre',colnames(filtered)))], 
                   barcode=filtered$barcode, eid=filtered$guide,eseq=NULL)

design_chx <- data.frame(intcpt = 1,
                         post = grepl("Post", colnames(mpra_nmd)),
                         turb = grepl('-1', colnames(mpra_nmd))) 

fit_chx <- mpralm(object = mpra_chx, design = design_chx,
                  aggregate = "sum", normalize = TRUE,
                  model_type = "indep_groups", plot = TRUE)

chx_tab <- topTable(fit_chx, coef = 2, number = Inf)
chx_tab$sgid = sapply(strsplit(rownames(chx_tab),'_'),'[[',1)
chx_tab$gene = sgd[match(chx_tab$sgid, sgd$name),'gene']
chx_tab$gene = ifelse(!(nzchar(chx_tab$gene)),chx_tab$sgid,chx_tab$gene)

#compress low p val for plotting
chx_tab$adj.P.Val =ifelse(chx_tab$adj.P.Val<1e-50,1e-50,chx_tab$adj.P.Val)

#find guides that were active in the nmd screen
chx_tab$nmd_gene = nmd_tab$up[match(rownames(chx_tab),rownames(nmd_tab))]

chx_plot = ggplot() +
  geom_point_rast(data = subset(chx_tab,!(chx_tab$nmd_gene)), aes(x = logFC, y = -log10(adj.P.Val)),color = '#bdbdbd', raster.dpi = 100) +
  geom_point_rast(data = subset(chx_tab,chx_tab$nmd_gene), aes(x = logFC, y = -log10(adj.P.Val)),color = '#f03b20', raster.dpi = 100) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,50)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = 'none',
        panel.border = element_blank())

if(!(file.exists('nmd_CHX_volcano.pdf'))){
  ggsave('nmd_CHX_volcano.pdf',plot = chx_plot, height = 5, width = 5, units = 'in')
}

