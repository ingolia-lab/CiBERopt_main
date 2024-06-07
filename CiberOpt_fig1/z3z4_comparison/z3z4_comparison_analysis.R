library(mpra)
library(ggplot2)
library(ggrastr)

setwd('~/CiBERopt_main//CiBERopt_counts/z3z4_counts/')

f = list.files(pattern="*.txt") 

tables=  lapply(f, function(x){
  f = read.table(x, header=F,sep = '')
  colnames(f) = c(x,'bcs')
  #need to load from position 1
  f = f[1:length(f[,1]),]
  return(f)
})

counts= Reduce(function(...) merge(..., all = TRUE, by = "bcs"), tables)
counts[is.na(counts)] <- 0
colnames(counts)[2:ncol(counts)] = sapply(strsplit(colnames(counts)[2:ncol(counts)],'_'),'[[',2)

setwd('~/CiBERopt_main/CiBERopt_counts/')
if (!file.exists("SGD_features.tab")) {
  sgdfile <- 'https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab'
  sgd <- download.file(sgdfile, destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

setwd('~/CiBERopt_main/CiberOpt_fig1/z3z4_comparison/')

#remove unassigned barcodes 
counts = counts[2:nrow(counts),]
#add the counts together from the separate sequencing run
counts = cbind(counts[,1],as.data.frame(colsum(counts[,2:ncol(counts)],group=colnames(counts)[2:ncol(counts)])))
colnames(counts)[1] = 'bcs'
counts$guide = paste0(sapply(strsplit(counts$bcs, '_'),'[[',1),'_',sapply(strsplit(counts$bcs, '_'),'[[',2))
counts$barcode = sapply(strsplit(counts$bcs, '_'),'[[',3)

#use pre samples as cutoff
filtered = counts[counts$`dna-pre-1`> 5 & 
                    counts$`dna-pre-2` > 5 &
                    counts$`z-cfp-pre-1` > 5 &
                    counts$`z-cfp-pre-2` > 5 &
                    counts$`z-yfp-pre-2` > 5 &
                    counts$`z-yfp-pre-2` > 5,]

#analyze the RNA to RNA counts, pre-vs-post comparison
mpra_RNA = MPRASet(DNA=filtered[,grepl('cfp',colnames(filtered))], 
                   RNA=filtered[,grepl('yfp',colnames(filtered))], 
                   barcode=NULL, eid=filtered$guide,eseq=filtered$barcode)

designRNA = data.frame(intcpt = 1,
                       post = grepl("post", colnames(mpra_RNA)),
                       turb = grepl('-1', colnames(mpra_RNA)))

fit_RNA <- mpralm(object = mpra_RNA, design = designRNA,
                  aggregate = "sum", normalize = TRUE,
                  model_type = "indep_groups", plot = TRUE)

rna_tab <- topTable(fit_RNA, coef = 2, number = Inf)
rna_tab$sgid = sapply(strsplit(rownames(rna_tab),'_'),'[[',1)
rna_tab$gene = sgd[match(rna_tab$sgid, sgd$name),'gene']


#analyze the RNA to DNA counts, pre-vs-post comparison
mpra_DNA = MPRASet(DNA=filtered[,grepl('dna',colnames(filtered))], 
                   RNA=filtered[,grepl('yfp',colnames(filtered))], 
                   barcode=NULL, eid=filtered$guide,eseq=filtered$barcode)

designDNA = data.frame(intcpt = 1,
                       post = grepl("post", colnames(mpra_DNA)),
                       turb = grepl('-1', colnames(mpra_DNA)))

fit_DNA <- mpralm(object = mpra_DNA, design = designRNA,
                  aggregate = "sum", normalize = TRUE,
                  model_type = "indep_groups", plot = TRUE)

dna_tab <- topTable(fit_DNA, coef = 2, number = Inf)
dna_tab$sgid = sapply(strsplit(rownames(dna_tab),'_'),'[[',1)
dna_tab$gene = sgd[match(dna_tab$sgid, sgd$name),'gene']



#plotting 
rna_tab$sigdiff = ifelse((rna_tab$logFC < -1 | rna_tab$logFC > 1) & rna_tab$adj.P.Val < 0.01, T, F)
dna_tab$sigdiff = ifelse((dna_tab$logFC < -1 | dna_tab$logFC > 1) & dna_tab$adj.P.Val < 0.01, T, F)

#if p-vals > 1e-50, max at 1e-50
dna_tab$adj.P.Val = ifelse(dna_tab$adj.P.Val < 1e-50, 1e-50, dna_tab$adj.P.Val)


rna_plot = ggplot() + 
  geom_point_rast(data= subset(rna_tab,!(rna_tab$sigdiff)), aes(x=logFC, y=-log10(adj.P.Val)), size = 2,raster.dpi = 600, col='#bdbdbd') + 
  geom_point_rast(data= subset(rna_tab,rna_tab$sigdiff), aes(x=logFC, y=-log10(adj.P.Val)), size = 2,raster.dpi = 600, col='#f03b20') +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01),limits = c(0,50)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=0.75,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none",
        panel.border = element_blank())
  

dna_plot = ggplot() + 
  geom_point_rast(data= subset(dna_tab,!(dna_tab$sigdiff)), aes(x=logFC, y=-log10(adj.P.Val)), size = 2,raster.dpi = 600, col='#bdbdbd') + 
  geom_point_rast(data= subset(dna_tab,dna_tab$sigdiff), aes(x=logFC, y=-log10(adj.P.Val)), size = 2,raster.dpi = 600, col='#f03b20') +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed', size=1) +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed',size=1) + 
  theme_classic() + xlab('Log2 Fold Change') + ylab('-log10(Padj)') + 
  scale_y_continuous(expand=c(0.01,0.01),limits = c(0,50)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio=0.75,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none",
        panel.border = element_blank())

if(!(file.exists('rna_ciber.pdf'))){
  ggsave('rna_ciber.pdf',rna_plot, height = 6, width = 6,units='in')
}
if(!(file.exists('dna_ciber.pdf'))){
  ggsave('dna_ciber.pdf',dna_plot, height = 6, width = 6, units = 'in')
}


####For supp, look at GO terms that change in the dnaplot
#check GO terms for DNA
sig_up_dna = unique(dna_tab[dna_tab$logFC > 1 & dna_tab$adj.P.Val < 0.01, 'gene'])
sig_down_dna = unique(dna_tab[dna_tab$logFC < -1 & dna_tab$adj.P.Val < 0.01, 'gene'])
all_dna = unique(dna_tab$gene)

setwd('~/CiBERopt_main/CiBERsupp/ciberOpt_supp2/DNAgoTerms/')

#for analysis, used fisher exact test and bonferroni correction in pantherDB. download results. 
write.table(x = sig_up_dna,'sig_up_dna_ids.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = sig_down_dna,'sig_down_dna_ids.txt',sep = '',quote = F,row.names = F, col.names = F)
write.table(x = all_dna,'all_dna_ids.txt',sep = '',quote = F,row.names = F, col.names = F)

  

