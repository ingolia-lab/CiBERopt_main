library(ggplot2)

setwd("~/CiBERopt_main/CiBERsupp/ciberOpt_supp4/nmd_goTerms/")

go = read.csv("nmd_GOterms_up.csv")
colnames(go) = go[7,]
go = go[8:nrow(go),]

go$Padj = as.numeric(go$`nmd_up.txt (P-value)`)
go$found = as.numeric(go$`nmd_up.txt (9)`)
go$fold_enrichment = go$`nmd_up.txt (fold Enrichment)`
#some fold enrichment >100, cut off at 100
go$fold_enrichment[1] = 100
go$fold_enrichment= as.numeric(go$fold_enrichment)
go$GO_process <- factor(go$`GO biological process complete`, 
                        levels = go$`GO biological process complete`[order(go$fold_enrichment)])

full_terms = ggplot(go, aes(x = fold_enrichment, y =GO_process, color=Padj,size=found)) + 
  geom_point() + labs(x = "Fold Enrichment", y = "GO Term") + 
  guides(size = guide_legend(title = "Genes Found")) + 
  theme_minimal()  +  theme(axis.text.y = element_text(size = 8))

if(!file.exists('nmd_goTerms.pdf')){
  ggsave('nmd_goTerms.pdf',plot = full_terms, height = 4, width = 8, units = 'in')
}
