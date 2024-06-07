library(ggplot2)
library(stringr)

setwd("~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/DNAgoTerms/")

go = read.csv("GOterms_DNA.csv")
colnames(go) = go[7,]
go = go[8:nrow(go),]

###generate full table for supplement
go$Padj = as.numeric(go$`sig_up_dna_ids.txt (P-value)`)
go$found = as.numeric(go$`sig_up_dna_ids.txt (429)`)
go$fold_enrichment = as.numeric(go$`sig_up_dna_ids.txt (fold Enrichment)`)
go$GO_process <- factor(go$`GO biological process complete`, 
                              levels = go$`GO biological process complete`[order(go$fold_enrichment)])



full_terms = ggplot(go, aes(x = fold_enrichment, y =GO_process, color=Padj,size=found)) + 
  geom_point() + labs(x = "Fold Enrichment", y = "GO Term") + 
  guides(size = guide_legend(title = "Genes Found")) + 
  theme_minimal()  +  theme(axis.text.y = element_text(size = 8))



