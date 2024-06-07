library(ggplot2)
library(ggpubr)
setwd("~/CiBERopt_main/CiberOpt_fig4/nmd_libraryQPCR/")

cq = read.csv("nmd_lib_qPCR.csv")

#average technical replicates
df = aggregate(Ct..dR.~strain+target+condition+bio_rep, cq, mean)
#set row with condition+strain
df$label = paste(df$strain, df$condition, df$bio_rep, sep='_')
df = as.data.frame.matrix(xtabs(Ct..dR. ~ label + target, df))
df$cq_diff = df$cit-df$scarlet
df$diff_amount = 2^-(df$cq_diff)
df$condition = sapply(strsplit(rownames(df),'_'),'[[',2)

df$condition = factor(df$condition, levels=c('pre','post','chx'))

lib_qpcr = ggplot(df, aes(x=condition, y=diff_amount, fill='black')) + 
  geom_bar(position='dodge', stat='summary',fun.y='mean', fill='#2171b5')+ 
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') + 
  geom_point(aes(x = condition), shape = 16, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) + 
  theme_classic() + ylab('YFP PTC mRNA Amount') + xlab('') + 
  scale_y_continuous(expand=c(0.0,0.0),limits = c(0,1.3)) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=18),
        axis.title.y=element_text(size=20),
        legend.position = "none",
        panel.border = element_blank())

#evaluate statistics on the cq_dff, all CHX < 0.05
pvals = compare_means(cq_diff~condition, df, method='t.test')

if(!(file.exists('nmd_lib_qpcr.pdf'))){
  ggsave('nmd_lib_qpcr.pdf',lib_qpcr, height=3,width=2.5)
}
