library(ggplot2)
library(ggpubr)
setwd("~/CiBERopt_main/CiberOpt_fig3/CL1_doa10_qpcr/")

cq = read.csv("wt_cl1_doa10_qpcr.csv")
cq = aggregate(Ct..dR.~name, cq,mean)
cq$biorep = paste0(sapply(strsplit(cq$name,'_'),'[[',1),'_',sapply(strsplit(cq$name,'_'),'[[',3))
cq$target = paste0(sapply(strsplit(cq$name,'_'),'[[',2))

df = as.data.frame.matrix(t(xtabs(Ct..dR.~target + biorep,cq)))
df$cq_diff = df$YFP-df$CFP
df$amount = 2^-df$cq_diff
df$cond =  sapply(strsplit(rownames(df),'_'),'[[',1)

#manually inspect p values, label in illustrator because it's easier
#perform on unexponentiated cq diff
pvals = compare_means(cq_diff~cond,df,method='t.test',paired = F)
#manually inspect
df$cond = factor(df$cond,levels=c('wt','cl1','doa10'))

p = ggplot(df, aes(x=cond,y=amount, fill='black')) +
  geom_bar(position='dodge', stat='summary',fun.y='mean')+ 
  scale_fill_manual(values='#2171b5') + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9), size=2) + 
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') +
  scale_y_continuous(expand=c(0,0),limits =c(0,0.12)) + ylab('Relative reporter expression') + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=18),
        axis.title.y=element_text(size=20),
        legend.position = "none",
        panel.border = element_blank())

if(!(file.exists('cl1_qpcr.pdf'))){
  ggsave('cl1_qpcr.pdf',p,height = 4,width = 3, units='in')
}

