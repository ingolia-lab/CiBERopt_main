library(ggplot2)

setwd("~/CiBERopt_main/CiberOpt_fig3/cl1_guide_validation/")

cq = read.csv("cl1_guides_qpcr.csv")
#collapse technical replicates
cq = aggregate(Ct..dR. ~ label, cq, mean)
cq$target = sapply(strsplit(cq$label,'_'),'[[',3)
cq$condition = paste0(sapply(strsplit(cq$label,'_'),'[[',1),'_',sapply(strsplit(cq$label,'_'),'[[',2))

df = as.data.frame.matrix(t(xtabs(Ct..dR. ~ target + condition,cq)))

df$condition = sapply(strsplit(rownames(df),'_'),'[[',1)
df$diff = 2^-(df$YFP-df$CFP)

df$condition = factor(df$condition, levels=c('empty','ubc6','rpt4','cdc48'))

qpcr_plot = ggplot(df, aes(x=condition, y=diff, fill='black')) + 
  geom_bar(position='dodge', stat='summary',fun.y='mean', fill='#2171b5')+ 
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') + 
  geom_point(aes(x = condition), shape = 16, size = 2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) + 
  theme_classic() + ylab('Relative Reporter Expression') + xlab('') + 
  scale_y_continuous(expand=c(0.0,0.0),limits=c(0,93)) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=18),
        axis.title.y=element_text(size=20),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('cl1_validation.pdf')){
  ggsave('cl1_validation.pdf',qpcr_plot, height=3,width=2.5, units = 'in')
}
