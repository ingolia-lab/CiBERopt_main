library(ggplot2)

setwd('~/CiBERopt_main/CiberOpt_fig2/bxb1_efficiency/')
counts = read.csv("transformation_efficiency.csv", stringsAsFactors=FALSE, row.names=1)
counts$date = sub("[A-Z]$", "", row.names(counts))
counts$effic = counts$Count * counts$Dilution / counts$fmol

bp = ggplot(data= counts, aes(x=Plasmid,y=effic,fill='#2171b5'))+ 
  geom_bar(position='dodge', stat='summary',fun.y='mean')+ 
  scale_fill_manual(values='#2171b5') + 
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') + 
  geom_point(aes(x = Plasmid), shape = 16, size = 2, 
             position = position_jitterdodge()) + 
  scale_y_continuous(expand=c(0.0,0.0),limits = c(1,30000),trans = 'log10') +
  ylab('Efficiency') + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=18),
        axis.title.y=element_text(size=20),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('bxb1_eff.pdf')){
  ggsave('bxb1_eff.pdf',bp, width = 4.5, height=6, unit='in')
}

  
  

