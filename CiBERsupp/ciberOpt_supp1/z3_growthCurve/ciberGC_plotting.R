
library(ggplot2)
library(tidyr)
setwd("~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/z3_growthCurve/")

gc =read.csv('hormoneGC.csv')
gc$Time..s. = as.numeric(substr(gc$Time..s.,1,(nchar(gc$Time..s.)-1)))

gc = gather(gc, key = "Name", value = "OD", -Time..s.)
gc$group = sapply(strsplit(gc$Name,'_'),'[[',1)

blues= c('#deebf7','#9ecae1','#3182bd')
reds = c('#fee0d2','#fc9272','#de2d26')


grow_plot = ggplot(gc, aes(x=Time..s./3600,y=OD, color=Name)) + 
  geom_line(size=1.5)  + 
  facet_wrap(~group,ncol=1) + 
  scale_color_manual(values=rep(c(blues,reds),2)) + theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size=14),
        legend.position = 'none',
        axis.title = element_text(size=16),
        panel.border = element_blank())

if(!file.exists('prog_gc.pdf')){
  ggsave('prog_gc.pdf',grow_plot,height=5,width=5,units='in')
}
  



