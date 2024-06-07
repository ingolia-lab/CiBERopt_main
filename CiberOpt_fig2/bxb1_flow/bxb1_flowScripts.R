library(ggplot2)
library(ggrastr)
library(ggpointdensity)
library(RColorBrewer)
library(tidyr)

setwd('~/CiBERopt_main/')
load(file = 'gateFlow.R')
setwd("~/CiBERopt_main/CiberOpt_fig2/bxb1_flow/")

file_names = list.files(pattern ='fcs')
#didn't collect FSC.H, so do the gating only of fsc.a/ssc.a. Collected on different instrument so absolute vals diff
list_tables <- lapply(file_names, 
                      function(x){data.frame(exprs(read.FCS(x,transformation = FALSE)))
                      })

get_ssc <- function(df, fsca_low=60000, fsca_high=100000,rat = 1.4){
  #ssc/fsca gated within 2 sigma of a linear fit and extreme high/lows trimmer
  ssc_lm = lm(df$SSC.A~df$FSC.A)
  df$good_ssc = ifelse(df$SSC.A < ssc_lm$coefficients[2]* df$FSC.A + ssc_lm$coefficients[1]+(2*sigma(ssc_lm)) &
                         df$SSC.A > ssc_lm$coefficients[2]* df$FSC.A + ssc_lm$coefficients[1]-(2*sigma(ssc_lm)) &
                         df$FSC.A > fsca_low & df$FSC.A < fsca_high, TRUE, FALSE)
  return(df)
}

df = lapply(list_tables,get_ssc)

#after applying gating
gated_tables = lapply(df, function(x){x[x$good_ssc == TRUE & x$FITC.A > 0 & x$PE.Texas.Red.A>0,  ]})
dna = c('Bxb1','Bxb1','CEN','CEN')
media = c('select','unsel','select','unsel')
for(i in 1:length(dna)){
  gated_tables[[i]]$dna = dna[i]
  gated_tables[[i]]$media = media[i]
  gated_tables[[i]]$name = paste0(dna[i],'_',media[i])
}

gated_tables = do.call(rbind,gated_tables)

cols = brewer.pal(9,'Blues')

flow_2d_plot = ggplot(gated_tables, aes(x=FITC.A,y=PE.Texas.Red.A),) + 
  geom_pointdensity(size=0.5, show.legend = FALSE) + 
  scale_alpha(range=c(0.1,1)) + 
  scale_color_gradientn(colors = cols[4:9]) + 
  geom_hline(yintercept = 2000, color='#bdbdbd',size=1.2) + 
  geom_vline(xintercept = 1000, color='#bdbdbd',size=1.2) + 
  scale_y_log10() + scale_x_log10() + 
  coord_fixed() + 
  xlab('Citrine') + ylab('mScarlet-i') + 
  facet_wrap(~name) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

flow_2d_plot = rasterize(flow_2d_plot, layers='Point',dpi=600)

#for the selected fluorescence, compare citrine distributions. Take yeast that express only YFP (none of the doubles)
select_yfp = gated_tables[gated_tables$FITC.A > 1000 & 
                            gated_tables$PE.Texas.Red.A < 2000 & 
                            gated_tables$media =='select',]

dens = ggplot(select_yfp, aes(x=FITC.A,color=dna,fill=dna)) + 
  geom_density(alpha=0.5) + scale_x_log10() + scale_y_continuous(expand = c(0,0)) + 
  scale_color_manual(values=c('#e41a1c','#377eb8')) +
  scale_fill_manual(values=c('#e41a1c','#377eb8')) + 
  theme_classic() + ylab('') + xlab('Citrine (AU)') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size=20)
        )

if(!(file.exists('bxb1_flow.pdf'))){
  ggsave('bxb1_flow.pdf',flow_2d_plot, height = 6, width = 6, units = 'in')
}

if(!(file.exists('bxb1_dens.pdf'))){
  ggsave('bxb1_dens.pdf',dens, height = 4, width = 6, units = 'in')
}

