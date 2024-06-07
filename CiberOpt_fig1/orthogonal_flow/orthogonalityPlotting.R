library(ggplot2)
library(flowCore)
library(minpack.lm)

setwd('~/CiBERopt_main/')
load(file = 'gateFlow.R')
setwd("~/CiBERopt_main/CiberOpt_fig1/orthogonal_flow/")


file_names = list.files(pattern ='fcs')

list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE & x$FITC.A > 0 ,  ]})


df = data.frame()

df = data.frame('strain' = sapply(strsplit(file_names,'_'),'[[',3),
                'yfp' = sapply(gated_tables,function(x){mean(x[,'FITC.A'])}))

df = aggregate(yfp ~ strain, df, mean)

df$tf = substr(df$strain,1,2)
df$promoter = substr(df$strain,3,4)

#normalize data to max of each TF
df$yfp = ifelse(df$tf == 'z3', df$yfp/max(df[df$tf=='z3','yfp']),df$yfp/max(df[df$tf=='z4','yfp']))

df$tf <- factor(df$tf, levels=c('z3','z4'))
df$promoter <- factor(df$promoter, levels=c('p4','p3'))


ortho_plot = ggplot(df, aes(tf, promoter, fill= yfp)) + 
  geom_tile() + scale_fill_gradientn(colors = c("#bdbdbd", "#8856a7")) + 
  geom_tile(color = "black", size = 0.5) + 
  theme_classic() + xlab('Transcription Factor') + ylab('Promoter') + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        panel.border = element_blank())  

if(!(file.exists('orthoPlot.pdf'))){
  ggsave('orthoPlot.pdf',ortho_plot, height = 5, width = 5, units='in')
}

