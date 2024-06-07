library(ggplot2)
library(flowCore)
library(minpack.lm)

setwd('~/CiBERopt_main/')
load(file = 'gateFlow.R')
setwd("~/CiBERopt_main/CiberOpt_fig1/hormoneTitration_z3z4/")


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
                'hormone' = sapply(strsplit(file_names,'_'),'[[',4),
                'yfp' = sapply(gated_tables,function(x){mean(x[,'FITC.A'])}))

#add hormone concentrations
df$hormone = rep(c(0,0,100,100,12.5,12.5,200,200,25,25,3.125,3.125,50,50,6.25,6.25),2)

#find model fits here
z3_df = df[grepl('z3',df$strain),]
z4_df = df[grepl('z4',df$strain),]

z3_fit = nlsLM(yfp~ A*hormone/(hormone + k),data = z3_df, start=list(A=25000, k =25))
z4_fit = nlsLM(yfp~ A*hormone/(hormone + k),data = z4_df, start=list(A=25000, k =25))

kds = paste(paste0('Z3 = ', round(coefficients(z3_fit)[2],2)), 
              paste0('Z4 = ', round(coefficients(z4_fit)[2],2)), sep='\n')

titration = ggplot(df, aes(x=hormone, y=yfp/1000, color=strain)) + geom_point(size=4) + 
  stat_smooth(method = 'nls', formula = y~a*(x/(x+k)),
              method.args = list(start=c(a=25, k=25)), se=FALSE,size=1) + 
  scale_color_manual(values=c('#e41a1c','#377eb8')) + 
  scale_y_continuous(expand=c(0.01,0.01)) +
  theme_classic() + xlab('[Progesterone] (nM)') + ylab('YFP*1000') + 
  annotate(geom="text", x=100, y=8, label=kds) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        legend.position="none",
        axis.ticks.length = unit(0.1,'cm'),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        panel.border = element_blank())


if(!(file.exists('titrationZ3Z4.pdf'))){
  ggsave('titrationZ3Z4.pdf',titration,width = 5, height = 3, units='in')
}

##########################
######for supp figs#######
##########################
#gating scheme for fwd scatter Height/Area
fsc_plot = ggplot(list_tables[[1]], aes(FSC.A/1000,FSC.H/1000, color = good_ah)) +
  rasterize(geom_point(),dpi=600) + scale_color_manual(values = c('#bdbdbd','#8856a7')) + 
  theme_classic() + ylab('FSC-H x1000') + xlab('FSC-A x1000') + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        panel.border = element_blank())

#gating scheme for fwd/side scatter area
ssc_plot = ggplot(list_tables[[1]], aes(FSC.A/1000,SSC.A/1000, color = good_ssc)) +
  rasterize(geom_point(),dpi=600) + scale_color_manual(values = c('#bdbdbd','#8856a7')) + 
  theme_classic() + ylab('FSC-H x1000') + xlab('SSC-A x1000') + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        panel.border = element_blank())


#density plot for induced/uninduced, representative 
data = gated_tables[c(1,7,17,23)]
ids = c('z3_zero','z3_pro','z4_zero','z4_pro')
for(i in seq_along(data)){
  data[[i]]$id = ids[i]
}

data = do.call(rbind, data)
data$group = ifelse(grepl('z3', data$id),'z3','z4')

data$id = factor(data$id, levels=ids)


prog_plots = ggplot(data, aes(x=FITC.A/1000,..scaled..)) + geom_density(aes(fill=id),alpha=0.5) +  
  scale_fill_manual(values=c('#bdbdbd', '#8856a7', '#bdbdbd', '#8856a7')) +
  scale_color_manual(values=c('#bdbdbd', '#8856a7', '#bdbdbd', '#8856a7')) + 
  facet_wrap(~group, ncol=1) + theme_classic() + xlab('YFP x1000') + ylab('Density') + 
  scale_y_continuous(expand=c(0,0)) + scale_x_log10(limits=c(0.01,1000), expand=c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none",
        panel.border = element_blank())


setwd("~/CiBERopt_main/CiBERsupp/ciberOpt_supp1/")
if(!(file.exists('fsc_plot.pdf'))){
  ggsave('fsc_plot.pdf',fsc_plot,height=5,width=5,units='in')
}
if(!(file.exists('ssc_plot.pdf'))){
  ggsave('ssc_plot.pdf',ssc_plot,height=5,width=5,units='in')
}
if(!(file.exists('progesterone_addition.pdf'))){
  ggsave('progesterone_addition.pdf',prog_plots,height=5,width=5,units='in')
}