library(tidyverse)

read.peri <- function(indata) {
  indata <- indata %>% group_by(gene, seqlen, dis5, frame5, rep) %>%
    summarise(reads=sum(seqsnum)) %>% ungroup()
  indata$rep <- paste0('r', indata$rep)
  inreps <- unique(indata$rep)
  indata <- indata %>% spread(rep, reads)
  indata[is.na(indata)] <- 0
  indata$rpm <- rpm.mean(indata, inreps)
  return(indata)
}

plot.peri <- function(indata, intitle, seqrange=c(-Inf, Inf)) {
  indata <- indata[indata$seqlen >= seqrange[1] & indata$seqlen <= seqrange[2], ]
  return(
    ggplot(indata, aes(x=as.factor(seqlen), weight=rpm, fill=as.factor(frame5)))+
      geom_bar(position='dodge')+theme_classic()+
      scale_fill_manual(values=c('#6699CC', '#D6EAF8', '#D6DBDF'), 
                        name='Frame')+
      theme(legend.position=c(.05, .8),
            legend.background=element_blank(),
            legend.text=element_text(size=8),
            legend.title=element_text(size=10),
            legend.key.size=unit(1, 'line'),
            axis.title=element_text(size=12),
            axis.text=element_text(size=10))+
      xlab('Fragment Length')+
      ylab('Normalised RPM')+
      labs(title=intitle)
    )
}

disome.ori <- read.table('disome_ori.txt', header=T, sep='\t', stringsAsFactors=F)
disome.peri <- read.peri(disome.ori[disome.ori$rep != 0, ])
rm(disome.ori)

mono.ori <- read.table('monosome_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.peri <- read.peri(mono.ori)
rm(mono.ori)

mrna.ori <- read.table('mRNA_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mrna.peri <- read.peri(mrna.ori)
rm(mrna.ori)

di.3at.ori <- read.table('di_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
di.3at.peri <- read.peri(di.3at.ori)
rm(di.3at.ori)

mono.3at.ori <- read.table('mono_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.3at.peri <- read.peri(mono.3at.ori)
rm(mono.3at.ori)

mrna.3at.ori <- read.table('mRNA_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mrna.3at.peri <- read.peri(mrna.3at.ori)
rm(mrna.3at.ori)

{
  pdf('P1_periodicity.pdf', width=8, height=3)
  plot(plot.peri(disome.peri, 'Disome', c(35, 70)))
  plot(plot.peri(mono.peri, 'Monosome'))
  plot(plot.peri(mrna.peri, 'mRNA'))
  plot(plot.peri(di.3at.peri, 'Disome, 3-AT', c(35, 70)))
  plot(plot.peri(mono.3at.peri, 'Monosome, 3-AT'))
  plot(plot.peri(mrna.3at.peri, 'mRNA, 3-AT'))
  dev.off()
}