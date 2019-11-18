library(tidyverse)
library(parallel)

##
dis5.sum <- mapply(function(indata, tag){
  indata <- peak.sum(indata, groupcol='dis5')
  repcol <- colnames(indata)[grep('r\\d+', colnames(indata), perl=T)]
  print(repcol)
  indata$rpm <- rpm.mean(indata, repcol)
  indata$tag <- tag
  indata[, repcol] <- NULL
  indata <- indata %>% group_by(gene) %>%
    mutate(rpm.nor=rpm/sum(rpm))
  indata <- as.data.frame(indata)
  return(indata)
}, list(disome.data, mono.data), c('disome', 'mono'), SIMPLIFY=F)



## read disome-seq and monosome-seq data ####
disome.peak <- read.table('./disome_peak.txt', header=T, sep='\t', stringsAsFactors=F)
disome.perc <- disome.peak[, c('gene', 'asite', 'read', 'rpm')]
disome.perc <- merge(disome.perc, cds.pep[, c('gene', 'len')], all.x=T, by='gene')
disome.perc$tosc <- disome.perc$asite - (disome.perc$len - 3)
disome.perc$asite.perc <- disome.perc$asite / (disome.perc$len - 3)
disome.perc$sc <- ifelse(disome.perc$asite.perc == 1, 1, 0)
# disome.perc$len <- NULL
disome.perc <- disome.perc %>% group_by(gene) %>%
  mutate(rpm.nor=rpm/sum(rpm)) %>% ungroup() %>%
  as.data.frame()

mono.peak <- read.table('./mono_peak.txt', header=T, sep='\t', stringsAsFactors=F)
mono.perc <- mono.peak[, c('gene', 'asite', 'read', 'rpm')]
mono.perc <- merge(mono.perc, cds.pep[, c('gene', 'len')], all.x=T, by='gene')
mono.perc$tosc <- mono.perc$asite - (mono.perc$len - 3)
mono.perc$asite.perc <- mono.perc$asite / (mono.perc$len - 3)
mono.perc$sc <- ifelse(mono.perc$asite.perc == 1, 1, 0)
mono.perc <- mono.perc %>% group_by(gene) %>%
  mutate(rpm.nor=rpm/sum(rpm)) %>% ungroup() %>%
  as.data.frame()


dis5.sum <- do.call('rbind', dis5.sum)
dis5.sum <- merge(dis5.sum, cds.pep[, c('gene', 'len')], all.x=T, by='gene')
dis5.sum$tosc <- dis5.sum$dis5 - (dis5.sum$len - 3)
dis5.sum$len <- NULL
head(dis5.sum)


plot.dis <- function(indata, inx, iny, low=-Inf, up=Inf, inweight='rpm.nor', indis='tosc') {
  oplot <- ggplot(indata[indata[, indis] >= low & indata[, indis] <= up, ],
                  aes_string(x=indis, weight=inweight, color='tag'))+
    geom_freqpoly(binwidth=1, show.legend=F)+theme_classic()+
    scale_color_manual(values=c("#FF9900", "#3399FF"))+
    # facet_wrap(~tag, scales='free_y', ncol=1, labeller=labeller(tag=c(disome='Disome', mono='Monosome')))+
    facet_wrap(~tag, scales='free_y', ncol=1)+
    geom_text(aes(label=lab), size=3,
              data=data.frame(lab=c('Disome', 'Monosome'),
                              tag=c('disome', 'mono')),
              x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, inherit.aes = F)+
    theme(strip.background=element_blank(),
          strip.text.x=element_blank(),
          axis.title=element_text(size=9),
          axis.text=element_text(size=8),
          axis.line=element_line(size=0.4),
          axis.ticks=element_line(size=0.4))+
    xlab(inx)+
    ylab(iny)
  if (indis == 'tosc'){
    oplot <- oplot+
      geom_vline(xintercept=c(-15, -45, -48, -75, -78),
                 linetype='dashed', color='#66cc00')
  }
  return(oplot)
}

plot.rela.dis <- function(indata1, indata2, intag1, intag2, inweight, iny){
  indata1$inweight <- indata1[, inweight] / (length(unique(indata1$gene))/100)
  indata2$inweight <- indata2[, inweight] / (length(unique(indata2$gene))/100)
  indata <- rbind(data.frame(indata1[, c('gene', 'inweight', 'asite.perc')], tag=intag1),
                  data.frame(indata2[, c('gene', 'inweight', 'asite.perc')], tag=intag2))

  return(
    ggplot(indata %>%
             filter(asite.perc>-0.1 & asite.perc<1.1),
           aes(x=asite.perc, weight=inweight, color=tag))+
      geom_freqpoly(binwidth=0.01, show.legend=F)+theme_classic()+
      geom_vline(xintercept=0, linetype='dashed', color='#006633')+
      geom_vline(xintercept=1, linetype='dashed', color='#66cc00')+
      scale_color_manual(values=c("#FF9900", "#3399FF"))+
      facet_wrap(~tag, scales='free_y', ncol=1)+
      geom_text(aes(label=lab), size=3,
                data=data.frame(lab=c(intag1, intag2),
                                tag=c(intag1, intag2)),
                x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, inherit.aes = F)+
      theme(strip.background=element_blank(),
            strip.text.x=element_blank(),
            axis.title=element_text(size=9, color='black'),
            axis.text=element_text(size=8, color='black'),
            axis.line=element_line(size=0.4),
            axis.ticks=element_line(size=0.4))+
      xlab('Relative A-site on CDS')+
      ylab(iny)
  )
}

geno.rscu <- cal.rscu(500, as.character(unique(disome.perc$gene)))

{
  pdf('P1_rpm_distance.pdf', width=5, height=2.5, useDingbats=F)
  plot(
    plot.dis(dis5.sum, 'Distance from the first base of stop codon', 'Normalised RPM', -100, 10)
  )
  plot(
    plot.dis(dis5.sum, 'Distance from the first base of start codon', 'Normalised RPM', -50, 200, indis='dis5')
  )
  plot(
    plot.rela.dis(disome.perc, mono.perc, 'Disome', 'Monosome', 'rpm.nor', 'Normalised RPM')
  )
  dev.off()
  }

write.table(disome.perc, 'disome.perc.txt', sep='\t', quote=F, row.names=F)
write.table(mono.perc, 'mono.perc.txt', sep='\t', quote=F, row.names=F)


## stop codon percentage####
table(cds.sc$stop_codon)

ggplot(data=data.frame(sc=factor(c('TAA', 'TGA', 'TAG'), levels=c('TAG', 'TGA', 'TAA')),
                       num=c(2793, 1742, 1349)),
       aes(x=factor(1), y=num, fill=sc))+
  geom_bar(stat='identity', size=1)+coord_polar(theta='y')+
  theme_classic()+
  theme(legend.position='top',
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_text(size=12),
        axis.text.x=element_text(size=10),
        legend.background=element_blank(),
        legend.text=element_text(size=10),
        legend.key.size=unit(.8, 'cm'))+
  ylab('# of Yeast genes (5884 in Total)')+
  scale_y_continuous(breaks=c(0, 2793, 4535, 5884))+
  scale_fill_brewer(palette='Set2')

mrna.read.ordered.gene <- mrna.read[order(mrna.read$rpm, decreasing=T),]$gene
table(cds.sc[cds.sc$gene %in% mrna.read.ordered.gene[1:589],]$stop_codon)

{
  pdf('p2_stop_codon_percentage.pdf', width=7.5, height=3.6, useDingbats=F)
  multiplot(ggplot(data=data.frame(sc=factor(c('TAA', 'TGA', 'TAG'), levels=c('TAG', 'TGA', 'TAA')),
                                   num=c(2793, 1742, 1349)),
                   aes(x=factor(1), y=num, fill=sc))+
              geom_bar(stat='identity', size=1)+coord_polar(theta='y')+
              theme_classic()+
              theme(legend.position='top',
                    axis.line=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_text(size=12),
                    axis.text.x=element_text(size=10),
                    legend.background=element_blank(),
                    legend.text=element_text(size=10),
                    legend.key.size=unit(.8, 'cm'))+
              ylab('# of Yeast genes (5884 in Total)')+
              scale_y_continuous(breaks=c(0, 2793, 4535, 5884))+
              scale_fill_brewer(palette='Set2'),
            ggplot(data=data.frame(sc=factor(c('TAA', 'TGA', 'TAG'), levels=c('TAG', 'TGA', 'TAA')),
                                   num=c(385, 105, 99)),
                   aes(x=factor(1), y=num, fill=sc))+
              geom_bar(stat='identity', size=1)+coord_polar(theta='y')+
              theme_classic()+
              theme(legend.position='top',
                    axis.line=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_text(size=12),
                    axis.text.x=element_text(size=10),
                    legend.background=element_blank(),
                    legend.text=element_text(size=10),
                    legend.key.size=unit(.8, 'cm'))+
              ylab('# of top 1/10 expression Yeast genes (589 in Total)')+
              scale_y_continuous(breaks=c(0, 385, 385+105, 385+105+99))+
              scale_fill_brewer(palette='Set2'),
            cols=2)
  dev.off()
}
## the end ####


