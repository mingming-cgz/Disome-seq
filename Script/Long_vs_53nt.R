library(tidyverse)
library(data.table)

## rich medium ####
dishort <- disome.ori[disome.ori$seqlen %in% 53:54, ]

table(dishort$rep)
# 1    2 
# 5547 3693
dishort[, c('sequnum', 'frame3', 'putative', 'region', 'dupl')] <- NULL
rownames(dishort) <- NULL

dishort$rawlen <- dishort$seqlen
dishort$seqlen <- dishort$rawlen + 5

dishort$cor.pos <- correct.pos(dishort)
dishort$asite <- get.asite(dishort)

dishort.peak <- dishort %>% group_by(rep, gene, asite) %>%
  summarise(read=sum(seqsnum)) %>% ungroup()
dishort.peak$rep <- paste0('r', dishort.peak$rep)
dishort.peak <- dishort.peak %>% spread(rep, read)
dishort.peak[is.na(dishort.peak)] <- 0
dishort.peak$read <- rowSums(dishort.peak[, c('r1', 'r2')])

dishort.peak <- merge(dishort.peak, cds.len, by='gene', all.x=T)
dishort.peak$sc <- ifelse(dishort.peak$len-3 == dishort.peak$asite, 1, 0)
dishort.peak$sc[dishort.peak$asite > dishort.peak$len-3] <- -1
dishort.peak$peakshare <- ifelse(dishort.peak$r1>0 & dishort.peak$r2>0, 1, 0)
dishort.peak$rpm1 <- dishort.peak$r1 / (sum(dishort.peak$r1)/1e6)
dishort.peak$rpm2 <- dishort.peak$r2 / (sum(dishort.peak$r2)/1e6)
dishort.peak$rpm <- rowMeans(dishort.peak[, c('rpm1', 'rpm2')])
dishort.peak <- dishort.peak %>% group_by(gene) %>%
  mutate(rpm.nor=rpm/sum(rpm)) %>% as.data.frame()

dishort.peak$asite.perc <- dishort.peak$asite / (dishort.peak$len-3)
dishort.peak$len <- NULL

length(unique(dishort.peak$gene))
# [1] 1876

{
  pdf('P3_5354nt_asite_nor.pdf', width = 2.2, height = 1.2)
  plot(
    ggplot(dishort.peak, aes(x=asite.perc, weight=rpm.nor/(1876/100)))+
      geom_freqpoly(binwidth=0.01, color='#FF9900')+theme_classic()+
      # geom_vline(xintercept = 1)+
      theme(axis.title = element_text(size=8, color='black'),
            axis.text = element_text(size=8, color='black'),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none')+
      xlab('Relative position on CDS')+ylab('Density')+
      scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.1, 1.1),
                         labels = c('Start', '50%', 'Stop'))
  )
  dev.off()
}


dishort.asite <- dishort.peak[, c('gene', 'asite', 'r1', 'r2', 'read', 'sc', 'peakshare')]
dishort.asite$aaseq <- peak.seq(dishort.asite, 'exit', 4, asite=T)
dishort.asite$codonseq <- peak.seq(dishort.asite, 'exit', 4, aa=F, asite=T)
dishort.asite <- dishort.asite[dishort.asite$sc != -1, ]

#### asite
# bgmrna.asite.codon.read <- enrich.posi('p4', diasite.codon.read, mrnaasite.codon.read, overlapgene=F, codon=T)
dshortasite.codon.read <- count_site(dishort.asite, 'codonseq', 1, aa=F, count_read=T)
bgmrna.ds.asite.codon.read <- enrich.posi('p4', dshortasite.codon.read, mrnaasite.codon.read, overlapgene=F, codon=T)


vsrich58.rich53.46 <- merge(bgmrna.ds.asite.codon.read[, c('mer', 'genenum', 'oddsratio', 'pvalue'),],
                            bgmrna.asite.codon.read[, c('mer', 'genenum', 'oddsratio', 'pvalue', 'optimality'),],
                            by=c('mer'), all=T, suffixes=c('.rich53', '.rich58'))
vsrich58.rich53.46 <- vsrich58.rich53.46[vsrich58.rich53.46$mer != '   ',]


vsr58.r53.a46.cor <- cor.test(vsrich58.rich53.46[!vsrich58.rich53.46$mer %in% c('TAA', 'TGA', 'TAG'),]$oddsratio.rich53,
                              vsrich58.rich53.46[!vsrich58.rich53.46$mer %in% c('TAA', 'TGA', 'TAG'),]$oddsratio.rich58, 
                              method = 's')

{
  pdf('P3_rich53vsrich58.pdf', width = 2.4, height = 2.63, useDingbats = F)
  plot(
    ggplot(vsrich58.rich53.46, 
           aes(x=log2(oddsratio.rich58), y=log2(oddsratio.rich53),
               color=optimality))+
      geom_point(size=1)+theme_classic()+
      geom_text(data=vsrich58.rich53.46[vsrich58.rich53.46$mer %in% c('TAG', 'TAA', 'TGA'),],
                aes(label=mer))+
      geom_abline(slope = 1, intercept = 0, linetype='dashed', color='grey')+
      xlim(-5, 6)+ylim(-5, 6)+
      labs(title='Disome, rich 58nt vs rich 53nt, 46-48nt asite')+
      scale_color_manual(values=c('black', '#3399ff'))+
      theme(axis.title = element_text(size=8, color='black'),
            axis.text = element_text(size=8, color='black'),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none')+
      xlab('Log2(A-site pausing score), 58nt, rich')+
      ylab('Log2(A-site pausing score), 53nt, rich')+
      annotate('text', x=-3, y=4, 
               label=sprintf('rho = %.2f\nP = %.2e\nN = 61',
                             vsr58.r53.a46.cor$estimate, 
                             vsr58.r53.a46.cor$p.value))
  )
  dev.off()
}


dspaa.aa.read.rmsc <- count_site(dishort.asite[dishort.asite$sc == 0, ], 'aaseq', 1, aa=T, count_read=T)
bgmrna.dspaa.read <- enrich.posi('p3', dspaa.aa.read.rmsc, mrnaasite.aa.read.rmsc, overlapgene=F, codon=F)

vsrich58.rich53.psite <- merge(bgmrna.psite.aa.read.rmsc[, c('mer', 'genenum', 'oddsratio', 'pvalue'),],
                               bgmrna.dspaa.read[, c('mer', 'genenum', 'oddsratio', 'pvalue'),],
                            by=c('mer'), all=T, suffixes=c('.rich58', '.rich53'))
vsrich58.rich53.psite <- vsrich58.rich53.psite[vsrich58.rich53.psite$mer != ' ',]
vsr58.r53.paa.cor <- cor.test(log2(vsrich58.rich53.psite$oddsratio.rich58), 
                              log2(vsrich58.rich53.psite$oddsratio.rich53), method = 's')


{
  pdf('P3_rich53vsrich58_psite.pdf', width = 2.4, height = 2.63, useDingbats = F)
  plot(
    ggplot(vsrich58.rich53.psite, 
           aes(log2(oddsratio.rich58), log2(oddsratio.rich53)))+
      theme_classic()+
      geom_point(size=1)+theme_classic()+
      geom_text(data=vsrich58.rich53.psite[vsrich58.rich53.psite$mer == 'P',],
                aes(label=mer))+
      geom_abline(slope = 1, intercept = 0, linetype='dashed', color='grey')+
      xlim(-0.83, 1.432)+ylim(-0.83, 1.432)+
      labs(title='Disome, rich 58nt vs rich 53nt, Psite')+
      scale_color_manual(values=c('black', '#3399ff'))+
      theme(axis.title = element_text(size=8, color='black'),
            axis.text = element_text(size=8, color='black'),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none')+
      xlab('Log2(P-site pausing score), 58nt, rich')+
      ylab('Log2(P-site pausing score), 53nt, rich')+
      annotate('text', x=-0.5, y=1.3, 
               label=sprintf('rho = %.2f\nP = %.2e\nN = 20',
                             vsr58.r53.paa.cor$estimate, 
                             vsr58.r53.paa.cor$p.value))
  )
  dev.off()
  }


## 3at ####
di3atshort <- di.3at.ori[di.3at.ori$seqlen %in% 53:54,]
table(di3atshort$rep)
# 1    2 
# 3331 3028


di3atshort[, c('sequnum', 'frame3', 'putative', 'region', 'dupl')] <- NULL
rownames(di3atshort) <- NULL

di3atshort$rawlen <- di3atshort$seqlen
di3atshort$seqlen <- di3atshort$rawlen + 5

di3atshort$cor.pos <- correct.pos(di3atshort)
di3atshort$asite <- get.asite(di3atshort)

di3atshort.peak <- di3atshort %>% group_by(rep, gene, asite) %>%
  summarise(read=sum(seqsnum)) %>% ungroup()
di3atshort.peak$rep <- paste0('r', di3atshort.peak$rep)
di3atshort.peak <- di3atshort.peak %>% spread(rep, read)
di3atshort.peak[is.na(di3atshort.peak)] <- 0
di3atshort.peak$read <- rowSums(di3atshort.peak[, c('r1', 'r2')])

di3atshort.peak <- merge(di3atshort.peak, cds.len, by='gene', all.x=T)
di3atshort.peak$sc <- ifelse(di3atshort.peak$len-3 == di3atshort.peak$asite, 1, 0)
di3atshort.peak$sc[di3atshort.peak$asite > di3atshort.peak$len-3] <- -1
di3atshort.peak$peakshare <- ifelse(di3atshort.peak$r1>0 & di3atshort.peak$r2>0, 1, 0)

di3atshort.peak$len <- NULL

di3atshort.asite <- di3atshort.peak
di3atshort.asite$aaseq <- peak.seq(di3atshort.asite, 'exit', 4, asite=T)
di3atshort.asite$codonseq <- peak.seq(di3atshort.asite, 'exit', 4, aa=F, asite=T)
di3atshort.asite <- di3atshort.asite[di3atshort.asite$sc != -1, ]

ds3at.acodon.read <- count_site(di3atshort.asite, 'codonseq', 1, aa=F, count_read=T)
bgmrna.ds3at.acodon.read <- enrich.posi('p4', ds3at.acodon.read, mrnaasite.3at.codon.read, overlapgene=F, codon=T)


vs3at.5358 <- merge(bgmrna.ds3at.acodon.read[, c('mer', 'genenum', 'oddsratio', 'pvalue', 'optimality'),],
                    bgmrna.3at.asite.codon.read[, c('mer', 'genenum', 'oddsratio', 'pvalue'),],
                    by=c('mer'), all=T, suffixes=c('.53', '.58'))
vs3at.5358 <- vs3at.5358[vs3at.5358$mer != '   ',]


{
  pdf('P3_3AT53vs3AT58.pdf', width = 2.4, height = 2.65, useDingbats = F)
  plot(
    ggplot(vs3at.5358, aes(x=log2(oddsratio.58), y=log2(oddsratio.53),
                           color=optimality))+
      geom_point(size=1)+theme_classic()+
      geom_text(data=vs3at.5358[vs3at.5358$mer %in% c('CAT', 'CAC', 'TAG', 'TAA', 'TGA'),],
                aes(label=mer))+
      geom_abline(slope = 1, intercept = 0, linetype='dashed', color='grey')+
      xlim(-5, 6)+ylim(-5, 6)+labs(title='Disome, 3at 58nt vs 53nt')+
      # theme(axis.title = element_text(size=8, color='black'),
      # axis.text = element_text(size=8, color='black'),
      # axis.line = element_line(color='black'))+
      xlab('Log2(A-site pausing score), 58nt, 3-AT')+
      ylab('Log2(A-site pausing score), 53nt, 3-AT')+
      scale_color_manual(values=c('black', '#3399ff'))+
      theme(axis.title = element_text(size=8, color='black'),
            axis.text = element_text(size=8, color='black'),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none')
  )
  dev.off()
}


## the end ####
