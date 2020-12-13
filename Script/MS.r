library(tidyverse)
library(data.table)

hmld.protg.raw <- read.table('./silac/HMLDI_2/proteinGroups.txt', header=T, sep='\t', stringsAsFactors=F)
lmhd.protg.raw <- read.table('./silac/LMHDI_2/proteinGroups.txt', header=T, sep='\t', stringsAsFactors=F)

geno.des <- read.table('~/projects/reference/annotation/s288c_R64-1-1_gene_des_SGD.txt',
                       header=T, sep='\t', stringsAsFactors=F, quote='', comment.char='')

# please don't use ratio h/l normalised
# Ratio H/L normalized -- Normalized ratio between two medium and light label partners.
# The median of ratio sub-populations was shifted to 1.


hmld.protg <- hmld.protg.raw[hmld.protg.raw$Potential.contaminant != '+' & hmld.protg.raw$Reverse != '+', ]
lmhd.protg <- lmhd.protg.raw[lmhd.protg.raw$Potential.contaminant != '+' & lmhd.protg.raw$Reverse != '+', ]

hmld.protg <- hmld.protg[, c('Protein.IDs', 'Majority.protein.IDs', 'Q.value', 'Score', 'Ratio.H.L', 'Intensity',
                             'Intensity.L', 'Intensity.H', 'Number.of.proteins', 'id', 'Peptides')]
lmhd.protg <- lmhd.protg[, c('Protein.IDs', 'Majority.protein.IDs', 'Q.value', 'Score', 'Ratio.H.L', 'Intensity',
                             'Intensity.L', 'Intensity.H', 'Number.of.proteins', 'id', 'Peptides')]

## removing Retrotransposon
## more Retrotransposon in disome

hmld.protg <- hmld.protg[hmld.protg$Number.of.proteins<10,]
lmhd.protg <- lmhd.protg[lmhd.protg$Number.of.proteins<10,]

hmld.protg <- hmld.protg %>% separate_rows(Majority.protein.IDs, sep=';')
lmhd.protg <- lmhd.protg %>% separate_rows(Majority.protein.IDs, sep=';')

hmld.protg$uni <- ifelse(hmld.protg$Number.of.proteins == 1, 1, 0)
lmhd.protg$uni <- ifelse(lmhd.protg$Number.of.proteins == 1, 1, 0)

ratio <- merge(hmld.protg[, c('Majority.protein.IDs', 'Ratio.H.L', 'uni', 'Q.value', 'Score')],
               lmhd.protg[, c('Majority.protein.IDs', 'Ratio.H.L', 'uni', 'Q.value', 'Score')],
               by='Majority.protein.IDs', all=T, suffixes=c('_hmld', '_lmhd'))

ratio$divsmomo_lmhd <- ratio$Ratio.H.L_lmhd
ratio$divsmomo_hmld <- 1/ratio$Ratio.H.L_hmld


hmld.protg$disome <- ifelse(hmld.protg$Majority.protein.IDs %in% unique(disome.peak$gene), 1, 0)
lmhd.protg$disome <- ifelse(lmhd.protg$Majority.protein.IDs %in% unique(disome.peak$gene), 1, 0)

hmld.protg <- merge(hmld.protg, geno.des, by.x='Majority.protein.IDs', by.y='gene', all.x=T)
lmhd.protg <- merge(lmhd.protg, geno.des, by.x='Majority.protein.IDs', by.y='gene', all.x=T)

hmld.protg <- merge(hmld.protg, expr[, c('gene', 'mono.den', 'di.den')],
                    by.x='Majority.protein.IDs', by.y='gene', all.x=T)

lmhd.protg <- merge(lmhd.protg, expr[, c('gene', 'mono.den', 'di.den')],
                    by.x='Majority.protein.IDs', by.y='gene', all.x=T)


lmhd.protg$genefunc <- ifelse(grepl('(chapero)|(hsp)|(Hsp)|(heat)|(fold)', lmhd.protg$des),
                               'chaperone', 'other')
hmld.protg$genefunc <- ifelse(grepl('(chapero)|(hsp)|(Hsp)|(heat)|(fold)', hmld.protg$des),
                              'chaperone', 'other')

lmhd.protg[lmhd.protg$genename == 'PMA1',]$genefunc <- 'other'
hmld.protg[hmld.protg$genename == 'PMA1',]$genefunc <- 'other'

lmhd.protg$genefunc <- ifelse(grepl('(RPL)|(RPS)', lmhd.protg$genename),
                              'ribosomal', lmhd.protg$genefunc)
hmld.protg$genefunc <- ifelse(grepl('(RPL)|(RPS)', hmld.protg$genename),
                              'ribosomal', hmld.protg$genefunc)

write.table(hmld.protg, './silac/hmld.protgroup.txt', quote=F, row.names=F, sep='\t')
write.table(lmhd.protg, './silac/lmhd.protgroup.txt', quote=F, row.names=F, sep='\t')

######
goterm <- read.table('~/projects/reference/annotation/go_slim_mapping.tab',
                     header=F, sep='\t', stringsAsFactors=F)
colnames(goterm) <- c('gene', 'genename', 'sgdid', 'go_aspect', 'go_slim_term', 'goid', 'feature_type')
goterm.p <- goterm[goterm$go_aspect == 'P',]
silac.genes <- unique(c(hmld.protg$Majority.protein.IDs, lmhd.protg$Majority.protein.IDs))
goterm.silac <- goterm.p[goterm.p$gene %in% silac.genes,]
goterm.silac$goid <- gsub('GO:', '', goterm.silac$goid)
goterm.silac <- goterm.silac %>% group_by(gene, genename, sgdid, go_aspect, feature_type) %>%
  summarise(go_slim_term=paste(go_slim_term, collapse=';')) %>% ungroup() %>% as.data.frame()

lmhd.protg.go <- na.omit(lmhd.protg[lmhd.protg$Ratio.H.L > 1,
                                    c('Majority.protein.IDs', 'Ratio.H.L', 'genename', 'des', 'genefunc')])
hmld.protg.go <- na.omit(hmld.protg[hmld.protg$Ratio.H.L < 1,
                                    c('Majority.protein.IDs', 'Ratio.H.L', 'genename', 'des', 'genefunc')])
hmld.protg.go$Ratio.H.L <- 1/hmld.protg.go$Ratio.H.L

colnames(lmhd.protg.go) <- c('gene', 'DIvsMO', 'genename', 'des', 'genefunc')
colnames(hmld.protg.go) <- c('gene', 'DIvsMO', 'genename', 'des', 'genefunc')

lmhd.protg.go <- merge(lmhd.protg.go, goterm.silac[, c('gene', 'go_slim_term')], by='gene', all.x=T)
hmld.protg.go <- merge(hmld.protg.go, goterm.silac[, c('gene', 'go_slim_term')], by='gene', all.x=T)

write.table(lmhd.protg.go, './silac/lmhd.protg.go.txt', quote=F, row.names=F, sep='\t')
write.table(hmld.protg.go, './silac/hmld.protg.go.txt', quote=F, row.names=F, sep='\t')
######

plot.silac.inten <- function(indata, monointen, diinten, incap){
  # require(scales)
  indata[, monointen] <- log2(indata[, monointen])
  indata[, diinten] <- log2(indata[, diinten])
  outplot <- ggplot(indata[indata$genefunc=='other',], aes_string(x=monointen, y=diinten))+
    geom_point(alpha=0.5, color='#999999')+
    geom_point(data=indata[indata$genefunc=='ribosomal',], aes_string(x=monointen, y=diinten),
               alpha=0.7, color='#e57fb0')+#'#0066cc'
    geom_point(data=indata[indata$genefunc=='chaperone',], aes_string(x=monointen, y=diinten),
               alpha=0.7, color='#4ca74a')+#'#FF3300'4ca74a
    geom_point(data=indata[indata$genename=='SSB1',], aes_string(x=monointen, y=diinten),
               alpha=0.7, color='#4ca74a')+ # #00ad45
    geom_point(data=indata[indata$genename=='FAS2',], aes_string(x=monointen, y=diinten),
               alpha=0.7, color='#8f4b99')+ #00ad45
    geom_text(data=indata[indata$genename=='SSB1',], aes_string(x=monointen, y=diinten, label='genename'),
              color='black', hjust=1.2)+
    geom_text(data=indata[indata$genename=='FAS2',], aes_string(x=monointen, y=diinten, label='genename'),
              color='black', hjust=1.2)+
    theme_classic()+
    xlim(16, 32)+ylim(16, 32)+
    geom_abline(slope=1, intercept=0, linetype='dashed',
                color='red', alpha=0.5)+
    # scale_x_continuous(trans=log2_trans())+
    # scale_y_continuous(trans=log2_trans())+
    xlab('log2(Protein intensity in monosome fraction)')+
    ylab('log2(Protein intensity in disome fraction)')+
    theme(legend.position='top',
          axis.title=element_text(size=8, color='black'),
          axis.text=element_text(size=8, color='black'))+
    labs(caption=sprintf('(Heavy label, %s)', incap))
    # scale_color_brewer(palette='Paired', direction=-1)
    # scale_color_manual(name='Gene function',
    #                    values=c('#0066cc', '#999999'))
  return(outplot)
}


{
  pdf('P7_silac_chaperone_nona.pdf', width=4.8, height=2.6, useDingbats=F)
  multiplot(plot.silac.inten(lmhd.protg[!is.na(lmhd.protg$Ratio.H.L),], 'Intensity.L', 'Intensity.H', 'disome'),
            plot.silac.inten(hmld.protg[!is.na(hmld.protg$Ratio.H.L),], 'Intensity.H', 'Intensity.L', 'monosome'),
            cols=2)
  dev.off()
}

# RQT4 / YKR023W
# HEL2 / YDR266C
# SLH1 / YGR271W
# RQT4 / YKR023W


lmhd.protg[lmhd.protg$Majority.protein.IDs %in% c('YKR023W', 'YDR266C', 'YGR271W', 'YKR023W'),]
hmld.protg[hmld.protg$Majority.protein.IDs %in% c('YKR023W', 'YDR266C', 'YGR271W', 'YKR023W'),]

##mh.test####
silac.com <- hmld.protg[, c('Majority.protein.IDs', 'Ratio.H.L', 'genefunc')]
silac.com$Ratio.H.L <- 1/silac.com$Ratio.H.L
# silac.com <- rbind(lmhd.protg[, c('Majority.protein.IDs', 'Ratio.H.L', 'genefunc')],
#                    silac.com)
silac.com <- merge(silac.com,
                   lmhd.protg[, c('Majority.protein.IDs', 'Ratio.H.L', 'genefunc')],
                   by='Majority.protein.IDs', all=T, suffixes=c('.hmld', '.lmhd'))
silac.com <- na.omit(silac.com)
write.table(silac.com, 'silac/silac.com.txt', sep='\t', quote=F, row.names=F)



lmhd.pep <- read.table('./silac/LMHDI_2/peptides.txt', header=T, sep='\t', stringsAsFactors=F)
lmhd.pep.posi <- lmhd.pep[lmhd.pep$Reverse != '+' & lmhd.pep$Potential.contaminant != '+',
                          c('Proteins', 'Leading.razor.protein', 'Start.position', 'End.position', 'Sequence')]
lmhd.pep.posi$cds.start <- (lmhd.pep.posi$Start.position-1)*3
lmhd.pep.posi$cds.end <- lmhd.pep.posi$End.position*3-1


cds.kr <- sapply(cds.pep$pep, function(x){
  which(strsplit(x, '')[[1]] %in% c('R', 'K'))
}, USE.NAMES=F)

cds.kr <- data.frame(gene=rep(cds.pep$gene, lapply(cds.kr, length) %>% unlist(use.names=F)),
                     aaposi=unlist(cds.kr))
cds.kr$gene <- as.character(cds.kr$gene)
cds.kr <- rbind(cds.kr, data.frame(gene=unique(cds.kr$gene),
                                   aaposi=rep(0, length(unique(cds.kr$gene)))))
cds.kr <- cds.kr[order(cds.kr$gene, cds.kr$aaposi, decreasing=F), ]
rownames(cds.kr) <- NULL

cds.kr$gene.start <- c('', cds.kr$gene[1:length(cds.kr$gene)-1])
cds.kr$aa.start <- c(NaN, cds.kr$aaposi[1:length(cds.kr$aaposi)-1])
cds.kr <- cds.kr[-1, c('gene.start', 'aa.start', 'gene', 'aaposi')]
cds.kr <- merge(cds.kr, cds.pep[, c('gene', 'len')],
                by.x='gene.start', by.y='gene')
cds.kr$len <- cds.kr$len/3-1
cds.kr$aaposi <- ifelse(cds.kr$gene.start != cds.kr$gene,
                        cds.kr$len, cds.kr$aaposi)
cds.kr[, c('gene', 'len')] <- NULL
cds.kr$aa.start <- cds.kr$aa.start+1
colnames(cds.kr) <- c('gene', 'aa.start', 'aa.end')
cds.kr$aa.len <- cds.kr$aa.end-cds.kr$aa.start+1
cds.kr$nt.start <- (cds.kr$aa.start-1)*3



lmhd.pep.cds <- lmhd.pep.posi[, c('Leading.razor.protein', 'cds.start')]
colnames(lmhd.pep.cds) <- c('gene', 'nt.start')
lmhd.pep.cds$tag <- 'nascent'
lmhd.pep.cds <- rbind(lmhd.pep.cds,
                      cbind(cds.kr[cds.kr$gene %in% unique(lmhd.pep.cds$gene), c('gene', 'nt.start')],
                            tag='cds'))

pep.odds <- lmhd.pep.cds
pep.odds <- merge(pep.odds, cds.pep[, c('gene', 'len')], by='gene')
pep.odds$len <- pep.odds$len-3
pep.odds$region <- ifelse(pep.odds$nt.start < pep.odds$len/2, 5, 3)
pep.odds$region <- paste0('p', pep.odds$region)
pep.odds <- pep.odds %>% group_by(gene, region, tag) %>%
  summarise(sites=length(nt.start)) %>% ungroup() %>%
  as.data.frame()
pep.odds$grp <- paste(pep.odds$region, pep.odds$tag, sep='.')
pep.odds <- pep.odds[, c('gene', 'sites', 'grp')] %>%
  spread(grp, sites) %>% ungroup() %>% as.data.frame()
pep.odds[is.na(pep.odds)] <- 0

pep.odds[, c('odds.ratio', 'p.value')] <- mapply(function(w,x,y,z){
  out <- fisher.test(matrix(c(w,x,y,z), byrow=T, nrow=2))
  return(c(out$estimate, out$p.value))
}, pep.odds$p5.nascent, pep.odds$p5.cds,
pep.odds$p3.nascent, pep.odds$p3.cds) %>% t()

pep.odds$or.type <- ifelse(pep.odds$odds.ratio>=1, 'en', 'dep')
pep.odds <- merge(pep.odds, goterm.silac[, c('gene', 'go_slim_term')],
                  by='gene', all.x=T)
pep.odds <- merge(pep.odds, geno.des[, c('gene', 'genename')],
                  by='gene', all.x=T)

pep.odds$metagrp <- ifelse(grepl('metabolic', pep.odds$go_slim_term), 'metabolic', 'other')
pep.odds$chapgrp <- ifelse(pep.odds$gene %in% unique(lmhd.protg[lmhd.protg$genefunc == 'chaperone',]$Majority.protein.IDs),
                           'chaperone', 'other')
pep.odds$ribogrp <- ifelse(grepl('(RPL)|(RPS)', pep.odds$genename), 'ribosomal', 'other')


pep.conti <- as.vector(t(pep.odds[, c('p5.nascent', 'p3.nascent', 'p5.cds', 'p3.cds')]))
pep.conti <- array(pep.conti, dim=c(2, 2, nrow(pep.odds)))
pepmg <- mantelhaen.test(pep.conti, exact=F)

pep.conti.meta <- as.vector(t(pep.odds[pep.odds$metagrp == 'metabolic', c('p5.nascent', 'p3.nascent', 'p5.cds', 'p3.cds')]))
pep.conti.meta <- array(pep.conti.meta, dim=c(2, 2, nrow(pep.odds[pep.odds$metagrp == 'metabolic',])))
pepmg.meta <- mantelhaen.test(pep.conti.meta, exact=F)

pep.conti.chap <- as.vector(t(pep.odds[pep.odds$chapgrp == 'chaperone', c('p5.nascent', 'p3.nascent', 'p5.cds', 'p3.cds')]))
pep.conti.chap <- array(pep.conti.chap, dim=c(2, 2, nrow(pep.odds[pep.odds$chapgrp == 'chaperone',])))
pepmg.chap <- mantelhaen.test(pep.conti.chap, exact=F)

pep.conti.ribo <- as.vector(t(pep.odds[pep.odds$ribogrp == 'ribosomal', c('p5.nascent', 'p3.nascent', 'p5.cds', 'p3.cds')]))
pep.conti.ribo <- array(pep.conti.ribo, dim=c(2, 2, nrow(pep.odds[pep.odds$ribogrp == 'ribosomal',])))
pepmg.ribo <- mantelhaen.test(pep.conti.ribo, exact=F)


############
head(expr)
lmhd.protg

silactest <- lmhd.protg
expr.silac <- expr

mseq <- merge(expr.silac[, c('gene', 'di.den', 'mono.den', 'rpm.mono', 'rpm.di')],
              silactest[!is.na(silactest$Ratio.H.L),
                        c('Majority.protein.IDs', 'Intensity.H', 
                          'Intensity.L', 'tag', 'genefunc')],
              by.y='Majority.protein.IDs',
              by.x='gene')

goterm.silac.grp <- goterm.silac
goterm.silac.grp$metabolic <- ifelse(grepl('metabolic', goterm.silac.grp$go_slim_term),
                                     'metabolic', 'other')

mseq$gotag <- mseq$genefunc
mseq$gotag <- ifelse(mseq$gene %in% goterm.silac.grp[goterm.silac.grp$metabolic == 'metabolic', ]$gene,
                     'metabolic', mseq$gotag)
# mseq$gotag <- ifelse(mseq$gene %in% silactest[silactest$genefunc == 'ribosomal', ]$Majority.protein.IDs,
#                      'ribosomal', mseq$gotag)
mseq$gotag <- ifelse(mseq$gene %in% mseq[mseq$genefunc %in% c('ribosomal', 'chaperone'), ]$gene,
                     'ribosomal', mseq$gotag)


library(lmodel2)
ribocor <- cor.test(log2(mseq[mseq$gotag == 'ribosomal',]$Intensity.H),
                   log2(mseq[mseq$gotag == 'ribosomal',]$rpm.di), method = 's')
riboma <- lmodel2(log2(rpm.di)~log2(Intensity.H),
                  mseq[mseq$gotag == 'ribosomal',])

othercor <- cor.test(log2(mseq[mseq$gotag == 'metabolic',]$Intensity.H),
                     log2(mseq[mseq$gotag == 'metabolic',]$rpm.di), method = 's')
otherma <- lmodel2(log2(rpm.di)~log2(Intensity.H),
                   mseq[mseq$gotag == 'metabolic',])


{
  pdf('P7_Ms_vs_seq.pdf', width = 5.5, height = 2.7, useDingbats = F)
  multiplot(
    ggplot(mseq, aes(x=log2(Intensity.L), y=log2(Intensity.H), color=gotag))+
      geom_point()+theme_classic()+
      xlim(16, 32)+ylim(16, 32)+
      scale_color_manual(values = c('#3591aa', '#999999', '#e57fb0'))+
      theme(axis.title = element_text(color='black', size=8),
            axis.text = element_text(color='black', size=8),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none')+
      geom_abline(slope = 1, intercept = 0, color='grey', linetype='dashed'),
    ggplot(mseq[mseq$gotag != 'other', ], 
           aes(x=log2(Intensity.H), y=log2(rpm.di), color=gotag))+
      geom_point()+theme_classic()+
      scale_color_manual(values = c('#3591aa', '#e57fb0'))+
      geom_abline(intercept = riboma$regression.results[2,2],
                  slope = riboma$regression.results[2,3], color='#e57fb0')+
      geom_abline(intercept = otherma$regression.results[2,2],
                  slope = otherma$regression.results[2,3], color='#3591aa')+
      annotate('text', x=23, y=4, color='#e57fb0',
               label=sprintf('rho=%.2f P=%.2f N = 140', ribocor$estimate, ribocor$p.value))+
      annotate('text', x=23, y=3, color='#3591aa',
               label=sprintf('rho=%.2f P=%.2e N = 25', othercor$estimate, othercor$p.value))+
      theme(axis.title = element_text(color='black', size=8),
            axis.text = element_text(color='black', size=8),
            axis.line = element_line(color='black'),
            axis.ticks = element_line(color='black'),
            legend.position = 'none'), cols = 2
  )
  dev.off()
}


## the end ####
