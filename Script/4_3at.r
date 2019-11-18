library(tidyverse)
library(parallel)
library(grid)
suppressMessages(library(IRanges))

correct.pos <- function(indata){
  outpos <- (1-indata$frame5)%%3-1 + indata$dis5
  return(outpos)
}

getseq <- function(starts, ends, ingene, aa){
  if (aa) seqcol <- 'pep' else seqcol <- 'cds'
  geneseq <- cds.pep[cds.pep$gene == ingene, seqcol]
  blk <- ''
  blke <- ''
  seqlen <- ends - starts + 1
  if (starts<=0) {
    if (ends<=0) {
      starts <- 0
      ends <- 0
      blk <- strrep(' ', seqlen)
    }
    else {
      blk <- strrep(' ', abs(starts)+1)
      starts <- 1
    }
  }
  if (ends > nchar(geneseq)){
    genelen <- nchar(geneseq)
    # print(geneseq)
    if (starts > genelen){
      starts <- 0
      ends <- 0
      blke <- strrep(' ', seqlen)
    }
    else {
      blke <- strrep(' ', ends-genelen)
      ends <- genelen
    }
  }
  
  return(paste0(blk, substr(geneseq, starts, ends), blke))
}

kmers <- function(indata, incol, mer, aa, count_read, readcol='read'){
  # mer for codon or amino acid, not nucleotide
  if (aa) seqstep <- 1 else seqstep <- 3 
  countlist <- mcmapply(function(x){
    inseqlen <- nchar(x)-(mer*seqstep-1)
    mers <- substring(x,
                      seq(1, inseqlen, seqstep),
                      seq(1+mer*seqstep-1, nchar(x), seqstep))
    # print(mers)
    return(mers)
  },as.character(indata[, incol]), SIMPLIFY=F, mc.cores=4L)
  names(countlist) <- NULL
  countlist <- do.call(rbind, countlist)
  countlist <- as.data.frame(countlist)
  posinum <- paste0('p', 1:ncol(countlist))
  colnames(countlist) <- posinum
  if (count_read){
    countlist <- cbind(indata[, c('gene', readcol)], countlist)
  }
  else {
    countlist <- cbind(data.frame(gene=indata$gene, read=1), countlist)
  }
  # print(countlist)
  countlist <- countlist %>% gather_('posi', 'mer', posinum)
  return(countlist)
}

peak.sum <- function(indata){
  indata$cor.pos <- correct.pos(indata)
  indata$asite <- get.asite(indata)
  indata <- indata %>% group_by(rep, gene, asite) %>%
    summarise(read=sum(seqsnum)) %>% ungroup()
  indata$rep <- paste0('r', indata$rep)
  repcols <- unique(indata$rep)
  indata <- indata %>% spread(rep, read)
  indata[is.na(indata)] <- 0
  indata$read <- rowSums(indata[, repcols])
  return(indata)
}

peak.rpm <- function(inlib){
  librpm <- inlib / (sum(inlib) * 1e-6)
  return(librpm)
}

rpm.mean <- function(indata, inrow){
  outdata <- mapply(function(x) peak.rpm(x),
                    indata[, inrow])
  # print(head(outdata))
  # print(ncol(outdata))
  outdata <- rowMeans(outdata)
  return(outdata)
}

## disome-3at ####
di.3at.ori <- read.table('di_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
di.3at.raw <- di.3at.ori[di.3at.ori$seqlen %in% c(58, 59, 61, 62), ]
di.3at.read <- data.read(di.3at.ori)
plot.read.cor(di.3at.read, 'r1', 'r2', '3at.disome')
rm(di.3at.ori)
di.3at.raw[, c('sequnum', 'frame3', 'putative', 'region', 'dupl')] <- NULL
rownames(di.3at.raw) <- NULL
di.3at.raw$cor.pos <- correct.pos(di.3at.raw)
di.3at.raw$lengroup <- di.3at.raw$seqlen %/% 3 *3

d3at.site.sum <- di.3at.raw %>% group_by(gene, lengroup, cor.pos, rep) %>%
  summarise(read=sum(seqsnum)) %>% ungroup()
d3at.site.sum$rep <- paste0('r', d3at.site.sum$rep)
d3at.site.sum <- d3at.site.sum %>% spread(rep, read)
d3at.site.sum[is.na(d3at.site.sum)] <- 0
d3at.site.sum$read <- rowSums(d3at.site.sum[, c('r1', 'r2')])
d3at.site.sum <- d3at.site.sum %>% as.data.frame()
d3at.site.sum$codonseq <- mapply(function(st, ed, ing, seqtype) getseq(st, ed, ing, seqtype),
                              d3at.site.sum$cor.pos+1, d3at.site.sum$cor.pos+d3at.site.sum$lengroup+6,
                              d3at.site.sum$gene, F)
d3at.site.sum$aaseq <- mapply(function(st, seqlen, ing, seqtype) {
  ed <- (st + seqlen+6)/3
  st <- st / 3 + 1
  getseq(st, ed, ing, seqtype)
}, d3at.site.sum$cor.pos, d3at.site.sum$lengroup, d3at.site.sum$gene, T)
d3at.site.sum$rpm <- rpm.mean(d3at.site.sum, c('r1', 'r2'))

## monosome-3at ####
mono.3at.ori <- read.table('mono_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.3at.raw <- mono.3at.ori[mono.3at.ori$seqlen %in% c(28,29), ]
mono.3at.read <- data.read(mono.3at.raw)
rm(mono.3at.ori)
mono.3at.raw[, c('sequnum', 'frame3', 'putative', 'region', 'dupl')] <- NULL
rownames(mono.3at.raw) <- NULL
mono.3at.raw$cor.pos <- correct.pos(mono.3at.raw)
mono.3at.raw$lengroup <- mono.3at.raw$seqlen %/% 3 *3
m3at.site.sum <- mono.3at.raw %>% group_by(gene, lengroup, cor.pos, rep) %>%
  summarise(read=sum(seqsnum)) %>% ungroup()
m3at.site.sum$rep <- paste0('r', m3at.site.sum$rep)
m3at.site.sum <- m3at.site.sum %>% spread(rep, read)
m3at.site.sum[is.na(m3at.site.sum)] <- 0
m3at.site.sum$read <- rowSums(m3at.site.sum[, c('r0', 'r1')])
m3at.site.sum <- m3at.site.sum %>% as.data.frame()
m3at.site.sum$codonseq <- mapply(function(st, ed, ing, seqtype) getseq(st, ed, ing, seqtype),
                                 m3at.site.sum$cor.pos+1, m3at.site.sum$cor.pos+m3at.site.sum$lengroup+6,
                                 m3at.site.sum$gene, F)
m3at.site.sum$aaseq <- mapply(function(st, seqlen, ing, seqtype) {
  ed <- (st + seqlen+6)/3
  st <- st / 3 + 1
  getseq(st, ed, ing, seqtype)
}, m3at.site.sum$cor.pos, m3at.site.sum$lengroup, m3at.site.sum$gene, T)
m3at.site.sum$rpm <- rpm.mean(m3at.site.sum, c('r0', 'r1'))
######

kmers.frame <- function(inframe, indata, mer=1, aa=F, count_read=T, readcol='rpm'){
  incol <- 'codonseq'
  # print(head(indata))
  indata[, incol] <- substring(indata[, incol], 1+inframe, nchar(indata[, incol]))
  # print(head(indata))
  return(
    kmers(indata, 'codonseq', mer, aa=aa, count_read=count_read, readcol=readcol)
  )
}

######
# A different way to calculate the A-site of the 3at-disome data
# The A-sites of both 58~59 nt and 61~62 nt are at the 16th codon (start from 0)
head(d3at.site.sum)
d3at.site.sum$asite <- d3at.site.sum$cor.pos + 15*3
d3at.site.sum$rpm <- rpm.mean(d3at.site.sum, c('r1', 'r2'))

d3at.sites <- d3at.site.sum %>% group_by(gene, asite) %>%
  summarise(rpm=sum(rpm)) %>% ungroup()
d3at.sites <- d3at.sites %>% group_by(gene) %>%
  mutate(zscore=(rpm-mean(rpm))/sd(rpm)) %>%
  as.data.frame()

## removing genes with only one read / peak
# d3at.sites <- na.omit(d3at.sites)

disome.sites <- disome.peak[, c('gene', 'asite', 'rpm')]
disome.sites <- disome.sites %>% group_by(gene) %>%
  mutate(zscore=(rpm-mean(rpm))/sd(rpm)) %>%
  as.data.frame()

mono.sites <- mono.peak[, c('gene', 'asite', 'rpm')]
mono.sites <- mono.sites %>% group_by(gene) %>%
  mutate(zscore=(rpm-mean(rpm))/sd(rpm)) %>%
  as.data.frame()

cds.len <- read.table('~/projects/reference/annotation/fasta_cds_pep.txt', header=T, sep='\t',
                      colClasses=c('character', rep('NULL', 2), 'numeric'))

multiplots <- function(datalist, coln){
  plots <- length(datalist)
  rown <- ceiling(plots / coln)
  rowf <- 1/rown
  colf <- 1/coln
  # print(paste(plots, rown, rowf, colf))
  for (i in 1:plots){
    plotrow <- ceiling(i / coln) - 1
    plotcol <- (i-1) %% coln
    print(paste(i, plotrow, plotcol))
    print(paste((colf/2+(plotcol)*colf), (rowf/2+(plotrow)*rowf)))
    print(
      datalist[[i]],
      vp=viewport(1/coln-0.01, 1/rown-0.01,
                  x=(colf/2+(plotcol)*colf), y=(rowf/2+(plotrow)*rowf))
    )
  }
}

######
mono.3at.ori <- read.table('mono_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.3at.peri <- read.peri(mono.3at.ori) %>% group_by(gene) %>%
  mutate(rpm.nor=rpm/sum(rpm)) %>% ungroup()
rm(mono.3at.ori)

di.3at.ori <- read.table('di_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
di.3at.peri <- read.peri(di.3at.ori) %>% group_by(gene) %>%
  mutate(rpm.nor=rpm/sum(rpm)) %>% ungroup()
rm(di.3at.ori)

cds.cdn.posi <- function(codon, aa){
  cdnmer <- function(x){return(substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3)))}
  # aamer <- function(x){return(strsplit(x, '')[[1]])}
  # aacol <- 'pep'
  count.codon <- function(mer){
    outdata <- sapply(cds.pep[, 'cds'], function(x) {which(cdnmer(x) == mer)}, USE.NAMES=F)
    outdata <- data.frame(cdnposi=unlist(outdata),
                          gene=rep(cds.pep$gene, lapply(outdata, length) %>% unlist(use.names=F)))
    outdata$cdnposi <- (outdata$cdnposi-1)*3
    outdata$codon <- mer
    return(outdata)
  }
  outdata <- do.call(rbind, lapply(codon, function(x) count.codon(x)))
  outdata$aa <- aa
  return(outdata)
}

cds.his <- cds.cdn.posi(c('CAT', 'CAC'), 'H')
# cds.ala <- cds.cdn.posi(c('GCT', 'GCC', 'GCA', 'GCG'), 'A')
# cds.leu <- cds.cdn.posi(c('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'), 'L')

around.aa <- function(aaposi, dis5posi, aarange){
  colnames(aaposi)[1] <- 'cdnposi'
  require(IRanges)
  
  aaposi <- with(aaposi,
                  RangedData(IRanges(start=cdnposi-aarange, end=cdnposi+aarange),
                             space=gene, cdnposi=cdnposi, codon=codon, aa=aa))
  dis5posi <- with(dis5posi,
                    RangedData(IRanges(start=dis5, end=dis5),
                               space=gene, dis5=dis5, rpm=rpm, rpm.nor=rpm.nor))
  fo <- findOverlaps(aaposi, dis5posi)
  mfo <- as.matrix(fo)
  aaId <- mfo[, 1]
  dis5Id <- mfo[, 2]
  
  outdata <- data.frame(gene=space(aaposi)[aaId], codon=aaposi$codon[aaId], aa=aaposi$aa[aaId],
                        posiaa=aaposi$cdnposi[aaId], dis5=dis5posi$dis5[dis5Id],
                        rpm=dis5posi$rpm[dis5Id], rpm.nor=dis5posi$rpm.nor[dis5Id])
  
  outdata$diff <- outdata$dis5 - outdata$posiaa
  return(outdata)
}

di.3at.his <- around.aa(cds.his, di.3at.peri[di.3at.peri$seqlen %in% c(58, 59, 61, 62),], 120)
mono.3at.his <- around.aa(cds.his, mono.3at.peri[mono.3at.peri$seqlen %in% c(28, 29), ], 120)

# di.3at.ala <- around.aa(cds.ala, di.3at.peri[di.3at.peri$seqlen %in% c(58, 59, 61, 62),], 120)
# mono.3at.ala <- around.aa(cds.ala, mono.3at.peri[mono.3at.peri$seqlen %in% c(28, 29), ], 120)

# di.3at.leu <- around.aa(cds.leu, di.3at.peri[di.3at.peri$seqlen %in% c(58, 59, 61, 62),], 120)
# mono.3at.leu <- around.aa(cds.leu, mono.3at.peri[mono.3at.peri$seqlen %in% c(28, 29), ], 120)

plot.disttoaa <- function(indata, datatype, xls=-120, xle=120, inweight='rpm.nor', incolor=NULL){
  # require(scales)
  require(RColorBrewer)
  return(
    ggplot(indata, aes_string(x='diff', weight=inweight, color=incolor))+
      geom_freqpoly(binwidth=1, alpha=0.4, size=0.5)+theme_classic()+
      xlab('Distance from Codon (nt)')+
      ylab('Normalised RPF')+
      scale_x_continuous(breaks=seq(xls, xle, 30), limits=c(xls, xle))+
      # scale_color_brewer(palette='Greens', direction=1)+
      scale_color_manual(values=brewer.pal(5, 'Set1'))+
      # scale_y_continuous(trans=log2_trans())+
      theme(axis.title=element_text(size=8),
            axis.text=element_text(size=8),
            legend.title=element_text(size=8),
            legend.text=element_text(size=8),
            legend.position='top')+
      annotate('text', label=datatype, x=-Inf, y=Inf, vjust=1.5, hjust=-0.1)
  )
}

{
  pdf('P2_his_distance_-60to30_rpm.pdf', width=5, height=4)
  multiplot(plot.disttoaa(di.3at.his, xls=-60, xle=30, datatype='Disome', incolor='codon', inweight='rpm'),
            plot.disttoaa(mono.3at.his, xls=-60, xle=30, datatype='Monosome', incolor='codon', inweight='rpm'),
            plot.disttoaa(di.3at.ala, xls=-60, xle=30, datatype='Disome', incolor='codon', inweight='rpm'),
            plot.disttoaa(mono.3at.ala, xls=-60, xle=30, datatype='Monosome', incolor='codon', inweight='rpm'),
            cols=2, byrow=F)
  dev.off()
}

## replication ####
mrna.3at.ori <- read.table('mRNA_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mrna.3at.read <- data.read(mrna.3at.ori)
rm(mrna.3at.ori)

mono.3at.ori <- read.table('mono_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.3at.read <- data.read(mono.3at.ori)
rm(mono.3at.ori)

di.3at.ori <- read.table('di_3at_ori.txt', header=T, sep='\t', stringsAsFactors=F)
di.3at.read <- data.read(di.3at.ori)
rm(di.3at.ori)

{
  pdf('P0_replication_3at.pdf', width=13, height=4, useDingbats=F)
  multiplot(plot.read.cor(di.3at.read, 'r1', 'r2', 'Disome'),
            plot.read.cor(mono.3at.read, 'r0', 'r1', 'Monosome'),
            plot.read.cor(mrna.3at.read, 'r0', 'r1', 'mRNA'),
            cols=3)
  dev.off()
}

di.mono.3at.cor <- merge(di.3at.read[, c('gene', 'rpm')],
                         mono.3at.read[, c('gene', 'rpm')],
                     by='gene', suffixes=c('.di', '.mono'))

{
  pdf('P0_replication_monovsdi_3at.pdf', width=5.3, height=5, useDingbats=FALSE)
  plot(plot.read.cor(di.mono.3at.cor, 'rpm.di', 'rpm.mono', 'mix'))
  dev.off()
}

## the end ####





























