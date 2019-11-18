library(tidyverse)
library(parallel)

## functions ####
# functions from p1.r
correct.pos <- function(indata){
  outpos <- (1-indata$frame5)%%3-1 + indata$dis5
  return(outpos)
}

get.asite <- function(indata){
  asite <- (indata$seqlen %/% 3)*3 + indata$cor.pos - 12
  return(asite)
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

peak.seq <- function(indata, part, aanum, aa=T, asite=F, seqafter=F){
  if (!(part %in% c('exit', 'read'))) stop('Wrong parts')
  if (seqafter){
    indata$asite <- indata$asite+(-3*as.numeric(!asite)+1) # including the asite position
    if (aa) indata$asite <- (indata$asite+2)/3 else aanum <- aanum*3
    indata$start <- indata$asite + aanum - 1
    out <- mapply(function(x, y, z, seqtype) getseq(x, y, z, seqtype),
                  indata$asite, indata$start, indata$gene, aa)
  }else{
    if (asite) indata$asite <- indata$asite+3 # including the asite position
    if (aa) indata$asite <- indata$asite/3 else aanum <- aanum*3
    indata$start <- indata$asite - aanum + 1
    out <- mapply(function(x, y, z, seqtype) getseq(x, y, z, seqtype),
                  indata$start, indata$asite, indata$gene, aa)
  }
  return(out)
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

count_kmer <- function(indata, incol, mer, aa=T, count_read=T){
  # mer for codon or amino acid, not nucleotide
  countlist <- kmers(indata, incol, mer, aa=aa, count_read=count_read)
  countlist <- countlist %>% group_by(gene, mer) %>%
    summarise(num=sum(read)) %>% ungroup()
  countlist <- countlist %>% group_by(gene) %>%
    mutate(notnum=sum(num)-num) %>% ungroup()
  return(countlist)
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

enrich <- function(indata, bgdata, overlapgene=F, mccores=4L, allcom=F, aa=F){
  start_time <- Sys.time()
  if (overlapgene){
    bgdata <- bgdata[bgdata$gene %in% unique(indata$gene),]
  }
  
  if (allcodons){
    inmer <- ifelse(codon.tb$Codon[1:64] )
  }
  
  entable <- merge(indata, bgdata, all=T, by=c('gene', 'mer'), suffixes=c('.in', '.bg'))
  entable[is.na(entable)] <- 0
  ensum <- entable %>% mutate(innum=num.in+notnum.in, bgnum=num.bg+notnum.bg)
  ensum <- merge(unique(ensum[ensum$innum !=0 , c('gene', 'innum')]),
                 unique(ensum[ensum$bgnum !=0 , c('gene', 'bgnum')]), by='gene', all=T)
  ensum[is.na(ensum)] <- 0
  entable <- merge(entable, ensum, by=c('gene'), all.x=T)
  entable$notnum.in <- ifelse(entable$notnum.in == 0 | entable$innum == 0, entable$innum, entable$notnum.in)
  entable$notnum.bg <- ifelse(entable$notnum.bg == 0 | entable$bgnum == 0, entable$bgnum, entable$notnum.bg)
  entable <- entable[, -c(7:8)]
  # return(entable)
  outdata <- entable %>% group_by(mer) %>%
    summarise(genenum=length(mer))
  print(head(outdata))
  outdata <- cbind(outdata, mcmapply(function(x) {
    mer.mhtest(x, entable)
  }, outdata$mer) %>% t() %>% as.data.frame(), mc.cores=mccores)
  rownames(outdata) <- NULL
  colnames(outdata)[3:5] <- c('oddsratio', 'pvalue', 'method')
  # print(head(outdata))
  outdata[c('oddsratio', 'pvalue')] <- sapply(outdata[c('oddsratio', 'pvalue')], 
                                   function(x) as.character(x) %>% as.numeric(x))
  print(Sys.time()- start_time)
  return(outdata)
}

mer.mhtest <- function(mer, indata){
  indata <- indata[indata$mer == mer, ]
  if (nrow(indata) < 2){
    # print(mer)
    out <- indata[, mhcols] %>% unlist() %>% matrix(nrow=2) %>% fisher.test()
    return(c(out$estimate, out$p.value, 'fisher'))
  }
  else{
    conti.table <- as.vector(t(indata[, mhcols]))
    conti.table <- array(conti.table, dim=c(2, 2, nrow(indata)))
    # print(c(mer, r))
    out <- mantelhaen.test(conti.table, exact=F)
    return(c(as.numeric(out$estimate), as.numeric(out$p.value), 'mantelhaen'))
  }
}
mhcols <- c('num.in', 'notnum.in', 'num.bg', 'notnum.bg')

## read data####
cds.pep <- read.table('~/projects/reference/annotation/fasta_cds_pep.txt',
                      header=T, sep='\t', stringsAsFactors=F)
cds.pep$pep <- gsub('O', '*', cds.pep$pep)

disome.peak.kmer <- read.table('./disome.peak.txt',header=T, sep='\t', stringsAsFactors=F,
                               colClasses=c('character', rep('numeric', 6), rep('NULL', 2)))
disome.peak.kmer <- disome.peak.kmer[disome.peak.kmer$asite >= 0,]
disome.peak.kmer$read <- rowSums(disome.peak.kmer[, c('r1', 'r2')])
disome.peak.kmer$aaseq <- peak.seq(disome.peak.kmer, 'exit', 21, asite=F)
disome.peak.kmer$codonseq <- peak.seq(disome.peak.kmer, 'exit', 21, aa=F)

### monosome
mono.peak <- read.table('./mono.peak.txt',header=T, sep='\t', stringsAsFactors=F,
                               colClasses=c('character', rep('numeric', 6), rep('NULL', 2)))
mono.peak.kmer <- mono.peak[mono.peak$asite >= 0,]
# mono.peak.kmer$read <- rowSums(mono.peak.kmer[, c('r1', 'r2', 'r0')])
mono.peak.kmer$aaseq <- peak.seq(mono.peak.kmer, 'exit', 21, asite=F)
mono.peak.kmer$codonseq <- peak.seq(mono.peak.kmer, 'exit', 21, aa=F)
###

## mRNA ##
mrna.ori <- read.table('mRNA_ori.txt', header=T, sep='\t',
                       colClasses=c('character', rep('numeric', 6), rep('NULL', 4), 'numeric', 'character'))
mrna.data <- mrna.ori[mrna.ori$seqlen %in% 28:29, ]
rm(mrna.ori)
mrna.peak <- peak.sum(mrna.data)
mrna.peak <- merge(mrna.peak, cds.pep[, c('gene', 'len')], all.x=T, by='gene')
mrna.peak$sc <- ifelse(mrna.peak$len-3 == mrna.peak$asite, 1, 0)
mrna.peak$sc[mrna.peak$asite > mrna.peak$len-3] <- -1
mrna.peak$len <- NULL
mrna.peak$rpm <- rpm.mean(mrna.peak, c('r0', 'r1', 'r2'))
# mrna.peak.kmer <- mrna.peak[mrna.peak$sc == 0 & mrna.peak$asite >= 0,]
mrna.peak.kmer <- mrna.peak[mrna.peak$asite >= 0,]
mrna.peak.kmer$aaseq <- peak.seq(mrna.peak.kmer, 'exit', 21)
mrna.peak.kmer$codonseq <- peak.seq(mrna.peak.kmer, 'exit', 21, aa=F)
mrna.peak.kmer$rpm <- rpm.mean(mrna.peak.kmer[, c('r0', 'r1', 'r2')])

{
  start_time <- Sys.time()
  di3mer.codon.rmsc <- count_kmer(disome.peak.kmer[disome.peak.kmer$sc == 0, ], 'codonseq', 3, aa=F, count_read=F)
  di3mer.aa.rmsc <- count_kmer(disome.peak.kmer[disome.peak.kmer$sc == 0, ], 'aaseq', 3, aa=T, count_read=F)
  mono3mer.codon.rmsc <- count_kmer(mono.peak.kmer[mono.peak.kmer$sc == 0, ], 'codonseq', 3, aa=F, count_read=F)
  mono3mer.aa.rmsc <- count_kmer(mono.peak.kmer[mono.peak.kmer$sc == 0, ], 'aaseq', 3, aa=T, count_read=F)
  mrna3mer.codon.rmsc <- count_kmer(mrna.peak.kmer[mrna.peak.kmer$sc == 0, ], 'codonseq', 3, aa=F, count_read=F)
  mrna3mer.aa.rmsc <- count_kmer(mrna.peak.kmer[mrna.peak.kmer$sc == 0, ], 'aaseq', 3, aa=T, count_read=F)
  
  bgmrna.3mer.aa.rmsc <- enrich(di3mer.aa.rmsc, mrna3mer.aa.rmsc, overlapgene=F)
  bgmrna.3mer.codon.rmsc <- enrich(di3mer.codon.rmsc, mrna3mer.codon.rmsc, overlapgene=F)
  monobgmrna.3mer.aa.rmsc <- enrich(mono3mer.aa.rmsc, mrna3mer.aa.rmsc, overlapgene=F)
  monobgmrna.3mer.codon.rmsc <- enrich(mono3mer.codon.rmsc, mrna3mer.codon.rmsc, overlapgene=F)
}


codon.tb <- read.table('~/codons.txt', header=T, sep='\t', stringsAsFactors=F)
codon.tb[codon.tb$Letter == 'O',]$Letter <- '*'
codon.tb[codon.tb$Codon == 'TAA',]$optimality <- '+'

kmer_read <- function(indata, incol, mer, aa=T){
  # mer for codon or amino acid, not nucleotide
  if (aa) seqstep <- 1 else seqstep <- 3
  countlist <- mcmapply(function(x, y){
    inseqlen <- nchar(x)-(mer*seqstep-1)
    mers <- substring(x,
                      seq(1, inseqlen, seqstep),
                      seq(1+mer*seqstep-1, nchar(x), seqstep))
    # print(mers)
    return(mers)
  }, indata[, incol], SIMPLIFY=F, mc.cores=4L)
  names(countlist) <- NULL
  countlist <- do.call(rbind, countlist)
  countlist <- as.data.frame(countlist)
  posinum <- paste0('p', 1:ncol(countlist))
  colnames(countlist) <- posinum
  # print(head(indata))
  countlist <- cbind(indata[, c('gene', 'read', 'peakshare')], countlist)
  countlist[is.na(countlist)] <- 0
  # print(countlist)
  countlist <- countlist %>% gather_('posi', 'mer', posinum)
  countlist <- countlist %>% group_by(mer) %>%
    summarise(readnum=sum(read), peaknum=length(read),
              peakshare=sum(peakshare), peakshareperc=sum(peakshare)/length(read)) %>% ungroup()
  return(countlist)
}

##a-site####
disome.asite <- disome.peak.kmer[, c('gene', 'asite', 'r1', 'r2', 'read', 'sc', 'peakshare')]
disome.asite$aaseq <- peak.seq(disome.asite, 'exit', 4, asite=T)
disome.asite$codonseq <- peak.seq(disome.asite, 'exit', 4, aa=F, asite=T)

mrna.asite <- mrna.peak.kmer[, c('gene', 'asite', 'r1', 'r2', 'read', 'sc')]
mrna.asite$aaseq <- peak.seq(mrna.asite, 'exit', 4, asite=T)
mrna.asite$codonseq <- peak.seq(mrna.asite, 'exit', 4, aa=F, asite=T)

mono.asite <- mono.peak.kmer[, c('gene', 'asite', 'r0', 'r1', 'r2', 'read', 'sc')]
mono.asite$aaseq <- peak.seq(mono.asite, 'exit', 4, asite=T)
mono.asite$codonseq <- peak.seq(mono.asite, 'exit', 4, aa=F, asite=T)

count_site <- function(indata, incol, mer, aa=T, count_read=T){
  # mer for codon or amino acid, not nucleotide
  # indata <- disome.asite
  # countlist <- kmers(indata, 'codonseq', 1, F, T)
  countlist <- kmers(indata, incol, mer, aa=aa, count_read=count_read)
  countlist <- countlist %>% group_by(gene, mer, posi) %>%
    summarise(num=sum(read)) %>% ungroup()
  countlist <- countlist %>% group_by(gene, posi) %>%
    mutate(notnum=sum(num)-num) %>% ungroup()
  return(countlist)
}

# rm stop codon
{
  start_time <- Sys.time()
  diasite.codon.read.rmsc <- count_site(disome.asite[disome.asite$sc == 0, ], 'codonseq', 1, aa=F, count_read=T)
  diasite.aa.read.rmsc <- count_site(disome.asite[disome.asite$sc == 0, ], 'aaseq', 1, aa=T, count_read=T)
  monoasite.codon.read.rmsc <- count_site(mono.asite[mono.asite$sc == 0, ], 'codonseq', 1, aa=F, count_read=T)
  monoasite.aa.read.rmsc <- count_site(mono.asite[mono.asite$sc == 0, ], 'aaseq', 1, aa=T, count_read=T)
  mrnaasite.codon.read.rmsc <- count_site(mrna.asite[mrna.asite$sc == 0, ], 'codonseq', 1, aa=F, count_read=T)
  mrnaasite.aa.read.rmsc <- count_site(mrna.asite[mrna.asite$sc == 0, ], 'aaseq', 1, aa=T, count_read=T)
  print(Sys.time()- start_time)
}


enrich.posi <- function(posi, indata, bgdata, overlapgene, codon){
  indata <- indata[indata$posi == posi, -3]
  bgdata <- bgdata[bgdata$posi == posi, -3]
  outdata <- enrich(indata, bgdata, overlapgene)
  
  if (codon){
    outdata <- merge(outdata, codon.tb[, c('Codon', 'Letter', 'optimality')],
                     by.x='mer', by.y='Codon', all.x=T)
  }
  return(outdata)
}

### disome
bgmrna.asite.codon.read <- enrich.posi('p4', diasite.codon.read, mrnaasite.codon.read, overlapgene=F, codon=T)
bgmrna.psite.aa.read.rmsc <- enrich.posi('p3', diasite.aa.read.rmsc, mrnaasite.aa.read.rmsc, overlapgene=F, codon=F)

### monosome
monobgmrna.asite.codon.read <- enrich.posi('p4', monoasite.codon.read, mrnaasite.codon.read, overlapgene=F, codon=T)
monobgmrna.psite.aa.read.rmsc <- enrich.posi('p3', monoasite.aa.read.rmsc, mrnaasite.aa.read.rmsc, overlapgene=F, codon=F)


##plots####

plot.mer.enrich <- function(indata, ingrep='H{3}|K{3}|R{3}',
                              colors=c('#6633CC', '#FF3300', '#FF33CC'),
                              # colors=c('#E78C62', '#6D83BD', '#6D83BD'),
                              # twotype=F,
                              mertext=T, cutoff=F, rmlegd=T){
  indata <- indata[!grepl(' ', indata$mer),]
  if (cutoff){
    indata$tag <- ifelse((indata$pvalue < cutoff) | (grepl(ingrep[1], indata$mer)),
    indata$mer, 'Other')
  }else{
    indata$tag <- ifelse(grepl(ingrep[1], indata$mer), 'type1', 'Other')
  }
  if (length(ingrep)>1){
    indata$tag <- ifelse(grepl(ingrep[2], indata$mer),
                         'type2', indata$tag)
  }
  print(table(indata$tag))
  oplot <- ggplot(indata[indata$tag == 'Other',], aes(x=log2(oddsratio), y=-log10(pvalue)))+
    geom_point(shape=1, alpha=0.8, color='#CCCCCC', size=1.5)+theme_classic()+
    geom_vline(xintercept=0, linetype='dashed', color='#3399FF')+
    theme(legend.position=c(.25,.85),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          axis.title=element_text(size=9),
          axis.text=element_text(size=8),
          axis.line=element_line(size=0.4),
          axis.ticks=element_line(size=0.4))+
    scale_x_continuous(limits=c(log2(0.10730160), log2(5.687468)),
                       breaks=seq(-3, 2, 1))+
    scale_y_continuous(limits=c(0, 30),
                       breaks=seq(0, 30, 10))+
    xlab(expression(log[2]*' (Odds Ratio)'))+
    ylab(expression(-log[10]~italic(P)))+
    geom_point(data=indata[indata$tag != 'Other',], alpha=0.8, size=1.5,
               aes(x=log2(oddsratio), y=-log10(pvalue), color=tag))
  # print(length(table(indata$tag))-1)
  if (length(table(indata$tag))-1<= length(colors)){
    oplot <- oplot + scale_colour_manual(values=colors,
                                         name='3-mers')
  }
  if (mertext){
    oplot <- oplot+
      geom_text(data=indata[indata$tag != 'Other',],
                aes(x=log2(oddsratio)-log2(1.2), y=-log10(pvalue), label=mer), size=2)
      # geom_text(data=indata[indata$tag != 'Other' | -log10(indata$pvalue)>10,],
      #           aes(x=log2(oddsratio)-log2(1.2), y=-log10(pvalue), label=mer), size=2)+
  }
  if (rmlegd){
    oplot <- oplot+theme(legend.position='none')
  }else{
    oplot <- oplot+
      theme(legend.position='top')+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
  }
  plot(oplot)
}

{
  pdf('P3_3mer_enrichment.pdf', width=4, height=3, useDingbats=F)
  plot.mer.enrich(bgmrna.3mer.aa.rmsc, c('(KKK)|(RRR)|(HHH)|(RKR)|(RKK)|(KRR)', '(QQQ)|(GGG)'),
                    c('#EA5514', '#6C84BA'), mertext=T)
  # plot.asite.enrich(bgmrna.3mer.aa.rmsc, '(KKK)|(RRR)|(HHH)|(RKR)|(RKK)|(KRR)|(QQQ)|(GGG)|(RGG)', c('#EA5514', '#6C84BA'),
  #                   twotype=T, mertext=T)
  plot.mer.enrich(monobgmrna.3mer.aa.rmsc[monobgmrna.3mer.aa.rmsc$method != 'fisher',],
                    c('(KKK)|(RRR)|(HHH)|(RKR)|(RKK)|(KRR)', '(QQQ)|(GGG)'),
                    c('#EA5514', '#6C84BA'), mertext=T)
  dev.off()
}

######

plot.asite.codon.enrich <- function(indata, site, codon=T, merfill='Letter', vline=c(-1, 1), flip=T){
  if (codon) {
    # indata$Letter <- factor(indata$Letter)
    # indata$optimality <- factor(indata$optimality)
    indata$mer <- factor(indata$mer, levels=indata[order(indata$oddsratio, decreasing=T),]$mer)
  }
  if (merfill=='stopcodon'){
    indata[, merfill] <- ifelse(indata$mer %in% c('TAA', 'TAG', 'TGA'), indata$mer, 'other')
  }
  # print(indata$mer)
  outplot <- ggplot(indata,
                    # aes(y=log2(oddsratio), x=mer))+
                    aes(y=oddsratio, x=mer))+
    geom_bar(stat='identity', aes_string(fill=merfill))+
    theme_classic()+
    scale_x_discrete(limits=levels(indata$mer))+
    scale_y_continuous(breaks=seq(0, 50, 10),
                       limits=c(0, 57))+
    # ylim(0, 57)+
    # geom_hline(yintercept=vline, linetype='dashed', color='#99ccff')+
    theme(axis.title=element_text(size=8, color='black'),
          axis.text.x=element_text(size=8, color='black', angle=-90),
          axis.text.y=element_text(size=8, color='black'),
          axis.line=element_line(size=0.4),
          axis.ticks=element_line(size=0.4),
          legend.position='top',#c(.8, .9),
          legend.background=element_blank(),
          legend.text=element_text(size=8),
          legend.title=element_text(size=8),
          legend.key.size=unit(.6, 'cm'),
          plot.caption=element_text(size=8))+
    labs(caption=sprintf('(%s-site)', site))+
    # ylab(expression(log[2]*'(Odds Ratio)'))
    ylab('A-site pausing score')
  if (merfill == 'stopcodon'){
    outplot <- outplot+
      scale_fill_manual(values=c("#8B9CC4", "#64BA9F", '#F08961', '#898989'),
                        name='stopcodon',
                        labels=c('A', 'B', 'C', 'D'),
                        guide='none')
      
  }
  if (merfill == 'optimality') {
    outplot <- outplot +
      scale_fill_manual(values=c("#FF9900", "#3399FF"),
                        name='Optimality',
                        labels=c('Unpreferred', 'Preferred'),
                        guide='none')
  }
  if (codon) {
    outplot <- outplot+
      # geom_text(aes(y=max(log2(oddsratio))*1.1, x=mer, label=optimality, color=optimality), size=3)+
      xlab('Codons')
  } else{
    outplot <- outplot + xlab('Animo Acids')
  }
  if (flip) outplot <- outplot+coord_flip()
  plot(outplot)
}

{
  pdf('P3_asite_codon_enrichment_opt_sc.pdf', width=2.5, height=6)
  plot.asite.codon.enrich(bgmrna.asite.codon.read[bgmrna.asite.codon.read$mer != '   ',], site='di A', merfill='stopcodon')
  plot.asite.codon.enrich(monobgmrna.asite.codon.read[monobgmrna.asite.codon.read$mer != '   ',], site='mono A', merfill='stopcodon')
  dev.off()
}

plot.psite.aa.vol <- function(indata, tlog=T){
  require(RColorBrewer)
  indata <- indata[!(indata$mer %in% c(' ', '*')),]
  indata$tag <- ifelse(indata$mer %in% c('P', 'G', 'N', 'K', 'C', 'R'), indata$mer, 'other')

  if (tlog) {
    indata$oddsratio <- log2(indata$oddsratio)
  }

  outplot <- ggplot(indata[indata$tag == 'other',], aes(x=oddsratio, y=-log10(pvalue), color=tag))+
    geom_point(color='#595757')+theme_classic()+
    # scale_color_manual(values=brewer.pal(8, 'Dark2'))+
    geom_point(data=indata[indata$tag != 'other',])+
    scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#E6AB02", "#A6761D"))+
    theme(axis.text=element_text(size=8, color='black'),
          axis.title=element_text(size=8, color='black'),
          legend.position='none')+
    ylab('-Log10P')+
    xlab('P-site pausing score')
    
    if (tlog){
     outplot <- outplot+
       geom_text(data=indata[indata$mer %in% c('P', 'G', 'N', 'K', 'C', 'R'), ],
                 aes(label=mer, x=oddsratio, y=-log10(pvalue), hjust=-1))
    } else {
      outplot <- outplot+
        geom_text(data=indata[indata$mer %in% c('P', 'G', 'N', 'K', 'C', 'R'), ],
                  aes(label=mer, x=oddsratio, y=-log10(pvalue), hjust=-1))+
        scale_x_continuous(limits=c(0.5, 2.59632603),
                           breaks=seq(0.5, 2.5, 0.5))+
        scale_y_continuous(limits=c(0, 186.0619),
                           breaks=seq(0, 160, 40))
    }
    
    # scale_color_brewer(palette='Set2', direction=-1)+
  return(outplot)
}

{
  pdf('P3_psite_aa_enrichment_vol_rmsc.pdf', width=2, height=1.8, useDingbats=F)
  plot(plot.psite.aa.vol(bgmrna.psite.aa.read.rmsc))
  plot(plot.psite.aa.vol(bgmrna.psite.aa.read.rmsc, tlog=F))
  plot(plot.psite.aa.vol(monobgmrna.psite.aa.read.rmsc, tlog=F))
  dev.off()
}

## helix ####
aa.ddG <- read.table('./helix_propensity.txt', header=T, sep='\t')
helix.ddG <- bgmrna.psite.aa.read.rmsc[bgmrna.psite.aa.read.rmsc$mer != ' ',]
helix.ddG <- merge(bgmrna.psite.aa.read.rmsc, aa.ddG, by.x='mer', by.y='aa')


# mono.helix.ddG <- monobgmrna.psite.aa.read.rmsc[monobgmrna.psite.aa.read.rmsc$mer != ' ',]
mono.helix.ddG <- merge(monobgmrna.psite.aa.read.rmsc, aa.ddG, by.x='mer', by.y='aa')
# ggplot(helix.ddG, aes(x=log2(oddsratio), y=ddG))+
#   geom_point()+theme_classic()+
#   geom_text(data=helix.ddG[helix.ddG$mer%in% c('R', 'K'),],
#             aes(x=log2(oddsratio), y=ddG, hjust=2, label=mer))+
#   annotate('text', y=Inf, x=-Inf, vjust=1, hjust=-0.3,
#            label=sprintf('cor=%.2f, P=%.2e', 0.5551833, 0.01105))

cor.test(helix.ddG$oddsratio, helix.ddG$ddG)
cor.test(log2(helix.ddG$oddsratio), helix.ddG$ddG)
cor.test(helix.ddG$oddsratio, helix.ddG$ddG, method='s')
cor.test(mono.helix.ddG$oddsratio, mono.helix.ddG$ddG)

{
  pdf('P3_ddG_oddsratio.pdf', width=2.8, height=2.8, useDingbats=F)
  plot(
    ggplot(helix.ddG, aes(x=oddsratio, y=ddG))+
      geom_point()+theme_classic()+
      geom_text(data=helix.ddG[helix.ddG$mer%in% c('R', 'P', 'G', 'N', 'K', 'C'),],
                aes(x=oddsratio+0.1, y=ddG, label=mer))+
      annotate('text', y=Inf, x=-Inf,vjust=1, hjust=-0.3,
               label=sprintf('cor=%.2f, P=%.2e', 0.6209966, 0.003477))+
      xlab('P-site pausing score')+
      theme(axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))
  )
  plot(
    ggplot(mono.helix.ddG, aes(x=oddsratio, y=ddG))+
      geom_point()+theme_classic()+
      geom_text(data=mono.helix.ddG[mono.helix.ddG$mer%in% c('R', 'P', 'G', 'N', 'K', 'C'),],
                aes(x=oddsratio+0.1, y=ddG, label=mer))+
      annotate('text', y=Inf, x=-Inf,vjust=1, hjust=-0.3,
               label=sprintf('cor=%.2f, P=%.2f', 0.195113, 0.4097))+
      xlab('P-site pausing score')+
      theme(axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))+
      scale_x_continuous(limits=c(0.5, 2.5963260),
                         breaks=seq(0.5, 2.5, 0.5))
  )
  plot(
    ggplot(helix.ddG, aes(x=log2(oddsratio), y=ddG))+
      geom_point()+theme_classic()+
      geom_text(data=helix.ddG[helix.ddG$mer%in% c('R', 'P', 'G', 'N', 'K', 'C'),],
                aes(x=log2(oddsratio), y=ddG, hjust=2, label=mer))+
      annotate('text', y=Inf, x=-Inf, vjust=1, hjust=-0.3,
               label=sprintf('cor=%.2f, P=%.2e', 0.5551833, 0.01105))+
      xlab('P-site pausing score')+
      theme(axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))
  )
  dev.off()
}

## the end####




































