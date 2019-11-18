library(tidyverse)
suppressMessages(library(GenomicAlignments))
suppressMessages(library(IRanges))

prot.domain <- read.table('./protein_domain/S288C_secondary_structure.txt',
                          header=T, sep='\t', stringsAsFactors=F)
# start and end position in amino acid
prot.domain$start <- as.numeric(prot.domain$start)
prot.domain$end <- as.numeric(prot.domain$end)
prot.domain <- na.omit(prot.domain)

prot.mid <- prot.domain[prot.domain$domain %in%
                          c('helix', 'strand', 'turn'), ]
prot.mid$uniprot_accession <- NULL
prot.mid <- prot.mid[order(prot.mid$gene, prot.mid$start),]
rownames(prot.mid) <- NULL
prot.mid$aalen <- prot.mid$end - prot.mid$start + 1

define.helix.gap <- function(indomain, helixlencutoff){
  # indomain <- prot.mid
  # helixlencutoff <- 4
  indomain <- indomain[indomain$domain == 'helix', ]
  colsindomain <- colnames(indomain)
  defind.gap <- function(indata){
    indata$start.follow <- indata$end+1
    indata <- cbind(indata,
                      rbind(indata[-1, c('gene', 'domain', 'start')],
                            c(NaN, NaN, NaN)))
    colnames(indata)[7:9] <- c('gene.follow', 'domain.follow', 'start.next.domain')
    indata$end.follow <- indata$start.next.domain - 1
    indata$gap <- ifelse(indata$start.follow != indata$start.next.domain, 1, 0)
    return(indata)
  }
  indomain <- defind.gap(indomain)
  indomain$com <- ifelse(indomain$gap == 1, 1, 0)
  continue.helix.head <- which(indomain$gap == 0)[!(which(indomain$gap == 0) %in% (which(indomain$gap == 0)+1))]
  indomain$com[continue.helix.head] <- 1
  continue.helix.ht <- indomain[!(indomain$gap == 0 & indomain$com == 0), ]
  continue.helix.ht <- continue.helix.ht[order(continue.helix.ht$gene, continue.helix.ht$start, continue.helix.ht$end),]
  rownames(continue.helix.ht) <- NULL
  com.h <- which(continue.helix.ht$gap == 0 & continue.helix.ht$com == 1)
  
  continue.helix.remain <- continue.helix.ht[-c(com.ht, (com.ht+1)),]
  continue.helix.ht <- cbind(continue.helix.ht[com.h, 1:4],
                             continue.helix.ht[com.h+1, 5:12])
  indomain <- rbind(continue.helix.ht, continue.helix.remain)
  rm(continue.helix.head, continue.helix.ht, com.h, continue.helix.remain)
  
  indomain <- indomain[order(indomain$gene, indomain$start),]
  indomain$aalen <- indomain$end - indomain$start + 1
  indomain <- indomain[indomain$aalen>helixlencutoff,]
  
  indomain <- indomain[, colsindomain]
  indomain <- defind.gap(indomain)
  
  out.helix <- merge(indomain, cds.pep[, c('gene', 'len')],
                      by='gene', all.x=T)
  out.helix$len <- out.helix$len / 3 - 1
  colnames(out.helix)[grep('len', colnames(out.helix))] <- 'pep.len'
  out.helix$end.follow <- ifelse(out.helix$gene != out.helix$gene.follow,
                                 out.helix$pep.len, out.helix$end.follow)
  out.helix$end.follow <- ifelse(out.helix$end == out.helix$pep.len,
                                 out.helix$start.follow, out.helix$end.follow)
  out.helix <- out.helix[out.helix$end <= out.helix$pep.len,]
  out.helix <- out.helix[, c('gene', 'domain', 'start', 'end', 'start.follow', 'end.follow')]
  out.helix$gaplen <- out.helix$end.follow - out.helix$start.follow + 1
  out.helix$end.follow <- ifelse(out.helix$gaplen > 30, out.helix$start.follow+30, out.helix$end.follow)
  out.helix$gaplen <- NULL
  
  out.helix$domain.g <- 'gap'
  tmp <- out.helix[, c('gene', 'domain.g', 'start.follow', 'end.follow')]
  colnames(tmp) <- c('gene', 'domain', 'start', 'end')
  out.helix <- rbind(out.helix[, c('gene', 'domain', 'start', 'end')], tmp)
  rm(tmp)
  
  out.helix <- na.omit(out.helix)
  out.helix$start_nt <- (out.helix$start-1)*3
  out.helix$end_nt <- out.helix$end*3-1
  
  out.helix <- out.helix[order(out.helix$gene, out.helix$start_nt, out.helix$end_nt),]
  rownames(out.helix) <- NULL
  return(out.helix)
}

helix.range <- function(indata){
  outdata <- with(indata, RangedData(IRanges(start=start_nt, end=end_nt),
                                      space=gene, domain=domain, gene=gene))#,
                                      # no=no, norev=norev))
  return(outdata)
}

read.in.helix <- function(inread, inhelix, aaposical=0){
  aaposical <- aaposical*3
  inread$asite <- inread$asite - aaposical
  # print(head(inread))
  datarange <- with(inread,
                    RangedData(IRanges(start=asite, width=1),
                               space=gene, rpm=rpm, gene=gene))
  # print(head(datarange))
  fo <- findOverlaps(datarange, inhelix, type='within')
  mfo <- as.matrix(fo)
  # print(nrow(mfo))
  
  helixID <- mfo[, 2]
  dataID <- mfo[, 1]
  
  minfo <- data.frame(helixgene=inhelix$gene[helixID], helixstart_nt=start(inhelix)[helixID],
                      helixend_nt=end(inhelix)[helixID],
                      # helixno=inhelix$no[helixID],helixnorev=inhelix$norev[helixID],
                      helixdomain=inhelix$domain[helixID],
                      #seqgene=datarange$gene[dataID],
                      seqposi=start(datarange)[dataID], seqrpm=datarange$rpm[dataID])
  minfo$tohelixend <- ifelse(minfo$helixdomain == 'gap',
                             minfo$seqposi - minfo$helixstart_nt,
                             minfo$seqposi-(minfo$helixend_nt+1))
  minfo$helixlen <- minfo$helixend_nt - minfo$helixstart_nt + 1
  return(minfo)
}


helix.gt4 <- define.helix.gap(prot.mid, 4)
helix.gt4.range <- helix.range(helix.gt4)

helix.disome.asite <- read.in.helix(disome.peak[disome.peak$sc != 1,], helix.gt4.range)
ggplot(helix.disome.asite, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+#xlim(-50, 50)+
  labs(caption='(disome)')

helix.mono.asite <- read.in.helix(mono.peak[mono.peak$sc != 1,], helix.gt4.range)
ggplot(helix.mono.asite, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+#xlim(-50, 50)+
  labs(caption='(monosome)')

helix.mrna.asite <- read.in.helix(mrna.peak[mrna.peak$sc != 1,], helix.gt4.range)
ggplot(helix.mrna.asite, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+#xlim(-50, 50)+
  labs(caption='(mRNA)')

disome.peak.pro <- disome.asite[grepl('\\w{2}P\\w{1}', disome.asite$aaseq),]
disome.peak.pro <- merge(disome.peak.pro, disome.peak[, c('gene', 'asite', 'rpm')],
                         by=c('gene', 'asite'))

mono.peak.pro <- mono.asite[grepl('\\w{2}P\\w{1}', mono.asite$aaseq),]
mono.peak.pro <- merge(mono.peak.pro, mono.peak[, c('gene', 'asite', 'rpm')],
                         by=c('gene', 'asite'))


helix.disome.pro.asite <- read.in.helix(disome.peak.pro[disome.peak.pro$sc != 1,], helix.gt4.range)
ggplot(helix.disome.pro.asite, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+#xlim(-50, 50)+
  labs(caption='(disome, pro)')

helix.mono.pro.asite <- read.in.helix(mono.peak.pro[mono.peak.pro$sc != 1,], helix.gt4.range)
ggplot(helix.mono.pro.asite, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+#xlim(-50, 50)+
  labs(caption='(disome, pro)')

pro.helix.gene <- unique(helix.disome.pro.asite[helix.disome.pro.asite$helixdomain == 'helix',]$helixgene)


write.table(helix.disome.pro.asite, 'helix.disome.pro.asite.txt', sep='\t', quote=F, row.names=F)
write.table(helix.disome.pro.asite[!(helix.disome.pro.asite$helixgene %in% pro.helix.gene),],
            'helix.disome.pro.rmhelixgenes.txt', sep='\t', quote=F, row.names=F)

helix.cds.pro <- read.in.helix(cds.pro.foralign, helix.gt4.range)
ggplot(helix.cds.pro, aes(x=tohelixend, weight=seqrpm))+
  geom_freqpoly(binwidth=3)+theme_classic()+xlim(-50, 50)+
  labs(caption='(genome, pro)')

######
{
  pdf('helix.test.pdf', width=18, height=6)
  multiplot(
    ggplot(helix.disome.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+#xlim(-50, 50)+
      labs(caption='(disome)'),
    ggplot(helix.disome.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+xlim(-50, 50)+
      labs(caption='(disome)'),
    ggplot(helix.mono.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+#xlim(-50, 50)+
      labs(caption='(monosome)'),
    ggplot(helix.mono.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+xlim(-50, 50)+
      labs(caption='(monosome)'),
    ggplot(helix.mrna.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+#xlim(-50, 50)+
      labs(caption='(mRNA)'),
    ggplot(helix.mrna.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+xlim(-50, 50)+
      labs(caption='(mRNA)'),
    ggplot(helix.disome.pro.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+#xlim(-50, 50)+
      labs(caption='(disome, pro)'),
    ggplot(helix.disome.pro.asite, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+xlim(-50, 50)+
      labs(caption='(disome, pro)'),
    ggplot(helix.cds.pro, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+#xlim(-50, 50)+
      labs(caption='(genome, pro)'),
    ggplot(helix.cds.pro, aes(x=tohelixend, weight=seqrpm))+
      geom_freqpoly(binwidth=1)+theme_classic()+xlim(-50, 50)+
      labs(caption='(genome, pro)'),
    
    cols=5, byrow=F
  )
  dev.off()
}
######
## proline out####

{
  pdf('helix.disome.pro.asite.pdf', width=3.4, height=2)
  plot(
    ggplot(helix.disome.pro.asite, aes(x=tohelixend/3, weight=seqrpm))+
      geom_bar()+theme_classic()+#xlim(-17, 17)+
      labs(caption='(disome, pro)')+
      theme(axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))+
      scale_x_continuous(breaks=c(-16, -8, 0, 8, 16),
                         limits=c(-17, 17))+
      xlab('Position relative to\nthe last amino acid of an a-helix')+
      ylab('Disome footprint\nabundance (RPM)')
  )
  plot(
    ggplot(helix.mono.pro.asite, aes(x=tohelixend/3, weight=seqrpm))+
      geom_bar()+theme_classic()+#xlim(-17, 17)+
      labs(caption='(mono, pro)')+
      theme(axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))+
      scale_x_continuous(breaks=c(-16, -8, 0, 8, 16),
                         limits=c(-17, 17))+
      xlab('Position relative to\nthe last amino acid of an a-helix')+
      ylab('Monosome footprint\nabundance (RPM)')
  )
  dev.off()
}



######




######



























































## the end####