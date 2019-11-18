#! /usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2018-08-27 15:31:06
# @Author  : Yanming Chen 
# @Version : 1.0.0

start_time <- Sys.time()
library(optparse)
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(IRanges))

optslist <- list(
    make_option(c('-f', '--file'), type='character', default=NULL,
                help='your input bam file'),
    make_option(c('-o', '--odir'), type='character', default=NULL,
                help='your output directory'),
    make_option(c('-n', '--nobarcode'), action='store_false', default=NULL,
                help='no barcode tag'),
    make_option(c('-s', '--se'), action='store_false', default=NULL,
                help='single end tag')
)

opt_parser <- OptionParser(option_list=optslist)
opt <- parse_args(opt_parser)

if(is.null(opt$file)){
    print_help(opt_parser)
    stop('Please input your .bam of the ', call.=FALSE)
}

if (is.null(opt$se)) {se <- TRUE} else {se <- opt$se}
if (is.null(opt$nobarcode)) {nobc <- TRUE} else {nobc <- opt$nobarcode}


readBam2TSSRangedDataFilter <- function(bamfile, barcode=T, filtunique=T, paired=T){
    
    
    bamfile <- file.path(bamfile, 'accepted_hits.bam')
    # filter with tag for the unique mapping or not
    if(filtunique){
        gbam <- readGAlignments(bamfile, use.names=T, param=ScanBamParam(what=c('mapq', 'flag'), tag='NH'))
    }else{
        gbam <- readGAlignments(bamfile, use.names=T, param=ScanBamParam(what=c('mapq', 'flag')))
    }


    if(filtunique){
        fqId <- mcols(gbam)$mapq < 30 | mcols(gbam)$NH != 1
    }
    else{
        fqId <- mcols(gbam)$mapq < 30
    }
    
    keepID = !fqId

    if (paired) {
       pairId <-  mcols(gbam)$flag == 147 | mcols(gbam)$flag == 163
       keepID = keepID & !pairId
    }
    
    # keepID = !fqId & !pairId

    print(paste0(sum(!keepID), '(', sum(!keepID)/length(gbam)*100, '%)',
                 'are filtered out in the quality control'))
    print(paste0(sum(keepID), '(', sum(keepID)/length(gbam)*100, '%)',
                 'are kept in the quality control'))
    
    gfbam <- gbam[keepID]
    pos <- ifelse(as.character(gfbam@strand) == '+', start(gfbam), end(gfbam))
    flen <- width(gfbam)
    uniIds <- paste(as.character(gfbam@seqnames),
                    as.character(gfbam@strand),
                    as.character(flen),
                    pos, sep='_')
    if (barcode){
        barSeqs <- unlist(lapply(strsplit(names(gfbam), ','), function(x) x[2]))
        rbls <- split(barSeqs, uniIds)
        unum <- unlist(lapply(rbls, function(x) length(unique(x))))
        snum <- unlist(lapply(rbls, length))
    }else {
        rbls <- split(1:length(gfbam), uniIds)
        unum <- unlist(lapply(rbls, length))
        snum <- unum
    }

    tssl <- strsplit(names(rbls), '_')
    tchr <- unlist(lapply(tssl, function(x) x[1]))
    tstrand <- unlist(lapply(tssl, function(x) x[2]))
    tflen <- unlist(lapply(tssl, function(x) x[3])) %>% as.numeric()
    tpos <- unlist(lapply(tssl, function(x) x[4])) %>% as.numeric()

    outdata <- RangedData(IRanges(start=tpos, width=1), strand=tstrand,
                          space=paste(tchr, tstrand, sep='_'),
                          chr=tchr, flen=tflen, unum=unum, snum=snum)
    # print(head(outdata))
    return(outdata)
}

ovORFsREADs <- function(extCDS, inreads, outfile){

    # extCDS -- query; inreads -- subject
    fo <- findOverlaps(extCDS, inreads)
    mfo <- as.matrix(fo)

    # sls <- split(mfo[, 2], mfo[, 1])
    genoId <- mfo[, 1]
    seqId <- mfo[, 2]
    # mfo <- cbind(mfo, 
    #              duplicate=ifelse(duplicated(seqId,fromLast=T)|
    #                               duplicated(seqId,fromLast=F), 1, 0))

    minfo <- data.frame(efir=extCDS$first[genoId], elast=extCDS$last[genoId],
                        gId=genoId, rId=seqId,
                        strand=extCDS$strand[genoId], gene=extCDS$geneId[genoId],
                        rank=extCDS$rank[genoId], rdis5=extCDS$rdis5[genoId],
                        rdis3=extCDS$rdis3[genoId], isFir=extCDS$isFir[genoId],
                        isLast=extCDS$isLast[genoId], seqchr=inreads$chr[seqId],
                        seqpos=start(inreads)[seqId], seqlen=inreads$flen[seqId],
                        sequnum=inreads$unum[seqId], seqsnum=inreads$snum[seqId])
    
    minfo$idis5 <- with(minfo,ifelse(strand=="+",seqpos-efir,elast-seqpos))
	minfo$idis3 <- with(minfo,ifelse(strand=="+",seqpos-elast, efir-seqpos))

    minfo$dis5 <- with(minfo,idis5+rdis5)
	minfo$dis3 <- with(minfo,idis3-rdis3)

    # minfo$inORF <- ifelse(with(minfo, dis5>=0 & dis3 <=0), 1, 0)
    minfo$frame5 <- with(minfo,dis5 %% 3)
    minfo$frame3 <- with(minfo,(dis3+2) %% 3)
    minfo$putative <- ifelse(minfo$gene %in% putative, 1, 0)

    minfo$region <- 0
    minfo$region[(minfo$idis5>=0&minfo$idis3<=0)] <- 'orf'
    minfo$region[(minfo$dis5<0&minfo$isFir==1)] <- '5putr'
    minfo$region[(minfo$dis3>0&minfo$isLast==1)] <- '3putr'
    # remove wrong alignment because of the peudoovelapping in annotation (+100)
    minfo <- minfo[!(minfo$region==0 & minfo$rId %in% minfo[minfo$region!=0,'rId']),]
    minfo[minfo$region==0,]$region <- 'intron' # ignore intron overlapping region    

    minfo$dupl <- ifelse(duplicated(minfo$rId,fromLast=T)| duplicated(minfo$rId,fromLast=F), 1, 0)
    
    write.table(minfo, outfile[1], sep='\t', quote=F, row.names=F)

    # if reads duplicatedly mapped to ORF and putative gene, remove putative gene
    minfo <- minfo[!(minfo$dupl==1 & minfo$putative == 1), ]
    write.table(minfo, outfile[2], sep='\t', quote=F, row.names=F)
    print(sum(minfo$seqsnum))
}

print(paste0('paired end: ', se, '; barcode: ', nobc))
home <- path.expand('~')
dubious <- read.table(paste0(home,'/projects/reference/annotation/dubious_R64-1-1.txt'),
                      header=F, sep='\t')$V1 %>% as.character
# only "Y" and "Q" genes in the annotation
anno <- readRDS(paste0(home,'/projects/reference/annotation/Saccharomyces_cerevisiae.R64-1-1.89.rds'))
# remove dubious gene
anno <- anno[!anno$geneId %in% dubious,]
putative <- read.table(paste0(home,'/projects/reference/annotation/putative_dubious_R64-1-1.txt'),
                      header=F, sep='\t')$V1 %>% as.character

outf <- unlist(strsplit(opt$file, '/'))
# print(outf)
outf <- gsub('.bam', '', outf[length(outf)])
outf <- c(file.path(opt$odir, paste0(outf, '_puta.out')),
          file.path(opt$odir, paste0(outf, '.out')))
# print(outf)

reads <- readBam2TSSRangedDataFilter(opt$file, barcode=nobc, paired=se)
ovORFsREADs(anno, reads, outf)
print(Sys.time()- start_time)