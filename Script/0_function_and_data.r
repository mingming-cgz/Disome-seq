library(tidyverse)

## readtable ####

read.data <- function(indir, ptn){
  incol <- c(rep('NULL', 5), rep('character', 1),
             rep('NULL', 7), rep('numeric', 3),
             rep('NULL', 2), rep('numeric', 5),
             'character', 'numeric')
  dirpath <- file.path('~/projects/disome/manuscript', indir, 'samout')
  #   print(dirpath)
  outdata <- data.frame()
  inf <- dir(dirpath)[!grepl('_puta', dir(dirpath)) & grepl(ptn, dir(dirpath))]
  # print(inf)
  for (f in inf){
    tag <- strsplit(f, '[_\\.]')[[1]]
    r <- gsub('\\D+-', '', tag[1])
    # print(tag)
    tag <- paste(indir, tag[2], sep='-')
    f <- file.path(dirpath, f)
    print(f)
    indata <- read.table(f, header=T, sep='\t', colClasses=incol)
    indata$rep <- as.numeric(substr(r, nchar(r), nchar(r)))-1
    indata$tag <- tag
    # print(head(indata))
    outdata <- rbind(outdata, indata)
  }
  return(outdata)
}

read.dir <- function(indir, outtag, ptn='\\.out'){
  outframe <- data.frame()
  for (i in indir){
    tmp <- read.data(i, ptn)
    outframe <- rbind(outframe, tmp)
  }
  outframe <- outframe[outframe$region != 'intron',]
  outtag <- paste0(outtag, '_ori.txt')
  write.table(outframe, outtag, sep='\t', row.names=F, quote=F)
}

read.dir(c('201710', '20171107', '20171116'), 'disome')
read.dir(c('monosome'), 'monosome')
read.dir(c('mRNA'), 'mRNA')
read.dir(c('di_3at'), 'di_3at')
read.dir(c('mono_3at'), 'mono_3at')
read.dir(c('mRNA_3at'), 'mRNA_3at')

###### 
correct.pos <- function(indata){
  outpos <- (1-indata$frame5)%%3-1 + indata$dis5
  return(outpos)
}

get.asite <- function(indata){
  asite <- (indata$seqlen %/% 3)*3 + indata$cor.pos - 12
  return(asite)
}

peak.sum <- function(indata, groupcol='asite'){
  if (sum(colnames(indata) == 'asite')==0){
    indata$cor.pos <- correct.pos(indata)
    indata$asite <- get.asite(indata) 
  }
  indata <- indata %>% group_by_at(vars(c('rep', 'gene', groupcol))) %>%
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

######
## correlation between replications ####

data.read <- function(inori, picklen=c(-Inf, Inf), inrep=NULL){
  inori <- inori[inori$seqlen >= picklen[1] & inori$seqlen <= picklen[2], ]
  inori <- inori %>% group_by(gene, rep) %>%
    summarise(reads=sum(seqsnum)) %>% ungroup()
  if (is.null(inrep)) {
    inori$rep <- paste0('r', inori$rep) 
  }
  inrep <- unique(inori$rep)
  inori <- inori %>% spread(rep, reads)
  inori[is.na(inori)] <- 0
  inori$read <- rowMeans(inori[, inrep])
  inori$rpm <- rpm.mean(inori, inrep)
  return(inori)
}

disome.ori <- read.table('disome_ori.txt', header=T, sep='\t', stringsAsFactors=F)
disome.ori <- disome.ori[disome.ori$rep != 0, ]
disome.read <- data.read(disome.ori, c(55, 65))
cor.test(disome.read$r1, disome.read$r2)
rm(disome.ori)

mono.ori <- read.table('monosome_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mono.read <- data.read(mono.ori, c(20,30))
mono.rpkm <- mono.read[, c('gene', 'rpm')]
mono.rpkm <- merge(mono.rpkm, cds.len, by='gene', all.x=T)
mono.rpkm$rpkm <- mono.rpkm$rpm / (mono.rpkm$len * 1e-3)
mono.rpkm$len <- NULL

mrna.ori <- read.table('mRNA_ori.txt', header=T, sep='\t', stringsAsFactors=F)
mrna.read <- data.read(mrna.ori)
rm(mrna.ori)


library(scales)
multiplot <- function(..., plotlist=NULL, cols=1, rows=NULL, layout=NULL, byrow=T, naposi=NULL){
  require(grid)
  plots <- c(list(...), plotlist)
  plotnum <- length(plots)
  # print(plots[1])
  if(is.null(layout)){
    if (is.null(rows)) rows <- ceiling(plotnum/cols)
    layout <- seq(1, cols*rows)
    if (!is.null(naposi)) {
      layout[naposi] <- NA
    }
    layout <- matrix(layout, ncol=cols, byrow=byrow)
    print(layout)
    layoutnum <- sort(as.numeric(layout[!is.na(layout)]), decreasing=F)
    print(layoutnum)
  }
  
  if (length(layoutnum) != plotnum) stop('Grid number is not equal to the plot number')
  
  if (plotnum==1){
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:plotnum){
      matchidx <- which(layout == layoutnum[i], arr.ind = TRUE)
      plot(plots[[i]], vp = viewport(layout.pos.row=matchidx[1],
                                      layout.pos.col=matchidx[2]))
    }
  }
}

plot.read.cor <- function(indata, row1, row2, datatag, takelog=T){
  # print(str(indata[, row1]))
  indata <- as.data.frame(indata)
  incor <- cor.test(indata[, row1], indata[, row2])
  if (takelog) {
    indata[, row1] <- log2(indata[, row1])
    indata[, row2] <- log2(indata[, row2])
  }
  datatag <- sprintf('# of %s Reads in Replication %s,',
                     datatag, gsub('r', '', c(row1, row2)))
  outplot <- ggplot(indata, aes_string(row1, row2))+
    geom_point()+theme_classic()+
    xlab(bquote(.(datatag[1])~log[2]))+
    ylab(bquote(.(datatag[2])~log[2]))+
    theme(axis.title=element_text(size=8, color='black'),
          axis.text=element_text(size=8, color='black'),
          axis.line=element_line(size=.6),
          axis.ticks=element_line(size=.6))+
    annotate('text', x=-Inf, y=Inf, vjust=1, hjust=-.05,
             label=sprintf('cor=%.4f    P=%.4e', incor$estimate, incor$p.value),
             size=5)
    
  return(outplot)
}

##the end####
























