exFHsibrel_GRM <- function(grm,pedigree,genolinkped,typefam,exclfamsize=1,outname){
  grm <- read.table(paste(grm),stringsAsFactors=F)[,1:3]
  cat('\n... Genomic relationships imported ...\n')
  colnames(grm) <- c('IID1','IID2','G')
  pedig <- read.table(paste(pedigree),stringsAsFactors=F)[,1:3]
  pedig <- na.omit(pedig)
  cat('... Pedigree imported ...\n')
  colnames(pedig) <- c('IID','Sire','Dam')
  genoanim <- read.table(paste(genolinkped),stringsAsFactors=F)[,1:2]
  cat('... IDs linking pedigree to genomic relationships imported ...\n')
  colnames(genoanim) <- c('IID','seqID')
  genoanimped <- merge(pedig,genoanim,by='IID')
  
  if(typefam=='FS'){
    genoanimped$famID <- paste(genoanimped$Sire,genoanimped$Dam,sep='_')
    genoanimped$famIDseqID <- unclass(as.factor(genoanimped$famID))
    cat('\n... Full-sib families created using pedigree information (sire_dam) ...\n')
    usefam <- sort(table(genoanimped$famIDseqID))
    allfam <- length(usefam)
    usefam <- sort(as.numeric(as.vector(names(usefam[usefam>exclfamsize]))))
    cat('... ',length(usefam),'Full-sib families (sire_dam) obtained ...\n')
    if((allfam-length(usefam))<=exclfamsize){
      cat('... ',allfam-length(usefam),' family(s) had only ',exclfamsize,' offspring, thus discarded ...\n')} else {
        cat('... ',allfam-length(usefam),' families had only one offspring, thus discarded ...\n')
      }
    cat('\n... Extracting full-sib relationships per family ...\n')
    genoanimped <- genoanimped[genoanimped$famIDseqID %in% usefam,]
    iterchecks <- round(length(usefam)/5,digits=0)
    FSgrm <- data.frame(IID1=character(),IID2=character(),G=numeric(),seqfamID=character())
    for(f in 1:length(usefam)){
      fullsibsg <- genoanimped[genoanimped$famIDseqID %in% usefam[f],]
      fullsibsg <- expand.grid(as.vector(fullsibsg$seqID),as.vector(fullsibsg$seqID),stringsAsFactors=F)
      fullsibsg <- fullsibsg[which(fullsibsg$Var1!=fullsibsg$Var2),c(2,1)]
      fullsibsg <- fullsibsg[order(fullsibsg$Var2,fullsibsg$Var1),]
      fullsibsg <- fullsibsg[which(fullsibsg$Var2<fullsibsg$Var1),]
      for(fsize in 1:nrow(fullsibsg)){
        fullsibgrm <- grm[which((grm$IID1==fullsibsg$Var2[fsize] & grm$IID2==fullsibsg$Var1[fsize]) | 
                                  (grm$IID2==fullsibsg$Var2[fsize] & grm$IID1==fullsibsg$Var1[fsize])),]
        fullsibgrm$famIDseqID <- usefam[f]
        FSgrm <- rbind.data.frame(FSgrm,fullsibgrm,stringsAsFactors=F)
      }
      if(f %% iterchecks==0){cat('... fullsib family ...',f,' (out of ',length(usefam),') ...\n')}
      rm(fullsibgrm,fullsibsg)
    }
    colnames(FSgrm) <- c('IID1','IID2','G','seqfamID')
    
    plotname <-paste(outname,'.tiff',sep='')
    dpi=300;width.cm<-18;height.cm<-13;width.in<-width.cm/2.54;height.in<-height.cm/2.54
    tiff(file=plotname,width=width.in*dpi,height=height.in*dpi,pointsize=10,units="px",res=dpi)
    layout(matrix(1:2,1,2))
    par(mar=c(4,4,0.1,0.1))
    xlim=range(FSgrm$G)+c(-0.20,0.20)
    hist(FSgrm$G,breaks=100,col=sample(colors(),1,replace=F),xlim=xlim,
         xlab='Full-sib genomic relationships',main='')
    abline(v=mean(FSgrm$G),col='red',lwd=2)
    abline(v=0.5,col='green',lwd=1)
    legend('topleft',legend=c(paste('Average G (',round(mean(FSgrm$G),3),')',sep=''),
                              'Pedigre-based expectation (= 0.500)'),
           col=c('red','green'),lty=1,lwd=1,cex=0.6,bty='n')
    par(mar=c(4,3.8,0.1,0.5))
    plot(density(FSgrm$G),col=sample(colors(),1,replace=F),xlim=xlim,
         xlab='Full-sib genomic relationships',main='',type='l',lwd=1.5)
    points(density(FSgrm$G,adjust=2),col=sample(colors(),1,replace=F),xlim=xlim,
           xlab='Full-sib genomic relationships',main='',type='l',lty=2)
    dev.off()
    write.csv(FSgrm,paste(outname,'_FSrelationships.csv',sep=''),quote=F,row.names=F)
    return(FSgrm)
  } else if (typefam=='HS'){
    genoanimped$famID <- paste(genoanimped$Sire,sep='_')
    genoanimped$famIDFS <- paste(genoanimped$Sire,genoanimped$Dam,sep='_')
    genoanimped$famIDFSseqID <- unclass(as.factor(genoanimped$famIDFS))
    genoanimped$famIDseqID <- unclass(as.factor(genoanimped$famID))
    genoanimped <- genoanimped[order(genoanimped$famIDFSseqID),]
    cat('\n... families created using pedigree information ...\n')
    usefam <- sort(table(genoanimped$famIDseqID))
    allfam <- length(usefam)
    usefam <- sort(as.numeric(as.vector(names(usefam[usefam>exclfamsize]))))
    if((allfam-length(usefam))<=exclfamsize){
      cat('... ',allfam-length(usefam),' family(s) had only ',exclfamsize,' offspring, thus discarded ...\n')} else {
        cat('... ',allfam-length(usefam),' families had only one offspring, thus discarded ...\n')
      }
    cat('\n... Extracting familial relationships ...\n')
    
    genoanimped <- genoanimped[genoanimped$famIDseqID %in% usefam,]
    iterchecks <- round(length(usefam)/5,digits=0)
    
    HSgrm <- data.frame(IID1=character(),IID2=character(),G=numeric(),seqfamID=character())
    for(f in 1:length(usefam)){
      halfsibsg <- genoanimped[genoanimped$famIDseqID %in% usefam[f],]
      halfsibsgID <- halfsibsg[,c('seqID','famIDFSseqID')]
      if(length(unique(halfsibsg$famIDFSseqID))>1){
        halfsibsg <- expand.grid(as.vector(halfsibsg$seqID),as.vector(halfsibsg$seqID),stringsAsFactors=F)
        halfsibsg <- halfsibsg[which(halfsibsg$Var1!=halfsibsg$Var2),c(2,1)]
        halfsibsg <- halfsibsg[order(halfsibsg$Var2,halfsibsg$Var1),]
        halfsibsg <- halfsibsg[which(halfsibsg$Var2<halfsibsg$Var1),]
        halfsibsg <- merge(halfsibsg,halfsibsgID,by.x=1,by.y=1)
        halfsibsg <- merge(halfsibsg,halfsibsgID,by.x=2,by.y=1)
        halfsibsg <- halfsibsg[order(halfsibsg$famIDFSseqID.x,halfsibsg$famIDFSseqID.y),]
        halfsibsg <- halfsibsg[which(halfsibsg$famIDFSseqID.x!=halfsibsg$famIDFSseqID.y),c(2,1)]
        halfsibsg <- halfsibsg[order(halfsibsg$Var2,halfsibsg$Var1),]
        
        for(fsize in 1:nrow(halfsibsg)){
          halfsibgrm <- grm[which((grm$IID1==halfsibsg$Var2[fsize] & grm$IID2==halfsibsg$Var1[fsize]) | 
                                    (grm$IID2==halfsibsg$Var2[fsize] & grm$IID1==halfsibsg$Var1[fsize])),]
          halfsibgrm$famIDseqID <- usefam[f]
          HSgrm <- rbind.data.frame(HSgrm,halfsibgrm,stringsAsFactors=F)
        }
        if(f %% iterchecks==0){cat('... family ...',f,' (out of ',length(usefam),') familial checks ...\n')}
        rm(halfsibgrm,halfsibsg)
      }
    }
    cat('... ',length(unique(HSgrm$famIDseqID)),'Half-sib families (sire families) obtained ...\n')
    
    plotname <-paste(outname,'.tiff',sep='')
    dpi=300;width.cm<-18;height.cm<-13;width.in<-width.cm/2.54;height.in<-height.cm/2.54
    tiff(file=plotname,width=width.in*dpi,height=height.in*dpi,pointsize=10,units="px",res=dpi)
    layout(matrix(1:2,1,2))
    par(mar=c(4,4,0.1,0.1))
    xlim=range(HSgrm$G)+c(-0.20,0.20)
    hist(HSgrm$G,breaks=100,col=sample(colors(),1,replace=F),xlim=xlim,
         xlab='Half-sib genomic relationships',main='')
    abline(v=mean(HSgrm$G),col='red',lwd=2)
    abline(v=0.25,col='green',lwd=1)
    legend('topleft',legend=c(paste('Average G (',round(mean(HSgrm$G),3),')',sep=''),
                              'Pedigre-based expectation (= 0.250)'),
           col=c('red','green'),lty=1,lwd=1,cex=0.6,bty='n')
    par(mar=c(4,3.8,0.1,0.5))
    plot(density(HSgrm$G),col=sample(colors(),1,replace=F),xlim=xlim,
         xlab='Half-sib genomic relationships',main='',type='l',lwd=1.5)
    points(density(HSgrm$G,adjust=2),col=sample(colors(),1,replace=F),xlim=xlim,
           xlab='Half-sib genomic relationships',main='',type='l',lty=2)
    dev.off()
    write.csv(HSgrm,paste(outname,'_FSrelationships.csv',sep=''),quote=F,row.names=F)
    return(HSgrm)
  }
}
