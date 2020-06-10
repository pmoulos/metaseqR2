read2count <- function(targets,annotation,fileType=targets$type,
    transLevel="gene",utrOpts=list(frac=1,minLength=300,downstream=50),
    interFeature=FALSE,rc=NULL) {
    if (missing(targets))
        stopwrap("You must provide the targets argument!")
    if (missing(annotation))
        stopwrap("You must provide an annotation GenomicRanges or data frame!")
    if (!requireNamespace("GenomicRanges"))
        stopwrap("The Bioconductor package GenomicRanges is required to ",
            "proceed!")
    if (fileType=="bed" && !requireNamespace("rtracklayer"))
        stopwrap("The Bioconductor package rtracklayer is required to process ",
            "BED files!")
    if (fileType %in% c("sam","bam")) {
        if (!requireNamespace("Rsamtools"))
            stopwrap("The Bioconductor package Rsamtools is required to ",
                "process BAM files!")
        if (!requireNamespace("GenomeInfoDb"))
            stopwrap("The Bioconductor package GenomeInfoDb is required to ",
                "continue with BAM files!")
    }
    if (!is.list(targets) && file.exists(targets))
        targets <- readTargets(targets)
    else if (!is.list(targets) && !file.exists(targets))
        stopwrap("You must provide a targets list or a valid targets file!")
    
    # Convert annotation to GRanges. If there are no sufficient columns to
    # create GRanges, then it will fail anyway. All the reducing/merging
    # procedures, would have already taken place.
    if (!is(annotation,"GenomicRanges"))
        annotationGr <- GRanges(annotation)
    else
        annotationGr <- annotation
    
    # Determine internal count type
    if (length(grep("MET",annotation$transcript_id[1])) > 0
        || length(grep("MEU",annotation$transcript_id[1])) > 0)
        countType <- "utr"
    else if (length(grep("MEX",annotation$exon_id[1])) > 0
        || length(grep("MTE",annotation$exon_id[1])) > 0)
        countType <- "exon"
    else
        countType <- "gene"
    
    # If the count type is "exon", we must reduce overlapping exons belonging to
    # multiple transcripts, so as to avoid inflating the final read count when
    # summing all exons. Included here mostly for backwards compatibility, until
    # new version stabilizes.
    backList <- .backPreprocessAnnotation(annotation,annotationGr,transLevel)
    annotationGr <- backList$annotationGr
    mergedAnnotation <- backList$mergedAnnotation
    
    if (countType == "utr" && !is.null(utrOpts)) {
        disp("Resizing transcript 3' UTRs...")
        w <- width(annotationGr)
        # 1. Percentage of 3' UTR only (upstream of end), ensuring min length
        nw <- utrOpts$frac*width(annotationGr)
        nw <- ifelse(nw < utrOpts$minLength,utrOpts$minLength,nw)
        annotationGr <- resize(annotationGr,width=nw,fix="end")
        # 2. Extend 3' UTR downstream
        annotationGr <- resize(annotationGr,
            width=width(annotationGr) + utrOpts$downstream,fix="start")
    }
    
    # Continue
    filesList <- targets$files
    sampleNames <- unlist(lapply(filesList,names),use.names=FALSE)
    sampleFiles <- unlist(filesList,use.names=FALSE)
    names(sampleFiles) <- sampleNames
    if (!is.null(targets$paired)) {
        paired <- unlist(targets$paired,use.names=FALSE)
        names(paired) <- sampleNames
    }
    else
        paired <- NULL
    if (!is.null(targets$stranded)) {
        stranded <- unlist(targets$stranded,use.names=FALSE)
        names(stranded) <- sampleNames
    }
    else
        stranded <- NULL
        
    # Initialize libsize
    libsize <- vector("list",length(sampleNames))
    names(libsize) <- sampleNames
    
    if (fileType=="bed") {
        # Initialize counts
        counts <- matrix(0,nrow=length(annotationGr),ncol=length(sampleNames))
        rownames(counts) <- names(annotationGr)
        colnames(counts) <- sampleNames
        
        retVal <- cmclapply(sampleNames,function(n,sampleFiles) {
            disp("Reading bed file ",basename(sampleFiles[n]),
                " for sample with name ",n,". This might take some time...")
            bed <- import.bed(sampleFiles[n],trackLine=FALSE)
            disp("  Checking for chromosomes not present in the annotation...")
            bed <- bed[which(!is.na(match(as(seqnames(bed),"character"),
                seqlevels(annotationGr))))]
            libsize <- length(bed)
            if (length(bed) > 0) {
                disp("  Counting reads overlapping with given annotation...")
                counts <- countOverlaps(annotationGr,bed)
            }
            else
                warnwrap(paste("No reads left after annotation chromosome ",
                    "presence check for sample ",n,sep=""))
            gc(verbose=FALSE)
            return(list(counts=counts,libsize=libsize))
        },sampleFiles,rc=rc)
        disp("  Finished counting!")
    }
    else if (fileType %in% c("sam","bam")) {
        if (fileType=="sam")
            sampleFiles <- .convertSam(sampleFiles)
        
        # Adjust seqlevelsStyle locally, only for counting, the finally reported
        # style is UCSC as several downstream steps depend on it
        annotationGr <- .adjustSeqlevelsStyle(annotationGr,sampleFiles[1])
        
        # We must now compare the Seqinfos of annotation and a BAM file to check
        # whether they match. If they do not match, then for now, we drop the
        # Seqinfo of annotation. In a future release, it will be converted to
        # the BAMs Seqinfo derived from the header.
        annotationGr <- .adjustOrDropSeqlevels(annotationGr,sampleFiles[1])
        
        # Initialize counts, after fixing annotationGr
        counts <- matrix(0,nrow=length(annotationGr),ncol=length(sampleNames))
        rownames(counts) <- names(annotationGr)
        colnames(counts) <- sampleNames
        
        retVal <- cmclapply(sampleNames,function(n,sampleFiles,paired,
            stranded) {
            disp("Reading bam file ",basename(sampleFiles[n])," for sample ",
                "with name ",n,". This might take some time...")
            if (!is.null(paired)) {
                p <- tolower(paired[n])
                if (p=="single") {
                    singleEnd <- TRUE
                    fragments <- FALSE
                    asMates <- FALSE
                }
                else if (p=="paired") {
                    singleEnd <- FALSE
                    fragments <- FALSE
                    asMates <- TRUE
                }
                else if (p=="mixed") {
                    singleEnd <- FALSE
                    fragments <- TRUE
                    asMates <- TRUE
                }
                else {
                    warnwrap("Information regarding single- or paired-end ",
                        "reads is not correctly provided! Assuming single...")
                    singleEnd <- TRUE
                    fragments <- FALSE
                    asMates <- FALSE
                }
            }
            else {
                singleEnd <- TRUE
                fragments <- FALSE
                asMates <- FALSE
            }
            if (!is.null(stranded)) {
                s <- tolower(stranded[n])
                if (s %in% c("forward","reverse"))
                    ignoreStrand <- FALSE
                else if (s=="no")
                    ignoreStrand <- TRUE
                else {
                    warnwrap("Information regarding strandedness of the reads ",
                        "is not correctly provided! Assuming unstranded...")
                    ignoreStrand <- TRUE
                }
            }
            else
                ignoreStrand <- TRUE
            # Check remoteness
            if (length(grep("^(http|ftp)",sampleFiles[n],perl=TRUE))>=1) {
                reads <- as(readGAlignments(file=sampleFiles[n]),"GRanges")
                libsize <- length(reads)
                isRemote <- TRUE
            }
            else {
                reads <- BamFile(sampleFiles[n],asMates=asMates)
                libsize <- countBam(reads,
                param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE)))$records
                isRemote <- FALSE
            }
            if (libsize>0) {
                disp("  Counting reads overlapping with given annotation...")
                if (singleEnd & !fragments)
                    disp("    ...for single-end reads...")
                else if (!singleEnd & !fragments)
                    disp("    ...for paired-end reads...")
                else if (!singleEnd & fragments)
                    disp("    ...for mixed single- and paired-end reads...")
                if (ignoreStrand)
                    disp("    ...ignoring strandedness...")
                else {
                    disp("    ...assuming ",s," sequenced reads...")
                    if (s=="reverse")
                        strand(annotationGr) <- ifelse(strand(
                            annotationGr)=="+","-","+")
                }
                if (isRemote)
                    disp("    ...for remote BAM file... might take longer...")
                counts <- summarizeOverlaps(annotationGr,reads,
                    singleEnd=singleEnd,fragments=fragments,
                    ignore.strand=ignoreStrand,inter.feature=interFeature)
                counts <- assays(counts)$counts
            }
            else
                warnwrap(paste("No reads left after annotation chromosome ",
                    "presence check for sample ",n,sep=""))
            gc(verbose=FALSE)
            return(list(counts=counts,libsize=libsize))
        },sampleFiles,paired,stranded,rc=rc)
        disp("  Finished counting!")
    }
    for (i in seq_along(retVal)) {
        counts[,i] <- retVal[[i]]$counts
        libsize[[i]] <- retVal[[i]]$libsize
    }
    
    return(list(counts=counts,libsize=libsize,mergedann=mergedAnnotation))
}

readTargets <- function(input,path=NULL) {
    if (missing(input) || !file.exists(input))
        stopwrap("The targets file should be a valid existing text file!")
    tab <- read.delim(input,strip.white=TRUE)
    samples <- as.character(tab[,1])
    conditions <- unique(as.character(tab[,3]))
    rawfiles <- as.character(tab[,2])
    if (!is.null(path)) {
        tmp <- dirname(rawfiles) # Test if there is already a path
        if (any(tmp=="."))
            rawfiles <- file.path(path,basename(rawfiles))
    }
    # Check if they exist!!!
    for (f in rawfiles) {
        if (!file.exists(f))
            stopwrap("Raw reads input file ",f," does not exist! Please check!")
    }
    if (length(samples) != length(unique(samples)))
        stopwrap("Sample names must be unique for each sample!")
    if (length(rawfiles) != length(unique(rawfiles)))
        stopwrap("File names must be unique for each sample!")
    sampleList <- vector("list",length(conditions))
    names(sampleList) <- conditions
    for (n in conditions)
        sampleList[[n]] <- samples[which(as.character(tab[,3])==n)]
    fileList <- vector("list",length(conditions))
    names(fileList) <- conditions
    for (n in conditions) {
        fileList[[n]] <- rawfiles[which(as.character(tab[,3])==n)]
        names(fileList[[n]]) <- samples[which(as.character(tab[,3])==n)]
    }
    if (ncol(tab)>3) { # Has info about single- or paired-end reads / strand
        if (ncol(tab)==4) { # Stranded or paired
            whats <- tolower(as.character(tab[,4]))
            if (!all(whats %in% c("yes","no","forward","reverse",
                "single","paired")))
                stopwrap("Unknown options for paired-end reads and/or ",
                    "strandedness in targets file.")
            what <- whats[1]
            if (what %in% c("single","paired")) {
                hasPairedInfo <- TRUE
                hasStrandedInfo <- FALSE
            }
            else {
                if (what %in% c("yes","no")) {
                    .deprecatedWarning("readTargets")
                    tmp <- as.character(tab[,4])
                    tmp[tmp=="yes"] <- "forward"
                    tab[,4] <- tmp
                    hasPairedInfo <- FALSE
                    hasStrandedInfo <- TRUE
                }
                if (what %in% c("forward","reverse","no")) {
                    hasPairedInfo <- FALSE
                    hasStrandedInfo <- TRUE
                }
            }
        }
        if (ncol(tab)==5) { # Both
            whatsPaired <- tolower(as.character(tab[,4]))
            if (!all(whatsPaired %in% c("single","paired","mixed")))
                stopwrap("Unknown option for type of reads (single, paired, ",
                    "mixed) in targets file.")
            whatsStrand <- tolower(as.character(tab[,5]))
            if (!all(whatsStrand %in% c("yes","no","forward","reverse")))
                stopwrap("Unknown option for read strandedness in targets file")
            if (any(whatsStrand=="yes")) {
                .deprecatedWarning("readTargets")
                tmp <- as.character(tab[,5])
                tmp[tmp=="yes"] <- "forward"
                tab[,5] <- tmp
            }
            hasPairedInfo <- TRUE
            hasStrandedInfo <- TRUE
        }
        if (hasPairedInfo && !hasStrandedInfo) {
            pairedList <- vector("list",length(conditions))
            names(pairedList) <- conditions
            for (n in conditions) {
                pairedList[[n]] <- character(length(sampleList[[n]]))
                names(pairedList[[n]]) <- sampleList[[n]]
                for (nn in names(pairedList[[n]]))
                    pairedList[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
        else
            pairedList <- NULL
        if (hasStrandedInfo && !hasPairedInfo) {
            strandedList <- vector("list",length(conditions))
            names(strandedList) <- conditions
            for (n in conditions) {
                strandedList[[n]] <- character(length(sampleList[[n]]))
                names(strandedList[[n]]) <- sampleList[[n]]
                for (nn in names(strandedList[[n]]))
                    strandedList[[n]][nn] <- 
                        as.character(tab[which(as.character(tab[,1])==nn),4])
            }
        }
        else
            strandedList <- NULL
        if (hasStrandedInfo && hasPairedInfo) {
            strandedList <- vector("list",length(conditions))
            names(strandedList) <- conditions
            for (n in conditions) {
                strandedList[[n]] <- character(length(sampleList[[n]]))
                names(strandedList[[n]]) <- sampleList[[n]]
                for (nn in names(strandedList[[n]]))
                    strandedList[[n]][nn] <- as.character(tab[which(
                        as.character(tab[,1])==nn),5])
            }
            pairedList <- vector("list",length(conditions))
            names(pairedList) <- conditions
            for (n in conditions) {
                pairedList[[n]] <- character(length(sampleList[[n]]))
                names(pairedList[[n]]) <- sampleList[[n]]
                for (nn in names(pairedList[[n]]))
                    pairedList[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
    }
    else
        pairedList <- strandedList <- NULL
    # Guess file type based on only one of them
    tmp <- fileList[[1]][1]
    if (length(grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "bam"
    else if (length(grep("\\.sam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "sam"
    else if (length(grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0))
        type <- "bed"
    else
        type <- NULL
    return(list(samples=sampleList,files=fileList,paired=pairedList,
        stranded=strandedList,type=type))
}

.writeTargets <- function(targets,outfile) {
    out <- data.frame(
        samplename=unlist(targets$samples,use.names=FALSE),
        filename=unlist(targets$files,use.names=FALSE),
        condition=rep(names(targets$samples),lengths(targets$samples))
    )
    if (!is.null(targets$paired))
        out$paired <- unlist(targets$paired,use.names=FALSE)
    if (!is.null(targets$stranded))
        out$stranded <- unlist(targets$stranded,use.names=FALSE)
    write.table(out,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
    return(out)
}

.convertSam <- function(sampleFiles) {
    for (n in names(sampleFiles)) {
        dest <- file.path(dirname(sampleFiles[n]),n)
        disp("Converting sam file ",basename(sampleFiles[n]),
            " to bam file ",basename(dest),"...")
        asBam(file=sampleFiles[n],destination=dest,overwrite=TRUE)
        sampleFiles[n] <- paste(dest,"bam",sep=".")
    }
    return(sampleFiles)
}

.adjustOrDropSeqlevels <- function(gr,f) {
    ciAnn <- .chromInfoFromSeqinfo(seqinfo(gr))
    ciBam <- .chromInfoFromBAM(f)
    rows <- intersect(rownames(ciAnn),rownames(ciBam))
    norows <- setdiff(rownames(ciAnn),rownames(ciBam))
    
    if (length(norows) > 0)
        gr <- dropSeqlevels(gr,value=norows,pruning.mode="coarse")
    
    if (length(rows) > 0) {
        ciAnn <- ciAnn[rows,,drop=FALSE]
        ciBam <- ciBam[rows,,drop=FALSE]
        if (!identical(ciAnn,ciBam)) # Different lengths, drop from annotation
            seqlengths(gr) <- rep(NA,length(seqlengths(gr)))
    }
    
    return(gr)
}

.adjustSeqlevelsStyle <- function(gr,f) {
    b <- BamFile(f)
    
    # Check seqlevels styles of the loaded genome for local counting use and
    # change the style of the local annotation.
    sb <- seqlevelsStyle(b)
    sg <- seqlevelsStyle(gr)
    if (!any(sg %in% sb))
        seqlevelsStyle(gr) <- sb[1]
    
    return(gr)
}

.backPreprocessAnnotation <- function(annotation,annotationGr,transLevel,
    rc=NULL) {
    if (length(grep("exon",names(mcols(annotationGr)))) > 0) {
        # countType is exon
        if (length(grep("MEX",annotation$exon_id[1]))) # Retrieved from previous
            mergedAnnotation <- annotation
        else {
            if (transLevel=="gene") {
                disp("Merging exons to create unique gene models...")
                annotationGr <- .reduceExonsOld(annotationGr)
            }
            #else if (transLevel=="transcript") {
            #   disp("Merging exons to create unique gene models...")
            #   annotationGr <- reduceExonsTranscript(annotationGr,rc=rc)
            #}
            #mergedAnnotation <- as.data.frame(annotationGr) # Bug?
            mergedAnnotation <- data.frame(
                chromosome=as.character(seqnames(annotationGr)),
                start=start(annotationGr),
                end=end(annotationGr),
                exon_id=if (!is.null(annotationGr$exon_id))
                    as.character(annotationGr$exon_id) else
                    as.character(annotationGr$name),
                gene_id=if (!is.null(annotationGr$gene_id))
                    as.character(annotationGr$gene_id) else
                    as.character(annotationGr$name),
                strand=as.character(strand(annotationGr)),
                gene_name=if (!is.null(annotationGr$gene_name))
                    as.character(annotationGr$gene_name) else 
                    if (!is.null(annotationGr$symbol))
                    as.character(annotationGr$name) else NULL,
                biotype=if (!is.null(annotationGr$biotype))
                    as.character(annotationGr$biotype) else NULL
            )
            rownames(mergedAnnotation) <- 
                as.character(mergedAnnotation$exon_id)
            names(annotationGr) <- annotationGr$exon_id
        }
    }
    else if (length(grep("transcript",names(mcols(annotationGr))))>0) {
        # countType may be utr
        if (length(grep("MET",annotation$transcript_id[1]))
            || length(grep("MEU",annotation$transcript_id[1])))
            # Retrieved from previous
            mergedAnnotation <- annotation
        else {
            if (transLevel=="gene") {
                disp("Merging transcript 3' UTRs to create unique ",
                    "gene models...")
                annotationGr <- .reduceTranscriptsOld(annotationGr,rc=rc)
            }
            if (transLevel=="transcript") {
                disp("Merging transcript 3' UTRs to create unique ",
                    "transcript models...")
                annotationGr <- 
                    .reduceTranscriptsUtrOld(annotationGr,rc=rc)
            }
            
            #mergedAnnotation <- as.data.frame(annotationGr) # Bug?
            mergedAnnotation <- data.frame(
                chromosome=as.character(seqnames(annotationGr)),
                start=start(annotationGr),
                end=end(annotationGr),
                transcript_id=if (!is.null(annotationGr$transcript_id))
                    as.character(annotationGr$transcript_id) else
                    as.character(annotationGr$name),
                gene_id=if (!is.null(annotationGr$gene_id))
                    as.character(annotationGr$gene_id) else
                    as.character(annotationGr$name),
                strand=as.character(strand(annotationGr)),
                gene_name=if (!is.null(annotationGr$gene_name))
                    as.character(annotationGr$gene_name) else 
                    if (!is.null(annotationGr$symbol))
                    as.character(annotationGr$name) else NULL,
                biotype=if (!is.null(annotationGr$biotype))
                    as.character(annotationGr$biotype) else NULL
            )
            rownames(mergedAnnotation) <- 
                as.character(mergedAnnotation$transcript_id)
            names(annotationGr) <- annotationGr$transcript_id
        }
        interFeature = FALSE # Quant-Seq
    }
    else
        mergedAnnotation <- NULL
    return(list(annotationGr=annotationGr,mergedAnnotation=mergedAnnotation))
}

.reduceExonsOld <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            exon_id=paste(x,"MEX",seq_len(n),sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    len <- unlist(cmclapply(redList,function(x) {
        return(sum(width(x)))
    },rc=rc))
    names(len) <- names(redList)
    return(list(model=do.call("c",redList),length=len))
}

.reduceTranscriptsOld <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MET",seq_len(n),sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    return(do.call("c",redList))
}

.reduceTranscriptsUtrOld <- function(gr,rc=NULL) {
    trans <- unique(as.character(gr$transcript_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(trans,function(x,a,g,b) {
        tmp <- a[a$transcript_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MEU",seq_len(n),sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    return(do.call("c",redList))
}

################################################################################

#read2countOld <- function(targets,annotation,fileType=targets$type,
#    transLevel="gene",utrFlank=500,interFeature=FALSE,rc=NULL) {
#    if (missing(targets))
#        stopwrap("You must provide the targets argument!")
#    if (missing(annotation))
#        stopwrap("You must provide an annotation data frame!")
#    if (!require(GenomicRanges))
#        stopwrap("The Bioconductor package GenomicRanges is required to ",
#            "proceed!")
#    if (fileType=="bed" && !require(rtracklayer))
#        stopwrap("The Bioconductor package rtracklayer is required to ",
#            "process BED files!")
#    if (fileType %in% c("sam","bam")) {
#        if (!require(Rsamtools))
#            stopwrap("The Bioconductor package Rsamtools is required to ",
#                "process BAM files!")
#    }
#    if (!is.list(targets) && file.exists(targets))
#        targets <- readTargets(targets)
#    else if (!is.list(targets) && !file.exists(targets))
#        stopwrap("You must provide a targets list or a valid targets file!")
    
#    # Convert annotation to GRanges
#    annotationGr <- makeGRangesFromDataFrame(
#        df=annotation,
#        keep.extra.columns=TRUE,
#        seqnames.field="chromosome"
#    )
#    # annotationGr <- GRanges(annotation) also works for our structure
    
#    # If the count type is "exon", we must reduce overlapping exons belonging 
#    # to multiple transcripts, so as to avoid inflating the final read count 
#    # when summing all exons
#    if (length(grep("exon",colnames(annotation)))>0) { # countType is exon
#        # Retrieved from previous
#        if (length(grep("MEX",annotation$exon_id[1]))) 
#            mergedAnnotation <- annotation
#        else {
#            if (transLevel=="gene") {
#                disp("Merging exons to create unique gene models...")
#                annotationGr <- reduceExons(annotationGr,rc=rc)
#            }
#            #else if (transLevel=="transcript") {
#            #   disp("Merging exons to create unique gene models...")
#            #   annotationGr <- reduceExonsTranscript(annotationGr,rc=rc)
#            #}
#            #mergedAnnotation <- as.data.frame(annotationGr) # Bug?
#            mergedAnnotation <- data.frame(
#                chromosome=as.character(seqnames(annotationGr)),
#                start=start(annotationGr),
#                end=end(annotationGr),
#                exon_id=if (!is.null(annotationGr$exon_id))
#                    as.character(annotationGr$exon_id) else
#                    as.character(annotationGr$name),
#                gene_id=if (!is.null(annotationGr$gene_id))
#                    as.character(annotationGr$gene_id) else
#                    as.character(annotationGr$name),
#                strand=as.character(strand(annotationGr)),
#                gene_name=if (!is.null(annotationGr$gene_name))
#                    as.character(annotationGr$gene_name) else 
#                    if (!is.null(annotationGr$symbol))
#                    as.character(annotationGr$name) else NULL,
#                biotype=if (!is.null(annotationGr$biotype))
#                    as.character(annotationGr$biotype) else NULL
#            )
#            rownames(mergedAnnotation) <- 
#                as.character(mergedAnnotation$exon_id)
#        }
#    }
#    else if (length(grep("transcript",colnames(annotation)))>0) {
#        # countType may be utr
#        if (length(grep("MET",annotation$transcript_id[1]))
#            || length(grep("MEU",annotation$transcript_id[1]))) 
#            # Retrieved from previous
#            mergedAnnotation <- annotation
#        else {
#            if (transLevel=="gene") {
#                disp("Merging transcript 3' UTRs to create unique ",
#                    "gene models...")
#                annotationGr <- reduceTranscriptsUtr(annotationGr,rc=rc)
#            }
#            if (transLevel=="transcript") {
#                disp("Merging transcript 3' UTRs to create unique ",
#                    "transcript models...")
#                annotationGr <- 
#                    reduceTranscriptsUtrTranscript(annotationGr,rc=rc)
#            }
#            if (utrFlank > 0) {
#                disp("Flanking merged transcript 3' UTRs per ",utrFlank,
#                    "bp...")
#                w <- width(annotationGr)
#                annotationGr <- promoters(annotationGr,upstream=utrFlank,
#                    downstream=0)
#                annotationGr <- resize(annotationGr,width=w+2*utrFlank)
#            }
#            #mergedAnnotation <- as.data.frame(annotationGr) # Bug?
#            mergedAnnotation <- data.frame(
#                chromosome=as.character(seqnames(annotationGr)),
#                start=start(annotationGr),
#                end=end(annotationGr),
#                transcript_id=if (!is.null(annotationGr$transcript_id))
#                    as.character(annotationGr$transcript_id) else
#                    as.character(annotationGr$name),
#                gene_id=if (!is.null(annotationGr$gene_id))
#                    as.character(annotationGr$gene_id) else
#                    as.character(annotationGr$name),
#                strand=as.character(strand(annotationGr)),
#                gene_name=if (!is.null(annotationGr$gene_name))
#                    as.character(annotationGr$gene_name) else 
#                    if (!is.null(annotationGr$symbol))
#                    as.character(annotationGr$name) else NULL,
#                biotype=if (!is.null(annotationGr$biotype))
#                    as.character(annotationGr$biotype) else NULL
#            )
#            rownames(mergedAnnotation) <- 
#                as.character(mergedAnnotation$transcript_id)
#        }
#        interFeature = FALSE # Quant-Seq
#    }
#    else
#        mergedAnnotation <- NULL
#    # Continue
#    filesList <- targets$files
#    sampleNames <- unlist(lapply(filesList,names),use.names=FALSE)
#    sampleFiles <- unlist(filesList,use.names=FALSE)
#    names(sampleFiles) <- sampleNames
#    if (!is.null(targets$paired)) {
#        paired <- unlist(targets$paired,use.names=FALSE)
#        names(paired) <- sampleNames
#    }
#    else
#        paired <- NULL
#    if (!is.null(targets$stranded)) {
#        stranded <- unlist(targets$stranded,use.names=FALSE)
#        names(stranded) <- sampleNames
#    }
#    else
#        stranded <- NULL
#    counts <- matrix(0,nrow=length(annotationGr),ncol=length(sampleNames))
#    if (length(grep("exon",colnames(annotation)))>0)
#        rownames(counts) <- as.character(annotationGr$exon_id)
#    else if (length(grep("transcript",colnames(annotation)))>0)
#        rownames(counts) <- as.character(annotationGr$transcript_id)
#    else
#        rownames(counts) <- as.character(annotationGr$gene_id)
#    colnames(counts) <- sampleNames
#    libsize <- vector("list",length(sampleNames))
#    names(libsize) <- sampleNames
#    if (fileType=="bed") {
#        retVal <- cmclapply(sampleNames,function(n,sampleFiles) {
#            disp("Reading bed file ",basename(sampleFiles[n]),
#                " for sample with name ",n,". This might take some time...")
#            bed <- import.bed(sampleFiles[n],trackLine=FALSE)
#            disp("  Checking for chromosomes not present in the annotation...")
#            bed <- bed[which(!is.na(match(as(seqnames(bed),"character"),
#                seqlevels(annotationGr))))]
#            libsize <- length(bed)
#            if (length(bed)>0) {
#                disp("  Counting reads overlapping with given annotation...")
#                counts <- countOverlaps(annotationGr,bed)
#            }
#            else
#                warnwrap(paste("No reads left after annotation chromosome ",
#                    "presence check for sample ",n,sep=""))
#            gc(verbose=FALSE)
#            return(list(counts=counts,libsize=libsize))
#        },sampleFiles,rc=rc)
#    }
#    else if (fileType %in% c("sam","bam")) {
#        if (fileType=="sam") {
#            for (n in sampleNames) {
#                dest <- file.path(dirname(sampleFiles[n]),n)
#                disp("Converting sam file ",basename(sampleFiles[n]),
#                    " to bam file ",basename(dest),"...")
#                asBam(file=sampleFiles[n],destination=dest,overwrite=TRUE)
#                sampleFiles[n] <- paste(dest,"bam",sep=".")
#            }
#        }
#        retVal <- cmclapply(sampleNames,function(n,sampleFiles,paired,
#            stranded) {
#            disp("Reading bam file ",basename(sampleFiles[n])," for sample ",
#                "with name ",n,". This might take some time...")
#            if (!is.null(paired)) {
#                p <- tolower(paired[n])
#                if (p=="single") {
#                    singleEnd <- TRUE
#                    fragments <- FALSE
#                    asMates <- FALSE
#                }
#                else if (p=="paired") {
#                    singleEnd <- FALSE
#                    fragments <- FALSE
#                    asMates <- TRUE
#                }
#                else if (p=="mixed") {
#                    singleEnd <- FALSE
#                    fragments <- TRUE
#                    asMates <- TRUE
#                }
#                else {
#                    warnwrap("Information regarding single- or paired-end ",
#                        "reads is not correctly provided! Assuming single...")
#                    singleEnd <- TRUE
#                    fragments <- FALSE
#                    asMates <- FALSE
#                }
#            }
#            else {
#                singleEnd <- TRUE
#                fragments <- FALSE
#                asMates <- FALSE
#            }
#            if (!is.null(stranded)) {
#                s <- tolower(stranded[n])
#                if (s %in% c("forward","reverse"))
#                    ignoreStrand <- FALSE
#                else if (s=="no")
#                    ignoreStrand <- TRUE
#                else {
#                   warnwrap("Information regarding strandedness of the reads ",
#                        "is not correctly provided! Assuming unstranded...")
#                    ignoreStrand <- TRUE
#                }
#            }
#            else
#                ignoreStrand <- TRUE
#            # Check remoteness
#            if (length(grep("^(http|ftp)",sampleFiles[n],perl=TRUE))>=1) {
#                reads <- as(readGAlignments(file=sampleFiles[n]),"GRanges")
#                libsize <- length(reads)
#                isRemote <- TRUE
#            }
#            else {
#                reads <- BamFile(sampleFiles[n],asMates=asMates)
#                libsize <- countBam(reads,
#                param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE)))$records
#                isRemote <- FALSE
#            }
#            if (libsize>0) {
#                disp("  Counting reads overlapping with given annotation...")
#                if (singleEnd & !fragments)
#                    disp("    ...for single-end reads...")
#                else if (!singleEnd & !fragments)
#                    disp("    ...for paired-end reads...")
#                else if (!singleEnd & fragments)
#                    disp("    ...for mixed single- and paired-end reads...")
#                if (ignoreStrand)
#                    disp("    ...ignoring strandedness...")
#                else {
#                    disp("    ...assuming ",s," sequenced reads...")
#                    if (s=="reverse")
#                        strand(annotationGr) <- ifelse(strand(
#                            annotationGr)=="+","-","+")
#                }
#                if (isRemote)
#                    disp("    ...for remote BAM file... might take longer...")
#                counts <- summarizeOverlaps(annotationGr,reads,
#                    singleEnd=singleEnd,fragments=fragments,
#                    ignore.strand=ignoreStrand,inter.feature=interFeature)
#                counts <- assays(counts)$counts
#            }
#            else
#                warnwrap(paste("No reads left after annotation chromosome ",
#                    "presence check for sample ",n,sep=""))
#            gc(verbose=FALSE)
#            return(list(counts=counts,libsize=libsize))
#        },sampleFiles,paired,stranded,rc=rc)
#    }
#    for (i in seq_along(retVal)) {
#        counts[,i] <- retVal[[i]]$counts
#        libsize[[i]] <- retVal[[i]]$libsize
#    }
    
#    return(list(counts=counts,libsize=libsize,mergedann=mergedAnnotation))
#}
