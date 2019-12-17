buildExport <- function(geneData,rawGeneCounts,normGeneCounts,flags,
    sampleList,cnt,statistics,
    rawList,normList,
    pMat=if (is(geneData,"GenomicRanges")) 
        matrix(NA,length(geneData),length(statistics))
        else matrix(NA,nrow(geneData),length(statistics)),
    adjpMat=if (is(geneData,"GenomicRanges"))
        matrix(NA,length(geneData),length(statistics))
        else matrix(NA,nrow(geneData),length(statistics)),
    sumP=if (is(geneData,"GenomicRanges")) rep(NA,length(geneData))
        else rep(NA,nrow(geneData)),
    adjSumP=if (is(geneData,"GenomicRanges")) rep(NA,length(geneData))
        else rep(NA,nrow(geneData)),
    exportWhat=c("annotation","p_value","adj_p_value","meta_p_value",
        "adj_meta_p_value","fold_change","stats","counts","flags"),
    exportScale=c("natural","log2","log10","rpgm","vst"),
    exportValues=c("raw","normalized"),
    exportStats=c("mean","median","sd","mad","cv","rcv"),
    logOffset=1,report=TRUE) {
    if (is(geneData,"GenomicRanges")) {
        geneData <- as.data.frame(geneData)
        gci <- grep("gc_content",colnames(geneData))
        nmi <- grep("gene_name",colnames(geneData))
        bti <- grep("biotype",colnames(geneData))
        baseInd <- c(1:3,6)
        if (length(gci) > 0)
            baseInd <- c(baseInd,gci,5)
        else
            baseInd <- c(baseInd,5)
        if (length(nmi) > 0)
            baseInd <- c(baseInd,nmi)
        if (length(bti) > 0)
            baseInd <- c(baseInd,bti)
        geneData <- geneData[,baseInd]
    }
        
    if (is.null(colnames(pMat)))
        colnames(pMat) <- statistics
    if (is.null(adjpMat))
        adjpMat=matrix(NA,nrow(geneData),length(statistics))
    if (is.null(colnames(adjpMat)))
        colnames(adjpMat) <- statistics
    
    export <- data.frame(row.names=rownames(geneData))
    if (report)
        exportHtml <- as.matrix(export)
    theNames <- character(0)
    if ("annotation" %in% exportWhat) {
        disp("      binding annotation...")
        export <- cbind(export,geneData)
        if (report)
            exportHtml <- cbind(exportHtml,.makeHtmlCells(geneData,
                type="text"))
        theNames <- c(theNames,colnames(geneData))
    }
    if ("p_value" %in% exportWhat) {
        disp("      binding p-values...")
        export <- cbind(export,pMat)
        if (report) 
            exportHtml <- cbind(exportHtml,.makeHtmlCells(pMat))
        theNames <- c(theNames,paste("p-value_",colnames(pMat),sep=""))
    }
    if ("adj_p_value" %in% exportWhat) {
        disp("      binding FDRs...")
        export <- cbind(export,adjpMat)
        if (report) 
            exportHtml <- cbind(exportHtml,.makeHtmlCells(adjpMat))
        theNames <- c(theNames,paste("FDR_",colnames(adjpMat),sep=""))
    }
    if ("meta_p_value" %in% exportWhat && length(statistics)>1) { 
        # Otherwise, it does not exist
        disp("      binding meta p-values...")
        export <- cbind(export,sumP)
        if (report) 
            exportHtml <- cbind(exportHtml,.makeHtmlCells(sumP))
        theNames <- c(theNames,paste("meta_p-value_",cnt,sep=""))
    }
    if ("adj_meta_p_value" %in% exportWhat && length(statistics)>1) {
        disp("      binding adjusted meta p-values...")
        export <- cbind(export,adjSumP)
        if (report) 
            exportHtml <- cbind(exportHtml,.makeHtmlCells(adjSumP))
        theNames <- c(theNames,paste("meta_FDR_",cnt,sep=""))
    }
    if ("fold_change" %in% exportWhat) {
        if ("normalized" %in% exportValues) {
            tmp <- makeFoldChange(cnt,sampleList,normGeneCounts,logOffset)
            if (!("natural" %in% exportScale || "natural" %in% exportScale)) {
                # In case only "rpgm" output requested, at least log2 fold
                # changes must be exported if requested
                disp("      binding log2 normalized fold changes...")
                export <- cbind(export,log2(tmp))
                if (report) 
                    exportHtml <- cbind(exportHtml,.makeHtmlCells(log2(tmp)))
                theNames <- c(theNames,paste("log2_normalized_fold_change_",
                    colnames(tmp),sep=""))
            }
            else {
                if ("natural" %in% exportScale) {
                disp("      binding natural normalized fold changes...")
                    export <- cbind(export,tmp)
                    if (report) 
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,
                        paste("natural_normalized_fold_change_",colnames(tmp),
                            sep=""))
                }
                if ("log2" %in% exportScale) {
                    disp("      binding log2 normalized fold changes...")
                    export <- cbind(export,log2(tmp))
                    if (report) 
                        exportHtml <- cbind(exportHtml,
                            .makeHtmlCells(log2(tmp)))
                    theNames <- c(theNames,paste("log2_normalized_fold_change_",
                        colnames(tmp),sep=""))
                }
            }
        }
        if ("raw" %in% exportValues) {
            tmp <- makeFoldChange(cnt,sampleList,rawGeneCounts,logOffset)
            if (!("natural" %in% exportScale || "natural" %in% exportScale)) {
                disp("      binding log2 raw fold changes...")
                export <- cbind(export,log2(tmp))
                if (report) 
                    exportHtml <- cbind(exportHtml,.makeHtmlCells(log2(tmp)))
                theNames <- c(theNames,paste("log2_raw_fold_change_",
                    colnames(tmp),sep=""))
            }
            else {
                if ("natural" %in% exportScale) {
                    disp("      binding natural raw fold changes...")
                    export <- cbind(export,tmp)
                    if (report) exportHtml <- cbind(exportHtml,
                        .makeHtmlCells(tmp))
                    theNames <- c(theNames,paste("natural_raw_fold_change_",
                        colnames(tmp),sep=""))
                }
                if ("log2" %in% exportScale) {
                    disp("      binding log2 raw fold changes...")
                    export <- cbind(export,log2(tmp))
                    if (report) 
                        exportHtml <- cbind(exportHtml,
                            .makeHtmlCells(log2(tmp)))
                    theNames <- c(theNames,paste("log2_raw_fold_change_",
                        colnames(tmp),sep=""))
                }
            }
            
        }
    }
    if ("stats" %in% exportWhat) {
        conds <- strsplit(cnt,"_vs_")[[1]]
        for (cond in conds) {
            if ("normalized" %in% exportValues) {
                if ("mean" %in% exportStats) {
                    disp("      binding normalized mean counts...")
                    tmp <- makeStat(sampleList[[cond]],normList,"mean",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report) 
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_mean_counts_",cond,sep=""))
                }
                if ("median" %in% exportStats) {
                    disp("      binding normalized median counts...")
                    tmp <- makeStat(sampleList[[cond]],normList,"median",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report) 
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_median_counts_",cond,sep=""))
                }
                if ("sd" %in% exportStats) {
                    disp("      binding normalized count sds...")
                    tmp <- makeStat(sampleList[[cond]],normList,"sd",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report) 
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_sd_counts_",cond,sep=""))
                }
                if ("mad" %in% exportStats) {
                    disp("      binding normalized count MADs...")
                    tmp <- makeStat(sampleList[[cond]],normList,"mad",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_mad_counts_",cond,sep=""))
                }
                if ("cv" %in% exportStats) {
                    disp("      binding normalized count CVs...")
                    tmp <- makeStat(sampleList[[cond]],normList,"cv",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_cv_counts_",cond,sep=""))
                }
                if ("rcv" %in% exportStats) {
                    disp("      binding normalized counts RCVs...")
                    tmp <- makeStat(sampleList[[cond]],normList,"rcv",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_normalized_rcv_counts_",cond,sep=""))
                }
            }
            if ("raw" %in% exportValues) {
                if ("mean" %in% exportStats) {
                    disp("      binding raw mean counts...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"mean",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_mean_counts_",cond,sep=""))
                }
                if ("median" %in% exportStats) {
                    disp("      binding raw median counts...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"median",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_median_counts_",cond,sep=""))
                }
                if ("sd" %in% exportStats) {
                    disp("      binding raw counts sds...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"sd",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_sd_counts_",cond,sep=""))
                }
                if ("mad" %in% exportStats) {
                    disp("      binding raw counts MADs...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"mad",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_mad_counts_",cond,sep=""))
                }
                if ("cv" %in% exportStats) {
                    disp("      binding raw counts CVs...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"cv",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_cv_counts_",cond,sep=""))
                }
                if ("rcv" %in% exportStats) {
                    disp("      binding raw counts RCVs...")
                    tmp <- makeStat(sampleList[[cond]],rawList,"rcv",
                        exportScale)
                    export <- cbind(export,tmp)
                    if (report)
                        exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                    theNames <- c(theNames,paste(colnames(tmp),
                        "_raw_rcv_counts_",cond,sep=""))
                }
            }
        }
    }
    if ("counts" %in% exportWhat) {
        conds <- strsplit(cnt,"_vs_")[[1]]
        for (cond in conds) {
            if ("normalized" %in% exportValues) {                    
                disp("      binding all normalized counts for ",cond,"...")
                tmp <- makeMatrix(sampleList[[cond]],normList,exportScale)
                export <- cbind(export,tmp)
                if (report)
                    exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                part1 <- rep(paste(exportScale,"_normalized_counts_",sep=""),
                    each=length(sampleList[[cond]]))
                part2 <- paste(part1,colnames(tmp),sep="")
                theNames <- c(theNames,part2)
            }
            if ("raw" %in% exportValues) {
                disp("      binding all raw counts for ",cond,"...")
                tmp <- makeMatrix(sampleList[[cond]],rawList,exportScale)
                export <- cbind(export,tmp)
                if (report)
                    exportHtml <- cbind(exportHtml,.makeHtmlCells(tmp))
                part1 <- rep(paste(exportScale,"_raw_counts_",sep=""),
                    each=length(sampleList[[cond]]))
                part2 <- paste(part1,colnames(tmp),sep="")
                theNames <- c(theNames,part2)
            }
        }
    }
    if ("flags" %in% exportWhat && !is.null(flags)) {
        disp("      binding filtering flags...")
        export <- cbind(export,as.data.frame(flags))
        if (report)
            exportHtml <- cbind(exportHtml,.makeHtmlCells(flags))
        theNames <- c(theNames,colnames(flags))
    }
    names(export) <- theNames

    if (!report)
        exportHtml <- NULL

    return (list(textTable=export,htmlTable=exportHtml,headers=theNames))
}

createSignalTracks <- function(targets,type=c("bigwig","wig"),org,urlBase=NULL,
    stranded=FALSE,normTo=1e+9,exportPath=".") {
    if (!requireNamespace(rtracklayer))
        stopwrap("Bioconductor package rtracklayer is required to build ",
            "tracks!")
    
    if (!is.list(targets)) {
        if (file.exists(targets))
            targets <- readTargets(targets)
        else
            stopwrap("targets must be either the result of readTargets ",
                "function or a valid targets file!")
    }
    
    # Predefined condition colors
    # Track colors - to be recycled in case of many classes for each dataset
    posBaseColours <- c("#B40000","#00B400","#0000B4","#B45200","#9B59B6",
        "#21BCBF","#BC4800","#135C34","#838F00","#4900B5")
    negBaseColours <- c("#FF7575","#79FF79","#8484FF","#FFB77C","#EBB7FF",
        "#63E5E7","#FFB88B","#5BFFA5","#F4FF75","#C69FFF")
    
    # Group colors
    posCol <- rep(posBaseColours,length.out=length(targets$samples))
    posCol <- rep(posCol,lengths(targets$samples))
    names(posCol) <- unlist(targets$samples,use.names=FALSE)
    if (stranded) {
        negCol <- rep(negBaseColours,length.out=length(targets$samples))
        negCol <- rep(negCol,lengths(targets$samples))
        names(negCol) <- names(posCol)
    }
    else
        negCol <- NULL
    
    # Read in BAMs and make Rles
    message("Creating tracks...")
    bams <- unlist(targets$files,use.names=FALSE)
    
    # Get seqinfo form BAM
    preSf <- metaseqR2:::.chromInfoFromBAM(bams[1])
    vchrs <- metaseqR2:::getValidChrs(org)
    preSf <- preSf[vchrs,,drop=FALSE]
    sf <- metaseqR2:::.chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)
    
    message("Reading BAM files to Rle...")
    if (stranded) {
        message("  stranded mode")
    }
    else {
        message("  unstranded mode")
        bw <- .createUnstrandedSignalTracks(bams,org,sf,normTo,vchrs,urlBase,
            exportPath,posCol)
    }
}

.createUnstrandedSignalTracks <- function(bams,org,sf,normTo,vchrs,urlBase,
    exportPath,posCol) {
    
    rles <- metaseqR2:::cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b))))
        cov <- coverage(reads)
        cov <- cov[v]
        seqinfo(cov) <- s
        return(cov)
    },vchrs,sf)
    names(rles) <- names(posCol)
    
    # Calculate normalization factors
    rawSums <- sapply(lapply(rles,sum),sum)
    normFacs <- normTo/rawSums
    names(normFacs) <- names(rles)
    
    # Normalize (should be quick)
    nrles <- metaseqR2:::cmclapply(names(rles),function(n,R,N) {
        for (chr in names(R[[n]]))
            R[[n]][[chr]] <- R[[n]][[chr]]*N[n]
        return(R[[n]])
    },rles,normFacs)
    names(nrles) <- names(rles)
    
    #for (n in names(nrles)) {
    #   export.bw(nrles[[n]])
    #}
}




