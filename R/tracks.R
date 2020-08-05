createSignalTracks <- function(targets,org,urlBase=NULL,stranded=FALSE,
    normTo=1e+9,exportPath=".",hubInfo=list(name="MyHub",shortLabel="My hub",
    longLabel="My hub",email="someone@example.com"),fasta=NULL,gtf=NULL,
    forceHub=FALSE,overwrite=FALSE,rc=NULL) {
    if (!requireNamespace("rtracklayer"))
        stopwrap("Bioconductor package rtracklayer is required to build ",
            "tracks!")
    
    if (!is.null(fasta) && !requireNamespace("Biostrings"))
        stopwrap("Bioconductor package Biostrings is required to read ",
            "FASTA files!")

    if (!is.list(targets)) {
        if (file.exists(targets))
            targets <- readTargets(targets)
        else
            stopwrap("targets must be either the result of readTargets ",
                "function or a valid targets file!")
    }

    # Read in BAMs and make Rles
    message("Creating tracks...")
    if (stranded) {
        message("  stranded mode")
        .createTrackHub(targets,org,normTo,stranded=TRUE,urlBase,
            exportPath,hubInfo,fasta,gtf,overwrite,rc=rc)
        return(paste0(urlBase,"/hub.txt"))
        #return(.createStrandedSignalTracks(targets,org,normTo,urlBase,
        #    exportPath,hubInfo,overwrite,rc=rc))
    }
    else {
        message("  unstranded mode")
        if (forceHub) {
            .createTrackHub(targets,org,normTo,stranded=FALSE,urlBase,
                exportPath,hubInfo,fasta,gtf,overwrite,rc=rc)
            return(paste0(urlBase,"/hub.txt"))
        }
        else {
            trackLines <- .createUnstrandedSignalTracks(targets,org,normTo,
                urlBase,exportPath,overwrite,asLines=TRUE,rc=rc)
            tracksFile <- file.path(exportPath,"tracks.txt")
            writeLines(trackLines,tracksFile)
            tracksLink <- paste0(urlBase,"/tracks.txt")
            return(tracksLink)
        }
        
    }
}

.createTrackHub <- function(targets,org,normTo,stranded=FALSE,urlBase,
    exportPath,hubInfo,fasta=NULL,gtf=NULL,overwrite=FALSE,rc=NULL) {
    supOrg <- org %in% getSupportedOrganisms()
    if (supOrg)
        org <- getUcscOrganism(org)
    trackDbPath <- file.path(exportPath,org)
    
    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    
    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    if (supOrg) {
        vchrs <- getValidChrs(org)
        preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
    }
    sf <- .chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)
    
    # Start creating the hub
    groups <- rep(names(targets$files),lengths(targets$files))
    if (is.null(groups))
        groups <- rep("Signal",length(bams))
    
    # First the dir structure
    if (!dir.exists(trackDbPath))
        dir.create(trackDbPath,recursive=TRUE)
    
    # Has a custom FASTA file been provided?
    has2bit <- FALSE
    if (!supOrg && !is.null(fasta) && is.character(fasta) 
        && file.exists(fasta)) {
        has2bit <- TRUE
        message("Reading file ",fasta)
        ss <- readDNAStringSet(fasta)
        names(ss) <- 
            vapply(strsplit(names(ss)," "),function(x) x[1],character(1))
        message("Writing file ",file.path(trackDbPath,paste0(org,".2bit")))
        export.2bit(ss,file.path(trackDbPath,paste0(org,".2bit")))
    }
    
    # Has a custom GTF file been provided? More complex...
    hasGtf <- FALSE
    if (!supOrg && !is.null(gtf) && is.character(gtf)
        && file.exists(gtf)) {
        hasGtf <- TRUE
        message("Reading file ",gtf," and creating bigBed")
        bb <- .externalGtfToBigBed(gtf,preSf)
        file.copy(bb,file.path(trackDbPath,paste0(org,".bigBed")),
            overwrite=TRUE)
    }
    
    # Write the genomes file
    gf <- .makeTrackhubGenomesFile(org,bams[1],has2bit)
    writeLines(paste(names(gf)," ",gf,sep=""),
        file.path(exportPath,"genomes.txt"))
        
    # Write the groups file
    if (file.exists(file.path(trackDbPath,"groups.txt")))
        unlink(file.path(trackDbPath,"groups.txt"))
    grp <- .makeTrackhubGroupsFile(targets,hasGtf)
    for (g in grp) {
        grm <- paste(paste(names(g)," ",g,sep=""),collapse="\n")
        aGroup <- c(grm,"\n\n")
        cat(aGroup,file=file.path(trackDbPath,"groups.txt"),append=TRUE)
    }
    
    # Create the bigWigs and write the trackDb file
    # Clear previous file if exists, otherwise append will cause problems
    if (file.exists(file.path(trackDbPath,"trackDb.txt")))
        unlink(file.path(trackDbPath,"trackDb.txt"))
    if (stranded) {
        # BigWig attributes
        pretdb <- .createStrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,overwrite=overwrite,rc=rc)
        for (p in pretdb) {
            tracks <- p$tracks
            p$tracks <- NULL
            meta <- paste(paste(names(p)," ",p,sep=""),collapse="\n")
            post <- paste(names(tracks$positive)," ",tracks$positive,sep="")
            post <- paste(paste0("    ",post),sep="",collapse="\n")
            negt <- paste(names(tracks$negative)," ",tracks$negative,sep="")
            negt <- paste(paste0("    ",negt),sep="",collapse="\n")
            aTrack <- c(meta,"\n\n",post,"\n\n",negt,"\n\n")
            cat(aTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",
                append=TRUE)
        }
    }
    else {
        pretdb <- .createUnstrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,overwrite=overwrite,asLines=FALSE,rc=rc)
        for (p in pretdb) {
            meta <- paste(paste(names(p)," ",p,sep=""),collapse="\n")
            aTrack <- c(meta,"\n\n")
            cat(aTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",
                append=TRUE)
        }
    }
    
    # Is there an annotation bigBed track?
    if (hasGtf) {
        annTrack <- list(
            track=paste0(org,"_genes"),
            type="bigBed",
            shortLabel=paste0(org," genes"),
            longLabel=paste0(org," gene models"),
            boxedCfg="on",
            color="184,0,212",
            exonArrows="on",
            visibility="dense",
            group="annotation",
            bigDataUrl=paste0(urlBase,"/",org,"/",paste0(org,".bigBed"))
        )
        annTrack <- paste(paste(names(annTrack)," ",annTrack,sep=""),
            collapse="\n")
        annTrack <- c(annTrack,"\n\n")
        cat(annTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",
            append=TRUE)
    }
    
    # Finally, write and return the hub file
    hub <- .makeTrackhubEntrypoint(hubInfo)
    
    writeLines(paste(names(hub)," ",hub,sep=""),
        file.path(exportPath,"hub.txt"))
}

.createStrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,hubInfo,fasta=NULL,gtf=NULL,overwrite=FALSE,rc=NULL) {
    # First check if bigWig files already exist
    supOrg <- org %in% getSupportedOrganisms()
    if (supOrg)
        org <- getUcscOrganism(org)
    
    #trackDbPath <- file.path(exportPath,getUcscOrganism(org))
    trackDbPath <- file.path(exportPath,org)
    posCheck <- file.path(trackDbPath,
        paste0(unlist(targets$samples,use.names=FALSE),"_plus.bigWig"))
    negCheck <- file.path(trackDbPath,
        paste0(unlist(targets$samples,use.names=FALSE),"_minus.bigWig"))
    posEx <- file.exists(posCheck)
    negEx <- file.exists(negCheck)
    if (all(posEx) && all(negEx) && !overwrite) {
        message("Requested tracks have already been created! ",
            "Use overwrite=TRUE to re-create them.")
        return("")
    }

    # Positive and negative color options
    posBaseColours <- .getPosBaseColors()
    negBaseColours <- .getNegBaseColors()
    posCol <- rep(posBaseColours,length.out=length(targets$samples))
    posCol <- rep(posCol,lengths(targets$samples))
    negCol <- rep(negBaseColours,length.out=length(targets$samples))
    negCol <- rep(negCol,lengths(targets$samples))
    names(posCol) <- names(negCol) <- unlist(targets$samples,use.names=FALSE)

    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    groups <- rep(names(targets$files),lengths(targets$files))
    if (is.null(groups))
        groups <- rep("Signal",length(bams))

    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    if (supOrg) {
        vchrs <- getValidChrs(org)
        preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
    }
    else
        vchrs <- rownames(preSf)
    sf <- .chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)

    # Get coverage and assign seqinfo for bigwig
    message("Reading positive strand reads from BAM files to Rle...")
    pbg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b,
            param=ScanBamParam(scanBamFlag(isMinusStrand=FALSE))))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(IRanges::slice(cov,1),"GRanges")
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(pbg) <- names(posCol)

    message("Reading negative strand reads from BAM files to Rle...")
    nbg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b,
            param=ScanBamParam(scanBamFlag(isMinusStrand=TRUE))))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(IRanges::slice(cov,1),"GRanges")
        # Inverse the coverage
        gr$score <- -gr$score
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(nbg) <- names(negCol)

    # Calculate normalization factors
    rawPosSums <- vapply(pbg,function(x) sum(x$score),numeric(1))
    rawNegSums <- vapply(nbg,function(x) sum(x$score),numeric(1))

    rat <- rawPosSums/-rawNegSums
    posNormTo <- rat*0.5*normTo
    negNormTo <- normTo - posNormTo

    posNormFacs <- 0.5*posNormTo/rawPosSums
    negNormFacs <- -0.5*negNormTo/rawNegSums

    names(posNormFacs) <- names(pbg)
    names(negNormFacs) <- names(nbg)

    # Normalize (should be quick)
    message("Normalizing positive strand...")
    npbg <- cmclapply(names(pbg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },pbg,posNormFacs)
    names(npbg) <- names(pbg)
    
    message("Normalizing negative strand...")
    nnbg <- cmclapply(names(nbg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },nbg,negNormFacs)
    names(nnbg) <- names(nbg)
    
    # Start creating the hub
    # First the dir structure
    if (!dir.exists(trackDbPath))
        dir.create(trackDbPath,recursive=TRUE)
        
    # Put the bigWig files there
    message("Exporting bigWig files...")
    posBwFiles <- file.path(trackDbPath,paste0(names(npbg),"_plus.bigWig"))
    names(posBwFiles) <- names(npbg)
    for (n in names(npbg)) {
        message("  exporting + strand for sample ",n," to ",posBwFiles[n])
        export.bw(npbg[[n]],posBwFiles[n])
    }
    negBwFiles <- file.path(trackDbPath,paste0(names(nnbg),"_minus.bigWig"))
    names(negBwFiles) <- names(nnbg)
    for (n in names(nnbg)) {
        message("  exporting - strand for sample ",n," to ",negBwFiles[n])
        export.bw(nnbg[[n]],negBwFiles[n])
    }
    
    ## Write the genomes file
    #gf <- .makeTrackhubGenomesFile(org)
    #writeLines(paste(names(gf)," ",gf,sep=""),
    #    file.path(exportPath,"genomes.txt"))
    
    ## Write the hub file
    #hub <- .makeTrackhubEntrypoint(hubInfo)
    #writeLines(paste(names(hub)," ",hub,sep=""),
    #    file.path(exportPath,"hub.txt"))
    
    # Write the trackDb file...
    pretdb <- .makeTrackhubMultiwig(names(posBwFiles),posBwFiles,negBwFiles,
        posCol,negCol,urlBase,org,groups)
    
    return(pretdb)
    
    ## Clear previous file if exists, otherwise append will cause problems
    #if (file.exists(file.path(trackDbPath,"trackDb.txt")))
    #    unlink(file.path(trackDbPath,"trackDb.txt"))
    #for (p in pretdb) {
    #    tracks <- p$tracks
    #    p$tracks <- NULL
    #    meta <- paste(paste(names(p)," ",p,sep=""),collapse="\n")
    #    post <- paste(names(tracks$positive)," ",tracks$positive,sep="")
    #    post <- paste(paste0("    ",post),sep="",collapse="\n")
    #    negt <- paste(names(tracks$negative)," ",tracks$negative,sep="")
    #    negt <- paste(paste0("    ",negt),sep="",collapse="\n")
    #    aTrack <- c(meta,"\n\n",post,"\n\n",negt,"\n\n")
    #    cat(aTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",
    #        append=TRUE)
    #}
    
    ## Generate the hub link
    #return(paste0(urlBase,"/hub.txt"))
}

.createUnstrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,overwrite,asLines=TRUE,rc=NULL) {
    supOrg <- org %in% getSupportedOrganisms()
    if (supOrg)
        org <- getUcscOrganism(org)
    if (!asLines) # Belong to a trackhub, org is attached to the exportPath
        exportPath <- file.path(exportPath,org)
    
    # First check if bigWig files already exist
    bCheck <- file.path(exportPath,paste0(unlist(targets$samples,
        use.names=FALSE),".bigWig"))
    bEx <- file.exists(bCheck)
    if (all(bEx) && !overwrite) {
        message("Requested tracks have already been created! ",
            "Use overwrite=TRUE to re-create them.")
        return("")
    }
    
    # Create color scheme
    posBaseColours <- .getPosBaseColors()
    posCol <- rep(posBaseColours,length.out=length(targets$samples))
    posCol <- rep(posCol,lengths(targets$samples))
    names(posCol) <- unlist(targets$samples,use.names=FALSE)
    
    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    groups <- rep(names(targets$files),lengths(targets$files))
    if (is.null(groups))
        groups <- rep("Signal",length(bams))
    
    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    if (supOrg) {
        vchrs <- getValidChrs(org)
        preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
    }
    else
        vchrs <- rownames(preSf)
    sf <- .chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)
    
    # Get standed coverage and assign seqinfo for bigwig
    message("Reading BAM files to Rle...")
    bg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(IRanges::slice(cov,1),"GRanges")
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(bg) <- names(posCol)
    
    # Calculate normalization factors
    rawSums <- vapply(bg,function(x) sum(x$score),numeric(1))
    normFacs <- normTo/rawSums
    names(normFacs) <- names(bg)
    
    # Normalize (should be quick)
    message("Normalizing...")
    nbg <- cmclapply(names(bg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },bg,normFacs)
    names(nbg) <- names(bg)
    
    message("Exporting bigWig files...")
    bwFiles <- file.path(exportPath,paste0(names(nbg),".bigWig"))
    names(bwFiles) <- names(nbg)
    for (n in names(nbg)) {
        message("  exporting for sample ",n," to ",bwFiles[n])
        export.bw(nbg[[n]],bwFiles[n])
    }
    
    if (asLines) { # For simple usage
        message("Creating track lines...")
        opts <- vector("list",length(posCol))
        trackLines <- character(length(bwFiles))
        for (i in seq_along(bwFiles)) {
            opts[[i]]$type <- "bigWig"
            opts[[i]]$name <- sub("^([^.]*).*","\\1",basename(bwFiles[i]))
            opts[[i]]$name <- gsub(" ","_",opts[[i]]$name)
            opts[[i]]$description <- paste0("\"",opts[[i]]$name," signal\"")
            opts[[i]]$color <- paste(t(col2rgb(posCol[i])),collapse=",")
            opts[[i]]$visibility <- "full"
            opts[[i]]$maxHeightPixels <- "128:64:16"
            opts[[i]]$bigDataUrl <- paste0(urlBase,"/",basename(bwFiles[i]))
            popts <- paste(names(opts[[i]]),"=",opts[[i]],sep="")
            trackLines[i] <- paste(popts,collapse=" ")
            trackLines[i] <- paste("track",trackLines[i])
        }
        return(trackLines)
    }
    else { # For trackhub usage
        opts <- .makeTrackhubSinglewig(names(bwFiles),bwFiles,posCol,urlBase,
            org,groups)
        return(opts)
    }
    
    #tracksFile <- file.path(exportPath,"tracks.txt")
    #writeLines(trackLines,tracksFile)
    #tracksLink <- paste0(urlBase,"/tracks.txt")
    
    #return(tracksLink)
}

.makeTrackhubSinglewig <- function(tnames,files,cols,url,org,groups) {
    if (org %in% getSupportedOrganisms())
        org <- getUcscOrganism(org)
    opts <- vector("list",length(cols))
    for (i in seq_along(files)) {
        opts[[i]]$track <- tolower(sub("^([^.]*).*","\\1",basename(files[i])))
        opts[[i]]$type <- "bigWig"
        opts[[i]]$shortLabel <- sub("^([^.]*).*","\\1",basename(files[i]))
        opts[[i]]$shortLabel <- gsub(" ","_",opts[[i]]$shortLabel)
        opts[[i]]$longLabel <- paste0(opts[[i]]$shortLabel," signal")
        opts[[i]]$color <- paste(t(col2rgb(cols[i])),collapse=",")
        opts[[i]]$visibility <- "full"
        opts[[i]]$maxHeightPixels <- "128:64:16"
        opts[[i]]$autoScale <- "on"
        opts[[i]]$group <- tolower(groups[i])
        opts[[i]]$bigDataUrl <- paste0(url,"/",org,"/",basename(files[i]))
    }
    return(opts)
}

.makeTrackhubMultiwig <- function(tnames,pfiles,nfiles,pcols,ncols,
    url,org,groups) {
    if (org %in% getSupportedOrganisms())
        org <- getUcscOrganism(org)
    opts <- vector("list",length(pcols))
    for (i in seq_along(pfiles)) {
        opts[[i]]$track <- paste0(tnames[i],"_stranded")
        opts[[i]]$type <- "bigWig"
        opts[[i]]$container <- "multiWig"
        opts[[i]]$aggregate <- "transparentOverlay"
        opts[[i]]$showSubtrackColorOnUi <- "on"
        opts[[i]]$shortLabel <- tnames[i]
        opts[[i]]$longLabel <- paste0(tnames[i]," signal")
        opts[[i]]$boxedCfg <- "on"
        opts[[i]]$autoScale <- "on"
        opts[[i]]$visibility <- "full"
        opts[[i]]$maxHeightPixels <- "128:64:16"
        opts[[i]]$group <- tolower(groups[i])
        opts[[i]]$tracks <- list(
            positive=list(
                track=paste0(tnames[i],"_plus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(pcols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",org,"/",paste0(tnames[i],
                    "_plus.bigWig"))
            ),
            negative=list(
                track=paste0(tnames[i],"_minus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(ncols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",org,"/",paste0(tnames[i],
                    "_minus.bigWig"))
            )
        )
    }
    return(opts)
}

.makeTrackhubGroupsFile <- function(targets,hasGtf=FALSE) {
    groups <- names(targets$files)
    if (is.null(groups))
        groups <- "Signal"
    
    if (hasGtf)
        grp <- vector("list",length(groups)+1)
    else
        grp <- vector("list",length(groups))
        
    for (i in seq_len(length(groups))) {
        grp[[i]]$name <- tolower(gsub(" ","_",groups[i]))
        grp[[i]]$label <- groups[i]
        grp[[i]]$priority <- i
        grp[[i]]$defaultIsClosed <- 0
    }
    
    n <- length(grp)
    if (hasGtf) {
        grp[[n]]$name <- "annotation"
        grp[[n]]$label <- "Annotation"
        grp[[n]]$priority <- n
        grp[[n]]$defaultIsClosed <- 0
    }
    
    return(grp)
}

.makeTrackhubGenomesFile <- function(org,b=NULL,has2bit=FALSE) {
    supOrg <- FALSE
    if (org %in% getSupportedOrganisms()) {
        supOrg <- TRUE
        org <- getUcscOrganism(org)
    }
    gf <- list(
        genome=org,
        groups=paste0(org,"/groups.txt"),
        trackDb=paste0(org,"/trackDb.txt")
    )
    if (!supOrg) {
        ci <- .chromInfoFromBAM(b)
        chr <- rownames(ci)[1]
        len <- as.integer(ci[1,1])
        middle <- round(len/2)
        start <- sample.int(middle,1)
        end <- sample((middle+1):(len-1),1)
        gf$organism <- org
        gf$description <- paste0("Hub for ",org)
        gf$defaultPos <- paste0(chr,":",start,"-",end)
    }
    if (has2bit)
        gf$twoBitPath <- paste0(org,"/",org,".2bit")
    return(gf)
}

.makeTrackhubEntrypoint <- function(h) {
    return(list(
        hub=h$name,
        shortLabel=h$shortLabel,
        longLabel=h$longLabel,
        genomesFile="genomes.txt",
        email=h$email
    ))
}

.getPosBaseColors <- function() {
    return(c("#B40000","#00B400","#0000B4","#B45200","#9B59B6","#21BCBF",
        "#BC4800","#135C34","#838F00","#4900B5"))
}

.getNegBaseColors <- function() {
    return(c("#FF7575","#79FF79","#8484FF","#FFB77C","#EBB7FF","#63E5E7",
        "#FFB88B","#5BFFA5","#F4FF75","#C69FFF"))
}

.externalGtfToBigBed <- function(gtf,ci) {
    # 1. Get the required tools
    # 1.1. gtfToGenePred
    gtfToGenePred <- file.path(tempdir(),"gtfToGenePred")
    if (!file.exists(file.path(tempdir(),"gtfToGenePred"))) {
        message("  Retrieving gtfToGenePred tool")
        download.file(
        "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred",
            gtfToGenePred,quiet=TRUE
        )
        system(paste("chmod 775",gtfToGenePred))
    }
    
    # 1.2. genePredToBigGenePred
    genePredToBigGenePred <- file.path(tempdir(),"genePredToBigGenePred")
    if (!file.exists(file.path(tempdir(),"genePredToBigGenePred"))) {
        message("  Retrieving genePredToBigGenePred tool")
        download.file(paste0("http://hgdownload.soe.ucsc.edu/admin/exe/",
            "linux.x86_64/genePredToBigGenePred"),genePredToBigGenePred,
            quiet=TRUE
        )
        system(paste("chmod 775",genePredToBigGenePred))
    }
    
    # 1.3. bedToBigBed
    bedToBigBed <- file.path(tempdir(),"bedToBigBed")
    if (!file.exists(file.path(tempdir(),"bedToBigBed"))) {
        message("  Retrieving bedToBigBed tool")
        download.file(
        "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed",
            bedToBigBed,quiet=TRUE
        )
        system(paste("chmod 775",bedToBigBed))
    }
    
    # 1.4. bigGenePred.as
    bigGenePred.as <- file.path(tempdir(),"bigGenePred.as")
    if (!file.exists(file.path(tempdir(),"bigGenePred.as"))) {
        message("  Retrieving bigGenePred.as definition")
        download.file(
            "https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as",
            bigGenePred.as,quiet=TRUE
        )
        system(paste("chmod 775",bigGenePred.as))
    }
    
    # 2. Run gtfToGenePred
    #gtfFile <- file.path(tempdir(),basename(gtf))
    message("  Converting ",basename(gtf)," to genePred")
    tmpGenePred <- file.path(tempdir(),paste(format(Sys.time(),"%Y%m%d%H%M%S"),
        "tgp",sep="."))
    command1 <- paste(gtfToGenePred,"-genePredExt",gtf,tmpGenePred)
    message("Executing: ",command1)
    system(command1)
    
    # 3. Run genePredToBigGenePred
    message("  Converting ",basename(tmpGenePred)," to bigGenePred")
    tmpBigGenePred <- file.path(tempdir(),paste(format(Sys.time(),
        "%Y%m%d%H%M%S"),"tbgp",sep="."))
    command2 <- paste(genePredToBigGenePred,tmpGenePred,tmpBigGenePred)
    message("Executing: ",command2)
    system(command2)
    
    # 4. Write the chromosome info file
    chromInfo <- file.path(tempdir(),"chromInfo.txt")
    write.table(ci,file=chromInfo,sep="\t",col.names=FALSE,quote=FALSE)
    
    # 5. Run sort
    message("  Sorting ",tmpBigGenePred)
    tmpSort <- file.path(tempdir(),paste(format(Sys.time(),
        "%Y%m%d%H%M%S"),"txt",sep="."))
    command3 <- paste("sort -k1,1 -k2n,2",tmpBigGenePred,">",tmpSort)
    message("Executing: ",command3)
    system(command3)
    
    # 6. Run bedToBigBed
    tmpBigBed <- file.path(tempdir(),paste(format(Sys.time(),
        "%Y%m%d%H%M%S"),"bigBed",sep="."))
    message("  Converting sorted ",basename(tmpSort)," to bigBed")
    command4 <- paste("bedToBigBed -type=bed12+8 -tab -as=",bigGenePred.as,
        " ",tmpSort," ",chromInfo," ",tmpBigBed,sep="")
    message("Executing: ",command4)
    system(command4)
    
    return(tmpBigBed)
}
