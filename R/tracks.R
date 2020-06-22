createSignalTracks <- function(targets,org,urlBase=NULL,stranded=FALSE,
    normTo=1e+9,exportPath=".",hubInfo=list(name="MyHub",shortLabel="My hub",
    longLabel="My hub",email="someone@example.com"),overwrite=FALSE,rc=NULL) {
    if (!requireNamespace("rtracklayer"))
        stopwrap("Bioconductor package rtracklayer is required to build ",
            "tracks!")

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
        return(.createStrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,hubInfo,overwrite,rc=rc))
    }
    else {
        message("  unstranded mode")
        return(.createUnstrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,overwrite,rc=rc))
    }
}

.createStrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,hubInfo,overwrite,rc=NULL) {
    # First check if bigWig files already exist
    trackDbPath <- file.path(exportPath,getUcscOrganism(org))
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

    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    vchrs <- getValidChrs(org)
    preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
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
        gr <- as(slice(cov,1),"GRanges")
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
        gr <- as(slice(cov,1),"GRanges")
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
    
    # Write the genomes file
    gf <- .makeTrackhubGenomesFile(org)
    writeLines(paste(names(gf)," ",gf,sep=""),
        file.path(exportPath,"genomes.txt"))
    
    # Write the hub file
    hub <- .makeTrackhubEntrypoint(hubInfo)
    writeLines(paste(names(hub)," ",hub,sep=""),
        file.path(exportPath,"hub.txt"))
    
    # Write the trackDb file...
    pretdb <- .makeTrackhubMultiwig(names(posBwFiles),posBwFiles,negBwFiles,
        posCol,negCol,urlBase,org)
    # Clear previous file if exists, otherwise append will cause problems
    if (file.exists(file.path(trackDbPath,"trackDb.txt")))
        unlink(file.path(trackDbPath,"trackDb.txt"))
    for (p in pretdb) {
        tracks <- p$tracks
        p$tracks <- NULL
        meta <- paste(paste(names(p)," ",p,sep=""),collapse="\n")
        post <- paste(names(tracks$positive)," ",tracks$positive,sep="")
        post <- paste(paste0("    ",post),sep="",collapse="\n")
        negt <- paste(names(tracks$negative)," ",tracks$negative,sep="")
        negt <- paste(paste0("    ",negt),sep="",collapse="\n")
        aTrack <- c(meta,"\n\n",post,"\n\n",negt,"\n\n")
        cat(aTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",append=TRUE)
    }
    
    # Generate the hub link
    return(paste0(urlBase,"/hub.txt"))
}

.createUnstrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,overwrite,rc=NULL) {
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
    
    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    vchrs <- getValidChrs(org)
    preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
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
        gr <- as(slice(cov,1),"GRanges")
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
    
    message("Creating track lines...")
    opts <- .makeUnstrandedTrackOpts(posCol)
    trackLines <- character(length(bwFiles))
    for (i in seq_along(bwFiles)) {
        opts[[i]]$name <- sub("^([^.]*).*","\\1",basename(bwFiles[i]))
        opts[[i]]$name <- gsub(" ","_",opts[[i]]$name)
        opts[[i]]$bigDataUrl <- paste0(urlBase,"/",basename(bwFiles[i]))
        popts <- paste(names(opts[[i]]),"=",opts[[i]],sep="")
        trackLines[i] <- paste(popts,collapse=" ")
        trackLines[i] <- paste("track",trackLines[i])
    }
    
    tracksFile <- file.path(exportPath,"tracks.txt")
    writeLines(trackLines,tracksFile)
    tracksLink <- paste0(urlBase,"/tracks.txt")
    return(tracksLink)
}

.makeUnstrandedTrackOpts <- function(cols) {
    opts <- vector("list",length(cols))
    for (i in seq_len(length(cols))) {
        opts[[i]]$type <- "bigWig"
        opts[[i]]$color <- paste(t(col2rgb(cols[i])),collapse=",")
        opts[[i]]$visibility <- "full"
        opts[[i]]$maxHeightPixels <- "128:64:16"
    }
    return(opts)
}

.makeTrackhubMultiwig <- function(tnames,pfiles,nfiles,pcols,ncols,
    url,org) {
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
        opts[[i]]$tracks <- list(
            positive=list(
                track=paste0(tnames[i],"_plus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(pcols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",getUcscOrganism(org),"/",
                    paste0(tnames[i],"_plus.bigWig"))
            ),
            negative=list(
                track=paste0(tnames[i],"_minus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(ncols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",getUcscOrganism(org),"/",
                    paste0(tnames[i],"_minus.bigWig"))
            )
        )
    }
    return(opts)
}

.makeTrackhubGenomesFile <- function(org) {
    return(list(
        genome=getUcscOrganism(org),
        trackDb=paste0(org,"/trackDb.txt")
    ))
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
