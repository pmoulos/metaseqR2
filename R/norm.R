normalizeEdaseq <- function(geneCounts,sampleList,normArgs=NULL,
    geneData=NULL,output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","edaseq")
    if (!is.matrix(geneCounts))
        geneCounts <- as.matrix(geneCounts)
    if (!is.null(geneData) && is(geneData,"GenomicRanges")) {
        gl <- NULL
        if (!is.null(attr(geneData,"geneLength")))
            gl <- attr(geneData,"geneLength")
        geneData <- as.data.frame(geneData)
        geneData <- geneData[,c(1,2,3,6,7,5,8,9)]
        if (!is.null(gl))
            attr(geneData,"geneLength") <- gl
    }
    if (!is.null(geneData) && is.null(attr(geneData,"geneLength")))
        attr(geneData,"geneLength") <- rep(1,nrow(geneCounts))
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    if (is.null(geneData)) {
        seqGenes <- newSeqExpressionSet(
            geneCounts,
            phenoData=AnnotatedDataFrame(
                data.frame(
                    conditions=classes,
                    row.names=colnames(geneCounts)
                )
            ),
            featureData=AnnotatedDataFrame(
                data.frame(
                    length=rep(1,nrow(geneCounts)),
                    row.names=rownames(geneCounts)
                )
            )
        )
        seqGenes <- betweenLaneNormalization(withinLaneNormalization(seqGenes,
            "length",which=normArgs$within.which),
            which=normArgs$between.which)
    }
    else {
        seqGenes <- newSeqExpressionSet(
            geneCounts,
            phenoData=AnnotatedDataFrame(
                data.frame(
                    conditions=classes,
                    row.names=colnames(geneCounts)
                )
            ),
            featureData=AnnotatedDataFrame(
                data.frame(
                    gc=geneData$gc_content,
                    length=attr(geneData,"geneLength"),
                    row.names=if (is.data.frame(geneData)) rownames(geneData)
                        else names(geneData)
                )
            )
        )
        seqGenes <- betweenLaneNormalization(withinLaneNormalization(seqGenes,
            "gc",which=normArgs$within.which),which=normArgs$between.which)
    }
    if (output=="matrix")
        return(normCounts(seqGenes)) # Class: matrix
    else if (output=="native")
        return(seqGenes) # Class: SeqExpressionSet
}

normalizeDeseq <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","deseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    cds <- newCountDataSet(geneCounts,classes)
    #cds <- DESeq::estimateSizeFactors(cds,locfunc=normArgs$locfunc)
    cds <- estimateSizeFactors(cds,locfunc=normArgs$locfunc)
    if (output=="native")
        return(cds) # Class: .CountDataSet
    else if (output=="matrix")
        #return(round(DESeq::counts(cds,normalized=TRUE))) # Class: matrix
        return(round(counts(cds,normalized=TRUE))) # Class: matrix
}

normalizeDeseq2 <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","deseq2")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    conditions=as.factor(classes)                 
    colData=DataFrame(conditions)
    design= as.formula(c("~", names(colData[1]))) 
    dds <- DESeqDataSetFromMatrix(geneCounts,colData,design=design,
        tidy=normArgs$tidy)
    dds <- DESeq2::estimateSizeFactors(dds,type=normArgs$type,
        locfunc=normArgs$locfunc)
    if (output=="native")
        return(dds) # Class: DESeqDataSet
    else if (output=="matrix")
        return(round(DESeq2::counts(dds,normalized=TRUE))) # Class: matrix
}

normalizeEdger <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","edger")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    dge <- DGEList(counts=geneCounts,group=classes)
    dge <- calcNormFactors(dge,method=normArgs$method,
        refColumn=normArgs$refColumn,logratioTrim=normArgs$logratioTrim,
        sumTrim=normArgs$sumTrim,doWeighting=normArgs$doWeighting,
        Acutoff=normArgs$Acutoff,p=normArgs$p)
    if (output=="native")
        return(dge) # Class: DGEList
    else if (output=="matrix") {
        scl <- dge$samples$lib.size * dge$samples$norm.factors
        return(round(t(t(dge$counts)/scl)*mean(scl)))
        #return(round(dge$pseudo.counts)) # Class: matrix
    }
}

normalizeNoiseq <- function(geneCounts,sampleList,normArgs=NULL,
    geneData=NULL,logOffset=1,output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","noiseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    if (!is.null(geneData)) {
        if (is(geneData,"GenomicRanges")) {
            gl <- NULL
            if (!is.null(attr(geneData,"geneLength")))
                gl <- attr(geneData,"geneLength")
            else
                gl <- width(geneData)
            geneData <- as.data.frame(geneData)
            geneData <- geneData[,c(1,2,3,6,7,5,8,9)]
            if (!is.null(gl))
                attr(geneData,"geneLength") <- gl
        }
        else {
            if (is.null(attr(geneData,"geneLength"))) {
                gl <- geneData$end - geneData$start + 1
                attr(geneData,"geneLength") <- gl
            }
        }
    }
    if (is.null(geneData)) {
        #nsObj <- NOISeq::readData(
        nsObj <- readData(
            data=geneCounts,
            factors=data.frame(class=classes)
        )
    }
    else {
        gcContent <- geneData$gc_content
        geneLength <- attr(geneData,"geneLength")
        biotype <- as.character(geneData$biotype)
        names(gcContent) <- names(biotype) <- names(geneLength) <- 
            if (is.data.frame(geneData)) rownames(geneData) else names(geneData)
        #nsObj <- NOISeq::readData(
        nsObj <- readData(
            data=geneCounts,
            length=geneLength,
            gc=gcContent,
            chromosome=geneData[,c(1,2,3)],
            factors=data.frame(class=classes),
            biotype=biotype
        )
    }
    normArgs$k <- logOffset # Set the zero fixing constant
    switch(normArgs$method,
        rpkm = {
            #M <- NOISeq::rpkm(assayData(nsObj)$exprs,long=normArgs$long,
            #    k=normArgs$k,lc=normArgs$lc)
            M <- noirpkm(assayData(nsObj)$exprs,long=normArgs$long)
        },
        uqua = {
            M <- noiuqua(assayData(nsObj)$exprs,long=normArgs$long,
                k=normArgs$k,lc=normArgs$lc)
        },
        tmm = {
            M <- noitmm(assayData(nsObj)$exprs,long=normArgs$long,
                k=normArgs$k,lc=normArgs$lc,refColumn=normArgs$refColumn,
                logratioTrim=normArgs$logratioTrim,sumTrim=normArgs$sumTrim,
                doWeighting=normArgs$doWeighting,Acutoff=normArgs$Acutoff)
        }
    )
    if (output=="native") {
        if (is.null(geneData))
            #return(NOISeq::readData(
            return(readData(
                data=M,
                factors=data.frame(class=classes)
        )) # Class: CD
        else    
            #return(NOISeq::readData(
            return(readData(
                data=M,
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )) # Class: CD
    }
    else if (output=="matrix")
        return(as.matrix(round(M))) # Class: matrix
}

normalizeNbpseq <- function(geneCounts,sampleList,normArgs=NULL,
    libsizeList=NULL,output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","nbpseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    if (is.null(libsizeList)) {
        libsizeList <- vector("list",length(classes))
        names(libsizeList) <- unlist(sampleList,use.names=FALSE)
        for (n in names(libsizeList))
            libsizeList[[n]] <- sum(geneCounts[,n])
    }
    libSizes <- unlist(libsizeList)
    norm.factors <- estimate.norm.factors(geneCounts,lib.sizes=libSizes,
        method=normArgs$method)
    #if (normArgs$main.method=="nbpseq")
    #    nbData <- prepare.nbData(geneCounts,libSizes=libSizes,
    #    norm.factors=norm.factors)
    #else if (normArgs$main.method=="nbsmyth")
        nbData <- prepare.nbp(geneCounts,classes,lib.sizes=libSizes,
            norm.factors=norm.factors,thinning=normArgs$thinning)
    if (output=="native")
        return(nbData) # Class: list or nbp
    else if (output=="matrix") {
        #if (normArgs$main.method=="nbpseq") {
        #    normCounts <- matrix(0,nrow(geneCounts),ncol(geneCounts))
        #    for (i in seq_len(ncol(geneCounts)))
        #        normCounts[,i] <- norm.factors[i]*geneCounts[,i]
        #}
        #else if (normArgs$main.method=="nbsmyth") 
            normCounts <- nbData$pseudo.counts
        return(as.matrix(round(normCounts))) # Class: matrix
    }
}

normalizeAbsseq <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","absseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    #classes <- asClassVector(sampleList)
    #groups= as.factor(classes)
    abs <- ABSDataSet(geneCounts,normMethod= normArgs$normMethod)  
    abs <- normalFactors(abs)
    excounts(abs) <- ABSSeq::counts(abs,norm=TRUE)
    if (output=="native")
        return(abs) # Class: ABSDataSet
    else if (output=="matrix")
        return(as.matrix(round(excounts(abs)))) # Class: matrix
}

normalizeDss <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","dss")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    
    # make the design
    classes <- asClassVector(sampleList)
    design <- as.data.frame(classes)
    colnames(design) <- "designs"
    
    # newSeqCountSet takes only matrices as gene.counts
    if (is(geneCounts,"data.frame")) {
        geneCountsTmp <- geneCounts
        samCols <- match(unlist(sampleList),colnames(geneCountsTmp))
        samCols <- samCols[which(!is.na(samCols))]
        geneCountsTmp <- as.matrix(geneCounts[,samCols,drop=FALSE])
    }
    if (is(geneCounts,"matrix"))
        geneCountsTmp <- geneCounts
    
    seqD <- newSeqCountSet(geneCountsTmp,design) # create the class
    # estimate normalization factors
    seqD <- estNormFactors(seqD,method=normArgs$method)
    
    if (output=="native")
        return(seqD) # Class: SeqCountSet
    else if (output=="matrix") {
        # I cannot extract normalized counts from DSS and thus I do this small 
        # dribble in order to output a matrix
        theDesign <- data.frame(condition=classes,row.names=colnames(seqD)) 
        cds <- newCountDataSet(as.matrix(round(assayData(seqD)$exprs)),
            theDesign$condition)
        #DESeq::sizeFactors(cds) <- normalizationFactor(seqD)
        sizeFactors(cds) <- normalizationFactor(seqD)
        #counts <- as.matrix(DESeq::counts(cds,normalized=TRUE))
        counts <- as.matrix(counts(cds,normalized=TRUE))
        return(as.matrix(round(counts))) # Class: matrix
    }
}
