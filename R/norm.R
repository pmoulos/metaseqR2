#' Normalization based on the EDASeq package
#'
#' This function is a wrapper over EDASeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of EDASeq normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"edaseq")} for an example and how
#' you can modify it.
#' @param geneData an optional annotation data frame (such the ones produced by
#' \code{getAnnotation}) which contains the GC content for each gene and from
#' which the gene lengths can be inferred by chromosome coordinates.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the EDASeq native S4
#' object (SeqExpressionSet). In the latter case it should be handled with suitable
#' EDASeq methods.
#' @return A matrix or a SeqExpressionSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' lengths <- round(1000*runif(nrow(dataMatrix)))
#' starts <- round(1000*runif(nrow(dataMatrix)))
#' ends <- starts + lengths
#' gc <- runif(nrow(dataMatrix))
#' geneData <- data.frame(
#'   chromosome=c(rep("chr1",nrow(dataMatrix)/2),rep("chr2",nrow(dataMatrix)/2)),
#'   start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc
#' )
#' norm.dataMatrix <- normalizeEdaseq(dataMatrix,sampleList,geneData=geneData)
#' diagplotBoxplot(norm.dataMatrix,sampleList)      
#'}
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
		geneData <- geneData[,c(1:3,6,7,5,8,9)]
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
        return(exprs(seqGenes)) # Class: matrix
    else if (output=="native")
        return(seqGenes) # Class: SeqExpressionSet
}

#' Normalization based on the DESeq package
#'
#' This function is a wrapper over DESeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of DESeq normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"deseq")} for an example and how you
#' can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the DESeq native S4
#' object (CountDataSet). In the latter case it should be handled with suitable
#' DESeq methods.
#' @return A matrix or a CountDataSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' norm.dataMatrix <- normalizeDeseq(dataMatrix,sampleList)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
normalizeDeseq <- function(geneCounts,sampleList,normArgs=NULL,
    output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","deseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    cds <- newCountDataSet(geneCounts,classes)
    cds <- DESeq::estimateSizeFactors(cds,locfunc=normArgs$locfunc)
    if (output=="native")
        return(cds) # Class: CountDataSet
    else if (output=="matrix")
        return(round(DESeq::counts(cds,normalized=TRUE))) # Class: matrix
}

#' Normalization based on the DESeq2 package
#'
#' This function is a wrapper over DESeq2 normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of DESeq2 normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"deseq2")} for an example and how you
#' can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the DESeq2 native S4
#' object (DESeqDataSet). In the latter case it should be handled with suitable
#' DESeq2 methods.
#' @return A matrix or a DESeqDataSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' norm.dataMatrix <- normalizeDeseq2(dataMatrix,sampleList)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
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

#' Normalization based on the edgeR package
#'
#' This function is a wrapper over edgeR normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of edgeR normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"edger")} for an example and how
#' you can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the edgeR native S4
#' object (DGEList). In the latter case it should be handled with suitable edgeR
#' methods.
#' @return A matrix or a DGEList with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' norm.dataMatrix <- normalizeEdger(dataMatrix,sampleList)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
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

#' Normalization based on the NOISeq package
#'
#' This function is a wrapper over NOISeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of NOISeq normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"noiseq")} for an example and how
#' you can modify it.
#' @param geneData an optional annotation data frame (such the ones produced by
#' \code{getAnnotation} which contains the GC content for each gene and from which
#' the gene lengths can be inferred by chromosome coordinates.
#' @param logOffset an offset to use to avoid infinity in logarithmic data
#' transformations.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the NOISeq native S4
#' object (SeqExpressionSet). In the latter case it should be handled with suitable
#' NOISeq methods.
#' @return A matrix with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' lengths <- round(1000*runif(nrow(dataMatrix)))
#' starts <- round(1000*runif(nrow(dataMatrix)))
#' ends <- starts + lengths
#' gc=runif(nrow(dataMatrix)),
#' geneData <- data.frame(
#'   chromosome=c(rep("chr1",nrow(dataMatrix)/2),rep("chr2",nrow(dataMatrix)/2)),
#'   start=starts,end=ends,gene_id=rownames(dataMatrix),gc_content=gc
#' )
#' norm.dataMatrix <- normalizeNoiseq(dataMatrix,sampleList,geneData)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
normalizeNoiseq <- function(geneCounts,sampleList,normArgs=NULL,
    geneData=NULL,logOffset=1,output=c("matrix","native")) {
    if (is.null(normArgs))
        normArgs <- getDefaults("normalization","noiseq")
    output <- tolower(output[1])
    checkTextArgs("output",output,c("matrix","native"))
    classes <- asClassVector(sampleList)
    if (!is.null(geneData) && is(geneData,"GenomicRanges")) {
		gl <- NULL
		if (!is.null(attr(geneData,"geneLength")))
			gl <- attr(geneData,"geneLength")
		else
			gl <- width(geneData)
		geneData <- as.data.frame(geneData)
		geneData <- geneData[,c(1:3,6,7,5,8,9)]
		if (!is.null(gl))
			attr(geneData,"geneLength") <- gl
	}
	else {
		if (is.null(attr(geneData,"geneLength"))) {
			gl <- geneData$end - geneData$start + 1
			attr(geneData,"geneLength") <- gl
		}
	}

    if (is.null(geneData)) {
        nsObj <- NOISeq::readData(
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
        nsObj <- NOISeq::readData(
            data=geneCounts,
            length=geneLength,
            gc=gcContent,
            chromosome=geneData[,1:3],
            factors=data.frame(class=classes),
            biotype=biotype
        )
    }
    normArgs$k <- logOffset # Set the zero fixing constant
    switch(normArgs$method,
        rpkm = {
            #M <- NOISeq::rpkm(assayData(nsObj)$exprs,long=normArgs$long,
            #    k=normArgs$k,lc=normArgs$lc)
            M <- rpkm(assayData(nsObj)$exprs,geneLength=normArgs$long)
        },
        uqua = {
            M <- NOISeq::uqua(assayData(nsObj)$exprs,long=normArgs$long,
                k=normArgs$k,lc=normArgs$lc)
        },
        tmm = {
            M <- NOISeq::tmm(assayData(nsObj)$exprs,long=normArgs$long,
                k=normArgs$k,lc=normArgs$lc,refColumn=normArgs$refColumn,
                logratioTrim=normArgs$logratioTrim,sumTrim=normArgs$sumTrim,
                doWeighting=normArgs$doWeighting,Acutoff=normArgs$Acutoff)
        }
    )
    if (output=="native") {
        if (is.null(geneData))
            return(NOISeq::readData(
            data=M,
            factors=data.frame(class=classes)
        )) # Class: CD
        else    
            return(NOISeq::readData(
                data=M,
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )) # Class: CD
    }
    else if (output=="matrix")
        return(as.matrix(round(M))) # Class: matrix
}

#' Normalization based on the NBPSeq package
#'
#' This function is a wrapper over DESeq normalization. It accepts a matrix of gene
#' counts (e.g. produced by importing an externally generated table of counts to
#' the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of NBPSeq normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"nbpseq")} for an example and how
#' you can modify it.
#' @param libsizeList an optional named list where names represent samples (MUST
#' be the same as the samples in \code{sampleList}) and members are the library
#' sizes (the sequencing depth) for each sample. If not provided, the default is
#' the column sums of the \code{geneCounts} matrix.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the NBPSeq native S4
#' object (a specific list). In the latter case it should be handled with suitable
#' NBPSeq methods.
#' @return A matrix with normalized counts or a list with the normalized counts
#' and other NBPSeq specific parameters.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' norm.dataMatrix <- normalizeNbpseq(dataMatrix,sampleList)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
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
    norm.factors <- estimate.norm.factors(geneCounts,libSizes=libSizes,
        method=normArgs$method)
    #if (normArgs$main.method=="nbpseq")
    #    nbData <- prepare.nbData(geneCounts,libSizes=libSizes,
    #    norm.factors=norm.factors)
    #else if (normArgs$main.method=="nbsmyth")
        nbData <- prepare.nbp(geneCounts,classes,libSizes=libSizes,
            norm.factors=norm.factors,thinning=normArgs$thinning)
    if (output=="native")
        return(nbData) # Class: list or nbp
    else if (output=="matrix") {
        #if (normArgs$main.method=="nbpseq") {
        #    normCounts <- matrix(0,nrow(geneCounts),ncol(geneCounts))
        #    for (i in 1:ncol(geneCounts))
        #        normCounts[,i] <- norm.factors[i]*geneCounts[,i]
        #}
        #else if (normArgs$main.method=="nbsmyth") 
            normCounts <- nbData$pseudo.counts
        return(as.matrix(round(normCounts))) # Class: matrix
    }
}

#' Normalization based on the ABSSeq package
#'
#' This function is a wrapper over ABSSeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param geneCounts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param normArgs a list of ABSSeq normalization parameters. See the result of
#' \code{getDefaults("normalization",} \code{"absseq")} for an example and how you
#' can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the ABSSeq native S4
#' object (ABSDataSet). In the latter case it should be handled with suitable
#' ABSSeq methods.
#' @return A matrix or a ABSDataSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplotBoxplot(dataMatrix,sampleList)
#'
#' norm.dataMatrix <- normalizeAbsseq(dataMatrix,sampleList)
#' diagplotBoxplot(norm.dataMatrix,sampleList)
#'}
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

#' Normalization based on the DSS package
#'
#' This function is a wrapper over DSS normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of DSS normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"DSS")} for an example and how you
#' can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the DSS native S4
#' object (ABSDataSet). In the latter case it should be handled with suitable
#' DSS methods.
#' @return A matrix or a ABSDataSet with normalized counts.
#' @author Dionysios Fanidis
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' dataMatrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(dataMatrix,sample.list)
#'
#' norm.dataMatrix <- normalize.absseq(dataMatrix,sample.list)
#' diagplot.boxplot(norm.dataMatrix,sample.list)
#'}
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
    if(class(geneCounts) == "data.frame") {
        geneCountsTmp <- geneCounts
        allCols <- 1:ncol(geneCountsTmp)
        samCols <- match(unlist(sampleList),colnames(geneCountsTmp))
        samCols <- samCols[which(!is.na(samCols))]
        geneCountsTmp <- geneCounts[,samCols]
        geneCountsTmp= as.matrix(geneCountsTmp)
    }
    if(class(geneCounts) == "matrix")
        geneCountsTmp <- geneCounts
    
    seqD <- newSeqCountSet(geneCountsTmp,design) # create the class
    seqD <- estNormFactors(seqD,method=normArgs$method) # estimate normalization factors
    
    if (output=="native")
        return(seqD) # Class: SeqCountSet
    else if (output=="matrix") {
        # I cannot extract normalized counts from DSS and thus I do this small 
        # dribble in order to output a matrix
        theDesign <- data.frame(condition=classes,row.names=colnames(seqD)) 
        cds <- newCountDataSet(as.matrix(round(assayData(seqD)$exprs)),
            theDesign$condition)
        DESeq::sizeFactors(cds) <- normalizationFactor(seqD)
        counts <- as.matrix(DESeq::counts(cds,normalized=TRUE))
        return(as.matrix(round(counts))) # Class: matrix
    }
}
