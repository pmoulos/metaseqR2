#' Filter gene expression based on exon counts
#'
#' This function performs the gene expression filtering based on exon read counts 
# and a set of exon filter rules. For more details see the main help pages of 
#' \code{\link{metaseqr}}.
#'
#' @param theCounts a named list created with the \code{\link{construct.gene.model}} 
#' function. See its help page for details.
#' @param geneData an annotation data frame usually obtained with 
#' \code{\link{getAnnotation}} containing the unique gene accession identifiers.
#' @param sampleList the list containing condition names and the samples under 
#' each condition.
#' @param exonFilters a named list with exon filters and their parameters. See 
#' the main help page of \code{\link{metaseqr}} for details.
#' @param restrict.cores in case of parallel execution of several subfunctions, 
#' the fraction of the available cores to use. In some cases if all available 
#' cores are used (\code{restrict.cores=1} and the system does not have sufficient 
#' RAM, the running machine might significantly slow down.
#' @return a named list with two members. The first member (\code{result} is a 
#' named list whose names are the exon filter names and its members are the filtered 
#' rownames of \code{geneData}. The second member is a matrix of binary flags 
#' (0 for non-filtered, 1 for filtered) for each gene. The rownames of the flag 
#' matrix correspond to gene ids.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' data("hg19.exonData",package="metaseqR")
#' exonCounts <- hg19.exonCounts
#' geneData <- getAnnotation("hg19","gene")
#' sampleList <- sampleList.hg19
#' exonFilters <- getDefaults("exonFilter")
#' theCounts <- construct.gene.model(exonCounts,sampleList,geneData)
#' filter.results <- filterExons(theCounts,geneData,sampleList,exonFilters)
#'}
filterExons <- function(theCounts,geneData,sampleList,exonFilters,rc=0.8) {
    exonFilterResult <- vector("list",length(exonFilters))
    names(exonFilterResult) <- names(exonFilters)
    flags <- matrix(0,nrow(geneData),1)
    rownames(flags) <- rownames(geneData)
    colnames(flags) <- c("MAE")
    theGenes <- as.character(geneData$gene_id)
    if (!is.null(exonFilters)) {
        for (xf in names(exonFilters)) {
            disp("Applying exon filter ",xf,"...")
            switch(xf,
                minActiveExons = {
                    pass <- vector("list",length(unlist(sampleList)))
                    names(pass) <- names(theCounts)
                    for (n in names(pass)) {
                        disp("  Checking read presence in exons for ",n,"...")
                        pass[[n]] <- theGenes
                        names(pass[[n]]) <- theGenes
                        pass[[n]] <- cmclapply(theCounts[[n]],function(x,f) {
                            if (length(x$count) == 1)
                                if (x$count[1]!=0)
                                    return(FALSE)
                                else
                                    return(TRUE)
                            else if (length(x$count) > 1 &&
                                length(x$count) <= f$exonsPerGene)
                                if (length(which(x$count!=0)) >= f$minExons)
                                    return(FALSE)
                                else
                                    return(TRUE)
                            else
                                if (length(which(x$count!=0)) >=
                                    ceiling(length(x$count)*f$frac))
                                    return(FALSE)
                                else
                                    return(TRUE)
                            },exonFilters$minActiveExons,rc=rc)
                        pass[[n]] <- do.call("c",pass[[n]])
                    }
                    passMatrix <- do.call("cbind",pass)
                    exonFilterResult[[xf]] <- theGenes[which(apply(
                        passMatrix,1,function(x) return(all(x))))]
                    flags[exonFilterResult[[xf]],"MAE"] <- 1
                }
                # More to come...
                # TODO: Write more rules based in exons
            )
        }
    }
    return(list(result=exonFilterResult,flags=flags))
}

#' Filter gene expression based on gene counts
#'
#' This function performs the gene expression filtering based on gene read counts 
#' and a set of gene filter rules. For more details see the main help pages of 
#' \code{\link{metaseqr}}.
#'
#' @param geneCounts a matrix of gene counts, preferably after the normalization 
#' procedure.
#' @param geneData an annotation data frame usually obtained with 
#' \code{\link{getAnnotation}} containing the unique gene accession identifiers.
#' @param geneFilters a named list with gene filters and their parameters. See 
#' the main help page of \code{\link{metaseqr}} for details.
#' @return a named list with three members. The first member (\code{result} is a 
#' named list whose names are the gene filter names and its members are the 
#' filtered rownames of \code{geneData}. The second member (\code{cutoff} is a 
#' named list whose names are the gene filter names and its members are the cutoff 
#' values corresponding to each filter. The third member is a matrix of binary
#' flags (0 for non-filtered, 1 for filtered) for each gene. The rownames of the 
#' flag matrix correspond to gene ids.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' data("mm9.geneData",package="metaseqR")
#' geneCounts <- mm9.geneCounts
#' sampleList <- mm9.sampleList
#' geneCounts <- normalizeEdger(as.matrix(geneCounts[,9:12]),sampleList)
#' geneData <- getAnnotation("mm9","gene")
#' geneFilters <- getDefaults("geneFilter","mm9")
#' filter.results <- filterGenes(geneCounts,geneData,geneFilters)
#'}
filterGenes <- function(geneCounts,geneData,geneFilters,sampleList) {
    geneFilterResult <- geneFilterCutoff <- vector("list",
        length(geneFilters))
    names(geneFilterResult) <- names(geneFilterCutoff) <-
        names(geneFilters)
    flags <- matrix(0,nrow(geneCounts),9)
    rownames(flags) <- rownames(geneCounts)
    colnames(flags) <- c("LN","AR","MD","MN","QN","KN","CM","BT","PR")
    for (gf in names(geneFilters)) {
        disp("Applying gene filter ",gf,"...")
        switch(gf,
            length = { # This is real gene length independently of exons
                if (!is.null(geneFilters$length)) {
                    geneFilterResult$length <- rownames(geneData)[which(
                        geneData$end - geneData$start <
                            geneFilters$length$length
                    )]
                    geneFilterCutoff$length <- geneFilters$length$length
                    flags[intersect(geneFilterResult$length,
                        rownames(geneCounts)),"LN"] <- 1
                    disp("  Threshold below which ignored: ",
                        geneFilters$length$length)
                }
                else
                    geneFilterCutoff$length <- NULL
            },
            avgReads = {
                if (!is.null(geneFilters$avgReads)) {
                    len <- attr(geneData,"geneLength")
                    if (!is.null(len))
                        len <- len[rownames(geneCounts)]
                    else {
                        gg <- geneData[rownames(geneCounts),]
                        len <- gg$end - gg$start
                    }
                    avgMat <- sweep(geneCounts,1,
                        len/geneFilters$avgReads$averagePerBp,"/")
                    qua <- max(apply(avgMat,2,quantile,
                        geneFilters$avgReads$quantile))
                    geneFilterResult$avgReads <- rownames(geneData)[which(
                        apply(avgMat,1,filterLow,qua))]
                    geneFilterCutoff$avgReads <- qua
                    flags[intersect(geneFilterResult$avgReads,
                        rownames(geneCounts)),"AR"] <- 1
                    disp("  Threshold below which ignored: ",qua)
                }
                else
                    geneFilterCutoff$avgReads <- NULL
            },
            expression = {
                if (!is.null(geneFilters$expression)) {
                    if (!is.null(geneFilters$expression$median) && 
                        geneFilters$expression$median) {
                        md <- median(geneCounts)
                        theDeadMedian <- rownames(geneCounts)[which(
                            apply(geneCounts,1,filterLow,md))]
                        disp("  Threshold below which ignored: ",md)
                    }
                    else
                        theDeadMedian <- md <- NULL
                    if (!is.null(geneFilters$expression$mean) &&
                        geneFilters$expression$mean) {
                        mn <- mean(geneCounts)
                        theDeadMean <- rownames(geneCounts)[which(apply(
                            geneCounts,1,filterLow,mn))]
                        disp("  Threshold below which ignored: ",mn)
                    }
                    else
                        theDeadMean <- mn <- NULL
                    if (!is.null(geneFilters$expression$quantile) &&
                        !is.na(geneFilters$expression$quantile)) {
                        qu <- quantile(geneCounts,
                            geneFilters$expression$quantile)
                        theDeadQuantile <- rownames(geneCounts)[which(
                            apply(geneCounts,1,filterLow,qu))]
                        disp("  Threshold below which ignored: ",qu)
                    }
                    else
                        theDeadQuantile <- qu <- NULL
                    if (!is.null(geneFilters$expression$known) &&
                        !is.na(geneFilters$expression$known)) {
                        # Think about the case of embedded
                        bioCut <- match(geneFilters$expression$known,
                            geneData$gene_name) 
                        bioCut <- bioCut[-which(is.na(bioCut))]
                        bioCutCounts <- as.vector(geneCounts[bioCut,])
                        theBioCut <- quantile(bioCutCounts,0.9)
                        theDeadKnown <- rownames(geneCounts)[which(apply(
                            geneCounts,1,filterLow,theBioCut))]
                        disp("  Threshold below which ignored: ",theBioCut)
                    }
                    else
                        theDeadKnown <- theBioCut <- NULL
                    if (!is.null(geneFilters$expression$custom) &&
                        !is.na(geneFilters$expression$custom)) {
                        # For future use
                        theDeadCustom <- NULL
                    }
                    else
                        theDeadCustom <- NULL
                    # Derive one common expression filter
                    geneFilterResult$expression$median <- theDeadMedian
                    geneFilterResult$expression$mean <- theDeadMean
                    geneFilterResult$expression$quantile <-
                        theDeadQuantile
                    geneFilterResult$expression$known <- theDeadKnown
                    geneFilterResult$expression$custom <- theDeadCustom
                    geneFilterCutoff$expression$median <- md
                    geneFilterCutoff$expression$mean <- mn
                    geneFilterCutoff$expression$quantile <- qu
                    geneFilterCutoff$expression$known <- theBioCut
                    geneFilterCutoff$expression$custom <- NULL
                    if (!is.null(theDeadMedian)) flags[theDeadMedian,
                        "MD"] <- 1
                    if (!is.null(theDeadMean)) flags[theDeadMean,"MN"] <- 1
                    if (!is.null(theDeadQuantile)) flags[theDeadQuantile,
                        "QN"] <- 1
                    if (!is.null(theDeadKnown)) flags[theDeadKnown,
                        "KN"] <- 1
                    if (!is.null(theDeadCustom)) flags[theDeadCustom,
                        "CM"] <- 1
                    #theDead <- list(theDeadMedian,theDeadMean,
                    #    theDeadQuantile,theDead.known,theDead.custom)
                    #geneFilterResult$expression <- Reduce("union",theDead)
                }
            },
            biotype = {
                if (!is.null(geneFilters$biotype)) {
                    filterOut <- names(which(unlist(geneFilters$biotype)))
                    # Necessary hack because of R naming system
                    if (length(grep("three_prime_overlapping_ncrna",
                        filterOut))>0) 
                        filterOut <- sub("three_prime_overlapping_ncrna",
                            "3prime_overlapping_ncrna",filterOut)
                    filterInd <- vector("list",length(filterOut))
                    names(filterInd) <- filterOut
                    for (bt in filterOut)
                        filterInd[[bt]] <- rownames(geneData)[which(
                            geneData$biotype==bt)]
                    geneFilterResult$biotype <- Reduce("union",filterInd)
                    geneFilterCutoff$biotype <- paste(filterOut,
                        collapse=", ")
                    disp("  Biotypes ignored: ",paste(filterOut,
                        collapse=", "))
                }
                else
                    geneFilterResult$biotype <- NULL
                if (!is.null(geneFilterResult$biotype) && 
                    length(geneFilterResult$biotype)>0) 
                    flags[geneFilterResult$biotype,"BT"] <- 1
            },
            presence = {
                if (!is.null(geneFilters$presence)) {
                    frac <- geneFilters$presence$frac
                    minc <- geneFilters$presence$minCount
                    pc <- geneFilters$presence$perCondition
                    
                    if (pc) {
                        npTmp <- namTmp <- vector("list",length(sampleList))
                        names(npTmp) <- names(namTmp) <- names(sampleList)
                        for (n in names(sampleList)) {
                            tmp <- geneCounts[,sampleList[[n]]]
                            nreq <- ceiling(frac*ncol(tmp))
                            npTmp[[n]] <- apply(tmp,1,function(x,m,n) {
                                w <- which(x>=m)
                                if (length(w)>0 && length(w)>n)
                                    return(FALSE)
                                return(TRUE)
                            },minc,nreq)
                            namTmp[[n]] <- rownames(tmp[which(npTmp[[n]]),,
                                drop=FALSE])
                        }
                        theDeadPresence <- Reduce("union",namTmp)
                    }
                    else {
                        nreq <- ceiling(frac*ncol(geneCounts))
                        notPresent <- apply(geneCounts,1,function(x,m,n) {
                            w <- which(x>=m)
                            if (length(w)>0 && length(w)>n)
                                return(FALSE)
                            return(TRUE)
                        },minc,nreq)
                        theDeadPresence <- 
                            rownames(geneCounts[which(notPresent),,
                                drop=FALSE])
                    }
                    
                    if (length(theDeadPresence)>0) {
                        geneFilterResult$presence <- theDeadPresence
                        flags[theDeadPresence,"PR"] <- 1
                        geneFilterCutoff$presence <- nreq
                        disp("  Threshold below which ignored: ",nreq)
                    }
                    else
                        geneFilterResult$presence <- 
                            geneFilterCutoff$presence <- NULL
                }
                else
                    geneFilterResult$presence <- NULL
            }
        )
    }
    return(list(result=geneFilterResult,cutoff=geneFilterCutoff,
        flags=flags))
}
