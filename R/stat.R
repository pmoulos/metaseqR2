statDeseq <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs) && !is(object,".CountDataSet"))
        statArgs <- getDefaults("statistics","deseq")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    classes <- asClassVector(sampleList)
    theDesign <- data.frame(condition=classes,row.names=colnames(object))
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    # Check if there is no replication anywhere
    if (all(vapply(sampleList,function(x) ifelse(length(x)==1,TRUE,FALSE),
        logical(1)))) {
        warnwrap("No replication detected! There is a possibility that ",
            "DESeq will fail to estimate dispersions...")
        methodDisp <- "blind"
        sharingModeDisp <- "fit-only"
        fitTypeDisp <- "local"
    }
    else {
        methodDisp <- statArgs$method
        sharingModeDisp <- statArgs$sharingMode
        fitTypeDisp <- statArgs$fitType
    }
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            cds <- object
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp,fitType=fitTypeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp,fitType=fitTypeDisp)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            #cds <- newCountDataSet(round(DESeq::counts(object,
            #   normalized=TRUE)),
            #    theDesign$condition)
            cds <- newCountDataSet(round(counts(object,normalized=TRUE)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            cds <- newCountDataSet(round(t(t(object$counts)/scl)*mean(scl)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            cds <- newCountDataSet(object,theDesign$condition)
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        nbData = { # Has been normalized with NBPSeq and main method was 
            # "nbpseq"
            cds <- newCountDataSet(as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),theDesign$condition)
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        nbp = { # Has been normalized with NBPSeq and main method was 
            # "nbsmyth"...
            cds <- newCountDataSet(as.matrix(round(object$pseudo.counts)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            cds <- newCountDataSet(as.matrix(round(excounts(object)),
                theDesign$condition))
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        },
        SeqCountSet = { # Has been normalized with DSS
            cds <- newCountDataSet(as.matrix(round(assayData(object)$exprs),
                theDesign$condition))
            #DESeq::sizeFactors(cds) <- normalizationFactor(object)
            sizeFactors(cds) <- normalizationFactor(object)
            #cds <- DESeq::estimateDispersions(cds,method=methodDisp,
            #    sharingMode=sharingModeDisp)
            cds <- estimateDispersions(cds,method=methodDisp,
                sharingMode=sharingModeDisp)
        }
    )
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        cons <- unique(unlist(con))
        if (length(con)==2) {
            res <- nbinomTest(cds,cons[1],cons[2])
            p[[conName]] <- res$pval
        }
        else {
            #cind <- match(cons,theDesign$condition)
            #if (any(is.na(cind)))
            #    cind <- cind[which(!is.na(cind))]
            cc <- names(unlist(con))
            cdsTmp <- cds[,cc]
            fit0 <- fitNbinomGLMs(cdsTmp,count~1)
            fit1 <- fitNbinomGLMs(cdsTmp,count~condition)
            p[[conName]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[conName]]) <- rownames(object)
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statDeseq2 <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","deseq2")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
        collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    
    classes <- asClassVector(sampleList)
    # Transform and rename classes it into a factor so as to use it to construct
    # a colData DataFrame
    conditions <- as.factor(classes)                 
    # colData constructed based on the sampleList given by the user 
    colData <- DataFrame(conditions)
    # Design constructed based on the colData made in the above line
    design <- as.formula(c("~", names(colData[1]))) 
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)

    switch(class(object)[1],
        DESeqDataSet = { # Has been normalized with DESeq2
            dds <- object
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType, 
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)                   
        },
        .CountDataSet = { # Has been normalized with DESeq
            #countData <- round(DESeq::counts(object, normalized=TRUE))
            countData <- round(counts(object, normalized=TRUE))
            dds <- DESeqDataSetFromMatrix(countData,colData,design=design,
                tidy=statArgs$tidy) 
            DESeq2::sizeFactors(dds) <- rep(1,ncol(countData))
            dds= DESeq2::estimateDispersions(dds,fitType=statArgs$fitType, 
                    maxit=statArgs$maxit,quiet=statArgs$quiet, 
                    modelMatrix=statArgs$modelMatrix)
        },
        DGEList = { # Has been normalized with edgeR              
            dds <- DESeqDataSetFromMatrix(object$counts,colData,
                design=design,tidy=statArgs$tidy)
            # edgeR has already normalized for library size. So what is done
            # is to pass size factors from edgeR to sizeFactors of DESeq2
            DESeq2::sizeFactors(dds) <- object$samples$norm.factors 
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType, 
                maxit=statArgs$maxit,quiet=statArgs$quiet, 
                modelMatrix=statArgs$modelMatrix)
        },
        SeqExpressionSet = { # Has been normalized with EDASeq
            # A matrix with the normalized counts (both within- and 
            # between-lane)
            countData <- normCounts(object) 
            dds <- DESeqDataSetFromMatrix(countData,colData,design=design,
                tidy=statArgs$tidy)
            # Because EDAseq has done library size norm, we set all 
            # sizeFactors equal to 1
            DESeq2::sizeFactors(dds) <- rep(1,ncol(dds)) 
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq or ABSSeq
            # I take the matrix ouput of these tools that contains 
            # normalized counts
            dds <- DESeqDataSetFromMatrix(object,colData,design=design,
                tidy=statArgs$tidy)
            # Because EDAseq has done library size norm (and NOISeq too, 
            # and EdgeR too), we set all sizeFactors equal 1
            DESeq2::sizeFactors(dds) <- rep(1,ncol(dds))
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)
        },
        nbData = { # Normalized with NBPSeq and main method was "nbpseq".
            # Using sweep to multiply raw counts with the respective 
            # norm.factors -> normCounts generated
            dds <- DESeqDataSetFromMatrix(
                round(sweep(object$counts,2,object$norm.factors,"*")),
                colData,design=design,tidy=statArgs$tidy) 
            DESeq2::sizeFactors(dds)= rep(1,ncol(dds))
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)
        },
        nbp = { # Normalized with NBPSeq and main method was "nbsmyth"; 
            #"nbsmyth" refers to prepare.nbp. 
            # We take the pseudocounts of the nbp object because they are the 
            # normalized ones for the lib.size
            dds <- DESeqDataSetFromMatrix(round(object$pseudo.counts),
                colData,design=design,tidy=statArgs$tidy) 
            DESeq2::sizeFactors(dds)= rep(1,ncol(dds))
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            dds <- DESeqDataSetFromMatrix(round(excounts(object)),colData,
                design=design,tidy=statArgs$tidy)
            DESeq2::sizeFactors(dds) <- rep(1,ncol(dds))
            dds <- DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
            modelMatrix=statArgs$modelMatrix)
        },
        SeqCountSet = { # Has been normalized with DSS
            dds <- DESeqDataSetFromMatrix(round(assayData(object)$exprs),
                colData,design=design,tidy=statArgs$tidy)
            DESeq2::sizeFactors(dds)=normalizationFactor(object)
            dds <-  DESeq2::estimateDispersions(dds,fitType=statArgs$fitType,
                maxit=statArgs$maxit,quiet=statArgs$quiet,
                modelMatrix=statArgs$modelMatrix)
        }
    )
    
    for (conName in names(contrastList)) {
        disp(" Contrast: ",conName)
        con <- contrastList[[conName]]
        cons <- unique(unlist(con))
        # nbinomWaldTest performs ALL possible pairwise comparisons between the 
        # levels of the factor specified in the dds's design.
        dds <- nbinomWaldTest(dds,betaPrior=statArgs$betaPrior,
            modelMatrix=statArgs$modelMatrix,betaTol=statArgs$betaTol, 
            maxit=statArgs$maxit,useOptim=statArgs$useOptim,
            quiet=statArgs$quiet,useT=statArgs$useT, 
            useQR=statArgs$useQR)
        # Then, one has just to ask for the results of the comparison he wants.
        if (length(con)==2) {
            # Because in conditions factor may be >2 levels, in a pairwise 
            # comparison I have to specify which results to show
            res <- DESeq2::results(dds,contrast=c("conditions",cons[[1]],
                cons[[2]]))
            p[[conName]] <- res$pvalue
        }
        else { # DE analysis of >2 levels in conditions factor
            cc <- names(unlist(con))
            # This assignment, except from keeping only the samples we want to 
            # compare, it also keeps only the respective levels under 
            # conditions!
            ddsTmp <- dds[,cc] 
            ddsTmp <- nbinomLRT(ddsTmp,full=design,reduced=~1,
                betaTol=statArgs$betaTol,maxit=statArgs$maxit,
                useOptim=statArgs$useOptim,quiet=statArgs$quiet, 
                useQR=statArgs$useQR)
            # Using the LRT to compare all the levels inside a factor 
            # (like an ANOVA)
            res <- DESeq2::results(ddsTmp)
            p[[conName]] <- res$pvalue
        }
        names(p[[conName]]) <- rownames(object)
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    
    return(p)
}                 

statEdger <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","edger")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    classes <- asClassVector(sampleList)
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #dge <- DGEList(counts=DESeq::counts(object,normalized=TRUE),
            #    group=classes)
            dge <- DGEList(counts=counts(object,normalized=TRUE),group=classes)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            dge <- DGEList(counts=DESeq2::counts(object,normalized=TRUE),
                group=classes)
        },  
        DGEList = { # Has been normalized with edgeR
            dge <- object    
        },
        matrix = { # Has been normalized with EDASeq or NOISeq or ABSSeq
            dge <- DGEList(object,group=classes)
        },
        nbData = { # Has been normalized with NBPSeq and method was "nbpseq"
            dge <- DGEList(counts=as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),group=classes)
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            dge <- DGEList(counts=as.matrix(round(object$pseudo.counts)),
                group=classes)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            dge <- DGEList(counts=as.matrix(round(excounts(object))),
                group=classes)
        },
        SeqCountSet = { # Has been normalized with DSS
            dge <- DGEList(counts=as.matrix(round(assayData(object)$exprs)),
                norm.factors=normalizationFactor(object),group=classes)
        }
    )
    # Dispersion estimate step
    # Check if there is no replication anywhere
    repli = TRUE
    if (all(vapply(sampleList,function(x) ifelse(length(x)==1,TRUE,FALSE),
        logical(1)))) {
        warnwrap("No replication when testing with edgeR! Consider using ",
            "another statistical test or just performing empirical analysis. ",
            "Setting bcv to 0.2...")
        repli <- FALSE
        bcv <- 0.2
    }
    if (repli) {
        if (statArgs$main.method=="classic") {
            dge <- estimateCommonDisp(dge,rowsum.filter=statArgs$rowsum.filter)
            dge <- estimateTagwiseDisp(dge,prior.df=statArgs$prior.df,
                trend=statArgs$trend,span=statArgs$span,
                    method=statArgs$tag.method,
                grid.length=statArgs$grid.length,
                    grid.range=statArgs$grid.range)
        }
        else if (statArgs$main.method=="glm") {
            design <- model.matrix(~0+classes,data=dge$samples)
            dge <- estimateGLMCommonDisp(dge,design=design,
                offset=statArgs$offset,
                method=statArgs$glm.method,subset=statArgs$subset,
                AveLogCPM=statArgs$AveLogCPM)
            dge <- estimateGLMTrendedDisp(dge,design=design,
                offset=statArgs$offset,
                method=statArgs$trend.method,AveLogCPM=statArgs$AveLogCPM)
            dge <- estimateGLMTagwiseDisp(dge,design=design,
                offset=statArgs$offset,
                dispersion=statArgs$dispersion,prior.df=statArgs$prior.df,
                span=statArgs$span,AveLogCPM=statArgs$AveLogCPM)
        }
    }
    # Actual statistical test
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        if (length(con)==2) {
            if (repli) {
                if (statArgs$main.method=="classic") {
                    res <- exactTest(dge,pair=unique(unlist(con)))
                }
                else if (statArgs$main.method=="glm") {
                    s <- unlist(con)
                    us <- unique(s)
                    ms <- match(names(s),rownames(dge$samples))
                    if (any(is.na(ms)))
                        ms <- ms[which(!is.na(ms))]
                    design <- model.matrix(~0+s,data=dge$samples[ms,])
                    colnames(design) <- us
                    #fit <- glmFit(dge[,ms],design=design,
                    #    offset=statArgs$offset,
                    #    weights=statArgs$weights,lib.size=statArgs$lib.size,
                    #    prior.count=statArgs$prior.count,
                    #    start=statArgs$start,method=statArgs$method)
                    fit <- glmFit(dge[,ms],design=design,
                        prior.count=statArgs$prior.count,
                        start=statArgs$start,method=statArgs$method)
                    co <- makeContrasts(paste(us[2],us[1],sep="-"),
                        levels=design)
                    lrt <- glmLRT(fit,contrast=co)
                    res <- topTags(lrt,n=nrow(dge))
                }
            }
            else {
                if (statArgs$main.method=="classic") {
                    res <- exactTest(dge,pair=unique(unlist(con)),
                        dispersion=bcv^2)
                }
                else if (statArgs$main.method=="glm") {
                    s <- unlist(con)
                    us <- unique(s)
                    ms <- match(names(s),rownames(dge$samples))
                    if (any(is.na(ms)))
                        ms <- ms[which(!is.na(ms))]
                    design <- model.matrix(~0+s,data=dge$samples[ms,])
                    colnames(design) <- us
                    #fit <- glmFit(dge[,ms],design=design,
                    #    offset=statArgs$offset,weights=statArgs$weights,
                    #    lib.size=statArgs$lib.size,
                    #    prior.count=statArgs$prior.count,
                    #    start=statArgs$start,
                    #    method=statArgs$method,dispersion=bcv^2)
                    fit <- glmFit(dge[,ms],design=design,dispersion=bcv^2,
                        prior.count=statArgs$prior.count,start=statArgs$start,
                        method=statArgs$method)
                    co <- makeContrasts(paste(us[2],us[1],sep="-"),
                        levels=design)
                    lrt <- glmLRT(fit,contrast=co)
                    res <- topTags(lrt,n=nrow(dge))
                }
            }
        }
        else { # GLM only
            s <- unlist(con)
            us <- unique(s)
            #design <- model.matrix(~0+s,data=dge$samples) # Ouch!
            ms <- match(names(s),rownames(dge$samples))
            if (any(is.na(ms)))
                ms <- ms[which(!is.na(ms))]
            design <- model.matrix(~s,data=dge$samples[ms,])
            if (repli)
                #fit <- glmFit(dge[,ms],design=design,offset=statArgs$offset,
                #    weights=statArgs$weights,lib.size=statArgs$lib.size,
                #    prior.count=statArgs$prior.count,start=statArgs$start,
                #    method=statArgs$method)
                fit <- glmFit(dge[,ms],design=design,start=statArgs$start,
                    prior.count=statArgs$prior.count)
            else
                #fit <- glmFit(dge[,ms],design=design,offset=statArgs$offset,
                #    weights=statArgs$weights,lib.size=statArgs$lib.size,
                #    prior.count=statArgs$prior.count,start=statArgs$start,
                #    method=statArgs$method,dispersion=bcv^2)
                #    dispersion=bcv^2)
                fit <- glmFit(dge[,ms],design=design,dispersion=bcv^2,
                    prior.count=statArgs$prior.count,start=statArgs$start)
                #lrt <- glmLRT(fit,coef=2:ncol(fit$design))
                lrt <- glmLRT(fit,coef=2:ncol(fit$design))
                res <- topTags(lrt,n=nrow(dge))
        }
        p[[conName]] <- res$table[,"PValue"]
        names(p[[conName]]) <- rownames(res$table)
        p[[conName]] <- p[[conName]][rownames(dge)]
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statLimma <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","limma")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    classes <- asClassVector(sampleList)
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #dge <- DGEList(DESeq::counts(object,normalized=TRUE),group=classes)
            dge <- DGEList(counts(object,normalized=TRUE),group=classes)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            dge <- DGEList(DESeq2::counts(object,normalized=TRUE),group=classes)
        },
        DGEList = { # Has been normalized with edgeR
            dge <- object
        },
        matrix = { # Has been normalized with EDASeq or NOISeq or EDASeq
            dge <- DGEList(object,group=classes)
        },
        nbData = { # Has been normalized with NBPSeq and method was "nbpseq"
            dge <- DGEList(counts=as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),group=classes)
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            dge <- DGEList(counts=as.matrix(round(object$pseudo.counts)),
                group=classes)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            dge <- DGEList(counts=as.matrix(round(excounts(object))),
                group=classes)
        },
        SeqCountSet = { # Has been normalized with DSS
            dge <- DGEList(counts=as.matrix(round(assayData(object)$exprs)),
                norm.factors=normalizationFactor(object),group=classes)
        }
    )
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        s <- unlist(con)
        us <- unique(s)
        ms <- match(names(s),rownames(dge$samples))
        if (any(is.na(ms)))
            ms <- ms[which(!is.na(ms))]
        if (length(con)==2) {
            design <- model.matrix(~0+s,data=dge$samples[ms,])
            colnames(design) <- us
            vom <- voom(dge[,ms],design,
                normalize.method=statArgs$normalize.method)
            fit <- lmFit(vom,design)
            fit <- eBayes(fit)
            co <- makeContrasts(contrasts=paste(us[2],us[1],sep="-"),
                levels=design)
            fit <- eBayes(contrasts.fit(fit,co))
            p[[conName]] <- fit$p.value[,1]
        }
        else {
            design <- model.matrix(~s,data=dge$samples[ms,])
            vom <- voom(dge[,ms],design,
                normalize.method=statArgs$normalize.method)
            fit <- lmFit(vom,design)
            fit <- eBayes(fit)
            res <- topTable(fit,coef=2:ncol(fit$design),number=nrow(vom))
            p[[conName]] <- res[,"P.Value"]
            names(p[[conName]]) <- rownames(res)
            p[[conName]] <- p[[conName]][rownames(dge)]
        }
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statNoiseq <- function(object,sampleList,contrastList=NULL,statArgs=NULL,
    geneData=NULL,logOffset=1) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","noiseq")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    if (!is.null(geneData) && is(geneData,"GenomicRanges")) {
        gl <- NULL
        if (!is.null(attr(geneData,"geneLength")))
            gl <- attr(geneData,"geneLength")
        geneData <- as.data.frame(geneData)
        geneData <- geneData[,c(1,2,3,6,7,5,8,9)]
        if (!is.null(gl))
            attr(geneData,"geneLength") <- gl
    }
    if (is.null(geneData)) {
        gcContent <- NULL
        chromosome <- NULL
        biotype <- NULL
        geneLength <- NULL
    }
    else {
        gcContent <- geneData$gc_content
        if (is.null(geneData$gc_content))
            gcContent <- rep(0.5,nrow(geneData))
        biotype <- as.character(geneData$biotype)
        names(gcContent) <- names(biotype) <- rownames(geneData)
        if (is.null(attr(geneData,"geneLength")))
            geneLength <- NULL
        else {
            geneLength <- attr(geneData,"geneLength")
            names(geneLength) <- rownames(geneData)
        }
    }
    classes <- asClassVector(sampleList)
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                #data=DESeq::counts(object,normalized=TRUE),
                data=counts(object,normalized=TRUE),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=DESeq2::counts(object,normalized=TRUE),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            dm <- round(t(t(object$counts)/scl)*mean(scl))            
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=dm,
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        ExpressionSet = { # Has been normalized with NOISeq
            nsObj <- object
        },
        matrix = { # Has been normalized with EDASeq
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=object,
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        nbData = { # Has been normalized with NBPSeq and main method "nbpseq"
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=as.matrix(round(sweep(object$counts,2,
                    object$norm.factors,"*"))),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=as.matrix(round(object$pseudo.counts)),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                data=as.matrix(round(excounts(object))),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        SeqCountSet = { # Has been normalized with DSS; 
            # Because NOISeq class does not have any slot for norm.factors and 
            # because DSS doesn't have any function to return the normalized 
            # counts on their own, I will transform the SeqCountSet into a cds 
            # object and then I will get the normalized counts
            theDesign <- data.frame(condition=classes,
                row.names=colnames(object))
            cds <- newCountDataSet(as.matrix(round(assayData(object)$exprs)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- normalizationFactor(object)
            sizeFactors(cds) <- normalizationFactor(object)
            #nsObj <- NOISeq::readData(
            nsObj <- readData(
                #data=DESeq::counts(cds,normalized=TRUE),
                data=counts(cds,normalized=TRUE),
                length=geneLength,
                gc=gcContent,
                chromosome=geneData[,c(1,2,3)],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        }
    )
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        if (length(con)==2) {
            statArgs$conditions=unique(unlist(con))
            if (any(vapply(sampleList,function(x) ifelse(length(x)==1,
                TRUE,FALSE),logical(1))))
                # At least one condition does not have replicates
                statArgs$replicates <- "no" 
            if (statArgs$replicates %in% c("technical","no"))
                res <- noiseq(nsObj,k=logOffset,norm="n",
                    replicates=statArgs$replicates,factor=statArgs$factor,
                    conditions=statArgs$conditions,pnr=statArgs$pnr,
                    nss=statArgs$nss,v=statArgs$v,lc=statArgs$lc)
            else
                res <- noiseqbio(nsObj,k=logOffset,norm="n",
                    nclust=statArgs$nclust,factor=statArgs$factor,
                    lc=statArgs$lc,conditions=statArgs$conditions,
                    r=statArgs$r,adj=statArgs$adj,a0per=statArgs$a0per,
                    cpm=statArgs$cpm,
                    filter=statArgs$filter,depth=statArgs$depth,
                    cv.cutoff=statArgs$cv.cutoff)
            # Beware! This is not the classical p-value!
            p[[conName]] <- 1 - res@results[[1]]$prob
        }
        else {
            warnwrap(paste("NOISeq differential expression algorithm does not ",
                "support multi-factor designs (with more than two conditions ",
                "to be compared)! Switching to DESeq for this comparison: ",
                conName))
            M <- assayData(nsObj)$exprs
            cds <- newCountDataSet(round(M),data.frame(condition=unlist(con),
                row.names=names(unlist(con))))
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method="blind",
            #    sharingMode="fit-only")
            cds <- estimateDispersions(cds,method="blind",
                sharingMode="fit-only")
            fit0 <- fitNbinomGLMs(cds,count~1)
            fit1 <- fitNbinomGLMs(cds,count~condition)
            p[[conName]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[conName]]) <- rownames(nsObj)
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statBayseq <- function(object,sampleList,contrastList=NULL,statArgs=NULL,
    libsizeList=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","bayseq")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    classes <- asClassVector(sampleList)
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #bayesData <- DESeq::counts(object,normalized=TRUE)
            bayesData <- counts(object,normalized=TRUE)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            bayesData <- DESeq2::counts(object,normalized=TRUE)
        },  
        DGEList = { # Has been normalized with edgeR
            scl <- object$samples$lib.size * object$samples$norm.factors
            bayesData <- round(t(t(object$counts)/scl)*mean(scl))
        },
        matrix = { # Has been normalized with EDASeq or NOISeq or ABSSeq
            bayesData <- object
        },
        nbData = {
            bayesData <- as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*")))
        },
        nbp = {
            bayesData <- as.matrix(round(object$pseudo.counts))
        },
        ABSDataSet = {
            bayesData <- as.matrix(round(excounts(object)))
        },
        SeqCountSet = { # Again the same trick I did back in NOISeq
            theDesign <- data.frame(condition=classes,
                row.names=colnames(object))
            cds <- newCountDataSet(as.matrix(round(assayData(object)$exprs)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- normalizationFactor(object)
            sizeFactors(cds) <- normalizationFactor(object)
            #bayesData <- as.matrix(DESeq::counts(cds,normalized=TRUE))
            bayesData <- as.matrix(counts(cds,normalized=TRUE))
        }
    )
    CD <- new("countData",data=bayesData,replicates=classes)
    CD@annotation <- data.frame(name=rownames(bayesData))
    if (is.null(libsizeList))
        baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
    else
        baySeq::libsizes(CD) <- unlist(libsizeList)
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        #cd <- CD[,names(unlist(con))]
        cd <- CD[,match(names(unlist(con)),colnames(CD@data))]
        if (length(con)==2)
            baySeq::groups(cd) <- list(NDE=rep(1,length(unlist(con))),
                DE=c(rep(1,length(con[[1]])),rep(2,length(con[[2]]))))
        else
            baySeq::groups(cd) <- list(NDE=rep(1,length(unlist(con))),
                DE=unlist(con,use.names=FALSE)) # Maybe this will not work
        baySeq::replicates(cd) <- as.factor(classes[names(unlist(con))])
        cd <- baySeq::getPriors.NB(cd,samplesize=statArgs$samplesize,
            samplingSubset=statArgs$samplingSubset,
            equalDispersions=statArgs$equalDispersions,
            estimation=statArgs$estimation,zeroML=statArgs$zeroML,
            consensus=statArgs$consensus,cl=statArgs$cl)
        cd <- baySeq::getLikelihoods(cd,pET=statArgs$pET,
            marginalise=statArgs$marginalise,subset=statArgs$subset,
            priorSubset=statArgs$priorSubset,bootStraps=statArgs$bootStraps,
            conv=statArgs$conv,nullData=statArgs$nullData,
            returnAll=statArgs$returnAll,returnPD=statArgs$returnPD,
            discardSampling=statArgs$discardSampling,cl=statArgs$cl)
        tmp <- baySeq::topCounts(cd,group="DE",number=nrow(cd))
        #p[[conName]] <- 1 - as.numeric(tmp[,"Likelihood"])
        p[[conName]] <- 1 - as.numeric(tmp[,"likes"])
        names(p[[conName]]) <- as.character(tmp$name)
        p[[conName]] <- p[[conName]][rownames(CD@data)]
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statNbpseq <- function(object,sampleList,contrastList=NULL,statArgs=NULL,
    libsizeList=NULL) {
    if (is.null(statArgs) && !is(object,"list"))
        statArgs <- getDefaults("statistics","nbpseq")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    classes <- asClassVector(sampleList)
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #counts <- round(DESeq::counts(object,normalized=TRUE))
            counts <- round(counts(object,normalized=TRUE))
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            counts <- round(DESeq2::counts(object,normalized=TRUE))
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            counts <- round(t(t(object$counts)/scl)*mean(scl))
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq or ABSSeq
            counts <- object
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        nbData = { # Has been normalized with NBPSeq
            object$counts <- as.matrix(object$counts)
            nbData <- object
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(nbData$counts[,n])
            }
            libSizes <- unlist(libsizeList)
            nbData$pseudo.lib.sizes=rep(1e+7,dim(object$counts)[2])
        },
        nbp = { # Same...
            object$pseudo.counts <- as.matrix(object$pseudo.counts)
            nbData <- object
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(nbData$counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            counts <- as.matrix(round(excounts(object)))
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        },
        SeqCountSet = { # Again the trick of NOISeq
            theDesign <- data.frame(condition=classes,
                row.names=colnames(object)) 
            cds <- newCountDataSet(as.matrix(round(assayData(object)$exprs)),
                theDesign$condition)
            #DESeq::sizeFactors(cds) <- normalizationFactor(object)
            sizeFactors(cds) <- normalizationFactor(object)
            #counts <- as.matrix(DESeq::counts(cds,normalized=TRUE))
            counts <- as.matrix(counts(cds,normalized=TRUE))
            if (is.null(libsizeList)) {
                libsizeList <- vector("list",length(classes))
                names(libsizeList) <- unlist(sampleList,use.names=FALSE)
                for (n in names(libsizeList))
                    libsizeList[[n]] <- sum(counts[,n])
            }
            libSizes <- unlist(libsizeList)
        }
    )
    # To avoid repeating the following chunk in the above
    if ((!is(object,"list") && !is(object,"nbp")) || is(object,"DGEList")) {
        #if (statArgs$main.method=="nbpseq") {
        #    nbData <- list(
        #        counts=as.matrix(counts),
        #        lib.sizes=libSizes,
        #        norm.factors=rep(1,dim(counts)[2]),
        #        eff.lib.sizes=libSizes*rep(1,dim(counts)[2]),
        #        rel.frequencies=as.matrix(sweep(counts,2,
        #            libSizes*rep(1,dim(counts)[2]),"/")),
        #        tags=matrix(row.names(counts),dim(counts)[1],1)
        #    )
        #}
        #else if (statArgs$main.method=="nbsmyth") {
        #    nbData <- new("nbp",list(
        #        counts=as.matrix(counts),
        #        lib.sizes=libSizes,
        #        grp.ids=classes,
        #        eff.lib.sizes=libSizes*rep(1,dim(counts)[2]),
        #        pseudo.counts=as.matrix(counts),
        #        pseudo.lib.sizes=colSums(as.matrix(counts))*rep(1,
        #            dim(counts)[2])
        #    ))
            nbData <- list(
                counts=as.matrix(counts),
                lib.sizes=libSizes,
                grp.ids=classes,
                eff.lib.sizes=libSizes*rep(1,dim(counts)[2]),
                pseudo.counts=as.matrix(counts),
                #pseudo.lib.sizes=colSums(as.matrix(counts)) *
                #    rep(1,dim(counts)[2])
                pseudo.lib.sizes=rep(1e+7,dim(counts)[2])
            )
            class(nbData) <- "nbp"
        #}
    }
    for (conName in names(contrastList)) {
        disp("  Contrast: ", conName)
        con <- contrastList[[conName]]
        cons <- unique(unlist(con))
        if (length(con)==2) {
            #if (statArgs$main.method=="nbpseq") {
            #    dispersions <- estimate.dispersion(nbData,
            #    model.matrix(~classes),
            #        method=statArgs$method$nbpseq)
            #    res <- test.coefficient(nbData,dispersion=dispersions,
            #        x=model.matrix(~classes),beta0=c(NA,0),
            #        tests=statArgs$tests,
            #        alternative=statArgs$alternative,print.level=1)
            #    #res <- nb.glm.test(nbData$counts,x=model.matrix(~classes),
            #        beta0=c(NA,0),libSizes=libSizes,
            #    #    dispersion.method=statArgs$method$nbpseq,
            #         tests=statArgs$tests)
            #    p[[conName]] <- res[[statArgs$tests]]$p.values
            #    #p[[conName]] <- res$test[[statArgs$tests]]$p.values
            #}
            #else if (statArgs$main.method=="nbsmyth") {
                obj <- suppressWarnings(estimate.disp(nbData,
                    model=statArgs$model$nbsmyth,print.level=0))
                obj <- exact.nb.test(obj,cons[1],cons[2],print.level=0)
                p[[conName]] <- obj$p.values
            #}
        }
        else {
            warnwrap(paste("NBPSeq differential expression algorithm does not ",
                "support ANOVA-like designs with more than two conditions to ",
                "be compared! Switching to DESeq for this comparison:",
                conName))
            cds <- newCountDataSet(nbData$counts,
                data.frame(condition=unlist(con),
                row.names=names(unlist(con))))
            #DESeq::sizeFactors(cds) <- rep(1,ncol(cds))
            sizeFactors(cds) <- rep(1,ncol(cds))
            #cds <- DESeq::estimateDispersions(cds,method="blind",
            #    sharingMode="fit-only")
            cds <- estimateDispersions(cds,method="blind",
                sharingMode="fit-only")
            fit0 <- fitNbinomGLMs(cds,count~1)
            fit1 <- fitNbinomGLMs(cds,count~condition)
            p[[conName]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[conName]]) <- rownames(nbData$counts)
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}

statAbsseq <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","absseq")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
        collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)

    # In statAbsseq "switch()" is inside "for" loop. Actually, there are two 
    # different "switch()" sections inside the loop. This happens because this 
    # tool takes only 2-level factors for analysis. There is a feature where we
    # can give a factor with more than 2 levels, but then the tool-specific
    # object must be created without the "groups" argument (see tools 
    # reference manual and vignette) and the function called for differential 
    # expression differs from that of the 2-level factors. This is why we place 
    # switch inside the "for" loop, in order to give the user the freedom to use
    # the tool to check between 2 or more levels of a factor during the same 
    # analysis.
    
    for(conName in names(contrastList)){
    disp(" Contrast: ",conName)
    con <- contrastList[[conName]]
    cons= unique(unlist(con))
        if(length(con)==2){
            # Checking if we have a pairwise comparison, because then we must 
            # resize a sampleList to only the 2 levels of interest and pass them
            # in the groups argument of ABSDataSet().
            sampleListTmp <- sampleList[cons]
            # In any other case, of length(groups)!=2 it will through an ERROR
            classes <- asClassVector(sampleListTmp)
            # Sample description as wanted by ABSDataSet for pairwise 
            # comparisons
            groupsTmp <- as.factor(classes)

            # Check for replication
            # We cannot do anything for the case that there is no replication, 
            # but to call separately callParameterWithoutReplicates (which by
            # the way has no other args except from the ABSDataSet) instead of 
            # callParameter and then continue with callDEs
            # All the above are automatically done by ABSSeq()
            # However, ABSSeqlm() cannot run without providing replicates (see
            # further down for another comment)
            switch(class(object)[1],
                ABSDataSet = {
                    abs <- ABSDataSet(excounts(object)[,names(unlist(con))],
                        groupsTmp,normMethod="user",
                        # "user" so as to pass my own sizeFactors
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                .CountDataSet = { # Has been normalized with DESeq
                    #countsTmp <- DESeq::counts(object)
                    countsTmp <- counts(object)
                    abs <- ABSDataSet(countsTmp[,names(unlist(con))],
                        # "user" so as to pass my own sizeFactors
                        groupsTmp,normMethod="user",
                        sizeFactor=sizeFactors(object)[names(unlist(con))], 
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                        # Here I used the DESeq's raw counts and size factors
                        # because that way I got a little better p-values
                },
                DESeqDataSet = { # Has been normalized with DESeq2
                    countsTmp <- DESeq2::counts(object)
                    abs <- ABSDataSet(countsTmp[,names(unlist(con))],
                        # "user" so as to pass my own sizeFactors
                        groupsTmp,normMethod= "user",
                        sizeFactor=sizeFactors(object)[names(unlist(con))], 
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC) 
                        # Here I used the DESeq's raw counts and size factors
                        #  because that way I got a little better p-values 
                },
                DGEList = { # Has been normalized with edgeR                
                    abs <- ABSDataSet(object$counts[,names(unlist(con))],
                        groupsTmp,normMethod="user", 
                        sizeFactor=object$samples[names(unlist(con)),3],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                SeqExpressionSet = { # Has been normalized with EDASeq
                    abs <- ABSDataSet(normCounts(object)[,
                        names(unlist(con))],groupsTmp,normMethod= "user", 
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },  
                matrix = { # Has been normalized with EDASeq or NOISeq or Edger
                    abs <- ABSDataSet(object[,names(unlist(con))],groupsTmp,
                        normMethod= "user",
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC) 
                        # Because EDAseq has done library size norm (and 
                        # NOISeq too), we set all sizeFactors equal to 1
                },              
                nbData = { # Has been normalized with NBPSeq and main method was
                    # "nbpseq"; 
                    # "nbpseq" refers to prepare.nbData that seems to have 
                    # stopped running in metaseqR (see norm.R lines 366-368)
                    normCounts=round(sweep(object$counts,2,
                        object$norm.factors,"*"))
                    abs <- ABSDataSet(object[,names(unlist(con))],groupsTmp,
                        normMethod="user", 
                        sizeFactor= rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                nbp = { # Has been normalized with NBPSeq and main method was 
                    #"nbsmyth"; 
                    #"nbsmyth" refers to prepare.nbp (see norm.R lines 370-371) 
                    abs <- ABSDataSet(round(object$pseudo.counts[,
                        names(unlist(con))]),groupsTmp,normMethod= "user", 
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC) 
                        # I take the pseudocounts of the nbp object because 
                        # they are the normalized ones for the lib.size
                },
                SeqCountSet = {
                    # isolate normalization factors so as to give them names and
                    # be able to subset them at will
                    nF <- normalizationFactor(object) 
                    names(nF) <- names(unlist(sampleList))
                    abs <- ABSDataSet(round(
                        assayData(object)$exprs[,names(unlist(con))]),
                        groupsTmp,normMethod= "user", 
                        sizeFactor=nF[names(unlist(con))],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                }
            )
        
            abs <- ABSSeq(abs,adjmethod=statArgs$adjmethod,
                replaceOutliers=statArgs$replaceOutliers, 
                useaFold=statArgs$useaFold,quiet= statArgs$quiet) 
                # It linearly calls normalFactors(),callParameter(),callDEs(); 
                # normalFactors() will use those given.
            res <- ABSSeq::results(abs)
            p[[conName]] <- res[,"pvalue"]
            names(p[[conName]]) <- rownames(res) 
            # ATTENTION: p[[con.names]] had to be named by the subsetted count 
            # table of each contrast and not at the end. If not thus there was a
            # problem with line 2975 of main.R
        }
        else { 
            # Here we build the ABSDataSet object WITHOUT the groups argument so
            # as to perform ANOVA. Otherwise, if I try to build a groups
            # argument with more than 2 levels to run the ANOVA, it will through
            # an ERROR 
            sampleListTmp <- sampleList[cons]
            classes <- asClassVector(sampleListTmp)
            groupsTmp <- as.factor(classes)
            design <- model.matrix(~0+groupsTmp)

            switch(class(object)[1],
                ABSDataSet = {
                    abs <- ABSDataSet(excounts(object)[,names(unlist(con))],
                        normMethod="user", 
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                .CountDataSet = { # Has been normalized with DESeq
                    #countsTmp <- DESeq::counts(object)
                    countsTmp <- counts(object)
                    abs <- ABSDataSet(countsTmp[,names(unlist(con))],
                        normMethod="user", 
                        sizeFactor=sizeFactors(object)[names(unlist(con))],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                DESeqDataSet = { # Has been normalized with DESeq2
                        countsTmp= DESeq2::counts(object)
                    abs <- ABSDataSet(countsTmp[,names(unlist(con))],
                        normMethod="user", 
                        sizeFactor=sizeFactors(object)[names(unlist(con))],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                        # Here we use the DESeq's raw counts and size factors
                        # because that way I got a little better p-values 
                },
                DGEList = { # Has been normalized with edgeR                
                    abs <- ABSDataSet(object$counts[,names(unlist(con))],
                        normMethod="user",
                        sizeFactor=object$samples[names(unlist(con)),3],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                SeqExpressionSet = { # Has been normalized with EDASeq
                    abs <- ABSDataSet(
                        EDASeq::normCounts(object)[,names(unlist(con))],
                        normMethod="user",
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)              
                },  
                matrix = { # Has been normalized with EDASeq or NOISeq
                    abs <- ABSDataSet(object[,names(unlist(con))],
                        normMethod="user", 
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC) 
                        # Because EDAseq has done library size norm (and 
                        # NOISeq too), we set all sizeFactors equal to 1
                },
                nbData = { 
                    # Has been normalized with NBPSeq and main method was 
                    # "nbpseq"; 
                    # "nbpseq" refers to prepare.nbData that seems to have 
                    # stopped running in metaseqR (see norm.R lines 366-368)
                    normCounts <- round(sweep(object$counts,2,
                        object$norm.factors,"*"))
                    abs <- ABSDataSet(object[,names(unlist(con))],
                        normMethod="user",
                        sizeFactor=rep(1,length(unlist(con))),
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                },
                nbp = { # Has been normalized with NBPSeq and main method was
                    # "nbsmyth"; 
                    # "nbsmyth" refers to prepare.nbp 
                    # (see norm.R lines 370-371) 
                    abs <- ABSDataSet(round(
                        object$pseudo.counts[,names(unlist(con))]),
                            normMethod="user", 
                            sizeFactor=rep(1,length(unlist(con))),
                            paired=statArgs$paired,
                            minDispersion=statArgs$minDispersion, 
                            minRates=statArgs$minRates,
                            maxRates=statArgs$maxRates, 
                            LevelstoNormFC=statArgs$LevelstoNormFC) 
                            # I take the pseudocounts of the nbp object because
                            # they are the normalized ones for the lib.size
                },
                SeqCountSet = {
                    # isolate normalization factors so as to give them names and
                    # be able to subset them at will
                    nF <- normalizationFactor(object) 
                    names(nF) <- names(unlist(sampleList))
                    abs <- ABSDataSet(round(assayData(
                        object)$exprs[,names(unlist(con))]),
                        normMethod="user", 
                        sizeFactor=nF[names(unlist(con))],
                        paired=statArgs$paired,
                        minDispersion=statArgs$minDispersion, 
                        minRates=statArgs$minRates,
                        maxRates=statArgs$maxRates, 
                        LevelstoNormFC=statArgs$LevelstoNormFC)
                }
            )   
            abs <- ABSSeqlm(abs,design=design,condA=colnames(design),
                lmodel=statArgs$lmodel,preval=statArgs$preval,
                qforkappa=statArgs$qforkappa,adjmethod=statArgs$adjmethod,
                scale=statArgs$scale,quiet=statArgs$quiet) 
                # It performs ANOVA because I have passed the design to condA 
                # and condB is NULL (default;I don't let the user to change it).
                # normalFactors() and aFoldcomplexDesign() are called in row by 
                # this function
            p[[conName]] <- abs[,"pvalue"]
            names(p[[conName]]) <- rownames(abs)
        }           
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)       
}

statDss <- function(object,sampleList,contrastList=NULL,statArgs=NULL) {
    if (is.null(statArgs))
        statArgs <- getDefaults("statistics","dss")
    if (is.null(contrastList))
        contrastList <- makeContrastList(paste(names(sampleList)[c(1,2)],
            collapse="_vs_"),sampleList)
    if (!is.list(contrastList))
        contrastList <- makeContrastList(contrastList,sampleList)
    
    classes <- asClassVector(sampleList)
    design <- as.data.frame(classes)
    # if colname is not designs then it throughs an error
    colnames(design) <- "designs" 
    p <- vector("list",length(contrastList))
    names(p) <- names(contrastList)
    
    switch(class(object)[1],
        .CountDataSet = { # Has been normalized with DESeq
            #countData <- round(DESeq::counts(object, normalized=TRUE))
            countData <- round(counts(object, normalized=TRUE))
            seqData <- newSeqCountSet(countData,design,
            normalizationFactor=rep(1,ncol(countData)))
            # estimating dispersions
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        DESeqDataSet = { # Has been normalized with DESeq2
            countData <- round(DESeq2::counts(object,normalized=TRUE))
            seqData <- newSeqCountSet(countData,design,
            normalizationFactor=rep(1,ncol(countData)))
            # estimating dispersions
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        DGEList = { # Has been normalized with edgeR              
            seqData <- newSeqCountSet(object$counts,design,
                normalizationFactor=object$samples$norm.factors)
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        SeqExpressionSet = { # Has been normalized with EDASeq
            seqData <- newSeqCountSet(normCounts(object),design,
                normalizationFactor=rep(1,ncol(normCounts))) 
            # Getting the normalized counts and passing in normFactors the 
            # value 1
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            seqData <- newSeqCountSet(object,design,
            normalizationFactor=rep(1,ncol(object)))
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        nbData ={ # Has been normalized with NBPSeq and method was "nbpseq".
            #seqData <- newSeqCountSet(round(sweep(object$counts,2,
            #   object$norm.factors,"*")),design,
            #   normalizationFactor=rep(1,ncol(object$counts)))
            seqData <- newSeqCountSet(normCounts(object),design,
                normalizationFactor=rep(1,ncol(object)))
            # Getting the normalized counts and passing in normFactors 
            # the value 1
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        }, 
        nbp = { # Has been normalized with NBPSeq and method was "nbsmyth".
            seqData <- newSeqCountSet(round(object$pseudo.counts),design,
                normalizationFactor=rep(1,ncol(object$pseudo.counts)))
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        ABSDataSet = { # Has been normalized with ABSSeq
            seqData <- newSeqCountSet(round(excounts(object)),design,
                normalizationFactor=rep(1,ncol(object$pseudo.counts)))
            seqData <- estDispersion(seqData,trend=statArgs$trend)
        },
        SeqCountSet = {
            seqData <- object
            seqData <- estDispersion(seqData,trend=statArgs$trend)       
        }
    )

    for (conName in names(contrastList)) {
        disp(" Contrast: ",conName)
        # It assignes all element values of contrastList to con (number of 
        # conditions)
        con <- contrastList[[conName]]
        # It unlists con and picks the unique values
        cons <- unique(unlist(con)) 

        if(length(con)==2) { # DE analysis of 2 conditions
            # We choose from all levels the two that will be used for the 
            # pairwise comparison: "cons[1],cons[2]"
            res <- waldTest(seqData,cons[1],cons[2],
                equal.var=statArgs$equal.var)
            # Order according to initial object! waldTest returns results sorted
            # in decreasing p-values...
            res <- res[rownames(object),,drop=FALSE]
            p[[conName]] <- res$pval
        }
        else { # DE analysis of >2 levels in conditions factor
            warnwrap(paste("DSS differential expression algorithm does not ",
                "support multi-level designs (with more than two levels in a ",
                "factor to be compared)! Switching to DESeq. ",
                "Comparison: ",conName))
            statArgs <- getDefaults("statistics","deseq")
            colData <- DataFrame(design)
            colnames(colData) <- "conditions"
            designTmp <- as.formula(c("~",names(colData[1])))
            
            # The same chunck with the norm.R script to get a cds file
            theDesign <- data.frame(condition=classes,
                row.names=colnames(seqData)) 
            cds <- newCountDataSet(as.matrix(round(assayData(seqData)$exprs)),
                conditions=theDesign$condition)
            # retrieve DSS-calculated normalizationFactors
            #DESeq::sizeFactors(cds) <- normalizationFactor(seqData)     
            sizeFactors(cds) <- normalizationFactor(seqData)     
            #cds <- DESeq::estimateDispersions(cds,method=statArgs$method,
            #        sharingMode=statArgs$sharingMode)
            cds <- estimateDispersions(cds,method=statArgs$method,
                    sharingMode=statArgs$sharingMode)

            cc <- names(unlist(con))
            cdsTmp <- cds[,cc]
            fit0 <- fitNbinomGLMs(cdsTmp,count~1)
            fit1 <- fitNbinomGLMs(cdsTmp,count~condition)
            p[[conName]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[conName]]) <- rownames(object)
        p[[conName]][which(is.na(p[[conName]]))] <- 1
    }
    return(p)
}                  
