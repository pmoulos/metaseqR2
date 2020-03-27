estimateAufcWeights <- function(counts,normalization,statistics,nsim=10,
    N=10000,samples=c(3,3),ndeg=c(500,500),top=500,modelOrg="mm9",fcBasis=1.5,
    drawFpc=FALSE,rc=NULL,...) {
    if (!requireNamespace("zoo"))
        stopwrap("R pacakage zoo is required in order to estimate AUFC ",
            "weights!")

    if (ncol(counts)<4)
        stopwrap("Cannot estimate AUFC weights with an initial dataset with ",
            "less than 4 samples!")
    else if (ncol(counts)>=4 && ncol(counts)<10) {
        reind <- sample(1:ncol(counts),20,replace=TRUE)
        counts <- counts[,reind]
    }
    parList <- estimateSimParams(counts,...)

    disp("Running simulations... This procedure requires time... Please ",
        "wait...")
    simResults <- cmclapply(1:nsim,function(x,normalization,statistics,N,
        parList,samples,ndeg,fcBasis,modelOrg) {
        D <- makeSimDataSd(N=N,param=parList,samples=samples,ndeg=ndeg,
            fcBasis=fcBasis,modelOrg=modelOrg)
        dd <- D$simdata
        
        if (!is.null(modelOrg)) {
            tmp <- metaseqr2(
                counts=dd,
                sampleList=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                embedCols=list(
                    idCol=4,
                    gcCol=5,
                    nameCol=7,
                    btCol=8
                ),
                org=modelOrg,
                countType="gene",
                normalization=normalization,
                statistics=statistics,
                metaP="simes",
                figFormat="png",
                preset="all_basic",
                exportWhere=tempdir(),
                restrictCores=rc,
                qcPlots=NULL,
                report=FALSE,
                runLog=FALSE,
                outList=TRUE
            )
        }
        else {
            tmp <- metaseqr2(
                counts=dd,
                sampleList=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                embedCols=list(
                    idCol=4,
                    gcCol=5,
                    nameCol=7,
                    btCol=8
                ),
                countType="gene",
                normalization=normalization,
                statistics=statistics,
                metaP="simes",
                figFormat="png",
                preset="all_basic",
                exportWhere=tempdir(),
                qcPlots=NULL,
                report=FALSE,
                runLog=FALSE,
                outList=TRUE
            )
        }

        # Retrieve several p-values
        pList <- vector("list",length(statistics))
        for (s in statistics) {
            field <- paste("p-value",s,sep="_")
            pList[[s]] <- tmp$data[[1]][,field]
            names(pList[[s]]) <- rownames(tmp$data[[1]])
        }
        pMatrix <- do.call("cbind",pList)
        return(list(simdata=D,pvalues=pMatrix))
    },normalization,statistics,N,parList,samples,ndeg,fcBasis,modelOrg,rc=rc)

    disp("Estimating AUFC weights... Please wait...")
    fpcObj <- cmclapply(simResults,function(x) {
        trueDe <- x$simdata$truedeg
        names(trueDe) <- rownames(x$simdata$simdata)
        pMatrix <- x$pvalues
        trueDe <- trueDe[rownames(pMatrix)]
        fdc <- diagplotFtd(trueDe,pMatrix,type="fpc",draw=FALSE)
    },rc=rc)
    avgFpc <- diagplotAvgFtd(fpcObj,draw=drawFpc)

    x <- 1:top
    aufc <- apply(avgFpc$avgFtdr$means[1:top,],2,function(x,i) {
        return(sum(diff(i)*rollmean(x,2)))
    },x)
    weightAufc <- (sum(aufc)/aufc)/sum(sum(aufc)/aufc)
    return(weightAufc)
}

makeSimDataTcc <- function(...) {
    if (suppressWarnings(!requireNamespace("TCC")))
        stopwrap("Bioconductor package TCC is required to create ",
            "simulated data!")
    #tcc <- simulateReadCounts(Ngene=Ngene,PDEG=PDEG,DEG.assign=DEG.assign,
    #    DEG.foldchange=DEG.foldchange,replicates=replicates)
    tcc <- TCC::simulateReadCounts(...)
    n <- nrow(tcc$count)
    # Now we have to simulate annotation
    chromosome <- paste("chr",1+round(20*runif(n)),sep="")
    start <- 1 + round(1e+6*runif(n))
    end <- start + 250 + round(1e+6*runif(n))
    gene_id <- gene_name <- rownames(tcc$count)
    gc_content <- runif(n)
    strand <- sample(c("+","-"),n,replace=TRUE)
    biotype <- sample(paste("biotype",1:10),n,replace=TRUE)
    simData <- data.frame(
        chromosome=chromosome,
        start=start,
        end=end,
        gene_id=gene_id,
        gc_content=gc_content,
        strand=strand,
        gene_name=gene_name,
        biotype=biotype
    )
    simData <- cbind(simData,tcc$count)
    return(list(simdata=simData,simparam=tcc$simulation))
}

makeSimDataSd <- function(N,param,samples=c(5,5),ndeg=rep(round(0.1*N),2),
    fcBasis=1.5,libsizeRange=c(0.7,1.4),libsizeMag=1e+7,modelOrg=NULL,
    simLengthBias=FALSE) {
    if (!is.null(modelOrg)) {
        modelOrg <- tolower(modelOrg)
        checkTextArgs("modelOrg",modelOrg,c("hg18","hg19","mm9","mm10",
            "rno5","dm3","rn5","rn6","danrer7","pantro4","tair10"),
            multiarg=FALSE)
        ann <- getAnnotation(modelOrg,"gene")
        realGc <- as.numeric(ann$gc_content)
        realStart <- as.numeric(ann$start)
        realEnd <- as.numeric(ann$end)
        realStrand <- as.character(ann$strand)
    }
    muHat <- param$muHat
    phiHat <- param$phiHat

    if (simLengthBias) {
        sind <- sort(muHat,index.return=TRUE)$ix
        muHat <- muHat[sind]
        phiHat <- phiHat[sind]
        if (length(muHat)>=N)
            ii <- sort(sample(1:length(muHat),N))
        else
            ii <- sort(sample(1:length(muHat),N,replace=TRUE))
    }
    else {
        if (length(muHat)>=N)
            ii <- sample(1:length(muHat),N)
        else
            ii <- sample(1:length(muHat),N,replace=TRUE)
    }

    s1 <- samples[1]
    s2 <- samples[2]
    L1 <- round(libsizeMag*runif(s1,min=libsizeRange[1],
        max=libsizeRange[2]))
    L2 <- round(libsizeMag*runif(s2,min=libsizeRange[1],
        max=libsizeRange[2]))

    lambda1 <- do.call("cbind",rep(list(muHat[ii]),s1))
    mu1 <- sweep(lambda1,2,L1/sum(lambda1[,1]),"*")
    sim1 <- matrix(0,N,s1)
    for (j in 1:s1)
        sim1[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu1[,j])

    v <- numeric(N)
    if (sum(ndeg)>0) {
        iUpdown <- sample(1:length(v),sum(ndeg))
        regDir <- rep(c(1,-1),c(ndeg[1],ndeg[2]))
        v[iUpdown] <- regDir
        lambda2 <- ((fcBasis + rexp(N))^v)*lambda1
        mu2 <- sweep(lambda2,2,L2/sum(lambda2[,1]),"*")
        sim2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim2[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu2[,j])
    }
    else {
        lambda2 <- lambda1
        mu2 <- sweep(lambda2,2,L2/sum(lambda2[,1]),"*")
        sim2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim2[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu2[,j])
    }

    # Now we have to simulate annotation
    chromosome <- paste("chr",1+round(20*runif(N)),sep="")
    gene_id <- gene_name <- paste("gene",1:N,sep="_")
    if (!is.null(modelOrg)) {
        if (length(realGc)>=N)
            sampleInd <- sample(1:length(realGc),N)            
        else
            sampleInd <- sample(1:length(realGc),N,replace=TRUE)
        gc_content <- realGc[sampleInd]
        start <- realStart[sampleInd]
        end <- realEnd[sampleInd]
        strand <- realStrand[sampleInd]
        if (simLengthBias) {
            lenix <- sort(end-start,index.return=TRUE)$ix
            start <- start[lenix]
            end <- end[lenix]
            gc_content <- gc_content[lenix]
            strand <- strand[lenix]
        }
    }
    else {
        gc_content <- runif(N)
        start <- 1 + round(1e+6*runif(N))
        end <- start + 250 + round(1e+6*runif(N))
        strand <- sample(c("+","-"),N,replace=TRUE)
        if (simLengthBias) {
            lenix <- sort(end-start,index.return=TRUE)$ix
            start <- start[lenix]
            end <- end[lenix]
            gc_content <- gc_content[lenix]
            strand <- strand[lenix]
        }
    }
    biotype <- sample(paste("biotype",1:10),N,replace=TRUE)
    simData <- data.frame(
        chromosome=chromosome,
        start=start,
        end=end,
        gene_id=gene_id,
        gc_content=gc_content,
        strand=strand,
        gene_name=gene_name,
        biotype=biotype
    )
    colnames(sim1) <- paste("G1_rep",1:s1,sep="")
    colnames(sim2) <- paste("G2_rep",1:s2,sep="")
    rownames(sim1) <- rownames(sim2) <- names(v) <- gene_id

    return(list(simdata=cbind(simData,sim1,sim2),truedeg=v))
}

estimateSimParams <- function(realCounts,libsizeGt=3e+6,rowmeansGt=5,
    eps=1e-11,rc=NULL,draw=FALSE) {
    if (is.data.frame(realCounts))
        mat <- as.matrix(realCounts)
    else if (is.matrix(realCounts))
        mat <- realCounts
    else if (file.exists(realCounts)) {
        realData <- read.delim(realCounts,row.names=1)
        mat <- as.matrix(realData)
    }
    else
        stopwrap("The input count data must be either a file, a matrix or a ",
            "data frame!")
    
    lowLib <- which(apply(mat,2,sum)<libsizeGt)
    if (length(lowLib)==ncol(mat))
        stopwrap("Cannot estimate simulation parameters as the library sizes ",
            "are too small! Try lowering the value of the libsizeGt ",
            "parameter...")
    if (length(lowLib)>0)
        mat <- mat[,-lowLib]
    disp("Downsampling counts...")
    dmat <- downsampleCounts(mat)
    lowCo <- which(apply(dmat,1,
        function(x) if (mean(x)<5) TRUE else FALSE))
    if (length(lowCo)>0)
        dmat <- dmat[-lowCo,]
    muHat <- apply(dmat,1,mean)
    disp("Estimating initial dispersion population...")
    phiEst <- apply(dmat,1,function(x) {
        m <- mean(x)
        v <- var(x)
        phi <- (v-m)/m^2
        return(phi)
    })
    phiInd <- which(phiEst>0)
    phiEst <- phiEst[phiInd]
    dmat <- dmat[phiInd,]
    disp("Estimating dispersions using log-likelihood...\n")
    init <- cmclapply(seq_along(1:nrow(dmat)),function(i,d,p) {
        list(y=d[i,],h=p[i])
    },dmat,phiEst,rc=rc)
    phiHat <- unlist(cmclapply(init,function(x,eps) {
        optimize(mlfo,c(x$h-1e-2,x$h+1e-2),y=x$y,tol=eps)$minimum
    },eps,rc=rc))
    if (draw) {
        dev.new()
        plot(log10(muHat[phiInd]),log10(phiHat),col="blue",pch=20,cex=0.5,
            xlab="",ylab="")
        title(xlab="mean",ylab="dispesion",font=2,cex=0.9)
        grid()
    }
    return(list(muHat=muHat[phiInd],phiHat=phiHat))
}

downsampleCounts <- function(counts) {
    libSizes <- apply(counts,2,sum)
    targetSize <- min(libSizes)
    toRemove <- libSizes-targetSize
    ii <- which(toRemove>0)
    dcounts <- counts
    for (i in ii) {
        tmp <- round(toRemove[i]*(counts[,i]/sum(counts[,i])))
        victimSize <- sum(tmp)
        if (victimSize>toRemove[i]) {
            dif <- victimSize - toRemove[i]
            #victims <- sample(1:length(tmp),dif)
            victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
            tmp[victims] <- tmp[victims] - 1
        }
        else if (victimSize<toRemove[i]) {
            dif <- toRemove[i] - victimSize
            #victims <- sample(1:length(tmp),dif)
            victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
            tmp[victims] <- tmp[victims] + 1
        }
        dcounts[,i] <- dcounts[,i] - tmp
    }
    return(dcounts)
}

mlfo <- function(phi,y) {
    N <- length(y)
    mu <- mean(y)
    -(sum(lgamma(y+1/phi)) - N*lgamma(1/phi) - sum(lgamma(y+1)) + 
        sum(y*log(mu*phi/(1+mu*phi))) - (N/phi)*log(1+mu*phi))
}

makePermutation <- function(counts,sampleList,contrast,repl=FALSE) {
    cnts <- strsplit(contrast,"_vs_")[[1]]
    virtualContrast <- paste(paste("VirtCond",1:length(cnts),sep=""),
        collapse="_vs_")
    vitualSampleList <- vector("list",length(sampleList))
    names(vitualSampleList) <- paste("VirtCond",1:length(sampleList),sep="")
    # Avoid the extreme case of returning a vector with all samples the same
    if (repl) {
        resample <- rep(1,ncol(counts))
        while(length(unique(resample))==1)
            resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
    }
    else
        resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
    virtualCounts <- counts[,resample]
    samples <- paste("VirtSamp",1:ncol(counts),sep="")
    colnames(virtualCounts) <- samples
    nsample <- sapply(sampleList,length)
    virtualSamples <- split(samples,rep(1:length(nsample),nsample))
    names(virtualSamples) <- names(vitualSampleList)
    for (n in names(vitualSampleList))
        vitualSampleList[[n]] <- virtualSamples[[n]]
    return(list(counts=virtualCounts,sampleList=vitualSampleList,
        contrast=virtualContrast))
}

calcOtr <- function(truth,p,sig=0.05) {
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    sigGenes <- trueIsects <- missed <- vector("list",ncol(pmat))
    names(sigGenes) <- names(trueIsects) <- names(missed) <- colnames(pmat)
    for (n in colnames(pmat)) {
        sigGenes[[n]] <- names(which(pmat[,n]<sig))
        trueIsects[[n]] <- intersect(sigGenes[[n]],names(which(truth!=0)))
        missed[[n]] <- setdiff(names(which(truth!=0)),trueIsects[[n]])
    }
    result <- data.frame(
        P=sapply(sigGenes,length),
        TP=sapply(trueIsects,length),
        FN=sapply(missed,length)
    )
    result$FP <- result$P - result$TP
    otr <- result$TP/(result$FP+result$FN)
    names(otr) <- rownames(result)
    return(list(result=result,otr=otr))
}

calcF1Score <- function(truth,p,sig=0.05) {
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    sigGenes <- trueIsects <- missed <- vector("list",ncol(pmat))
    names(sigGenes) <- names(trueIsects) <- names(missed) <- colnames(pmat)
    for (n in colnames(pmat)) {
        sigGenes[[n]] <- names(which(pmat[,n]<sig))
        trueIsects[[n]] <- intersect(sigGenes[[n]],names(which(truth!=0)))
        missed[[n]] <- setdiff(names(which(truth!=0)),trueIsects[[n]])
    }
    result <- data.frame(
        P=sapply(sigGenes,length),
        TP=sapply(trueIsects,length),
        FN=sapply(missed,length)
    )
    result$FP <- result$P - result$TP
    f1 <- 2*result$TP/(2*result$TP+result$FP+result$FN)
    names(f1) <- rownames(result)
    return(list(result=result,f1=f1))
}
