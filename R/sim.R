#' Estimate AUFC weights
#'
#' This function automatically estimates weights for the \code{"weight"} and
#' \code{"dperm_weight"} options of metaseqR for combining p-values from multiple
#' statistical tests. It creates simulated dataset based on real data and then
#' performs statistical analysis with metaseqR several times in order to derive
#' False Discovery Curves. Then, the average areas under the false discovery curves
#' are used to construct weights for each algorithm, according to its performance
#' when using simulated data.
#'
#' @param counts the real raw counts table from which the simulation parameters
#' will be estimated. It must not be normalized and must contain only integer
#' counts, without any other annotation elements and unique gene identifiers as
#' the rownames attribute.
#' @param normalization same as \code{normalization} in \code{link{metaseqr}}.
#' @param statistics same as \code{statistics} in \code{link{metaseqr}}.
#' @param nsim the number of simulations to perform to estimate the weights. It
#' default to 10.
#' @param N the number of genes to produce. See \code{link{makeSimDataSd}}.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fcBasis the minimum fold-change for deregulation.
#' @param top the top \code{top} best ranked (according to p-value) to use, to
#' calculate area under the false discovery curve.
#' @param modelOrg the organism from which the data are derived. It must be one
#' of \code{\link{metaseqr}} supported organisms.
#' @param seed a number to be used as seed for reproducible simulations. Defaults
#' to \code{NULL} (NOT reproducible results!).
#' @param drawFpc draw the averaged false discovery curves? Default to \code{FALSE}.
#' @param multic whether to run in parallel (if package \code{parallel} is present
#' or not.
#' @param ... Further arguments to be passed to \code{\link{estimateSimParams}}.
#' @value A vector of weights to be used in \link{\code{metaseqr}} with the
#' \code{weights} option.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data("mm9.geneData",package="metaseqR")
#' multic <- checkParallel(0.8)
#' weights <- estimateAufcWeights(
#'   counts=as.matrix(mm9.geneCounts[,9:12]),
#'   normalization="edaseq",
#'   statistics=c("deseq","edger"),
#'   nsim=3,N=100,ndeg=c(10,10),top=10,modelOrg="mm9",
#'   seed=100,multic=multic,libsizeGt=1e+5
#' )
#'}
estimateAufcWeights <- function(counts,normalization,statistics,nsim=10,
    N=10000,samples=c(3,3),ndeg=c(500,500),top=500,modelOrg="mm9",fcBasis=1.5,
    seed=NULL,drawFpc=FALSE,rc=NULL,...) {
    if (!require(zoo))
        stopwrap("R pacakage zoo is required in order to estimate AUFC ",
            "weights!")

    if (is.null(seed)) {
        seedStart <- round(100*runif(1))
        seedEnd <- seedStart + nsim - 1
        seed <- as.list(seedStart:seedEnd)
    }
    else {
        set.seed(seed)
        seedStart <- round(100*runif(1))
        seedEnd <- seedStart + nsim - 1
        seed <- as.list(seedStart:seedEnd)
    }
    
    if (ncol(counts)<4)
        stopwrap("Cannot estimate AUFC weights with an initial dataset with ",
            "less than 4 samples!")
    else if (ncol(counts)>=4 && ncol(counts)<10) {
        set.seed(seedStart)
        reind <- sample(1:ncol(counts),20,replace=TRUE)
        counts <- counts[,reind]
    }
    parList <- estimateSimParams(counts,...)

    disp("Running simulations... This procedure requires time... Please ",
        "wait...")
    simResults <- cmclapply(seed,function(x,normalization,statistics,N,
        parList,samples,ndeg,fcBasis,modelOrg) {
        D <- makeSimDataSd(N=N,param=parList,samples=samples,ndeg=ndeg,
            fcBasis=fcBasis,modelOrg=modelOrg,seed=x)
        dd <- D$simdata
        
        if (!is.null(modelOrg)) {
            tmp <- metaseqr(
                counts=dd,
                sampleList=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                idCol=4,
                gcCol=5,
                nameCol=7,
                btCol=8,
                org=modelOrg,
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
        else {
            tmp <- metaseqr(
                counts=dd,
                sampleList=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                idCol=4,
                gcCol=5,
                nameCol=7,
                btCol=8,
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

#' Create simulated counts using TCC package
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' \code{simulateReadCounts} function from the Bioconductor
#' package TCC and it adds simulated annoation elements. For further information
#' please consult the TCC package documentation. Note that the produced data are
#' based in an Arabidopsis dataset.
#'
#' @param ... parameters to the \code{simulateReadCounts} function.
#' @return A list with the following members: \code{simdata} holding the simulated
#' dataset complying with metaseqr requirements, and \code{simparam} holding the
#' simulation parameters (see TCC documentation).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' dd <- make.simData(Ngene=10000,PDEG=0.2,DEG.assign=c(0.9,0.1),
#'   DEG.foldchange=c(5,5),replicates=c(3,3))
#' head(dd$simdata)
#'}
makeSimDataTcc <- function(...) {
    if (suppressWarnings(!require(TCC)))
        stopwrap("Bioconductor package TCC is required to create simulated data!")
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

#' Create simulated counts using the Soneson-Delorenzi method
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' method presented in (Soneson and Delorenzi, BMC Bioinformatics, 2013). For the
#' time being, it creates only simulated datasets with two conditions.
#'
#' @param N the number of genes to produce.
#' @param param a named list with negative binomial parameter sets to sample from.
#' The first member is the mean parameter to sample from (\code{muHat}} and the
#' second the dispersion (\code{phiHat}). This list can be created with the
#' \code{\link{estimateSimParams}} function.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fcBasis the minimum fold-change for deregulation.
#' @param libsizeRange a vector with 2 numbers (generally small, see the default),
#' as they are multiplied with \code{libsizeMag}.
#' These numbers control the library sized of the synthetic data to be produced.
#' @param libsizeMag a (big) number to multiply the \code{libsizeRange} to
#' produce library sizes.
#' @param modelOrg the organism from which the real data are derived from. It
#' must be one of the supported organisms (see the main \code{\link{metaseqr}}
#' help page). It is used to sample real values for GC content.
#' @param simLengthBias a boolean to instruct the simulator to create genes
#' whose read counts is proportional to their length. This is achieved by sorting
#' in increasing order the mean parameter of the negative binomial distribution
#' (and the dispersion according to the mean) which will cause an increasing gene
#' count length with the sampling. The sampled lengths are also sorted so that in
#' the final gene list, shorter genes have less counts as compared to the longer
#' ones. The default is FALSE.
#' @param seed a seed to use with random number generation for reproducibility.
#' @return A named list with two members. The first member (\code{simdata})
#' contains the synthetic dataset 
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # File "bottomly_read_counts.txt" from the ReCount database
#' download.file("http://bowtie-bio.sourceforge.net/recount/countTables/bottomly_count_table.txt",
#'   destfile="~/bottomly_count_table.txt")
#' N <- 10000
#' parList <- estimateSimParams("~/bottomly_read_counts.txt")
#' sim <- makeSimDataSd(N,parList)
#' synth.data <- sim$simdata
#' true.deg <- which(sim$truedeg!=0)
#'}
makeSimDataSd <- function(N,param,samples=c(5,5),ndeg=rep(round(0.1*N),2),
    fcBasis=1.5,libsizeRange=c(0.7,1.4),libsizeMag=1e+7,modelOrg=NULL,
    simLengthBias=FALSE,seed=NULL) {
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

    if (!is.null(seed)) 
		set.seed(seed)
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
    if (!is.null(seed))
		set.seed(seed)
    L1 <- round(libsizeMag*runif(s1,min=libsizeRange[1],
        max=libsizeRange[2]))
    if (!is.null(seed)) set.seed(2*seed)
    L2 <- round(libsizeMag*runif(s2,min=libsizeRange[1],
        max=libsizeRange[2]))

    lambda1 <- do.call("cbind",rep(list(muHat[ii]),s1))
    mu1 <- sweep(lambda1,2,L1/sum(lambda1[,1]),"*")
    sim1 <- matrix(0,N,s1)
    for (j in 1:s1) {
        if (!is.null(seed)) set.seed(seed+j)
        sim1[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu1[,j])
    }

    v <- numeric(N)
    if (sum(ndeg)>0) {
        if (!is.null(seed)) set.seed(seed)
        iUpdown <- sample(1:length(v),sum(ndeg))
        regDir <- rep(c(1,-1),c(ndeg[1],ndeg[2]))
        v[iUpdown] <- regDir
        if (!is.null(seed)) set.seed(seed+19051980)
        lambda2 <- ((fcBasis + rexp(N))^v)*lambda1
        mu2 <- sweep(lambda2,2,L2/sum(lambda2[,1]),"*")
        sim2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim2[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu2[,j])
    }
    else {
        if (!is.null(seed)) set.seed(seed+19051980)
        lambda2 <- lambda1
        mu2 <- sweep(lambda2,2,L2/sum(lambda2[,1]),"*")
        sim2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim2[,j] <- rnbinom(N,size=1/phiHat[ii],mu=mu2[,j])
    }

    # Now we have to simulate annotation
    if (!is.null(seed))
		set.seed(seed)
    chromosome <- paste("chr",1+round(20*runif(N)),sep="")
    gene_id <- gene_name <- paste("gene",1:N,sep="_")
    if (!is.null(modelOrg)) {
        if (!is.null(seed)) set.seed(seed)
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
        if (!is.null(seed))
			set.seed(seed)
        gc_content <- runif(N)
        if (!is.null(seed))
			set.seed(seed)
        start <- 1 + round(1e+6*runif(N))
        if (!is.null(seed))
			set.seed(seed)
        end <- start + 250 + round(1e+6*runif(N))
        if (!is.null(seed))
			set.seed(seed)
        strand <- sample(c("+","-"),N,replace=TRUE)
        if (simLengthBias) {
            lenix <- sort(end-start,index.return=TRUE)$ix
            start <- start[lenix]
            end <- end[lenix]
            gc_content <- gc_content[lenix]
            strand <- strand[lenix]
        }
    }
    if (!is.null(seed))
		set.seed(seed)
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

#' Estimate negative binomial parameters from real data
#'
#' This function reads a read counts table containing real RNA-Seq data (preferebly
#' with more than 20 samples so as to get as much accurate as possible estimations)
#' and calculates a population of count means and dispersion parameters which can
#' be used to simulate an RNA-Seq dataset with synthetic genes by drawing from a
#' negative binomial distribution. This function works in the same way as described
#' in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al., BMC
#' Genomics, 2012).
#'
#' @param realCounts a text tab-delimited file with real RNA-Seq data. The file
#' should strictly contain a unique gene name (e.g. Ensembl accession) in the
#' first column and all other columns should contain read counts for each gene.
#' Each column must be named with a unique sample identifier. See examples in the
#' ReCount database \link{http://bowtie-bio.sourceforge.net/recount/}.
#' @param libsizeGt a library size below which samples are excluded from parameter
#' estimation (default: 3000000).
#' @param rowmeansGt a row means (mean counts over samples for each gene) below
#' which genes are excluded from parameter estimation (default: 5).
#' @param eps the tolerance for the convergence of \code{\link{optimize}} function.
#' Defaults to 1e-11.
#' @param restrict.cores in case of parallel optimization, the fraction of the
#' available cores to use.
#' @param seed a seed to use with random number generation for reproducibility.
#' @param draw boolean to determine whether to plot the estimated simulation
#' parameters (mean and dispersion) or not. Defaults to \code{FALSE} (do not draw
#' a mean-dispersion scatterplot).
#' @return A named list with two members: \code{muHat} which contains negative
#' binomial mean estimates and \code{phiHat} which contains dispersion.
#' estimates
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # Dowload locally the file "bottomly_count_table.txt" from the ReCount datbase
#' download.file("http://bowtie-bio.sourceforge.net/recount/countTables/bottomly_count_table.txt",
#'   destfile="~/bottomly_count_table.txt")
#' # Estimate simulation parameters
#' parList <- estimateSimParams("~/bottomly_count_table.txt")
#'}
estimateSimParams <- function(realCounts,libsizeGt=3e+6,rowmeansGt=5,
    eps=1e-11,rc=NULL,seed=42,draw=FALSE) {
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
    dmat <- downsampleCounts(mat,seed)
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

#' Downsample read counts
#'
#' This function downsamples the library sizes of a read counts table to the lowest
#' library size, according to the methdology used in  (Soneson and Delorenzi,
#' BMC Bioinformatics, 2013).
#'
#' @param counts the read counts table which is subjected to downsampling.
#' @param seed random seed for reproducible downsampling.
#' @return The downsampled counts matrix.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # Dowload locally the file "bottomly_read_counts.txt" from
#' # the ReCount database
#' download.file(paste("http://bowtie-bio.sourceforge.net/",
#'   "recount/countTables/bottomly_count_table.txt",sep=""),
#'  destfile="~/bottomly_count_table.txt")
#' M <- as.matrix(read.delim("~/bottomly_count_table.txt",row.names=1))
#' D <- downsampleCounts(M)
#'}
downsampleCounts <- function(counts,seed=42) {
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

#' MLE dispersion estimate
#'
#' MLE function used to estimate negative binomial dispersions from real RNA-Seq
#' data, as in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al.,
#' BMC Genomics, 2012). Internal use.
#'
#' @param phi the parameter to be optimized.
#' @param y count samples used to perform the optimization.
#' @return objective function value.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
mlfo <- function(phi,y) {
    N <- length(y)
    mu <- mean(y)
    -(sum(lgamma(y+1/phi)) - N*lgamma(1/phi) - sum(lgamma(y+1)) + 
        sum(y*log(mu*phi/(1+mu*phi))) - (N/phi)*log(1+mu*phi))
}

#' Create counts matrix permutations
#'
#' This function creates a permuted read counts matrix based on the \code{contrast}
#' argument (to define new virtual contrasts of the same number) and on the
#' \code{sampleList} to derive the number of samples for each virtual condition.
#' It is a helper for the \code{\link{metaPerm}} function.
#'
#' @param counts the gene read counts matrix.
#' @param sampleList the list containing condition names and the samples under
#' each condition.
#' @param contrast the contrasts vector. See the main \code{\link{metaseqr}} help
#' page.
#' @param repl the same as the replace argument in \code{\link{sample}} function.
#' @return A list with three members: the matrix of permuted per sample read counts,
#' the virtual sample list and the virtual contrast to be used with the \code{stat.*}
#' functions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data("mm9.geneData",package="metaseqR")
#' per <- makePermutation(mm9.geneCounts,sampleList.mm9,"e14.5_vs_adult_8_weeks")
#'}
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

#' Calculate the ratio TP/(FP+FN)
#'
#' This function calculates the ratio of True Positives to the sum of False
#' Positives and False Negatives given a matrix of p-values (one for each
#' statistical test used) and a vector of ground truth (DE or non-DE). This
#' function serves as a method evaluation helper.
#'
#' @param truth the ground truth differential expression vector. It should contain
#' only zero and non-zero elements, with zero denoting non-differentially expressed
#' genes and non-zero, differentially expressed genes. Such a vector can be obtained
#' for example by using the \code{\link{makeSimDataSd}} function, which creates
#' simulated RNA-Seq read counts based on real data. It MUST be named with gene
#' names, the same as in \code{p}.
#' @param p a p-value matrix whose rows correspond to each element in the
#' \code{truth} vector. If the matrix has a \code{colnames} attribute, a legend
#' will be added to the plot using these names, else a set of column names will
#' be auto-generated. \code{p} can also be a list or a data frame. In any case,
#' each row (or element) MUST be named with gene names (the same as in \code{truth}).
#' @param sig a significance level (0 < \code{sig} <=1).
#' @return A named list with two members. The first member is a data frame with
#' the numbers used to calculate the TP/(FP+FN) ratio and the second member is
#' the ratio TP/(FP+FN) for each statistical test.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p1 <- 0.001*matrix(runif(300),100,3)
#' p2 <- matrix(runif(300),100,3)
#' p <- rbind(p1,p2)
#' rownames(p) <- paste("gene",1:200,sep="_")
#' colnames(p) <- paste("method",1:3,sep="_")
#' truth <- c(rep(1,40),rep(-1,40),rep(0,10),rep(1,10),rep(2,10),rep(0,80))
#' names(truth) <- rownames(p)
#' otr <- calcOtr(truth,p)
#'}
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

#' Calculate the F1-score
#'
#' This function calculates the F1 score (2*(precision*recall/precision+racall)
#' or 2*TP/(2*TP+FP+FN) given a matrix of p-values (one for each statistical test
#' used) and a vector of ground truth (DE or non-DE). This function serves as a
#' method evaluation helper.
#'
#' @param truth the ground truth differential expression vector. It should contain
#' only zero and non-zero elements, with zero denoting non-differentially expressed
#' genes and non-zero, differentially expressed genes. Such a vector can be obtained
#' for example by using the \code{\link{makeSimDataSd}} function, which creates
#' simulated RNA-Seq read counts based on real data. It MUST be named with gene
#' names, the same as in \code{p}.
#' @param p a p-value matrix whose rows correspond to each element in the
#' \code{truth} vector. If the matrix has a \code{colnames} attribute, a legend
#' will be added to the plot using these names, else a set of column names will
#' be auto-generated. \code{p} can also be a list or a data frame. In any case,
#' each row (or element) MUST be named with gene names (the same as in \code{truth}).
#' @param sig a significance level (0 < \code{sig} <=1).
#' @return A named list with two members. The first member is a data frame with
#' the numbers used to calculate the F1-score and the second member is the
#' F1-score for each statistical test.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p1 <- 0.001*matrix(runif(300),100,3)
#' p2 <- matrix(runif(300),100,3)
#' p <- rbind(p1,p2)
#' rownames(p) <- paste("gene",1:200,sep="_")
#' colnames(p) <- paste("method",1:3,sep="_")
#' truth <- c(rep(1,40),rep(-1,40),rep(0,10),rep(1,10),rep(2,10),rep(0,80))
#' names(truth) <- rownames(p)
#' f1 <- calcF1Score(truth,p)
#'}
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
