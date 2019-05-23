#' Meta-analysis using several RNA-Seq statistics
#'
#' This function calculates the combined p-values when multiple statistical algorithms
#' are applied to the input dataset. It is a helper and it requires very specific 
#' arguments so it should not be used individually.
#'
#' @param cpList a named list whose names are the contrasts requested from metaseqR.
#' Each member is a p-value matrix whose colnames are the names of the statistical
#' tests applied to the data. See the main \code{\link{metaseqr}} help page.
#' @param metaP the p-value combination method to use. See the main 
#' \code{\link{metaseqr}} help page.
#' @param counts the normalized and possibly filtered read counts matrix. See the
#' main \code{\link{metaseqr}} help page.
#' @param sampleList the list containing condition names and the samples under
#' each condition. See the main \code{\link{metaseqr}} help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main 
#' \code{\link{metaseqr}} help page.
#' @param statArgs the parameters for each statistical argument. See the main
#' \code{\link{metaseqr}} help page.
#' @param libsizeList a list with library sizes. See the main \code{\link{metaseqr}}
#' and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm.
#' @param reprod create reproducible permutations. Ideally one would want to create
#' the same set of indices for a given dataset so as to create reproducible p-values.
#' If \code{reprod=TRUE}, a fixed seed is used by \code{metaPerm} for all the
#' datasets analyzed with \code{metaseqr}. If \code{reprod=FALSE}, then the
#' p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod}
#' can be a numeric vector of seeds with the same length as \code{nperm} so that
#' the user can supply his/her own seeds.
#' @param multic use multiple cores to execute the premutations. This is an
#' external parameter and implies the existence of parallel package in the execution
#' environment. See the main \code{\link{metaseqr}} help page.
#' @return A named list with combined p-values. The names are the contrasts and
#' the list members are combined p-value vectors, one for each contrast.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
metaTest <- function(cpList,metaP=c("simes","bonferroni","fisher",
    "dperm_min","dperm_max","dperm_weight","fperm","whitlock","minp","maxp",
    "weight","pandora","none"),counts,sampleList,statistics,statArgs,
    libsizeList,nperm=10000,weight=rep(1/length(statistics),
    length(statistics)),reprod=TRUE,rc=NULL) {
    checkTextArgs("metaP",metaP,c("simes","bonferroni","fisher","dperm_min",
        "dperm_max","dperm_weight","fperm","whitlock","minp","maxp","weight",
        "pandora","none"))
    contrast <- names(cpList)
    disp("Performing meta-analysis with ",metaP)
    if (metaP=="pandora")
        metaP <- "weight"
    switch(metaP,
        fisher = {
            sumpList <- cmclapply(cpList,function(x) {
                tmp <- fisherMethod(x,p.corr="none",
                    zeroSub=.Machine$double.xmin)
                rp <- tmp$p.value
                names(rp) <- rownames(x)
                return(rp)
            },rc=rc)
        },
        fperm = {
            sumpList <- cmclapply(cpList,function(x) {
                if (multic)
                    tmp <- fisherMethodPerm(x,p.corr="none",B=nperm,
                        mc.cores=getOption("cores"),zeroSub=1e-32)
                else
                    tmp <- fisherMethodPerm(x,p.corr="none",B=nperm,
                        zeroSub=.Machine$double.xmin)
                return(tmp$p.value)
            },rc=rc)
        },
        whitlock = {
            sumpList <- cmclapply(cpList,function(x) 
                return(apply(x,1,combine.test,method="z.transform")),rc=rc)
        },
        simes = {
            sumpList <- cmclapply(cpList,function(x) {
                return(apply(x,1,combineSimes))
            },rc=rc)
        },
        bonferroni = {
            sumpList <- cmclapply(cpList,function(x) {
                return(apply(x,1,combineBonferroni))
            },rc=rc)
        },
        minp = {
            sumpList <- cmclapply(cpList,function(x) {
                return(apply(x,1,combineMinp))
            },rc=rc)
        },
        maxp = {
            sumpList <- cmclapply(cpList,function(x) {
                return(apply(x,1,combineMaxp))
            },rc=rc)
        },
        weight = {
            sumpList <- wapply(cpList,function(x) {
                return(apply(x,1,combineWeight,weight))
            },rc=rc)
        },
        dperm_min = {
            sumpList <- vector("list",length(cpList))
            names(sumpList) <- names(cpList)
            conl <- as.list(contrast)
            names(conl) <- contrast
            tempPList <- cmclapply(conl,metaPerm,
                counts=counts,sampleList=sampleList,
                statistics=statistics,statArgs=statArgs,
                libsizeList=libsizeList,
                nperm=nperm,weight=weight,
                select="min",reprod=reprod,
                rc=rc)
            originalPList <- cmclapply(cpList,function(x,m,w=NULL) {
                x[which(is.na(x))] <- 1
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) return(prod(p^w)),
                            w))
                    }
                )
            },"min",rc=rc)
            for (cc in names(originalPList))
            {
                pc <- cbind(tempPList[[cc]],originalPList[[cc]])
                ly <- ncol(pc)
                sumpList[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
        },
        dperm_max = {
            sumpList <- vector("list",length(cpList))
            names(sumpList) <- names(cpList)
            conl <- as.list(contrast)
            names(conl) <- contrast
            tempPList <- cmclapply(conl,metaPerm,
                counts=counts,sampleList=sampleList,
                statistics=statistics,statArgs=statArgs,
                libsizeList=libsizeList,
                nperm=nperm,weight=weight,
                select="max",reprod=reprod,
                rc=rc)
            originalPList <- cmclapply(cpList,function(x,m,w=NULL) {
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) return(prod(p^w)),
                            w))
                    }
                )
            },"max",rc=rc)
            for (cc in names(originalPList)) {
                pc <- cbind(tempPList[[cc]],originalPList[[cc]])
                ly <- ncol(pc)
                sumpList[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
            #assign("perm.list",tempPList,envir=.GlobalEnv)
            #assign("oList",originalPList,envir=.GlobalEnv)
        },
        dperm_weight = {
            sumpList <- vector("list",length(cpList))
            names(sumpList) <- names(cpList)
            conl <- as.list(contrast)
            names(conl) <- contrast
            tempPList <- cmclapply(conl,metaPerm,
                counts=counts,sampleList=sampleList,
                statistics=statistics,statArgs=statArgs,
                libsizeList=libsizeList,
                nperm=nperm,weight=weight,
                select="weight",reprod=reprod,
                rc=rc)
            originalPList <- cmclapply(cpList,function(x,m,w=NULL) {
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) {return(prod(p^w))},
                            w))
                    }
                )
            },"weight",weight,rc=rc)
            for (cc in names(originalPList)) {
                pc <- cbind(tempPList[[cc]],originalPList[[cc]])
                ly <- ncol(pc)
                sumpList[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
            #assign("perm.list",tempPList,envir=.GlobalEnv)
            #assign("oList",originalPList,envir=.GlobalEnv)
        },
        none = {
            # A default value must be there to use with volcanos, we say the one
            # of the first statistic in order of input
            sumpList <- cmclapply(cpList,function(x) return(x[,1]),rc=rc)
        }
    )
    return(sumpList)
}

#' Permutation tests for meta-analysis
#'
#' This function performs permutation tests in order to derive a meta p-value by
#' combining several of the statistical algorithms of metaseqr. This is probably
#' the most accurate way of combining multiple statistical algorithms for RNA-Seq
#' data, as this issue is different from the classic interpretation of the term
#' "meta-analysis" which implies the application of the same statistical test on
#' different datasets treating the same subject/experiment. For other methods, see
#' also the main \code{\link{metaseqr}} help page. You should keep in mind that
#' the permutation procedure can take a long time, even when executed in parallel.
#'
#' @param counts a normalized read counts table, one row for each gene, one column
#' for each sample.
#' @param sampleList the list containing condition names and the samples under
#' each condition. See the main \code{\link{metaseqr}} help page.
#' @param contrast the contrasts to be tested by each statistical algorithm. See
#' the main \code{\link{metaseqr}} help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main
#' \code{\link{metaseqr}} help page.
#' @param statArgs the parameters for each statistical algorithm. See the main
#' \code{\link{metaseqr}} help page.
#' @param libsizeList a list with library sizes. See the main \code{\link{metaseqr}}
#' and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm.
#' @param select how to select the initial vector of p-values. It can be \code{"min"}
#' to select the minimum p-value for each gene (more conservative), \code{"max"}
#' to select the maximum p-value for each gene (less conservative), \code{"weight"}
#' to apply the weights to the p-value vector for each gene and derive a weighted
#' p-value.
#' @param replace same as the \code{replace} argument in the \code{\link{sample}}
#' function. Implies bootstraping or simple resampling without replacement. It can
#' also be \code{"auto"}, to determine bootstraping or not with the following rule:
#' if \code{ncol(counts)<=6} \code{replace=FALSE else} \code{replace=TRUE}. This
#' protects from the case of having zero variability across resampled conditions.
#' In such cases, most statistical tests would crash.
#' @param multic use multiple cores to execute the premutations. This is an 
#' external parameter and implies the existence of parallel package in the
#' execution environment. See the main \code{\link{metaseqr}} help page.
#' @param reprod create reproducible permutations. Ideally one would want to
#' create the same set of indices for a given dataset so as to create reproducible
#' p-values. If \code{reprod=TRUE}, a fixed seed is used by \code{metaPerm} for
#' all the datasets analyzed with \code{metaseqr}. If \code{reprod=FALSE}, then
#' the p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod} can
#' be a numeric vector of seeds with the same length as \code{nperm} so that the
#' user can supply his/her own seeds.
#' @return A vector of meta p-values
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
metaPerm <- function(contrast,counts,sampleList,statistics,statArgs,
    libsizeList,nperm=10000,weight=rep(1/ncol(counts),ncol(counts)),
    select=c("min","max","weight"),replace="auto",reprod=TRUE,rc=NULL) {
    checkTextArgs("select",select,c("min","max","weight"))
    if (replace=="auto") {
        if (ncol(counts)<=6)
            replace=FALSE
        else
            replace=TRUE
    }
    # We will construct relist in a way so that we can assign seeds for random
    # number generation and track progress at the same time
    if (is.logical(reprod)) {
        relist <- vector("list",nperm)
        if (reprod) {
            relist <- cmclapply(seq_along(relist),function(i) {
				return(list(seed=i,prog=i))
			},rc=rc)
        }
        else
            relist <- cmclapply(seq_along(relist),function(i) {
				return(list(seed=round(1e+6*runif(1)),prog=i))
			},rc=rc)
    }
    else if (is.numeric(reprod)) {
        if (length(reprod) != nperm)
            stopwrap("When reprod is numeric, it must have the same length as ",
                "nperm!")
        relist <- cmclapply(seq_along(reprod),function(i) {
            return(list(seed=reprod[i],prog=i))
        },rc=rc)
    }
    else
        stopwrap("reprod must be either a logical or a numeric vector!")
    disp("  Resampling procedure started...")
    # In this case, we must not use wapply as we want to be able to track progress
    # through mc.preschedule...
    if (!is.null(rc))
        pp <- mclapply(relist,metaWorker,counts,sampleList,contrast,
            statistics,replace,statArgs,libsizeList,select,weight,
            mc.preschedule=FALSE,mc.cores=getOption("cores"))
    else
        pp <- lapply(relist,metaWorker,counts,sampleList,contrast,statistics,
            replace,statArgs,libsizeList,select,weight)
    disp("  Resampling procedure ended...")
    return(do.call("cbind",pp))
}

#' Permutation tests helper
#'
#' This function performs the statistical test for each permutation. Internal use
#' only.
#'
#' @param x a virtual list with the random seed and the permutation index.
#' @param co the counts matrix.
#' @param sl the sample list.
#' @param cnt the contrast name.
#' @param s the statistical algorithms.
#' @param sa the parameters for each statistical algorithm.
#' @param ll a list with library sizes.
#' @param r same as the \code{replace} argument in the \code{\link{sample}} function.
#' @param el min, max or weight.
#' @param w the weights when \code{el="weight"}.
#' @return A matrix of p-values.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
metaWorker <- function(x,co,sl,cnt,s,r,sa,ll,el,w) {
    set.seed(x$seed)
    disp("    running permutation #",x$prog)
    pl <- makePermutation(co,sl,cnt,r)
    ppmat <- matrix(NA,nrow(co),length(s))
    colnames(ppmat) <- s
    for (alg in s) {
        #disp("      running permutation tests with: ",alg)
        tcl <- makeContrastList(pl$contrast,pl$sampleList)
        switch(alg,
            deseq = {
                pList <- suppressMessages(statDeseq(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            deseq2 = {
                pList <- suppressMessages(statDeseq2(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            edger = {
                pList <- suppressMessages(statEdger(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            noiseq = {
                pList <- suppressMessages(statNoiseq(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            bayseq = {
                pList <- suppressMessages(statBayseq(pl$counts,pl$sampleList,
                    tcl,sa[[alg]],ll))
            },
            limma = {
                pList <- suppressMessages(statLimma(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            nbpseq = {
                pList <- suppressMessages(statNbpseq(pl$counts,pl$sampleList,
                    tcl,sa[[alg]],ll))
            },
            absseq = {
                pList <- suppressMessages(statAbsseq(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            },
            dss = {
                pList <- suppressMessages(statDss(pl$counts,pl$sampleList,
                    tcl,sa[[alg]]))
            }
        )
        ppmat[,alg] <- as.numeric(pList[[1]])
    }
    ppmat[which(is.na(ppmat))] <- 1
    switch(el,
        min = {
            pIter <- apply(ppmat,1,min)
        },
        max = {
            pIter <- apply(ppmat,1,max)
        },
        weight = {
            pIter <- apply(ppmat,1,function(p,w) return(prod(p^w)),w)
        }
    )
    return(pIter)
}

#' Combine p-values with Simes' method
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using the Simes' method (see reference in the main \code{\link{metasqr}}
#' help page or in the vignette).
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combineSimes(p)
#'}
combineSimes <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    m <- length(p)
    y <- sort(p)
    s <- min(m*(y/(1:m)))
    return(min(c(s,1)))
}

#' Combine p-values with Bonferroni's method
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using the Bonferroni's method (see reference in the main
#' \code{\link{metasqr}} help page or in the vignette).
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combineBonferroni(p)
#'}
combineBonferroni <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    b <- length(p)*min(p)
    return(min(c(1,b)))
}

#' Combine p-values using weights
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using p-value weights.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @param w a weights vector, must sum to 1.
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combineWeight(p,w=c(0.2,0.5,0.3))
#'}
combineWeight <- function(p,w) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    return(prod(p^w))
}

#' Combine p-values using the minimum p-value
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR by taking the minimum p-value.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.min(p)
#'}
combineMinp <- function(p) { return(min(p)) }

#' Combine p-values using the maximum p-value
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR by taking the maximum p-value.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.max(p)
#'}
combineMaxp <- function(p) { return(max(p)) }

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisherMethod <- function(pvals,method=c("fisher"),p.corr=c("bonferroni","BH",
    "none"),zeroSub=0.00001,na.rm=FALSE,mc.cores=NULL) {
    stopifnot(method %in% c("fisher"))
    stopifnot(p.corr %in% c("none","bonferroni","BH"))
    stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
    stopifnot(zeroSub>=0 & zeroSub<=1 || length(zeroSub)!=1)
    if(is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
    ##substitute p-values of 0
    pvals[pvals == 0] <- zeroSub
    if(is.null(mc.cores)) {
        fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisherSum,
            zeroSub=zeroSub,na.rm=na.rm)))
    } 
    else {
        fisher.sums <- parallel::mclapply(1:nrow(pvals), function(i) {
            fisherSum(pvals[i,],zeroSub=zeroSub,na.rm=na.rm)
        }, mc.cores=mc.cores)
        fisher.sums <- data.frame(do.call(rbind,fisher.sums))
    }
    
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
      bonferroni = p.adjust(fisher.sums$p.value,"bonferroni"),
      BH = p.adjust(fisher.sums$p.value,"BH"),
      none = fisher.sums$p.value
    )
    return(fisher.sums)
}

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisherMethodPerm <- function(pvals,p.corr=c("bonferroni","BH","none"),
    zeroSub=0.00001,B=10000,mc.cores=NULL,blinker=1000) {
    stopifnot(is.na(blinker) || blinker>0)
    stopifnot(p.corr %in% c("none","bonferroni","BH"))
    stopifnot(all(pvals>=0,na.rm=TRUE) & all(pvals<=1,na.rm=TRUE))
    stopifnot(zeroSub>=0 & zeroSub<=1 || length(zeroSub)!=1)
    if(is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr)!=1,"BH",p.corr)
    pvals[pvals==0] <- zeroSub
  
    res.perm <- lapply(1:nrow(pvals),function(i) {
        if(!is.na(blinker) & i%%blinker==0)
        message("=", appendLF=FALSE)
        ##which studies contribute to S (don't have a NA in row i)
        good.p <- which(!is.na(pvals[i,]))
        S.obs= fisherSum(pvals[i,good.p], na.rm=FALSE)
        if(is.null(mc.cores)) {
            S.rand <- unlist(lapply(1:B, function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- sapply(good.p, function(pc){
                sample(na.exclude(pvals[,pc]),1)
            })
            fisherSum(myp)$S
        }))
        } else {
        S.rand <- unlist(parallel::mclapply(1:B, function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- sapply(good.p, function(pc) {
                sample(na.exclude(pvals[,pc]),1)
            })
            fisherSum(myp)$S
            }, mc.cores=mc.cores))
        }
        p.value <- sum(S.rand>=S.obs$S)/B
        data.frame(S=S.obs$S, num.p=S.obs$num.p, p.value=p.value)
    })
    res.perm <- data.frame(do.call(rbind, res.perm))
  
    if(!is.na(blinker) && blinker>0)
        message()
    ## rownames(res.perm) <- rownames(pvals)
    res.perm$p.adj <- switch(p.corr,
      bonferroni = p.adjust(res.perm$p.value,"bonferroni"),
      BH = p.adjust(res.perm$p.value,"BH"),
      none = res.perm$p.value)
    return(res.perm)
}

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisherSum <- function(p,zeroSub=0.00001,na.rm=FALSE) {
    if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
        stop("You provided bad p-values")
    stopifnot(zeroSub>=0 & zeroSub<=1 || length(zeroSub)!=1)
    p[p==0] <- zeroSub
    if (na.rm)
        p <- p[!is.na(p)]
    S = -2*sum(log(p))
    res <- data.frame(S=S,num.p=length(p))
    return(res)
}
