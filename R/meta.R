metaTest <- function(cpList,metaP=c("simes","bonferroni","fisher",
    "dperm_min","dperm_max","dperm_weight","fperm","whitlock","minp","maxp",
    "weight","pandora","none"),counts,sampleList,statistics,statArgs,
    libsizeList,nperm=10000,weight=rep(1/length(statistics),
    length(statistics)),rc=NULL) {
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
                if (!is.null(rc))
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
            sumpList <- cmclapply(cpList,function(x) {
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
                select="min",rc=rc)
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
            for (cc in names(originalPList)) {
                pc <- cbind(tempPList[[cc]],originalPList[[cc]])
                ly <- ncol(pc)
                sumpList[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[seq_len(m-1)]<y[m]))/(m-1)),ly)
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
                select="max",rc=rc)
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
                    return(length(which(y[seq_len(m-1)]<y[m]))/(m-1)),ly)
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
                select="weight",rc=rc)
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
                    return(length(which(y[seq_len(m-1)]<y[m]))/(m-1)),ly)
            }
        },
        none = {
            # A default value must be there to use with volcanos, we say the one
            # of the first statistic in order of input
            sumpList <- cmclapply(cpList,function(x) return(x[,1]),rc=rc)
        }
    )
    return(sumpList)
}

metaPerm <- function(contrast,counts,sampleList,statistics,statArgs,
    libsizeList,nperm=10000,weight=rep(1/ncol(counts),ncol(counts)),
    select=c("min","max","weight"),replace="auto",rc=NULL) {
    checkTextArgs("select",select,c("min","max","weight"))
    if (replace=="auto") {
        if (ncol(counts)<=6)
            replace=FALSE
        else
            replace=TRUE
    }
    # We will construct relist in a way so that we can assign seeds for random
    # number generation and track progress at the same time
    relist <- vector("list",nperm)
    relist <- cmclapply(seq_along(relist),function(i) {
        return(list(seed=round(1e+6*runif(1)),prog=i))
    },rc=rc)
    disp("  Resampling procedure started...")
    # In this case, we must not use cmclapply as we want to be able to track 
    # progress through mc.preschedule...
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

metaWorker <- function(x,co,sl,cnt,s,r,sa,ll,el,w) {
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

combineSimes <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    m <- length(p)
    y <- sort(p)
    s <- min(m*(y/(seq_len(m))))
    return(min(c(s,1)))
}

combineBonferroni <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    b <- length(p)*min(p)
    return(min(c(1,b)))
}

combineWeight <- function(p,w) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    return(prod(p^w))
}

combineMinp <- function(p) { return(min(p)) }

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
        fisher.sums <- parallel::mclapply(seq_len(nrow(pvals)), function(i) {
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
  
    resPerm <- lapply(seq_len(nrow(pvals)),function(i) {
        if(!is.na(blinker) & i%%blinker==0)
        message("=", appendLF=FALSE)
        ##which studies contribute to S (don't have a NA in row i)
        good.p <- which(!is.na(pvals[i,]))
        S.obs= fisherSum(pvals[i,good.p], na.rm=FALSE)
        if(is.null(mc.cores)) {
            Srand <- unlist(lapply(seq_len(B), function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- vapply(good.p, function(pc){
                sample(na.exclude(pvals[,pc]),1)
            },numeric(1))
            fisherSum(myp)$S
        }))
        } else {
        Srand <- unlist(parallel::mclapply(seq_len(B), function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- vapply(good.p, function(pc) {
                sample(na.exclude(pvals[,pc]),1)
            },numeric(1))
            fisherSum(myp)$S
            }, mc.cores=mc.cores))
        }
        p.value <- sum(Srand>=S.obs$S)/B
        data.frame(S=S.obs$S, num.p=S.obs$num.p, p.value=p.value)
    })
    resPerm <- data.frame(do.call(rbind, resPerm))
  
    if(!is.na(blinker) && blinker>0)
        message()
    ## rownames(resPerm) <- rownames(pvals)
    resPerm$p.adj <- switch(p.corr,
      bonferroni = p.adjust(resPerm$p.value,"bonferroni"),
      BH = p.adjust(resPerm$p.value,"BH"),
      none = resPerm$p.value)
    return(resPerm)
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
