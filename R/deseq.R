################################################################################

load(system.file("extra/scvBiasCorrectionFits.rda",package="metaseqR2"))

# class_and_slots.R

setClass(".CountDataSet",
    contains="eSet",
    representation=representation(
        fitInfo="environment",
        dispTable="character",
        multivariateConditions="logical"
    ),
    prototype=prototype(new("VersionedBiobase",
        versions=c(classVersion("eSet"),.CountDataSet="1.1.0")
    ))
)

newCountDataSet <- function(countData,conditions,sizeFactors=NULL,
    phenoData=NULL,featureData=NULL) {
    countData <- as.matrix(countData)
    if (any(round(countData) != countData))
        stop("The countData is not integer.")
    mode(countData) <- "integer"

    if (is.null(sizeFactors))
        sizeFactors <- rep(NA_real_,ncol(countData))
    else
        warning("The 'sizeFactor' argument is deprecated. Use ",
            "'estimateSizeFactors'.")
    if (is.null(phenoData))
        phenoData <- annotatedDataFrameFrom(countData,byrow=FALSE)
    if (is.null(featureData))
        featureData <- annotatedDataFrameFrom(countData,byrow=TRUE)
        
    phenoData$`sizeFactor` <- sizeFactors
    varMetadata(phenoData)["sizeFactor","labelDescription"] <-
        "size factor (relative estimate of sequencing depth)"
    
    if (is(conditions,"matrix"))
        conditions <- as.data.frame(conditions)
    
    if (is(conditions,"data.frame") || is(conditions,"AnnotatedDataFrame")) {
        stopifnot(nrow(conditions) == ncol(countData))
        conditions <- as(conditions,"AnnotatedDataFrame")
        dimLabels(conditions) <- dimLabels(phenoData)
        rownames(pData(conditions)) <- rownames(pData(phenoData))
        phenoData <- combine(phenoData,conditions)
        multivariateConditions <- TRUE
        rvft <- c(`_all`=NA_character_)
    } 
    else {
        conditions <- factor(conditions)
        stopifnot(length(conditions) == ncol(countData))
        phenoData$`condition` <- factor(conditions)
        varMetadata(phenoData)["condition","labelDescription"] <-
            "experimental condition,treatment or phenotype"
        multivariateConditions <- FALSE
        rvft <- rep(NA_character_,length(levels(conditions)))
    }
    
    cds <- new(".CountDataSet",
        assayData=assayDataNew("environment",counts=countData),
        phenoData=phenoData,
        featureData=featureData,
        multivariateConditions=multivariateConditions,
        fitInfo=new.env(hash=TRUE),
        dispTable=rvft)
                
    cds
}

setValidity(".CountDataSet",function(object) {
    if (length(object@multivariateConditions) != 1)
        return("multivariateConditions is not scalar.")
    if (!"sizeFactor"  %in% names(pData(object)))
        return("phenoData does not contain a 'sizeFactor' columns.")
    if (!is(pData(object)$`sizeFactor`,"numeric"))
        return("The 'sizeFactor' column in phenoData is not numeric.")
    if (!object@multivariateConditions) {
        if (!"condition" %in% names(pData(object)))
            return("phenoData does not contain a 'condition' columns.")
        if (!is(pData(object)$`condition`,"factor"))
            return("The 'condition' column in phenoData is not a factor.")
    }
    if (!is.integer(counts(object)))
        return("the count data is not in integer mode")
    if (any(counts(object) < 0))
        return("the count data contains negative values")
    TRUE
})

setMethod("counts",signature(object=".CountDataSet"),
    function(object,normalized=FALSE) {
    if(!normalized)
        assayData(object)[["counts"]] 
    else {
        if(any(is.na(sizeFactors(object))))
            stop("Please first calculate size factors or set ",
                "normalized=FALSE")
        else
            t(t(assayData(object)[["counts"]])/sizeFactors(object))
    }
})

setReplaceMethod("counts",signature(object=".CountDataSet",value="matrix"),
    function(object,value) {
        assayData(object)[["counts"]] <- value
        validObject(object)
        object
})

setMethod("sizeFactors",signature(object=".CountDataSet"),function(object) {
    sf <- pData(object)$sizeFactor
    names(sf) <- colnames(counts(object))
    sf
})

setReplaceMethod("sizeFactors",signature(object=".CountDataSet",
    value="numeric"),
    function(object,value) {
        pData(object)$sizeFactor <- value
        validObject(object)
        object
})

setMethod("conditions",signature(object=".CountDataSet"),function(object,...) {
    if (length(list(...))!=0)
        warning("in conditions: Ignoring second and/or further arguments.")
    if (object@multivariateConditions 
        && !("condition" %in% colnames(pData(object))))
    stop("Could not find 'condition' column in pData.")
    conds <- pData(object)$`condition`
    names(conds) <- colnames(counts(object))
    conds
})

setReplaceMethod("conditions",signature(object=".CountDataSet"),
    function(object,value) {
        if (object@multivariateConditions)
            stop("The 'conditions<-' accessor is only for simple ",
                "single-factor conditions,but you have specified multivariate ",
                "conditions. Access them via 'pData<-'.")
        pData(object)$`condition` <- factor(value)
        validObject(object)
        object
})

setMethod("dispTable",signature(object=".CountDataSet"),function(object) {
    object@dispTable
})

setReplaceMethod("dispTable",signature(object=".CountDataSet"),
    function(object,value) {
        object@dispTable <- value
        validObject(object)
        object
})

fitInfo <- function(cds,name=NULL) {
    stopifnot(is(cds,".CountDataSet"))
    if (length(ls(cds@fitInfo)) == 0)
        stop("No fits available. Call 'estimateDispersions' first.")
    if (length(ls(cds@fitInfo)) > 1 && is.null(name))
        stop("More than one fitInfo object available. Specify by name. (See ",
            "'ls(cds@fitInfo)' for a list.)")
    if (length(ls(cds@fitInfo)) == 1 && is.null(name))
        name=ls(cds@fitInfo)[1]
    cds@fitInfo[[name]]
}

################################################################################

################################################################################

# core.R

estimateSizeFactorsForMatrix <- function(counts,locfunc=median) {
    loggeomeans <- rowMeans(log(counts))
    apply(counts,2,function(cnts)
        exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans)])))
}

getBaseMeansAndVariances <- function(counts,sizeFactors) {
    # Divides the counts by sizeFactors and calculates the estimates for
    # base means and variances for each gene.
    data.frame(
        baseMean=rowMeans(t(t(counts)/sizeFactors)),
        baseVar=rowVars(t(t(counts)/sizeFactors))
    )
}

modelMatrixToConditionFactor <- function(modelMatrix) {
    nr <- nrow(modelMatrix)
    if (nr<2)
        stop("nrow(modelMatrix) must be >=2.")
    mmconds <- seq_len(nr)
    for (i in 2:nr) {
        for (j in seq_len(i-1)) {
            if (all(modelMatrix[i,] == modelMatrix[j,])) {
                mmconds[i] <- mmconds[j]
                break
            }
        }
    }
    factor(as.integer(factor(mmconds)))
}

getBaseMeansAndPooledVariances <- function(counts,sizeFactors,conditions) {
    basecounts <- t(t(counts)/sizeFactors)
    replicated_sample <- conditions %in% names(which(table(conditions)>1))
    df <- sum(replicated_sample) - length(unique(conditions[replicated_sample]))
    data.frame(
        baseMean=rowMeans(basecounts),
        baseVar=rowSums(
            vapply(tapply((seq_len(ncol(counts)))[replicated_sample],
                factor(conditions[replicated_sample]),function(cols)
                rowSums((basecounts[,cols] - rowMeans(basecounts[,cols]))^2)),
                identity,numeric(1))
        )/df
    )
}

parametricDispersionFit <- function(means,disps) {
    coefs <- c(.1,1)
    iter <- 0
    while(TRUE) {
        residuals <- disps/(coefs[1] + coefs[2]/means)
        good <- which((residuals > 1e-4) & (residuals < 15))
        fit <- glm(disps[good] ~ I(1/means[good]),
            family=Gamma(link="identity"),start=coefs)
        oldcoefs <- coefs
        coefs <- coefficients(fit)
        if (!all(coefs > 0))
            stop("Parametric dispersion fit failed. Try a local fit and/or a ",
            "pooled estimation. (See '?estimateDispersions')")
        if (sum(log(coefs/oldcoefs)^2) < 1e-6)
            break
        iter <- iter + 1
        if (iter > 10) {
            warning("Dispersion fit did not converge.")
            break 
        }
    }

    names(coefs) <- c("asymptDisp","extraPois")
    ans <- function(q) {
        coefs[1] + coefs[2]/q
    }
    attr(ans,"coefficients") <- coefs
    ans
}

estimateAndFitDispersionsFromBaseMeansAndVariances <- function(means,variances,
    sizeFactors,fitType=c("parametric","local"),locfit_extra_args=list(),
    lp_extra_args=list(),adjustForBias=TRUE) {
    fitType <- match.arg(fitType)

    xim <- mean(1/sizeFactors)
    dispsAll <- (variances - xim*means)/means^2

    variances <- variances[means > 0]
    disps <- dispsAll[means > 0]
    means <- means[means > 0]

    if (adjustForBias)
        disps <- adjustScvForBias(disps,length(sizeFactors))

    if (fitType == "local") {
        fit <- do.call("locfit",c(
            list(
                variances ~ do.call("lp",c(list(log(means)),lp_extra_args)),
                family="gamma"
            ),
            locfit_extra_args
        ))

        rm(means)
        rm(variances)

        if (adjustForBias)
            ans <- function(q) {
                adjustScvForBias(
                    pmax((safepredict(fit,log(q)) - xim*q)/q^2,1e-8),
                    length(sizeFactors)
                )
            }
        else
            ans <- function(q) {
                pmax((safepredict(fit,log(q)) - xim*q)/q^2,1e-8)
            }

        # Note: The 'pmax' construct above serves to limit the overdispersion 
        # to a minimum of 10^-8,which should be indistinguishable from 0 but 
        # ensures numerical stability.
    } else if (fitType == "parametric") {
        ans <- parametricDispersionFit(means,disps)
    } else
        stop("Unknown fitType.")

    attr(ans,"fitType") <- fitType
    list(disps=dispsAll,dispFunc=ans)
}

profileLogLikelihood <- function(disp,mm,y,muhat) {
    # calculate the log likelihood:
    if (length(disp) != length(y))
        disp <- rep(disp,length(y))

    ll <- sum(vapply(seq(along=y),function(i)
        dnbinom(y[i],mu=muhat[i],size=1/disp[i],log=TRUE),numeric(1)))

    # transform the residuals,i.e.,y - muhat,to the linear
    # predictor scale by multiplying them with the derivative
    # of the link function,i.e.,by 1/muhat,and add this to the
    # linear predictors,log(muhat),to get the predictors that
    # are used in the IWLS regression
    z <- log(muhat) + (y - muhat)/muhat

    # the variance function of the NB is as follows
    v0 <- muhat + disp*muhat^2

    # transform the variance vector to linear predictor scale by
    # multiplying with the squared derivative of the link to
    # get the (reciprocal) weights for the IWLS
    w <- 1/((1/muhat)^2*v0)

    # All we need from the IRLS run is the QR decomposition of
    # its matrix
    qrres <- qr(mm*sqrt(w))

    # from it,we extract we leverages and calculate the Cox-Reid
    # term:
    cr <- sum(log(abs(diag(qrres$qr)[seq_len(qrres$rank)])))

    # return the profile log likelihood:
    ll - cr
}

estimateAndFitDispersionsWithCoxReid <- function(counts,modelFormula,modelFrame,
    sizeFactors,fitType=c("parametric","local"),locfit_extra_args=list(),
    lp_extra_args=list(),initialGuess=.1) {
    if (as.character(modelFormula[2]) == "count")
        modelFormula <- modelFormula[-2]
    mm <- model.matrix(modelFormula,modelFrame)
    disps <- apply(counts,1,function(y) {
        fit <- try(
            glm.fit(mm,y,family=MASS::negative.binomial(initialGuess),
                offset=log(sizeFactors)),
            silent=TRUE
        )
        if (inherits(fit,"try-error"))
            NA
        else {
            if (df.residual(fit) == 0)
                stop("No residual degrees of freedom. Most likely the design ",
                    "is lacking sufficient replication.")
            exp(
                optimize(
                    function(logalpha)
                        profileLogLikelihood(exp(logalpha),mm,y,
                        fitted.values(fit)),log(c(1e-11,1e5)),
                        maximum=TRUE
                )$maximum
            )
        } 
    })

    means <- colMeans(t(counts)/sizeFactors)
    xim <- mean(1/sizeFactors)

    if (fitType == "local") {
        fit <- do.call("locfit",c(
            list(
                disps[means>0] ~ do.call("lp",c(list(log(means[means>0])),
                lp_extra_args)),
                family="gamma"),
            locfit_extra_args
        ))

        rm(means)

        ans <- function(q)
            pmax((safepredict(fit,log(q)) - xim*q)/q^2,1e-8)
            
        # Note: The 'pmax' construct above serves to limit the overdispersion 
        # to a minimum of 10^-8,which should be indistinguishable from 0 but 
        # ensures numerical stability.

    } else if (fitType == "parametric")
        ans <- parametricDispersionFit(means,disps)
    else
        stop("Unkknown fitType.")

    attr(ans,"fitType") <- fitType
    list(disps=disps,dispFunc=ans)
}

safepredict <- function(fit,x) {
    # A wrapper around predict to avoid the issue that predict.locfit cannot
    # propagate NAs and NaNs properly.

    res <- rep.int(NA_real_,length(x))
    res[is.finite(x)] <- predict(fit,x[is.finite(x)])
    res
}

nbinomTestForMatrices <- function(countsA,countsB,sizeFactorsA,sizeFactorsB,
    dispsA,dispsB) {

    kAs <- rowSums(cbind(countsA))
    kBs <- rowSums(cbind(countsB))

    mus <- rowMeans(cbind(
        t(t(countsA)/sizeFactorsA),
        t(t(countsB)/sizeFactorsB)))

    fullVarsA <- pmax(mus*sum(sizeFactorsA) + dispsA*mus^2*sum(sizeFactorsA^2),
        mus*sum(sizeFactorsA)*(1+1e-8))
    fullVarsB <- pmax(mus*sum(sizeFactorsB) + dispsB*mus^2*sum(sizeFactorsB^2),
        mus*sum(sizeFactorsB)*(1+1e-8))

    sumDispsA <- (fullVarsA - mus*sum(sizeFactorsA))/(mus*sum(sizeFactorsA))^2
    sumDispsB <- (fullVarsB - mus*sum(sizeFactorsB))/(mus*sum(sizeFactorsB))^2

    vapply(seq(along=kAs),function(i) {

        if (kAs[i] == 0 & kBs[i] == 0)
            return(NA)

        # probability of all possible counts sums with the same total count:
        ks <- 0 : (kAs[i] + kBs[i])
        ps <- dnbinom(ks,mu=mus[i]*sum(sizeFactorsA),size=1/sumDispsA[i]) *
            dnbinom(kAs[i] + kBs[i] - ks,mu=mus[i]*sum(sizeFactorsB),
            size=1/sumDispsB[i])

        # probability of observed count sums:
        pobs <- dnbinom(kAs[i],mu=mus[i]*sum(sizeFactorsA),size=1/sumDispsA[i])*
            dnbinom(kBs[i],mu=mus[i]*sum(sizeFactorsB),size=1/sumDispsB[i])

        stopifnot(pobs == ps[kAs[i]+1])
        if (kAs[i]*sum(sizeFactorsB) < kBs[i]*sum(sizeFactorsA))
            numer <- ps[seq_len((kAs[i]+1))]
        else
            numer <- ps[(kAs[i]+1) : length(ps)]
        min(1,2*sum(numer)/sum(ps))
    },numeric(1))
}

# Note: The following function is never called; it is here only for
# documentation purposes,as it has been used to produce the data object
# scvBiasCorrectionFits,which is stored in the file
# inst/scvBiasCorrectionFits.rda,gets loadewd by the line after this
# function and is used by the function adjustScvForBias
#
# To do: the correct place for this type of documentation is a vignette
#
prepareScvBiasCorrectionFits <- function(maxnrepl=15,mu=100000,ngenes=10000,
    true_raw_scv=c(seq(0,2,length.out=100)[-1],seq(2,10,length.out=20)[-1]))
    lapply(2:maxnrepl,function(m) {
        est_raw_scv <- vapply(true_raw_scv,function(alpha) {
            k <- matrix(rnbinom(ngenes*m,mu=mu,size=1/alpha),ncol=m)
            k <- k[rowSums(k)>0,]
            mean(rowVars(k)/rowMeans(k)^2) },numeric(1))
        locfit(true_raw_scv ~ lp(est_raw_scv,nn=.2))
})

adjustScvForBias <- function(scv,nsamples) {
    stopifnot(nsamples > 1)
    if (nsamples - 1 > length(scvBiasCorrectionFits))
        scv
    else
        ifelse(scv > .02,
            pmax(safepredict(scvBiasCorrectionFits[[nsamples-1]],scv),
            1e-8*scv),scv)
}

nbkd.sf <- function(r,sf) {
    fam <- list(
        family=sprintf("nbkd,r=%g",r),
        link="log_sf",
        linkfun=function(mu) { log(mu/sf) },
        linkinv=function(eta) { pmax(sf*exp(eta),.Machine$double.eps) },
        mu.eta=function (eta) { pmax(sf*exp(eta),.Machine$double.eps) },
        variance=function(mu) { mu + mu^2/r },
        dev.resids=function(y,mu,wt) {
            2*wt*(ifelse(y > 0,y*log(y/mu),0) +
                (r + y)*log((r+mu)/(r+y)))
        },
        initialize=expression({
            n <- rep.int(1,nobs) # What is n?
            mustart <- y + 0.1
        }),
        valid.mu <- function(mu) {all(mu > 0)},
        valid.eta <- function(eta) {TRUE},
        simulate <- NA
    )

    class(fam) <- "family"
    fam 
}

fitNbinomGLMsForMatrix <- function(counts,sizeFactors,rawScv,modelFormula,
    modelFrame,quiet=FALSE,reportLog2=TRUE,glmControl=list()) {
    stopifnot(length(sizeFactors) == ncol(counts))
    stopifnot(length(rawScv) == nrow(counts))
    stopifnot(nrow(modelFrame) == ncol(counts))

    stopifnot(is(modelFormula,"formula"))
    if (as.character(modelFormula[[1]]) != "~")
        stop("Formula does not have a '~' as top-level operator.")
    if (as.character(modelFormula[[2]]) != "count")
        stop("Left-hand side of model formula must be 'count'.")

    goodRows <- is.finite(rawScv) & rowSums(counts) > 0

    modelMatrix <- model.matrix(modelFormula[c(1,3)],modelFrame)
    #res <- t(sapply(which(goodRows),function(i) {
    #    if (!quiet & i %% 1000 == 0)
    #        cat('.')
    #    nbfam <- nbkd.sf(1/rawScv[i],sizeFactors)
    #    fit <- try(
    #        glm.fit(modelMatrix,counts[i,],family=nbfam,control=glmControl),
    #        silent=TRUE)
    #    if (!inherits(fit,"try-error"))
    #        c(
    #            coefficients(fit),
    #            deviance=deviance(fit),
    #            df.residual=fit$df.residual,
    #            converged=fit$converged
    #        )
    #    else {
    #        coefs <- rep(NA,ncol(modelMatrix))
    #        names(coefs)  <- colnames(modelMatrix)
    #        #warning(as.character(fit))
    #        c(coefs,deviance=NA,df.residual=NA,converged=FALSE)
    #     } 
    #}))
    res <- t(vapply(which(goodRows),function(i) {
        if (!quiet & i %% 1000 == 0)
            cat('.')
        nbfam <- nbkd.sf(1/rawScv[i],sizeFactors)
        fit <- try(
            glm.fit(modelMatrix,counts[i,],family=nbfam,control=glmControl),
            silent=TRUE)
        if (!inherits(fit,"try-error"))
            c(
                coefficients(fit),
                deviance=deviance(fit),
                df.residual=fit$df.residual,
                converged=fit$converged
            )
        else {
            coefs <- rep(NA,ncol(modelMatrix))
            names(coefs)  <- colnames(modelMatrix)
            #warning(as.character(fit))
            c(coefs,deviance=NA,df.residual=NA,converged=FALSE)
        }
    },vector("list",ncol(modelMatrix)+3)))

    if (!quiet)
        cat("\n")

    df.residual <- max(na.omit(res[,"df.residual"]))
    if (!all(na.omit(res[,"df.residual"] == df.residual))) {
        res[res[,"df.residual"] != df.residual,"deviance"] <- NA
        warning("Some deviances set to NA due to reduction in degrees of ",
            "freedom.")
    }

    # Put in the NAs
    res2 <- data.frame(
        row.names=row.names(counts),
        apply(res,2,function(col) {
            a <- rep(NA_real_,nrow(counts))
            a[goodRows] <- col
            a 
        }))
    colnames(res2) <- colnames(res)

    if (reportLog2)
        res2[,seq_len(ncol(res2)-3)] <- res2[,seq_len(ncol(res2)-3)]/log(2)

    res2$converged <- as.logical(res2$converged)

    res2 <- res2[,colnames(res) != "df.residual"]
    attr(res2,"df.residual") <- df.residual
    res2
}

################################################################################

################################################################################

# methods.R

setMethod("estimateSizeFactors",signature(object=".CountDataSet"),
    function(object,locfunc=median,...) {
        if (length(list(...)) != 0)
            warning("in estimateSizeFactors: Ignoring extra argument(s).")
    sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object),locfunc)
    object
})

setMethod("estimateDispersions",signature(object=".CountDataSet"),
    function(object,method=c("pooled","pooled-CR","per-condition","blind"),
        sharingMode=c("maximum","fit-only","gene-est-only"),
        fitType=c("parametric","local"),
        locfit_extra_args=list(),lp_extra_args=list(),
        modelFrame=NULL,modelFormula=count ~ condition,...) {
    
    stopifnot(is(object,".CountDataSet"))
    if (any(is.na(sizeFactors(object))))
        stop("NAs found in size factors. Have you called already ",
            "'estimateSizeFactors'?")
    method <- match.arg(method)
    sharingMode <- match.arg(sharingMode)
    fitType <- match.arg(fitType)
    if (length(list(...)) != 0)
        warning("in estimateDispersions: Ignoring extra argument(s).")
    if (object@multivariateConditions 
        && !method %in% c("blind","pooled","pooled-CR"))
        stop("You have specified multivariate conditions (i.e.,passed a data ",
            "frame with conditions). In this case,you cannot use method ",
            "'per-condition'.")
    if (sharingMode == "gene-est-only")
        warning("in estimateDispersions: sharingMode=='gene-est-only' will ",
            "cause inflated numbers of false positives unless you have many ",
            "replicates.")
    
    # Remove results from previous fits
    fData(object) <- fData(object)[,!colnames(fData(object)) %in% paste("disp",
        object@dispTable,sep="_"),drop=FALSE]
    object@dispTable <- character()
    object@fitInfo=new.env(hash=TRUE)

    if (method == "blind") {
        bmv <- getBaseMeansAndVariances(counts(object),sizeFactors(object))
        dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances(
            bmv$baseMean,bmv$baseVar,sizeFactors(object),fitType,
            locfit_extra_args,lp_extra_args)
        object@fitInfo[["blind"]] <- list(
            perGeneDispEsts=dispsAndFunc$disps,
            dispFunc=dispsAndFunc$dispFunc,
            fittedDispEsts=dispsAndFunc$dispFunc(bmv$baseMean),
            df=ncol(counts(object)) - 1,
            sharingMode=sharingMode 
        )

        if (object@multivariateConditions)
            dispTable(object) <- c("_all"="blind")
        else {
            a <- rep("blind",length(levels(conditions(object))))
            names(a) <- levels(conditions(object))
            object@dispTable <- a 
        }
    } 
    else if (method == "per-condition") {
        replicated <- names(which(tapply(conditions(object),
            conditions(object),length) > 1))
        if (length(replicated) < 1)
            stop("None of your conditions is replicated. Use method='blind' ",
                "to estimate across conditions,or 'pooled-CR',if you have ",
                "crossed factors.")
        nonreplicated <- names(which(tapply(conditions(object),
            conditions(object),length) == 1))
        overall_basemeans <- rowMeans(counts(object,normalized=TRUE))

        for(cond in replicated) {
            cols <- conditions(object)==cond
            bmv <- getBaseMeansAndVariances(counts(object)[,cols],
                sizeFactors(object)[cols])
            dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances(
                bmv$baseMean,bmv$baseVar,sizeFactors(object)[cols],fitType,
                locfit_extra_args,lp_extra_args)
            object@fitInfo[[cond]] <- list(
                perGeneDispEsts=dispsAndFunc$disps,
                dispFunc=dispsAndFunc$dispFunc,
                fittedDispEsts=dispsAndFunc$dispFunc(overall_basemeans),     
                df=sum(cols)-1,
                sharingMode=sharingMode
            ) 
        }
        object@dispTable <- vapply(levels(conditions(object)),function(cond)
            ifelse(cond %in% replicated,cond,"max"),character(1))
    } 
    else if (method == "pooled" || method == "pooled-CR") {
        if (method == "pooled") {
            if (object@multivariateConditions) {
                if (is.null(modelFrame))
                    modelFrame <- 
                        pData(object)[,colnames(pData(object)) != "sizeFactor"]
                conds <- modelMatrixToConditionFactor(modelFrame) 
            }
            else
                conds <- conditions(object)
            if (!any(duplicated(conds)))
                stop("None of your conditions is replicated. Use ",
                    "method='blind' to estimate across conditions,or ",
                    "'pooled-CR',if you have crossed factors.")
            bmv <- getBaseMeansAndPooledVariances(counts(object),
                sizeFactors(object),conds)
            baseMeans <- bmv$baseMean
            dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances(
                bmv$baseMean,bmv$baseVar,sizeFactors(object),fitType,
                locfit_extra_args,lp_extra_args)
            df <- ncol(counts(object)) - length(unique(conds))
        } 
        else {  # method == "pooled-CR"
            if (is.null(modelFrame))
                modelFrame <- 
                    pData(object)[,colnames(pData(object)) != "sizeFactor",
                        drop=FALSE]
            
            baseMeans <- rowMeans(counts(object,normalized=TRUE))
            dispsAndFunc <- estimateAndFitDispersionsWithCoxReid(counts(object),
                modelFormula,modelFrame,sizeFactors(object),fitType,
                locfit_extra_args,lp_extra_args)
            df <- NA
        }

        object@fitInfo[["pooled"]] <- list(
            perGeneDispEsts=dispsAndFunc$disps,
            dispFunc=dispsAndFunc$dispFunc,
            fittedDispEsts=dispsAndFunc$dispFunc(baseMeans),
            df=df,
            sharingMode=sharingMode
        )

        dt <- if (object@multivariateConditions) 
            c("_all"="pooled") 
        else 
            character(0)
        
        if ("condition" %in% colnames(pData(object))) {
            a <- rep("pooled",length(levels(conditions(object))))
            names(a) <- levels(conditions(object))
            dt=c(dt,a)
        }
        dispTable(object) <- dt
        
    } 
    else
        stop(sprintf("Invalid method '%s'.",method))

    for (n in ls(object@fitInfo))
        fData(object)[[paste("disp",n,sep="_")]] <-
            switch(sharingMode,
                `fit-only` = object@fitInfo[[n]]$fittedDispEsts,
                `gene-est-only` = {
                    a <- object@fitInfo[[n]]$perGeneDispEsts
                    a[is.nan(a)] <- 0
                    pmax(a,1e-8) 
                },
                `maximum` = pmax(object@fitInfo[[n]]$fittedDispEsts,
                    object@fitInfo[[n]]$perGeneDispEsts,na.rm=TRUE),
                stop(sprintf("Invalid sharingMode '%s'.",sharingMode))
            ) ## switch
    
    if ("max" %in% object@dispTable)
        fData(object)[["disp_max"]] <- do.call(pmax,
            c(fData(object)[,colnames(fData(object)) %in% paste("disp",
                object@dispTable,sep="_"),drop=FALSE],na.rm=TRUE))

    validObject(object)
    object
})

nbinomTest <- function(cds,condA,condB,pvals_only=FALSE,eps=NULL) {
    stopifnot(is(cds,".CountDataSet"))
    if (all(is.na(dispTable(cds))))
        stop("Call 'estimateDispersions' first.")

    if (dispTable(cds)[condA] == "blind" || dispTable(cds)[condB] == "blind") {
        if (fitInfo(cds,"blind")$sharingMode != "fit-only")
            warning('You have used \'method="blind"\' in estimateDispersion ",
                "without also setting \'sharingMode="fit-only"\'. This will ",
                "not yield useful results.')
    }

    stopifnot(condA %in% levels(conditions(cds)))
    stopifnot(condB %in% levels(conditions(cds)))
    if (!is.null(eps))
        warning("The 'eps' argument is defunct and hence ignored.")

    colA <- conditions(cds)==condA
    colB <- conditions(cds)==condB

    bmv <- getBaseMeansAndVariances(counts(cds)[,colA|colB],
        sizeFactors(cds)[colA|colB])

    rawScvA <- fData(cds)[,paste("disp",dispTable(cds)[condA],sep="_")]
    rawScvB <- fData(cds)[,paste("disp",dispTable(cds)[condB],sep="_")]

    pval <- nbinomTestForMatrices(
        counts(cds)[,colA],
        counts(cds)[,colB],
        sizeFactors(cds)[colA],
        sizeFactors(cds)[colB],
        rawScvA,
        rawScvB
    )

    if (pvals_only)
        pval
    else {
        bmvA <- getBaseMeansAndVariances(counts(cds)[,colA],
            sizeFactors(cds)[colA])
        bmvB <- getBaseMeansAndVariances(counts(cds)[,colB],
            sizeFactors(cds)[colB])
        data.frame(
            id=rownames(counts(cds)),
            baseMean=bmv$baseMean,
            baseMeanA=bmvA$baseMean,
            baseMeanB=bmvB$baseMean,
            foldChange=bmvB$baseMean/bmvA$baseMean,
            log2FoldChange=log2(bmvB$baseMean/bmvA$baseMean),
            pval=pval,
            padj=p.adjust(pval,method="BH"),
            stringsAsFactors=FALSE
        ) 
    }
}

getVarianceStabilizedData <- function(cds) {
    stopifnot(is(cds,".CountDataSet"))
    if ("blind" %in% ls(cds@fitInfo))
        fitInfo <- cds@fitInfo[["blind"]]
    else if ("pooled" %in% ls(cds@fitInfo))
        fitInfo <- cds@fitInfo[["pooled"]]
    else
        stop("Use 'estimateDispersions' with 'method=\"blind\"' ",
            "(or \"pooled\") before calling 'getVarianceStabilizedData'")
    ncounts <- t(t(counts(cds))/sizeFactors(cds))
    if (attr(fitInfo$dispFunc,"fitType") == "parametric") {
        coefs <- attr(fitInfo$dispFunc,"coefficients")
        vst <- function(q)
            log((1 + coefs["extraPois"] + 2*coefs["asymptDisp"]*q +
                2*sqrt(coefs["asymptDisp"]*q*(1 + coefs["extraPois"] + 
                coefs["asymptDisp"]*q)))
                /(4*coefs["asymptDisp"]))/log(2)
        vst(ncounts)
    }
    else {
        # non-parametric fit -> numerical integration
        xg <- sinh(seq(asinh(0),asinh(max(ncounts)),length.out=1000))[-1]
        xim <- mean(1/sizeFactors(cds))
        baseVarsAtGrid <- fitInfo$dispFunc(xg)*xg^2 + xim*xg
        integrand <- 1/sqrt(baseVarsAtGrid)
        splf <- splinefun(
            asinh((xg[-1] + xg[-length(xg)])/2),
            cumsum((xg[-1] - xg[-length(xg)]) *
                (integrand[-1] + integrand[-length(integrand)])/2)
        )
        h1 <- quantile(rowMeans(ncounts),.95)
        h2 <- quantile(rowMeans(ncounts),.999)
        eta <- (log2(h2) - log2(h1))/(splf(asinh(h2)) - splf(asinh(h1)))
        xi <- log2(h1) - eta*splf(asinh(h1))
        tc <- vapply(colnames(counts(cds)),function(clm)
            eta*splf(asinh(ncounts[,clm])) + xi,numeric(1))
        rownames(tc) <- rownames(counts(cds))
        tc
    }
}

varianceStabilizingTransformation <- function (cds) {
    new("ExpressionSet",
        exprs=getVarianceStabilizedData(cds),
        phenoData=phenoData(cds),
        featureData=featureData(cds),
        experimentData=experimentData(cds),
        annotation=annotation(cds),
        protocolData=protocolData(cds)
    )
}

fitNbinomGLMs <- function(cds,modelFormula,glmControl=list()) {
    stopifnot(is(cds,".CountDataSet"))
    if ("disp_pooled" %in% colnames(fData(cds)))
        disps <- fData(cds)$disp_pooled
    else if ("disp_blind" %in% colnames(fData(cds))) {
        if (fitInfo(cds,"blind")$sharingMode != "fit-only")
            warning('You have used \'method="blind"\' in estimateDispersion ",
                "without also setting \'sharingMode="fit-only"\'. This will ",
                "not yield useful results.')
        disps <- fData(cds)$disp_blind
    } 
    else
        stop("Call 'estimateDispersions' with 'method=\"pooled\"' ",
            "'(or 'blind') first.")

    fitNbinomGLMsForMatrix(counts(cds),sizeFactors(cds),disps,modelFormula,
        pData(cds),glmControl=glmControl)
}

nbinomGLMTest <- function(resFull,resReduced)
    1 - pchisq(resReduced$deviance - resFull$deviance,
    attr(resReduced,"df.residual") - attr(resFull,"df.residual"))

################################################################################
