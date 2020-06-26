################################################################################

# classes.R

setClass("Biodetection",representation(dat="list"))
setClass("CD",representation(dat="list"))
setClass("CountsBio",representation(dat="list"))
setClass("GCbias",representation(dat="list"))
setClass("lengthbias",representation(dat="list"))
setClass("Saturation",representation(dat="list"))
setClass("PCA",representation(dat="list"))

setGeneric("explo.plot",function(object,...) standardGeneric("explo.plot"))

setMethod("explo.plot","Biodetection",
    function(object,samples=c(1,2),plottype=c("persample","comparison"),
        toplot="protein_coding",...)
        biodetection.plot(object@dat,samples=samples,plottype=plottype,
            toplot=toplot,...))

#setMethod("explo.plot","CD",function(object,samples=NULL,...) 
#    cd.plot(object@dat,samples=samples,...))

setMethod("explo.plot","CountsBio",
    function(object,samples=c(1,2),toplot="global",
        plottype=c("barplot","boxplot"),...)
        countsbio.plot(object@dat,samples,toplot,plottype,...))

#setMethod("explo.plot","GCbias",
#    function(object,samples=NULL,toplot="global",...)
#    GC.plot(object@dat,samples=samples,toplot=toplot,...))

#setMethod("explo.plot","lengthbias",
#    function(object,samples=NULL,toplot="global",...)
#    length.plot(object@dat,samples=samples,toplot=toplot,...))

#setMethod("explo.plot","Saturation",
#    function(object,samples=NULL,toplot=1,yleftlim=NULL,yrightlim=NULL,...)
#    saturation.plot(object@dat,samples=samples,toplot=toplot,yleftlim=yleftlim,
#        yrightlim=yrightlim,...))

#setMethod("explo.plot","PCA",
#    function(object,samples=c(1,2),plottype="scores",factor=NULL)
#    PCA.plot(object@dat,samples=samples,plottype=plottype,factor=factor))


# Show methods for exploration objects
setMethod("show","Biodetection",function(object) {
    cat("\n Reference genome: \n==========\n")
    names(dimnames(object@dat$genome))=NULL
    print(object@dat$genome)
    for (i in seq_len(length(object@dat$biotables))) {
        cat("\n",names(object@dat$biotables)[i],"\n==========\n")
        print(object@dat$biotables[[i]])
    }
})

setMethod("show","CD",function(object) {
    cat("\n Confidence intervals for median of M to compare each sample ",
        "to reference:\n=======\n")
    print(object@dat$DiagnosticTest)
    cat("\n Reference sample is:\n=======\n")
    print(object@dat$refColumn)
})

setMethod("show","CountsBio",function(object) {
    cat("\n Summary: \n============\n")              
    print(object@dat$summary[[1]])            
})

setMethod("show","GCbias",function(object) {
    x <- object@dat$RegressionModels
    for (i in seq_len(length(x))) {
        print(names(x)[i])
        print(summary(x[[i]]))
    }
})

setMethod("show","lengthbias",function(object) {
    x <- object@dat$RegressionModels
    for (i in seq_len(length(x))) {
        print(names(x)[i])
        print(summary(x[[i]]))
    }
})

setMethod("show","Saturation",function(object) {
    x <- dat2save(object)
    cat("\n Number of detected features at each sequencing ",
        "depth: \n============\n")
    for (i in seq_len(length(x))) {
        print(names(x)[i])
        print(x[[i]])
    }            
})

setMethod("show","PCA",function(object) {
    x <- object$result$var.exp
    x <- round(x*100,4)
    colnames(x) <- c("%Var","Cum %Var")
    rownames(x) <- paste("PC",seq_len(nrow(x)))
    cat("\n Percentage of total variance explained by each component: ",
        "\n============\n")
    print(x)
})


# Coercion methods for exploration objects
setGeneric("dat2save",function(object) standardGeneric("dat2save"))
setMethod("dat2save","Biodetection",function(object) object@dat)
setMethod("dat2save","CD",function(object) object@dat)
setMethod("dat2save","CountsBio",function(object) object@dat$summary)
setMethod("dat2save","GCbias",function(object) object@dat$data2plot)
setMethod("dat2save","lengthbias",function(object) object@dat$data2plot)
setMethod("dat2save","Saturation",function(object) {
    muestras=vector("list",length=length(object@dat$depth))
    names(muestras)=names(object@dat$depth)
    for (i in seq_len(length(muestras))) {
        muestras[[i]] <- object@dat$depth[[i]]
        for (j in seq_len(length(object@dat$saturation))) {
            muestras[[i]] <- cbind(muestras[[i]],
                object@dat$saturation[[j]][[i]])
        }
        colnames(muestras[[i]]) <- c("depth",names(object@dat$saturation))
    }

    muestras
})

setMethod("dat2save","PCA",function(object) object@dat$result)

setClass("myInfo",representation(method="character",k="numeric",
    lc="numeric",factor="vector",v="numeric",nss="numeric",pnr="numeric",
    comparison="vector",replicates="character"))
setClass("Output",representation(results="list"),contains="myInfo")

setValidity("Output",function(object) {
    if (!(is.character(object@method))) {
        return(paste("Method must be a string"))
    } else if (!(is.numeric(object@k))) {
        return(paste("k must be numeric"))
    } else if (!(is.numeric(object@lc))) {
        return(paste("lc must be numeric"))
    } else if (!(is.vector(object@factor))) {
        return(paste("Factor must be a vector of strings"))
    } else if (!(is.numeric(object@v))) {
        return(paste("v must be numeric"))
    } else if (!(is.numeric(object@nss))) {
        return(paste("nss must be numeric"))
    } else if (!(is.numeric(object@pnr))) {
        return(paste("pnr must be numeric"))
    } else if (!(is.vector(object@comparison))) {
        return(paste("Comparison must be a vector of strings"))
    } else if (!(is.list(object@results))) {
        return(paste("Results must be a list of data.frames"))
    } else {
        return(TRUE)
    }
})

Output <- function (data,method,k,lc,factor,v,nss,pnr,comparison,replicates) {
    new("Output",results=data,method=method,k=k,lc=lc,factor=factor,v=v,nss=nss,
        pnr=pnr,comparison=comparison,replicates=replicates)
}

setMethod("show","Output",function(object) {
    if (object@method == "n")
        object@method="none"
    for (i in seq_len(length(object@results))) {
        cat("\nSummary",i,"\n=========\n")
        cat("\nYou are comparing",object@comparison[i],"from",object@factor[i],
            "\n\n")
        print(head(object@results[[i]][order(object@results[[i]][,5],
            decreasing=TRUE),]))
    }
    cat("\nNormalization\n")
    cat("\tmethod:",object@method,"\n")
    cat("\tk:",object@k,"\n")
    cat("\tlc:",object@lc,"\n")

    # Simulated samples
    if (object@replicates == "no") {
        cat("\nYou are working with simulated replicates:\n")
        cat("\tpnr:",object@pnr,"\n")
        cat("\tnss:",object@nss,"\n")
        cat("\tv:",object@v,"\n")
    }
    # With biological or technical replicates
    else {
        cat("\nYou are working with",object@replicates,"replicates\n")
    }           
})

################################################################################

################################################################################

# allMD.R

allMD <- function (input,factor,conditions,k=0.5,replicates,norm="rpkm",pnr=0.2,
    nss=5,v=0.02,lc=0) {
    # Check if the factor introduced is already defined
    # If the factor introduced is defined and has more than 2 conditions,it will
    # check if the conditions specified are defined too
    condition_fac <- FALSE
    condition_lev <- FALSE

    datos1 <- datos2 <- matrix()

    for (i in colnames(pData(input))) {
        if (factor == i) {
            condition_fac <- TRUE
            if (!is.factor(pData(input)[,i])) 
                pData(input)[,i]=as.factor(pData(input)[,i])
            if (length(levels(pData(input)[,i])) == 2) {
                if (!is.null(assayData(input)$exprs)) {
                    datos1 <- 
                        assayData(input)$exprs[,which(
                            pData(input)[,i]==levels(pData(input)[,i])[1]),
                            drop=FALSE]
                    datos2 <- 
                        assayData(input)$exprs[,which(
                            pData(input)[,i] ==levels(pData(input)[,i])[2]),
                            drop=FALSE]
                } 
                else {
                    datos1 <- 
                        assayData(input)$counts[,which(
                            pData(input)[,i]==levels(pData(input)[,i])[1]),
                            drop=FALSE]
                    datos2 <- assayData(input)$counts[,which(
                        pData(input)[,i]==levels(pData(input)[,i])[2]),
                            drop=FALSE]
                }

                # Define the comparison string
                comparison <- paste(levels(pData(input)[,i])[1],
                    levels(pData(input)[,i])[2],sep=" - ")

                condition_lev <- TRUE
            }

            else {
                if (is.null(conditions))
                    stop("Error. You must specify which conditions you wish ",
                        "to compare when the factor has two or more ",
                        "conditions.\n")
                if (length(conditions) != 2)
                    stop("Error. The argument conditions must contain the 2 ",
                        "conditions you wish to compare.")
            
                l <- conditions %in% pData(input)[,i]
                # If they are defined,they will be TRUE
                if (l[1] == TRUE && l[2] == TRUE) {
                    if (!is.null(assayData(input)$exprs)) {
                        datos1 <- assayData(input)$exprs[,which(
                            pData(input)[,i] == conditions[1]),drop=FALSE]
                        datos2 <- assayData(input)$exprs[,which(
                            pData(input)[,i] == conditions[2]),drop=FALSE]
                    }
                    else {
                        datos1 <- assayData(input)$counts[,which(
                            pData(input)[,i] == conditions[1]),drop=FALSE]
                        datos2 <- assayData(input)$counts[,which(
                            pData(input)[,i] == conditions[2]),drop=FALSE]
                    }
                    # Define the comparison string
                    comparison <- paste(conditions[1],conditions[2],sep=" - ")
                    condition_lev <- TRUE
                }
            }
        }
    }

    if (!condition_fac)
        stop("The factor you have written does not correspond with any of the ",
            "ones you have defined.")
            
    if (!condition_lev)
        stop("The conditions you have written don't exist in the factor ",
            "specified.\n")

    # Correction to make it work when there are simulated samples
    if (replicates == "no")
        replicates <- "technical"

    n1 <- ncol(as.matrix(datos1))
    n2 <- ncol(as.matrix(datos2))

    if (norm == "n") { # no normalization
        datos1 <- round(datos1,100)
        datos2 <- round(datos2,100)
    }

    if (is.null(k)) {
        m1 <- min(datos1[noceros(datos1,num=FALSE)],na.rm=TRUE)
        m2 <- min(datos2[noceros(datos2,num=FALSE)],na.rm=TRUE)
        mm <- min(m1,m2)
        k <- mm/2    
    } 

    # Total counts for each gene:
    suma1 <- rowSums(as.matrix(datos1))
    suma2 <- rowSums(as.matrix(datos2))

    # All genes
    todos <- rownames(as.matrix(datos1))

    # Genes with counts in any condition
    concounts <- names(which(suma1+suma2 > 0))

    long <- 1000
    g.sinL <- NULL

    if (!is.null(featureData(input)@data$Length)) {
        g.sinL <- names(which(is.na(featureData(input)@data$Length)))
        if (any(!is.na(featureData(input)@data$Length)) == TRUE) 
        long <- featureData(input)@data[concounts,"Length"]
    }
  
    if (replicates == "technical") { ### technical replicates
        suma1 <- suma1[concounts]
        suma2 <- suma2[concounts]    

        #----------------------------------------------------------------------#
        # Normalization of counts for each condition (aggregating replicates)

        if (norm == "rpkm") {      # RPKM
            suma1.norm <- noirpkm(suma1,long=long,k=k,lc=lc)
            suma2.norm <- noirpkm(suma2,long=long,k=k,lc=lc)
        }

        if (norm == "uqua") {
            suma.norm <- noiuqua(cbind(suma1,suma2),long=long,lc=lc,k=k)
            suma1.norm <- as.matrix(suma.norm[ ,1])
            suma2.norm <- as.matrix(suma.norm[ ,2])
        }

        if (norm == "tmm") {
            suma.norm <- noitmm(as.matrix(cbind(suma1,suma2)),long=long,lc=lc,
                k=k)
            suma1.norm <- as.matrix(suma.norm[ ,1])
            suma2.norm <- as.matrix(suma.norm[ ,2])
        }
    }

    #-------------------------------------------------------------------------#

    ## Noise distribution

    if ((n1+n2)>2) {   # with real samples
        datitos <- cbind(datos1,datos2)
        datitos <- datitos[concounts,]

        gens.sin0 <- setdiff(concounts,g.sinL)

        if (norm == "n") {       # no normalization
            datitos.0 <- sinceros(datitos,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "rpkm") {      # RPKM
            datitos.0 <- noirpkm(datitos,long=long,k=k,lc=lc)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "uqua") {      # Upper Quartile
            datitos.0 <- noiuqua(datitos,long=long,lc=lc,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "tmm") {
            datitos.0 <- noitmm(datitos,long=long,lc=lc,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        datos1.norm <- datitos.norm[ ,seq_len(n1)]
        datos2.norm <- datitos.norm[ ,(n1+1):(n1+n2)]

        if (n1 > 1) {
            MD1 <- MD(dat=datos1.norm)
        } 
        else { 
            MD1 <- NULL 
        }

        if (n2 > 1) {
            MD2 <- MD(dat=datos2.norm)
        } 
        else { 
            MD2 <- NULL 
        }
    } 
    else { # with simulated samples
        if (nss == 0) {
            nss <- 5
        }
        
        datos.sim <- sim.samples(counts1=sinceros(suma1,k=k),
            counts2=sinceros(suma2,k=k),pnr=pnr,nss=nss,v=v)

        nn <- vapply(datos.sim,ncol,numeric(1))
        dat.sim.norm <- vector("list",length=2)
        datitos <- cbind(datos.sim[[1]],datos.sim[[2]])
        rownames(datitos)=names(suma1)

        sumita <- rowSums(datitos)
        g.sin0 <- names(which(sumita > 0))
        gens.sin0 <- setdiff(g.sin0,g.sinL)

        if (norm == "n") {       # no normalization
            datitos.0 <- sinceros(datitos,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "rpkm") {      # RPKM
            datitos.0 <- noirpkm(datitos,long=long,k=k,lc=lc)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "uqua") {      # Upper Quartile
            datitos.0 <- noiuqua(datitos,long=long,lc=lc,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        if (norm == "tmm") {
            datitos.0 <- noitmm(datitos,long=long,lc=lc,k=k)
            datitos.norm <- datitos.0[gens.sin0,]
        }

        dat.sim.norm[[1]] <- datitos.norm[,seq_len(nn[1])]
        dat.sim.norm[[2]] <- datitos.norm[,(nn[1]+1):sum(nn)]

        MD1 <- MD(dat=dat.sim.norm[[1]])
        MD2 <- MD(dat=dat.sim.norm[[2]])
    }

    Mr <- c(as.numeric(MD1$M),as.numeric(MD2$M))
    Dr <- c(as.numeric(MD1$D),as.numeric(MD2$D))

    #-------------------------------------------------------------------------#

    ## M and D for different experimental conditions

    if (replicates == "technical" & norm != "n") {
        MDs <- MD(dat=cbind(suma1.norm,suma2.norm))

        lev1 <- suma1.norm[,1]
        lev1 <- lev1[todos]
        
        lev2 <- suma2.norm[,1]
        lev2 <- lev2[todos]
    } 
    else {
        if ((n1+n1) == 2) {
            datos1.norm <- sinceros(as.matrix(datos1)[concounts,],k=k)
            datos2.norm <- sinceros(as.matrix(datos2)[concounts,],k=k)      
        }  

        resum1.norm <- rowMeans(as.matrix(datos1.norm))
        resum2.norm <- rowMeans(as.matrix(datos2.norm))

        lev1 <- resum1.norm[todos]
        lev2 <- resum2.norm[todos]

        MDs <- MD(dat=cbind(resum1.norm,resum2.norm))
    }

    ## Completing M and D
    names(lev1) <- names(lev2) <- todos

    Ms <- as.numeric(MDs$M)
    names(Ms) <- rownames(MDs$M)
    Ms <- Ms[todos]
    names(Ms) <- todos

    Ds <- as.numeric(MDs$D)
    names(Ds) <- rownames(MDs$D)
    Ds <- Ds[todos]
      names(Ds) <- todos

    ## Results
    list("k"=k,"comp"=comparison,"Level1"=lev1,"Level2"=lev2,"Ms"=Ms,"Ds"=Ds,
        "Mn"=Mr,"Dn"=Dr)
}

allMDbio <- function(input,factor,conditions,k=0.5,norm="rpkm",lc=1,r=10,
    a0per=0.9,nclust=15,filter=1,depth=NULL,cv.cutoff=0,cpm=1) {

    # Check if the factor introduced is already defined
    # If the factor introduced is defined and has more than 2 conditions,
    # it will check if the conditions specified are defined too
    condition_fac <- FALSE
    condition_lev <- FALSE

    datos1 <- datos2 <- matrix()

    for (i in colnames(pData(input))) {
        if (factor == i) {
            condition_fac <- TRUE

            if (!is.factor(pData(input)[,i])) 
                pData(input)[,i] <- as.factor(pData(input)[,i])

            if (length(levels(pData(input)[,i])) == 2) {
                if (!is.null(assayData(input)$exprs)) {
                    datos1 <- assayData(input)$exprs[,which(
                        pData(input)[,i]==levels(pData(input)[,i])[1])]
                    datos2 <- assayData(input)$exprs[,which(
                        pData(input)[,i] ==levels(pData(input)[,i])[2])]
                } 
                else {
                    datos1 <- assayData(input)$counts[,which(
                        pData(input)[,i] ==levels(pData(input)[,i])[1])]
                    datos2 <- assayData(input)$counts[,which(
                        pData(input)[,i] ==levels(pData(input)[,i])[2])]
                }

            # Define the comparison string
            comparison <- paste(levels(pData(input)[,i])[1],
                levels(pData(input)[,i])[2],sep=" - ")

            condition_lev <- TRUE

            if (!((ncol(datos1) > 1) && (ncol(datos2) > 1)))
                stop("Error. NOISeqBIO needs at least 2 biological replicates ",
                    "per condition.\n")
            }
            else {
                if (is.null(conditions))
                    stop("Error. You must specify which conditions you wish ",
                        "to compare when the factor has two or more ",
                        "conditions.\n")
                if (length(conditions) != 2)
                    stop("Error. The argument conditions must contain the 2 ",
                        "conditions you wish to compare.")
            
                l <- conditions %in% pData(input)[,i]
                # If they are defined,they will be TRUE
                if (l[1] == TRUE && l[2] == TRUE) {
                    if (!is.null(assayData(input)$exprs)) {
                        datos1 <- assayData(input)$exprs[,which(
                            pData(input)[,i] == conditions[1])]
                        datos2 <- assayData(input)$exprs[,which(
                            pData(input)[,i] == conditions[2])]
                    }
                    else {
                        datos1 <- assayData(input)$counts[,which(
                            pData(input)[,i] == conditions[1])]
                        datos2 <- assayData(input)$counts[,which(
                            pData(input)[,i] == conditions[2])]
                    }
                    
                    # Define the comparison string
                    comparison <- paste(conditions[1],conditions[2],sep=" - ")
                    condition_lev <- TRUE
                }
            }
        }
    }

    if (condition_fac == FALSE)
        stop("The factor specified does not correspond with any of the ",
            "ones you have defined.") 
    if (condition_lev == FALSE)
        stop("The conditions specified don't exist for the factor specified.\n")

    #--------------------------------------------------------------------------#

    # Number of observations within each condition 
    n1 <- ncol(as.matrix(datos1))
    n2 <- ncol(as.matrix(datos2))
    if (max(n1,n2) == 1) 
        stop("There is only one replicate per condition. Please,use ",
            "NOISeq instead of NOISeqBIO.\n")

    # Rounding off data
    if (norm == "n") {      # no normalization
        datos1 <- round(datos1,10)
        datos2 <- round(datos2,10)
    }

    # Computing k
    if (is.null(k)) {
        m1 <- min(datos1[noceros(datos1,num=FALSE)],na.rm=TRUE)
        m2 <- min(datos2[noceros(datos2,num=FALSE)],na.rm=TRUE)
        k <- min(m1,m2)/2          
    }

    # Total counts for each gene:
    suma1 <- rowSums(as.matrix(datos1))
    suma2 <- rowSums(as.matrix(datos2))

    # Genes with counts in any condition
    concounts <- names(which(suma1+suma2 > 0))

    # All genes
    todos <- rownames(as.matrix(datos1))

    # Gene length
    long <- 1000
    g.sinL <- NULL  # genes with no length defined

    if (!is.null(featureData(input)@data$Length)) {
        g.sinL <- names(which(is.na(featureData(input)@data$Length)))
        if (any(!is.na(featureData(input)@data$Length)) == TRUE) 
            long <- featureData(input)@data[concounts,"Length"]
    }

    # Genes with counts and with length
    gens.sin0 <- setdiff(concounts,g.sinL)

    # cond1 and cond2 in the same matrix
    datitos <- cbind(datos1,datos2)  
    datitos <- datitos[concounts,]  # selecting only genes with counts

    # Sequencing depth when filtering method=3
    if (filter == 3 && is.null(depth)) 
        depth <- colSums(datitos)

    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#

    ## Normalization  
    if (norm == "n") {       # no normalization
        datitos.0 <- sinceros(datitos,k=k)
        datitos.norm <- datitos.0[gens.sin0,]
    }

    if (norm == "rpkm") {      # RPKM
        datitos.0 <- noirpkm(datitos,long=long,k=k,lc=lc)
        datitos.norm <- datitos.0[gens.sin0,]
    }
  
    if (norm == "uqua") {      # Upper Quartile
        datitos.0 <- noiuqua(datitos,long=long,lc=lc,k=k)
        datitos.norm <- datitos.0[gens.sin0,] 
    }

    if (norm == "tmm") {
        datitos.0 <- noitmm(datitos,long=long,lc=lc,k=k)
        datitos.norm <- datitos.0[gens.sin0,]       
    }

    #-------------------------------------------------------------------------#
  
    ## Filtering out low count features
    if (filter != 0) {
        datos.filt <- filtered.data(dataset=datitos.norm,
            factor=c(rep("cond1",n1),rep("cond2",n2)),norm=TRUE,depth=depth,
            method=filter,cv.cutoff=cv.cutoff,cpm=cpm)
    } 
    else {
        datos.filt <- datitos.norm
    }
  
    datos1.filt <- datos.filt[,seq_len(n1)]
    datos2.filt <- datos.filt[,(n1+1):(n1+n2)]   
  
    #-------------------------------------------------------------------------#
  
    ## Noise distribution
  
    Zr <- NULL
  
    if (n1+n2 <= 9) {  # sharing information within clusters
        Zr=share.info(mydata=datos.filt,n1=n1,n2=n2,r=r,nclust=nclust)
    } 
    else {   # r permutations
        for (i in seq_len(r)) {
            print(paste("r =",i))
            mipermu <- sample(seq_len(n1+n2))
            mipermu <- datos.filt[,mipermu]
            mean1 <- rowMeans(mipermu[,seq_len(n1)])
            mean2 <- rowMeans(mipermu[,(n1+1):(n1+n2)])
            sd1 <- apply(mipermu[,seq_len(n1)],1,sd)
            sd2 <- apply(mipermu[,(n1+1):(n1+n2)],1,sd)
            
            myparam <- list("n"=c(n1,n2),"sd"=cbind(sd1,sd2))
            MDperm <- MDbio(dat=cbind(mean1,mean2),param=myparam,a0per=a0per)
            Zr <- cbind(Zr,myDfunction(mydif=MDperm$D,myrat=MDperm$M,stat=1,
                coef=0.5))
        }
    }

    #-------------------------------------------------------------------------#
  
    ## Z-score for different experimental conditions (SIGNAL)

    mean1 <- rowMeans(as.matrix(datos1.filt))
    mean2 <- rowMeans(as.matrix(datos2.filt))        

    sd1 <- apply(as.matrix(datos1.filt),1,sd)
    sd2 <- apply(as.matrix(datos2.filt),1,sd)

    myparam=list("n"=c(n1,n2),"sd"=cbind(sd1,sd2))

    MDs <- MDbio(dat=cbind(mean1,mean2),param=myparam,a0per=a0per)
    
    Zs <- myDfunction(mydif=MDs$D,myrat=MDs$M,stat=1,coef=0.5)
  
    #-------------------------------------------------------------------------#

    ## Completing M and D (in signal)

    lev1 <- mean1[todos]
    lev2 <- mean2[todos]
    names(lev1) <- names(lev2) <- todos

    Zs <- as.numeric(Zs)
    names(Zs) <- rownames(MDs$M)

    Zs <- Zs[todos]
    names(Zs) <- todos

    ## Computing Zn
    Zn <- as.numeric(Zr)

    #-------------------------------------------------------------------------#
  
    ## Results
    list("k"=k,"comp"=comparison,"Level1"=lev1,"Level2"=lev2,"Zs"=Zs,"Zn"=Zn)
}

## Function to summarize difference and ratio information (D and D0)
myDfunction <- function (mydif,myrat,stat,coef) {
    if (stat == 1) { # linear combination of difference and ratio
        myDvalues=coef*mydif + (1-coef)*myrat
    }
    if (stat == 2) { # distance to origin from (ratio,difference)
        myDvalues=sign(mydif) * sqrt((mydif)^2 + (myrat)^2)  
    }
    myDvalues
}

################################################################################

################################################################################

# ARSyNcomponents.R

# This program selects the number of components that explain more than the 
# Variability% If Variability="average" the number of components will be those 
# that explain more than the average variation of the principal components
# For residuals model the number of components selected are 
# beta*average-variability.

ARSyNcomponents <- function(asca=asca,Variability=0.75,beta=2) {
    MODEL<-asca[-length(asca)]
    M<-length(MODEL)-1
    output<-NULL

    if (Variability=="average") {
        for (i in seq_len(M)) { 
            lim <- 1/rankMatrix(MODEL[[i]]$X)[1]    
            t <- table(MODEL[[i]]$var.exp[,1]>lim)
            if (length(t)==1) {
                t[2]=0
            }
            t <- t[2]
            names(t)<-names(MODEL)[i]
            output<-c(output,t)
        }
    }

    if (Variability!="average") {
        lim<-Variability
        for (i in seq_len(M)) { 
            t <- which(MODEL[[i]]$var.exp[,2]>lim)[1]
            names(t) <- names(MODEL)[i]
            output <- c(output,t)
        }
    }
    
    ### Residuals model
    
    i <- M+1
    lim <- beta*1/rankMatrix(MODEL[[i]]$X)[1]   
    t <- table(MODEL[[i]]$var.exp[,1]>lim)
    if (length(t)==1) {
        t[2]=0
    }
    t <- t[2]
    names(t)<-names(MODEL)[i]
    output <- c(output,t)

    output
}

################################################################################

################################################################################

# ARSyNmodel.R

ARSyNmodel <- function(X=X,Factors=2,Designa=Designa,Designb=Designb,
    Designc=Designc,Variability="average",Join=TRUE,Interaction=TRUE,beta=2) {

    if (Factors==1) {
        Fac0 <- c(1,2)
        names(Fac0) <- c("Model.a","Model.res")
        asca0 <- ASCA.1f(X=X,Designa=Designa,Fac=Fac0)
        Fac <- ARSyNcomponents(asca0,Variability=Variability,beta=beta)
        for (i in seq_len(length(Fac))) {
            Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
        }
        asca <- ASCA.1f(X=X,Designa=Designa,Fac=Fac0)
    }



    if (Factors==2) {
        Fac0 <- c(1,2,2,2)
        names(Fac0) <- c("Model.a","Model.b","Model.ab","Model.res")
        if (Join[1]) {
            names(Fac0)[3] <- c("Model.bab")
        }
        asca0 <- ASCA.2f(X=X,Designa=Designa,Designb=Designb,Fac=Fac0,Join=Join,
            Interaction=Interaction)
        Fac <- ARSyNcomponents(asca0,Variability=Variability,beta=beta)
        for (i in seq_len(length(Fac))) {
            Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
        }
        asca <- ASCA.2f(X=X,Designa=Designa,Designb=Designb,Fac=Fac0,Join=Join,
            Interaction=Interaction)
    }

    if (Factors==3) {
        Fac0 <- c(0,2,2,2,2,2,2,2)
        names(Fac0) <- c("Model.a","Model.b","Model.c","Model.ab","Model.ac",
            "Model.bc","Model.abc","Model.res")
        if (Join[1]) {
            names(Fac0)[4]<-c("Model.bab")
        }
        if (Join[2]) {
            names(Fac0)[5]<-c("Model.cac")
        }
        asca0 <- ASCA.3f(X=X,Designa=Designa,Designb=Designb,Designc=Designc,
            Fac=Fac0,Join=Join,Interaction=Interaction)
        Fac<-ARSyNcomponents(asca0,Variability=Variability,beta=beta)
        for (i in seq_len(length(Fac))) {
            Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
        }
        asca <- ASCA.3f(X=X,Designa=Designa,Designb=Designb,Designc=Designc,
            Fac=Fac0,Join=Join,Interaction=Interaction)
    }

    Input <- asca$Input
    asca$Input <- c(Input,Factors)
    names(asca$Input) <- c(names(Input),"Factors")

    output <- asca
    output
}

################################################################################

################################################################################

# ARSyNseq.R

ARSyNseq <- function(data,factor=NULL,batch=FALSE,norm="rpkm",logtransf=FALSE,
    Variability=0.75,beta=2) {

    Join <- TRUE
    Interaction <- TRUE

    dat <- as.matrix(assayData(data)$exprs)
    long <- featureData(data)@data[rownames(dat),"Length"]
    if (is.null(long)) 
        long <- 1000

    if (norm == "rpkm") {
        dat <- noirpkm(dat,long=long,k=0,lc=1)
    }

    if (norm == "uqua") {
        dat <- noiuqua(dat,long=long,lc=1,k=0)
    }

    if (norm == "tmm") {
        dat <- noitmm(dat,long=long,lc=1,k=0)      
    }

    if (!logtransf)   dat <- log(dat + 1)

    X <- t(dat)         #conditions x genes

    #----------------------------------

    if (is.null(factor)) {
        Covariates <- t(pData(data))
        Num.factors <- nrow(Covariates)
        labels.factors <- rownames(Covariates)
        Design <- list(NULL,NULL,NULL)
        for (i in seq_len(Num.factors)) {
            x <- as.character(Covariates[i,])
            Design[[i]] <- make.ASCA.design(x)
        } 
    }
    else {
        Covariates <- pData(data)[,factor]
        Num.factors <- 1
        labels.factors <- factor
        Design <- list(NULL,NULL,NULL)
        x <- as.character(Covariates)
        Design[[1]] <- make.ASCA.design(x)
    }

    ####################################
    ### --- Execute ASCAmodel 
    ####################################

    my.asca <- ARSyNmodel(Factors=Num.factors,X=X,Designa=Design[[1]],
        Designb=Design[[2]],Designc=Design[[3]],Join=Join,
        Interaction=Interaction,Variability=Variability,beta=beta)

    #################################### 
    ### --- Writing filtered matrix 
    ####################################

    X.filtered <- X
    M <- length(my.asca)-1

    if (!batch) {
        X.filtered <- X.filtered-my.asca[[M]]$TP
    }

    if(batch) {
        X.filtered <- X.filtered-my.asca[[1]]$TP 
    }

    data.filtered <- t(X.filtered)

    if (!logtransf)
        data.filtered <- exp(data.filtered)+1
    
    exprs(data) <- data.filtered

    return(data)
}

################################################################################

# ASCA1f.R

ASCA.1f <- function(X=X,Designa=Designa,Designb=NULL,Designc=NULL,Fac=c(1,2),
    Join=NULL,Interaction=NULL) {

    n <- ncol(X)
    p <- nrow(X)
    I <- ncol(Designa)

    Faca <- Fac[1] # number components Model a (time)
    Facres <- Fac[2] # number components Residues

    #----------------------- Calculate Overall Mean ----------------------------

    offset <- apply(X,2,mean)
    Xoff <- X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))

    #-----------------------  PART I: Submodel a (TIME) ------------------------

    Model.a<-ASCAfun1(Xoff,Designa,Faca)
    Xres<-Xoff-Model.a$X

    #------------------------Collecting models ---------------------------------

    models <- ls(pattern="Model")
    output <- vector(mode="list")
    Xres <- Xoff
    for (i in seq_len(length(models))) {
        mymodel <- get(models[i], envir=environment())
        output <- c(output, list(mymodel))
        Xres <- Xres - mymodel$X
        rm(mymodel)
        gc()
    }
    names(output) <- models

    #------------------------- PART III: Submodel res --------------------------

    Model.res <- ASCAfun.res(Xres,Facres)

    Model.res <- list(Model.res)
    names(Model.res) <- c("Model.res")
    output <- c(output,Model.res)

    #------------------------- Add Input data to the Output --------------------

    Input <- list(X,Designa,Designb,Designc,Fac,Join,Interaction)
    names(Input)<-c("X", "Designa", "Designb", "Designc", "Fac", "Join",
        "Interaction")
    Input <- list(Input)
    names(Input) <- "Input"
    output <- c(output,Input)
    
    output
}

################################################################################

################################################################################

# ASCA2f.R

ASCA.2f <- function(X=X,Designa=Designa,Designb=Designb,Designc=NULL,
    Fac=c(1,2,2,2),Join = TRUE,Interaction=TRUE) {

    n <- ncol(X)
    p <- nrow(X)
    I <- ncol(Designa)
    J <- ncol(Designb)

    Faca <- Fac[1] # number components Model a (time)
    Facb <- Fac[2] # number components Model b  (second factor)
    Facab <- Fac[3] # number components Model ab (interaction) 
    Facres <- Fac[4] # number components Residues

    #----------------------- Calculate Overall Mean ----------------------------

    offset <- apply(X,2,mean)
    Xoff <- X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))

    #-----------------------  PART I: Submodel a (TIME) ------------------------

    Model.a <- ASCAfun1(Xoff,Designa,Faca)
    Xres <- Xoff-Model.a$X

    #-------------------------- PART II: Submodel b and ab ---------------------

    if (!Join) {
        Model.b<-ASCAfun1(Xoff,Designb,Facb)
        if (Interaction) {
            Model.ab<-ASCAfun2(Xoff,Designa,Designb,Facab)
        }
    }

    if (Join) {
        Model.bab <- ASCAfun12(Xoff,Designa,Designb,Facab)
    }

    #------------------------Collecting models ---------------------------------

    models <- ls(pattern="Model")
    output <- vector(mode="list")
    Xres <- Xoff
    for (i in seq_len(length(models))) {
        mymodel <- get(models[i],envir=environment())
        output <- c(output,list(mymodel))
        Xres <- Xres - mymodel$X
        rm(mymodel)
        gc()    
    }

    names(output) <- models

    #------------------------- PART III: Submodel res --------------------------

    Model.res <- ASCAfun.res(Xres,Facres)

    Model.res <- list(Model.res)
    names(Model.res) <- c("Model.res")
    output <- c(output,Model.res)

    #------------------------- Add Input data to the Output --------------------

    Input <- list(X,Designa,Designb,Designc,Fac,Join,Interaction,Xoff)
    names(Input) <- c("X","Designa","Designb","Designc","Fac","Join",
        "Interaction","Xoff")
    Input <- list(Input)
    names(Input) <- "Input"
    output <- c(output,Input)

    output
}

################################################################################

################################################################################

# ASCA3f.R

ASCA.3f <- function(X=X,Designa=Designa,Designb=Designb,Designc=Designc,
    Fac=c(1,2,2,2,2,2,2,2),Join=c(TRUE,TRUE),
    Interaction=c(TRUE,TRUE,TRUE,TRUE)) {

    n <- ncol(X)
    p <- nrow(X)
    I <- ncol(Designa)
    J <- ncol(Designb)
    K <- ncol(Designc)

    Faca <- Fac[1]# number components Model a (time)
    Facb <- Fac[2] # number components Model b  (second factor)
    Facc <- Fac[3] # number components Model c  (third factor)
    Facab <- Fac[4] # number components Model ab (interaction) 
    Facac <- Fac[5] # number components Model ac (interaction) 
    Facbc <- Fac[6] # number components Model bc (interaction) 
    Facabc <- Fac[7]  # number components Model abc (interaction) 
    Facres <- Fac[8]
    
    #----------------------- Calculate Overall Mean ----------------------------

    offset <- apply(X,2,mean)
    Xoff <- X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))

    #-----------------------  PART I: Submodel a (TIME) ------------------------

    Model.a <- ASCAfun1(Xoff,Designa,Faca)
    Xres <- Xoff-Model.a$X

    #-------------------------- PART II.1: Submodel b.ab------------------------
    if(!Join[1]) {
        Model.b <- ASCAfun1(Xoff,Designb,Facb)
        if (Interaction[1]) {
            Model.ab <- ASCAfun2(Xoff,Designa,Designb,Facab)
        }
    }

    if(Join[1]) {
        Model.bab<-ASCAfun12(Xoff,Designa,Designb,Facab)
    }
    
    #-------------------------- PART II.2: Submodel (c.ac) ---------------------
    if(!Join[2]) {
        Model.c <- ASCAfun1(Xoff,Designc,Facc)
        if (Interaction[2]) {
            Model.ac<-ASCAfun2(Xoff,Designa,Designc,Facac)
        }
    }

    if(Join[2]) {
        Model.cac<-ASCAfun12(Xoff,Designa,Designc,Facac)
    }
    
    #-------------------------- PART II.3: Submodel (bc) -----------------------

    if (Interaction[3]) {
        Model.bc<-ASCAfun2(Xoff,Designb,Designc,Facbc)
    }
    
    #-------------------------- PART II.4: Submodel (abc) ----------------------

    if (Interaction[4]) {
        Model.abc<-ASCAfun.triple(Xoff,Designa,Designb,Designc,Facabc)
    }
    #------------------------Collecting models ---------------------------------

    models <- ls(pattern="Model")
    output <- vector(mode="list")
    Xres <- Xoff
    for (i in seq_len(length(models))) {
        mymodel <- get(models[i], envir=environment())
        output <- c(output, list(mymodel))
        Xres <- Xres - mymodel$X
        rm(mymodel)
        gc()
    }
    names(output) <- models

    #------------------------- PART III: Submodel res --------------------------

    Model.res <- ASCAfun.res(Xres,Facres)

    Model.res <- list(Model.res)
    names(Model.res) <- c("Model.res")
    output <- c(output,Model.res)

    #------------------------- Add Input data to the Output --------------------

    Input <- list(X, Designa, Designb, Designc, Fac, Join,Interaction,Xoff)
    names(Input) <- c("X","Designa","Designb","Designc","Fac","Join",
        "Interaction","Xoff")
    Input <- list(Input)
    names(Input) <- "Input"
    output<-c(output,Input)

    output
}

################################################################################

################################################################################

# ASCAfun1.R

ASCAfun1 <- function(X,Design,Fac) {
    n <- ncol(X) # number of genes
    I <- ncol(Design) # number of levels in the factor

    NK <- NULL
    XK <- matrix(NA,nrow=I,ncol=n)

    for (i in seq_len(I)) {
        sub <- X[Design[,i]==1,]
        if (is.null(nrow(sub))) { #when there isn't replicates
            NK[i] <- 1
            XK[i,] <- sub 
        }
        else {
        NK[i]<-nrow(sub)
        XK[i,]<-apply(sub,2,mean) }
    }
    NK<-sqrt(NK)

    # Weigh the data of the Submodel with the corresponding number of 
    # measurement occasions

    XKw <- NK*XK

    PCA <- PCA.GENES(XKw)
    scw <- PCA$scores[,seq_len(Fac)]
    ld <- PCA$loadings[,seq_len(Fac)]
    ssq <- PCA$var.exp
    if (Fac==1) {
        scw <- as.matrix(scw)
        ld <- as.matrix(ld)
    }      
    if (Fac==0) {
        scw<-as.matrix(rep(0,I))
        ld<-as.matrix(rep(0,n))
    }        
    
    # Re-weigth the scores
    sc <- scw/NK
    XKrec <- sc%*%t(ld)

    Xa <- NULL
    TPa <- NULL
    for (i in seq_len(nrow(X))) {
        position<-which(Design[i,]==1)
        Xa<-rbind(Xa,XK[position,])
        TPa<-rbind(TPa,XKrec[position,])
    }

    Ea <- Xa-TPa

    #Leverage & SPE 
    leverage<-apply(ld^2,1,sum)
    SPE<-apply(Ea^2,2,sum)

    output <- list(XK,sc,ld,ssq,Xa,TPa,Ea,leverage,SPE)
    names(output) <- c("data","scores","loadings","var.exp","X","TP","E",
        "leverage","SPE")
    
    output
}

################################################################################

################################################################################

# ASCAfun2.R

ASCAfun2 <- function(X,Desa,Desb,Fac) {
    n <- ncol(X) # number of genes
    I <- ncol(Desa) # number of levels in the factor TIME
    J <- ncol(Desb) # number of levels in the other factor

    XK1 <- matrix(NA,nrow=I,ncol=n)

    for (i in seq_len(I)) {
        sub <- X[Desa[,i]==1,]
        if (is.null(nrow(sub))) { #when there isn't replicates
            XK1[i,] <- sub 
        }
        else {
            XK1[i,] <- apply(sub,2,mean) 
        }
    }

    XK2 <- matrix(NA,nrow=J,ncol=n)

    for (j in seq_len(J)) {
        sub <- X[Desb[,j]==1,]
        if (is.null(nrow(sub))) { #when there isn't replicates
            XK2[j,] <- sub 
        }
        else {
            XK2[j,]<-apply(sub,2,mean) 
        }
    }

    NK <- matrix(NA,nrow=I,ncol=J)
    XK <- matrix(NA,nrow=I*J,ncol=n)

    k=1
    for (j in seq_len(J)){
        for (i in seq_len(I)) {
            sub <- X[(Desa[,i]+Desb[,j])==2,]
            if (is.null(nrow(sub))) {  #when there isn't replicates
                NK[i,j] <- 1
                XK[k,] <- sub-XK1[i,]-XK2[j,]
            }
            else {
                NK[i,j] <- sqrt(nrow(sub))
                XK[k,] <- apply(sub,2,mean)-XK1[i,]-XK2[j,] 
            }
            k <- k+1
        }
    }

    XKw <- XK*(as.numeric(NK))

    PCA <- PCA.GENES(XKw)
    scw <- PCA$scores[,seq_len(Fac)]
    ld <- PCA$loadings[,seq_len(Fac)]
    ssq <- PCA$var.exp
    if (Fac==1) {
        scw <- as.matrix(scw)
        ld <- as.matrix(ld) 
    }
    if (Fac==0) {
        scw <- as.matrix(rep(0,I*J))
        ld <- as.matrix(rep(0,n))
    }
    
    # Re-weigth the scores
    sc <- scw/(as.numeric(NK))
    
    XKrec <- sc%*%t(ld)

    Xab <- NULL
    TPab <- NULL

    for (i in seq_len(nrow(X))) {
        position1 <- which(Desa[i,]==1)
        position2 <- which(Desb[i,]==1)
        Xab <- rbind(Xab,XK[I*(position2-1)+position1,])
        TPab <- rbind(TPab,XKrec[I*(position2-1)+position1,])
    }
    
    Eab <- Xab-TPab

    leverage <- apply(ld^2,1,sum)
    SPE <- apply(Eab^2,2,sum)

    output <- list(XK,sc,ld,ssq,Xab,TPab,Eab,leverage,SPE)
    names(output)<-c("data","scores","loadings","var.exp","X","TP","E",
        "leverage","SPE")
    
    output
}

################################################################################

################################################################################

# ASCAfun12.R

ASCAfun12 <- function (X,Desa,Desb,Fac) {
    n <- ncol(X) # number of genes
    I <- ncol(Desa) # number of levels in the factor TIME
    J <- ncol(Desb) # number of levels in the other factor

    XK1 <- matrix(NA,nrow=I,ncol=n)

    for (i in seq_len(I)) {
         sub <- X[Desa[,i]==1,]
         XK1[i,] <- apply(sub,2,mean)
    }
    
    NK <- matrix(NA,nrow=I,ncol=J)
    XK <- matrix(NA,nrow=I*J,ncol=n)

    k <- 1
    for (j in seq_len(J)) {
        for (i in seq_len(I)) {
            sub <- X[(Desa[,i]+Desb[,j])==2,]
            if (is.null(nrow(sub))) {  #when there isn't replicates
                NK[i,j] <- 1
                XK[k,] <- sub-XK1[i,]
            }
            else {
                NK[i,j] <- sqrt(nrow(sub))
                XK[k,] <- apply(sub,2,mean)-XK1[i,] 
            }
            k <- k+1
        }
    }

    XKw <- XK*(as.numeric(NK))

    PCA <- PCA.GENES(XKw)
    scw <- PCA$scores[,seq_len(Fac)]
    ld <- PCA$loadings[,seq_len(Fac)]
    ssq <- PCA$var.exp
      
    if (Fac==1) {
        scw <- as.matrix(scw)
        ld <- as.matrix(ld) 
    }
    if (Fac==0) {
        scw <- as.matrix(rep(0,I*J))
        ld <- as.matrix(rep(0,n))
    }
    
    # Re-weigth the scores
    sc <- scw/(as.numeric(NK))
    
    XKrec <- sc%*%t(ld)

    Xab <- NULL
    TPab <- NULL
    for (i in seq_len(nrow(X))){
        position1 <- which(Desa[i,]==1)
        position2 <- which(Desb[i,]==1)
        Xab <- rbind(Xab,XK[I*(position2-1)+position1,])
        TPab <- rbind(TPab,XKrec[I*(position2-1)+position1,])
    }
    Eab <- Xab-TPab

    #Leverage & SPE 
    leverage <- apply(ld^2,1,sum)
    SPE <- apply(Eab^2,2,sum)
    output <- list(XK,sc,ld,ssq,Xab,TPab,Eab,leverage,SPE)
    names(output) <- c("data","scores","loadings","var.exp","X","TP","E",
        "leverage","SPE")
    
    output
}

################################################################################

################################################################################

# ASCAfunres.R

ASCAfun.res <- function (X,Fac) {
    PCA <- PCA.GENES(X)
    sc <- PCA$scores[,seq_len(Fac)]
    ld <- PCA$loadings[,seq_len(Fac)]
    ssq <- PCA$var.exp
    if (Fac==1) {
        sc <- as.matrix(sc)
        ld <- as.matrix(ld)
    }
    TPres<-sc%*%t(ld)
    if (Fac==0) {
        sc <- 0
        ld <- 0
        TPres<-matrix(0,nrow(X),ncol(X))
    }      
    Eres<-X-TPres
    output <- list(sc,ld,ssq,X,TPres,Eres)
    names(output) <- c("scores","loadings","var.exp","X","TP","E")
    output
}

################################################################################

################################################################################

# ASCAfun-triple.R

ASCAfun.triple <- function(X,Desa,Desb,Desc,Fac) {

    n <- ncol(X) # number of genes
    I <- ncol(Desa) # number of levels in the factor TIME
    J <- ncol(Desb) # number of levels in the other factor
    H <- ncol(Desc) # number of levels in the other factor

    #Matrices con medias efectos individuales
    XK1 <- matrix(NA,nrow=I,ncol=n)
    for (i in seq_len(I)) {
        sub <- X[Desa[,i]==1,]
        XK1[i,] <- apply(sub,2,mean)
    }

    XK2 <- matrix(NA,nrow=J,ncol=n)
    for (j in seq_len(J)) {
        sub <- X[Desb[,j]==1,]
        XK2[j,] <- apply(sub,2,mean)
    }

    XK3 <- matrix(NA,nrow=H,ncol=n)
    for (h in seq_len(H)) {
        sub <- X[Desc[,h]==1,]
        XK3[h,] <- apply(sub,2,mean)
    }

    #Matrices con medias de efectos simples

    XK12 <- matrix(NA,nrow=I*J,ncol=n)
    k <- 1
    for (j in seq_len(J)) {
        for (i in seq_len(I)) {
            sub <- X[(Desa[,i]+Desb[,j])==2,]
            XK12[k,] <- apply(sub,2,mean)
            k <- k+1
        }
    }

    XK13 <- matrix(NA,nrow=I*H,ncol=n)
    k <- 1
    for (h in seq_len(H)) {
        for (i in seq_len(I)) {
            sub <- X[(Desa[,i]+Desc[,h])==2,]
            XK13[k,] <- apply(sub,2,mean)
            k <- k+1
        }
    }

    XK23 <- matrix(NA,nrow=J*H,ncol=n)
    k <- 1
    for (h in seq_len(H)){
        for (j in seq_len(J)){
            sub <- X[(Desb[,j]+Desc[,h])==2,]
            XK23[k,] <- apply(sub,2,mean)
            k <- k+1
        }
    }

    NK <- matrix(NA,nrow=I,ncol=J*H)
    XK <- matrix(NA,nrow=I*J*H,ncol=n)
    k <- 1

    for (h in seq_len(H)) {
        for (j in seq_len(J)) {
            for (i in seq_len(I)) {
                sub <- as.matrix(rbind(X[(Desa[,i]+Desb[,j]+Desc[,h])==3,]))
                NK[i,(h-1)*J+j] <- sqrt(nrow(sub))
                XK[k,] <- apply(sub,2,mean)+XK1[i,]+
                    XK2[j,]+XK3[h,]-XK12[(j-1)*I+i,]-XK13[(h-1)*I+i,]-
                    XK23[(h-1)*J+j,]
                k <- k+1
            }
        }
    }

    XKw <- XK*(as.numeric(NK))

    PCA <- PCA.GENES(XKw)
    scw <- PCA$scores[,seq_len(Fac)]
    ld <- PCA$loadings[,seq_len(Fac)]
    ssq <- PCA$var.exp
    if (Fac==1) {
        scw <- as.matrix(scw)
        ld <- as.matrix(ld) 
    }
    if (Fac==0) {
        scw <- as.matrix(rep(0,I*J*H))
        ld <- as.matrix(rep(0,n))
    }
    # Re-weigth the scores
    sc <- scw/(as.numeric(NK))
    
    XKrec <- sc%*%t(ld)

    Xabc <- NULL
    TPabc <- NULL

    for (i in seq_len(nrow(X))){
        position1 <- which(Desa[i,]==1)
        position2 <- which(Desb[i,]==1)
        position3 <- which(Desc[i,]==1)
        Xabc <- rbind(Xabc,XK[I*(position2-1)+I*J*(position3-1)+position1,])
        TPabc <- rbind(TPabc,XKrec[I*(position2-1)+
            I*J*(position3-1)+position1,])
    }
    Eabc <- Xabc-TPabc

    #Leverage & SPE 
    leverage <- apply(ld^2,1,sum)
    SPE <- apply(Eabc^2,2,sum)

    output <- list(XK,sc,ld,ssq,Xabc,TPabc,Eabc,leverage,SPE)
    names(output) <- c("data","scores","loadings","var.exp","X","TP","E",
        "leverage","SPE")
    
    output
}

################################################################################

################################################################################

# auxiliar.R

busca <- function (x,S) {
    which(S[,1] == x[1] & S[,2] == x[2])
}

int.mult <- function(lista,todos=NULL) {
    if (is.null(todos)) {
        todos <- unlist(lista)
    }
    comunes <- todos
    for (i in seq_len(length(lista))) {
        comunes <- intersect(comunes,lista[[i]])
    }
    comunes
}

n.menor <- function (x,S1,S2) {
    length(which(S1 <= x[1] &  S2 <= x[2]))
}

noceros <- function (x,num=TRUE,k=0) {
    nn <- length(which(x > k))
    if (num) {
        nn 
    } 
    else {
        if(nn > 0) { 
            which(x > k) 
        } 
        else { 
            NULL 
        }
    }
}

sinceros <- function (datos,k) {
    datos <- as.matrix(datos)
    datos0 <- as.matrix(datos)
    
    if (is.null(k)) {
        mini0 <- min(datos[noceros(datos,num=FALSE,k=0)])
        kc <- mini0/2
        datos0[datos0 == 0] <- kc
    } 
    else {
        datos0[datos0 == 0] <- k
    }
    datos0
}

sim.samples <- function(counts1,counts2=NULL,pnr=1,nss=5,v=0.02) {
    seqdep <- c(sum(counts1),sum(counts2))
    num.reads1 <- (pnr + c(-v,v))*seqdep[1]   
    muestras <- vector("list")
    muestras$c1 <- NULL
    for (s in seq_len(nss)) {
        tama <- round(runif(1,num.reads1[1],num.reads1[2]),0)
        muestras$c1 <- cbind(muestras$c1,rmultinom(1,size=tama,prob=counts1))
    }
    
    if (!is.null(counts2)) {
        num.reads2 <- (pnr + c(-v,v))*seqdep[2]
        muestras$c2 <- NULL
        for (s in seq_len(nss)) {
            tama <- round(runif(1,num.reads2[1],num.reads2[2]),0)
            muestras$c2 <- cbind(muestras$c2,
                rmultinom(1,size=tama,prob=counts2))
        }
    }
    muestras
}

ranking <- function(results) {
    M <- results$M
    D <- results$D
    prob <- results$prob
    
    ## Changing NA by 0
    M[is.na(M)] <- 0
    D[is.na(D)] <- 0
    prob[is.na(prob)] <- 0
    
    
    ## Ranking
    ranking5 <- sqrt(M*M + D*D)*sign(M)
    
    ## Ranking results
    theranking <- data.frame(rownames(results),ranking5)
    rownames(theranking) <- NULL
    colnames(theranking) <- c("ID","statistic")
    
    theranking
}

plot.y2 <- function(x,yright,yleft,yrightlim=range(yright,na.rm=TRUE),
    yleftlim=range(yleft,na.rm=TRUE),xlim=range(x,na.rm=TRUE),xlab=NULL,
    yylab=c("",""),lwd=c(2,2),pch=c(1,2),col=c(1,2),type=c("o","o"),
    linky=TRUE,smooth=0,bg=c("white","white"),lwds=1,length=10,...,x2=NULL,
    yright2=NULL,yleft2=NULL,col2=c(3,4)) {
  
    ## Plotting RIGHT axis data
    plot(x,yright,axes=FALSE,ylab="",xlab=xlab,ylim=yrightlim,xlim=xlim,
        pch=pch[1],type=type[1],lwd=lwd[1],col=col[1],...)
  
    axis(4,pretty(yrightlim,length),col=1,col.axis=1)

    if (is.null(yright2) == FALSE) {
        points(x2,yright2,type=type[1],pch=pch[1],lwd=lwd[1],col=col2[1],...)
    }
  
    if (smooth != 0) 
        lines(supsmu(x,yright,span=smooth),col=col[1],lwd=lwds,...)
  
    if (yylab[1]=="") {
        mtext(deparse(substitute(yright)),side=4,outer=FALSE,line=2,col=1,...)
    } 
    else {
        mtext(yylab[1],side=4,outer=FALSE,line=2,col=1,...)
    }

    par(new=TRUE)

    ## Plotting LEFT axis data
    plot(x,yleft,axes=FALSE,ylab="" ,xlab=xlab,ylim=yleftlim,xlim=xlim,bg=bg[1],
        pch=pch[2],type=type[2],lwd=lwd[2],col=col[2],...)
    box()
    axis(2,pretty(yleftlim,length),col=1,col.axis=1)

    if (is.null(yleft2) == FALSE) {
        points(x2,yleft2,type=type[2],pch=pch[2],bg=bg[2],lwd=lwd[2],
            col=col2[2],...)
    }
  
    if (smooth != 0) 
        lines(supsmu(x,yleft,span=smooth),col=col[2],lwd=lwds,...)
  
    if(yylab[2] == "") {
        mtext(deparse(substitute(yleft)),side=2,outer=FALSE,line=2,col=1,...)
    } 
    else {
        mtext(yylab[2],side=2,outer=FALSE,line=2,col=1,...)
    }

    ## X-axis
    axis(1,at=pretty(xlim,length))
}

logscaling <- function(data,base=2,k=1) {
    logmaximo <- round(max(data,na.rm=TRUE),0)
    numceros <- nchar(logmaximo)-1
    etiquetas <- c(0,10^(seq_len(numceros)))
    donde <- log(etiquetas + k,base=base)
    data <- log(data + k,base=base)
    list("data"=data,"at"=donde,"labels"=etiquetas)
}

miscolores <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]

################################################################################

################################################################################

# biodetection.plot.R

biodetection.dat <- function(input,factor=NULL,k=0) {
    if (inherits(input,"eSet") == FALSE)
        stop("Error. The input data must be an eSet object\n")

    if (any(!is.na(featureData(input)@data$Biotype)) == FALSE)
        stop ("No biological classification was provided.\nPlease run ",
            "addData() function to add this information\n")

    if (!is.null(assayData(input)$exprs)) {
        dat <- as.matrix(assayData(input)$exprs)
        mysamples <- colnames(assayData(input)$exprs)
    } 
    else {
        dat <- as.matrix(assayData(input)$counts)
        mysamples <- colnames(assayData(input)$counts)
    }
    
    numgenes <- nrow(dat)

    if (is.null(factor)) {  # per sample  
        cat("Biotypes detection is to be computed for:\n")
        print(colnames(dat))
        biotablas <- vector("list",length=NCOL(dat))
        names(biotablas) <- colnames(dat)
    } 
    else {  # per condition
        mifactor <- pData(input)[,factor]
        niveles <- levels(mifactor)
        cat("Biotypes detection is to be computed for:\n")
        print(niveles)
        biotablas <- vector("list",length=length(niveles))
        names(biotablas) <- niveles
        dat <- vapply(niveles,
            function (k) rowSums(as.matrix(dat[,grep(k,mifactor)])),numeric(1))
    }

    infobio <- as.character(featureData(input)@data$Biotype)

    genome <- 100*table(infobio)/sum(table(infobio))
    ordre <- order(genome,decreasing=TRUE)

    for (i in seq_len(length(biotablas))) {
        detect <- dat[,i] > k
        perdet1 <- genome*table(infobio,detect)[names(genome),"TRUE"]/
            table(infobio)[names(genome)]
        perdet2 <- 100*table(infobio,detect)[names(genome),"TRUE"] /
            sum(table(infobio,detect)[,"TRUE"])
        biotablas[[i]] <- as.matrix(rbind(perdet1[ordre],perdet2[ordre]))
        rownames(biotablas[[i]]) <- c("detectionVSgenome","detectionVSsample")
    }

    mybiotable <- list("genome"=genome[ordre],"biotables"=biotablas,
        "genomesize"=numgenes)

    mybiotable
}

biodetection.plot <- function(dat,samples=c(1,2),
    plottype=c("persample","comparison"),toplot="protein_coding",
    toreport=FALSE,...) {

    mypar=par(no.readonly=TRUE)
    plottype=match.arg(plottype)

    if (length(samples) > 2) {
        stop("ERROR: This function cannot generate plots for more than 2 ",
            "samples.\nPlease,use it as many times as needed to generate the ",
            "plots for all your samples.\n")
    }
    
    if (is.numeric(samples)) 
        samples=names(dat$biotables)[samples]
    biotable1 <- rbind(dat$genome,dat$biotables[[samples[1]]],
        rep(0,length(dat$genome)))

    # Computing ylim for left and right axis         
    if (ncol(biotable1) >= 3) {
        ymaxL <- ceiling(max(biotable1[,c(1,2,3)],na.rm=TRUE))
        ymaxR <- max(biotable1[,-c(1,2,3)],na.rm=TRUE)                
    } 
    else {
        ymaxL <- ceiling(max(biotable1,na.rm=TRUE))
        ymaxR=0 
    }

    if (length(samples) == 2) {
        biotable2 <- rbind(dat$genome,dat$biotables[[samples[2]]],
            rep(0,length(dat$genome)))
        if (ncol(biotable2) >= 3) {
            ymax2 <- ceiling(max(biotable2[,c(1,2,3)],na.rm=TRUE))
            ymax2sin <- max(biotable2[,-c(1,2,3)],na.rm=TRUE)
            ymaxR <- ceiling(max(ymaxR,ymax2sin))
        } 
        else {
            ymax2 <- ceiling(max(biotable2,na.rm=TRUE))            
        }
        ymaxL=max(ymaxL,ymax2)     
    }

    # Rescaling biotables (datos2)
    if (length(samples) == 2) {                
        if (ncol(biotable2) >= 3) 
            biotable2[,-c(1,2,3)] <- biotable2[,-c(1,2,3)]*ymaxL/ymaxR
    }

    # Rescaling biotables (datos1)    
    if (ncol(biotable1) >= 3) 
        biotable1[,-c(1,2,3)] <- biotable1[,-c(1,2,3)]*ymaxL/ymaxR

    ## PLOTS

    if (length(samples) == 1) {     # Plot (1 sample) - 2 scales                
        par(mar=c(11,4,2,2))
        barplot(biotable1[c(1,3),],main=samples[1],xlab=NULL,ylab="%features",
            axis.lty=1,legend=FALSE,beside=TRUE,col=c("grey",2),las=2,
            ylim=c(0,ymaxL),border=c("grey",2))
        barplot(biotable1[c(2,4),],main=samples[1],xlab=NULL,ylab="%features",
            axis.lty=1,legend=FALSE,beside=TRUE,col=c(2,1),las=2,density=30,
            ylim=c(0,ymaxL),border=2,add=TRUE)

        # if number of biotypes >= 3 so we have left and right axis
        if (ymaxR > 0) {
            axis(side=4,at=pretty(c(0,ymaxL),n=5),
                labels=round(pretty(c(0,ymaxL),n=5)*ymaxR/ymaxL,1))
            abline(v=9.5,col=3,lwd=2,lty=2)                
        }        
        
        legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",2,2),
            density=c(NA,30,NA),border=c("grey",2,2),
            legend=c("% in genome","detected","% in sample"))
    } 
    else {     # Plot (2 samples)
        par(mar=c(11,4,2,2))
        if (plottype == "persample") {
            barplot(biotable1[c(1,3),],main=samples[1],xlab=NULL,
                ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c("grey",2),las=2,ylim=c(0,ymaxL),border=c("grey",2))
            barplot(biotable1[c(2,4),],main=samples[1],xlab=NULL,
                ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c(2,1),las=2,density=30,ylim=c(0,ymaxL),border=2,add=TRUE)
            
            # if number of biotypes >= 3 so we have left and right axis
            if (ymaxR > 0) {
                axis(side=4,at=pretty(c(0,ymaxL),n=5),
                    labels=round(pretty(c(0,ymaxL),n=5)*ymaxR/ymaxL,1))
                abline(v=9.5,col=3,lwd=2,lty=2)
            }
            
            legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",2,2),
                density=c(NA,30,NA),border=c("grey",2,2),
                legend=c("% in genome","detected","% in sample"))

            # Datos2                
            barplot(biotable2[c(1,3),],main=samples[2],xlab=NULL,
                ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c("grey",4),las=2,ylim=c(0,ymaxL),border=c("grey",4))
            barplot(biotable2[c(2,4),],main=samples[2],xlab=NULL,
                ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,col=c(4,1),
                las=2,density=30,ylim=c(0,ymaxL),border=4,add=TRUE)
            
            # if number of biotypes >= 3 so we have left and right axis
            if (ymaxR > 0) {
                axis(side=4,at=pretty(c(0,ymaxL),n=5),
                    labels=round(pretty(c(0,ymaxL),n=5)*ymaxR/ymaxL,1))
                abline(v=9.5,col=3,lwd=2,lty=2)
            }

            legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",4,4),
                density=c(NA,30,NA),border=c("grey",4,4),
                legend=c("% in genome","detected","% in sample"))
        }
        
        # A plot comparing two samples with regard to genome and for % in sample
        if (plottype == "comparison") {
            lefttable <- rbind(100*dat$biotables[[samples[1]]][1,]/dat$genome,
                100*dat$biotables[[samples[2]]][1,]/dat$genome)
            righttable <- rbind(dat$biotables[[samples[1]]][2,],
                dat$biotables[[samples[2]]][2,])

            if (length(toplot) > 1) {
                toplot <- toplot[1]
                print("WARNING: More than one biotype was provided, the ",
                    "proportion test will only by applied to the first ",
                    "biotype.")
            }
            
            if ((toplot != 1) && (toplot != "global")) {
                numgenes <- dat$genomesize
                myx <- round(righttable[,toplot]*numgenes/100,0)
                mytest <- prop.test(x=myx,n=rep(numgenes,2),
                    alternative="two.sided")
                if (is.numeric(toplot)) 
                    toplot <- colnames(righttable)[toplot]
            }
            
            asumar <- colSums(righttable)
            asumar <- which(asumar < 0.25)
            if (length(asumar) > 1) {
                righttable <- cbind(righttable[,-asumar],
                    rowSums(righttable[,asumar]))
                colnames(righttable)[ncol(righttable)] <- "Others"
            }

            # Detection in the genome
            bbb <- barplot(lefttable,main="Biotype detection over genome total",
                xlab=NULL,ylab="% detected features",axis.lty=1,legend=FALSE,
                cex.names=0.8,beside=TRUE,col=c(2,4),las=2,density=80,
                border=c(2,4),ylim=c(0,100))
            bbb <- colSums(bbb)/2
            lines(bbb,dat$genome,pch=20,type="o",lwd=2)

            # %detection in the sample                
            barplot(righttable,main="Relative biotype abundance in sample",
                xlab=NULL,ylab="Relative % biotypes",axis.lty=1,legend=FALSE,
                beside=TRUE,col=c(2,4),las=2,border=c(2,4))

            legend(x="topright",bty="n",horiz=FALSE,pch=c(15,15,20),
                lwd=c(NA,NA,1),legend=c(samples,"% in genome"),col=c(2,4,1))

            if ((toplot != 1) && (toplot != "global")) {
                print(paste("Percentage of",toplot,"biotype in each sample:"))
                names(mytest$estimate) <- samples
                print(round(mytest$estimate*100,4))
                print(paste("Confidence interval at 95% for the difference of",
                    "percentages:",samples[1],"-",samples[2]))
                print(round(mytest$conf.int[c(1,2)]*100,4))
                if (mytest$p.value < 0.05) {
                    print(paste("The percentage of this biotype is",
                        "significantly DIFFERENT for these two samples",
                        "(p-value =",signif(mytest$p.value,4),")."))
                } 
                else {
                    print(paste("The percentage of this biotype is NOT",
                        "significantly different for these two samples",
                        "(p-value =",signif(mytest$p.value,4),")."))
                }
            }
        }
    }
    
    # Reset with the default values
    if (!toreport) 
        par(mypar)
}

################################################################################

################################################################################

# countsbio.plot.R

countsbio.dat <- function (input,biotypes=NULL,factor=NULL,norm=FALSE)  {
    if (inherits(input,"eSet") == FALSE)
        stop("Error. You must give an eSet object\n")

    if (!is.null(assayData(input)$exprs))
        datos <- assayData(input)$exprs
    else
        datos <- assayData(input)$counts

    depth <- round(colSums(datos)/10^6,1)
    names(depth) <- colnames(datos)
    ceros <- which(rowSums(datos) == 0)
    hayceros <- length(ceros) > 0
        
    if (hayceros) {
        warning("Warning:",length(ceros),"features with 0 counts in all ",
            "samples are to be removed for this analysis.")
        datos0=datos[-ceros,,drop=FALSE]
    } 
    else { 
        datos0=datos
    }
    
    nsam <- NCOL(datos)
    
    if (nsam == 1) {
        datos <- as.matrix(datos)
        datos0 <- as.matrix(datos0)
    }

    # Per condition
    if (is.null(factor)) {    # per sample    
        print("Count distributions are to be computed for:")
        print(colnames(datos))
    } 
    else {    # per condition
        mifactor <- pData(input)[,factor]
        niveles <- levels(mifactor)
        print("Counts per million distributions are to be computed for:")
        print(niveles)
        if (norm) {
            datos <- vapply(niveles,function (k) { 
                rowMeans(as.matrix(datos[,grep(k,mifactor)]))
            },numeric(1))
            datos0 <- vapply(niveles,function (k) { 
                rowMeans(as.matrix(datos0[,grep(k,mifactor)]))
            },numeric(1))
        } 
        else {
            datos <- vapply(niveles,function (k) { 
                10^6*rowMeans(t(t(datos[,grep(k,mifactor)])/colSums(as.matrix(
                    datos[,grep(k,mifactor)]))))
            },numeric(1))
            datos0 <- vapply(niveles,function (k) { 
                10^6*rowMeans(t(t(datos0[,grep(k,mifactor)])/colSums(as.matrix(
                    datos0[,grep(k,mifactor)]))))
            },numeric(1))            
        }
        colnames(datos) <- colnames(datos0) <- niveles
        depth <- vapply(niveles,function (k) {
            paste(range(depth[grep(k,mifactor)]),collapse="-")
        },character(1))
    }

    # Biotypes
    if (!is.null(featureData(input)$Biotype)) {
        if (hayceros) {
            infobio0 <- as.character(featureData(input)$Biotype)[-ceros] 
        } 
        else { 
            infobio0 <- as.character(featureData(input)$Biotype) 
        }
        infobio <- as.character(featureData(input)$Biotype)
    } 
    else { 
        infobio0 <- NULL 
        infobio <- NULL 
    }

    if (!is.null(infobio)) {        
        if(is.null(biotypes)) {
            biotypes <- unique(infobio)
            names(biotypes) <- biotypes        
        }     
        # which genes belong to each biotype
        biog <- lapply(biotypes,function(x) { 
            which(is.element(infobio0,x)) 
        })
        names(biog) <- biotypes
        bionum <- c(NROW(datos0),vapply(biog,length,integer(1)))
        names(bionum) <- c("global",names(biotypes)) 
        bio0 <- which(bionum == 0)
        if (length(bio0) > 0) 
            bionum <- bionum[-bio0]
    } 
    else { 
        biotypes <- NULL
        bionum <- NULL 
    }
    
    # Create the summary matrix information
    if (is.null(bionum)) {
        resumen <- vector("list",length=1)
        names(resumen) <- "global"
    } 
    else {
        resumen <- vector("list",length=length(bionum))
        names(resumen) <- names(bionum)
    }

    cuentas <- c(0,1,2,5,10)

    if (is.null(factor)) {
        if (norm) {
            datosCPM <- datos
        } 
        else {
            datosCPM=10^6 * t(t(datos)/colSums(as.matrix(datos)))
        }
    } else { datosCPM=datos }

    for (i in seq_len(length(resumen))) {
        if (i == 1) {
            datosR=datosCPM             
        } else {             
            if(!is.null(infobio)) {
                datosR <- datosCPM[which(infobio == names(resumen)[i]),]
                if (!is(datosR,"matrix")) {
                    datosR <- t(as.matrix(datosR))
                }
            }
        }

        nfeatures <- nrow(datosR)
        datosR <- datosR[which(rowSums(datosR) > 0),]
        if (!is(datosR,"matrix")) {
            datosR <- t(as.matrix(datosR)) 
        }

        myglobal <- NULL
        mypersample <- NULL

        for (kk in seq_len(length(cuentas))) {
            mypersample <- rbind(mypersample,apply(datosR,2,function (x) { 
                length(which(x > cuentas[kk])) 
            }))
            myglobal <- c(myglobal,sum(apply(datosR,1,function (x) { 
                max(x) > cuentas[kk] 
            })))
        } 

        mypersample <- round(100*mypersample/nfeatures,1)
        mypersample <- rbind(mypersample,depth)
        rownames(mypersample) <- seq_len(nrow(mypersample))
        myglobal <- c(round(100*myglobal/nfeatures,1),nfeatures)

        resumen[[i]] <- data.frame(c(paste("CPM >",cuentas),"depth"),
            mypersample,"total"=myglobal)
        colnames(resumen[[i]])[1] <- names(resumen)[i]
        colnames(resumen[[i]])[2:(ncol(resumen[[i]])-1)] <- colnames(datosR)
    }

    ## results
    cosas <- list("result"=datos0,"bionum"=bionum,"biotypes"=infobio0,
        "summary"=resumen)
    
    cosas
}

countsbio.plot <- function(dat,samples=c(1,2),toplot="global",
    plottype=c("barplot","boxplot"),toreport=FALSE,...) {

    # dat: Data coming from countsbio.dat function
    # samples: Samples to be plotted. If NULL,all samples are plotted.
    # toplot: Name of biotype (including "global") to be plotted.

    mypar <- par(no.readonly=TRUE)
    plottype <- match.arg(plottype)

    ## Preparing data
    if (is.null(samples)) { 
        if (NCOL(dat$result) == 1) {
            samples=1
        } else {
            samples <- seq_len(NCOL(dat$result)) 
        }        
    }
    if (is.numeric(toplot)) 
        toplot <- names(dat$summary)[toplot]
    if (is.numeric(samples) && !is.null(colnames(dat$result))) 
        samples <- colnames(dat$result)[samples]

    if (plottype == "barplot") {
        if ((exists("ylab") && !is.character(ylab)) || !exists("ylab")) 
            ylab <- ""

        datos <- dat$summary[[toplot]]
        mytotal <- as.numeric(datos[,"total"])
        datos <- as.matrix(datos[,samples])
        rownames(datos) <- as.character(dat$summary[[toplot]][,1])

        par(mar=c(6,4,4,2))
        barplot(as.numeric(datos[1,]),col=miscolores[1],las=2,main="",ylab="",
            density=70,ylim=c(0,100),cex.axis=0.8,names.arg="",...)
        for (i in 2:(length(mytotal)-2)) {
            barplot(as.numeric(datos[i,]),col=miscolores[i],las=2,main="",
                ylab="",add=TRUE,density=70,ylim=c(0,100),cex.axis=0.8,
                names.arg="",...)
        }
        
        bp <- barplot(as.numeric(datos[(length(mytotal)-1),]),
            col=miscolores[(length(mytotal)-1)],las=2,
            main=paste(toupper(toplot)," (",mytotal[length(mytotal)],")",
            sep=""),ylab="Sensitivity (%)",add=TRUE,names.arg=colnames(datos),
            cex.axis=0.8,density=70,ylim=c(0,100),cex.names=0.8,...)
        
        for (j in seq_len((length(mytotal)-1))) 
            abline(h=mytotal[j],col=miscolores[j],lwd=2)
        if (length(samples) <= 10) {
            mtext(side=3,text=datos["depth",],adj=0.5,at=bp,cex=0.8)
        } 
        else {
            mtext(side=3,text=datos["depth",],at=bp,cex=0.7,las=2)
        }
        legend("top",rownames(datos)[-length(mytotal)],fill=miscolores,
            density=70,bty="n",ncol=3)
        par(mar=c(5,4,4,4) + 0.1)
    }

    if (plottype == "boxplot") {
        conteos <- as.matrix(dat$result[,samples])
        if (is.numeric(samples)) 
            colnames(conteos) <- colnames(dat$result)[samples]
        else 
            colnames(conteos) <- samples
        num <- dat$bionum[toplot]
        if (is.null(num)) {
            if (toplot == "global") {
                num=nrow(conteos)
            } else {
                num=0
            }
        } 
        infobio <- dat$biotypes

        if (num == 0 && toplot != "global") 
            stop("Error: No data available. Please,change toplot parameter.")

        if ((exists("ylab") && !is.character(ylab)) || !exists("ylab"))
            ylab <- "Expression values"

        ## Plots            
        if (length(samples) == 1) {
            escala <- logscaling(conteos,base=2)
            if (is.null(infobio)) {
                boxplot(escala$data,col=miscolores[1],ylab=ylab,main="",
                    yaxt="n",...)
            } 
            else {
                par(mar=c(10,4,4,2))    
                boxplot(escala$data ~ infobio,col=miscolores,ylab=ylab,
                    main=colnames(conteos),las=2,cex.axis=0.8,cex.lab=0.9,
                    yaxt="n",...)
                cuantos <- dat$bionum[-1]
                cuantos <- cuantos[sort(names(cuantos))]
                mtext(cuantos,3,at=seq_len(length(cuantos)),cex=0.6,las=2)
            }
        } 
        else {
            if (toplot != "global") 
                conteos=conteos[which(infobio == toplot),]
            
            escala <- logscaling(conteos,base=2)
            main <- paste(toupper(toplot)," (",num,")",sep="")

            par(mar=c(6,4,2,2))    
            boxplot(escala$data,col=miscolores,ylab=ylab,main=main,las=2,
                cex.lab=0.9,cex.axis=0.8,yaxt="n",...)
        }        
        axis(side=2,at=escala$at,labels=escala$labels)
    }

    if (!toreport)
        par(mypar)
}

################################################################################

################################################################################

# dat.R

dat <- function(input,type = c("biodetection","cd","countsbio","GCbias",
    "lengthbias","saturation","PCA"),k=0,ndepth=6,factor=NULL,norm=FALSE,
    refColumn=1,logtransf = FALSE) {

    type <- match.arg(type)
    if (type == "biodetection") {
        output <- new("Biodetection",
            dat=biodetection.dat(input,factor=factor,k=k))        
    }
    
    #if (type == "cd") {
    #    output <- new("CD",dat=cd.dat(input,norm=norm,refColumn=refColumn))
    #}

    if (type == "countsbio") {
        output <- new("CountsBio",
            dat=countsbio.dat(input,factor=factor,norm=norm))
    }

    #if (type == "GCbias") {
    #    output <- new("GCbias",
    #        dat=GC.dat(input,factor=factor,norm=norm))
    #}

    #if (type == "lengthbias") {
    #    output <- new("lengthbias",
    #        dat=length.dat(input,factor=factor,norm=norm))
    #}
    
    if (type == "saturation") {
        output <- new("Saturation",dat=saturation.dat(input,k=k,ndepth=ndepth))
    }
    
    #if (type == "PCA") {
    #    output=new("PCA",dat=PCA.dat(input,norm=norm,logtransf=logtransf))
    #}
    
    output
}

################################################################################

################################################################################

# DE.plot.R

DE.plot <- function(output,q=NULL,graphic=c("MD","expr","chrom","distr"),
    pch=20,cex=0.5,col=1,pch.sel=1,cex.sel=0.6,col.sel=2,log.scale=TRUE,
    chromosomes=NULL,join=FALSE,...) {
    mypar <- par(no.readonly=TRUE)

    if (!is(output,"Output"))
        stop("Error. Output argument must contain an object generated by ",
            "noiseq or noiseqbio functions.\n")

    graphic <- match.arg(graphic)    
    noiseqbio <- "theta" %in% colnames(output@results[[1]])[c(1,2,3,4)]
    
    ## MD plot
    if (graphic == "MD") {
        if (noiseqbio) {   
            M <- output@results[[1]][,"log2FC"]
            D <- abs(output@results[[1]][,1] - output@results[[1]][,2])+1
            names(M) <- names(D) <- rownames(output@results[[1]])
            
            plot(M,D,pch=pch,xlab="M",ylab="D",cex=cex,col=col,log="y",...)
            
            if(!is.null(q)) {
                mySelection <- rownames(degenes(output,q))
                points(M[mySelection],D[mySelection],col=col.sel,pch=pch.sel,
                    cex=cex.sel)
            }
        }
        else {
            plot(output@results[[1]][,"M"],(1+output@results[[1]][,"D"]),
                pch=pch,xlab="M",ylab="D",cex=cex,col=col,log="y",...)
            if (!is.null(q)) {
                mySelection <- rownames(degenes(output,q))
                points(output@results[[1]][mySelection,"M"],
                    output@results[[1]][mySelection,"D"]+1,col=col.sel,
                    pch=pch.sel,cex=cex.sel)
            }
        }
    }
    
    ## Expression plot: Condition1 vs Condition2
    else if (graphic == "expr") {
        data <- cbind(output@results[[1]][1],output@results[[1]][2])
        rownames(data) <- rownames(output@results[[1]])
        colnames(data) <- colnames(output@results[[1]][c(1,2)])

        if (log.scale) {
            k <- min(data,na.rm=TRUE)
            escala <- logscaling(data,base=2,k=k)
            
            plot(escala$data,pch=pch,cex=cex,col=col,yaxt="n",xaxt="n",...)
            axis(side=1,at=escala$at,labels=escala$labels)
            axis(side=2,at=escala$at,labels=escala$labels)
            
            if(!is.null(q)) {
                mySelection <- rownames(degenes(output,q))    
                points(escala$data[mySelection,],col=col.sel,pch=pch.sel,
                    cex=cex.sel)
            }
        } 
        else {
            plot(data,pch=pch,cex=cex,col=col,...)            
            if(!is.null(q)) {
                mySelection <- rownames(degenes(output,q))
                points(data[mySelection,],col=col.sel,pch=pch.sel,cex=cex.sel)
            }
        }
    }
    
    ## MANHATTAN PLOT
    else if (graphic == "chrom") {
        mydata <- data.frame(as.character(output@results[[1]][,"Chrom"]),
            output@results[[1]][,c("GeneStart","GeneEnd")],
            output@results[[1]][,c(1,2)])
        mydata <- na.omit(mydata)
        colnames(mydata) <- c("chr","start","end",colnames(mydata)[-c(1,2,3)])
        
        if (is.null(chromosomes)) {    # todos los cromosomas
            chromosomes <- unique(mydata$chr)
            chromosomes <- sort(chromosomes)
        }
        
        print("REMEMBER. You are plotting these chromosomes and in this order:")
        print(chromosomes)
        
        # logarithmic scale
        if (log.scale) {
            if (min(mydata[,-c(1,2,3)],na.rm=TRUE) < 1) {
                kk <- -min(mydata[,-c(1,2,3)],na.rm=TRUE)+1
            } 
            else { 
                kk <- 0 
            }
            mydata[,-c(1,2,3)] <- log2(mydata[,-c(1,2,3)]+kk)
        }
        
        # Selecting chromosomes and ordering positions
        ordenat <- NULL
        for (cromo in chromosomes) {
            myselec <- mydata[mydata[,"chr"] == cromo,]
            myselec <- myselec[order(myselec[,"start"]),]
            ordenat <- rbind(ordenat,myselec)
        }
        
        sel.ord <- NULL
        
        if (!is.null(q)) {
            # up-regulated in first condition
            cond1 <- rownames(degenes(output,q,M="up"))
            # up-regulated in second condition
            cond2 <- rownames(degenes(output,q,M="down"))
            cond1 <- (rownames(ordenat) %in% cond1)*1
            cond2 <- (rownames(ordenat) %in% cond2)*1
            sel.ord <- cbind(cond1,cond2)
            rownames(sel.ord) <- rownames(ordenat)
        }
        
        chr.long <- aggregate(ordenat[,"end"],
            by=list(as.character(ordenat[,"chr"])),max)
        chr.long <- chr.long[match(chromosomes,chr.long[,1]),]
        
        if (join) { # si todos los cromosomas van en el mismo plot
            total.long <- sum(as.numeric(chr.long$x))
            chr.start <- cumsum(c(1,chr.long$x[-length(chr.long$x)]))
            names(chr.start) <- chr.long[,1]
            plot(c(1,total.long),c(-max(ordenat[,5]),max(ordenat[,4])),
                     type="n",xlab="",ylab="Expression data",xaxt="n")
            axis(side=1,at=chr.start,labels=chr.long[,1],font=2)
            abline(h=0,lty=2,lwd=0.5)        
            
            for (ch in chromosomes) {
                dat.chr <- ordenat[which(ordenat[,"chr"] == ch),]
                dat.chr[,c("start","end")] <- dat.chr[,c("start","end")] +
                    chr.start[ch] - 1
                rect(xleft=dat.chr[,"start"],ybottom=0,xright=dat.chr[,"end"],
                    ytop=dat.chr[,4],col="grey",border=NA)
                rect(xleft=dat.chr[,"start"],ybottom=-dat.chr[,5],
                    xright=dat.chr[,"end"],ytop=0,col="grey",border=NA)

                if (!is.null(q)) {                    
                    aux <- which(rownames(sel.ord) %in% rownames(dat.chr))
                    sel.chr1 <- dat.chr[,4]*sel.ord[aux,1]
                    sel.chr2 <- -dat.chr[,5]*sel.ord[aux,2]
                    rect(xleft=dat.chr[,"start"],ybottom=0,
                        xright=dat.chr[,"end"],ytop=sel.chr1,col=2,border=NA)
                    rect(xleft=dat.chr[,"start"],ybottom=sel.chr2,
                        xright=dat.chr[,"end"],ytop=0,col=3,border=NA)
                }

                segments(x0=dat.chr[,"start"],y0=0,x1=dat.chr[,"end"],y1=0,
                    col=4,lwd=0.5) # annotated genes
            }
            text(c(1,1),0.9*c(max(ordenat[,4]),-max(ordenat[,5])),
                colnames(mydata)[4:5],font=3,adj=0,col=2:3)
            # a plot for each chromosome
        } 
        else {        
            num.chr <- length(chromosomes)
            k <- 20
            long.prop <- round((chr.long$x / min(chr.long$x))*k,0)
            while (max(long.prop)*num.chr > 500 | max(long.prop) > 50) {
                k <- k-1
                long.prop <- round((chr.long$x / min(chr.long$x))*k,0)
            }
            forlayout <- matrix(0,num.chr,max(long.prop))
            for (i in seq_len(num.chr)) {
                forlayout[i,seq_len(long.prop[i])] <- i
            }
            layout(forlayout)
            
            if (num.chr > 1)
                par(mar=c(2,4.1,0.1,0.1))
            miylim <- c(-max(ordenat[,5],na.rm=TRUE),
                max(ordenat[,4],na.rm=TRUE))
            
            for (i in seq_len(num.chr)) {
                apintar <- ordenat[which(ordenat[,"chr"] == chromosomes[i]),2:5]
                plot(c(1,chr.long[i,2]),miylim,type="n",xlab="",
                    ylab=chromosomes[i],xaxt="n",font.lab=2,cex.lab=1.3)
                rect(xleft=apintar[,"start"],ybottom=0,xright=apintar[,"end"],
                    ytop=apintar[,3],col="grey",border=NA)
                rect(xleft=apintar[,"start"],ybottom=-apintar[,4],
                    xright=apintar[,"end"],ytop=0,col="grey",border=NA)
                
                if (!is.null(q)) {
                    aux <- which(rownames(sel.ord) %in% rownames(apintar))
                    asel <- sel.ord[aux,]*apintar[,3:4]
                    rect(xleft=apintar[,"start"],ybottom=0,
                        xright=apintar[,"end"],ytop=asel[,1],col=2,border=NA)
                    rect(xleft=apintar[,"start"],ybottom=-asel[,2],
                        xright=apintar[,"end"],ytop=0,col=3,border=NA)
                }            
                abline(h=0,lty=2,lwd=0.5)
                segments(x0=apintar[,"start"],y0=0,x1=apintar[,"end"],y1=0,
                    col=4,lwd=0.5) # annotated genes
                etiq <- mypretty(c(1,chr.long[i,2]),n=10)
                axis(side=1,at=etiq,labels=etiq)
                text(c(1,1),0.9*miylim,colnames(mydata)[5:4],font=3,adj=0,
                    col=3:2)
            }
        }
    }
    
    ## DEG distribution across biotypes/chromosomes
    else if (graphic == "distr") {        
        if(!is.null(q)) { # Computing DEG
            mySelection <- rownames(degenes(output,q))
            detect <- rownames(output@results[[1]]) %in% mySelection
        }
        else {
            stop("You must specify a valid value for q\n")
        }

        if ("Chrom" %in% colnames(output@results[[1]])) {
            numplot <- 1
            infobio <- output@results[[1]][,"Chrom"]
            genome <- 100*table(infobio)/sum(table(infobio))
            ordre <- order(genome,decreasing=TRUE)            
            perdet1 <- genome*table(infobio,detect)[names(
                genome),2]/table(infobio)[names(genome)]
            perdet2 <- 100*table(infobio,detect)[names(
                genome),2]/sum(table(infobio,detect)[,2])
            ceros <- rep(0,length(genome))
            biotable1 <- as.matrix(rbind(genome[ordre],perdet1[ordre],
                perdet2[ordre],ceros))
            rownames(biotable1) <- c("genome","degVSgenome","deg","ceros")
            if (!is.null(chromosomes))
                biotable1=biotable1[,chromosomes]
            ymaxL1 <- ceiling(max(biotable1,na.rm=TRUE))
        } 
        else { 
            numplot <- 0 
        }
        
        if ("Biotype" %in% colnames(output@results[[1]])) {
            numplot <- c(numplot,1)
            infobio <- output@results[[1]][,"Biotype"]
            genome <- 100*table(infobio)/sum(table(infobio))
            ordre <- order(genome,decreasing=TRUE)            
            perdet1 <- genome*table(infobio,detect)[names(
                genome),2]/table(infobio)[names(genome)]
            perdet2 <- 100*table(infobio,detect)[names(
                genome),2]/sum(table(infobio,detect)[,2])
            ceros <- rep(0,length(genome))
            biotable2 <- as.matrix(rbind(genome[ordre],perdet1[ordre],
                perdet2[ordre],ceros))
            rownames(biotable2) <- c("genome","degVSgenome","deg","ceros")
            higher2 <- which(biotable2[1,] > 2)
            lower2 <- which(biotable2[1,] <= 2)
            if (length(higher2) > 0) {       
                ymaxL2 <- ceiling(max(biotable2[,higher2],na.rm=TRUE)) 
                if (length(lower2) > 0) {
                    ymaxR2 <- ceiling(max(biotable2[,lower2],na.rm=TRUE))
                    biotable2[,lower2] <- biotable2[,lower2]*ymaxL2/ymaxR2
                }
                else { 
                    ymaxR2 <- ymaxL2 
                }
                
            } 
            else { 
                ymaxR2 <- ceiling(max(biotable2[,lower2],na.rm=TRUE))
                ymaxL2 <- ymaxR2                             
            }
        } 
        else { 
            numplot <- c(numplot,0) 
        }
        
        # Plot
        if (sum(numplot) == 0) 
            stop("Biotype or chromosome information is needed for this plot\n")
        
        if (sum(numplot) == 1) {    # 1 Plot 
            par(mar=c(10,4,2,2))
            if (numplot[1] == 1) {
                barplot(biotable1[c(1,3),],
                    main="DEG distribution across chromosomes",
                    xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,
                    beside=TRUE,col=c("grey",2),las=2,ylim=c(0,ymaxL1),
                    border=c("grey",2))
                barplot(biotable1[c(2,4),],
                    main="DEG distribution across chromosomes",
                    xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,
                    beside=TRUE,col=c(2,1),las=2,density=30,ylim=c(0,ymaxL1),
                    xborder=2,add=TRUE)
                legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",2,2),
                    density=c(NA,30,NA),border=c("grey",2,2),
                    legend=c("%Chrom in genome","%DEG in Chrom",
                    "%Chrom in DEG"))
            }
            
            if (numplot[2] == 1) {
                barplot(biotable2[c(1,3),],
                    main="DEG distribution across biotypes",
                    xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,
                    beside=TRUE,col=c("grey",4),las=2,ylim=c(0,ymaxL2),
                    border=c("grey",4))
                barplot(biotable2[c(2,4),],
                    main="DEG distribution across biotypes",
                    xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,
                    beside=TRUE,col=c(4,1),las=2,density=30,ylim=c(0,ymaxL2),
                    border=4,add=TRUE)
                axis(side=4,at=pretty(c(0,ymaxL2),n=5),
                    labels=round(pretty(c(0,ymaxL2),n=5)*ymaxR2/ymaxL2,1))
                if (ymaxR2 != ymaxL2) {
                    abline(v=3*length(higher2) + 0.5,col=3,lwd=2,lty=2)
                }
                legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",4,4),
                    density=c(NA,30,NA),border=c("grey",4,4),
                    legend=c("%Biotype in genome","%DEG in Biotype",
                    "%Biotype in DEG"))
            }
        }
        
        if (sum(numplot) == 2) {    # 2 Plots 
            par(mar=c(10,4,2,2),mfrow=c(1,2))
            
            # Chromosomes
            barplot(biotable1[c(1,3),],
                main="DEG distribution across chromosomes",
                xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c("grey",2),las=2,ylim=c(0,ymaxL1),border=c("grey",2))
            barplot(biotable1[c(2,4),],
                main="DEG distribution across chromosomes",xlab=NULL,
                ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,col=c(2,1),
                las=2,density=30,ylim=c(0,ymaxL1),border=2,add=TRUE)
            legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",2,2),
                density=c(NA,30,NA),border=c("grey",2,2),
                legend=c("%Chrom in genome","%DEG in Chrom","%Chrom in DEG"))
            
            # Biotypes
            barplot(biotable2[c(1,3),],main="DEG distribution across biotypes",
                xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c("grey",4),las=2,ylim=c(0,ymaxL2),border=c("grey",4))
            barplot(biotable2[c(2,4),],main="DEG distribution across biotypes",
                xlab=NULL,ylab="%features",axis.lty=1,legend=FALSE,beside=TRUE,
                col=c(4,1),las=2,density=30,ylim=c(0,ymaxL2),border=4,add=TRUE)
            axis(side=4,at=pretty(c(0,ymaxL2),n=5),
                labels=round(pretty(c(0,ymaxL2),n=5)*ymaxR2/ymaxL2,1))
            if (ymaxR2 != ymaxL2) {
                abline(v=3*length(higher2) + 0.5,col=3,lwd=2,lty=2)
            }        
            legend(x="topright",bty="n",horiz=FALSE,fill=c("grey",4,4),
                density=c(NA,30,NA),border=c("grey",4,4),
                legend=c("%Biotype in genome","%DEG in Biotype",
                    "%Biotype in DEG"))
        }
    }
    par(mypar)
}

## Auxiliar function
mypretty <- function(x,n=5) {    
    mywidth <- diff(x)
    mybin <- ceiling(mywidth/(n-1))
    mylabels0 <- x[1] + mybin*(0:(n-1))
    ndig <- nchar(as.character(mylabels0))
    mylabels <- mylabels0
    k <- 0
    for (i in 2:(n-1)) {
        mylabels[i] <- (round(mylabels0[i]/10^(ndig[i]-1),k))*10^(ndig[i]-1)
        while (mylabels[i] == mylabels[i-1]) {
            k <- k+1
            mylabels[i] <- (round(mylabels0[i]/10^(ndig[i]-1),k))*10^(ndig[i]-1)
        }
    }
    mylabels
}

################################################################################

################################################################################

# degenes.R

degenes <- function (object,q=0.95,M=NULL) {
    if (!is(object,"Output"))
        stop("You must give the object returned by the noiseq function\n")

    x <- object@results[[1]]
    noiseqbio="theta" %in% colnames(x)[seq_len(4)]

    if (noiseqbio) {
        y <- na.omit(x[c("theta","prob")])
        colnames(y)[1]="M"
    } else {
        y <- na.omit(x[c("M","D","prob")])
    }

    if (is.null(M)) {
        losdeg <- y[y[,"prob"] > q,]
        print(paste(dim(losdeg)[1],"differentially expressed features"))

    }
    else if (M == "up") {
        estos <- y[y[,"M"] > 0,]
        losdeg <- estos[estos[,"prob"] > q,]
        print(paste(dim(losdeg)[1],
            "differentially expressed features (up in first condition)"))
    } 
    else if (M == "down") {
        estos <- y[y[,"M"] < 0,]
        losdeg <- estos[estos[,"prob"] > q,]
        print(paste(dim(losdeg)[1],
            "differentially expressed features (down in first condition)"))
    }
    else {
        stop("ERROR! Value for parameter M is not valid. Please,choose ",
            "between NULL,'up' or 'down'")
    }
    
    # Restore the object with the same "results" structure
    losdeg <- x[rownames(losdeg),]
    losdeg[order(losdeg[,"prob"],decreasing=TRUE),]
}

################################################################################

################################################################################

# fewrepliactes.R

share.info = function (mydata, n1, n2, r, nclust)   {
    # clustering data by k-means algorithm
    # 1.a) Normalized data
    gc()
    cl <- suppressWarnings(kmeans(mydata,nclust,nstart=25,iter.max=nclust + 30))

    cat("...k-means clustering done\n")
    cat(paste("Size of",nclust,"clusters:\n"))
    print(cl$size)
                        
    # Creating pseudo-data
    cluster.data <- lapply(seq_len(nclust),function (k) { 
        mydata[cl$cluster == k,]
    })
        
    # Resampling
    npermu <- cl$size*r
    npermu <- vapply(npermu,function (x) min(x,1000),numeric(1))

    cat("Resampling cluster...")

    myres=vector("list",length=nclust)

    for (i in seq_len(nclust)) {
        print(i)

        #if (cl$size[i] > 1) {     # OPTION 2.A
        # OPTION 2.C: small clusters
        if (cl$size[i] > 1 && cl$size[i] < 1000) {
            #myres[[i]]=t(sapply(seq_len(npermu[i]),function (j) {
            #    permu=sample(cluster.data[[i]])
            #    nn1=n1*cl$size[i]
            #    nn2=n2*cl$size[i]
            #    mean1=mean(permu[seq_len(nn1)])
            #    mean2=mean(permu[(nn1+1):(nn1+nn2)])
            #    sd1=sd(as.numeric(permu[seq_len(nn1)]))
            #    sd2=sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
            #    data.frame("M"=log2(mean1/mean2),"D"=mean1-mean2,
            #        "M.sd"=sqrt(sd1^2/(mean1^2*log(2)^2*nn1) + 
            #        sd2^2/(mean2^2*log(2)^2*nn2)),
            #        "D.sd"=sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))
            #}))
            myres[[i]]=t(vapply(seq_len(npermu[i]),function (j) {
                permu=sample(cluster.data[[i]])
                nn1=n1*cl$size[i]
                nn2=n2*cl$size[i]
                mean1=mean(permu[seq_len(nn1)])
                mean2=mean(permu[(nn1+1):(nn1+nn2)])
                sd1=sd(as.numeric(permu[seq_len(nn1)]))
                sd2=sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
                data.frame("M"=log2(mean1/mean2),"D"=mean1-mean2,
                    "M.sd"=sqrt(sd1^2/(mean1^2*log(2)^2*nn1) + 
                    sd2^2/(mean2^2*log(2)^2*nn2)),
                    "D.sd"=sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))
            },vector("list",4)))
        }
        
        if (cl$size[i] >= 1000) {    # OPTION 2.C & 2.D: big clusters
            # Option 2.D: clustering big clusters
            cl2 <- kmeans(cluster.data[[i]],nclust,nstart=25,
                iter.max=nclust + 20)
            cat(paste("Size of",nclust,"subclusters of cluster:",i)); cat("\n")
            print(cl2$size)
            subcluster.data <- lapply(seq_len(nclust),function (k) { 
                cluster.data[[i]][cl2$cluster == k,] 
            })
            npermu2 <- cl2$size*r
            npermu2 <- vapply(npermu2,function (x) min(x,1000),numeric(1))
            ## modified to reduce the number of permutations to be done
            myres2 <- vector("list",length=nclust)
            
            for (h in seq_len(nclust)) {
                if (cl2$size[h] > 1) {
                    #myres2[[h]] <- t(sapply(seq_len(npermu2[h]),function (j) {
                    #    permu=sample(subcluster.data[[h]])
                    #    nn1 <- n1*cl2$size[h]
                    #    nn2 <- n2*cl2$size[h]             
                    #    mean1 <- mean(permu[seq_len(nn1)])
                    #    mean2 <- mean(permu[(nn1+1):(nn1+nn2)])
                    #    sd1 <- sd(as.numeric(permu[seq_len(nn1)]))
                    #    sd2 <- sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
                    #    data.frame("M"=log2(mean1/mean2),"D"=mean1-mean2,
                    #        "M.sd"=sqrt(sd1^2/(mean1^2*log(2)^2*nn1) +
                    #        sd2^2/(mean2^2*log(2)^2*nn2)),
                    #        "D.sd"=sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))
                    #}))
                    myres2[[h]] <- t(vapply(seq_len(npermu2[h]),function (j) {
                        permu=sample(subcluster.data[[h]])
                        nn1 <- n1*cl2$size[h]
                        nn2 <- n2*cl2$size[h]             
                        mean1 <- mean(permu[seq_len(nn1)])
                        mean2 <- mean(permu[(nn1+1):(nn1+nn2)])
                        sd1 <- sd(as.numeric(permu[seq_len(nn1)]))
                        sd2 <- sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
                        data.frame("M"=log2(mean1/mean2),"D"=mean1-mean2,
                            "M.sd"=sqrt(sd1^2/(mean1^2*log(2)^2*nn1) +
                            sd2^2/(mean2^2*log(2)^2*nn2)),
                            "D.sd"=sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))
                    },vector("list",4)))
                }
            }
            myres[[i]] <- do.call("rbind",myres2)
        }
    }                                        

    # Computing Zr for noise
    cat("Computing Z for noise...\n")

    # 4.A) a0: Global for all R*G permutations
    myres <- do.call("rbind",myres)
    
    a0.M <- quantile(as.numeric(myres[,"M.sd"]),probs=0.9,na.rm=TRUE)
    a0.D <- quantile(as.numeric(myres[,"D.sd"]),probs=0.9,na.rm=TRUE)
    
    M <- as.numeric(myres[,"M"])/(a0.M + as.numeric(myres[,"M.sd"]))
    D <- as.numeric(myres[,"D"])/(a0.D + as.numeric(myres[,"D.sd"]))
    (M + D)/2
}

################################################################################

# filter.low.counts.R

CV = function(data) { 
    100*sd(data,na.rm=TRUE)/mean(data,na.rm=TRUE)
}

filtered.data <- function(dataset,factor,norm=TRUE,depth=NULL,method=1,
    cv.cutoff=100,cpm=1,p.adj="fdr") {    
    dataset0=dataset[rowSums(dataset) > 0,]
    dataset=dataset0

    if ((method == 3) && (norm)) {
        if (is.null(depth)) {
            stop("ERROR: Sequencing depth for each column in dataset must ",
                "be provided.\n")
        }
        dataset <- t(t(dataset0)/(colSums(dataset0)/depth))
    }
    
    if ((method < 3) && (!norm)) {
        dataset <- 10^6*t(t(dataset0)/colSums(dataset0))
    }

    grupos <- unique(factor)
    cumple <- NULL

    cat("Filtering out low count features...\n")

    for (gg in grupos) {        
        datos <- as.matrix(dataset[,grep(gg,factor)])

        if (method == 1) {
            if (ncol(datos) == 1) {
                cumplecond=(datos > cpm)
            } 
            else {
                cumplecond <- 
                    (apply(datos,1,CV) < cv.cutoff)*(rowMeans(datos) > cpm)
                cumplecond[which(is.na(cumplecond) == TRUE)] <- 0 
            }
            cumple <- cbind(cumple,cumplecond)
        }

        if (method == 2) {
            if (ncol(datos) == 1) 
                stop("ERROR: At least 2 replicates per condition are ",
                    "required to apply this method.")
            mytest <- apply(datos,1,function (x) { 
                suppressWarnings(wilcox.test(x,alternative="greater",
                    conf.int=FALSE,mu=0))$"p.value" 
            })
            mytest <- p.adjust(mytest,method=p.adj)
            cumple <- cbind(cumple,1*(mytest < 0.05))    
        }
        
        if (method == 3) {                     
            p0=cpm/10^6                                
            mytest <- apply(datos,1,function (x) 
                suppressWarnings(prop.test(sum(x),n=sum(datos),p=p0,
                    alternative="greater"))$"p.value")
            mytest <- p.adjust(mytest,method=p.adj)
            cumple <- cbind(cumple,1*(mytest < 0.05))                     
        }
    }

    cumple <- which(rowSums(as.matrix(cumple)) >= 1) 

    cat(paste(length(cumple),"features are to be kept for differential ",
        "expression analysis with filtering method",method))
    cat("\n")

    dataset0[cumple,]
}

################################################################################

################################################################################

# makeASCAdesign.R

make.ASCA.design <- function(x) {
    x <- as.factor(x)
    levels <- unique(x)
    n <- length(x)
    p <- length(levels)
    Design <- matrix(0,nrow=n,ncol=p)
    for (i in seq_len(n)){
        for (j in seq_len(p)){
            if (x[i]==levels[j]) {
                Design[i,j]=1
            }
        }
    }
    colnames(Design)<-levels
    output<-Design
    output
}

################################################################################

################################################################################

# MD.R

MD <- function (dat=dat,selec=c(seq_len(nrow(dat)))) {
    pares <- as.matrix(combn(ncol(dat), 2))
    if (NCOL(pares) > 30) {  
        sub30 <- sample(seq_len(NCOL(pares)),size=30,replace=FALSE)
        pares <- pares[,sub30]
    }
    mm <- NULL
    dd <- NULL
    for (i in seq_len(ncol(pares))) {
        a <- dat[selec,pares[1,i]]
        b <- dat[selec,pares[2,i]]
        mm <- cbind(mm,log(a/b,2))
        dd <- cbind(dd,abs(a-b))
    }
    list("M"=mm,"D"=dd)
}

MDbio <- function (dat=dat,selec= c(seq_len(nrow(dat))),param=NULL,a0per=0.9) {
    pares <- as.matrix(combn(ncol(dat), 2))
    mm <- NULL
    dd <- NULL
    for (i in seq_len(ncol(pares))) {
        a <- dat[selec,pares[1,i]]
        b <- dat[selec,pares[2,i]]
        mm <- cbind(mm, log(a/b, 2))
        dd <- cbind(dd, (a-b))  
    }

    ## Correcting (M,D) 
    sd.M <- sqrt(param$sd[,1]^2/(dat[,1]^2*log(2)^2*param$n[1]) + 
        param$sd[,2]^2/(dat[,2]^2*log(2)^2*param$n[2]))

    sd.D = sqrt(param$sd[,1]^2/param$n[1] + param$sd[,2]^2/param$n[2])

    if(is.null(a0per)) {
        a0.M = a0.D = 0
    } 
    else {
        if (a0per == "B") {
            B <- 100
            a0.M <- B*max(sd.M,na.rm=TRUE)
            a0.D <- B*max(sd.D,na.rm=TRUE)
      
        } 
        else {
            a0per <- as.numeric(a0per)
            a0.M <- quantile(sd.M, probs = a0per, na.rm = TRUE)
            a0.D <- quantile(sd.D, probs = a0per, na.rm = TRUE)
        }
    }

    mm <- mm / (a0.M + sd.M)
    dd <- dd / (a0.D + sd.D)

    # Results
    list("M"=mm,"D"=dd)
}

################################################################################

################################################################################

# noiseq.R

noiseq <- function (input,k=0.5,norm=c("rpkm","uqua","tmm","n"),
    replicates=c("technical","biological","no"),factor=NULL,conditions=NULL,
    pnr=0.2,nss=5,v=0.02,lc=0) {
    if (inherits(input,"eSet") == FALSE)
        stop("Error. You must give an object generated by the readData ",
            "function\n")

    if (is.null(factor))
        stop("Error. You must specify the factor to know which conditions ",
            "you wish to compare. Please,give the argument 'factor'.\n")

    replicates <- match.arg(replicates)
    if (replicates == "biological") 
        print("WARNING: Your experiment has biological replicates. You ",
            "should consider using NOISeqBIO instead of NOISeq.")
    norm <- match.arg(norm)
    
    # (M,D) for signal and noise
    print("Computing (M,D) values...")
    miMD <- allMD(input,factor,k=k,replicates=replicates,norm=norm,
        conditions=conditions,pnr=pnr,nss=nss,v=v,lc=lc)

    ## Probability of differential expression

    print("Computing probability of differential expression...")
    prob.concounts <- probdeg(miMD$Ms,miMD$Ds,miMD$Mn,miMD$Dn)$prob

    if (!is.null(assayData(input)$exprs))
        todos <- rownames(as.matrix(assayData(input)$exprs))
    else
        todos <- rownames(as.matrix(assayData(input)$counts))
    
    prob <- prob.concounts[todos]
    names(prob) <- todos

    ## Results
    resultat <- data.frame("level_1"=miMD$Level1,"level_2"=miMD$Level2,
        "M"=miMD$Ms,"D"=miMD$Ds,"prob"=prob)
    rownames(resultat) <- todos

    # We change the name of the conditions to "name_mean"
    colnames(resultat)[1] <- 
        paste(unlist(strsplit(miMD$comp," "))[1],"mean",sep="_")
    colnames(resultat)[2] <- 
        paste(unlist(strsplit(miMD$comp," "))[3],"mean",sep="_")

    resultat <- data.frame(resultat,"ranking"=ranking(resultat)$statistic)
    if (!is.null(featureData(input)@data$Length))
        resultat <- data.frame(resultat,
            "Length"=as.numeric(as.character(featureData(input)@data[todos,
            "Length"])),stringsAsFactors=FALSE)
    if (!is.null(featureData(input)@data$GC))
        resultat <- data.frame(resultat,"GC"=as.numeric(as.character(
            featureData(input)@data[todos,"GC"])),stringsAsFactors=FALSE)
    if (!is.null(featureData(input)@data$Chromosome))
        resultat <- data.frame(resultat,"Chrom"=as.character(
            featureData(input)@data$Chromosome),"GeneStart"=as.numeric(
            as.character(featureData(input)@data$GeneStart)),
            "GeneEnd"=as.numeric(as.character(featureData(input)@data$GeneEnd)),
            stringsAsFactors=FALSE)
    if (!is.null(featureData(input)@data$Biotype))
        resultat <- data.frame(resultat,"Biotype"=as.character(
            featureData(input)@data[todos,"Biotype"]),stringsAsFactors=FALSE)
    #resultat[order(resultat[,5],decreasing=TRUE),]
    
    Output(data=list(resultat),method=norm,k=miMD$k,lc=lc,factor=factor,v=v,
        nss=nss,pnr=pnr,comparison=miMD$comp,replicates=replicates)
}

################################################################################

################################################################################

# noiseqbio.R

noiseqbio <- function (input,k=0.5,norm=c("rpkm","uqua","tmm","n"),nclust=15,
    plot=FALSE,factor=NULL,conditions=NULL,lc=0,r=50,adj=1.5,a0per=0.9,filter=1,
    depth=NULL,cv.cutoff=500,cpm=1) {
    if (inherits(input,"eSet") == FALSE)
        stop("ERROR: You must give an object generated by the readData ",
            "function\n")

    if (is.null(factor))
        stop("ERROR: You must specify the factor to know which conditions ",
            "you wish to compare. Please,give the argument 'factor'.\n")

    if (min(table(input@phenoData@data[,factor])) < 2)
        stop("ERROR: To run NOISeqBIO at least two replicates per condition ",
            "are needed. Please, run NOISeq if there are not enough ",
            "replicates in your experiment.\n")

    norm <- match.arg(norm)

    # Z-scores for signal and noise
    cat("Computing Z values...\n")
    miMD <- allMDbio(input,factor,k=k,norm=norm,conditions=conditions,lc=lc,
        r=r,a0per=a0per,nclust=nclust,filter=filter,depth=depth,
        cv.cutoff=cv.cutoff,cpm=cpm)

    ##### Probability of differential expression
    cat("Computing probability of differential expression...\n")

    ## KDE estimators of f0 and f    
    desde <- min(c(miMD$Zs,miMD$Zn),na.rm=TRUE)
    hasta <- max(c(miMD$Zs,miMD$Zn),na.rm=TRUE)

    fdens <- density(miMD$Zs,adjust=adj,n=5000,from=desde,to=hasta,na.rm=TRUE)
    f <- approxfun(fdens)

    f0dens <- density(miMD$Zn,adjust=adj,n=5000,from=desde,to=hasta,na.rm=TRUE)
    f0 <- approxfun(f0dens)

    if (f0(0)/f(0) < 1) 
        print("WARNING: f0(0)/f(0) < 1 => FP with Z~0 will be detected.")

    f0.f <- f0(miMD$Zs)/f(miMD$Zs)

    ## f0,f plot
    if (plot) {
        plot(f0dens,lwd=2,col=4,main=paste("r=",r,"; a0per=",a0per,sep=""),
            xlab="theta")
        lines(fdens,lwd=2,col=1)
        legend("topleft",c("f","f0"),col=c(1,4),lwd=2)
    }

    ## ESTIMATION of p0
    p0 <- min(1/f0.f,na.rm=TRUE)
    p0=min(p0,1)

    cat(paste("p0 =",p0)); cat("\n")

    ## PROBABILITY of DIFFERENTIAL EXPRESSION
    myprob <- 1 - p0*f0.f
    if (min(myprob,na.rm=TRUE) < 0) {
        print(summary(f0.f))
    }
    names(myprob) <- names(miMD$Zs)

    cat("Probability\n")
    print(summary(myprob))

    if (!is.null(assayData(input)$exprs))
        todos <- rownames(as.matrix(assayData(input)$exprs))
    else
        todos <- rownames(as.matrix(assayData(input)$counts))

    myprob <- myprob[todos]
    names(myprob) <- todos

    ## Results
    resultat <- data.frame("level_1"=miMD$Level1,"level_2"=miMD$Level2,
        "theta"=miMD$Zs,"prob"=myprob,"log2FC"=log2(miMD$Level1/miMD$Level2))
    rownames(resultat) <- todos

    colnames(resultat)[1] <- 
        paste(unlist(strsplit(miMD$comp," "))[1],"mean",sep="_")
    colnames(resultat)[2] <- 
        paste(unlist(strsplit(miMD$comp," "))[3],"mean",sep="_")

    # resultat <- data.frame(resultat,"ranking"=ranking(resultat)$statistic)
    if (!is.null(featureData(input)@data$Length))
        resultat <- data.frame(resultat,"length"=as.numeric(as.character(
            featureData(input)@data[todos,"Length"])),stringsAsFactors=FALSE)
    if (!is.null(featureData(input)@data$GC))
        resultat <- data.frame(resultat,"GC"=as.numeric(as.character(
            featureData(input)@data[todos,"GC"])),stringsAsFactors=FALSE)
    if (!is.null(featureData(input)@data$Chromosome))
        resultat <- data.frame(resultat,"Chrom"=(as.character(
            featureData(input)@data$Chromosome)),
            "GeneStart"=as.numeric(as.character(featureData(
            input)@data$GeneStart)),"GeneEnd"=as.numeric(as.character(
            featureData(input)@data$GeneEnd)),stringsAsFactors=FALSE)

    if (!is.null(featureData(input)@data$Biotype))
        resultat <- data.frame(resultat,"Biotype"=as.character(
            featureData(input)@data[todos,"Biotype"]),stringsAsFactors=FALSE)
    
    Output(data=list(resultat),method=norm,k=miMD$k,lc=lc,factor=factor,
        comparison=miMD$comp,replicates="biological",v=0,nss=0,pnr=0)
}

################################################################################

################################################################################

# normalization.R

noitmm <- function (datos,long=1000,lc=0,k=0,refColumn=1,logratioTrim=.3,
    sumTrim=0.05,doWeighting=TRUE,Acutoff=-1e10) {
    if (!is.null(ncol(long))) {
        mynames <- long[,1]
        long <- long[,2]
        names(long) <- mynames
    }

    L <- (long/1000)^lc
    datos <- datos/L

    total <- colSums(as.matrix(datos))

    datos0 <- sinceros(datos,k)

    if (ncol(as.matrix(datos)) > 1) {
        fk <- .calcNormFactors(as.matrix(datos),refColumn=refColumn,
            method="TMM",logratioTrim=logratioTrim,sumTrim=sumTrim,
            doWeighting=doWeighting,Acutoff=Acutoff)

        fk <- fk * (total/mean(total))
        datos.norm <- t(t(datos0)/fk)
    } 
    else {
        datos.norm <- datos0/L
    }
    na.omit(datos.norm)    
    
}

noirpkm <- function (datos,long=1000,lc=1,k=0) {
    if (!is.null(ncol(long))) {
        mynames=long[,1]
        long=long[,2]
        names(long)=mynames
    }

    total <- colSums(as.matrix(datos))
    datos0 <- sinceros(datos,k)                                     
    datos.norm <- (t(t(datos0)/total)*10^6)/((long/1000)^lc)
    na.omit(datos.norm)     
}

noiuqua <- function (datos,long=1000,lc=0,k=0) {
    if (!is.null(ncol(long))) {
        mynames=long[,1]
        long=long[,2]
        names(long)=mynames
    }

    L <- (long/1000)^lc
    datos=datos/L

    datos0 <- sinceros(datos,k)

    if (ncol(as.matrix(datos)) > 1) {
        sumatot <- rowSums(datos)
        supertot <- sum(sumatot)
        counts0 <- which(sumatot == 0)

        if (length(counts0) > 0) {
            datitos <- datos[-counts0,]
        } else {
            datitos <- datos
        }

        q3 <- apply(datitos,2,quantile,probs=0.75)
        d <- q3*supertot/sum(q3)

        datos.norm <- t(t(datos0)/d)*10^6

    } else {

        datos.norm <- datos0/L

    }
    
    na.omit(datos.norm)
}

## Taken from the edgeR package with minor modifications
.calcNormFactors <- function(object,method=c("TMM","quantile"),refColumn=NULL,
    logratioTrim=.3,sumTrim=0.05,doWeighting=TRUE,Acutoff=-1e10,quantile=0.75) {
    method <- match.arg(method)
    if( is.matrix(object) ) {
        if(is.null(refColumn))
            refColumn <- 1
        data <- object
        libsize <- colSums(data)
    } 
    else {
        stop(".calcNormFactors() only operates on 'matrix' objects")
    }

    f <- switch(method,
        TMM = apply(data,2,.calcFactorWeighted,ref=data[,refColumn],
            logratioTrim=logratioTrim,sumTrim=sumTrim,doWeighting=doWeighting,
            Acutoff=Acutoff),quantile=.calcFactorQuantile(data,libsize,
            q=quantile))
    f <- f/exp(mean(log(f)))

    return(f)
}

.calcFactorQuantile <- function(data,lib.size,q=0.75) {
    y <- t(t(data)/lib.size)
    f <- apply(y,2,function(x) quantile(x,p=q))
    f/exp(mean(log(f)))
}

.calcFactorWeighted <- function(obs,ref,logratioTrim=.3,sumTrim=0.05,
    doWeighting=TRUE,Acutoff=-1e10) {
    if( all(obs==ref) )
        return(1)

    obs <- as.numeric(obs)
    ref <- as.numeric(ref)

    nO <- sum(obs)
    nR <- sum(ref)
    # log ratio of expression,accounting for library size
    logR <- log2((obs/nO)/(ref/nR)) 
    absE <- (log2(obs/nO) + log2(ref/nR))/2 # absolute expression
    v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref # estimated asymptotic variance

    # remove infinite values,cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]

    # taken from the original mean() function
    n <- sum(fin)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    
    #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    # a fix from leonardo ivan almonacid cardenas,since rank() can return
    # non-integer values when there are a lot of ties
    keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS 
        & rank(absE)<=hiS)
    
    if (doWeighting) 
        2^(sum(logR[keep]/v[keep],na.rm=TRUE) / sum(1/v[keep],na.rm=TRUE))
    else
        2^(mean(logR[keep],na.rm=TRUE))
}

################################################################################

################################################################################

# PCA.GENES.R

PCA.GENES <- function(X) {
    n <- ncol(X)
    p <- nrow(X)
    offset <- apply(X,2,mean)
    Xoff <- X-(cbind(matrix(1,p,1))%*%rbind(offset))

    eigen <- eigen(Xoff%*%t(Xoff)/(p-1))
    var <- cbind(eigen$values/sum(eigen$values),
        cumsum(eigen$values/sum(eigen$values)))

    loadings2 <- eigen$vectors
    scores2 <- t(Xoff)%*%loadings2

    normas2 <- sqrt(apply(scores2^2,2,sum))

    scores1 <- loadings2%*%diag(normas2)
    loadings1 <- scores2%*%diag(1/normas2)

    output <- list(eigen,var,scores1,loadings1)
    names(output) <- c("eigen","var.exp","scores","loadings")
    output
}

################################################################################

################################################################################

# probdeg.R

probdeg <- function (Mg,Dg,Mn,Dn,prec=2) {
    tot <- length(Mn)   # number of points in noise distribution
    gens <- names(Mg)
    Mruido <- abs(round(Mn,prec))
    Druido <- round(Dn,prec)
    Mgen <- abs(round(Mg,prec))
    Dgen <- round(Dg,prec)
    MDgen <- na.omit(cbind(Mgen, Dgen))
    MDunic <- unique(MDgen)
    Nres <- apply(MDunic,1,n.menor,S1=Mruido,S2=Druido)
    lugares <- apply(MDgen,1,busca,S=MDunic)
    Nconj <- Nres[lugares]
    names(Nconj) <- names(lugares)
    laprob <- Nconj / tot
    laprob <- laprob[gens]
    names(laprob) <- gens
    Nconj <- Nconj[gens]
    names(Nconj) <- gens
    laprob <- list("prob"=laprob,"numDE"=Nconj,"numNOISE"=tot)
    laprob
}

################################################################################

################################################################################

# readData.R

readData <- function(data=NULL,factors=NULL,length=NULL,biotype=NULL,
    chromosome=NULL,gc=NULL) {
    if (is.null(data))
        stop("Expression information must be provided to the readData function")

    if (is.null(factors))
        stop("Condition information must be provided to the readData funcion")

    if (!is.null(length) && !is.vector(length) && !is.data.frame(length)
        && !is.matrix(length))
        stop( "The length info should be a vector or a data.frame/matrix.")

    if (!is.null(gc) && !is.vector(gc) && !is.data.frame(gc) && !is.matrix(gc))
        stop( "The GC content info should be a vector or a data.frame/matrix.")
    
    if (!is.null(chromosome) && ncol(chromosome) != 3)
        stop( "The chromosome object should be a matrix or data.frame with ",
            "3 columns: chromosome,start position and end position.")

    if (!is.null(biotype) && !is.vector(biotype) && !is.data.frame(biotype) 
        && !is.matrix(biotype))
        stop( "The biotype info should be a vector or a data.frame/matrix.")

    countData <- as.matrix(data)

    rowNames <- rownames(countData)
    
    if (nrow(factors) == ncol(countData)) {
        rownames(factors) <- colnames(countData)
    }
    else {
        stop ("Number of rows in factors must be equal to number of ",
            "columns in data.\n")
    }

    pheno <- AnnotatedDataFrame(data=as.data.frame(factors))
    input <- ExpressionSet(assayData=countData,phenoData=pheno)

    if (!is.null(length))
        input <- addData(data=input,length=length)

    if (!is.null(gc))
        input <- addData(data=input,gc=gc)

    if (!is.null(biotype))
        input <- addData(data=input,biotype=biotype)

    if (!is.null(chromosome))
        input <- addData(data=input,chromosome=chromosome)

    input
}

addData <- function(data,length=NULL,biotype=NULL,chromosome=NULL,factors=NULL,
    gc=NULL) {
    if (!inherits(data,"eSet"))
        stop("Error. You must give an eSet object.")

    if (!is.null(length) && !is.vector(length) && !is.data.frame(length) 
        && !is.matrix(length))
        stop( "The length info should be a vector or a data.frame/matrix.")

    if (!is.null(gc) && !is.vector(gc) && !is.data.frame(gc) && !is.matrix(gc))
        stop( "The GC content info should be a vector or a data.frame/matrix.")

    if (!is.null(biotype) && !is.vector(biotype) && !is.data.frame(biotype) 
        && !is.matrix(biotype))
        stop( "The biotype info should be a vector or a data.frame/matrix.")
    
    if (!is.null(chromosome) && ncol(chromosome) != 3)
        stop( "The chromosome object should be a matrix or data.frame with ",
            "3 columns: chromosome,start position and end position.")

    if (!is.null(assayData(data)$exprs))
        rowNames <- rownames(assayData(data)$exprs)
    else
        rowNames <- rownames(assayData(data)$counts)

    # If exists length
    if (!is.null(length)) {
        Length <- rep(NA,length(rowNames))
        names(Length) <- rowNames
        if (is.vector(length)) {
            Length[rowNames] <- as.numeric(as.character(length[rowNames]))
        } 
        else if (is.data.frame(length) || is.matrix(length)) {
            if (ncol(length) == 2) {
                # We assume that the feature names are in the first column and 
                # the length in the second
                rownames(length) <- length[,1]
                Length[rowNames] <- as.numeric(as.character(length[rowNames,2]))
            } 
            else if (ncol(length) == 1) {
                # We assume that the length are in the first column and the 
                # feature names in the rownames
                Length[rowNames] <- as.numeric(as.character( 
                    length[rowNames,1]))
            } 
            else {
                stop( "The length matrix/data.frame contains more columns ",
                    "than expected.")
            }
        }
        featureData(data)@data <- cbind(featureData(data)@data,Length)
    }

    # If exists gc
    if (!is.null(gc)) {
        GC <- rep(NA,length(rowNames))
        names(GC) <- rowNames
        if (is.vector(gc)) {
            GC[rowNames] <- as.numeric(as.character(gc[rowNames]))
        } 
        else if (is.data.frame(gc) || is.matrix(gc)) {
            if (ncol(gc) == 2) {
                # We assume that the feature names are in the first column and 
                # the GC content in the second
                rownames(gc) <- gc[,1]
                GC[rowNames] <- as.numeric(as.character(gc[rowNames,2]))
            } 
            else if (ncol(gc) == 1) {
                # We assume that the GC contents are in the first column 
                # and the feature names in the rownames
                GC[rowNames] <- as.numeric(as.character(gc[rowNames,1]))
            } 
            else {
                stop( "The GC matrix/data.frame contains more columns ",
                    "than expected.")
            }
        }
        featureData(data)@data <- cbind(featureData(data)@data,GC)
    }
    
    # If exists biotype
    if (!is.null(biotype)) {
        Biotype <- rep(NA,length(rowNames))
        names(Biotype) <- rowNames
        if (is.vector(biotype)) {
            Biotype[rowNames] <- as.character(biotype[rowNames])
        } 
        else if (is.data.frame(biotype) || is.matrix(biotype)) {
            if (ncol(biotype) == 2) {
                # We assume that the feature names are in the first column and 
                # the biotypes in the second
                rownames(biotype) <- biotype[,1]
                Biotype[rowNames] <- as.character( biotype[rowNames,2])
            } 
            else if (ncol(biotype) == 1) {
                # We assume that the biotypes are in the first column and 
                # the feature names in the rownames
                Biotype[rowNames] <- as.character( biotype[rowNames,1]    )
            } 
            else {
                stop("The biotype matrix/data.frame contains more ",
                    "columns than expected.")
            }
        }        
        featureData(data)@data <- cbind(featureData(data)@data,Biotype)
        featureData(data)@data$Biotype <- 
            as.character(featureData(data)@data$Biotype)
    }
    
    # If exists chromosome
    if (!is.null(chromosome)) {
        Chromosome <- GeneStart <- GeneEnd <- rep(NA,length(rowNames))
        names(Chromosome) <- names(GeneStart) <- names(GeneEnd) <- rowNames
        Chromosome[rowNames] <- as.character(chromosome[rowNames,1])
        GeneStart[rowNames] <- as.numeric(as.character(chromosome[rowNames,2]))
        GeneEnd[rowNames] <- as.numeric(as.character(chromosome[rowNames,3]))
        featureData(data)@data <- cbind(featureData(data)@data,Chromosome,
            GeneStart,GeneEnd)
    }

    # If exists new factors
    if (!is.null(factors))
        phenoData(data)@data <- cbind(phenoData(data)@data,factors)
    
    data
}

################################################################################

################################################################################

# saturation.plot.R

#### SATURATION PLOTS

## Data for saturation plot (with or without biotypes)

saturation.dat <- function (input,k=0,biotypes=NULL,ndepth=6) {
    if (inherits(input,"eSet") == FALSE)
        stop("Error. You must give an eSet object\n")

    if (!is.null(assayData(input)$exprs))
        datos <- assayData(input)$exprs
    else
        datos <- assayData(input)$counts
    
    if (!is.null(featureData(input)$Biotype))
        infobio <- featureData(input)$Biotype        
    else
        infobio <- NULL

    nsam <- NCOL(datos)

    if (!is.null(infobio)) {        
        if(is.null(biotypes)) {
            biotypes <- unique(infobio)
            names(biotypes) <- biotypes                 
        }     
    } else { biotypes=NULL }     

    satura <- vector("list",length=length(biotypes)+1)
    names(satura) <- c("global",names(biotypes))

    ndepth1 <- ceiling(ndepth/2)

    datos <- round(datos,0)
    datos0 <- datos + 0.2
    
    # Random subsamples for each sample
    submuestras <- seq.depth <- vector("list",length=nsam)
    names(submuestras) <- names(seq.depth) <- colnames(datos)

    for (n in seq_len(nsam)) { # simulating subsamples for each sample
        total <- sum(datos[,n]) # total counts in sample n
        # simulation for each depth and real depth
        varias <- vector("list",length=ndepth+1)

        for (i in seq_len(ndepth1)) { # simulating depths < real depth
            muestra <- rmultinom(10,size=round(total/(ndepth1+1),0),
                prob=datos[,n])
            if (i == 1)
                varias[[i]] <- muestra 
            else
                varias[[i]] <- varias[[i-1]] + muestra
        }

        varias[[ndepth1+1]] <- as.matrix(datos[,n])

        for (i in (ndepth1+2):(ndepth+1)) { # simulating depths < real depth
            muestra <- rmultinom(10,size=round(total/(ndepth1+1),0),
                prob=datos0[,n])
            if (i == ndepth1+2)
                varias[[i]] <- matrix(varias[[i-1]],ncol=10,
                    nrow=nrow(varias[[i-1]])) + muestra
            else
                varias[[i]] <- varias[[i-1]] + muestra
        }

        submuestras[[n]] <- varias
        seq.depth[[n]] <- c(round(total/(ndepth1+1),0)*seq_len(ndepth1),total,
            round(total/(ndepth1+1),0)*((ndepth1+2):(ndepth+1)))
    }

    # Global saturation
    satura[[1]] <- vector("list",length=nsam)
    names(satura[[1]]) <- colnames(datos)

    for (n in seq_len(nsam)) { # for each sample
        satura[[1]][[n]] <- vapply(submuestras[[n]],function(x) { 
            mean(apply(x,2,noceros,k=k)) 
        },numeric(1))
    }

    # Per biotypes
    if (!is.null(infobio)) { # if biotypes available
        biog <- lapply(biotypes,function(x) { 
            which(is.element(infobio,x)) 
        })
        names(biog) <- names(biotypes)
        for (j in 2:length(satura)) { # for each biotype
            satura[[j]] <- vector("list",length=nsam)
            names(satura[[j]]) <- colnames(datos)
            for (n in seq_len(nsam)) { # for each sample
                conbio <- lapply(submuestras[[n]],function(x) { 
                    as.matrix(x[biog[[j-1]],]) 
                })
                satura[[j]][[n]] <- vapply(conbio,function(x) { 
                    mean(apply(x,2,noceros,k=k)) 
                },numeric(1))
            }
        }
    } 
    else
        biog <- NULL

    # computing detection increasing per million reads
    newdet <- vector("list",length=1+length(biotypes))
    names(newdet) <- c("global",names(biotypes))
    
    for (j in seq_len(length(newdet))) {
        newdet[[j]] <- vector("list",length=nsam)
        names(newdet[[j]]) <- colnames(datos)
        for (n in seq_len(nsam)) {
            puntos <- data.frame("x"=seq.depth[[n]],"y"=satura[[j]][[n]])
            pendi <- NULL
            for(i in 2:nrow(puntos))                
                pendi <- c(pendi,
                    (puntos$y[i]-puntos$y[i-1])/(puntos$x[i]-puntos$x[i-1]))
            newdet[[j]][[n]] <- c(NA,pendi*1000000)
        }
    }

    bionum <- c(NROW(datos),vapply(biog,length,integer(1)))
    names(bionum) <- c("global",names(biog))
    
    # Results at real sequencing depth
    real <- vector("list",length=length(satura))
    names(real) <- names(satura)
    realdepth <- vapply(seq.depth,function (x) x[ndepth1+1],numeric(1))/10^6
    for (i in seq_len(length(real))) {
        real[[i]] <- data.frame("depth"=realdepth,
            "detec"=vapply(satura[[i]],function (x) x[ndepth1+1],numeric(1)))
        rownames(real[[i]]) <- colnames(datos)
    }
    
    # Results
    satura <- list("saturation"=satura,"bionum"=bionum,"depth"=seq.depth,
        "newdet"=newdet,"real"=real)
    
    satura
}

################################################################################
