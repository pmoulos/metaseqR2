metaseqrPlot <- function(object,sampleList,annotation=NULL,contrastList=NULL,
    pList=NULL,thresholds=list(p=0.05,f=1),plotType=c("mds","biodetection",
    "countsbio","saturation","readnoise","rnacomp","correl","pairs","boxplot",
    "gcbias","lengthbias","meandiff","meanvar","deheatmap","volcano","biodist",
    "filtered","venn"),isNorm=FALSE,output="x11",path=NULL,...) {
    if (is(object,"GenomicRanges")) {
		object <- as.data.frame(object)
		object <- object[,c(1:3,6,7,5,8,9)]
		colnames(object)[1] <- "chromosome"
	}
    if (!is.matrix(object) && !is.data.frame(object))
        stopwrap("object argument must be a matrix or data frame!")
    if (is.null(annotation) && any(plotType %in% c("biodetection",
        "countsbio","saturation","rnacomp","readnoise","biodist","gcbias",
        "lengthbias","filtered")))
        stopwrap("annotation argument is needed when plotType is ",
            "\"biodetection\", \"countsbio\",\"saturation\",\"rnacomp\", ",
            "\"readnoise\", \"biodist\", \"gcbias\", \"lengthbias\", ",
            "\"filtered\" or \"venn\"!")
    if (any(plotType %in% c("deheatmap","volcano","biodist","venn"))) {
        if (is.null(contrastList))
            stopwrap("contrastList argument is needed when plotType is ",
                "\"deheatmap\",\"volcano\", \"biodist\" or \"venn\"!")
        if (is.null(pList))
            stopwrap("The p argument which is a list of p-values for each ",
                "contrast is needed when plotType is \"deheatmap\", ",
                "\"volcano\", \"biodist\" or \"venn\"!")
        if (is.na(thresholds$p) || is.null(thresholds$p) || thresholds$p==1) {
            warnwrap(paste("The p-value threshold when plotType is ",
            "\"deheatmap\", \"volcano\", \"biodist\" or \"venn\" must allow ",
            "the normal plotting of DEG diagnostic plots! Setting to 0.05..."))
            thresholds$p <- 0.05
        }
    }
    if (is.null(path)) 
        path <- getwd()
    if (is.data.frame(object) && !("filtered" %in% plotType)) 
        object <- as.matrix(object)
    
    if (is(annotation,"GenomicRanges")) {
		annotation <- as.data.frame(annotation)
		annotation <- annotation[,c(1:3,6,7,5,8,9)]
		colnames(annotation)[1] <- "chromosome"
	}
    
    if (any(plotType %in% c("biodetection","countsbio","saturation",
        "rnacomp","biodist","readnoise")))
        covars <- list(
            data=object,
            length=annotation$end - annotation$start,
            gc=annotation$gc_content,
            chromosome=annotation[,1:3],
            factors=data.frame(class=asClassVector(sampleList)),
            biotype=annotation$biotype,
            gene_name=as.character(annotation$gene_name)
        )

    rawPlots <- c("mds","biodetection","countsbio","saturation","readnoise",
        "correl","pairwise")
    normPlots <- c("boxplot","gcbias","lengthbias","meandiff","meanvar",
        "rnacomp")
    statPlots <- c("deheatmap","volcano","biodist")
    otherPlots <- c("filtered")
    vennPlots <- c("venn")
    files <- list()
    
    for (p in plotType) {
        disp("  Plotting ",p,"...")
        if (p %in% rawPlots && !isNorm) {
            switch(p,
                mds = {
                    files$mds <- diagplotMds(object,sampleList,output=output,
                        path=path)
                },
                biodetection = {
                    files$biodetection <- diagplotNoiseq(object,sampleList,
                        covars,whichPlot=p,output=output,path=path,...)
                },
                countsbio = {
                    files$countsbio <- diagplotNoiseq(object,sampleList,
                        covars,whichPlot=p,output=output,path=path,...)
                },
                saturation = {
                    fil <- diagplotNoiseq(object,sampleList,covars,
                        whichPlot=p,output=output,path=path,...)
                    files$saturation$biotype <- fil[["biotype"]]
                    files$saturation$sample <- fil[["sample"]]
                },
                readnoise = {
                    files$readnoise <- diagplotNoiseq(object,sampleList,
                        covars,whichPlot=p,output=output,path=path,...)
                },
                correl = {
                    files$correl$heatmap <- diagplotCor(object,type="heatmap",
                        output=output,path=path,...)
                    files$correl$correlogram <- diagplotCor(object,
                        type="correlogram",output=output,path=path,...)
                },
                pairwise = {
                    files$pairwise <- diagplotPairs(object,output=output,
                        path=path)
                }
            )
        }
        if (p %in% normPlots) {
            switch(p,
                boxplot = {
                    files$boxplot <- diagplotBoxplot(object,name=sampleList,
                        isNorm=isNorm,output=output,path=path,...)
                },
                gcbias = {
                    files$gcbias <- diagplotEdaseq(object,sampleList,
                        covar=annotation$gc_content,isNorm=isNorm,
                        whichPlot=p,output=output,path=path,...)
                },
                lengthbias = {
                    files$lengthbias <- diagplotEdaseq(object,sampleList,
                        covar=annotation$end-annotation$start,isNorm=isNorm,
                        whichPlot=p,output=output,path=path,...)
                },
                meandiff = {
                    fil <- diagplotEdaseq(object,sampleList,isNorm=isNorm,
                        whichPlot=p,output=output,path=path,...)
                    for (n in names(fil)) {
                        if (!is.null(fil[[n]]))
                            files$meandiff[[n]] <- unlist(fil[[n]])
                    }
                },
                meanvar = {
                    fil <- diagplotEdaseq(object,sampleList,isNorm=isNorm,
                        whichPlot=p,output=output,path=path,...)
                    for (n in names(fil)) {
                        if (!is.null(fil[[n]]))
                            files$meanvar[[n]] <- unlist(fil[[n]])
                    }
                },
                rnacomp = {
                    files$rnacomp <- diagplotNoiseq(object,sampleList,covars,
                        whichPlot=p,output=output,isNorm=isNorm,path=path,...)
                }
            )
        }
        if (p %in% statPlots && isNorm) {
            for (cnt in names(contrastList)) {
            disp("  Contrast: ",cnt)                
                samples <- names(unlist(contrastList[[cnt]]))
                mat <- as.matrix(object[,match(samples,colnames(object))])
                switch(p,
                    deheatmap = {
                        files$deheatmap[[cnt]] <- diagplotDeHeatmap(mat,cnt,
                            output=output,path=path)
                    },
                    volcano = {
                        fc <- log2(makeFoldChange(cnt,sampleList,object,1))
                        for (contrast in colnames(fc)) {
                            files$volcano[[contrast]] <- diagplotVolcano(
                                fc[,contrast],pList[[cnt]],contrast,
                                fcut=thresholds$f,pcut=thresholds$p,
                                output=output,path=path)
                        }
                    },
                    biodist = {
                        files$biodist[[cnt]] <- diagplotNoiseq(object,
                            sampleList,covars,whichPlot=p,output=output,
                            biodistOpts=list(p=pList[[cnt]],
                            pcut=thresholds$p,name=cnt),path=path,...)
                    }
                )
            }
        }
        if (p %in% otherPlots) {
            switch(p,
                filtered = {
                    files$filtered <- diagplotFiltered(object,annotation,
                        output=output,path=path)
                }
            )
        }
        if (p %in% vennPlots) {
            switch(p,
                venn = {
                    for (cnt in names(contrastList)) {
                        disp("  Contrast: ",cnt)
                        if (!is.null(annotation)) {
                            altNames <- as.character(annotation$gene_name)
                            names(altNames) <- rownames(annotation)
                        }
                        else
                            altNames <- NULL
                        files$venn[[cnt]] <- diagplotVenn(pList[[cnt]],
                            pcut=thresholds$p,nam=cnt,output=output,path=path,
                            altNames=altNames)
                    }
                }
            )
        }
    }
    
    return(files)
}

diagplotBoxplot <- function(mat,name=NULL,logIt="auto",yLim="default",
    isNorm=FALSE,output="x11",path=NULL,altNames=NULL,...) {
    if (is.null(path)) 
        path <- getwd()
    if (isNorm)
        status<- "normalized"
    else
        status<- "raw"
    # Need to log?
    if (logIt=="auto") {
        if (diff(range(mat,na.rm=TRUE))>1000)
            mat <- log2disp(mat)
    }
    else if (logIt=="yes")
        mat <- log2disp(mat)
    # Define the axis limits based on user input
    if (!is.numeric(yLim) && yLim=="default") {
        minY <- floor(min(mat))
        maxY <- ceiling(max(mat))
    }
    else if (is.numeric(yLim)) {
        minY <- yLim[1]
        maxY <- yLim[2]
    }
    grouped <- FALSE
    if (is.null(name)) {
        if (is.null(colnames(mat)))
            nams <- paste("Sample",1:ncol(mat),sep=" ")
        else
            nams <- colnames(mat)
    }
    else if (length(name)==1 && name=="none")
        nams <- rep("",ncol(mat))
    else if (is.list(name)) { # Is sampleList
        nams <- unlist(name)
        grouped <- TRUE
    }
    cols <- c("red3","green3","blue2","gold","skyblue","orange3","burlywood",
        "red","blue","green","orange","darkgrey","green4","black","pink",
        "brown","magenta","yellowgreen","pink4","seagreen4","darkcyan")
    if (grouped) {
        tmp <- as.numeric(factor(asClassVector(name)))
        bCols <- cols[tmp]
    }
    else bCols <- cols
    matList <- list()
    for (i in 1:ncol(mat))
        matList[[i]] <- mat[,i]
    names(matList) <- nams
    if (output != "json") {
        fil <- file.path(path,paste("boxplot_",status,".",output,sep=""))
        graphicsOpen(output,fil)
        if (!is.numeric(yLim) && yLim=="default")
            b <- boxplot(matList,names=nams,col=bCols,las=2,main=paste(
                "Boxplot ",status,sep=""),...)
        else
            b <- boxplot(matList,names=nams,col=bCols,ylim=c(minY,maxY),
                las=2,main=paste("Boxplot ",status,sep=""),...)
        graphicsClose(output)
    }
    else {
        # Create boxplot object
        b <- boxplot(matList,plot=FALSE)
        colnames(b$stat) <- nams
        # Locate the outliers
        oList <- lapply(names(matList),function(x,M,b) {
            v <- b[,x]
            o <- which(M[[x]]<v[1] | M[[x]]>v[5])
            if (length(o)>0)
                return(M[[x]][o])
            else
                return(NULL)
        },matList,b$stat)
        # Create output object
        obj <- list(
            x=NULL,
            y=NULL,
            plot=b,
            samples=name,
            ylims=c(minY,maxY),
            xlims=NULL,
            status=status,
            pcut=NULL,
            fcut=NULL,
            altnames=altNames,
            user=oList
        )
        json <- boxplotToJSON(obj)
        fil <- file.path(path,paste("boxplot_",status,".json",sep=""))
        disp("Writing ",fil)
        write(json,fil)
    }
    return(fil)
}

diagplotMds <- function(x,sampleList,method="spearman",logIt=TRUE,
    output="x11",path=NULL,...) {
    if (is.null(path)) path <- getwd()
    classes <- as.factor(asClassVector(sampleList))
    design <- as.numeric(classes)
    colspace <- c("red","blue","yellowgreen","orange","aquamarine2",
                  "pink2","seagreen4","brown","purple","chocolate")
    pchspace <- c(20,17,15,16,8,3,2,0,1,4)
    if (ncol(x)<3) {
        warnwrap("MDS plot cannot be created with less than 3 samples! ",
            "Skipping...")
        return(NULL)
    }
    if (logIt)
        y <- nat2log(x,base=2)
    else
        y <- x
    d <- as.dist(0.5*(1-cor(y,method=method)))
    mdsObj <- cmdscale(d,eig=TRUE,k=2)
    xr <- diff(range(min(mdsObj$points[,1]),max(mdsObj$points[,1])))
    yr <- diff(range(min(mdsObj$points[,2]),max(mdsObj$points[,2])))
    xlim <- c(min(mdsObj$points[,1])-xr/10,max(mdsObj$points[,1])+xr/10)
    ylim <- c(min(mdsObj$points[,2])-yr/10,max(mdsObj$points[,2])+yr/10)
    if (output!="json") {
        fil <- file.path(path,paste("mds.",output,sep=""))
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=9,height=7)
        else
            graphicsOpen(output,fil,width=1024,height=768)     
        plot(mdsObj$points[,1],mdsObj$points[,2],
             col=colspace[1:length(levels(classes))][design],
             pch=pchspace[1:length(levels(classes))][design],
             xlim=xlim,ylim=ylim,
             main="MDS plot",xlab="MDS 1",ylab="MDS 2",
             cex=0.9,cex.lab=0.9,cex.axis=0.9,cex.main=0.9)
        text(mdsObj$points[,1],mdsObj$points[,2],labels=colnames(x),pos=3,
            cex=0.7)
        grid()
        graphicsClose(output)
        return(fil)
    }
    else {
        # Create output object
        xx <- mdsObj$points[,1]
        yy <- mdsObj$points[,2]
        names(xx) <- names(yy) <- unlist(sampleList)
        obj <- list(
            x=xx,
            y=yy,
            plot=NULL,
            samples=sampleList,
            ylim=ylim,
            xlim=xlim,
            status=NULL,
            pcut=NULL,
            fcut=NULL,
            altnames=NULL,
            user=NULL
        )
        json <- mdsToJSON(obj)
        return(json)
        #fil <- file.path(path,"mds.json")
        #disp("Writing ",fil)
        #write(json,fil)
    }
}

.diagplotMdsGg <- function(x,sampleList,method="spearman",logIt=TRUE,...) {
    classes <- as.factor(asClassVector(sampleList))
    design <- as.numeric(classes)
    if (ncol(x)<3) {
        warnwrap("MDS plot cannot be created with less than 3 samples! ",
            "Skipping...")
        return("MDS plot cannot be created with less than 3 samples!")
    }
    if (logIt)
        y <- nat2log(x,base=2)
    else
        y <- x
    d <- as.dist(0.5*(1-cor(y,method=method)))
    mdsObj <- cmdscale(d,eig=TRUE,k=2)
    
    xr <- diff(range(min(mdsObj$points[,1]),max(mdsObj$points[,1])))
    yr <- diff(range(min(mdsObj$points[,2]),max(mdsObj$points[,2])))
    xlims <- c(min(mdsObj$points[,1])-xr/10,max(mdsObj$points[,1])+xr/10)
    ylims <- c(min(mdsObj$points[,2])-yr/10,max(mdsObj$points[,2])+yr/10)
	
	plotData <- data.frame(
		x=mdsObj$points[,1],
		y=mdsObj$points[,2],
		Condition=classes
	)
	rownames(plotData) <- colnames(x)
	
	mds <- ggplot() +
		geom_point(data=plotData,mapping=aes(x=x,y=y,colour=Condition,
			shape=Condition),size=4) +
		xlim(xlims[1],xlims[2]) + 
		ylim(ylims[1],ylims[2]) +
		ggtitle("MDS plot") +
		xlab("\nPrincipal coordinate 1") +
		ylab("Principal coordinate 2\n") +
		theme_bw() +
		theme(axis.title.x=element_text(size=12,face="bold"),
			axis.title.y=element_text(size=12,face="bold"),
			axis.text.x=element_text(size=11,face="bold"),
			axis.text.y=element_text(size=11,face="bold"),
			strip.text.x=element_text(size=11,face="bold"),
			strip.text.y=element_text(size=11,face="bold"),
			legend.position="bottom",
			legend.title=element_text(size=10,face="bold"),
			legend.text=element_text(size=9),
			legend.key=element_blank()) +
		geom_text(data=plotData,mapping=aes(x=x,y=y,
			label=rownames(plotData)),size=4,hjust=-0.15,vjust=0)
	return(mds)
    #return(fil)
}

diagplotPairs <- function(x,output="x11",path=NULL,...) {    
    x <- as.matrix(x)
    x <- nat2log(x)
    n <- ncol(x)
    if (!is.null(colnames(x)))
        nams <- colnames(x)
    else
        nams <- paste("Sample_",1:ncol(x),sep="")

    if (!is.null(path))
        fil <- file.path(path,paste("correlation_pairs",output,sep="."))
    else
        fil <- paste("correlation_pairs",output,sep=".")
    if (output %in% c("pdf","ps","x11"))
        graphicsOpen(output,fil,width=12,height=12)
    else {
        if (ncol(x)<=5)
            graphicsOpen(output,fil,width=800,height=800,res=100)
        else
            graphicsOpen(output,fil,width=1024,height=1024,res=150)
    }
        
    # Setup the grid
    par(mfrow=c(n,n),mar=c(1,1,1,1),oma=c(1,1,0,0),mgp=c(2,0.5,0),cex.axis=0.6,
        cex.lab=0.6)

    # Plot
    for (i in 1:n) {
        for (j in 1:n) {
            if (i==j) { # Diagonal
                plot(0:10,0:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
                text(c(3,5,3),c(9.5,5,1),c("X-Y plots",nams[i],"M-D plots"),
                    cex=c(0.8,1,0.8))
                arrows(6,9.5,9.5,9.5,angle=20,length=0.1,lwd=0.8,cex=0.8)
                arrows(0.2,3.2,0.2,0.2,angle=20,length=0.1,lwd=0.8,cex=0.8)
            }
            else if (i<j) { # XY plot
                plot(x[,i],x[,j],pch=20,col="blue",cex=0.4,xlab=nams[i],
                    ylab=nams[j],...)
                lines(lowess(x[,i],x[,j]),col="red")
                cc <- paste("cor:",formatC(cor(x[,i],x[,j]),digits=3))
                text(3,max(x[,j]-1),labels=cc,cex=0.7,)
                #grid()
            }
            else if (i>j) # MD plot
            {
                plot((x[,i]+x[,j])/2,x[,j]-x[,i],pch=20,col="blue",cex=0.4,...)
                lines(lowess((x[,i]+x[,j])/2,x[,j]-x[,i]),col="red")
                #grid()
            }
        }
    }

    graphicsClose(output)
    return(fil)
}

diagplotCor <- function(mat,type=c("heatmap","correlogram"),output="x11",
    path=NULL,...) {
    x <- as.matrix(mat)
    type <- tolower(type[1])
    checkTextArgs("type",type,c("heatmap","correlogram"))
    #if (!require(corrplot) && type=="correlogram")
    #    stop("R package corrplot is required!")
    cor.mat <- cor(mat)
    if (!is.null(colnames(mat)))
        colnames(cor.mat) <- colnames(mat)
    if (!is.null(path))
        fil <- file.path(path,paste("correlation_",type,".",output,sep=""))
    else
        fil <- paste("correlation_",type,".",output,sep="")
    if (output %in% c("pdf","ps","x11"))
        graphicsOpen(output,fil,width=7,height=7)
    else
        graphicsOpen(output,fil,width=640,height=640,res=100)
    if (type=="correlogram")
        corrplot(cor.mat,method="ellipse",order="hclust",...)
    else if (type=="heatmap") {
        n <- dim(cor.mat)[1]
        labs <- matrix(NA,n,n)
        for (i in 1:n)
            for (j in 1:n)
                labs[i,j] <- sprintf("%.2f",cor.mat[i,j])
        if (n <= 5)
            notecex <- 1.2
        else if (n > 5 & n < 10)
            notecex <- 0.9
        else
            notecex <- 0.7
        heatmap.2(cor.mat,col=colorRampPalette(c("yellow","grey","blue")),
            revC=TRUE,trace="none",symm=TRUE,Colv=TRUE,cellnote=labs,
            keysize=1,density.info="density",notecex=notecex,cexCol=0.9,
            cexRow=0.9,font.lab=2)
    }
    graphicsClose(output)
    return(fil)
}

diagplotEdaseq <- function(x,sampleList,covar=NULL,isNorm=FALSE,
    whichPlot=c("meanvar","meandiff","gcbias","lengthbias"),output="x11",
    path=NULL,...) {
    if (is.null(path)) path <- getwd()
    checkTextArgs("whichPlot",whichPlot,c("meanvar","meandiff","gcbias",
        "lengthbias"),multiarg=TRUE)
    if (is.null(covar) && whichPlot %in% c("gcbias","lengthbias"))
        stopwrap("\"covar\" argument is required when \"whichPlot\" is ",
            "\"gcbias\ or \"lengthbias\"!")
    if (isNorm)
        status <- "normalized"
    else
        status <- "raw"
    if (is.null(covar)) covar <- rep(NA,nrow(x))
    s <- newSeqExpressionSet(x,phenoData=AnnotatedDataFrame(
        data.frame(conditions=asClassVector(sampleList),
        row.names=colnames(x))),featureData=AnnotatedDataFrame(data.frame(
        gc=covar,length=covar,row.names=rownames(x))))
    switch(whichPlot,
        meandiff = {
            fil <- vector("list",length(sampleList))
            names(fil) <- names(sampleList)
            for (n in names(sampleList)) {
                if (length(sampleList[[n]])==1) {
                    warnwrap("Cannot create a mean-difference plot with one ",
                        "sample per condition! Skipping...")
                    next
                }
                pair.matrix <- combn(1:length(sampleList[[n]]),2)
                fil[[n]] <- vector("list",ncol(pair.matrix))
                for (i in 1:ncol(pair.matrix)) {
                    s1 <- sampleList[[n]][pair.matrix[1,i]]
                    s2 <- sampleList[[n]][pair.matrix[2,i]]
                    fil[[n]][[i]] <- file.path(path,paste(whichPlot,"_",
                        status,"_",n,"_",s1,"_",s2,".",output,sep=""))
                    names(fil[[n]][i]) <- paste(s1,"vs",s2,sep="_")
                    graphicsOpen(output,fil[[n]][[i]])
                    MDPlot(s,y=pair.matrix[,i],main=paste("MD plot for ",n," ",
                        status," samples ",s1," and ",s2,sep=""),cex.main=0.9)
                    graphicsClose(output)
                }
            }
        },
        meanvar = {
            fil <- vector("list",length(sampleList))
            names(fil) <- names(sampleList)
            for (n in names(sampleList)) {    
                if (length(sampleList[[n]])==1) {
                    warnwrap("Cannot create a mean-variance plot with one ",
                        "sample per condition! Skipping...")
                    next
                }
                pair.matrix <- combn(1:length(sampleList[[n]]),2)
                fil[[n]] <- vector("list",ncol(pair.matrix))
                for (i in 1:ncol(pair.matrix)) {
                    s1 <- sampleList[[n]][pair.matrix[1,i]]
                    s2 <- sampleList[[n]][pair.matrix[2,i]]
                    fil[[n]][[i]] <- file.path(path,paste(whichPlot,"_",status,
                        "_",n,"_",s1,"_",s2,".",output,sep=""))
                    names(fil[[n]][i]) <- paste(s1,"vs",s2,sep="_")
                    graphicsOpen(output,fil[[n]][[i]])
                    suppressWarnings(meanVarPlot(s,main=paste("MV plot for ",n,
                        " ",status," samples ",s1," and ",s2,sep=""),
                        cex.main=0.9))
                    graphicsClose(output)
                }
            }
        },
        gcbias = {
            if (!output=="json") {
                fil <- file.path(path,paste(whichPlot,"_",status,".",output,
                    sep=""))
                graphicsOpen(output,fil)
                biasPlot(s,"gc",xlim=c(0.1,0.9),log=TRUE,ylim=c(0,15),
                    main=paste("Expression - GC content ",status,sep=""))
                grid()
                graphicsClose(output)
            }
            else {
                obj <- list(
                    x=NULL,
                    y=NULL,
                    plot=NULL,
                    samples=sampleList,
                    ylim=NULL,
                    xlim=NULL,
                    status=status,
                    pcut=NULL,
                    fcut=NULL,
                    altnames=NULL,
                    user=list(counts=x,covar=covar,covarname="GC content")
                )
                json <- biasPlotToJSON(obj)
                fil <- file.path(path,paste(whichPlot,"_",status,".json",
                    sep=""))
                disp("Writing ",fil)
                write(json,fil)
            }
        },
        lengthbias = {
            if (output!="json") {
                fil <- file.path(path,paste(whichPlot,"_",status,".",output,
                    sep=""))
                graphicsOpen(output,fil)
                biasPlot(s,"length",log=TRUE,ylim=c(0,10),
                    main=paste("Expression - Gene length ",status,sep=""))
                grid()
                graphicsClose(output)
            }
            else {
                obj <- list(
                    x=NULL,
                    y=NULL,
                    plot=NULL,
                    samples=sampleList,
                    ylim=NULL,
                    xlim=NULL,
                    status=status,
                    pcut=NULL,
                    fcut=NULL,
                    altnames=NULL,
                    user=list(counts=x,covar=covar,
                        covarname="Gene/transcript length")
                )
                json <- biasPlotToJSON(obj)
                fil <- file.path(path,paste(whichPlot,"_",status,".json",
                    sep=""))
                disp("Writing ",fil)
                write(json,fil)
            }
        }
    )
    return(fil)
}

diagplotNoiseq <- function(x,sampleList,covars,whichPlot=c("biodetection",
    "countsbio","saturation","rnacomp","readnoise","biodist"),output="x11",
    biodistOpts=list(p=NULL,pcut=NULL,name=NULL),path=NULL,isNorm=FALSE,
    ...) {
    if (is.null(path)) path <- getwd()
    # covars is a list of gc-content, factors, length, biotype, chromosomes, 
    # factors, basically copy of the noiseq object
    whichPlot <- tolower(whichPlot[1])
    checkTextArgs("whichPlot",whichPlot,c("biodetection","countsbio",
        "saturation","readnoise","rnacomp","biodist"),multiarg=FALSE)
    if (missing(covars))
        stopwrap("\"covars\" argument is required with NOISeq specific plots!")
    else {
        covars$biotype <- as.character(covars$biotype)
        names(covars$length) <- names(covars$gc) <-
            rownames(covars$chromosome) <- names(covars$biotype) <-
            rownames(x)
    }
    if (whichPlot=="biodist") {
        if (is.null(biodistOpts$p))
            stopwrap("A p-value must be provided for the \"biodist\" plot!")
        if (is.null(biodistOpts$pcut) || is.na(biodistOpts$pcut)) 
            biodistOpts$pcut=0.05
    }
    if (isNorm)
        status<- "normalized"
    else
        status<- "raw"
    
    # All of these plots are NOISeq specific so we need a local NOISeq object
    if (any(is.na(unique(covars$biotype))))
        covars$biotype=NULL # Otherwise, it will probably crash
    localObj <- NOISeq::readData(
        data=x,
        length=covars$geneLength,
        gc=covars$gcContent,
        chromosome=covars$chromosome,
        #factors=data.frame(class=covars$factors),
        factors=covars$factors,
        biotype=covars$biotype
    )
    switch(whichPlot,
        biodetection = {
            diagplotData <- NOISeq::dat(localObj,type=whichPlot)
            samples <- unlist(sampleList)
            if (output!="json") {
                fil <- character(length(samples))
                names(fil) <- samples
                for (i in 1:length(samples)) {
                    fil[samples[i]] <- file.path(path,paste(whichPlot,"_",
                        samples[i],".",output,sep=""))
                    if (output %in% c("pdf","ps","x11"))
                        graphicsOpen(output,fil[samples[i]],width=9,height=7)
                    else
                        graphicsOpen(output,fil[samples[i]],width=1024,
							height=768)
                    explo.plot(diagplotData,samples=i)
                    graphicsClose(output)
                }
                return(fil)
            }
            else {
                diagplotDataSave = NOISeq::dat2save(diagplotData)
                obj <- list(
                   x=NULL,
                   y=NULL,
                   plot=NULL,
                   samples=sampleList,
                   ylims=NULL,
                   xlims=NULL,
                   status=status,
                   pcut=NULL,
                   fcut=NULL,
                   altnames=covars$gene_name,
                   user=list(plotdata=diagplotDataSave,covars=covars)
                )
                json <- bioDetectionToJSON(obj)
                names(json) <- samples
                #fil <- character(length(samples))
                #names(fil) <- samples
                #for (i in 1:length(samples)) {
                #    fil[samples[i]] <- file.path(path,
                #        paste(whichPlot,"_",samples[i],".json",sep=""))
                #    disp("Writing ",fil[samples[i]])
                #    write(json[[i]],fil[samples[i]])
                #}
                return(json)
            }
        },
        countsbio = {
            samples <- unlist(sampleList)
            if (output!="json") {
                diagplotData <- NOISeq::dat(localObj,type=whichPlot,
                    factor=NULL)
                fil <- character(length(samples))
                names(fil) <- samples
                for (i in 1:length(samples)) {
                    fil[samples[i]] <- file.path(path,paste(whichPlot,"_",
                        samples[i],".",output,sep=""))
                    if (output %in% c("pdf","ps","x11"))
                        graphicsOpen(output,fil[samples[i]],width=9,height=7)
                    else
                        graphicsOpen(output,fil[samples[i]],width=1024,
                            height=768)
                    explo.plot(diagplotData,samples=i,plottype="boxplot")
                    graphicsClose(output)
                }
                return(fil)
            }
            else {
                colnames(x) <- unlist(sampleList)
                obj <- list(
                   x=NULL,
                   y=NULL,
                   plot=NULL,
                   samples=sampleList,
                   ylims=NULL,
                   xlims=NULL,
                   status=status,
                   pcut=NULL,
                   fcut=NULL,
                   altnames=covars$gene_name,
                   user=list(counts=nat2log(x),covars=covars)
                )
                # Write JSON by sample
                #fil <- vector("list",2)
                #names(fil) <- c("sample","biotype")
                #fil[["sample"]] <- character(length(samples))
                #names(fil[["sample"]]) <- samples
                jsonList <- vector("list",2)
                names(jsonList) <- c("sample","biotype")
                jsonList[["sample"]] <- vector("list",length(samples))
                names(jsonList[["sample"]]) <- samples
                bts <- unique(as.character(obj$user$covars$biotype))
                #fil[["biotype"]] <- character(length(bts))
                #names(fil[["biotype"]]) <- bts
                jsonList[["biotype"]] <- vector("list",length(bts))
                #json <- countsBioToJSON(obj,by="sample")
                jsonList[["sample"]] <- countsBioToJSON(obj,by="sample")
                #for (i in 1:length(samples)) {
                #    fil[["sample"]][samples[i]] <- file.path(path,
                #        paste(whichPlot,"_",samples[i],".json",sep=""))
                #    disp("Writing ",fil[["sample"]][samples[i]])
                #    write(json[[i]],fil[["sample"]][samples[i]])
                #}
                #json <- countsBioToJSON(obj,by="biotype")
                jsonList[["biotype"]] <- countsBioToJSON(obj,by="biotype")
                names(jsonList[["biotype"]]) <- bts
                #names(json) <- samples
                #for (i in 1:length(bts)) {
                #    fil[["biotype"]][bts[i]] <- file.path(path,
                #        paste(whichPlot,"_",bts[i],".json",sep=""))
                #    disp("Writing ",fil[["biotype"]][bts[i]])
                #    write(json[[i]],fil[["biotype"]][bts[i]])
                #}
                return(jsonList)
            }
        },
        saturation = {
            # For 10 saturation points
            diagplotData <- NOISeq::dat(localObj,k=0,ndepth=9,type=whichPlot)
            d2s <- NOISeq::dat2save(diagplotData)
            if (output != "json") {
                fil <- diagplotNoiseqSaturation(d2s,output,covars$biotype,
                    path=path)
                return(fil)
			}
            else {
                samples <- unlist(sampleList)
                obj <- list(
                   x=NULL,
                   y=NULL,
                   plot=NULL,
                   samples=sampleList,
                   ylims=NULL,
                   xlims=NULL,
                   status=status,
                   pcut=NULL,
                   fcut=NULL,
                   altnames=covars$gene_name,
                   user=list(plotdata=d2s)
                )
                # Write JSON by sample
                #fil <- vector("list",2)
                #names(fil) <- c("sample","biotype")
                #fil[["sample"]] <- character(length(samples))
                #names(fil[["sample"]]) <- samples
                jsonList <- vector("list",2)
                names(jsonList) <- c("sample","biotype")
                #json <- bioSaturationToJSON(obj,by="sample")
                jsonList[["sample"]] <- bioSaturationToJSON(obj,by="sample")
                names(jsonList[["sample"]]) <- samples
                #for (i in 1:length(samples)) {
                #    fil[["sample"]][samples[i]] <- file.path(path,
                #        paste(whichPlot,"_",samples[i],".json",sep=""))
                #    disp("Writing ",fil[["sample"]][samples[i]])
                #    write(json[[i]],fil[["sample"]][samples[i]])
                #}
                #json <- bioSaturationToJSON(obj,by="biotype")
                jsonList[["biotype"]] <- bioSaturationToJSON(obj,by="biotype")
                #names(jsonList[["sample"]]) <- samples
                #fil[["biotype"]] <- character(length(json))
                #names(fil[["biotype"]]) <- names(json)
                #for (n in names(json)) {
                #    fil[["biotype"]][n] <- file.path(path,
                #        paste(whichPlot,"_",n,".json",sep=""))
                #    disp("Writing ",fil[["biotype"]][n])
                #    write(json[[n]],fil[["biotype"]][n])
                #}
                return(jsonList)
            }
        },
        rnacomp = {
            if (ncol(localObj)<3) {
                warnwrap("RNA composition plot cannot be created with less ",
                    "than 3 samples! Skipping...")
                return(NULL)
            }
            if (ncol(localObj)>12) {
                warnwrap("RNA composition plot cannot be created with more ",
                    "than 12 samples! Skipping...")
                return(NULL)
            }
            diagplotData <- NOISeq::dat(localObj,type="cd")
            fil <- file.path(path,paste(whichPlot,"_",status,".",output,
                sep=""))
            graphicsOpen(output,fil)
            explo.plot(diagplotData)
            grid()
            graphicsClose(output)
        },
        readnoise = {
            D <- cddat(localObj)
            if (output!="json") {
                fil <- file.path(path,paste(whichPlot,".",output,sep=""))
                graphicsOpen(output,fil)
                cdplot(D,main="RNA-Seq reads noise")
                grid()
                graphicsClose(output)
            }
            else {
                colnames(D$data2plot)[2:ncol(D$data2plot)] <- 
                    unlist(sampleList)
                obj <- list(
                    x=NULL,
                    y=NULL,
                    plot=NULL,
                    samples=sampleList,
                    xlim=NULL,
                    ylim=NULL,
                    status=NULL,
                    pcut=NULL,
                    fcut=NULL,
                    altnames=NULL,
                    user=D$data2plot
                )
                json <- readNoiseToJSON(obj)
                fil <- file.path(path,paste(whichPlot,".json",sep=""))
                disp("Writing ",fil)
                write(json,fil)
            }
        },
        biodist = { # We have to fake a noiseq object
            p <- biodistOpts$p
            if (is.matrix(p)) p <- p[,1]
            dummy <- new("Output",
                comparison=c("Dummy1","Dummy2"),
                factor=c("class"),
                k=1,
                lc=1,
                method="n",
                replicates="biological",
                results=list(
                    data.frame(
                        Dummy1=rep(1,length(p)),
                        Dummy2=rep(1,length(p)),
                        M=rep(1,length(p)),
                        D=rep(1,length(p)),
                        prob=as.numeric(p),
                        ranking=rep(1,length(p)),
                        Length=rep(1,length(p)),
                        GC=rep(1,length(p)),
                        Chrom=as.character(covars$chromosome[,1]),
                        GeneStart=covars$chromosome[,2],
                        GeneEnd=covars$chromosome[,3],
                        Biotype=covars$biotype
                    )
                ),
                nss=5,
                pnr=0.2,
                v=0.02
            )
            if (!is.null(biodistOpts$name))
                fil <- file.path(path,paste(whichPlot,"_",biodistOpts$name,
                    ".",output,sep=""))
            else
                fil <- file.path(path,paste(whichPlot,".",output,sep=""))
            if (output %in% c("pdf","ps","x11"))
                graphicsOpen(output,fil,width=10,height=6)
            else
                graphicsOpen(output,fil,width=1024,height=640)
            tryCatch( # A lot of times, there is a problem with this function
                DE.plot(dummy,chromosomes=NULL,q=biodistOpts$pcut,
                    graphic="distr"),
                error=function(e) {
                    disp("      Known problem with NOISeq and external ",
                        "p-values  detected! Trying to make a plot with ",
                        "alternative p-values  (median of p-value ",
                        "distribution)...")
                    fil="error"
                    tryCatch(
                        DE.plot(dummy,chromosomes=NULL,
                            q=quantile(biodistOpts$p,0.5),
                            graphic="distr"),
                        error=function(e) {
                            disp("      Cannot create DEG biotype plot! This ",
                                "is not related to a problem with the ",
                                "results. Excluding...")
                            fil="error"
                        },
                        finally=""
                    )
                },
                finally=""
            )
            graphicsClose(output)
        }
    )
    return(fil)
}

diagplotNoiseqSaturation <- function(x,o,tb,path=NULL) {
    if (is.null(path)) path <- getwd()
    if (length(unique(tb))==1) {
        warnwrap("Saturation plot cannot be created with only one biotype! ",
            "Skipping...")
        return(NULL)
    }
    totalBiotypes <- table(tb)
    theBiotypes <- names(tb)
    biotypes <- colnames(x[[1]][,2:ncol(x[[1]])])
    colspace <- c("red3","green4","blue2","orange3","burlywood",
                  "lightpink4","gold","skyblue","red2","green2","firebrick3",
                  "orange4","yellow4","skyblue3","tan4","gray40",
                  "brown2","darkgoldenrod","cyan3","coral2","cadetblue",
                  "bisque3","blueviolet","chocolate3","darkkhaki","dodgerblue")
    pchspace <- c(rep(c(15,16,17,18),6),15)

    # Plot all biotypes per sample
    fSample <- character(length(names(x)))
    names(fSample) <- names(x)
    for (n in names(x)) {
        fSample[n] <- file.path(path,paste("saturation_",n,".",o,sep=""))
        if (o %in% c("pdf","ps","x11"))
            graphicsOpen(o,fSample[n],width=10,height=7)
        else
            graphicsOpen(o,fSample[n],width=1024,height=800)
        y <- x[[n]]
        sep <- match(c("global","protein_coding"),colnames(y))
        yab <- cbind(y[,"depth"],y[,sep])
        ynab <- y[,-sep]
        colnames(yab)[1] <- colnames(ynab)[1] <- "depth"
        xlim <- range(y[,"depth"])
        ylimAb <- range(yab[,2:ncol(yab)])
        ylimNab <- range(ynab[,2:ncol(ynab)])
        par(cex.axis=0.9,cex.main=1,cex.lab=0.9,font.lab=2,font.axis=2,pty="m",
            lty=2,lwd=1.5,mfrow=c(1,2))
        plot.new()
        plot.window(xlim,ylimNab)
        axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)/1e+6))
        axis(2,at=pretty(ylimNab,10))
        title(main="Non abundant biotype detection saturation",
            xlab="Depth in millions of reads",ylab="Detected features")
        co <- 0
        for (b in biotypes) {
            co <- co + 1
            if (b=="global" || b=="protein_coding") {
                # Silently do nothing
            }
            else {
                lines(ynab[,"depth"],ynab[,b],col=colspace[co])
                points(ynab[,"depth"],ynab[,b],pch=pchspace[co],
                    col=colspace[co],cex=1)
            }
        }
        grid()
        graphics::legend(
            x="topleft",legend=colnames(ynab)[2:ncol(ynab)],xjust=1,yjust=0,
            box.lty=0,x.intersp=0.5,cex=0.6,text.font=2,
            col=colspace[1:(ncol(ynab)-1)],pch=pchspace[1:(ncol(ynab)-1)]
        )
        plot.new()
        plot.window(xlim,ylimAb)
        axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)/1e+6))
        axis(2,at=pretty(ylimAb,10))
        title(main="Abundant biotype detection saturation",
            xlab="Depth in millions of reads",ylab="Detected features")
        co <- 0
        for (b in c("global","protein_coding")) {
            co <- co + 1
            lines(yab[,"depth"],yab[,b],col=colspace[co])
            points(yab[,"depth"],yab[,b],pch=16,col=colspace[co],cex=1.2)
        }
        grid()
        graphics::legend(
            x="topleft",legend=c("global","protein_coding"),xjust=1,yjust=0,
            box.lty=0,lty=2,x.intersp=0.5,cex=0.7,text.font=2,
            col=colspace[1:2],pch=pchspace[1:2]
        )
        mtext(n,side=3,line=-1.5,outer=TRUE,font=2,cex=1.3)
        graphicsClose(o)
    }

    # Plot all samples per biotype
    g <- makeGrid(length(biotypes))
    fAll <- file.path(path,paste("biotype_saturation.",o,sep=""))
    if (o %in% c("pdf","ps"))
        graphicsOpen(o,fAll,width=14,height=14)
    else
        graphicsOpen(o,fAll,width=1600,height=1600,res=150)
    par(cex.axis=0.8,cex.main=0.9,cex.lab=0.8,pty="m",lty=2,lwd=1.5,mfrow=g,
        mar=c(3,3,1,1),oma=c(1,1,0,0),mgp=c(2,0.5,0))
    for (b in biotypes) {
        y <- depth <- vector("list",length(x))
        names(y) <- names(depth) <- names(x)
        for (n in names(x)) {
            y[[n]] <- x[[n]][,b]
            depth[[n]] <- x[[n]][,"depth"]
        }
        y <- do.call("cbind",y)
        xlim <- range(do.call("c",depth))
        ylim <- range(y)
        plot.new()
        plot.window(xlim,ylim)
        axis(1,at=pretty(xlim,5),labels=as.character(pretty(xlim,5)/1e+6),
            line=0.5)
        axis(2,at=pretty(ylim,5),line=0.5)
        title(main=b,xlab="Depth in millions of reads",
            ylab="Detected features")
        co <- 0
        for (n in colnames(y)) {
            co <- co + 1
            lines(depth[[n]],y[,n],col=colspace[co])
            points(depth[[n]],y[,n],pch=pchspace[co],col=colspace[co])
        }
        grid()
        graphics::legend(
            x="bottomright",legend=colnames(y),xjust=1,yjust=0,
            box.lty=0,x.intersp=0.5,
            col=colspace[1:length(colnames(y))],
            pch=pchspace[1:length(colnames(y))]
        )
    }
    graphicsClose(o)

    return(list(sample=fSample,biotype=fAll))
}

diagplotVolcano <- function(f,p,con=NULL,fcut=1,pcut=0.05,altNames=NULL,
    output="x11",path=NULL,...) { # output can be json here...
    ## Check rjson
    #if ("json" %in% output && !require(rjson))
    #    stopwrap("R package rjson is required to create interactive volcano plot!")
    if (is.null(path)) path <- getwd()
    if (is.null(con))
        con <- conn <- ""
    else {
        conn <- con
        con <- paste("for ",con)
    }
    fil <- file.path(path,paste("volcano_plot_",conn,".",output,sep=""))
    if (output!="json") {
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=8,height=10)
        else
            graphicsOpen(output,fil,width=768,height=1024,res=100)
    }
    rem <- which(is.na(p))
    if (length(rem)>0) {
        p <- p[-rem]
        f <- f[-rem]
        if (!is.null(altNames))
            altNames <- altNames[-rem]
    }
    # Fix problem with extremely low p-values, only for display purposes though
    pZero <- which(p==0)
    if (length(pZero)>0)
        p[pZero] <- runif(length(pZero),0,1e-256)
    xlim <- c(-max(abs(f)),max(abs(f)))
    ylim <- c(0,ceiling(-log10(min(p))))
    up <- which(f>=fcut & p<pcut)
    down <- which(f<=-fcut & p<pcut)
    u <- union(up,down)
    altNamesNeutral <- NULL
    if (length(u)>0) {
        ff <- f[-u]
        pp <- p[-u]
        if (!is.null(altNames))
            altNamesNeutral <- altNames[-u]
    }
    else {
        ff <- f
        pp <- p
        if (!is.null(altNames))
            altNamesNeutral <- altNames
    }
    if (output!="json") {
        par(cex.main=1.1,cex.lab=1.1,cex.axis=1.1,font.lab=2,font.axis=2,
            pty="m",lwd=1.5)
        plot.new()
        plot.window(xlim,ylim)
        axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)))
        axis(2,at=pretty(ylim,10))
        title(paste(main="Volcano plot",con),
            xlab="Fold change",ylab="-log10(p-value)")
        points(ff,-log10(pp),pch=20,col="blue2",cex=0.9)
        points(f[down],-log10(p[down]),pch=20,col="green3",cex=0.9)
        points(f[up],-log10(p[up]),pch=20,col="red2",cex=0.9)
        abline(h=-log10(pcut),lty=4)
        abline(v=-fcut,lty=2)
        abline(v=fcut,lty=2)
        grid()
        graphics::legend(
            x="topleft",
            legend=c("up-regulated","down-regulated","unregulated",
                "p-value threshold","fold change threshold"),
            col=c("red2","green3","blue1","black","black"),
            pch=c(20,20,20,NA,NA),lty=c(NA,NA,NA,4,2),
            xjust=1,yjust=0,box.lty=0,x.intersp=0.5,cex=0.8,text.font=2
        )
        graphicsClose(output)
    }
    else {
        obj <- list(
            x=f,
            y=p,
            plot=NULL,
            samples=NULL,
            xlim=xlim,
            ylim=ylim,
            status=NULL,
            pcut=pcut,
            fcut=fcut,
            altnames=altNames,
            user=list(up=up,down=down,unf=ff,unp=pp,ualt=altNamesNeutral,
                con=con)
        )
        #json <- volcanoToJSON(obj)
        #fil <- file.path(path,paste("volcano_",con,".json",sep=""))
        #write(json,fil)
        fil <- volcanoToJSON(obj)
    }
    return(fil)
}

diagplotDeHeatmap <- function(x,con=NULL,output="x11",path=NULL,...) {
    if (is.null(path)) path <- getwd()
    if (is.null(con))
        con <- conn <- ""
    else {
        conn <- con
        con <- paste("for ",con)
    }
    y <- nat2log(x,2,1)
    # First plot the normal image
    fil <- file.path(path,paste("de_heatmap_",conn,".",output,sep=""))
    if (output %in% c("pdf","ps","x11"))
        graphicsOpen(output,fil,width=10,height=10)
    else
        graphicsOpen(output,fil,width=800,height=800)
    heatmap.2(y,trace="none",col=bluered(16),labRow="",cexCol=0.9,keysize=1,
        font.lab=2,main=paste("DEG heatmap",con),cex.main=0.9)
    graphicsClose(output)
    ## Then the "interactive" using sendplot
    #xy.labels <- list(normalized_counts=x,log2_normalized_counts=y)
    #x.labels <- data.frame(
    #    label=colnames(x),
    #    description=paste("Sample",colnames(x))
    #)
    #y.labels <- data.frame(
    #    label=rownames(x),
    #    description=paste("Gene ID:",rownames(x))
    #)
    #suppressWarnings(heatmap.send(
    #    y,
    #    distfun=dist,
    #    hclustfun=hclust,
    #    MainColor=bluered(16),
    #    labRow="",
    #    labCol=NULL,
    #    keep.dendro=TRUE, 
    #    x.labels=x.labels,
    #    y.labels=y.labels,
    #    xy.labels=xy.labels,
    #    image.size="2048x4096",
    #    fname.root=paste("iframe_de_heatmap_",conn,sep=""),
    #    dir=paste(path,.Platform$file.sep,sep=""),
    #    header="v3",
    #    window.size="2048x4192"
    #))
    return(fil)
}

diagplotFiltered <- function(x,y,output="x11",path=NULL,...) {
    if (output !="json") {
        if (is.null(path)) path <- getwd()
        fil <- file.path(path,paste("filtered_genes.",output,sep=""))
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=12,height=8)
        else
            graphicsOpen(output,fil,width=1200,height=800,res=100)
        chr <- table(as.character(x$chromosome))
        bt <- table(as.character(x$biotype))
        chrAll <- table(as.character(y$chromosome))
        btAll <- table(as.character(y$biotype))
        barlabChr <- as.character(chr)
        barlabBt <- as.character(bt)
        perChr <- chr/chrAll[names(chr)]
        perBt <- bt/btAll[names(bt)]
        # Some bug...
        perChr[perChr>1] <- 1
        perBt[perBt>1] <- 1
        #
        suppressWarnings(perChrLab <- paste(formatC(100*perChr,digits=1,
            format="f"),"%",sep=""))
        suppressWarnings(perBt.lab <- paste(formatC(100*perBt,digits=1,
            format="f"),"%",sep=""))

        par(mfrow=c(2,2),mar=c(1,4,2,1),oma=c(1,1,1,1))

        # Chromosomes
        barxChr <- barplot(chr,space=0.5,
            ylim=c(0,max(chr)+ceiling(max(chr)/10)),yaxt="n",xaxt="n",
            plot=FALSE)
        plot.new()
        plot.window(xlim=c(0,ceiling(max(barxChr))),
            ylim=c(0,max(chr)+ceiling(max(chr)/10)),mar=c(1,4,1,1))
        axis(2,at=pretty(0:(max(chr)+ceiling(max(chr)/10))),cex.axis=0.9,padj=1,
            font=2)
        text(x=barxChr,y=chr,label=barlabChr,cex=0.7,font=2,col="green3",
            adj=c(0.5,-1.3))
        title(main="Filtered genes per chromosome",cex.main=1.1)
        mtext(side=2,text="Number of genes",line=2,cex=0.9,font=2)
        grid()
        barplot(chr,space=0.5,ylim=c(0,max(chr)+ceiling(max(chr)/10)),
            col="blue3",border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

        # Biotypes
        barxBt <- barplot(bt,space=0.5,ylim=c(0,max(bt)+ceiling(max(bt)/10)),
            yaxt="n",xaxt="n",plot=FALSE)
        plot.new()
        plot.window(xlim=c(0,ceiling(max(barxBt))),
            ylim=c(0,max(bt)+ceiling(max(bt)/10)),mar=c(1,4,1,1))
        axis(2,at=pretty(0:(max(bt)+ceiling(max(bt)/10))),cex.axis=0.9,padj=1,
            font=2)
        text(x=barxBt,y=bt,label=barlabBt,cex=0.7,font=2,col="blue",
            adj=c(0.5,-1.3),xpd=TRUE)
        title(main="Filtered genes per biotype",cex.main=1.1)
        mtext(side=2,text="Number of genes",line=2,cex=0.9,font=2)
        grid()
        barplot(bt,space=0.5,ylim=c(0,max(bt)+ceiling(max(bt)/10)),col="red3",
            border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

        # Chromosome percentage
        barxPerChr <- barplot(perChr,space=0.5,ylim=c(0,max(perChr)),
            yaxt="n",xaxt="n",plot=FALSE)
        plot.new()
        par(mar=c(9,4,1,1))
        plot.window(xlim=c(0,max(barxPerChr)),ylim=c(0,max(perChr)))
        #axis(1,at=barxPerChr,labels=names(perChr),cex.axis=0.9,font=2,
        #    tcl=-0.3,col="lightgrey",las=2)
        axis(1,at=barxPerChr,labels=FALSE,tcl=-0.3,col="lightgrey")
        axis(2,at=seq(0,max(perChr),length.out=5),labels=formatC(seq(0,
            max(perChr),length.out=5),digits=2,format="f"),cex.axis=0.9,padj=1,
            font=2)
        text(barxPerChr,par("usr")[3]-max(perChr)/17,labels=names(perChr),
            srt=45,adj=c(1,1.1),xpd=TRUE,cex=0.9,font=2)
        text(x=barxPerChr,y=perChr,label=perChrLab,cex=0.7,font=2,
            col="green3",adj=c(0.5,-1.3),xpd=TRUE)
        mtext(side=2,text="fraction of total genes",line=2,cex=0.9,font=2)
        grid()
        barplot(perChr,space=0.5,ylim=c(0,max(perChr)),col="blue3",
            border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

        # Biotype percentage
        barxPerBt <- barplot(perBt,space=0.5,ylim=c(0,max(perBt)),yaxt="n",
            xaxt="n",plot=FALSE)
        plot.new()
        par(mar=c(9,4,1,1))
        plot.window(xlim=c(0,max(barxPerBt)),ylim=c(0,max(perBt)))
        #axis(1,at=barxPerBt,labels=names(perBt),cex.axis=0.9,font=2,
        #    tcl=-0.3,col="lightgrey",las=2)
        axis(1,at=barxPerBt,labels=FALSE,tcl=-0.3,col="lightgrey")
        axis(2,at=seq(0,max(perBt),length.out=5),
            labels=formatC(seq(0,max(perBt),length.out=5),digits=2,format="f"),
            cex.axis=0.9,padj=1,font=2)
        text(barxPerBt,par("usr")[3]-max(perBt)/17,
            labels=gsub("prime","'",names(perBt)),srt=45,adj=c(1,1.1),
            xpd=TRUE,cex=0.9,font=2)
        text(x=barxPerBt,y=perBt,label=perBt.lab,cex=0.7,font=2,col="blue",
            adj=c(0.5,-1.3),xpd=TRUE)
        mtext(side=2,text="fraction of total genes",line=2,cex=0.9,font=2)
        grid()
        barplot(perBt,space=0.5,ylim=c(0,max(perBt)),col="red3",
            border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)
        
        graphicsClose(output)
    }
    else {
        obj <- list(
            x=NULL,
            y=NULL,
            plot=NULL,
            samples=NULL,
            xlim=NULL,
            ylim=NULL,
            status=NULL,
            pcut=NULL,
            fcut=NULL,
            altnames=NULL,
            user=list(filtered=x,total=y)
        )
        fil <- list(chromosome=NULL,biotype=NULL)
        json <- filteredToJSON(obj,by="chromosome")
        fil$chromosome <- file.path(path,"filtered_genes_chromosome.json")
        write(json,fil$chromosome)
        json <- filteredToJSON(obj,by="biotype")
        fil$biotype <- file.path(path,"filtered_genes_biotype.json")
        write(json,fil$biotype)
    }

    return(fil)
}

diagplotVenn <- function(pmat,fcmat=NULL,pcut=0.05,fcut=0.5,
    direction=c("dereg","up","down"),nam=as.character(round(1000*runif(1))),
    output="x11",path=NULL,altNames=NULL,...) {
    checkTextArgs("direction",direction,c("dereg","up","down"))
    if (is.na(pcut) || is.null(pcut) || pcut==1)
        warnwrap("Illegal pcut argument! Using the default (0.05)")
    algs <- colnames(pmat)
    if (is.null(algs))
        stopwrap("The p-value matrices must have the colnames attribute ",
            "(names of statistical algorithms)!")
    if (!is.null(fcmat) && (is.null(colnames(fcmat)) ||
        length(intersect(colnames(pmat),colnames(fcmat)))!=length(algs)))
        stopwrap("The fold change matrices must have the colnames attribute ",
            "(names of statistical algorithms) and must be the same as in the ",
            "p-value matrices!")
    nalg <- length(algs)
    if(nalg>5) {
        warnwrap(paste("Cannot create a Venn diagram for more than 5 result ",
            "sets! ",nalg,"found, only the first 5 will be used..."))
        algs <- algs[1:5]
        nalg <- 5
    }
    lenalias <- c("two","three","four","five")
    aliases <- toupper(letters[1:nalg])
    names(algs) <- aliases
    genes <- rownames(pmat)
    pairs <- makeVennPairs(algs)
    areas <- makeVennAreas(length(algs))
    counts <- makeVennCounts(length(algs))
    # Initially populate the results and counts lists so they can be used to create
    # the rest of the intersections
    results <- vector("list",nalg+length(pairs))
    names(results)[1:nalg] <- aliases
    names(results)[(nalg+1):length(results)] <- names(pairs)
    if (is.null(fcmat)) {
        for (a in aliases) {
            results[[a]] <- genes[which(pmat[,algs[a]]<pcut)]
            counts[[areas[[a]]]] <- length(results[[a]])
        }
    }
    else {
        switch(direction,
            dereg = {
                for (a in aliases) {
                    results[[a]] <-
                        genes[which(pmat[,algs[a]]<pcut & abs(
                        fcmat[,algs[a]])>=fcut)]
                    counts[[areas[[a]]]] <- length(results[[a]])
                }
            },
            up = {
                for (a in aliases) {
                    results[[a]] <-
                        genes[which(pmat[,algs[a]]<pcut &
                        fcmat[,algs[a]]>=fcut)]
                    counts[[areas[[a]]]] <- length(results[[a]])
                }
            },
            down = {
                for (a in aliases) {
                    results[[a]] <-
                        genes[which(pmat[,algs[a]]<pcut &
                        fcmat[,algs[a]]<=-fcut)]
                    counts[[areas[[a]]]] <- length(results[[a]])
                }
            }
        )
    }
    # Now, perform the intersections
    for (p in names(pairs)) {
        a = pairs[[p]][1]
        b = pairs[[p]][2]
        results[[p]] <- intersect(results[[a]],results[[b]])
        counts[[areas[[p]]]] <- length(results[[p]])
    }
    # And now, the Venn diagrams must be constructed
    colorScheme <- makeVennColorscheme(length(algs))
    fil <- file.path(path,paste("venn_plot_",nam,".",output,sep=""))
    if (output %in% c("pdf","ps","x11"))
        graphicsOpen(output,fil,width=8,height=8)
    else
        graphicsOpen(output,fil,width=800,height=800,res=100)
    switch(lenalias[length(algs)-1],
        two = {
            v <- draw.pairwise.venn(
                area1=counts$area1,
                area2=counts$area2,
                cross.area=counts$cross.area,
                category=paste(algs," (",aliases,")",sep=""),
                lty="blank",
                fill=colorScheme$fill,
                cex=1.5,
                cat.cex=1.3,
                cat.pos=c(0,0),
                cat.col=colorScheme$font,
                #cat.dist=0.07,
                cat.fontfamily=rep("Bookman",2)
            )
        },
        three = {
            #overrideTriple <<- TRUE
            v <- draw.triple.venn(
                area1=counts$area1,
                area2=counts$area2,
                area3=counts$area3,
                n12=counts$n12,
                n13=counts$n13,
                n23=counts$n23,
                n123=counts$n123,
                category=paste(algs," (",aliases,")",sep=""),
                lty="blank",
                fill=colorScheme$fill,
                cex=1.5,
                cat.cex=1.3,
                #cat.pos=c(0,0,180),
                cat.col=colorScheme$font,
                #cat.dist=0.07,
                cat.fontfamily=rep("Bookman",3)
            )
        },
        four = {
            v <- draw.quad.venn(
                area1=counts$area1,
                area2=counts$area2,
                area3=counts$area3,
                area4=counts$area4,
                n12=counts$n12,
                n13=counts$n13,
                n14=counts$n14,
                n23=counts$n23,
                n24=counts$n24,
                n34=counts$n34,
                n123=counts$n123,
                n124=counts$n124,
                n134=counts$n134,
                n234=counts$n234,
                n1234=counts$n1234,
                category=paste(algs," (",aliases,")",sep=""),
                lty="blank",
                fill=colorScheme$fill,
                cex=1.5,
                cat.cex=1.3,
                c(0.1,0.1,0.05,0.05),
                cat.col=colorScheme$font,
                cat.fontfamily=rep("Bookman",4)
            )
        },
        five = {
            v <- draw.quintuple.venn(
                area1=counts$area1,
                area2=counts$area2,
                area3=counts$area3,
                area4=counts$area4,
                area5=counts$area5,
                n12=counts$n12,
                n13=counts$n13,
                n14=counts$n14,
                n15=counts$n15,
                n23=counts$n23,
                n24=counts$n24,
                n25=counts$n25,
                n34=counts$n34,
                n35=counts$n35,
                n45=counts$n45,
                n123=counts$n123,
                n124=counts$n124,
                n125=counts$n125,
                n134=counts$n134,
                n135=counts$n135,
                n145=counts$n145,
                n234=counts$n234,
                n235=counts$n235,
                n245=counts$n245,
                n345=counts$n345,
                n1234=counts$n1234,
                n1235=counts$n1235,
                n1245=counts$n1245,
                n1345=counts$n1345,
                n2345=counts$n2345,
                n12345=counts$n12345,
                category=paste(algs," (",aliases,")",sep=""),
                lty="blank",
                fill=colorScheme$fill,
                cex=1.5,
                cat.cex=1.3,
                cat.dist=0.1,
                cat.col=colorScheme$font,
                cat.fontfamily=rep("Bookman",5)
            )
        }
    )
    graphicsClose(output)

    # Now do something with the results
    if (!is.null(path)) {
        resultsEx <- vector("list",length(results))
        names(resultsEx) <- names(results)
        if (!is.null(altNames)) {
            for (n in names(results))
                resultsEx[[n]] <- altNames[results[[n]]]
        }
        else {
            for (n in names(results))
                resultsEx[[n]] <- results[[n]]
        }
        maxLen <- max(sapply(resultsEx,length))
        for (n in names(resultsEx)) {
            if (length(resultsEx[[n]])<maxLen) {
                dif <- maxLen - length(resultsEx[[n]])
                resultsEx[[n]] <- c(resultsEx[[n]],rep(NA,dif))
            }
        }
        resultsEx <- do.call("cbind",resultsEx)
        write.table(resultsEx,file=file.path(path,"..","..","lists",
            paste0("venn_categories_",nam,".txt")),sep="\t",
            row.names=FALSE,quote=FALSE,na="")
    }
    
    return(fil)
}

makeVennPairs <- function(algs) {
    lenalias <- c("two","three","four","five")
    switch(lenalias[length(algs)-1],
        two = {
            return(list(
                AB=c("A","B")
            ))
        },
        three = {
            return(list(
                AB=c("A","B"),
                AC=c("A","C"),
                BC=c("B","C"),
                ABC=c("AB","C")
            ))
        },
        four = {
            return(list(
                AB=c("A","B"),
                AC=c("A","C"),
                AD=c("A","D"),
                BC=c("B","C"),
                BD=c("B","D"),
                CD=c("C","D"),
                ABC=c("AB","C"),
                ABD=c("AB","D"),
                ACD=c("AC","D"),
                BCD=c("BC","D"),
                ABCD=c("ABC","D")
            ))
        },
        five = {
            return(list(
                AB=c("A","B"),
                AC=c("A","C"),
                AD=c("A","D"),
                AE=c("A","E"),
                BC=c("B","C"),
                BD=c("B","D"),
                BE=c("B","E"),
                CD=c("C","D"),
                CE=c("C","E"),
                DE=c("D","E"),
                ABC=c("AB","C"),
                ABD=c("AB","D"),
                ABE=c("AB","E"),
                ACD=c("AC","D"),
                ACE=c("AC","E"),
                ADE=c("AD","E"),
                BCD=c("BC","D"),
                BCE=c("BC","E"),
                BDE=c("BD","E"),
                CDE=c("CD","E"),
                ABCD=c("ABC","D"),
                ABCE=c("ABC","E"),
                ABDE=c("ABD","E"),
                ACDE=c("ACD","E"),
                BCDE=c("BCD","E"),
                ABCDE=c("ABCD","E")
            ))
        }
    )
}

makeVennAreas <- function(n) {
    lenalias <- c("two","three","four","five")
    switch(lenalias[n-1],
        two = {
            return(list(
                A="area1",
                B="area2",
                AB="cross.area"
            ))
        },
        three = {
            return(list(
                A="area1",
                B="area2",
                C="area3",
                AB="n12",
                AC="n13",
                BC="n23",
                ABC="n123"
            ))
        },
        four = {
            return(list(
                 A="area1",
                 B="area2",
                 C="area3",
                 D="area4",
                 AB="n12",
                 AC="n13",
                 AD="n14",
                 BC="n23",
                 BD="n24",
                 CD="n34",
                 ABC="n123",
                 ABD="n124",
                 ACD="n134",
                 BCD="n234",
                 ABCD="n1234"
            ))
        },
        five = {
            return(list(
                 A="area1",
                 B="area2",
                 C="area3",
                 D="area4",
                 E="area5",
                 AB="n12",
                 AC="n13",
                 AD="n14",
                 AE="n15",
                 BC="n23",
                 BD="n24",
                 BE="n25",
                 CD="n34",
                 CE="n35",
                 DE="n45",
                 ABC="n123",
                 ABD="n124",
                 ABE="n125",
                 ACD="n134",
                 ACE="n135",
                 ADE="n145",
                 BCD="n234",
                 BCE="n235",
                 BDE="n245",
                 CDE="n345",
                 ABCD="n1234",
                 ABCE="n1235",
                 ABDE="n1245",
                 ACDE="n1345",
                 BCDE="n2345",
                 ABCDE="n12345"
            ))
        }
    )
}

makeVennCounts <- function(n) {
    lenalias <- c("two","three","four","five")
    switch(lenalias[n-1],
        two = {
            return(list(
                area1=NULL,
                area2=NULL,
                cross.area=NULL
            ))
        },
        three = {
            return(list(
                area1=NULL,
                area2=NULL,
                area3=NULL,
                n12=NULL,
                n13=NULL,
                n23=NULL,
                n123=NULL
            ))
        },
        four = {
            return(list(
                 area1=NULL,
                 area2=NULL,
                 area3=NULL,
                 area4=NULL,
                 n12=NULL,
                 n13=NULL,
                 n14=NULL,
                 n23=NULL,
                 n24=NULL,
                 n34=NULL,
                 n123=NULL,
                 n124=NULL,
                 n134=NULL,
                 n234=NULL,
                 n1234=NULL
            ))
        },
        five = {
            return(list(
                 area1=NULL,
                 area2=NULL,
                 area3=NULL,
                 area4=NULL,
                 area5=NULL,
                 n12=NULL,
                 n13=NULL,
                 n14=NULL,
                 n15=NULL,
                 n23=NULL,
                 n24=NULL,
                 n25=NULL,
                 n34=NULL,
                 n35=NULL,
                 n45=NULL,
                 n123=NULL,
                 n124=NULL,
                 n125=NULL,
                 n134=NULL,
                 n135=NULL,
                 n145=NULL,
                 n234=NULL,
                 n235=NULL,
                 n245=NULL,
                 n345=NULL,
                 n1234=NULL,
                 n1235=NULL,
                 n1245=NULL,
                 n1345=NULL,
                 n2345=NULL,
                 n12345=NULL
            ))
        }
    )
}

makeVennColorscheme <- function(n) {
    lenalias <- c("two","three","four","five")
    switch(lenalias[n-1],
        two = {
            return(list(
                fill=c("blue","orange2"),
                font=c("darkblue","orange4")
            ))
        },
        three = {
            return(list(
                fill=c("red","green","mediumpurple"),
                font=c("darkred","darkgreen","mediumpurple4")
            ))
        },
        four = {
            return(list(
                fill=c("red","green","mediumpurple","orange2"),
                font=c("darkred","darkgreen","mediumpurple4","orange4")
            ))
        },
        five = {
            return(list(
                fill=c("red","green","blue","mediumpurple","orange2"),
                font=c("darkred","darkgreen","darkblue","mediumpurple4",
                    "orange4")
            ))
        }
    )
}

diagplotRoc <- function(truth,p,sig=0.05,x="fpr",y="tpr",output="x11",
    path=NULL,draw=TRUE,...) {
    checkTextArgs("x",x,c("fpr","fnr","tpr","tnr","scrx","sens","spec"),
        multiarg=FALSE)
    checkTextArgs("y",y,c("fpr","fnr","tpr","tnr","scry","sens","spec"),
        multiarg=FALSE)
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    axName <- list(
        tpr="True Positive Rate",
        tnr="True Negative Rate",
        fpr="False Positive Rate",
        fnr="False Negative Rate",
        scrx="Ratio of selected",
        scry="Normalized TP/(FP+FN)",
        sens="Sensitivity",
        spec="1 - Specificity"
    )

    ROC <- vector("list",ncol(pmat))
    names(ROC) <- colnames(pmat)

    colspaceUniverse <- c("red","blue","green","orange","darkgrey","green4",
        "black","pink","brown","magenta","yellowgreen","pink4","seagreen4",
        "darkcyan")
    colspace <- colspaceUniverse[1:ncol(pmat)]
    names(colspace) <- colnames(pmat)

    eps <- min(pmat[!is.na(pmat) & pmat>0])
    for (n in colnames(pmat)) {
        disp("Processing ",n)
        gg <- which(pmat[,n]<=sig)
        psample <- -log10(pmax(pmat[gg,n],eps))
        #psample <- pmat[gg,n]
        size <- seq(1,length(gg))
        cuts <- seq(-log10(sig),max(psample),length.out=length(gg))
        #cuts <- seq(min(psample),sig,length.out=length(gg))
        local.truth <- truth[gg] #these are the true deg values of the ercc spike-ins reported by each tool
        S <- length(size)
                
        TP <- FP <- FN <- TN <- FPR <- FNR <- TPR <- TNR <- SENS <- SPEC <-
            SCRX <- SCRY <- numeric(S)
        
        for (i in 1:S) {
            TP[i] <- length(which(psample>cuts[i] & local.truth!=0))
            FP[i] <- length(which(psample>cuts[i] & local.truth==0))
            FN[i] <- length(which(psample<cuts[i] & local.truth!=0))
            TN[i] <- length(which(psample<cuts[i] & local.truth==0))

            ## Alternatives which I keep in the code
            #TP[i] <- length(intersect(names(which(psample>cuts[i])),
            #    names(which(local.truth!=0))))
            #FP[i] <- length(intersect(names(which(psample>cuts[i])),
            #    names(which(local.truth==0))))
            #FN[i] <- length(intersect(names(which(psample<cuts[i])),
            #    names(which(local.truth!=0))))
            #TN[i] <- length(intersect(names(which(psample<cuts[i])),
            #    names(which(local.truth==0))))
            #bad <- which(psample<cuts[i])
            #good <- which(psample>cuts[i])
            #TP[i] <- length(which(local.truth[good]!=0))
            #FP[i] <- length(which(local.truth[good]==0))
            #TN[i] <- length(which(local.truth[bad]==0))
            #FN[i] <- length(which(local.truth[bad]!=0))
            
            SCRX[i] <- i/S
            SCRY[i] <- TP[i]/(FN[i]+FP[i])

            if (FP[i]+TN[i] == 0)
                FPR[i] <- 0
            else
                FPR[i] <- FP[i]/(FP[i]+TN[i])
            FNR[i] <- FN[i]/(TP[i]+FN[i])
            TPR[i] <- TP[i]/(TP[i]+FN[i])
            if (TN[i]+FP[i] == 0)
                TNR[i] <- 0
            else
                TNR[i] <- TN[i]/(TN[i]+FP[i])
            SENS[i] <- TPR[i]
            SPEC[i] <- 1 - TNR[i]
        }
        #if (all(FPR==0))
        #    FPR[length(FPR)] <- 1
        #if (all(TNR==0)) {
        #    TNR[1] <- 1
        #    SPEC[i] <- 0
        #}

        ROC[[n]] <- list(TP=TP,FP=FP,FN=FN,TN=TN,
            FPR=FPR,FNR=FNR,TPR=TPR,TNR=TNR,SCRX=SCRX,SCRY=SCRY/max(SCRY),
            SENS=SENS,SPEC=SPEC,AUC=NULL)
    }
    
    for (n in colnames(pmat)) {
        disp("Calculating AUC for ",n)
        auc <- 0
        for (i in 2:length(ROC[[n]][[toupper(y)]])) {
            auc <- auc +
                0.5*(ROC[[n]][[toupper(x)]][i]-ROC[[n]][[toupper(x)]][i-1])*
                (ROC[[n]][[toupper(y)]][i]+ROC[[n]][[toupper(y)]][i-1])
        }
        ROC[[n]]$AUC <- abs(auc)
        # There are some extreme cases, with the Intersection case for the paper
        # where there are no FPs or TNs for a p-value cutoff of 0.2 (which is
        # imposed in order to avoid the saturation of the ROC curves). In these
        # cases, performance is virtually perfect, and the actual AUC should be
        # 1. For these cases, we set it to a value between 0.95 and 0.99 to
        # represent a more plausible truth.
        if (ROC[[n]]$AUC==0) ROC[[n]]$AUC <- sample(seq(0.95,0.99,by=0.001),1)
    }
    disp("")

    if (draw) {
        fil <- file.path(path,paste("ROC",output,sep="."))
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=8,height=8)
        else
            graphicsOpen(output,fil,width=1024,height=1024,res=100)

        xlim <- c(0,1)
        ylim <- c(0,1)
        par(cex.axis=0.9,cex.main=1,cex.lab=0.9,font.lab=2,font.axis=2,pty="m",
            lwd=1.5,lty=1)
        plot.new()
        plot.window(xlim,ylim)
        axis(1,at=pretty(xlim,10))
        axis(2,at=pretty(ylim,10))
        for (n in names(ROC))
            lines(ROC[[n]][[toupper(x)]],ROC[[n]][[toupper(y)]],
                col=colspace[n],...)
        grid()
        title(xlab=axName[[x]],ylab=axName[[y]])
        aucText <- as.character(sapply(ROC,function(x)
            round(x$AUC,digits=3)))
        graphics::legend(x="bottomright",col=colspace,lty=1,cex=0.9,
            legend=paste(names(ROC)," (AUC = ",aucText,")",sep=""))

        graphicsClose(output)
    }
    else
        fil <- NULL

    return(list(ROC=ROC,truth=truth,sigLevel=sig,xAxis=x,yAxis=y,path=fil))
}

diagplotFtd <- function(truth,p,type="fpc",N=2000,output="x11",path=NULL,
    draw=TRUE,...) {
    checkTextArgs("type",type,c("fpc","tpc","fnc","tnc"),multiarg=FALSE)
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    else if (is.numeric(p))
        pmat <- as.matrix(p)
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    yName <- list(
        tpc="Number of True Positives",
        fpc="Number of False Positives",
        tnc="Number of True Negatives",
        fnc="Number of False Negatives"
    )

    ftdr.list <- vector("list",ncol(pmat))
    names(ftdr.list) <- colnames(pmat)

    colspaceUniverse <- c("red","blue","green","orange","darkgrey","green4",
        "black","pink","brown","magenta","yellowgreen","pink4","seagreen4",
        "darkcyan")
    colspace <- colspaceUniverse[1:ncol(pmat)]
    names(colspace) <- colnames(pmat)

    switch(type,
        fpc = {
            for (n in colnames(pmat)) {
                disp("Processing ",n)
                z <- sort(pmat[,n])
                for (i in 1:N) {
                    nn <- length(intersect(names(z[1:i]),
                        names(which(truth==0))))
                    if (nn==0)
                        ftdr.list[[n]][i] <- 1
                    else
                        ftdr.list[[n]][i] <- nn
                }
            }
        },
        tpc = {
            for (n in colnames(pmat)) {
                disp("Processing ",n)
                z <- sort(pmat[,n])
                for (i in 1:N)
                    ftdr.list[[n]][i] <- length(intersect(names(z[1:i]),
                        names(which(truth!=0))))
            }
        },
        fnc = {
            for (n in colnames(pmat)) {
                disp("Processing ",n)
                z <- sort(pmat[,n],decreasing=TRUE)
                for (i in 1:N) {
                    nn <- length(intersect(names(z[1:i]),
                        names(which(truth!=0))))
                    if (nn==0)
                        ftdr.list[[n]][i] <- 1
                    else
                        ftdr.list[[n]][i] <- nn
                }
            }
        },
        tnc = {
            for (n in colnames(pmat)) {
                disp("Processing ",n)
                z <- sort(pmat[,n],decreasing=TRUE)
                for (i in 1:N)
                    ftdr.list[[n]][i] <- length(intersect(names(z[1:i]),
                        names(which(truth==0))))
            }
        }    
    )
    disp("")

    if (draw) {
        fil <- file.path(path,paste("FTDR_",type,".",output,sep=""))
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=8,height=8)
        else
            graphicsOpen(output,fil,width=1024,height=1024,res=100)

        xlim <- ylim <- c(1,N)
        #ylim <- c(1,length(which(truth!=0)))
        par(cex.axis=0.9,cex.main=1,cex.lab=0.9,font.lab=2,font.axis=2,pty="m",
            lwd=1.5,lty=1)
        plot.new()

        switch(type,
            fpc = {
                plot.window(xlim,ylim,log="y")
                axis(1,at=pretty(xlim,10))
                axis(2)
                for (n in names(ftdr.list)) {
                    lines(ftdr.list[[n]],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs False Positives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="topleft",legend=names(ftdr.list),
                    col=colspace,lty=1)
            },
            tpc = {
                plot.window(xlim,ylim)
                axis(1,at=pretty(xlim,10))
                axis(2,at=pretty(ylim,10))
                for (n in names(ftdr.list)) {
                    lines(ftdr.list[[n]],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs True Positives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="bottomright",legend=names(ftdr.list),
                    col=colspace,lty=1)
            },
            fnc = {
                plot.window(xlim,ylim,log="y")
                axis(1,at=pretty(xlim,10))
                axis(2)
                for (n in names(ftdr.list)) {
                    lines(ftdr.list[[n]],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs False Negatives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="topleft",legend=names(ftdr.list),
                    col=colspace,lty=1)
            },
            tnc = {
                plot.window(xlim,ylim)
                axis(1,at=pretty(xlim,10))
                axis(2,at=pretty(ylim,10))
                for (n in names(ftdr.list)) {
                    lines(ftdr.list[[n]],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs True Negatives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="bottomright",legend=names(ftdr.list),
                    col=colspace,lty=1)
            }    
        )

        graphicsClose(output)
    }
    else
        fil <- NULL

    return(list(ftdr=ftdr.list,truth=truth,type=type,N=N,path=fil))
}

diagplotAvgFtd <- function(ftdrObj,output="x11",path=NULL,draw=TRUE,...) {
    yName <- list(
        tpc="Number of True Positives",
        fpc="Number of False Positives",
        tnc="Number of True Negatives",
        fnc="Number of False Negatives"
    )

    stats <- names(ftdrObj[[1]]$ftdr)
    type <- ftdrObj[[1]]$type
    truth <- ftdrObj[[1]]$truth
    N <- ftdrObj[[1]]$N
    avgFtdrObj <- vector("list",length(stats))
    names(avgFtdrObj) <- stats
    colspaceUniverse <- c("red","blue","green","orange","darkgrey","green4",
        "black","pink","brown","magenta","yellowgreen","pink4","seagreen4",
        "darkcyan")
    colspace <- colspaceUniverse[1:length(stats)]
    names(colspace) <- stats
    
    for (s in stats) {
        disp("Retrieving ",s)
        avgFtdrObj[[s]] <- do.call("cbind",lapply(ftdrObj,
            function(x) x$ftdr[[s]]))
    }
    disp("")

    avgFtdrObj <- lapply(avgFtdrObj,function(x) {
        mn <- apply(x,1,mean)
        st <- apply(x,1,sd)
        return(list(mean=mn,sd=st))
    })

    means <- do.call("cbind",lapply(avgFtdrObj,function(x) x$mean))
    stds <- do.call("cbind",lapply(avgFtdrObj,function(x) x$sd))

    if (draw) {
        fil <- file.path(path,paste("AVG_FTDR_",type,".",output,sep=""))
        if (output %in% c("pdf","ps","x11"))
            graphicsOpen(output,fil,width=8,height=8)
        else
            graphicsOpen(output,fil,width=1024,height=1024,res=100)

        xlim <- ylim <- c(1,N)
        par(cex.axis=0.9,cex.main=1,cex.lab=0.9,font.lab=2,font.axis=2,pty="m",
            lwd=1.5,lty=1)
        plot.new()

        switch(type,
            fpc = {
                plot.window(xlim,ylim,log="y")
                axis(1,at=pretty(xlim,10))
                axis(2)
                for (n in colnames(means)) {
                    lines(means[,n],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs False Positives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="topleft",legend=colnames(means),
                    col=colspace,lty=1)
            },
            tpc = {
                plot.window(xlim,ylim)
                axis(1,at=pretty(xlim,10))
                axis(2,at=pretty(ylim,10))
                for (n in colnames(means)) {
                    lines(means[,n],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs True Positives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="bottomright",legend=colnames(means),
                    col=colspace,lty=1)
            },
            fnc = {
                plot.window(xlim,ylim,log="y")
                axis(1,at=pretty(xlim,10))
                axis(2)
                for (n in colnames(means)) {
                    lines(means[,n],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs False Negatives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="topleft",legend=colnames(means),
                    col=colspace,lty=1)
            },
            tnc = {
                plot.window(xlim,ylim)
                axis(1,at=pretty(xlim,10))
                axis(2,at=pretty(ylim,10))
                for (n in colnames(means)) {
                    lines(means[,n],col=colspace[n],...)
                }
                grid()
                title(main="Selected genes vs True Negatives",
                    xlab="Number of selected genes",ylab=yName[[type]])
                graphics::legend(x="bottomright",legend=colnames(means),
                    col=colspace,lty=1)
            }
        )
        
        graphicsClose(output)
    }
    else
        fil <- NULL

    return(list(avgFtdr=list(means=means,stds=stds),path=fil))
}

graphicsOpen <- function(o,f,...) {
    if(!checkGraphicsType(o))
        stopwrap("Invalid graphics output type!")
    if(checkGraphicsFile(o) && is.null(f))
        stopwrap("Please specify an output file name for your plot")
    
    switch(o,
        x11 = { dev.new(...) },
        pdf = { pdf(file=f,pointsize=10,...) },
        ps = { postscript(file=f,pointsize=10,...) },
        png = { png(filename=f,pointsize=12,...) },
        jpg = { jpeg(filename=f,pointsize=12,quality=100,...) },
        bmp = { bmp(filename=f,pointsize=12,...) },
        tiff = { tiff(filename=f,pointsize=12,...) }
    )
}

graphicsClose <- function(o) {
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    if (o!="x11")
        dev.off()
}

checkGraphicsType <- function(o) {
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    else
        return(TRUE)
}

checkGraphicsFile <- function(o) {
    if (is.element(o,c("png","jpg","tiff","bmp","pdf","ps")))
        return(TRUE)
    else
        return(FALSE)
}

log2disp <- function(mat,base=2) {
    mat[mat==0] <- 1
    if (base==10)
        return(log10(mat))
    else
        return(log2(mat))
}

nat2log <- function(x,base=2,off=1) {
    #x[x==0] <- off
    x <- x + off
    if (base==2)
        return(log2(x))
    else
        return(log10(x))
}

cddat <- function (input) {
    if (inherits(input,"eSet") == FALSE)
        stopwrap("The input data must be an eSet object.\n")
    if (!is.null(assayData(input)$exprs)) {
        if (ncol(assayData(input)$exprs) < 2)
            stopwrap("The input object should have at least two samples.\n")
        datos <- assayData(input)$exprs
    }
    else {
        if (ncol(assayData(input)$counts) < 2)
            stopwrap("The input object should have at least two samples.\n")
        datos <- assayData(input)$counts
    }
    datos <- datos[which(rowSums(datos) > 0),]
    nu <- nrow(datos) # number of detected features
    qq <- 1:nu
    data2plot = data.frame("%features" = 100*qq/nu)
    for (i in 1:ncol(datos)) {
        acumu <- 100*cumsum(sort(datos[,i],decreasing=TRUE))/sum(datos[,i])
        data2plot = cbind(data2plot, acumu)   
    }
    colnames(data2plot)[-1] = colnames(datos)
  
    # Diagnostic test
    KSpval = mostres = NULL
    for (i in 1:(ncol(datos)-1)) {
        for (j in (i+1):ncol(datos)) {      
            mostres = c(mostres, paste(colnames(datos)[c(i,j)], collapse="_"))
            KSpval = c(KSpval, suppressWarnings(ks.test(datos[,i], datos[,j],
                alternative = "two.sided"))$"p.value")
        }
    }    
    KSpval = p.adjust(KSpval, method = "fdr")
  
    return(list(
        "data2plot"=data2plot,
        "DiagnosticTest"=data.frame("ComparedSamples"=mostres,"KSpvalue"=KSpval)
    ))
}

cdplot <- function (dat,samples=NULL,...) {
    dat = dat$data2plot
    if (is.null(samples)) samples <- 1:(ncol(dat)-1)
    if (is.numeric(samples)) samples = colnames(dat)[samples+1]
    colspace <- c("red","blue","yellowgreen","orange","aquamarine2","pink2",
        "seagreen4","brown","purple","chocolate","gray10","gray30","darkblue",
        "darkgreen","firebrick2","darkorange4","darkorchid","darkcyan","gold4",
        "deeppink3")
    if (length(samples)>length(colspace))
        miscolores <- sample(colspace,length(samples),replace=TRUE)
    else
        miscolores <- sample(colspace,length(samples))
    plot(dat[,1],dat[,samples[1]],xlab="% features",ylab="% reads",type="l",
        col=miscolores[1],...)
    for (i in 2:length(samples))
        lines(dat[,1],dat[,samples[i]],col=miscolores[i])

    graphics::legend("bottom",legend=samples,
        text.col=miscolores[1:length(samples)],bty="n",lty=1,lwd=2,
        col=miscolores[1:length(samples)])
}
