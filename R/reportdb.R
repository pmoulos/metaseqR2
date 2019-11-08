.makeReportEnv <- function(e) {
	re <- new.env(parent=globalenv())
	
    re$offlineReport <- e$offlineReport
    re$reportDb <- e$reportDb
    re$org <- e$org
    re$refdb <- e$refdb
    re$utrOpts <- e$utrOpts
    re$qcPlots <- e$qcPlots
    re$sampleList <- e$sampleList
    re$geneData <- e$geneData
    re$geneCounts <- e$geneCounts
    re$geneDataExpr <- e$geneDataExpr
    re$normGenesExpr <- e$normGenesExpr
    re$fromRaw <- e$fromRaw
    re$transLevel <- e$transLevel
    re$countType <- e$countType
    re$fileType <- e$fileType
    re$geneFilters <- e$geneFilters
    re$exonFilters <- e$exonFilters
    re$normalization <- e$normalization
    re$statistics <- e$statistics
    re$metaP <- e$metaP
    re$pcut <- e$pcut
    re$contrast <- e$contrast
    re$adjustMethod <- e$adjustMethod
    re$sumpList <- e$sumpList
    re$logOffset <- e$logOffset
    re$countsName <- e$countsName
    re$excludeList <- e$excludeList
    re$libsizeList <- e$libsizeList
    re$reportMessages <- e$reportMessages
    re$preset <- e$preset
    re$whenApplyFilter <- e$whenApplyFilter
    re$normArgs <- e$normArgs
    re$statArgs <- e$statArgs
    re$figFormat <- e$figFormat
    re$exportWhere <- e$exportWhere
    re$exportWhat <- e$exportWhat
    re$exportScale <- e$exportScale
    re$exportValues <- e$exportValues
    re$exportStats <- e$exportStats
    re$execTime <- e$execTime
    re$theZeros <- e$theZeros
    re$theDead <- e$theDead
    re$exonFilterResult <- e$exonFilterResult
    re$geneFilterResult <- e$geneFilterResult
    re$geneFilterCutoff <- e$geneFilterCutoff
    re$runLog <- e$runLog
    re$geneDataFiltered <- e$geneDataFiltered
    re$totalGeneData <- e$totalGeneData
    re$normGenes <- e$normGenes
    re$geneCountsZero <- e$geneCountsZero
    re$contrastList <- e$contrastList
    re$reportTop <- e$reportTop
    re$cpList <- e$cpList
    re$reportTables <- e$reportTables
    re$exportCountsTable <- e$exportCountsTable
    
    re$FUN_CALL <- e$FUN_CALL
    re$PROJECT_PATH <- e$PROJECT_PATH
    
    return(re)
}

.createSqlPlotDb <- function(obj) {
	sampleList <- obj$sampleList
	contrastList <- obj$contrastList
	qcPlots <- obj$qcPlots
	geneCounts <- obj$geneCounts
	geneData <- obj$geneData
	geneDataExpr <- obj$geneDataExpr
	geneDataFiltered <- obj$geneDataFiltered
	totalGeneData <- obj$totalGeneData
	normGenes <- obj$normGenes
	normGenesExpr <- obj$normGenesExpr
	cpList <- obj$cpList
	sumpList <- obj$sumpList
	pcut <- obj$pcut
	metaP <- obj$metaP
	PROJECT_PATH <- obj$PROJECT_PATH
	
	samples <- unlist(sampleList)
	nsa <- length(samples)
	
	covarsRaw <- covarsStat <- NULL
	if (any(qcPlots %in% c("biodetection","countsbio","saturation",
		"rnacomp","readnoise"))) {
		covarsRaw <- list(
			data=geneCounts,
			length=width(geneData),
			gc=as.numeric(geneData$gc_content),
			chromosome=data.frame(
				chromosome=as.character(seqnames(geneData)),
				start=start(geneData),
				end=end(geneData)
			),
			factors=data.frame(class=asClassVector(sampleList)),
			biotype=as.character(geneData$biotype),
			gene_name=as.character(geneData$gene_name)
		)
	}

	if ("biodist" %in% qcPlots) {
		covarsStat <- list(
			data=normGenesExpr,
			length=width(geneDataExpr),
			gc=as.numeric(geneDataExpr$gc_content),
			chromosome=data.frame(
				chromosome=as.character(seqnames(geneDataExpr)),
				start=start(geneDataExpr),
				end=end(geneDataExpr)
			),
			factors=data.frame(class=asClassVector(sampleList)),
			biotype=as.character(geneDataExpr$biotype),
			gene_name=as.character(geneDataExpr$gene_name)
		)
	}

	rdb <- file.path(PROJECT_PATH$data,"reportdb.sqlite")
	disp("Opening plot database in ",rdb)
	con <- dbConnect(dbDriver("SQLite"),dbname=rdb)
	rs <- .initReportDbTables(con)

	if ("mds" %in% qcPlots) {
		disp("Importing mds...")
		json <- diagplotMds(geneCounts,sampleList,output="json")
		.dbImportPlot(con,"MDS","mds","generic",json)
	}
	if ("biodetection" %in% qcPlots) {
		disp("  Importing biodetection...")
		cap <- capture.output({
			json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
				whichPlot="biodetection",output="json")
		})
		for (s in samples) {
			disp("    ",s)
			name <- paste("biodetection",s,sep="_")
			.dbImportPlot(con,name,"biodetection","generic",json[[s]])
		}
	}
	if ("countsbio" %in% qcPlots) {
		disp("  Importing countsbio...")
		cap <- capture.output({
			jsonList <- diagplotNoiseq(geneCounts,sampleList,
				covars=covarsRaw,whichPlot="countsbio",output="json")
		})
		for (s in samples) {
			disp("    ",s)
			name <- paste("countsbio",s,sep="_")
			.dbImportPlot(con,name,"countsbio","sample",
				jsonList[["sample"]][[s]])
		}
		for (b in names(jsonList[["biotype"]])) {
			disp("    ",b)
			name <- paste("countsbio",b,sep="_")
			.dbImportPlot(con,name,"countsbio","biotype",
				jsonList[["biotype"]][[b]])
		}
	}
	if ("saturation" %in% qcPlots) {
		disp("  Importing saturation...")
		jsonList <- diagplotNoiseq(geneCounts,sampleList,
			covars=covarsRaw,whichPlot="saturation",output="json")
		for (s in unlist(sampleList)) {
			disp("    ",s)
			name <- paste("saturation",s,sep="_")
			.dbImportPlot(con,name,"saturation","sample",
				jsonList[["sample"]][[s]])
		}
		for (b in names(jsonList[["biotype"]])) {
			disp("    ",b)
			name <- paste("saturation",b,sep="_")
			.dbImportPlot(con,name,"saturation","biotype",
				jsonList[["biotype"]][[b]])
		}
	}
	if ("readnoise" %in% qcPlots) {
		disp("  Importing readnoise...")
		json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
			whichPlot="readnoise",output="json")
		.dbImportPlot(con,"ReadNoise","readnoise","generic",json)
	}
	if ("pairwise" %in% qcPlots) {
		disp("  Importing pairwise...")
		jsonList <- diagplotPairs(geneCounts,output="json")
		for (name in names(jsonList$xy)) {
			disp("    ",name)
			.dbImportPlot(con,name,"pairwise","xy",jsonList$xy[[name]])
			.dbImportPlot(con,name,"pairwise","md",jsonList$md[[name]])
		}
	}
	if ("filtered" %in% qcPlots) {
		disp("  Importing filtered...")
		jsonList <- diagplotFiltered(geneDataFiltered,totalGeneData,
			output="json")
		.dbImportPlot(con,"filtered_chromosome","filtered","chromosome",
			jsonList[["chromosome"]])
		.dbImportPlot(con,"filtered_biotype","filtered","biotype",
			jsonList[["biotype"]])
	}

	if ("boxplot" %in% qcPlots) {
		disp("  Importing boxplot...")
		jsonUnorm <- diagplotBoxplot(geneCounts,name=sampleList,
			isNorm=FALSE,output="json")
		jsonNorm <- diagplotBoxplot(normGenes,name=sampleList,
			isNorm=FALSE,output="json")
		.dbImportPlot(con,"Boxplot","boxplot","unorm",jsonUnorm)
		.dbImportPlot(con,"Boxplot","boxplot","norm",jsonNorm)
	}
	if ("gcbias" %in% qcPlots) {
		disp("  Importing gcbias...")
		covar <- as.numeric(geneData$gc_content)
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,covar=covar,
			isNorm=FALSE,whichPlot="gcbias",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,covar=covar,
			isNorm=TRUE,whichPlot="gcbias",output="json")
		.dbImportPlot(con,"GCBias","gcbias","unorm",jsonUnorm)
		.dbImportPlot(con,"GCBias","gcbias","norm",jsonNorm)
	}

	if ("lengthbias" %in% qcPlots) {
		disp("  Importing lengthbias...")
		covar <- width(geneData)
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,covar=covar,
			isNorm=FALSE,whichPlot="lengthbias",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,covar=covar,
			isNorm=TRUE,whichPlot="lengthbias",output="json")
		.dbImportPlot(con,"LengthBias","lengthbias","unorm",jsonUnorm)
		.dbImportPlot(con,"LengthBias","lengthbias","norm",jsonNorm)
	}

	if ("meandiff" %in% qcPlots) {
		disp("  Importing meandif...")
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,isNorm=FALSE,
			whichPlot="meandiff",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,isNorm=TRUE,
			whichPlot="meandiff",output="json")
		for (n in names(jsonUnorm)) {
			for (s in names(jsonUnorm[[n]])) {
				disp("    ",s)
				.dbImportPlot(con,s,"meandiff","unorm",
					jsonUnorm[[n]][[s]])
				.dbImportPlot(con,s,"meandiff","norm",
					jsonNorm[[n]][[s]])
			}
		}
	}
	if ("meanvar" %in% qcPlots) {
		disp("  Importing meanvar...")
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,isNorm=FALSE,
			whichPlot="meanvar",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,isNorm=TRUE,
			whichPlot="meanvar",output="json")
		.dbImportPlot(con,"MeanVar","meanvar","unorm",jsonUnorm)
		.dbImportPlot(con,"MeanVar","meanvar","orm",jsonNorm)
	}
	if ("rnacomp" %in% qcPlots) {
		disp("  Importing rnacomp...")
		cap1 <- capture.output({
			jsonUnorm <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
				whichPlot="rnacomp",isNorm=FALSE,output="json")
		})
		cap2 <- capture.output({
			jsonNorm <- diagplotNoiseq(normGenes,sampleList,covars=covarsRaw,
				whichPlot="rnacomp",isNorm=TRUE,output="json")
		})
		.dbImportPlot(con,"RnaComp","rnacomp","unorm",jsonUnorm)
		.dbImportPlot(con,"RnaComp","rnacomp","norm",jsonNorm)
	}
	if ("volcano" %in% qcPlots) {
		disp("  Importing volcano")
		nn <- names(contrastList)
		for (n in nn) {
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			for (contr in colnames(fc)) {
				disp("    ",n," ",contr)
				json <- diagplotVolcano(fc[,contr],sumpList[[n]],contr,
					altNames=geneDataExpr$gene_name,output="json")
				.dbImportPlot(con,paste("volcano",contr,sep="_"),
					"volcano","generic",json)
			}
		}
	}
	if ("mastat" %in% qcPlots) {
		disp("  Importing mastat")
		nn <- names(contrastList)
		for (n in nn) {
			m <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			a <- makeA(n,sampleList,normGenesExpr,1)
			for (contr in colnames(m)) {
				disp("    ",n," ",contr)
				json <- diagplotMa(m[,contr],a[,contr],sumpList[[n]],
					contr,altNames=geneDataExpr$gene_name,output="json")
				.dbImportPlot(con,paste("mastat",contr,sep="_"),
					"mastat","generic",json)
			}
		}
	}
	if ("biodist" %in% qcPlots) {
		disp("  Importing biodist")
		nn <- names(contrastList)
		for (n in nn) {
			disp("    ",n)
			cap <- capture.output({
				json <- diagplotNoiseq(normGenesExpr,sampleList,
					covars=covarsStat,whichPlot="biodist",
					biodistOpts=list(p=cpList[[cnt]],pcut=pcut,name=cnt),
					output="json")
			})
			.dbImportPlot(con,paste("biodist",n,sep="_"),"biodist",
				"chromosome",json$chromosome)
			.dbImportPlot(con,paste("biodist",n,sep="_"),"biodist",
				"biotype",json$biotype)
		}
	}
	if ("statvenn" %in% qcPlots) {
		disp("  Importing statvenn")
		nn <- names(contrastList)
		geneNames <- as.character(geneDataExpr$gene_name)
		names(geneNames) <- names(geneDataExpr)
		for (n in nn) {
			disp("    ",n)
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			# To narrow We take always the first condition versus the reference
			# to narrow down a bit the list of DE genes
			fc <- fc[,1,drop=FALSE]
			fcmat <- matrix(data=apply(fc,2,function(x) {
				rep(x,length(colnames(cpList[[n]])))
			}),ncol=ncol(fc)*ncol(cpList[[n]]))
			colnames(fcmat) <- colnames(cpList[[n]])
			# Deregulated
			json <- makeJVennStatData(cpList[[n]],fcmat=fcmat,pcut=pcut,
				direction="dereg",altNames=geneNames[rownames(cpList[[n]])])
			.dbImportPlot(con,paste("statvenn",n,"dereg",sep="_"),"statvenn",
				"stat_dereg",toJSON(json,auto_unbox=TRUE,null="null"))
			# Up
			json <- makeJVennStatData(cpList[[n]],fcmat=fcmat,pcut=pcut,
				direction="up",altNames=geneNames[rownames(cpList[[n]])])
			.dbImportPlot(con,paste("statvenn",n,"up",sep="_"),"statvenn",
				"stat_up",toJSON(json,auto_unbox=TRUE,null="null"))
			# Down
			json <- makeJVennStatData(cpList[[n]],fcmat=fcmat,pcut=pcut,
				direction="down",altNames=geneNames[rownames(cpList[[n]])])
			.dbImportPlot(con,paste("statvenn",n,"down",sep="_"),"statvenn",
				"stat_down",toJSON(json,auto_unbox=TRUE,null="null"))
		}
	}
	
	if ("foldvenn" %in% qcPlots) {
		disp("  Importing foldvenn")
		nn <- names(contrastList)
		
		# We need to construct a matrix with p-values for each contrast and
		# each statistical algorithm
		if (ncol(cpList[[1]]) > 1) {
			N <- names(cpList)
			prePmat <- lapply(N,function(n,C,S) {
				many <- C[[n]]
				one <- S[[n]]
				mat <- cbind(many,one)
				colnames(mat)[ncol(mat)] <- metaP
				return(mat)
			},cpList,sumpList)
			names(prePmat) <- N
		}
		else
			prePmat <- sumpList
		
		# Now we need to reformat prePmat to a list names according to each
		# algorithm
		algs <- colnames(prePmat[[1]])
		pmat <- lapply(algs,function(a,P) {
			algMat <- matrix(NA,nrow(P[[1]]),length(P))
			colnames(algMat) <- names(P)
			rownames(algMat) <- rownames(P[[1]])
			for (cnt in names(P))
				algMat[,cnt] <- prePmat[[cnt]][,a]
			return(algMat)
		},prePmat)
		
		# We need to construct a matrix with fold changes for each contrast
		fList <- vector("list",length(nn))
		names(fList) <- nn
		for (n in nn) {
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			fc <- fc[,1,drop=FALSE]
			fList[[n]] <- fc
		}
		fcmat <- do.call("cbind",fList)
		fcmat <- fcmat[rownames(pmat),]
		
		# Based on pmat and fcmat, we calculate Venns
		geneNames <- geneDataExpr$gene_name
		names(geneNames) <- names(geneDataExpr)
		for (n in names(pmat)) {
			disp("    ",n)
			jsonList[[n]] <- list(dereg=NULL,up=NULL,down=NULL)
			for (d in c("dereg","up","down"))
				json <- makeJVennFoldData(pmat=pmat,fcmat=fcmat,pcut=pcut,
					direction=d,altNames=geneNames[rownames(pmat[[n]])])
				.dbImportPlot(con,paste("foldvenn",n,d,sep="_"),"foldvenn",
					paste0("fold_",d),toJSON(json,auto_unbox=TRUE,null="null"))
		}
	}
	
	if ("deregulogram" %in% qcPlots) {
		disp("  Importing deregulogram")
		cntPairs <- combn(names(contrastList),2)
		for (i in 1:ncol(cntPairs)) {
			disp("    ",cntPairs[1,i]," and ",cntPairs[2,i])
			namc <- paste0(cntPairs[1,i],"__",cntPairs[2,i])
			fmat <- cbind(
				log2(makeFoldChange(cntPairs[1,i],sampleList,
					normGenesExpr,1))[,1,drop=FALSE],
				log2(makeFoldChange(cntPairs[2,i],sampleList,
					normGenesExpr,1))[,1,drop=FALSE]
			)
			pmat <- do.call("cbind",sumpList[c(cntPairs[1,i],
				cntPairs[2,i])])
			colnames(pmat) <- colnames(fmat)
			json <- diagplotDeregulogram(fmat,pmat,pcut=pcut,fcut=1,
				output="json")
			.dbImportPlot(con,paste("deregulogram",namc,sep="_"),"deregulogram",
				"generic",json)
		}
	}

	# Close SQLite connection
	dbDisconnect(con)
}

.createDexiePlotDb <- function(obj) {
	sampleList <- obj$sampleList
	contrastList <- obj$contrastList
	qcPlots <- obj$qcPlots
	geneCounts <- obj$geneCounts
	geneData <- obj$geneData
	geneDataExpr <- obj$geneDataExpr
	geneDataFiltered <- obj$geneDataFiltered
	totalGeneData <- obj$totalGeneData
	normGenes <- obj$normGenes
	normGenesExpr <- obj$normGenesExpr
	cpList <- obj$cpList
	sumpList <- obj$sumpList
	pcut <- obj$pcut
	metaP <- obj$metaP
	PROJECT_PATH <- obj$PROJECT_PATH
	
	samples <- unlist(sampleList)
	nsa <- length(samples)
	
	covarsRaw <- covarsStat <- NULL
	if (any(qcPlots %in% c("biodetection","countsbio","saturation",
		"rnacomp","readnoise"))) {
		covarsRaw <- list(
			data=geneCounts,
			length=width(geneData),
			gc=as.numeric(geneData$gc_content),
			chromosome=data.frame(
				chromosome=as.character(seqnames(geneData)),
				start=start(geneData),
				end=end(geneData)
			),
			factors=data.frame(class=asClassVector(sampleList)),
			biotype=as.character(geneData$biotype),
			gene_name=as.character(geneData$gene_name)
		)
	}

	if ("biodist" %in% qcPlots) {
		covarsStat <- list(
			data=normGenesExpr,
			length=width(geneDataExpr),
			gc=as.numeric(geneDataExpr$gc_content),
			chromosome=data.frame(
				chromosome=as.character(seqnames(geneDataExpr)),
				start=start(geneDataExpr),
				end=end(geneDataExpr)
			),
			factors=data.frame(class=asClassVector(sampleList)),
			biotype=as.character(geneDataExpr$biotype),
			gene_name=as.character(geneDataExpr$gene_name)
		)
	}

	plots <- list()
	listIndex <- 0
	
	if ("mds" %in% qcPlots) {
		disp("  Importing mds...")
		json <- diagplotMds(geneCounts,sampleList,output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="MDS",
			type="mds",
			subtype="generic",
			json=fromJSON(json)
		)
	}
	
	if ("biodetection" %in% qcPlots) {
		disp("  Importing biodetection...")
		cap <- capture.output({
			json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
				whichPlot="biodetection",output="json")
		})
		for (i in 1:length(json)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("biodetection",samples[i],sep="_"),
				type="biodetection",
				subtype="generic",
				json=fromJSON(json[[samples[i]]])
			)
		}
	}
	
	if ("countsbio" %in% qcPlots) {
		disp("  Importing countsbio...")
		cap <- capture.output({
			jsonList <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
				whichPlot="countsbio",output="json")
		})
		for (i in 1:length(jsonList$sample)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("countsbio",samples[i],sep="_"),
				type="countsbio",
				subtype="sample",
				json=fromJSON(jsonList$sample[[samples[i]]])
			)
		}	
		biotypes <- names(jsonList$biotype)
		for (i in 1:length(jsonList$biotype)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("countsbio",biotypes[i],sep="_"),
				type="countsbio",
				subtype="biotype",
				json=fromJSON(jsonList$biotype[[biotypes[i]]])
			)
		}
	}
	
	if ("saturation" %in% qcPlots) {
		disp("  Importing saturation...")
		jsonList <- diagplotNoiseq(geneCounts,sampleList,
			covars=covarsRaw,whichPlot="saturation",output="json")
		for (i in 1:length(jsonList$sample)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("saturation",samples[i],sep="_"),
				type="saturation",
				subtype="sample",
				json=fromJSON(jsonList$sample[[samples[i]]])
			)
		}
		biotypes <- names(jsonList$biotype)
		for (i in 1:length(jsonList$biotype)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("saturation",biotypes[i],sep="_"),
				type="saturation",
				subtype="biotype",
				json=fromJSON(jsonList$biotype[[biotypes[i]]])
			)
		}
	}
	
	if ("readnoise" %in% qcPlots) {
		disp("  Importing readnoise...")
		json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
			whichPlot="readnoise",output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="ReadNoise",
			type="readnoise",
			subtype="generic",
			json=fromJSON(json)
		)
	}
	
	if ("pairwise" %in% qcPlots) {
		disp("  Importing pairwise...")
		jsonList <- diagplotPairs(geneCounts,output="json")
		nams <- names(jsonList$xy)
		for (i in 1:length(jsonList$xy)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=nams[i],
				type="pairwise",
				subtype="xy",
				json=fromJSON(jsonList$xy[[nams[i]]])
			)
		}
		for (i in 1:length(jsonList$md)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=nams[i],
				type="pairwise",
				subtype="md",
				json=fromJSON(jsonList$md[[nams[i]]])
			)
		}
	}
	
	if ("filtered" %in% qcPlots) {
		disp("  Importing filtered...")
		jsonList <- diagplotFiltered(geneDataFiltered,totalGeneData,
			output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="filtered_chromosome",
			type="filtered",
			subtype="chromosome",
			json=fromJSON(jsonList[["chromosome"]])
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="filtered_biotype",
			type="filtered",
			subtype="biotype",
			json=fromJSON(jsonList[["biotype"]])
		)
	}

	if ("boxplot" %in% qcPlots) {
		disp("  Importing boxplot...")
		jsonUnorm <- diagplotBoxplot(geneCounts,name=sampleList,
			isNorm=FALSE,output="json")
		jsonNorm <- diagplotBoxplot(normGenes,name=sampleList,
			isNorm=FALSE,output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="Boxplot",
			type="boxplot",
			subtype="unorm",
			json=fromJSON(jsonUnorm)
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="Boxplot",
			type="boxplot",
			subtype="norm",
			json=fromJSON(jsonNorm)
		)
	}
	
	if ("gcbias" %in% qcPlots) {
		disp("  Importing gcbias...")
		covar <- as.numeric(geneData$gc_content)
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,covar=covar,
			isNorm=FALSE,whichPlot="gcbias",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,covar=covar,
			isNorm=TRUE,whichPlot="gcbias",output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="GCBias",
			type="gcbias",
			subtype="unorm",
			json=fromJSON(jsonUnorm)
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="GCBias",
			type="gcbias",
			subtype="norm",
			json=fromJSON(jsonNorm)
		)
	}

	if ("lengthbias" %in% qcPlots) {
		disp("  Importing lengthbias...")
		covar <- width(geneData)
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,covar=covar,
			isNorm=FALSE,whichPlot="lengthbias",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,covar=covar,
			isNorm=TRUE,whichPlot="lengthbias",output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="LengthBias",
			type="lengthbias",
			subtype="unorm",
			json=fromJSON(jsonUnorm)
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="LengthBias",
			type="lengthbias",
			subtype="norm",
			json=fromJSON(jsonNorm)
		)
	}

	if ("meandiff" %in% qcPlots) {
		disp("  Importing meandif...")
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,isNorm=FALSE,
			whichPlot="meandiff",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,isNorm=TRUE,
			whichPlot="meandiff",output="json")
		np <- sum(lengths(jsonUnorm))
		for (i in 1:length(jsonUnorm)) {
			for (j in 1:length(jsonUnorm[[i]])) {
				listIndex <- listIndex + 1
				plots[[listIndex]] <- list(
					name=names(jsonUnorm[[i]])[j],
					type="meandiff",
					subtype="unorm",
					json=fromJSON(jsonUnorm[[i]][[j]])
				)
			}
		}
		for (i in 1:length(jsonNorm)) {
			for (j in 1:length(jsonNorm[[i]])) {
				listIndex <- listIndex + 1
				plots[[listIndex]] <- list(
					name=names(jsonNorm[[i]])[j],
					type="meandiff",
					subtype="norm",
					json=fromJSON(jsonNorm[[i]][[j]])
				)
			}
		}
	}
	
	if ("meanvar" %in% qcPlots) {
		disp("  Importing meanvar...")
		jsonUnorm <- diagplotEdaseq(geneCounts,sampleList,isNorm=FALSE,
			whichPlot="meanvar",output="json")
		jsonNorm <- diagplotEdaseq(normGenes,sampleList,isNorm=TRUE,
			whichPlot="meanvar",output="json")
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="MeanVar",
			type="meanvar",
			subtype="unorm",
			json=fromJSON(jsonUnorm)
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="MeanVar",
			type="meanvar",
			subtype="norm",
			json=fromJSON(jsonNorm)
		)
	}
	
	if ("rnacomp" %in% qcPlots) {
		disp("  Importing rnacomp...")
		cap1 <- capture.output({
			jsonUnorm <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
				whichPlot="rnacomp",isNorm=FALSE,output="json")
		})
		cap2 <- capture.output({
			jsonNorm <- diagplotNoiseq(normGenes,sampleList,covars=covarsRaw,
				whichPlot="rnacomp",isNorm=TRUE,output="json")
		})
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="RnaComp",
			type="rnacomp",
			subtype="unorm",
			json=fromJSON(jsonUnorm)
		)
		listIndex <- listIndex + 1
		plots[[listIndex]] <- list(
			name="RnaComp",
			type="rnacomp",
			subtype="norm",
			json=fromJSON(jsonNorm)
		)
	}
	
	if ("volcano" %in% qcPlots) {
		disp("  Importing volcano")
		nn <- names(contrastList)
		json <- list()
		namc <- character(0)
		index <- 0
		for (n in nn) {
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			for (contr in colnames(fc)) {
				disp("    ",n," ",contr)
				index <- index + 1
				namc <- c(namc,contr)
				json[[index]] <- diagplotVolcano(fc[,contr],sumpList[[n]],contr,
					altNames=geneDataExpr$gene_name,output="json")
			}
		}
		names(json) <- namc
		for (i in 1:length(json)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("volcano",namc[i],sep="_"),
				type="volcano",
				subtype="generic",
				json=fromJSON(json[[i]])
			)
		}
	}
	
	if ("mastat" %in% qcPlots) {
		disp("  Importing mastat")
		nn <- names(contrastList)
		json <- list()
		namc <- character(0)
		index <- 0
		for (n in nn) {
			m <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			a <- makeA(n,sampleList,normGenesExpr,1)
			for (contr in colnames(m)) {
				disp("    ",n," ",contr)
				index <- index + 1
				namc <- c(namc,contr)
				json[[index]] <- diagplotMa(m[,contr],a[,contr],sumpList[[n]],
					contr,altNames=geneDataExpr$gene_name,output="json")
			}
		}
		names(json) <- namc
		for (i in 1:length(json)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("mastat",namc[i],sep="_"),
				type="mastat",
				subtype="generic",
				json=fromJSON(json[[i]])
			)
		}
	}
	
	if ("biodist" %in% qcPlots) {
		disp("  Importing biodist")
		nn <- names(contrastList)
		jsonList <- vector("list",length(nn))
		names(jsonList) <- nn
		for (n in nn) {
			disp("    ",n)
			cap <- capture.output({
				jsonList[[n]] <- diagplotNoiseq(normGenesExpr,sampleList,
					covars=covarsStat,whichPlot="biodist",
					biodistOpts=list(p=cpList[[n]],pcut=pcut,name=n),
					output="json")
			})
		}
		for (i in 1:length(jsonList)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("biodist",nn[i],sep="_"),
				type="biodist",
				subtype="chromosome",
				json=fromJSON(jsonList[[i]]$chromosome)
			)
		}
		for (i in 1:length(jsonList)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("biodist",nn[i],sep="_"),
				type="biodist",
				subtype="biotype",
				json=fromJSON(jsonList[[i]]$biotype)
			)
		}
	}
	
	if ("statvenn" %in% qcPlots) {
		disp("  Importing statvenn")
		nn <- names(contrastList)
		jsonList <- vector("list",length(nn))
		names(jsonList) <- nn
		geneNames <- geneDataExpr$gene_name
		names(geneNames) <- names(geneDataExpr)
		for (n in nn) {
			disp("    ",n)
			jsonList[[n]] <- list(dereg=NULL,up=NULL,down=NULL)
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			# To narrow We take always the first condition versus the reference
			# to narrow down a bit the list of DE genes
			fc <- fc[,1,drop=FALSE]
			fcmat <- matrix(data=apply(fc,2,function(x) {
				rep(x,length(colnames(cpList[[n]])))
			}),ncol=ncol(fc)*ncol(cpList[[n]]))
			colnames(fcmat) <- colnames(cpList[[n]])
			jsonList[[n]]$dereg <- makeJVennStatData(cpList[[n]],fcmat=fcmat,
				pcut=pcut,direction="dereg",
				altNames=geneNames[rownames(cpList[[n]])])
			jsonList[[n]]$up <- makeJVennStatData(cpList[[n]],fcmat=fcmat,
				pcut=pcut,direction="up",
				altNames=geneNames[rownames(cpList[[n]])])
			jsonList[[n]]$down <- makeJVennStatData(cpList[[n]],fcmat=fcmat,
				pcut=pcut,direction="down",
				altNames=geneNames[rownames(cpList[[n]])])
		}
		for (i in 1:length(jsonList)) {
			for (d in names(jsonList[[i]])) {
				listIndex <- listIndex + 1
				plots[[listIndex]] <- list(
					name=paste("statvenn",nn[i],d,sep="_"),
					type="statvenn",
					subtype=paste0("stat_",d),
					json=jsonList[[i]][[d]]
				)
			}
		}
	}
	
	if ("foldvenn" %in% qcPlots) {
		disp("  Importing foldvenn")
		nn <- names(contrastList)
		
		# We need to construct a matrix with p-values for each contrast and
		# each statistical algorithm
		if (ncol(cpList[[1]]) > 1) {
			N <- names(cpList)
			prePmat <- lapply(N,function(n,C,S) {
				many <- C[[n]]
				one <- S[[n]]
				mat <- cbind(many,one)
				colnames(mat)[ncol(mat)] <- metaP
				return(mat)
			},cpList,sumpList)
			names(prePmat) <- N
		}
		else
			prePmat <- sumpList
		
		# Now we need to reformat prePmat to a list names according to each
		# algorithm
		algs <- colnames(prePmat[[1]])
		pmat <- lapply(algs,function(a,P) {
			algMat <- matrix(NA,nrow(P[[1]]),length(P))
			colnames(algMat) <- names(P)
			rownames(algMat) <- rownames(P[[1]])
			for (cnt in names(P))
				algMat[,cnt] <- prePmat[[cnt]][,a]
			return(algMat)
		},prePmat)
		names(pmat) <- algs
		
		# We need to construct a matrix with fold changes for each contrast
		fList <- vector("list",length(nn))
		names(fList) <- nn
		for (n in nn) {
			fc <- log2(makeFoldChange(n,sampleList,normGenesExpr,1))
			fc <- fc[,1,drop=FALSE]
			fList[[n]] <- fc
		}
		fcmat <- do.call("cbind",fList)
		fcmat <- fcmat[rownames(pmat[[1]]),]
		
		# Based on pmat and fcmat, we calculate Venns
		jsonList <- vector("list",length(pmat))
		names(jsonList) <- names(pmat)
		geneNames <- geneDataExpr$gene_name
		names(geneNames) <- names(geneDataExpr)
		for (n in names(jsonList)) {
			disp("    ",n)
			jsonList[[n]] <- list(dereg=NULL,up=NULL,down=NULL)
			for (d in c("dereg","up","down"))
				jsonList[[n]][[d]] <- makeJVennFoldData(pmat=pmat[[n]],
					fcmat=fcmat,pcut=pcut,direction=d,
					altNames=geneNames[rownames(pmat[[n]])])
		}
		for (n in names(jsonList)) {
			for (d in names(jsonList[[n]])) {
				listIndex <- listIndex + 1
				plots[[listIndex]] <- list(
					name=paste("foldvenn",n,d,sep="_"),
					type="foldvenn",
					subtype=paste0("fold_",d),
					json=jsonList[[n]][[d]]
				)
			}
		}
	}
	
	if ("deregulogram" %in% qcPlots) {
		disp("  Importing deregulogram")
		cntPairs <- combn(names(contrastList),2)
		json <- vector("list",ncol(cntPairs))
		namc <- character(ncol(cntPairs))
		counter <- 0
		for (i in 1:ncol(cntPairs)) {
			counter <- counter + 1
			disp("    ",cntPairs[1,i]," and ",cntPairs[2,i])
			namc[counter] <- paste0(cntPairs[1,i],"__",cntPairs[2,i])
			fmat <- cbind(
				log2(makeFoldChange(cntPairs[1,i],sampleList,
					normGenesExpr,1))[,1,drop=FALSE],
				log2(makeFoldChange(cntPairs[2,i],sampleList,
					normGenesExpr,1))[,1,drop=FALSE]
			)
			pmat <- do.call("cbind",sumpList[c(cntPairs[1,i],
				cntPairs[2,i])])
			colnames(pmat) <- colnames(fmat)
			json[[counter]] <- diagplotDeregulogram(fmat,pmat,fcut=1,pcut=pcut,
				output="json")
		}
		for (i in 1:length(json)) {
			listIndex <- listIndex + 1
			plots[[listIndex]] <- list(
				name=paste("deregulogram",namc[i],sep="_"),
				type="deregulogram",
				subtype="generic",
				json=fromJSON(json[[i]])
			)
		}
	}
	
	disp("Writing plot database in ",file.path(PROJECT_PATH$data,"reportdb.js"))
	jsonPlots <- toJSON(plots,auto_unbox=TRUE,null="null",pretty=TRUE)
	code <- paste0("var plotData = ",jsonPlots)
	cat(code,file=file.path(PROJECT_PATH$data,"reportdb.js"),sep="\n")
}

.initReportDbTables <- function(con) {
	queries <- .reportDbTblDef()
	rs <- dbSendQuery(con,queries[[1]])
	if (dbHasCompleted(rs))
		dbClearResult(rs)
	for (n in names(queries)) {
		rs <- dbSendStatement(con,queries[[n]])
		if (dbHasCompleted(rs))
			dbClearResult(rs)
	}
}

.reportDbTblDef <- function() {
	return(list(
		enable_fkey="PRAGMA foreign_keys=1;",
		plot=paste(
			"CREATE TABLE IF NOT EXISTS plot (",
			"_id INTEGER PRIMARY KEY AUTOINCREMENT,",
			"name TEXT,",
			"type TEXT,",
			"subtype TEXT,",
			"json TEXT",
			");"
		),
		clear="DELETE FROM plot;"
	))
}

.dbImportPlot <- function(con,name,type,subtype=NULL,json) {
	if (is.null(subtype))
		query <- paste("INSERT INTO plot (name, type, subtype, json) ",
			"VALUES (","'",name,"', ","'",type,"', NULL, :j)",sep="")
	else
		query <- paste("INSERT INTO plot (name, type, subtype, json) ",
			"VALUES (","'",name,"', ","'",type,"', ","'",
			subtype,"', :j)",sep="")
	#print(query)
	nr <- dbExecute(con,query,params=list(j=json))
	return(nr)
}

.downloadJsLibs <- function(prPath,reportDb) {
	disp("Downloading required JavaScript libraries...")
	if (!file.exists(file.path(prPath$js,"pace.min.js")))
		download.file(paste0("https://raw.github.com/HubSpot/pace/",
			"v1.0.0/pace.min.js"),
			file.path(prPath$js,"pace.min.js"))
	if (!file.exists(file.path(prPath$js,"highcharts.js")))
		download.file("https://code.highcharts.com/highcharts.js",
			file.path(prPath$js,"highcharts.js"))
	if (!file.exists(file.path(prPath$js,"highcharts-more.js")))
		download.file("https://code.highcharts.com/highcharts-more.js",
			file.path(prPath$js,"highcharts-more.js"))
	if (!file.exists(file.path(prPath$js,"exporting.js")))
		download.file("https://code.highcharts.com/modules/exporting.js",
			file.path(prPath$js,"exporting.js"))
	if (!file.exists(file.path(prPath$js,"offline-exporting.js")))
		download.file(
			"https://code.highcharts.com/modules/offline-exporting.js",
			file.path(prPath$js,"offline-exporting.js"))
	if (!file.exists(file.path(prPath$js,"export-data.js")))
		download.file(
			"https://code.highcharts.com/modules/export-data.js",
			file.path(prPath$js,"export-data.js"))
	if (!file.exists(file.path(prPath$js,"canvas2svg.js")))
		download.file(
			"http://jvenn.toulouse.inra.fr/app/js/canvas2svg.js",
			file.path(prPath$js,"canvas2svg.js"))
	if (!file.exists(file.path(prPath$js,"jvenn.min.js")))
		download.file(
			"http://jvenn.toulouse.inra.fr/app/js/jvenn.min.js",
			file.path(prPath$js,"jvenn.min.js"))
	if (reportDb == "sqlite") {
		if (!file.exists(file.path(prPath$js,"sql.js")))
			download.file(
		"https://cdnjs.cloudflare.com/ajax/libs/sql.js/0.5.0/js/sql.js",
		file.path(prPath$js,"sql.js"))
	}
	else if (reportDb == "dexie") {
		if (!file.exists(file.path(prPath$js,"dexie.min.js")))
			download.file(
				"https://unpkg.com/dexie@2.0.4/dist/dexie.min.js",
				file.path(prPath$js,"dexie.min.js"))
	}
}

.formatForReport <- function(x,o=NULL,r=NULL,v=NULL) {
	# Replace seqnames name
	colnames(x)[1] <- "chromosome"
	
	# Link to UCSC genome browser with position
	if (!is.null(o)) {
		x$chromosome <- 
			paste('<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=',
				getUcscOrganism(o),'&position=',paste(x$chromosome,":",
				x$start,"-",x$end,sep=""),'" target="_blank">',x$chromosome,
				'</a>',sep="")
	}
	
	# Link to gene pages - only for Ensembl, RefSeq, not custom GTFs! UCSC does
	# not offer a related page.
	if (!is.null(r)) {
		if (r == "ensembl") {
			h <- getHost(o,v)
			if ("gene_id" %in% names(x)) # transLevel = "gene"
				x$gene_id <- paste('<a href="',h,'/Gene/Summary?g=',x$gene_id,
					'" target="_blank">',x$gene_id,'</a>',sep="")
			else if ("transcript_id" %in% names(x)) {# transLevel = "transcript"
				# Stub
			}
			else if ("exon_id" %in% names(x)) { # transLevel = "exon"
				# Stub
			}
		}
		else if (r == "refseq") {
			if ("gene_id" %in% names(x)) # transLevel = "gene"
				x$gene_id <- paste('<a href="https://www.ncbi.nlm.nih.gov/',
					'nuccore/',x$gene_id,'" target="_blank">',x$gene_id,'</a>',
					sep="")
			else if ("transcript_id" %in% names(x)) { # transLevel = "transcript"
				# Stub
			}
			else if ("exon_id" %in% names(x)) { # transLevel = "exon"
				# Stub
			}
		}
	}
	
	# Round rpgm and fold change
	tmpi1 <- grep("log2_normalized_fold_change_",colnames(x))
	tmpi2 <- grep("rpgm_normalized_mean_counts_",colnames(x))
	tmpi <- c(tmpi1,tmpi2)
	if (length(tmpi) > 0)
		x[,tmpi] <- round(x[,tmpi],digits=4)
	
	# Round p-values and FDR
	tmpi1 <- grep("p-value",colnames(x))
	tmpi2 <- grep("FDR",colnames(x))
	tmpi <- c(tmpi1,tmpi2)
	if (length(tmpi) > 0)
		x[,tmpi] <- round(x[,tmpi],digits=6)
	
	colnames(x) <- gsub("log2_normalized_fold_change_","",colnames(x))
	colnames(x) <- gsub("rpgm_normalized_mean_counts_","",colnames(x))
	tmpi <- grep("p-value",colnames(x))
	if (length(tmpi) > 0)
		colnames(x)[tmpi] <- "p-value"
	tmpi <- grep("FDR",colnames(x))
	if (length(tmpi) > 0)
		colnames(x)[tmpi] <- "FDR"
	
	return(x)
}
