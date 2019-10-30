metaseqr2 <- function(
    counts,
    sampleList,
    excludeList=NULL,
    fileType=c("auto","sam","bam","bed"),
    path=NULL,
    contrast=NULL,
    libsizeList=NULL,
    embedCols=list(
		idCol=4,
		gcCol=NA,
		nameCol=NA,
		btCol=NA
    ),
    annotation=NULL,
    org=c("hg18","hg19","hg38","mm9","mm10","rn5","rn6","dm3","dm6",
        "danrer7","pantro4","susscr3","tair10","equcab2"),
    refdb=c("ensembl","ucsc","refseq"),
    version="auto",
    transLevel=c("gene","transcript","exon"),
    countType=c("gene","exon","utr"),
    utrOpts=list(
		frac=1,
		minLength=300,
		downstream=50
    ),
    exonFilters=list(
        minActiveExons=list(
            exonsPerGene=5,
            minExons=2,
            frac=1/5
        )
    ),
    geneFilters=list(
        length=list(
            length=500
        ),
        avgReads=list(
            averagePerBp=100,
            quantile=0.25
        ),
        expression=list(
            median=TRUE,
            mean=FALSE,
            quantile=NA,
            known=NA,
            custom=NA
        ),
        biotype=getDefaults("biotypeFilter",org[1]),
        presence=list(
            frac=0.25,
            minCount=10,
            perCondition=FALSE
        )
    ),
    whenApplyFilter=c("postnorm","prenorm"),
    normalization=c("deseq","deseq2","edaseq","edger","noiseq","nbpseq",
        "absseq","dss","each","none"),
    normArgs=NULL,
    statistics=c("deseq","deseq2","edger","noiseq","bayseq","limma","nbpseq",
        "absseq","dss"),
    statArgs=NULL,
    adjustMethod=sort(c(p.adjust.methods,"qvalue")),
    metaP=if (length(statistics)>1) c("simes","bonferroni","fisher",
        "dperm_min","dperm_max","dperm_weight","fperm","whitlock","minp","maxp",
        "weight","pandora","none") else "none",
    weight=rep(1/length(statistics),length(statistics)),
    nperm=10000,
    reprod=TRUE,
    pcut=NA,
    logOffset=1,
    preset=NULL, # An analysis strictness preset
    qcPlots=c(
        "mds","biodetection","countsbio","saturation","readnoise","filtered",
        "correl","pairwise","boxplot","gcbias","lengthbias","meandiff",
        "meanvar","rnacomp","deheatmap","volcano","biodist","mastat","venn"
    ),
    figFormat=c("png","jpg","tiff","bmp","pdf","ps"),
    outList=FALSE,
    exportWhere=NA, # An output directory for the project
    exportWhat=c("annotation","p_value","adj_p_value","meta_p_value",
        "adj_meta_p_value","fold_change","stats","counts","flags"),
    exportScale=c("natural","log2","log10","vst","rpgm"),
    exportValues=c("raw","normalized"),
    exportStats=c("mean","median","sd","mad","cv","rcv"),
    exportCountsTable=FALSE,
    restrictCores=0.6,
    report=TRUE,
    reportTop=0.1,
    reportTemplate="default",
    saveGeneModel=TRUE,
    verbose=TRUE,
    runLog=TRUE,
    reportDb=c("sqlite","dexie"),
    localDb=file.path(system.file(package="metaseqR2"),"annotation.sqlite"),
    progressFun=NULL,
    offlineReport=TRUE,
    ...
) {
	# Save function call for report
	FUN_CALL <- deparse(sys.call())
	
    # Check for older argument names and adjust
    args <- as.list(match.call())
    backCheckedArgs <- .backwardsConvertArgs(args)
    if (backCheckedArgs$mixDetected)
        stop("You must provide either old metaseqR arguments (e.g sample.list)",
            " or new metaseqR2 arguments (e.g. sampleList)! Not both!")
    if (backCheckedArgs$backDetected) {
        argMap <- .backwardsMapOld2New()
        # Check what happens with the new argument embedCols
        if (!is.null(backCheckedArgs$args$embedCols)) {
			embedCols <- backCheckedArgs$args$embedCols
			backCheckedArgs$args$embedCols <- NULL
		}
        for (n in names(backCheckedArgs$args))
            assign(argMap[[n]],eval(backCheckedArgs$args[[n]]))
        # There are some nested arguments that also need correction
        if (!is.null(backCheckedArgs$args$gene.filters)) {
            if (!is.null(geneFilters$avg.reads)) {
                geneFilters$avgReads <- geneFilters$avg.reads
                if (!is.null(geneFilters$avg.reads$average.per.bp)) {
                    geneFilters$avgReads$averagePerBp <-
                        geneFilters$avg.reads$average.per.bp
                    geneFilters$avg.reads$average.per.bp <- NULL
                }
                geneFilters$avg.reads <- NULL
            }
            if (!is.null(geneFilters$presence)) {
                if (!is.null(geneFilters$presence$min.count)) {
                    geneFilters$presence$minCount <-
                        geneFilters$presence$min.count
                    geneFilters$presence$min.count <- NULL
                }
                if (!is.null(geneFilters$presence$per.condition)) {
                    geneFilters$presence$perCondition <-
                        geneFilters$presence$per.condition
                    geneFilters$presence$per.condition <- NULL
                }
            }
        }
        if (!is.null(backCheckedArgs$args$exon.filters)) {
            if (!is.null(exonFilters$min.active.exons)) {
                exonFilters$minActiveExons <- exonFilters$min.active.exons
                if (!is.null(exonFilters$min.active.exons$exons.per.gene)) {
                    exonFilters$minActiveExons$exonsPerGene <-
                        exonFilters$min.active.exons$exons.per.gene
                    exonFilters$min.active.exons$exons.per.gene <- NULL
                }
                if (!is.null(exonFilters$min.active.exons$min.exons)) {
                    exonFilters$minActiveExons$minExons <-
                        exonFilters$min.active.exons$min.exons
                    exonFilters$min.active.exons$min.exons <- NULL
                }
                exonFilters$min.active.exons <- NULL
            }
        }
    }
    else
        # Check if there are any mispelled or invalid parameters and throw a
        # warning
        checkMainArgs(args)    
    
    # Check essential arguments
    fromRaw <- fromPrevious <- FALSE
    if (missing(counts) && (missing(sampleList) || is.list(sampleList)))
        stop("You must provide a file with genomic region (gene, exon, etc.) ",
            "counts or an input targets file to create input from! If the ",
            "counts file is missing, sampleList cannot be missing or it must ",
            "be a targets file with at least three columns! See the ",
            "readTargets function. counts may also be a gene model list ",
            "(see the documentation)")
    if (!missing(counts) && !missing(sampleList) && is.character(counts) 
        && file.exists(counts) && length(grep(".RData$",counts))>0) {
        warning("When restoring a previous analysis, sampleList argument is ",
            "not necessary! Ignoring...",immediate.=TRUE)
        fromPrevious <- TRUE
        disp("Restoring previous analysis from ",basename(counts))
        tmpEnv <- .backwardsCompatibility(counts)
        sampleList <- tmpEnv$sampleList
        countType <- tmpEnv$countType
    }
    if (!missing(counts) && missing(sampleList) && is.character(counts) 
        && file.exists(counts) && length(grep(".RData$",counts))>0) {
        # Time to load previous analysis if existing
        fromPrevious <- TRUE
        message("Restoring previous analysis from ",basename(counts))
        tmpEnv <- .backwardsCompatibility(counts)
        sampleList <- tmpEnv$sampleList
        countType <- tmpEnv$countType
    }
    if (missing(sampleList) && !fromPrevious || (!is.list(sampleList) &&
        !file.exists(sampleList)))
        stop("You must provide a list with condition names and sample names ",
            "(same as in the counts file) or an input file to create the ",
            "sample list from!")
    if (!missing(sampleList) && !is.list(sampleList) 
        && file.exists(sampleList) && !missing(counts) && !fromPrevious)
        sampleList <- makeSampleList(sampleList)
    if (!missing(sampleList) && !is.list(sampleList) 
        && file.exists(sampleList) && missing(counts)) {
        counts <- NULL
        theList <- readTargets(sampleList,path=path)
        sampleList <- theList$samples
        fileList <- theList$files
        if (tolower(fileType[1])=="auto")
            fileType <- theList$type
        if (is.null(fileType))
            stop(paste("The type of the input files could not be recognized!",
                "Please specify (BAM or BED)..."))
        fromRaw <- TRUE
    }
    
    ## If report requested, RSQLite must be present
    #if (report && !requireNamespace(RSQLite))
	#	stop("R package RSQLite is required to build metaseqR reports!")
	if (report && !requireNamespace("pander"))
		stop("R package pander is required to build metaseqR2 reports!")
    
    # Initialize environmental variables
    HOME <- system.file(package="metaseqR2")
    TEMPLATE <- HOME
    # Globalize the project's verbosity and logger
    if (fromRaw)
        PROJECT_PATH <- makeProjectPath(exportWhere)
    else
        PROJECT_PATH <- makeProjectPath(exportWhere,counts)
    assign("VERBOSE",verbose,envir=metaEnv)
    if (runLog)
        logger <- create.logger(logfile=file.path(PROJECT_PATH$logs,
            "metaseqr_run.log"),level=2,logformat="%d %c %m")
    else
        logger <- NULL
    assign("LOGGER",logger,envir=metaEnv)
    
    # Check if sample names match in file/df and list, otherwise meaningless to 
    # proceed
    if (!fromRaw && !fromPrevious) {
        if (!is.data.frame(counts) && !is.list(counts)) {
            if (file.exists(counts)) {
                aline <- read.delim(counts,nrows=5) # Read the 1st lines
                aline <- colnames(aline)
            }
            else
                stopwrap("The counts file you provided does not exist!")
        }
        else if (is.data.frame(counts))
            aline <- colnames(counts)
        else if (is.list(counts))
            aline <- names(counts)
        samples <- unlist(sampleList,use.names=FALSE)
        if (length(which(!is.na(match(samples,aline)))) != length(samples))
            stopwrap("The sample names provided in the counts file/list do ",
                "not match with those of the sampleList!")
    }
    
    # If exclude list given, check that it's a subset of sampleList, otherwise
    # just ignore excludeList
    if (!is.null(excludeList) && !is.na(excludeList)) {
        sl <- unlist(sampleList)
        el <- unlist(excludeList)
        if (length(intersect(sl,el)) != length(el)) {
            warnwrap("Some samples in excludeList do not match those in the ",
                "initial sampleList! Ignoring...",now=TRUE)
            excludeList <- NULL
        }
    }

    fileType <- tolower(fileType[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    version <- version[1]
    transLevel <- tolower(transLevel[1])
    countType <- tolower(countType[1])
    whenApplyFilter <- tolower(whenApplyFilter[1])
    normalization <- tolower(normalization[1])
    adjustMethod <- adjustMethod[1]
    metaP <- tolower(metaP[1])
    statistics <- tolower(statistics)
    figFormat <- tolower(figFormat)
    if (!is.null(qcPlots)) 
        qcPlots <- tolower(qcPlots)
    exportWhat <- tolower(exportWhat)
    exportScale <- tolower(exportScale)
    exportValues <- tolower(exportValues)
    exportStats <- tolower(exportStats)
    reportDb <- tolower(reportDb[1])
    if (!is.null(preset)) 
        preset <- tolower(preset[1])
    
    if (fromRaw)
        countsName <- "imported sam/bam/bed files"
    else {
        if (!is.data.frame(counts) && !is.null(counts) && !is.list(counts)) {
            checkFileArgs("counts",counts)
            if (fromPrevious)
                countsName <- "previously stored project"
            else
                countsName <- basename(counts)
        }
        else if (is.list(counts) && !is.data.frame(counts)) {
            countsName <- "previously stored gene model"
        }
        else {
            countsName <- "imported custom data frame"
        }
    }
	
	# annotation cannot be "embedded" when countType is other than "gene" or
	# transLevel is other than "gene"
	if (!is.null(annotation) && !is.list(annotation)
		&& is.character(annotation) && annotation == "embedded" 
		&& countType != "gene")
		stopwrap("The annotation argument cannot be \"embedded\" when the ",
			"countType argument is other than \"gene\"!")

	# annotation must be a list to be fed to buildCustomAnnotation
	if (is.list(annotation)) {
		# members are checked by buildCustomAnnotation if required
		# We only need to check that the gtfFile exists here
		if (!("gtf" %in% names(annotation)))
			stopwrap("A gtf field must be provided with an existing GTF file ",
				"when providing a list with custom annotation elements!")
		if ("gtf" %in% names(annotation) && is.character(annotation$gtf)
			&& !file.exists(annotation$gtf))
			stopwrap("An existing GTF file must be provided when providing a ",
				"list with custom annotation elements!")
	}
	
    if (is.list(counts) && !is.data.frame(counts)
		&& (countType=="exon" || countType=="utr")
        && annotation=="embedded") {
        warnwrap("annotation cannot be \"embedded\" when importing a stored ",
            "gene model! Setting to NULL...")
        #annotation <- "download"
        annotation <- NULL
    }

    if (metaP %in% c("weight","pandora","dperm_weight") && 
        abs(1-sum(weight))>1e-5)
        stopwrap("The weights given for p-value combination should sum to 1!")

    checkTextArgs("fileType",fileType,c("auto","sam","bam","bed"),
        multiarg=FALSE)
    if (!is.null(annotation) && !is.list(annotation) 
		&& !file.exists(annotation))
		checkTextArgs("annotation",annotation,"embedded",multiarg=FALSE)
	if (is.character(localDb) && file.exists(localDb)) {
		if (!.userOrg(org,localDb) && is.null(annotation))
			checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5",
				"rn6","dm3","dm6","danrer7","pantro4","susscr3","tair10",
				"equcab2"),multiarg=FALSE)
	}
	else {
		# So only some annotations can be fetched on-the-fly
		if (is.null(annotation)) 
			checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5",
				"rn6","dm3","dm6","danrer7","pantro4","susscr3","tair10",
				"equcab2"),multiarg=FALSE)
	}
	if (is.character(localDb) && file.exists(localDb)) {
		if (!.userRefdb(refdb,localDb) && is.null(annotation))
			checkTextArgs("refdb",refdb,c("ensembl","ucsc","refseq"),
				multiarg=FALSE)
	}
	else {
		# So only some annotations can be fetched on-the-fly
		if (is.null(annotation)) 
			checkTextArgs("refdb",refdb,c("ensembl","ucsc","refseq"),
				multiarg=FALSE)
	}
    checkTextArgs("transLevel",transLevel,c("gene","transcript","exon"),
        multiarg=FALSE)
    checkTextArgs("countType",countType,c("gene","exon","utr"),
        multiarg=FALSE)
    checkTextArgs("whenApplyFilter",whenApplyFilter,c("postnorm",
        "prenorm"),multiarg=FALSE)
    checkTextArgs("normalization",normalization,c("edaseq","deseq","deseq2",
        "edger","noiseq","nbpseq","absseq","dss","each","none"),multiarg=FALSE)
    checkTextArgs("statistics",statistics,c("deseq","deseq2","edger","noiseq",
        "bayseq","limma","nbpseq","absseq","dss"),multiarg=TRUE)
    checkTextArgs("metaP",metaP,c("simes","bonferroni","fisher","dperm_min",
        "dperm_max","dperm_weight","fperm","whitlock","minp","maxp","weight",
        "pandora","none"),multiarg=FALSE)
    checkTextArgs("figFormat",figFormat,c("png","jpg","tiff","bmp","pdf",
        "ps"),multiarg=TRUE)
    checkTextArgs("exportWhat",exportWhat,c("annotation","p_value",
        "adj_p_value","meta_p_value","adj_meta_p_value","fold_change","stats",
        "counts","flags"),multiarg=TRUE)
    checkTextArgs("exportScale",exportScale,c("natural","log2","log10",
        "rpgm","vst"),multiarg=TRUE)
    checkTextArgs("exportValues",exportValues,c("raw","normalized"),
        multiarg=TRUE)
    checkTextArgs("exportStats",exportStats,c("mean","median","sd","mad",
        "cv","rcv"),multiarg=TRUE)
    checkTextArgs("reportDb",reportDb,c("sqlite","dexie"))
    if (!is.null(preset))
        checkTextArgs("preset",preset,c("all_basic","all_normal","all_full",
            "medium_basic","medium_normal","medium_full","strict_basic",
            "strict_normal","strict_full"),multiarg=FALSE)
    if (!is.null(qcPlots))
        checkTextArgs("qcPlots",qcPlots,c("mds","biodetection","countsbio",
            "saturation","readnoise","correl","pairwise","boxplot","gcbias",
            "lengthbias","meandiff","meanvar","rnacomp","deheatmap","volcano",
            "biodist","filtered","mastat","venn"),multiarg=TRUE)
    if (!is.na(restrictCores)) checkNumArgs("restrictCores",restrictCores,
        "numeric",c(0,1),"botheq")
    if (!is.na(pcut)) 
		checkNumArgs("pcut",pcut,"numeric",c(0,1),"botheq")
    if (!is.null(embedCols$gcCol) && !is.na(embedCols$gcCol)) 
		checkNumArgs("embedCols$gcCol",embedCols$gcCol,"numeric",0,"gt")
    if (!is.null(embedCols$nameCol) && !is.na(embedCols$nameCol)) 
		checkNumArgs("embedCols$nameCol",embedCols$nameCol,"numeric",0,"gt")
    if (!is.null(embedCols$btCol) && !is.na(embedCols$btCol)) 
		checkNumArgs("embedCols$btCol",embedCols$btCol,"numeric",0,"gt")
    if (!is.na(logOffset)) checkNumArgs("logOffset",logOffset,"numeric",0,"gt")
    checkNumArgs("nperm",nperm,"numeric",10,"gt")
    if (!is.null(reportTop))
        checkNumArgs("reportTop",reportTop,"numeric",c(0,1),"both")
    if (!is.null(contrast)) {
        checkContrastFormat(contrast,sampleList)
        contrast <- unique(contrast)
    }
    if ("bayseq" %in% statistics) libsizeList <- checkLibsize(libsizeList,
        sampleList)
    # Check the genomic database version argument
    if (is.character(version)) {
        version <- tolower(version)
        checkTextArgs("version",version,c("auto"),multiarg=FALSE)
    }
    else
        checkNumArgs("version",version,"numeric")
    
    # Check utrOpts if countType is "utr"
    if (countType == "utr") {
		utrOptsDef <- getDefaults("utrOpts")
		if (!is.null(utrOpts$frac)) {
			checkNumArgs("utrOpts$frac",utrOpts$frac,"numeric",c(0,1),"both")
			utrOptsDef$frac <- utrOpts$frac
		}
		if (!is.null(utrOpts$minLength)) {
			checkNumArgs("utrOpts$minLength",utrOpts$minLength,"numeric",0,
				"gt")
			utrOptsDef$minLength <- utrOpts$minLength
		}
		if (!is.null(utrOpts$downstream)) {
			checkNumArgs("utrOpts$downstream",utrOpts$downstream,"numeric",0,
				"gte")
			utrOptsDef$downstream <- utrOpts$downstream
		}
		utrOpts <- utrOptsDef
	}
    
    # Check main functionality packages
    checkPackages(metaP,qcPlots)
    # Check the case of embedded annotation, not given gc and gene name columns
    # Checks about countType have been performed before
    if (!is.null(annotation) && annotation == "embedded") {
        if (is.na(embedCols$gcCol) && countType=="gene" 
			&& normalization=="edaseq")
            stopwrap("The column that contains the gene GC content ",
                "(\"embedCols$gcCol\") argument is required when ",
                "\"annotation\" is \"embedded\"!")
                
        if (is.na(embedCols$nameCol) && !is.na(geneFilters$expression$known)) {
            warnwrap("The column that contains the HUGO gene symbols ",
                "(\"embedCols$nameCol\") is missing with embedded annotation! ",
                "Gene name expression filter will not be available...")
            geneFilters$expression$known=NA
            if ("volcano" %in% qcPlots)
                warnwrap("The column that contains the HUGO gene symbols ",
                    "(\"embedCols$nameCol\") is missing with embedded ",
                    "annotation! Interactive volcano plots will not contain ",
                    "gene names...")
        }
        
        if (is.na(embedCols$btCol) && countType=="gene") {
            warnwrap("The column that contains the gene biotypes ",
				"(\"embedCols$btCol\") is missing with embedded annotation! ",
                "Biotype filters and certain plots will not be available...")
            geneFilters$biotype=NULL
            toRemove <- match(c("biodetection","countsbio","saturation",
                "biodist","filtered","rnacomp","readnoise"),qcPlots)
            noMatch <- which(is.na(toRemove))
            if (length(noMatch)>0)
                toRemove <- toRemove[-noMatch]
            if (length(toRemove)>0)
                qcPlots <- qcPlots[-toRemove]
        }
    }
    if (org %in% c("hg18","hg38","mm10") && (refdb %in% c("ucsc","refseq"))) {
        warnwrap("Gene/exon biotypes cannot be retrieved when organism is ",
            "\"hg18\", \"hg38\", \"mm10\" and annotation database is ",
            "\"ucsc\" or \"refseq\"!\nBiotype filters and certain plots will ",
            "not be available...")
        geneFilters$biotype=NULL
        toRemove <- match(c("biodetection","countsbio","saturation",
            "biodist","filtered"),qcPlots)
        noMatch <- which(is.na(toRemove))
        if (length(noMatch)>0)
            toRemove <- toRemove[-noMatch]
        if (length(toRemove)>0)
            qcPlots <- qcPlots[-toRemove]
    }
    
    # Check if drawing a Venn diagram is possible
    if ("venn" %in% qcPlots && length(statistics)==1) {
        warnwrap("The creation of a Venn diagram is possible only when more ",
            "than one statistical algorithms are used (meta-analysis)! ",
            "Removing from figures list...")
        toRemove <- match("venn",qcPlots)
        noMatch <- which(is.na(toRemove))
        if (length(noMatch)>0)
            toRemove <- toRemove[-noMatch]
        if (length(toRemove)>0)
            qcPlots <- qcPlots[-toRemove]
    }
    
    # Check additional input arguments for normalization and statistics
    algArgs <- validateAlgArgs(normalization,statistics,normArgs,statArgs)
    normArgs <- algArgs$normArgs
    statArgs <- algArgs$statArgs
    
    # Override settings if a preset is given
    if (!is.null(preset)) {
        presetOpts <- getPresetOpts(preset,org)
        exonFilters <- presetOpts$exonFilters
        geneFilters <- presetOpts$geneFilters
        pcut <- presetOpts$pcut
        exportWhat <- presetOpts$exportWhat
        exportScale <- presetOpts$exportScale
        exportValues <- presetOpts$exportValues
        exportStats <- presetOpts$exportStats
    }

    if (report) {
        reportMessages <- makeReportMessages("en")
        if (!is.null(qcPlots) && !("png" %in% figFormat)) {
            warnwrap("png format is required in order to build a report! ",
                "Adding to figure output formats...")
            figFormat <- c(figFormat,"png")
        }
    }
    
    # Check for replicates. DESeq2 cannot perform without them and ABSSeq cannot
    # perform complex model analysis if not provided with replicates.
    if ("deseq2" %in% statistics && any(lengths(sampleList) < 2)) {	
		stopwrap("DESeq2 does not support analysis without ",
			"replication. Please choose another algorithm or ",
			"provide replicates.")
    }
    
    if (!is.null(contrast) && "absseq" %in% statistics) {        
		contrastListRepCheck <- makeContrastList(contrast,sampleList)
		for (conName in names(contrastListRepCheck)) {
			con <- contrastListRepCheck[[conName]]
			cons <- unique(unlist(con))
			# Checking if we have more than 2 conditions to compare 
			# simultaneously
			if (length(cons)>2) {
				if (any(lengths(sampleList) < 2))
					# Checking if there are no replicates
					stopwrap("ABSSeq cannot perform complex design",
						"analysis without replicates.")
			}
		}
    }
    
    # Display initialization report
    TB <- Sys.time()
    disp(strftime(Sys.time()),": Data processing started...\n")
    ############################################################################
    if ("dss" %in% statistics && any(lengths(sampleList)) < 2) {
        message=c(
            "\n=================================================\n", 
            "DSS will be run without replicates.\n",
            "Results must be interpreted with caution.\n",
            "===================================================\n"
        )
        disp(message)
    }
    
    disp("Read counts file: ",countsName)
    disp("Conditions: ",paste(names(sampleList),collapse=", "))
    disp("Samples to include: ",paste(unlist(sampleList),collapse=", "))
    if (!is.null(excludeList) && !is.na(excludeList))
        disp("Samples to exclude: ",paste(unlist(excludeList),collapse=", "))
    else
        disp("Samples to exclude: none")
    disp("Requested contrasts: ",paste(contrast,collapse=", "))
    if (!is.null(libsizeList)) {
        disp("Library sizes: ")
        for (n in names(libsizeList))
            disp("  ",paste(n,libsizeList[[n]],sep=": "))
    }
    disp("Annotation: ",annotation)
    disp("Organism: ",org)
    disp("Reference source: ",refdb)
    disp("Count type: ",countType)
    if (countType == "utr")
        disp("3' UTR flanking: ",utrFlank)
    if (!is.null(preset))
        disp("Analysis preset: ",preset)
    disp("Transcriptional level: ",transLevel)
    if (!is.null(exonFilters)) {
        disp("Exon filters: ",paste(names(exonFilters),collapse=", "))
        for (ef in names(exonFilters)) {
            disp("  ",ef,": ")
            for (efp in names(exonFilters[[ef]])) {
                if (length(exonFilters[[ef]][[efp]])==1 && 
                    is.function(exonFilters[[ef]][[efp]]))
                    print(exonFilters[[ef]][[efp]])
                else if (length(exonFilters[[ef]][[efp]])==1)
                    disp("    ",paste(efp,exonFilters[[ef]][[efp]],sep=": "))
                else if (length(exonFilters[[ef]][[efp]])>1)
                    disp("    ",paste(efp,paste(exonFilters[[ef]][[efp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Exon filters: none applied")
    if (!is.null(geneFilters)) {
        disp("Gene filters: ",paste(names(geneFilters),collapse=", "))
        for (gf in names(geneFilters)) {
            disp("  ",gf,": ")
            for (gfp in names(geneFilters[[gf]])) {
                if (length(geneFilters[[gf]][[gfp]])==1 && 
                    is.function(geneFilters[[gf]][[gfp]]))
                    print(geneFilters[[gf]][[gfp]])
                else if (length(geneFilters[[gf]][[gfp]])==1)
                    disp("    ",paste(gfp,geneFilters[[gf]][[gfp]],sep=": "))
                else if (length(geneFilters[[gf]][[gfp]])>1)
                    disp("    ",paste(gfp,paste(geneFilters[[gf]][[gfp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Gene filters: none applied")
    disp("Filter application: ",whenApplyFilter)
    disp("Normalization algorithm: ",normalization)
    if (!is.null(normArgs)) {
        disp("Normalization arguments: ")
        for (na in names(normArgs)) {
            if (length(normArgs[[na]])==1 && is.function(normArgs[[na]])) {
                disp("  ",na,": ")
                disp(as.character(substitute(normArgs[[na]])))
            }
            else if (length(normArgs[[na]])==1)
                disp("  ",paste(na,normArgs[[na]],sep=": "))
            else if (length(normArgs[[na]])>1)
                disp("  ",paste(na,paste(normArgs[[na]],collapse=", "),
                    sep=": "))
        }
    }
    disp("Statistical algorithm: ",paste(statistics,collapse=", "))
    if (!is.null(statArgs)) {
        disp("Statistical arguments: ")
        for (sa in names(statArgs)) {
            if (length(statArgs[[sa]])==1 && is.function(statArgs[[sa]])) {
                disp("  ",sa,": ")
                disp(as.character(substitute(statArgs[[na]])))
            }
            else if (length(statArgs[[sa]])==1)
                disp("  ",paste(sa,statArgs[[sa]],sep=": "))
            else if (length(statArgs[[sa]])>1)
                disp("  ",paste(sa,paste(statArgs[[sa]],collapse=", "),
                    sep=": "))
        }
    }
    disp("Meta-analysis method: ",metaP)
    disp("Multiple testing correction: ",adjustMethod)
    if (!is.na(pcut)) 
        disp("p-value threshold: ",pcut)
    disp("Logarithmic transformation offset: ",logOffset)
    if (!is.null(preset)) 
        disp("Analysis preset: ",preset)
    disp("Quality control plots: ",paste(qcPlots,collapse=", "))
    disp("Figure format: ",paste(figFormat,collapse=", "))
    if (!is.na(exportWhere)) 
        disp("Output directory: ",exportWhere)
    disp("Output data: ",paste(exportWhat,collapse=", "))
    disp("Output scale(s): ",paste(exportScale,collapse=", "))
    disp("Output values: ",paste(exportValues,collapse=", "))
    if ("stats" %in% exportWhat)
        disp("Output statistics: ",paste(exportStats,collapse=", "),"\n")
        
    if (is.function(progressFun)) {
        text <- paste("Starting the analysis...")
        progressFun(detail=text)
    }
    ############################################################################
	
	# Somewhere here we must load and/or construct annotation
	
	# geneData are always required, unless annotation is embedded, where we
	# need only geneData as embedded annotation is not allowed anywhere else
	disp("Loading gene annotation...")
	if (!is.null(annotation) && annotation == "embedded") {
		# The following should work if annotation elements are arranged 
		# in MeV-like data style
		if (!is.data.frame(counts)) {
			disp("Reading counts file ",countsName,"...")
			geneCounts <- read.delim(counts)
		}
		else
			geneCounts <- counts
		rownames(geneCounts) <- as.character(geneCounts[,embedCols$idCol])
		allCols <- 1:ncol(geneCounts)
		samCols <- match(unlist(sampleList),colnames(geneCounts))
		samCols <- samCols[which(!is.na(samCols))]
		annCols <- allCols[-samCols]
		geneData <- geneCounts[,annCols]
		geneCounts <- geneCounts[,samCols]		
		colnames(geneData)[embedCols$idCol] <- "gene_id"
		if (!is.na(embedCols$gcCol)) {
			colnames(geneData)[embedCols$gcCol] <- "gc_content"
			if (max(geneData$gc_content<=1)) # Is already divided
				geneData$gc_content = 100*geneData$gc_content
		}
		if (!is.na(embedCols$nameCol)) 
			colnames(geneData)[embedCols$nameCol] <- "gene_name"
		if (!is.na(embedCols$btCol))
			colnames(geneData)[embedCols$btCol] <- "biotype"
		geneData <- GRanges(geneData)
		names(geneData) <- geneData$gene_id
	}
	else {
		geneData <- tryCatch(loadAnnotation(org,refdb,level=transLevel,
			type="gene",version=version,db=localDb,rc=restrictCores),
			error=function(e) {
				# Not found and user based
				if (!is.null(annotation)) { # Has been checked
					gtfFile <- annotation$gtf
					metadata <- annotation
					metadata$gtf <- NULL
					geneData <- importCustomAnnotation(gtfFile,metadata,
						"gene","gene")
				}
				else
					stop("Please provide an existing organism or a list ",
						"with annotation metadata and GTF file!")
			},finally="")
	}
	
	# Now start examining additional required data per case...
    if (countType=="exon") {
        # We need to load the exon annotation and see what counts have been
        # provided
        if (!fromPrevious) {
			# Load exon annotation
			disp("Loading exon annotation...")
			exonData <- tryCatch(loadAnnotation(org,refdb,level=transLevel,
				type=countType,version=version,db=localDb,summarized=TRUE,
				rc=restrictCores),
			error=function(e) {
				# Not found and user based
				if (!is.null(annotation)) { # Has been checked
					gtfFile <- annotation$gtf
					metadata <- annotation
					metadata$gtf <- NULL
					exonData <- importCustomAnnotation(gtfFile,metadata,
						"gene","exon")
				}
				else
					stop("Please provide an existing organism or a list ",
						"with annotation metadata and GTF file!")
			},finally="")
			
            # Load/read counts
            if (!is.null(counts)) { # Provided
                if (!is.data.frame(counts) && !is.list(counts)) {
                    disp("Reading counts file ",countsName,"...")
                    exonCounts <- read.delim(counts,row.names=1)
                }
                else { # Already a data frame as input
					if (is.character(counts[,1])) {
						exonCounts <- as.matrix(counts[,-1])
						rownames(exonCounts) <- as.character(counts[,1])
					}
                    else { # Should be named!
						if (is.null(rownames(counts)))
							stopwrap("A counts data frame as input should ",
								"have rownames!")
						exonCounts <- counts
					}
                }
            }
            else { # Coming from read2count
                if (fromRaw) { # Double check
                    r2c <- read2count(theList,exonData,fileType,
                        utrFlank,rc=restrictCores)
                    exonCounts <- r2c$counts
                    # Merged exon data! - Not needed anymore
                    #exonData <- r2c$mergedann
                    if (is.null(libsizeList))
                        libsizeList <- r2c$libsize
                    if (exportCountsTable) {
						exonDataExp <- as.data.frame(exonData)
						exonDataExp <- exonDataExp[,c(1:3,6,7,5,8,9)]
						names(exonDataExp)[1] <- "chromosome"
						disp("Exporting raw read counts table to ",
                            file.path(PROJECT_PATH[["lists"]],
                                "raw_counts_table.txt.gz"))
                        resFile <- file.path(PROJECT_PATH[["lists"]],
                            "raw_counts_table.txt.gz")
                        gzfh <- gzfile(resFile,"w")
                        write.table(cbind(
                            exonDataExp[rownames(exonCounts),],
                            exonCounts
						),gzfh,sep="\t",row.names=FALSE,quote=FALSE)
                        close(gzfh)
                    }
                }
            }
            exonCounts <- exonCounts[,unlist(sampleList,use.names=FALSE)]

            # Get the exon counts per gene model
            disp("Checking chromosomes in exon counts and gene annotation...")
            geneData <- .reduceGeneData(exonData[rownames(exonCounts)],geneData)
            disp("Processing exons...")
            theCounts <- constructGeneModel(exonCounts,exonData,type="exon",
				rc=restrictCores)

            if (saveGeneModel) {
                disp("Saving gene model to ",file.path(PROJECT_PATH[["data"]],
                    "gene_model.RData"))
                save(theCounts,exonData,geneData,sampleList,transLevel,
                    countType,file=file.path(PROJECT_PATH$data,
					"gene_model.RData"),compress=TRUE)
            }
        }
        else {
            theCounts <- tmpEnv$theCounts
            exonData <- tmpEnv$exonData
            geneData <- tmpEnv$geneData
        }

        # Exclude any samples not wanted (when e.g. restoring a previous project
        # and having determined that some samples are of bad quality
        if (!is.null(excludeList) && !is.na(excludeList)) {
            for (n in names(excludeList)) {
                sampleList[[n]] <- setdiff(sampleList[[n]],
                    excludeList[[n]])
                if (length(sampleList[[n]])==0) # Removed whole condition
                    sampleList[n] <- NULL
            }
            theCounts <- theCounts[unlist(sampleList)]
        }

        # Apply exon filters
        if (!is.null(exonFilters)) {
			exonFilterOut <- filterExons(theCounts,geneData,sampleList,
                exonFilters)
            exonFilterResult <- exonFilterOut$result
            exonFilterFlags <- exonFilterOut$flags
        }
        else
            exonFilterResult <- exonFilterFlags <- NULL
        
        disp("Summarizing count data...")
        theGeneCounts <- theExonLengths <- vector("list",
            length(unlist(sampleList)))
        names(theGeneCounts) <- names(theExonLengths) <- names(theCounts)
        for (n in names(theGeneCounts))
            theGeneCounts[[n]] <- sapply(theCounts[[n]],sum)
        geneCounts <- do.call("cbind",theGeneCounts)
        # Based on the sum of their exon lengths
        lengthList <- attr(theCounts,"lengthList")
        geneLength <- sapply(lengthList,sum)
        # Could also be
        # geneLength <- attr(exonData,"activeLength")
        
        # In case there are small differences between annotation data and  
        # external file, due to e.g. slightly different Ensembl versions
        geneData <- geneData[rownames(geneCounts)]
        totalGeneData <- geneData # We need this for some total stats
    }
    
    if (countType=="utr") {
        # We need to load the utr annotation and see what counts have been
        # provided
        if (!fromPrevious) {
			disp("Loading 3' UTR annotation...")
			transcriptData <- tryCatch(loadAnnotation(org,refdb,
				level=transLevel,type=countType,version=version,
				summarized=TRUE,rc=restrictCores),error=function(e) {
				# Not found and user based
				if (!is.null(annotation)) { # Has been checked
					gtfFile <- annotation$gtf
					metadata <- annotation
					metadata$gtf <- NULL
					transcriptData <- importCustomAnnotation(gtfFile,metadata,
						transLevel,countType)
				}
				else
					stop("Please provide an existing organism or a list with ",
						"annotation metadata and GTF file!")
			},finally="")
			
            # Load/read counts
            if (!is.null(counts)) { # Otherwise coming ready from read2count
                if (!is.data.frame(counts) && !is.list(counts)) {
                    disp("Reading counts file ",countsName,"...")
                    transcriptCounts <- read.delim(counts)
                }
                else { # Already a data frame as input
                    if (is.character(counts[,1])) {
						transcriptCounts <- as.matrix(counts[,-1])
						rownames(transcriptCounts) <- as.character(counts[,1])
					}
                    else { # Should be named!
						if (is.null(rownames(counts)))
							stopwrap("A counts data frame as input should ",
								"have rownames!")
						transcriptCounts <- counts
					}
				}
            }
            else { # Coming from read2count
                if (fromRaw) { # Double check
                    r2c <- read2count(theList,transcriptData,fileType,
                        transLevel,utrFlank,rc=restrictCores)
                    transcriptCounts <- r2c$counts
                    if (is.null(libsizeList))
                        libsizeList <- r2c$libsize
                    if (exportCountsTable) {
						transcriptDataExp <- as.data.frame(transcriptData)
						transcriptDataExp <- 
							transcriptDataExp[,c(1:3,6,7,5,8,9)]
						names(transcriptDataExp)[1] <- "chromosome"
                        disp("Exporting raw read counts table to ",
                            file.path(PROJECT_PATH[["lists"]],
                            "raw_counts_table.txt.gz"))
                        resFile <- file.path(PROJECT_PATH[["lists"]],
                            "raw_counts_table.txt.gz")
                        gzfh <- gzfile(resFile,"w")
                        write.table(cbind(
							transcriptDataExp[rownames(transcriptCounts),],
							transcriptCounts
						),gzfh,sep="\t",row.names=FALSE,quote=FALSE)
                        close(gzfh)
                    }
                }
            }
            transcriptCounts <- transcriptCounts[,unlist(sampleList,
				use.names=FALSE)]

            # Get the transcript counts per gene model
            disp("Checking chromosomes in transcript counts and gene ",
                "annotation...")
            geneData <- 
                .reduceGeneData(transcriptData[rownames(transcriptCounts)],
                    geneData)
            disp("Processing transcripts...")
            theCounts <- constructGeneModel(transcriptCounts,transcriptData,
				type="utr",rc=restrictCores)
			
            if (saveGeneModel) {
                disp("Saving gene model to ",file.path(PROJECT_PATH[["data"]],
                    "gene_model.RData"))
                save(theCounts,transcriptData,geneData,sampleList,countType,
					transLevel,file=file.path(PROJECT_PATH$data,
					"gene_model.RData"),compress=TRUE)
            }
        }
        else {
            theCounts <- tmpEnv$theCounts
            transcriptData <- tmpEnv$transcriptData
            geneData <- tmpEnv$geneData
        }

        # Exclude any samples not wanted (when e.g. restoring a previous project
        # and having determined that some samples are of bad quality
        if (!is.null(excludeList) && !is.na(excludeList)) {
            for (n in names(excludeList)) {
                sampleList[[n]] <- setdiff(sampleList[[n]],
                    excludeList[[n]])
                if (length(sampleList[[n]])==0) # Removed whole condition
                    sampleList[n] <- NULL
            }
            theCounts <- theCounts[unlist(sampleList)]
        }

        disp("Summarizing count data...")
        theGeneCounts <- theTranscriptLengths <- vector("list",
            length(unlist(sampleList)))
        names(theGeneCounts) <- names(theTranscriptLengths) <- 
            names(theCounts)
        for (n in names(theGeneCounts))
			theGeneCounts[[n]] <- sapply(theCounts[[n]],sum)
        geneCounts <- do.call("cbind",theGeneCounts)
        # Based on the sum of their transcript lengths
        lengthList <- attr(theCounts,"lengthList")
        geneLength <- sapply(lengthList,sum)
        
        # In case there are small differences between annotation data and 
        # external file, due to e.g. slightly different Ensembl versions
        geneData <- geneData[rownames(geneCounts),]
        totalGeneData <- geneData # We need this for some total stats
        
        # No exon filters have been applied, we are at "utr" level
        exonFilterResult <- NULL
    }
    
    else if (countType=="gene") {
		# geneData has already been loaded and also geneCounts in the case of
		# embedded annotation
        if (!fromPrevious) {
            # Load/read counts
            if (!is.null(counts)  && !is.list(counts)) {
                if (!is.data.frame(counts)) { # Else it's already here
                    disp("Reading counts file ",countsName,"...")
                    geneCounts <- read.delim(counts)
                }
                else { # Already a data frame as input
					if (is.character(counts[,1])) {
						geneCounts <- as.matrix(counts[,-1])
						rownames(geneCounts) <- as.character(counts[,1])
					}
					else { # Should be named!
						if (is.null(rownames(counts)))
							stopwrap("A counts data frame as input should ",
								"have rownames!")
						geneCounts <- counts
					}
				}
            }
            else { # Coming from read2count
                if (fromRaw) { # Double check
                    r2c <- read2count(theList,geneData,fileType,
                        utrFlank,rc=restrictCores)
                    geneCounts <- r2c$counts
                    if (is.null(libsizeList))
                        libsizeList <- r2c$libsize
                    if (exportCountsTable) {
						geneDataExp <- as.data.frame(geneData)
						geneDataExp <- geneDataExp[,c(1:3,6,7,5,8,9)]
						names(geneDataExp)[1] <- "chromosome"
                        disp("Exporting raw read counts table to ",
                            file.path(PROJECT_PATH[["lists"]],
                            "raw_counts_table.txt.gz"))
                        resFile <- file.path(PROJECT_PATH[["lists"]],
                            "raw_counts_table.txt.gz")
                        gzfh <- gzfile(resFile,"w")
                        write.table(cbind(
							geneDataExp[rownames(geneCounts),],
                            geneCounts
                        ),gzfh,sep="\t",row.names=FALSE,quote=FALSE)
                        close(gzfh)
                    }
                }
            }
        }
        else {
            geneCounts <- tmpEnv$geneCounts
            geneData <- tmpEnv$geneData
        }
        
        totalGeneData <- geneData # We need this for some total stats
        exonFilterResult <- NULL

        geneData <- geneData[rownames(geneCounts)]
        geneLength <- width(geneData)
        names(geneLength) <- names(geneData)
        
        # Exclude any samples not wanted (when e.g. restoring a previous project
        # and having determined that some samples are of bad quality
        if (!is.null(excludeList) && !is.na(excludeList)) {
            for (n in names(excludeList)) {
                sampleList[[n]] <- setdiff(sampleList[[n]],
                    excludeList[[n]])
                if (length(sampleList[[n]])==0) # Removed whole condition
                    sampleList[n] <- NULL
            }
            geneCounts <- geneCounts[,unlist(sampleList,use.names=FALSE)]
        }
        
        if (saveGeneModel) {
            disp("Saving gene model to ",file.path(PROJECT_PATH[["data"]],
                "gene_model.RData"))
            save(geneCounts,geneData,sampleList,countType,
                file=file.path(PROJECT_PATH$data,"gene_model.RData"),
                compress=TRUE)
        }
    }

    # Transform GC-content and biotype - should never be required in the new
    # annotation strategy
    if (is.null(geneData$gc_content))
        geneData$gc_content <- rep(0.5,length(geneData))
    if (is.null(geneData$biotype))
        geneData$biotype <- rep("gene",length(geneData))
    names(geneLength) <- rownames(geneCounts)
    attr(geneData,"geneLength") <- geneLength

    ############################################################################
    # BEGIN FILTERING SECTION
    ############################################################################
    
    if (is.function(progressFun)) {
        text <- paste("Filtering...")
        progressFun(detail=text)
    }
    
    # GC bias is NOT alleviated if we do not remove the zeros!!!
    disp("Removing genes with zero counts in all samples...")
    theZeros <- which(apply(geneCounts,1,filterLow,0))
    if (length(theZeros) > 0) {
        # Store the filtered, maybe we do some stats
        geneCountsZero <- geneCounts[theZeros,,drop=FALSE]
        geneDataZero <- geneData[theZeros]
        attr(geneDataZero,"geneLength") <- geneLength[theZeros]
        theZeroNames <- names(geneData)[theZeros]
        # Then remove
        geneCounts <- geneCounts[-theZeros,,drop=FALSE]
        geneData <- geneData[-theZeros]
        attr(geneData,"geneLength") <- geneLength[-theZeros]
    }
    else
        geneCountsZero <- geneDataZero <- theZeroNames <- NULL

    # Store un-normalized gene counts for export purposes
    geneCountsUnnorm <- geneCounts

    # Apply filtering prior to normalization if desired
    if (whenApplyFilter=="prenorm") {
        # However, a first round of normalization has to be performed in
        # order to get proper expression filters
        disp("Prefiltering normalization with: ",normalization)
        switch(normalization,
            edaseq = {
                tempGenes <- normalizeEdaseq(geneCounts,sampleList,
                    normArgs,geneData,output="matrix")
            },
            deseq = {
                tempGenes <- normalizeDeseq(geneCounts,sampleList,normArgs,
                    output="matrix")
            },
            deseq2 = {
                tempGenes <- normalizeDeseq2(geneCounts,sampleList,normArgs,
                    output="matrix")
            },   
            edger = {
                tempGenes <- normalizeEdger(geneCounts,sampleList,normArgs,
                    output="matrix")
            },
            noiseq = {  
                tempGenes <- normalizeNoiseq(geneCounts,sampleList,
                    normArgs,geneData,logOffset,output="matrix")
            },
            nbpseq = {
                tempGenes <- normalizeNbpseq(geneCounts,sampleList,
                    normArgs,libsizeList,output="matrix")
            },
            absseq = {
                tempGenes <- normalizeAbsseq(geneCounts,sampleList,normArgs,
                    output="matrix")
            },
			dss = {
				tempGenes <- normalizeDss(geneCounts,sampleList,normArgs,
					output="native")
			},
            none = {
                # In case some external normalization is applied (e.g. equal
                # read counts from all samples)
                tempGenes <- geneCounts
            }
        )
        
        # Now filter
        if (!is.null(geneFilters)) {
            geneFilterOut <- filterGenes(tempGenes,geneData,geneFilters,
                sampleList)
            geneFilterResult <- geneFilterOut$result
            geneFilterCutoff <- geneFilterOut$cutoff
            geneFilterFlags <- geneFilterOut$flags
        }
        else
            geneFilterResult <- geneFilterCutoff <- geneFilterFlags <- NULL

        # Unify the filters and filter
        theDeadGenes <- list(
            geneFilterResult$expression$median,
            geneFilterResult$expression$mean,
            geneFilterResult$expression$quantile,
            geneFilterResult$expression$known,
            geneFilterResult$expression$custom
        )
        theDead <- unique(unlist(c(geneFilterResult,exonFilterResult)))
        # Some genes filtered by zero, were present in exon filters, not yet applied
        if (countType=="exon")
            theDead <- setdiff(theDead,theZeroNames)
        
        # All method specific objects are row-index subsettable
        if (length(theDead)>0) {
            # Store the filtered for later export or some stats
            geneCountsDead <- geneCounts[theDead,]
            geneCountsUnnorm <- geneCountsUnnorm[theDead,]
            geneDataDead <- geneData[theDead]
            attr(geneDataDead,"geneLength") <- attr(geneData,
                "geneLength")[theDead]
            # Now filter
            theDeadInd <- match(theDead,rownames(geneCounts))
            geneCountsExpr <- geneCounts[-theDeadInd,]
            geneDataExpr <- geneData[-theDeadInd]
            attr(geneDataExpr,"geneLength") <- attr(geneData,
                "geneLength")[-theDeadInd]
        }
        else {
            geneCountsExpr <- geneCounts
            geneDataExpr <- geneData
            geneCountsDead <- geneDataDead <- geneCountsUnnorm <- NULL
        }
        
        if (is.function(progressFun)) {
            text <- paste("Normalizing...")
            progressFun(detail=text)
        }
        
        disp("Normalizing with: ",normalization)
        switch(normalization,
            edaseq = {
                normGenes <- normalizeEdaseq(geneCountsExpr,sampleList,normArgs,
					geneDataExpr,output="matrix")
            },
            deseq = {
                normGenes <- normalizeDeseq(geneCountsExpr,sampleList,normArgs,
                    output="native")
            },
            deseq2 = {
                normGenes <- normalizeDeseq2(geneCountsExpr,sampleList,normArgs,
                    output="native")
            },
            edger = {
                normGenes <- normalizeEdger(geneCountsExpr,sampleList,normArgs,
                    output="native")
            },
            noiseq = {
                normGenes <- normalizeNoiseq(geneCountsExpr,sampleList,normArgs,
                    geneDataExpr,logOffset,output="matrix")
            },
            nbpseq = {
                normGenes <- normalizeNbpseq(geneCountsExpr,sampleList,normArgs,
                    libsizeList,output="native")
            },
            absseq = {
                normGenes <- normalizeAbsseq(geneCountsExpr,sampleList,normArgs,
                    output="native")
            },
			dss = {
                normGenes <- normalizeDss(geneCountsExpr,sampleList,normArgs,
                    output="native")
            },   
            none = {
                normGenes <- geneCountsExpr
            }
        )
        normGenesExpr <- normGenes
    }
    else if (whenApplyFilter=="postnorm") {
        if (is.function(progressFun)) {
            text <- paste("Normalizing...")
            progressFun(detail=text)
        }
        
     # Apply filtering after normalization if desired (default)
     disp("Normalizing with: ",normalization)   
        switch(normalization,
            edaseq = {
                normGenes <- normalizeEdaseq(geneCounts,sampleList,
                    normArgs,geneData,output="matrix")
            },
            deseq = {
                normGenes <- normalizeDeseq(geneCounts,sampleList,normArgs,
                    output="native")
            },
            deseq2 = {
                normGenes <- normalizeDeseq2(geneCounts,sampleList,
                    normArgs,output="native")
            },
            edger = {
                normGenes <- normalizeEdger(geneCounts,sampleList,normArgs,
                    output="native")
            },
            noiseq = {
                normGenes <- normalizeNoiseq(geneCounts,sampleList,
                    normArgs,geneData,logOffset,output="matrix")
            },
            nbpseq = {
                normGenes <- normalizeNbpseq(geneCounts,sampleList,
                    normArgs,libsizeList,output="native")
            },
            absseq = {
                normGenes <- normalizeAbsseq(geneCounts,sampleList,
                    normArgs,output="native")
            },
			dss = {
                normGenes <- normalizeDss(geneCounts,sampleList,normArgs,
                    output="native")
            },
            none = {
                normGenes <- geneCounts
            }
        )
        
        switch(class(normGenes), 
            CountDataSet = { # Has been normalized with DESeq
                tempMatrix <- round(DESeq::counts(normGenes,normalized=TRUE))
            },
            DESeqDataSet = { # Has been normalized with DESeq2
                tempMatrix <- round(DESeq2::counts(normGenes,normalized=TRUE))
            },
            DGEList = { # Has been normalized with edgeR
                # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
                scl <- normGenes$samples$lib.size *
                    normGenes$samples$norm.factors
                tempMatrix <- round(t(t(normGenes$counts)/scl)*mean(scl))
            },
            matrix = { # Has been normalized with EDASeq or NOISeq or nothing
                tempMatrix <- normGenes
            },
            data.frame = { # Has been normalized with or nothing
                tempMatrix <- as.matrix(normGenes)
            },
            nbData = { # Has been normalized with NBPSeq and main was "nbpseq"
                tempMatrix <- as.matrix(round(sweep(normGenes$counts,2,
                    normGenes$norm.factors,"*")))
            },
            nbp = { # Has been normalized with NBPSeq and main was "nbsmyth"
                 tempMatrix <- as.matrix(round(normGenes$pseudo.counts))
            },
            ABSDataSet = { # Has been normalized with ABSSeq
                 tempMatrix <- as.matrix(round(excounts(normGenes)))
            },
			SeqCountSet = { # Has been normalized with DSS
				# Dribble for taking a mtx out of SeqCountSet class
				classes <- asClassVector(sampleList)
				theDesign <- data.frame(condition=classes,
					row.names=colnames(normGenes)) 
            	cds <- newCountDataSet(as.matrix(round(
					assayData(normGenes)$exprs)),theDesign$condition)
            	DESeq::sizeFactors(cds) <- normalizationFactor(normGenes)
				tempMatrix <- as.matrix(round(DESeq::counts(cds,
					normalized=TRUE)))
			}
        )

        # Implement gene filters after normalization
        if (!is.null(geneFilters)) {
            geneFilterOut <- filterGenes(tempMatrix,geneData,geneFilters,
                sampleList)
            geneFilterResult <- geneFilterOut$result
            geneFilterCutoff <- geneFilterOut$cutoff
            geneFilterFlags <- geneFilterOut$flags
        }
        else
            geneFilterResult <- geneFilterCutoff <-
                geneFilterFlags <- NULL

        # Unify the filters and filter
        theDeadGenes <- list(
            geneFilterResult$expression$median,
            geneFilterResult$expression$mean,
            geneFilterResult$expression$quantile,
            geneFilterResult$expression$known,
            geneFilterResult$expression$custom
        )
        #when running maually exonFilterResult=NULL
        theDead <- unique(unlist(c(geneFilterResult,exonFilterResult)))
        # Some genes filtered by zero, were present in exon filters, not yet applied
        if (countType=="exon")
            theDead <- setdiff(theDead,theZeroNames)
        
        # All method specific objects are row-index subsettable
        if (length(theDead) > 0) {
            # Store the filtered for later export or some stats
            geneCountsDead <- tempMatrix[theDead,]
            geneCountsUnnorm <- geneCountsUnnorm[theDead,]
            geneDataDead <- geneData[theDead]
            attr(geneDataDead,"geneLength") <- attr(geneData,
                "geneLength")[theDead]
            # Now filter
            theDeadInd <- match(theDead,rownames(tempMatrix))  
            switch(class(normGenes),                           
                CountDataSet = { # Has been normalized with DESeq
                    normGenesExpr <- normGenes[-theDeadInd,]
                },
                DESeqDataSet = { # Has been normalized with DESeq2
                    normGenesExpr <- normGenes[-theDeadInd,]
                },
                DGEList = { # edgeR bug???
                    normGenesExpr <- normGenes[-theDeadInd,]
                    normGenesExpr$AveLogCPM <-
                        normGenesExpr$AveLogCPM[-theDeadInd]
                },
                matrix = { # Has been normalized with EDASeq or NOISeq
                    normGenesExpr <- normGenes[-theDeadInd,]
                },
                data.frame = { # Has been normalized with EDASeq or NOISeq
                    normGenesExpr <- as.matrix(normGenes[-theDeadInd,])
                },
                nbData = { # Has been normalized NBPSeq, main.method="nbpseq"
                    normGenesExpr <- normGenes
                    normGenesExpr$counts <-
                        as.matrix(normGenesExpr$counts[-theDeadInd,])
                    normGenesExpr$rel.frequencies <- 
                        normGenesExpr$rel.frequencies[-theDeadInd,]
                    normGenesExpr$tags <-
                        as.matrix(normGenesExpr$tags[-theDeadInd,])
                },
                nbp = {
                    normGenesExpr <- normGenes
                    normGenesExpr$counts <-
                        as.matrix(normGenesExpr$counts[-theDeadInd,])
                    normGenesExpr$pseudo.counts <- 
                        as.matrix(normGenesExpr$pseudo.counts[-theDeadInd,])
                    normGenesExpr$pseudo.libSizes <- 
                        colSums(as.matrix(normGenesExpr$pseudo.counts))*
                            rep(1,dim(normGenesExpr$counts)[2])
                },
                ABSDataSet = { 
                    normGenesExpr <- normGenes
                    ABSSeq::counts(normGenesExpr) <-
                        ABSSeq::counts(normGenesExpr)[-theDeadInd,]
                    excounts(normGenesExpr) <-
                        excounts(normGenesExpr)[-theDeadInd,]
                },
				SeqCountSet = { 
					normGenesExpr <- normGenes
					# Subset raw counts table
					normGenesExpr <- normGenesExpr[-theDeadInd,] 
					# All the other tools here return a matrix. DSS must do the 
					# same. --> otherwise error in rbind of line 3192
					classes <- asClassVector(sampleList)
				    theDesign <- data.frame(condition=classes,
						row.names=colnames(normGenesExpr)) 
					cds <- newCountDataSet(as.matrix(round(
						assayData(normGenesExpr)$exprs)),theDesign$condition)
            	    DESeq::sizeFactors(cds) <- 
						normalizationFactor(normGenesExpr)
				 	normGenesExpr <- as.matrix(round(DESeq::counts(cds,
						normalized=TRUE)))
				}
            )
            geneCountsExpr <- geneCounts[rownames(normGenesExpr),]
            geneDataExpr <- geneData[-theDeadInd]
            attr(geneDataExpr,"geneLength") <-
                attr(geneData,"geneLength")[-theDeadInd]
            
        }
        else {
            normGenesExpr <- normGenes
            geneCountsExpr <- geneCounts
            geneDataExpr <- geneData
            geneCountsDead <- geneDataDead <- geneCountsUnnorm <- NULL
        }
    }
    
    # Store the final filtered, maybe we do some stats
    geneDataFiltered <- c(geneDataZero,geneDataDead)
    if (!is.null(geneDataFiltered) && length(geneDataFiltered) > 0) {
        disp(length(geneDataFiltered)," genes filtered out")
        if (!is.null(geneDataZero) && length(geneDataZero) > 0)
            attr(geneDataFiltered,"geneLength") <- c(attr(geneDataZero,
                "geneLength"),attr(geneDataDead,"geneLength"))
        else
            attr(geneDataFiltered,"geneLength") <-
                attr(geneDataDead,"geneLength")
    }
    if (!is.null(geneFilters) || !is.null(exonFilters))
        disp(length(geneDataExpr)," genes remain after filtering")

    ############################################################################
    # END FILTERING SECTION
    ############################################################################

    # There is a small case that no genes are left after filtering...
    if(any(dim(normGenesExpr)==0))
        stopwrap("No genes left after gene and/or exon filtering! Try again ",
            "with no filtering or less strict filter rules...")

    if (is.function(progressFun)) {
        text <- paste("Statistical testing...")
        progressFun(detail=text)
    }
    
    # Run the statistical test, normGenes is always a method-specific object,
    # handled in the metaseqr.stat.R stat.* functions
    cpList <- vector("list",length(contrast))
    names(cpList) <- contrast
    contrastList <- makeContrastList(contrast,sampleList)
    for (n in names(cpList)) {
        cpList[[n]] <- vector("list",length(statistics))
        names(cpList[[n]]) <- statistics
    }
    for (alg in statistics) {
        disp("Running statistical tests with: ",alg)    
        switch(alg,
            deseq = {
                pList <- statDeseq(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            deseq2 = {
                pList <- statDeseq2(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },   
            edger = {
                pList <- statEdger(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            noiseq = {
                pList <- statNoiseq(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]],geneDataExpr,logOffset)
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            bayseq = {
                pList <- statBayseq(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]],libsizeList)
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            limma = {
                pList <- statLimma(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            nbpseq = {
                pList <- statNbpseq(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]],libsizeList)
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            absseq = {
                pList <- statAbsseq(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            },
            dss = {
                pList <- statDss(normGenesExpr,sampleList,contrastList,
                    statArgs[[alg]])
                # Order pList genes as in normGenesExpr, because dss outputs 
                # gene names in a different order
                pList[[1]] <- pList[[1]][rownames(normGenesExpr)]
                if (!is.na(pcut)) {
                    for (con in names(contrastList))
                        disp("  Contrast ",con,": found ",
                            length(which(pList[[con]]<=pcut))," genes")
                }
            }
        )
        for (n in names(pList))
            cpList[[n]][[alg]] <- pList[[n]]
    }
    for (n in names(cpList))
        cpList[[n]] <- do.call("cbind",cpList[[n]])

    # Create the adjusted p-value matrices (if needed)
    if ("adj_p_value" %in% exportWhat) {
        adjCpList <- cmclapply(cpList,
            function(x,a) return(apply(x,2,p.adjust,a)),adjustMethod,
            rc=restrictCores)
        for (n in names(cpList)) {
            noi <- grep("noiseq",colnames(cpList[[n]]))
            if (length(noi)>0) {
                # DESeq has not run in this case, FDR cannot be calculated
                if (length(strsplit(n,"_vs_")[[1]])==2)
                    adjCpList[[n]][,noi] <- rep(NA,nrow(cpList[[n]]))
            }
        }
    }
    else
        adjCpList <- NULL

    # At this point, all method-specific objects must become matrices for  
    # exporting and plotting
    switch(class(normGenesExpr),
        CountDataSet = { # Has been processed with DESeq
            normGenes <- round(DESeq::counts(normGenes,normalized=TRUE))
            normGenesExpr <- round(DESeq::counts(normGenesExpr,normalized=TRUE))
        },
        DESeqDataSet = { # Has been processed with DESeq2
            normGenes <- round(DESeq2::counts(normGenes,normalized=TRUE))
            normGenesExpr <- round(DESeq2::counts(normGenesExpr,
                normalized=TRUE))
        },
        DGEList = { # Has been processed with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            sclR <- normGenes$samples$lib.size*normGenes$samples$norm.factors
            normGenes <- round(t(t(normGenes$counts)/sclR)*mean(sclR))
            sclN <- normGenesExpr$samples$lib.size *
                normGenesExpr$samples$norm.factors
            normGenesExpr <- round(t(t(normGenesExpr$counts)/sclN) *
                mean(sclN))
        },
        nbData = {
            normGenes <- as.matrix(round(sweep(normGenes$counts,2,
                normGenes$norm.factors,"*")))
            normGenesExpr <- as.matrix(round(sweep(normGenesExpr$counts,2,
                normGenes$norm.factors,"*")))
        },
        nbp = {
            normGenes <- as.matrix(round(normGenes$pseudo.counts))
            normGenesExpr <- as.matrix(round(normGenesExpr$pseudo.counts))
        },
        ABSDataSet = {
            normGenes <- as.matrix(round(excounts(normGenes)))
            normGenesExpr <- as.matrix(round(excounts(normGenesExpr)))
        },
		SeqCountSet = {
			# Way to take normalized counts out of DSS as a matrix
            classes <- asClassVector(sampleList)
			theDesign <- data.frame(condition=classes,
				row.names=colnames(normGenes)) 
        	cds <- newCountDataSet(as.matrix(round(assayData(normGenes)$exprs)),
               	 theDesign$condition)
       		DESeq::sizeFactors(cds) <- normalizationFactor(normGenes)
			normGenes <- as.matrix(round(DESeq::counts(cds,normalized=TRUE)))
           
			cdsExpr <- newCountDataSet(as.matrix(round(assayData(
				normGenesExpr)$exprs)),theDesign$condition)
       		DESeq::sizeFactors(cdsExpr) <- normalizationFactor(normGenesExpr)
			normGenesExpr <- as.matrix(round(DESeq::counts(cdsExpr,
				normalized=TRUE)))
        }
        # We don't need the matrix case
    )

    # Now that everything is a matrix, export the normalized counts if asked
    if (exportCountsTable) {
		geneDataExprExp <- as.data.frame(geneDataExpr)
		geneDataExprExp <- geneDataExprExp[,c(1:3,6,7,5,8,9)]
		colnames(geneDataExprExp)[1] <- "chromosome"
		geneDataFilteredExp <- as.data.frame(geneDataFiltered)
		geneDataFilteredExp <- geneDataFilteredExp[,c(1:3,6,7,5,8,9)]
		colnames(geneDataFilteredExp)[1] <- "chromosome"
        disp("Exporting and compressing normalized read counts table to ",
            file.path(PROJECT_PATH[["lists"]],"normalized_counts_table.txt"))
        expo <- cbind(
            rbind(geneDataExprExp,geneDataFilteredExp),
            rbind(normGenesExpr,geneCountsZero,geneCountsDead)
        )
        resFile <- file.path(PROJECT_PATH[["lists"]],
            "normalized_counts_table.txt.gz")
        gzfh <- gzfile(resFile,"w")
        write.table(expo,gzfh,sep="\t",row.names=FALSE,quote=FALSE)
        close(gzfh)
    }

    # Calculate meta-statistics, if more than one statistical algorithm has been used
    if (length(statistics)>1) {
        sumpList <- metaTest(
            cpList=cpList,
            metaP=metaP,
            counts=normGenesExpr,
            sampleList=sampleList,
            statistics=statistics,
            statArgs=statArgs,
            libsizeList=libsizeList,
            nperm=nperm,
            weight=weight,
            reprod=reprod,
            rc=restrictCores
        )
    }
    # We assign the p-values from the only statistic used to sumpList in order
    # to use it for stat plots
    else
        sumpList <- cpList
    # Useless for one statistics but just for safety
    if ("adj_meta_p_value" %in% exportWhat) 
        adjSumpList <- cmclapply(sumpList,
            function(x,a) return(p.adjust(x,a)),adjustMethod,rc=restrictCores)
    else
        adjSumpList <- NULL

    ############################################################################ 
    #                            EXPORT SECTION
    ############################################################################

    if (is.function(progressFun)) {
        text <- paste("Exporting...")
        progressFun(detail=text)
    }

    # Bind all the flags
    if (countType=="gene")
        flags <- geneFilterFlags
    else if (countType=="utr")
        flags <- geneFilterFlags
    else if (countType=="exon") {
        if (!is.null(exonFilterFlags)) {
            flags <- cbind(geneFilterFlags,
                as.matrix(exonFilterFlags[rownames(geneFilterFlags),]))
            nams <- c(colnames(geneFilterFlags),colnames(exonFilterFlags))
            rownames(flags) <- rownames(geneFilterFlags)
            colnames(flags) <- nams
        }
        else
            flags <- geneFilterFlags
    }
    
    disp("Building output files...")
    if (outList) 
        out <- makeExportList(contrast) 
    else 
        out <- NULL
    if (report) 
        html <- makeExportList(contrast) 
    else 
        html <- NULL
    if ("rpgm" %in% exportScale)
        fa <- attr(geneData,"geneLength")
    else
        fa <- NULL
    if ("normalized" %in% exportValues) {
        fac <- fa[rownames(normGenesExpr)]
        normList <- makeTransformation(normGenesExpr,exportScale,fac,logOffset)
    }
    else
        normList <- NULL
    if ("raw" %in% exportValues) {
        fac <- fa[rownames(geneCountsExpr)]
        rawList <- makeTransformation(geneCountsExpr,exportScale,fac,logOffset)
    }
    else
        rawList <- NULL
    if ("flags" %in% exportWhat)
        goodFlags <- flags[rownames(normGenesExpr),]
    else
        goodFlags <- NULL

    if (!is.null(geneCountsZero) || !is.null(geneCountsDead)) {
        geneCountsFiltered <- rbind(geneCountsZero,geneCountsDead)
        geneCountsUnnormFiltered <- rbind(geneCountsZero,
            geneCountsUnnorm)
        if ("normalized" %in% exportValues) {
            fac <- fa[rownames(geneCountsFiltered)]
            normListFiltered <- makeTransformation(geneCountsFiltered,
                exportScale,fac,logOffset)
        }
        else
            normListFiltered <- NULL
        if ("raw" %in% exportValues) {
            fac <- fa[rownames(geneCountsUnnormFiltered)]
            rawListFiltered <- makeTransformation(
                geneCountsUnnormFiltered,exportScale,fac,logOffset)
        }
        else
            rawListFiltered <- NULL
        if ("flags" %in% exportWhat && !is.null(flags))
            allFlags <- rbind(
                matrix(1,nrow(geneCountsZero),ncol(flags)),
                    flags[rownames(geneCountsDead),]
            )
        else
            allFlags <- NULL
    }
    else {
        geneCountsFiltered <- NULL
        geneCountsUnnormFiltered <- NULL
        allFlags <- NULL
    }
    
    reportTables <- vector("list",length(contrast))
	names(reportTables) <- contrast
    
    counter <- 1
    for (cnt in contrast) {
        disp("  Contrast: ",cnt)
        disp("    Adding non-filtered data...")
        theExport <- buildExport(
            geneData=geneDataExpr,
            rawGeneCounts=geneCountsExpr,
            normGeneCounts=normGenesExpr,
            flags=goodFlags,
            sampleList=sampleList,
            cnt=cnt,
            statistics=statistics,
            rawList=rawList,
            normList=normList,
            pMat=cpList[[cnt]],
            adjpMat=adjCpList[[cnt]],
            sumP=sumpList[[cnt]],
            adjSumP=adjSumpList[[cnt]],
            exportWhat=exportWhat,
            exportScale=exportScale,
            exportValues=exportValues,
            exportStats=exportStats,
            logOffset=logOffset,
            #report=report
            report=FALSE
        )
        
        # If report requested, build a more condensed summary table, while the
        # complete tables are available for download
        if (report) { 
			if (length(statistics) > 1)
				ew <- c("annotation","meta_p_value","adj_meta_p_value",
					"fold_change","stats")
			else
				ew <- c("annotation","p_value","adj_p_value","fold_change",
					"stats")
			esc <- "rpgm"
			ev <- "normalized"
			est <- "mean"
			
			faR <- attr(geneDataExpr,"geneLength")			
			facR <- faR[rownames(normGenesExpr)]
			normListR <- makeTransformation(normGenesExpr,esc,facR,logOffset)
			
			disp("    Adding report data...")
			reportTables[[cnt]] <- buildExport(
				geneData=geneDataExpr,
				rawGeneCounts=geneCountsExpr,
				normGeneCounts=normGenesExpr,
				flags=goodFlags,
				sampleList=sampleList,
				cnt=cnt,
				statistics=statistics,
				rawList=NULL,
				normList=normListR,
				pMat=cpList[[cnt]],
				adjpMat=adjCpList[[cnt]],
				sumP=sumpList[[cnt]],
				adjSumP=adjSumpList[[cnt]],
				exportWhat=ew,
				exportScale=esc,
				exportValues=ev,
				exportStats=est,
				logOffset=logOffset,
				report=FALSE
			)$textTable
			
			if (!is.null(reportTop)) {
				topi <- ceiling(reportTop*nrow(reportTables[[cnt]]))
				reportTables[[cnt]] <- reportTables[[cnt]][1:topi,,drop=FALSE]
			}
		}

        # Adjust the export based on what statistics have been done and a 
        # possible p-value cutoff
        export <- theExport$textTable
        colnames(export)[1] <- "chromosome"
        if (report)
            exportHtml <- theExport$htmlTable
        if (!is.na(pcut)) {
            if (length(statistics)>1) {
                switch(metaP,
                    fisher = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    fperm = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    whitlock = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    dperm_min = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    dperm_max = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    dperm_weight = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    minp = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    maxp = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    weight = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    pandora = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    simes = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    },
                    none = {
                        cutInd <- which(sumpList[[cnt]]<=pcut)
                    }
                )
                pp <- sumpList[[cnt]][cutInd]
                export <- export[cutInd,]
                export <- export[order(pp),]
                if (report) {
                    exportHtml <- exportHtml[cutInd,]
                    exportHtml <- exportHtml[order(pp),]
                }
            }
            else {
                cutInd <- which(sumpList[[cnt]]<=pcut)
                pp <- sumpList[[cnt]][cutInd,]
                export <- export[cutInd,]
                export <- export[order(pp),]
                if (report) {
                    exportHtml <- exportHtml[cutInd,]
                    exportHtml <- exportHtml[order(pp),]
                }
            }
        }
        else {
            pp <- sumpList[[cnt]]
            export <- export[order(pp),]
            if (report)
                exportHtml <- exportHtml[order(pp),]
        }

        # Final safety trigger
        naInd <- grep("NA",rownames(export))
        if (length(naInd)>0) {
            export <- export[-naInd,]
            if (report) 
                exportHtml <- exportHtml[-naInd,]
        }

        resFile <- file.path(PROJECT_PATH[["lists"]],
            paste("metaseqr_sig_out_",cnt,".txt.gz",sep=""))
        disp("    Writing output...")
        gzfh <- gzfile(resFile,"w")
        write.table(export,gzfh,quote=FALSE,row.names=FALSE,sep="\t")
        close(gzfh)
        if (outList)
            out[[cnt]] <- export

        if (!is.null(geneCountsZero) || !is.null(geneCountsDead)) {
            disp("    Adding filtered data...")
            theExportFiltered <- buildExport(
                geneData=geneDataFiltered,
                rawGeneCounts=geneCountsUnnormFiltered,
                normGeneCounts=geneCountsFiltered,
                flags=allFlags,
                sampleList=sampleList,
                cnt=cnt,
                statistics=statistics,
                rawList=rawListFiltered,
                normList=normListFiltered,
                exportWhat=exportWhat,
                exportScale=exportScale,
                exportValues=exportValues,
                exportStats=exportStats,
                logOffset=logOffset,
                report=FALSE
            )

            # Now we should be having theExport and theExportFiltered. We do not 
            # generate html output for filtered or total results just a 
            # compressed text file. We thus have to append theExport$textTable
            # and theExportFiltered$htmlTable before writing the final output...
            exportAll <- rbind(theExport$textTable,theExportFiltered$textTable)
            # ...and order them somehow... alphabetically according to row
            # names, as the annotation might not have been bundled...
            exportAll <- exportAll[order(rownames(exportAll)),]
            
            # Here, both filtered and unfiltered genes are passed to the output
            # list.
            if (outList)
                out[[cnt]] <- exportAll   
            
            resFile <- file.path(PROJECT_PATH[["lists"]],paste(
                "metaseqr_all_out_",cnt,".txt.gz",sep=""))
            disp("    Writing output...")
            gzfh <- gzfile(resFile,"w")
            write.table(exportAll,gzfh,quote=FALSE,row.names=FALSE,sep="\t")
            close(gzfh)
        }
    }

    ############################################################################
    # END EXPORT SECTION
    ############################################################################

    ############################################################################
    # BEGIN PLOTTING SECTION
    ############################################################################
    
    if (is.function(progressFun)) {
        text <- paste("Plotting...")
        progressFun(detail=text)
    }
    
    # Check if we have more than 6 samples in total, pairwise plots are not
    # meaningful
    if (length(unlist(sampleList)) > 6 && "pairwise" %in% qcPlots) {
		warnwrap("Pairwise sample comparison plot becomes indistinguishable ",
			"for more than 6 samples! Removing from plots...")
		qcPlots <- qcPlots[-which(qcPlots == "pairwise")]
	}
	
    if (!is.null(qcPlots)) {
        disp("Creating quality control graphs...")
        plots <- list(
            raw=c("mds","biodetection","countsbio","saturation","readnoise",
                "correl","pairwise"),
            norm=c("boxplot","gcbias","lengthbias","meandiff","meanvar",
                "rnacomp"),
            stat=c("deheatmap","volcano","mastat","biodist"),
            other=c("filtered"),
            venn=c("venn")
        )
        figRaw <- figUnorm <- figNorm <- figStat <- figOther <- figVenn <- 
            vector("list",length(figFormat))
        names(figRaw) <- names(figUnorm) <- names(figNorm) <-
            names(figStat) <- names(figOther) <- names(figVenn) <-
            figFormat
        for (fig in figFormat) {
            disp("Plotting in ",fig," format...")
            figRaw[[fig]] <- metaseqrPlot(geneCounts,sampleList,
                annotation=geneData,plotType=intersect(qcPlots,plots$raw),
                isNorm=FALSE,output=fig,path=PROJECT_PATH$qc)
            
            figUnorm[[fig]] <- metaseqrPlot(geneCounts,sampleList,
                annotation=geneData,plotType=intersect(qcPlots,plots$norm),
                isNorm=FALSE,output=fig,path=PROJECT_PATH$normalization)
            
            if (whenApplyFilter=="prenorm") # The annotation dimensions change...
                figNorm[[fig]] <- metaseqrPlot(normGenes,sampleList,
                    annotation=geneDataExpr,plotType=intersect(qcPlots,
                    plots$norm),isNorm=TRUE,output=fig,
                    path=PROJECT_PATH$normalization) 
            else if (whenApplyFilter=="postnorm")
                figNorm[[fig]] <- metaseqrPlot(normGenes,sampleList,
                    annotation=geneData,plotType=intersect(qcPlots,
                    plots$norm),isNorm=TRUE,output=fig,
                    path=PROJECT_PATH$normalization)
            
            figStat[[fig]] <- metaseqrPlot(normGenesExpr,sampleList,
                annotation=geneDataExpr,contrastList=contrastList,
                pList=sumpList,thresholds=list(p=pcut,f=1),
                plotType=intersect(qcPlots,plots$stat),isNorm=TRUE,
                output=fig,path=PROJECT_PATH$statistics)
            if (!is.null(geneDataFiltered))
                figOther[[fig]] <- metaseqrPlot(geneDataFiltered,
                    sampleList,annotation=totalGeneData,
                    plotType=intersect(qcPlots,plots$other),isNorm=FALSE,
                    output=fig,path=PROJECT_PATH$qc)
            else 
                figOther[[fig]] <- NULL
            
            if ("venn" %in% qcPlots)
                figVenn[[fig]] <- metaseqrPlot(normGenesExpr,
                    sampleList,annotation=geneDataExpr,
                    contrastList=contrastList,pList=cpList,
                    thresholds=list(p=pcut,f=1),
                    plotType=intersect(qcPlots,plots$venn),
                    output=fig,path=PROJECT_PATH$statistics)
        }
        
        ########################################################################
        #assign("sampleList",sampleList,envir=parent.frame())
		#assign("geneCounts",geneCounts,envir=parent.frame())
		#assign("normGenes",normGenes,envir=parent.frame())
		assign("normGenesExpr",normGenes,envir=parent.frame())
		assign("sumpList",sumpList,envir=parent.frame())
		assign("contrastList",contrastList,envir=parent.frame())
		#assign("geneData",geneData,envir=parent.frame())
		########################################################################
    }

    ############################################################################
    # END PLOTTING SECTION
    ############################################################################

    ############################################################################
    # BEGIN REPORTING SECTION
    ############################################################################
    
    if (report) {
		# Help data
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
		
		if (reportDb == "sqlite") {
			samples <- unlist(sampleList)
			nsa <- length(samples)
			
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
				json <- diagplotNoiseq(geneCounts,sampleList,covars=covarsRaw,
					whichPlot="biodetection",output="json")
				for (s in samples) {
					disp("    ",s)
					name <- paste("biodetection",s,sep="_")
					.dbImportPlot(con,name,"biodetection","generic",json[[s]])
				}
			}
			if ("countsbio" %in% qcPlots) {
				disp("  Importing countsbio...")
				jsonList <- diagplotNoiseq(geneCounts,sampleList,
					covars=covarsRaw,whichPlot="countsbio",output="json")
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
				.dbImportPlot(con,"ReadNoise","readnoise",NULL,json)
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
						#nam <- paste(n,s,sep="_")
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
				jsonUnorm <- diagplotNoiseq(geneCounts,sampleList,
					covars=covarsRaw,whichPlot="rnacomp",isNorm=FALSE,
					output="json")
				jsonNorm <- diagplotNoiseq(normGenes,sampleList,
					covars=covarsRaw,whichPlot="rnacomp",isNorm=TRUE,
					output="json")
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
					json <- diagplotNoiseq(normGenesExpr,sampleList,
						covars=covarsStat,whichPlot="biodist",
						biodistOpts=list(p=cpList[[cnt]],pcut=pcut,name=cnt),
						output="json")
					.dbImportPlot(con,paste("biodist",n,sep="_"),"biodist",
						"chromosome",json$chromosome)
					.dbImportPlot(con,paste("biodist",n,sep="_"),"biodist",
						"biotype",json$biotype)
				}
			}
			if ("venn" %in% qcPlots) {
				disp("  Importing venn")
				nn <- names(contrastList)
				geneNames <- as.character(geneDataExpr$gene_name)
				names(geneNames) <- names(geneDataExpr)
				for (n in nn) {
					disp("    ",n)
					json <- makeJVennData(cpList[[n]],pcut=pcut,
						altNames=geneNames[rownames(cpList[[n]])])
					.dbImportPlot(con,paste("venn",n,sep="_"),"venn","generic",
						toJSON(json,auto_unbox=TRUE,null="null"))
				}
			}
			
			# Close SQLite connection
			dbDisconnect(con)
		}
		
		disp("Creating HTML report...")
		
		if (!is.null(qcPlots)) {
            # First create zip archives of the figures
            disp("Compressing figures...")
            zipfiles <- file.path(PROJECT_PATH$plots,paste("metaseqr_figures_",
                figFormat,".zip",sep=""))
            names(zipfiles) <- figFormat
            for (f in figFormat) {
                files <- c(
                    dir(PROJECT_PATH$qc,pattern=paste(".",f,sep=""),
                        full.names=TRUE),
                    dir(PROJECT_PATH$normalization,pattern=paste(".",f,sep=""),
                        full.names=TRUE),
                    dir(PROJECT_PATH$statistics,pattern=paste(".",f,sep=""),
                        full.names=TRUE)
                )
                zip(zipfiles[f],files)
            }
            # Then create the final figure variables which brew will find...
            figRaw <- figRaw[["png"]]
            figUnorm <- figUnorm[["png"]]
            figNorm <- figNorm[["png"]]
            figStat <- figStat[["png"]]
            figOther <- figOther[["png"]]
            figVenn <- figVenn[["png"]]
        }
        else
            figRaw <- figUnorm <- figNorm <- figStat <- figOther <-
                figVenn <- NULL
		
		## Then see what is going of if default report changed
        #if (tolower(reportTemplate)=="default") {
        #    if (exists("TEMPLATE")) {
        #        reportTemplate=list(
        #            rmd=file.path(TEMPLATE,"metaseqr2_report.Rmd"),
        #            loader=file.path(TEMPLATE,"dna_loader.gif")
        #        )
        #    }
        #    else
        #        reportTemplate=list(rmd=NULL,loader=NULL)
        #}
        #if (!is.null(reportTemplate$rmd)) {
        #    if (file.exists(reportTemplate$rmd)) {
        #        template <- reportTemplate$rmd
        #        hasTemplate <- TRUE
        #    }
        #    else {
        #        warnwrap(paste("The template file",reportTemplate$rmd,
        #            "was not ","found! The HTML report will NOT be generated."))
        #        hasTemplate <- FALSE
        #    }
        #}
        #else {
        #    warnwrap(paste("The report option was enabled but no template ",
        #        "file is provided! The HTML report will NOT be generated."))
        #    hasTemplate <- FALSE
        #}
        #if (!is.null(reportTemplate$loader)) {
        #    if (file.exists(reportTemplate$loader))
        #        file.copy(from=reportTemplate$loader,to=PROJECT_PATH$media)
        #    else
        #        warnwrap(paste("The report logo image",reportTemplate$loader,
        #            "was not found!"))
        #}
        #else
        #    warnwrap(paste("The report loader image was not provided!"))
		
		# Here we must download all required libraries and put them in the js
		# folder of the report to make it available offline
		if (offlineReport) {
			disp("Downloading required JavaScript libraries...")
			if (!file.exists(file.path(PROJECT_PATH$js,"pace.js")))
				download.file(paste0("https://raw.github.com/HubSpot/pace/",
					"v1.0.0/pace.min.js"),
					file.path(PROJECT_PATH$js,"pace.min.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"highcharts.js")))
				download.file("https://code.highcharts.com/highcharts.js",
					file.path(PROJECT_PATH$js,"highcharts.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"highcharts-more.js")))
				download.file("https://code.highcharts.com/highcharts-more.js",
					file.path(PROJECT_PATH$js,"highcharts-more.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"exporting.js")))
				download.file("https://code.highcharts.com/modules/exporting.js",
					file.path(PROJECT_PATH$js,"exporting.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"offline-exporting.js")))
				download.file(
					"https://code.highcharts.com/modules/offline-exporting.js",
					file.path(PROJECT_PATH$js,"offline-exporting.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"export-data.js")))
				download.file(
					"https://code.highcharts.com/modules/export-data.js",
					file.path(PROJECT_PATH$js,"export-data.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"canvas2svg.js")))
				download.file(
					"http://jvenn.toulouse.inra.fr/app/js/canvas2svg.js",
					file.path(PROJECT_PATH$js,"canvas2svg.js"))
			if (!file.exists(file.path(PROJECT_PATH$js,"jvenn.min.js")))
				download.file(
					"http://jvenn.toulouse.inra.fr/app/js/jvenn.min.js",
					file.path(PROJECT_PATH$js,"jvenn.min.js"))
			if (reportDb == "sqlite") {
				if (!file.exists(file.path(PROJECT_PATH$js,"sql.js")))
					download.file(
				"https://cdnjs.cloudflare.com/ajax/libs/sql.js/0.5.0/js/sql.js",
				file.path(PROJECT_PATH$js,"sql.js"))
			}
			else if (reportDb == "dexie") {
				if (!file.exists(file.path(PROJECT_PATH$js,"dexie.min.js")))
					download.file(
						"https://unpkg.com/dexie@2.0.4/dist/dexie.min.js",
						file.path(PROJECT_PATH$js,"dexie.min.js"))
			}
		}
		
		#if (hasTemplate) {
			execTime <- elap2human(TB)
            TEMP <- environment()
            
            REPORT_ENV <- .makeReportEnv(environment())
                        
            #file.copy(file.path(TEMPLATE,"metaseqr2_report.Rmd"),
            file.copy("/media/raid/software/metaseqR2-local/inst/metaseqr2_report.Rmd",
			#file.copy("C:/software/metaseqR2-local/inst/metaseqr2_report.Rmd",
				file.path(PROJECT_PATH$main,"metaseqr2_report.Rmd"),
				overwrite=TRUE)
			invisible(knitr::knit_meta(class=NULL,clean=TRUE))
			render(
				input=file.path(PROJECT_PATH$main,"metaseqr2_report.Rmd"),
				output_file="index.html",
				output_dir=PROJECT_PATH$main,
				#output_format="html_document",
				#envir=new.env(parent=globalenv()),
				envir=TEMP#,
				#encoding="UTF-8"
			)
			gc(verbose=FALSE)
			# Remove the Rmd file after rendering the report
			unlink(file.path(PROJECT_PATH$main,"metaseqr2_report.Rmd"))
			
			####################################################################
			# Replace a jQuery incompatibility bug line.. Hopefully to remove...
			L <- readLines(file.path(PROJECT_PATH$main,"index.html"))
			ln <- grep(paste("$(\".menu\").find(\"li[data-target=\" +",
				"window.page + \"]\").trigger(\"click\");"),L,fixed=TRUE)
			if (length(ln) == 1)
				L[ln] <- paste("$(\".menu\").find(\"li[data-target='\" +",
				"window.page + \"']\").trigger(\"click\");")
			cat(L,file=file.path(PROJECT_PATH$main,"index.html"),sep="\n")
			####################################################################
        #}
    }

    ############################################################################
    # END REPORTING SECTION
    ############################################################################
    
    disp("\n",strftime(Sys.time()),": Data processing finished!\n")
    execTime <- elap2human(TB)
    disp("\n","Total processing time: ",execTime,"\n\n")
                             
    if (outList) {
        tmp <- c(geneDataExpr,geneDataFiltered)
        a <- c(attr(geneDataExpr,"geneLength"),
            attr(geneDataFiltered,"geneLength"))
        names(a) <- rownames(tmp)
        attr(tmp,"geneLength") <- a
        for (n in names(cpList)) {
            if (!is.null(geneDataFiltered)) {
                filler <- matrix(NA,length(geneDataFiltered),ncol(cpList[[n]]))
                rownames(filler) <- names(geneDataFiltered)
                colnames(filler) <- colnames(cpList[[n]])
            }
            else
                filler <- NULL
            cpList[[n]] <- rbind(cpList[[n]],filler)
            cpList[[n]] <- cpList[[n]][rownames(tmp),,drop=FALSE]
           }
        if (!is.null(adjCpList)) {
            for (n in names(adjCpList)) {
                if (!is.null(geneDataFiltered)) {
                    filler <- matrix(NA,length(geneDataFiltered),
                        ncol(adjCpList[[n]]))
                    rownames(filler) <- names(geneDataFiltered)
                    colnames(filler) <- colnames(cpList[[n]])
                }
                else
                    filler <- NULL
                adjCpList[[n]] <- rbind(adjCpList[[n]],filler)
                adjCpList[[n]] <- adjCpList[[n]][rownames(tmp),,drop=FALSE]
            }
        }
        if (!is.null(sumpList)) {
            for (n in names(sumpList)) {
                if (!is.null(geneDataFiltered)) {
                    filler <- rep(NA,length(geneDataFiltered))
                    names(filler) <- names(geneDataFiltered)
                }
                else
                    filler <- NULL
                if (is.matrix(sumpList[[n]])) {
                    # The following if-else was added because some runs had 
                    # filler=NULL while others do not. Sometimes this created a 
                    # problem and the simulations stopped. Now fixed.
                    if (is.null(filler)) {
                        sumpList[[n]] <- rbind(sumpList[[n]],filler)
                        sumpList[[n]] <- sumpList[[n]][rownames(tmp),,
                            drop=FALSE]
                    }
                    else {
                        filler=as.matrix(filler)
                        sumpList[[n]] <- rbind(sumpList[[n]],filler)
                        sumpList[[n]] <- sumpList[[n]][rownames(tmp),,
                            drop=FALSE]
                    }
                }
                else {
                    sumpList[[n]] <- c(sumpList[[n]],filler)
                    sumpList[[n]] <- sumpList[[n]][rownames(tmp)]
                }
           }
        }
        if (!is.null(adjSumpList)) {
           for (n in names(adjSumpList)) {
               if (!is.null(geneDataFiltered)) {
                   filler <- rep(NA,length(geneDataFiltered))
                   names(filler) <- names(geneDataFiltered)
               }
               else
                   filler <- NULL
               adjSumpList[[n]] <- c(adjSumpList[[n]],filler)
               adjSumpList[[n]] <- adjSumpList[[n]][rownames(tmp)]
           }
        }

        complete <- list(
            #call=as.list(match.call()),
            params=list(
                sampleList=sampleList,
                excludeList=excludeList,
                fileType=fileType,
                path=path,
                contrast=contrast,
                libsizeList=libsizeList,
                embedCols=embedCols,
                annotation=annotation,
                org=org,
                refdb=refdb,
                countType=countType,
                exonFilters=exonFilters,
                geneFilters=geneFilters,
                whenApplyFilter=whenApplyFilter,
                normalization=normalization,
                normArgs=normArgs,
                statistics=statistics,
                statArgs=statArgs,
                adjustMethod=adjustMethod,
                metaP=metaP,
                weight=weight,
                nperm=nperm,
                reprod=reprod,
                pcut=pcut,
                logOffset=logOffset,
                preset=preset,
                qcPlots=qcPlots,
                figFormat=figFormat,
                outList=outList,
                exportWhere=exportWhere,
                exportWhat=exportWhat,
                exportScale=exportScale,
                exportValues=exportValues,
                exportStats=exportStats,
                exportCountsTable=exportCountsTable,
                restrictCores=restrictCores,
                report=report,
                reportTop=reportTop,
                reportTemplate=reportTemplate,
                saveGeneModel=saveGeneModel,
                verbose=verbose,
                runLog=runLog
            ),
            filterCutoffs=list(
                exonFilter=list(
                    minActiveExons=NULL
                ),
                geneFilter=list(
                    length=geneFilters$length$length,
                    avgReads=if (is.null(geneFilterCutoff)) NULL else
                        round(geneFilterCutoff$avgReads,digits=5),
                    expression=list(
                        median=geneFilterCutoff$expression$median,
                        mean=geneFilterCutoff$expression$mean,
                        quantile=geneFilterCutoff$expression$quantile,
                        known=geneFilterCutoff$expression$known,
                        custom=geneFilterCutoff$expression$custom
                    ),
                    biotype=if (is.null(geneFilters$biotype)) NULL else
                        paste(names(geneFilters$biotype)[which(
                            unlist(geneFilters$biotype))],collapse=", ")
                ),
                zeroFiltered=length(theZeros),
                exonFiltered=length(unique(unlist(exonFilterResult))),
                geneFiltered=length(unique(unlist(geneFilterResult))),
                totalFiltered=length(theZeros)+length(theDead)
            ),
            geneData=tmp,
            rawCounts=rbind(geneCountsExpr,geneCountsUnnormFiltered),
            normCounts=rbind(normGenesExpr,geneCountsFiltered),
            flags=rbind(goodFlags,allFlags),
            sampleList=sampleList,
            contrastList=contrast,
            pValue=cpList,
            fdr=adjCpList,
            metaPValue=sumpList,
            metaFdr=adjSumpList
        )

        return(list(data=out,html=html,complete=complete))
    }
    
} # End metaseqr2


constructGeneModel <- function(countData,annoData,type,rc=NULL) {
	# countData is a matrix, must have exact rownames as annoData GRanges
	if (!all(rownames(countData) == names(annoData)))
		stop("The rownames of the provided counts are not the same (or in ",
			" the same order) as the provided annotation GRanges!")
	
	# Messaging...
	msg <- "provided count regions"
	if (!missing(type)) {
		if (type == "exon")
			msg <- "exons"
		else if (type == "utr")
			msg <- "transcripts (UTR regions)"
	}
	
	# Create splitting factor
	theGenes <- factor(annoData$gene_id,levels=unique(annoData$gene_id))
	
	# Split the counts vectors...
	theCounts <- cmclapply(colnames(countData),function(n,f,M) {
		disp("  Separating ",msg," per gene for ",n,"...")
        co <- M[,n]
        names(co) <- rownames(M)
        return(split(co,f))
	},theGenes,countData,rc=rc)
	names(theCounts) <- colnames(countData)
	
	# and the width of the merged exons
	w <- width(annoData)
	if (type == "exon")
		names(w) <- annoData$exon_id
	else if (type == "utr")
		names(w) <- annoData$transcript_id
	lengthList <- split(w,theGenes)
	attr(theCounts,"lengthList") <- lengthList
	# The exon lengths information is essentially redundant in the new version,
	# unless we need it for some QC distribution later on...
	
    return(theCounts)
}

.reduceGeneData <- function(exonData,geneData) {
    exonChrs <- unique(as.character(seqnames(exonData)))
    geneChrs <- unique(as.character(seqnames(geneData)))
    if (length(exonChrs) != length(geneChrs)) {
        m <- match(seqnames(geneData),exonChrs)
        geneData <- geneData[which(!is.na(m))]
    }
    return(geneData)
}

.userOrg <- function(org,db=NULL) {
	ua <- getUserAnnotations(db)
	if (nrow(ua) == 0)
		return(FALSE)
	orgs <- unique(ua$organism)
	return(org %in% orgs)
}

.userRefdb <- function(refdb,db=NULL) {
	us <- getUserAnnotations(db)
	if (nrow(us) == 0)
		return(FALSE)
	refdbs <- unique(us$source)
	return(refdb %in% refdbs)
}

.backwardsCompatibility <- function(dataFile) {
    tmpEnv <- new.env()
    load(dataFile,tmpEnv)
    if (!is.null(tmpEnv$the.counts)) { # Old file
        # Copy to new
        tmpEnv$theCounts <- tmpEnv$the.counts
        tmpEnv$countsType <- tmpEnv$count.type
        tmpEnv$sampleList <- tmpEnv$sampleList
        tmpEnv$geneData <- tmpEnv$gene.data
        if (!is.null(tmpEnv$transcript.data))
            tmpEnv$transcriptData <- tmpEnv$transcript.data
        if (!is.null(tmpEnv$exon.data))
            tmpEnv$exonData <- tmpEnv$exon.data
        
        # Delete old
        if (exists("the.counts",envir=tmpEnv)) {
			rm(the.counts,count.type,sample.list,gene.data,envir=tmpEnv)
			if (!is.null(tmpEnv$transcript.data))
				rm(transcript.data,envir=tmpEnv)
			if (!is.null(tmpEnv$exon.data))
				rm(exon.data,envir=tmpEnv)
		}        
    }
    return(tmpEnv)
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
			x$gene_id <- paste('<a href="',h,'/Gene/Summary?g=',x$gene_id,
				'" target="_blank">',x$gene_id,'</a>',sep="")
		}
		else if (r == "refseq") {
			x$gene_id <- paste('<a href="https://www.ncbi.nlm.nih.gov/nuccore/',
				x$gene_id,'" target="_blank">',x$gene_id,'</a>',sep="")
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

#############################################################################

# SQL(ite) backed-up ready (for later shiny app)

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

.makeReportEnv <- function(e) {
	re <- new.env(parent=globalenv())
    re$offlineReport <- e$offlineReport
    re$reportDb <- e$reportDb
    re$qcPlots <- e$qcPlots
    re$sampleList <- e$sampleList
    re$geneData <- e$geneData
    re$geneCounts <- e$geneCounts
    re$geneData <- e$geneDataExpr
    re$geneCounts <- e$normGenesExpr
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
    re$contrastList <- e$contrastList
    re$reportTop <- e$reportTop
    re$cpList <- e$cpList
    re$reportTables <- e$reportTables
    re$exportCountsTable <- e$exportCountsTable
    
    re$FUN_CALL <- e$FUN_CALL
    re$PROJECT_PATH <- e$PROJECT_PATH
    
    return(re)
}
