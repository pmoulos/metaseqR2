checkMainArgs <- function(mainArgs) {
    inArgs <- names(mainArgs)[-1] # 1st member name of calling function
    validArgs <- .getValidArgs()
    invalid <- setdiff(inArgs,validArgs)
    if (length(invalid) > 0) {
        for (i in 1:length(invalid))
            warnwrap("Unknown input argument to metaseqr pipeline: ",invalid[i],
                " ...Ignoring...",now=TRUE)
    }
}

checkTextArgs <- function(argName,argValue,argList,multiarg=FALSE) {
    if (multiarg) {
        argValue <- tolower(argValue)
        if (!all(argValue %in% argList))
            stopwrap("\"",argName,"\""," parameter must be one or more of ",
                paste(paste("\"",argList,sep=""),collapse="\", "),"\"!")
    }
    else {
        argSave <- argValue[1]
        argValue <- tolower(argValue[1])
        # An exception must be added for annotation because it can be an  
        # external file too
        #if (argName=="annotation") { 
        #    if (!(argValue %in% argList) && !file.exists(argSave))
        #        stopwrap("\"",argName,"\""," parameter must be one of ",
        #            paste(paste("\"",argList,sep=""),collapse="\", "),
        #                "\" or an existing file!")
        #}
        #else {
            if (!(argValue %in% argList))
                stopwrap("\"",argName,"\""," parameter must be one of ",
                    paste(paste("\"",argList,sep=""),collapse="\", "),"\"!")
        #}
    }
}

checkNumArgs <- function(argName,argValue,argType,argBounds,direction) {
    switch(argType,
        numeric = {
            if (!is.numeric(argValue))
                stopwrap("\"",argName,"\"",
                    " parameter must be a numeric value!")
            if (!missing(argBounds)) {
                switch(direction,
                    both = {
                        if (argValue<=argBounds[1] ||
                            argValue>=argBounds[2])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric ","value larger than or equal to ",
                                argBounds[1]," and smaller than or equal to ",
                                argBounds[2],"!")
                    },
                    botheq = {
                        if (argValue<argBounds[1] || argValue>argBounds[2])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric value larger than ",argBounds[1],
                                " and smaller than ",argBounds[2],"!")
                    },
                    gt = {
                        if (argValue<=argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric value greater than ",argBounds[1],"!")
                    },
                    lt = {
                        if (argValue>=argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric value lower than ",argBounds[1],"!")
                    },
                    gte = {
                        if (argValue<argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric value greater than or equal to ",
                                argBounds[1],"!")
                    },
                    lte = {
                        if (argValue>argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be a ",
                                "numeric value lower than or equal to ",
                                argBounds[1],"!")
                    }
                )
            }
        },
        integer = {
            if (!is.integer(argValue))
                stopwrap("\"",argName,"\""," parameter must be an integer!")
            if (!missing(argBounds)) {
                switch(direction,
                    both = {
                        if (argValue<=argBounds[1] ||
                            argValue>=argBounds[2])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer larger than or equal to ",
                                argBounds[1]," and smaller than or equal to ",
                                argBounds[2],"!")
                    },
                    botheq = {
                        if (argValue<argBounds[1] || argValue>argBounds[2])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer larger than or equal to ",
                                argBounds[1]," and smaller than or equal to ",
                                argBounds[2],"!")
                    },
                    gt = {
                        if (argValue<=argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer greater than ",argBounds[1],"!")
                    },
                    lt = {
                        if (argValue>=argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer lower than ",argBounds[1],"!")
                    },
                    gte = {
                        if (argValue<argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer greater than or equal to ",
                                argBounds[1],"!")
                    },
                    lte = {
                        if (argValue>argBounds[1])
                            stopwrap("\"",argName,"\""," parameter must be ",
                                "an integer lower than or equal to ",
                                argBounds[1],"!")
                    }
                )
            }
        }
    )
}

checkFileArgs <- function(argName,argValue) {
    if (!file.exists(argValue))
        stopwrap("\"",argName,"\""," parameter must be an existing file!")
}

checkContrastFormat <- function(cnt,sampleList) {
    # This function will break cnt and check that all contrast counter parts are 
    # members of the names of the sampleList and contain the string "_vs_" as 
    # many times as the names of the sampleList minus 1. If satisfied return 
    # TRUE else error.
    cnts <- strsplit(cnt,"_vs_")
    #if (length(unique(unlist(cnts))) != length(names(sampleList)))
    if (!all(unique(unlist(cnts)) %in% names(sampleList)))
        stopwrap("Condition names in sample list and contrast list do not ",
            "match! Check if the contrasts follow the appropriate format (e.g.",
            " \"_vs_\" separating contrasting conditions...")
    if (length(unique(cnt))!=length(cnt))
        warnwrap("Duplicates found in the contrasts list! Duplicates will be ",
            "ignored...")
}

checkLibsize <- function(libsizeList,sampleList) {
    if (!is.null(libsizeList)) {
        if (length(intersect(names(libsizeList),unlist(sampleList,
            use.names=FALSE)))!=length(unlist(sampleList,
            use.names=FALSE))) {
            warnwrap("Sample names in \"libsizeList\" and \"sampleList\" do ",
                "not match! Library sizes will be estimated from count data...")
            return(NULL)
        }
        else return(libsizeList)
    }
    else
        return(NULL)
}

checkPackages <- function(m,p) {
    # Check meta-analysis packages
    if (m=="whitlock" && !require(survcomp))
        stopwrap("Bioconductor package survcomp is required for \"whitlock\" ",
            "p-value meta analysis!")
    if ("venn" %in% p && !require(VennDiagram))
        stopwrap("R package VennDiagram is required for some of the selected ",
            "QC plots!")
}

.backwardsConvertArgs <- function(args) {
    args <- args[-1] # 1st member name of calling function
    inArgs <- names(args)
    oldArgs <- .getBackwardsValidArgs()
    newArgs <- .getValidArgs()
    old <- intersect(inArgs,oldArgs)
    
    # Detect mixed arguments by excluding the common and then checking the rest
    # A mix can be detected if we exclude the common ones (e.g. 'normalization')
    # from both newArgs and oldArgs and intersect the rest
    common <- intersect(newArgs,oldArgs)
    restNew <- setdiff(newArgs,common)
    restOld <- setdiff(oldArgs,common)
    if (any(inArgs %in% restNew) && any(inArgs %in% restOld))
        return(list(backDetected=FALSE,mixDetected=TRUE,args=NULL))
    
    # Continue, as the mix has been covered. The call with old arguments 
    # contains also the common
    oldInCall <- setdiff(intersect(inArgs,oldArgs),common)
    if (length(oldInCall) > 0) {
        backArgs <- list()
        backArgs$backDetected <- TRUE
        backArgs$mixDetected <- FALSE
        backArgs$args <- args[oldInCall]
        
        # If old call contains id.col, gc.col, name.col or bt.col, these must be
        # grouped. We define anyway and fill them as needed.
        backArgs$args$embedCols <- list(
			idCol=NA,
			gcCol=NA,
			nameCol=NA,
			btCol=NA
        )
        if ("id.col" %in% oldInCall) {
			backArgs$args$embedCols$idCol <- args$id.col
			backArgs$args$id.col <- NULL
		}
		if ("gc.col" %in% oldInCall) {
			backArgs$args$embedCols$gcCol <- args$gc.col
			backArgs$args$gc.col <- NULL
		}
		if ("name.col" %in% oldInCall) {
			backArgs$args$embedCols$nameCol <- args$name.col
			backArgs$args$name.col <- NULL
		}
		if ("bt.col" %in% oldInCall) {
			backArgs$args$embedCols$btCol <- args$bt.col
			backArgs$args$bt.col <- NULL
		}
        
        return(backArgs)
    }
    else
        return(list(backDetected=FALSE,mixDetected=FALSE,args=NULL))
}

.backwardsMapOld2New <- function() {
    oldArgNames <- .getBackwardsValidArgs()
    newArgNames <- .getValidArgs()
    
    # The main differences are: i) there is no gene.file in the new args,
    # ii) the idCol, gcCol, nameCol and btCol arguments are within embedCols
    map <- list()
    
    # 1. Remove "id.col","gc.col","name.col","bt.col" from old arguments so as
    # to group them. Remove also gene.file.
    r1 <- match(c("id.col","gc.col","name.col","bt.col","gene.file"),
		oldArgNames)
    oldArgNames <- oldArgNames[-r1]
    # 2. Remove embedCols from new arguments as they can't be mapped directly
    r2 <- match("embedCols",newArgNames)
    newArgNames <- newArgNames[-r2]
    
    # We can now map
    map[oldArgNames] <- newArgNames
    return(map)
}

.getValidArgs <- function() {
    return(c(
        "counts","sampleList","excludeList","fileType","path","contrast",
        "libsizeList","embedCols","annotation","org","transLevel","countType",
        "utrFlank","exonFilters","geneFilters","whenApplyFilter",
        "normalization","normArgs","statistics","statArgs","adjustMethod",
        "metaP","weight","nperm","reprod","pcut","logOffset","preset","qcPlots",
        "figFormat","outList","exportWhere","exportWhat","exportScale",
        "exportValues","exportStats","exportCountsTable","restrictCores",
        "report","refdb","reportTop","reportTemplate","verbose","runLog",
        "saveGeneModel","version","localDb"
    ))
}

.getBackwardsValidArgs <- function() {
    return(c(
        "counts","sample.list","exclude.list","file.type","path","contrast",
        "libsize.list","id.col","gc.col","name.col","bt.col","annotation",
        "gene.file","org","trans.level","count.type","utr.flank","exon.filters",
        "gene.filters","when.apply.filter","normalization","norm.args",
        "statistics","stat.args","adjust.method","meta.p","weight","nperm",
        "reprod","pcut","log.offset","preset","qc.plots","fig.format",
        "out.list","export.where","export.what","export.scale","export.values",
        "export.stats","export.counts.table","restrict.cores","report","refdb",
        "report.top","report.template","verbose","run.log","save.gene.model",
        "version","local.db.home"
    ))
}
