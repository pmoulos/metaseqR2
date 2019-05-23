#' Main argument validator
#'
#' Checks if the arguments passed to \code{\link{metaseqr}} are valid and throws
#' a warning about the invalid ones (which are ignored anyway because of the
#' \code{...} in \code{\link{metaseqr}}. However, for this reason this function
#' is useful as some important parameter faults might go unnoticed in the beginning
#' and cause a failure afterwards. Internal use.
#' 
#' @param mainArgs a list of parameters with which metaseqr is called (essentially,
#' the output of \code{\link{match.call}}.
#' @author Panagiotis Moulos
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

#' Text argument validator
#'
#' Checks if one or more given textual argument(s) is/are member(s) of a list of 
#' correct arguments. It's a more package-specific function similar to 
#' \code{\link{match.arg}}. Mostly for internal use.
#' 
#' @param argName the name of the argument that is checked (for display purposes).
#' @param argValue the value(s) of the argument to be checked.
#' @param argList a vector of valid argument values for \code{argValue} to be 
#' matched against.
#' @param multiarg a logical scalar indicating whether \code{argName} accepts 
#' multiple arguments or not. In that case, all of the values in \code{argValue} 
#' are checked against \code{argList}.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' checkTextArgs("count.type",count.type,c("gene","exon"),multiarg=FALSE)
#' checkTextArgs("statistics",statistics,c("deseq","edger","noiseq","bayseq",
#'  "limma"), multiarg=TRUE)
#'}
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
        # An exception must be added for annotation because it can be an external 
        # file too
        if (argName=="annotation") { 
            if (!(argValue %in% argList) && !file.exists(argSave))
                stopwrap("\"",argName,"\""," parameter must be one of ",
                    paste(paste("\"",argList,sep=""),collapse="\", "),
                        "\" or an existing file!")
        }
        else {
            if (!(argValue %in% argList))
                stopwrap("\"",argName,"\""," parameter must be one of ",
                    paste(paste("\"",argList,sep=""),collapse="\", "),"\"!")
        }
    }
}

#' Numeric argument validator
#'
#' Checks if one or more given numeric argument(s) satisfy several rules concerning 
#' numeric arguments, e.g. proper bounds or proper format (e.g. it must be a number 
#' and not a character). Mostly for internal use.
#' 
#' @param argName the name of the argument that is checked (for display purposes).
#' @param argValue the value(s) of the argument to be checked.
#' @param argType either the string \code{"numeric"} to denote generic double-like 
#' R numerics or \code{"integer"} for integer values.
#' @param argBounds a numeric or a vector with 2 elements, restraining 
#' \code{argValue} to be within the bounds defined by the input vector or e.g. 
#' larger (smaller) than the numeric value. See examples.
#' @param direction a string denoting to which direction the \code{argValue} 
#' should be compared with \code{argBounds}. For example, \code{"both"} should 
#' be given with a two element vector against which, \code{argValue} will be 
#' checked to see whether it is smaller than the low boundary or larger than the 
#' higher boundary. In that case, the function will throw an error. The direction 
#' parameter can be one of: \code{"both"} (described above), \code{"botheq"} (as 
#' above, but the \code{arg.val} is also checked for equality -closed intervals), 
#' \code{"gt"} or \code{"gte"} (check whether \code{arg.val} is smaller or smaller 
#' than or equal to the first value of \code{argBounds}), \code{"lt"} or \code{"lte"} 
#' (check whether \code{arg.val} is larger or larger than or equal to the first 
#' value of \code{argBounds}).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' pcut <- 1.2 # A probability cannot be larger than 1! It will throw an error!
#' checkNumArgs("pcut",pcut,"numeric",c(0,1),"botheq")
#' pcut <- 0.05 # Pass
#' checkNumArgs("pcut",pcut,"numeric",c(0,1),"botheq")
#' gc.col <- 3.4 # A column in a file cannot be real! It will throw an error!
#' checkNumArgs("gc.col",gc.col,"integer",0,"gt")
#' gc.col <- 5 # Pass
#' checkNumArgs("gc.col",gc.col,"integer",0,"gt")
#'}
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

#' File argument validator
#'
#' Checks if a file exists for specific arguments requiring a file input. Internal 
#' use only.
#'
#' @param argName argument name to display in a possible error.
#' @param argValue the filename to check.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # OK
#' checkFileArgs("file",system.file("metaseqr_report.html",package="metaseqR"))
#' # Error!
#' checkFileArgs("file",system.file("metaseqr_report.htm",package="metaseqR"))
#'}
checkFileArgs <- function(argName,argValue) {
    if (!file.exists(argValue))
        stopwrap("\"",argName,"\""," parameter must be an existing file!")
}

#' Contrast validator
#'
#' Checks if the contrast vector follows the specified format. Internal use only.
#'
#' @param cnt contrasts vector.
#' @param sampleList the input sample list.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' dontrun{
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' cnt <- c("A_vs_B") # Will work
#' #cnt <- c("A_vs_C") ## Will throw error!
#' checkContrastFormat(cnt,sampleList)
#}
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

#' Library size validator
#'
#' Checks the names of the supplied library sizes. Internal use only.
#'
#' @param libsizeList the samples-names library size list.
#' @param sampleList the input sample list.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' sampleList <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' libsizeList.1 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,B2=1.3e+6,B3=1.5e+6)
#' libsizeList.2 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,B2=1.3e+6)
#' checkLibsize(libsizeList.1,sampleList) # Will work
#' #checkLibsize(libsizeList.2,sampleList) # Will throw error!
#'}
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

#' Required packages validator
#'
#' Checks if all the required packages are present according to metaseqr input 
#' options. Internal use only.
#'
#' @param m meta-analysis method
#' @param p qc plot types
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' checkPackages(c("simes","whitlock"),c("gcbias","correl"))
#}
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
       return(backArgs)
    }
    else
        return(list(backDetected=FALSE,mixDetected=FALSE,args=NULL))
}

.backwardsMapOld2New <- function() {
    oldArgNames <- .getBackwardsValidArgs()
    newArgNames <- .getValidArgs()
    map <- list()
    map[oldArgNames] <- newArgNames
    return(map)
}

.getValidArgs <- function() {
    return(c(
        "counts","sampleList","excludeList","fileType","path","contrast",
        "libsizeList","idCol","gcCol","nameCol","btCol","annotation","geneFile",
        "org","transLevel","countType","utrFlank","exonFilters","geneFilters",
        "whenApplyFilter","normalization","normArgs","statistics","statArgs",
        "adjustMethod","metaP","weight","nperm","reprod","pcut","logOffset",
        "preset","qcPlots","figFormat","outList","exportWhere","exportWhat",
        "exportScale","exportValues","exportStats","exportCountsTable",
        "restrictCores","report","refdb","reportTop","reportTemplate",
        "verbose","runLog","saveGeneModel","version","localDbHome"
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
