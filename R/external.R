.myCreateFilterXMLchunk <- function(filterChunk,mart) {
    individualFilters <- vapply(names(filterChunk),
        FUN=function(filter,values,mart) {
            # if the filter exists and is boolean we do this
            if (filter %in% listFilters(mart,what="name") && 
                grepl('boolean',filterType(filter=filter,mart=mart))) {
                if (!is.logical(values[[filter]])) 
                    stop("'",filter, "' is a boolean filter and needs a ",
                        "corresponding logical value of TRUE or FALSE to ",
                        "indicate if the query should retrieve all data that ",
                        "fulfill the boolean or alternatively that all data ", 
                        "that not fulfill the requirement should be retrieved.",
                        call. = FALSE)
                val <- ifelse(values[[filter]],yes=0,no=1)
                val <- paste0("\" excluded = \"",val,"\" ")
                
            } else { 
                # otherwise the filter isn't boolean, or doesn't exist
                if(is.numeric(values[[filter]])) 
                    values[[filter]] <- as.integer(values[[filter]])
                val <- paste0(values[[filter]], collapse = ",")
                # convert " to ' to avoid truncating the query string
                val <- gsub(x=val,pattern="\"",replacement="'",
                    fixed=TRUE)
                val <- paste0('" value = "',val,'" ')
            }
            filterXML <- paste0("<Filter name = \"",filter,val,"/>")
            return(filterXML)
        },FUN.VALUE=character(1),filterChunk,mart,USE.NAMES=FALSE)
    filterXML <- paste0(individualFilters,collapse="")
    return(filterXML)
}

.mySubmitQueryXML <- function (host, query) {
    httr::set_config(httr::config(http_version = 1L))
    res <- httr::POST(url = host, body = list(query = query),
        httr::timeout(2000))
    if (httr::status_code(res) == 302) {
        host <- stringr::str_match(string=res$all_headers[[1]]$headers$location,
            pattern = "//([a-zA-Z./]+)\\??;?redirectsrc")[, 2]
        res <- httr::POST(url = host, body = list(query = query),
            config = list(httr::timeout(2000)))
    }
    return(suppressMessages(httr::content(res)))
}

.myGenerateFilterXML <- function(filters="",values,mart) {
    # return empty string if no filter specified & this isn't ensembl
    # specifying no filter is generally bad, as it will get 'everything'
    # and we might encounter the time out problem
    if (filters[1] == "") {
        return("")
    }
    # if we have multiple filters, the values must be specified as a list.
    if (length(filters) > 1 && !is(values,"list")) {
        stop("If using multiple filters, the 'value' has to be a list.\n",
            "For example, a valid list for 'value' could be: ",
            "list(affyid=c('1939_at','1000_at'),chromosome='16')\n",
            "Here we select on Affymetrix identifier and chromosome, only ",
            "results that pass both filters will be returned");
    }
    # it's easy to not realise you're passing a data frame here, so check
    if (is.data.frame(values) && ncol(values == 1)) {
        values <- values[,1]
    }
    if (!is.list(values)) {
        values <- list(values)
    }
    names(values) <- filters
    values <- .mySplitValues(list(values),maxChunkSize=1000)
    filterXML_list <- lapply(values,.myCreateFilterXMLchunk,mart)
    return(filterXML_list)
}

.mySplitValues <- function(valuesList,maxChunkSize=500) {
    vLength <- vapply(valuesList[[1]],FUN=length,FUN.VALUE=integer(1))
    
    if(all(vLength <= maxChunkSize)) {
        return(valuesList)
    } else {
        # pick the next filter to split
        vIdx <- min(which(vLength > maxChunkSize))
        nchunks <- (vLength[vIdx] %/% maxChunkSize) + 1
        splitIdx <- rep(seq_len(nchunks),
            each=ceiling(vLength[vIdx]/nchunks))[seq_len(vLength[vIdx])]
        
        # a new list we will populate with the chunks
        tmpList <- list()
        for(i in seq_len(nchunks)) {
            for(j in seq_along(valuesList)) {
                listIdx <- ((i - 1) * length(valuesList)) + j
                tmpList[[listIdx]] <- valuesList[[j]]
                tmpList[[listIdx]][[vIdx]] <-
                    tmpList[[listIdx]][[vIdx]][which(splitIdx == i)]
            }
        }
        # recursively call the function to process next filter
        valuesList <- .mySplitValues(tmpList,maxChunkSize=maxChunkSize)
    }
    return(valuesList)
}

.myGetBM <- function(attributes,filters="",values="",mart,curl=NULL,
    checkFilters=TRUE,verbose=FALSE,uniqueRows=TRUE,bmHeader=FALSE,quote="\"") {
    .myMartCheck(mart)
    if(missing(attributes))
        stop("Argument 'attributes' must be specified.")
    
    if (is.list(filters) && !missing(values))
        warning("Argument 'values' should not be used when argument 'filters'",
            "is a list and will be ignored.")
    if (is.list(filters) && is.null(names(filters)))
        stop("Argument 'filters' must be a named list when sent as a list.")
    if (!is.list(filters) && all(filters != "") && missing(values))
        stop("Argument 'values' must be specified.")
    if (length(filters) > 0 && length(values) == 0)
        stop("Values argument contains no data.")
    
    if (is.list(filters)) {
        values <- filters
        filters <- names(filters)
    }
    
    if (!is(uniqueRows,"logical"))
        stop("Argument 'uniqueRows' must be a logical value, so either TRUE ",
            "or FALSE")
    
    # force the query to return the 'english text' header names with the result
    # we use these later to match and order attribute/column names    
    callHeader <- TRUE
    xmlQuery <- paste0("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query>",
        #"<Query  virtualSchemaName = '",biomaRt:::martVSchema(mart),
        "<Query  virtualSchemaName = '",mart@vschema,
        "' uniqueRows = '",as.numeric(uniqueRows),
        "' count = '0' datasetConfigVersion = '0.6' header='",
        as.numeric(callHeader),"' requestid= 'biomaRt'> <Dataset name = '",
        #biomaRt:::martDataset(mart),"'>")
        mart@dataset,"'>")
    
    # checking the Attributes
    invalid <- !(attributes %in% listAttributes(mart, what="name"))
    if(any(invalid))
        stop(paste("Invalid attribute(s):", paste(attributes[invalid],
            collapse=", "),"\nPlease use the function 'listAttributes' to get ",
            "valid attribute names"))
    
    # attribute are ok lets add them to the query
    attributeXML = paste("<Attribute name = '",attributes, "'/>",collapse="",
        sep="")
    
    #checking the filters
    if (filters[1] != "" && checkFilters) {
        invalid <- !(filters %in% listFilters(mart, what="name"))
        if(any(invalid))
            stop(paste("Invalid filters(s):", paste(filters[invalid],
                collapse=", "),"\nPlease use the function 'listFilters' to ",
                "get valid filter names"))
    }
    
    # filterXML is a list containing filters with reduced numbers of values
    # to meet the 500 value limit in BioMart queries
    filterXmlList <- .myGenerateFilterXML(filters,values,mart)
    
    resultList <- list()
    
    # we submit a query for each chunk of the filter list
    for (i in seq_along(filterXmlList)) {
        filterXML <- filterXmlList[[i]]
        fullXmlQuery <- paste(xmlQuery,attributeXML,filterXML,
            "</Dataset></Query>",sep="")
        
        if(verbose)
            message(fullXmlQuery)
        
        # we choose a separator based on whether '?redirect=no' is present
        #sep <- ifelse(grepl(x=biomaRt:::martHost(mart),
        sep <- ifelse(grepl(x=mart@host,
            pattern=".+\\?.+"), "&", "?")
        
        #postRes <- .mySubmitQueryXML(host=paste0(biomaRt:::martHost(mart),sep),
        postRes <- .mySubmitQueryXML(host=paste0(mart@host,sep),
            query=fullXmlQuery)
        
        if (verbose) {
            writeLines("#################\nResults from server:")
            print(postRes)
        }
        
        if (!(is.character(postRes) && (length(postRes)==1L)))
            stop("The query to the BioMart webservice returned an invalid ",
                "result: biomaRt expected a character string of length 1.\n",
                "Please report this on the support site at", 
                "http://support.bioconductor.org")
        
        if (gsub("\n","",postRes,fixed=TRUE,useBytes=TRUE)== "") { 
            # meaning an empty result
            result <- as.data.frame(matrix("",ncol=length(attributes),nrow=0),
                stringsAsFactors=FALSE)
        } else {
            if (length(grep("^Query ERROR",postRes)) > 0L)
                stop(postRes)
            
            # convert the serialized table into a dataframe
            con <- textConnection(postRes)
            result <- read.table(con,sep="\t",header=callHeader,quote=quote,
                comment.char="",check.names=FALSE,stringsAsFactors=FALSE)
            if(verbose) {
                writeLines("#################\nParsed results:")
                print(result)
            }
            close(con)
            
            if (!(is(result,"data.frame") 
                && (ncol(result)==length(attributes)))) {
                print(head(result))
                stop("The query to the BioMart webservice returned an invalid ",
                    "result: the number of columns in the result table does ",
                    "not equal the number of attributes in the query.\n",
                    "Please report this on the support site at",
                    "http://support.bioconductor.org")
            }
        }
        
        resultList[[i]] <- .mySetResultColNames(result,mart=mart,
            attributes=attributes,bmHeader=bmHeader)
    }
    
    # collate results
    result <- do.call('rbind',resultList)
    return(result)
}

.myMartCheck = function(mart,biomart=NULL) {
    if (missing(mart) || !is(mart,'Mart')) {
        stop("You must provide a valid Mart object. To create a Mart ",
            "object use the function: useMart. Check ?useMart for more ",
            "information.")
    }
    if(!is.null(biomart)){
        martcheck <- mart@biomart
        bmok <- FALSE
        for (k in seq_len(length(biomart))) {
            if(martcheck[1] == biomart[k]) {
                bmok <- TRUE
            }
        }
        if (!bmok) {
            stop(paste("This function only works when used with the ",
                biomart," BioMart.",sep=""))
        }
    }
    if (mart@dataset == "") {
        stop("No dataset selected, please select a dataset first.\nYou can ",
            "see the available datasets by using the listDatasets function\n",
            "See ?listDatasets for more information. Then you should create\n",
            "the Mart object by using the useMart function.\nSee ?useMart ",
            "for more information")
    }
}

.mySetResultColNames <- function(result,mart,attributes,bmHeader=FALSE) {
    # get all avaialble sttributes and 
    # filter only for the ones we've actually asked for
    att <- listAttributes(mart,what=c("name","description"))
    att <- att[which(att[,'name'] %in% attributes),]
    if (length(which(duplicated(att[,'description']))) > 
        length(which(duplicated(att)))) {
        warning("Cannot unambiguously match attribute names\nIgnoring ",
            "bmHeader argument and using biomart description field")
        return(result)
    }
    
    resultNames = colnames(result)
    # match the returned column names with the attribute names
    matches <- match(resultNames,att[,2],NA)
    if (any(is.na(matches))) {
        warning("Problems assigning column names. Currently using ",
            "the biomart description field. You may wish to set these ",
            "manually.")
        return(result)
    }
    # if we want to use the attribute names we specified, do this, 
    # otherwise we use the header returned with the query
    if (!bmHeader) {
        colnames(result) <- att[matches,1]
    }
    # now put things in the order we actually asked for the attributes in
    result <- result[,match(att[matches,1],attributes),drop=FALSE]
    
    return(result)
}

.GENE_TYPES <- c("gene","pseudogene","transposable_element_gene")
.TX_TYPES <- c("transcript","pseudogenic_transcript","primary_transcript",
    "mRNA","ncRNA","rRNA","snoRNA","snRNA","tRNA","tmRNA","miRNA",
    "miRNA_primary_transcript","RNase_P_RNA","RNase_MRP_RNA","SRP_RNA",
    "misc_RNA","antisense_RNA", "antisense","lnc_RNA","antisense_lncRNA",
    "transcript_region","pseudogenic_tRNA","scRNA","guide_RNA",
    "telomerase_RNA","vault_RNA","Y_RNA")
.EXON_TYPES <- c("exon","pseudogenic_exon","coding_exon",
    "five_prime_coding_exon","interior_coding_exon",
    "three_prime_coding_exon","exon_of_single_exon_gene","interior_exon",
    "noncoding_exon","five_prime_noncoding_exon",
    "three_prime_noncoding_exon")
.CDS_TYPES <- c("CDS","transposable_element_CDS","CDS_predicted",
    "edited_CDS")
.STOP_CODON_TYPES <- "stop_codon"
.GFF_FEATURE_TYPES <- c(.GENE_TYPES,.TX_TYPES,.EXON_TYPES,.CDS_TYPES,
    .STOP_CODON_TYPES)
