buildAnnotationDatabase <- function(organisms,sources,
    #db=file.path(path.expand("~"),".metaseqR","annotation.sqlite"),
    db=file.path(system.file(package="metaseqR2"),"annotation.sqlite"),
    forceDownload=TRUE,rc=NULL) {
    if (missing(organisms))
        organisms <- getSupportedOrganisms()
    if (missing(sources))
        sources <- getSupportedRefDbs()
    
    orgIsList <- FALSE
    if (!is.list(organisms) && is.character(organisms)
        && "ensembl" %in% sources)
        warning("When ensembl is in the annotation sources to download, it is ",
            "advised to provide organisms as\na named list with names the ",
            "requested organisms and list members the Ensembl versions.\n",
            "Otherwise, the latest Ensembl version for each organism will be ",
            "used.",immediate.=TRUE)
    if (is.list(organisms)) {
        if (is.null(names(organisms)))
            stop("When organisms is a list, it must be named!")
        orgList <- organisms
        organisms <- names(organisms)
        orgIsList <- TRUE
    }
    
    checkTextArgs("organisms",organisms,getSupportedOrganisms(),multiarg=TRUE)
    checkTextArgs("sources",sources,getSupportedRefDbs(),multiarg=TRUE)
    
    if (!requireNamespace("GenomeInfoDb"))
        stop("R package GenomeInfoDb is required to construct annotation ",
            "stores!")
    
    # Check database path
    if (!dir.exists(dirname(db)))
        dir.create(dirname(db),recursive=TRUE,mode="0755")
    
    # Initialize or open the annotation SQLite datatabase
    message("Opening metaseqR2 SQLite database ",db)
    con <- initDatabase(db)
    
    for (s in sources) {
        for (o in organisms) {
            # Retrieving genome info. We will be inserting the seqinfo to the
            # sqlite multiple times for the sake of simplicity regarding later
            # deletion by foreign key cascade.
            message("Retrieving genome information for ",o," from ",s)
            sf <- getSeqInfo(o)
            
            # Now, we must derive versioning. For Ensembl it's obvious as it has
            # official versions. For UCSC we need to keep a date track. Then 
            # inside metaseqR, if date not provided, it will automatically
            # detect the latest one (as in Ensembl with versions, if not
            # provided).
            if (s == "ensembl") {
                if (orgIsList)
                    vs <- orgList[[o]]
                else {
                    vss <- getUcscToEnsembl(o)
                    vs <- vss[length(vss)]
                }
                vs <- .validateEnsemblVersions(o,vs)
                if (is.null(vs)) { # Something went wrong... Get latest again
                    vss <- getUcscToEnsembl(o)
                    vs <- vss[length(vss)]
                }
            }
            else if (s %in% getSupportedUcscDbs())
                vs <- format(Sys.Date(),"%Y%m%d")
                        
            for (v in vs) {
                #storePath <- file.path(home,s,o,v)
                #if (!dir.exists(storePath))
                #    dir.create(storePath,recursive=TRUE,mode="0755")
                    
                # Retrieve gene annotations
                if (.annotationExists(con,o,s,v,"gene") && !forceDownload)
                #if (file.exists(file.path(storePath,"gene.rda")) 
                #    && !forceDownload)
                    message("Gene annotation for ",o," from ",s," version ",v,
                        " has already been created and will be skipped.\nIf ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving gene annotation for ",o," from ",s,
                        " version ",v)
                    ann <- getAnnotation(o,"gene",refdb=s,ver=v,rc=rc)
                    #gene <- makeGRangesFromDataFrame(
                    #    df=ann,
                    #    seqinfo=sf,
                    #    keep.extra.columns=TRUE,
                    #    seqnames.field="chromosome"
                    #)
                    #save(gene,file=file.path(storePath,"gene.rda"),
                    #    compress=TRUE)
                    
                    # First drop if previously exists
                    nr <- .dropAnnotation(con,o,s,v,"gene")
                    # Then insert to the contents table so as to get the content
                    # id to attach in the annotation table
                    nr <- .insertContent(con,o,s,v,"gene")
                    nid <- .annotationExists(con,o,s,v,"gene",out="id")
                    # If something happens, the whole procedure will break
                    # anyway
                    # Add content_id
                    ann$content_id <- rep(nid,nrow(ann))
                    sfGene <- sf
                    sfGene$content_id <- rep(nid,nrow(sfGene))
                    # Write genes and seqinfo
                    dbWriteTable(con,"gene",ann,row.names=FALSE,append=TRUE)
                    dbWriteTable(con,"seqinfo",sfGene,row.names=FALSE,
                        append=TRUE)
                }
                
                # Retrieve transcript annotations
                if (.annotationExists(con,o,s,v,"transcript") && !forceDownload)
                    message("Transcript annotation for ",o," from ",s,
                        " version ",v," has already been created and will be ",
                        "skipped.\nIf you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    message("Retrieving transcript annotation for ",o,
                        " from ",s," version ",v)
                    ann <- getAnnotation(o,"transcript",refdb=s,ver=v,rc=rc)
                    nr <- .dropAnnotation(con,o,s,v,"transcript")
                    nr <- .insertContent(con,o,s,v,"transcript")
                    nid <- .annotationExists(con,o,s,v,"transcript",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfTranscript <- sf
                    sfTranscript$content_id <- rep(nid,nrow(sfTranscript))
                    dbWriteTable(con,"transcript",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfTranscript,row.names=FALSE,
                        append=TRUE)
                }
                
                # Then summarize the transcripts and write again with type 
                # sum_transcript
                if (.annotationExists(con,o,s,v,"summarized_transcript") 
                    && !forceDownload)
                    message("Summarized transcript annotation for ",o," from ",
                        s," version ",v," has already been created and will ",
                        "be skipped.\nIf you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!.annotationExists(con,o,s,v,"transcript")) 
                        stop("Transcript annotation for ",o," from ",s,
                            " version ",v," is required in order to build ",
                            "predefined merged transcript regions for read ",
                            "counting.\nPlease rerun the ",
                            "buildAnnotationDatabase function with ",
                            "appropriate parameters.")
                    annGr <- .loadPrebuiltAnnotation(con,o,s,v,"transcript",
                        "gene")
                    message("Merging transcripts for ",o," from ",s,
                        " version ",v)
                    #annGr <- reduceTranscripts(annGr)
                    #ann <- as.data.frame(annGr)
                    annList <- reduceTranscripts(annGr)
                    ann <- as.data.frame(annList$model)
                    ann <- ann[,c(1,2,3,6,7,5,8,9)]
                    names(ann)[1] <- "chromosome"
                    ann$chromosome <- as.character(ann$chromosome)
                    ann <- ann[order(ann$chromosome,ann$start),]
                    
                    nr <- .dropAnnotation(con,o,s,v,"summarized_transcript")
                    nr <- .insertContent(con,o,s,v,"summarized_transcript")
                    nid <- .annotationExists(con,o,s,v,"summarized_transcript",
                        out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfSumTranscript <- sf
                    sfSumTranscript$content_id <- rep(nid,nrow(sfSumTranscript))
                    dbWriteTable(con,"summarized_transcript",ann,
                        row.names=FALSE,append=TRUE)
                    dbWriteTable(con,"seqinfo",sfSumTranscript,row.names=FALSE,
                        append=TRUE)
                }
                
                # Retrieve 3' UTR annotations
                if (.annotationExists(con,o,s,v,"utr") && !forceDownload)
                    message("3' UTR annotation for ",o," from ",s," version ",
                        v," has already been created and will be skipped.\nIf ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving 3' UTR annotation for ",o,
                        " from ",s," version ",v)
                    ann <- getAnnotation(o,"utr",refdb=s,ver=v,rc=rc)
                    nr <- .dropAnnotation(con,o,s,v,"utr")
                    nr <- .insertContent(con,o,s,v,"utr")
                    nid <- .annotationExists(con,o,s,v,"utr",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfUtr <- sf
                    sfUtr$content_id <- rep(nid,nrow(sfUtr))
                    dbWriteTable(con,"utr",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfUtr,row.names=FALSE,
                        append=TRUE)
                }
                
                # Then summarize the 3'utrs per gene and write again with type 
                # sum_utr
                if (.annotationExists(con,o,s,v,"summarized_3utr")
                    && !forceDownload)
                    message("Summarized 3' UTR annotation for ",o," from ",
                        s," version ",v," has already been created and will ",
                        "be skipped.\nIf you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!.annotationExists(con,o,s,v,"utr")) 
                        stop("3' UTR annotation for ",o," from ",s," version ",
                            vs," is required in order to build predefined ",
                            "merged 3' UTR regions for read counting.\nPlease ",
                            "rerun the buildAnnotationDatabase function with ",
                            "appropriate parameters.")
                    annGr <- .loadPrebuiltAnnotation(con,o,s,v,"gene","utr")
                    message("Merging gene 3' UTRs for ",o," from ",s,
                        " version ",v)
                    #annGr <- reduceTranscripts(annGr)
                    #ann <- as.data.frame(annGr)
                    annList <- reduceTranscripts(annGr)
                    ann <- as.data.frame(annList$model)
                    ann <- ann[,c(1,2,3,6,7,5,8,9)]
                    names(ann)[1] <- "chromosome"
                    ann$chromosome <- as.character(ann$chromosome)
                    ann <- ann[order(ann$chromosome,ann$start),]
                    nr <- .dropAnnotation(con,o,s,v,"summarized_3utr")
                    nr <- .insertContent(con,o,s,v,"summarized_3utr")
                    nid <- .annotationExists(con,o,s,v,"summarized_3utr",
                        out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfSumUtr <- sf
                    sfSumUtr$content_id <- rep(nid,nrow(sfSumUtr))
                    dbWriteTable(con,"summarized_3utr",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfSumUtr,row.names=FALSE,
                        append=TRUE)
                    
                    activeLength <- annList$length
                    nr <- .dropAnnotation(con,o,s,v,"active_utr_length")
                    nr <- .insertContent(con,o,s,v,"active_utr_length")
                    nid <- .annotationExists(con,o,s,v,"active_utr_length",
                        out="id")
                    active <- data.frame(
                        name=names(activeLength),
                        length=activeLength,
                        content_id=rep(nid,length(activeLength))
                    )
                    dbWriteTable(con,"active_utr_length",active,row.names=FALSE,
                        append=TRUE)
                }
                
                # Then summarize the 3'utrs per transcript and write again with 
                # type sum_utr
                if (.annotationExists(con,o,s,v,"summarized_3utr_transcript")
                    && !forceDownload)
                    message("Summarized 3' UTR annotation per transcript for ",
                        o," from ",s," version ",v," has already been created ",
                        "and will be skipped.\nIf you wish to recreate it ",
                        "choose forceDownload = TRUE.")
                else {
                    if (!.annotationExists(con,o,s,v,"utr"))
                        stop("3' UTR annotation per transcript for ",o," from ",
                            s," version ",v," is required in order to build ",
                            "predefined merged 3'UTR regions for read ",
                            "counting.\nPlease rerun the ",
                            "buildAnnotationDatabase function with ",
                            "appropriate parameters.")
                    annGr <- 
                        .loadPrebuiltAnnotation(con,o,s,v,"transcript","utr")
                    message("Merging transcript 3' UTRs for ",o," from ",s,
                        " version ",v)
                    #annGr <- reduceTranscriptsUtr(annGr)
                    #ann <- as.data.frame(annGr)
                    annList <- reduceTranscriptsUtr(annGr)
                    ann <- as.data.frame(annList$model)
                    ann <- ann[,c(1,2,3,6,7,5,8,9)]
                    names(ann)[1] <- "chromosome"
                    ann$chromosome <- as.character(ann$chromosome)
                    ann <- ann[order(ann$chromosome,ann$start),]
                    nr <- 
                        .dropAnnotation(con,o,s,v,"summarized_3utr_transcript")
                    nr <- .insertContent(con,o,s,v,"summarized_3utr_transcript")
                    nid <- .annotationExists(con,o,s,v,
                        "summarized_3utr_transcript",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfSumUtrTranscript <- sf
                    sfSumUtrTranscript$content_id <- 
                        rep(nid,nrow(sfSumUtrTranscript))
                    dbWriteTable(con,"summarized_3utr_transcript",ann,
                        row.names=FALSE,append=TRUE)
                    dbWriteTable(con,"seqinfo",sfSumUtrTranscript,
                        row.names=FALSE,append=TRUE)
                    
                    activeLength <- annList$length
                    nr <- .dropAnnotation(con,o,s,v,"active_trans_utr_length")
                    nr <- .insertContent(con,o,s,v,"active_trans_utr_length")
                    nid <- .annotationExists(con,o,s,v,
                        "active_trans_utr_length",out="id")
                    active <- data.frame(
                        name=names(activeLength),
                        length=activeLength,
                        content_id=rep(nid,length(activeLength))
                    )
                    dbWriteTable(con,"active_trans_utr_length",active,
                        row.names=FALSE,append=TRUE)
                }
                
                # Retrieve exon annotations
                if (.annotationExists(con,o,s,v,"exon") && !forceDownload)
                    message("Exon annotation for ",o," from ",s," version ",v,
                        " has already been created and will be skipped.\nIf ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving exon annotation for ",o," from ",s,
                        " version ",v)
                    ann <- getAnnotation(o,"exon",refdb=s,ver=v,rc=rc)
                    nr <- .dropAnnotation(con,o,s,v,"exon")
                    nr <- .insertContent(con,o,s,v,"exon")
                    nid <- .annotationExists(con,o,s,v,"exon",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfExon <- sf
                    sfExon$content_id <- rep(nid,nrow(sfExon))
                    dbWriteTable(con,"exon",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfExon,row.names=FALSE,
                        append=TRUE)
                }
                
                # Retrieve extended annotations
                if (.annotationExists(con,o,s,v,"transexon") && !forceDownload)
                    message("Extended exon annotation for ",o," from ",s,
                        " version ",v," has already been created and will be ",
                        "skipped.\nIf you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    message("Retrieving extended exon annotation for ",o,
                        " from ",s," version ",v)
                    ann <- getAnnotation(o,"transexon",refdb=s,ver=v,rc=rc)
                    nr <- .dropAnnotation(con,o,s,v,"transexon")
                    nr <- .insertContent(con,o,s,v,"transexon")
                    nid <- .annotationExists(con,o,s,v,"transexon",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfTrExon <- sf
                    sfTrExon$content_id <- rep(nid,nrow(sfExon))
                    dbWriteTable(con,"transexon",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfTrExon,row.names=FALSE,
                        append=TRUE)
                }
                
                # Then summarize the exons per gene and write again with type 
                # summarized_exon
                if (.annotationExists(con,o,s,v,"summarized_exon")
                    && !forceDownload)
                    message("Summarized exon annotation for ",o," from ",s,
                        " version ",v," has already been created and will be ",
                        "skipped.\nIf you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!.annotationExists(con,o,s,v,"exon")) 
                        stop("Exon annotation for ",o," from ",s," version ",v,
                            " is required in order to build predefined merged ",
                            "exon regions for RNA-Seq (exon) coverage ",
                            "calculations.\nPlease rerun the ",
                            "buildAnnotationDatabase function with ",
                            "appropriate parameters.")
                            
                    annGr <- .loadPrebuiltAnnotation(con,o,s,v,"exon","exon")
                    message("Merging exons for ",o," from ",s," version ",v)
                    annList <- reduceExons(annGr)
                    ann <- as.data.frame(annList$model)
                    ann <- ann[,c(1,2,3,6,7,5,8,9)]
                    names(ann)[1] <- "chromosome"
                    ann$chromosome <- as.character(ann$chromosome)
                    ann <- ann[order(ann$chromosome,ann$start),]
                    nr <- .dropAnnotation(con,o,s,v,"summarized_exon")
                    nr <- .insertContent(con,o,s,v,"summarized_exon")
                    nid <- .annotationExists(con,o,s,v,"summarized_exon",
                        out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfSumExon <- sf
                    sfSumExon$content_id <- rep(nid,nrow(sfSumExon))
                    dbWriteTable(con,"summarized_exon",ann,row.names=FALSE,
                        append=TRUE)
                    dbWriteTable(con,"seqinfo",sfSumExon,row.names=FALSE,
                        append=TRUE)        
                            
                    activeLength <- annList$length
                    nr <- .dropAnnotation(con,o,s,v,"active_length")
                    nr <- .insertContent(con,o,s,v,"active_length")
                    nid <- .annotationExists(con,o,s,v,"active_length",
                        out="id")
                    active <- data.frame(
                        name=names(activeLength),
                        length=activeLength,
                        content_id=rep(nid,length(activeLength))
                    )
                    dbWriteTable(con,"active_length",active,row.names=FALSE,
                        append=TRUE)
                }
                
                # Then summarize the exons per transcript and write again with 
                # type summarized_transcript_exon
                if (.annotationExists(con,o,s,v,"summarized_transcript_exon")
                    && !forceDownload)
                    message("Summarized exon annotation per transcript for ",o,
                        " from ",s," version ",v," has already been created ",
                        "and will be skipped.\nIf you wish to recreate it ",
                        "choose forceDownload = TRUE.")
                else {
                    if (!.annotationExists(con,o,s,v,"transexon")) 
                        stop("Extended exon annotation for ",o," from ",s,
                            " version ",v," is required in order to build ",
                            "summarized exon annotation per transcript\n",
                            "Please rerun the buildAnnotationDatabase ",
                            "function with appropriate parameters.")
                
                    annGr <- .loadPrebuiltAnnotation(con,o,s,v,"transcript",
                        "exon")
                    message("Merging exons for ",o," from ",s," version ",v)
                    annList <- reduceTranscriptsExons(annGr)
                    ann <- as.data.frame(annList$model)
                    ann <- ann[,c(1,2,3,6,7,5,8,9)]
                    names(ann)[1] <- "chromosome"
                    ann$chromosome <- as.character(ann$chromosome)
                    ann <- ann[order(ann$chromosome,ann$start),]
                    nr <- .dropAnnotation(con,o,s,v,
                        "summarized_transcript_exon")
                    nr <- .insertContent(con,o,s,v,"summarized_transcript_exon")
                    nid <- .annotationExists(con,o,s,v,
                        "summarized_transcript_exon",out="id")
                    ann$content_id <- rep(nid,nrow(ann))
                    sfSumTrExon <- sf
                    sfSumTrExon$content_id <- rep(nid,nrow(sfSumTrExon))
                    dbWriteTable(con,"summarized_transcript_exon",ann,
                        row.names=FALSE,append=TRUE)
                    dbWriteTable(con,"seqinfo",sfSumTrExon,row.names=FALSE,
                        append=TRUE)        
                            
                    activeLength <- annList$length
                    nr <- .dropAnnotation(con,o,s,v,"active_trans_exon_length")
                    nr <- .insertContent(con,o,s,v,"active_trans_exon_length")
                    nid <- .annotationExists(con,o,s,v,
                        "active_trans_exon_length",out="id")
                    active <- data.frame(
                        name=names(activeLength),
                        length=activeLength,
                        content_id=rep(nid,length(activeLength))
                    )
                    dbWriteTable(con,"active_trans_exon_length",active,
                        row.names=FALSE,append=TRUE)
                }
            }
        }
    }
    
    dbDisconnect(con)
}

# GTF only!
buildCustomAnnotation <- function(gtfFile,metadata,
    db=file.path(system.file(package="metaseqR2"),"annotation.sqlite"),
    #db=file.path(path.expand("~"),".metaseqR","annotation.sqlite"),
    rewrite=TRUE) {
    # Check metadata
    if (is.null(metadata$organism))
        stop("An organism name must be provided with metadata!")
    if (is.null(metadata$source)) {
        warning("A source should be provided with metadata! Using 'inhouse'...",
            immediate.=TRUE)
        metadata$source <- "inhouse"
    }
    if (is.null(metadata$version)) {
        warning("A version should be provided with metadata! Using today...",
            immediate.=TRUE)
        metadata$version <- format(Sys.Date(),"%Y%m%d")
    }
    if (is.null(metadata$chromInfo)) {
        warning("Chromosomal lengths should be provided with metadata! ",
            "Only chromosome names will be available... ",immediate.=TRUE)
        metadata$chromInfo <- NULL
    }
    else {
        str <- metadata$chromInfo
        if (is.character(str)) {
            out <- tryCatch(open(Rsamtools::BamFile(str)),error=function(e) e)
            if (inherits(out,"error")) # Not a BAM file, try to read.delim
                metadata$chromInfo <- read.delim(str,row.names=1)
            else
                metadata$chromInfo <- .chromInfoFromBAM(str)
        }
        if (!is.data.frame(metadata$chromInfo))
            stop("metadata$chromInfo must be a data frame!")
    }
    
    # Check database path
    if (!dir.exists(dirname(db)))
        dir.create(dirname(db),recursive=TRUE,mode="0755")
    
    # Initialize or open the annotation SQLite datatabase
    message("Opening metaseqR2 SQLite database ",db)
    con <- initDatabase(db)
    
    parsed <- parseCustomGtf(gtfFile)
    
    s <- tolower(metadata$source)
    o <- tolower(metadata$organism)
    v <- metadata$version
        
    # Retrieve gene annotations
    if (.annotationExists(con,o,s,v,"gene") && !rewrite)
        message("Gene annotation for ",o," from ",s," version ",v," has ",
            "already been created and will be skipped.\nIf you wish to ",
            "recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving gene annotation for ",o," from ",s," version ",v,
            " from ",gtfFile)
        ann <- annotationFromCustomGtf(parsed,level="gene",type="gene",
            asdf=TRUE)
        nr <- .dropAnnotation(con,o,s,v,"gene")
        nr <- .insertContent(con,o,s,v,"gene",1)
        nid <- .annotationExists(con,o,s,v,"gene",out="id")
        ann$content_id <- rep(nid,nrow(ann))
        sfGene <- .chromInfoToSeqInfoDf(metadata$chromInfo)
        sfGene$content_id <- rep(nid,nrow(sfGene))
        dbWriteTable(con,"gene",ann,row.names=FALSE,append=TRUE)
        dbWriteTable(con,"seqinfo",sfGene,row.names=FALSE,append=TRUE)
    }
    
    # Retrieve transcript annotations
    if (.annotationExists(con,o,s,v,"transcript") && !rewrite)
        message("Transcript annotation for ",o," from ",s," version ",v," has ",
            "already been created and will be skipped.\nIf you wish to ",
            "recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving transcript annotation for ",o," from ",s,
            " version ",v)
        ann <- annotationFromCustomGtf(parsed,level="transcript",type="gene",
            asdf=TRUE)
        nr <- .dropAnnotation(con,o,s,v,"transcript")
        nr <- .insertContent(con,o,s,v,"transcript",1)
        nid <- .annotationExists(con,o,s,v,"transcript",out="id")
        ann$content_id <- rep(nid,nrow(ann))
        sfTranscript <- .chromInfoToSeqInfoDf(metadata$chromInfo)
        sfTranscript$content_id <- rep(nid,nrow(sfTranscript))
        dbWriteTable(con,"transcript",ann,row.names=FALSE,append=TRUE)
        dbWriteTable(con,"seqinfo",sfTranscript,row.names=FALSE,append=TRUE)
    }
    
    # Then summarize the transcripts and write again with type sum_transcript
    if (.annotationExists(con,o,s,v,"summarized_transcript") && !rewrite)
        message("Summarized transcript annotation for ",o," from ",s,
            " version ",v," has already been created and will be skipped.\nIf ",
            "you wish to recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving summarized transcript annotation for ",o," from ",s,
            " version ",v)
        ann <- annotationFromCustomGtf(parsed,level="transcript",type="gene",
            summarized=TRUE,asdf=TRUE)
        nr <- .dropAnnotation(con,o,s,v,"summarized_transcript")
        nr <- .insertContent(con,o,s,v,"summarized_transcript",1)
        nid <- .annotationExists(con,o,s,v,"summarized_transcript",out="id")
        ann$content_id <- rep(nid,nrow(ann))
        sfSumTranscript <- .chromInfoToSeqInfoDf(metadata$chromInfo)
        sfSumTranscript$content_id <- rep(nid,nrow(sfSumTranscript))
        dbWriteTable(con,"summarized_transcript",ann,row.names=FALSE,
            append=TRUE)
        dbWriteTable(con,"seqinfo",sfSumTranscript,row.names=FALSE,append=TRUE)
    }
    
    # Retrieve 3' UTR annotations
    if (.annotationExists(con,o,s,v,"utr") && !rewrite)
        message("3' UTR annotation for ",o," from ",s," version ",v," has ",
            "already been created and will be skipped.\nIf you wish to ",
            "recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving 3' UTR annotation for ",o," from ",s," version ",v)
        ann <- annotationFromCustomGtf(parsed,level="gene",type="utr",asdf=TRUE)
        if (nrow(ann) > 0) {
            nr <- .dropAnnotation(con,o,s,v,"utr")
            nr <- .insertContent(con,o,s,v,"utr",1)
            nid <- .annotationExists(con,o,s,v,"utr",out="id")
            ann$content_id <- rep(nid,nrow(ann))
            sfUtr <- .chromInfoToSeqInfoDf(metadata$chromInfo)
            sfUtr$content_id <- rep(nid,nrow(sfUtr))
            dbWriteTable(con,"utr",ann,row.names=FALSE,append=TRUE)
            dbWriteTable(con,"seqinfo",sfUtr,row.names=FALSE,append=TRUE)
        }
        else
            message("3' UTR annotation for ",o," from ",s," version ",v," is ",
                "not available in the provided GTF file.")
    }
    
    # Then summarize the 3'UTRs and write again with type sum_transcript
    if (.annotationExists(con,o,s,v,"summarized_3utr") && !rewrite)
        message("Summarized 3' UTR annotation for ",o," from ",s," version ",v,
            " has already been created and will be skipped.\nIf you wish to ",
            "recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving summarized 3' UTR annotation per gene for ",o,
            " from ",s," version ",v)
        ann <- annotationFromCustomGtf(parsed,level="gene",type="utr",
            summarized=TRUE,asdf=TRUE)
        
        if (nrow(ann) > 0) {
            activeLength <- attr(ann,"activeLength")
            
            nr <- .dropAnnotation(con,o,s,v,"summarized_3utr")
            nr <- .insertContent(con,o,s,v,"summarized_3utr",1)
            nid <- .annotationExists(con,o,s,v,"summarized_3utr",out="id")
            ann$content_id <- rep(nid,nrow(ann))
            sfSumUtr <- .chromInfoToSeqInfoDf(metadata$chromInfo)
            sfSumUtr$content_id <- rep(nid,nrow(sfSumUtr))
            dbWriteTable(con,"summarized_3utr",ann,row.names=FALSE,append=TRUE)
            dbWriteTable(con,"seqinfo",sfSumUtr,row.names=FALSE,append=TRUE)
            
            nr <- .dropAnnotation(con,o,s,v,"active_utr_length")
            nr <- .insertContent(con,o,s,v,"active_utr_length",1)
            nid <- .annotationExists(con,o,s,v,"active_utr_length",out="id")
            active <- data.frame(
                name=names(activeLength),
                length=activeLength,
                content_id=rep(nid,length(activeLength))
            )
            dbWriteTable(con,"active_utr_length",active,row.names=FALSE,
                append=TRUE)
        }
        else
            message("3' UTR annotation for ",o," from ",s," version ",v," is ",
                "not available in the provided GTF file.")
    }
    
    # Then summarize the 3'utrs per transcript and write again with 
    # type sum_utr
    if (.annotationExists(con,o,s,v,"summarized_3utr_transcript") && !rewrite)
        message("Summarized 3' UTR annotation per transcript for ",o," from ",s,
            " version ",v," has already been created and will be skipped.\nIf ",
            "you wish to recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving summarized 3' UTR annotation per transcript for ",o,
            " from ",s," version ",v)
        ann <- annotationFromCustomGtf(parsed,level="transcript",type="utr",
            summarized=TRUE,asdf=TRUE)
        
        if (nrow(ann) > 0) {
            activeLength <- attr(ann,"activeLength")
            nr <- .dropAnnotation(con,o,s,v,"summarized_3utr_transcript")
            nr <- .insertContent(con,o,s,v,"summarized_3utr_transcript",1)
            nid <- .annotationExists(con,o,s,v,"summarized_3utr_transcript",
                out="id")
            ann$content_id <- rep(nid,nrow(ann))
            sfSumUtrTranscript <- .chromInfoToSeqInfoDf(metadata$chromInfo)
            sfSumUtrTranscript$content_id <- rep(nid,nrow(sfSumUtrTranscript))
            dbWriteTable(con,"summarized_3utr_transcript",ann,row.names=FALSE,
                append=TRUE)
            dbWriteTable(con,"seqinfo",sfSumUtrTranscript,row.names=FALSE,
                append=TRUE)
            
            nr <- .dropAnnotation(con,o,s,v,"active_trans_utr_length")
            nr <- .insertContent(con,o,s,v,"active_trans_utr_length",1)
            nid <- .annotationExists(con,o,s,v,"active_trans_utr_length",
                out="id")
            active <- data.frame(
                name=names(activeLength),
                length=activeLength,
                content_id=rep(nid,length(activeLength))
            )
            dbWriteTable(con,"active_trans_utr_length",active,row.names=FALSE,
                append=TRUE)
        }
        else
            message("3' UTR annotation for ",o," from ",s," version ",v," is ",
                "not available in the provided GTF file.")
    }
    
    # Retrieve exon annotations
    if (.annotationExists(con,o,s,v,"exon") && !rewrite)
        message("Exon annotation for ",o," from ",s," version ",v," has ",
            "already been created and will be skipped.\nIf you wish to ",
            "recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving exon annotation for ",o," from ",s," version ",v)
        ann <- annotationFromCustomGtf(parsed,level="exon",type="exon",
            asdf=TRUE)
        nr <- .dropAnnotation(con,o,s,v,"exon")
        nr <- .insertContent(con,o,s,v,"exon",1)
        nid <- .annotationExists(con,o,s,v,"exon",out="id")
        ann$content_id <- rep(nid,nrow(ann))
        sfExon <- .chromInfoToSeqInfoDf(metadata$chromInfo)
        sfExon$content_id <- rep(nid,nrow(sfExon))
        dbWriteTable(con,"exon",ann,row.names=FALSE,append=TRUE)
        dbWriteTable(con,"seqinfo",sfExon,row.names=FALSE,append=TRUE)
    }
    
    # Then summarize the exons and write again with type sum_exon
    if (.annotationExists(con,o,s,v,"summarized_exon") && !rewrite)
        message("Summarized exon annotation for ",o," from ",s," version ",v,
            " has already been created and will be skipped.\nIf you wish to ",
            "to recreate it choose rewrite = TRUE.")
    else {
        message("Retrieving summarized exon annotation for ",o," from ",s,
            " version ",v)
        ann <- annotationFromCustomGtf(parsed,level="gene",type="exon",
            summarized=TRUE,asdf=TRUE)
        activeLength <- attr(ann,"activeLength")
        
        nr <- .dropAnnotation(con,o,s,v,"summarized_exon")
        nr <- .insertContent(con,o,s,v,"summarized_exon",1)
        nid <- .annotationExists(con,o,s,v,"summarized_exon",out="id")
        ann$content_id <- rep(nid,nrow(ann))
        sfSumExon <- .chromInfoToSeqInfoDf(metadata$chromInfo)
        sfSumExon$content_id <- rep(nid,nrow(sfSumExon))
        dbWriteTable(con,"summarized_exon",ann,row.names=FALSE,append=TRUE)
        dbWriteTable(con,"seqinfo",sfSumExon,row.names=FALSE,append=TRUE)

        nr <- .dropAnnotation(con,o,s,v,"active_length")
        nr <- .insertContent(con,o,s,v,"active_length",1)
        nid <- .annotationExists(con,o,s,v,"active_length",out="id")
        active <- data.frame(
            name=names(activeLength),
            length=activeLength,
            content_id=rep(nid,length(activeLength))
        )
        dbWriteTable(con,"active_length",active,row.names=FALSE,append=TRUE)
    }
}

.chromInfoFromBAM <- function(bam) {
    # Danger of including non-canonical chromosomes
    b <- scanBamHeader(bam)
    ci <- as.data.frame(b[[bam]]$targets)
    names(ci) <- "length"
    return(ci)
}

.chromInfoFromSeqinfo <- function(sf) {
    sf <- as.data.frame(sf)
    sf <- sf[,1,drop=FALSE]
    names(sf) <- "length"
    return(sf)
}

# Load annotation must be capable of reading custom annotation files if imported
# with metaseqR facilities
loadAnnotation <- function(genome,refdb,level=c("gene","transcript","exon"),
    type=c("gene","exon","utr"),version="auto",
    db=file.path(system.file(package="metaseqR2"),"annotation.sqlite"),
    summarized=FALSE,asdf=FALSE,rc=NULL) {
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to load annotation from database!")
    
    genome <- tolower(genome[1])
    refdb <- tolower(refdb[1])
    level <- level[1]
    checkTextArgs("level",level,c("gene","transcript","exon"),multiarg=FALSE)
    type <- type[1]
    checkTextArgs("type",type,c("gene","exon","utr"),multiarg=FALSE)
    if (version != "auto")
        checkNumArgs("version",version,"numeric")
    
    # Check if local storage has been set
    onTheFly <- FALSE
    if (file.exists(db)) {
        # Open the connection
        drv <- dbDriver("SQLite")
        con <- dbConnect(drv,dbname=db)
        
        # Is the general resource (organism, source) installed?
        if (!.annotationExists(con,genome,refdb)) {
            warning("The requested annotation does not seem to exist in the ",
                "database! It will be loaded on the fly.\nConsider importing ",
                "it by using buildAnnotationDatabase.")
            onTheFly <- TRUE
        }
        
        # If main source exists, decide on version
        if (!onTheFly) {
            if (version != "auto") {
                if (!.annotationExists(con,genome,refdb,version,
                    .annotationTypeFromInputArgs(level,type))) {
                    warning("The requested annotation version does not seem ",
                        "to exist! Have you run buildAnnotationDatabase or ",
                        "possibly mispelled? Will use newest existing version.",
                        immediate.=TRUE)
                    version <- "auto"
                }
            }
            if (version == "auto") {
                # Check if annotation exists, has been performed before
                vers <- .installedVersions(con,genome,refdb)
                vers <- sort(vers,decreasing=TRUE)
                version <- vers[1]
            }
            ann <- .loadPrebuiltAnnotation(con,genome,refdb,version,level,type,
                summarized)
            dbDisconnect(con)
            
            if (asdf) {
                a <- attr(ann,"activeLength")
                ann <- as.data.frame(unname(ann))
                ann <- ann[,c(1,2,3,6,7,5,8,9)]
                names(ann)[1] <- "chromosome"
                if (!is.null(a))
                    attr(ann,"activeLength") <- a
                return(ann)
            }
            else
                return(ann)
        }
    }
    else
        onTheFly <- TRUE
        
    if (onTheFly) {
        if (genome %in% getSupportedOrganisms()
            && refdb %in% getSupportedRefDbs()) {
            message("Getting latest annotation on the fly for ",genome," from ",
                refdb)
            ann <- .loadAnnotationOnTheFly(genome,refdb,level,type,rc)
            if (asdf) {
                a <- attr(ann,"activeLength")
                ann <- as.data.frame(unname(ann))
                ann <- ann[,c(1,2,3,6,7,5,8,9)]
                names(ann)[1] <- "chromosome"
                if (!is.null(a))
                    attr(ann,"activeLength") <- a
                return(ann)
            }
            else
                return(ann)
        }
        else {
            stop("genome and refdb not in supported automatically download ",
                "annotation options. Please use importCustomAnnotation.")
        }
    }
}

.loadPrebuiltAnnotation <- function(con,genome,refdb,version,level,type,
    summarized=FALSE) {
    metaType <- .annotationTypeFromInputArgs(level,type,summarized)
    cid <- .annotationExists(con,genome,refdb,version,metaType,out="id")
    if (metaType == "summarized_exon")
        tName <- "active_length"
    else if (metaType == "summarized_3utr")
        tName <- "active_utr_length"
    else if (metaType == "summarized_3utr_transcript")
        tName <- "active_trans_utr_length"
    else if (metaType == "summarized_transcript_exon")
        tName <- "active_trans_exon_length"
    cida <- .annotationExists(con,genome,refdb,version,"active_length",out="id")
        
    querySet <- .makeAnnotationQuerySet(metaType,cid,cida)
    
    preAnn <- dbGetQuery(con,querySet$main)
    preAnn$`_id` <- NULL
    preAnn$content_id <- NULL
    ann <- GRanges(preAnn)
    seqlevels(ann) <- unique(preAnn$chromosome)
    
    preSf <- dbGetQuery(con,querySet$seqinfo)
    preSf$`_id` <- NULL
    preSf$content_id <- NULL
    rownames(preSf) <- as.character(preSf[,1])
    if (genome %in% getSupportedOrganisms())
        sfg <- getUcscOrganism(genome)
    else
        sfg <- genome
    sf <- Seqinfo(seqnames=preSf[,1],seqlengths=preSf[,2],
        isCircular=rep(FALSE,nrow(preSf)),genome=sfg)
    
    if (length(sf) > 0) {
        if (length(seqlevels(ann)) != length(seqlevels(sf)))
            # If a subset, this is enough
            seqinfo(ann) <- sf[intersect(seqlevels(ann),seqlevels(sf))]
        else if (!all(seqlevels(ann) == seqlevels(sf))) {
            # Must also be sorted in the same way
            seqlevels(ann) <- seqlevels(sf)
            seqinfo(ann) <- sf
        }
        else
            seqinfo(ann) <- sf
    }
    
    ann <- .nameAnnotationFromMetaType(ann,metaType)
    
    preActive <- NULL
    if (!is.null(querySet$active)) {
        preActive <- dbGetQuery(con,querySet$active)
        preActive$`_id` <- NULL
        preActive$content_id <- NULL
        active <- as.integer(preActive$length)
        names(active) <- as.character(preActive$name)
        # FIXME: Will have to take care later of this
        #active <- active[names(ann)]
        attr(ann,"activeLength") <- active
    }
    
    return(ann)
}

.loadAnnotationOnTheFly <- function(genome,refdb,level,type,rc=NULL) {
    message("Retrieving genome information for ",genome," from ",refdb)
    sf <- getSeqInfo(genome,asSeqinfo=TRUE)
    
    switch(level,
        gene = {
            switch(type,
                gene = {
                    message("Retrieving latest gene annotation for ",genome,
                        " from ",refdb)
                    ann <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
                    annGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    names(annGr) <- as.character(annGr$gene_id)
                },
                exon = {
                    message("Retrieving latest exon annotation for ",genome,
                        " from ",refdb)
                    ann <- getAnnotation(genome,"exon",refdb=refdb,rc=rc)
                    tmpGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    message("Merging exons for ",genome," from ",refdb)
                    annList <- reduceExons(tmpGr)
                    annGr <- annList$model
                    names(annGr) <- as.character(annGr$exon_id)
                    activeLength <- annList$length
                    names(activeLength) <- unique(annGr$gene_id)
                    attr(annGr,"activeLength") <- activeLength
                },
                utr = {
                    message("Retrieving latest 3' UTR annotation for ",genome,
                        " from ",refdb)
                    ann <- getAnnotation(genome,"utr",refdb=refdb,rc=rc)
                    tmpGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    message("Merging 3' UTRs for ",genome," from ",refdb)
                    #annGr <- reduceTranscripts(tmpGr)
                    #names(annGr) <- as.character(annGr$transcript_id)
                    annList <- reduceExons(tmpGr)
                    annGr <- annList$model
                    names(annGr) <- as.character(annGr$transcript_id)
                    activeLength <- annList$length
                    names(activeLength) <- unique(annGr$gene_id)
                    attr(annGr,"activeLength") <- activeLength
                }
            )
        },
        transcript = {
            switch(type,
                gene = {
                    message("Retrieving latest transcript annotation for ",
                        genome," from ",refdb)
                    ann <- getAnnotation(genome,"transcript",refdb=refdb,rc=rc)
                    annGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    names(annGr) <- as.character(annGr$transcript_id)
                },
                exon = {
                    # Stub
                    # TODO: "summarized_exon_by_transcript.rda"
                },
                utr = {
                    message("Retrieving latest 3' UTR annotation for ",genome,
                        " from ",refdb)
                    ann <- getAnnotation(genome,"utr",refdb=refdb,rc=rc)
                    annGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    message("Merging 3' UTRs for ",genome," from ",refdb)
                    #annGr <- reduceTranscriptsUtr(annGr)
                    #names(annGr) <- as.character(annGr$transcript_id)
                    annList <- reduceTranscriptsUtr(tmpGr)
                    annGr <- annList$model
                    names(annGr) <- as.character(annGr$transcript_id)
                    activeLength <- annList$length
                    names(activeLength) <- unique(annGr$transcript_id)
                    attr(annGr,"activeLength") <- activeLength
                }
            )
        },
        exon = {
            switch(type,
                exon = {
                    message("Retrieving latest exon annotation for ",genome,
                        " from ",refdb)
                    ann <- getAnnotation(genome,"exon",refdb=refdb,rc=rc)
                    annGr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    names(annGr) <- as.character(annGr$exon_id)
                }
            )
        }
    )
    
    return(annGr)
}

.nameAnnotationFromMetaType <- function(ann,type) {
    switch(type,
        gene = {
            names(ann) <- as.character(ann$gene_id)
        },
        summarized_exon = {
            names(ann) <- as.character(ann$exon_id)
        },
        exon = {
            names(ann) <- as.character(ann$exon_id)
        },
        summarized_3utr = {
            names(ann) <- as.character(ann$transcript_id)
        },
        utr = {
            names(ann) <- as.character(ann$transcript_id)
        },
        summarized_transcript = {
            names(ann) <- as.character(ann$transcript_id)
        },
        transcript = {
            names(ann) <- as.character(ann$transcript_id)
        },
        summarized_3utr_transcript = {
            names(ann) <- as.character(ann$transcript_id)
        }
    )
    return(ann)
}

.annotationTypeFromInputArgs <- function(level,type,summarized=FALSE) {
    switch(level,
        gene = {
            switch(type,
                gene = {
                    return("gene")
                },
                exon = {
                    if (summarized)
                        return("summarized_exon")
                    else
                        return("exon")
                },
                utr = {
                    if (summarized)
                        return("summarized_3utr")
                    else
                        return("utr")
                }
            )
        },
        transcript = {
            switch(type,
                gene = {
                    if (summarized)
                        return("summarized_transcript")
                    else
                        return("transcript")
                },
                exon = {
                    if (summarized)
                        return("summarized_transcript_exon")
                    else
                        return("transexon")
                },
                utr = {
                    if (summarized)
                        return("summarized_3utr_transcript")
                    else
                        return("utr")
                }
            )
        },
        exon = {
            switch(type,
                gene = {
                    if (summarized)
                        return("summarized_exon")
                    else
                        return("exon")
                },
                exon = {
                    if (summarized)
                        return("summarized_exon")
                    else
                        return("exon")
                },
                utr = {
                    return("utr")
                }
            )
        }
    )
}

importCustomAnnotation <- function(gtfFile,metadata,
    level=c("gene","transcript","exon"),type=c("gene","exon","utr")) {
    level <- tolower(level[1])
    type <- tolower(type[1])
    # Check metadata
    if (is.null(metadata$organism)) {
        tmpOrg <- paste("species",format(Sys.Date(),"%Y%m%d"),sep="_")
        warning("An organism name must be provided with metadata for ",
            "reporting purposes! Using ",tmpOrg,immediate.=TRUE)
        metadata$organism <- tmpOrg
    }

    if (is.null(metadata$source)) {
        tmpSource <- paste("source",format(Sys.Date(),"%Y%m%d"),sep="_")
        warning("A source should be provided with metadata for reporting ",
            "purposes ! Using ",tmpSource,immediate.=TRUE)
        metadata$source <- tmpSource
    }
    if (is.null(metadata$version)) {
        tmpVer <- paste("version",format(Sys.Date(),"%Y%m%d"),sep="_")
        warning("A version should be provided with metadata for reporting ",
            "purposes! Using ",tmpVer,immediate.=TRUE)
        metadata$version <- tmpVer
    }
    if (is.null(metadata$chromInfo)) {
        warning("Chromosomal lengths should be provided with metadata! ",
            "Only chromosome names will be available... ",immediate.=TRUE)
        metadata$chromInfo <- NULL
    }
    else {
        str <- metadata$chromInfo
        if (is.character(str) && file.exists(str)) {
            out <- tryCatch(open(Rsamtools::BamFile(str)),error=function(e) e)
            if (inherits(out,"error")) # Not a BAM file, try to read.delim
                metadata$chromInfo <- read.delim(str,row.names=1)
            else
                metadata$chromInfo <- .chromInfoFromBAM(str)
        }
        if (!is.data.frame(metadata$chromInfo))
            stop("metadata$chromInfo must be a data frame!")
    }
    
    # For display meta information
    s <- metadata$source
    o <- metadata$organism
    v <- metadata$version
    
    # Parse the GTF file... If something wrong, it will crash here
    parsed <- parseCustomGtf(gtfFile)
    
    switch(level,
        gene = {
            switch(type,
                gene = {
                    message("Retrieving gene annotation for ",o," from ",s,
                        " version ",v," from ",gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="gene",
                        type="gene")
                },
                exon = {
                    message("Retrieving summarized exon annotation for ",o,
                        " from ",s," version ",v," from ",gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="gene",
                        type="exon",summarized=TRUE)
                },
                utr = {
                    message("Retrieving summarized 3' UTR annotation per gene ",
                        "for ",o," from ",s," version ",v," from ",gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="gene",
                        type="utr",summarized=TRUE)
                }
            )
        },
        transcript = {
            switch(type,
                gene = {
                    message("Retrieving transcript annotation for ",o," from ",
                        s," version ",v," from ",gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="transcript",
                        type="gene")
                },
                exon = {
                    # Stub
                    # TODO: "summarized_exon_by_transcript.rda"
                },
                utr = {
                    message("Retrieving summarized 3' UTR annotation per ",
                        "transcript for ",o," from ",s," version ",v," from ",
                        gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="transcript",
                        type="utr",summarized=TRUE)
                }
            )
        },
        exon = {
            switch(type,
                exon = {
                    message("Retrieving exon annotation for ",o," from ",s,
                        " version ",v," from ",gtfFile)
                    annGr <- annotationFromCustomGtf(parsed,level="exon",
                        type="exon")
                }
            )
        }
    )
    
    return(annGr)
}

.validateDbCon <- function(obj) {
    # obj can be an already opened connection or the db file or missing. In the
    # latter case, the function looks at the package filesystem location
    if (is.null(obj)) {
        db <- file.path(system.file(package="metaseqR"),"annotation.sqlite")
        drv <- dbDriver("SQLite")
        con <- tryCatch(dbConnect(drv,dbname=db),error=function(e) {
            message("Caught error: ",e)
            stop("Have you constructed the metaseqR annotation database?")
        },finally="")
    }
    if (file.exists(obj)) {
        drv <- dbDriver("SQLite")
        con <- tryCatch(dbConnect(drv,dbname=obj),error=function(e) {
            message("Caught error: ",e)
            stop("Is obj an SQLite database?")
        },finally="")
    }
    else if (is(obj,"SQLiteConnection"))
        con <- obj
    return(con)
}

getInstalledAnnotations <- function(obj=NULL) {
    con <- .validateDbCon(obj)
    content <- .browseContent(con)
    dbDisconnect(con)   
    return(content[,-1])
}

getUserAnnotations <- function(obj=NULL) {
    con <- .validateDbCon(obj)
    userContent <- .browseUserContent(con)
    dbDisconnect(con)
    return(userContent[,-1])
}

correctTranscripts <- function(ann) {
    rownames(ann) <- paste("T",seq_len(nrow(ann)),sep="_")
    len <- ann[,3] - ann[,2]
    len <- len[-which(is.na(len))]
    len[len==0] <- 1
    defUtrLen <- round(2^mean(log2(len)))
    nas <- which(is.na(ann$start))
    annNa <- ann[nas,]
    annNa$start <- annNa$tstart
    annNa$end <- annNa$tend
    tmp <- makeGRangesFromDataFrame(df=annNa)
    tmp <- flank(resize(tmp,width=1,fix="end"),start=FALSE,width=defUtrLen)
    ann[names(tmp),"start"] <- start(tmp)
    ann[names(tmp),"end"] <- end(tmp)
    return(ann)
}

reduceTranscripts <- function(gr) {
    # Get the gene ids to use as split factor
    gene <- unique(as.character(gr$gene_id))
    
    # Get the GRanges metadata to create a map for later
    meta <- elementMetadata(gr)
    if (is.null(gr$gene_name))
        meta$gene_name <- meta$gene_id
    if (is.null(gr$biotype))
        meta$biotype <- rep("gene",nrow(meta))
    
    # There will be a transcript_id
    ir <- which(names(meta) == "transcript_id")
    map <- meta[,-ir]
    
    # Remove duplicates and name the map
    d <- which(duplicated(map))
    if (length(d) > 0)
        map <- map[-d,]
    # In some cases, transcripts are recorded with multiple biotypes...
    if (length(map$gene_id) != length(unique(map$gene_id))) {
        mapp <- map
        ir <- which(names(mapp) == "biotype")
        mapp <- mapp[,-ir]
        d <- which(duplicated(mapp))
        if (length(d) > 0)
            map <- map[-d,]
    }
    rownames(map) <- map$gene_id
    map <- map[gene,]
    
    # Gene names and biotypes for later reconstruction
    gn <- as.character(map$gene_name)
    bt <- as.character(map$biotype)
    
    # Split the initial GRanges to apply rest operations
    grList <- split(gr,gr$gene_id)
    # Re-order for common reference with the map
    grList <- grList[gene]
    # Now, reduce to merge overlaping exons
    grList <- reduce(grList)
    
    # Start the reconstruction by getting new lengths and create names
    lens <- lengths(grList)
    inds <- unlist(lapply(lens,function(j) return(seq_len(j))),use.names=FALSE)
    grNew <- unname(unlist(grList))
    gene_id <- rep(gene,lens)
    transcript_id <- paste(gene_id,"MET",inds,sep="_")
    gene_name <- rep(gn,lens)
    biotype <- rep(bt,lens)
    newMeta <- DataFrame(
        transcript_id=transcript_id,
        gene_id=gene_id,
        gene_name=gene_name,
        biotype=biotype
    )
    mcols(grNew) <- newMeta
    
    # grNew is the GRanges to return. In order to get the activeLength, we split
    # again per gene_id in a temp variable
    tmp <- split(grNew,grNew$gene_id)
    tmp <- tmp[gene]
    len <- vapply(width(tmp),sum,integer(1))
    
    #return(grNew)
    return(list(model=grNew,length=len))
}

reduceTranscriptsUtr <- function(gr) {
    # Get the transcript ids to use ordering element
    trans <- unique(as.character(gr$transcript_id))
    
    # Get the GRanges metadata to create a map for later
    meta <- elementMetadata(gr)
    if (is.null(gr$gene_name))
        meta$gene_name <- meta$gene_id
    if (is.null(gr$biotype))
        meta$biotype <- rep("gene",nrow(meta))
    
    # Make the map and remove duplicates
    map <- meta
    d <- which(duplicated(map))
    if (length(d) > 0)
        map <- map[-d,,drop=FALSE]
    rownames(map) <- map$transcript_id
    map <- map[trans,]
    
    # Gene ids, names and biotypes for later reconstruction
    gi <- as.character(map$gene_id)
    gn <- as.character(map$gene_name)
    bt <- as.character(map$biotype)
    
    # Split the initial GRanges to apply rest operations
    grList <- split(gr,gr$transcript_id)
    # Re-order for common reference with the map
    grList <- grList[trans]
    # Now, reduce to merge overlaping exons
    grList <- reduce(grList)
    
    # Start the reconstruction by getting new lengths and create names
    lens <- lengths(grList)
    inds <- unlist(lapply(lens,function(j) return(seq_len(j))),use.names=FALSE)
    grNew <- unname(unlist(grList))
    transcript_id <- paste(rep(trans,lens),"MEU",inds,sep="_")
    gene_id <- rep(gi,lens)
    gene_name <- rep(gn,lens)
    biotype <- rep(bt,lens)
    newMeta <- DataFrame(
        transcript_id=transcript_id,
        #gene_id=gene_id,
        gene_id=rep(trans,lens),
        gene_name=gene_name,
        biotype=biotype
    )
    mcols(grNew) <- newMeta
    
    # grNew is the GRanges to return. In order to get the activeLength, we split
    # again per gene_id in a temp variable
    tmp <- split(grNew,grNew$gene_id)
    tmp <- tmp[trans]
    len <- vapply(width(tmp),sum,integer(1))
    
    #return(grNew)
    return(list(model=grNew,length=len))
}

reduceExons <- function(gr) {
    # Get the gene ids to use as split factor
    gene <- unique(as.character(gr$gene_id))
    
    # Get the GRanges metadata to create a map for later
    meta <- elementMetadata(gr)
    if (is.null(gr$gene_name))
        meta$gene_name <- meta$gene_id
    if (is.null(gr$biotype))
        meta$biotype <- rep("gene",nrow(meta))
    
    # There will be an exon_id
    ir <- which(names(meta) == "exon_id")
    map <- meta[,-ir]
    
    # Remove duplicates and name the map
    d <- which(duplicated(map))
    if (length(d) > 0)
        map <- map[-d,]
    rownames(map) <- map$gene_id
    map <- map[gene,]
    
    # Gene names and biotypes for later reconstruction
    gn <- as.character(map$gene_name)
    bt <- as.character(map$biotype)
    
    # Split the initial GRanges to apply rest operations
    grList <- split(gr,gr$gene_id)
    # Re-order for common reference with the map
    grList <- grList[gene]
    # Now, reduce to merge overlaping exons
    grList <- reduce(grList)
    
    # Start the reconstruction by getting new lengths and create names
    lens <- lengths(grList)
    inds <- unlist(lapply(lens,function(j) return(seq_len(j))),use.names=FALSE)
    grNew <- unname(unlist(grList))
    gene_id <- rep(gene,lens)
    exon_id <- paste(gene_id,"MEX",inds,sep="_")
    gene_name <- rep(gn,lens)
    biotype <- rep(bt,lens)
    newMeta <- DataFrame(
        exon_id=exon_id,
        gene_id=gene_id,
        gene_name=gene_name,
        biotype=biotype
    )
    mcols(grNew) <- newMeta
    
    # grNew is the GRanges to return. In order to get the activeLength, we split
    # again per gene_id in a temp variable
    tmp <- split(grNew,grNew$gene_id)
    tmp <- tmp[gene]
    len <- vapply(width(tmp),sum,integer(1))
    
    return(list(model=grNew,length=len))
}

reduceTranscriptsExons <- function(gr) {
    # Get the transcript ids to use ordering element
    trans <- unique(as.character(gr$transcript_id))
    
    # Get the GRanges metadata to create a map for later
    meta <- elementMetadata(gr)
    if (is.null(gr$gene_name))
        meta$gene_name <- meta$gene_id
    if (is.null(gr$biotype))
        meta$biotype <- rep("gene",nrow(meta))
        
    # There will be an exon_id
    ir <- which(names(meta) == "exon_id")
    map <- meta[,-ir]
    
    # Make the map and remove duplicates
    d <- which(duplicated(map))
    if (length(d) > 0)
        map <- map[-d,,drop=FALSE]
    rownames(map) <- map$transcript_id
    map <- map[trans,]
    
    # Gene ids, names and biotypes for later reconstruction
    #gi <- as.character(map$gene_id)
    gn <- as.character(map$gene_name)
    bt <- as.character(map$biotype)
    
    # Split the initial GRanges to apply rest operations
    grList <- split(gr,gr$transcript_id)
    # Re-order for common reference with the map
    grList <- grList[trans]
    # Now, reduce to merge overlaping exons
    grList <- reduce(grList)
    
    # Start the reconstruction by getting new lengths and create names
    lens <- lengths(grList)
    inds <- unlist(lapply(lens,function(j) return(seq_len(j))),use.names=FALSE)
    grNew <- unname(unlist(grList))
    exon_id <- paste(rep(trans,lens),"MTE",inds,sep="_")
    #transcript_id <- rep(gi,lens)
    #gene_id <- rep(gi,lens)
    gene_name <- rep(gn,lens)
    biotype <- rep(bt,lens)
    newMeta <- DataFrame(
        exon_id=exon_id,
        #gene_id=gene_id,
        transcript_id=rep(trans,lens),
        gene_name=gene_name,
        biotype=biotype
    )
    mcols(grNew) <- newMeta
    
    # grNew is the GRanges to return. In order to get the activeLength, we split
    # again per gene_id in a temp variable
    tmp <- split(grNew,grNew$transcript_id)
    tmp <- tmp[trans]
    len <- vapply(width(tmp),sum,integer(1))
    
    #return(grNew)
    return(list(model=grNew,length=len))
}

getAnnotation <- function(org,type,refdb="ensembl",ver=NULL,rc=NULL) {
    org <- tolower(org)
    switch(refdb,
        ensembl = { return(getEnsemblAnnotation(org,type,ver)) },
        ucsc = { return(getUcscAnnotation(org,type,refdb,rc=rc)) },
        refseq = { return(getUcscAnnotation(org,type,refdb,rc=rc)) }
    )
}

getEnsemblAnnotation <- function(org,type,ver=NULL) {
    if (org=="tair10")
        dat <- "plants_mart"
    else
        dat <- "ENSEMBL_MART_ENSEMBL"
    host <- getHost(org,ver)
    message("Using Ensembl host ",host)
    mart <- useMart(biomart=dat,host=host,dataset=getDataset(org))

    chrsExp <- paste("^",getValidChrs(org),"$",sep="",collapse="|")
    if (type=="gene") {
        bm <- tryCatch(
            getBM(attributes=.getGeneAttributes(org),mart=mart),
            error=function(e) {
                message("Caught error: ",e)
                message("This error is most probably related to httr package ",
                    "internals! Using fallback solution...")
                .myGetBM(attributes=.getGeneAttributes(org),mart=mart)
            },
            finally=""
        )
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$start_position,
            end=bm$end_position,
            gene_id=bm$ensembl_gene_id,
            gc_content=if (org %in% 
                c("hg18","hg19","mm9","rn5","dm3","danrer7")) 
                bm$percentage_gc_content else bm$percentage_gene_gc_content,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","hg19","mm9")) bm$external_gene_id 
                else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- ann$gene_id
    }
    else if (type=="transcript") {
        bm <- tryCatch(
            getBM(attributes=.getTranscriptAttributes(org),mart=mart),
            error=function(e) {
                message("Caught error: ",e)
                message("This error is most probably related to httr package ",
                    "internals! Using fallback solution...")
                .myGetBM(attributes=.getTranscriptAttributes(org),mart=mart)
            },
            finally=""
        )
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$transcript_start,
            end=bm$transcript_end,
            transcript_id=bm$ensembl_transcript_id,
            gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","hg19","mm9")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- as.character(ann$transcript_id)
    }
    else if (type=="utr") {
        bm <- tryCatch(
            getBM(attributes=.getTranscriptUtrAttributes(org),mart=mart),
            error=function(e) {
                message("Caught error: ",e)
                message("This error is most probably related to httr package ",
                    "internals! Using fallback solution...")
                .myGetBM(attributes=.getTranscriptUtrAttributes(org),mart=mart)
            },
            finally=""
        )
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$`3_utr_start`,
            end=bm$`3_utr_end`,
            tstart=bm$transcript_start,
            tend=bm$transcript_end,
            transcript_id=bm$ensembl_transcript_id,
            gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","hg19","mm9","tair10")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        ann <- correctTranscripts(ann)
        ann <- ann[,c("chromosome","start","end","transcript_id","gene_id",
            "strand","gene_name","biotype")]
    }
    else if (type=="exon") {
        bm <- tryCatch(
            getBM(attributes=.getExonAttributes(org),mart=mart),
            error=function(e) {
                message("Caught error: ",e)
                message("This error is most probably related to httr package ",
                    "internals! Using fallback solution...")
                .myGetBM(attributes=.getExonAttributes(org),mart=mart)
            },
            finally=""
        )
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$exon_chrom_start,
            end=bm$exon_chrom_end,
            exon_id=bm$ensembl_exon_id,
            gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","hg19","mm9")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- ann$exon_id
    }
    else if (type=="transexon") {
        bm <- tryCatch(
            getBM(attributes=.getTranscriptExonAttributes(org),mart=mart),
            error=function(e) {
                message("Caught error: ",e)
                message("This error is most probably related to httr package ",
                    "internals! Using fallback solution...")
                .myGetBM(attributes=.getTranscriptExonAttributes(org),mart=mart)
            },
            finally=""
        )
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$exon_chrom_start,
            end=bm$exon_chrom_end,
            exon_id=bm$ensembl_exon_id,
            transcript_id=bm$ensembl_transcript_id,
            #gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","hg19","mm9")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
    }
    ann <- ann[order(ann$chromosome,ann$start),]
    ann <- ann[grep(chrsExp,ann$chromosome),]
    ann$chromosome <- as.character(ann$chromosome)
    
    return(ann)
}

getUcscAnnotation <- function(org,type,refdb="ucsc",chunkSize=500,rc=NULL) {
    if (!requireNamespace("RMySQL")) {
        rmysqlPresent <- FALSE
        warning("R package RMySQL is not present! Annotation will be ",
            "retrieved by downloading temporary files from UCSC and the usage
            of a temporary SQLite database...",immediate.=TRUE)
    }
    else
        rmysqlPresent <- TRUE
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to use annotation from UCSC!")
    
    if (org=="tair10") {
        warnwrap("Arabidopsis thaliana genome is not supported by UCSC Genome ",
            "Browser database! Switching to Ensembl...")
        return(getEnsemblAnnotation("tair10",type))
    }
    
    validChrs <- getValidChrs(org)
    chrsExp <- paste("^",paste(validChrs,collapse="$|^"),"$",sep="")

    dbOrg <- getUcscOrganism(org)
    if (rmysqlPresent) {
        # The UTR case is handled later
        if (type != "utr") {
            dbCreds <- .getUcscCredentials()
            drv <- dbDriver("MySQL")
            con <- dbConnect(drv,user=dbCreds[2],password=NULL,dbname=dbOrg,
                host=dbCreds[1])
            type2 <- ifelse(type=="transexon","exon",type)
            query <- getUcscQuery(org,type2,refdb)
            rawAnn <- dbGetQuery(con,query)
            dbDisconnect(con)
        }
    }
    else {
        # This should return the same data frame as the db query
        #if (type == "transcript")
        #    tmpSqlite <- getUcscDbl(org,"gene",refdb)
        #else if (type %in% c("gene","exon"))
        #    tmpSqlite <- getUcscDbl(org,type,refdb)
        
        tmpSqlite <- getUcscDbl(org,refdb)
        
        if (type != "utr") { # Otherwise diect download is used
            drv <- dbDriver("SQLite")
            con <- dbConnect(drv,dbname=tmpSqlite)
            query <- getUcscQuery(org,type,refdb)
            rawAnn <- dbGetQuery(con,query)
            dbDisconnect(con)
        }
    }
    if (type=="gene") {
        ann <- rawAnn
        ann <- ann[grep(chrsExp,ann$chromosome,perl=TRUE),]
        ann$chromosome <- as.character(ann$chromosome)
        rownames(ann) <- ann$gene_id
        #rownames(ann) <- ann$transcript_id
        #tmpAnn <- rawAnn
        #tmpAnn <- tmpAnn[grep(chrsExp,tmpAnn$chromosome,perl=TRUE),]
        #tmpAnn$chromosome <- as.character(tmpAnn$chromosome)
        #rownames(tmpAnn) <- tmpAnn$transcript_id
        ## Split the UCSC transcripts per gene name
        ##tmp <- split(tmpAnn,tmpAnn$gene_name)
        # Retrieve the longest transcript (as per Ensembl convention)
        #ann <- do.call("rbind",cmclapply(tmp,function(x) {
        #    size <- x$end - x$start
        #    selected <- which(size == max(size))
        #    return(x[selected[1],,drop=FALSE])
        #},rc=rc))
        #names(ann)[4] <- "gene_id"
        gcContent <- getGcContent(ann,org)
        ann$gc_content <- gcContent
        # rownames are leftover from splitting above
        rownames(ann) <- ann$gene_id
    }
    else if (type=="transcript") {
        ann <- rawAnn
        ann <- ann[grep(chrsExp,ann$chromosome,perl=TRUE),]
        d <- which(duplicated(ann))
        if (length(d) > 0)
            ann <- ann[-d,]
        ann$chromosome <- as.character(ann$chromosome)
        ann$transcript_id <- as.character(ann$transcript_id)
        ann$gene_id <- as.character(ann$transcript_id)
        # There are still duplicated transcript ids
        d <- which(duplicated(ann$transcript_id))
        iter <- 1
        while(length(d) > 0) {
            ann$transcript_id[d] <- paste(ann$transcript_id[d],iter,sep="_")
            iter <- iter + 1
            d <- which(duplicated(ann$transcript_id))
        }
        rownames(ann) <- ann$transcript_id
        ann <- ann[,c(1,2,3,4,8,5:7)]
    }
    else if (type=="exon" || type=="transexon") {
        rawAnn <- rawAnn[grep(chrsExp,rawAnn$chromosome,perl=TRUE),]
        exList <- cmclapply(as.list(seq_len(nrow(rawAnn))),function(x,d) {
            r <- d[x,,drop=FALSE]
            starts <- as.numeric(strsplit(r[,"start"],",")[[1]])
            ends <- as.numeric(strsplit(r[,"end"],",")[[1]])
            nexons <- length(starts)
            ret <- data.frame(
                rep(r[,"chromosome"],nexons),
                starts,ends,
                paste(r[,"exon_id"],"_e",seq_len(nexons),sep=""),
                rep(r[,"strand"],nexons),
                rep(r[,"gene_id"],nexons),
                rep(r[,"gene_name"],nexons),
                rep(r[,"biotype"],nexons)
            )
            names(ret) <- names(r)
            rownames(ret) <- ret$exon_id
            return(ret)
        },rawAnn,rc=rc)
        
        # For some reason rbind takes ages for large datasets... We have to 
        # split in chunks of 1000
        N <- length(exList)
        mo <- N%%chunkSize
        if (mo == 0) {
            fl <- N/chunkSize
            fac <- factor(rep(seq_len(fl),each=chunkSize))
        }
        else {
            fl <- (N-mo)/chunkSize
            fac <- factor(c(rep(seq_len(fl),each=chunkSize),rep(fl,mo)))
        }
        exChunkedList <- split(exList,fac)
        # Merge the chunks
        tmp <- cmclapply(names(exChunkedList),function(n,X) {
            message("Binding chunk ",n,"...")
            return(do.call("rbind",X[[n]]))
        },exChunkedList,rc=rc)
        # Final merge
        message("Binding all chunks...")
        tmpAnn <- do.call("rbind",tmp)
        
        ann <- data.frame(
            chromosome=as.character(tmpAnn$chromosome),
            start=tmpAnn$start,
            end=tmpAnn$end,
            exon_id=as.character(tmpAnn$exon_id),
            gene_id=as.character(tmpAnn$gene_id),
            strand=as.character(tmpAnn$strand),
            gene_name=as.character(tmpAnn$gene_name),
            biotype=as.character(tmpAnn$biotype)
        )
        rownames(ann) <- ann$exon_id
        
        if (type == "transexon")
            names(ann)[5] <- "transcript_id"
    }
    else if (type=="utr") {
        # We are already supposed to have the necessary elements to construct
        # the required data frame. Should be something like
        # 1. Read the GTF as TxDb
        # 2. Import the GTF with rtracklayer
        # 3. Call the 3UTR function of TxDb on (1)
        # 4. Construct a map to add gene_name and other things (essentially
        #    the biotype if we make it, if not we have to connect it from the
        #    respective gene annotation during the pipeline execution)
        # 5. Add the mapped elements to the GRanges/data.frame
        # 6. Return a data.frame
        #
        # All the process is streamlined in getUcscUtr function
        
        preAnn <- getUcscUtr(org,refdb)
        preAnn <- as.data.frame(preAnn)
        preAnn <- preAnn[,c(1,2,3,6,7,5,8,9)]
        preAnn <- preAnn[grep(chrsExp,preAnn$seqnames,perl=TRUE),]
        names(preAnn) <- c("chromosome","start","end","transcript_id","gene_id",
            "strand","gene_name","biotype")
        # preAnn now has exon names as rownames... OK...
        #rownames(preAnn) <- paste("T",seq_len(nrow(preAnn)),sep="_")
        ann <- preAnn
    }
    
    ann <- ann[order(ann$chromosome,ann$start),]
    return(ann)
}

getUcscUtr <- function(org,refdb="ucsc") {
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    checkTextArgs("refdb",refdb,getSupportedUcscDbs())
    
    if (!requireNamespace("GenomicFeatures"))
        stop("Bioconductor package GenomicFeatures is required!")
    
    httpBase <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/",
        getUcscOrganism(org),"/database/",sep="")
    tableUtr <- getUcscTableNameUtr(org,refdb) # Need one table
    message("  Retrieving table ",tableUtr," for 3'UTR annotation generation ",
        "from ",refdb," for ",org)
    download.file(paste(httpBase,tableUtr,sep=""),file.path(tempdir(),
        paste(tableUtr,".txt.gz",sep="")),quiet=TRUE)
    
    # Do the conversion stuff. AS there is no easy way to check if genePredToGtf
    # exists in the system, we should download it on the fly (once for the 
    # session). If no Linux machine, then problem.
    genePredToGtf <- file.path(tempdir(),"genePredToGtf")
    if (!file.exists(file.path(tempdir(),"genePredToGtf"))) {
        message("  Retrieving genePredToGtf tool")
        download.file(
        "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf",
            genePredToGtf,quiet=TRUE
        )
        system(paste("chmod 775",genePredToGtf))
    }
    
    # Then do the conversion... No need for windows case as Kent tools do not
    # work in Windows anyway
    gtfFile <- file.path(tempdir(),paste(tableUtr,".gtf",sep=""))
    message("  Converting ",tableUtr," to GTF")
    tmpName <- file.path(tempdir(),paste(format(Sys.time(),"%Y%m%d%H%M%S"),
        "tgtf",sep="."))
    commandUcsc <- paste(
        "zcat ",file.path(tempdir(),paste(tableUtr,".txt.gz",sep=""))," | ",
        "cut -f1-10 - | ",genePredToGtf," file stdin ",tmpName," -source=",
        tableUtr," -utr && grep -vP '\t\\.\t\\.\t' ",tmpName," > ",gtfFile,
        sep=""
    )
    commandRefSeq <- paste(
        "zcat ",file.path(tempdir(),paste(tableUtr,".txt.gz",sep=""))," | ",
        "cut -f2- | ",genePredToGtf," file stdin ",tmpName," -source=",
        tableUtr," -utr && grep -vP '\t\\.\t\\.\t' ",tmpName," > ",gtfFile,
        sep=""
    )
    command <- commandRefSeq
    # There is an exception for organisms that do not exist in UCSC databases
    # so we must use RefSeq
    ucscUnsup <- c("rn5","rn6","dm3","dm6","danrer7","danrer10","danrer11",
        "pantro4","pantro5","susscr3","susscr11","equcab2")
    if (refdb == "ucsc" && !(org %in% ucscUnsup))
        command <- commandUcsc
    
    # Run the command and process the data
    message("Executing: ",command)
    system(command)
    parsed <- parseCustomGtf(gtfFile)
    return(.makeGeneUtrFromTxDb(parsed$txdb,parsed$map,asdf=FALSE))
}

getGcContent <- function(ann,org) {
    if (missing(ann))
        stop("A valid annotation data frame must be provided in order to ",
            "retrieve GC-content.")
    org <- tolower(org[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    # Convert annotation to GRanges
    message("Converting annotation to GenomicRanges object...")
    annGr <- tryCatch(GRanges(ann),error=function(e) {
        if (packageVersion("GenomicRanges")<1.14)
            return(GRanges(
                seqnames=Rle(ann[,1]),
                ranges=IRanges(start=ann[,2],end=ann[,3]),
                strand=Rle(ann[,6]),
                name=as.character(ann[,4])
            ))
        else
            return(makeGRangesFromDataFrame(
                df=ann,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome"
            ))
    },finaly="")
    
    bsg <- loadBsGenome(org)
    if (is(bsg,"BSgenome")) {
        message("Getting DNA sequences...")
        seqs <- tryCatch(getSeq(bsg,names=annGr),error=function(e) {
            warning("Caught error ",e,immediate.=TRUE)
            message("Cannot get ",org," sequences! Returning NA GC content...")
            return(NA)
        },finally="")
        if (!is.na(seqs[1])) {
            message("Getting GC content...")
            freqMatrix <- alphabetFrequency(seqs,as.prob=TRUE,baseOnly=TRUE)
            gcContent <- apply(freqMatrix,1,function(x) round(100*sum(x[2:3]),
                digits=2))
        }
        else
            gcContent <- rep(NA,nrow(ann))
    }
    else
        gcContent <- rep(NA,nrow(ann))
    names(gcContent) <- as.character(ann[,4])
    return(gcContent)
}

getSeqInfo <- function(org,asSeqinfo=FALSE) {
    sf <- tryCatch(.chromInfoWrapperGID(getUcscOrganism(org)),
        error=function(e) {
            message("GenomeInfoDb chrom info fetch mechanism failed with the ",
                "following error: ")
            message(e)
            message("")
            message("Trying a direct download...")
            getChromInfo(getUcscOrganism(org))
        },finally="")
    rownames(sf) <- as.character(sf[,1])
    sf <- sf[getValidChrs(org),]
    if (asSeqinfo)
        return(Seqinfo(seqnames=sf[,1],seqlengths=sf[,2],
            isCircular=sf[,3],genome=getUcscOrganism(org)))
    else
        return(data.frame(chromosome=sf[,1],length=as.integer(sf[,2])))
}

getUcscOrganism <- function(org) {
    switch(org,
        hg18 = { return("hg18") },
        hg19 = { return("hg19") },
        hg38 = { return("hg38") },
        mm9 = { return("mm9") },
        mm10 = { return("mm10") },
        rn5 = { return("rn5") },
        rn6 = { return("rn6") },
        dm3 = { return("dm3") },
        dm6 = { return("dm6") },
        danrer7 = { return("danRer7") },
        danrer10 = { return("danRer10") },
        danrer11 = { return("danRer11") },
        pantro4 = { return("panTro4") },
        pantro5 = { return("panTro5") },
        susscr3 = { return("susScr3") },
        susscr11 = { return("susScr11") },
        equcab2 = { return("equCab2") },
        tair10 = { return("TAIR10") }
    )
}

getBsOrganism <- function(org) {
    switch(org,
        hg18 = {
            return("BSgenome.Hsapiens.UCSC.hg18")
        },
        hg19 = {
            return("BSgenome.Hsapiens.UCSC.hg19")
        },
        hg38 = {
            return("BSgenome.Hsapiens.UCSC.hg38")
        },
        mm9 = {
            return("BSgenome.Mmusculus.UCSC.mm9")
        },
        mm10 = {
            return("BSgenome.Mmusculus.UCSC.mm10")
        },
        rn5 = {
            return("BSgenome.Rnorvegicus.UCSC.rn5")
        },
        rn6 = {
            return("BSgenome.Rnorvegicus.UCSC.rn6")
        },
        dm3 = {
            return("BSgenome.Dmelanogaster.UCSC.dm3")
        },
        dm6 = {
            return("BSgenome.Dmelanogaster.UCSC.dm6")
        },
        danrer7 = {
            return("BSgenome.Drerio.UCSC.danRer7")
        },
        danrer10 = {
            return("BSgenome.Drerio.UCSC.danRer10")
        },
        danrer11 = {
            warning("danRer11 is not supported by BSgenome! Please use Ensembl",
                " as annotation source if GC content is important.",
                immediate.=TRUE)
            return(NA)
            #return("BSgenome.Drerio.UCSC.danRer11")
        },
        pantro4 = {
            warning("panTro4 is not supported by BSgenome! Please use Ensembl ",
                "as annotation source if GC content is important.",
                immediate.=TRUE)
            return(NA)
        },
        pantro5 = {
            return("BSgenome.Ptroglodytes.UCSC.panTro5")
        },
        susscr3 = {
            return("BSgenome.Sscrofa.UCSC.susScr3")
        },
        susscr11 = {
            warning("susScr11 is not supported by BSgenome! Please use ",
                "Ensembl as annotation source if GC content is important.",
                immediate.=TRUE)
            return(NA)
        },
        equcab2 = {
            warning("equCab2 is not supported by BSgenome! Please use Ensembl ",
                "as annotation source if GC content is important.",
                immediate.=TRUE)
            return(NA)
        },
        tair10 = {
            warning("TAIR10 is not supported by BSgenome! Please use Ensembl ",
                "as annotation source if GC content is important.",
                immediate.=TRUE)
            return(NA)
        }
    )
}

loadBsGenome <- function(org) {
    if (!requireNamespace("BiocManager"))
        stop("The Bioconductor package BiocManager is required to ",
            "proceed!")
    if (!requireNamespace("BSgenome"))
        stop("The Bioconductor package BSgenome is required to proceed!")
    bsOrg <- getBsOrganism(org)
    if (!is.na(bsOrg)) {
        if (bsOrg %in% BSgenome::installed.genomes())
            bsObj <- BSgenome::getBSgenome(getUcscOrganism(org))
        else {
            BiocManager::install(bsOrg,update=FALSE,ask=FALSE)
            bsObj <- BSgenome::getBSgenome(getUcscOrganism(org))
        }
        return(bsObj)
    }
    else
        return(NA)
}

getChromInfo <- function(org,
    goldenPath="http://hgdownload.cse.ucsc.edu/goldenPath/") {
    download.file(paste(goldenPath,org,"/database/chromInfo.txt.gz",sep=""),
        file.path(tempdir(),"chromInfo.txt.gz"),quiet=TRUE)
    chromInfo <- read.delim(file.path(tempdir(),"chromInfo.txt.gz"),
        header=FALSE)
    chromInfo <- chromInfo[,1,2]
    chromInfo[,1] <- as.character(chromInfo[,1])
    chromInfo$V3 <- rep(FALSE,nrow(chromInfo))
    m <- grep("M",chromInfo[,1])
    if (length(m) > 0)
        chromInfo$V3[m] <- TRUE
    return(chromInfo)
}

getHost <- function(org,ver=NULL) {
    if (!requireNamespace("biomaRt"))
        stop("The Bioconductor package biomaRt is required to proceed!")
    
    org <- tolower(org[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    if (!is.null(ver) 
        && (!is.numeric(ver) || is.na(suppressWarnings(as.numeric(ver)))))
        stop("ver must be numeric or coercible to numeric if not NULL!")
        
    if (org == "tair10")
        return("plants.ensembl.org")
    
    aver <- getUcscToEnsembl(org)
    if (!is.null(ver) && !(ver %in% aver)) {
        warning("Version ",ver," not available/existing for ",org,"! Will ",
            "use the latest available version...",immediate.=TRUE)
        ver <- NULL
    }
    
    if (is.null(ver)) {
        u2e <- ucscToEnsembl()
        vers <- u2e[[org]]
        ver <- vers[length(vers)]
    }
    
    ea <- biomaRt::listEnsemblArchives()
    i <- grep(as.character(ver),ea[,"version"])
    if (length(i) > 0) {
        if (ea[i,"current_release"] == "*")
            return("http://www.ensembl.org")
        else
            return(ea[i,"url"])
    }
    else {
        warning("Version ",ver," not found in biomaRt archives for ",org,"! ",
            "Will use the latest available version for ",org,"...",
            immediate.=TRUE)
        u2e <- ucscToEnsembl()
        vers <- u2e[[org]]
        for (v in rev(vers)) {
            newi <- grep(as.character(v),ea[,"version"])
            if (length(newi) > 0)
                return(ea[newi,"url"])
        }
    }
    # If everything fails...
    return(NULL)
}

getUcscToEnsembl <- function(org) {
    u2e <- ucscToEnsembl()
    return(u2e[[org]])
}

checkUcscToEnsembl <- function(org,ver) {
    u2e <- getUcscToEnsembl()
    return(ver %in% u2e[[org]])
}

ucscToEnsembl <- function() {
    return(list(
        hg18=67,
        hg19=74:75,
        hg38=76:100,
        mm9=67,
        mm10=74:100,
        rn5=74:79,
        rn6=80:100,
        dm3=c(67,74:78),
        dm6=79:100,
        danrer7=c(67,74:79),
        danrer10=80:91,
        danrer11=92:100,
        pantro4=c(67,74:90),
        pantro5=91:100,
        #pantro6=,
        susscr3=c(67,74:89),
        susscr11=90:100,
        equcab2=c(67,74:100)
    ))
}

getAltHost <- function(org) {
    .Deprecated("getHost")
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("uswest.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("uswest.ensembl.org") },
        rn5 = { return("uswest.ensembl.org") },
        dm3 = { return("uswest.ensembl.org") },
        dm6 = { return("uswest.ensembl.org") },
        danrer7 = { return("uswest.ensembl.org") },
        danrer10 = { return("uswest.ensembl.org") },
        danrer11 = { return("uswest.ensembl.org") },
        pantro4 = { return("uswest.ensembl.org") },
        pantro5 = { return("uswest.ensembl.org") },
        #pantro6 = { return("uswest.ensembl.org") },
        susscr3 = { return("uswest.ensembl.org") },
        susscr11 = { return("www.ensembl.org") }
    )
}

getDataset <- function(org) {
    switch(org,
        hg18 = { return("hsapiens_gene_ensembl") },
        hg19 = { return("hsapiens_gene_ensembl") },
        hg38 = { return("hsapiens_gene_ensembl") },
        mm9 = { return("mmusculus_gene_ensembl") },
        mm10 = { return("mmusculus_gene_ensembl") },
        rn5 = { return("rnorvegicus_gene_ensembl") },
        rn6 = { return("rnorvegicus_gene_ensembl") },
        dm3 = { return("dmelanogaster_gene_ensembl") },
        dm6 = { return("dmelanogaster_gene_ensembl") },
        danrer7 = { return("drerio_gene_ensembl") },
        danrer10 = { return("drerio_gene_ensembl") },
        danrer11 = { return("drerio_gene_ensembl") },
        pantro4 = { return("ptroglodytes_gene_ensembl") },
        pantro5 = { return("ptroglodytes_gene_ensembl") },
        #pantro6 = { return("ptroglodytes_gene_ensembl") },
        susscr3 = { return("sscrofa_gene_ensembl") },
        susscr11 = { return("sscrofa_gene_ensembl") },
        equcab2 = { return("ecaballus_gene_ensembl") },
        tair10 = { return("athaliana_eg_gene") }
    )
}

getValidChrs <- function(org) {
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        rn6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        dm6 = {
            return(c(
                "chr2L","chr2R","chr3L","chr3R","chr4","chrX"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        pantro5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        pantro6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        susscr11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        equcab2 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr26","chr27","chr28","chr29","chr3","chr30",
                "chr31","chr4","chr5","chr6","chr7","chr8","chr9","chrX"#,"chrY"
            ))
        },
        tair10 = {
            return(c(
                "chr1","chr2","chr3","chr4","chr5"
            ))
        }
    )
}

getValidChrsWithMit <- function(org) {
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY","chrM"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY","chrM"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY","chrM"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY","chrM"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY","chrM"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrM"
            ))
        },
        rn6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrM"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        dm6 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY",
                "chrM"
            ))
        },
        pantro5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY",
                "chrM"
            ))
        },
        pantro6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY",
                "chrM"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY","chrM"
            ))
        },
        susscr11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY","chrM"
            ))
        },
        equcab2 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr26","chr27","chr28","chr29","chr3","chr30",
                "chr31","chr4","chr5","chr6","chr7","chr8","chr9","chrX",
                #"chrY",
                "chrM"
            ))
        },
        tair10 = {
            return(c(
                "chr1","chr2","chr3","chr4","chr5"
            ))
        }
    )
}

.getGeneAttributes <- function(org) {
    if (org %in% c("hg18","hg19","mm9","tair10"))
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_id",
            "gene_biotype"
        ))
    else if (org %in% c("rn5","danrer7","dm3"))
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_name",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gene_gc_content",
            "strand",
            "external_gene_name",
            "gene_biotype"
        ))
}

.getTranscriptAttributes <- function(org) {
    if (org %in% c("hg18","hg19","mm9","tair10"))
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

.getTranscriptUtrAttributes <- function(org) {
    if (org %in% c("hg18","hg19","mm9","tair10"))
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "3_utr_start",
            "3_utr_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "3_utr_start",
            "3_utr_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

.getExonAttributes <- function(org) {
    if (org %in% c("hg18","hg19","mm9","tair10"))
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

.getTranscriptExonAttributes <- function(org) {
    if (org %in% c("hg18","hg19","mm9","tair10"))
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_transcript_id",
            #"ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_transcript_id",
            #"ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

getBiotypes <- function(org) {
    if (!(org %in% getSupportedOrganisms()))
        return(NULL)
    switch(org,
        hg18 = {
            return(c("unprocessed_pseudogene","pseudogene","miRNA",
                "retrotransposed","protein_coding","processed_pseudogene",
                "snRNA","snRNA_pseudogene","Mt_tRNA_pseudogene",
                "miRNA_pseudogene","misc_RNA","tRNA_pseudogene","snoRNA",
                "scRNA_pseudogene","rRNA_pseudogene","snoRNA_pseudogene","rRNA",
                "misc_RNA_pseudogene","IG_V_gene","IG_D_gene","IG_J_gene",
                "IG_C_gene","IG_pseudogene","scRNA"))
        },
        hg19 = {
            return(c("pseudogene","lincRNA","protein_coding","antisense",
                "processed_transcript","snRNA","sense_intronic","miRNA",
                "misc_RNA","snoRNA","rRNA","polymorphic_pseudogene",
                "sense_overlapping","3prime_overlapping_ncrna","TR_V_gene",
                "TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene",
                "TR_J_pseudogene","IG_C_gene","IG_C_pseudogene","IG_J_gene",
                "IG_J_pseudogene","IG_D_gene","IG_V_gene","IG_V_pseudogene"))
        },
        hg38 = {
            return(c("protein_coding","polymorphic_pseudogene","lincRNA",
                "unprocessed_pseudogene","processed_pseudogene","antisense",
                "processed_transcript","transcribed_unprocessed_pseudogene",
                "sense_intronic","unitary_pseudogene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","sense_overlapping",
                "transcribed_processed_pseudogene","miRNA","snRNA","misc_RNA",
                "rRNA","snoRNA","IG_J_pseudogene","IG_J_gene","IG_D_gene",
                "3prime_overlapping_ncrna","IG_C_gene","IG_C_pseudogene",
                "pseudogene","TR_V_pseudogene","Mt_tRNA","Mt_rRNA",
                "translated_processed_pseudogene","TR_J_gene","TR_C_gene",
                "TR_D_gene","TR_J_pseudogene","LRG_gene"))
        },
        mm9 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "lincRNA","snoRNA","processed_transcript","misc_RNA","rRNA",
                "sense_overlapping","sense_intronic","polymorphic_pseudogene",
                "non_coding","3prime_overlapping_ncrna","IG_C_gene",
                "IG_J_gene","IG_D_gene","IG_V_gene","ncrna_host"))
        },
        mm10 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "snoRNA","lincRNA","processed_transcript","misc_RNA","rRNA",
                "sense_intronic","sense_overlapping","polymorphic_pseudogene",
                "IG_C_gene","IG_J_gene","IG_D_gene","IG_LV_gene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","TR_V_pseudogene",
                "3prime_overlapping_ncrna"))
        },
        dm3 = {
            return(c("protein_coding","ncRNA","snoRNA","pre_miRNA","pseudogene",
                "snRNA","tRNA","rRNA"))
        },
        dm6 = {
            return(c("protein_coding","ncRNA","snoRNA","pre_miRNA","pseudogene",
                "snRNA","tRNA","rRNA"))
        },
        rn5 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","misc_RNA"))
        },
        rn6 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","misc_RNA"))
        },
        danrer7 = {
            return(c("antisense","protein_coding","miRNA","snoRNA","rRNA",
                "lincRNA","processed_transcript","snRNA","pseudogene",
                "sense_intronic","misc_RNA","polymorphic_pseudogene",
                "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                "non_coding","sense_overlapping"
            ))
        },
        danrer10 = {
            return(c("antisense","protein_coding","miRNA","snoRNA","rRNA",
                "lincRNA","processed_transcript","snRNA","pseudogene",
                "sense_intronic","misc_RNA","polymorphic_pseudogene",
                "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                "non_coding","sense_overlapping"
            ))
        },
        danrer11 = {
            return(c("antisense","protein_coding","miRNA","snoRNA","rRNA",
                "lincRNA","processed_transcript","snRNA","pseudogene",
                "sense_intronic","misc_RNA","polymorphic_pseudogene",
                "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                "non_coding","sense_overlapping"
            ))
        },
        pantro4 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","snRNA","snoRNA","misc_RNA"))
        },
        pantro5 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","snRNA","snoRNA","misc_RNA"))
        },
        pantro6 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","snRNA","snoRNA","misc_RNA"))
        },
        susscr3 = {
            return(c("antisense","protein_coding","lincRNA","pseudogene",
                "processed_transcript","miRNA","rRNA","snRNA","snoRNA",
                "misc_RNA","non_coding","IG_C_gene","IG_J_gene",
                "IG_V_gene","IG_V_pseudogene"))
        },
        equcab2 = {
            return(c("miRNA","misc_RNA","protein_coding","pseudogene","rRNA",
                "processed_pseudogene","snoRNA","snRNA"))
        },
        tair10 = {
            return(c("miRNA","ncRNA","protein_coding","pseudogene","rRNA",
                "snoRNA","snRNA","transposable_element","tRNA"))
        }
    )
}

getSupportedRefDbs <- function() {
    return(c("ensembl","ucsc","refseq"))
}

getSupportedOrganisms <- function() {
    return(c("hg18","hg19","hg38","mm9","mm10","rn5","rn6","dm3","dm6",
        "danrer7","danrer10","danrer11","pantro4","pantro5","susscr3",
        #"pantro6",
        "susscr11","equcab2","tair10"))
}

getSupportedUcscDbs <- function() {
    return(c("ucsc","refseq"))
}

importCustomGtf <- function(gtfFile,level=c("gene","transcript","exon"),
    type=c("gene","exon","utr")) {
    if (!requireNamespace("GenomicFeatures"))
        stop("Bioconductor package GenomicFeatures is required to build ",
            "custom annotation!")
        
    # Some argument checking
    level <- level[1]
    type <- type[1]
    checkTextArgs("level",level,c("gene","transcript","exon"),multiarg=FALSE)
    checkTextArgs("type",type,c("gene","exon","utr"),multiarg=FALSE)
    
    # Import the GTF with rtracklayer to create a map of available metadata
    message("  Importing GTF ",gtfFile," as GTF to make id map")
    desiredColumns <- c("type","gene_id","transcript_id","exon_id",
        "gene_name","gene_biotype")
    gr <- import(gtfFile,format="gtf",colnames=desiredColumns,
        feature.type=.GFF_FEATURE_TYPES)
    grdf <- as.data.frame(gr)
    grdf <- grdf[grdf$type=="exon",]

    # Need to map gene ids to names and biotypes. We need a collapsed 
    # structure exon_id - transcript_id - gene_id - gane_name - biotype
    message("  Making id map")
    map <- .makeIdMap(grdf)
    
    message("  Importing GTF ",gtfFile," as TxDb")
    txdb <- suppressWarnings(makeTxDbFromGFF(gtfFile))
    
    switch(level,
        gene = {
            switch(type,
                gene = {
                    return(.makeGeneGeneFromTxDb(txdb,map))
                },
                exon = {
                    return(.makeGeneExonFromTxDb(txdb,map))
                },
                utr = {
                    return(.makeGeneUtrFromTxDb(txdb,map))
                }
            )
        },
        transcript = {
            switch(type,
                gene = {
                    return(.makeTranscriptGeneFromTxDb(txdb,map))
                },
                exon = {
                    # Stub
                },
                utr = {
                    return(.makeTranscriptUtrFromTxDb(txdb,map))
                }
            )
        },
        exon = {
            switch(type,
                exon = {
                    return(.makeExonExonFromTxDb(txdb,map))
                }
            )
        }
    )
}

# parsed <- parseCustomGtf(gffFile)
parseCustomGtf <- function(gtfFile) {
    if (!requireNamespace("GenomicFeatures"))
        stop("Bioconductor package GenomicFeatures is required to parse ",
            "custom GTF file!")
    
    # Import the GTF with rtracklayer to create a map of available metadata
    message("  Importing GTF ",gtfFile," as GTF to make id map")
    desiredColumns <- c("type","gene_id","transcript_id","exon_id",
        "gene_name","gene_biotype")
    # Let it recognize automatically GFF/GTF
    gr <- import(gtfFile,colnames=desiredColumns,
        feature.type=.GFF_FEATURE_TYPES)
    grdf <- as.data.frame(gr)
    grdf <- grdf[grdf$type=="exon",]

    # Need to map gene ids to names and biotypes. We need a collapsed 
    # structure exon_id - transcript_id - gene_id - gane_name - biotype
    message("  Making id map")
    map <- .makeIdMap(grdf)
    
    message("  Importing GTF ",gtfFile," as TxDb")
    txdb <- suppressWarnings(makeTxDbFromGFF(gtfFile))
    
    return(list(txdb=txdb,map=map))
}

annotationFromCustomGtf <- function(parsed,level=c("gene","transcript","exon"),
    type=c("gene","exon","utr"),summarized=FALSE,asdf=FALSE) {
    # Some argument checking
    if (!is.logical(summarized))
        stop("summarized must be TRUE or FALSE")
    level <- level[1]
    type <- type[1]
    checkTextArgs("level",level,c("gene","transcript","exon"),multiarg=FALSE)
    checkTextArgs("type",type,c("gene","exon","utr"),multiarg=FALSE)
    
    txdb <- parsed$txdb
    map <- parsed$map
    
    switch(level,
        gene = {
            switch(type,
                gene = {
                    return(.makeGeneGeneFromTxDb(txdb,map,asdf))
                },
                exon = {
                    if (summarized)
                        return(.makeSumGeneExonFromTxDb(txdb,map,asdf))
                    else
                        return(.makeGeneExonFromTxDb(txdb,map,asdf))
                },
                utr = {
                    if (summarized)
                        return(.makeSumGeneUtrFromTxDb(txdb,map,asdf))
                    else
                        return(.makeGeneUtrFromTxDb(txdb,map,asdf))
                }
            )
        },
        transcript = {
            switch(type,
                gene = {
                    if (summarized)
                        return(.makeSumTranscriptGeneFromTxDb(txdb,map,asdf))
                    else
                        return(.makeTranscriptGeneFromTxDb(txdb,map,asdf))
                },
                exon = {
                    # Stub
                },
                utr = {
                    if (summarized)
                        return(.makeSumTranscriptUtrFromTxDb(txdb,map,asdf))
                    else
                        return(.makeTranscriptUtrFromTxDb(txdb,map,asdf))
                }
            )
        },
        exon = {
            switch(type,
                exon = {
                    return(.makeExonExonFromTxDb(txdb,map,asdf))
                }
            )
        }
    )
}

# ann <- .makeGeneGeneFromTxDb(txdb,map,FALSE)
.makeGeneGeneFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- genes(txdb)
        # We need to add gc_content, gene_name, biotype from the map
        # Remove the exon_id and transcript_id columns from the map
        # for the gene case
        ir <- c(which(names(map)=="exon_id"),
            which(names(map)=="transcript_id"))
        if (length(ir) > 0)
            smap <- map[,-ir,drop=FALSE]
        dgt <- which(duplicated(smap))
        if (length(dgt) > 0)
            smap <- smap[-dgt,]
        
        # Add metadata from the map
        rownames(smap) <- smap$gene_id
        smap <- smap[names(gr),,drop=FALSE]
        gr$gene_name <- smap$gene_name
        gr$biotype <- smap$biotype
    }
    else {
        gr <- genes(txdb,columns=c("gene_id","tx_type"))
        gr$gene_id <- as.character(gr$gene_id)
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(gr$tx_type)
        gr$tx_type <- NULL
    }
    
    gr$gc_content = rep(50,length(gr))
    ann <- as.data.frame(gr)
    ann <- ann[,c(1:3,6,9,5,7,8)]
    names(ann)[1] <- "chromosome"
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

# ann <- .makeGeneExonFromTxDb(txdb,map,FALSE)
.makeGeneExonFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- exons(txdb,columns="exon_name")
        
        # There are cases where exons are unnamed and just the structure 
        useMap <- TRUE
        if (any(is.na(gr$exon_name)))
            useMap <- FALSE
        
        if (useMap) {
            names(gr) <- gr$exon_name
        
            # We need to add gene_id, gene_name, biotype from the map
            # Remove the transcript_id column from the map for the gene
            # case
            ir <- which(names(map)=="transcript_id")
            if (length(ir) > 0)
                smap <- map[,-ir,drop=FALSE]
            dgt <- which(duplicated(smap))
            if (length(dgt) > 0)
                smap <- smap[-dgt,]
            
            # Add metadata from the map
            rownames(smap) <- smap$exon_id
            smap <- smap[names(gr),,drop=FALSE]
            gr$gene_id <- smap$gene_id
            gr$gene_name <- smap$gene_name
            gr$biotype <- smap$biotype
        }
        else {
            # In order to create the exon annotation, we need to manually 
            # overlap and assign to genes
            ge <- genes(txdb,columns=c("gene_id","tx_type"))
            ov <- findOverlaps(gr,ge)
            exonPos <- queryHits(ov)
            genePos <- subjectHits(ov)
            levs <- unique(genePos)
            
            dup <- which(duplicated(exonPos))
            if (length(dup) > 0) {
                exonPos <- exonPos[-dup]
                genePos <- genePos[-dup]
            }
            
            S <- split(exonPos,factor(genePos,levels=levs))
            gr$gene_id <- gr$gene_name <- rep(names(ge)[levs],lengths(S))
            gr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],
                lengths(S))
            #gr$biotype <- rep("gene",length(gr))
            if (is(gr$gene_id,"CharacterList"))
                gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
            else
                gtmp <- gr$gene_id
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
            gr$exon_name <- names(gr)
        }
    }
    else {
        gr <- exons(txdb,columns=c("exon_name","gene_id","tx_type"))
        if (is(gr$gene_id,"CharacterList"))
            gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
        else
            gtmp <- gr$gene_id
        if (any(is.na(gr$exon_name)))            
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
        else
            names(gr) <- gr$exon_name
        gr$gene_id <- gtmp
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(lapply(gr$tx_type,function(x) x[1]))
        gr$tx_type <- NULL
    }
    
    ann <- as.data.frame(unname(gr))
    ann <- ann[,c(1,2,3,6,7,5,8,9)]
    names(ann)[c(1,4)] <- c("chromosome","exon_id")
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    if (length(unique(ann$exon_id)) == length(ann$exon_id))
        # Should always be TRUE
        rownames(ann) <- as.character(ann$exon_id)
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

# ann <- .makeSumGeneExonFromTxDb(txdb,map,FALSE)
.makeSumGeneExonFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- exons(txdb,columns="exon_name")
        
        useMap <- TRUE
        if (any(is.na(gr$exon_name)))
            useMap <- FALSE
        
        if (useMap) {
            names(gr) <- gr$exon_name
            
            # We need to add gene_id, gene_name, biotype from the map
            # Remove the transcript_id column from the map for the gene
            # case
            ir <- which(names(map)=="transcript_id")
            if (length(ir) > 0)
                smap <- map[,-ir,drop=FALSE]
            dgt <- which(duplicated(smap))
            if (length(dgt) > 0)
                smap <- smap[-dgt,]
            
            # Add metadata from the map
            rownames(smap) <- smap$exon_id
            smap <- smap[names(gr),,drop=FALSE]
            gr$gene_id <- smap$gene_id
            gr$gene_name <- smap$gene_name
            gr$biotype <- smap$biotype
        }
        else {
            # In order to create the exon annotation, we need to manually 
            # overlap and assign to genes
            ge <- genes(txdb,columns=c("gene_id","tx_type"))
            ov <- findOverlaps(gr,ge)
            exonPos <- queryHits(ov)
            genePos <- subjectHits(ov)
            levs <- unique(genePos)
            
            dup <- which(duplicated(exonPos))
            if (length(dup) > 0) {
                exonPos <- exonPos[-dup]
                genePos <- genePos[-dup]
            }
            
            S <- split(exonPos,factor(genePos,levels=levs))
            gr$gene_id <- gr$gene_name <- rep(names(ge)[levs],lengths(S))
            gr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],
                lengths(S))
            if (is(gr$gene_id,"CharacterList"))
                gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
            else
                gtmp <- gr$gene_id
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
            gr$exon_name <- names(gr)
        }
    }
    else {
        gr <- exons(txdb,columns=c("exon_name","gene_id","tx_type"))
        if (is(gr$gene_id,"CharacterList"))
            gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
        else
            gtmp <- gr$gene_id
        if (any(is.na(gr$exon_name)))            
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
        else
            names(gr) <- gr$exon_name
        gr$gene_id <- gtmp
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(lapply(gr$tx_type,function(x) x[1]))
        gr$tx_type <- NULL
    }
    
    ann <- as.data.frame(unname(gr))
    ann <- ann[,c(1,2,3,6,8,5,7,9)]
    names(ann)[c(1,4)] <- c("chromosome","exon_id")
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    message("  summarizing exons per gene for imported GTF")
    annList <- reduceExons(GRanges(ann))
    sexon <- annList$model
    names(sexon) <- as.character(sexon$exon_id)
    activeLength <- annList$length
    names(activeLength) <- unique(sexon$gene_id)
    attr(sexon,"activeLength") <- activeLength
    
    if (asdf) {
        eann <- as.data.frame(sexon)
        eann <- eann[,c(1,2,3,6,8,5,7,9)]
        names(eann)[c(1,4)] <- c("chromosome","exon_id")
        attr(eann,"activeLength") <- activeLength
        return(eann)
    }
    else
        return(sexon)
}

# ann <- .makeGeneUtrFromTxDb(txdb,map,FALSE)
.makeGeneUtrFromTxDb <- function(txdb,map,asdf) {
    utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
    utrGr <- unlist(utrList)
    
    # There may be no sufficient information to create UTRs
    if (length(utrGr) == 0) {
        warning("No UTR information was found in the provided GTF file! Will ",
            "return an empty object...",immediate.=TRUE)
        ann <- data.frame(chromosome=1,start=1,end=1,
            transcript_id=1,gene_id=1,strand=1,gene_name=1,
            biotype=1,row.names=1)
        # Strange but required to be compatible with GRanges
        ann <- ann[-1,]
        if (asdf)
            return(ann)
        else
            return(GRanges(ann))
    }
    
    utrGr$transcript_id <- names(utrGr)
    utrTmp <- as.data.frame(unname(utrGr))
    keep <- c("seqnames","start","end","transcript_id",
        "exon_rank","strand","exon_name")
    utr <- utrTmp[,keep]
    
    useMap <- TRUE
    if (any(is.na(utr$exon_name)) || is.null(map))
        # Impossible to take info from map...
        useMap <- FALSE
    else
        rownames(utr) <- paste(utr$exon_name,utr$transcript_id,sep="_")
        
    if (useMap) {
        rownames(map) <- paste(map$exon_id,map$transcript_id,sep="_")

        # Different case with map here
        smap <- map[rownames(utr),]
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the gene case
        ir <- which(names(smap)=="exon_id")
        if (length(ir) > 0)
            smap <- smap[,-ir,drop=FALSE]
        
        # Add metadata from the map
        utr$gene_id <- smap$gene_id
        utr$gene_name <- smap$gene_name
        utr$biotype <- smap$biotype
    }
    else {
        # Add metadata from the utr frame
        utr$gene_id <- utr$transcript_id
        utr$gene_name <- utr$gene_id
        
        # To get biotypes
        ge <- genes(txdb,columns=c("gene_id","tx_type"))
        ov <- findOverlaps(utrGr,ge)
        utrPos <- queryHits(ov)
        genePos <- subjectHits(ov)
        levs <- unique(genePos)
        
        dup <- which(duplicated(utrPos))
        if (length(dup) > 0) {
            utrPos <- utrPos[-dup]
            genePos <- genePos[-dup]
        }
        
        S <- split(utrPos,factor(genePos,levels=levs))
        utr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],lengths(S))
        #utr$biotype <- rep("gene",nrow(utr))
    }
    
    ann <- utr[,c(1,2,3,4,8,6,9,10)]
    names(ann)[1] <- "chromosome"
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

# ann <- .makeSumGeneUtrFromTxDb(txdb,map,FALSE)
.makeSumGeneUtrFromTxDb <- function(txdb,map,asdf) {
    utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
    utrGr <- unlist(utrList)
    
    if (length(utrGr) == 0) {
        warning("No UTR information was found in the provided GTF file! Will ",
            "return an empty object...",immediate.=TRUE)
        ann <- data.frame(chromosome=1,start=1,end=1,
            transcript_id=1,gene_id=1,strand=1,gene_name=1,
            biotype=1,row.names=1)
        # Strange but required to be compatible with GRanges
        ann <- ann[-1,]
        if (asdf)
            return(ann)
        else
            return(GRanges(ann))
    }
    
    utrGr$transcript_id <- names(utrGr)
    utrTmp <- as.data.frame(unname(utrGr))
    keep <- c("seqnames","start","end","transcript_id",
        "exon_rank","strand","exon_name")
    utr <- utrTmp[,keep]
    
    useMap <- TRUE
    if (any(is.na(utr$exon_name)) || is.null(map))
        useMap <- FALSE
    else
        rownames(utr) <- paste(utr$exon_name,utr$transcript_id,sep="_")
    
    if (useMap) {
        rownames(map) <- paste(map$exon_id,map$transcript_id,sep="_")
    
        # Different case with map here
        smap <- map[rownames(utr),]
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the gene case
        ir <- which(names(smap)=="exon_id")
        if (length(ir) > 0)
            smap <- smap[,-ir,drop=FALSE]
        
        # Add metadata from the map
        utr$gene_id <- smap$gene_id
        utr$gene_name <- smap$gene_name
        utr$biotype <- smap$biotype
    }
    else {
        # Add metadata from the utr frame
        utr$gene_id <- utr$transcript_id
        utr$gene_name <- utr$gene_id
        
        # To get biotypes
        ge <- genes(txdb,columns=c("gene_id","tx_type"))
        ov <- findOverlaps(utrGr,ge)
        utrPos <- queryHits(ov)
        genePos <- subjectHits(ov)
        levs <- unique(genePos)
        
        dup <- which(duplicated(utrPos))
        if (length(dup) > 0) {
            utrPos <- utrPos[-dup]
            genePos <- genePos[-dup]
        }
        
        S <- split(utrPos,factor(genePos,levels=levs))
        utr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],lengths(S))
    }
    ann <- utr[,c(1,2,3,4,8,6,9,10)]
    names(ann)[1] <- "chromosome"
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    message("  summarizing UTRs per gene for imported GTF")
    #s3utr <- reduceTranscripts(GRanges(ann))
    annList <- reduceTranscripts(GRanges(ann))
    s3utr <- annList$model
    names(s3utr) <- as.character(s3utr$transcript_id)
    activeLength <- annList$length
    names(activeLength) <- unique(as.character(s3utr$gene_id))
    
    if (asdf) {
        sann <- as.data.frame(s3utr)
        sann <- sann[,c(1,2,3,6,8,5,7,9)]
        if (useMap)
            names(sann)[c(1,4)] <- c("chromosome","gene_id")
        else {
            sann <- sann[,c(1,2,3,7,5,6,4,8)]
            names(sann)[c(1,4,7)] <- c("chromosome","gene_id","transcript_id")
        }
        attr(sann,"activeLength") <- activeLength
        return(sann)
    }
    else
        return(s3utr)
}

# ann <- .makeTranscriptGeneFromTxDb(txdb,map,FALSE)
.makeTranscriptGeneFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- transcripts(txdb,columns="tx_name")
        names(gr) <- gr$tx_name
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the transcript 
        # case
        ir <- which(names(map)=="exon_id")
        if (length(ir) > 0)
            smap <- map[,-ir,drop=FALSE]
        dgt <- which(duplicated(smap))
        if (length(dgt) > 0)
            smap <- smap[-dgt,]
        
        # Add metadata from the map
        rownames(smap) <- smap$transcript_id
        smap <- smap[names(gr),,drop=FALSE]
        gr$gene_id <- smap$gene_name
        gr$gene_name <- smap$gene_name
        gr$biotype <- smap$biotype
    }
    else {
        gr <- transcripts(txdb,columns=c("tx_name","gene_id","tx_type"))
        if (is(gr$gene_id,"CharacterList"))
            gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
        else
            gtmp <- as.character(gr$gene_id)
        gr$gene_id <- gtmp
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(lapply(gr$tx_type,function(x) x[1]))
        gr$tx_type <- NULL
    }
    
    ann <- as.data.frame(gr)
    ann <- ann[,c(1,2,3,6,7,5,8,9)]
    names(ann)[c(1,4)] <- c("chromosome","transcript_id")
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

# ann <- .makeSumTranscriptGeneFromTxDb(txdb,map,FALSE)
.makeSumTranscriptGeneFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- transcripts(txdb,columns="tx_name")
        names(gr) <- gr$tx_name
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the transcript 
        # case
        ir <- which(names(map)=="exon_id")
        if (length(ir) > 0)
            smap <- map[,-ir,drop=FALSE]
        dgt <- which(duplicated(smap))
        if (length(dgt) > 0)
            smap <- smap[-dgt,]
        
        # Add metadata from the map
        rownames(smap) <- smap$transcript_id
        smap <- smap[names(gr),,drop=FALSE]
        gr$gene_id <- smap$gene_name
        gr$gene_name <- smap$gene_name
        gr$biotype <- smap$biotype
    }
    else {
        gr <- transcripts(txdb,columns=c("tx_name","gene_id","tx_type"))
        if (is(gr$gene_id,"CharacterList"))
            gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
        else
            gtmp <- as.character(gr$gene_id)
        gr$gene_id <- gtmp
        #gr$gene_id <- as.character(gr$gene_id)
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(lapply(gr$tx_type,function(x) x[1]))
        gr$tx_type <- NULL
    }
    
    ann <- as.data.frame(gr)
    ann <- ann[,c(1,2,3,6,7,5,8,9)]
    names(ann)[c(1,4)] <- c("chromosome","transcript_id")
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    #stranscript <- reduceTranscripts(GRanges(ann))
    annList <- reduceTranscripts(GRanges(ann))
    stranscript <- annList$model
    names(stranscript) <- as.character(stranscript$transcript_id)
    activeLength <- annList$length
    #names(activeLength) <- as.character(stranscript$transcript_id)
    names(activeLength) <- unique(stranscript$gene_id)
    
    if (asdf) {
        sann <- as.data.frame(stranscript)
        sann <- sann[,c(1,2,3,6,8,5,7,9)]
        names(sann)[c(1,4)] <- c("chromosome","transcript_id")
        attr(sann,"activeLength") <- activeLength
        return(sann)
    }
    else
        return(stranscript)
}

# ann <- .makeTranscriptExonFromTxDb(txdb,map,FALSE)
.makeTranscriptExonFromTxDb <- function(txdb,map,asdf) {
    # Stub
    # TODO: "summarized_exon_by_transcript.rda"
}

# ann <- .makeTranscriptUtrFromTxDb(txdb,map,FALSE)
.makeTranscriptUtrFromTxDb <- function(txdb,map,asdf) {
    utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
    utrGr <- unlist(utrList)
    
    if (length(utrGr) == 0) {
        warning("No UTR information was found in the provided GTF file! Will ",
            "return an empty object...",immediate.=TRUE)
        ann <- data.frame(chromosome=1,start=1,end=1,
            transcript_id=1,gene_id=1,strand=1,gene_name=1,
            biotype=1,row.names=1)
        # Strange but required to be compatible with GRanges
        ann <- ann[-1,]
        if (asdf)
            return(ann)
        else
            return(GRanges(ann))
    }
    
    utrGr$transcript_id <- names(utrGr)
    utrTmp <- as.data.frame(unname(utrGr))
    keep <- c("seqnames","start","end","transcript_id",
        "exon_rank","strand","exon_name")
    utr <- utrTmp[,keep]
    
    useMap <- TRUE
    if (any(is.na(utr$exon_name)) || is.null(map))
        useMap <- FALSE
    else
        rownames(utr) <- paste(utr$exon_name,utr$transcript_id,sep="_")
    
    if (useMap) {
        rownames(map) <- paste(map$exon_id,map$transcript_id,sep="_")
    
        # Different case with map here
        smap <- map[rownames(utr),]
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the gene case
        ir <- which(names(map)=="exon_id")
        if (length(ir) > 0)
            smap <- smap[,-ir,drop=FALSE]
        
        # Add metadata from the map
        #rownames(smap) <- smap$transcript_id
        #smap <- smap[rownames(utr),]
        utr$gene_id <- smap$gene_id
        utr$gene_name <- smap$gene_name
        utr$biotype <- smap$biotype
    }
    else {
        # Add metadata from the utr frame
        utr$gene_id <- utr$transcript_id
        utr$gene_name <- utr$gene_id
        
        # To get biotypes
        ge <- genes(txdb,columns=c("gene_id","tx_type"))
        ov <- findOverlaps(utrGr,ge)
        utrPos <- queryHits(ov)
        genePos <- subjectHits(ov)
        levs <- unique(genePos)
        
        dup <- which(duplicated(utrPos))
        if (length(dup) > 0) {
            utrPos <- utrPos[-dup]
            genePos <- genePos[-dup]
        }
        
        S <- split(utrPos,factor(genePos,levels=levs))
        utr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],lengths(S))
        #utr$biotype <- rep("gene",nrow(utr))
    }
    
    #ann <- utr[,c(1,2,3,4,7,6,8,9)]
    ann <- utr[,c(1,2,3,4,8,6,9,10)]
    names(ann)[1] <- "chromosome"
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

# ann <- .makeSumTranscriptUtrFromTxDb(txdb,map,FALSE)
.makeSumTranscriptUtrFromTxDb <- function(txdb,map,asdf) {
    utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
    utrGr <- unlist(utrList)
    
    if (length(utrGr) == 0) {
        warning("No UTR information was found in the provided GTF file! Will ",
            "return an empty object...",immediate.=TRUE)
        ann <- data.frame(chromosome=1,start=1,end=1,
            transcript_id=1,gene_id=1,strand=1,gene_name=1,
            biotype=1,row.names=1)
        # Strange but required to be compatible with GRanges
        ann <- ann[-1,]
        if (asdf)
            return(ann)
        else
            return(GRanges(ann))
    }
    
    utrGr$transcript_id <- names(utrGr)
    utrTmp <- as.data.frame(unname(utrGr))
    keep <- c("seqnames","start","end","transcript_id",
        "exon_rank","strand","exon_name")
    utr <- utrTmp[,keep]
    
    useMap <- TRUE
    if (any(is.na(utr$exon_name)) || is.null(map))
        useMap <- FALSE
    else
        rownames(utr) <- paste(utr$exon_name,utr$transcript_id,sep="_")
    
    if (useMap) {
        rownames(map) <- paste(map$exon_id,map$transcript_id,sep="_")
        
        # Different case with map here
        smap <- map[rownames(utr),]
        
        # We need to add gene_id, gene_name, biotype from the map
        # Remove the exon_id column from the map for the gene case
        ir <- which(names(map)=="exon_id")
        if (length(ir) > 0)
            smap <- smap[,-ir,drop=FALSE]
        
        # Add metadata from the map
        utr$gene_id <- smap$gene_id
        utr$gene_name <- smap$gene_name
        utr$biotype <- smap$biotype
    }
    else {
        # Add metadata from the utr frame
        utr$gene_id <- utr$transcript_id
        utr$gene_name <- utr$gene_id
        
        # To get biotypes
        ge <- genes(txdb,columns=c("gene_id","tx_type"))
        ov <- findOverlaps(utrGr,ge)
        utrPos <- queryHits(ov)
        genePos <- subjectHits(ov)
        levs <- unique(genePos)
        
        dup <- which(duplicated(utrPos))
        if (length(dup) > 0) {
            utrPos <- utrPos[-dup]
            genePos <- genePos[-dup]
        }
        
        S <- split(utrPos,factor(genePos,levels=levs))
        utr$biotype <- rep(unlist(ge$tx_type,use.names=FALSE)[levs],lengths(S))
    }
    
    #ann <- utr[,c(1,2,3,4,7,6,8,9)]
    ann <- utr[,c(1,2,3,4,8,6,9,10)]
    names(ann)[1] <- "chromosome"
    ann$chromosome <- as.character(ann$chromosome)
    ann <- ann[order(ann$chromosome,ann$start),]
    
    message("  summarizing UTRs per gene for imported GTF")
    #s3utrTranscript <- reduceTranscriptsUtr(GRanges(ann))
    annList <- reduceTranscriptsUtr(GRanges(ann))
    s3utrTranscript <- annList$model
    names(s3utrTranscript) <- 
        as.character(s3utrTranscript$transcript_id)
    activeLength <- annList$length
    #names(activeLength) <- as.character(s3utrTranscript$transcript_id)
    # Remember, gene_id in UTR's per transcript is the transcript
    names(activeLength) <- unique(as.character(s3utrTranscript$gene_id))
    
    if (asdf) {
        sann <- as.data.frame(s3utrTranscript)
        sann <- sann[,c(1,2,3,6,8,5,7,9)]
        names(sann)[c(1,4)] <- c("chromosome","transcript_id")
        attr(sann,"activeLength") <- activeLength
        return(sann)
    }
    else
        return(s3utrTranscript)
}

# ann <- .makeExonExonFromTxDb(txdb,map,FALSE)
.makeExonExonFromTxDb <- function(txdb,map,asdf) {
    if (!is.null(map)) {
        gr <- exons(txdb,columns="exon_name")
        
        # There are cases where exons are unnamed and just the structure 
        useMap <- TRUE
        if (any(is.na(gr$exon_name)))
            useMap <- FALSE
        
        if (useMap) {
            names(gr) <- gr$exon_name
            # We need to add gene_id, gene_name, biotype from the map
            # Remove the transcript_id column from the map for the gene
            # case
            ir <- which(names(map)=="transcript_id")
            if (length(ir) > 0)
                smap <- map[,-ir,drop=FALSE]
            dgt <- which(duplicated(smap))
            if (length(dgt) > 0)
                smap <- smap[-dgt,]
            
            # Add metadata from the map
            rownames(smap) <- smap$exon_id
            smap <- smap[names(gr),,drop=FALSE]
            gr$gene_id <- smap$gene_id
            gr$gene_name <- smap$gene_name
            gr$biotype <- smap$biotype
            
            ann <- as.data.frame(unname(gr))
            #ann <- ann[,c(1,2,3,6,8,5,7,9)]
            ann <- ann[,c(1,2,3,6,7,5,8,9)]
            names(ann)[c(1,4)] <- c("chromosome","exon_id")
            ann$chromosome <- as.character(ann$chromosome)
            ord <- order(ann$chromosome,ann$start)
            ann <- ann[ord,]
            smap <- smap[ord,]
            if (length(unique(ann$exon_id)) == length(ann$exon_id))
                rownames(ann) <- as.character(ann$exon_id)
            else
                rownames(ann) <- rownames(smap)
        }
        else {
            # In order to create the exon annotation, we need to manually 
            # overlap and assign to genes
            ge <- genes(txdb)
            ov <- findOverlaps(gr,ge)
            exonPos <- queryHits(ov)
            genePos <- subjectHits(ov)
            levs <- unique(genePos)
            
            dup <- which(duplicated(exonPos))
            if (length(dup) > 0) {
                exonPos <- exonPos[-dup]
                genePos <- genePos[-dup]
            }
            
            S <- split(exonPos,factor(genePos,levels=levs))
            gr$gene_id <- gr$gene_name <- rep(names(ge)[levs],lengths(S))
            gr$biotype <- rep("gene",length(gr))
            if (is(gr$gene_id,"CharacterList"))
                gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
            else
                gtmp <- gr$gene_id
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
            gr$exon_name <- names(gr)
            
            ann <- as.data.frame(unname(gr))
            ann <- ann[,c(1,2,3,6,7,5,8,9)]
            names(ann)[c(1,4)] <- c("chromosome","exon_id")
            ann$chromosome <- as.character(ann$chromosome)
            ord <- order(ann$chromosome,ann$start)
            ann <- ann[ord,]
            if (length(unique(ann$exon_id)) == length(ann$exon_id))
                # Should always be TRUE
                rownames(ann) <- as.character(ann$exon_id)
        }
    }
    else {
        gr <- exons(txdb,columns=c("exon_name","gene_id","tx_type"))
        if (is(gr$gene_id,"CharacterList"))
            gtmp <- unlist(lapply(gr$gene_id,function(x) x[1]))
        else
            gtmp <- gr$gene_id
        if (any(is.na(gr$exon_name)))            
            names(gr) <- paste(seqnames(gr),":",start(gr),"-",end(gr),"_",
                gtmp,sep="")
        else
            names(gr) <- gr$exon_name
        gr$gene_id <- gtmp
        gr$gene_name <- gr$gene_id
        gr$biotype <- unlist(lapply(gr$tx_type,function(x) x[1]))
        gr$tx_type <- NULL
        ann <- as.data.frame(unname(gr))
        ann <- ann[,c(1,2,3,6,7,5,8,9)]
        names(ann)[c(1,4)] <- c("chromosome","exon_id")
        ann$chromosome <- as.character(ann$chromosome)
        ord <- order(ann$chromosome,ann$start)
        ann <- ann[ord,]
        if (length(unique(ann$exon_id)) == length(ann$exon_id))
            # Should always be TRUE
            rownames(ann) <- as.character(ann$exon_id)
    }
    
    if (asdf)
        return(ann)
    else
        return(GRanges(ann))
}

.makeIdMap <- function(grdf) {
    # Some basic names must be present in order to make a map
    if (!any(c("exon_id","transcript_id","gene_id") %in% names(grdf))
        || any(is.na(grdf$exon_id))
        || any(is.na(grdf$transcript_id))
        || any(is.na(grdf$gene_id)))
        return(NULL)

    #hasGeneName <- hasBiotype <- FALSE
    if (!all(is.na(grdf$gene_name)) && !all(is.na(grdf$gene_biotype))) {
        map <- grdf[,c("exon_id","transcript_id","gene_id","gene_name",
            "gene_biotype")]
        names(map)[5] <- "biotype"
        #hasGeneName <- hasBiotype <- TRUE
    }
    else if (all(is.na(grdf$gene_name)) && !all(is.na(grdf$gene_biotype))) {
        map <- grdf[,c("exon_id","transcript_id","gene_id","gene_id",
            "gene_biotype")]
        names(map)[4:5] <- c("gene_name","biotype")
        #hasBiotype <- TRUE
    }
    else if (!all(is.na(grdf$gene_name)) && all(is.na(grdf$gene_biotype))) {
        map <- grdf[,c("exon_id","transcript_id","gene_id","gene_name")]
        map$gene_biotype <- rep("gene",nrow(map))
        names(map)[5] <- "biotype"
        #hasGeneName <- TRUE
    }
    else {
        map <- grdf[,c("exon_id","transcript_id","gene_id","gene_id")]
        map$gene_biotype <- rep("gene",nrow(map))
        names(map)[4:5] <- c("gene_name","biotype")
    }
    nag <- which(is.na(map$gene_name))
    if (length(nag) > 0)
        map$gene_name[nag] <- map$gene_id[nag]
    nab <- which(is.na(map$biotype))
    if (length(nab) > 0)
        map$biotype[nag] <- "gene"
    return(map)
}

initDatabase <- function(db) {
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to build the annotation database!")
    if (missing(db))
        stop("A database file must be provided!")
    drv <- dbDriver("SQLite")
    if (file.exists(db))
        # The database has been created at least with the tables defined
        con <- dbConnect(drv,dbname=db)
    else {
        # Create database file and define tables
        con <- dbConnect(drv,dbname=db)
        rs <- .initTables(con)
    }
    return(con)
}

.initTables <- function(con) {
    queries <- .localTblDef()
    rs <- dbSendQuery(con,queries[[1]])
    if (dbHasCompleted(rs))
        dbClearResult(rs)
    for (n in names(queries)) {
        rs <- dbSendStatement(con,queries[[n]])
        if (dbHasCompleted(rs))
            dbClearResult(rs)
    }
}

.chromInfoToSeqInfoDf <- function(ci,o="custom",circ=FALSE,asSeqinfo=FALSE) {
    if (asSeqinfo)
        return(Seqinfo(seqnames=rownames(ci),seqlengths=ci[,1],
            isCircular=rep(circ,nrow(ci)),genome=o))
    else
        return(data.frame(chromosome=rownames(ci),length=as.integer(ci[,1])))
}

.chromInfoWrapperGID <- function(o) {
    if (packageVersion("GenomeInfoDb")>=1.23)
        return(GenomeInfoDb::getChromInfoFromUCSC(o))
    else
        return(GenomeInfoDb::fetchExtendedChromInfoFromUCSC(o))
}

.validateEnsemblVersions <- function(o,v) {
    ea <- biomaRt::listEnsemblArchives()
    found <- as.character(v) %in% ea[,"version"]
    if (any(!found)) {
        nf <- which(!found)
        warning("Version(s) ",paste0(v[nf],collapse=", ")," not found in ",
            "biomaRt archives for current organism!\nWill use the latest ",
            "available version(s)...",immediate.=TRUE)
        newv <- v[-nf]
        if (length(newv) == 0) {
            warning("No Ensembl versions left after validation! Will return ",
                "only the latest available...",immediate.=TRUE)
            u2e <- ucscToEnsembl()
            vers <- u2e[[o]]
            for (newv in rev(vers)) {
                newi <- grep(as.character(newv),ea[,"version"])
                if (length(newi) > 0)
                    return(newv)
            }
        }
        return(newv)
    }
    return(v)
}
