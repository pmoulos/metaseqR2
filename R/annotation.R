buildAnnotationStore <- function(organisms,sources,
    home=file.path(path.expand("~"),".metaseqR"),forceDownload=TRUE,rc=NULL) {

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
    
    for (s in sources) {
        for (o in organisms) {
            # Retrieving genome info
            message("Retrieving genome information for ",o," from ",s)
            sf <- getSeqInfo(o)
            
            # Now, we must build the directory structure. For Ensembl it's
            # obvious as it has official versions. For UCSC we need to keep a
            # date track. Then inside metaseqR, if date not provided, it will
            # automatically detect the latest one (as in Ensembl with versions,
            # if not provided).
            if (s == "ensembl") {
                if (orgIsList)
                    vs <- orgList[[o]]
                else {
                    vss <- getUcscToEnsembl(o)
                    vs <- vss[length(vss)]
                }
            }
            else if (s %in% getSupportedUcscDbs())
                vs <- format(Sys.Date(),"%Y%m%d")
                        
            for (v in vs) {
                storePath <- file.path(home,s,o,v)
                if (!dir.exists(storePath))
                    dir.create(storePath,recursive=TRUE,mode="0755")
                    
                # Retrieve gene annotations
                if (file.exists(file.path(storePath,"gene.rda")) 
                    && !forceDownload)
                    message("Gene annotation for ",o," from ",s," version ",v,
                        " has already been created and will be skipped. If ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving gene annotation for ",o," from ",s,
                        " version ",v)
                    ann <- getAnnotation(o,"gene",refdb=s,ver=v,rc=rc)
                    gene <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    save(gene,file=file.path(storePath,"gene.rda"),
                        compress=TRUE)
                }
                
                # Retrieve transcript annotations
                if (file.exists(file.path(storePath,"transcript.rda"))
                    && !forceDownload)
                    message("Transcript annotation for ",o," from ",s,
                        " version ",v," has already been created and will be ",
                        "skipped. If you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    message("Retrieving transcript annotation for ",o,
                        " from ",s," version ",v)
                    ann <- getAnnotation(o,"transcript",refdb=s,ver=v,rc=rc)
                    transcript <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    save(transcript,file=file.path(storePath,"transcript.rda"),
                        compress=TRUE)
                    message("Splitting transcripts per gene for ",o," from ",s,
                        " version ",v)
                    transcriptByGene <- split(transcript,transcript$gene_id)
                    save(transcriptByGene,file=file.path(storePath,
                        "transcript_by_gene.rda"),compress=TRUE)
                }
                
                # Then summarize the transcripts and write again with type 
                # sum_transcript
                if (file.exists(file.path(storePath,
                    "summarized_transcript.rda")) && !forceDownload)
                    message("Summarized transcript annotation for ",o," from ",
                        s," version ",v," has already been created and will ",
                        "be skipped. If you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!file.exists(file.path(storePath,"transcript.rda")))
                        stop("Transcript annotation for ",o," from ",s,
                            " version ",v," is required in order to build ",
                            "predefined merged exon regions for read counting.",
                            " Please rerun the buildAnnotationStore function ",
                            "with appropriate parameters.")
                    tr <- load(file.path(storePath,"transcript.rda"))
                    annGr <- transcript
                    message("Merging transcripts for ",o," from ",s,
                        " version ",v)
                    stranscript <- reduceTranscripts(annGr,rc=rc)
                    names(stranscript) <- 
                        as.character(stranscript$transcript_id)
                    save(stranscript,file=file.path(storePath,
                        "summarized_transcript.rda"),compress=TRUE)
                    message("Splitting transcripts per gene for ",o," from ",s,
                        " version ",v)
                    stranscriptByGene <- split(stranscript,stranscript$gene_id)
                    save(stranscriptByGene,file=file.path(storePath,
                        "summarized_transcript_by_gene.rda"),compress=TRUE)
                }
                
                # Retrieve 3' UTR annotations
                if (file.exists(file.path(storePath,"3utr.rda"))
                    && !forceDownload)
                    message("3' UTR annotation for ",o," from ",s," version ",
                        v," has already been created and will be skipped. If ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving 3' UTR annotation for ",o,
                        " from ",s," version ",v)
                    ann <- getAnnotation(o,"utr",refdb=s,ver=v,rc=rc)
                    utr <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    save(utr,file=file.path(storePath,"3utr.rda"),
                        compress=TRUE)
                    message("Splitting 3' UTR per gene for ",o," from ",s,
                        " version ",v)
                    utrByGene <- split(utr,utr$gene_id)
                    save(utrByGene,file=file.path(storePath,
                        "3utr_by_gene.rda"),compress=TRUE)
                }
                
                # Then summarize the 3'utrs per gene and write again with type 
                # sum_utr
                if (file.exists(file.path(storePath,
                    "summarized_3utr.rda")) && !forceDownload)
                    message("Summarized 3utr annotation for ",o," from ",
                        s," version ",v," has already been created and will ",
                        "be skipped. If you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!file.exists(file.path(storePath,"3utr.rda")))
                        stop("3' UTR annotation for ",o," from ",s," version ",
                            vs," is required in order to build predefined ",
                            "merged exon regions for read counting. Please ",
                            "rerun the buildAnnotationStore function with ",
                            "appropriate parameters.")
                    tr <- load(file.path(storePath,"3utr.rda"))
                    annGr <- utr
                    message("Merging 3utrs for ",o," from ",s," version ",v)
                    s3utr <- reduceTranscripts(annGr,rc=rc)
                    names(s3utr) <- as.character(s3utr$gene_id)
                    save(s3utr,file=file.path(storePath,"summarized_3utr.rda"),
                        compress=TRUE)
                    message("Splitting 3' UTRs per gene for ",o," from ",s,
                        " version ",v)
                    s3utrByGene <- split(s3utr,s3utr$gene_id)
                    save(s3utrByGene,file=file.path(storePath,
                        "summarized_3utr_by_gene.rda"),compress=TRUE)
                }
                
                # Then summarize the 3'utrs per transcript and write again with 
                # type sum_utr
                if (file.exists(file.path(storePath,
                    "summarized_3utr_transcript.rda")) && !forceDownload)
                    message("Summarized 3utr annotation per transcript for ",o,
                        " from ",s," version ",v," has already been created  ",
                        "and will be skipped. If you wish to recreate it ",
                        "choose forceDownload = TRUE.")
                else {
                    if (!file.exists(file.path(storePath,"3utr.rda")))
                        stop("3' UTR annotation per transcript for ",o," from ",
                            s," version ",v," is required in order to build ",
                            "predefined merged exon regions for read ",
                            "counting. Please rerun the buildAnnotationStore ",
                            "function with appropriate parameters.")
                    tr <- load(file.path(storePath,"3utr.rda"))
                    annGr <- utr
                    message("Merging 3utrs for ",o," from ",s," version ",v)
                    s3utrTranscript <- reduceTranscriptsUtr(annGr,rc=rc)
                    names(s3utrTranscript) <- 
                        as.character(s3utrTranscript$transcript_id)
                    save(s3utrTrascript,file=file.path(storePath,
                        "summarized_3utr_transcript.rda"),compress=TRUE)
                    message("Splitting 3' UTRs per gene for ",o," from ",s,
                        " version ",v)
                    s3utrTranscriptByTranscript <- 
                        split(s3utrTranscript,s3utrTranscript$transcript_id)
                    save(s3utrTranscriptByTranscript,file=file.path(storePath,
                        "summarized_3utr_transcript_by_transcript.rda"),
                        compress=TRUE)
                }
                
                # Retrieve exon annotations
                if (file.exists(file.path(storePath,"exon.rda")) 
                    && !forceDownload)
                    message("Exon annotation for ",o," from ",s," version ",v,
                        " has already been created and will be skipped. If ",
                        "you wish to recreate it choose forceDownload = TRUE.")
                else {
                    message("Retrieving exon annotation for ",o," from ",s,
                        " version ",v)
                    ann <- getAnnotation(o,"exon",refdb=s,ver=v,rc=rc)
                    exon <- makeGRangesFromDataFrame(
                        df=ann,
                        seqinfo=sf,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    save(exon,file=file.path(storePath,"exon.rda"),
                        compress=TRUE)
                    message("Splitting exons per gene for ",o," from ",s,
                        " version ",v)
                    exonByGene <- split(exon,exon$gene_id)
                    save(exonByGene,file=file.path(storePath,
                        "exon_by_gene.rda"),compress=TRUE)
                }
                
                # Then summarize the exons and write again with type sum_exon
                if (file.exists(file.path(storePath,"summarized_exon.rda")) 
                    && !forceDownload)
                    message("Summarized exon annotation for ",o," from ",s,
                        " version ",v," has already been created and will be ",
                        "skipped. If you wish to recreate it choose ",
                        "forceDownload = TRUE.")
                else {
                    if (!file.exists(file.path(storePath,"exon.rda")))
                        stop("Exon annotation for ",o," from ",s," version ",v,
                            " is required in order to build predefined merged ",
                            "exon regions for RNA-Seq (exon) coverage ",
                            "calculations. Please rerun the ",
                            "buildAnnotationStore function with appropriate ",
                            "parameters.")
                    ex <- load(file.path(storePath,"exon.rda"))
                    annGr <- unlist(exon,use.names=FALSE)
                    message("Merging exons for ",o," from ",s," version ",v)
                    annList <- reduceExons(annGr,rc=rc)
                    sexon <- annList$model
                    names(sexon) <- as.character(sexon$exon_id)
                    activeLength <- annList$length
                    names(activeLength) <- unique(sexon$gene_id)
                    save(sexon,activeLength,
                        file=file.path(storePath,"summarized_exon.rda"),
                        compress=TRUE)
                    message("Splitting summarized exons per gene for ",o,
                        " from ",s," version ",v)
                    sexonByGene <- split(sexon,sexon$gene_id)
                    save(sexonByGene,activeLength,file=file.path(storePath,
                        "summarized_exon_by_gene.rda"),compress=TRUE)
                }
            }
        }
    }
}

# There should be another function for importing and registering annotation
importAnnotation <- function() {
	# TODO: Lots of checks!
	# 1. According to type, some columns MUST exist
	# 2. Is it a GTF? If yes we must deploy the custom GTF reading 
	if (file.exists(genome)) {
		genome <- read.delim(genome) # Must be bed like
		rownames(genome) <- as.character(genome[,4])
		genomeRanges <- makeGRangesFromDataFrame(
			df=genome,
			keep.extra.columns=TRUE
		)
	}
}

loadAnnotation <- function(genome=getSupportedOrganisms(),
	refdb=getSupportedRefDbs,level=c("gene","transcript","exon"),
	type=c("gene","exon","utr"),version="auto",
	home=file.path(path.expand("~"),".metaseqR"),rc=NULL) {
	
	home <- "/media/raid/tmp/metaseqR2/test_ann"
	
	# Firstly do some argument checking
	genome <- genome[1]
	checkTextArgs("genome",genome,c("hg18","hg19","hg38","mm9","mm10","rn5",
        "rn6","dm3","dm6","danrer7","pantro4","susscr3","tair10","equcab2",
        "custom"),multiarg=FALSE)
    refdb <- refdb[1]
    checkTextArgs("refdb",refdb,c("ensembl","ucsc","refseq"),multiarg=FALSE)
    level <- level[1]
    checkTextArgs("level",level,c("gene","transcript","exon"),multiarg=FALSE)
    type <- type[1]
    checkTextArgs("type",type,c("gene","exon","utr"),multiarg=FALSE)
    
	# Check if local storage has been set
	storageHome <- file.path(home,refdb,genome)
	if (dir.exists(storageHome)) {
		# Decide on version
		if (version != "auto") {
			storageHomeVersion <- file.path(storageHome,version)
			if (!dir.exists(storageHomeVersion)) {
				warning("The annotation directory ",storageHomeVersion,
					" does not seem to exist! Have you run ",
					"buildAnnotationStorage or possibly mispelled? Will ",
					"use newest existing version.",immediate.=TRUE)
				version <- "auto"
			}
		}
		if (version == "auto") {
			vers <- suppressWarnings(as.numeric(dir(storageHome)))
			if (any(is.na(vers))) {
				of <- vers[which(is.na(vers))]
				stop("Corrupted annotation storage directory ",
					"contents! ->",of,"<- Either remove offending ",
					"files/directories or rebuild.")
			}
			vers <- sort(vers,decreasing=TRUE)
			version <- vers[1]
			storageHomeVersion <- file.path(storageHome,version)
		}
		
		switch(level,
			gene = {
				switch(type,
					gene = {
						g <- load(file.path(storageHomeVersion,"gene.rda"))
						return(gene)
					},
					exon = {
						e <- load(file.path(storageHomeVersion,
							"summarized_exon_by_gene.rda"))
						return(sexon)
					},
					utr = {
						u <- load(file.path(storageHomeVersion,
							"summarized_3utr_by_gene.rda"))
						return(s3utrByGene)
					}
				)
			},
			transcript = {
				switch(type,
					gene = {
						g <- load(file.path(storageHomeVersion,
							"transcript.rda"))
						return(transcript)
					},
					exon = {
						# Stub
						# TODO: "summarized_exon_by_transcript.rda"
					},
					utr = {
						u <- load(file.path(storageHomeVersion,
							"summarized_3utr_transcript_by_transcript.rda"))
						return(s3utrTranscriptByTranscript)
					}
				)
			},
			exon = {
				switch(type,
					exon = {
						e <- load(file.path(storageHomeVersion,
							"summarized_exon.rda"))
					}
				)
			}
		)
	}
	else { # On the fly latest
		message("Getting latest annotation on the fly for ",genome," from ",
			refdb)
		sf <- getSeqInfo(genome)
		
		switch(level,
			gene = {
				switch(type,
					gene = {
						ann <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
						annGr <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
					},
					exon = {
						ann <- getAnnotation(genome,"exon",refdb=refdb,rc=rc)
						annGr <- unlist(exon,use.names=FALSE)
						message("Merging exons for ",genome," from ",refdb)
						annList <- reduceExons(annGr,rc=rc)
						sexon <- annList$model
						names(sexon) <- as.character(sexon$exon_id)
						activeLength <- annList$length
						names(activeLength) <- unique(sexon$gene_id)
						annGr <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
						attr(annGr,"activeLength") <- activeLength
						#annGr <- split(annGr,annGr$gene_id)
						#attr(annGr,"activeLength") <- activeLength
					},
					utr = {
						ann <- getAnnotation(genome,"utr",refdb=refdb,rc=rc)
						utr <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
						message("Merging 3utrs for ",genme," from ",refdb)
						annGr <- reduceTranscripts(utr,rc=rc)
						names(annGr) <- as.character(annGr$gene_id)
						#annGr <- split(annGr,annGr$gene_id)
					}
				)
			},
			transcript = {
				switch(type,
					gene = {
						ann <- getAnnotation(o,"transcript",refdb=s,ver=v,rc=rc)
						anGr <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
					},
					exon = {
						# Stub
						# TODO: "summarized_exon_by_transcript.rda"
					},
					utr = {
						ann <- getAnnotation(genome,"utr",refdb=refdb,rc=rc)
						annGr <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
						message("Merging 3utrs for ",genome," from ",refdb)
						annGr <- reduceTranscriptsUtr(annGr,rc=rc)
						names(annGr) <- as.character(annGr$transcript_id)
						#annGr <- split(annGr,annGr$transcript_id)
					}
				)
			},
			exon = {
				switch(type,
					exon = {
						ann <- getAnnotation(o,"exon",refdb=s,ver=v,rc=rc)
						exon <- makeGRangesFromDataFrame(
							df=ann,
							seqinfo=sf,
							keep.extra.columns=TRUE,
							seqnames.field="chromosome"
						)
						annGr <- unlist(exon,use.names=FALSE)
						message("Merging exons for ",genome," from ",refdb)
						annList <- reduceExons(annGr,rc=rc)
						annGr <- annList$model
						names(annGr) <- as.character(annGr$exon_id)
						activeLength <- annList$length
						names(activeLength) <- unique(annGr$gene_id)
						attr(annGr,"activeLength") <- activeLength
					}
				)
			}
		)
		
		return(annGr)
	}
}

correctTranscripts <- function(ann) {
    rownames(ann) <- paste("T",1:nrow(ann),sep="_")
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

reduceTranscripts <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MET",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    return(do.call("c",redList))
}

reduceTranscriptsUtr <- function(gr,rc=NULL) {
    trans <- unique(as.character(gr$transcript_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(trans,function(x,a,g,b) {
        tmp <- a[a$transcript_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MEU",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    return(do.call("c",redList))
}

reduceExons <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    redList <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            exon_id=paste(x,"MEX",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    len <- unlist(cmclapply(redList,function(x) {
        return(sum(width(x)))
    },rc=rc))
    names(len) <- names(redList)
    return(list(model=do.call("c",redList),length=len))
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
        bm <- getBM(attributes=getGeneAttributes(org),mart=mart)
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
        bm <- getBM(attributes=getTranscriptAttributes(org),mart=mart)
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
        bm <- getBM(attributes=getTranscriptUtrAttributes(org),mart=mart)
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
        bm <- getBM(attributes=getExonAttributes(org),mart=mart)        
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
        dbCreds <- .getUcscCredentials()
        drv <- dbDriver("MySQL")
        con <- dbConnect(drv,user=dbCreds[2],password=NULL,dbname=dbOrg,
            host=dbCreds[1])
        if (type == "utr") {
            # Proceed to GTF generation workflow
            rawAnn <- getUcscUtr(org,refdb)
        }
        else
            query <- getUcscQuery(org,type,refdb)
        rawAnn <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    else {
        # This should return the same data frame as the db query
        if (type == "transcript")
            tmpSqlite <- getUcscDbl(org,"gene",refdb)
        else if (type == "utr") {
            # Proceed to GTF generation workflow
            # rawAnn <- getUcscUtr(org,refdb)
        }
        else
            tmpSqlite <- getUcscDbl(org,type,refdb)
        drv <- dbDriver("SQLite")
        con <- dbConnect(drv,dbname=tmpSqlite)
        query <- getUcscQuery(org,type,refdb)
        rawAnn <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    if (type=="gene") {
        tmpAnn <- rawAnn
        tmpAnn <- tmpAnn[grep(chrsExp,tmpAnn$chromosome,perl=TRUE),]
        tmpAnn$chromosome <- as.character(tmpAnn$chromosome)
        rownames(tmpAnn) <- tmpAnn$transcript_id
        # Split the UCSC transcripts per gene name
        tmp <- split(tmpAnn,tmpAnn$gene_name)
        # Retrieve the longest transcript (as per Ensembl convention)
        ann <- do.call("rbind",cmclapply(tmp,function(x) {
            size <- x$end - x$start
            selected <- which(size == max(size))
            return(x[selected[1],,drop=FALSE])
        },rc=rc))
        names(ann)[4] <- "gene_id"
        gcContent <- getGcContent(ann,org)
        ann$gc_content <- gcContent
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
        ann <- ann[,c(1:4,8,5:7)]
    }
    else if (type=="exon") {
        rawAnn <- rawAnn[grep(chrsExp,rawAnn$chromosome,perl=TRUE),]
        exList <- cmclapply(as.list(1:nrow(rawAnn)),function(x,d,s) {
            r <- d[x,]
            starts <- as.numeric(strsplit(r[,"start"],",")[[1]])
            ends <- as.numeric(strsplit(r[,"end"],",")[[1]])
            nexons <- length(starts)
            ret <- data.frame(
                rep(r[,"chromosome"],nexons),
                starts,ends,
                paste(r[,"exon_id"],"_e",1:nexons,sep=""),
                rep(r[,"strand"],nexons),
                rep(r[,"transcript_id"],nexons),
                rep(r[,"gene_name"],nexons),
                rep(r[,"biotype"],nexons)
            )
            names(ret) <- names(r)
            rownames(ret) <- ret$exon_id
            return(ret)
        },rawAnn,validChrs,rc=rc)
        
        # For some reason rbind takes ages for large datasets... We have to 
        # split in chunks of 1000
        N <- length(exList)
        mo <- N%%chunkSize
        if (mo == 0) {
            fl <- N/chunkSize
            fac <- factor(rep(1:fl,each=chunkSize))
        }
        else {
            fl <- (N-mo)/chunkSize
            fac <- factor(c(rep(1:fl,each=chunkSize),rep(fl,mo)))
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
            gene_id=as.character(tmpAnn$transcript_id),
            strand=as.character(tmpAnn$strand),
            gene_name=as.character(tmpAnn$gene_name),
            biotype=as.character(tmpAnn$biotype)
        )
        rownames(ann) <- ann$exon_id
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
        
        # First get the UTRs
        preAnn <- getUcscUtr(org,refdb)
        preAnn <- preAnn[grep(chrsExp,preAnn$seqnames,perl=TRUE),]
        
        # Then try and attach gene names and biotypes...
        helperAnn <- loadAnnotation(org,"gene",refdb=refdb)
        map <- elementMetadata(helperAnn)
        map <- map[,c(1,3,4)]
        map <- as.data.frame(map[,c(1,3)],row.names=as.character(map[,2]))
        mapp <- map[preAnn$gene_id,]
        preAnn <- cbind(preAnn,mapp)
        preAnn <- preAnn[,c(1:4,6,8,7,9)]
        names(preAnn) <- c("chromosome","start","end","transcript_id","strand",
			"gene_id","gene_name","biotype")
		rownames(preAnn) <- paste("T",1:nrow(preAnn),sep="_")
		ann <- preAnn
    }
    
    ann <- ann[order(ann$chromosome,ann$start),]
    return(ann)
}

getGcContent <- function(ann,org) {
    if (missing(ann))
        stop("A valid annotation data frame must be provided in order to ",
            "retrieve GC-content.")
    org <- tolower(org[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    # Convert annotation to GRanges
    message("Converting annotation to GenomicRanges object...")
    if (packageVersion("GenomicRanges")<1.14)
        annGr <- GRanges(
            seqnames=Rle(ann[,1]),
            ranges=IRanges(start=ann[,2],end=ann[,3]),
            strand=Rle(ann[,6]),
            name=as.character(ann[,4])
        )
    else
        annGr <- makeGRangesFromDataFrame(
            df=ann,
            keep.extra.columns=TRUE,
            seqnames.field="chromosome"
        )
    bsg <- loadBsGenome(org)
    if (is(bsg,"BSgenome")) {
        message("Getting DNA sequences...")
        seqs <- getSeq(bsg,names=annGr)
        message("Getting GC content...")
        freqMatrix <- alphabetFrequency(seqs,as.prob=TRUE,baseOnly=TRUE)
        gcContent <- apply(freqMatrix,1,function(x) round(100*sum(x[2:3]),
            digits=2))
    }
    else
        gcContent <- rep(NA,nrow(ann))
    names(gcContent) <- as.character(ann[,4])
    return(gcContent)
}

getSeqInfo <- function(org) {
	sf <- tryCatch(GenomeInfoDb::fetchExtendedChromInfoFromUCSC(
		getUcscOrganism(org)),error=function(e) {
			message("GenomeInfoDb::fetchExtendedChromInfoFromUCSC ",
				"failed with the following error: ")
			message(e)
			message("")
			message("Trying a direct download...")
			getChromInfo(getUcscOrganism(org))
		},finally="")
	rownames(sf) <- as.character(sf[,1])
	sf <- sf[getValidChrs(org),]
	sf <- Seqinfo(seqnames=sf[,1],seqlengths=sf[,2],
		isCircular=sf[,3],genome=getUcscOrganism(org))
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
        dm6 = { return("dm3") },
        danrer7 = { return("danRer7") },
        danrer10 = { return("danRer10") },
        pantro4 = { return("panTro4") },
        pantro5 = { return("panTro5") },
        susscr3 = { return("susScr3") },
        susscr11 = { return("susScr3") },
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
            #warning("danRer7 is not supported by BSgenome! Please use Ensembl ",
            #    "as annotation source if GC content is important.",
            #    immediate.=TRUE)
            #return(NA)
            # Is default for refseq for some reason...
            return("BSgenome.Drerio.UCSC.danRer7")
        },
        danrer10 = {
            #warning("danRer10 is not supported by BSgenome! Please use Ensembl",
            #    " as annotation source if GC content is important.",
            #    immediate.=TRUE)
            #return(NA)
            #  Is default for refseq for some reason...
            return("BSgenome.Drerio.UCSC.danRer10")
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
    if (!requireNamespace("BiocInstaller"))
        stop("The Bioconductor package BiocInstaller is required to ",
            "proceed!")
    if (!requireNamespace("BSgenome"))
        stop("The Bioconductor package BSgenome is required to ",
            "proceed!")
    bsOrg <- getBsOrganism(org)
    if (!is.na(bsOrg)) {
        if (bsOrg %in% BSgenome::installed.genomes())
            bsObj <- getBSgenome(getUcscOrganism(org))
        else {
            BiocInstaller::biocLite(bsOrg)
            bsObj <- getBSgenome(getUcscOrganism(org))
        }
        return(bsObj)
    }
    else
        return(NA)
}

getChromInfo <- function(org,
    goldenPath="http://hgdownload.cse.ucsc.edu/goldenPath/") {
    download.file(paste(goldenPath,org,"/database/chromInfo.txt.gz",sep=""),
        file.path(tempdir(),"chromInfo.txt.gz"))
    chromInfo <- read.delim(file.path(tempdir(),"chromInfo.txt.gz"),
        header=FALSE)
    chromInfo <- chromInfo[,1:2]
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
    if (length(i) > 0)
        return(ea[i,"url"])
    else
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
        hg38=76:96,
        mm9=67,
        mm10=74:96,
        rn5=74:79,
        rn6=80:96,
        dm3=c(67,74:78),
        dm6=79:96,
        danrer7=c(67,74:79),
        danrer10=80:91,
        danrer11=92:96,
        pantro4=c(67,74:90),
        pantro5=91:96,
        #pantro6=,
        susscr3=c(67,74:89),
        susscr11=90:96,
        equcab2=c(67,74:96)
    ))
}

getHostOld <- function(org) {
    .Deprecated("getHost")
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("www.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("www.ensembl.org") },
        rn5 = { return("grch37.ensembl.org") },
        rn6 = { return("www.ensembl.org") },
        dm3 = { return("grch37.ensembl.org") },
        dm6 = { return("www.ensembl.org") },
        danrer7 = { return("grch37.ensembl.org") },
        danrer10 = { return("www.ensembl.org") },
        danrer11 = { return("www.ensembl.org") },
        pantro4 = { return("grch37.ensembl.org") },
        pantro5 = { return("www.ensembl.org") },
        #pantro6 = { return("www.ensembl.org") },
        susscr3 = { return("grch37.ensembl.org") },
        susscr11 = { return("www.ensembl.org") }
    )
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
                "chr31","chr4","chr5","chr6","chr7","chr8","chr9","chrX",#"chrY",
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

getGeneAttributes <- function(org) {
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

getTranscriptAttributes <- function(org) {
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

getTranscriptUtrAttributes <- function(org) {
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

getExonAttributes <- function(org) {
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
        "danrer7","danrer10","danrer11","pantro4","pantro5","susscr3",#"pantro6",
        "susscr11","equcab2","tair10"))
}

getSupportedUcscDbs <- function() {
    return(c("ucsc","refseq"))
}

#organisms <- list(
#	hg18=67,
#	hg19=75,
#	hg38=95:96,
#	mm9=67,
#	mm10=95:96,
#	rn5=79,
#	rn6=95:96,
#	dm3=78,
#	dm6=95:96,
#	danrer7=79,
#	danrer10=91,
#	danrer11=95:96,
#	pantro4=90,
#	pantro5=95:96,
#	susscr3=89,
#	susscr11=95:96,
#	equcab2=95:96
#)
#sources <- c("ensembl","ucsc","refseq")
