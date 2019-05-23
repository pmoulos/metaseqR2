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
    if (refdb == "ucsc") {
        command <- paste(
            "zcat ",file.path(tempdir(),paste(tableUtr,".txt.gz",sep=""))," | ",
            "cut -f1-10 - | ",genePredToGtf," file stdin ",gtfFile," -source=",
            tableUtr," -utr",sep=""
        )
        system(command)
        
        # Then make TxDb through makeTxDbFromGFF and transform to data frame
        # to return later
        message("  Importing GTF ",tableUtr," as TxDb")
        txdb <- makeTxDbFromGFF(gtfFile)
        utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
        utrGr <- unlist(utrList)
        utrGr$transcript_id <- names(utrGr)
        utrTmp <- as.data.frame(unname(utrGr))
        keep <- c("seqnames","start","end","transcript_id","exon_rank","strand")
        utr <- utrTmp[,keep]
        utr$gene_id <- utr$transcript_id
    }
    else if (refdb == "refseq") {
        command <- paste(
            "zcat ",file.path(tempdir(),paste(tableUtr,".txt.gz",sep=""))," | ",
            "cut -f2- | ",genePredToGtf," file stdin ",gtfFile," -source=",
            tableUtr," -utr",sep=""
        )
        system(command)
        
        # Then make TxDb through makeTxDbFromGFF and transform to data frame
        # to return later
        message("  Importing GTF ",tableUtr," as TxDb")
        txdb <- makeTxDbFromGFF(gtfFile)
        utrList <- threeUTRsByTranscript(txdb,use.names=TRUE)
        utrGr <- unlist(utrList)
        utrGr$transcript_id <- names(utrGr)
        utrTmp <- as.data.frame(unname(utrGr))
        keep <- c("seqnames","start","end","transcript_id","exon_rank","strand")
        utr <- utrTmp[,keep]
        
		# Import again the GTF with rtracklayer to create the names map
		message("  Importing GTF ",tableUtr," as GTF to make id map")
		desiredColumns <- c("type","gene_id","transcript_id","exon_id",
		 "gene_name","gene_biotype")
		gr <- import(gtfFile,format="gtf",colnames=desiredColumns,
			feature.type=GenomicFeatures:::GFF_FEATURE_TYPES)
		grdf <- as.data.frame(gr)
		grdf <- grdf[grdf$type=="exon",]
		map <- grdf[,c("transcript_id","gene_id")]
		
		# Remove total duplicates (e.g. exons mapping to multiple gns, trns)
		dgt <- which(duplicated(map))
		if (length(dgt) > 0)
			map <- map[-dgt,]
		geneIds <- map$gene_id
		names(geneIds) <- map$transcript_id
		utr$gene_id <- geneIds[utr$transcript_id]
    }
    
    return(utr)
}

#~ importCustomGtf <- function() {
#~      ########################################################################
#~      # The folowing chunk belongs to a more generic custom GTF importing
#~      # function. UCSC GTFs do not have gene names and biotypes
        
#~      # Import again the GTF with rtracklayer to create the names map
#~      message("  Importing GTF ",tableUtr," as GTF to make id map")
#~      desiredColumns <- c("type","gene_id","transcript_id","exon_id",
#~          "gene_name","gene_biotype")
#~      gr <- import(gtfFile,format="gtf",colnames=desiredColumns,
#~          feature.type=GenomicFeatures:::GFF_FEATURE_TYPES)
#~      grdf <- as.data.frame(gr)
#~      grdf <- grdf[grdf$type=="exon",]

#~      # Need to map gene ids to names and biotypes. We need a collapsed 
#~      # structure exon_id - transcript_id - gene_id - gane_name - biotype
#~      message("  Making ",tableUtr," id map")
#~      hasGeneName <- hasBiotype <- FALSE
#~      if (!all(is.na(grdf$gene_name)) && !all(is.na(grdf$gene_biotype))) {
#~          map <- grdf[,c("exon_id","transcript_id","gene_id","gene_name",
#~              "gene_biotype")]
#~          names(map)[5] <- "biotype"
#~          hasGeneName <- hasBiotype <- TRUE
#~      }
#~      else if (all(is.na(grdf$gene_name)) && !all(is.na(grdf$gene_biotype))) {
#~          map <- grdf[,c("exon_id","transcript_id","gene_id","gene_name")]
#~          map$biotype <- rep("gene",nrow(map))
#~          hasGeneName <- TRUE
#~      }
#~      else if (!all(is.na(grdf$gene_name)) && all(is.na(grdf$gene_biotype))) {
#~          map <- grdf[,c("exon_id","transcript_id","gene_id","gene_id",
#~              "gene_biotype")]
#~          names(map)[3:5] <- c("gene_id","gene_name","biotype")
#~          hasBiotype <- TRUE
#~      }
#~      else {
#~          map <- grdf[,c("exon_id","transcript_id","gene_id","gene_id",
#~              "gene_id")]
#~          map$biotype <- rep("gene",nrow(map))
#~          names(map)[3:4] <- c("gene_id","gene_name")
#~      }

#~      # The map is not of very much use if there are no gene names and
#~      # biotypes... The only true usage is gene_id to plug into utrs...
#~      # If has both gene name and biotype, then ok, otherwise we will have to
#~      # download or load annotations on the fly.
#~      # Let's  use it for map gene_id to transcript_id at first
#~      submap <- map[,c("transcript_id","gene_id")]
#~      # Remove total duplicates (e.g. exons mapping to multiple gns, trns)
#~      dgt <- which(duplicated(submap))
#~      if (length(dgt) > 0)
#~          submap <- submap[-dgt,]
#~      geneIds <- submap$gene_id
#~      names(geneIds) <- submap$transcript_id
#~      utr$gene_id <- geneIds[utr$transcript_id]
#~      ########################################################################
#~ }

getUcscTableNameUtr <- function(org,refdb) {
    switch(refdb,
        ucsc = {
            switch(org,
                hg18 = {
                    return("knownGene")
                },
                hg19 = {
                    return("knownGene")
                },
                hg38 = {
                    return("knownGene")
                },
                mm9 = {
                    return("knownGene")
                },
                mm10 = {
                    return("knownGene")
                },
                rn5 = {
                    return("mgcGenes")
                },
                rn6 = {
                    return("mgcGenes")
                },
                dm3 = {
                    return("flyBaseGene")
                },
                dm6 = {
                    warning("No UCSC Genome annotation for Drosophila ",
                        "melanogaster v6! Will use RefSeq instead...",
                        immediate.=TRUE)
                    return("refGene")
                },
                danrer7 = {
                    return("mgcGenes")
                },
                danrer10 = {
                    return("mgcGenes")
                    
                },
                danrer11 = {
                    warning("No UCSC Genome annotation for Danio rerio v11! ",
                        "Will use RefSeq instead...",immediate.=TRUE)
                    return("refGene")
                },
                pantro4 = {
                    warning("No UCSC Genome annotation for Pan ",
                        "troglodytes! Will use RefSeq instead...",
                        immediate.=TRUE)
                    return("refGene")
                },
                pantro5 = {
                    warning("No UCSC Genome annotation for Pan ",
                        "troglodytes! Will use RefSeq instead...",
                        immediate.=TRUE)
                    return("refGene")
                },
                susscr3 = {
                    warning("No UCSC Genome annotation for Sus ",
                        "scrofa v3! Will use RefSeq instead...",
                        immediate.=TRUE)
                        return("refGene")
                },
                susscr11 = {
                    warning("No UCSC Genome annotation for Sus ",
                        "scrofa v11! Will use RefSeq instead...",
                        immediate.=TRUE)
                    return("refGene")
                },
                equcab2 = {
                    warning("No UCSC Genome annotation for Equus ",
                        "caballus v2! Will use RefSeq instead...",
                        immediate.=TRUE)
                    return("refGene")
                }
            )
        },
        refseq = {
            return("refGene")
        }
    )
}

#' Download annotation from UCSC servers, according to organism and source
#'
#' Directly downloads UCSC and RefSeq annotation files from UCSC servers to be
#' used with metaseqR. This functionality is used when the package RMySQL is not
#' available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A data frame with annotation elements.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.file <- get.ucsc.file("hg18","gene","ucsc")
#'}
getUcscDbl <- function(org,refdb="ucsc") {
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    checkTextArgs("refdb",refdb,getSupportedUcscDbs())
    
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to use annotation from UCSC!")

    httpBase <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/",
        getUcscOrganism(org),"/database/",sep="")
    tableDefs <- getUcscTabledef(org,refdb,"fields")
    fileList <- vector("list",length(tableDefs))
    names(fileList) <- names(tableDefs)
    for (n in names(fileList))
        fileList[[n]] <- paste(httpBase,n,".txt.gz",sep="")
        
    # Fill the fields for each table
    drv <- dbDriver("SQLite")
    dbTmp <- tempfile()
    con <- dbConnect(drv,dbname=dbTmp)
    message("  Retrieving tables for temporary SQLite ",refdb," ",org,
        " subset database")
    for (n in names(fileList)) {
        message("    Retrieving table ",n)
        download.file(fileList[[n]],file.path(tempdir(),
            paste(n,".txt.gz",sep="")),quiet=TRUE)
        sqlDf <- read.delim(file.path(tempdir(),paste(n,".txt.gz",sep="")),
            row.names=NULL,header=FALSE,strip.white=TRUE)
        names(sqlDf) <- tableDefs[[n]]
        dbWriteTable(con,n,sqlDf,row.names=FALSE)
    }
    dbDisconnect(con)
    return(dbTmp)
}

#' Get SQLite UCSC table defintions, according to organism and source
#'
#' Creates a list of UCSC Genome Browser database tables and their SQLite
#' definitions with the purpose of creating a temporary SQLite database to be 
#' used used with metaseqR. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A list with SQLite table definitions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.tabledefs <- get.ucsc.tables("hg18","gene","ucsc")
#'}
getUcscTabledef <- function(org,refdb="ucsc",what="queries") {
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    what <- tolower(what[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    checkTextArgs("refdb",refdb,getSupportedUcscDbs())
    checkTextArgs("what",what,c("queries","fields"))
    switch(refdb,
		ucsc = {
			return(.getUcscTabledefUcsc(org,what))
		},
		refseq = {
			return(.getUcscTabledefRefseq(org,what))
		}
	)
}

#' Create SQLite UCSC table template defintions
#'
#' Returns an SQLIte table template defintion, according to  UCSC Genome Browser 
#' database table schemas. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines. Internal use only.
#'
#' @param tab name of UCSC database table.
#' @param what \code{"queries"} for SQLite table definitions or \code{"fields"}
#' for table column names.
#' @return An SQLite table definition.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.table.tmpl <- get.ucsc.tbl.tpl("knownCanonical")
#'}
getUcscTblTpl <- function(tab,what="queries") {
    if (what=="queries") {
        switch(tab,
            knownCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`knownCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownGene = {
                return(paste(
                    "CREATE TABLE",
                    "`knownGene` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`proteinID` TEXT NOT NULL DEFAULT '',",
                    "`alignID` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            refFlat = {
                return(paste("CREATE TABLE",
                    "`refFlat` (",
                    "`geneName` TEXT NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            refGene = {
                return(paste("CREATE TABLE",
                    "`refGene` (",
                    "`bin` UNSIGNED INTEGER NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`score` INTEGER NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    "`name2` TEXT NOT NULL,",
                    "`cdsStartStat` TEXT NOT NULL,",
                    "`cdsEndStat` TEXT NOT NULL,",
                    "`exonFrames` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            knownToEnsembl = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToEnsembl` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            ensemblSource = {
                return(paste(
                    "CREATE TABLE",
                    "`ensemblSource` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`source` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            mgcGenes = {
                return(paste(
                    "CREATE TABLE `mgcGenes` (",
                    "`bin` UNSIGNED INTEGER NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`score` INTEGER DEFAULT NULL,",
                    "`name2` TEXT NOT NULL,",
                    "`cdsStartStat` TEXT NOT NULL,",
                    "`cdsEndStat` TEXT NOT NULL,",
                    "`exonFrames` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            ensemblToGeneName = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToGeneName` (",
                    "`name` TEXT NOT NULL,",
                    "`value` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER unsigned NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            flyBaseGene = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseGene` (",
                    "`bin` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            }
        )
    }
    else if (what=="fields") {
        switch(tab,
            knownCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                "transcript","protein"))
            },
            knownGene = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","proteinID",
                    "alignID"))
            },
            knownToRefSeq = {
                return(c("name","value"))
            },
            refFlat = {
                return(c("geneName","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            refGene = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score",
                "name2","cdsStartStat","cdsEndStat","exonFrames"))
            },
            knownToEnsembl = {
                return(c("name","value"))
            },
            ensemblSource = {
                return(c("name","source"))
            },
            mgcGenes = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
                    "score","name2","cdsStartStat","cdsEndStat","exonFrames"
                ))
            },
            ensemblToGeneName = {
                return(c("name","value"))
            },
            flyBaseCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                    "transcript","protein"))
            },
            flyBaseGene = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            flyBaseToRefSeq = {
                return(c("name","value"))
            }
        )
    }
}

#' Return queries for the UCSC Genome Browser database, according to organism and
#' source
#'
#' Returns an SQL query to be used with a connection to the UCSC Genome Browser
#' database and fetch metaseqR supported organism annotations. This query is
#' constructed based on the data source and data type to be returned.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A valid SQL query.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.query <- get.ucsc.query("hg18","gene","ucsc")
#'}
getUcscQuery <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("type",type,c("gene","exon","transcript"))
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    return(.getUcscQueryUcscGene(org))
                },
                refseq = {
                    return(.getUcscQueryRefseqGene(org))
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    return(.getUcscQueryUcscExon(org))
                },
                refseq = {
                    return(.getUcscQueryRefseqExon(org))                    
                }
            )
        },
        transcript = {
            switch(refdb,
                ucsc = {
                    return(.getUcscQueryUcscTranscript(org))
                },
                refseq = {
                    return(.getUcscQueryRefseqTranscript(org))
                }
            )
        }
    )
}

.getUcscTabledefUcsc <- function(org,what="queries") {
    switch(org,
        hg18 = {
            return(list(
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownGene=getUcscTblTpl("knownGene",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        hg19 = {
            return(list(
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownGene=getUcscTblTpl("knownGene",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        hg38 = {
            return(list(
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownGene=getUcscTblTpl("knownGene",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        mm9 = {
            return(list(
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownGene=getUcscTblTpl("knownGene",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        mm10 = {
            return(list(
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownGene=getUcscTblTpl("knownGene",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        rn5 = {
            return(list(
                mgcGenes=getUcscTblTpl("mgcGenes",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        rn6 = {
            return(list(
                mgcGenes=getUcscTblTpl("mgcGenes",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        dm3 = {
            return(list(
                flyBaseCanonical=
                    getUcscTblTpl("flyBaseCanonical",what),
                flyBaseGene=
                    getUcscTblTpl("flyBaseGene",what),
                flyBaseToRefSeq=
                    getUcscTblTpl("flyBaseToRefSeq",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what),
                refFlat=getUcscTblTpl("refFlat",what)
            ))
        },
        dm6 = {
            warning("No UCSC Genome annotation for Drosophila ",
                "melanogaster v6! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        danrer7 = {
            return(list(
                mgcGenes=getUcscTblTpl("mgcGenes",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        danrer10 = {
            return(list(
                mgcGenes=getUcscTblTpl("mgcGenes",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        pantro4 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        pantro5 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        #pantro6 = {
        #    warning("No UCSC Genome annotation for Pan ",
        #        "troglodytes! Will use RefSeq instead...",
        #        immediate.=TRUE)
        #    return(list(
        #        refFlat=getUcscTblTpl("refFlat",what),
        #        ensemblToGeneName=
        #            getUcscTblTpl("ensemblToGeneName",what),
        #        ensemblSource=
        #            getUcscTblTpl("ensemblSource",what)
        #    ))
        #},
        susscr3 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v3! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        susscr11 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v11! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        equcab2 = {
            warning("No UCSC Genome annotation for Equus ",
                "caballus v2! Will use RefSeq instead...",
                immediate.=TRUE)
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        }
    )
}

.getUcscTabledefRefseq <- function(org,what="queries") {
    switch(org,
        hg18 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what)
            ))
        },
        hg19 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        hg38 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what)
            ))
        },
        mm9 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        mm10 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                knownToRefSeq=
                    getUcscTblTpl("knownToRefSeq",what),
                knownCanonical=
                    getUcscTblTpl("knownCanonical",what),
                knownToEnsembl=
                    getUcscTblTpl("knownToEnsembl",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        rn5 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        rn6 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        dm3 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        dm6 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        danrer7 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        danrer10 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        pantro4 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        pantro5 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        #pantro6 = {
        #    return(list(
        #        refFlat=getUcscTblTpl("refFlat",what),
        #        ensemblToGeneName=
        #            getUcscTblTpl("ensemblToGeneName",what),
        #        ensemblSource=
        #            getUcscTblTpl("ensemblSource",what)
        #    ))
        #},
        susscr3 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        susscr11 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        },
        equcab2 = {
            return(list(
                refFlat=getUcscTblTpl("refFlat",what),
                ensemblToGeneName=
                    getUcscTblTpl("ensemblToGeneName",what),
                ensemblSource=
                    getUcscTblTpl("ensemblSource",what)
            ))
        }
    )
}

.getUcscQueryUcscGene <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg19 = {
            return(paste("SELECT knownCanonical.chrom AS `chromosome`,",
				"`chromStart` AS `start`,",
				"`chromEnd` AS `end`,",
				"`transcript` AS `gene_id`,",
				"0 AS `gc_content`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownCanonical` INNER JOIN `knownGene`",
				"ON knownCanonical.transcript=knownGene.name",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY `gene_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should have been like hg19 but it's like hg18
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT knownCanonical.chrom AS `chromosome`,",
				"`chromStart` AS `start`,",
				"`chromEnd` AS `end`,",
				"`transcript` AS `gene_id`,",
				"0 AS `gc_content`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownCanonical` INNER JOIN `knownGene`",
				"ON knownCanonical.transcript=knownGene.name",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY `gene_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT knownCanonical.chrom AS `chromosome`,",
				"`chromStart` AS `start`,",
				"`chromEnd` AS `end`,",
				"`transcript` AS `gene_id`,",
				"0 AS `gc_content`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownCanonical` INNER JOIN `knownGene`", 
				"ON knownCanonical.transcript=knownGene.name",
				"INNER JOIN `knownToRefSeq`", 
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY `gene_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        rn6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm3 = {
            return(paste("SELECT flyBaseCanonical.chrom AS `chromosome`,",
				"`chromStart` AS `start`,",
				"`chromEnd` AS `end`,",
				"`transcript` AS `gene_id`,",
				"0 AS `gc_content`,",
				"flyBaseGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `flyBaseCanonical` INNER JOIN `flyBaseGene`", 
				"ON flyBaseCanonical.transcript=flyBaseGene.name",
				"INNER JOIN `flyBaseToRefSeq`",
				"ON flyBaseCanonical.transcript=flyBaseToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON flyBaseToRefSeq.value=refFlat.name",
				"INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm6 = {
            warning("No UCSC Genome annotation for Drosophila ",
                "melanogaster v6! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`", 
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer7 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
			"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
			"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
			"mgcGenes.chrom AS `chromosome`,",
			"`txStart` AS `start`,",
			"`txEnd` AS `end`,",
			"mgcGenes.name AS `gene_id`,",
			"0 AS `gc_content`,",
			"mgcGenes.strand AS `strand`,",
			"`name2` AS `gene_name`,",
			"`source` AS `biotype`",
			"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
			"ON mgcGenes.name2=ensemblToGeneName.value",
			"INNER JOIN `ensemblSource`",
			"ON ensemblToGeneName.name=ensemblSource.name",
			"GROUP BY `gene_name`",
			"ORDER BY `chromosome`,`start`)"))
        },
        danrer10 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro4 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v4! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro5 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v5! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr3 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v3! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr11 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v11! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        equcab2 = {
            warning("No UCSC Genome annotation for Equus ",
                "caballus v2! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        }
    )
}

.getUcscQueryRefseqGene <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT  refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg19 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should be the same as hg19 but is as hg18
            return(paste("SELECT  refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        rn6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm3 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer7 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer10 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro4 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro5 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr3 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        equcab2 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`gene_id`,",
				"`gc_content`,`strand`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        }
    )
}

.getUcscQueryUcscExon <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `transcript_id`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`" ,
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`" ,
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY `exon_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        hg19 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.exonStarts AS `start`,",
				"knownGene.exonEnds AS `end`,",
				"knownGene.name AS `exon_id`,",
				"knownGene.strand AS `strand`,",
				"`transcript` AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownCanonical`", 
				"ON knownGene.name=knownCanonical.transcript",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`", 
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should be the same as hg19 but is as hg18
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.exonStarts AS `start`,",
				"knownGene.exonEnds AS `end`,",
				"knownGene.name AS `exon_id`,",
				"knownGene.strand AS `strand`,",
				"`transcript` AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownCanonical`", 
				"ON knownGene.name=knownCanonical.transcript",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.exonStarts AS `start`,",
				"knownGene.exonEnds AS `end`,",
				"knownGene.name AS `exon_id`,",
				"knownGene.strand AS `strand`,",
				"`transcript` AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownCanonical`",
				"ON knownGene.name=knownCanonical.transcript",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.exonStarts AS `start`,",
				"knownGene.exonEnds AS `end`,",
				"knownGene.name AS `exon_id`,",
				"knownGene.strand AS `strand`,",
				"`transcript` AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownCanonical`",
				"ON knownGene.name=knownCanonical.transcript",
				"INNER JOIN `knownToRefSeq`",
				"ON knownCanonical.transcript=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`exonStarts` AS `start`,",
				"`exonEnds` AS `end`,",
				"mgcGenes.name AS `exon_id`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `gene_id`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        rn6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`exonStarts` AS `start`,",
				"`exonEnds` AS `end`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `exon_id`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `gene_id`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`" ,
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm3 = {
            return(paste("SELECT flyBaseCanonical.chrom AS `chromosome`,",
				"flyBaseGene.exonStarts AS `start`,",
				"flyBaseGene.exonEnds AS `end`,",
				"`transcript` AS `exon_id`,",
				"flyBaseGene.strand AS `strand`,",
				"`transcript` AS `transcript_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `flyBaseCanonical` INNER JOIN `flyBaseGene` ON",
				"flyBaseCanonical.transcript=flyBaseGene.name",
				"INNER JOIN `flyBaseToRefSeq`",
				"ON flyBaseCanonical.transcript=flyBaseToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON flyBaseToRefSeq.value=refFlat.name",
				"INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm6 = {
            warning("No UCSC Genome annotation for Drosophila ",
                "melanogaster v6! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer7 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`exonStarts` AS `start`,",
				"`exonEnds` AS `end`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `exon_id`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `gene_id`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer10 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`exonStarts` AS `start`,",
				"`exonEnds` AS `end`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `exon_id`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `gene_id`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"mgcGenes.chrom AS `chromosome`,",
				"`exonStarts` AS `start`,",
				"`exonEnds` AS `end`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `exon_id`,",
				"mgcGenes.strand AS `strand`,",
				"mgcGenes.name AS `gene_id`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro4 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v4! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro5 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v5! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr3 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v3! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr11 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v11! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        equcab2 = {
            warning("No UCSC Genome annotation for Equus ",
                "caballus v11! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        }
    )
}

.getUcscQueryRefseqExon <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT  refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds  AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg19 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds  AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`", 
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should be the same as hg19 but is as hg18
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds  AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`", 
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds  AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds  AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`", 
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"GROUP BY refFlat.name",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        rn6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm3 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        dm6 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer7 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer10 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        danrer11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro4 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        pantro5 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr3 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        susscr11 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        },
        equcab2 = {
            return(paste("SELECT `chromosome`,`start`,`end`,`exon_id`,",
				"`strand`,`gene_id`,`gene_name`,`biotype` FROM",
				"(SELECT MAX(`txEnd` - `txStart`) AS `width`,",
				"refFlat.chrom AS `chromosome`,",
				"refFlat.exonStarts AS `start`,",
				"refFlat.exonEnds AS `end`,",
				"refFlat.name AS `exon_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name", 
				"GROUP BY `gene_name`",
				"ORDER BY `chromosome`,`start`)"))
        }
    )
}

.getUcscQueryUcscTranscript <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.txStart AS `start`,",
				"knownGene.txEnd AS `end`,",
				"knownGene.name AS `transcript_id`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownToRefSeq`", 
				"ON knownGene.name=knownToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg19 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.txStart AS `start`,",
				"knownGene.txEnd AS `end`,",
				"knownGene.name AS `transcript_id`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownToRefSeq`", 
				"ON knownGene.name=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownGene.name=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat` ON", 
				"knownToRefSeq.value=refFlat.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should be the same as hg19 but is like hg18
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.txStart AS `start`,",
				"knownGene.txEnd AS `end`,",
				"knownGene.name AS `transcript_id`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownToRefSeq`",
				"ON knownGene.name=knownToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON knownToRefSeq.value=refFlat.name",
				"GROUP BY knownGene.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.txStart AS `start`,",
				"knownGene.txEnd AS `end`,",
				"knownGene.name AS `transcript_id`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownToRefSeq`",
				"ON knownGene.name=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownGene.name=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat` ON",
				"knownToRefSeq.value=refFlat.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT knownGene.chrom AS `chromosome`,",
				"knownGene.txStart AS `start`,",
				"knownGene.txEnd AS `end`,",
				"knownGene.name AS `transcript_id`,",
				"knownGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `knownGene` INNER JOIN `knownToRefSeq`", 
				"ON knownGene.name=knownToRefSeq.name",
				"INNER JOIN `knownToEnsembl`",
				"ON knownGene.name=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"INNER JOIN `refFlat` ON",
				"knownToRefSeq.value=refFlat.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `transcript_id`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        rn6 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `transcript_id`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`",
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm3 = {
            return(paste("SELECT flyBaseGene.chrom AS `chromosome`,",
				"flyBaseGene.txStart AS `start`,",
				"flyBaseGene.txEnd AS `end`,",
				"flyBaseGene.name AS `transcript_id`,",
				"flyBaseGene.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `flyBaseGene` INNER JOIN `flyBaseToRefSeq`",
				"ON flyBaseGene.name=flyBaseToRefSeq.name",
				"INNER JOIN `refFlat`",
				"ON flyBaseToRefSeq.value=refFlat.name",
				"INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm6 = {
            warning("No UCSC Genome annotation for Drosophila ",
                "melanogaster v6! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer7 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `transcript_id`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer10 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `transcript_id`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer11 = {
            return(paste("SELECT mgcGenes.chrom AS `chromosome`,",
				"`txStart` AS `start`,",
				"`txEnd` AS `end`,",
				"mgcGenes.name AS `transcript_id`,",
				"mgcGenes.strand AS `strand`,",
				"`name2` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `mgcGenes` INNER JOIN `ensemblToGeneName`", 
				"ON mgcGenes.name2=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        pantro4 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v4! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        pantro5 = {
            warning("No UCSC Genome annotation for Pan ",
                "troglodytes v5! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        susscr3 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v3! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        susscr11 = {
            warning("No UCSC Genome annotation for Sus ",
                "scrofa v11! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        equcab2 = {
            warning("No UCSC Genome annotation for Equus ",
                "caballus v2! Will use RefSeq instead...",
                immediate.=TRUE)
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`",
				"ON ensemblToGeneName.value=refFlat.geneName",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        }
    )
}

.getUcscQueryRefseqTranscript <- function(org) {
    switch(org,
        hg18 = {
            return(paste("SELECT  refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"refFlat.name AS `gene_id`",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"ORDER BY `chromosome`, `start`"))
        },
        hg19 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"ORDER BY `chromosome`,`start`"))
        },
        hg38 = {
			# Should be the same as hg19 but is as hg18
            return(paste("SELECT  refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `gene_id`,",
				"0 AS `gc_content`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"'NA' AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"ORDER BY `chromosome`,`start`"))
        },
        mm9 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"ORDER BY `chromosome`,`start`"))
        },
        mm10 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `knownToRefSeq`",
				"ON refFlat.name=knownToRefSeq.value",
				"INNER JOIN `knownCanonical`",
				"ON knownToRefSeq.name=knownCanonical.transcript",
				"INNER JOIN `knownToEnsembl`",
				"ON knownCanonical.transcript=knownToEnsembl.name",
				"INNER JOIN `ensemblSource`",
				"ON knownToEnsembl.value=ensemblSource.name",
				"ORDER BY `chromosome`,`start`"))
        },
        rn5 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        rn6 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm3 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        dm6 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer7 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer10 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        danrer11 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        pantro4 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        pantro5 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        susscr3 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        susscr11 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        },
        equcab2 = {
            return(paste("SELECT refFlat.chrom AS `chromosome`,",
				"refFlat.txStart AS `start`,",
				"refFlat.txEnd AS `end`,",
				"refFlat.name AS `transcript_id`,",
				"refFlat.strand AS `strand`,",
				"`geneName` AS `gene_name`,",
				"`source` AS `biotype`",
				"FROM `refFlat` INNER JOIN `ensemblToGeneName`", 
				"ON refFlat.geneName=ensemblToGeneName.value",
				"INNER JOIN `ensemblSource`",
				"ON ensemblToGeneName.name=ensemblSource.name",
				"GROUP BY `transcript_id`",
				"ORDER BY `chromosome`,`start`"))
        }
    )
}

#' Return host, username and password for UCSC Genome Browser database
#'
#' Returns a character vector with a hostname, username and password to connect
#' to the UCSC Genome Browser database to retrieve annotation. Internal use.
#'
#' @return A named character vector.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.creds <- get.ucsc.credentials()
#'}
.getUcscCredentials <- function() {
    return(c(
        host="genome-mysql.cse.ucsc.edu",
        user="genome",
        password=""
    ))
}

..getAllUcsc <- function(to) {
	if (missing(to))
		stop("Please provide a path to put the sqlite databases!")
	if (!dir.exists(to))
		dir.create(to,recursive=TRUE)
		
	for (refdb in getSupportedUcscDbs())
		for (org in getSupportedOrganisms())
			..getUcscSqlite(org,refdb,to)
}

..getUcscSqlite <- function(org,refdb,to) {
	if (missing(to))
		stop("Please provide a path to put the sqlite databases!")
	if (!dir.exists(to))
		dir.create(to,recursive=TRUE)
	
	message("RETRIEVING ",org," FROM ",refdb)
	dbTmp <- getUcscDbl(org=org,refdb=refdb)
	name <- paste(org,"_",refdb,".sqlite",sep="")
	file.copy(dbTmp,file.path(to,name),recursive=TRUE)
	return()
}
