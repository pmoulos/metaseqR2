# The new function with the new weights.
getWeights <- function(org=c("human","chimpanzee","mouse","fruitfly",
    "arabidopsis","rat")) {
    org <- tolower(org)
    checkTextArgs("org",org,c("human","chimpanzee","mouse","fruitfly",
        "arabidopsis","rat"))
    switch(org,
        human = {
            return(c(
                deseq=0.03077317,
                deseq2=0.16791844,
                edger=0.07695718,
                limma=0.18600075,
                nbpseq=0.03315111,
                noiseq=0.05058633,
                bayseq=0.15290457,
                absseq=0.12843044,
                dss=0.17327802
            ))
        },
        chimpanzee = {
            return(c(
                deseq=0.03345859,
                deseq2=0.17691534,
                edger=0.08214923,
                limma=0.18127409,
                nbpseq=0.03209271,
                noiseq=0.05261594,
                bayseq=0.15328300,
                absseq=0.13449179,
                dss=0.15371930
            ))
        },
        mouse = {
            return(c(
                deseq=0.02707875,
                deseq2=0.16434559,
                edger=0.12614496,
                limma=0.17038977,
                nbpseq=0.02614409,
                noiseq=0.02064747,
                bayseq=0.15508511,
                absseq=0.14800242,
                dss=0.16216184
            ))
        },
        fruitfly = {
            return(c(
                deseq=0.006526282,
                deseq2=0.172079991,
                edger=0.064114940,
                limma=0.193692581,
                nbpseq=0.005907387,
                noiseq=0.019310076,
                bayseq=0.177483603,
                absseq=0.166022212,
                dss=0.194862929
            ))
        },
        arabidopsis = {
            return(c(
                deseq=0.02745255,
                deseq2=0.18604304,
                edger=0.05530561,
                limma=0.20806899,
                nbpseq=0.02567801,
                noiseq=0.04377844,
                bayseq=0.16104719,
                absseq=0.09621622,
                dss=0.19640995
            ))
        },
        rat = {
            return(c(
                deseq=0.03409202,
                deseq2=0.15203180,
                edger=0.14542665,
                limma=0.15327172,
                nbpseq=0.03501810,
                noiseq=0.01954866,
                bayseq=0.15587551,
                absseq=0.15289613,
                dss=0.15183941
            ))
        }
    )
}

getDefaults <- function(what,method=NULL) {
    if (what %in% c("normalization","statistics") && is.null(method))
        stopwrap("The method argument must be provided when what is ",
            "\"normalization\" or \"statistics\"!")
    switch(what,
        normalization = {
            switch(method,
                edaseq = {
                    return(list(within.which="loess",between.which="full"))
                },
                deseq = {
                    return(list(locfunc=median))
                },
                deseq2 = {
                    return(list(
                        tidy=FALSE, # DESeqDataSet object creation args
                        type="ratio",locfunc=stats::median
                    ))
                },                
                edger = {
                    return(list(
                        method="TMM",refColumn=NULL,logratioTrim=0.3,
                        sumTrim=0.05,doWeighting=TRUE,Acutoff=-1e10,p=0.75
                    ))
                },
                noiseq = {
                    return(list(
                        method="tmm", # which normalization
                        long=1000,lc=1,k=1, # common arguments
                        refColumn=1,logratioTrim=0.3,sumTrim=0.05,
                        # TMM normalization arguments
                        doWeighting=TRUE,Acutoff=-1e+10 
                    ))
                },
                nbpseq = {
                    return(list(main.method="nbsmyth",method="AH2010",
                        thinning=TRUE))
                },
                absseq = {
                    # method for estimating sizeFactors
                    return(list(normMethod="qtotal"))
                },
                dss = {
                    return(list(method="lr"))                        
                }
            )
        },           
        statistics = {
            switch(method,
                deseq = {
                    return(list(method="blind",sharingMode="fit-only",
                        fitType="local"))
                },
                deseq2 = {
                    # If sth is changed from here, change it also from the 
                    # switch beggining at line 804
                    return(list(      
                        # DESeqDataSetFromMatrix default args
                        tidy=FALSE,
                        # sizeFactors default args won't be passed because user 
                        # doesn't have to change anything at this point
                        # estimateDispersions default args
                        fitType= "parametric", maxit= 100, quiet= FALSE, 
                        modelMatrix= NULL, 
                        #minmu= 0.5, #removed because it caused 
                        # "Error ...bla-bla... unused argument minmu= 0.5"
                        #nbinomWaldTest default args
                        betaPrior= FALSE,
                        #betaPriorVar, modelMatrixType, 
                        betaTol=1e-08,useOptim=TRUE,
                        useT=FALSE,useQR=TRUE,
                        #results default args
                        lfcThreshold=0,altHypothesis="greaterAbs",
                        independentFiltering=TRUE, alpha=0.1,pAdjustMethod="BH",
                        format= "DataFrame",addMLE= FALSE,parallel=FALSE#,
                        #BPPARAM=bpparam()
                        #nbinomLRT default args are all passed above
                    ))
                },  
                edger = {
                    return(list(
                        main.method="classic", # classic or glm fit
                        rowsum.filter=5,prior.df=10,
                        trend="movingave",span=NULL,
                        # classic estimateCommonDisp arguments
                        tag.method="grid",grid.length=11,grid.range=c(-6,6),
                        # classic estimateTagwiseDisp arguments
                        offset=NULL,glm.method="CoxReid",subset=10000,
                        # glm estimateGLMCommonDisp and estimateGLMTrendedDisp 
                        # arguments
                        AveLogCPM=NULL,trend.method="auto",
                        # glm estimateGLMTagwiseDisp arguments
                        dispersion=NULL,offset=NULL,weights=NULL, 
                        # glmFit arguments
                        lib.size=NULL,prior.count=0.125,start=NULL,
                        method="auto",test="chisq",
                        # glmLRT arguments
                        abundance.trend=TRUE,robust=FALSE,
                        winsor.tail.p=c(0.05,0.1) # glmLFTest arguments
                    ))
                },
                noiseq = {
                    return(list(
                        k=0.5,norm="n",replicates="biological",
                        factor="class",conditions=NULL,pnr=0.2,
                        nss=5,v=0.02,lc=1, 
                        # noiseq general and specific arguments
                        nclust=15,r=100,adj=1.5,
                        a0per=0.9,filter=0,depth=NULL,        
                        cv.cutoff=500,cpm=1 # noiseqbio specific arguments
                        
                    ))
                },
                bayseq = {
                    return(list(samplesize=10000,samplingSubset=NULL,
                        equalDispersions=TRUE,estimation="QL",zeroML=FALSE,
                        consensus=FALSE,moderate=TRUE,pET="BIC",
                        marginalise=FALSE,subset=NULL,priorSubset=NULL,
                        bootStraps=1,conv=1e-4,nullData=FALSE,returnAll=FALSE,
                        returnPD=FALSE,discardSampling=FALSE,cl=NULL))
                },
                limma = {
                    return(list(normalize.method="none"))
                },
                nbpseq = {
                    return(list(
                        main.method="nbsmyth",
                        model=list(nbpseq="log-linear-rel-mean",nbsmyth="NBP"),
                        tests="HOA",
                        alternative="two.sided"
                    ))
                },
                absseq = {
                    # if sth is changed from here, change it also from the
                    # switch beggining at line 802
                   return(list(
                        #ABSDataSet default args
                        paired = FALSE, minDispersion = NULL, minRates = 0.1,
                        maxRates = 0.3, LevelstoNormFC = 100,
                        #ABSSeq default args
                        adjmethod = "BH", replaceOutliers = TRUE, 
                        useaFold = FALSE, quiet = FALSE,
                        #ABSSeqlm default args
                        lmodel = TRUE, preval = 0.05,
                        qforkappa = 0, scale = FALSE
                    ))
                },
                dss = {
                    # if sth is changed from here, change it also from the 
                    # switch beggining at line 802
                    return(list(
                        #estDispersion default args
                        trend=FALSE,
                        #waldTest default args
                        equal.var=FALSE
                    ))
                }
            )
        },
        utrOpts = {
            return(list(
                frac=1,
                minLength=300,
                downstram=50
            ))
        },
        geneFilter = {
            return(list(
                length=list(
                    length=500
                ),
                avgReads=list(
                    average.per.bp=100,
                    quantile=0.75
                ),
                expression=list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA,
                    custom=NA
                ),
                biotype=getDefaults("biotypeFilter",method[1]),
                presence=list(
                    frac=0.25,
                    min.count=10,
                    per.condition=FALSE
                )
            ))
        },
        exonFilter = {
            return(list(
                mnrpx=list(
                    exons.per.gene=5,
                    min.exons=2,
                    frac=1/5
                )
            ))
        },
        biotypeFilter = {
            switch(method,
                hg18 = {
                    return(list(
                        unprocessed_pseudogene=TRUE,
                        pseudogene=FALSE,
                        miRNA=FALSE,
                        retrotransposed=FALSE,
                        protein_coding=FALSE,
                        processed_pseudogene=FALSE,
                        snRNA=FALSE,
                        snRNA_pseudogene=TRUE,
                        Mt_tRNA_pseudogene=TRUE,
                        miRNA_pseudogene=TRUE,
                        misc_RNA=FALSE,
                        tRNA_pseudogene=TRUE,
                        snoRNA=FALSE,
                        scRNA_pseudogene=TRUE,
                        rRNA_pseudogene=TRUE,
                        snoRNA_pseudogene=TRUE,
                        rRNA=TRUE,
                        misc_RNA_pseudogene=TRUE,
                        IG_V_gene=FALSE,
                        IG_D_gene=FALSE,
                        IG_J_gene=FALSE,
                        IG_C_gene=FALSE,
                        IG_pseudogene=TRUE,
                        scRNA=FALSE
                    ))
                },
                hg19 = {
                    return(list(
                        pseudogene=FALSE,
                        lincRNA=FALSE,
                        protein_coding=FALSE,
                        antisense=FALSE,
                        processed_transcript=FALSE,
                        snRNA=FALSE,
                        sense_intronic=FALSE,
                        miRNA=FALSE,
                        misc_RNA=FALSE,
                        snoRNA=FALSE,
                        rRNA=TRUE,
                        polymorphic_pseudogene=FALSE,
                        sense_overlapping=FALSE,
                        three_prime_overlapping_ncrna=FALSE,
                        TR_V_gene=FALSE,
                        TR_V_pseudogene=TRUE,
                        TR_D_gene=FALSE,
                        TR_J_gene=FALSE,
                        TR_C_gene=FALSE,
                        TR_J_pseudogene=TRUE,
                        IG_C_gene=FALSE,
                        IG_C_pseudogene=TRUE,
                        IG_J_gene=FALSE,
                        IG_J_pseudogene=TRUE,
                        IG_D_gene=FALSE,
                        IG_V_gene=FALSE,
                        IG_V_pseudogene=TRUE
                    ))
                },
                hg38 = {
                    return(list(
                        protein_coding=FALSE,
                        polymorphic_pseudogene=FALSE,
                        lincRNA=FALSE,
                        unprocessed_pseudogene=TRUE,
                        processed_pseudogene=FALSE,
                        antisense=FALSE,
                        processed_transcript=FALSE,
                        transcribed_unprocessed_pseudogene=FALSE,
                        sense_intronic=FALSE,
                        unitary_pseudogene=TRUE,
                        IG_V_gene=FALSE,
                        IG_V_pseudogene=TRUE,
                        TR_V_gene=FALSE,
                        sense_overlapping=FALSE,
                        transcribed_processed_pseudogene=FALSE,
                        miRNA=FALSE,
                        snRNA=FALSE,
                        misc_RNA=FALSE,
                        rRNA=TRUE,
                        snoRNA=FALSE,
                        IG_J_pseudogene=TRUE,
                        IG_J_gene=FALSE,
                        IG_D_gene=FALSE,
                        three_prime_overlapping_ncrna=FALSE,
                        IG_C_gene=FALSE,
                        IG_C_pseudogene=TRUE,
                        pseudogene=TRUE,
                        TR_V_pseudogene=TRUE,
                        Mt_tRNA=TRUE,
                        Mt_rRNA=TRUE,
                        translated_processed_pseudogene=FALSE,
                        TR_J_gene=FALSE,
                        TR_C_gene=FALSE,
                        TR_D_gene=FALSE,
                        TR_J_pseudogene=TRUE,
                        LRG_gene=FALSE
                    ))
                },
                mm9 = {
                    return(list(
                        pseudogene=FALSE,
                        snRNA=FALSE,
                        protein_coding=FALSE,
                        antisense=FALSE,
                        miRNA=FALSE,
                        lincRNA=FALSE,
                        snoRNA=FALSE,
                        processed_transcript=FALSE,
                        misc_RNA=FALSE,
                        rRNA=TRUE,
                        sense_overlapping=FALSE,
                        sense_intronic=FALSE,
                        polymorphic_pseudogene=FALSE,
                        non_coding=FALSE,
                        three_prime_overlapping_ncrna=FALSE,
                        IG_C_gene=FALSE,
                        IG_J_gene=FALSE,
                        IG_D_gene=FALSE,
                        IG_V_gene=FALSE,
                        ncrna_host=FALSE
                    ))
                },
                mm10 = {
                    return(list(
                        pseudogene=FALSE,
                        snRNA=FALSE,
                        protein_coding=FALSE,
                        antisense=FALSE,
                        miRNA=FALSE,
                        snoRNA=FALSE,
                        lincRNA=FALSE,
                        processed_transcript=FALSE,
                        misc_RNA=FALSE,
                        rRNA=TRUE,
                        sense_intronic=FALSE,
                        sense_overlapping=FALSE,
                        polymorphic_pseudogene=FALSE,
                        IG_C_gene=FALSE,
                        IG_J_gene=FALSE,
                        IG_D_gene=FALSE,
                        IG_LV_gene=FALSE,
                        IG_V_gene=FALSE,
                        IG_V_pseudogene=TRUE,
                        TR_V_gene=FALSE,
                        TR_V_pseudogene=TRUE,
                        three_prime_overlapping_ncrna=FALSE
                    ))
                },
                dm3 = {
                    return(list(
                        protein_coding=FALSE,
                        ncRNA=FALSE,
                        snoRNA=FALSE,
                        pre_miRNA=FALSE,
                        pseudogene=FALSE,
                        snRNA=FALSE,
                        tRNA=FALSE,
                        rRNA=TRUE
                    ))
                },
                rn5 = {
                    return(list(
                        protein_coding=FALSE,
                        pseudogene=FALSE,
                        processed_pseudogene=FALSE,
                        miRNA=FALSE,
                        rRNA=TRUE,
                        misc_RNA=FALSE
                    ))
                },
                rn6 = {
                    return(list(
                        antisense=FALSE , 
                        lincRNA=FALSE, 
                        miRNA=FALSE, 
                        misc_RNA=FALSE, 
                        Mt_rRNA=TRUE, 
                        Mt_tRNA=TRUE, 
                        processed_pseudogene=FALSE, 
                        processed_transcript=FALSE, 
                        protein_coding=FALSE, 
                        pseudogene=FALSE, 
                        ribozyme=FALSE, 
                        rRNA=TRUE, 
                        scaRNA=FALSE, 
                        sense_intronic=FALSE, 
                        snoRNA=FALSE, 
                        snRNA=FALSE, 
                        sRNA=FALSE, 
                        TEC=FALSE, 
                        transcribed_processed_pseudogene=FALSE, 
                        transcribed_unprocessed_pseudogene=FALSE, 
                        unprocessed_pseudogene=TRUE
                    ))
                },
                danrer7 = {
                    return(list(
                        antisense=FALSE,
                        protein_coding=FALSE,
                        miRNA=FALSE,
                        snoRNA=FALSE,
                        rRNA=TRUE,
                        lincRNA=FALSE,
                        processed_transcript=FALSE,
                        snRNA=FALSE,
                        pseudogene=FALSE,
                        sense_intronic=FALSE,
                        misc_RNA=FALSE,
                        polymorphic_pseudogene=FALSE,
                        IG_V_pseudogene=TRUE,
                        IG_C_pseudogene=TRUE,
                        IG_J_pseudogene=TRUE,
                        non_coding=FALSE,
                        sense_overlapping=FALSE
                    ))
                },
                pantro4 = {
                    return(list(
                        protein_coding=FALSE,
                        pseudogene=FALSE,
                        processed_pseudogene=FALSE,
                        miRNA=FALSE,
                        rRNA=TRUE,
                        snRNA=FALSE,
                        snoRNA=FALSE,
                        misc_RNA=FALSE
                    ))
                },
                susscr3 = {
                    return(list(
                        antisense=FALSE,
                        protein_coding=FALSE,
                        lincRNA=FALSE,
                        pseudogene=FALSE,
                        processed_transcript=FALSE,
                        miRNA=FALSE,
                        rRNA=TRUE,
                        snRNA=FALSE,
                        snoRNA=FALSE,
                        misc_RNA=FALSE,
                        non_coding=FALSE,
                        IG_C_gene=FALSE,
                        IG_J_gene=FALSE,
                        IG_V_gene=FALSE,
                        IG_V_pseudogene=TRUE
                    ))
                },
                tair10 = {
                    return(list(
                        miRNA=FALSE,
                        ncRNA=FALSE,
                        protein_coding=FALSE,
                        pseudogene=FALSE,
                        rRNA=TRUE,
                        snoRNA=FALSE,
                        snRNA=FALSE,
                        transposable_element=FALSE,
                        tRNA=FALSE
                    ))
                },
                equcab2 = {
                    return(list(
                        miRNA=FALSE,
                        misc_RNA=FALSE,
                        protein_coding=FALSE,
                        pseudogene=FALSE,
                        processed_pseudogene=FALSE,
                        rRNA=TRUE,
                        snoRNA=FALSE,
                        snRNA=FALSE
                    ))
                }
            )
        }
    )
}

validateAlgArgs <- function(normalization,statistics,normArgs,statArgs) {
    if (normalization=="each") {
        if (!is.null(normArgs)) {
            for (s in statistics) {
                if (!is.null(normArgs[[s]])) {
                    switch(s,
                        deseq = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("normalization",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        deseq2 = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("normalization",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        edger = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        limma = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics","edger",tmp)
                            normArgs[[s]] <- getDefaults("normalization",
                                "edger")
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        nbpseq = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        noiseq = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        bayseq = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics","edger",tmp)
                            normArgs[[s]] <- getDefaults("normalization",
                                "edger")
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        },
                        absseq = {
                            #the arguments' list passed by the user
                            tmp <- normArgs[[s]]
                            #it validates the arguments given
                            tmp <- validateListArgs("statistics",s,tmp) 
                            #it gets defaults and thus ignores everything else
                            normArgs[[s]] <- getDefaults("normalization", 
                                s) 
                            if (length(tmp)>0)
                                # it changes default values of valid arguments 
                                # to user-specified 
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp) 
                        },
                        dss = {
                            tmp <- normArgs[[s]]
                            tmp <- validateListArgs("statistics",s,tmp)
                            normArgs[[s]] <- getDefaults("normalization",s)
                            if (length(tmp)>0)
                                normArgs[[s]] <- setArg(normArgs[[s]],tmp)
                        }                        
                    )
                }
                else {
                    switch(s,
                        deseq = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "deseq")
                        },
                        deseq2 = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "deseq2")
                        },
                        edger = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "edger")
                        },
                        limma = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "edger")
                        },
                        nbpseq = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "nbpseq")
                        },
                        noiseq = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "noiseq")
                        },
                        bayseq = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "edger")
                        },    
                        absseq = {
                            normArgs[[s]] <- getDefaults(normalization,
                                "absseq")
                        },
                        dss = {
                            normArgs[[s]] <- getDefaults(normalization,"dss")
                        }
                    )
                }
            }
        }
        else {
            normArgs <- vector("list",length(statistics))
            names(normArgs) <- statistics
            for (s in statistics) {
                switch(s,
                    deseq = {
                        normArgs[[s]] <- getDefaults(normalization,"deseq")
                    },
                    deseq2 = {
                        normArgs[[s]] <- getDefaults(normalization,"deseq2")
                    },    
                    edger = {
                        normArgs[[s]] <- getDefaults(normalization,"edger")
                    },
                    limma = {
                        normArgs[[s]] <- getDefaults(normalization,"edger")
                    },
                    nbpseq = {
                        normArgs[[s]] <- getDefaults(normalization,"nbpseq")
                    },
                    noiseq = {
                        normArgs[[s]] <- getDefaults(normalization,"noiseq")
                    },
                    bayseq = {
                        normArgs[[s]] <- getDefaults(normalization,"edger")
                    },
                    absseq = {
                        normArgs[[s]] <- getDefaults(normalization,"absseq")
                    },   
                    dss = {
                        norm.args[[s]] <- getDefaults(normalization,"dss")
                    }
                )
            }
        }
    }
    else {
        if (!is.null(normArgs))
        {
            tmp <- normArgs
            tmp <- validateListArgs("normalization",normalization,tmp)
            normArgs <- getDefaults("normalization",normalization)
            if (length(tmp)>0)
                normArgs <- setArg(normArgs,tmp)
        }
        else
            normArgs <- getDefaults("normalization",normalization)
    }
    for (s in statistics) {
        if (!is.null(statArgs[[s]])) {
            tmp <- statArgs[[s]]
            tmp <- validateListArgs("statistics",s,tmp)
            statArgs[[s]] <- getDefaults("statistics",s)
            if (length(tmp)>0)
                statArgs[[s]] <- setArg(statArgs[[s]],tmp)
        }
        else
            statArgs[[s]] <- getDefaults("statistics",s)
    }
    return(list(normArgs=normArgs,statArgs=statArgs))
}

validateListArgs <- function(what,method=NULL,argList) {
    what <- tolower(what)
    checkTextArgs("what",what,c("normalization","statistics","geneFilter",
        "exonFilter","biotype.filter"))
    if (what %in% c("normalization","statistics") && is.null(method))
        stopwrap("The method argument must be provided when what is ",
            "\"normalization\" or \"statistics\"!")
    switch(what,
        normalization = {
            switch(method,
                edaseq = { 
                    valid <- names(argList) %in% c("within.which",
                        "between.which")
                    not.valid <- which(!valid)
                },
                deseq = {
                    valid <- names(argList) %in% c("locfunc")
                    not.valid <- which(!valid) 
                },
                deseq2 = {
                    valid <- names(argList) %in% c("tidy","type","locfunc")
                    not.valid <- which(!valid) 
                },
                edger = {
                    valid <- names(argList) %in% c("method","refColumn",
                        "logratioTrim","sumTrim","doWeighting","Acutoff","p")
                    not.valid <- which(!valid)
                },
                noiseq = {
                    valid <- names(argList) %in% c("method","long","lc","k",
                        "refColumn","logratioTrim","sumTrim","doWeighting",
                        "Acutoff")
                    not.valid <- which(!valid)
                },
                nbpseq = {
                    valid <- names(argList) %in% c("main.method","method",
                        "thinning")
                    not.valid <- which(!valid)
                },
                absseq = {
                    valid <- names(argList) %in% c("normMethod")
                    not.valid <- which(!valid)
                },
                dss = {
                    valid <- names(argList) %in% c("method")
                    notValid <- which(!valid)
                }
            )
            if (length(not.valid)>0) {
                warnwrap(paste("The following",method,what,"argument names",
                    "are invalid and will be ignored:",
                    paste(names(argList)[not.valid],collapse=", ")))
                argList[not.valid] <- NULL
            }
            return(argList)
        },
        statistics = {
            switch(method,
                deseq = {
                    valid <- names(argList) %in% c("method","sharingMode",
                        "fitType")
                    not.valid <- which(!valid)
                },
                deseq2 = {
                    valid <- names(argList) %in% c("tidy","fitType","maxit",
                        "quiet","modelMatrix","minmu","betaPrior",
                        #"betaPriorVar","modelMatrixType", #removed because they
                        # are also removed in get defaults
                        "betaTol","useOptim","useT","df","useQR","lfcThreshold",
                        "altHypothesis","independentFiltering","alpha",
                        "pAdjustMethod","format","addMLE","parallel","BPPARAM")
                    not.valid <- which(!valid)
                },
                edger = {
                    valid <- names(argList) %in% c("main.method",
                        "rowsum.filter","prior.df","trend","span","tag.method",
                        "grid.length","grid.range","offset","glm.method",
                        "subset","AveLogCPM","trend.method","dispersion",
                        "offset","weights","lib.size","prior.count","start",
                        "method","abundance.trend","robust",
                        "winsor.tail.p")
                    not.valid <- which(!valid)
                },
                noiseq = {
                    valid <- names(argList) %in% c("k","norm","replicates",
                        "factor","conditions","pnr","nss","v","lc","nclust","r",
                        "adj","a0per","filter","depth","cv.cutoff","cpm")
                    not.valid <- which(!valid)
                },
                bayseq = {
                    valid <- names(argList) %in% c("samplesize",
                        "samplingSubset","equalDispersions","estimation",
                        "zeroML","consensus","moderate","pET","marginalise",
                        "subset","priorSubset","bootStraps","conv","nullData",
                        "returnAll","returnPD","discardSampling","cl")
                    not.valid <- which(!valid)
                },
                limma = {
                    valid <- names(argList) %in% c("normalize.method")
                    not.valid <- which(!valid)
                },
                nbpseq = {
                    valid <- names(argList) %in% c("main.method","method",
                        "tests","alternative")
                    not.valid <- which(!valid)
                },
                absseq = {
                    valid <- names(argList) %in% c("paired","minDispersion",
                        "minRates","maxRates","LevelstoNormFC","adjmethod",
                        "replaceOutliers","useaFold","quiet","lmodel","preval",
                        "qforkappa","scale")
                    not.valid <- which(!valid)
                },
                dss = {
                    valid <- names(argList) %in% c("trend","equal.var")
                    not.valid <- which(!valid)
                }
            )
            if (length(not.valid)>0) {
                warnwrap(paste("The following",method,what,"argument names",
                    "are invalid and will be ignored:",
                    paste(names(argList)[not.valid],collapse=", ")))
                argList[not.valid] <- NULL
            }
            return(argList)
        },
        geneFilter = {
            valid.1 <- names(argList) %in% c("length","avgReads","expression",
                "biotype","presence")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warnwrap(paste("The following",method,what,"argument names",
                    "are invalid and will be ignored:",
                    paste(names(argList)[not.valid.1],collapse=", ")))
                argList[not.valid.1] <- NULL
            }
            if (length(argList)>0) {
                for (n in names(argList)) {
                    switch(n,
                        length = {
                            valid.2 <- names(argList[[n]]) %in% c("length")
                            not.valid.2 <- which(!valid.2)
                        },
                        avgReads = {
                            valid.2 <- names(argList[[n]]) %in%
                                c("average.per.bp","quantile")
                            not.valid.2 <- which(!valid.2)
                        },
                        expression = {
                            valid.2 <- names(argList[[n]]) %in% c("median",
                                "mean","quantile","known","custom")
                            not.valid.2 <- which(!valid.2)
                        },
                        presence = {
                            valid.2 <- names(argList[[n]]) %in% c("frac",
                                "min.count","per.condition")
                            not.valid.2 <- which(!valid.2)
                        }
                    )
                    if (length(not.valid.2)>0) {
                        warnwrap(paste("The following",method,what,
                            "sub-argument names are invalid and will be",
                            "ignored:",paste(names(argList[[n]])[not.valid.2],
                            collapse=", ")))
                        argList[[n]][not.valid.2] <- NULL
                    }
                }
            }
            return(argList)
        },
        exonFilter = {
            valid.1 <- names(argList) %in% c("mnrpx")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warnwrap(paste("The following",method,what,"argument names",
                    "are invalid and will be ignored:",
                    paste(names(argList)[not.valid.1],collapse=", ")))
                argList[not.valid.1] <- NULL
            }
            if (length(argList)>0) {
                for (n in names(argList)) {
                    switch(n,
                        mnrpx = {
                            valid.2 <- names(argList[[n]]) %in%
                                c("exons.per.gene","min.exons","frac")
                            not.valid.2 <- which(!valid.2)
                        }
                    )
                    if (length(not.valid.2)>0) {
                        warnwrap(paste("The following",method,what,
                            "sub-argument names are invalid and will be",
                            "ignored:",paste(names(argList[[n]])[not.valid.2],
                            collapse=", ")))
                        argList[[n]][not.valid.2] <- NULL
                    }
                }
            }
            return(argList)
        },
        biotypeFilter = {
            switch(method,
                hg18 = {
                    valid <- names(argList) %in% c("unprocessed_pseudogene",
                        "pseudogene","miRNA","retrotransposed","protein_coding",
                        "processed_pseudogene","snRNA","snRNA_pseudogene",
                        "Mt_tRNA_pseudogene","miRNA_pseudogene","misc_RNA",
                        "tRNA_pseudogene","snoRNA","scRNA_pseudogene",
                        "rRNA_pseudogene","snoRNA_pseudogene","rRNA",
                        "misc_RNA_pseudogene","IG_V_gene","IG_D_gene",
                        "IG_J_gene","IG_C_gene","IG_pseudogene","scRNA")
                    not.valid <- which(!valid)
                },
                hg19 = {
                    valid <- names(argList) %in% c("pseudogene","lincRNA",
                        "protein_coding","antisense","processed_transcript",
                        "snRNA","sense_intronic","miRNA","misc_RNA","snoRNA",
                        "rRNA","polymorphic_pseudogene","sense_overlapping",
                        "three_prime_overlapping_ncrna","TR_V_gene",
                        "TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene",
                        "TR_J_pseudogene","IG_C_gene","IG_C_pseudogene",
                        "IG_J_gene","IG_J_pseudogene","IG_D_gene","IG_V_gene",
                        "IG_V_pseudogene")
                    not.valid <- which(!valid)
                },
                mm9 = {
                    valid <- names(argList) %in% c("pseudogene","snRNA",
                        "protein_coding","antisense","miRNA","lincRNA","snoRNA",
                        "processed_transcript","misc_RNA","rRNA",
                        "sense_overlapping","sense_intronic",
                        "polymorphic_pseudogene","non_coding",
                        "three_prime_overlapping_ncrna","IG_C_gene","IG_J_gene",
                        "IG_D_gene","IG_V_gene","ncrna_host")
                    not.valid <- which(!valid)
                },
                mm10 = {
                    valid <- names(argList) %in% c("pseudogene","snRNA",
                        "protein_coding","antisense","miRNA","snoRNA","lincRNA",
                        "processed_transcript","misc_RNA","rRNA",
                        "sense_intronic","sense_overlapping",
                        "polymorphic_pseudogene","IG_C_gene","IG_J_gene",
                        "IG_D_gene","IG_LV_gene","IG_V_gene","IG_V_pseudogene",
                        "TR_V_gene","TR_V_pseudogene",
                        "three_prime_overlapping_ncrna")
                    not.valid <- which(!valid)
                },
                dm3 = {
                    valid <- names(argList) %in% c("protein_coding","ncRNA",
                        "snoRNA","pre_miRNA","pseudogene","snRNA","tRNA","rRNA")
                    not.valid <- which(!valid)
                },
                rn5 = {
                    valid <- names(argList) %in% c("protein_coding",
                        "pseudogene","processed_pseudogene","miRNA","rRNA",
                        "misc_RNA")
                    not.valid <- which(!valid)
                },
                rn6 = {
                    valid <- names(argList) %in% c("antisense","lincRNA",
                        "miRNA","misc_RNA","Mt_rRNA","Mt_tRNA",
                        "processed_pseudogene","processed_transcript",
                        "protein_coding","pseudogene","ribozyme","rRNA",
                        "scaRNA","sense_intronic","snoRNA","snRNA","sRNA",
                        "TEC","transcribed_processed_pseudogene",
                        "transcribed_unprocessed_pseudogene",
                        "unprocessed_pseudogene")
                    not.valid <- which(!valid)
                },
                danrer7 = {
                    valid <- names(argList) %in% c("antisense",
                        "protein_coding","miRNA","snoRNA","rRNA","lincRNA",
                        "processed_transcript","snRNA","pseudogene",
                        "sense_intronic","misc_RNA","polymorphic_pseudogene",
                        "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                        "non_coding","sense_overlapping")
                    not.valid <- which(!valid)
                },
                pantro4 = {
                    valid <- names(argList) %in% c("protein_coding",
                        "pseudogene","processed_pseudogene","miRNA","rRNA",
                        "snRNA","snoRNA","misc_RNA")
                    not.valid <- which(!valid)
                },
                susscr3 = {
                    valid <- names(argList) %in% c("antisense",
                        "protein_coding","lincRNA","pseudogene",
                        "processed_transcript","miRNA","rRNA","snRNA","snoRNA",
                        "misc_RNA","non_coding","IG_C_gene","IG_J_gene",
                        "IG_V_gene","IG_V_pseudogene")
                    not.valid <- which(!valid)
                },
                tair10 = {
                    valid <- names(argList) %in% c("miRNA","ncRNA",
                        "protein_coding","pseudogene","rRNA","snoRNA",
                        "snRNA","transposable_element","tRNA")
                    not.valid <- which(!valid)
                },
                equcab2 = {
                    valid <- names(argList) %in% c("miRNA","misc_RNA",
                        "protein_coding","pseudogene","processed_pseudogene",
                        "rRNA","snoRNA","snRNA")
                }
            )
            if (length(not.valid)>0) {
                warnwrap(paste("The following",method,what,"argument names",
                    "are invalid and will be ignored:",
                    paste(names(argList)[not.valid],collapse=", ")))
                argList[not.valid] <- NULL
            }
            return(argList)
        }
    )
}

getStrictBiofilter <- function(org) {
    switch(org,
        hg18 = {
            return(list(
                unprocessed_pseudogene=TRUE,
                pseudogene=TRUE,
                miRNA=FALSE,
                retrotransposed=FALSE,
                protein_coding=FALSE,
                processed_pseudogene=TRUE,
                snRNA=FALSE,
                snRNA_pseudogene=TRUE,
                Mt_tRNA_pseudogene=TRUE,
                miRNA_pseudogene=TRUE,
                misc_RNA=TRUE,
                tRNA_pseudogene=TRUE,
                snoRNA=TRUE,
                scRNA_pseudogene=TRUE,
                rRNA_pseudogene=TRUE,
                snoRNA_pseudogene=TRUE,
                rRNA=TRUE,
                misc_RNA_pseudogene=TRUE,
                IG_V_gene=FALSE,
                IG_D_gene=FALSE,
                IG_J_gene=FALSE,
                IG_C_gene=FALSE,
                IG_pseudogene=TRUE,
                scRNA=FALSE
            ))
        },
        hg19 = {
            return(list(
                pseudogene=TRUE,
                lincRNA=FALSE,
                protein_coding=FALSE,
                antisense=FALSE,
                processed_transcript=FALSE,
                snRNA=FALSE,
                sense_intronic=FALSE,
                miRNA=FALSE,
                misc_RNA=FALSE,
                snoRNA=TRUE,
                rRNA=TRUE,
                polymorphic_pseudogene=TRUE,
                sense_overlapping=FALSE,
                three_prime_overlapping_ncrna=FALSE,
                TR_V_gene=FALSE,
                TR_V_pseudogene=TRUE,
                TR_D_gene=FALSE,
                TR_J_gene=FALSE,
                TR_C_gene=FALSE,
                TR_J_pseudogene=TRUE,
                IG_C_gene=FALSE,
                IG_C_pseudogene=TRUE,
                IG_J_gene=FALSE,
                IG_J_pseudogene=TRUE,
                IG_D_gene=FALSE,
                IG_V_gene=FALSE,
                IG_V_pseudogene=TRUE
            ))
        },
        hg38 = {
            return(list(
                protein_coding=FALSE,
                polymorphic_pseudogene=TRUE,
                lincRNA=FALSE,
                unprocessed_pseudogene=TRUE,
                processed_pseudogene=TRUE,
                antisense=FALSE,
                processed_transcript=FALSE,
                transcribed_unprocessed_pseudogene=TRUE,
                sense_intronic=FALSE,
                unitary_pseudogene=TRUE,
                IG_V_gene=FALSE,
                IG_V_pseudogene=TRUE,
                TR_V_gene=FALSE,
                sense_overlapping=FALSE,
                transcribed_processed_pseudogene=TRUE,
                miRNA=FALSE,
                snRNA=FALSE,
                misc_RNA=FALSE,
                rRNA=TRUE,
                snoRNA=TRUE,
                IG_J_pseudogene=TRUE,
                IG_J_gene=FALSE,
                IG_D_gene=FALSE,
                three_prime_overlapping_ncrna=FALSE,
                IG_C_gene=FALSE,
                IG_C_pseudogene=TRUE,
                pseudogene=TRUE,
                TR_V_pseudogene=TRUE,
                Mt_tRNA=TRUE,
                Mt_rRNA=TRUE,
                translated_processed_pseudogene=TRUE,
                TR_J_gene=FALSE,
                TR_C_gene=FALSE,
                TR_D_gene=FALSE,
                TR_J_pseudogene=TRUE,
                LRG_gene=FALSE
            ))
        },
        mm9 = {
            return(list(
                pseudogene=TRUE,
                snRNA=FALSE,
                protein_coding=FALSE,
                antisense=FALSE,
                miRNA=FALSE,
                lincRNA=FALSE,
                snoRNA=TRUE,
                processed_transcript=FALSE,
                misc_RNA=TRUE,
                rRNA=TRUE,
                sense_overlapping=FALSE,
                sense_intronic=FALSE,
                polymorphic_pseudogene=TRUE,
                non_coding=FALSE,
                three_prime_overlapping_ncrna=FALSE,
                IG_C_gene=FALSE,
                IG_J_gene=FALSE,
                IG_D_gene=FALSE,
                IG_V_gene=FALSE,
                ncrna_host=FALSE
            ))
        },
        mm10 = {
            return(list(
                pseudogene=TRUE,
                snRNA=FALSE,
                protein_coding=FALSE,
                antisense=FALSE,
                miRNA=FALSE,
                snoRNA=TRUE,
                lincRNA=FALSE,
                processed_transcript=FALSE,
                misc_RNA=TRUE,
                rRNA=TRUE,
                sense_intronic=FALSE,
                sense_overlapping=FALSE,
                polymorphic_pseudogene=TRUE,
                IG_C_gene=FALSE,
                IG_J_gene=FALSE,
                IG_D_gene=FALSE,
                IG_LV_gene=FALSE,
                IG_V_gene=FALSE,
                IG_V_pseudogene=TRUE,
                TR_V_gene=FALSE,
                TR_V_pseudogene=TRUE,
                three_prime_overlapping_ncrna=FALSE
            ))
        },
        dm3 = {
            return(list(
                protein_coding=FALSE,
                ncRNA=FALSE,
                snoRNA=TRUE,
                pre_miRNA=FALSE,
                pseudogene=TRUE,
                snRNA=FALSE,
                tRNA=FALSE,
                rRNA=TRUE
            ))
        },
        rn5 = {
            return(list(
                protein_coding=FALSE,
                pseudogene=TRUE,
                processed_pseudogene=FALSE,
                miRNA=FALSE,
                rRNA=TRUE,
                misc_RNA=TRUE
            ))
        },
        rn6 = {
            return(list(
                antisense=FALSE, 
                lincRNA=FALSE, 
                miRNA=FALSE, 
                misc_RNA=TRUE, 
                Mt_rRNA=TRUE, 
                Mt_tRNA=TRUE, 
                processed_pseudogene=FALSE, 
                processed_transcript=FALSE, 
                protein_coding=FALSE, 
                pseudogene=FALSE, 
                ribozyme=FALSE, 
                rRNA=TRUE, 
                scaRNA=FALSE, 
                sense_intronic=FALSE, 
                snoRNA=TRUE, 
                snRNA=FALSE, 
                sRNA=FALSE, 
                TEC=FALSE, 
                transcribed_processed_pseudogene=FALSE, 
                transcribed_unprocessed_pseudogene=FALSE, 
                unprocessed_pseudogene=TRUE
            ))
        },
        danrer7 = {
            return(list(
                antisense=FALSE,
                protein_coding=FALSE,
                miRNA=FALSE,
                snoRNA=TRUE,
                rRNA=TRUE,
                lincRNA=FALSE,
                processed_transcript=FALSE,
                snRNA=FALSE,
                pseudogene=TRUE,
                sense_intronic=FALSE,
                misc_RNA=TRUE,
                polymorphic_pseudogene=TRUE,
                IG_V_pseudogene=TRUE,
                IG_C_pseudogene=TRUE,
                IG_J_pseudogene=TRUE,
                non_coding=FALSE,
                sense_overlapping=FALSE
            ))
        },
        pantro4 = {
            return(list(
                protein_coding=FALSE,
                pseudogene=TRUE,
                processed_pseudogene=TRUE,
                miRNA=FALSE,
                rRNA=TRUE,
                snRNA=TRUE,
                snoRNA=TRUE,
                misc_RNA=TRUE
            ))
        },
        susscr3 = {
            return(list(
                antisense=FALSE,
                protein_coding=FALSE,
                lincRNA=FALSE,
                pseudogene=TRUE,
                processed_transcript=FALSE,
                miRNA=FALSE,
                rRNA=TRUE,
                snRNA=TRUE,
                snoRNA=TRUE,
                misc_RNA=TRUE,
                non_coding=FALSE,
                IG_C_gene=TRUE,
                IG_J_gene=TRUE,
                IG_V_gene=TRUE,
                IG_V_pseudogene=FALSE
            ))
        },
        tair10 = {
            return(list(
                miRNA=FALSE,
                ncRNA=FALSE,
                protein_coding=FALSE,
                pseudogene=TRUE,
                rRNA=TRUE,
                snoRNA=TRUE,
                snRNA=TRUE,
                transposable_element=FALSE,
                tRNA=TRUE
            ))
        },
        equcab2 = {
            return(list(
                miRNA=FALSE,
                misc_RNA=TRUE,
                protein_coding=FALSE,
                pseudogene=FALSE,
                processed_pseudogene=FALSE,
                rRNA=TRUE,
                snoRNA=TRUE,
                snRNA=TRUE
            ))
        }
    )
}

getPresetOpts <- function(preset,org) {
    # Override filter rules and maybe normArgs and statArgs
    switch(preset,
        all_basic = {
            exonFilters <- NULL
            geneFilters <- NULL
            pcut <- NA
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean")
        },
        all_normal = {
            exonFilters <- NULL
            geneFilters <- NULL
            pcut <- NA
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean","sd","cv")
        },
        all_full = {
            exonFilters <- NULL
            geneFilters <- NULL
            pcut <- NA
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts","flags")
            exportScale <- c("natural","log2","log10","vst")
            exportValues <- c("raw","normalized")
            exportStats <- c("mean","median","sd","mad","cv","rcv")
        },
        medium_basic = {
            exonFilters <- list(
                minActiveExons=list(
                    exonsPerGene=5,
                    minExons=2,
                    frac=1/5
                )
            )
            geneFilters <- list(
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
                biotype=getDefaults("biotypeFilter",org[1])
            )
            pcut <- 0.05
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean")
        },
        medium_normal = {
            exonFilters <- list(
                minActiveExons=list(
                    exonsPerGene=5,
                    minExons=2,
                    frac=1/5
                )
            )
            geneFilters <- list(
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
                biotype=getDefaults("biotypeFilter",org[1])
            )
            pcut <- 0.05
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean","sd","cv")
        },
        medium_full = {
            exonFilters <- list(
                minActiveExons=list(
                    exonsPerGene=5,
                    minExons=2,
                    frac=1/5
                )
            )
            geneFilters <- list(
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
                biotype=getDefaults("biotypeFilter",org[1])
            )
            pcut <- 0.05
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts","flags")
            exportScale <- c("natural","log2","log10","vst")
            exportValues <- c("raw","normalized")
            exportStats <- c("mean","median","sd","mad","cv","rcv")
        },
        strict_basic = {
            exonFilters=list(
                minActiveExons=list(
                    exons_per_gene=4,
                    min_exons=2,
                    frac=1/4
                )
            )
            geneFilters=list(
                length=list(
                    length=750
                ),
                avgReads=list(
                    averagePerBp=100,
                    quantile=0.5
                ),
                expression=list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA,
                    custom=NA
                ),
                biotype=getStrictBiofilter(org[1])
            )
            pcut <- 0.01
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean")
        },
        strict_normal = {
            exonFilters=list(
                minActiveExons=list(
                    exonsPerGene=4,
                    minExons=2,
                    frac=1/4
                )
            )
            geneFilters=list(
                length=list(
                    length=750
                ),
                avgReads=list(
                    averagePerBp=100,
                    quantile=0.5
                ),
                expression=list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA,
                    custom=NA
                ),
                biotype=getStrictBiofilter(org[1])
            )
            pcut <- 0.01
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts")
            exportScale <- c("natural","log2")
            exportValues <- c("normalized")
            exportStats <- c("mean","sd","cv")
        },
        strict_full = {
            exonFilters=list(
                minActiveExons=list(
                    exonsPerGene=4,
                    minExons=2,
                    frac=1/4
                )
            )
            geneFilters=list(
                length=list(
                    length=750
                ),
                avgReads=list(
                    averagePerBp=100,
                    quantile=0.5
                ),
                expression=list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA,
                    custom=NA
                ),
                biotype=getStrictBiofilter(org[1])
            )
            pcut <- 0.01
            exportWhat <- c("annotation","p_value","adj_p_value",
                "meta_p_value","adj_meta_p_value","fold_change","stats",
                "counts","flags")
            exportScale <- c("natural","log2","log10","vst")
            exportValues <- c("raw","normalized")
            exportStats <- c("mean","median","sd","mad","cv","rcv")
        }
    )
    presetOpts <- list(
        exonFilters=exonFilters,
        geneFilters=geneFilters,
        pcut=pcut,
        exportWhat=exportWhat,
        exportScale=exportScale,
        exportValues=exportValues,
        exportStats=exportStats
    )
    return(presetOpts)
}

makeFoldChange <- function(contrast,sampleList,dataMatrix,logOffset=1) {
    conds <- strsplit(contrast,"_vs_")[[1]]
    foldMat <- matrix(0,nrow(dataMatrix),length(conds)-1)
    for (i in seq_len(length(conds)-1)) { # Last condition is ALWAYS reference
    #for (i in 2:length(conds)) { # First condition is ALWAYS reference
        samplesNom <- sampleList[[conds[i]]]
        samplesDenom <- sampleList[[conds[length(conds)]]]
        #samplesDenom <- sampleList[[conds[1]]]
        nom <- dataMatrix[,match(samplesNom,colnames(dataMatrix)),drop=FALSE]
        denom <- dataMatrix[,match(samplesDenom,colnames(dataMatrix)),
            drop=FALSE]
        if (!is.matrix(nom)) 
            nom <- as.matrix(nom) # Cover the case with no replicates...
        if (!is.matrix(denom)) 
            denom <- as.matrix(denom)
        meanNom <- apply(nom,1,mean)
        meanDenom <- apply(denom,1,mean)
        #meanNom <- ifelse(meanNom==0,logOffset,meanNom)
        if (any(meanNom==0)) 
            meanNom <- meanNom + logOffset
        #meanDenom <- ifelse(meanDenom==0,logOffset,meanDenom)
        if (any(meanDenom==0)) 
            meanDenom <- meanDenom + logOffset
        #foldMat[,i-1] <- meanNom/meanDenom
        foldMat[,i] <- meanNom/meanDenom
    }
    rownames(foldMat) <- rownames(dataMatrix)
    #colnames(foldMat) <- paste(conds[1],"_vs_",conds[2:length(conds)],sep="")
    colnames(foldMat) <- paste(conds[seq_len(length(conds)-1)],"_vs_",
        conds[length(conds)],sep="")
    return(foldMat)
}

makeA <- function(contrast,sampleList,dataMatrix,logOffset=1) {
    conds <- strsplit(contrast,"_vs_")[[1]]
    aMat <- matrix(0,nrow(dataMatrix),length(conds)-1)
    for (i in seq_len(length(conds)-1)) { # Last condition is ALWAYS reference
    #for (i in 2:length(conds)) { # First condition is ALWAYS reference
        samplesTrt <- sampleList[[conds[i]]]
        samplesCnt <- sampleList[[conds[length(conds)]]]
        #samplesCnt <- sampleList[[conds[1]]]
        trt <- dataMatrix[,match(samplesTrt,colnames(dataMatrix)),drop=FALSE]
        cnt <- dataMatrix[,match(samplesCnt,colnames(dataMatrix)),drop=FALSE]
        meanTrt <- apply(trt,1,mean)
        meanCnt <- apply(cnt,1,mean)
        if (any(meanTrt==0)) 
            meanTrt <- meanTrt + logOffset
        if (any(meanCnt==0)) 
            meanCnt <- meanCnt + logOffset
        aMat[,i] <- 0.5*(log2(meanTrt)+log2(meanCnt))
        #aMat[,i-1] <- 0.5*(log2(meanTrt)+log2(meanCnt))
    }
    rownames(aMat) <- rownames(dataMatrix)
    colnames(aMat) <- paste(conds[seq_len(length(conds)-1)],"_vs_",
        conds[length(conds)],sep="")
    #colnames(aMat) <- paste(conds[2:length(conds)],"_vs_",conds[1],sep="")
    return(aMat)
}

makeAvgExpression <- function(contrast,sampleList,dataMatrix,logOffset=1) {
    conds <- strsplit(contrast,"_vs_")[[1]]
    aMat <- matrix(0,nrow(dataMatrix),length(conds)-1)
    for (i in seq_len(length(conds)-1)) { # Last condition is ALWAYS reference
    #for (i in 2:length(conds)) { # First condition is ALWAYS reference
        samplesNom <- sampleList[[conds[i]]]
        #samplesDenom <- sampleList[[conds[1]]]
        samplesDenom <- sampleList[[conds[length(conds)]]]
        nom <- dataMatrix[,match(samplesNom,colnames(dataMatrix)),drop=FALSE]
        denom <- dataMatrix[,match(samplesDenom,colnames(dataMatrix)),
            drop=FALSE]
        if (!is.matrix(nom)) # Cover the case with no replicates...
            nom <- as.matrix(nom) 
        if (!is.matrix(denom)) 
            denom <- as.matrix(denom)
        meanNom <- apply(nom,1,mean)
        meanDenom <- apply(denom,1,mean)
        if (any(meanNom==0)) 
            meanNom <- meanNom + logOffset
        if (any(meanDenom==0)) 
            meanDenom <- meanDenom + logOffset
        aMat[,i] <- 0.5*(log2(meanNom)+log2(meanDenom))
        #aMat[,i-1] <- 0.5*(log2(meanNom)+log2(meanDenom))
    }
    rownames(aMat) <- rownames(dataMatrix)
    #colnames(aMat) <- paste(conds[1],"_vs_",conds[2:length(conds)],sep="")
    colnames(aMat) <- paste(conds[seq_len(length(conds)-1)],"_vs_",
        conds[length(conds)],sep="")
    return(aMat)
}

.makeHtmlCells <- function(mat,type="numeric",digits=3) {
    if (type=="numeric")
        tmp <- format(mat,digits=digits)
        #tmp <- formatC(mat,digits=digits,format="f")
    else
        tmp <- mat
    if (!is.matrix(tmp)) tmp <- as.matrix(tmp)
    tmp <- apply(tmp,c(1,2),function(x) paste("<td>",x,"</td>",sep=""))
    return(tmp)
}

.makeHtmlRows <- function(mat) {
    tmp <- apply(mat,1,paste,collapse="")
    tmp <- paste("<tr>",tmp,"</tr>",sep="")
    return(tmp)
}

.makeHtmlHeader <- function(h) {
    tmp <- paste("<th>",h,"</th>",sep="")
    tmp <- paste(tmp,collapse="")
    tmp <- paste("<tr>",tmp,"</tr>",sep="")
    return(tmp)
}

.makeHtmlBody <- function(mat) {
    tmp <- paste(mat,collapse="")
    return(tmp)
}

.makeHtmlTable <- function(b,h=NULL,id=NULL) {
    if (!is.null(id))
        html <- paste("<table id=\"",id,"\" class=\"datatable\">",sep="")
    else
        html <- "<table class=\"datatable\">"
    if (!is.null(h))
        html <- paste(html,"<thead>",h,"</thead>",sep="")
    html <- paste(html,"<tbody>",b,"</tbody></table>",sep="")
    return(html)
}    

makeTransformation <- function(dataMatrix,exportScale,
    scf=NULL,logOffset=1) {
    mat <- vector("list",length(exportScale))
    names(mat) <- exportScale
    if (!is.matrix(dataMatrix)) dataMatrix <- as.matrix(dataMatrix)
    if (is.null(scf)) scf <- rep(1,nrow(dataMatrix))
    for (scl in exportScale) {
        switch(scl,
            natural = {
                mat[[scl]] <- dataMatrix
            },
            log2 = {
                mat[[scl]] <- nat2log(dataMatrix,base=2,off=logOffset)
            },
            log10 = {
                mat[[scl]] <- nat2log(dataMatrix,base=10,off=logOffset)
            },
            vst = {
                fit <- vsn2(dataMatrix,verbose=FALSE)
                mat[[scl]] <- predict(fit,newdata=dataMatrix)
            },
            rpgm = {
                mat[[scl]] <- dataMatrix
                for (i in seq_len(ncol(dataMatrix)))
                    mat[[scl]][,i] <- dataMatrix[,i]/scf
            }
        )
    }
    return(mat)
}

makeStat <- function(samples,dataList,stat,exportScale) {
    statResult <- vector("list",length(exportScale))
    names(statResult) <- exportScale
    for (scl in exportScale) {
        statData <- dataList[[scl]][,match(samples,colnames(dataList[[scl]]))]
        if (!is.matrix(statData))
            statData <- as.matrix(statData)
        switch(stat,
            mean = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(round(mean(x)))
                    else
                        return(mean(x))
                },scl)
            },
            median = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(round(median(x)))
                    else
                        return(median(x))
                },scl)
            },
            sd = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(ceiling(sd(x)))
                    else
                        return(sd(x))
                },scl)
            },
            mad = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(ceiling(mad(x)))
                    else
                        return(mad(x))
                },scl)
            },
            cv = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(ceiling(sd(x))/round(mean(x)))
                    else
                        return(sd(x)/mean(x))
                },scl)
            },
            rcv = {
                statResult[[scl]] <- apply(statData,1,function(x,s) {
                    if (s=="natural")
                        return(ceiling(mad(x))/round(median(x))) 
                    else
                        return(mad(x)/median(x))
                },scl)
            }
        )
    }
    return(do.call("cbind",statResult))
}

makeMatrix <- function(samples,dataList,exportScale="natural") {
    mat <- vector("list",length(exportScale))
    names(mat) <- exportScale
    for (scl in exportScale) {
        mat.data <- dataList[[scl]][,match(samples,
            colnames(dataList[[scl]]))]
        if (!is.matrix(mat.data)) {
            mat.data <- as.matrix(mat.data)
            colnames(mat.data) <- samples
        }
        mat[[scl]] <- mat.data
    }
    return(do.call("cbind",mat))
}

makeContrastList <- function(contrast,sampleList) {
    # Construction
    contrastList <- vector("list",length(contrast))
    names(contrastList) <- contrast
    # First break the contrast vector
    cnts <- strsplit(contrast,"_vs_")
    names(cnts) <- names(contrastList)
    # Create list members
    for (n in names(contrastList)) {
        contrastList[[n]] <- vector("list",length(cnts[[n]]))
        for (i in seq_len(length(cnts[[n]]))) {
            contrastList[[n]][[i]] <- rep(cnts[[n]][i],
                length(sampleList[[cnts[[n]][i]]]))
            names(contrastList[[n]][[i]]) <- sampleList[[cnts[[n]][[i]]]]
        }
    }
    return(contrastList)
}

makeSampleList <- function(input,type=c("simple","targets")) {
    if (missing(input) || !file.exists(input))
        stopwrap("File to make sample list from should be a valid existing ",
            "text file!")
    type <- tolower(type[1])
    checkTextArgs("type",type,c("simple","targets"),multiarg=FALSE)
    tab <- read.delim(input)
    samples <- as.character(tab[,1])
    if (type=="simple") {
        ii <- 2
        conditions <- unique(as.character(tab[,ii]))
    }
    else if (type=="targets") {
        ii <- 3
        conditions <- unique(as.character(tab[,ii]))
    }
    if (length(samples) != length(unique(samples)))
        stopwrap("Sample names must be unique for each sample!")
    sampleList <- vector("list",length(conditions))
    names(sampleList) <- conditions
    for (n in conditions)
        sampleList[[n]] <- samples[which(as.character(tab[,ii])==n)]
    return(sampleList)
}

makeProjectPath <- function(path,f=NULL) {
    if (is.na(path) || is.null(path)) {
        if (!is.data.frame(f) && !is.null(f) && !is.list(f) && file.exists(f)) {
            if (length(grep(".RData$",f))>0)
                mainPath <- file.path(getwd(),paste("metaseqr_result_",
                    format(Sys.time(),format="%Y%m%d%H%M%S"),sep=""))
            else
                mainPath <- file.path(dirname(f),paste("metaseqr_result_",
                    format(Sys.time(),format="%Y%m%d%H%M%S"),sep=""))
        }
        else
            mainPath <- file.path(getwd(),paste("metaseqr_result_",
                format(Sys.time(),format="%Y%m%d%H%M%S"),sep=""))
        projectPath <- makePathStruct(mainPath)
    }
    else {
        success <- tryCatch(
            if (!file.exists(path)) dir.create(path,recursive=TRUE) else TRUE,
            error=function(e) {
                disp("Cannot create ",path,"! Is it a valid system path? Is ",
                    "there a write permissions problem? Reverting to ",
                    "automatic creation...")
                return(FALSE)
            },
            finally=""
        )
        if (success)
            projectPath <- makePathStruct(path)
        else
            projectPath <- makeProjectPath(NA,f)
    }
    return(projectPath)
}

makePathStruct <- function(mainPath) {
    projectPath <- list(
        main=mainPath,
        js=file.path(mainPath,"js"),
        media=file.path(mainPath,"media"),
        data=file.path(mainPath,"data"),
        logs=file.path(mainPath,"logs"),
        lists=file.path(mainPath,"lists"),
        tracks=file.path(mainPath,"tracks"),
        plots=file.path(mainPath,"plots"),
        qc=file.path(mainPath,"plots","qc"),
        normalization=file.path(mainPath,"plots","normalization"),
        statistics=file.path(mainPath,"plots","statistics")
    )
    for (p in names(projectPath))
        if (!file.exists(projectPath[[p]]))
            dir.create(projectPath[[p]],recursive=TRUE)
    return(projectPath)
}

makeExportList <- function(con) {
    f <- vector("list",length(con))
    names(f) <- con
    return(f)
}

makeGrid <- function(n) {
    m <- 0
    while (n > m*m)
        m <- m+1
    if (n < m*m) {
        k <- m-1
        if (n > m*k)
            k <- k+1
        else {
            while (n > m*k)
                k=k-1
        }
    }
    else
        k <- m
    return(c(m,k))
}

makeReportMessages <- function(lang) {
    switch(lang,
        en = {
            messages <- list(
                org=list(
                    hg18=paste("human (<em>Homo sapiens</em>),",
                        "genome version alias hg18"),
                    hg19=paste("human (<em>Homo sapiens</em>),",
                        "genome version alias hg19"),
                    hg38=paste("human (<em>Homo sapiens</em>),",
                        "genome version alias hg38"),
                    mm9=paste("mouse (<em>Mus musculus</em>),",
                        "genome version alias mm9"),
                    mm10=paste("mouse (<em>Mus musculus</em>),",
                        "genome version alias mm10"),
                    rn5=paste("rat (<em>Rattus norvegicus</em>),",
                        "genome version  alias rn5"),
                    rn6=paste("rat (<em>Rattus norvegicus</em>),",
                        "genome version alias rn6"),
                    dm3=paste("fruitfly (<em>Drosophila melanogaster</em>),",
                        "genome version alias dm3"),
                    dm6=paste("fruitfly (<em>Drosophila melanogaster</em>),",
                        "genome version alias dm6"),
                    danrer7=paste("zebrafish (<em>Danio rerio</em>),",
                        "genome version alias danRer7"),
                    danrer10=paste("zebrafish (<em>Danio rerio</em>),",
                        "genome version alias danRer10"),
                    danrer11=paste("zebrafish (<em>Danio rerio</em>),",
                        "genome version alias danRer11"),
                    pantro4=paste("chimpanzee (<em>Pan troglodytes</em>),",
                        "genome version alias panTro4"),
                    pantro5=paste("chimpanzee (<em>Pan troglodytes</em>),",
                        "genome version alias panTro5"),
                    susscr3=paste("pig (<em>Sus scrofa</em>),",
                        "genome version alias susScr3"),
                    susscr11=paste("pig (<em>Sus scrofa</em>),",
                        "genome version alias susScr11"),
                    equcab2=paste("horse (<em>Equus cabalus</em>),",
                        "genome version alias equcab2"),
                    tair10=paste("arabidopsis (<em>Arabidobsis thaliana</em>)",
                        ",","genome version alias TAIR10")
                ),
                refdb=list(
                    ensembl="Ensembl genomes",
                    ucsc="UCSC genomes database",
                    refseq="RefSeq database",
                    other="User provided"
                ),
                whenfilter=list(
                    prenorm="before normalization",
                    postnorm="after normalization"
                ),
                norm=list(
                    edaseq="EDASeq",
                    deseq="DESeq",
                    deseq2="DESeq2",
                    edger="edgeR",
                    noiseq="NOISeq",
                    nbpseq="NBPSeq",
                    absseq="ABSSeq",
                    dss="DSS",
                    each="the same as the corresponding statistical test"
                ),
                stat=list(
                    deseq="DESeq",
                    deseq2="DESeq2",
                    edger="edgeR",
                    noiseq="NOISeq",
                    bayseq="baySeq",
                    limma="limma",
                    nbpseq="NBPSeq",
                    absseq="ABSSeq",
                    dss="DSS"
                ),
                meta=list(
                    intersection="intersection of individual results",
                    union="union of individual results",
                    fisher="Fisher's method",
                    fperm="Fisher's method with permutations",
                    dperm_min=paste("samples permutation based method with",
                        "minimum p-values"),
                    dperm_max=paste("samples permutation based method with",
                        "maximum p-values"),
                    dperm_weight=paste("samples permutation based method with",
                        "weighted p-values"),
                    minp="minimum p-value across results",
                    maxp="maximum p-value across results",
                    weight="PANDORA weighted p-value across results",
                    pandora="PANDORA weighted p-value across results",
                    simes="Simes correction and combination method",
                    whitlock=paste("Whitlock's Z-transformation method",
                        "(Bioconductor package survcomp)"),
                    none=paste("no meta-analysis, reported p-values from the",
                        "first supplied statistical algorithm")
                ),
                adjust=list(
                    holm="Holm FWER",
                    hochberg="Hochberg DFR",
                    hommel="Hommel FWER",
                    bonferroni="Bonferroni FWER",
                    bh="Benjamini-Hochberg FDR",
                    by="Benjamini-Yekutiely FDR",
                    fdr="Benjamini-Hochberg FDR",
                    none="no multiple test correction",
                    qvalue="Storey-Tibshirani FDR"
                ),
                plots=list(
                    mds="multidimensional scaling",
                    biodetection="biotype detection",
                    countsbio="biotype counts",
                    saturation="sample and biotype saturation",
                    rnacomp="RNA composition",
                    boxplot="boxplots",
                    gcbias="GC-content bias",
                    lengthbias="transcript length bias",
                    meandiff="mean-difference plot",
                    meanvar="mean-variance plot",
                    deheatmap="DEG heatmap",
                    volcano="volcano plot",
                    biodist="DEG biotype detection",
                    filtered="filtered biotypes",
                    correl="correlation heatmap and correlogram",
                    pairwise="pairwise scatterplots between samples",
                    venn="Venn diagrams"
                ),
                export=list(
                    annotation="Annotation",
                    p_value="p-value",
                    meta_p_value="Combined p-value",
                    adj_p_value="Adjusted p-value (FDR)",
                    adj_meta_p_value="Adjusted combined p-value (FDR)",
                    fold_change="Fold change",
                    stats="Statistics",
                    counts="Read counts",
                    natural="Natural scale",
                    log2="log2 scale",
                    log10="log10 scale",
                    vst="Variance stabilization transformation",
                    rpgm="Reads per Gene Model",
                    raw="Raw values",
                    normalized="Normalized values",
                    mean="Mean",
                    median="Median",
                    sd="Standard deviation",
                    mad="Median Absolute Deviation (MAD)",
                    cv="Coefficient of Variation",
                    rcv="Robust Coefficient of Variation"
                ),
                preset=list(
                    all_basic=paste("use all genes and export all genes and",
                        "basic annotation and statistics elements"),
                    all_normal=paste("use all genes and export all genes and",
                        "normal annotation and statistics elements"),
                    all_full=paste("use all genes and export all genes and all",
                        "available annotation and statistics elements"),
                    medium_basic=paste("apply a medium set of filters and",
                        "export statistically significant genes and basic",
                        "annotation and statistics elements"),
                    medium_normal=paste("apply a medium set of filters and",
                        "export statistically significant genes and normal",
                        "annotation and statistics elements"),
                    medium_full=paste("apply a medium set of filters and",
                        "export statistically significant genes and all",
                        "available annotation and statistics elements"),
                    strict_basic=paste("apply a strict set of filters and",
                        "export statistically significant genes and basic",
                        "annotation and statistics elements"),
                    strict_normal=paste("apply a medium set of filters and",
                        "export statistically significant genes and normal",
                        "annotation and statistics elements"),
                    strict_full=paste("apply a medium set of filters and",
                        "export statistically significant genes and all",
                        "available annotation and statistics elements")
                ),
                explain=list(
                    mds=paste(
                "Multidimensional Scaling (MDS) plots constitute a means",
                "of visualizing the level of similarity of individual cases",
                "of a dataset. It is similar to Principal Component Analysis",
                "(PCA), but instead of using the covariance matrix to find",
                "similarities between cases, MDS uses absolute distance",
                "metrics such as the classical Euclidean distance. Because",
                "of the relative linear relations between sequencing samples,",
                "it provides a more realistic clustering of samples. MDS",
                "serves quality control and it can be interpreted as follows:",
                "when the distance between samples of the same biological",
                "condition in the MDS space is small, this is an indication",
                "of high correlation and reproducibility between them. When",
                "this distance is larger or heterogeneous (e.g. the 3rd",
                "sample of a triplicate set is further from the other 2),",
                "this constitutes an indication of low correlation and",
                "reproducibility between samples. It can help exclude poor",
                "samples from further analysis.",collapse=" "
                    ),
                    biodetection=paste(
                "The biotype detection bar diagrams are a set of quality",
                "control charts that show the percentage of each biotype",
                "in the genome (i.e. in the whole set of features provided,",
                "for example, protein coding genes, non coding RNAs or",
                "pseudogenes) in red bars, the proportion of which has been",
                "detected in a sample before normalization and after",
                "basic filtering by removing features with zero counts in",
                "green bars, and the percentage of each biotype within the",
                "sample in blue bars. The difference between red bars and",
                "blue bars is that red bars show the percentage of a",
                "feature in the genome while blue bars show the percentage",
                "in the sample. Thus, the blue bars may sometimes be higher",
                "than the green bars because certain features (e.g. protein",
                "coding genes) may be detected within a sample with a higher",
                "proportion relative to their presence in the genome, as",
                "compared with other features. For example, while the",
                "percentage of protein coding genes in the whole genome is",
                "already higher than other biotypes, this percentage is",
                "expected to be even higher in an RNA-Seq experiment where one",
                "expects protein-coding genes to exhibit greater abundance.",
                "The vertical line separates the most abundant (yellow band)",
                "biotypes (on the left-hand side, corresponding to the",
                "left axis scale) from the rest (on the right-hand side,",
                "corresponding to the right axis scale, red band). Otherwise,",
                "lower abundance biotypes would be indistinguishable.",
                "Unexpected outcomes in this quality control chart (e.g.",
                "very low detection of protein coding genes) would signify",
                "possible low quality of a sample.",collapse=" "
                    ),
                    countsbio=paste(
                "Biotype detection counts boxplots are a set of quality",
                "control charts that depict both the biological classification",
                "for the detected features and the actual distribution of",
                "the read counts for each biological type. The boxplot",
                "comprises a means of summarizing the read counts distribution",
                "of a sample in the form of a bar with extending lines,",
                "as a commonly used way of graphically presenting groups of",
                "numerical data. A boxplot also indicates which observations,",
                "if any, might be considered outliers and is able to visually",
                "show different types of populations, without making any",
                "assumptions about the underlying statistical distribution.",
                "The spacing between the different parts of the box help",
                "indicate variance, skewness and identify outliers. The",
                "thick bar inside the colored box is the median of the",
                "observations while the box extends over the Interquartile",
                "Range of the observations. The whiskers extend up (down)",
                "to +/-1.5xIQR. Unexpected outcomes (e.g. protein coding",
                "read count distribution similar to pseudogene read count",
                "distribution) indicates poor sample quality.",collapse=" "
                    ),
                    saturation=paste(
                "Read and biotype saturation plots are a set of quality",
                "control charts that depict the read count saturation",
                "levels at several sequencing depths. Thus, they comprise",
                "a means of assessing whether the sequencing depth of an",
                "RNA-Seq experiment is sufficient in order to detect the",
                "biological features under investigation. These quality",
                "control charts are separated in two subgroups: the first",
                "(read saturation per biotype for all samples)",
                "is a set of plots, one for each biological feature (e.g.",
                "protein coding, pseudogene, lincRNA, etc.), that depict",
                "the number of detected features in different sequencing",
                "depths and for all samples in the same plot. The second",
                "subgroup (read saturation per sample for all biotypes)",
                "is a set of plots similar to the above, but with",
                "one pair of plots with two panels for each sample,",
                "presenting all biological features. The left panel depicts",
                "the saturation levels for the less abundatnt features,",
                "while the right panel, the saturation for the more abundant",
                "features, as placing them all together would make the",
                "less abundant features indistinguishable. All the saturation",
                "plots should be interpreted as follows: if the read counts",
                "for a biotype tend to be saturated, the respective curve",
                "should tend to reach a plateau at higher depths. Otherwise,",
                "more sequencing is needed for the specific biotype.",
                        collapse=" "
                    ),
                    readnoise=paste(
                "The read noise plots depict the percentage of biological",
                "features detected when subsampling the total number of",
                "reads. Very steep curves in read noise plots indicate",
                "that although the sequencing depth reaches its maximum,",
                "a relatively small percentage of total features is detected,",
                "indicating that the level of background noise is relatively",
                "high. Less steep RNA composition curves, indicate less noise.",
                "When a sample's curve deviate from the rest, it may",
                "indicate lower or higher quality, depending on the curves",
                "of the rest of the samples.",collapse=" "
                    ),
                    correl=paste(
                "Sample correlation plots depict the accordance of",
                "RNA-Seq samples, as this is manifested through the",
                "read counts table used with the metaseqR2 pipeline, with",
                "representations that both use the correlation matrix",
                "(a matrix which depicts all the pairwise correlations",
                "between each pair of samples) of the read counts matrix.",
                "The correlation representation is a clustered heatmap which",
                "depicts the correlations of samples as color-scaled",
                "images and the hierarchical clustering tree depicts the",
                "grouping of the samples according to their correlation.",
                "If samples from the same group not being clustered together",
                "provides an indication that there might be a quality",
                "problem with the dataset.",collapse=" "
                    ),
                    pairwise=paste(
                "Pairwise comparison plots are split in two parts:",
                "the upper part consists of a simple scatterplot for",
                "all pairwise sample comparisons. It is a simple measure",
                "of between sample correlation using all the available",
                "data points instead of only the correlation matrix. The",
                "lower part consists of mean-difference plots for all",
                "pairwise sample comparisons. A mean-difference plot (or",
                "a Bland-Altman plots) is a method of data plotting used",
                "in analyzing the agreement between two different",
                "assays/variables. In this graphical method the differences",
                "(or alternatively the ratios) between the two variables",
                "are plotted against the averages of the two. Such a plot",
                "is useful, for example, for analyzing data with strong",
                "correlation between the x and y axes, when the (x,y) dots on",
                "the plot are close to the diagonal x=y. In this case, the",
                "value of the transformed variable X is approximately the same",
                "as x and y and the variable Y shows the difference between",
                "x and y. In both represantations, irregular shapes of the",
                "red smoother lines are an indication of poor correlation",
                "between samples or of other systematic bias sources,",
                "which is usually corrected through data normalization.",
                collapse=" "
                    ),
                    rnacomp=paste(
                "The RNA composition plots depict differences in the",
                "distributions of reads in the same biological features",
                "across samples. The following is taken from the NOISeq",
                "vignette: <em>'...when two samples have different RNA",
                "composition, the distribution of sequencing reads across",
                "the features is different in such a way that although",
                "a feature had the same number of read counts in both",
                "samples, it would not mean that it was equally expressed",
                "in both... To check if this bias is present in the data,",
                "the RNA composition plot and the correponding diagnostic",
                "test can be used. In this case, each sample s is compared",
                "to the reference sample r (which can be arbitrarily",
                "chosen). To do that, M values are computed as",
                "log2(counts_sample = counts_reference). If no bias is",
                "present, it should be expected that the median of M",
                "values for each comparison is 0. Otherwise, it would be",
                "indicating that expression levels in one of the samples",
                "tend to be higher than in the other, and this could lead",
                "to false discoveries when computing differencial expression.",
                "Confidence intervals for the M median are also computed by",
                "bootstrapping. If value 0 does not fall inside the interval,",
                "it means that the deviation of the sample with regard",
                "to the reference sample is statistically significant.",
                "Therefore, a normalization procedure is required.'</em>",
                collapse=" "
                    ),
                boxplot=paste(
                "The boxplot comprises a means of summarizing the read",
                "counts distribution of a sample in the form of a bar",
                "with extending lines, as a commonly used way of",
                "graphically presenting groups of numerical data. A",
                "boxplot also indicates which observations, if any, might",
                "be considered outliers and is able to visually show",
                "different types of populations, without making any",
                "assumptions about the underlying statistical distribution.",
                "The spacings between the different parts of the box help",
                "indicate variance, skewness and identify outliers. The",
                "thick bar inside the colored box is the median of the",
                "observations while the box extends over the Interquartile",
                "Range of the observations. The whiskers extend up (down)",
                "to +/-1.5xIQR. Similar boxplots indicate good quality of",
                "normalization. If boxplots remain dissimilar after",
                "normalization, another normalization algorithm may have to be",
                "examined. The un-normalized boxplots show the need for data",
                "normalization in order for the data from different",
                "samples to follow the same underlying distribution and",
                "statistical testing to become possible.",collapse=" "
                ),
                    gcbias=paste(
                "The GC-content bias plot is a quality control chart that",
                "shows the possible dependence of the read counts (in log2",
                "scale) under a gene to the GC content percentage of that",
                "gene. In order for the statistical tests to be able to",
                "detect statistical significance which occurs due to real",
                "biological effects and not through other systematic biases",
                "present in the data (e.g. possible GC-content bias),",
                "the latter should be accounted for by the applied",
                "normalization algorithm. Although the tests are performed",
                "for each gene across biological conditions one could assume",
                "that the GC content does not represent a bias, as it is the",
                "same for the tested gene across samples and conditions.",
                "However, Risso et al. (2011) showed that GC-content could",
                "could have an impact on the statistical testing procedure.",
                "The GC-content bias plot depicts the dependence of the",
                "read counts to the GC content before and after normalization.",
                "The smoothing lines for each sample, should be as 'straight'",
                "as possible after normalization. In addition, if the",
                "smoothing lines differ significantly between biological",
                "conditions, this would constitute a possible quality warning.",
                collapse=" "
                    ),
                    lengthbias=paste(
                "The gene/transcript length bias plot is a quality control",
                "chart that shows the possible dependence of read counts (in",
                "log2 scale) under a gene to the length of that gene (whole",
                "gene or sum of exons depending on the analysis). In order",
                "for the statistical tests to be able to detect statistical",
                "significance which occurs due to real biological effects",
                "and not by other systematic biases present in the data",
                "(e.g. possible length bias), the latter should be accounted",
                "for by the applied normalization algorithm. Although the",
                "tests are performed for each gene across biological",
                "conditions, one could assume that the gene length does not",
                "represent a bias, as it is the same for the tested gene",
                "across samples and conditions. However, it has been shown in",
                "several studies that gene length could have an impact on the",
                "statistical testing procedure. The length bias plot",
                "depicts the dependence of the read counts to the",
                "gene/transcript length before and after normalization.",
                "The smoothing lines for each sample, should be as 'straight'",
                "as possible after normalization. In addition, if the",
                "smoothing lines differ significantly between biological",
                "conditions, this would constitute a possible quality warning.",
                collapse=" "
                    ),
                    meandiff=paste(
                "A mean-difference plot (or a Bland-Altman plot) is a",
                "method of data plotting used in analyzing the agreement",
                "between two different assays/variables. In this graphical",
                "method the differences (or alternatively the ratios)",
                "between the two variables are plotted against the averages",
                "of the two. Such a plot is useful, for example, for analyzing",
                "data with strong correlation between the x and y axes, when",
                "the (x,y) dots on the plot are close to the diagonal x=y.",
                "In this case, the value of the transformed variable X is",
                "approximately the same as x and y and variable Y shows the",
                "difference between x and y. When the data cloud in a mean",
                "difference plot is centered around the horizontal zero line,",
                "this is an indication of good data quality and good",
                "normalization results. On the other hand, when the data",
                "cloud deviates from the center line or has a 'banana'",
                "shape, this constitutes an indication of systematic biases",
                "present in the data and that either the chosen normalization",
                "algorithm has not worked well, or that data are not",
                "normalized. The smoothing curve that traverses the data",
                "(red curve) summarizes the above trends.",collapse=" "
                    ),
                    meanvar=paste(
                "The mean-variance plot comprises a graphical means of",
                "displaying a possible relationship between the means of",
                "gene expression (counts) values and their variances",
                "across replicates of a gene expression experiment. Thus",
                "data can be inspected for possible overdispersion (greater",
                "variability in a dataset than would be expected based on",
                "a given simple statistical model). In such plots for",
                "RNA-Seq data, overdispersion is usually manifested as",
                "increasing variance with increasing gene expression",
                "(counts) and it is summarized through a smoothing curve",
                "(red curve). The following is taken from the EDASeq package",
                "vignette: '<em>...although the Poisson distribution",
                "is a natural and simple way to model count data, it has",
                "the limitation of assuming equality of the mean and",
                "variance. For this reason, the negative binomial",
                "distribution has been proposed as an alternative when the",
                "data show over-dispersion...'</em> If overdispersion is",
                "not present, the data cloud is expected to be evenly",
                "scattered around the smoothing curve.",collapse=" "
                    ),
                    deheatmap=paste(
                    "Differentially Expressed Genes (DEGs) heatmaps depict",
                    "how well samples from different conditions cluster",
                    "together according to their expression values after",
                    "normalization and statistical testing, for each requested",
                    "statistical contrast. If samples from the same biological",
                    "condition do not cluster together, this would constitute",
                    "a warning sign regarding the quality of the samples. In",
                    "addition, DEG heatmaps provide an initial view of",
                    "possible clusters of co-expressed genes.",collapse=" "
                    ),
                    volcano=paste(
                "A volcano plot is a scatterplot that is often used when",
                "analyzing high-throughput -omics data (e.g. microarray",
                "data, RNA-Seq data) to give an overview of interesting",
                "genes. The log2 fold change is plotted on the x-axis and",
                "the negative log10 p-value is plotted on the y-axis. A",
                "volcano plot combines the results of a statistical test",
                "(aka, p-values) with the magnitude of the change enabling",
                "quick visual identification of those genes that display",
                "large-magnitude changes and that are also statistically",
                "significant. The horizontal dashed line sets the threshold",
                "for statistical significance, while the vertical dashed",
                "lines set the thresholds for biological significance. It",
                "should be noted that volcano plots become harder to",
                "interpret when using more than one statistical algorithm",
                "and performing meta-analysis. This happens because the genes",
                "that have stronger evidence of being differentially",
                "expressed obtain lower p-values while the rest either",
                "remain at similar levels or obtain higher p-values.",
                "The result is a 'warped' volcano plot, with two",
                "main data clouds: one in the upper part of the plot, and",
                "one in the lower part of the plot. You can always zoom in",
                "when using interactive mode (the default).",collapse=" "
                    ),
                 mastat=paste(
                 "A mean-difference (or MA) plot with overlaid statistical",
                 "information (p-value and fold change thresholds manifested",
                 "as points with different colors) is a very useful graphic",
                 "that enables the visualization of the results of",
                 "differential expression analysis. It differs from the",
                 "volcano plot regarding what is displayed in the axes system.",
                 "While a volcano plot displays the fold change (x-axis)",
                 "versus the statistical significance (y-axis), an MA plot",
                 "with statistical scores depicts average expression over the",
                 "biological conditions that are compared (x-axis) versus the",
                 "fold change of the comparison. Statistical significance",
                 "categorization is added as point coloring and statistical",
                 "significance is indicated only by different colors and not",
                 "by the position to the axes system as in the volcano plot.",
                 "This plot is useful when it is of little interest how",
                 "statistically significant a gene/transcript is (we are",
                 "interested only in the fact that it is) but someone is",
                 "interested in actual expression and fold change values",
                 "instead.",collapse=" "
                 ),
                deregulogram=paste(
                "The de-regulogram is a scatterplot of fold changes between ",
                "two different contrasts. It depicts whether the DEGs between",
                "the two selected contrasts follow a concordant or discordant",
                "regulation pattern. For each (common) DEG, the x-axis and",
                "y-axis represent the log<sub>2</sub> fold change of the two",
                "contrasts. The location of each point along the four",
                "quartiles can directly show its regulation pattern in the two",
                "comparisons. Therefore, the dots localized in the second or",
                "the fourth quartile, illustrate DEGs with a common regulation",
                "pattern, while those localized in the first or third quartile",
                "represent DEGs with opposite patterns of regulation."
                 ),
                    biodist=paste(
                "The chromosome and biotype distributions bar diagram for",
                "Differentially Expressed Genes (DEGs) is split in two",
                "panels: i)in the upper panel DEGs are distributed per",
                "chromosome and the percentage of each chromosome in the",
                "genome is presented in red bars, the percentage of DEGs",
                "in each chromosome is presented in green bars and the",
                "percentage of certain chromosomes in the distribution of",
                "DEGs is presented in blue bars; ii)in the lower panel,",
                "DEGs are distributed per biotype and the percentage of",
                "each biotype in the genome (i.e. in the whole set of",
                "features provided, for example, protein coding genes, non",
                "coding RNAs or pseudogenes) is presented in red bars,",
                "the percentage of DEGs in each biotype is presented in",
                "green bars and the percentage of each biotype in DEGs is",
                "presented in blue lines. The vertical line separates the most",
                "abundant biotypes (on the left-hand side, corresponding to",
                "the left axis scale), from the rest(on the right-hand side,",
                "corresponding to the right axis scale). Otherwise, the lower",
                "abundance, biotypes would be indistinguishable.",collapse=" "
                    ),
                    filtered=paste(
                "The chromosome and biotype distribution of filtered genes",
                "is a quality control chart with two rows and four panels:",
                "on the left panel of the first row, the bar chart depicts",
                "the numbers of filtered genes per chromosome (actual numbers",
                "shown above the bars). On the right panel of the first row,",
                "the bar chart depicts the numbers of filtered genes per",
                "biotype (actual numbers shown above the bars). On the left",
                "panel of the second row, the bar chart depicts the fraction",
                "of filtered genes to the total genes per chromosome",
                "(actual percentages shown above the bars). On the right",
                "panel of the second row, the bar chart depicts the fraction",
                "of the filtered genes to the total genes per biotype",
                "(actual percentages shown above the bars). This plot",
                "should indicate possible quality problems when for example",
                "the filtered genes for a specific chromosome (or the",
                "fraction) is much higher than the rest. Generally,",
                "the fractions per chromosome should be uniform and the",
                "fractions per biotype should be proportional to the biotype",
                "fraction relative to the genome.",collapse=" "
                    ),
                    statvenn=paste(
                    "Venn diagrams are an intuitive way of presenting",
                    "overlaps between lists, based on the overlap of basic",
                    "geometrical shapes. The numbers of overlapping genes per",
                    "statistical algorithm are shown in the different areas",
                    "of the Venn diagrams, one for each contrast. Apart from",
                    "a p-value cutoff, a fold change threshold of 0.5 in",
                    "log<sub>2</sub> scale is applied for each contrast. For",
                    "multi-condition contrasts, the first condition is used to",
                    "calculate the fold change against the reference.",
                    collapse=" "
                    ),
                    foldvenn=paste(
                    "Venn diagrams are an intuitive way of presenting",
                    "overlaps between lists, based on the overlap of basic",
                    "geometrical shapes. The numbers of overlapping genes per",
                    "statistical contrast are shown in the different areas",
                    "of the Venn diagrams, one for each contrast. Apart from",
                    "a p-value cutoff, a fold change threshold of 0.5 in",
                    "log<sub>2</sub> scale is applied for each contrast. For",
                    "multi-condition contrasts, the first condition is used to",
                    "calculate the fold change against the reference.",
                    collapse=" "
                    )
                ),
                references=list(
                    main=paste("Moulos, P., Hatzis, P. (2015). Systematic",
                        "integration of RNA-Seq statistical algorithms for",
                        "accurate detection of differential gene expression",
                        "patterns. Nucleic Acids Research 43(4), e25."),
                    filein=list(
                        sam=paste("Statham, A.L., Strbenac, D., Coolen, M.W.,",
                            "Stirzaker, C., Clark, S.J., Robinson, M.D.",
                            " (2010). Repitools: an R package for the analysis",
                            "of enrichment-based epigenomic data.",
                            "Bioinformatics 26(13), 1662-1663."),
                        bam=paste("Statham, A.L., Strbenac, D., Coolen, M.W.,",
                            "Stirzaker, C., Clark, S.J., Robinson, M.D. (2010)",
                            "Repitools: an R package for the analysis of",
                            "enrichment-based epigenomic data. Bioinformatics",
                            "26(13), 1662-1663."),
                        bed=paste("Lawrence, M., Gentleman, R., Carey, V.",
                            "(2009). rtracklayer: an R package for interfacing",
                            "with genome browsers. Bioinformatics 25(14),",
                            "1841-1842.")
                    ),
                    norm=list(
                        edaseq=paste("Risso, D., Schwartz, K., Sherlock, G.,",
                            "and Dudoit, S. (2011). GC-content normalization",
                            "for RNA-Seq data. BMC Bioinformatics 12, 480."),
                        deseq=paste("Anders, S., and Huber, W. (2010).",
                            "Differential expression analysis for sequence",
                            "count data. Genome Biol 11, R106."),
                        deseq2=paste("Love, M.I., Huber, W., Anders, S.",
                            "Moderated estimation of fold change and",
                            "dispersion for RNA-seq data with DESeq2.",
                            "Genome Biology 15(12):550 (2014)"), 
                        edger=paste("Robinson, M.D., McCarthy, D.J., and",
                            "Smyth, G.K. (2010). edgeR: a Bioconductor package",
                            "for differential expression analysis of digital",
                            "gene expression data. Bioinformatics 26,",
                            "139-140."),
                        noiseq=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        nbpseq=paste("Di, Y, Schafer, D., Cumbie, J.S., and",
                            "Chang, J.H. (2011). The NBP Negative Binomial",
                            "Model for Assessing Differential Gene Expression",
                            "from RNA-Seq. Statistical Applications in",
                            "Genetics and Molecular Biology 10(1), 1-28."),
                        absseq=paste("Wentao Yang, Philip Rosenstiel and",
                             "Hinrich Schulenburg: ABSSeq: a new RNA-Seq",
                             "analysis method based on modelling absolute",
                             "expression differences. BMC Genomics 2016;",
                             "17: 541"),
                        dss=paste("Hao Wu, Chi Wang, Zhijin Wu (2013):",
                             "A new shrinkage estimator for dispersion",
                             "improves differential expression detection in",
                             "RNA-seq data. Biostatistics, 14(2):232-43."),
                        none=NULL
                    ),
                    stat=list(
                        deseq=paste("Anders, S., and Huber, W. (2010).",
                            "Differential expression analysis for sequence",
                            "count data. Genome Biol 11, R106."),
                        deseq2=paste("Love, M.I., Huber, W., Anders, S.",
                            "Moderated estimation of fold change and",
                            "dispersion for RNA-seq data with DESeq2.",
                            "Genome Biology 15(12):550 (2014)"),
                        edger=paste("Robinson, M.D., McCarthy, D.J., and",
                            "Smyth, G.K. (2010). edgeR: a Bioconductor package",
                            "for differential expression analysis of digital",
                            "gene expression data. Bioinformatics 26,",
                            "139-140."),
                        noiseq=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        limma=paste("Smyth, G. (2005). Limma: linear models",
                            "for microarray data. In Bioinformatics and",
                            "Computational Biology Solutions using R and",
                            "Bioconductor, G. R., C. V., D. S., I. R., and",
                            "H. W., eds. (New York, Springer), pp. 397-420."),
                        bayseq=paste("Hardcastle, T.J., and Kelly, K.A.",
                            "(2010). baySeq: empirical Bayesian methods for",
                            "identifying differential expression in sequence",
                            "count data. BMC Bioinformatics 11, 422."),
                        nbpseq=paste("Di, Y, Schafer, D., Cumbie, J.S., and",
                            "Chang, J.H. (2011). The NBP Negative Binomial",
                            "Model for Assessing Differential Gene Expression",
                            "from RNA-Seq. Statistical Applications in",
                            "Genetics and Molecular Biology 10(1), 1-28."),
                        ebseq=paste("Leng, N., Dawson, J.A., Thomson, J.A.,",
                            "Ruotti, V., Rissman, A.I., Smits, B.M., Haag,",
                            "J.D., Gould, M.N., Stewart, R.M., and",
                            "Kendziorski, C. (2013). EBSeq: an empirical",
                            "Bayes hierarchical model for inference in",
                            "RNA-seq experiments. Bioinformatics 29, ",
                            "1035-1043"),
                        absseq=paste("Wentao Yang, Philip Rosenstiel and",
                             "Hinrich Schulenburg: ABSSeq: a new RNA-Seq",
                             "analysis method based on modelling absolute",
                             "expression differences BMC Genomics 2016; 17:",
                             "541"),
                        dss=paste("Hao Wu, Chi Wang, Zhijin Wu (2013):",
                             "A new shrinkage estimator for dispersion",
                             "improves differential expression detection in",
                             "RNA-seq data. Biostatistics, 14(2):232-43.",
                             "doi:10.1093/biostatistics/kxs033")
                    ),
                    meta=list(
                        fisher=paste("Fisher, R.A. (1932). Statistical",
                            "Methods for Research Workers (Edinburgh, Oliver",
                            "and Boyd)."),
                        fperm=paste("Fisher, R.A. (1932). Statistical",
                            "Methods for Research Workers (Edinburgh, Oliver",
                            "and Boyd)."),
                        whitlock=c(
                            paste("Whitlock, M.C. (2005). Combining",
                                "probability from independent tests:",
                                "the weighted Z-method is superior to Fisher's",
                                "approach. J Evol Biol 18, 1368-1373."),
                            paste("Schroder, M.S., Culhane, A.C., Quackenbush,",
                                "J., and Haibe-Kains, B. (2011). survcomp:",
                                "an R/Bioconductor package for performance",
                                "assessment and comparison of survival",
                                "models. Bioinformatics 27, 3206-3208.")
                        ),
                        weight=paste("Genovese, C.R., Roeder, K., Wasserman,",
                            "L. (2006). False discovery control with p-value",
                            "weighting. Biometrika 93 (3): 509-524."),
                        simes=paste("Simes, R. J. (1986). An improved",
                            "Bonferroni procedure for multiple tests of",
                            "significance. Biometrika 73 (3): 751-754."),
                        none=NULL
                    ),
                    multiple=list(
                        BH=paste("Benjamini, Y., and Hochberg, Y. (1995). ",
                            "Controlling the False Discovery Rate: A Practical",
                            "and Powerful Approach to Multiple Testing.",
                            "Journal of the Royal Statistical Society Series",
                            "B (Methodological) 57, 289-300."),
                        fdr=paste("Benjamini, Y., and Hochberg, Y. (1995). ",
                            "Controlling the False Discovery Rate: A Practical",
                            "and Powerful Approach to Multiple Testing.",
                            "Journal of the Royal Statistical Society Series",
                            "B (Methodological) 57, 289-300."),
                        BY=paste("Benjamini, Y., and Yekutieli, D. (2001). The",
                            "control of the false discovery rate in multiple",
                            "testing under dependency. Annals of Statistics",
                            "26, 1165-1188."),
                        bonferroni=paste("Shaffer, J.P. (1995). Multiple",
                            "hypothesis testing. Annual Review of",
                            "Psychology 46, 561-576."),
                        holm=paste("Holm, S. (1979). A simple sequentially",
                            "rejective multiple test procedure. Scandinavian",
                            "Journal of Statistics 6, 65-70."),
                        hommel=paste("Hommel, G. (1988). A stagewise rejective",
                            "multiple test procedure based on a modified",
                            "Bonferroni test. Biometrika 75, 383-386."),
                        hochberg=paste("Hochberg, Y. (1988). A sharper",
                            "Bonferroni procedure for multiple tests of",
                            "significance. Biometrika 75, 800-803."),
                        qvalue=paste("Storey, J.D., and Tibshirani, R. (2003).",
                            "Statistical significance for genomewide studies.",
                            "Proc Natl Acad Sci U S A 100, 9440-9445.")
                    ),
                    figure=list(
                        mds=paste("Planet, E., Attolini, C.S., Reina, O.,",
                            "Flores, O., and Rossell, D. (2012). htSeqTools:",
                            "high-throughput sequencing quality control,",
                            "processing and visualization in R. Bioinformatics",
                            "28, 589-590."),
                        biodetection=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        countsbio=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        saturation=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        readnoise=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        gcbias=paste("Risso, D., Schwartz, K., Sherlock, G.,",
                            "and Dudoit, S. (2011). GC-content normalization",
                            "for RNA-Seq data. BMC Bioinformatics 12, 480."),
                        lengthbias=paste("Risso, D., Schwartz, K., Sherlock,",
                            "G., and Dudoit, S. (2011). GC-content",
                            "normalization for RNA-Seq data.",
                            "BMC Bioinformatics 12, 480."),
                        meandiff=paste("Risso, D., Schwartz, K., Sherlock, G.,",
                            "and Dudoit, S. (2011). GC-content normalization",
                            "for RNA-Seq data. BMC Bioinformatics 12, 480."),
                        meanvar=paste("Risso, D., Schwartz, K., Sherlock, G.,",
                            "and Dudoit, S. (2011). GC-content normalization",
                            "for RNA-Seq data. BMC Bioinformatics 12, 480."),
                        rnacomp=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        biodist=paste("Tarazona, S., Garcia-Alcalde, F.,",
                            "Dopazo, J., Ferrer, A., and Conesa, A. (2011).",
                            "Differential expression in RNA-seq: a matter of",
                            "depth. Genome Res 21, 2213-2223."),
                        statvenn=paste("Chen, H., and Boutros, P.C. (2011).",
                            "VennDiagram: a package for the generation of",
                            "highly-customizable Venn and Euler diagrams in R.",
                            "BMC Bioinformatics 12, 35."),
                        foldvenn=paste("Chen, H., and Boutros, P.C. (2011).",
                            "VennDiagram: a package for the generation of",
                            "highly-customizable Venn and Euler diagrams in R.",
                            "BMC Bioinformatics 12, 35."),
                        filtered=NULL
                    )
                )
            )
        }
    )
    return(messages)
}

#makeHighchartsPoints <- function(x,y,a=NULL) {
#    if (length(x)>0) {
#        n <- names(x)
#        x <- unname(x)
#        y <- unname(y)
#        stru <- vector("list",length(x))
#        if (is.null(a)) {
#            for (i in seq_len(length(x)))
#                stru[[i]] <- list(
#                    x=round(x[i],digits=3),
#                    y=round(y[i],digits=3),
#                    name=n[i]
#                )
#        }
#        else {
#            for (i in seq_len(length(x)))
#                stru[[i]] <- list(
#                    x=round(x[i],digits=3),
#                    y=round(y[i],digits=3),
#                    name=n[i],
#                    alt_name=a[i]
#                )
#        }
#    }
#    else
#        stru <- list(x=NULL,y=NULL,name=NULL,alt_name=NULL)
#    return(stru)
#}

makeHighchartsPoints <- function(x,y,a=NULL,p=NULL,simple=FALSE) {
    if (length(x)>0) {
        n <- names(x)
        x <- unname(x)
        y <- unname(y)
        if (!is.null(a))
            a <- unname(a)
        if (!is.null(p))
            p <- unname(p)
        #stru <- vector("list",length(x))
        if (simple) {
            return(lapply(seq_along(x),function(i,x,y) {
                return(c(x[i],y[i]))
           },x,y))
        }
        if (is.null(a) && is.null(p)) {
            stru <- lapply(seq_along(x),function(i,x,y,n) {
                return(list(
                    x=round(x[i],digits=3),
                    y=round(y[i],digits=3),
                    name=n[i]
                ))
           },x,y,n)
        }
        else if (!is.null(a) && is.null(p)) {
            stru <- lapply(seq_along(x),function(i,x,y,n,a) {
                return(list(
                    x=round(x[i],digits=3),
                    y=round(y[i],digits=3),
                    name=n[i],
                    alt_name=a[i]
                ))
           },x,y,n,a)
        }
        else if (is.null(a) && !is.null(p)) {
            stru <- lapply(seq_along(x),function(i,x,y,n,p) {
                return(list(
                    x=round(x[i],digits=3),
                    y=round(y[i],digits=3),
                    name=n[i],
                    sig=round(p[i],digits=3)
                ))
           },x,y,n,p)
        }
        else if (!is.null(a) && !is.null(p)) {
            stru <- lapply(seq_along(x),function(i,x,y,n,a,p) {
                return(list(
                    x=round(x[i],digits=3),
                    y=round(y[i],digits=3),
                    name=n[i],
                    alt_name=a[i],
                    sig=round(p[i],digits=3)
                ))
           },x,y,n,a,p)
        }
    }
    else
        stru <- list(x=NULL,y=NULL,name=NULL,alt_name=NULL,sig=NULL)
    return(stru)
}

asClassVector <- function(sampleList) {
    classes <- vector("list",length(sampleList))
    names(classes) <- names(sampleList)
    for (n in names(sampleList))
        classes[[n]] <- rep(n,times=length(sampleList[[n]]))
    classes <- unlist(classes,use.names=FALSE)
    names(classes) <- unlist(sampleList,use.names=FALSE)
    return(classes)
}

getArg <- function(argList,argName) {
    return(argList[argName])
}

setArg <- function(argList,argName,argValue=NULL) {
    if (is.list(argName))
        argList[names(argName)] <- argName
    else if (is.character(argName)) {
        tmp <- vector("list",length(argName))
        names(tmp) <- argName
        i <- 0
        for (n in argName) {
            i <- i + 1
            tmp[[n]] <- argValue[i]
        }
        argList[argName] <- tmp
    }
    return(argList)
}

wpAdjust <- function(p,m) {
    if (m=="qvalue")
        return(qvalue(p))
    else
        return(p.adjust(p,method=m))
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            else 
                m <- FALSE
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

filterLow <- function(x,f) { return(all(x<=f)) }

filterHigh <- function(x,f) { return(all(x>=f)) }

disp <- function(...) {
    #verbose <- get("VERBOSE",envir=metaEnv)
    verbose <- TRUE
    if (!is.null(verbose) && verbose) {
        message("\n",...,appendLF=FALSE)
    }
    logger <- get("LOGGER",envir=metaEnv)
    levalias <- c("one","two","three","four","five")
    if (!is.null(logger)) {
        switch(levalias[level(logger)],
            one = { debug(logger,paste0(...)) },
            two = { info(logger,gsub("\\n","",paste0(...))) },
            three = { warn(logger,gsub("\\n","",paste0(...))) },
            four = { error(logger,gsub("\\n","",paste0(...))) },
            five = { fatal(logger,gsub("\\n","",paste0(...))) }
        )
    }
}

stopwrap <- function(...,t="fatal") {
    logger <- get("LOGGER",envir=metaEnv)
    if (!is.null(logger)) {
        if (t=="fatal")
            fatal(logger,gsub("\\n","",paste0(...)))
        else
            error(logger,gsub("\\n","",paste0(...)))
    }
    stop(paste0(...))
}

warnwrap <- function(...,now=FALSE) {
    logger <- get("LOGGER",envir=metaEnv)
    if (!is.null(logger))
        warn(logger,gsub("\\n","",paste0(...)))
    if (now)
        warning(paste0(...),call.=FALSE,immediate.=TRUE)
    else
        warning(paste0(...),call.=FALSE)
}

exampleCountData <- function(ngenes) {
    if (missing(ngenes))
        ngenes <- 10000
    q0 <- rexp(ngenes,rate=1/250)
    is_DE <- runif(ngenes) < 0.3
    lfc <- rnorm(ngenes,sd=2)
    q0A <- ifelse(is_DE,q0*2^(lfc/2),q0)
    q0B <- ifelse(is_DE,q0*2^(-lfc/2),q0)
    true_sf <- c(1.0, 1.3,0.7,0.9,1.6)
    conds <- c("A","A","B","B","B")
    m <- t(vapply(seq_len(ngenes),function(i)
        vapply(seq_len(5),function(j)
            rnbinom(1,mu=true_sf[j]*ifelse(conds[j]=="A",q0A[i],q0B[i] ),
                size=1/0.2),numeric(1)),numeric(5)))
    colnames(m) <- c("A1","A2","B1","B2","B3")
    rownames(m) <- paste("gene",seq_len(ngenes),
        ifelse(is_DE,"T","F"),sep="_")
    return(m)
}

elap2human <- function(start.time) {
    start.time <- as.POSIXct(start.time)
    dt <- difftime(Sys.time(),start.time,units="secs")
    ndt <- as.numeric(dt)
    if (ndt<60)
        format(.POSIXct(dt,tz="GMT"),"%S seconds")
    else if (ndt>=60 && ndt<3600)
        format(.POSIXct(dt,tz="GMT"),"%M minutes %S seconds")
    else if (ndt>=3600 && ndt<86400)
        format(.POSIXct(dt,tz="GMT"),"%H hours %M minutes %S seconds")
    else if (ndt>=86400)
        format(.POSIXct(dt,tz="GMT"),"%d days %H hours %M minutes %S seconds")
}

.deprecatedWarning <- function(func) {
    switch(func,
        readTargets = {
            warnwrap("\"yes\" and \"no\" for read strandedness have been ",
                "deprecated. Please use \"forward\", \"forward\" or \"no\". ",
                "Replacing \"yes\" with \"forward\"...")
        }
    )
}

## The old function. It has NULLs in the new fields that are filled in
## getWeights2 (see furhter down).
#..getWeightsOld <- function(org=c("human","chimpanzee","mouse","fruitfly",
#    "arabidopsis","rat")) {
#    org <- tolower(org)
#    checkTextArgs("org",org,c("human","chimpanzee","mouse","fruitfly",
#        "arabidopsis","rat"))
#    switch(org,
#        human = {
#            return(c(
#                deseq=0.05772458,
#                deseq2=NULL,
#                edger=0.14321672,
#                limma=0.34516089,
#                nbpseq=0.06108182,
#                noiseq=0.11595169,
#                bayseq=0.27686431,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        chimpanzee = {
#            return(c(
#                deseq=0.06026782,
#                deseq2=NULL,
#                edger=0.14964358,
#                limma=0.33500306,
#                nbpseq=0.05814585,
#                noiseq=0.11337043,
#                bayseq=0.28356925,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        mouse = {
#            return(c(
#                deseq=0.05257695,
#                deseq2=NULL,
#                edger=0.24161354,
#                limma=0.29957277,
#                nbpseq=0.04914485,
#                noiseq=0.06847809,
#                bayseq=0.28861381,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        fruitfly = {
#            return(c(
#                deseq=0.01430269,
#                deseq2=NULL,
#                edger=0.12923339,
#                limma=0.38315685,
#                nbpseq=0.01265952,
#                noiseq=0.06778537,
#                bayseq=0.39286218,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        arabidopsis = {
#            return(c(
#                deseq=0.04926122,
#                deseq2=NULL,
#                edger=0.10130858,
#                limma=0.40842011,
#                nbpseq=0.04596652,
#                noiseq=0.09336509,
#                bayseq=0.30167848,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        chimp = {
#            return(c(
#                deseq=NULL,
#                deseq2=NULL,
#                edger=NULL,
#                limma=NULL,
#                nbpseq=NULL,
#                noiseq=NULL,
#                bayseq=NULL,
#                absseq=NULL,
#                dss=NULL
#            ))
#        },
#        rat = {
#            return(c(
#                deseq=NULL,
#                deseq2=NULL,
#                edger=NULL,
#                limma=NULL,
#                nbpseq=NULL,
#                noiseq=NULL,
#                bayseq=NULL,
#                absseq=NULL,
#                dss=NULL
#            ))
#        }
#    )
#}
