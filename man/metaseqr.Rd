\name{metaseqr}
\alias{metaseqr}
\alias{metaseqr.main}
\title{The main metaseqr pipeline}
\usage{
    metaseqr(counts, sampleList, excludeList = NULL,
        fileType = c("auto", "sam", "bam", "bed"),
        path = NULL, contrast = NULL, libsizeList = NULL,
        id.col = 4, gc.col = NA, name.col = NA, bt.col = NA,
        annotation = c("download", "embedded"), gene.file = NULL,
        org = c("hg18", "hg19", "hg38", "mm9", "mm10", "rn5", "rn6",
            "dm3", "dm6", "danrer7", "pantro4", "susscr3", "tair10", 
            "equcab2","custom"),
        refdb = c("ensembl", "ucsc", "refseq"),
        trans.level = c("gene", "transcript", "exon"),
        count.type = c("gene", "exon","utr"),
        utr.flank = 500,
        exon.filters = list(min.active.exons = list(exons.per.gene = 5, 
                min.exons = 2, frac = 1/5)),
        gene.filters = list(length = list(length = 500), 
                avg.reads = list(average.per.bp = 100, quantile = 0.25), 
                expression = list(median = TRUE, mean = FALSE, quantile = NA, 
                    known = NA, custom = NA), 
                biotype = get.defaults("biotype.filter", org[1]),
                presence = list(frac = 0.25, min.count = 10,
					per.condition = FALSE)),
        when.apply.filter = c("postnorm", "prenorm"),
        normalization = c("deseq", "edaseq", "edger", "noiseq", "nbpseq", 
                "each", "none"),
        norm.args = NULL,
        statistics = c("deseq", "edger", "noiseq", "bayseq", "limma", 
                "nbpseq"),
        stat.args = NULL,
        adjust.method = sort(c(p.adjust.methods, "qvalue")),
        meta.p = if (length(statistics) > 1) c("simes", "bonferroni", "fisher", 
                "dperm.min", "dperm.max", "dperm.weight", "fperm", "whitlock", 
                "minp", "maxp", "weight", "pandora", "none") else "none",
        weight = rep(1/length(statistics), length(statistics)),
        nperm = 10000, reprod=TRUE, pcut = NA, log.offset = 1,
        preset = NULL,
        qc.plots = c("mds", "biodetection", "countsbio", "saturation", 
                "readnoise", "filtered", "correl", "pairwise", "boxplot", 
                "gcbias", "lengthbias", "meandiff", "meanvar", "rnacomp", 
                "deheatmap", "volcano", "biodist"),
        fig.format = c("png", "jpg", "tiff", "bmp", "pdf", "ps"),
        out.list = FALSE, export.where = NA,
        export.what = c("annotation", "p.value", "adj.p.value", 
                "meta.p.value", "adj.meta.p.value", "fold.change", 
                "stats", "counts","flags"),
        export.scale = c("natural", "log2", "log10", "vst", "rpgm"),
        export.values = c("raw", "normalized"),
        export.stats = c("mean", "median", "sd", "mad", "cv", 
                "rcv"),
        export.counts.table = FALSE,
        restrict.cores = 0.6, report = TRUE, report.top = 0.1,
        report.template = "default", save.gene.model = TRUE,
        verbose = TRUE, run.log = TRUE, progress.fun = NULL, ...)
}
\arguments{
    \item{counts}{a text tab-delimited file containing gene
    or exon counts in one of the following formats: i) the
    first column contains unique gene or exon identifiers and
    the rest of the columns contain the read counts for each
    sample. Thus the first cell of each row is a gene or exon
    accession and the rest are integers representing the
    counts for that accession. In that case, the
    \code{annotation} parameter should strictly be
    \code{"download"} or an external file in proper format.
    ii) The first n columns should contain gene or exon
    annotation elements like chromosomal locations, gene
    accessions, exon accessions, GC content etc. In that
    case, the \code{annotation} parameter can also be
    \code{"embedded"}. The ideal embedded annotation contains
    8 columns, chromosome, gene or exon start, gene or exon
    end, gene or exon accession, GC-content (fraction or
    percentage), strand, HUGO gene symbol and gene biotype
    (e.g. "protein_coding" or "ncRNA"). When the
    \code{annotation} parameter is "embedded", certain of
    these features are mandatory (co-ordinates and
    accessions). If they are not present, the pipeline will
    not run. If additional elements are not present (e.g. GC
    content or biotypes), certain features of metaseqr will
    not be available. For example, EDASeq normalization will
    not be performed based on a GC content covariate but
    based on gene length which is not what the authors of
    EDASeq suggest. If biotypes are not present, a lot of
    diagnostic plots will not be available. If the HUGO gene
    symbols are missing, the final annotation will contain
    only gene accessions and thus be less comprehensible.
    Generally, it's best to set the \code{annotation}
    parameter to \code{"download"} to ensure the most
    comprehensible results. Counts can be a data frame 
    satisfying the above conditions. It is a data frame 
    by default when \code{read2count} is used. counts can 
    also be an .RData file (output of \code{\link{save}} 
    function which contains static input elements (list 
    containing the gene model (exon counts for each gene 
    constructed by the \code{\link{construct.gene.model}} 
    function, gene and exon annotation to avoid 
    re-downloading and/or gene counts depending on 
    \code{count.type}). This kind of input facilitates the 
    re-analysis of the same experiment, using different 
    filtering, normalization and statistical algorithms. 
    Finally, counts can be a list representing the gene model 
    (exon counts for each gene) constructed by the 
    \code{\link{construct.gene.model}} function (provided for 
    backwards compatibility). This .RData file can be generated 
    by setting \code{save.gene.model=TRUE} when performing data 
    analysis for the first time.}

    \item{sample.list}{a list containing condition names and
    the samples under each condition. It should have the
    format \code{sample.list <-}
    \code{list(ConditionA=c("Sample_A1",} \code{"Sample_A2",
    "Sample_A3"),} \code{ConditionB=c("Sample_B1",
    "Sample_B2"),} \code{ConditionC=c("Sample_C1",
    "Sample_C2"))}. The names of the samples in list members
    MUST match the column names containing the read counts in
    the counts file. If they do not match, the pipeline will
    either crash or at best, ignore several of your samples.
    Alternative, \code{sample.list} can be a small
    tab-delimited file structured as follows: the first line
    of the external tab delimited file should contain column
    names (names are not important). The first column MUST
    contain UNIQUE sample names and the second column MUST
    contain the biological condition where each of the
    samples in the first column should belong to. In this
    case, the function \code{\link{make.sample.list}} is
    used. If the \code{counts} argument is missing, the
    \code{sample.list} argument MUST be a targets text
    tab-delimited file which contains the sample names, the
    BAM/BED file names and the biological conditions/groups
    for each sample/file. The file should be text
    tab-delimited and structured as follows: the first line
    of the external tab delimited file should contain column
    names (names are not important). The first column MUST
    contain UNIQUE sample names. The second column MUST
    contain the raw BAM/BED files WITH their full path.
    Alternatively, the \code{path} argument should be
    provided (see below). The third column MUST contain the
    biological condition where each of the samples in the
    first column should belong to.}
    
    \item{exclude.list}{ a list of samples to exclude, in the 
    same (list) format as \code{sample.list} above.}

    \item{path}{an optional path where all the BED/BAM files
    are placed, to be prepended to the BAM/BED file names in
    the targets file. If not given and if the files in the
    second column of the targets file do not contain a path
    to a directory, the current directory is assumed to be
    the BAM/BED file container.}

    \item{file.type}{the type of raw input files. It can be
    \code{"auto"} for auto-guessing, \code{"bed"} for BED
    files, \code{"sam"} for SAM files or \code{"bam"} for BAM
    files.}

    \item{contrast}{a character vector of contrasts to be
    tested in the statistical testing step(s) of the metaseqr
    pipeline. Each element of contrast should STRICTLY have the
    format "ConditionA_vs_ConditionB_vs_...". A valid example
    based on the \code{sample.list} above is \code{contrast
    <- c("ConditionA_vs_ConditionB",}
    \code{"ConditionA_vs_ConditionC",}
    \code{"ConditionA_vs_ConditionB_vs_ConditionC")}. The
    first element of pairwise contrasts (e.g. "ConditionA"
    above) MUST be the control condition or any reference
    that ConditionB is checked against. metaseqr uses this
    convention to properly calculate fold changes. If it's
    NULL, a contrast between the first two members of the
    \code{sample.list} will be auto-generated.}

    \item{libsize.list}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sample.list}) and members are the library sizes
    (the sequencing depth) for each sample. For example
    \code{libsize.list <- list(Sample_A1=32456913,}
    \code{Sample_A2=4346818)}.}

    \item{id.col}{an integer denoting the column number in
    the file (or data frame) provided with the counts
    argument, where the unique gene or exon accessions are.
    Default to \code{4} which is the standard feature name
    column in a BED file.}

    \item{gc.col}{an integer denoting the column number in
    the file (or data frame) provided with the \code{counts}
    argument, where each gene's GC content is given. If not
    provided, GC content normalization provided by EDASeq
    will not be available.}

    \item{name.col}{an integer denoting the column number in
    the file (or data frame) provided with the counts
    argument, where the HUGO gene symbols are given. If not
    provided, it will not be available when reporting
    results. In addition, the \code{"known"} gene filter will
    not be available.}

    \item{bt.col}{an integer denoting the column number in
    the file (or data frame) provided with the counts
    argument, where the gene biotypes are given. If not
    provided, the \code{"biodetection"}, \code{"countsbio"},
    \code{"saturation"}, \code{"filtered"} and
    \code{"biodist"} plots will not be available.}

    \item{annotation}{instructs metaseqr where to find the
    annotation for the given counts file. It can be one of i)
    \code{"download"} (default) for automatic downloading of
    the annotation for the organism specified by the org
    parameter (using biomaRt), ii) \code{"embedded"} if the
    annotation elements are embedded in the read counts file
    or iv) a file specified by the user which should be as
    similar as possible to the \code{"download"} case, in
    terms of column structure.}
    
    \item{gene.file}{an external gene annotation file 
    required when \code{annotation="embedded"}, 
    \code{count.type="exon"} or \code{count.type="utr"}
    and \code{org="custom"}. See the result of 
    \code{\link{get.annotation}} on the format of the 
    external gene file.}

    \item{org}{the supported organisms by metaseqr. These can
    be, for human genomes \code{"hg18"}, \code{"hg19"} or 
    \code{"hg38"},  for mouse genomes \code{"mm9"}, \code{"mm10"}, 
    for rat genomes \code{"rn5"} or \code{"rn6"}, for drosophila 
    genome \code{"dm3"} or \code{"dm6"}, for zebrafish genome 
    \code{"danrer7"}, for chimpanzee genome \code{"pantro4"}, 
    for pig genome \code{"susscr3"}, for Arabidopsis 
    thaliana genome \code{"tair10"} and for Equus caballus genome
    \code{"equcab2"}. Finally, \code{"custom"} will instruct
    metaseqR to completely ignore the \code{org} argument and
    depend solely on annotation file provided by the user.}

    \item{refdb}{the reference annotation repository from
    which to retrieve annotation elements to use with
    metaseqr. It can be one of \code{"ensembl"} (default),
    \code{"ucsc"} or \code{"refseq"}.}
	
	\item{trans.level}{perform differential expression 
	analysis at which transcriptional unit, can be one of
	\code{"gene"} (default), \code{"transcript"} for 
	reporting differential expression at the transcript
	level or \code{"exon"} for exon level.}
	
    \item{count.type}{the type of reads inside the counts
    file. It can be one of \code{"gene"}, \code{"exon"} or
    \code{"utr"} for quant seq (lexogen) protocol. This is
    a very important and mandatory parameter as it defines
    the course of the workflow.}
    
    \item{utr.flank}{the number of base pairs to flank the
    3' UTR of transcripts when analyzing Quant-Seq data.}

    \item{exon.filters}{a named list whose names are the
    names of the supported exon filters and its members the
    filter parameters. See section "Exon filters" below for
    details.}

    \item{gene.filters}{a named list whose names are the
    names of the supported gene filters and its members the
    filter parameters. See section "Gene filters" below for
    details.}

    \item{when.apply.filter}{a character string determining
    when to apply the exon and/or gene filters, relative to
    normalization. It can be \code{"prenorm"} to apply apply
    the filters and exclude genes from further processing
    before normalization, or \code{"postnorm"} to apply the
    filters after normalization (default). In the case of
    \code{when.apply.filter="prenorm"}, a first normalization
    round is applied to a copy of the gene counts matrix in
    order to derive the proper normalized values that will
    constitute the several expression-based filtering
    cutoffs.}

    \item{normalization}{the normalization algorithm to be
    applied on the count data. It can be one of
    \code{"edaseq"} for EDASeq normalization, \code{"deseq"}
    for the normalization algorithm (individual options 
    specified by the \code{norm.args} argument) in the DESq 
    package (default), \code{"edger"} for the
    normalization algorithms present in the edgeR package
    (specified by the \code{norm.args} argument),
    \code{"noiseq"} for the normalization algorithms present
    in the NOISeq package (specified by the \code{norm.args}
    argument), \code{"nbpseq"} for the normalization
    algorithms present in the NBPSeq package (specified by
    the \code{norm.args} argument) or \code{"none"} to not
    normalize the data (highly unrecommended). It can also 
    be \code{"each"} where in this case, the normalization 
    applied will be specific to each statistical test used
    (i.e. the normalization method bundled with each package 
    and used in its examples and documentation). The last
    choice is for future use!}

    \item{norm.args}{a named list whose names are the names
    of the normalization algorithm parameters and its members
    parameter values. See section "Normalization parameters"
    below for details. Leave \code{NULL} for the defaults of
    \code{normalization}. If \code{normalization="each"}, it 
    must be a named list of lists, where each sub-list 
    contains normalization parameters specific to each 
    statistical test to be used. The last choice is for 
    future use!}

    \item{statistics}{one or more statistical analyses to be
    performed by the metaseqr pipeline.It can be one or more
    of \code{"deseq"} (default) to conduct statistical
    test(s) implemented in the DESeq package, \code{"edger"}
    to conduct statistical test(s) implemented in the edgeR
    package, \code{"limma"} to conduct the RNA-Seq version of
    statistical test(s) implemented in the limma package,
    \code{"noiseq"} to conduct statistical test(s)
    implemented in the NOISeq package, \code{"bayseq"} to
    conduct statistical test(s) implemented in the baySeq
    package and \code{"nbpseq"} to conduct statistical
    test(s) implemented in the NBPSeq package. In any case
    individual algorithm parameters are controlled by the
    contents of the \code{stat.args} list.}

    \item{stat.args}{a named list whose names are the names
    of the statistical algorithms used in the pipeline. Each
    member is another named list whose names are the
    algorithm parameters and its members are the parameter
    values. See section "Statistics parameters" below for
    details. Leave \code{NULL} for the defaults of
    \code{statistics}.}

    \item{adjust.method}{the multiple testing p-value
    adjustment method. It can be one of
    \code{\link{p.adjust.methods}} or \code{"qvalue"} from
    the qvalue Bioconductor package. Defaults to \code{"BH"}
    for Benjamini-Hochberg correction.}

    \item{meta.p}{the meta-analysis method to combine
    p-values from multiple statistical tests . It can be 
    one of \code{"simes"} (default), \code{"bonferroni"}, 
    \code{"minp"}, \code{"maxp"}, \code{"weight"}, \code{"pandora"}, 
    \code{"dperm.min"}, \code{"dperm.max"}, \code{"dperm.weight"}, 
    \code{"fisher"}, \code{"fperm"}, \code{"whitlock"} or 
    \code{"none"}. For the \code{"fisher"} and \code{"fperm"} 
    methods, see the documentation of the R package MADAM. For 
    the \code{"whitlock"} method, see the documentation of the 
    survcomp Bioconductor package. With the \code{"maxp"} 
    option, the final p-value is the maximum p-value out of 
    those returned by each statistical test. This is 
    equivalent to an "intersection" of the results derived 
    from each algorithm so as to have a final list with the 
    common genes returned by all statistical tests. Similarly, 
    when \code{meta.p="minp"}, is equivalent to a "union" of 
    the results derived from each algorithm so as to have a 
    final list with all the genes returned by all statistical 
    tests. The latter can be used as a very lose statistical 
    threshold to aggregate results from all methods regardless 
    of their False Positive Rate. With the \code{"simes"} 
    option, the method proposed by Simes (Simes, R. J., 1986)
    is used. With the \code{"dperm.min"}, \code{"dperm.max"}, 
    \code{"dperm.weight"} options, a permutation procedure is 
    initialed, where \code{nperm} permutations are performed 
    across the samples of the normalized counts matrix, 
    producing \code{nperm} permuted instances of the initital 
    dataset. Then, all the chosen statistical tests are 
    re-executed for each permutation. The final p-value is
    the number of times that the p-value of the permuted 
    datasets is smaller than the original dataset. 
    The p-value of the original dataset is created based on
    the choice of one of \code{dperm.min}, \code{dperm.max} 
    or \code{dperm.weight} options. In case of 
    \code{dperm.min}, the intial p-value vector is consists 
    of the minimum p-value resulted from the applied 
    statistical tests for each gene. The maximum p-value 
    is used with the \code{dperm.max} option. With the 
    \code{dperm.weight} option, the \code{weight} 
    weighting vector for each statistical test is used to 
    weight each p-value according to the power of 
    statistical tests (some might work better for a 
    specific dataset). Be careful as the permutation 
    procedure usually requires a lot of time. However, 
    it should be the most accurate. This method will NOT 
    work when there are no replicated samples across 
    biological conditions. In that case, use 
    \code{meta.p="simes"} instead. Finally, there are the 
    \code{"minp"}, \code{"maxp"} and \code{"weight"}
    options which correspond to the latter three methods 
    but without permutations. Generally, permutations 
    would be accurate to use when the experiment includes
    >5 samples per condition (or even better 7-10) which 
    is rather rare in RNA-Seq experiments. Finally, 
    \code{"pandora"} is the same as \code{"weight"} and is
    added to be in accordance with the metaseqR paper.}

    \item{weight}{a vector of weights with the same length as
    the \code{statistics} vector containing a weight for each
    statistical test. It should sum to 1. \strong{Use with
    caution with the} \code{dperm.weight} \strong{parameter!
    Theoretical background is not yet} \strong{solid and only
    experience shows improved results!}}

    \item{nperm}{the number of permutations performed to
    derive the meta p-value when \code{meta.p="fperm"} or
    \code{meta.p="dperm"}. It defaults to 10000.}

    \item{reprod}{create reproducible permutations when 
    \code{meta.p="dperm.min"}, \code{meta.p="dperm.max"} 
    or \code{meta.p="dperm.weight"}. Ideally one would 
    want to create the same set of indices for a given 
    dataset so as to create reproducible p-values. If 
    \code{reprod=TRUE}, a fixed seed is used by 
    \code{meta.perm} for all the datasets analyzed with 
    \code{metaseqr}. If \code{reprod=FALSE}, then the 
    p-values will not be reproducible, although statistical 
    significance is not expected to change for a large 
    number of resambling. Finally, \code{reprod} can be 
    a numeric vector of seeds with the same length as 
    \code{nperm} so that the user can supply his/her 
    own seeds.}

    \item{pcut}{a p-value cutoff for exporting differentially
    genes, default is to export all the non-filtered genes.}

    \item{log.offset}{an offset to be added to values during
    logarithmic transformations in order to avoid Infinity
    (default is \code{1}).}

    \item{preset}{an analysis strictness preset.
    \code{preset} can be one of \code{"all.basic"},
    \code{"all.normal"}, \code{"all.full"},
    \code{"medium.basic"}, \code{"medium.normal"},
    \code{"medium.full"}, \code{"strict.basic"},
    \code{"strict.normal"} or \code{"strict.full"}, each of
    which control the strictness of the analysis and the
    amount of data to be exported. For an explanation of the
    presets, see the section "Presets" below.}

    \item{qc.plots}{a set of diagnostic plots to show/create.
    It can be one or more of \code{"mds"},
    \code{"biodetection"}, \code{"rnacomp"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"readnoise"}, \code{"filtered"}, \code{"boxplot"},
    \code{"gcbias"}, \code{"lengthbias"}, \code{"meandiff"},
    \code{"meanvar"}, \code{"deheatmap"}, \code{"volcano"},
    \code{"biodist"}, \code{"venn"}. The \code{"mds"} stands
    for Mutlti-Dimensional Scaling and it creates a PCA-like
    plot but using the MDS dimensionality reduction instead.
    It has been succesfully used for NGS data (e.g. see the
    package htSeqTools) and it shows how well samples from
    the same condition cluster together. For
    \code{"biodetection"}, \code{"countsbio"},
    \code{"saturation"}, \code{"rnacomp"},
    \code{"readnoise"}, \code{"biodist"} see the vignette of
    NOISeq package. The \code{"saturation"} case has been
    rewritten in order to display more samples in a more
    simple way. See the help page of
    \code{\link{diagplot.noiseq.saturation}}. In addition,
    the \code{"readnoise"} plots represent an older version
    or the RNA composition plot included in older versions of
    NOISeq. For \code{"gcbias"}, \code{"lengthbias"},
    \code{"meandiff"}, \code{"meanvar"} see the vignette of
    EDASeq package. \code{"lenghtbias"} is similar to
    \code{"gcbias"} but using the gene length instead of the
    GC content as covariate. The \code{"boxplot"} option
    draws boxplots of log2 transformed gene counts. The
    \code{"filtered"} option draws a 4-panel figure with the
    filtered genes per chromosome and per biotype, as
    absolute numbers and as fractions of the genome. See also
    the help page of \code{\link{diagplot.filtered}}. The
    \code{"deheatmap"} option performs hierarchical
    clustering and draws a heatmap of differentially
    expressed genes. In the context of diagnostic plots, it's
    useful to see if samples from the same groups cluster
    together after statistical testing. The \code{"volcano"}
    option draws a volcano plot for each contrast and if a
    report is requested, an interactive volcano plot is
    presented in the HTML report. The \code{"venn"} option
    will draw an up to 5-way Venn diagram depicting the
    common and specific to each statistical algorithm genes
    and for each contrast, when meta-analysis is performed.
    The \code{"correl"} option creates two correlation
    graphs: the first one is a correlation heatmap (a
    correlation matrix which depicts all the pairwise
    correlations between each pair of samples in the counts
    matrix is drawn as a clustered heatmap) and the second
    one is a correlogram plot, which summarizes the
    correlation matrix in the form of ellipses (for an
    explanation please see the vignette/documentation of the
    R package corrplot. Set \code{qc.plots=NULL} if you don't
    want any diagnostic plots created.}

    \item{fig.format}{the format of the output diagnostic
    plots. It can be one or more of \code{"png"},
    \code{"jpg"}, \code{"tiff"}, \code{"bmp"}, \code{"pdf"},
    \code{"ps"}. The native format \code{"x11"} (for direct
    display) is not provided as an option as it may not
    render the proper display of some diagnostic plots in
    some devices.}

    \item{out.list}{a logical controlling whether to export a
    list with the results in the running environment.}

    \item{export.where}{an output directory for the project
    results (report, lists, diagnostic plots etc.)}

    \item{export.what}{the content of the final lists. It can
    be one or more of \code{"annotation"}, to bind the
    annoation elements for each gene, \code{"p.value"}, to
    bind the p-values of each method, \code{"adj.p.value"},
    to bind the multiple testing adjusted p-values,
    \code{"meta.p.value"}, to bind the combined p-value from
    the meta-analysis, \code{"adj.meta.p.value"}, to bind the
    corrected combined p-value from the meta-analysis,
    \code{"fold.change"}, to bind the fold changes of each
    requested contrast, \code{"stats"}, to bind several
    statistics calclulated on raw and normalized counts (see
    the \code{export.stats} argument), \code{"counts"}, to
    bind the raw and normalized counts for each sample.}

    \item{export.scale}{export values from one or more
    transformations applied to the data. It can be one or
    more of \code{"natural"}, \code{"log2"}, \code{"log10"},
    \code{"vst"} (Variance Stabilizing Transormation, see the
    documentation of DESeq package) and \code{"rpgm"} which
    is ratio of mapped reads per gene model (either the gene
    length or the sum of exon lengths, depending on 
    \code{count.type} argument). Note that this is not RPKM
    as reads are already normalized for library size using 
    one of the supported normalization methods. Also, 
    \code{"rpgm"} might be misleading when \code{normalization} 
    is other than \code{"deseq"}.}

    \item{export.values}{It can be one or more of
    \code{"raw"} to export raw values (counts etc.) and
    \code{"normalized"} to export normalized counts.}

    \item{export.stats}{calculate and export several
    statistics on raw and normalized counts, condition-wise.
    It can be one or more of \code{"mean"}, \code{"median"},
    \code{"sd"}, \code{"mad"}, \code{"cv"} for the
    Coefficient of Variation, \code{"rcv"} for a robust
    version of CV where the median and the MAD are used
    instead of the mean and the standard deviation.}

    \item{export.counts.table}{exports also the calculated 
    read counts table when input is read from bam files 
    and exports also the normalized count table in all 
    cases. Defaults to \code{FALSE}.}

    \item{restrict.cores}{in case of parallel execution of
    several subfunctions, the fraction of the available cores
    to use. In some cases if all available cores are used
    (\code{restrict.cores=1} and the system does not have
    sufficient RAM, the pipeline running machine might
    significantly slow down.}

    \item{report}{a logical value controlling whether to
    produce a summary report or not. Defaults to
    \code{TRUE}.}

    \item{report.top}{a fraction of top statistically 
    significant genes to append to the HTML report. This 
    helps in keeping the size of the report as small as 
    possible, as appending the total gene list might 
    create a huge HTML file. Users can always retrieve 
    the whole gene lists from the report links. Defaults 
    to \code{0.1} (top 10% of statistically significant 
    genes). Set to \code{NA} or \code{NULL} to append all 
    the statistically significant genes to the HTML report.}

    \item{report.template}{an HTML template to use for the
    report. Do not change this unless you know what you are
    doing.}

    \item{save.gene.model}{in case of exon analysis, a list
    with exon counts for each gene will be saved to the file
    \code{export.where/data/gene_model.RData}. This file can
    be used as input to metaseqR for exon count based 
    analysis, in order to avoid the time consuming step of 
    assembling the counts for each gene from its exons}

    \item{verbose}{print informative messages during
    execution? Defaults to \code{TRUE}.}

    \item{run.log}{write a log file of the \code{metaseqr}
    run using package log4r. Defaults to \code{TRUE}. The
    filename will be auto-generated.}
    
    \item{progress.fun}{a function which updates a 
    \code{Progress} object from shiny. This function must
    accept a \code{detail} argument. See 
    http://shiny.rstudio.com/articles/progress.html}

    \item{...}{further arguments that may be passed to
    plotting functions, related to \code{\link{par}}.}
}
\value{
    If \code{out.list} is \code{TRUE}, a named list whose
    length is the same as the number of requested contrasts.
    Each list member is named according to the corresponding
    contrast and contains a data frame of differentially
    expressed genes for that contrast. The contents of the
    data frame are defined by the \code{export.what,
    export.scale, export.stats, export.values} parameters. If
    \code{report} is \code{TRUE}, the output list contains
    two main elements. The first is described above (the
    analysis results) and the second contains the same
    results but in HTML formatted tables.
}
\description{
    This function is the main metaseqr workhorse and
    implements the main metaseqr workflow which performs data
    read, filtering, normalization and statistical selection,
    creates diagnostic plots and exports the results and a
    report if requested. The metaseqr function is responsible
    for assembling all the steps of the metaseqr pipeline
    which i) reads the input gene or exon read count table
    ii) performs prelimininary filtering of data by removing
    chrM and other non-essential information for a typical
    differential gene expression analysis as well as a
    preliminary expression filtering based on the exon
    counts, if an exon read count file is provided. iii)
    performs data normalization with one of currently widely
    used algorithms, including EDASeq (Risso et al., 2011),
    DESeq (Anders and Huber, 2010), edgeR (Robinson et al.,
    2010), NOISeq (Tarazona et al., 2012) or no normalization
    iv) performs a second stage of filtering based on the
    normalized gene expression according to several gene
    filters v) performs statistical testing with one or more
    of currently widely used algorithms, including DESeq
    (Anders and Huber, 2010), edgeR (Robinson et al., 2010),
    NOISeq (Tarazona et al., 2012), limma (Smyth et al.,
    2005) for RNA-Seq data, baySeq (Hardcastle et al., 2012)
    vi) in the case of multiple statistical testing
    algorithms, performs meta-analysis using one of five
    available methods (see the meta.p argument) vii) exports
    the resulting differentially expressed gene list in text
    tab-delimited format viii) creates a set of diagnostic
    plots either available in the aforementioned packages or
    metaseqr specific ones and ix) creates a comprehensive
    HTML report which summarizes the run information, the
    results and the diagnostic plots. Certain diagnostic
    plots (e.g. the volcano plot) can be interactive with the
    use of the external Highcharts
    (http://www.highcharts.com) JavaScript library for
    interactive graphs. Although the inputs to the metaseqr
    workflow are many, in practice, setting only very few of
    them and accepting the defaults as the rest can result in
    quite comprehensible results for mainstream organisms
    like mouse, human, fly and rat.
}
\note{
    Please note that currently only gene and exon annotation
    from Ensembl (http://www.ensembl.org), UCSC and RefSeq 
    are supported. Thus, the unique gene or exon ids in the 
    counts files should correspond to valid Ensembl, UCSC or 
    RefSeq gene or exon accessions for the organism of interest. 
    If you are not sure about the source of your counts file or 
    do not know how to produce it, it's better to start from the 
    original BAM/BED files (metaseqr will use the
    \code{\link{read2count}} function to create a counts
    file). Keep in mind that in the case of BED files, the
    performance will be significantly lower and the overall
    running time significanlty higher as the R functions
    which are used to read BED files to proper structures
    (GenomicRanges) and calculate the counts are quite slow.
    An alternative way is maybe the easyRNASeq package
    (Delhomme et al, 2012). The \code{\link{read2count}}
    function does not use this package but rather makes use
    of standard Bioconductor functions to handle NGS data. If
    you wish to work outside R, you can work with other
    popular read counters such as the HTSeq read counter
    (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html).
    Please also note that in the current version, the members
    of the \code{gene.filters} and \code{exon.filters} lists
    are not checked for validity so be careful to supply with
    correct names otherwise the pipeline will crash or at the
    best case scenario, will ignore the filters. Also note
    that when you are supplying metaseqr wtih an exon counts
    table, gene annotation is always downloaded so please be
    sure to have a working internet connection. In addition
    to the above, if you have a multiple core system, be very
    careful on how you are using the \code{restrict.cores}
    argument and generally how many cores you are using with
    scripts purely written in R. The analysis with exon read
    data can very easily cause memory problems, so unless you
    have more than 64Gb of RAM available, consider setting
    restrict.cores to something like 0.2 when working with
    exon data. Finally, if you do not wish to download the
    same annotation again and again when performing multiple
    analyses, it is best to use the
    \code{\link{get.annotation}} function to download and
    store the resulting data frames in local files and then
    use these files with the \code{annotation} option.

    Please note that the \strong{meta-analysis} feature
    provided by metaseqr does not satisfy the strict definition 
    of "meta-analysis", which is the combination of multiple 
    similar datasets under the same statistical methodology. 
    Instead it is the use of mulitple statistical tests applied 
    to the same data. For the Simes method, please consult also
    "Simes, R. J. (1986). "An improved Bonferroni procedure
    for multiple tests of significance". Biometrika 73 (3):
    751–754."
}
\section{Exon filters}{
    The exon filters are a set of filters which are applied
    after the gene models are assembled from the read counts
    of individual exons and before the gene expression is
    summarized from the exons belonging to each gene. These
    filters can be applied when the input read counts file
    contains exon reads. It is not applicable when the input
    file already contains gene counts. Such filters can be
    for example "accept genes where all the exons contain
    more than x reads" or "accept genes where there is read
    presence in at least m/n exons, n being the total exons
    of the gene". Such filters are NOT meant for detecting
    differential splicing as also the whole metaseqr
    pipeline, thus they should not be used in that context.
    The \code{exon.filters} argument is a named list of
    filters, where the names are the filter names and the
    members are the filter parameters (named lists with
    parameter name, parameter value). See the usage of the
    \code{metaseqr} function for an example of how these
    lists are structured. The supported exon filters in the
    current version are: i) \code{min.active.exons} which
    implements a filter for demanding m out of n exons of a
    gene to have a certain read presence with parameters
    \code{exons.per.gene}, \code{min.exons} and \code{frac}.
    The filter is described as follows: if a gene has up to
    \code{exons.per.gene} exons, then read presence is
    required in at least \code{min.exons} of them, else read
    presence is required in a \code{frac} fraction of the
    total exons. With the default values, the filter
    instructs that if a gene has up to 5 exons, read presence
    is required in at least 2, else in at least 20% of the
    exons, in order to be accepted. More filters will be
    implemented in future versions and users are encouraged
    to propose exon filter ideas to the author by mail. See
    \code{metaseqr} usage for the defaults. Set
    \code{exon.filters=NULL} to not apply any exon filtering.
}

\section{Gene filters}{
    The gene filters are a set of filters applied to gene
    expression as this is manifested through the read
    presence on each gene and are preferably applied after
    normalization. These filters can be applied both when the
    input file or data frame contains exon read counts and
    gene read counts. Such filter can be for example "accept
    all genes above a certain count threshold" or "accept all
    genes with expression above the median of the normalized
    counts distribution" or "accept all with length above a
    certain threshold in kb" or "exclude the 'pseudogene'
    biotype from further analysis". The supported gene
    filters in the current version, which have the same
    structure as the exon filters (named list of lists with
    filter names, parameter names and parameter arguments)
    are: i) \code{length} which implements a length filter
    where genes are accepted for further analysis if they are
    above \code{length} (its parameter) kb. ii)
    \code{avg.reads} which implements a filter where a gene
    is accepted for further analysis if it has more average
    reads than the \code{quantile} of the average count
    distribution per \code{average.per.bp} base pairs. In
    summary, the reads of each gene are averaged per
    \code{average.per.bp} based on each gene's length (in
    case of exons, input the "gene's length" is the sum of
    the lengths of exons) and the \code{quantile} quantile of
    the average counts distribution is calculated for each
    sample. Genes passing the filter should have an average
    read count larger than the maximum of the vector of the
    quantiles calculated above. iii) \code{expression} which
    implements a filter based on the overall expression of a
    gene. The parameters of this filter are: \code{median},
    where genes below the median of the overall count
    distribution are not accepted for further analysis (this
    filter has been used to distinguish between "expressed"
    and "not expressed" genes in several cases, e.g. (Mokry
    et al., NAR, 2011) with a logical as value, \code{mean}
    which is the same as \code{median} but using the mean,
    \code{quantile} which is the same as the previous two but
    using a specific quantile of the total counts
    distribution, \code{known}, where in this case, a set of
    known not-expressed genes in the system under
    investigation are used to estimate an expression cutoff.
    This can be quite useful, as the genes are filtered based
    on a "true biological" cutoff instead of a statistical
    cutoff. The value of this filter is a character vector of
    HUGO gene symbols (MUST be contained in the annotation,
    thus it's better to use \code{annotation="download"})
    whose counts are used to build a "null" expression
    distribution. The 90th quantile of this distribution is
    then the expression cutoff. This filter can be combined
    with any other filter. Be careful with gene names as they
    are case sensitive and must match exactly ("Pten" is
    different from "PTEN"!). iv) \code{biotype} where in this
    case, genes with a certain biotype (MUST be contained in
    the annotation, thus it's better to use
    \code{annotation="download"}) are excluded from the
    analysis. This filter is a named list of logical, where
    names are the biotypes in each genome and values are
    \code{TRUE} or \code{FALSE}. If the biotype should be
    excluded, the value should be \code{TRUE} else
    \code{FALSE}. See the result of
    \code{get.defaults("biotype.filter","hg19")} for an
    example. Finally, in future versions there will be
    support for user-defined filters in the form of a
    function. v) \code{presence} where in this case, a gene
    is further considered for statistical testing if 
    \code{frac} (x100 for a percentage value) have more
    than \code{min.count} reads across all samples
    (\code{per.condition=FALSE}) or across the samples
    of each condition (\code{per.condition=TRUE}).
}

\section{Normalization parameters}{
    The normalization parameters are passed again as a named
    list where the names of the members are the normalization
    parameter names and the values are the normalization
    parameter values. You should check the documentation of
    the packages EDASeq, DESeq, edgeR, NOISeq and NBPSeq for
    the parameter names and parameter values. There are a few
    exceptions in parameter names: in case of
    \code{normalization="edaseq"} the only parameter names
    are \code{within.which} and \code{between.which},
    controlling the withing lane/sample and between
    lanes/samples normalization algorithm. In the case 
    of \code{normalization="nbpseq"}, there is one
    additional parameter called \code{main.method} which can
    take the calues \code{"nbpseq"} or \code{"nbsmyth"}.
    These values correspond to the two different workflows
    available in the NBPSeq package. Please, consult the
    NBPSeq package documentation for further details. For the
    rest of the algorithms, the parameter names are the same
    as the names used in the respective packages. For
    examples, please use the \code{\link{get.defaults}}
    function.
}

\section{Statistics parameters}{
    The statistics parameters as passed to statistical
    algorithms in metaseqr, exactly with the same way as the
    normalization parametes above. In this case, there is one
    more layer in list nesting. Thus, \code{stat.args} is a
    named list whose names are the names the algorithms used
    (see the \code{statistics} parameter). Each member is
    another named list,with parameters to be used for each
    statistical algorithm. Again, the names of the member
    lists are parameter names and the values of the member
    lists are parameter values. You should check the
    documentations of DESeq, edgeR, NOISeq, baySeq, limma and
    NBPSeq for these parameters. There are a few exceptions
    in parameter names: In case of \code{statistics="edger"},
    apart from the rest of the edgeR statistical testing
    arguments, there is the argument \code{main.method} which
    can be either \code{"classic"} or \code{"glm"}, again
    defining whether the binomial test or GLMs will be used
    for statistical testing. For examples, please use the
    \code{\link{get.defaults}} function. When
    \code{statistics="nbpseq"}, apart from the rest arguments
    of the NBPSeq functions \code{estimate.disp} and
    \code{estimate.dispersion}, there is the argument
    \code{main.method} which can be \code{"nbpseq"} or
    \code{"nbsmyth"}. This argument determines the parameters
    to be used by the \code{estimate.dispersion} function or
    by the \code{estimate.disp} function to estimate RNA-Seq
    count dispersions. The difference between the two is that
    they constitute different starting points for the two
    workflows in the package NBPSeq. The first worklfow (with
    \code{main.method="nbpseq"} and the
    \code{estimate.dispersion} function is NBPSeq package
    specific, while the second (with
    \code{main.method="nbsmyth"} and the \code{estimate.disp}
    function is similar to the workflow of the edgeR package.
    For additional information regarding the statistical
    testing in NBPSeq, please consult the documentation of
    the NBPSeq package. \strong{Additinally, please note that
    there is currently a problem with the NBPSeq package and
    the workflow that is specific to the NBPSeq package. The
    problem has to do with function exporting as there are
    certain functions which are not recognized from the
    package internally. For this reason and until it is
    fixed, only the Smyth workflow will be available with the
    NBPSeq package (thus}
    \code{stat.args$main.method="nbpseq"} \strong{ will not
    be available)!}
}

\section{Presets}{
    The analysis presets are a set of keywords (only one can
    be used) that predefine some of the parameters of the
    metaseqr pipeline. For the time being they are quite
    simple and they control i) the strictness of filtering
    and statistical thresholding with three basic levels
    ("all", "medium", "strict") and ii) the data columns that
    are exported, again in three basic ways ("basic",
    "normal", "full") controlling the amount of data to be
    exported. These keywords can be combined with a dot in
    the middle (e.g. \code{"all.basic"} to define an analysis
    preset. When using analysis presets, the following
    argumentsof metaseqr are overriden: \code{exon.filters},
    \code{gene.filters}, \code{pcut}, \code{export.what},
    \code{export.scale}, \code{export.values},
    \code{exon.stats}. If you want to explicitly control the
    above arguments, the \code{preset} argument should be set
    to \code{NULL} (default). Following is a synopsis of the
    different presets and the values of the arguments they
    moderate: \itemize{ \item \code{"all.basic"}: use all
    genes (do not filter) and export all genes and basic
    annotation and statistics elements. In this case, the
    above described arguments become: \itemize{ \item
    \code{exon.filters=NULL} \item \code{gene.filters=NULL}
    \item \code{pcut=1} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change")} \item
    \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean")} } \item
    \code{"all.normal"}: use all genes (do not filter) and
    export all genes and normal annotation and statistics
    elements. In this case, the above described arguments
    become: \itemize{ \item \code{exon.filters=NULL} \item
    \code{gene.filters=NULL} \item \code{pcut=1} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean","sd","cv")} } In this case,
    the above described arguments become: \itemize{ \item
    \code{exon.filters=NULL} \item \code{gene.filters=NULL}
    \item \code{pcut=1} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item
    \code{export.scale=c("natural","log2","log10","vst")}
    \item \code{export.values=c("raw","normalized")} \item
    \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
    } \item \code{"medium.basic"}: apply a medium set of
    filters and and export statistically significant genes
    and basic annotation and statistics elements. In this
    case, the above described arguments become: \itemize{
    \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,min.exons=2,frac=1/5))} 
    \item \code{gene.filters=list(length=list(length=500),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.25),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.05} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change")} \item
    \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean")} } \item
    \code{"medium.normal"}: apply a medium set of filters and
    and export statistically significant genes and normal
    annotation and statistics elements. In this case, the
    above described arguments become: \itemize{ \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,min.exons=2,frac=1/5))} 
    \item \code{gene.filters=list(length=list(length=500),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.25),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.05} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean","sd","cv")} } and statistics
    elements. In this case, the above described arguments
    become: \itemize{ \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,min.exons=2,frac=1/5))} 
    \item \code{gene.filters=list(length=list(length=500),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.25),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.05} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item
    \code{export.scale=c("natural","log2","log10","vst")}
    \item \code{export.values=c("raw","normalized")} \item
    \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
    } \item \code{"strict.basic"}: apply a strict set of
    filters and and export statistically significant genes
    and basic annotation and statistics elements. In this
    case, the above described arguments become: \itemize{
    \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,min.exons=2,frac=1/4))} 
    \item \code{gene.filters=list(length=list(length=750),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.5),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.01} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change")} \item
    \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean")} } \item
    \code{"strict.normal"}: apply a strict set of filters and
    and export statistically significant genes and normal
    annotation and statistics elements. In this case, the
    above described arguments become: \itemize{ \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,min.exons=2,frac=1/4))} 
    \item \code{gene.filters=list(length=list(length=750),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.5),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.01} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item \code{export.scale=c("natural","log2")} \item
    \code{export.values=c("normalized")} \item
    \code{export.stats=c("mean","sd","cv")} } and statistics
    elements. In this case, the above described arguments
    become: \itemize{ \item
    \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,min.exons=2,frac=1/4))} 
    \item \code{gene.filters=list(length=list(length=750),} 
    \code{avg.reads=list(average.per.bp=100,quantile=0.5),} 
    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,known=NA,custom=NA),} 
    \code{biotype=get.defaults("biotype.filter",org[1]))}
    \item \code{pcut=0.01} \item
    \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",} 
    \code{"adj.meta.p.value","fold.change","stats","counts")}
    \item
    \code{export.scale=c("natural","log2","log10","vst")}
    \item \code{export.values=c("raw","normalized")} \item
    \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
    } }
}
\examples{
\donttest{
# An example pipeline with exon counts
data("hg19.exon.data",package="metaseqR")
metaseqr(
 counts=hg19.exon.counts,
 sample.list=list(normal="normal",paracancerous="paracancerous",cancerous="cancerous"),
 contrast=c("normal_vs_paracancerous","normal_vs_cancerous",
        "normal_vs_paracancerous_vs_cancerous"),
 libsize.list=libsize.list.hg19,
 id.col=4,
 annotation="download",
 org="hg19",
 count.type="exon",
 normalization="edaseq",
 statistics="deseq",
 pcut=0.05,
 qc.plots=c("mds", "biodetection", "countsbio", "saturation", "rnacomp", 
        "boxplot", "gcbias", "lengthbias", "meandiff", "readnoise","meanvar", 
        "readnoise", "deheatmap", "volcano", "biodist", "filtered"),
 fig.format=c("png","pdf"),
 export.what=c("annotation","p.value","adj.p.value","fold.change","stats",
        "counts"),
 export.scale=c("natural","log2","log10","vst"),
 export.values=c("raw","normalized"),
 export.stats=c("mean","median","sd","mad","cv","rcv"),
 restrict.cores=0.8,
 gene.filters=list(
     length=list(
         length=500
     ),
     avg.reads=list(
         average.per.bp=100,
         quantile=0.25
     ),
     expression=list(
         median=TRUE,
         mean=FALSE
     ),
     biotype=get.defaults("biotype.filter","hg18")
 )
)

# An example pipeline with gene counts
data("mm9.gene.data",package="metaseqR")
result <- metaseqr(
 counts=mm9.gene.counts,
 sample.list=list(e14.5=c("e14.5_1","e14.5_2"), adult_8_weeks=c("a8w_1","a8w_2")),
 contrast=c("e14.5_vs_adult_8_weeks"),
 libsize.list=libsize.list.mm9,
 annotation="download",
 org="mm9",
 count.type="gene",
 normalization="edger",
 statistics=c("deseq","edger","noiseq"),
 meta.p="fisher",
 pcut=0.05,
 fig.format=c("png","pdf"),
 export.what=c("annotation","p.value","meta.p.value","adj.meta.p.value",
        "fold.change"),
 export.scale=c("natural","log2"),
 export.values="normalized",
 export.stats=c("mean","sd","cv"),
 export.where=getwd(),
 restrict.cores=0.8,
 gene.filters=list(
     length=list(
         length=500
     ),
     avg.reads=list(
            average.per.bp=100,
            quantile=0.25
     ),
     expression=list(
            median=TRUE,
            mean=FALSE,
            quantile=NA,
            known=NA,
            custom=NA
     ),
     biotype=get.defaults("biotype.filter","mm9")
 ),
 out.list=TRUE
)
head(result$data[["e14.5_vs_adult_8_weeks"]])
}
}
\author{
    Panagiotis Moulos
}

