#' @docType data
#' @name hg19.exon.counts
#' @title Human RNA-Seq data with three conditions, three samples
#' @description This data set contains RNA-Seq exon read counts for 3 chromosomes.
#' Data are derived from three colon tissue types (normal, paracancerous, cancerous).
#' It contains three coditions (normal, paracancerous, cancerous) with one replicate
#' each (three samples in total). It also contains a predefined \code{sample.list}
#' and \code{libsize.list} named \code{sample.list.hg19} and \code{libsize.list.hg19}.
#' Data were downloaded from GEO (GSE33782) and the corresponding reference is
#' Wu et al., Transcriptome profiling of the cancer, adjacent non-tumor and distant
#' normal tissues from a colorectal cancer patient by deep sequencing. PLoS One
#' 2012, 7(8), e41001.
#' @usage hg19.exon.counts
#' @format a \code{data.frame} with exon read counts and some embedded annotation,
#' one row per exon.
#' @source GEO (http://www.ncbi.nlm.nih.gov/geo/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name sample.list.hg19
#' @title Human RNA-Seq data with three conditions, three samples
#' @description The sample list for \code{hg19.exon.counts}. See the data set
#' description.
#' @usage sample.list.hg19
#' @format a named \code{list} with condition and sample names.
#' @source GEO (http://www.ncbi.nlm.nih.gov/geo/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name libsize.list.hg19
#' @title Human RNA-Seq data with three conditions, three samples
#' @description The library size list for \code{hg19.exon.counts}. See the data
#' set description.
#' @usage libsize.list.hg19
#' @format a named \code{list} with library sizes.
#' @source GEO (http://www.ncbi.nlm.nih.gov/geo/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm9.gene.counts
#' @title mouse RNA-Seq data with two conditions, four samples
#' @description This data set contains RNA-Seq gene read counts for 3 chromosomes.
#' The data were downloaded from the ENCODE public repository and are derived
#' from the study of Mortazavi et al., 2008 (Mortazavi A, Williams BA, McCue K, 
#' Schaeffer L, Wold B. Mapping and quantifying mammalian transcriptomes by RNA-Seq.
#' Nat Methods. 2008 Jul;5(7):621-8). In their experiment, the authors studied
#' among others genes expression at two developmental stages of mouse liver cells.
#' It has two conditions-developmental stages (e14.5,
#' adult_8_weeks) and four samples (e14.5_1, e14.5_2, a8w_1, a8w_2). It also
#' contains a predefined \code{sample.list} and \code{libsize.list}
#' named \code{sample.list.mm9} and \code{libsize.list.mm9}.
#' @usage mm9.gene.counts
#' @format a \code{data.frame} with gene read counts and some embedded annotation,
#' one row per gene.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name sample.list.mm9
#' @title Mouse RNA-Seq data with two conditions, four samples
#' @description The sample list for \code{mm9.gene.counts}. See the data set
#' description.
#' @usage sample.list.mm9
#' @format a named \code{list} with condition and sample names.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name libsize.list.mm9
#' @title Mouse RNA-Seq data with two conditions, four samples
#' @description The library size list for \code{mm9.gene.counts}. See the data set
#' description.
#' @usage libsize.list.mm9
#' @format a named \code{list} with library sizes.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
