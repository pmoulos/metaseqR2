#' @docType data
#' @name mm9GeneCounts
#' @title mouse RNA-Seq data with two conditions, four samples
#' @description This data set contains RNA-Seq gene read counts for 3 
#' chromosomes.
#' The data were downloaded from the ENCODE public repository and are derived
#' from the study of Mortazavi et al., 2008 (Mortazavi A, Williams BA, McCue K, 
#' Schaeffer L, Wold B. Mapping and quantifying mammalian transcriptomes by 
#' RNA-Seq.
#' Nat Methods. 2008 Jul;5(7):621-8). In their experiment, the authors studied
#' among others genes expression at two developmental stages of mouse liver 
#' cells.
#' It has two conditions-developmental stages (e14.5,
#' adult_8_weeks) and four samples (e14.5_1, e14.5_2, a8w_1, a8w_2). It also
#' contains a predefined \code{sampleList} and \code{libsizeList}
#' named \code{sampleListMm9} and \code{libsizeListMm9}.
#' @usage mm9.gene.counts
#' @format a \code{data.frame} with gene read counts and some embedded 
#' annotation, one row per gene.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name sampleListMm9
#' @title Mouse RNA-Seq data with two conditions, four samples
#' @description The sample list for \code{mm9GeneCounts}. See the data set
#' description.
#' @usage sampleListMm9
#' @format a named \code{list} with condition and sample names.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name libsizeListMm9
#' @title Mouse RNA-Seq data with two conditions, four samples
#' @description The library size list for \code{mm9GeneCounts}. See the data set
#' description.
#' @usage libsizeListMm9
#' @format a named \code{list} with library sizes.
#' @source ENCODE (http://genome.ucsc.edu/encode/)
#' @author Panagiotis Moulos
NULL
