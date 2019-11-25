# metaseqR2-local

An R package for the analysis, meta-analysis and result reporting of RNA-Seq gene expression data - Next Generation!

## Installation of the pre-Bioconductor release

```
git clone https://github.com/pmoulos/metaseqR2-local.git
mkdir metaseqR2-build
rsync -avr --exclude=README.md --exclude=.git --exclude=.gitignore  \
    ./metaseqR2-local/ ./metaseqR2-build/metaseqR2
cd ./metaseqR2-build
R CMD build ./metaseqR2
```

This will take some time to build the vignettes. If you do not need them:

```
R CMD build --no-build-vignettes ./metaseqR2
```

And then install

```
R CMD INSTALL metaseqR2_0.0.1.tar.gz
```

Please report any issues [here](https://github.com/pmoulos/metaseqR2-local/issues). 

## List of required packages

metaseqR2 would benefit from the existence of all the following packages:

 * ABSSeq
 * baySeq
 * Biobase
 * BiocGenerics
 * BiocManager
 * BiocParallel
 * BiocStyle
 * biomaRt
 * Biostrings
 * BSgenome
 * corrplot
 * DESeq
 * DESeq2
 * DSS
 * DT
 * EDASeq
 * edgeR
 * GenomeInfoDb
 * GenomicAlignments
 * GenomicFeatures
 * GenomicRanges
 * gplots
 * graphics
 * grDevices
 * heatmaply
 * htmltools
 * httr
 * IRanges
 * jsonlite
 * knitr
 * limma
 * log4r
 * magrittr
 * methods
 * NBPSeq
 * NOISeq
 * pander
 * parallel
 * qvalue
 * rmarkdown
 * rmdformats
 * RMySQL
 * Rsamtools
 * RSQLite
 * rtracklayer
 * RUnit
 * S4Vectors
 * stats
 * stringr
 * SummarizedExperiment
 * survcomp
 * TCC
 * utils
 * VennDiagram
 * vsn
 * zoo

A recent version of [Pandoc](https://pandoc.org/) is also required, ideally
above 2.0.
