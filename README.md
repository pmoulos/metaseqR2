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

## metaseqR2 annotation database

If you do not wish to build annotation databases on your own using the
```buildAnnotationDatabase``` function, you can find complete pre-built 
annotation SQLite databases [here](https://drive.google.com/drive/folders/15lOY9PBggCcaoohO_0rQTvExXenqah55?usp=sharing). New versions will be constructed from time to time, most probably
whenever a new Ensembl release comes live.

The prebuilt annotations contain:

* Annotations for all supported organsisms for all types of metaseqR2 analyses.
* For every supported organism:
  + For the latest version of each genome, the latest two Ensembl required
  annotations.
  + For all other versions of each genome, the latest Ensembl required
  annotations supporting that particular version.
  + UCSC and RefSeq annotations as fetched in the day of the build (denoted
  by the folder name in the above link).
  
The SQLite database must be placed in ```system.file(package="metaseqR2")``` and
named ```annotation.sqlite```, that is
```file.path(system.file(package="metaseqR2"),"annotation.sqlite")```. Otherwise
you will have to provide your desired location in each ```metaseqr2``` call.
Alternatively, on-the-fly download is still supported but is inneficient.

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
* harmonicmeanp
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
* splines
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
