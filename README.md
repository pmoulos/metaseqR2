<!-- badges: start -->
  ![Bioconductor build](http://www.bioconductor.org/shields/build/devel/bioc/metaseqR2.svg)
  ![Bioconductor platforms](http://www.bioconductor.org/shields/availability/3.12/metaseqR2.svg)
  ![Bioconductor dependencies](http://www.bioconductor.org/shields/dependencies/devel/metaseqR2.svg)
  </br>
  ![GitHub](https://img.shields.io/github/license/pmoulos/metaseqR2)
  ![GitHub repo size](https://img.shields.io/github/repo-size/pmoulos/metaseqR2)
  ![GitHub issues](https://img.shields.io/github/issues/pmoulos/metaseqR2)
<!-- badges: end -->

# metaseqR2

An R package for the analysis, meta-analysis and result reporting of RNA-Seq 
gene expression data - Next Generation!

## Citation

metaseqR2 along with further research regarding the abilities of the PANDORA 
algorithm was published in [Briefings in Bioinformatics](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa156/5890504). If you use metaseqR2 in your research, please cite:

> Dionysios Fanidis, Panagiotis Moulos: **Integrative, normalization-insusceptible statistical analysis of RNA-Seq data, with improved differential expression and unbiased downstream functional analysis**, *Briefings in Bioinformatics*, 2020, bbaa156, DOI: [10.1093/bib/bbaa156](https://doi.org/10.1093/bib/bbaa156)

## Installation from Bioconductor

```
if (!requireNamespace("BiocManager",quietly=TRUE))
    install.packages("BiocManager")

library(BiocManager)

BiocManager::install("metaseqR2")

# or for development version to be installed
# BiocManager::install("metaseqR2",version="devel")
```

## Installation from GitHub

Use with caution as the latest version may be unstable, although typical
Bioconductor checks are executed before each push.

```
if (!requireNamespace("devtools",quietly=TRUE))
    install.packages("devtools")

library(devtools)
install_github("pmoulos/metaseqR2")
```

## Installation from source

The same things apply regarding stability.

```
git clone https://github.com/pmoulos/metaseqR2.git
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
R CMD INSTALL ./metaseqR2_x.y.z.tar.gz
```

Please report any issues [here](https://github.com/pmoulos/metaseqR2-local/issues). 

## metaseqR2 annotation database

If you do not wish to build annotation databases on your own using the
```buildAnnotationDatabase``` function, you can find complete pre-built 
annotation SQLite databases
[here](https://drive.google.com/drive/folders/15lOY9PBggCcaoohO_0rQTvExXenqah55?usp=sharing). 
New versions will be constructed from time to time, most probably whenever a new
Ensembl release comes live.

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
* DESeq2
* DSS
* DT
* EDASeq
* edgeR
* harmonicmeanp
* genefilter
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
* Matrix
* methods
* NBPSeq
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
