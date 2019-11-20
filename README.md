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
