test_metaseqr <- function() {
    data("mm9.gene.data",package="metaseqR")
    exDir <- tempdir()

    result1 <- metaseqr(
        counts=mm9.gene.counts,
        sampleList=sample.list.mm9,
        contrast=c("e14.5_vs_adult_8_weeks"),
        libsizeList=libsize.list.mm9,
        annotation="download",
        org="mm9",
        countType="gene",
        normalization="edger",
        statistics=c("edger","limma"),
        metaP="simes",
        preset="medium.basic",
        qcPlots="mds",
        figFormat="png",
        exportWhere=exDir,
        outList=TRUE
    )
    checkTrue(file.exists(file.path(exDir,"index.html")))
    checkTrue(file.exists(file.path(exDir,"plots","qc","mds.png")))
    checkTrue(file.exists(file.path(exDir,"lists")))
    checkTrue(nrow(result1[[1]][[1]])>0)
    checkEqualsNumeric(ncol(result1[[1]][[1]]),16)
}
