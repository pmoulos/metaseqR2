test_metaseqr <- function() {
    data("mm9GeneData",package="metaseqR2")
    exDir <- tempdir()

    result1 <- metaseqr2(
        counts=mm9GeneCounts,
        sampleList=sampleListMm9,
        contrast=c("adult_8_weeks_vs_e14.5"),
        libsizeList=libsizeListMm9,
        annotation="embedded",
        org="mm9",
        countType="gene",
        normalization="edger",
        statistics=c("edger","limma"),
        metaP="simes",
        preset="medium_basic",
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
