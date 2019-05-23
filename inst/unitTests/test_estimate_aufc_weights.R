test_estimate_aufc_weights <- function() {
    data("mm9.gene.data",package="metaseqR")
    rc <- NULL
    weights <- estimateAufcWeights(
       counts=as.matrix(mm9.gene.counts[,9:12]),
       normalization="edaseq",
       statistics=c("edger","limma"),
       nsim=1,N=10,ndeg=c(2,2),top=4,modelOrg="mm9",
       seed=42,rc=rc,libsizeGt=1e+5
    )
    checkEqualsNumeric(weights,c(0.5384615,0.4615385),tolerance=1e-5)
    checkEqualsNumeric(sum(weights),1,tolerance=1e-9)
}
