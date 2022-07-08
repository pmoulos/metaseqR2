test_estimate_aufc_weights <- function() {
    data("mm9GeneData",package="metaseqR2")
    set.seed(42)
    weights <- estimateAufcWeights(
       counts=as.matrix(mm9GeneCounts[,9:12]),
       normalization="edaseq",
       statistics=c("edger","limma"),
       nsim=1,N=10,ndeg=c(2,2),top=4,modelOrg="mm10",
       rc=0.01,libsizeGt=1e+5
    )
    #checkEqualsNumeric(weights,c(0.5384615,0.4615385),tolerance=1e-5)
    checkEqualsNumeric(weights,c(0.5,0.5),tolerance=1e-5)
    checkEqualsNumeric(sum(weights),1,tolerance=1e-9)
}
