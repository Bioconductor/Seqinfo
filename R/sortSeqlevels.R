### =========================================================================
### sortSeqlevels()
### -------------------------------------------------------------------------


setGeneric("sortSeqlevels", signature="x",
    function(x, X.is.sexchrom=NA) standardGeneric("sortSeqlevels")
)

setMethod("sortSeqlevels", "character",
    function(x, X.is.sexchrom=NA)
    {
        x[order(rankSeqlevels(x, X.is.sexchrom=X.is.sexchrom))]
    }
)

setMethod("sortSeqlevels", "ANY",
    function(x, X.is.sexchrom=NA)
    {
        seqlevels(x) <- sortSeqlevels(seqlevels(x),
                                      X.is.sexchrom=X.is.sexchrom)
        x
    }
)

