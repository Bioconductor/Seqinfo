### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter
###

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
        standardGeneric("seqinfo<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames() getter and setter
###

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", signature="x",
    function(x, value) standardGeneric("seqnames<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels() getter and setter
###

setGeneric("seqlevels", function(x) standardGeneric("seqlevels"))

### Default "seqlevels" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlevels", "ANY", function(x) seqlevels(seqinfo(x)))

setGeneric("seqlevels<-", signature="x",
    function(x,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
        standardGeneric("seqlevels<-")
)

### NOT exported but used in GenomicFeatures.
### Return -3L for "renaming" mode, -2L for "strict subsetting" mode (no new
### seqlevels added), -1L for "extended subsetting" mode (new seqlevels added),
### or an integer vector containing the mapping from the new to the old
### seqlevels for "general" mode (i.e. a combination of renaming and/or
### subsetting). Note that the vector describing the "general" mode is
### guaranteed to contain no negative values.
getSeqlevelsReplacementMode <- function(new_seqlevels, old_seqlevels)
{
    if (!is.character(new_seqlevels)
     || anyNA(new_seqlevels)
     || anyDuplicated(new_seqlevels))
        stop(wmsg("the supplied 'seqlevels' must be a character vector ",
                  "with no NAs and no duplicates"))
    nsl_names <- names(new_seqlevels)
    if (!is.null(nsl_names)) {
        nonempty_names <- nsl_names[!(nsl_names %in% c(NA, ""))]
        if (any(duplicated(nonempty_names)) ||
            length(setdiff(nonempty_names, old_seqlevels)) != 0L)
            stop(wmsg("the names of the supplied 'seqlevels' contain ",
                      "duplicates or invalid sequence levels"))
        return(match(nsl_names, old_seqlevels))
    }
    if (all(new_seqlevels %in% old_seqlevels))
        return(-2L)
    if (length(new_seqlevels) != length(old_seqlevels))
        return(-1L)
    is_renamed <- new_seqlevels != old_seqlevels
    tmp <- intersect(new_seqlevels[is_renamed], old_seqlevels[is_renamed])
    if (length(tmp) != 0L)
        return(-1L)
    return(-3L)
}

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
    {
        ## Make the new Seqinfo object.
        x_seqinfo <- seqinfo(x)
        seqlevels(x_seqinfo) <- value
        ## Map the new sequence levels to the old ones.
        new2old <- getSeqlevelsReplacementMode(value, seqlevels(x))
        if (identical(new2old, -3L)) {
            ## "renaming" mode
            new2old <- seq_along(value)
        } else if (identical(new2old, -2L) || identical(new2old, -1L)) {
            ## "subsetting" mode
            new2old <- match(value, seqlevels(x))
        }
        ## Do the replacement.
        seqinfo(x, new2old=new2old, pruning.mode=pruning.mode) <- x_seqinfo
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels0() getter
###
### Currently applicable to TxDb objects only.
###

setGeneric("seqlevels0", function(x) standardGeneric("seqlevels0"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### restoreSeqlevels()
###
### Currently applies to TxDb only.
###

restoreSeqlevels <- function(x)
{
    seqlevels(x) <- seqlevels0(x)
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlengths() getter and setter
###

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

### Default "seqlengths" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", signature="x",
    function(x, value) standardGeneric("seqlengths<-")
)

### Default "seqlengths<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlengths", "ANY",
    function(x, value)
    {
        seqlengths(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCircular() getter and setter
###

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

### Default "isCircular" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", signature="x",
    function(x, value) standardGeneric("isCircular<-")
)

### Default "isCircular<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("isCircular", "ANY",
    function(x, value)
    {
        isCircular(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genome() getter and setter
###

setGeneric("genome", function(x) standardGeneric("genome"))

### Default "genome" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("genome", "ANY", function(x) genome(seqinfo(x)))

setGeneric("genome<-", signature="x",
    function(x, value) standardGeneric("genome<-")
)

### Default "genome<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("genome", "ANY",
    function(x, value)
    {
        genome(seqinfo(x)) <- value
        x
    }
)

