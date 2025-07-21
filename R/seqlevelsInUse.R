### =========================================================================
### seqlevelsInUse()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevelsInUse() getter
###

setGeneric("seqlevelsInUse", function(x) standardGeneric("seqlevelsInUse"))

### Covers GenomicRanges, SummarizedExperiment, GAlignments, and any object
### for which the seqnames are returned as a factor-Rle.
setMethod("seqlevelsInUse", "Vector",
    function(x)
    {
        f <- runValue(seqnames(x))
        levels(f)[tabulate(f, nbins=nlevels(f)) != 0L]
    }
)

### Covers GRangesList and GAlignmentsList objects.
setMethod("seqlevelsInUse", "CompressedList",
    function(x) seqlevelsInUse(x@unlistData)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A bunch of (non-exported) helper functions aimed at facilitating the
### implementation of seqinfo() setter methods
###

normarg_new2old <- function(new2old, new_N, old_N,
                            new_what="the supplied 'seqinfo'",
                            old_what="the current 'seqinfo'")
{
    if (!is.numeric(new2old))
        stop(wmsg("'new2old' must be NULL or an integer vector"))
    if (length(new2old) != new_N)
        stop(wmsg("when not NULL, 'new2old' must ",
                  "have the length of ", new_what))
    if (!is.integer(new2old))
        new2old <- as.integer(new2old)
    min_new2old <- suppressWarnings(min(new2old, na.rm=TRUE))
    if (min_new2old != Inf) {
        if (min_new2old < 1L)
            stop(wmsg("when not NULL, 'new2old' must ",
                      "contain positive values or NAs"))
        if (max(new2old, na.rm=TRUE) > old_N)
            stop(wmsg("'new2old' cannot contain values ",
                      "greater than the length of ", old_what))
    }
    if (any(duplicated(new2old) & !is.na(new2old)))
        stop(wmsg("non-NA values in 'new2old' must be unique"))
    new2old
}

### Validate and reverse the 'new2old' mapping.
.reverse_new2old <- function(new2old, new_N, old_N,
                             new_what="the supplied 'seqinfo'",
                             old_what="the current 'seqinfo'")
{
    if (is.null(new2old))
        return(NULL)
    new2old <- normarg_new2old(new2old, new_N, old_N, new_what, old_what)
    S4Vectors:::reverseIntegerInjection(new2old, old_N)
}       

### NOT exported but used in GenomicRanges, SummarizedExperiment, and
### GenomicAlignments, in their seqinfo() setter methods.
### The dangling seqlevels in 'x' are those seqlevels that the user wants to
### drop but are in use.
getDanglingSeqlevels <- function(x, new2old=NULL,
                            pruning.mode=c("error", "coarse", "fine", "tidy"),
                            new_seqlevels)
{
    pruning.mode <- match.arg(pruning.mode)
    if (!is.character(new_seqlevels) || any(is.na(new_seqlevels)))
        stop(wmsg("the supplied 'seqlevels' must be a character vector ",
                  "with no NAs"))
    if (is.null(new2old))
        return(character(0))
    new_N <- length(new_seqlevels)
    old_seqlevels <- seqlevels(x)
    old_N <- length(old_seqlevels)
    old2new <- .reverse_new2old(new2old, new_N, old_N,
                                new_what="the supplied 'seqlevels'",
                                old_what="the current 'seqlevels'")
    seqlevels_to_drop <- old_seqlevels[is.na(old2new)]
    seqlevels_in_use <- seqlevelsInUse(x)
    dangling_seqlevels <- intersect(seqlevels_to_drop, seqlevels_in_use)
    if (length(dangling_seqlevels) != 0L && pruning.mode == "error")
        stop(wmsg("The following seqlevels are to be dropped but are ",
                  "currently in use (i.e. have ranges on them): ",
                  paste(dangling_seqlevels, collapse = ", "), ".\n",
                  "Please use the 'pruning.mode' argument to control how ",
                  "to prune 'x', that is, how to remove the ranges in 'x' ",
                  "that are on these sequences. For example, do something ",
                  "like:"),
             "\n    seqlevels(x, pruning.mode=\"coarse\") <- new_seqlevels",
             "\n  or:",
             "\n    keepSeqlevels(x, new_seqlevels, pruning.mode=\"coarse\")",
             "\n  See ?seqinfo for a description of the pruning modes.")
    dangling_seqlevels
}

### Compute the new seqnames resulting from new seqlevels.
### Assumes that 'seqnames(x)' is a 'factor' Rle (which is true if 'x' is a
### GRanges or GAlignments object, but not if it's a GRangesList object),
### and returns a 'factor' Rle of the same length (and same runLength vector).
### Always used in the context of the seqinfo() setter i.e. 'new_seqlevels'
### comes from 'seqlevels(value)' where 'value' is the supplied Seqinfo object.
makeNewSeqnames <- function(x, new2old=NULL, new_seqlevels)
{
    ## Should never happen.
    stopifnot(is.character(new_seqlevels), all(!is.na(new_seqlevels)))
    new_N <- length(new_seqlevels)
    old_N <- length(seqlevels(x))
    x_seqnames <- seqnames(x)
    if (!is.null(new2old)) {
        old2new <- .reverse_new2old(new2old, new_N, old_N,
                                    new_what="the supplied Seqinfo object",
                                    old_what="seqinfo(x)")
        tmp <- runValue(x_seqnames)
        levels(tmp) <- new_seqlevels[old2new]
        runValue(x_seqnames) <- factor(as.character(tmp), levels=new_seqlevels)
        return(x_seqnames)
    }
    if (new_N >= old_N &&
        identical(new_seqlevels[seq_len(old_N)], seqlevels(x)))
    {
        levels(x_seqnames) <- new_seqlevels
        return(x_seqnames)
    }
    SEQLEVELS_ARE_NOT_THE_SAME <- c(
        "The seqlevels in the supplied Seqinfo object ",
        "are not the same as the seqlevels in 'x'. "
    )
    if (length(intersect(seqlevels(x), new_seqlevels)) == 0L)
        stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
                  "Please use the 'new2old' argument to specify the ",
                  "mapping between the formers and the latters."))
    if (setequal(seqlevels(x), new_seqlevels))
        stop(wmsg("The seqlevels in the supplied Seqinfo object ",
                  "are not in the same order as the seqlevels in 'x'. ",
                  "Please reorder the seqlevels in 'x' with:"),
             "\n\n",
             "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
             wmsg("before calling the 'seqinfo()' setter."),
             "\n  ",
             wmsg("For any more complicated mapping between the new ",
                  "and old seqlevels (e.g. for a mapping that will ",
                  "result in the renaming of some seqlevels in 'x'), ",
                  "please use the 'new2old' argument."))
    if (all(seqlevels(x) %in% new_seqlevels))
        stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
                  "To map them to the seqlevels of the same name in 'x', ",
                  "the easiest way is to propagate them to 'x' with:"),
             "\n\n",
             "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
             wmsg("before calling the 'seqinfo()' setter."),
             "\n  ",
             wmsg("For any more complicated mapping, please use ",
                  "the 'new2old' argument."))
    stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
              "To map them to the seqlevels of the same name in 'x', ",
              "the easiest way is to propagate them to 'x' with:"),
         "\n\n",
         "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
         wmsg("before calling the 'seqinfo()' setter. ",
              "Note that you might need to specify a pruning mode ",
              "(via the 'pruning.mode' argument) if this operation ",
              "will drop seqlevels that are in use in 'x'."),
         "\n  ",
         wmsg("For any more complicated mapping, please use ",
              "the 'new2old' argument."))
}

### Returns a logical vector of the same length as 'new_seqinfo' indicating
### whether the length or circularity flag of the corresponding sequence
### have changed. Assumes that the 'new2old' mapping is valid (see
### .reverse_new2old() function above for what this means exactly).
### NAs in 'new2old' are propagated to the result.
sequenceGeometryHasChanged <- function(new_seqinfo, old_seqinfo, new2old=NULL)
{
    ans_len <- length(new_seqinfo)
    if (is.null(new2old)) {
        idx1 <- idx2 <- seq_len(ans_len)
    } else {
        idx1 <- which(!is.na(new2old))
        idx2 <- new2old[idx1]
    }
    new_seqlengths <- seqlengths(new_seqinfo)[idx1]
    old_seqlengths <- seqlengths(old_seqinfo)[idx2]
    new_isCircular <- isCircular(new_seqinfo)[idx1]
    old_isCircular <- isCircular(old_seqinfo)[idx2]
    hasNotChanged <- function(x, y)
        (is.na(x) & is.na(y)) | (is.na(x) == is.na(y) & x == y)
    ans <- logical(ans_len)
    ans[] <- NA
    ans[idx1] <- !(hasNotChanged(new_seqlengths, old_seqlengths) &
                   hasNotChanged(new_isCircular, old_isCircular))
    ans
}

