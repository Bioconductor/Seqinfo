\name{sortSeqlevels}

\alias{sortSeqlevels}
\alias{sortSeqlevels,character-method}
\alias{sortSeqlevels,ANY-method}

\title{Sort the sequence levels of an object}

\description{
  A generic function and methods to sort the sequence levels of an
  object in "natural" order.
}

\usage{
sortSeqlevels(x, X.is.sexchrom=NA)
}

\arguments{
  \item{x}{
    Any object containing sequence information i.e. with a \code{seqinfo()}
    component.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome
    or to chromosome with Roman Numeral X. If \code{NA}, \code{sortSeqlevels}
    does its best to "guess".
  }
}

\value{
  The input object \code{x} with its seqlevels sorted.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \code{\link{seqlevels}} getter and setter.

    \item \code{\link{rankSeqlevels}}, on which \code{sortSeqlevels} is
          based.

    \item \link{Seqinfo} objects.

    \item \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
          objects in the \pkg{GenomicRanges} package.
  }
}

\examples{
sortSeqlevels(c("11", "Y", "1", "10", "9", "M", "2"))

seqlevels <- c("chrXI", "chrY", "chrI", "chrX", "chrIX", "chrM", "chrII")
sortSeqlevels(seqlevels)
sortSeqlevels(seqlevels, X.is.sexchrom=TRUE)
sortSeqlevels(seqlevels, X.is.sexchrom=FALSE)

seqlevels <- c("chr2RHet", "chr4", "chrUextra", "chrYHet",
               "chrM", "chrXHet", "chr2LHet", "chrU",
               "chr3L", "chr3R", "chr2R", "chrX")
sortSeqlevels(seqlevels)

library(GenomicRanges)
gr <- GRanges()
seqlevels(gr) <- seqlevels
sortSeqlevels(gr)
}

\keyword{methods}
