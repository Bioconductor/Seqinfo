\name{seqlevelsInUse}

\alias{seqlevelsInUse}
\alias{seqlevelsInUse,Vector-method}
\alias{seqlevelsInUse,CompressedList-method}

\title{Get the sequence levels in use}

\description{
  Get the sequence levels that are "used" by an object, that is,
  the sequence levels on which the object defines genomic ranges or
  features.
}

\usage{
seqlevelsInUse(x)
}

\arguments{
  \item{x}{
    Any object containing sequence information i.e. with a \code{seqinfo()}
    component.
  }
}

\value{
  The sequence levels in use in a character vector.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \code{\link{seqlevels}} getter and setter.

    \item \link{Seqinfo} objects.

    \item \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
          objects in the \pkg{GenomicRanges} package.
  }
}

\examples{
library(GenomicRanges)

gr <- GRanges(rep(c("chr2", "chr3", "chrM"), 2), IRanges(1:6, 10))

## Add new seqlevels:
seqlevels(gr) <- c("chr1", seqlevels(gr), "chr4")
seqlevels(gr)
seqlevelsInUse(gr)

## Drop all unused seqlevels:
seqlevels(gr) <- seqlevelsInUse(gr)
}

\keyword{methods}
