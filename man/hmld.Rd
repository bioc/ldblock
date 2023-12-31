% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmld.R
\name{hmld}
\alias{hmld}
\title{import hapmap LD data and create a structure for its management;
generates a sparse matrix representation of pairwise LD statistics and binds
metadata on variant name and position}
\usage{
hmld(hmgztxt, poptag, chrom, genome = "hg19", stat = "Dprime")
}
\arguments{
\item{hmgztxt}{name of gzipped text file as distributed at
\url{hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/}. It
will be processed by \code{\link{read.delim}}.}

\item{poptag}{heuristic tag identifying population}

\item{chrom}{heuristic tag for chromosome name}

\item{genome}{genome tag}

\item{stat}{statistic to use, "Dprime", "R2", and "LOD" are options}
}
\value{
instance of ldstruct class
}
\description{
import hapmap LD data and create a structure for its management;
generates a sparse matrix representation of pairwise LD statistics and binds
metadata on variant name and position
}
\examples{

getClass("ldstruct")
# see vignette

}
\keyword{models}
