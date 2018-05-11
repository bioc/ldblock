\name{hmld}
\alias{hmld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
import hapmap LD data and create a structure for its management
}
\description{
import hapmap LD data and create a structure for its management
}
\usage{
hmld(hmgztxt, poptag, chrom, genome = "hg19", stat = "Dprime")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmgztxt}{
name of gzipped text file as distributed
at 
\url{hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/}.
It will be processed by \code{\link{read.delim}}.
}
  \item{poptag}{
heuristic tag identifying population
}
  \item{chrom}{
heuristic tag for chromosome name
}
  \item{genome}{
genome tag
}
  \item{stat}{
statistic to use, "Dprime", "R2", and "LOD" are options
}
}
\details{
generates a sparse matrix representation of pairwise
LD statistics and binds metadata on variant name and
position
}
\value{
instance of ldstruct class
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
getClass("ldstruct")
# see vignette
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
