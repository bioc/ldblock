\name{downloadPopByChr}
\alias{downloadPopByChr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
download hapmap resource with LD estimates
}
\description{
download hapmap resource with LD estimates
}
\usage{
downloadPopByChr(chrname = "chr1", 
 popname = "CEU", 
 urlTemplate = "http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/ld_\%\%CHRN\%\%_\%\%POPN\%\%.txt.gz", 
 targfolder = Sys.getenv("LDBLOCK_TXTGZ_DIR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrname}{
UCSC format tag for chromosome
}
  \item{popname}{
hapmap three letter code for population, e.g. 'CEU'
}
  \item{urlTemplate}{
pattern for creating URL given chr and pop
}
  \item{targfolder}{
destination
}
}
\details{
delivers HapMap LD data to `targfolder`
}
\value{
just run for side effect of download.file
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
 downloadPopByChr()
 }
}
\keyword{ models }
