\name{expandSnpSet}
\alias{expandSnpSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Given a set of SNP identifiers, use LD to expand the set to
include linked loci
}
\description{
Given a set of SNP identifiers, use LD to expand the set to
include linked loci
}
\usage{
expandSnpSet(rsl, lb = 0.8, ldstruct, chrn = "chr17", popn = "CEU", txtgzfn = dir(system.file("hapmap", package = "ldblock"), full = TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rsl}{
input list -- SNPs not found in the LD structure are simply
returned along with those found, and the expansion list, all
combined in a vector
}
  \item{lb}{
lower bound on statistic used to retrieve loci in LD
}
  \item{ldstruct}{
instance of \code{\link[ldblock]{ldstruct-class}}
}
  \item{chrn}{
chromosome identifier
}
  \item{popn}{
population identifier (one of 'CEU', 'MEX', ...)
}
  \item{txtgzfn}{
path to gzipped hapmap file with LD information
}
}
\details{
direct use of elementwise arithmetic comparison
}
\value{
character vector
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
As of 2015, it appears that locus names are more informative
than addresses for determining SNP identity across resources.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  og = Sys.getenv("LDBLOCK_TXTGZ_DIR")
  on.exit( Sys.setenv("LDBLOCK_TXTGZ_DIR" = og ) )
  Sys.setenv("LDBLOCK_TXTGZ_DIR"=system.file("hapmap", package="ldblock"))
  ld17 = hmld(chr="chr17", pop="CEU")
  ee = expandSnpSet( ld17@allrs[1:10], ldstruct = ld17 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
