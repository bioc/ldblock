\name{ldByGene}
\alias{ldByGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
obtain LD statistics in region specified by a gene model
}
\description{
Obtain LD statistics in region specified by a gene model.
}
\usage{
ldByGene(sym = "MMP24", 
  vcf = system.file("vcf/c20exch.vcf.gz", package = "gQTLstats"), 
  flank = 1000, 
  vcfSLS = "NCBI", 
  genomeSLS = "hg19", 
  stats = "D.prime", depth = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym}{
A standard gene symbol for use with \code{\link{genemodel}}
}
  \item{vcf}{
Path to a tabix-indexed VCF file
}
  \item{flank}{
number of basepairs to flank gene model for search
}
  \item{vcfSLS}{
seqlevelsStyle (SLS) token for VCF; will be imposed on gene model
}
  \item{genomeSLS}{
character tag for genome, to be used with \code{\link[VariantAnnotation]{readVcf}}
}
  \item{stats}{
passed to \code{\link[snpStats]{ld}}
}
  \item{depth}{
passed to \code{\link[snpStats]{ld}} 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
sparse matrix representation of selected LD statistic,
as returned by \code{\link[snpStats]{ld}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (.Platform$OS.type == "unix") {
 ld1 = ldByGene(depth=150)
 image(ld1[1:200,1:200], col.reg=heat.colors(120), colorkey=TRUE,
   main="SNPs in MMP24 (chr20)") 
 }
}
\keyword{ models }
