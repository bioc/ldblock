\name{s3_1kg}
\alias{s3_1kg}
\alias{stack1kg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a URL referencing 1000 genomes content in AWS S3
}
\description{
Create a URL referencing 1000 genomes content in AWS S3.
}
\usage{
s3_1kg(chrnum, tag = "20130502", wrap = function(x) TabixFile(x), tmpl = NULL, dropchr = TRUE)
stack1kg(chrs=as.character(1:22))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrnum}{
a character string denoting a chromosome, such as '22'
}
  \item{chrs}{
a vector of chromosome names for extraction from 1000 genomes VCF collection
}
  \item{tag}{
a character string identifying the version, ignored if \code{tmpl}
is non-null; valid \code{tag} values are the default or "20101123"
}
  \item{wrap}{
The URL is returned after evaluating \code{wrap} on it; default is useful 
when Tabix indexing is to be used
}
  \item{tmpl}{
alternate template for full URL, useful if versions prior to 2010 are
of interest
}
  \item{dropchr}{
if TRUE \code{chrnum} will have 'chr' removed if present
}
  \item{index}{
a logical indicating if the vcf index files should be created
(for \code{stack1kg})
}
}
\details{
stack1kg produces a VcfStack instance with references to VCF
for 1000 genomes autosomal chrs.  S3-resident
VCF files with version "v5a.20130502" are used.
}
\value{
by default, a \code{\link[Rsamtools]{TabixFile}} instance
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
%%%  ~~who you are~~
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
s3_1kg("22")
\dontrun{
 require(VariantAnnotation)
 scanVcfHeader(s3_1kg("22"))
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
