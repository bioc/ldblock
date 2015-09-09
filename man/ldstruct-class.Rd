\name{ldstruct-class}
\Rdversion{1.1}
\docType{class}
\alias{ldstruct-class}
\alias{ldmat,ldstruct-method}
\alias{ldmat}

\title{Class \code{"ldstruct"}}
\description{
Manage information about LD statistics as reported by HapMap.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ldstruct", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ldmat}:}{Object of class \code{"dsCMatrix"} sparse
representation of statistics }
    \item{\code{chrom}:}{Object of class \code{"character"} chromosome
tag in UCSC format }
    \item{\code{genome}:}{Object of class \code{"character"} genome tag }
    \item{\code{allpos}:}{Object of class \code{"numeric"} coordinates }
    \item{\code{poptag}:}{Object of class \code{"character"} hapmap
founder population tag, 'CEU', 'MEX' etc. }
    \item{\code{statInUse}:}{Object of class \code{"character"} code
for statistic retrieved, one of 'Dprime', 'LOD', 'R2' }
    \item{\code{allrs}:}{Object of class \code{"character"} all
SNP identifiers, sometimes in affy format }
  }
}
\section{Methods}{
  \describe{
    \item{ldmat}{\code{signature(x = "ldstruct")}: extract
sparse matrix }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("ldstruct")
}
\keyword{classes}
