# using retrievals from
# hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/

#' container for LD data 
#' @importClassesFrom Matrix dsCMatrix
setClass("ldstruct", representation(ldmat = "dsCMatrix", chrom="character",
  genome="character", allpos="numeric", poptag="character", statInUse="character",
  allrs="character"))
setMethod("show", "ldstruct", function(object) {
  cat(paste0("ldstruct for population ", object@poptag, ", chrom ", object@chrom, "\n"))
  d = dim(object@ldmat)
  cat("dimensions:", d[1], "x", d[2], "; statistic is", object@statInUse, "\n")
  cat("object structure:\n")
  show(getSlots(class(object)))
})
setGeneric("ldmat", function(x)standardGeneric("ldmat"))

#' accessor for matrix component
#' @aliases ldmat,ldstruct-method
#' @param x instance of ldstruct
#' @export
setMethod("ldmat", "ldstruct", function(x) x@ldmat)



#' import hapmap LD data and create a structure for its management;
#' generates a sparse matrix representation of pairwise LD statistics and binds
#' metadata on variant name and position
#' @param hmgztxt name of gzipped text file as distributed at
#' \url{hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/}. It
#' will be processed by \code{\link{read.delim}}.
#' @param poptag heuristic tag identifying population
#' @param chrom heuristic tag for chromosome name
#' @param genome genome tag
#' @param stat statistic to use, "Dprime", "R2", and "LOD" are options
#' @return instance of ldstruct class
#' @keywords models
#' @examples
#' 
#' getClass("ldstruct")
#' # see vignette
#' 
#' @export hmld
hmld = function(hmgztxt, poptag, chrom, genome="hg19", stat="Dprime") {
 if (missing(hmgztxt)) {
    if (file.exists(fn <- paste0(Sys.getenv("LDBLOCK_TXTGZ_DIR"), "/",
             paste0("ld_", chrom, "_", poptag, ".txt.gz")))) hmgztxt = fn
    else stop("hmgztxt not supplied and no file in Sys.getenv('LDBLOCK_TXTGZ_DIR'); consider running downloadPopByChr()")
    }
 stopifnot(is.character(hmgztxt))
 stopifnot(length(grep("txt.gz$", hmgztxt))==1)
# requireNamespace("GenomicRanges")
 requireNamespace("Matrix")
 message(paste0("importing ", hmgztxt))
 lddf = read.delim(gzfile(hmgztxt), sep=" ", header=FALSE, stringsAsFactors=FALSE)
 message("done.")
 names(lddf) = c("pos1", "pos2", "pop", "rs1", "rs2", "Dprime",
   "R2", "LOD", "fbin")
 posvec = lddf$pos1
 names(posvec) = lddf$rs1
 extras = setdiff(lddf$rs2, lddf$rs1)
 if (length(extras)>0) {
   p2 = lddf$pos2
   names(p2) = lddf$rs2
   extrpos = p2[extras]
   posvec = c(posvec, extrpos)
   }
 urs = union(lddf$rs1, lddf$rs2)  # all rs numbers in play
 nur = length(urs)
 frs1 = factor(lddf$rs1, levels=urs)
 frs2 = factor(lddf$rs2, levels=urs)
 mm = Matrix::sparseMatrix(i=as.integer(frs1), j=as.integer(frs2), x=lddf[[stat]], 
    dims=c(nur,nur), dimnames=list(as.character(urs), as.character(urs)),
    symmetric=TRUE)
 diag(mm) = 1.0
 posvec = posvec[as.character(urs)]
 new("ldstruct", ldmat=mm, poptag=poptag, chrom=chrom, genome=genome,
    allrs = names(posvec), allpos=as.numeric(posvec), statInUse=stat)
}



#' download hapmap resource with LD estimates
#' 
#' download hapmap resource with LD estimates
#' 
#' delivers HapMap LD data to `targfolder`
#' 
#' @param chrname UCSC format tag for chromosome
#' @param popname hapmap three letter code for population, e.g. 'CEU'
#' @param urlTemplate pattern for creating URL given chr and pop
#' @param targfolder destination
#' @return just run for side effect of download.file
#' @keywords models
#' @examples
#' 
#' \dontrun{
#'  downloadPopByChr()
#'  }
#' 
#' @export downloadPopByChr
downloadPopByChr = function(chrname="chr1", popname="CEU",
   urlTemplate="http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/ld_%%CHRN%%_%%POPN%%.txt.gz", targfolder=Sys.getenv("LDBLOCK_TXTGZ_DIR")) {
   tmp = sub("%%CHRN%%", chrname, urlTemplate)
   tmp = sub("%%POPN%%", popname, tmp)
   download.file(tmp, paste0(targfolder, "/", basename(tmp)))
}



#' Given a set of SNP identifiers, use LD to expand the set to include linked
#' loci
#' 
#' Given a set of SNP identifiers, use LD to expand the set to include linked
#' loci
#' 
#' direct use of elementwise arithmetic comparison
#' 
#' @param rsl input list -- SNPs not found in the LD structure are simply
#' returned along with those found, and the expansion list, all combined in a
#' vector
#' @param lb lower bound on statistic used to retrieve loci in LD
#' @param ldstruct instance of \code{\link[ldblock]{ldstruct-class}}
#' @param chrn chromosome identifier
#' @param popn population identifier (one of 'CEU', 'MEX', ...)
#' @param txtgzfn path to gzipped hapmap file with LD information
#' @return character vector
#' @note As of 2015, it appears that locus names are more informative than
#' addresses for determining SNP identity across resources.
#' @keywords models
#' @examples
#' 
#'   og = Sys.getenv("LDBLOCK_TXTGZ_DIR")
#'   on.exit( Sys.setenv("LDBLOCK_TXTGZ_DIR" = og ) )
#'   Sys.setenv("LDBLOCK_TXTGZ_DIR"=system.file("hapmap", package="ldblock"))
#'   ld17 = hmld(chr="chr17", pop="CEU")
#'   ee = expandSnpSet( ld17@allrs[1:10], ldstruct = ld17 )
#' 
#' @export expandSnpSet
expandSnpSet = function(rsl, lb=.8, ldstruct, chrn="chr17", popn="CEU",
   txtgzfn = dir(system.file("hapmap", package="ldblock"), full.names=TRUE)) {
 if (missing(ldstruct)) curhm = hmld(txtgzfn, poptag=popn, chrom=chrn)@ldmat
 else curhm = ldstruct@ldmat
 bad = setdiff(rsl, rownames(curhm))
 if (length(bad)>0) warning(paste0("dropping ", length(bad), " rsn not matched in ld matrix"))
 rsl = setdiff(rsl, bad)
 hits = apply(curhm[rsl,], 1, function(x) which(x >= lb))
 sort(unique(c(unlist(hits), rsl, bad)))
}

