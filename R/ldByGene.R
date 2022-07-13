#' Obtain LD statistics in region specified by a gene model.
#' @param sym A standard gene symbol for use with \code{genemodel}
#' @param vcf Path to a tabix-indexed VCF file
#' @param flank number of basepairs to flank gene model for search
#' @param vcfSLS seqlevelsStyle (SLS) token for VCF; will be imposed on gene
#' model
#' @param genomeSLS character tag for genome, to be used with
#' \code{readVcf}
#' @param stats passed to \code{\link[snpStats]{ld}}
#' @param depth passed to \code{\link[snpStats]{ld}}
#' @return sparse matrix representation of selected LD statistic, as returned
#' by \code{\link[snpStats]{ld}}
#' @note Uses an internal function genemod4ldblock, that relies
#' on EnsDb.Hsapiens.v75 to get gene model.
#' @keywords models
#' @examples
#' if (interactive()) {  # there is a warning owing to non-SNV present
#' ld1 = ldByGene(depth=150)
#' image(ld1[1:200,1:200], col.reg=heat.colors(120), colorkey=TRUE,
#'  main="SNPs in MMP24 (chr20)") 
#' }
#' @export ldByGene
ldByGene = function(sym="MMP24", 
   vcf=system.file("vcf/c20exch.vcf.gz", package="ldblock"),
   flank=1000, vcfSLS = "NCBI", genomeSLS="hg19",
   stats = "D.prime", depth=10) {  # assumes Homo.sapiens
if (!requireNamespace("snpStats")) stop("install snpStats to use this function")
if (!requireNamespace("GenomeInfoDb")) stop("install GenomeInfoDb to use this function")
if (!requireNamespace("VariantAnnotation")) stop("install VariantAnnotation to use this function")
if (!requireNamespace("Rsamtools")) stop("install Rsamtools to use this function")
tf = Rsamtools::TabixFile(vcf)
mod = range(genemodel4ldblock(sym))
GenomeInfoDb::seqlevelsStyle(mod) = vcfSLS # to match vcf
toget = VariantAnnotation::ScanVcfParam(which=(mod+flank)) # a little boundary added
r = VariantAnnotation::readVcf(tf, param=toget, genome=genomeSLS)
gt = VariantAnnotation::genotypeToSnpMatrix(r)[[1]]
ldout = snpStats::ld(gt, stats=stats, depth=depth) # sparse matrix
ldout
}

genemodel4ldblock = function(sym,
   annoresource) {
 if (!requireNamespace("ensembldb")) stop("install ensembldb to use this function")
 if (missing(annoresource) & requireNamespace("EnsDb.Hsapiens.v75")) 
       annoresource = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
 else if (missing(annoresource)) stop("supply annoresource or install EnsDb.Hsapiens.v75 to use this function")
 gg = ensembldb::transcripts(annoresource, columns="gene_name")
 ind = grep(sym, gg$gene_name)
 if (length(ind)==0) stop("cannot resolve symbol in EnsDb.Hsapiens.v75")
 gg[ind]
}
