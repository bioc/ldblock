#' obtain LD statistics in region specified by a gene model
#' 
#' Obtain LD statistics in region specified by a gene model.
#' 
#' 
#' @param sym A standard gene symbol for use with \code{\link[erma]{genemodel}}
#' @param vcf Path to a tabix-indexed VCF file
#' @param flank number of basepairs to flank gene model for search
#' @param vcfSLS seqlevelsStyle (SLS) token for VCF; will be imposed on gene
#' model
#' @param genomeSLS character tag for genome, to be used with
#' \code{\link[VariantAnnotation]{readVcf}}
#' @param stats passed to \code{\link[snpStats]{ld}}
#' @param depth passed to \code{\link[snpStats]{ld}}
#' @return sparse matrix representation of selected LD statistic, as returned
#' by \code{\link[snpStats]{ld}}
#' @keywords models
#' @examples
#' 
#' ld1 = ldByGene(depth=150)
#' image(ld1[1:200,1:200], col.reg=heat.colors(120), colorkey=TRUE,
#'  main="SNPs in MMP24 (chr20)") 
#' 
#' @export ldByGene
ldByGene = function(sym="MMP24", 
   vcf=system.file("vcf/c20exch.vcf.gz", package="gQTLstats"),
   flank=1000, vcfSLS = "NCBI", genomeSLS="hg19",
   stats = "D.prime", depth=10) {  # assumes Homo.sapiens
tf = TabixFile(vcf)
mod = range(genemodel(sym))
seqlevelsStyle(mod) = vcfSLS # to match vcf
toget = ScanVcfParam(which=(mod+flank)) # a little boundary added
r = readVcf(tf, param=toget, genome=genomeSLS)
gt = genotypeToSnpMatrix(r)[[1]]
ldout = ld(gt, stats=stats, depth=depth) # sparse matrix
ldout
}

