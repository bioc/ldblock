#' Create a URL referencing 1000 genomes content in AWS S3.
#' 
#' stack1kg produces a VcfStack instance with references to VCF for 1000
#' genomes autosomal chrs.  S3-resident VCF files with version "v5a.20130502"
#' are used.
#' 
#' @param chrnum a character string denoting a chromosome, such as '22'
#' @param tag a character string identifying the version, ignored if
#' \code{tmpl} is non-null; valid \code{tag} values are the default or
#' "20101123"
#' @param wrap The URL is returned after evaluating \code{wrap} on it; default
#' is useful when Tabix indexing is to be used
#' @param tmpl alternate template for full URL, useful if versions prior to
#' 2010 are of interest
#' @param dropchr if TRUE \code{chrnum} will have 'chr' removed if present
#' @return by default, a \code{\link{TabixFile}} instance %%
#' @keywords models
#' @examples
#' s3_1kg("22")
#' \dontrun{
#'  require(VariantAnnotation)
#'  scanVcfHeader(s3_1kg("22"))
#'  }
#' @export
s3_1kg = function(chrnum, tag="20130502", wrap = function(x) TabixFile(x),
 tmpl=NULL, dropchr=TRUE) {

  if (dropchr) chrnum = gsub("chr", "", chrnum)
  tmpl2010 = "http://1000genomes.s3.amazonaws.com/release/20110521/ALL.chr%%NUM%%.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

  tmpl2013 = "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr%%NUM%%.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

  if (is.null(tmpl)) {
   if (tag == "20130502") tmpl = tmpl2013
   else if (tag == "20101123") tmpl = tmpl2010
   else stop("must supply valid tag or tmpl, see ?s3_1kg in ldblock package")
   }

  wrap(gsub("%%NUM%%", chrnum, tmpl))
}

# March 2019 -- observed some issues getting data from AWS S3
# EBI images seem newer so will use them

#' couple together a group of VCFs
#' @import GenomicFiles
#' @importFrom BiocGenerics path
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom Rsamtools "index<-"
#' @importFrom VariantAnnotation VcfFileList
#' @param chrs a vector of chromosome names for extraction from 1000 genomes
#' VCF collection
#' @param index logical telling whether VcfStack should attempt to create the
#' local index; for 1000 genomes, the tbi are in the cloud and will be used by
#' readVcf so FALSE is appropriate
#' @param useEBI  logical(1) defaults to TRUE ... use tabix-indexed vcf from EBI
#' @return VcfStack instance
#' @note The seqinfo component of returned stack will have NA for genome.
#' Please set it manually; for useEBI=TRUE this would be GRCh38.
#' @examples
#' if (interactive()) {
#'   st1 = stack1kg()
#'   st1
#'   }
#' @export
stack1kg = function(chrs=as.character(1:22), index=FALSE, useEBI=TRUE)
{
func = s3_1kg
if (useEBI) func = ebi_1kg
fs = tmp = sapply(chrs,function(x) path(func(x)))
names(tmp) = as.character(chrs)
tmp = VcfStack(tmp, seqinfo=Seqinfo(chrs), index=index)
fis = tmp@files
updf = lapply(fis, function(x) {index(x) = paste0(path(x), ".tbi"); x})
tmp@files = VcfFileList(updf)
tmp
}

ebi_1kg = function(chrnum, wrap = function(x) TabixFile(x),
 tmpl=NULL, dropchr=TRUE) {
  if (dropchr) chrnum = gsub("chr", "", chrnum)
template =  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr%s.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  wrap(sprintf(template, chrnum))
}

