#' Create a URL referencing 1000 genomes content in AWS S3.
#' stack1kg produces a VcfStack instance with references to VCF for 1000
#' genomes autosomal chrs.  S3-resident VCF files with version "v5a.20130502"
#' are used.
#' @importFrom methods getSlots new show
#' @importFrom utils download.file read.delim
#' @param chrnum a character string denoting a chromosome, such as '22'
#' @param tmpl alternate template for full URL, useful if versions prior to
#' 2010 are of interest
#' @param dropchr if TRUE \code{chrnum} will have 'chr' removed if present
#' @return by default, a TabixFile instance 
#' @note The "wrap" parameter has been removed.  A TabixFile structure will be returned.  The
#' tag parameter has been removed.  Supply a tmpl argument if you are not using 20130502 version.
#' @keywords models
#' @examples
#' requireNamespace("Rsamtools")
#' s3_1kg("22") # try scanVcfHeader from VariantAnnotation
#' @export
s3_1kg = function(chrnum, tmpl, dropchr=TRUE) {
  if (!requireNamespace("Rsamtools")) stop("install Rsamtools to use this function")
  if (dropchr) chrnum = gsub("chr", "", chrnum)
#
  tmpl = paste0("http://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr", chrnum, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
  Rsamtools::TabixFile(tmpl)
}

# March 2019 -- observed some issues getting data from AWS S3
# EBI images seem newer so will use them

#' couple together a group of VCFs
#' @importFrom BiocGenerics path
#' @param chrs a vector of chromosome names for extraction from 1000 genomes
#' VCF collection
#' @param index logical telling whether VcfStack should attempt to create the
#' local index; for 1000 genomes, the tbi are in the cloud and will be used by
#' readVcf so FALSE is appropriate
#' @param useEBI  logical(1) defaults to FALSE ... if TRUE, use tabix-indexed vcf from EBI,
#' but in July 2022 the EBI FTP site does not respond.  If FALSE, the AWS Open Data
#' access path is used
#' @return VcfStack instance
#' @note The seqinfo component of returned stack will have NA for genome.
#' Please set it manually; for useEBI=TRUE this would be GRCh38; very likely so
#' for useEBI=FALSE, but this should be checked.
#' @examples
#' if (interactive()) {
#'   st1 = stack1kg()
#'   st1
#'   }
#' @export
stack1kg = function(chrs=as.character(1:22), index=FALSE, useEBI=FALSE)
{
func = s3_1kg
if (!requireNamespace("GenomeInfoDb")) stop("please install GenomeInfoDb to use this package")
if (!requireNamespace("VariantAnnotation")) stop("please install VariantAnnotation to use this package")
if (!requireNamespace("Rsamtools")) stop("please install Rsamtools to use this package")
if (!requireNamespace("GenomicFiles")) stop("please install GenomicFiles to use this package")
if (useEBI) func = ebi_1kg
fs = tmp = sapply(chrs,function(x) BiocGenerics::path(func(x)))
names(tmp) = as.character(chrs)
tmp = GenomicFiles::VcfStack(tmp, seqinfo=GenomeInfoDb::Seqinfo(chrs), index=index)
fis = tmp@files
updf = lapply(fis, function(x) {Rsamtools::index(x) = paste0(path(x), ".tbi"); x})
tmp@files = VariantAnnotation::VcfFileList(updf)
tmp
}

ebi_1kg = function(chrnum, wrap, tmpl=NULL, dropchr=TRUE) { # = function(x) TabixFile(x),
 if (missing(wrap)) wrap = Rsamtools::TabixFile
  if (dropchr) chrnum = gsub("chr", "", chrnum)
template =  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr%s.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  wrap(sprintf(template, chrnum))
}

