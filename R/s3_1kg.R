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

stack1kg = function(chrs=as.character(1:22), index=FALSE)
{
fs = tmp = sapply(chrs,function(x) path(s3_1kg(x)))
names(tmp) = as.character(chrs)
tmp = VcfStack(tmp, seqinfo=Seqinfo(chrs), index=index)
fis = tmp@files
updf = lapply(fis, function(x) {index(x) = paste0(path(x), ".tbi"); x})
tmp@files = VcfFileList(updf)
t1 = TabixFile(fs[1])
seqinfo(tmp) = seqinfo(scanVcfHeader(t1))
tmp
}

