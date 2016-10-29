
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

