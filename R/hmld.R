# using retrievals from
# hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/

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
setMethod("ldmat", "ldstruct", function(x) x@ldmat)

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
 mm = sparseMatrix(i=as.integer(frs1), j=as.integer(frs2), x=lddf[[stat]], 
    dims=c(nur,nur), dimnames=list(as.character(urs), as.character(urs)),
    symmetric=TRUE)
 diag(mm) = 1.0
 posvec = posvec[as.character(urs)]
 new("ldstruct", ldmat=mm, poptag=poptag, chrom=chrom, genome=genome,
    allrs = names(posvec), allpos=as.numeric(posvec), statInUse=stat)
}

downloadPopByChr = function(chrname="chr1", popname="CEU",
   urlTemplate="http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/2009-02_phaseIII_r2/ld_%%CHRN%%_%%POPN%%.txt.gz", targfolder=Sys.getenv("LDBLOCK_TXTGZ_DIR")) {
   tmp = sub("%%CHRN%%", chrname, urlTemplate)
   tmp = sub("%%POPN%%", popname, tmp)
   download.file(tmp, paste0(targfolder, "/", basename(tmp)))
}

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

