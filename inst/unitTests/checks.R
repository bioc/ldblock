checkfe = function() {
  checkTrue(file.exists( system.file("hapmap/ld_chr17_CEU.txt.gz",
     package="ldblock") ) )
  }
checkfe()

checkImp = function() {
  Sys.setenv("LDBLOCK_TXTGZ_DIR"=system.file("hapmap", package="ldblock"))
  ld17 = hmld(chr="chr17", pop="CEU")
  checkTrue(is(ld17, "ldstruct"))
  checkTrue(is(ldmat(ld17), "dsCMatrix"))
  checkTrue(nrow(ldmat(ld17))==36621)
  checkTrue(ld17@statInUse == "Dprime")
}
checkImp()

checkExp = function() {
  Sys.setenv("LDBLOCK_TXTGZ_DIR"=system.file("hapmap", package="ldblock"))
  ld17 = hmld(chr="chr17", pop="CEU")
  ee = expandSnpSet( ld17@allrs[1:10], ldstruct = ld17 )
  checkTrue( length(ee) == 138 )
}
  
