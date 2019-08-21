#' use LDmat API from NCI LDlink service
#' @importFrom httr GET content
#' @param rsvec character vector of SNP ids
#' @param pop three letter code for HapMap population, defaults to CEU
#' @param type 'r2' or 'd', defaults to 'd' implying d-prime
#' @return data.frame
#' @param token the API token provided by NCI, defaults to value of environment variable LDLINK_TOKEN
#' @examples 
#' if (interactive()) ldmat(c("rs77749396","rs9303279","rs9303280","rs9303281"))
#' @export
ldmat = function(rsvec, pop="CEU", type="d", token=Sys.getenv("LDLINK_TOKEN")) {
  snpstr = paste0(rsvec, collapse="%0A")
  reqstr = "https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?snps=%s&pop=%s&r2_d=%s&token=%s"
  req = sprintf(reqstr, snpstr, pop, type, token)
  ans = try(httr::GET(req))
  if (!is.null(httr::content(ans)$error)) stop(content(ans)$error)
  if (!inherits(ans, "try-error"))
     return(read.delim(textConnection(httr::content(ans, "text")), sep="\t"))
  else return(ans)
}


# POST METHOD HAS LARGER CAPACITY BUT DOES NOT SEEM FASTER 
#curl -k -H "Content-Type: application/json" -X POST -d '{"snps": "rs3\nrs4", "pop": "CEU","r2_d": "d"}' 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?&token=faketoken123'
#
