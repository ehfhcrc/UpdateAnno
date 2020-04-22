#' Function to update annotation of a character vector of gene aliases
#'
#' @param aliases character vector of gene aliases
#' @import data.table
#' @export
#'
mapAlias2Symbol <- function(aliases){
  vec_dt <- data.table(aliases)
  setnames(vec_dt, "aliases", "ALIAS")
  vec_dt[hgncAlias2Symbol, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  return(vec_dt$SYMBOL)
}
