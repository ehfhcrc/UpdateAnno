#' Function to update annotation of a character vector of gene aliases
#'
#' @param vec character vector of gene aliases
#' @import data.table
#' @export
#'

# Main Method
updateAnno <- function(vec){
  vec_dt <- data.table(vec)
  setnames(vec_dt, "vec", "ALIAS")
  std <- data.table(select(org.Hs.eg.db,
                           keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
                           columns = c("ALIAS", "SYMBOL"),
                           keytype = "SYMBOL"))
  vec_dt[std, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  return(vec_dt$SYMBOL)
}
