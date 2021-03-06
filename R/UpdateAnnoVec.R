#' Function to update annotation of a character vector of gene aliases
#'
#' @param vec character vector of gene aliases
#' @import data.table
#' @importFrom AnnotationDbi select keys
#' @export
#'

# Main Method
updateAnno <- function(vec){
  vec_dt <- data.table(vec)
  setnames(vec_dt, "vec", "ALIAS")
  std <- suppressMessages(data.table(AnnotationDbi::select(org.Hs.eg.db,
                           keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"),
                           columns = c("ALIAS", "SYMBOL"),
                           keytype = "SYMBOL")))
  vec_dt[std, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  return(vec_dt$SYMBOL)
}
