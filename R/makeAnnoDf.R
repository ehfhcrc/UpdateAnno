#' Function to map probe ids to gene symbols from an annotation package
#'
#' @param annoPkg name of annotation package to use
#' @param outPath file path for tsv to write out for upload to ImmuneSpace
#' @export
#'
makeAnnoDf <- function(annoPkg, outPath = NULL){
  library(annoPkg, character.only = TRUE)
  prbLs <- paste0(gsub("\\.db", "", annoPkg), "ENTREZID")
  probes <- ls(eval(parse(text = prbLs)))
  symLs <- gsub("ENTREZID", "SYMBOL", prbLs)
  syms <- unlist(mget(probes, eval(parse(text = symLs))))
  syms <- syms[ !is.na(syms)] # No NAs or "" allowed in IS FAS
  res <- data.frame(Probe_ID = names(syms),
                    Gene_Symbol = syms,
                    stringsAsFactors = F)
  if( !is.null(outPath) ){
    write.table(res,
                file = outPath,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
  }
  return(res)
}
