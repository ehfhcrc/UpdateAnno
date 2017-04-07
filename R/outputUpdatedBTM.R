#' Function to output a gmt file of the latest updated btm data from library
#'
#' @export
#'

# Main Method
outputUpdatedBTM <- function(){
  newbtm <- lapply(names(updated_btm_list), function(x){
    GeneSet(setName = x, geneIds = updated_btm_list[[x]])
  })
  gsc <- GeneSetCollection(newbtm)
  path2newgmt <- "./updated_btm.gmt"
  toGmt(gsc, path2newgmt)
}
