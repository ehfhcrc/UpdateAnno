#' Original HIPC BTM list
#'
#' A dataset containing the geneSet pathways and geneIds for
#' a BTM file with the original date as October 08, 2013 generated
#' by the HIPC collaborators and stored in the HIPC google drive.
#'
#' @format A list of lists with charactar vector elements
#' \describe{
#'   \item{names}{geneSet pathway names}
#'   \item{elements}{character vector of gene symbols}
#'   ...
#' }
"orig_btm_list"

#' Updated HIPC BTM list
#'
#' A dataset containing the geneSet pathways and geneIds for
#' a BTM file that is updated using the most recent version of
#' the org.Hs.eq.db package each time the UpdateAnno package is
#' re-installed with build_vignette = TRUE.
#'
#' @format A list of lists with charactar vector elements
#' \describe{
#'   \item{names}{geneSet pathway names}
#'   \item{elements}{character vector of gene symbols}
#'   ...
#' }
"updated_btm_list"

#' Updated HIPC BTM data frame
#'
#' A dataset containing the geneSet pathways and geneIds for
#' a BTM file with the original date as October 08, 2013 generated
#' by the HIPC collaborators and stored in the HIPC google drive.
#'
#' @format A dataframe with pathways and geneIds
#' \describe{
#'   \item{pathway}{geneSet pathway names}
#'   \item{SYMBOL}{updated gene symbols}
#'   ...
#' }
"updated_btm_df"
