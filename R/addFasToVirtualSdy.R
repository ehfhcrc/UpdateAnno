#' Function to migrate a copy of FA / FAS from /Studies/ to a virtual study container
#'
#' @param ISserver immunespace server, either test or prod
#' @param virtualSdy subdirectory within HIPC container
#' @param fasGrep grep statement to use in subsetting FAS to add just a few new sets
#' @import Rlabkey
#' @export
#'

# Main Method
addAnnoToVirtualSdy <- function(ISserver, virtualSdy, fasGrep = NULL){

    # Convert server to baseUrl
    baseUrl <- ifelse( ISserver == "prod",
                     "https://www.immunespace.org",
                     "https://test.immunespace.org")

    # Path assumed to be from HIPC directory
    vSdyPath <- paste0("/HIPC/", virtualSdy)

    # retrieve FAS
    fas <- labkey.selectRows(baseUrl = baseUrl,
                             folderPath = "/Studies/",
                             schemaName = "microarray",
                             queryName = "FeatureAnnotationSet",
                             colNameOpt = "fieldname",
                             showHidden = T)
    if (!is.null(fasGrep)) {
      fas <- fas[ grep(fasGrep, fas$Name), ]
    }
    fas <- fas[ order(fas$RowId), ]
    toImport <- data.frame(fas, stringsAsFactors = F)
    toImport[is.na(toImport)] <- ""
    toImport <- toImport[ grep("RowId", colnames(toImport), invert = T)]

    # Get container ID
    container <- labkey.selectRows(baseUrl = baseUrl,
                                    folderPath = vSdyPath,
                                    schemaName = "core",
                                    queryName = "containers",
                                    colNameOpt = "fieldname",
                                    colSelect = "EntityId")
    toImport$Container <- container[1,1]

    # Import FAS to vSdy
    impFas <- labkey.importRows(baseUrl = baseUrl,
                                folderPath = vSdyPath,
                                schemaName = "microarray",
                                queryName = "FeatureAnnotationSet",
                                toImport = toImport)

    # Get new rowId mappings
    newFas <- labkey.selectRows(baseUrl = baseUrl,
                                folderPath = vSdyPath,
                                schemaName = "microarray",
                                queryName = "FeatureAnnotationSet",
                                colNameOpt = "fieldname",
                                showHidden = T)
    if (!is.null(fasGrep)) {
      newFas <- newFas[ grep(fasGrep, newFas$Name), ]
    }
    newFas <- newFas[ order(match(newFas$Name, fas$Name)), ]

    if (!all(fas$Name == newFas$Name)) {
      stop("Imported FAS is not the same as /Studies/. Please fix.")
    }

    rowMap <- data.frame(old = fas$RowId,
                         new = newFas$RowId,
                         stringsAsFactors = F)

    # Update the featureAnnotation and Import
    # NOTE: colSelect "all" option aka "*" used in order to get FASid
    fasIdFilt <- makeFilter(c('FeatureAnnotationSetId', 'IN', paste(fas$RowId, collapse = ';')))
    fa <- labkey.selectRows(baseUrl = baseUrl,
                            folderPath = "/Studies/",
                            schemaName = "microarray",
                            queryName = "FeatureAnnotation",
                            colNameOpt = "fieldname",
                            colSelect = "*",
                            colFilter = fasIdFilt,
                            showHidden = T)
    newFa <- fa[ , colnames(fa) %in% c("Container","FeatureAnnotationSetId", "FeatureId", "GeneSymbol")]
    newFa$Container <- unique(toImport$Container)
    newFa$FeatureAnnotationSetId <- rowMap$new[ match(newFa$FeatureAnnotationSetId, rowMap$old)]
    newFa[is.na(newFa)] <- ""
    doneFa <- labkey.importRows(baseUrl = baseUrl,
                                folderPath = vSdyPath,
                                schemaName = "microarray",
                                queryName = "FeatureAnnotation",
                                toImport = newFa)
}
