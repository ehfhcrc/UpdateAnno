#' Function to migrate a copy of FA / FAS from /Studies/ to a virtual study container
#'
#' @param ISserver immunespace server, either test or prod
#' @param virtualSdy subdirectory within HIPC container
#' @import Rlabkey
#' @export
#'

# Main Method
addAnnoToVirtualSdy <- function(ISserver, virtualSdy){

    # Convert server to baseUrl
    if(ISserver == "test"){
      baseUrl <- "https://test.immunespace.org"
    }else if(ISserver == "prod"){
      baseUrl <- "https://www.immunespace.org"
    }else{
      stop("server name not recognized. Must be `test` or `prod`.")
    }

    # Path assumed to be from HIPC directory
    vSdyPath <- paste0("/HIPC/", virtualSdy)

    # retrieve FAS
    fas <- labkey.selectRows(baseUrl = baseUrl,
                             folderPath = "/Studies/",
                             schemaName = "microarray",
                             queryName = "FeatureAnnotationSet",
                             colNameOpt = "fieldname",
                             showHidden = T)
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
    fas <- fas[ order(fas$Name), ]
    newFas <- newFas[ order(newFas$Name), ]

    if(!all.equal(fas$Name, newFas$Name)){
      stop("Imported FAS is not the same as /Studies/. Please fix.")
    }

    rowMap <- data.frame(old = fas$RowId,
                         new = newFas$RowId,
                         stringsAsFactors = F)

    # Update the featureAnnotation and Import
    # NOTE: colSelect "all" option aka "*" used in order to get FASid
    fa <- labkey.selectRows(baseUrl = baseUrl,
                            folderPath = "/Studies/",
                            schemaName = "microarray",
                            queryName = "FeatureAnnotation",
                            colNameOpt = "fieldname",
                            colSelect = "*",
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
