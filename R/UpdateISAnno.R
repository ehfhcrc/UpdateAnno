# Evan Henrich
# June 2017
# Purpose: To update the annotation throughout the ImmuneSpace website to the latest
# official gene symbols.

# HOW-TO:
# 1. ssh rsT or rsP
# 2. become ImmuneSpace user (able to rename files AND access LK with auth)
# 3. R (to open up R)
# 4. devtools::install_github("RGLab/UpdateAnno", build_vignette = T)
#    This ensures that the latest symbols are loaded into the data in the library by the vignette.
# 5. library(UpdateAnno)
# 6. runUpdateAnno(baseUrl = "https://test.immunespace.org" )

# FOR TESTING ONLY:
# Dependencies --------------------------------------------------------
# library(UpdateAnno)
# library(data.table)
# library(Rlabkey)
# library(org.Hs.eg.db)
# library(ImmuneSpaceR)
# library(limma)
# library(Biobase)
# library(gtools)

# NOTE: Main Script @ bottom of page
#' @import Rlabkey
#' @import org.Hs.eg.db
#' @import ImmuneSpaceR
#' @import limma
#' @import Biobase
#' @import gtools

#######################################################################
###        1. Update Microarray.FeatureAnnotation table             ###
#######################################################################
# NOTES: This will take care of GeneExpressionExplorer Module, which uses
# the Microarray.FeatureAnnotation table to populate a dropdown for selection
# of genes of interest.  GEE will look for the FASid that was given to the original
# FAS when it was uploaded and this is challenging to change, therefore
# it was decided to move the original to a new FAS called "myOriginalFasName_orig"
# that gets a new FASid.  Con$getGEMatrix() then looks for this new FASid when
# populating the probe level gene symbols with the arg: `annotation = "default"`.

#' @export updateFAS
updateFAS <- function(baseUrl){

  # vars ---------------------------------------------------
  folderPath <- "/Studies/"
  schemaName <- "Microarray"

  # helper fn ----------------------------------------------
  getAnno <- function(nm, currFAS, baseUrl){
    annoSetId <- currFAS$RowId[ currFAS$Name == nm ]
    fasQuery <- sprintf("SELECT * from FeatureAnnotation
                          where FeatureAnnotationSetId='%s';",
                        annoSetId)
    features <- labkey.executeSql(baseUrl = baseUrl,
                                  folderPath = folderPath,
                                  schemaName = schemaName,
                                  sql = fasQuery,
                                  colNameOpt = "fieldname",
                                  showHidden = T)
    toDrop <- c("Created", "CreatedBy", "Modified", "ModifiedBy")
    features <- features[ , !(colnames(features) %in% toDrop) ]
  }

  currFas <- function(baseUrl){
    res <- data.table(labkey.selectRows(baseUrl = baseUrl,
                                        folderPath = folderPath,
                                        schemaName = schemaName,
                                        queryName = "FeatureAnnotationSet",
                                        colNameOpt = "fieldname",
                                        showHidden = TRUE ))
  }


  # MAIN ------------------------------------------------------------
  # for each name in fas$name see if there is an updated version then work on the
  # updated version if there is one or create one if there is not. Exceptions are those
  # with "Do Not Update", "non-updatable" or similar string in the Comment column.
  currFAS <- currFas(baseUrl)
  fasNms <- currFAS$Name[ grep("[N|n]o[n|t|][-| ][u|U]pdat[e|able]", currFAS$Comment, invert = T) ]

  lapply(fasNms, FUN = function(nm){
    message(paste0("Updating: ", nm))
    currAnno <- getAnno(nm, currFAS, baseUrl)
    orNm <- paste0(nm, "_orig")
    if( orNm %in% currFAS$Name ){
      message("Updating FAS")
      # if orig is present means that update has been performed at least once.
      # Update the rows of the previously updated anno using the orig
      # as the base for mapping to ensure that updates of updates are avoided
      orAnno <- getAnno(orNm, currFAS, baseUrl)
      orAnno$GeneSymbol <- updateAnno(orAnno$GeneSymbol)
      currAnno$GeneSymbol <- orAnno$GeneSymbol[ match(currAnno$FeatureId, orAnno$FeatureId) ]
      currAnno[ is.na(currAnno) ] <- ""
      toUpdate <- data.frame(currAnno, stringsAsFactors = F)
      done <- labkey.updateRows(baseUrl = baseUrl,
                                folderPath = folderPath,
                                schemaName = schemaName,
                                queryName = "FeatureAnnotation",
                                toUpdate = toUpdate)
    }else{
      print("Creating New updated FAS")
      # Create featureAnnotationSet with "_orig" name
      toImport <- data.frame(currFAS[ currFAS$Name == nm, ])
      toImport <- toImport[ , !(colnames(toImport) == "RowId") ]
      toImport$Name <- orNm
      addFas <- labkey.importRows(baseUrl = baseUrl,
                                  folderPath = folderPath,
                                  schemaName = schemaName,
                                  queryName = "FeatureAnnotationSet",
                                  toImport = toImport)

      # check that new "_orig" set has been imported correctly
      nowFas <- currFas(baseUrl)
      if( !(toImport$Name[[1]] %in% nowFas$Name) ){
        stop("Original FAS (",toImport$Name[[1]],") not imported correctly")
      }

      # Prep for importRows by removing rowIds (will be given) and using new FASid
      orAnno <- currAnno[ , !(colnames(currAnno) == "RowId") ]
      orAnno$FeatureAnnotationSetId <- nowFas$RowId[ nowFas$Name == toImport$Name[[1]] ]
      orAnno[ is.na(orAnno) ] <- ""
      toImport <- data.frame(orAnno, stringsAsFactors = F)
      addFeatures <- labkey.importRows(baseUrl = baseUrl,
                                       folderPath = folderPath,
                                       schemaName = schemaName,
                                       queryName = "FeatureAnnotation",
                                       toImport = toImport)

      featureChk <- labkey.selectRows(baseUrl = baseUrl,
                                      folderPath = folderPath,
                                      schemaName = schemaName,
                                      queryName = "FeatureAnnotation",
                                      colFilter = makeFilter(c("FeatureAnnotationSetId",
                                                               "EQUALS",
                                                               unique(toImport$FeatureAnnotationSetId)))
                                      )
      if( all.equal(featureChk, toImport) ){
        # Now update the old fasId rows with new geneSymbols
        currAnno$GeneSymbol <- updateAnno(currAnno$GeneSymbol)
        currAnno[ is.na(currAnno) ] <- ""
        toUpdate <- data.frame(currAnno, stringsAsFactors = F)
        done <- labkey.updateRows(baseUrl = baseUrl,
                                  folderPath = folderPath,
                                  schemaName = schemaName,
                                  queryName = "FeatureAnnotation",
                                  toUpdate = toUpdate)
      }else{
        stop("Original FA not uploaded correctly to *_orig table")
      }


    }
  })
  return(TRUE)
}

#######################################################################
###            2. Update the expression matrices                    ###
#######################################################################
# if running by hand on rsT/rsP
# library(ImmuneSpaceR)
# library(Rlabkey)
# library(data.table)


#' @export updateEMs
updateEMs <- function(sdy, runsDF){
  print(paste0("working on study: ", sdy))
  # get file basenames
  dirPath <- file.path("/share/files/Studies",
                       sdy,
                       "@files/analysis/exprs_matrices")
  fls <- list.files(dirPath)
  tmp <- unique(unlist(strsplit(fls, split = ".tsv", fixed = TRUE)))
  baseNms <- tmp[ !(tmp %in% c(".summary",".orig",".summary.orig")) ]

  # go through each baseNm to update summary tsv
  sapply(baseNms, function(nm){
    print(paste0("working on baseNm: ", nm))
    baseFls <- fls[ grep(nm, fls) ]

    # Rename original summary file to tsv.summary.orig if necessary (first time only)
    if( !(paste0(nm, ".tsv.summary.orig") %in% baseFls) ){
      sumFl <- paste0(dirPath, "/", nm, ".tsv.summary")
      dmp <- file.rename( sumFl, paste0(sumFl, ".orig") )
    }

    # get probe-level original df and update annotation using only features for FASid
    # Reason for limiting to FASid and therefore needing executeSql instead of SelectRows
    # is that original FAS entries may have same probes mapped to different genes based
    # on changes in bioconductor libraries over time.  ExecuteSql is only way to get FASid
    # because it is a lookup even though it is in microarray.FeatureAnnotation.
    prbEM <- fread(file.path(dirPath, paste0(nm, ".tsv")))
    annoSetId <- runsDF$featureset[ runsDF$name == nm ]
    sqlStr <- sprintf("SELECT FeatureAnnotationSetId, FeatureId, GeneSymbol
                    from FeatureAnnotation
                    where FeatureAnnotationSetId='%s';", annoSetId)
    features <- labkey.executeSql(baseUrl = baseUrl,
                                  folderPath = "/Studies/",
                                  schemaName = "Microarray",
                                  sql = sqlStr,
                                  colNameOpt = "fieldname")
    prbEM[ , V1 := as.character(V1)] # for SDY80 where probes have integer vals
    prbEM <- prbEM[features, gene_symbol := GeneSymbol, on = c(V1 = "FeatureId")]

    # Summarize - lifted from Create-Matrix.R
    em <- prbEM[ !is.na(gene_symbol) & gene_symbol != "NA" ]
    sumEM <- em[ , lapply(.SD, mean), by="gene_symbol", .SDcols = grep("^BS", colnames(em)) ]

    write.table(sumEM, file = paste0(dirPath, "/", nm, ".tsv.summary"), sep = "\t")
  })
  return(TRUE)
}

#######################################################################
###            3. Update DGEA-Results                               ###
#######################################################################
# This table is used in the module Differential Gene Expression Analysis
# GEARres <- labkey.selectRows(baseUrl = baseUrl,
#                              folderPath = "/Studies/",
#                              schemaName = "gene_expression",
#                              queryName = "gene_expression_analysis_results",
#                              colNameOpt = "fieldname",
#                              showHidden = TRUE)

# This method is based on the DGEA.Rmd found in the DGEA module

#' @export updateGEAR
updateGEAR <- function(sdy, baseUrl, runsDF){
  print(paste0("working on study: ", sdy))
  infostring <- ""
  labkey.url.base <- baseUrl
  contrast <- c("study_time_collected", "study_time_collected_unit")
  con <- CreateConnection(sdy)
  con$GeneExpressionInputs()
  GEA_list <- vector("list")
  GEAR_list <- vector("list")

  runs <- runsDF$name[ runsDF$folder_name == sdy ]

  idx <- 1 # analysis accession key
  for(run in runs){
    print(paste0("working on run: ", run))
    EM <- con$getGEMatrix(run, outputType = "normalized", annotation = "latest") # note params!
    pd <- data.table(pData(EM))
    pd <- pd[, coef := do.call(paste, .SD), .SDcols = contrast]
    to_drop <- unique(pd[study_time_collected <= 0, coef])
    pd <- pd[coef %in% to_drop, coef := "baseline"]
    pd <- pd[, coef := factor(coef, levels = c("baseline",
                                               grep("baseline",
                                                    value = TRUE,
                                                    invert = TRUE,
                                                    mixedsort(unique(coef)))))]

    # pd$coefs are timepoints so if only 1 then can't do differential
    if( length(unique(pd$coef)) > 1 ){
      mm <- model.matrix(formula("~participant_id + coef"), pd)

      # Check if it's RNA-seq or microarrays
      if( max(exprs(EM)) > 100 ){ EM <- voom(EM) }
      fit <- lmFit(EM, mm)
      fit <- eBayes(fit)

      # Prep for coefficients work
      cm <- con$getDataset("cohort_membership")
      cm <- unique( cm[, list(cohort, arm_accession)] )
      coefs <- grep("^coef", colnames(mm), value = TRUE)

      for(coef in coefs){
        analysis_accession <- paste0("GEA", idx)
        TP <- gsub("coef", "", coef)
        arm_name <- unique(pData(EM)$cohort)
        arm_accession <- cm[cohort == arm_name, arm_accession]
        arm_name[ is.null(arm_name) ] <- NA
        description <- paste0("Differential expression in ", run, ", ", TP, " vs. baseline")

        GEA_list[[idx]] <- data.table(analysis_accession = analysis_accession,
                                      expression_matrix = run,
                                      arm_name = arm_name,
                                      arm_accession = arm_accession,
                                      coefficient = gsub("^coef", "", coef),
                                      description = description)

        tt <- data.table(topTable(fit, coef = coef, number = Inf))

        tt <- if( sum(tt$adj.P.Val < 0.02) < 100 ){
                tt[order(adj.P.Val)][1:min(nrow(tt), 100)]
              }else{
                tt[adj.P.Val < 0.02]
              }

        if(nrow(tt) > 0){
          tt[, c("analysis_accession", "coefficient") := list(analysis_accession, coef)]
          tt[, coefficient := gsub("coef","", coefficient) ]
          GEAR_list[[idx]] <- data.table(tt)
        }
      }
      idx <- idx + 1
    }else{
      print("lengths(coefs) < 2 therefore no differential")
    }
  }

  if( length(GEA_list) > 0 ){
    # Set HTTP call vars
    folderPath <- con$config$labkey.url.path
    schemaGE <- "gene_expression"
    queryGEA <- "gene_expression_analysis"
    queryRes <- "gene_expression_analysis_results"

    # delete old GEA
    currGEA <- labkey.selectRows(baseUrl = baseUrl,
                                 folderPath = folderPath,
                                 schemaName = schemaGE,
                                 queryName = queryGEA,
                                 colNameOpt = "rname",
                                 showHidden = T)
    if( nrow(currGEA) != 0 ){
      deleteGEA <- labkey.deleteRows(baseUrl = baseUrl,
                                     folderPath = folderPath,
                                     schemaName = schemaGE,
                                     queryName = queryGEA,
                                     toDelete = currGEA)
      if( deleteGEA$rowsAffected != nrow(currGEA) ){
        stop("currGEA not deleted correctly")
      }
    }

    # push newGEA b/c listings may be different in terms of idx than old
    newGEA <- rbindlist(GEA_list)
    doneGEA <- labkey.importRows(baseUrl = baseUrl,
                                 folderPath = folderPath,
                                 schemaName = schemaGE,
                                 queryName = queryGEA,
                                 toImport = newGEA)

    if( doneGEA$rowsAffected != nrow(newGEA) ){
      stop("newGEA not imported correctly")
    }

    # GEAR gets deleted and then new rows imported because
    # new mappings will be different and do not want to have leftovers
    if(length(GEAR_list) != 0){
      currGEAR <- labkey.selectRows(baseUrl = baseUrl,
                                    folderPath = folderPath,
                                    schemaName = schemaGE,
                                    queryName = queryRes,
                                    colNameOpt = "rname",
                                    showHidden = T)

      if( nrow(currGEAR) != 0 ){
        delGEAR <- labkey.deleteRows(baseUrl = baseUrl,
                                     folderPath = folderPath,
                                     schemaName = schemaGE,
                                     queryName = queryRes,
                                     toDelete = currGEAR)

        postDeleteGEAR <- labkey.selectRows(baseUrl = baseUrl,
                                            folderPath = folderPath,
                                            schemaName = schemaGE,
                                            queryName = queryRes,
                                            colNameOpt = "rname",
                                            showHidden = T)

        if( nrow(postDeleteGEAR) != 0){
          stop("not all GEAR deleted correctly")
        }
      }

      # Import new GEAR
      GEAR <- rbindlist(GEAR_list)
      GEAR[ is.na(GEAR) ] <- ""
      setnames(GEAR,
               c("FeatureId", "gene_symbol", "adj.P.Val", "AveExpr", "logFC", "P.Value", "t"),
               c("feature_id", "gene_symbol", "adj_p_val", "ave_expr", "log_fc", "p_value", "statistic"))
      toImport <- data.frame(GEAR, stringsAsFactors = F)
      resGEAR <- labkey.importRows(baseUrl = baseUrl,
                                   folderPath = folderPath,
                                   schemaName = schemaGE,
                                   queryName = queryRes,
                                   toImport = toImport)
    }
  }
  return(TRUE)
}

#######################################################################
###            4. Update GSEA Modules                               ###
#######################################################################
# gene set (gs) module lists data was originally in library(ImmuneSpaceData)
# The gs are now in library(UpdateAnno) and are loaded in the GSEA Rmd.
# Following is the code snippet that shows how this is done:

# switch(set,
#        `MSigDB c7` = {data(msigdb_immunologic_signatures);
#          gene_sets <- msigdb_immunologic_signatures;
#          func_rnames <- msig_rnames},
#        `Blood transcription` = {data(emory_blood_transcript_modules);
#          gene_sets <- emory_blood_transcript_modules;
#          func_rnames <- emory_rnames},
#        `G2 (Trial 8) Modules` = {data(chaussabel_modules);
#          gene_sets <- chaussabel_modules;
#          func_rnames <- baylor_rnames})

#######################################################################
###                         MAIN SCRIPT                             ###
#######################################################################

#' @export runUpdateAnno
runUpdateAnno <- function(baseUrl){
  folderpath <- "/Studies/"

  # Double-check whether to run on test or prod
  chk <- readline(prompt = paste0("You are running on ",
                                  baseUrl,
                                  ". Continue? [T/f] "))

  if( !(chk %in% c("", "t", "T")) ){ return("Operations ended.") }

  # Update the secondary / updated FeatureAnnotation set
  updateFAS(baseUrl)

  # Get studies with gene expression matrices
  runsDF <- labkey.selectRows(baseUrl = baseUrl,
                              folderPath = folderPath,
                              schemaName = "assay.ExpressionMatrix.matrix",
                              queryName = "Runs",
                              colNameOpt = "rname")
  sdys <- unique(runsDF$folder_name)
  lapply(sdys, updateEMs, runsDF = runsDF) # update flat files for summary only
  lapply(sdys, updateGEAR, baseUrl = baseUrl, runsDF = runsDF) # using con$getGEMatrix(outputType = "normalized", annotation = "latest")
}


