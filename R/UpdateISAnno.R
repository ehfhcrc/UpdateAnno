# Evan Henrich
# June 2017
# Purpose: To update the annotation throughout the ImmuneSpace website to the latest
# official gene symbols.

# HOW-TO:
# 1. ssh rsT or rsP
# 2. s (to become su)
# 3. R (to open up R)
# 4. devtools::install_github("ehfhcrc/UpdateAnno", build_vignette = T)
#    This ensures that the latest symbols are loaded into the data in the library by the vignette.
# 5. library(UpdateAnno)
# 5. runUpdateAnno()

# FOR TESTING ONLY:
# Dependencies --------------------------------------------------------
library(UpdateAnno)
library(data.table)
library(Rlabkey)
library(org.Hs.eg.db)
library(ImmuneSpaceR)

# NOTE: Main Script @ bottom of page

#######################################################################
###        1. Update Microarray.FeatureAnnotation table             ###
#######################################################################
# NOTES: This will take care of GeneExpressionExplore Module, which uses
# the Microarray.FeatureAnnotation table to populate a dropdown for selection
# of genes of interest.

# https://github.com/RGLab/LabKeyModules/blob/master/GeneExpressionExplorer/web/GeneExpressionExplorer.js#L233

updateFAS <- function(){
  currAnno <- data.table(labkey.selectRows(baseUrl = baseUrl,
                                           folderPath = "/Studies/",
                                           schemaName = "Microarray",
                                           queryName = "FeatureAnnotation",
                                           colNameOpt = "fieldname",
                                           showHidden = TRUE ))

  currAnno$GeneId <- updateAnno(currAnno$GeneSymbol) # uses org.Hs.eq.db under-the-hood

  # toUpdate <- ????

  # Ideally, push new gs back up to microarray.FeatureAnnotation so it's there!
  # updatedAnno <- labkey.updateRows(baseUrl = baseUrl,
  #                                  folderPath = "/Studies/",
  #                                  schemaName = "Microarray",
  #                                  queryName = "FeatureAnnotation",
  #                                  toUpdate = toUpdate)
}

#######################################################################
###            2. Update the expression matrices                    ###
#######################################################################
# allDirs <- list.dirs(path = "testEMs") # local on EH machine

updateEMs <- function(){
  # get list of study folders
  allDirs <- list.dirs(path = "/share/files/Studies/") # on rsT / rsP
  emDirs <- allDirs[ grep("*exprs_matrices", allDirs) ]

  # Get list of matrices
  runsDF <- labkey.selectRows(baseUrl = baseUrl,
                              folderPath = "/Studies/",
                              schemaName = "assay.ExpressionMatrix.matrix",
                              queryName = "Runs",
                              colNameOpt = "rname")

  # for each folder
  dmp <- lapply(emDirs, function(dir){
    fls <- list.files(dir)

    # get file basenames
    if( length(fls) > 0 ){
      tmp <- unique(unlist(strsplit(fls, split = ".tsv", fixed = TRUE)))
      baseNms <- tmp[ !(tmp %in% c(".summary",".orig",".summary.orig")) ]

      # go through each baseNm to update probe level and summary tsv
      sapply(baseNms, function(nm){
        baseFls <- fls[ grep(nm, fls) ]

        # Make tsv.orig if necessary (first time only)
        if( length(grep("*.orig", baseFls)) == 0){
          dmp <- sapply(baseFls, function(fl){ file.rename(file.path(dir,fl),
                                                           file.path(dir, paste0(fl, ".orig"))) })
        }

        # get probe-level original df and update annotation using only features for FASid
        # Reason for limiting to FASid and therefore needing executeSql instead of SelectRows
        # is that original FAS entries may have same probes mapped to different genes based
        # on changes in bioconductor libraries over time.  executeSql is only way to get FASid
        # because it is a lookup.
        prbEM <- fread(file.path(dir, paste0(nm, ".tsv.orig")))
        annoSetId <- runsDF$featureset[ runsDF$name == nm ]
        sqlStr <- sprintf("SELECT FeatureAnnotationSetId, FeatureId, GeneSymbol
                        from FeatureAnnotation
                        where FeatureAnnotationSetId='%s';", annoSetId)
        features <- labkey.executeSql(baseUrl = baseUrl,
                                      folderPath = "/Studies/",
                                      schemaName = "Microarray",
                                      sql = sqlStr,
                                      colNameOpt = "fieldname")
        prbEM <- prbEM[features, gene_symbol := GeneSymbol, on = c(V1 = "FeatureId")]

        # NOT taking NAs out here - does happen in summary process though
        write.table(prbEM, file = paste0(dir, "/", nm, ".tsv"), sep = "\t")

        # Summarize - lifted from Create-Matrix.R
        em <- prbEM[!is.na(gene_symbol) & gene_symbol != "NA"]
        sumEM <- em[ , lapply(.SD, mean), by="gene_symbol", .SDcols = grep("^BS", colnames(em)) ]

        write.table(sumEM, file = paste0(dir, "/", nm, ".tsv.summary"), sep = "\t")
      })
    }
  })
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

  # Any time that updateAnno is run, will need to blow away the table and redo via
  # https://github.com/RGLab/LabKeyModules/blob/master/DifferentialExpressionAnalysis/reports/schemas/assay.ExpressionMatrix.matrix/Runs/DifferentialExpressionAnalysis.Rmd
    # labkey.deleteRows()??

  # Not clear how the Rmd is run, but could loop over the code here to import new rows

# library(Rlabkey)
# library(ImmuneSpaceR)
# library(limma)
# library(Biobase)
# library(data.table)
# library(gtools)
#
# contrast <- c("study_time_collected", "study_time_collected_unit")
# study <- basename(labkey.url.path)
# i <- 1 #analysis accession key
# infostring <- ""
#
# Differential expression analysis
#
# Returns genes differentially expressed with an FDR of less than 20%, or top 100 lowest FDR.
#
# study: r study
# contrast: r contrast
#
# con <- CreateConnection(study)
# con$GeneExpressionInputs()
# getAllCoefs <- function(con){
#   runs <- con$data_cache$GE_matrices$name
#   if(is.null(runs)){
#     stop("This study does not have gene-expression data")
#   }
#   coefs <- con$data_cache$GE_inputs
#   coefs <- unique(coefs[, c("arm_name", contrast), with = FALSE])
#   coefs <- coefs[study_time_collected > 0]
#   return(coefs)
# }
#
#
# if(is.null(con$data_cache[["GE_matrices"]])){
#   nCoefs <- -1
#   infostring <- "There is no HIPCMatrix run in this study"
# } else{
#   runs <- con$data_cache$GE_matrices[, list(cohort, name)]
#   setnames(runs, "cohort", "arm_name")
#   allCoefs <- getAllCoefs(con)
#   allCoefs[, coefficient := do.call(paste, .SD), .SDcols = contrast]
#   existGEA <- data.table(labkey.selectRows(labkey.url.base, labkey.url.path, "gene_expression",
#                                            "gene_expression_analysis", colNameOpt = "rname"))
#   if(nrow(existGEA) > 0){
#     q1 <- quote(arm_name)
#     q2 <- quote(coefficient)
#     existGEA[, key := paste0(eval(q1), eval(q2))]
#     allCoefs <- allCoefs[, key := paste0(eval(q1), eval(q2))]
#     allCoefs <- allCoefs[!key %in% existGEA$key]
#     arm_todo <- unique(allCoefs$arm_name)
#     i <- max(as.numeric(gsub("^GEA", "", existGEA$analysis_accession))) + 1
#     runs <- runs[arm_name %in% arm_todo]
#   }
#   runs <- gsub(".tsv$", "", runs$name)
#   nCoefs <- nrow(allCoefs)
# }
#
# if(nCoefs == 0){
#   infostring <- "This analysis has already been run. You can visualize the results using the Data Explorer module."
# } else if(nCoefs > 0){
#   infostring <- paste("There will be", nCoefs, "new differential expression analysis.")
# }
#
# r infostring
#
# if(nCoefs == 0){
#   opts_chunk$set(eval = FALSE, echo = FALSE)
# }
#
# idx <- 1
# GEA_list <- vector("list")
# GEAR_list <- vector("list")
# #GEA_list <- vector("list", nCoefs)
# #GEAR_list <- vector("list", nCoefs)
# for(run in runs){
#   EM <- con$getGEMatrix(run)
#   pd <- data.table(pData(EM))
#   cm <- con$getDataset("cohort_membership")
#   cm <- unique(cm[, list(cohort, arm_accession)])
#   pd <- pd[, coef := do.call(paste, .SD), .SDcols = contrast]
#   to_drop <- unique(pd[study_time_collected <= 0, coef])
#   pd <- pd[coef %in% to_drop, coef := "baseline"]
#   pd <- pd[, coef := factor(coef,levels = c("baseline",
#                                             grep("baseline", value = TRUE, invert = TRUE, mixedsort(unique(coef)))))]
#   #mm <- model.matrix(formula("~ParticipantId + coef"), pd)
#   mm <- model.matrix(formula("~participant_id + coef"), pd)
#
#   # Check if it's RNA-seq or microarrays
#   if(max(exprs(EM))>100){
#     fit <- lmFit(voom(EM), mm)
#   } else{
#     fit <- lmFit(EM, mm)
#   }
#
#   fit <- eBayes(fit)
#
#   coefs <- grep("^coef", colnames(mm), value = TRUE)
#   for(coef in coefs){
#     analysis_accession <- paste0("GEA", i)
#     TP <- gsub("coef", "", coef)
#     arm_name <- unique(pData(EM)$cohort)
#     arm_accession <- cm[cohort == arm_name, arm_accession]
#     if(is.null(arm_name)){ arm_name <- NA }
#     description <- paste0("Differential expression in ", run, ", ", TP, " vs. baseline")
#
#     GEA_list[[idx]] <- data.table(analysis_accession = analysis_accession,
#                                   expression_matrix = run, arm_name = arm_name,
#                                   arm_accession = arm_accession,
#                                   coefficient = gsub("^coef", "", coef), description = description)
#     tt <- data.table(topTable(fit, coef = coef, number = Inf))
#     ttDE <- tt[adj.P.Val < 0.02]
#     if(nrow(ttDE) < 100){
#       tt <- tt[order(adj.P.Val)][1:min(nrow(tt), 100)]
#     } else{
#       tt <- ttDE
#     }
#     if(nrow(tt) > 0){
#       tt[, c("analysis_accession", "coefficient") := list(analysis_accession, coef)]
#       GEAR_list[[idx]] <- data.table(tt)
#     }
#     i <- i + 1
#     idx <- idx+1
#   }
# }
# GEA <- rbindlist(GEA_list)
# if(length(GEAR_list) == 0){
#   warning("No feature was found to be differentially expressed at any timepoint.")
#   opts_chunk$set(eval = FALSE, echo = FALSE)
# } else{
#   GEAR <- rbindlist(GEAR_list)
#   setnames(GEAR, c("FeatureId", "gene_symbol", "adj.P.Val", "AveExpr", "logFC", "P.Value", "t"), c("feature_id", "gene_symbol", "adj_p_val", "ave_expr", "log_fc", "p_value", "statistic"))
# }
#
# res_GEA <- labkey.importRows(labkey.url.base, labkey.url.path, "gene_expression",
#                              "gene_expression_analysis", toImport = GEA)
# res_GEAR <- labkey.importRows(labkey.url.base, labkey.url.path, "gene_expression",
#                               "gene_expression_analysis_results", toImport = GEAR)



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

# fast way to figure out PROD vs TEST
con <- CreateConnection("")
baseUrl <- con$config$labkey.url.base

runUpdateAnno <- function(){
  updateFAS()
  updateEMs()
  updateGEAR()
}


