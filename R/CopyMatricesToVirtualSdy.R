copyMatricesToVirtualSdy <- function(baseUrl, virtualSdy){
  library(Rlabkey)
  runs <- labkey.selectRows(baseUrl = baseUrl,
                            folderPath = "/Studies/",
                            schemaName = "assay.ExpressionMatrix.matrix",
                            queryName = "Runs",
                            showHidden = TRUE)

  # subset to those in vSdy
  vSdyList <- list()
  vSdyList[["IS2"]] <- c("SDY56, SDY61, SDY63, SDY67, SDY80, SDY180, SDY212, SDY224, SDY269, SDY270, SDY400, SDY404, SDY520, SDY984, SDY1260, SDY1264, SDY1276, SDY1289, SDY1293")
  vSdyList[["IS1"]]  <- c("SDY63","SDY67","SDY80","SDY212","SDY400","SDY404")
  studies <- vSdyList[[virtualSdy]]
  studies <- strsplit(studies, ", ")[[1]]
  runs <- runs[ runs$Study %in% studies, ]

  # for each run cp over the normalized tsv only // 08.09.18 - will change and be unnecessary
  root <- "/share/files/Studies/"
  mid <- "/@files/analysis/exprs_matrices/"

  vPath <- paste0("/share/files/HIPC/", virtualSdy, "/@files/analysis/exprs_matrices/")
  for(i in 1:nrow(runs)){
    x <- runs[i, ]$Name
    print(x)

    mxPathOnSdy <- paste0(root, x$Study, mid, x$Name, ".tsv")
    mxPathOnVirt <- paste0(vPath, "Run", x$`Row Id`, "/", x$Name, ".tsv")

    # the latest / normalized .tsv are not same?
    file.copy(from = mxPathOnSdy,
              to = mxPathOnVirt,
              overwrite = TRUE)

    # if IS1
    if(virtualSdy == "IS1"){
      file.copy(from = paste0(mxPathOnSdy, ".immsig"),
                to = paste0(mxPathOnVirt, ".immsig"),
                overwrite = FALSE)
    }


    # file.copy(from = paste0(mxPathOnSdy, ".raw"),
    #           to = paste0(mxPathOnVirt, ".raw"),
    #           overwrite = FALSE)
    #
    # file.copy(from = paste0(mxPathOnSdy, ".summary"),
    #           to = paste0(mxPathOnVirt, ".summary"),
    #           overwrite = FALSE)
    #
    # file.copy(from = paste0(mxPathOnSdy, ".summary.orig"),
    #           to = paste0(mxPathOnVirt, ".summary.orig"),
    #           overwrite = FALSE)
  }
}

