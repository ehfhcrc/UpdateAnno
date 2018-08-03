library(Rlabkey)

baseUrl <- "https://test.immunespace.org"
# get runs from /Studies
runs <- labkey.selectRows(baseUrl = baseUrl,
                          folderPath = "/Studies/",
                          schemaName = "assay.ExpressionMatrix.matrix",
                          queryName = "Runs",
                          showHidden = TRUE)

# subset to those in IS2
studies <- c("SDY56, SDY61, SDY63, SDY67, SDY80, SDY180, SDY212, SDY224, SDY269, SDY270, SDY400, SDY404, SDY520, SDY984, SDY1260, SDY1264, SDY1276, SDY1289, SDY1293")
studies <- strsplit(studies, ", ")[[1]]
runs <- runs[ runs$Study %in% studies, ]

# for each run cp over the .tsv.raw, .tsv.summary, .tsv.summary.orig files
root <- "/share/files/Studies/"
mid <- "/@files/analysis/exprs_matrices/"

vPath <- "/share/files/HIPC/IS2/@files/analysis/exprs_matrices/"
for(i in 1:nrow(runs)){
  x <- runs[i, ]

  mxPathOnSdy <- paste0(root, x$Study, mid, x$Name, ".tsv")
  mxPathOnVirt <- paste0(vPath, "Run", x$`Row Id`, "/", x$Name, ".tsv")

  file.copy(from = paste0(mxPathOnSdy, ".raw"),
            to = paste0(mxPathOnVirt, ".raw"))

  file.copy(from = paste0(mxPathOnSdy, ".summary"),
            to = paste0(mxPathOnVirt, ".summary"))

  file.copy(from = paste0(mxPathOnSdy, ".summary.orig"),
            to = paste0(mxPathOnVirt, ".summary.orig"))
}
