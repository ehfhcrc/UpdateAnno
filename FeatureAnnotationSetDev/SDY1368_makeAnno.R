# Evan Henrich
# October 2018
# Notes: SDY1368 uses the alternative probe id (just digits instead of ILMN_123123 format) for
# HumanHt_v4 from Illumina.  Therefore creating a customAnno.

library(data.table)

# SDY1368_EM.txt was just renamed from original supp file from GSE
# GSE86331_P98_2011DS75_Stanford_Batch16_2HT12_V4_BM_V4_0_R1_20Sept2011_Raw_Data_GX11.txt.gz
link <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86331&format=file&file=GSE86331%5FP98%5F2011DS75%5FStanford%5FBatch16%5F2HT12%5FV4%5FBM%5FV4%5F0%5FR1%5F20Sept2011%5FRaw%5FData%5FGX11%2Etxt%2Egz"
fl <- "FeatureAnnotationSetDev/SDY1368_EM.txt.gz"
download.file(link, fl)
GEOquery::gunzip(fl)
em <- fread(gsub("\\.gz", "", fl))

# Remove unmapped
res <- data.frame(Probe_ID = em$ProbeID,
                  Gene_Symbol = em$TargetID,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY1368_customAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# clean up
file.remove(gsub("\\.gz", "", fl))
