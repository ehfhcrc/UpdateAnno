# Evan Henrich
# October 2018
# Notes: SDY1368 uses the alternative probe id (just digits instead of ILMN_123123 format) for
# HumanHt_v4 from Illumina.  Therefore creating a customAnno.

library(data.table)

# get unigene from norm_exprs by cutting off suffix of probes
# SDY1368_EM.txt was just renamed from original supp file from GSE
# GSE86331_P98_2011DS75_Stanford_Batch16_2HT12_V4_BM_V4_0_R1_20Sept2011_Raw_Data_GX11.txt.gz
em <- data.table::fread("FeatureAnnotationSetDev/SDY1368_EM.txt")

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
