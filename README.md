UpdateAnno
==========

## Purpose
The `UpdateAnno` package is designed to ensure that the Gene Symbols used throughout the ImmuneSpace portal are consistent and up to date.  There are a number of places where Gene Symbols play an important role:
- Gene Expression Explorer - Pulls Gene Symbols from the FeatureAnnotation query
- Gene Set Enrichment Analysis - Pulls Gene Symbols from the static gene set module data objects in this package, e.g. chaussabel or emory.
- Differential Gene Expression Analysis - Uses the GeneExpressionAnalysisResults query to display differentially expressed genes.
 

## Inner Workings
During installation, the package's vignette `updateAnnoGMT.Rmd` also installs the latest version of the `org.Hs.eq.db` R package and uses the ALIAS-to-SYMBOL mapping within `org.Hs.eq.db` to update the `updated_btm_list` and `updated_btm_df` objects in the packages data folder.  Since the updating process depends on the vignette, it is important to ensure that the vignette is built during installation.

## Use
To install with the vignette: 
`devtools::install_github("RGLab/UpdateAnno", build_vignettes=T)`

After installation, a user may access these updated data objects by doing:
`library(updateAnno)`

This will include the data objects in the calling environment.

Then on `rsT` or `rsP`, the user should call the `runUpdateAnno()` function which will do the following:
- Update the FeatureAnnotation table with the latest gene symbols
- Update the gene expression matrices (flat files / tsv) to use the latest symbols in the summary versions
- Update the GeneExpressionAnalysisResults table with newly processed results.

The package also includes two methods that may assist the gene set enrichment analysis process:
 * `UpdateAnnoVec` - a function that takes a character vector with gene aliases and maps the aliases to official gene SYMBOLS from the `org.Hs.eq.db`.
 * `outputUpdatedBTM` - a function that writes out a gmt file with the most recently updated BTM object.
