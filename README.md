UpdateAnno
==========

## Purpose
The `UpdateAnno` package is designed to support gene set enrichment analyses by automatically updating the gene symbols in an internal data object holding a blood transcription module created for the ImmuneSignatures project by the HIPC collaborators.

## Inner Workings
During installation, the package's vignette `updateAnnoGMT.Rmd` also installs the latest version of the `org.Hs.eq.db` R package and uses the ALIAS-to-SYMBOL mapping within `org.Hs.eq.db` to update the `updated_btm_list` and `updated_btm_df` objects in the packages data folder.  Since the updating process depends on the vignette, it is important to ensure that the vignette is built during installation.

## Use
After installation, a user may access these updated data objects by doing:
`library(updateAnno)`

This will include the data objects in the calling environment.

The package also includes two methods that may assist the gene set enrichment analysis process:
 * `UpdateAnnoVec` - a function that takes a character vector with gene aliases and maps the aliases to official gene SYMBOLS from the `org.Hs.eq.db`.
 * `outputUpdatedBTM` - a function that writes out a gmt file with the most recently updated BTM object.