UpdateAnno
==========

## Purpose
The `UpdateAnno` package is designed to ensure that the Gene Symbols used throughout the ImmuneSpace portal are consistent and up to date.  There are a number of places where Gene Symbols play an important role:
- Gene Expression Explorer - Pulls Gene Symbols from the FeatureAnnotation query
- Gene Set Enrichment Analysis - Pulls Gene Symbols from the static gene set module data objects in this package, e.g. chaussabel or emory.
- Differential Gene Expression Analysis - Uses the GeneExpressionAnalysisResults query to display differentially expressed genes.
 

## Inner Workings
In the `data-raw` sub-directory, the `updateDataWithLatestHGNCMap.Rmd` Rmarkdown file is used to generate a data table from the latest HUGO Gene Nomenclature Committee dataset.  This resource is used instead of the NCBI database via biomaRt or the equivalent `org.Hs.eg.db` as those data sources contained mappings that were deemed incorrect (such as "ACTB" > "POTEF") during the ImmuneSignatures 2 project.

However, there are two edge-cases that must be handled:
1. An alias maps to itself as a symbol as well as other symbols.  In this case, we have selected to use the self-mapping and removed any other mappings.  Our rationale is that the other-mapping symbols are historic artifacts and no longer accurate.
2. An alias maps to multiple symbols that do not include itself.  In this case, we drop these aliases from the matrix because we do not have a good way to know which symbol is the most accurate. 

## Use
To install with the vignette: 
`devtools::install_github("RGLab/UpdateAnno")`

After installation, a user may access these updated data objects by doing:
`library(updateAnno)`

This will include the data objects in the calling environment.

Then on `rsT` or `rsP`, the user should call the `runUpdateAnno()` function which will do the following:
- Update the FeatureAnnotation table with the latest gene symbols
- Update the gene expression matrices (flat files / tsv) to use the latest symbols in the summary versions
- Update the GeneExpressionAnalysisResults table with newly processed results.

The package also includes methods that may assist in processes for the ImmuneSpace platform:
 * `mapAlias2Symbol` - Function takes a character vector with gene aliases and maps the aliases to official gene SYMBOLS from the HGNC.
 * `outputUpdatedBTM` - Function writes out a gmt file with the most recently updated BTM object.
 * `makeAnnoDF` - Function takes an R package of annotation and creates a FeatureAnnotationSet tsv suitable for upload to the ImmuneSpace platform.
 * `addFasToVirtualSdy` - Function adds feature annotation to a virtual study container
 * `copyMatricesToVirtualSdy` - Function copies matrices (flat files) to virtual study in /share folder
