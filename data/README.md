# Analytic input files

The repository does not redistribute source data. Retrieve the indicators from
the sources reported in the manuscript and Supplementary Tables 1 and 2, then
construct the following UTF-8 comma-separated files using the supplied
templates and `metadata/data_construction_protocol.csv`.

## `analysis_data.csv`

Contains the variables required for the final network DEA and original scales
used in reporting. `epa_inv` and `treatmort_inv` must equal
`max(x) - x + 1`. The workflow verifies both transformations.

## `context_data.csv`

Contains the five continuous contextual indicators and health-system type.
Health-system type is treated explicitly as nominal in the Gower dissimilarity.

Both files must contain the same 30 unique two-letter `dmu` identifiers. Blank
cells are not permitted in the final analytic files. Values described as
imputed in Supplementary Table 1 must be entered after applying the documented
mean-imputation rule. Exact numerical reproduction requires the same source
vintages, reference periods, harmonisation decisions, and imputed values used
for the article.
