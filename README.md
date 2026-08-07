# Two-stage network DEA and contextual analyses: replication code

Version 2.0.0 provides one continuous R workflow for the final input-oriented
VRS network DEA and the sensitivity, cluster, and exploratory contextual
analyses reported in the revised manuscript.

## Scope

The workflow:

1. validates the two analytic input tables and the reported inversions;
2. performs synthetic numerical checks;
3. estimates the bounded two-stage VRS network DEA with equal stage weights,
   a common endogenously projected undesirable link, and no secondary tie-break;
4. estimates CRS, alternative-stage-weight, targeted-variable-exclusion, and
   leave-one-country-out sensitivity specifications;
5. performs the single-linkage diagnostic, GDP Tukey-fence check, primary
   four-cluster complete-linkage analysis using Gower dissimilarity, candidate
   two- to six-cluster comparisons, and the all-country and no-HCS-type checks;
6. exports the primary dendrogram, cluster profiles, Kruskal-Wallis tests,
   fractional-logit predictive margins, robust Wald tests, Holm-adjusted
   contrasts, and average marginal effects.

## Data availability

Source data are not redistributed in this repository. Users must retrieve the
indicators from the sources and reference periods documented in the manuscript,
Supplementary Tables 1 and 2, and
`metadata/data_construction_protocol.csv`. The supplied templates define the
required file structure. The complete analytical workflow is reproducible;
exact numerical reproduction additionally requires reconstruction of the two
analytic input files from the documented sources.

Create these UTF-8 comma-separated files:

- `data/analysis_data.csv`
- `data/context_data.csv`

See `data/README.md`, the templates, and `metadata/variable_dictionary.csv`.

## Run the complete workflow

Install the two contributed packages once:

```r
install.packages(c("lpSolve", "cluster"))
```

Open R in the repository root and run:

```r
source("run_all.R")
```

Outputs are written to:

- `results/01_network_dea/`
- `results/01b_dea_sensitivity/`
- `results/02_cluster_context/`

To recreate only the primary dendrogram after supplying the context data, run:

```r
source("make_dendrogram.R")
```

The run stops if required input, transformation, solver, constraint, link,
VRS, or score-bound checks fail. Exact package versions and `sessionInfo()` are
written to `results/`.

## Expected numerical messages

The model deliberately applies no secondary tie-break. R may therefore warn
that some stage scores or projected links are not unique on the optimal face.
These warnings document alternative optimal solutions; they do not indicate a
failed solver or feasibility diagnostic. Numerical ranges are exported for
inspection.

## Tested environment

The integrated workflow was tested with R 4.4.1, `lpSolve` 5.6-23, and
`cluster` 2.1-6. Runtime versions recorded by each execution are authoritative.

## Citation

The software author and GitHub repository are recorded in `CITATION.cff`.
Add the associated article citation and the new version-specific Zenodo DOI
after they become available. Cite the version-specific Zenodo record when
referring to the exact submitted code.

## Licence

Unless otherwise stated, the source code and original documentation in this
repository are licensed under the MIT License. This licence does not apply to
third-party data or materials referenced by the analysis. No OECD, Eurostat,
European Commission, WHO, or other third-party data are redistributed here.
