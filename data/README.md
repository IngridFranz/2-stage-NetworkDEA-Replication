# Analysis data

The source data are not redistributed in this repository. They originate from
the OECD and the additional sources cited in the manuscript and may be subject
to their respective access and reuse conditions.

Create `analysis_data.csv` in this directory using the columns in
`input_template.csv`. Definitions are provided in
`../metadata/variable_dictionary.csv`. The file must contain one row per
healthcare system and no missing values.

The file must represent the harmonised model input after the transformations
and imputations reported in the manuscript:

- `epa_inv = max(epa_original) - epa_original + 1`;
- `treatmort_inv = max(treatmort_original) - treatmort_original + 1`;
- documented imputed values enter `beddays_imp` and `consult_imp`.

Do not commit `analysis_data.csv` unless redistribution is permitted by the
source licences and the journal's data-sharing policy.
