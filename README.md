# Two-stage network DEA replication code

This repository contains the code for the final primary model reported in the
accompanying manuscript: a bounded, input-oriented two-stage network DEA model
under variable returns to scale (VRS), with an endogenously projected
undesirable intermediate.

## Model

Stage 1 uses healthcare expenditure, hospital-bed density, GP density, and the
transformed digital-health indicator as discretionary radial inputs. Avoidable
diabetes-related hospitalisations form the undesirable intermediate. Stage 2
uses bed-days and physician consultations as discretionary radial inputs and
transformed treatable mortality and self-assessed health as final desirable
outcomes.

Both stages reproduce the same endogenous projected intermediate. It may remain
unchanged or improve relative to the observed value and is inherited by Stage 2
as a non-discretionary link. It is not radially contracted in Stage 2. The
objective is `0.5 * theta1 + 0.5 * theta2`, with separate VRS convexity
conditions and `0 <= theta1, theta2 <= 1`. No secondary tie-break is applied.

## Contents

- `run_analysis.R`: single entry point;
- `R/network_dea_model.R`: model, diagnostics, projections, and exports;
- `data/input_template.csv`: required data structure;
- `data/README.md`: data-preparation instructions;
- `metadata/variable_dictionary.csv`: roles and transformations;
- `results/reference/`: reference efficiency and diagnostic outputs;
- `LICENSE`: MIT licence for the code.

## Data availability

Source data are not redistributed. They are available from the OECD and the
other sources cited in the manuscript, subject to provider terms. Reproduction
requires construction of `data/analysis_data.csv` following the manuscript and
`data/README.md`.

## Requirements

- R 4.4.1
- `lpSolve` 5.6.23

Install the package with `install.packages("lpSolve")`.

## Run

1. Create `data/analysis_data.csv` using the supplied template and dictionary.
2. Set the repository root as the R working directory.
3. Run `source("run_analysis.R")`.

The script runs synthetic numerical checks, estimates all healthcare systems,
checks solver and constraint diagnostics, and writes outputs to `results/`.

Peer weights can differ across equivalent optimal solutions or software
versions without changing efficiency scores. No secondary tie-break is used.

## Citation

Please cite the accompanying article. Add the final bibliographic reference and
repository DOI here after acceptance and repository archiving.

## Licence

Code is released under the MIT License. Data remain subject to the terms of
their original providers.
