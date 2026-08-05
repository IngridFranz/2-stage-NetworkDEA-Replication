# Reproduce the final bounded, input-oriented VRS network DEA model.

if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop("Package `lpSolve` is required. Install it with install.packages('lpSolve').")
}

input_file <- file.path("data", "analysis_data.csv")
if (!file.exists(input_file)) {
  stop("Missing data/analysis_data.csv. See data/README.md.")
}

options(ndea.functions_only = TRUE)
source(file.path("R", "network_dea_model.R"))

analysis_data <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
synthetic_check <- run_endogenous_link_synthetic_test(tolerance = 1e-7)
if (!all(synthetic_check$tests$Passed)) stop("A synthetic numerical check failed.")

result_vrs <- run_endogenous_link_network_dea(
  data = analysis_data,
  rts = "vrs",
  stage_weights = c(0.5, 0.5),
  exclude_incomplete = FALSE,
  tolerance = 1e-7,
  optimal_face_tolerance = 1e-8,
  tie_break = "none"
)

export_endogenous_link_results(
  result_vrs,
  output_directory = "results",
  prefix = "endogenous_link_vrs"
)
utils::write.csv(synthetic_check$tests, file.path("results", "synthetic_numerical_checks.csv"), row.names = FALSE)

d <- result_vrs$diagnostic_results
passed <- all(d$Primary_Solver_Status == 0L) &&
  all(d$Selected_Solver_Status == 0L) &&
  all(d$Max_Constraint_Violation <= 1e-7) &&
  all(d$Score_Range_OK) &&
  all(!d$Theta1_Above_One) &&
  all(!d$Theta2_Above_One)
if (!passed) stop("The model completed, but at least one required diagnostic failed.")

message("Analysis completed successfully. Results are in results/.")
