# Complete reproduction: final VRS network DEA followed by cluster and
# exploratory contextual analyses. Run from the repository root.

options(stringsAsFactors = FALSE, ndea.functions_only = TRUE)

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
code_files <- c(
  "R/utils.R",
  "R/network_dea_model.R",
  "R/dea_sensitivity.R",
  "R/cluster_analysis.R",
  "R/create_dendrogram.R",
  "R/cluster_comparisons.R",
  "R/fractional_logit.R"
)
missing_code <- code_files[!file.exists(file.path(repo_root, code_files))]
if (length(missing_code)) {
  stop("Run run_all.R from the repository root. Missing: ",
       paste(missing_code, collapse = ", "))
}

required_packages <- c("lpSolve", "cluster")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install required package(s): install.packages(c(",
       paste(sprintf("'%s'", missing_packages), collapse = ", "), "))")
}

input_files <- c(
  analysis = file.path(repo_root, "data", "analysis_data.csv"),
  context = file.path(repo_root, "data", "context_data.csv")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Missing input file(s): ", paste(missing_inputs, collapse = ", "),
       ". See data/README.md and metadata/variable_dictionary.csv.")
}

output_root <- file.path(repo_root, "results")
overwrite <- identical(tolower(Sys.getenv("REPLICATION_OVERWRITE")), "true")
if (dir.exists(output_root) && length(list.files(output_root)) && !overwrite) {
  stop("results/ is not empty. Move it, or set REPLICATION_OVERWRITE=true.")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dea_output <- file.path(output_root, "01_network_dea")
dea_sensitivity_output <- file.path(output_root, "01b_dea_sensitivity")
context_output <- file.path(output_root, "02_cluster_context")
dir.create(dea_output, recursive = TRUE, showWarnings = FALSE)
dir.create(dea_sensitivity_output, recursive = TRUE, showWarnings = FALSE)
dir.create(context_output, recursive = TRUE, showWarnings = FALSE)

for (file in code_files) source(file.path(repo_root, file), local = FALSE)

analysis_data <- read_csv_strict(input_files[["analysis"]])
context_data <- read_csv_strict(input_files[["context"]])

required_analysis <- c(
  "dmu", "exp", "hospbeds", "GPs", "epa", "epa_inv", "diab",
  "beddays_imp", "consult_imp", "treatmort", "treatmort_inv", "SAH"
)
require_columns(analysis_data, required_analysis, "analysis_data.csv")
assert_complete_unique(analysis_data, "dmu", "analysis_data.csv")
if (!all(vapply(analysis_data[setdiff(required_analysis, "dmu")],
                is.numeric, logical(1)))) {
  stop("All analysis variables except dmu must be numeric.")
}
validate_context(context_data)
if (!identical(sort(analysis_data$dmu), sort(context_data$dmu))) {
  stop("DMU identifiers differ between analysis_data.csv and context_data.csv.")
}

# Verify both reported linear inversions exactly. The +1 constants and maxima
# are derived from the observed 30-country sample, as described in the paper.
tol <- 1e-7
expected_epa_inv <- max(analysis_data$epa) - analysis_data$epa + 1
expected_treatmort_inv <- max(analysis_data$treatmort) - analysis_data$treatmort + 1
if (max(abs(expected_epa_inv - analysis_data$epa_inv)) > tol) {
  stop("epa_inv does not equal max(epa) - epa + 1.")
}
if (max(abs(expected_treatmort_inv - analysis_data$treatmort_inv)) > tol) {
  stop("treatmort_inv does not equal max(treatmort) - treatmort + 1.")
}

# Synthetic test before empirical estimation.
synthetic <- run_endogenous_link_synthetic_test(tolerance = tol)
if (!all(synthetic$tests$Passed)) stop("A synthetic DEA check failed.")
write_result(synthetic$tests, dea_output, "00_synthetic_numerical_checks.csv")

# Final prespecified main model.
dea_result <- run_endogenous_link_network_dea(
  data = analysis_data,
  rts = "vrs",
  stage_weights = c(0.5, 0.5),
  exclude_incomplete = FALSE,
  tolerance = tol,
  optimal_face_tolerance = 1e-8,
  tie_break = "none"
)
export_endogenous_link_results(
  dea_result,
  output_directory = dea_output,
  prefix = "endogenous_link_vrs"
)

d <- dea_result$diagnostic_results
required_diagnostics_passed <-
  all(d$Primary_Solver_Status == 0L) &&
  all(d$Selected_Solver_Status == 0L) &&
  all(d$Max_Constraint_Violation <= tol) &&
  all(d$Score_Range_OK) &&
  all(!d$Theta1_Above_One) &&
  all(!d$Theta2_Above_One)
if (!required_diagnostics_passed) {
  stop("Network DEA completed, but at least one required diagnostic failed.")
}

# Prespecified DEA sensitivity and robustness analyses reported in the paper.
run_dea_sensitivity(
  data = analysis_data,
  baseline_fit = dea_result,
  output_dir = dea_sensitivity_output,
  tolerance = tol
)

# Pass DEA results directly to all downstream analyses.
efficiency <- dea_result$efficiency_results[c(
  "DMU", "Eff_Stage1", "Eff_Stage2", "Eff_Overall"
)]
names(efficiency) <- c("dmu", "eff_stage1", "eff_stage2", "eff_overall")

dea_reporting_variables <- analysis_data[c(
  "dmu", "exp", "epa", "GPs", "hospbeds", "diab", "beddays_imp",
  "consult_imp", "treatmort", "SAH"
)]
names(dea_reporting_variables) <- c(
  "dmu", "healthcare_expenditure", "digital_health_infrastructure",
  "gp_density", "hospital_bed_density", "avoidable_diabetes_admissions",
  "bed_days", "physician_consultations", "treatable_mortality", "sah"
)

cluster_results <- run_cluster_analysis(context_data, context_output)
create_primary_dendrogram(cluster_results, context_output)
run_cluster_comparisons(
  context = context_data,
  dea_variables = dea_reporting_variables,
  primary_assignments = cluster_results$primary_assignments,
  output_dir = context_output
)
run_fractional_logit(
  context = context_data,
  efficiency = efficiency,
  primary_assignments = cluster_results$primary_assignments,
  output_dir = context_output
)

software <- data.frame(
  Software = c("R", "lpSolve", "cluster"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(utils::packageVersion("lpSolve")),
    as.character(utils::packageVersion("cluster"))
  )
)
write_result(software, output_root, "98_software_versions.csv")
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_root, "99_sessionInfo.txt"))

manifest <- data.frame(
  Section = c("Network DEA", "DEA sensitivity", "Cluster/context"),
  Output_directory = c("results/01_network_dea",
                       "results/01b_dea_sensitivity",
                       "results/02_cluster_context"),
  Status = c("completed and diagnostics passed", "completed", "completed")
)
write_result(manifest, output_root, "97_run_manifest.csv")

cat("\nComplete reproduction succeeded.\nOutput: ",
    normalizePath(output_root, winslash = "/"), "\n", sep = "")
