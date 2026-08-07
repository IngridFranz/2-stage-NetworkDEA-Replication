dea_labels <- c(
  healthcare_expenditure = "Healthcare expenditure",
  digital_health_infrastructure = "Digital health infrastructure",
  gp_density = "GP density",
  hospital_bed_density = "Hospital bed density",
  avoidable_diabetes_admissions = "Avoidable diabetes-related admissions",
  bed_days = "Bed-days",
  physician_consultations = "Physician consultations",
  treatable_mortality = "Treatable mortality",
  sah = "Self-assessed health"
)

run_cluster_comparisons <- function(context, dea_variables, primary_assignments,
                                    output_dir) {
  x <- merge(primary_assignments[c("dmu", "cluster")], dea_variables,
             by = "dmu", all.x = TRUE, sort = FALSE)
  if (nrow(x) != 27L || anyNA(x)) stop("DEA-variable/cluster merge failed.")
  descriptives <- do.call(rbind, lapply(names(dea_labels), function(v) {
    do.call(rbind, lapply(levels(x$cluster), function(g) {
      values <- x[x$cluster == g, v]
      data.frame(Variable = unname(dea_labels[v]), Cluster = g, N = length(values),
                 Mean = mean(values), SD = stats::sd(values),
                 Median = stats::median(values), Minimum = min(values),
                 Maximum = max(values))
    }))
  }))
  tests <- do.call(rbind, lapply(names(dea_labels), function(v) {
    test <- stats::kruskal.test(x[[v]] ~ x$cluster)
    data.frame(Variable = unname(dea_labels[v]), N = nrow(x),
               Kruskal_Wallis_chi_square = unname(test$statistic),
               df = unname(test$parameter), P_value = test$p.value)
  }))
  tests$P_Holm_across_nine_variables <- stats::p.adjust(tests$P_value, "holm")
  write_result(descriptives, output_dir, "10_dea_variables_by_cluster.csv")
  write_result(tests, output_dir, "11_kruskal_wallis_dea_variables.csv")
  invisible(list(data = x, descriptives = descriptives, tests = tests))
}

