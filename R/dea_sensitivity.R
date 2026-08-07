frontier_flag <- function(x, tolerance = 1e-7) abs(x - 1) <= tolerance

safe_spearman <- function(x, y) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  unname(stats::cor(x, y, method = "spearman"))
}

fit_sensitivity_spec <- function(data, label, rts = "vrs",
                                 weights = c(.5, .5),
                                 stage1 = c("exp", "hospbeds", "GPs", "epa_inv"),
                                 stage2_outputs = c("treatmort_inv", "SAH"),
                                 tolerance = 1e-7) {
  fit <- withCallingHandlers(
    run_endogenous_link_network_dea(
      data = data, id_var = "dmu", stage1_inputs = stage1,
      link_var = "diab", stage2_inputs = c("beddays_imp", "consult_imp"),
      stage2_outputs = stage2_outputs, stage_weights = weights, rts = rts,
      exclude_incomplete = FALSE, tolerance = tolerance,
      optimal_face_tolerance = 1e-8, tie_break = "none"
    ),
    warning = function(w) {
      if (grepl("non-unique", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  fit$efficiency_results$Specification <- label
  fit$diagnostic_results$Specification <- label
  fit
}

compare_fit <- function(fit, baseline, tolerance = 1e-7) {
  b <- baseline[, c("DMU", "Eff_Stage1", "Eff_Stage2", "Eff_Overall")]
  x <- fit$efficiency_results[, c("DMU", "Eff_Stage1", "Eff_Stage2", "Eff_Overall")]
  z <- merge(b, x, by = "DMU", suffixes = c("_Baseline", "_Current"))
  metrics <- c("Stage1", "Stage2", "Overall")
  out <- do.call(rbind, lapply(metrics, function(m) {
    a <- z[[paste0("Eff_", m, "_Baseline")]]
    c <- z[[paste0("Eff_", m, "_Current")]]
    data.frame(
      Metric = m, N = length(a), Spearman = safe_spearman(a, c),
      Mean_baseline = mean(a), Mean_current = mean(c),
      Mean_absolute_difference = mean(abs(c - a)),
      Maximum_absolute_difference = max(abs(c - a)),
      Stable_frontier_status = mean(frontier_flag(a, tolerance) ==
                                      frontier_flag(c, tolerance)),
      Frontier_current = sum(frontier_flag(c, tolerance))
    )
  }))
  out$Specification <- unique(fit$efficiency_results$Specification)
  out[, c("Specification", setdiff(names(out), "Specification"))]
}

run_dea_sensitivity <- function(data, baseline_fit, output_dir,
                                tolerance = 1e-7) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  baseline <- baseline_fit$efficiency_results
  specifications <- list(
    CRS_50_50 = list(rts = "crs", weights = c(.5, .5)),
    VRS_75_25 = list(rts = "vrs", weights = c(.75, .25)),
    VRS_25_75 = list(rts = "vrs", weights = c(.25, .75)),
    VRS_without_digital_health = list(
      rts = "vrs", weights = c(.5, .5),
      stage1 = c("exp", "hospbeds", "GPs")),
    VRS_without_SAH = list(
      rts = "vrs", weights = c(.5, .5),
      stage2_outputs = "treatmort_inv"),
    VRS_without_digital_health_and_SAH = list(
      rts = "vrs", weights = c(.5, .5),
      stage1 = c("exp", "hospbeds", "GPs"),
      stage2_outputs = "treatmort_inv")
  )
  fits <- lapply(names(specifications), function(label) {
    s <- specifications[[label]]
    do.call(fit_sensitivity_spec, c(list(data = data, label = label,
      tolerance = tolerance), s))
  })
  names(fits) <- names(specifications)
  summaries <- do.call(rbind, lapply(fits, compare_fit,
                                     baseline = baseline,
                                     tolerance = tolerance))
  efficiency <- do.call(rbind, lapply(fits, `[[`, "efficiency_results"))
  diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostic_results"))
  write_result(summaries, output_dir, "01_specification_summary.csv")
  write_result(efficiency, output_dir, "02_all_efficiency_results.csv")
  write_result(diagnostics, output_dir, "03_all_diagnostics.csv")

  loo <- lapply(seq_len(nrow(data)), function(i) {
    omitted <- data$dmu[i]
    fit <- fit_sensitivity_spec(data[-i, , drop = FALSE],
      paste0("leave_out_", omitted), tolerance = tolerance)
    e <- fit$efficiency_results
    e$Omitted_DMU <- omitted
    e
  })
  loo_all <- do.call(rbind, loo)
  write_result(loo_all, output_dir, "04_leave_one_out_all_results.csv")
  loo_summary <- do.call(rbind, lapply(loo, function(e) {
    b <- baseline[baseline$DMU %in% e$DMU, ]
    b <- b[match(e$DMU, b$DMU), ]
    data.frame(
      Omitted_DMU = unique(e$Omitted_DMU), N = nrow(e),
      Spearman_Stage1 = safe_spearman(b$Eff_Stage1, e$Eff_Stage1),
      Spearman_Stage2 = safe_spearman(b$Eff_Stage2, e$Eff_Stage2),
      Spearman_Overall = safe_spearman(b$Eff_Overall, e$Eff_Overall),
      Maximum_absolute_difference_Stage1 = max(abs(b$Eff_Stage1-e$Eff_Stage1)),
      Maximum_absolute_difference_Stage2 = max(abs(b$Eff_Stage2-e$Eff_Stage2)),
      Maximum_absolute_difference_Overall = max(abs(b$Eff_Overall-e$Eff_Overall)),
      Stable_frontier_status_Overall = mean(frontier_flag(b$Eff_Overall,tolerance)==
                                               frontier_flag(e$Eff_Overall,tolerance))
    )
  }))
  write_result(loo_summary, output_dir, "05_leave_one_out_summary.csv")
  invisible(list(specifications = fits, leave_one_out = loo))
}
