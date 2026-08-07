efficiency_outcomes <- c(eff_stage1 = "Stage 1 efficiency",
                         eff_stage2 = "Stage 2 efficiency",
                         eff_overall = "Overall efficiency")

fit_fractional <- function(formula, data) {
  model <- stats::glm(formula, data = data,
                      family = stats::quasibinomial(link = "logit"))
  list(model = model, vcov = hc1_vcov(model))
}

factor_design <- function(model, variable, level, levels_all) {
  nd <- setNames(data.frame(factor(level, levels = levels_all)), variable)
  as.numeric(stats::model.matrix(stats::delete.response(stats::terms(model)), nd)[1, ])
}

factor_outputs <- function(fit, variable, data, outcome_label) {
  model <- fit$model; V <- fit$vcov; beta <- stats::coef(model)
  levels_all <- levels(data[[variable]])
  margins <- do.call(rbind, lapply(levels_all, function(g) {
    X <- factor_design(model, variable, g, levels_all)
    p <- stats::plogis(sum(X * beta)); grad <- p * (1 - p) * X
    se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
    data.frame(Outcome = outcome_label, Group = g, N = sum(data[[variable]] == g),
               Predicted_efficiency = p, Robust_SE = se,
               CI_lower = max(0, p - stats::qnorm(.975) * se),
               CI_upper = min(1, p + stats::qnorm(.975) * se))
  }))
  pairs <- utils::combn(levels_all, 2, simplify = FALSE)
  contrasts <- do.call(rbind, lapply(pairs, function(pair) {
    xa <- factor_design(model, variable, pair[1], levels_all)
    xb <- factor_design(model, variable, pair[2], levels_all)
    pa <- stats::plogis(sum(xa * beta)); pb <- stats::plogis(sum(xb * beta))
    grad <- pa * (1 - pa) * xa - pb * (1 - pb) * xb
    difference <- pa - pb
    se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
    p <- 2 * stats::pnorm(abs(difference / se), lower.tail = FALSE)
    data.frame(Outcome = outcome_label, Group_1 = pair[1], Group_2 = pair[2],
               Difference = difference, Robust_SE = se,
               CI_lower = difference - stats::qnorm(.975) * se,
               CI_upper = difference + stats::qnorm(.975) * se,
               P_unadjusted = p)
  }))
  contrasts$P_Holm <- stats::p.adjust(contrasts$P_unadjusted, "holm")
  global <- wald_test(model, V); global$Outcome <- outcome_label
  global$N <- stats::nobs(model)
  list(margins = margins, contrasts = contrasts,
       global = global[c("Outcome", "N", "Wald_chi_square", "df", "P_value")])
}

continuous_ame <- function(fit, variable, increment, outcome_label, data) {
  model <- fit$model; V <- fit$vcov; beta <- stats::coef(model)
  X <- stats::model.matrix(model)
  p <- stats::plogis(drop(X %*% beta))
  slope_index <- match(variable, names(beta))
  ame_unit <- mean(p * (1 - p) * beta[slope_index])
  gradient_rows <- p * (1 - p) * (1 - 2 * p) * beta[slope_index] * X
  gradient_rows[, slope_index] <- gradient_rows[, slope_index] + p * (1 - p)
  gradient <- colMeans(gradient_rows) * increment
  estimate <- ame_unit * increment
  se <- sqrt(as.numeric(t(gradient) %*% V %*% gradient))
  data.frame(Outcome = outcome_label, Indicator = variable, Increment = increment,
             N = nrow(data), AME = estimate, Robust_SE = se,
             CI_lower = estimate - stats::qnorm(.975) * se,
             CI_upper = estimate + stats::qnorm(.975) * se,
             P_value = 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE))
}

run_fractional_logit <- function(context, efficiency, primary_assignments,
                                 output_dir) {
  tol <- 1e-7
  for (v in names(efficiency_outcomes))
    efficiency[[v]] <- pmin(pmax(efficiency[[v]], 0), 1)
  full <- merge(context, efficiency, by = "dmu", all = FALSE, sort = FALSE)
  if (nrow(full) != 30L) stop("Context/efficiency merge failed.")

  cluster_data <- merge(primary_assignments[c("dmu", "cluster")], efficiency,
                        by = "dmu", all.x = TRUE, sort = FALSE)
  cluster_data$cluster <- factor(cluster_data$cluster)
  hcs_data <- full[full$hcs_type != 5, , drop = FALSE]
  hcs_data$hcs_type <- factor(hcs_data$hcs_type)

  efficiency_descriptives <- do.call(rbind, lapply(
    names(efficiency_outcomes), function(y) {
      do.call(rbind, lapply(levels(cluster_data$cluster), function(g) {
        values <- cluster_data[cluster_data$cluster == g, y]
        data.frame(Outcome = unname(efficiency_outcomes[y]), Cluster = g,
                   N = length(values), Mean = mean(values), SD = stats::sd(values),
                   Median = stats::median(values), Minimum = min(values),
                   Maximum = max(values))
      }))
    }))

  categorical_results <- function(data, variable) {
    outputs <- lapply(names(efficiency_outcomes), function(y) {
      fit <- fit_fractional(stats::reformulate(variable, y), data)
      factor_outputs(fit, variable, data, unname(efficiency_outcomes[y]))
    })
    list(margins = do.call(rbind, lapply(outputs, `[[`, "margins")),
         contrasts = do.call(rbind, lapply(outputs, `[[`, "contrasts")),
         global = do.call(rbind, lapply(outputs, `[[`, "global")))
  }
  cluster_models <- categorical_results(cluster_data, "cluster")
  hcs_models <- categorical_results(hcs_data, "hcs_type")

  increments <- c(prev_index = 10, gdp_cap = 10000, gini = .1,
                  rural_pop = 10, elderly_pop = 10)
  continuous <- do.call(rbind, lapply(names(efficiency_outcomes), function(y) {
    do.call(rbind, lapply(names(increments), function(v) {
      fit <- fit_fractional(stats::reformulate(v, y), full)
      continuous_ame(fit, v, increments[[v]], unname(efficiency_outcomes[y]), full)
    }))
  }))

  write_result(efficiency_descriptives, output_dir,
               "19_efficiency_descriptives_by_cluster.csv")
  write_result(cluster_models$margins, output_dir, "20_cluster_predictive_margins.csv")
  write_result(cluster_models$global, output_dir, "21_cluster_global_wald_tests.csv")
  write_result(cluster_models$contrasts, output_dir, "22_cluster_pairwise_contrasts.csv")
  write_result(hcs_models$margins, output_dir, "23_hcs_predictive_margins.csv")
  write_result(hcs_models$global, output_dir, "24_hcs_global_wald_tests.csv")
  write_result(hcs_models$contrasts, output_dir, "25_hcs_pairwise_contrasts.csv")
  write_result(continuous, output_dir, "26_continuous_indicator_average_marginal_effects.csv")
  invisible(list(cluster = cluster_models, hcs = hcs_models, continuous = continuous))
}
