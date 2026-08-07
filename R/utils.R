read_csv_strict <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8-BOM")
  if (!nrow(x)) stop("Input file has no data rows: ", path)
  x
}

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

assert_complete_unique <- function(x, id, label, expected_n = 30L) {
  if (nrow(x) != expected_n) stop(label, " must contain ", expected_n, " rows.")
  if (anyNA(x) || any(!stats::complete.cases(x))) stop(label, " contains missing values.")
  if (anyDuplicated(x[[id]])) stop(label, " contains duplicate identifiers.")
}

validate_context <- function(x) {
  required <- c("country", "dmu", "elderly_pop", "gini", "rural_pop",
                "gdp_cap", "prev_index", "hcs_type")
  require_columns(x, required, "context_data.csv")
  assert_complete_unique(x, "dmu", "context_data.csv")
  numeric_vars <- setdiff(required, c("country", "dmu", "hcs_type"))
  if (!all(vapply(x[numeric_vars], is.numeric, logical(1))))
    stop("All continuous context variables must be numeric.")
  if (length(unique(x$hcs_type)) != 5L)
    stop("Expected five health-system-type categories.")
}

validate_dea_variables <- function(x) {
  required <- c("dmu", "healthcare_expenditure", "digital_health_infrastructure",
                "gp_density", "hospital_bed_density",
                "avoidable_diabetes_admissions", "bed_days",
                "physician_consultations", "treatable_mortality", "sah")
  require_columns(x, required, "dea_variables.csv")
  assert_complete_unique(x, "dmu", "dea_variables.csv")
  if (!all(vapply(x[setdiff(required, "dmu")], is.numeric, logical(1))))
    stop("All DEA variables must be numeric and in their manuscript reporting scales.")
}

validate_efficiency <- function(x) {
  required <- c("dmu", "eff_stage1", "eff_stage2", "eff_overall")
  require_columns(x, required, "dea_efficiency.csv")
  assert_complete_unique(x, "dmu", "dea_efficiency.csv")
  tol <- 1e-7
  for (v in required[-1]) {
    if (!is.numeric(x[[v]]) || any(x[[v]] < -tol | x[[v]] > 1 + tol))
      stop(v, " must lie within [0,1], allowing only 1e-7 numerical tolerance.")
    x[[v]] <- pmin(pmax(x[[v]], 0), 1)
  }
  overall_check <- abs(x$eff_overall - 0.5 * (x$eff_stage1 + x$eff_stage2))
  if (max(overall_check) > tol) stop("Overall scores do not equal 0.5*(Stage 1 + Stage 2).")
  invisible(TRUE)
}

validate_country_sets <- function(...) {
  inputs <- list(...)
  sets <- lapply(inputs, function(x) sort(as.character(x$dmu)))
  if (!all(vapply(sets[-1], identical, logical(1), sets[[1]])))
    stop("DMU identifiers do not match across the three input files.")
}

write_result <- function(x, output_dir, filename) {
  utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}

choose_two <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sc <- sum(choose_two(tab)); sr <- sum(choose_two(rowSums(tab)))
  sk <- sum(choose_two(colSums(tab))); expected <- sr * sk / choose_two(n)
  maximum <- 0.5 * (sr + sk)
  if (maximum == expected) return(1)
  (sc - expected) / (maximum - expected)
}

hc1_vcov <- function(model) {
  X <- stats::model.matrix(model)
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted(model)
  w <- pmax(mu * (1 - mu), .Machine$double.eps)
  bread <- solve(crossprod(X, X * w))
  score <- X * as.numeric(y - mu)
  n <- nrow(X); p <- ncol(X)
  (n / (n - p)) * bread %*% crossprod(score) %*% bread
}

wald_test <- function(model, vcov, omit_intercept = TRUE) {
  indices <- seq_along(stats::coef(model))
  if (omit_intercept) indices <- indices[-1]
  b <- stats::coef(model)[indices]
  V <- vcov[indices, indices, drop = FALSE]
  chi2 <- as.numeric(t(b) %*% solve(V, b))
  data.frame(Wald_chi_square = chi2, df = length(indices),
             P_value = stats::pchisq(chi2, length(indices), lower.tail = FALSE))
}

