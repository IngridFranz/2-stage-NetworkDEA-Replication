cluster_continuous <- c("elderly_pop", "gini", "rural_pop", "gdp_cap", "prev_index")
candidate_k <- 2:6
excluded_high_gdp <- c("IE", "LU", "NO")

gower_distance <- function(x, include_hcs = TRUE) {
  vars <- cluster_continuous
  if (include_hcs) vars <- c(vars, "hcs_type")
  z <- x[vars]
  if (include_hcs) z$hcs_type <- factor(z$hcs_type)
  cluster::daisy(z, metric = "gower",
                 type = if (include_hcs) list(nominal = "hcs_type") else list())
}

within_dispersion <- function(D, idx) {
  if (length(idx) <= 1L) return(0)
  sum(D[idx, idx, drop = FALSE]^2) / (2 * length(idx))
}

total_within <- function(D, groups) {
  sum(vapply(split(seq_along(groups), groups),
             function(i) within_dispersion(D, i), numeric(1)))
}

distance_ch <- function(D, groups) {
  n <- length(groups); k <- length(unique(groups))
  total <- within_dispersion(D, seq_len(n)); within <- total_within(D, groups)
  ((total - within) / (k - 1)) / (within / (n - k))
}

duda_hart <- function(tree, D, k) {
  parent <- stats::cutree(tree, k); child <- stats::cutree(tree, k + 1L)
  split_parent <- unique(parent)[vapply(unique(parent), function(g)
    length(unique(child[parent == g])) > 1L, logical(1))]
  if (length(split_parent) != 1L) return(c(NA, NA, NA))
  members <- which(parent == split_parent)
  children <- split(members, child[members])
  je1 <- within_dispersion(D, members)
  je2 <- sum(vapply(children, function(i) within_dispersion(D, i), numeric(1)))
  c(je2 / je1, if (je2 > 0) ((je1 - je2) / je2) * (length(members) - 2) else NA,
    length(members))
}

evaluate_tree <- function(tree, distance, specification) {
  D <- as.matrix(distance)
  do.call(rbind, lapply(candidate_k, function(k) {
    groups <- stats::cutree(tree, k); sil <- cluster::silhouette(groups, distance)
    sizes <- table(groups); dh <- duda_hart(tree, D, k)
    data.frame(Specification = specification, Number_of_clusters = k,
      Distance_based_CH = distance_ch(D, groups),
      Duda_Hart_Je2_over_Je1 = dh[1], Pseudo_T_squared = dh[2],
      Duda_split_cluster_size = dh[3], Mean_silhouette = mean(sil[, 3]),
      Minimum_silhouette = min(sil[, 3]),
      Negative_silhouettes = sum(sil[, 3] < 0),
      Minimum_cluster_size = min(sizes), Maximum_cluster_size = max(sizes),
      Singleton_clusters = sum(sizes == 1L))
  }))
}

assignment_table <- function(x, tree, distance, specification) {
  do.call(rbind, lapply(candidate_k, function(k) {
    groups <- stats::cutree(tree, k); sil <- cluster::silhouette(groups, distance)
    data.frame(Specification = specification, Number_of_clusters = k,
      country = x$country, dmu = x$dmu, hcs_type = x$hcs_type,
      Cluster = groups, Silhouette_width = sil[, 3])
  }))
}

fit_cluster_spec <- function(x, include_hcs, label) {
  distance <- gower_distance(x, include_hcs)
  tree <- stats::hclust(distance, method = "complete")
  list(data = x, distance = distance, tree = tree,
       criteria = evaluate_tree(tree, distance, label),
       assignments = assignment_table(x, tree, distance, label))
}

save_cluster_pdf <- function(result, path, title) {
  grDevices::pdf(path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(result$tree, labels = result$data$country, main = title,
                 xlab = "", sub = "", hang = -1, cex = .7)
  for (k in candidate_k) graphics::plot(
    cluster::silhouette(stats::cutree(result$tree, k), result$distance),
    main = paste(title, "- silhouette, k =", k), border = NA)
}

run_cluster_analysis <- function(context, output_dir) {
  context$hcs_type <- factor(context$hcs_type)
  all_distance <- gower_distance(context, TRUE)
  single_tree <- stats::hclust(all_distance, method = "single")
  grDevices::pdf(file.path(output_dir, "01_single_linkage_dendrogram.pdf"), 11, 8.5)
  graphics::plot(single_tree, labels = context$country,
                 main = "Preliminary single-linkage diagnostic: all 30 countries",
                 xlab = "", sub = "", hang = -1, cex = .7)
  grDevices::dev.off()

  single_assignments <- do.call(rbind, lapply(candidate_k, function(k) {
    groups <- stats::cutree(single_tree, k)
    sizes <- table(groups)
    data.frame(Number_of_clusters = k, country = context$country,
               dmu = context$dmu, Cluster = groups,
               Cluster_size = as.integer(sizes[as.character(groups)]))
  }))
  D_single <- as.matrix(all_distance)
  diag(D_single) <- Inf
  nearest <- apply(D_single, 1, which.min)
  nearest_neighbours <- data.frame(
    country = context$country, dmu = context$dmu,
    nearest_country = context$country[nearest],
    nearest_dmu = context$dmu[nearest],
    nearest_Gower_distance = apply(D_single, 1, min)
  )
  write_result(single_assignments, output_dir,
               "01a_single_linkage_assignments.csv")
  write_result(nearest_neighbours, output_dir,
               "01b_single_linkage_nearest_neighbours.csv")

  q1 <- unname(stats::quantile(context$gdp_cap, .25))
  q3 <- unname(stats::quantile(context$gdp_cap, .75)); iqr <- q3 - q1
  fence <- q3 + 1.5 * iqr
  outlier_check <- data.frame(country = context$country, dmu = context$dmu,
    gdp_cap = context$gdp_cap, upper_Tukey_fence = fence,
    above_upper_fence = context$gdp_cap > fence)
  observed <- sort(outlier_check$dmu[outlier_check$above_upper_fence])
  if (!identical(observed, sort(excluded_high_gdp)))
    stop("GDP Tukey-fence countries are ", paste(observed, collapse = ", "),
         "; expected IE, LU, NO. Review the input data.")
  write_result(outlier_check, output_dir, "02_gdp_tukey_outlier_check.csv")

  restricted <- context[!context$dmu %in% excluded_high_gdp, , drop = FALSE]
  specs <- list(
    main_27_hcs = fit_cluster_spec(restricted, TRUE, "Main: 27 countries, nominal HCS type"),
    all_30_hcs = fit_cluster_spec(context, TRUE, "Sensitivity: all 30 countries, nominal HCS type"),
    restricted_27_no_hcs = fit_cluster_spec(restricted, FALSE, "Sensitivity: 27 countries, no HCS type"),
    all_30_no_hcs = fit_cluster_spec(context, FALSE, "Additional: all 30 countries, no HCS type")
  )
  criteria <- do.call(rbind, lapply(specs, `[[`, "criteria"))
  assignments <- do.call(rbind, lapply(specs, `[[`, "assignments"))
  write_result(criteria, output_dir, "03_cluster_selection_criteria.csv")
  write_result(assignments, output_dir, "04_country_assignments_silhouettes.csv")

  primary <- subset(specs$main_27_hcs$assignments, Number_of_clusters == 4)
  # Stable manuscript labels: order clusters by their presentation in Table 4,
  # without changing any country assignment or recomputing the hierarchy.
  primary$cluster <- factor(
    primary$Cluster,
    levels = 1:4,
    labels = c("A", "D", "B", "C")
  )
  primary <- primary[c("country", "dmu", "hcs_type", "cluster", "Silhouette_width")]
  write_result(primary, output_dir, "05_primary_four_cluster_assignments.csv")

  profile_data <- merge(primary[c("dmu", "cluster")], restricted,
                        by = "dmu", all.x = TRUE, sort = FALSE)
  context_profiles <- do.call(rbind, lapply(levels(primary$cluster), function(g) {
    z <- profile_data[profile_data$cluster == g, ]
    do.call(rbind, lapply(cluster_continuous, function(v) {
      data.frame(Cluster = g, Variable = v, N = nrow(z),
                 Mean = mean(z[[v]]), SD = stats::sd(z[[v]]),
                 Median = stats::median(z[[v]]), Minimum = min(z[[v]]),
                 Maximum = max(z[[v]]))
    }))
  }))
  hcs_composition <- as.data.frame(table(
    Cluster = profile_data$cluster, HCS_type = profile_data$hcs_type
  ))
  hcs_composition$Proportion_within_cluster <- ave(
    hcs_composition$Freq, hcs_composition$Cluster,
    FUN = function(v) v / sum(v)
  )
  write_result(context_profiles, output_dir,
               "05a_primary_cluster_context_profiles.csv")
  write_result(hcs_composition, output_dir,
               "05b_primary_cluster_hcs_composition.csv")

  comparison_pairs <- list(
    c("main_27_hcs", "all_30_hcs"),
    c("main_27_hcs", "restricted_27_no_hcs"),
    c("all_30_hcs", "all_30_no_hcs")
  )
  ari <- do.call(rbind, lapply(comparison_pairs, function(pair) {
    do.call(rbind, lapply(candidate_k, function(k) {
      a <- subset(specs[[pair[1]]]$assignments, Number_of_clusters == k,
                  c("dmu", "Cluster")); names(a)[2] <- "a"
      b <- subset(specs[[pair[2]]]$assignments, Number_of_clusters == k,
                  c("dmu", "Cluster")); names(b)[2] <- "b"
      z <- merge(a, b, by = "dmu")
      data.frame(Specification_A = pair[1], Specification_B = pair[2],
                 Number_of_clusters = k, N_common = nrow(z),
                 Adjusted_Rand_Index = adjusted_rand_index(z$a, z$b))
    }))
  }))
  write_result(ari, output_dir, "06_assignment_stability_ari.csv")

  save_cluster_pdf(specs$main_27_hcs,
    file.path(output_dir, "07_main_cluster_diagnostics.pdf"), "Main analysis")
  save_cluster_pdf(specs$all_30_hcs,
    file.path(output_dir, "08_sensitivity_all30_diagnostics.pdf"), "All 30 countries")
  save_cluster_pdf(specs$restricted_27_no_hcs,
    file.path(output_dir, "09_sensitivity_without_hcs_diagnostics.pdf"), "Without HCS type")
  list(primary_assignments = primary, specifications = specs)
}
