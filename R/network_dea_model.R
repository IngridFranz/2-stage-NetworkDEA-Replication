# =============================================================================
# TWO-STAGE NETWORK DEA WITH AN ENDOGENOUS UNDESIRABLE LINK
# =============================================================================
#
# Stage 1:
#   exp, hospbeds, GPs, epa_inv  ->  diab
#
# Stage 2:
#   projected diab + beddays_imp + consult_imp
#       -> treatmort_inv + SAH
#
# Core link constraints for evaluated DMU o:
#
#   sum_j lambda_j * diab_j = diab_target_o
#   sum_j mu_j     * diab_j = diab_target_o
#   diab_target_o <= diab_o
#
# Interpretation:
# - Stage 1 may benchmark the undesirable intermediate to an equal or lower
#   level than observed.
# - Stage 2 must inherit exactly that same projected level.
# - Stage 2 cannot independently or radially contract the link.
#
# Main model:
#   min w1 * theta1 + w2 * theta2
#   subject to 0 <= theta1 <= 1 and 0 <= theta2 <= 1
#
# The separate upper bounds ensure that both stage-specific radial input
# factors retain the conventional input-oriented DEA efficiency domain.
#
# The code also calculates ranges of theta1, theta2 and diab_target on the
# optimal overall-efficiency face. This identifies non-unique decompositions.
#
# Required package:
#   install.packages("lpSolve")
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Package and data checks
# -----------------------------------------------------------------------------

check_packages <- function() {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop(
      "Package `lpSolve` is required. Install it with:\n",
      "install.packages(\"lpSolve\")"
    )
  }

  invisible(TRUE)
}


validate_network_data <- function(
    data,
    id_var = "dmu",
    stage1_inputs = c("exp", "hospbeds", "GPs", "epa_inv"),
    link_var = "diab",
    stage2_inputs = c("beddays_imp", "consult_imp"),
    stage2_outputs = c("treatmort_inv", "SAH"),
    exclude_incomplete = FALSE
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  required_vars <- unique(c(
    id_var,
    stage1_inputs,
    link_var,
    stage2_inputs,
    stage2_outputs
  ))

  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0L) {
    stop(
      "Missing required variables: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  d <- data[, required_vars, drop = FALSE]
  d[[id_var]] <- as.character(d[[id_var]])

  if (anyNA(d[[id_var]]) || any(d[[id_var]] == "")) {
    stop("DMU identifiers must be complete and non-empty.")
  }

  if (anyDuplicated(d[[id_var]])) {
    duplicate_ids <- unique(d[[id_var]][duplicated(d[[id_var]])])

    stop(
      "Duplicate DMU identifiers: ",
      paste(duplicate_ids, collapse = ", ")
    )
  }

  numeric_vars <- setdiff(required_vars, id_var)

  non_numeric <- numeric_vars[
    !vapply(d[numeric_vars], is.numeric, logical(1))
  ]

  if (length(non_numeric) > 0L) {
    stop(
      "The following DEA variables are not numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }

  incomplete <- !stats::complete.cases(d)

  excluded <- data.frame(
    DMU = character(0),
    Reason = character(0),
    stringsAsFactors = FALSE
  )

  if (any(incomplete)) {
    excluded <- data.frame(
      DMU = d[[id_var]][incomplete],
      Reason = vapply(
        which(incomplete),
        function(i) {
          missing_here <- names(d)[is.na(d[i, ])]

          paste0(
            "Missing: ",
            paste(missing_here, collapse = ", ")
          )
        },
        character(1)
      ),
      stringsAsFactors = FALSE
    )

    if (!isTRUE(exclude_incomplete)) {
      stop(
        "Incomplete observations found. ",
        "Set `exclude_incomplete = TRUE` only if complete-case exclusion ",
        "is part of the prespecified analysis."
      )
    }

    d <- d[!incomplete, , drop = FALSE]
  }

  if (nrow(d) < 2L) {
    stop("At least two complete DMUs are required.")
  }

  for (v in numeric_vars) {
    x <- d[[v]]

    if (any(!is.finite(x))) {
      stop("Variable `", v, "` contains non-finite values.")
    }

    if (any(x < 0)) {
      stop("Variable `", v, "` contains negative values.")
    }
  }

  radial_inputs <- c(stage1_inputs, stage2_inputs)

  zero_radial_inputs <- apply(
    d[, radial_inputs, drop = FALSE],
    1,
    function(x) any(x <= 0)
  )

  if (any(zero_radial_inputs)) {
    stop(
      "Zero or negative values occur in radially contracted inputs for: ",
      paste(d[[id_var]][zero_radial_inputs], collapse = ", "),
      ". Do not add an arbitrary constant without methodological justification."
    )
  }

  if (any(d[[link_var]] <= 0)) {
    stop(
      "The undesirable link `",
      link_var,
      "` must be strictly positive."
    )
  }

  if (any(d[, stage2_outputs, drop = FALSE] <= 0)) {
    stop("All desirable Stage-2 outputs must be strictly positive.")
  }

  role_table <- data.frame(
    Variable = c(
      stage1_inputs,
      link_var,
      stage2_inputs,
      stage2_outputs
    ),
    Stage = c(
      rep("Stage 1", length(stage1_inputs)),
      "Stage 1 -> Stage 2",
      rep("Stage 2", length(stage2_inputs)),
      rep("Stage 2", length(stage2_outputs))
    ),
    Role = c(
      rep("Discretionary radial input", length(stage1_inputs)),
      "Endogenously projected undesirable intermediate; non-discretionary in Stage 2",
      rep("Discretionary radial input", length(stage2_inputs)),
      rep("Desirable final output", length(stage2_outputs))
    ),
    stringsAsFactors = FALSE
  )

  list(
    data = d,
    excluded = excluded,
    role_table = role_table,
    id_var = id_var,
    stage1_inputs = stage1_inputs,
    link_var = link_var,
    stage2_inputs = stage2_inputs,
    stage2_outputs = stage2_outputs
  )
}


prepare_matrices <- function(validated) {
  d <- validated$data

  x1 <- as.matrix(
    d[, validated$stage1_inputs, drop = FALSE]
  )

  z <- as.numeric(d[[validated$link_var]])

  x2 <- as.matrix(
    d[, validated$stage2_inputs, drop = FALSE]
  )

  y2 <- as.matrix(
    d[, validated$stage2_outputs, drop = FALSE]
  )

  storage.mode(x1) <- "double"
  storage.mode(z) <- "double"
  storage.mode(x2) <- "double"
  storage.mode(y2) <- "double"

  list(
    dmu = d[[validated$id_var]],
    x1 = x1,
    z = z,
    x2 = x2,
    y2 = y2
  )
}


# -----------------------------------------------------------------------------
# 2. LP helpers
# -----------------------------------------------------------------------------

add_constraint_row <- function(
    matrix_rows,
    directions,
    rhs_values,
    row,
    direction,
    rhs
) {
  list(
    matrix_rows = c(matrix_rows, list(row)),
    directions = c(directions, direction),
    rhs_values = c(rhs_values, rhs)
  )
}


solve_lp_model <- function(
    objective,
    direction,
    constraint_matrix,
    constraint_directions,
    constraint_rhs
) {
  fit <- lpSolve::lp(
    direction = direction,
    objective.in = objective,
    const.mat = constraint_matrix,
    const.dir = constraint_directions,
    const.rhs = constraint_rhs,
    all.int = FALSE,
    all.bin = FALSE,
    transpose.constraints = TRUE,
    compute.sens = 0
  )

  fit
}


append_constraint <- function(
    constraint_matrix,
    constraint_directions,
    constraint_rhs,
    row,
    direction,
    rhs
) {
  list(
    matrix = rbind(constraint_matrix, row),
    directions = c(constraint_directions, direction),
    rhs = c(constraint_rhs, rhs)
  )
}


optimize_on_face <- function(
    base_matrix,
    base_directions,
    base_rhs,
    primary_objective,
    primary_optimum,
    face_tolerance,
    secondary_objective,
    direction = c("min", "max"),
    extra_matrix = NULL,
    extra_directions = NULL,
    extra_rhs = NULL
) {
  direction <- match.arg(direction)

  face <- append_constraint(
    constraint_matrix = base_matrix,
    constraint_directions = base_directions,
    constraint_rhs = base_rhs,
    row = primary_objective,
    direction = "<=",
    rhs = primary_optimum + face_tolerance
  )

  A <- face$matrix
  dirs <- face$directions
  rhs <- face$rhs

  if (!is.null(extra_matrix)) {
    A <- rbind(A, extra_matrix)
    dirs <- c(dirs, extra_directions)
    rhs <- c(rhs, extra_rhs)
  }

  solve_lp_model(
    objective = secondary_objective,
    direction = direction,
    constraint_matrix = A,
    constraint_directions = dirs,
    constraint_rhs = rhs
  )
}


# -----------------------------------------------------------------------------
# 3. Build the DMU-specific joint model
# -----------------------------------------------------------------------------

build_joint_model <- function(
    x1,
    z,
    x2,
    y2,
    evaluated_index,
    stage_weights = c(0.5, 0.5),
    rts = c("vrs", "crs"),
    tolerance = 1e-7
) {
  rts <- match.arg(rts)

  if (
    length(stage_weights) != 2L ||
    any(!is.finite(stage_weights)) ||
    any(stage_weights <= 0) ||
    abs(sum(stage_weights) - 1) > tolerance
  ) {
    stop(
      "`stage_weights` must contain two positive values that sum to 1."
    )
  }

  n <- nrow(x1)
  m1 <- ncol(x1)
  m2 <- ncol(x2)
  s2 <- ncol(y2)
  o <- evaluated_index

  # Variable order:
  # theta1, theta2, z_target, lambda_1...lambda_n, mu_1...mu_n
  idx_theta1 <- 1L
  idx_theta2 <- 2L
  idx_z_target <- 3L
  idx_lambda <- 4L:(3L + n)
  idx_mu <- (4L + n):(3L + 2L * n)
  n_vars <- 3L + 2L * n

  primary_objective <- numeric(n_vars)
  primary_objective[idx_theta1] <- stage_weights[1]
  primary_objective[idx_theta2] <- stage_weights[2]

  rows <- list()
  dirs <- character(0)
  rhs <- numeric(0)

  add_row <- function(row, direction, value) {
    rows[[length(rows) + 1L]] <<- row
    dirs <<- c(dirs, direction)
    rhs <<- c(rhs, value)
  }

  # Conventional input-oriented efficiency domain for each stage.
  # lpSolve imposes non-negativity by default; these constraints add the
  # required separate upper bounds.
  row <- numeric(n_vars)
  row[idx_theta1] <- 1
  add_row(row, "<=", 1)

  row <- numeric(n_vars)
  row[idx_theta2] <- 1
  add_row(row, "<=", 1)

  # Stage-1 radial input constraints:
  # sum_j lambda_j x1_ji <= theta1 * x1_oi
  for (i in seq_len(m1)) {
    row <- numeric(n_vars)
    row[idx_theta1] <- -x1[o, i]
    row[idx_lambda] <- x1[, i]

    add_row(row, "<=", 0)
  }

  # Stage-1 link production:
  # sum_j lambda_j z_j = z_target
  row <- numeric(n_vars)
  row[idx_lambda] <- z
  row[idx_z_target] <- -1

  add_row(row, "=", 0)

  # One-sided undesirable-link improvement:
  # z_target <= z_o
  row <- numeric(n_vars)
  row[idx_z_target] <- 1

  add_row(row, "<=", z[o])

  # Stage-2 radial dedicated-input constraints:
  # sum_j mu_j x2_jt <= theta2 * x2_ot
  for (t in seq_len(m2)) {
    row <- numeric(n_vars)
    row[idx_theta2] <- -x2[o, t]
    row[idx_mu] <- x2[, t]

    add_row(row, "<=", 0)
  }

  # Stage-2 continuity / inherited link:
  # sum_j mu_j z_j = z_target
  row <- numeric(n_vars)
  row[idx_mu] <- z
  row[idx_z_target] <- -1

  add_row(row, "=", 0)

  # Desirable final outputs:
  # sum_j mu_j y2_jr >= y2_or
  for (r in seq_len(s2)) {
    row <- numeric(n_vars)
    row[idx_mu] <- y2[, r]

    add_row(row, ">=", y2[o, r])
  }

  # VRS convexity in each stage
  if (rts == "vrs") {
    row <- numeric(n_vars)
    row[idx_lambda] <- 1

    add_row(row, "=", 1)

    row <- numeric(n_vars)
    row[idx_mu] <- 1

    add_row(row, "=", 1)
  }

  list(
    objective = primary_objective,
    constraint_matrix = do.call(rbind, rows),
    constraint_directions = dirs,
    constraint_rhs = rhs,
    indices = list(
      theta1 = idx_theta1,
      theta2 = idx_theta2,
      z_target = idx_z_target,
      lambda = idx_lambda,
      mu = idx_mu
    ),
    n_vars = n_vars,
    rts = rts,
    stage_weights = stage_weights
  )
}


# -----------------------------------------------------------------------------
# 4. Solve one DMU and inspect uniqueness
# -----------------------------------------------------------------------------

solve_one_dmu <- function(
    x1,
    z,
    x2,
    y2,
    evaluated_index,
    dmu_names,
    stage_weights = c(0.5, 0.5),
    rts = c("vrs", "crs"),
    tolerance = 1e-7,
    optimal_face_tolerance = 1e-8,
    tie_break = c("none", "min_link")
) {
  rts <- match.arg(rts)
  tie_break <- match.arg(tie_break)

  o <- evaluated_index

  model <- build_joint_model(
    x1 = x1,
    z = z,
    x2 = x2,
    y2 = y2,
    evaluated_index = o,
    stage_weights = stage_weights,
    rts = rts,
    tolerance = tolerance
  )

  primary <- solve_lp_model(
    objective = model$objective,
    direction = "min",
    constraint_matrix = model$constraint_matrix,
    constraint_directions = model$constraint_directions,
    constraint_rhs = model$constraint_rhs
  )

  if (primary$status != 0L) {
    stop(
      "Primary LP failed for DMU `",
      dmu_names[o],
      "`; lpSolve status = ",
      primary$status
    )
  }

  primary_optimum <- primary$objval

  face_tolerance <- max(
    optimal_face_tolerance,
    abs(primary_optimum) * optimal_face_tolerance
  )

  idx <- model$indices

  # ---------------------------------------------------------------------------
  # Global ranges on the optimal overall-efficiency face
  # ---------------------------------------------------------------------------

  unit_theta1 <- numeric(model$n_vars)
  unit_theta1[idx$theta1] <- 1

  unit_theta2 <- numeric(model$n_vars)
  unit_theta2[idx$theta2] <- 1

  unit_link <- numeric(model$n_vars)
  unit_link[idx$z_target] <- 1

  theta1_min_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_theta1,
    direction = "min"
  )

  theta1_max_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_theta1,
    direction = "max"
  )

  theta2_min_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_theta2,
    direction = "min"
  )

  theta2_max_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_theta2,
    direction = "max"
  )

  link_min_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_link,
    direction = "min"
  )

  link_max_fit <- optimize_on_face(
    model$constraint_matrix,
    model$constraint_directions,
    model$constraint_rhs,
    model$objective,
    primary_optimum,
    face_tolerance,
    unit_link,
    direction = "max"
  )

  range_fits <- list(
    theta1_min_fit,
    theta1_max_fit,
    theta2_min_fit,
    theta2_max_fit,
    link_min_fit,
    link_max_fit
  )

  if (any(vapply(range_fits, function(x) x$status != 0L, logical(1)))) {
    stop(
      "At least one optimal-face diagnostic LP failed for DMU `",
      dmu_names[o],
      "`."
    )
  }

  theta1_min <- theta1_min_fit$solution[idx$theta1]
  theta1_max <- theta1_max_fit$solution[idx$theta1]
  theta2_min <- theta2_min_fit$solution[idx$theta2]
  theta2_max <- theta2_max_fit$solution[idx$theta2]
  link_min <- link_min_fit$solution[idx$z_target]
  link_max <- link_max_fit$solution[idx$z_target]

  # ---------------------------------------------------------------------------
  # Optional symmetric conceptual tie-break:
  # among overall-optimal solutions, choose the lowest undesirable link.
  #
  # This does not change the primary overall score.
  # ---------------------------------------------------------------------------

  selected <- primary

  if (tie_break == "min_link") {
    selected <- link_min_fit
  }

  sol <- selected$solution

  theta1 <- sol[idx$theta1]
  theta2 <- sol[idx$theta2]
  z_target <- sol[idx$z_target]
  lambda <- sol[idx$lambda]
  mu <- sol[idx$mu]

  overall <- sum(stage_weights * c(theta1, theta2))

  x1_peer <- colSums(x1 * lambda)
  x2_peer <- colSums(x2 * mu)
  y2_peer <- colSums(y2 * mu)

  z_stage1 <- sum(lambda * z)
  z_stage2 <- sum(mu * z)

  x1_radial <- theta1 * x1[o, ]
  x2_radial <- theta2 * x2[o, ]

  stage1_input_violation <- max(
    pmax(x1_peer - x1_radial, 0)
  )

  stage2_input_violation <- max(
    pmax(x2_peer - x2_radial, 0)
  )

  output_violation <- max(
    pmax(y2[o, ] - y2_peer, 0)
  )

  stage1_link_residual <- abs(z_stage1 - z_target)
  stage2_link_residual <- abs(z_stage2 - z_target)
  link_continuity_residual <- abs(z_stage1 - z_stage2)
  one_sided_link_violation <- max(z_target - z[o], 0)

  rts_residual <- if (rts == "vrs") {
    max(
      abs(sum(lambda) - 1),
      abs(sum(mu) - 1)
    )
  } else {
    NA_real_
  }

  objective_residual <- abs(overall - primary_optimum)

  max_constraint_violation <- max(
    c(
      stage1_input_violation,
      stage2_input_violation,
      output_violation,
      stage1_link_residual,
      stage2_link_residual,
      link_continuity_residual,
      one_sided_link_violation,
      if (is.na(rts_residual)) 0 else rts_residual
    )
  )

  score_range_ok <- all(
    is.finite(c(theta1, theta2, overall))
  ) &&
    theta1 >= -tolerance &&
    theta2 >= -tolerance &&
    theta1 <= 1 + tolerance &&
    theta2 <= 1 + tolerance &&
    overall >= -tolerance &&
    overall <= 1 + tolerance

  # These flags are retained as explicit numerical checks of the LP bounds.
  theta1_above_one <- theta1 > 1 + tolerance
  theta2_above_one <- theta2 > 1 + tolerance

  efficiency <- data.frame(
    DMU = dmu_names[o],
    RTS = toupper(rts),
    Eff_Stage1 = theta1,
    Eff_Stage2 = theta2,
    Eff_Overall = overall,
    Observed_Link = z[o],
    Projected_Link = z_target,
    Link_Reduction = z[o] - z_target,
    Link_Reduction_Percent = 100 * (z[o] - z_target) / z[o],
    Tie_Break = tie_break,
    stringsAsFactors = FALSE
  )

  diagnostics <- data.frame(
    DMU = dmu_names[o],
    RTS = toupper(rts),
    Primary_Solver_Status = primary$status,
    Selected_Solver_Status = selected$status,
    Primary_Overall_Optimum = primary_optimum,
    Objective_Residual = objective_residual,
    Stage1_Link_Residual = stage1_link_residual,
    Stage2_Link_Residual = stage2_link_residual,
    Link_Continuity_Residual = link_continuity_residual,
    One_Sided_Link_Violation = one_sided_link_violation,
    RTS_Residual = rts_residual,
    Max_Constraint_Violation = max_constraint_violation,
    Score_Range_OK = score_range_ok,
    Theta1_Above_One = theta1_above_one,
    Theta2_Above_One = theta2_above_one,
    Self_Weight_Stage1 = lambda[o],
    Self_Weight_Stage2 = mu[o],
    Active_Peers_Stage1 = sum(lambda > tolerance),
    Active_Peers_Stage2 = sum(mu > tolerance),
    Theta1_Min_Overall_Optimum = theta1_min,
    Theta1_Max_Overall_Optimum = theta1_max,
    Theta1_Range = theta1_max - theta1_min,
    Theta2_Min_Overall_Optimum = theta2_min,
    Theta2_Max_Overall_Optimum = theta2_max,
    Theta2_Range = theta2_max - theta2_min,
    Link_Min_Overall_Optimum = link_min,
    Link_Max_Overall_Optimum = link_max,
    Link_Range = link_max - link_min,
    Stage_Scores_Unique = (
      theta1_max - theta1_min <= tolerance &&
      theta2_max - theta2_min <= tolerance
    ),
    Link_Unique = link_max - link_min <= tolerance,
    stringsAsFactors = FALSE
  )

  targets <- rbind(
    data.frame(
      DMU = dmu_names[o],
      Stage = "Stage 1",
      Variable_Type = "Discretionary input",
      Variable = colnames(x1),
      Observed = as.numeric(x1[o, ]),
      Radial_Target = as.numeric(x1_radial),
      Peer_Projection = as.numeric(x1_peer),
      Residual_Slack = pmax(
        as.numeric(x1_radial - x1_peer),
        0
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      DMU = dmu_names[o],
      Stage = "Link",
      Variable_Type = "Undesirable intermediate",
      Variable = "diab",
      Observed = z[o],
      Radial_Target = z_target,
      Peer_Projection = z_stage1,
      Residual_Slack = z[o] - z_target,
      stringsAsFactors = FALSE
    ),
    data.frame(
      DMU = dmu_names[o],
      Stage = "Stage 2",
      Variable_Type = "Discretionary input",
      Variable = colnames(x2),
      Observed = as.numeric(x2[o, ]),
      Radial_Target = as.numeric(x2_radial),
      Peer_Projection = as.numeric(x2_peer),
      Residual_Slack = pmax(
        as.numeric(x2_radial - x2_peer),
        0
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      DMU = dmu_names[o],
      Stage = "Stage 2",
      Variable_Type = "Desirable output",
      Variable = colnames(y2),
      Observed = as.numeric(y2[o, ]),
      Radial_Target = as.numeric(y2[o, ]),
      Peer_Projection = as.numeric(y2_peer),
      Residual_Slack = pmax(
        as.numeric(y2_peer - y2[o, ]),
        0
      ),
      stringsAsFactors = FALSE
    )
  )

  peer_stage1 <- data.frame(
    Evaluated_DMU = dmu_names[o],
    Peer_DMU = dmu_names,
    Lambda = lambda,
    stringsAsFactors = FALSE
  )

  peer_stage2 <- data.frame(
    Evaluated_DMU = dmu_names[o],
    Peer_DMU = dmu_names,
    Mu = mu,
    stringsAsFactors = FALSE
  )

  list(
    efficiency = efficiency,
    diagnostics = diagnostics,
    targets = targets,
    peer_stage1 = peer_stage1,
    peer_stage2 = peer_stage2
  )
}


# -----------------------------------------------------------------------------
# 5. Run all DMUs
# -----------------------------------------------------------------------------

run_endogenous_link_network_dea <- function(
    data,
    id_var = "dmu",
    stage1_inputs = c("exp", "hospbeds", "GPs", "epa_inv"),
    link_var = "diab",
    stage2_inputs = c("beddays_imp", "consult_imp"),
    stage2_outputs = c("treatmort_inv", "SAH"),
    stage_weights = c(0.5, 0.5),
    rts = c("vrs", "crs"),
    exclude_incomplete = FALSE,
    tolerance = 1e-7,
    optimal_face_tolerance = 1e-8,
    tie_break = c("none", "min_link")
) {
  check_packages()

  rts <- match.arg(rts)
  tie_break <- match.arg(tie_break)

  validated <- validate_network_data(
    data = data,
    id_var = id_var,
    stage1_inputs = stage1_inputs,
    link_var = link_var,
    stage2_inputs = stage2_inputs,
    stage2_outputs = stage2_outputs,
    exclude_incomplete = exclude_incomplete
  )

  matrices <- prepare_matrices(validated)

  solved <- lapply(
    seq_len(nrow(matrices$x1)),
    function(o) {
      solve_one_dmu(
        x1 = matrices$x1,
        z = matrices$z,
        x2 = matrices$x2,
        y2 = matrices$y2,
        evaluated_index = o,
        dmu_names = matrices$dmu,
        stage_weights = stage_weights,
        rts = rts,
        tolerance = tolerance,
        optimal_face_tolerance = optimal_face_tolerance,
        tie_break = tie_break
      )
    }
  )

  efficiency_results <- do.call(
    rbind,
    lapply(solved, `[[`, "efficiency")
  )

  # Scores numerically equal to one must share frontier rank 1. Rounding is
  # used only for ranking; exported efficiency estimates remain unchanged.
  rank_value <- function(x) {
    ifelse(abs(x - 1) <= tolerance, 1, x)
  }

  efficiency_results$Rank_Stage1 <- rank(
    -rank_value(efficiency_results$Eff_Stage1),
    ties.method = "min"
  )

  efficiency_results$Rank_Stage2 <- rank(
    -rank_value(efficiency_results$Eff_Stage2),
    ties.method = "min"
  )

  efficiency_results$Rank_Overall <- rank(
    -rank_value(efficiency_results$Eff_Overall),
    ties.method = "min"
  )

  diagnostic_results <- do.call(
    rbind,
    lapply(solved, `[[`, "diagnostics")
  )

  target_results <- do.call(
    rbind,
    lapply(solved, `[[`, "targets")
  )

  peer_stage1 <- do.call(
    rbind,
    lapply(solved, `[[`, "peer_stage1")
  )

  peer_stage2 <- do.call(
    rbind,
    lapply(solved, `[[`, "peer_stage2")
  )

  if (any(
    diagnostic_results$Max_Constraint_Violation > tolerance,
    na.rm = TRUE
  )) {
    warning(
      "At least one DMU has a constraint violation above tolerance."
    )
  }

  if (any(
    diagnostic_results$Objective_Residual >
      max(tolerance, optimal_face_tolerance * 10),
    na.rm = TRUE
  )) {
    warning(
      "At least one selected solution is not sufficiently close ",
      "to the primary overall optimum."
    )
  }

  if (any(!diagnostic_results$Stage_Scores_Unique)) {
    warning(
      "At least one DMU has non-unique stage scores at the optimal ",
      "overall efficiency. Inspect the reported score ranges."
    )
  }

  if (any(!diagnostic_results$Link_Unique)) {
    warning(
      "At least one DMU has a non-unique projected link at the optimal ",
      "overall efficiency. No secondary tie-break was applied; ",
      "inspect the reported link range."
    )
  }

  model_metadata <- data.frame(
    Item = c(
      "Model",
      "Orientation",
      "Returns to scale",
      "Link treatment",
      "Stage-1 weight",
      "Stage-2 weight",
      "Primary objective",
      "Stage-factor bounds",
      "Tie-break",
      "Tolerance",
      "Optimal-face tolerance",
      "Number of DMUs"
    ),
    Value = c(
      "Two-stage network DEA with an endogenously projected undesirable link",
      "Input-oriented",
      toupper(rts),
      paste0(
        link_var,
        "_target is common to both stages and <= observed ",
        link_var
      ),
      stage_weights[1],
      stage_weights[2],
      "w1 * theta1 + w2 * theta2",
      "0 <= theta1 <= 1; 0 <= theta2 <= 1",
      tie_break,
      tolerance,
      optimal_face_tolerance,
      nrow(validated$data)
    ),
    stringsAsFactors = FALSE
  )

  result <- list(
    efficiency_results = efficiency_results,
    diagnostic_results = diagnostic_results,
    target_results = target_results,
    peer_weights = list(
      stage1 = peer_stage1,
      stage2 = peer_stage2
    ),
    model_metadata = model_metadata,
    variable_roles = validated$role_table,
    excluded_dmus = validated$excluded,
    analysis_data = validated$data
  )

  class(result) <- c(
    "endogenous_link_network_dea",
    class(result)
  )

  result
}


print.endogenous_link_network_dea <- function(
    x,
    digits = 6,
    ...
) {
  ordered <- x$efficiency_results[
    order(x$efficiency_results$Rank_Overall),
  ]

  print(
    ordered,
    digits = digits,
    row.names = FALSE,
    ...
  )

  invisible(x)
}


# -----------------------------------------------------------------------------
# 6. Synthetic numerical checks
# -----------------------------------------------------------------------------

run_endogenous_link_synthetic_test <- function(
    tolerance = 1e-7
) {
  synthetic <- data.frame(
    dmu = c("A", "B", "C", "D", "E"),
    exp = c(1.0, 1.2, 1.5, 1.1, 2.0),
    hospbeds = c(1.0, 1.2, 1.5, 1.1, 2.0),
    GPs = c(1.0, 1.2, 1.5, 1.1, 2.0),
    epa_inv = c(1.0, 1.2, 1.5, 1.1, 2.0),
    diab = c(4.0, 5.0, 6.0, 4.5, 7.0),
    beddays_imp = c(1.0, 1.2, 1.5, 1.1, 2.0),
    consult_imp = c(1.0, 1.2, 1.5, 1.1, 2.0),
    treatmort_inv = c(5.0, 5.1, 4.9, 5.2, 4.0),
    SAH = c(5.0, 5.1, 4.9, 5.2, 4.0),
    stringsAsFactors = FALSE
  )

  base <- run_endogenous_link_network_dea(
    data = synthetic,
    rts = "vrs",
    tolerance = tolerance,
    tie_break = "none"
  )

  scaled <- synthetic
  scaled$exp <- 1000 * scaled$exp

  scaled_fit <- run_endogenous_link_network_dea(
    data = scaled,
    rts = "vrs",
    tolerance = tolerance,
    tie_break = "none"
  )

  comparison <- merge(
    base$efficiency_results[
      c("DMU", "Eff_Overall")
    ],
    scaled_fit$efficiency_results[
      c("DMU", "Eff_Overall")
    ],
    by = "DMU",
    suffixes = c("_base", "_scaled")
  )

  scale_difference <- max(
    abs(
      comparison$Eff_Overall_base -
      comparison$Eff_Overall_scaled
    )
  )

  diagnostics <- base$diagnostic_results
  efficiency <- base$efficiency_results

  tests <- data.frame(
    Test = c(
      "All primary LPs solved",
      "Stage-1 link equations",
      "Stage-2 link equations",
      "Link continuity",
      "Projected link never exceeds observed link",
      "Maximum constraint violation",
      "Objective reconstruction",
      "Positive unit-rescaling invariance"
    ),
    Passed = c(
      all(diagnostics$Primary_Solver_Status == 0),
      all(diagnostics$Stage1_Link_Residual <= tolerance),
      all(diagnostics$Stage2_Link_Residual <= tolerance),
      all(diagnostics$Link_Continuity_Residual <= tolerance),
      all(
        efficiency$Projected_Link <=
          efficiency$Observed_Link + tolerance
      ),
      all(
        diagnostics$Max_Constraint_Violation <= tolerance
      ),
      all(
        diagnostics$Objective_Residual <=
          max(tolerance, 1e-7)
      ),
      scale_difference <= tolerance
    ),
    Maximum = c(
      max(diagnostics$Primary_Solver_Status),
      max(diagnostics$Stage1_Link_Residual),
      max(diagnostics$Stage2_Link_Residual),
      max(diagnostics$Link_Continuity_Residual),
      max(
        efficiency$Projected_Link -
          efficiency$Observed_Link
      ),
      max(diagnostics$Max_Constraint_Violation),
      max(diagnostics$Objective_Residual),
      scale_difference
    ),
    stringsAsFactors = FALSE
  )

  list(
    data = synthetic,
    result = base,
    scaled_result = scaled_fit,
    tests = tests
  )
}


# -----------------------------------------------------------------------------
# 7. Export
# -----------------------------------------------------------------------------

export_endogenous_link_results <- function(
    result,
    output_directory = "endogenous_link_ndea_bounded_output",
    prefix = "endogenous_link_vrs_bounded"
) {
  if (!inherits(result, "endogenous_link_network_dea")) {
    stop(
      "`result` must be returned by ",
      "run_endogenous_link_network_dea()."
    )
  }

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  objects <- list(
    efficiency = result$efficiency_results,
    diagnostics = result$diagnostic_results,
    targets = result$target_results,
    peers_stage1 = result$peer_weights$stage1,
    peers_stage2 = result$peer_weights$stage2,
    model_metadata = result$model_metadata,
    variable_roles = result$variable_roles,
    excluded_dmus = result$excluded_dmus,
    analysis_data = result$analysis_data
  )

  for (nm in names(objects)) {
    utils::write.csv(
      objects[[nm]],
      file.path(
        output_directory,
        paste0(prefix, "_", nm, ".csv")
      ),
      row.names = FALSE,
      na = ""
    )
  }

  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(
      output_directory,
      paste0(prefix, "_sessionInfo.txt")
    )
  )

  invisible(normalizePath(output_directory))
}


# =============================================================================
# 8. APPLICATION
# =============================================================================
#
# install.packages("lpSolve")
# source("NDEA_endogenous_UODI_link_VRS.R")
#
# First run the synthetic checks:
#
if (!isTRUE(getOption("ndea.functions_only", FALSE))) {
test <- run_endogenous_link_synthetic_test()
print(test$tests)
#
# Then run the empirical VRS model:
#
result_vrs <- run_endogenous_link_network_dea(
data = data,
 rts = "vrs",
 stage_weights = c(0.5, 0.5),
 exclude_incomplete = FALSE,
 tolerance = 1e-7,
 optimal_face_tolerance = 1e-8,
tie_break = "none"
)
#
# Main results:
#
result_vrs$efficiency_results
#
#Diagnostics for Germany and Turkey:
#
subset(
result_vrs$efficiency_results,
 DMU %in% c("DE", "TR")
)
#
 subset(
 result_vrs$diagnostic_results,
DMU %in% c("DE", "TR")
)
#
# Active Stage-2 peers:
#
subset(
 result_vrs$peer_weights$stage2,
 Evaluated_DMU %in% c("DE", "TR") &
 Mu > 1e-8
)
#
# Key uniqueness check:
#
subset(
result_vrs$diagnostic_results,
 !Stage_Scores_Unique | !Link_Unique,
 select = c(
 DMU,
 Theta1_Min_Overall_Optimum,
 Theta1_Max_Overall_Optimum,
 Theta2_Min_Overall_Optimum,
 Theta2_Max_Overall_Optimum,
 Link_Min_Overall_Optimum,
 Link_Max_Overall_Optimum
   )
 )

# Export:

 export_endogenous_link_results(result_vrs)
}
#
# =============================================================================
