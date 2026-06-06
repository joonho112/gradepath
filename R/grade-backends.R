# =============================================================================
# grade-backends.R -- grade-IP solver backends
# -----------------------------------------------------------------------------
# The default parity backend is a robust Gurobi cascade: in-process R binding,
# fenced R binding, gurobi_cl indicator LP, then gurobi_cl big-M LP. Explicit
# open backends are also supported: HiGHS uses the highs package directly,
# while GLPK and SYMPHONY route through ROI plugins.
# =============================================================================

.gp_gurobi_paths <- c(
  "in_process",
  "subprocess_R_binding",
  "gurobi_cl_indicator",
  "gurobi_cl_bigM"
)

# Explicit license-free backends. `highs` is a direct highs-package route;
# `glpk` and `symphony` are ROI plugin routes.
.gp_open_backends <- c("highs", "glpk", "symphony")

.gp_gurobi_cache <- new.env(parent = emptyenv())

.gradepath_abort_backend <- function(message, ..., call = FALSE) {
  .gradepath_abort(
    message,
    ...,
    class = c("gradepath_backend_error", "gradepath_error"),
    call = call
  )
}

.gradepath_abort_backend_status <- function(status,
                                            message,
                                            ...,
                                            class,
                                            call = FALSE) {
  .gp_status_abort(
    status,
    message,
    ...,
    class = unique(c(class, "gradepath_backend_error")),
    call = call
  )
}

.gradepath_abort_backend_unavailable <- function(message, ..., class = NULL, call = FALSE) {
  .gradepath_abort_backend_status(
    "SOLVER_BACKEND_UNAVAILABLE",
    message,
    ...,
    class = unique(c(class, "gradepath_backend_unavailable")),
    call = call
  )
}

.gradepath_abort_solver_infeasible <- function(message, ..., class = NULL, call = FALSE) {
  .gradepath_abort_backend_status(
    "SOLVER_INFEASIBLE",
    message,
    ...,
    class = unique(c(class, "gradepath_solver_infeasible")),
    call = call
  )
}

.gradepath_abort_solver_output_invalid <- function(message, ..., class = NULL, call = FALSE) {
  .gradepath_abort_backend_status(
    "SOLVER_OUTPUT_INVALID",
    message,
    ...,
    class = unique(c(class, "gradepath_solver_output_invalid")),
    call = call
  )
}

.gradepath_abort_solver_objective_mismatch <- function(message, ..., class = NULL, call = FALSE) {
  .gradepath_abort_backend_status(
    "SOLVER_OBJECTIVE_MISMATCH",
    message,
    ...,
    class = unique(c(class, "gradepath_solver_objective_mismatch")),
    call = call
  )
}

.gradepath_abort_solver_canonical_mismatch <- function(message, ..., class = NULL, call = FALSE) {
  .gradepath_abort_backend_status(
    "SOLVER_CANONICAL_MISMATCH",
    message,
    ...,
    class = unique(c(class, "gradepath_solver_canonical_mismatch")),
    call = call
  )
}

.gp_gurobi_abort_unavailable <- function(detect, failures = list(), call = FALSE) {
  detect <- .gp_gurobi_validate_detect(detect)
  failure_text <- ""
  if (length(failures) > 0L) {
    failure_text <- paste0(
      "\nAttempted path failures:\n",
      paste(
        sprintf("  - %s: %s", names(failures), unlist(failures, use.names = FALSE)),
        collapse = "\n"
      )
    )
  }

  .gp_status_abort(
    "SOLVER_BACKEND_UNAVAILABLE",
    paste0(
      "Gurobi is gradepath's strongly recommended default backend (KRW's ",
      "solver; by far the fastest on the grade integer program), but no stable ",
      "Gurobi invocation path is available.\n",
      "Detected: gurobi R package = %s; callr = %s; Matrix = %s; ",
      "gurobi_cl = %s%s; gurobi_cl smoke test = %s; ",
      "fenced R-binding smoke test = %s.%s\n",
      "Install Gurobi and activate a license first: ",
      "https://www.gurobi.com/academia/\n",
      "Then install the Gurobi R package from your Gurobi distribution and ",
      "ensure `gurobi_cl` is on PATH. As a license-free last resort (several ",
      "times slower on the real 97-firm solve) use ",
      "`gp_control(backend = \"highs\")`."
    ),
    .gp_bool_label(detect$gurobi),
    .gp_bool_label(detect$callr),
    .gp_bool_label(detect$Matrix),
    .gp_bool_label(detect$gurobi_cl),
    if (nzchar(detect$gurobi_cl_path)) paste0(" (", detect$gurobi_cl_path, ")") else "",
    .gp_bool_label(detect$gurobi_cl_smoke),
    .gp_bool_label(detect$binding_smoke),
    failure_text,
    class = c(
      "gradepath_gurobi_unavailable",
      "gradepath_backend_unavailable",
      "gradepath_backend_error"
    ),
    call = call
  )
}

.gp_bool_label <- function(x) {
  if (isTRUE(x)) {
    "yes"
  } else if (identical(x, FALSE)) {
    "no"
  } else {
    "not checked"
  }
}

.gp_system2 <- function(command,
                        args = character(),
                        stdout = "",
                        stderr = "",
                        timeout = 0) {
  system2(
    command = command,
    args = args,
    stdout = stdout,
    stderr = stderr,
    timeout = timeout
  )
}

.gp_gurobi_validate_detect <- function(detect) {
  required <- c(
    "gurobi",
    "callr",
    "Matrix",
    "jsonlite",
    "gurobi_cl",
    "gurobi_cl_path",
    "gurobi_cl_smoke",
    "binding_smoke"
  )
  if (!is.list(detect) || !identical(names(detect), required)) {
    .gradepath_abort_internal("Invalid Gurobi detection result.")
  }
  detect
}

.gp_gurobi_detect <- function(smoke = FALSE) {
  smoke <- .gradepath_validate_scalar_logical(smoke, "smoke")
  cli <- unname(Sys.which("gurobi_cl"))
  has_cli <- nzchar(cli)
  structure(
    list(
      gurobi = requireNamespace("gurobi", quietly = TRUE),
      callr = requireNamespace("callr", quietly = TRUE),
      Matrix = requireNamespace("Matrix", quietly = TRUE),
      jsonlite = requireNamespace("jsonlite", quietly = TRUE),
      gurobi_cl = has_cli,
      gurobi_cl_path = if (has_cli) cli else "",
      gurobi_cl_smoke = if (isTRUE(smoke)) {
        has_cli && .gp_gurobi_cl_smoke_ok()
      } else {
        NA
      },
      binding_smoke = if (isTRUE(smoke)) .gp_gurobi_binding_smoke_ok() else NA
    ),
    class = "gp_gurobi_detection"
  )
}

.gp_gurobi_cli_available <- function(detect) {
  detect <- .gp_gurobi_validate_detect(detect)
  isTRUE(detect$gurobi_cl) && isTRUE(detect$gurobi_cl_smoke)
}

.gp_gurobi_binding_smoke_ok <- function(timeout = 10) {
  timeout <- .gradepath_validate_scalar_numeric(
    timeout,
    "timeout",
    finite = TRUE,
    lower = 0,
    include_lower = FALSE
  )

  cache_key <- paste0("timeout_", timeout)
  if (exists(cache_key, envir = .gp_gurobi_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gp_gurobi_cache, inherits = FALSE))
  }

  ok <- FALSE
  if (requireNamespace("gurobi", quietly = TRUE) &&
      requireNamespace("callr", quietly = TRUE)) {
    res <- tryCatch(
      callr::r(
        function() {
          model <- list(
            A = matrix(1, nrow = 1L, ncol = 1L),
            obj = 1,
            sense = "<",
            rhs = 1,
            lb = 0,
            ub = 1,
            vtype = "B",
            modelsense = "min"
          )
          out <- gurobi::gurobi(model, list(OutputFlag = 0))
          list(ok = identical(out$status, "OPTIMAL"))
        },
        timeout = timeout,
        spinner = FALSE,
        package = FALSE
      ),
      error = function(e) list(ok = FALSE)
    )
    ok <- is.list(res) && isTRUE(res$ok)
  }

  assign(cache_key, ok, envir = .gp_gurobi_cache)
  ok
}

.gp_gurobi_cl_smoke_ok <- function(timeout = 10) {
  timeout <- .gradepath_validate_scalar_numeric(
    timeout,
    "timeout",
    finite = TRUE,
    lower = 0,
    include_lower = FALSE
  )

  cli <- unname(Sys.which("gurobi_cl"))
  if (!nzchar(cli)) {
    return(FALSE)
  }

  dir <- tempfile("gradepath-gurobi-cl-smoke-")
  dir.create(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  lp <- file.path(dir, "smoke.lp")
  json <- file.path(dir, "smoke.json")
  writeLines(
    c(
      "Minimize",
      " obj: x",
      "Subject To",
      " c1: x <= 1",
      "Bounds",
      " 0 <= x <= 1",
      "Binaries",
      " x",
      "End"
    ),
    lp,
    useBytes = TRUE
  )

  ok <- tryCatch(
    {
      out <- .gp_system2(
        cli,
        args = c("OutputFlag=0", "JSONSolDetail=1", paste0("ResultFile=", json), lp),
        stdout = TRUE,
        stderr = TRUE,
        timeout = timeout
      )
      status_code <- attr(out, "status")
      if (is.null(status_code)) {
        status_code <- 0L
      }
      identical(as.integer(status_code), 0L) && file.exists(json)
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  isTRUE(ok)
}

.gp_available_gurobi_paths <- function(smoke = TRUE) {
  detect <- .gp_gurobi_detect(smoke = smoke)
  paths <- character()
  if (isTRUE(detect$gurobi) &&
      isTRUE(detect$Matrix) &&
      isTRUE(detect$callr) &&
      isTRUE(detect$binding_smoke)) {
    paths <- c(paths, "in_process")
  }
  if (isTRUE(detect$gurobi) &&
      isTRUE(detect$Matrix) &&
      isTRUE(detect$callr) &&
      isTRUE(detect$binding_smoke)) {
    paths <- c(paths, "subprocess_R_binding")
  }
  if (.gp_gurobi_cli_available(detect)) {
    paths <- c(paths, "gurobi_cl_indicator", "gurobi_cl_bigM")
  }
  attr(paths, "detect") <- detect
  paths
}

.gp_gurobi_path_available <- function(path, detect) {
  path <- .gradepath_validate_scalar_character(
    path,
    "path",
    allowed = .gp_gurobi_paths
  )
  detect <- .gp_gurobi_validate_detect(detect)

  switch(
    path,
    in_process = isTRUE(detect$gurobi) &&
      isTRUE(detect$Matrix) &&
      isTRUE(detect$callr) &&
      isTRUE(detect$binding_smoke),
    subprocess_R_binding = isTRUE(detect$gurobi) &&
      isTRUE(detect$Matrix) &&
      isTRUE(detect$callr) &&
      isTRUE(detect$binding_smoke),
    gurobi_cl_indicator = .gp_gurobi_cli_available(detect),
    gurobi_cl_bigM = .gp_gurobi_cli_available(detect)
  )
}

.gp_open_backend_packages <- function(backend) {
  backend <- .gradepath_validate_scalar_character(
    backend,
    "backend",
    allowed = .gp_open_backends
  )

  switch(
    backend,
    highs = c("highs", "Matrix"),
    glpk = c("ROI", "ROI.plugin.glpk", "slam"),
    symphony = c("ROI", "ROI.plugin.symphony", "slam")
  )
}

.gp_open_backend_available <- function(backend) {
  packages <- .gp_open_backend_packages(backend)
  all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))
}

.gp_available_open_backends <- function() {
  .gp_open_backends[
    vapply(.gp_open_backends, .gp_open_backend_available, logical(1))
  ]
}

.gp_open_solver_controls <- function(control = NULL,
                                     solver_options = NULL,
                                     backend) {
  backend <- .gradepath_validate_scalar_character(
    backend,
    "backend",
    allowed = .gp_open_backends
  )
  if (!is.null(control) && !is.null(solver_options)) {
    .gradepath_abort_control(
      "Supply either `control` or `solver_options`, not both."
    )
  }
  if (!is.null(control)) {
    control <- validate_gp_control(control)
    solver_options <- control$solver_options
  }
  solver_options <- .gp_control_validate_solver_options(solver_options)
  names_options <- names(solver_options)

  controls <- list()
  consumed <- character()
  time_keys <- intersect(c("time_limit", "max_time"), names_options)
  gap_keys <- intersect(c("mip_gap", "mip_rel_gap"), names_options)

  if (length(time_keys) > 1L) {
    .gradepath_abort_control(
      "Do not supply multiple time-limit controls for backend `%s`: %s.",
      backend,
      paste(time_keys, collapse = ", ")
    )
  }
  if (length(gap_keys) > 1L) {
    .gradepath_abort_control(
      "Do not supply multiple MIP-gap controls for backend `%s`: %s.",
      backend,
      paste(gap_keys, collapse = ", ")
    )
  }

  if (length(time_keys) == 1L) {
    key <- time_keys[[1L]]
    value <- solver_options[[key]]
    consumed <- c(consumed, key)
    if (identical(backend, "highs")) {
      controls$time_limit <- value
    } else if (identical(backend, "glpk")) {
      controls$tm_limit <- 1000 * value
    } else {
      controls$time_limit <- value
    }
  }

  if (length(gap_keys) == 1L) {
    key <- gap_keys[[1L]]
    value <- solver_options[[key]]
    consumed <- c(consumed, key)
    if (identical(backend, "highs")) {
      controls$mip_rel_gap <- value
    } else if (identical(backend, "symphony")) {
      controls$gap_limit <- value
    } else {
      .gradepath_abort_control(
        "Backend `glpk` does not expose a MIP-gap control."
      )
    }
  }

  backend_specific <- solver_options[setdiff(names_options, consumed)]
  overlap <- intersect(names(controls), names(backend_specific))
  if (length(overlap) > 0L) {
    .gradepath_abort_control(
      "Backend `%s` received duplicate solver-control keys after mapping: %s.",
      backend,
      paste(overlap, collapse = ", ")
    )
  }

  defaults <- switch(
    backend,
    highs = list(output_flag = FALSE, log_to_console = FALSE),
    glpk = list(verbose = FALSE),
    symphony = list(verbosity = -2)
  )
  c(defaults, controls, backend_specific)
}

.gp_gurobi_params <- function(control = NULL, solver_options = NULL) {
  if (!is.null(control) && !is.null(solver_options)) {
    .gradepath_abort_control(
      "Supply either `control` or `solver_options`, not both."
    )
  }
  if (!is.null(control)) {
    control <- validate_gp_control(control)
    solver_options <- control$solver_options
  }
  solver_options <- .gp_control_validate_solver_options(solver_options)

  option_names <- names(solver_options)
  if (all(c("time_limit", "TimeLimit") %in% option_names)) {
    .gradepath_abort_control(
      "Do not supply both `time_limit` and Gurobi-native `TimeLimit`."
    )
  }
  if (all(c("mip_gap", "MIPGap") %in% option_names)) {
    .gradepath_abort_control(
      "Do not supply both `mip_gap` and Gurobi-native `MIPGap`."
    )
  }
  unsupported <- intersect(c("max_time", "mip_rel_gap"), option_names)
  if (length(unsupported) > 0L) {
    .gradepath_abort_control(
      "Gurobi does not use `%s`; use `time_limit` / `mip_gap` or Gurobi-native `TimeLimit` / `MIPGap`.",
      paste(unsupported, collapse = "`, `")
    )
  }

  params <- list(
    OutputFlag = 0,
    FeasibilityTol = 1e-9,
    IntFeasTol = 1e-9,
    OptimalityTol = 1e-9
  )
  for (nm in option_names) {
    value <- solver_options[[nm]]
    if (identical(nm, "time_limit")) {
      params$TimeLimit <- value
    } else if (identical(nm, "mip_gap")) {
      params$MIPGap <- value
    } else {
      params[[nm]] <- value
    }
  }

  params
}

.gp_gurobi_varnames <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  n <- problem$n_units
  d_names <- as.vector(t(outer(
    seq_len(n),
    seq_len(n),
    FUN = function(i, j) paste0("D_", i, "_", j)
  )))
  c(d_names, paste0("G_", seq_len(n)))
}

.gp_gurobi_validate_start <- function(problem, warm_start = NULL) {
  problem <- validate_gp_grade_problem(problem)
  if (is.null(warm_start)) {
    return(NULL)
  }
  warm_start <- .gradepath_validate_numeric_vector(warm_start, "warm_start")
  if (length(warm_start) != problem$constraint_ncol) {
    .gradepath_abort_backend(
      "`warm_start` must have length n^2 + n for the grade problem."
    )
  }
  warm_start
}

.gp_gurobi_pair_index <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  n <- problem$n_units
  pair_i <- rep(seq_len(n), each = n)
  pair_j <- rep(seq_len(n), times = n)
  off <- pair_i != pair_j
  data.frame(
    i = pair_i[off],
    j = pair_j[off],
    d = problem$d_indices[cbind(pair_i[off], pair_j[off])],
    stringsAsFactors = FALSE
  )
}

.gp_gurobi_indicator_gencon <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "The Gurobi R-binding paths require the suggested package `Matrix`."
    )
  }

  pairs <- .gp_gurobi_pair_index(problem)
  n_constraints <- 2L * nrow(pairs)
  gencon <- vector("list", n_constraints)

  for (k in seq_len(nrow(pairs))) {
    i <- pairs$i[[k]]
    j <- pairs$j[[k]]
    d <- pairs$d[[k]]
    gi <- problem$grade_indices[[i]]
    gj <- problem$grade_indices[[j]]
    gencon[[k]] <- list(
      binvar = d,
      binval = 1L,
      a = Matrix::sparseVector(
        i = c(gi, gj),
        x = c(1, -1),
        length = problem$constraint_ncol
      ),
      sense = ">",
      rhs = 1,
      name = paste0("ind_D_", i, "_", j, "_one")
    )
    gencon[[nrow(pairs) + k]] <- list(
      binvar = d,
      binval = 0L,
      a = Matrix::sparseVector(
        i = c(gj, gi),
        x = c(1, -1),
        length = problem$constraint_ncol
      ),
      sense = ">",
      rhs = 0,
      name = paste0("ind_D_", i, "_", j, "_zero")
    )
  }

  gencon
}

.gp_gurobi_model_indicator <- function(problem, warm_start = NULL) {
  problem <- validate_gp_grade_problem(problem)
  warm_start <- .gp_gurobi_validate_start(problem, warm_start)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "The Gurobi R-binding paths require the suggested package `Matrix`."
    )
  }

  model <- list(
    A = Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(0L, problem$constraint_ncol)
    ),
    obj = problem$objective,
    modelsense = "min",
    sense = character(),
    rhs = numeric(),
    vtype = problem$types,
    lb = problem$lower,
    ub = problem$upper,
    varnames = .gp_gurobi_varnames(problem),
    genconind = .gp_gurobi_indicator_gencon(problem)
  )
  if (!is.null(warm_start)) {
    model$start <- warm_start
  }
  model
}

.gp_gurobi_model_bigM <- function(problem, warm_start = NULL) {
  problem <- validate_gp_grade_problem(problem)
  warm_start <- .gp_gurobi_validate_start(problem, warm_start)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "The Gurobi R-binding paths require the suggested package `Matrix`."
    )
  }

  model <- list(
    A = Matrix::sparseMatrix(
      i = problem$constraint_i,
      j = problem$constraint_j,
      x = problem$constraint_x,
      dims = c(problem$constraint_nrow, problem$constraint_ncol)
    ),
    obj = problem$objective,
    modelsense = "min",
    sense = rep("<", problem$constraint_nrow),
    rhs = problem$rhs,
    vtype = problem$types,
    lb = problem$lower,
    ub = problem$upper,
    varnames = .gp_gurobi_varnames(problem),
    constrnames = paste0("bigM_", seq_len(problem$constraint_nrow))
  )
  if (!is.null(warm_start)) {
    model$start <- warm_start
  }
  model
}

.gp_gurobi_solve_inprocess <- function(problem, params, warm_start = NULL) {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable("The Gurobi R package is not installed.")
  }
  model <- .gp_gurobi_model_indicator(problem, warm_start = warm_start)
  gurobi::gurobi(model, params)
}

.gp_gurobi_callr_timeout <- function(params) {
  time_limit <- params$TimeLimit
  if (is.null(time_limit)) {
    return(300)
  }
  max(30, as.numeric(time_limit) + 30)
}

.gp_gurobi_cli_timeout <- function(params) {
  .gp_gurobi_callr_timeout(params)
}

.gp_gurobi_solve_subprocess <- function(problem, params, warm_start = NULL) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable("The `callr` package is not installed.")
  }
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable("The Gurobi R package is not installed.")
  }
  model <- .gp_gurobi_model_indicator(problem, warm_start = warm_start)
  callr::r(
    function(model, params) {
      if (!requireNamespace("gurobi", quietly = TRUE)) {
        stop("The Gurobi R package is not installed in the child process.")
      }
      gurobi::gurobi(model, params)
    },
    args = list(model = model, params = params),
    timeout = .gp_gurobi_callr_timeout(params),
    spinner = FALSE,
    package = FALSE
  )
}

.gp_gurobi_num <- function(x) {
  x <- as.numeric(x)
  x[abs(x) < 1e-15] <- 0
  sprintf("%.17g", x)
}

.gp_gurobi_lp_terms <- function(coefs, varnames, force_first = TRUE) {
  idx <- which(abs(coefs) > 1e-15)
  if (length(idx) == 0L) {
    if (isTRUE(force_first)) {
      return(paste0("0 ", varnames[[1L]]))
    }
    return("0")
  }
  vapply(
    seq_along(idx),
    function(k) {
      coef <- coefs[[idx[[k]]]]
      sign <- if (coef < 0) "-" else "+"
      body <- paste(.gp_gurobi_num(abs(coef)), varnames[[idx[[k]]]])
      if (k == 1L && coef >= 0 && !isTRUE(force_first)) {
        body
      } else {
        paste(sign, body)
      }
    },
    character(1)
  )
}

.gp_gurobi_write_lp_common <- function(problem, file, constraints) {
  problem <- validate_gp_grade_problem(problem)
  varnames <- .gp_gurobi_varnames(problem)
  objective_terms <- .gp_gurobi_lp_terms(problem$objective, varnames)
  d_names <- varnames[seq_len(problem$n_units * problem$n_units)]
  g_names <- varnames[problem$grade_indices]

  lines <- c(
    "Minimize",
    paste(" obj:", objective_terms[[1L]]),
    if (length(objective_terms) > 1L) paste0(" ", objective_terms[-1L]) else character(),
    "Subject To",
    constraints,
    "Bounds",
    paste0(
      " ",
      .gp_gurobi_num(problem$lower),
      " <= ",
      varnames,
      " <= ",
      .gp_gurobi_num(problem$upper)
    ),
    "Binaries",
    paste(" ", paste(d_names, collapse = " ")),
    "Generals",
    paste(" ", paste(g_names, collapse = " ")),
    "End"
  )
  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

.gp_gurobi_write_lp_indicator <- function(problem, file) {
  problem <- validate_gp_grade_problem(problem)
  pairs <- .gp_gurobi_pair_index(problem)
  varnames <- .gp_gurobi_varnames(problem)
  constraints <- " dummy_diagonal_projection: 0 D_1_1 <= 0"
  for (k in seq_len(nrow(pairs))) {
    i <- pairs$i[[k]]
    j <- pairs$j[[k]]
    d_name <- varnames[[pairs$d[[k]]]]
    gi <- varnames[[problem$grade_indices[[i]]]]
    gj <- varnames[[problem$grade_indices[[j]]]]
    constraints <- c(
      constraints,
      paste0(
        " ind_D_", i, "_", j, "_one: ",
        d_name,
        " = 1 -> ",
        gi,
        " - ",
        gj,
        " >= 1"
      ),
      paste0(
        " ind_D_", i, "_", j, "_zero: ",
        d_name,
        " = 0 -> ",
        gj,
        " - ",
        gi,
        " >= 0"
      )
    )
  }
  .gp_gurobi_write_lp_common(problem, file, constraints)
}

.gp_gurobi_write_lp_bigM <- function(problem, file) {
  problem <- validate_gp_grade_problem(problem)
  varnames <- .gp_gurobi_varnames(problem)
  constraints <- character(problem$constraint_nrow)
  for (row in seq_len(problem$constraint_nrow)) {
    idx <- which(problem$constraint_i == row)
    coefs <- numeric(problem$constraint_ncol)
    coefs[problem$constraint_j[idx]] <- problem$constraint_x[idx]
    constraints[[row]] <- paste0(
      " bigM_",
      row,
      ": ",
      paste(.gp_gurobi_lp_terms(coefs, varnames, force_first = FALSE), collapse = " "),
      " <= ",
      .gp_gurobi_num(problem$rhs[[row]])
    )
  }
  .gp_gurobi_write_lp_common(problem, file, constraints)
}

.gp_gurobi_write_prm <- function(params, file) {
  lines <- vapply(
    names(params),
    function(nm) paste(nm, as.character(params[[nm]])),
    character(1)
  )
  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

.gp_gurobi_cli_param_args <- function(params) {
  vapply(
    names(params),
    function(nm) paste0(nm, "=", as.character(params[[nm]])),
    character(1)
  )
}

.gp_gurobi_write_mst <- function(problem, warm_start, file) {
  problem <- validate_gp_grade_problem(problem)
  warm_start <- .gp_gurobi_validate_start(problem, warm_start)
  varnames <- .gp_gurobi_varnames(problem)
  idx <- which(is.finite(warm_start))
  lines <- c(
    "# MIP start",
    paste(varnames[idx], .gp_gurobi_num(warm_start[idx]))
  )
  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

.gp_gurobi_cli_timeout_text <- function(text) {
  text <- paste(as.character(text), collapse = "\n")
  grepl("timed out|timeout|elapsed time limit", text, ignore.case = TRUE)
}

.gp_gurobi_cli_output_text <- function(out) {
  text <- paste(as.character(out), collapse = "\n")
  if (nzchar(text)) text else "<no output>"
}

.gp_gurobi_cli_run <- function(cli, args, timeout) {
  warnings <- character()
  out <- tryCatch(
    withCallingHandlers(
      .gp_system2(
        cli,
        args = args,
        stdout = TRUE,
        stderr = TRUE,
        timeout = timeout
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(conditionMessage(e), status = 1L)
    }
  )

  status_code <- attr(out, "status")
  timeout_flag <- isTRUE(attr(out, "timeout")) ||
    .gp_gurobi_cli_timeout_text(c(as.character(out), warnings)) ||
    (!is.null(status_code) && identical(as.integer(status_code), 124L))
  if (is.null(status_code)) {
    status_code <- if (timeout_flag) {
      124L
    } else if (length(warnings) > 0L) {
      1L
    } else {
      0L
    }
  }
  structure(
    c(as.character(out), warnings),
    status = as.integer(status_code),
    timeout = timeout_flag
  )
}

.gp_gurobi_cli_status_label <- function(status) {
  if (is.null(status) || length(status) != 1L || is.na(status)) {
    return(NA_character_)
  }
  if (is.numeric(status)) {
    map <- c(
      `1` = "LOADED",
      `2` = "OPTIMAL",
      `3` = "INFEASIBLE",
      `4` = "INF_OR_UNBD",
      `5` = "UNBOUNDED",
      `6` = "CUTOFF",
      `7` = "ITERATION_LIMIT",
      `8` = "NODE_LIMIT",
      `9` = "TIME_LIMIT",
      `10` = "SOLUTION_LIMIT",
      `11` = "INTERRUPTED",
      `12` = "NUMERIC",
      `13` = "SUBOPTIMAL"
    )
    out <- unname(map[[as.character(as.integer(status))]])
    if (is.null(out)) NA_character_ else out
  } else {
    toupper(as.character(status))
  }
}

.gp_json_number <- function(text, key) {
  pattern <- paste0('"', key, '"\\s*:\\s*([-+0-9.eE]+)')
  hit <- regexec(pattern, text, perl = TRUE)
  value <- regmatches(text, hit)[[1L]]
  if (length(value) < 2L) {
    return(NA_real_)
  }
  as.numeric(value[[2L]])
}

.gp_gurobi_parse_json_base <- function(file) {
  text <- paste(readLines(file, warn = FALSE), collapse = "\n")
  vars <- list()
  var_hit <- gregexpr(
    '\\{\\s*"VarName"\\s*:\\s*"([^"]+)"\\s*,\\s*"X"\\s*:\\s*([-+0-9.eE]+)',
    text,
    perl = TRUE
  )
  var_match <- regmatches(text, var_hit)[[1L]]
  if (length(var_match) > 0L && !identical(var_match, character(0))) {
    vars <- lapply(var_match, function(x) {
      parts <- regmatches(
        x,
        regexec(
          '\\{\\s*"VarName"\\s*:\\s*"([^"]+)"\\s*,\\s*"X"\\s*:\\s*([-+0-9.eE]+)',
          x,
          perl = TRUE
        )
      )[[1L]]
      list(VarName = parts[[2L]], X = as.numeric(parts[[3L]]))
    })
  }
  list(
    SolutionInfo = list(
      Status = .gp_json_number(text, "Status"),
      ObjVal = .gp_json_number(text, "ObjVal"),
      ObjBound = .gp_json_number(text, "ObjBound"),
      MIPGap = .gp_json_number(text, "MIPGap"),
      Runtime = .gp_json_number(text, "Runtime")
    ),
    Vars = vars
  )
}

.gp_gurobi_parse_json <- function(file, varnames) {
  if (!file.exists(file)) {
    .gradepath_abort_solver_output_invalid("Gurobi JSON solution file was not written.")
  }
  sol <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(file, simplifyVector = FALSE)
  } else {
    .gp_gurobi_parse_json_base(file)
  }
  info <- sol$SolutionInfo
  primal <- stats::setNames(rep(0, length(varnames)), varnames)
  vars <- sol$Vars
  if (length(vars) > 0L) {
    for (v in vars) {
      nm <- v$VarName
      if (!is.null(nm) && nm %in% varnames) {
        primal[[nm]] <- as.numeric(v$X)
      }
    }
  }
  list(
    x = as.numeric(primal),
    objval = as.numeric(info$ObjVal),
    status = .gp_gurobi_cli_status_label(info$Status),
    objbound = if (!is.null(info$ObjBound)) as.numeric(info$ObjBound) else NA_real_,
    mipgap = if (!is.null(info$MIPGap)) as.numeric(info$MIPGap) else NA_real_,
    runtime = if (!is.null(info$Runtime)) as.numeric(info$Runtime) else NA_real_
  )
}

.gp_gurobi_parse_sol <- function(file, varnames) {
  if (!file.exists(file)) {
    .gradepath_abort_solver_output_invalid("Gurobi SOL solution file was not written.")
  }
  lines <- readLines(file, warn = FALSE)
  primal <- stats::setNames(rep(0, length(varnames)), varnames)
  obj <- NA_real_
  obj_line <- grep("^# Objective value =", lines, value = TRUE)
  if (length(obj_line) > 0L) {
    obj <- as.numeric(sub("^# Objective value =\\s*", "", obj_line[[1L]]))
  }
  value_lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  for (line in value_lines) {
    parts <- strsplit(trimws(line), "\\s+")[[1L]]
    if (length(parts) >= 2L && parts[[1L]] %in% varnames) {
      primal[[parts[[1L]]]] <- as.numeric(parts[[2L]])
    }
  }
  list(
    x = as.numeric(primal),
    objval = obj,
    status = NA_character_,
    objbound = NA_real_,
    mipgap = NA_real_,
    runtime = NA_real_
  )
}

.gp_gurobi_cli_solve <- function(problem,
                                 params,
                                 warm_start = NULL,
                                 encoding = c("indicator", "bigM")) {
  problem <- validate_gp_grade_problem(problem)
  encoding <- match.arg(encoding)
  warm_start <- .gp_gurobi_validate_start(problem, warm_start)
  cli <- unname(Sys.which("gurobi_cl"))
  if (!nzchar(cli)) {
    .gradepath_abort_backend_unavailable("`gurobi_cl` is not on PATH.")
  }

  dir <- tempfile("gradepath-gurobi-")
  dir.create(dir, recursive = TRUE)
  keep_files <- isTRUE(getOption("gradepath.keep_solver_files", FALSE))
  if (!keep_files) {
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  }

  lp_file <- file.path(dir, paste0("gradepath-", encoding, ".lp"))
  prm_file <- file.path(dir, "gradepath.prm")
  mst_file <- file.path(dir, "gradepath.mst")
  json_file <- file.path(dir, "gradepath.json")
  log_file <- file.path(dir, "gurobi_cl.log")

  if (identical(encoding, "indicator")) {
    .gp_gurobi_write_lp_indicator(problem, lp_file)
  } else {
    .gp_gurobi_write_lp_bigM(problem, lp_file)
  }
  # Keep an audit copy of the solver parameters, but pass them to gurobi_cl as
  # inline Name=value arguments so the actual invocation is visible in the
  # command line and matches the smoke-test path.
  .gp_gurobi_write_prm(params, prm_file)
  if (!is.null(warm_start)) {
    .gp_gurobi_write_mst(problem, warm_start, mst_file)
  }

  result_file <- json_file
  args <- c(
    .gp_gurobi_cli_param_args(params),
    "JSONSolDetail=1",
    paste0("ResultFile=", result_file),
    if (!is.null(warm_start)) paste0("InputFile=", mst_file) else character(),
    lp_file
  )
  process_timeout <- .gp_gurobi_cli_timeout(params)
  out <- .gp_gurobi_cli_run(cli, args = args, timeout = process_timeout)
  writeLines(as.character(out), log_file, useBytes = TRUE)
  status_code <- attr(out, "status")
  if (is.null(status_code)) {
    status_code <- 0L
  }
  if (isTRUE(attr(out, "timeout"))) {
    .gp_status_abort(
      "SOLVER_TIME_LIMIT",
      "`gurobi_cl` timed out after %s seconds on the %s model: %s",
      format(process_timeout),
      encoding,
      .gp_gurobi_cli_output_text(out),
      class = c("gradepath_gurobi_cli_timeout", "gradepath_backend_error")
    )
  }
  if (!identical(as.integer(status_code), 0L)) {
    .gradepath_abort(
      "`gurobi_cl` failed on the %s model: %s",
      encoding,
      .gp_gurobi_cli_output_text(out),
      class = c("gradepath_gurobi_cli_error", "gradepath_backend_error", "gradepath_error")
    )
  }
  if (!file.exists(json_file)) {
    .gradepath_abort(
      "`gurobi_cl` finished on the %s model but did not write the JSON solution file.",
      encoding,
      class = c("gradepath_gurobi_cli_missing_json", "gradepath_backend_error", "gradepath_error")
    )
  }

  parsed <- .gp_gurobi_parse_json(json_file, .gp_gurobi_varnames(problem))
  parsed$artifacts <- list(
    directory = if (keep_files) dir else NA_character_,
    model = basename(lp_file),
    params = basename(prm_file),
    params_style = "inline_cli_args_with_prm_audit_copy",
    start = if (!is.null(warm_start)) basename(mst_file) else NA_character_,
    solution = basename(result_file),
    log = basename(log_file)
  )
  parsed
}

.gp_gurobi_status <- function(solver_status, mipgap = NA_real_, params = list()) {
  solver_status <- .gp_gurobi_cli_status_label(solver_status)
  if (is.na(solver_status)) {
    return("suboptimal")
  }
  if (identical(solver_status, "OPTIMAL")) {
    # Gurobi may report OPTIMAL after satisfying a user-supplied MIPGap. Keep a
    # positive within-tolerance gap visible as `gap_reached`; acceptance-mode
    # fallback can then look for a strict zero-gap path.
    gap_limit <- params$MIPGap
    if (!is.null(gap_limit) &&
        is.finite(mipgap) &&
        mipgap > 1e-9 &&
        mipgap <= as.numeric(gap_limit) + 1e-12) {
      return("gap_reached")
    }
    return("optimal")
  }
  if (identical(solver_status, "TIME_LIMIT")) {
    return("time_limit")
  }
  if (solver_status %in% c("SUBOPTIMAL", "SOLUTION_LIMIT", "INTERRUPTED")) {
    return("suboptimal")
  }
  tolower(solver_status)
}

.gp_gurobi_scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L || anyNA(x[[1L]])) {
    NA_real_
  } else {
    as.numeric(x[[1L]])
  }
}

.gp_validate_solver_primal <- function(x) {
  tryCatch(
    .gradepath_validate_numeric_vector(x, "solver result x"),
    error = function(e) {
      if (inherits(e, "gradepath_error")) {
        .gradepath_abort_solver_output_invalid("%s", conditionMessage(e))
      }
      stop(e)
    }
  )
}

.gp_gurobi_normalize_result <- function(result,
                                        problem,
                                        params,
                                        path,
                                        encoding,
                                        warm_start = NULL,
                                        matrix_opt = NULL) {
  problem <- validate_gp_grade_problem(problem)
  path <- .gradepath_validate_scalar_character(path, "path", allowed = .gp_gurobi_paths)
  encoding <- .gradepath_validate_scalar_character(
    encoding,
    "encoding",
    allowed = .gp_grade_problem_encodings
  )

  solver_status <- .gp_gurobi_cli_status_label(result$status)
  if (solver_status %in% c("INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED")) {
    .gradepath_abort_solver_infeasible(
      "Gurobi reported an infeasible grade problem: %s.",
      solver_status
    )
  }

  primal <- .gp_validate_solver_primal(result$x)
  if (length(primal) != problem$constraint_ncol) {
    .gradepath_abort_solver_output_invalid(
      "Gurobi returned a primal vector with the wrong length."
    )
  }
  if (any(primal < problem$lower - 1e-6) ||
      any(primal > problem$upper + 1e-6)) {
    .gradepath_abort_solver_output_invalid(
      "Gurobi returned primal values outside problem bounds."
    )
  }
  objval <- .gp_gurobi_scalar_or_na(result$objval)
  if (!is.finite(objval)) {
    .gradepath_abort_solver_output_invalid(
      "Gurobi did not return a finite objective value."
    )
  }

  raw_from_primal <- sum(problem$objective * primal)
  if (!isTRUE(all.equal(raw_from_primal, objval, tolerance = 1e-6))) {
    .gradepath_abort_solver_objective_mismatch(
      "Gurobi objective mismatch: returned objval does not equal c'x."
    )
  }

  D <- tryCatch(
    .gradepath_grade_decision_matrix(problem, primal),
    error = function(e) {
      if (inherits(e, "gradepath_grade_error")) {
        .gradepath_abort_solver_output_invalid("%s", conditionMessage(e))
      }
      stop(e)
    }
  )
  assignment <- .gradepath_grade_assignment(D, problem$ids)
  canon <- 2 * objval
  if (!is.null(matrix_opt)) {
    canon_from_matrix <- .gradepath_grade_objective_value(
      matrix_opt,
      D,
      problem$lambda
    )
    if (!isTRUE(all.equal(canon_from_matrix, canon, tolerance = 1e-6))) {
      .gradepath_abort_solver_canonical_mismatch(
        "Canonical objective mismatch: grade objective is not 2 * Gurobi raw objval."
      )
    }
    canon <- canon_from_matrix
  }

  mipgap <- .gp_gurobi_scalar_or_na(result$mipgap)
  problem_hash <- gp_problem_hash(problem)
  problem_signature <- gp_problem_feasible_signature(problem)
  out <- list(
    primal = primal,
    objval = objval,
    status = .gp_gurobi_status(solver_status, mipgap, params),
    solver_meta = list(
      backend = "gurobi",
      path = path,
      encoding = encoding,
      problem_hash = problem_hash,
      problem_signature = problem_signature,
      solver_status = solver_status,
      objbound = .gp_gurobi_scalar_or_na(result$objbound),
      mipgap = mipgap,
      runtime = .gp_gurobi_scalar_or_na(result$runtime),
      params = params,
      timeout = identical(solver_status, "TIME_LIMIT"),
      warmstart = if (is.null(warm_start)) "none" else if (grepl("^gurobi_cl", path)) "mst_file" else "vector"
    ),
    validation = list(
      problem_hash = problem_hash,
      problem_signature = problem_signature,
      encoding = encoding,
      grades = assignment$grade,
      assignment = assignment,
      decision_matrix = D,
      raw_objective_from_primal = raw_from_primal,
      canon = canon,
      diagonal_zero = all(diag(D) == 0L)
    )
  )
  if (!is.null(result$artifacts)) {
    out$solver_meta$artifacts <- result$artifacts
  }
  structure(out, class = c("gp_gurobi_result", "list"))
}

.gp_open_status <- function(status, message = NULL) {
  status_text <- paste(
    toupper(as.character(c(status, message))),
    collapse = " "
  )
  if (grepl("OPTIMAL|GLP_OPT|TM_OPTIMAL", status_text)) {
    return("optimal")
  }
  if (grepl("TIME|TM_TIME|LIMIT", status_text)) {
    return("time_limit")
  }
  if (grepl("GAP", status_text)) {
    return("gap_reached")
  }
  if (grepl("FEASIBLE|SUBOPTIMAL|SOLUTION", status_text)) {
    return("suboptimal")
  }
  tolower(as.character(status[[1L]]))
}

.gp_validate_grade_matrix_opt <- function(problem, matrix_opt, name = "matrix_opt") {
  if (is.null(matrix_opt)) {
    return(NULL)
  }
  problem <- validate_gp_grade_problem(problem)
  matrix_opt <- .gradepath_validate_numeric_matrix(matrix_opt, name)
  if (!identical(dim(matrix_opt), c(problem$n_units, problem$n_units))) {
    .gradepath_abort_backend(
      "`%s` must be an n x n matrix aligned with the grade problem.",
      name
    )
  }
  if (any(abs(diag(matrix_opt)) > 1e-12)) {
    .gradepath_abort_backend(
      "`%s` must have projected zero diagonals before solving.",
      name
    )
  }
  matrix_opt
}

.gp_open_normalize_result <- function(result,
                                      problem,
                                      params,
                                      backend,
                                      path,
                                      encoding = "bigM",
                                      matrix_opt = NULL) {
  problem <- validate_gp_grade_problem(problem)
  backend <- .gradepath_validate_scalar_character(
    backend,
    "backend",
    allowed = c(.gp_open_backends, "alabama_dp")
  )
  path <- .gradepath_validate_scalar_character(path, "path")
  encoding <- .gradepath_validate_scalar_character(encoding, "encoding")

  solver_status <- as.character(result$status[[1L]])
  if (grepl("INFEAS|INF_OR_UNBD|UNBOUNDED", toupper(solver_status))) {
    .gradepath_abort_solver_infeasible(
      "Open backend `%s` reported an infeasible grade problem: %s.",
      backend,
      solver_status
    )
  }

  primal <- .gp_validate_solver_primal(result$x)
  if (length(primal) != problem$constraint_ncol) {
    .gradepath_abort_solver_output_invalid(
      "Open backend returned a primal vector with the wrong length."
    )
  }
  if (any(primal < problem$lower - 1e-6) ||
      any(primal > problem$upper + 1e-6)) {
    .gradepath_abort_solver_output_invalid(
      "Open backend returned primal values outside problem bounds."
    )
  }
  objval <- .gp_gurobi_scalar_or_na(result$objval)
  if (!is.finite(objval)) {
    .gradepath_abort_solver_output_invalid(
      "Open backend did not return a finite objective value."
    )
  }

  raw_from_primal <- sum(problem$objective * primal)
  if (!isTRUE(all.equal(raw_from_primal, objval, tolerance = 1e-6))) {
    .gradepath_abort_solver_objective_mismatch(
      "Open-backend objective mismatch: returned objval does not equal c'x."
    )
  }

  D <- tryCatch(
    .gradepath_grade_decision_matrix(problem, primal),
    error = function(e) {
      if (inherits(e, "gradepath_grade_error")) {
        .gradepath_abort_solver_output_invalid("%s", conditionMessage(e))
      }
      stop(e)
    }
  )
  assignment <- .gradepath_grade_assignment(D, problem$ids)
  canon <- 2 * objval
  if (!is.null(matrix_opt)) {
    canon_from_matrix <- .gradepath_grade_objective_value(
      matrix_opt,
      D,
      problem$lambda
    )
    if (!isTRUE(all.equal(canon_from_matrix, canon, tolerance = 1e-6))) {
      .gradepath_abort_solver_canonical_mismatch(
        "Canonical objective mismatch: grade objective is not 2 * raw objval."
      )
    }
    canon <- canon_from_matrix
  }

  solver_message <- if (is.null(result$message)) NA_character_ else as.character(result$message[[1L]])
  problem_hash <- gp_problem_hash(problem)
  problem_signature <- gp_problem_feasible_signature(problem)
  out <- list(
    primal = primal,
    objval = objval,
    status = .gp_open_status(solver_status, solver_message),
    solver_meta = list(
      backend = backend,
      path = path,
      encoding = encoding,
      problem_hash = problem_hash,
      problem_signature = problem_signature,
      solver_status = solver_status,
      solver_message = solver_message,
      objbound = .gp_gurobi_scalar_or_na(result$objbound),
      mipgap = .gp_gurobi_scalar_or_na(result$mipgap),
      runtime = .gp_gurobi_scalar_or_na(result$runtime),
      params = params,
      timeout = grepl("TIME|LIMIT", toupper(paste(solver_status, solver_message))),
      warmstart = "none"
    ),
    validation = list(
      problem_hash = problem_hash,
      problem_signature = problem_signature,
      encoding = encoding,
      grades = assignment$grade,
      assignment = assignment,
      decision_matrix = D,
      raw_objective_from_primal = raw_from_primal,
      canon = canon,
      diagonal_zero = all(diag(D) == 0L)
    )
  )
  structure(out, class = c("gp_open_backend_result", "list"))
}

gp_assert_same_problem <- function(ref, other, n = NULL, tol = 1e-9) {
  for (arg in c("ref", "other")) {
    value <- get(arg, inherits = FALSE)
    if (!is.list(value) || is.null(value$validation)) {
      .gradepath_abort_backend("`%s` must be a normalized solver result.", arg)
    }
    if (is.null(value$validation$problem_hash)) {
      .gradepath_abort_backend("`%s$validation$problem_hash` is missing.", arg)
    }
    if (is.null(value$validation$problem_signature)) {
      .gradepath_abort_backend("`%s$validation$problem_signature` is missing.", arg)
    }
    if (is.null(value$primal)) {
      .gradepath_abort_backend("`%s$primal` is missing.", arg)
    }
  }

  if (!identical(other$validation$problem_hash, ref$validation$problem_hash)) {
    .gradepath_abort_solver_canonical_mismatch(
      "Solver paths did not solve the same canonical grade problem."
    )
  }
  if (!identical(other$validation$problem_signature, ref$validation$problem_signature)) {
    .gradepath_abort_solver_canonical_mismatch(
      "Solver paths did not solve the same feasible grade problem."
    )
  }

  if (!is.null(n)) {
    n <- .gradepath_validate_integerish(n, "n", min = 1L)
    diag_pos <- gp_diag_positions(n)
    if (length(ref$primal) < max(diag_pos) ||
        length(other$primal) < max(diag_pos)) {
      .gradepath_abort_backend("Solver result primal vector is too short for `n`.")
    }
    if (any(abs(ref$primal[diag_pos]) > tol) ||
        any(abs(other$primal[diag_pos]) > tol)) {
      .gradepath_abort_solver_canonical_mismatch(
        "Solver paths did not preserve the projected zero diagonal."
      )
    }
  }

  invisible(TRUE)
}

gp_assert_same_answer <- function(ref, other, tol = 1e-9) {
  for (arg in c("ref", "other")) {
    value <- get(arg, inherits = FALSE)
    if (!is.list(value) || is.null(value$validation)) {
      .gradepath_abort_backend("`%s` must be a normalized solver result.", arg)
    }
    for (field in c("grades", "canon")) {
      if (is.null(value$validation[[field]])) {
        .gradepath_abort_backend(
          "`%s$validation$%s` is missing.",
          arg,
          field
        )
      }
    }
    if (is.null(value$objval)) {
      .gradepath_abort_backend("`%s$objval` is missing.", arg)
    }
  }

  if (!identical(as.integer(other$validation$grades),
                 as.integer(ref$validation$grades))) {
    .gradepath_abort_solver_canonical_mismatch("Solver paths returned different grade vectors.")
  }
  if (!isTRUE(all.equal(
    other$validation$canon,
    2 * other$objval,
    tolerance = tol
  ))) {
    .gradepath_abort_solver_canonical_mismatch(
      "Solver path does not satisfy the canonical 2x objective invariant."
    )
  }
  if (!isTRUE(all.equal(
    other$validation$canon,
    ref$validation$canon,
    tolerance = tol
  ))) {
    .gradepath_abort_solver_canonical_mismatch(
      "Solver paths returned different canonical objectives."
    )
  }

  invisible(TRUE)
}

.gp_gurobi_attempt <- function(path, status, reason = NA_character_) {
  data.frame(
    path = as.character(path),
    status = as.character(status),
    reason = as.character(reason),
    stringsAsFactors = FALSE
  )
}

.gp_gurobi_attempts_frame <- function(attempts) {
  if (length(attempts) == 0L) {
    return(data.frame(
      path = character(),
      status = character(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, attempts)
  row.names(out) <- NULL
  out
}

.gp_gurobi_attach_attempts <- function(result, attempts, acceptance_mode) {
  result$solver_meta$acceptance_mode <- isTRUE(acceptance_mode)
  result$solver_meta$attempts <- .gp_gurobi_attempts_frame(attempts)
  result
}

.gp_gurobi_attempt_status_from_error <- function(error) {
  if (inherits(error, "gp_status_solver_time_limit") ||
      inherits(error, "gradepath_gurobi_cli_timeout")) {
    return("time_limit")
  }
  if (inherits(error, "gp_status_solver_infeasible") ||
      inherits(error, "gradepath_solver_infeasible")) {
    return("infeasible")
  }
  "error"
}

# Robustness policy:
# - default mode is invocation fallback only. The first available path that
#   returns a normalized solver result is returned, even if its normalized status
#   is not acceptance-ready; downstream producer-status rules decide acceptance
#   readiness.
# - acceptance mode is solution-quality fallback. Non-acceptance-ready solver
#   results are recorded and later paths are attempted when available.
gp_gurobi_robust <- function(problem,
                             control = NULL,
                             warm_start = NULL,
                             matrix_opt = NULL,
                             force_path = NULL,
                             acceptance_mode = FALSE) {
  problem <- validate_gp_grade_problem(problem)
  control <- if (is.null(control)) gp_control() else validate_gp_control(control)
  acceptance_mode <- .gradepath_validate_scalar_logical(
    acceptance_mode,
    "acceptance_mode"
  )
  if (!identical(control$backend, "gurobi")) {
    .gradepath_abort_backend(
      "`control$backend` must be `gurobi` for `gp_gurobi_robust()`."
    )
  }
  warm_start <- .gp_gurobi_validate_start(problem, warm_start)
  matrix_opt <- .gp_validate_grade_matrix_opt(problem, matrix_opt)
  if (!is.null(force_path)) {
    # `force_path` accepts a character VECTOR (length >= 1) of valid Gurobi path
    # names, so a genuine multi-leg cascade can be forced (e.g. exercising a
    # two-path fall-through). A scalar still works because it is a length-1
    # vector. Every element must be a known Gurobi path.
    force_path <- .gradepath_validate_character_vector(
      force_path,
      "force_path"
    )
    if (!all(force_path %in% .gp_gurobi_paths)) {
      .gradepath_abort(
        "`force_path` must contain only Gurobi path names: %s.",
        paste(.gp_gurobi_paths, collapse = ", ")
      )
    }
  }

  params <- .gp_gurobi_params(control = control)
  detect <- .gp_gurobi_detect(smoke = TRUE)
  path_order <- if (is.null(force_path)) .gp_gurobi_paths else force_path
  failures <- list()
  attempts <- list()
  fallback_result <- NULL

  for (path in path_order) {
    if (!.gp_gurobi_path_available(path, detect)) {
      failures[[path]] <- "required Gurobi component(s) not available"
      attempts[[length(attempts) + 1L]] <- .gp_gurobi_attempt(
        path,
        "unavailable",
        failures[[path]]
      )
      next
    }

    result <- tryCatch(
      switch(
        path,
        in_process = .gp_gurobi_solve_inprocess(
          problem,
          params,
          warm_start = warm_start
        ),
        subprocess_R_binding = .gp_gurobi_solve_subprocess(
          problem,
          params,
          warm_start = warm_start
        ),
        gurobi_cl_indicator = .gp_gurobi_cli_solve(
          problem,
          params,
          warm_start = warm_start,
          encoding = "indicator"
        ),
        gurobi_cl_bigM = .gp_gurobi_cli_solve(
          problem,
          params,
          warm_start = warm_start,
          encoding = "bigM"
        )
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      failures[[path]] <- conditionMessage(result)
      attempts[[length(attempts) + 1L]] <- .gp_gurobi_attempt(
        path,
        .gp_gurobi_attempt_status_from_error(result),
        failures[[path]]
      )
      next
    }

    encoding <- if (identical(path, "gurobi_cl_bigM")) "bigM" else "indicator"
    normalized <- .gp_gurobi_normalize_result(
      result = result,
      problem = problem,
      params = params,
      path = path,
      encoding = encoding,
      warm_start = warm_start,
      matrix_opt = matrix_opt
    )
    attempt_reason <- NA_character_
    if (!.gp_status_acceptance_ready(normalized$status)) {
      attempt_reason <- if (isTRUE(acceptance_mode)) {
        "non_acceptance_ready_in_acceptance_mode"
      } else {
        "non_acceptance_ready_default_policy"
      }
    }
    attempts[[length(attempts) + 1L]] <- .gp_gurobi_attempt(
      path,
      normalized$status,
      attempt_reason
    )
    if (!isTRUE(acceptance_mode) ||
        .gp_status_acceptance_ready(normalized$status)) {
      return(.gp_gurobi_attach_attempts(
        normalized,
        attempts,
        acceptance_mode = acceptance_mode
      ))
    }
    if (is.null(fallback_result)) {
      fallback_result <- normalized
    }
  }

  if (!is.null(fallback_result)) {
    return(.gp_gurobi_attach_attempts(
      fallback_result,
      attempts,
      acceptance_mode = acceptance_mode
    ))
  }

  .gp_gurobi_abort_unavailable(detect, failures = failures)
}

.gp_roi_op <- function(problem) {
  problem <- validate_gp_grade_problem(problem)
  for (pkg in c("ROI", "slam")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      .gradepath_abort_backend_unavailable(
        "ROI open backends require the suggested package `%s`.",
        pkg
      )
    }
  }

  constraint_matrix <- slam::simple_triplet_matrix(
    i = problem$constraint_i,
    j = problem$constraint_j,
    v = problem$constraint_x,
    nrow = problem$constraint_nrow,
    ncol = problem$constraint_ncol
  )

  ROI::OP(
    objective = ROI::L_objective(problem$objective),
    constraints = ROI::L_constraint(
      L = constraint_matrix,
      dir = problem$direction,
      rhs = problem$rhs
    ),
    types = problem$types,
    maximum = FALSE,
    bounds = ROI::V_bound(
      li = seq_along(problem$lower),
      ui = seq_along(problem$upper),
      lb = problem$lower,
      ub = problem$upper
    )
  )
}

.gp_open_solution <- function(problem, backend, controls) {
  if (identical(backend, "highs")) {
    .gp_highs_solution(problem, controls)
  } else {
    .gp_roi_solution(problem, backend, controls)
  }
}

.gp_roi_result_status <- function(result) {
  status_code <- NA_integer_
  status_symbol <- NA_character_
  status_message <- NA_character_

  if (is.list(result$status)) {
    if (!is.null(result$status$code)) {
      status_code <- as.integer(result$status$code)
    }
    if (is.list(result$status$msg)) {
      if (!is.null(result$status$msg$symbol)) {
        status_symbol <- as.character(result$status$msg$symbol)
      }
      if (!is.null(result$status$msg$message)) {
        status_message <- as.character(result$status$msg$message)
      }
    }
  }

  list(
    code = status_code,
    symbol = status_symbol,
    message = status_message
  )
}

.gp_roi_solution <- function(problem, backend, controls) {
  backend <- .gradepath_validate_scalar_character(
    backend,
    "backend",
    allowed = c("glpk", "symphony")
  )
  plugin <- switch(
    backend,
    glpk = "ROI.plugin.glpk",
    symphony = "ROI.plugin.symphony"
  )
  if (!requireNamespace(plugin, quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "Backend `%s` requires the suggested package `%s`.",
      backend,
      plugin
    )
  }

  op <- .gp_roi_op(problem)
  result <- ROI::ROI_solve(op, solver = backend, control = controls)
  status <- .gp_roi_result_status(result)

  list(
    x = as.numeric(ROI::solution(result)),
    objval = as.numeric(ROI::solution(result, "objval")),
    status = if (is.na(status$symbol)) status$code else status$symbol,
    message = if (is.na(status$message)) paste("ROI code", status$code) else status$message,
    objbound = NA_real_,
    mipgap = NA_real_,
    runtime = NA_real_
  )
}

.gp_highs_types <- function(types) {
  ifelse(types == "B", "I", types)
}

.gp_highs_solution <- function(problem, controls) {
  problem <- validate_gp_grade_problem(problem)
  if (!requireNamespace("highs", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "Backend `highs` requires the suggested package `highs`."
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .gradepath_abort_backend_unavailable(
      "Backend `highs` requires the suggested package `Matrix`."
    )
  }

  A <- Matrix::sparseMatrix(
    i = problem$constraint_i,
    j = problem$constraint_j,
    x = problem$constraint_x,
    dims = c(problem$constraint_nrow, problem$constraint_ncol)
  )
  control <- do.call(highs::highs_control, controls)
  result <- highs::highs_solve(
    L = problem$objective,
    lower = problem$lower,
    upper = problem$upper,
    A = A,
    lhs = rep(-Inf, problem$constraint_nrow),
    rhs = problem$rhs,
    types = .gp_highs_types(problem$types),
    maximum = FALSE,
    control = control
  )

  info <- result$info
  list(
    x = as.numeric(result$primal_solution),
    objval = as.numeric(result$objective_value),
    status = result$status_message,
    message = result$status_message,
    objbound = if (!is.null(info$mip_dual_bound)) as.numeric(info$mip_dual_bound) else NA_real_,
    mipgap = if (!is.null(info$mip_gap)) as.numeric(info$mip_gap) else NA_real_,
    runtime = NA_real_
  )
}

gp_open_solve <- function(problem,
                          backend,
                          control = NULL,
                          matrix_opt = NULL) {
  problem <- validate_gp_grade_problem(problem)
  backend <- .gradepath_validate_scalar_character(
    backend,
    "backend",
    allowed = .gp_open_backends
  )
  control <- if (is.null(control)) gp_control(backend = backend) else validate_gp_control(control)
  if (!identical(control$backend, backend)) {
    .gradepath_abort_backend(
      "`control$backend` must match the requested open backend `%s`.",
      backend
    )
  }
  matrix_opt <- .gp_validate_grade_matrix_opt(problem, matrix_opt)
  if (!.gp_open_backend_available(backend)) {
    .gradepath_abort_backend_unavailable(
      "Requested open backend `%s` is unavailable; install: %s.",
      backend,
      paste(.gp_open_backend_packages(backend), collapse = ", ")
    )
  }

  controls <- .gp_open_solver_controls(control = control, backend = backend)
  raw <- .gp_open_solution(problem, backend, controls)

  .gp_open_normalize_result(
    result = raw,
    problem = problem,
    params = controls,
    backend = backend,
    path = if (identical(backend, "highs")) "highs" else paste0("roi_", backend),
    encoding = "bigM",
    matrix_opt = matrix_opt
  )
}

gp_alabama_dp <- function(problem, matrix_opt = NULL) {
  problem <- validate_gp_grade_problem(problem)
  if (abs(problem$lambda) > 1e-12) {
    .gradepath_abort_backend(
      "`gp_alabama_dp()` is an oracle only at lambda = 0; use a MILP backend for lambda > 0."
    )
  }
  matrix_opt <- .gp_validate_grade_matrix_opt(problem, matrix_opt)

  D <- matrix(
    0L,
    nrow = problem$n_units,
    ncol = problem$n_units,
    dimnames = list(problem$ids, problem$ids)
  )
  primal <- numeric(problem$constraint_ncol)
  primal[seq_len(problem$n_units * problem$n_units)] <- as.vector(t(D))
  primal[problem$grade_indices] <- 1

  .gp_open_normalize_result(
    result = list(
      x = primal,
      objval = 0,
      status = "ORACLE",
      message = "lambda = 0 pooled-grade oracle",
      objbound = 0,
      mipgap = 0,
      runtime = 0
    ),
    problem = problem,
    params = list(lambda = 0),
    backend = "alabama_dp",
    path = "alabama_dp_lambda0",
    encoding = "dp",
    matrix_opt = matrix_opt
  )
}
