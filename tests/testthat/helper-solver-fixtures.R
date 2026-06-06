step45_pairwise <- function(P, backend = "highs") {
  ids <- rownames(P)
  if (is.null(ids)) {
    ids <- paste0("u", seq_len(nrow(P)))
    dimnames(P) <- list(ids, ids)
  }
  new_gp_pairwise(
    ids = ids,
    matrix = P,
    power = 0L,
    cleanup = list(antisymmetry = TRUE, diagonal = 0.5, zero_floor = 1e-7),
    source = list(
      stage = "posterior",
      rule = "outer_product",
      assumption = "one_level_independence"
    ),
    control = gp_control(backend = backend),
    provenance = .gradepath_new_provenance(producer = "step45_fixture")
  )
}

step45_matrix_opt <- function(P) {
  out <- P
  diag(out) <- 0
  storage.mode(out) <- "double"
  out
}

step45_problem <- function(P, lambda = 1, backend = "highs") {
  pairwise <- step45_pairwise(P, backend = backend)
  list(
    P = P,
    matrix_opt = step45_matrix_opt(P),
    pairwise = pairwise,
    problem = gp_grade_problem(pairwise, lambda = lambda, encoding = "bigM")
  )
}

step45_strict3_matrix <- function() {
  ids <- paste0("u", 1:3)
  P <- matrix(
    c(
      0.50, 0.95, 0.97,
      0.05, 0.50, 0.92,
      0.03, 0.08, 0.50
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  P
}

step45_n9_matrix <- function() {
  ids <- paste0("u", 1:4)
  P <- matrix(
    c(
      0.50, 0.95, 0.97, 0.99,
      0.05, 0.50, 0.95, 0.97,
      0.03, 0.05, 0.50, 0.95,
      0.01, 0.03, 0.05, 0.50
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  P
}

step45_solver_candidates <- function() {
  candidates <- list()
  detect <- tryCatch(.gp_gurobi_detect(smoke = TRUE), error = function(e) NULL)
  if (!is.null(detect)) {
    for (path in c("in_process", "subprocess_R_binding")) {
      if (isTRUE(.gp_gurobi_path_available(path, detect))) {
        candidates[[paste0("gurobi:", path)]] <- list(
          family = "gurobi",
          path = path
        )
      }
    }
    if (isTRUE(.gp_gurobi_cl_smoke_ok())) {
      for (path in c("gurobi_cl_indicator", "gurobi_cl_bigM")) {
        if (isTRUE(.gp_gurobi_path_available(path, detect))) {
          candidates[[paste0("gurobi:", path)]] <- list(
            family = "gurobi",
            path = path
          )
        }
      }
    }
  }

  for (backend in .gp_available_open_backends()) {
    candidates[[paste0("open:", backend)]] <- list(
      family = "open",
      backend = backend
    )
  }

  candidates
}

step45_solver_results <- function(problem, matrix_opt) {
  candidates <- step45_solver_candidates()
  results <- list()
  failures <- list()

  for (name in names(candidates)) {
    candidate <- candidates[[name]]
    result <- tryCatch(
      {
        if (identical(candidate$family, "gurobi")) {
          gp_gurobi_robust(
            problem,
            control = gp_control(
              backend = "gurobi",
              solver_options = list(time_limit = 20)
            ),
            force_path = candidate$path,
            matrix_opt = matrix_opt
          )
        } else {
          gp_open_solve(
            problem,
            backend = candidate$backend,
            control = gp_control(
              backend = candidate$backend,
              solver_options = list(time_limit = 20)
            ),
            matrix_opt = matrix_opt
          )
        }
      },
      error = function(e) e
    )

    if (inherits(result, "error")) {
      failures[[name]] <- conditionMessage(result)
    } else {
      results[[name]] <- result
    }
  }

  attr(results, "candidates") <- names(candidates)
  attr(results, "failures") <- failures
  results
}

step45_skip_without_solver <- function(results, min_results = 1L) {
  failures <- attr(results, "failures", exact = TRUE)
  if (length(failures) > 0L) {
    testthat::fail(paste(
      "available solver path(s) failed:",
      paste(sprintf("%s: %s", names(failures), unlist(failures)), collapse = "; ")
    ))
  }
  testthat::skip_if(
    length(results) < min_results,
    sprintf("need at least %d available solver result(s)", min_results)
  )
}
