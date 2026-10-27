# Helpers for testing that change statistics of effects match the difference
# in target statistics before and after a move. Based on target_change_match()
# in R/2x_effects_test.R, but with state, cache and resource covariates as
# arguments instead of the hard-coded lililiState / lililiCache.

# the arguments of an effect with the origin of the test data, for effects
# that depend on an origin
with_origin <- function(effect, args) {
  if ("origin" %in% names(formals(get(effect, mode = "function"))) && is.null(args$origin)) {
    args$origin <- "origin"
  }
  args
}

# sum of target contributions over all origin-destination pairs, or over the
# destinations for effects without an origin
target_stat <- function(effect, args, state, cache) {
  f <- get(effect, mode = "function")
  args <- with_origin(effect, args)
  destinations <- seq_along(state[[state[[state$dep.var]]$nodeSet[1]]]$ids)
  origins <- if (is.null(args$origin)) NA else
    seq_along(state[[state[[args$origin]]$nodeSet[1]]]$ids)
  total <- 0
  for (i in origins) {
    for (j in destinations) {
      total <- total + do.call(f, c(list(
        dep.var = 1, state = state, cache = cache,
        i = i, j = j, edge = 1, update = 1,
        getTargetContribution = TRUE
      ), args))
    }
  }
  total
}

# move resource `edge` from its current destination to `j_new`, the same way
# the sampler does in simulateOneStep(): remove with the old cache, update the
# cache, then add with the updated cache
check_move <- function(effect, args, edge, j_new,
                       state = lililiState, cache = lililiCache,
                       resourceCovariates = "sex") {
  f <- get(effect, mode = "function")
  args <- with_origin(effect, args)
  dv <- state$dep.var
  i <- if (is.null(args$origin)) NA else state[[args$origin]]$data[edge]
  j <- state[[dv]]$data[edge]

  change <- do.call(f, c(list(
    dep.var = 1, state = state, cache = cache,
    i = i, j = j, edge = edge, update = -1,
    getTargetContribution = FALSE
  ), args))

  updatedCache <- cache
  updatedCache[[1]] <- updateWeightedCache(cache[[1]], edge, j, state, update = -1)
  # as in simulateOneStep(), the moving resource's destination is NA while
  # the change of adding it is calculated
  removedState <- state
  removedState[[dv]]$data[edge] <- NA

  change <- change + do.call(f, c(list(
    dep.var = 1, state = removedState, cache = updatedCache,
    i = i, j = j_new, edge = edge, update = 1,
    getTargetContribution = FALSE
  ), args))

  newState <- state
  newState[[dv]]$data[edge] <- j_new
  newCache <- createInternalCache(newState, resourceCovariates = resourceCovariates)

  c(
    i = i, j = j,
    change = change,
    before = target_stat(effect, args, state, cache),
    after = target_stat(effect, args, newState, newCache)
  )
}

# moves to test: random ones plus moves out of loops, into loops and no-ops.
# With REALM_FULL_EFFECT_TESTS=true, every resource is moved once.
sample_moves <- function(n = 30, seed = 1, state = lililiState) {
  dest <- state[[state$dep.var]]$data
  orig <- state$origin$data
  nNodes <- length(state[[state[[state$dep.var]]$nodeSet[1]]]$ids)
  set.seed(seed)

  if (identical(Sys.getenv("REALM_FULL_EFFECT_TESTS"), "true")) {
    random <- seq_along(dest)
  } else {
    random <- sample(length(dest), n)
  }
  loopEdges <- which(orig == dest)[1:3]
  nonLoopEdges <- which(orig != dest)
  intoLoop <- nonLoopEdges[1:3]
  noOp <- nonLoopEdges[4:5]

  data.frame(
    edge = c(random, loopEdges, intoLoop, noOp),
    j_new = c(
      sample(nNodes, length(random), replace = TRUE),
      vapply(loopEdges, function(e) sample(setdiff(seq_len(nNodes), orig[e]), 1), numeric(1)),
      orig[intoLoop],
      dest[noOp]
    )
  )
}

expect_change_matches <- function(effect, args = list(), moves = sample_moves(),
                                  state = lililiState, cache = lililiCache,
                                  resourceCovariates = "sex") {
  for (k in seq_len(nrow(moves))) {
    r <- check_move(effect, args, moves$edge[k], moves$j_new[k],
                    state, cache, resourceCovariates)
    # tolerance relative to the size of the target statistic, as changes can
    # be tiny (e.g. AC effects with large weights) while targets are large
    diff <- r[["after"]] - r[["before"]]
    tol <- 1e-8 * max(1, abs(r[["before"]]))
    expect(
      is.finite(r[["change"]]) && abs(r[["change"]] - diff) <= tol,
      sprintf(
        "%s: moving edge %d (%s -> %d) to %d gives change statistic %.10g, but target statistics differ by %.10g",
        effect, moves$edge[k], format(r[["i"]]), r[["j"]], moves$j_new[k], r[["change"]], diff
      )
    )
  }
}
