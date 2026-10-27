# Helpers for testing the C++ sampler and effects against their R versions

# The C++ code sums in long double, as R's sum(), colSums() and cumsum() do, so
# C++ and R values are identical. R built without long double (e.g. CRAN's
# 'noLD' check) sums in double, which changes the last digits of R's values;
# there, the values are compared with a small tolerance instead.
expect_same <- function(object, expected, ...) {
  if (capabilities("long.double")) {
    expect_identical(object, expected, ...)
  } else {
    expect_equal(object, expected, tolerance = 1e-10, ...)
  }
}

# the C++ and R samplers give the same chain for the same seed: same
# destinations, cache and random number state afterwards
compare_chains <- function(effects, parameters, n = 5000,
                           multinomial = FALSE, state = myState, seed = 42) {
  dv <- state$dep.var
  oldOptions <- options(realm.useCpp = FALSE)
  on.exit(options(oldOptions), add = TRUE)
  set.seed(seed)
  r <- simulateNSteps(dv, state, state$cache, effects, parameters,
                      multinomialProposal = multinomial, n = n)
  rNext <- runif(1)
  options(realm.useCpp = TRUE)
  set.seed(seed)
  cpp <- simulateChain(dv, state, state$cache, effects, parameters,
                       burnIn = n, iterations = 0, thinning = 0,
                       multinomialProposal = multinomial)
  cppNext <- runif(1)

  expect_identical(cpp$state[[dv]]$data, toListState(r$state, dv)[[dv]]$data)
  expect_equal(cpp$cache, toListCache(r$cache), tolerance = 0)
  expect_identical(cppNext, rNext)
}

# every effect of `effects` gives exactly the same change statistic and target
# contribution in C++ as in R (up to rounding without long double, see above),
# for random moves including moves into and out of loops (of the first origin)
expect_cpp_matches_r <- function(effects, state = myState, n = 100, seed = 3) {
  dv <- state$dep.var
  dest <- state[[dv]]$data
  orig <- state[[state$origins[1]]]$data
  nNodes <- length(state[[state[[dv]]$nodeSet[1]]]$ids)
  set.seed(seed)
  E <- sample(length(dest), n, TRUE)
  J <- sample(nNodes, n, TRUE)
  loopMoves <- seq_len(max(1, n %/% 10))
  J[loopMoves] <- orig[E[loopMoves]]
  fromLoops <- which(orig == dest)
  if (length(fromLoops) > 0) {
    E[loopMoves + length(loopMoves)] <- fromLoops[seq_along(loopMoves)]
  }
  U <- sample(c(-1, 1), n, TRUE)
  origins <- effectOrigins(effects)

  for (k in seq_along(effects$name)) {
    f <- effects$effectFormulas[[k]]
    for (m in 1:n) {
      i <- if (is.na(origins[k])) NA else state[[origins[k]]]$data[E[m]]
      expected <- c(
        change = f(dep.var = dv, state = state, cache = state$cache,
                   i = i, j = J[m], edge = E[m], update = U[m]),
        target = f(dep.var = dv, state = state, cache = state$cache,
                   i = i, j = J[m], edge = E[m], update = U[m],
                   getTargetContribution = TRUE)
      )
      expect_same(
        cppEffectValues(effects, k, state, state$cache, i, J[m], E[m], U[m]),
        expected,
        label = paste(effects$name[k], "edge", E[m], "to", J[m], "update", U[m])
      )
    }
  }
  # whole target statistics, which some C++ effects calculate differently
  old <- options(realm.useCpp = FALSE)
  rStatistics <- getNetworkStatistics(dv, state, state$cache, effects)
  options(old)
  expect_same(cppNetworkStatistics(dv, state, state$cache, effects), rStatistics)
}

# the example state with a continuous dyadic covariate (similarity of
# organisation sizes) and a 'same' matrix for region; in the global environment, as addEffect() finds the
# state by its name there
local({
  cppTestState <- myState
  cppTestState$size_similarity <- createNetwork(myState$size$sim,
                                                nodeSet = c("organisations", "organisations"))
  # region with the 'same' matrix used by same_covariate
  cppTestState$region <- createNodeVariable(myState$region$data, nodeSet = "organisations",
                                            addSame = TRUE)
  assign("cppTestState", cppTestState, envir = globalenv())
})
