# The sampler updates state and cache in place (see toMutableState and
# toMutableCache in R/1_hiddenFunctions.R). These tests check that the cache
# stays consistent with the destinations, that the user's objects are never
# changed, and that stored draws are independent copies.

inplace_effects <- createEffects(myState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(concentration_AC, origin = "origin", alpha = 4) |>
  addEffect(in_weights_AC, alpha = 4) |>
  addEffect(loops_resource_covar, origin = "origin", edge.attribute = "sex") |>
  addEffect(avoiding_dissimilar_covar_cont, origin = "origin", edge.attribute = "second_or")

test_that("cache matches a rebuilt cache after many in-place steps", {
  state <- lililiState
  cache <- createInternalCache(state, resourceCovariates = c("sex", "second_or"))
  effects <- inplace_effects
  dv <- state$dep.var

  for (multinomial in c(FALSE, TRUE)) {
    set.seed(11)
    res <- simulateNSteps(dv, state, cache, effects, c(0.5, 0.4, 0.3, 0.2, 0.05),
                          multinomialProposal = multinomial,
                          n = if (multinomial) 300 else 3000)
    simState <- toListState(res$state, dv)
    simCache <- toListCache(res$cache)

    expect_s3_class(simState[[dv]], "dependent.realm")
    expect_false(identical(simState[[dv]]$data, state[[dv]]$data))
    expect_equal(
      simCache,
      createInternalCache(simState, resourceCovariates = c("sex", "second_or")),
      tolerance = 1e-12
    )
  }
})

test_that("in-place update and restore match the list version", {
  state <- lililiState
  cache <- createInternalCache(state, resourceCovariates = c("sex", "second_or"))
  dv <- state$dep.var
  mutable <- toMutableCache(cache)

  for (edge in c(1, 50, 300)) {
    j <- state[[dv]]$data[edge]
    expected <- updateWeightedCache(cache[[dv]], edge, j, state = state, update = -1)
    cells <- getCacheCells(mutable[[dv]], edge, j, state)
    updateCacheInPlace(mutable[[dv]], edge, j, state = state, update = -1)
    expect_identical(toListCache(mutable)[[dv]], expected)

    restoreCacheCells(mutable[[dv]], edge, j, state, cells)
    expect_identical(toListCache(mutable)[[dv]], cache[[dv]])
  }
})

test_that("simulation does not change the user's state or cache", {
  state <- myState
  frozen <- unserialize(serialize(state, NULL))

  set.seed(3)
  simulateNSteps(state$dep.var, state, state$cache, myEffects, rep(0.2, 6),
                 n = 1000)
  expect_identical(state, frozen)

  invisible(capture.output(
    simulateMobilityNetworks(state, myEffects, rep(0.2, 6),
                             burnin = 200, thinning = 200, nSimulations = 2)
  ))
  expect_identical(state, frozen)
})

test_that("stored simulations are independent lists", {
  set.seed(4)
  invisible(capture.output(
    sims <- simulateMobilityNetworks(myState, myEffects, rep(0.2, 6),
                                     burnin = 200, thinning = 300, nSimulations = 3)
  ))
  dv <- myState$dep.var

  for (s in sims) {
    expect_type(s$cache[[dv]], "list")
    expect_s3_class(s$state[[dv]], "dependent.realm")
    expect_equal(s$cache, createInternalCache(
      s$state, resourceCovariates = names(myState$cache[[dv]]$destinations$resourceInflows)
    ))
  }
  expect_false(identical(sims[[1]]$state[[dv]]$data, sims[[2]]$state[[dv]]$data))
  expect_false(identical(sims[[2]]$state[[dv]]$data, sims[[3]]$state[[dv]]$data))
})

test_that("draws returned by simulateStatisticVectors are independent lists", {
  dv <- myState$dep.var
  set.seed(5)
  res <- simulateStatisticVectors(dv, myState, myState$cache, myEffects, rep(0.2, 6),
                                  burnIn = 100, iterations = 3, thinning = 300,
                                  returnDeps = TRUE)
  deps <- res[[2]]
  expect_length(deps, 3)
  for (k in seq_along(deps)) {
    expect_type(deps[[k]]$cache[[dv]], "list")
    expect_equal(
      unname(getNetworkStatistics(dv, deps[[k]]$state, deps[[k]]$cache, myEffects)),
      unname(res[[1]][k, ])
    )
  }
  expect_false(identical(deps[[1]]$state[[dv]]$data, deps[[3]]$state[[dv]]$data))
})
