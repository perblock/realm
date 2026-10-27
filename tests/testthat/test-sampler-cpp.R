# The C++ sampler (src/sampler.cpp) must give exactly the same chain as the R
# sampler for the same seed: same destinations, cache, statistics and random
# number state afterwards.

test_that("C++ and R chains are identical with the loops effect", {
  effects <- createEffects(myState) |> addEffect(loops, origin = "origin")
  for (theta in c(0, 1.5, -1)) {
    compare_chains(effects, theta)
  }
})

test_that("C++ statistics equal the R statistics of each draw", {
  effects <- createEffects(myState) |> addEffect(loops, origin = "origin")
  dv <- myState$dep.var
  old <- options(realm.useCpp = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(5)
  r <- simulateStatisticVectors(dv, myState, myState$cache, effects, 1,
                                burnIn = 200, iterations = 20, thinning = 200)
  options(realm.useCpp = TRUE)
  set.seed(5)
  cpp <- simulateChain(dv, myState, myState$cache, effects, 1,
                       burnIn = 200, iterations = 20, thinning = 200)
  expect_same(unname(cpp$statistics), unname(r))
})

test_that("models with effects written only in R use the R sampler", {
  # a user-written effect without a C++ version; addEffect() finds effects by
  # name in the global environment
  assign("myLoops", function(dep.var = 1, origin, state, cache, i, j, edge, update,
                             getTargetContribution = FALSE) {
    if (getTargetContribution) return(cache[[dep.var]][[origin]]$valuedNetwork[i, j] * (i == j))
    if (i == j) return(update)
    0
  }, envir = globalenv())
  on.exit(rm("myLoops", envir = globalenv()), add = TRUE)

  effects <- createEffects(myState) |> addEffect(loops, origin = "origin") |> addEffect(myLoops, origin = "origin")
  expect_null(cppEffectSpec(effects, myState))
  expect_false(useCppSampler(effects, myState))

  set.seed(1)
  res <- simulateStatisticVectors(myState$dep.var, myState, myState$cache, effects,
                                  c(0.5, 0.5), burnIn = 100, iterations = 3, thinning = 100)
  # both effects count loops
  expect_equal(res[, 1], res[, 2])
})

test_that("C++ effects give the same values as the R effects", {
  effects <- createEffects(myState) |>
    addEffect(loops, origin = "origin") |>
    addEffect(concentration_AC, origin = "origin", alpha = 4) |>
    addEffect(in_weights_AC, alpha = 2) |>
    addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion") |>
    addEffect(alter_covariate, node.attribute = "size") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "region", edge.attribute = "sex") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size",
              edge.attribute = "second_or") |>
    addEffect(loops_resource_covar, origin = "origin", edge.attribute = "second_or")
  expect_cpp_matches_r(effects)
})

test_that("C++ and R chains are identical for the example model", {
  compare_chains(myEffects, myResDN$estimates, n = 5000)
  compare_chains(myEffects, rep(0, 6), n = 5000)
})

test_that("C++ and R chains are identical with the multinomial proposal", {
  compare_chains(myEffects, myResDN$estimates, n = 1000, multinomial = TRUE)
  compare_chains(myEffects, rep(0, 6), n = 1000, multinomial = TRUE)
})

test_that("estimation gives the same results with the C++ and the R sampler", {
  skip_on_cran()
  alg <- createAlgorithm(myState, myEffects, burnInN1 = 300, thinningN1 = 100,
                         burnInN2 = 300, thinningN2 = 100, nsubN2 = 1,
                         initialIterationsN2 = 15, burnInN3 = 300, thinningN3 = 100,
                         iterationsN3 = 30)
  run <- function(useCpp) {
    old <- options(realm.useCpp = useCpp)
    on.exit(options(old))
    set.seed(8)
    estimateMobilityNetwork(myState, myEffects, alg)
  }
  expect_true(useCppSampler(myEffects, myState))
  resR <- run(FALSE)
  resCpp <- run(TRUE)
  expect_equal(resCpp$estimates, resR$estimates, tolerance = 1e-12)
  expect_equal(resCpp$covarianceMatrix, resR$covarianceMatrix, tolerance = 1e-10)
  expect_equal(resCpp$convergenceStatistics, resR$convergenceStatistics, tolerance = 1e-10)
})

test_that("the multinomial proposal works with a single effect", {
  effects <- createEffects(myState) |> addEffect(loops, origin = "origin")
  compare_chains(effects, 1, n = 500, multinomial = TRUE)
})

test_that("simulation needs one parameter per effect", {
  expect_error(
    simulateMobilityNetworks(myState, myEffects, parameters = rep(0.1, 5),
                             burnin = 10, thinning = 10,
                             nSimulations = 1),
    "one value per effect"
  )
})
