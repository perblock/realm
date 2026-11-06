# Starting values of the estimation, and models of only exogenous effects,
# for which the MPLE is the maximum likelihood estimate.

quickAlgorithm <- function(state, effects, ...) {
  createAlgorithm(state, effects, multinomialProposal = FALSE, burnInN1 = 200, thinningN1 = 50,
                  iterationsN1 = 40, burnInN2 = 200, thinningN2 = 50, nsubN2 = 1,
                  initialIterationsN2 = 10, burnInN3 = 200, thinningN3 = 50,
                  iterationsN3 = 20, ...)
}

# the exogenous model of the example data, a model with endogenous effects, and
# an exogenous model with a second origin of another nodeset, origin dummies
# and fixed effects
exoEffects <- createEffects(myState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion") |>
  addEffect(alter_covariate, node.attribute = "size") |>
  addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex") |>
  addEffect(loops_resource_covar, origin = "origin", edge.attribute = "sex")

local({
  regionOf <- orgRegion + 1
  regions <- createNodeSet(2)
  originRegion <- createOrigin(regionOf[myState$origin$data], nodeSet = c("regions", "people"))
  sameRegionRegion <- createNetwork(outer(1:2, regionOf, "==") * 1,
                                    nodeSet = c("regions", "organisations"))
  for (name in c("transfers", "origin", "people", "organisations", "size", "sex")) {
    assign(name, myState[[name]])
  }
  fixedState <- realmDataCreate(transfers, origin, originRegion, people, organisations, regions,
                                size, sex, sameRegionRegion, fixedEffectDummies = TRUE)
  assign("fixedState", fixedState, envir = globalenv())
})

fixedEffects <- createEffects(fixedState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(dyadic_covariate, origin = "originRegion", node.attribute = "sameRegionRegion") |>
  addEffect(origin_to_node_covar, origin = "originRegion", node.attribute = "size") |>
  addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex") |>
  addFixedEffects(fixedState)

mixedEffects <- exoEffects |> addEffect(concentration_AC, origin = "origin", alpha = 4)

test_that("exogenousEffects are exactly the effects whose statistics add up over resources", {
  skip_on_cran()
  ns <- asNamespace("realm")
  state <- myState
  state$sizeSame <- createNodeVariable(myState$size$data, nodeSet = "organisations", addSame = TRUE)
  dv <- state$dep.var
  y0 <- state[[dv]]
  origins <- state$origin$data
  sex <- state$sex$data
  destinations <- state$organisations$ids
  n <- length(y0$data)

  # with the state as simulated in round r, moves two resources that interact
  # if the effect has interactions: with the same origin or the same sex, to
  # the same destination, in half of the cases the origin of the first
  moves <- lapply(1:100, function(r) {
    set.seed(r)
    other <- sample(n, 60)
    first <- sample(n, 1)
    candidates <- setdiff(if (r %% 3 == 0) which(sex == sex[first]) else
      which(origins == origins[first]), first)
    if (length(candidates) == 0) candidates <- setdiff(seq_len(n), first)
    list(other = other, otherTo = sample(destinations, 60, TRUE),
         resources = c(first, candidates[sample.int(length(candidates), 1)]),
         to = if (r %% 2 == 0) origins[first] else sample(destinations, 1))
  })

  specification <- function(effect) {
    formalNames <- names(formals(get(effect, envir = ns)))
    args <- list(effect)
    if ("origin" %in% formalNames) args$origin <- "origin"
    if ("attribute.index" %in% formalNames) {
      args$attribute.index <- if (attributeRoles[[effect]] == "dyad") "sameRegion" else
        if (effect == "same_covariate") "sizeSame" else "size"
    }
    if ("resource.attribute.index" %in% formalNames) {
      args$resource.attribute.index <- if (grepl("_cont|additional_origin", effect)) "second_or" else "sex"
    }
    if ("dummy.value" %in% formalNames) args$dummy.value <- 2
    args
  }

  additive <- function(effect) {
    effects <- suppressWarnings(createEffectsObject(list(specification(effect)),
                                                    checkProcessState = state))
    statistic <- function(data) {
      s <- state
      s[[dv]]$data <- data
      cache <- createInternalCache(s, resourceCovariates = c("sex", "second_or"))
      as.numeric(getNetworkStatistics(dv, s, cache, effects))
    }
    for (m in moves) {
      y <- y0$data
      y[m$other] <- m$otherTo
      a <- y; a[m$resources[1]] <- m$to
      b <- y; b[m$resources[2]] <- m$to
      ab <- a; ab[m$resources[2]] <- m$to
      if (abs(statistic(ab) - statistic(b) - statistic(a) + statistic(y)) > 1e-9) {
        return(FALSE)
      }
    }
    TRUE
  }

  effectNames <- Filter(function(name) {
    f <- get(name, envir = ns)
    is.function(f) && "getTargetContribution" %in% names(formals(f))
  }, getNamespaceExports(ns))
  expect_setequal(Filter(additive, effectNames), exogenousEffects)
})

test_that("effects are classified by their name; unknown effects are endogenous", {
  expect_equal(isEndogenousEffect(mixedEffects), c(rep(FALSE, 5), TRUE))
  expect_false(any(isEndogenousEffect(fixedEffects)))
  expect_true(isEndogenousEffect(list(name = "my_own_effect origin")))
})

test_that("getStartingValues gives the values of each method", {
  dv <- myState$dep.var
  cache <- myState$cache
  frozen <- unserialize(serialize(myState, NULL))

  zero <- getStartingValues("zero", dv, myState, cache, mixedEffects)
  expect_equal(zero, setNames(rep(0, 6), mixedEffects$name))

  expect_equal(getStartingValues("full MPLE", dv, myState, cache, mixedEffects),
               getMPLE(dv, myState, cache, mixedEffects))

  partial <- getStartingValues("partial MPLE", dv, myState, cache, mixedEffects)
  expect_named(partial, mixedEffects$name)
  expect_equal(unname(partial[6]), 0)
  expect_equal(partial[1:5], getMPLE(dv, myState, cache, exoEffects))
  expect_false(isTRUE(all.equal(partial[1:5], getMPLE(dv, myState, cache, mixedEffects)[1:5])))

  # without endogenous effects, and with only endogenous effects
  expect_equal(getStartingValues("partial MPLE", dv, myState, cache, exoEffects),
               getMPLE(dv, myState, cache, exoEffects))
  onlyEndogenous <- createEffects(myState) |> addEffect(concentration_AC, origin = "origin", alpha = 4)
  expect_equal(unname(getStartingValues("partial MPLE", dv, myState, cache, onlyEndogenous)), 0)

  expect_identical(myState, frozen)
})

test_that("starting values are zero if the MPLE does not exist", {
  local_mocked_bindings(getMPLE = function(...) NULL)
  for (method in c("partial MPLE", "full MPLE")) {
    expect_output(values <- getStartingValues(method, myState$dep.var, myState, myState$cache,
                                              mixedEffects),
                  "not available; starting from zero")
    expect_equal(unname(values), rep(0, 6))
  }
})

test_that("an estimation uses the starting values of startingValues", {
  dv <- myState$dep.var
  for (method in c("RM", "GT")) {
    alg <- if (method == "RM") quickAlgorithm(myState, mixedEffects) else
      createAlgorithm(myState, mixedEffects, phase2Method = "GT", burnInN2 = 200, thinningN2 = 50,
                      iterationsGT = 30, finalIterationsGT = 60, maxIterationsGT = 2,
                      burnInN3 = 200, thinningN3 = 50, iterationsN3 = 40, multinomialProposal = FALSE)
    estimate <- function(...) {
      set.seed(11)
      estimateMobilityNetwork(myState, mixedEffects, alg, ...)$estimates
    }
    partial <- getStartingValues("partial MPLE", dv, myState, myState$cache, mixedEffects)
    full <- getMPLE(dv, myState, myState$cache, mixedEffects)

    expect_identical(estimate(), estimate(startingValues = "partial MPLE"))
    expect_identical(estimate(), estimate(initialParameters = partial))
    expect_identical(estimate(startingValues = "full MPLE"), estimate(initialParameters = full))
    expect_identical(estimate(startingValues = "zero"), estimate(initialParameters = rep(0, 6)))
    expect_false(identical(estimate(), estimate(startingValues = "full MPLE")))
    expect_false(identical(estimate(), estimate(startingValues = "zero")))
    # initialParameters come before startingValues
    expect_identical(estimate(startingValues = "zero", initialParameters = full),
                     estimate(initialParameters = full))
  }

  expect_error(estimateMobilityNetwork(myState, mixedEffects, quickAlgorithm(myState, mixedEffects),
                                       startingValues = "MPLE"),
               "should be one of")
})

test_that("prevAns comes before startingValues", {
  alg <- quickAlgorithm(myState, mixedEffects)
  set.seed(5)
  previous <- estimateMobilityNetwork(myState, mixedEffects, alg)
  estimate <- function(...) {
    set.seed(6)
    capture.output(res <- estimateMobilityNetwork(myState, mixedEffects, alg, prevAns = previous, ...))
    res$estimates
  }
  expect_identical(estimate(), estimate(startingValues = "zero"))
})

# the maximum likelihood estimate of a multinomial logit, from a design built
# from the raw data without realm's change statistics: the Poisson model of the
# counts of each destination for each group of resources with the same
# predictors. design(g, j) gives the predictors of destination j for the first
# resource of group g.
independentFit <- function(destinations, groups, nDestinations, design) {
  rows <- do.call(rbind, lapply(unique(groups), function(g) {
    first <- which(groups == g)[1]
    x <- t(vapply(seq_len(nDestinations), function(j) design(first, j), numeric(length(design(first, 1)))))
    data.frame(group = g, count = tabulate(destinations[groups == g], nDestinations), x)
  }))
  predictors <- setdiff(names(rows), c("group", "count"))
  model <- glm(reformulate(c("0", "factor(group)", predictors), "count"),
               family = poisson, data = rows, control = list(epsilon = 1e-12, maxit = 100))
  coefficients <- summary(model)$coefficients[predictors, ]
  list(estimates = unname(coefficients[, 1]), standardErrors = unname(coefficients[, 2]))
}

exoFit <- local({
  o <- myState$origin$data
  sex <- myState$sex$data
  size <- myState$size$data
  sameRegion <- myState$sameRegion$data
  independentFit(myState$transfers$data, paste(o, sex), length(size), function(e, j) {
    c(loops = o[e] == j, dyadic = sameRegion[o[e], j], alter = size[j],
      sexSize = sex[e] * size[j], sexLoops = (o[e] == j) * sex[e])
  })
})

test_that("the default returns the maximum likelihood estimates of an exogenous model", {
  alg <- quickAlgorithm(myState, exoEffects)
  expect_output(res <- estimateMobilityNetwork(myState, exoEffects, alg), "returning the maximum")

  expect_s3_class(res, "result.realm")
  expect_named(res$estimates, exoEffects$name)
  expect_equal(unname(res$estimates), exoFit$estimates, tolerance = 1e-6)
  expect_equal(unname(res$standardErrors), exoFit$standardErrors, tolerance = 1e-6)
  expect_equal(res$covarianceMatrix, t(res$covarianceMatrix))
  expect_true(all(eigen(res$covarianceMatrix)$values > 0))
  expect_equal(diag(res$covarianceMatrix), res$standardErrors^2, ignore_attr = TRUE)
  expect_lt(max(abs(res$convergenceStatistics)), 1e-6)
  expect_identical(res$algorithm, alg)
  expect_null(res$deps)
  expect_output(print(res), "StandardErrors|Standard")
})

test_that("the R pseudo-likelihood gives the same results for an exogenous model", {
  old <- options(realm.useCpp = FALSE)
  on.exit(options(old))
  capture.output(res <- estimateMobilityNetwork(myState, exoEffects, quickAlgorithm(myState, exoEffects)))
  expect_equal(unname(res$estimates), exoFit$estimates, tolerance = 1e-6)
  expect_equal(unname(res$standardErrors), exoFit$standardErrors, tolerance = 1e-6)
})

test_that("an exogenous model with two origins, origin dummies and fixed effects is correct", {
  o <- fixedState$origin$data
  r <- fixedState$originRegion$data
  sex <- fixedState$sex$data
  size <- fixedState$size$data
  sameRegion <- fixedState$sameRegionRegion$data
  dummies <- as.integer(sub(".*dummy", "", grep("dummy", fixedEffects$name, value = TRUE)))
  fit <- independentFit(fixedState$transfers$data, paste(o, sex), length(size), function(e, j) {
    c(loops = o[e] == j, dyadic = sameRegion[r[e], j], originSize = (r[e] == 2) * size[j],
      sexSize = sex[e] * size[j], d = (j == dummies) * 1)
  })
  capture.output(res <- estimateMobilityNetwork(fixedState, fixedEffects,
                                                quickAlgorithm(fixedState, fixedEffects)))
  expect_equal(unname(res$estimates), fit$estimates, tolerance = 1e-5)
  expect_equal(unname(res$standardErrors), fit$standardErrors, tolerance = 1e-5)
})

test_that("a simulation-based estimation agrees with the MPLE of an exogenous model", {
  skip_on_cran()
  capture.output(mple <- estimateMobilityNetwork(myState, exoEffects, quickAlgorithm(myState, exoEffects)))
  alg <- createAlgorithm(myState, exoEffects, multinomialProposal = FALSE, iterationsN3 = 2000)
  set.seed(21)
  simulated <- estimateMobilityNetwork(myState, exoEffects, alg, startingValues = "zero")
  expect_lt(max(abs(simulated$estimates - mple$estimates) / mple$standardErrors), 0.25)
  expect_lt(max(abs(simulated$standardErrors / mple$standardErrors - 1)), 0.1)
})

test_that("phases are skipped only with the default for exogenous models", {
  alg <- quickAlgorithm(myState, exoEffects)
  ran <- function(...) stop("a phase was run")

  local({
    local_mocked_bindings(runPhase1 = ran, runPhase2 = ran, runPhase3 = ran)
    expect_no_error(capture.output(estimateMobilityNetwork(myState, exoEffects, alg)))
    expect_error(estimateMobilityNetwork(myState, exoEffects, alg, startingValues = "full MPLE"),
                 "a phase was run")
    expect_error(estimateMobilityNetwork(myState, exoEffects, alg, startingValues = "zero"),
                 "a phase was run")
    expect_error(estimateMobilityNetwork(myState, exoEffects, alg, initialParameters = rep(0, 5)),
                 "a phase was run")
    expect_error(estimateMobilityNetwork(myState, mixedEffects, quickAlgorithm(myState, mixedEffects)),
                 "a phase was run")
  })

  # a result of the MPLE can be continued
  capture.output(mple <- estimateMobilityNetwork(myState, exoEffects, alg))
  capture.output(res <- estimateMobilityNetwork(myState, exoEffects, alg, prevAns = mple))
  expect_lt(max(abs(res$estimates - mple$estimates) / mple$standardErrors), 2)

  # without an MPLE, the phases run from zero
  local({
    local_mocked_bindings(getMPLE = function(...) NULL, getMPLEFit = function(...) NULL, runPhase1 = ran)
    expect_error(capture.output(estimateMobilityNetwork(myState, exoEffects, alg)), "a phase was run")
  })
})

test_that("returnDeps simulates only phase 3 at the MPLE of an exogenous model", {
  alg <- createAlgorithm(myState, exoEffects, multinomialProposal = FALSE, burnInN3 = 2000,
                         thinningN3 = 500, iterationsN3 = 300)
  ran <- function(...) stop("a phase was run")
  capture.output(plain <- estimateMobilityNetwork(myState, exoEffects, alg))
  set.seed(31)
  local_mocked_bindings(runPhase1 = ran, runPhase2 = ran)
  capture.output(res <- estimateMobilityNetwork(myState, exoEffects, alg, returnDeps = TRUE))

  expect_identical(res$estimates, plain$estimates)
  expect_identical(res$standardErrors, plain$standardErrors)
  expect_identical(res$covarianceMatrix, plain$covarianceMatrix)
  expect_length(res$deps, 300)
  expect_equal(dim(res$simulatedStatistics), c(300, 5))
  expect_equal(unname(res$observedStatistics),
               unname(getNetworkStatistics(myState$dep.var, myState, myState$cache, exoEffects)))
  # at the maximum likelihood estimate the expected statistics are the observed ones
  expect_lt(max(abs(res$convergenceStatistics)), 0.4)
  expect_false(any(res$convergenceStatistics == 0))

  expect_type(autoCorrelationTest(res), "double")
  candidate <- createEffects(myState) |> addEffect(concentration_AC, origin = "origin", alpha = 4)
  expect_s3_class(scoreTest(res, candidate), "scoretest.realm")
  expect_s3_class(gofMobilityNetwork(ans = res, gofFunction = getIndegree, lvls = 1:80),
                  "gof.stats.realm")
})

test_that("returnDeps works in parallel for an exogenous model", {
  skip_on_cran()
  skip_if_not(testthat::is_checking(), "workers need the installed package")
  alg <- createAlgorithm(myState, exoEffects, multinomialProposal = FALSE, burnInN3 = 500,
                         thinningN3 = 200, iterationsN3 = 40)
  set.seed(41)
  capture.output(res <- estimateMobilityNetwork(myState, exoEffects, alg, returnDeps = TRUE,
                                                parallel = TRUE, cpus = 2))
  expect_length(res$deps, 40)
  expect_equal(unname(res$estimates), exoFit$estimates, tolerance = 1e-6)
})
