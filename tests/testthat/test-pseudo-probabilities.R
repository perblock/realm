# getPseudoProbabilities gives the probability of each resource's observed
# destination given the destinations of the other resources.

test_that("getPseudoProbabilities is exported", {
  expect_true("getPseudoProbabilities" %in% getNamespaceExports("realm"))
})

test_that("with exogenous effects they are the probabilities of the multinomial logit", {
  o <- myState$origin$data
  d <- myState$transfers$data
  sex <- myState$sex$data
  size <- myState$size$data
  sameRegion <- myState$sameRegion$data
  nDestinations <- length(size)

  effects <- createEffects(myState) |>
    addEffect(loops, origin = "origin") |>
    addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion") |>
    addEffect(alter_covariate, node.attribute = "size") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex")
  capture.output(fit <- estimateMobilityNetwork(myState, effects, createAlgorithm(myState, effects)))
  p <- getPseudoProbabilities(fit, effects)

  # independently of realm: a Poisson glm on the raw data, one group per
  # combination of origin and sex
  groups <- paste(o, sex)
  rows <- do.call(rbind, lapply(unique(groups), function(g) {
    first <- which(groups == g)[1]
    data.frame(group = g, destination = seq_len(nDestinations),
               count = tabulate(d[groups == g], nDestinations), n = sum(groups == g),
               loops = as.numeric(o[first] == seq_len(nDestinations)),
               dyadic = sameRegion[o[first], ],
               alter = size,
               sexSize = sex[first] * size)
  }))
  model <- glm(count ~ 0 + factor(group) + loops + dyadic + alter + sexSize, family = poisson,
               data = rows, control = list(epsilon = 1e-12, maxit = 100))
  expected <- (fitted(model) / rows$n)[match(paste(groups, d), paste(rows$group, rows$destination))]

  expect_type(p, "double")
  expect_null(dim(p))
  expect_null(names(p))
  expect_length(p, length(d))
  expect_equal(p, unname(expected), tolerance = 1e-6)
})

test_that("with endogenous effects they give the log pseudo-likelihood", {
  p <- getPseudoProbabilities(myResDN, myEffects)
  expect_length(p, length(myState$transfers$data))
  expect_true(all(p > 0 & p < 1))
  evaluate <- mpleObjectiveCpp(myState$dep.var, myState, myState$cache, myEffects)
  expect_equal(sum(log(p)), evaluate(myResDN$estimates, derivatives = FALSE)$logPL)
})

test_that("getPseudoProbabilities rejects what is not a fit or does not fit the effects", {
  expect_error(getPseudoProbabilities(myState, myEffects),
               "ans must be the result of an estimation. Input that caused the error: an object of class 'processState.realm'",
               fixed = TRUE)

  reordered <- myEffects
  reordered$name <- myEffects$name[c(2, 1, 3:6)]
  reordered$effectFormulas <- myEffects$effectFormulas[c(2, 1, 3:6)]
  expect_error(getPseudoProbabilities(myResDN, reordered),
               paste0("in the same order. Input that caused the error: '", myEffects$name[2], "'"),
               fixed = TRUE)

  fewer <- createEffects(myState) |> addEffect(loops, origin = "origin")
  expect_error(getPseudoProbabilities(myResDN, fewer),
               "Input that caused the error: 'number of effects: 1 instead of 6'", fixed = TRUE)

  other <- myEffects
  other$name[6] <- "another effect"
  expect_error(getPseudoProbabilities(myResDN, other),
               "Input that caused the error: 'another effect'", fixed = TRUE)
  expect_error(getPseudoProbabilities(myResDN, list()), "in the same order")
})

# a fitted model with the given parameters
pseudoFit <- function(state, effects, theta) {
  structure(list(estimates = setNames(theta, effects$name), state = state), class = "result.realm")
}

# the probabilities in C++ and in R
bothRoutes <- function(fit, effects) {
  expect_true(useCppSampler(effects, fit$state))
  old <- options(realm.useCpp = FALSE)
  on.exit(options(old))
  r <- getPseudoProbabilities(fit, effects)
  options(realm.useCpp = TRUE)
  list(cpp = getPseudoProbabilities(fit, effects), r = r)
}

test_that("C++ and R give the same probabilities", {
  routes <- bothRoutes(myResDN, myEffects)
  expect_type(routes$cpp, "double")
  expect_null(attributes(routes$cpp))
  expect_equal(routes$cpp, routes$r, tolerance = 1e-12)
})

test_that("C++ and R agree for effects that differ between resources (grouped in C++)", {
  effects <- createEffects(myState) |>
    addEffect(loops, origin = "origin") |>
    addEffect(concentration_AC, origin = "origin", alpha = 4) |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "region", edge.attribute = "second_or") |>
    addEffect(loops_resource_covar, origin = "origin", edge.attribute = "second_or") |>
    addEffect(concentration_AC_resource_covar_bin, origin = "origin", alpha = 4,
              edge.attribute = "sex") |>
    addEffect(dyadic_covariate_resource_attribute, origin = "origin",
              node.attribute = "sameRegion", edge.attribute = "sex")
  set.seed(5)
  routes <- bothRoutes(pseudoFit(myState, effects, rnorm(7, 0, 0.3)), effects)
  expect_equal(routes$cpp, routes$r, tolerance = 1e-12)
  # resources that differ in a covariate of the model differ in probability
  expect_gt(length(unique(routes$cpp)), length(unique(paste(myState$origin$data, myState$transfers$data))))
})

test_that("C++ and R agree with large utilities", {
  effects <- createEffects(myState) |> addEffect(alter_covariate, node.attribute = "size")
  routes <- bothRoutes(pseudoFit(myState, effects, 200), effects)
  expect_false(anyNA(routes$cpp))
  expect_equal(routes$cpp, routes$r, tolerance = 1e-12)
})

test_that("resources and destinations not considered when sampling give NA or are no choice", {
  effects <- createEffects(myState) |>
    addEffect(loops, origin = "origin") |>
    addEffect(alter_covariate, node.attribute = "size")
  theta <- c(1.2, 0.3)
  full <- getMultinomialStatistics(myState, effects)
  utility <- as.vector(as.matrix(full[, -(1:3)]) %*% theta)
  d <- myState$transfers$data
  excludedDestination <- 3
  excludedResources <- c(1, 50, 700)

  state <- myState
  state$people$considerWhenSampling[excludedResources] <- FALSE
  state$organisations$considerWhenSampling[excludedDestination] <- FALSE
  routes <- bothRoutes(pseudoFit(state, effects, theta), effects)

  expected <- vapply(seq_along(d), function(e) {
    u <- utility[full$resource == e]
    if (e %in% excludedResources || d[e] == excludedDestination) return(NA_real_)
    exp(u[d[e]]) / sum(exp(u[-excludedDestination]))
  }, numeric(1))
  expect_equal(routes$cpp, expected, tolerance = 1e-12)
  expect_equal(routes$r, expected, tolerance = 1e-12)
  expect_true(all(is.na(routes$cpp[excludedResources])))
  expect_true(all(is.na(routes$cpp[d == excludedDestination])))
  expect_equal(sum(is.na(routes$cpp)),
               length(union(excludedResources, which(d == excludedDestination))))

  # only resources excluded: the same as the log pseudo-likelihood
  state <- myState
  state$people$considerWhenSampling[excludedResources] <- FALSE
  p <- getPseudoProbabilities(pseudoFit(state, effects, theta), effects)
  evaluate <- mpleObjectiveCpp(state$dep.var, state, state$cache, effects)
  expect_equal(sum(log(p), na.rm = TRUE), evaluate(theta, derivatives = FALSE)$logPL)
})

test_that("C++ and R both give NA for resources with a missing covariate of the model", {
  # createNodeVariable() refuses missing values; a state can still contain them
  effects <- createEffects(myState) |>
    addEffect(loops, origin = "origin") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex")
  state <- myState
  state$sex$data[c(5, 300)] <- NA
  routes <- bothRoutes(pseudoFit(state, effects, c(1, 0.2)), effects)
  expect_equal(which(is.na(routes$cpp)), c(5L, 300L))
  expect_equal(routes$cpp, routes$r, tolerance = 1e-12)
})

test_that("C++ and R agree with two origins, also of different nodesets", {
  local({
    regionOf <- orgRegion + 1
    regions <- createNodeSet(2)
    originRegion <- createOrigin(regionOf[myState$origin$data], nodeSet = c("regions", "people"))
    for (name in c("transfers", "origin", "people", "organisations", "size", "sex")) {
      assign(name, myState[[name]])
    }
    assign("pseudoOrigins",
           realmDataCreate(transfers, origin, originRegion, people, organisations, regions, size, sex),
           envir = globalenv())
  })
  on.exit(rm("pseudoOrigins", envir = globalenv()), add = TRUE)
  effects <- createEffects(pseudoOrigins) |>
    addEffect(loops, origin = "origin") |>
    addEffect(concentration_AC, origin = "origin", alpha = 4) |>
    addEffect(concentration_AC, origin = "originRegion", alpha = 24) |>
    addEffect(origin_to_node_covar, origin = "originRegion", node.attribute = "size") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size", edge.attribute = "sex")
  set.seed(6)
  routes <- bothRoutes(pseudoFit(pseudoOrigins, effects, rnorm(5, 0, 0.2)), effects)
  expect_equal(routes$cpp, routes$r, tolerance = 1e-12)
})

test_that("the R route is used for effects without a C++ version", {
  assign("myLoopsPseudo", function(dep.var = 1, origin, state, cache, i, j, edge, update,
                                   getTargetContribution = FALSE) {
    loops(dep.var, origin, state, cache, i, j, edge, update, getTargetContribution)
  }, envir = globalenv())
  on.exit(rm("myLoopsPseudo", envir = globalenv()), add = TRUE)
  own <- createEffects(myState) |> addEffect(myLoopsPseudo, origin = "origin")
  standard <- createEffects(myState) |> addEffect(loops, origin = "origin")
  expect_false(useCppSampler(own, myState))
  expect_equal(getPseudoProbabilities(pseudoFit(myState, own, 1), own),
               getPseudoProbabilities(pseudoFit(myState, standard, 1), standard),
               tolerance = 1e-12)
})

