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
