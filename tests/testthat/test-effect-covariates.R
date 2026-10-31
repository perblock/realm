# The covariates of effects must be covariates of the state of the kind and
# size the effect reads. Effects index their data directly, so a nodeset, the
# dependent variable, an origin or a covariate of the wrong kind would give
# wrong statistics, or read outside the data in C++.

# the example state with a 'same' matrix for size (for same_covariate)
covariateState <- myState
covariateState$sizeSame <- createNodeVariable(myState$size$data, nodeSet = "organisations",
                                              addSame = TRUE)

effectsWithCovariates <- local({
  ns <- asNamespace("realm")
  Filter(function(name) {
    f <- get(name, envir = ns)
    is.function(f) && "getTargetContribution" %in% names(formals(f)) &&
      any(c("attribute.index", "resource.attribute.index") %in% names(formals(f)))
  }, getNamespaceExports(ns))
})

# a correct specification of each of these effects, and the inputs of the
# wrong kind for each covariate argument
correctArgs <- function(effect) {
  formalNames <- names(formals(get(effect, envir = asNamespace("realm"))))
  args <- list()
  if ("origin" %in% formalNames) args$origin <- "origin"
  if ("attribute.index" %in% formalNames) {
    args$attribute.index <- if (attributeRoles[[effect]] == "dyad") "sameRegion" else
      if (effect == "same_covariate") "sizeSame" else "size"
  }
  if ("resource.attribute.index" %in% formalNames) {
    args$resource.attribute.index <- if (grepl("_cont|additional_origin", effect)) "second_or" else "sex"
  }
  args
}
wrongInputs <- list(
  attribute.index = function(effect) {
    other <- if (attributeRoles[[effect]] == "dyad") "size" else "sameRegion"
    c("organisations", "people", "sex", "transfers", "sizee", "origin", "dep.var", other)
  },
  resource.attribute.index = function(effect) {
    c("organisations", "people", "size", "sameRegion", "transfers", "origin", "sexx")
  }
)

test_that("effects accept covariates of the right kind", {
  for (effect in effectsWithCovariates) {
    expect_no_error(suppressWarnings(
      createEffectsObject(list(c(list(effect), correctArgs(effect))),
                          checkProcessState = covariateState)
    ))
  }
})

test_that("effects reject nodesets, the dependent variable, origins and covariates of the wrong kind", {
  for (effect in effectsWithCovariates) {
    args <- correctArgs(effect)
    for (arg in intersect(names(args), names(wrongInputs))) {
      for (value in wrongInputs[[arg]](effect)) {
        spec <- c(list(effect), args)
        spec[[arg]] <- value
        expect_error(suppressWarnings(createEffectsObject(list(spec), checkProcessState = covariateState)),
                     label = paste(effect, arg, "=", value))
      }
    }
  }
})

test_that("a misspelled covariate is named in the error", {
  expect_error(createEffectsObject(list(list("alter_covariate", attribute.index = "sizee")),
                                   checkProcessState = covariateState),
               paste("Effect alter_covariate: node.attribute must be a covariate that applies",
                     "to the nodeset 'organisations'. Input that caused the error: 'sizee'"),
               fixed = TRUE)
})

test_that("covariates are checked again when estimation and simulation start", {
  # effects created for a state, used with a state where the covariate is a
  # network
  effects <- createEffects(myState) |> addEffect(alter_covariate, node.attribute = "size")
  badState <- myState
  badState$size <- createNetwork(myState$sameRegion$data, nodeSet = c("organisations", "organisations"))
  alg <- createAlgorithm(badState, effects)
  expect_error(estimateMobilityNetwork(badState, effects, alg),
               "Input that caused the error: 'size'")
  expect_error(simulateMobilityNetworks(badState, effects, parameters = 0,
                                        burnin = 10, thinning = 10, nSimulations = 1),
               "Input that caused the error: 'size'")
})

test_that("the C++ effects check the size of their covariates", {
  # effects created without a state are not checked in R; the C++ effects
  # stop instead of reading outside the data
  for (spec in list(list("dyadic_covariate", origin = "origin", attribute.index = "organisations"),
                    list("dyadic_covariate", origin = "origin", attribute.index = "size"),
                    list("alter_covariate", attribute.index = "sameRegion"),
                    list("loops_node_covar", origin = "origin", attribute.index = "people"),
                    list("sim_covariate", origin = "origin", attribute.index = "size",
                         origin.attribute.index = "sex"))) {
    effects <- createEffectsObject(list(spec))
    expect_error(cppNetworkStatistics(myState$dep.var, myState, myState$cache, effects),
                 "the covariate of", label = paste(spec, collapse = " "))
  }
})

test_that("the elements of a state have unique names", {
  expect_error(createProcessState(list(transfers = myState$transfers, origin = myState$origin,
                                       people = myState$people,
                                       organisations = myState$organisations,
                                       organisations = myState$sameRegion),
                                  "transfers"),
               "must be unique; used more than once: organisations")
})
