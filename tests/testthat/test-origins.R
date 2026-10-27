# The dependent variable holds only the destinations; origins are separate
# objects, and effects that depend on an origin name it. A state can have
# several origins, also of another nodeset than the destinations (then the
# origin's network is rectangular and there are no loops).

# the example state with a second origin: the region (1 or 2) of the origin
# organisation, with a region covariate and a region x organisation network;
# in the global environment, as addEffect() finds the state by its name there
local({
  regionOf <- orgRegion + 1
  regions <- createNodeSet(2)
  originRegion <- createOrigin(regionOf[myState$origin$data], nodeSet = c("regions", "people"))
  regionSize <- createNodeVariable(c(5, 9), nodeSet = "regions", addSame = FALSE, addSim = FALSE)
  sameRegionRegion <- createNetwork(outer(1:2, regionOf, "==") * 1,
                                    nodeSet = c("regions", "organisations"))
  for (name in c("transfers", "origin", "people", "organisations", "sameRegion",
                 "size", "sex", "second_or")) {
    assign(name, myState[[name]])
  }
  # the region as a covariate of the regions and of the organisations
  regionId <- createNodeVariable(1:2, nodeSet = "regions", addSame = FALSE, addSim = FALSE)
  orgRegionId <- createNodeVariable(regionOf, nodeSet = "organisations")
  twoOrigins <- realmDataCreate(transfers, origin, originRegion, people, organisations,
                                regions, sameRegion, size, sex, second_or, regionSize,
                                sameRegionRegion, regionId, orgRegionId)
  assign("twoOrigins", twoOrigins, envir = globalenv())
})

twoOriginEffects <- createEffects(twoOrigins) |>
  addEffect(loops, origin = "origin") |>
  addEffect(concentration_GW, origin = "origin", alpha = 2) |>
  addEffect(concentration_AC, origin = "originRegion", alpha = 2) |>
  addEffect(dyadic_covariate, origin = "originRegion", node.attribute = "sameRegionRegion") |>
  addEffect(present_relations, origin = "originRegion") |>
  addEffect(concentration_prop_orig_cov, origin = "originRegion", node.attribute = "regionSize") |>
  addEffect(joining_similar_avoiding_dissimilar_covar_bin, origin = "originRegion",
            edge.attribute = "sex") |>
  addEffect(alter_covariate, node.attribute = "size")

test_that("a state can have several origins, also of another nodeset", {
  expect_equal(twoOrigins$origins, c("origin", "originRegion"))
  cache <- twoOrigins$cache[[twoOrigins$dep.var]]
  expect_equal(dim(cache$originRegion$valuedNetwork), c(2, 17))
  expect_false(cache$originRegion$sameNodeSet)
  expect_true(cache$origin$sameNodeSet)
  # all origins and the destinations count the same resources
  expect_equal(colSums(cache$originRegion$valuedNetwork), cache$destinations$inflow)
  expect_equal(colSums(cache$origin$valuedNetwork), cache$destinations$inflow)
  expect_output(print(twoOrigins), "origin: originRegion with 2 regions")
})

test_that("C++ effects with several origins give the same values as the R effects", {
  expect_cpp_matches_r(twoOriginEffects, state = twoOrigins)
})

test_that("C++ and R chains are identical with several origins", {
  parameters <- c(1, 0.5, 0.3, 0.2, -0.3, 0.1, 0.2, 0.05)
  compare_chains(twoOriginEffects, parameters, n = 3000, state = twoOrigins)
  compare_chains(twoOriginEffects, parameters, n = 300, state = twoOrigins, multinomial = TRUE)
})

test_that("the caches of all origins match a rebuilt cache after many steps", {
  dv <- twoOrigins$dep.var
  resourceCovariates <- names(twoOrigins$cache[[dv]]$destinations$resourceInflows)
  parameters <- c(1, 0.5, 0.3, 0.2, -0.3, 0.1, 0.2, 0.05)
  for (useCpp in c(TRUE, FALSE)) {
    old <- options(realm.useCpp = useCpp)
    set.seed(6)
    res <- simulateNSteps(dv, twoOrigins, twoOrigins$cache, twoOriginEffects, parameters,
                          n = if (useCpp) 5000 else 1000)
    options(old)
    state <- toListState(res$state, dv)
    expect_false(identical(state[[dv]]$data, twoOrigins[[dv]]$data))
    expect_equal(toListCache(res$cache),
                 createInternalCache(state, resourceCovariates = resourceCovariates),
                 tolerance = 1e-12)
  }
})

test_that("change statistics match target differences for a rectangular origin", {
  dv <- twoOrigins$dep.var
  cache <- twoOrigins$cache
  moves <- data.frame(edge = c(1, 20, 300, 500), j_new = c(3, 17, 1, 8))
  for (spec in list(
    list("concentration_AC", origin = "originRegion", alpha = 2),
    list("dyadic_covariate", origin = "originRegion", attribute.index = "sameRegionRegion"),
    list("present_relations", origin = "originRegion"),
    list("concentration_prop_orig_cov", origin = "originRegion", attribute.index = "regionSize"),
    list("concentration_rankGW", origin = "originRegion", lambda = 2)
  )) {
    expect_change_matches(spec[[1]], spec[-1], moves, state = twoOrigins, cache = cache,
                          resourceCovariates = names(cache[[dv]]$destinations$resourceInflows))
  }
})

# same_covariate and sim_covariate comparing a covariate of the regions of
# origin with one of the destination organisations
comparisonEffects <- createEffects(twoOrigins) |>
  addEffect(same_covariate, origin = "originRegion", node.attribute = "orgRegionId",
            origin.attribute = "regionId") |>
  addEffect(sim_covariate, origin = "originRegion", node.attribute = "size",
            origin.attribute = "regionSize") |>
  addEffect(same_covariate, origin = "origin", node.attribute = "orgRegionId",
            origin.attribute = "orgRegionId", exclude.loops = TRUE)

test_that("same and sim covariates compare covariates of the origin and the destinations", {
  expect_cpp_matches_r(comparisonEffects, state = twoOrigins)
  compare_chains(comparisonEffects, c(0.5, 0.3, -0.2), n = 3000, state = twoOrigins)
  compare_chains(comparisonEffects, c(0.5, 0.3, -0.2), n = 300, state = twoOrigins,
                 multinomial = TRUE)
  dv <- twoOrigins$dep.var
  moves <- data.frame(edge = c(1, 20, 300, 500), j_new = c(3, 17, 1, 8))
  for (spec in list(
    list("same_covariate", origin = "originRegion", attribute.index = "orgRegionId",
         origin.attribute.index = "regionId"),
    list("sim_covariate", origin = "originRegion", attribute.index = "size",
         origin.attribute.index = "regionSize")
  )) {
    expect_change_matches(spec[[1]], spec[-1], moves, state = twoOrigins, cache = twoOrigins$cache,
                          resourceCovariates = names(twoOrigins$cache[[dv]]$destinations$resourceInflows))
  }
  # same region of origin and destination, as the dyadic covariate
  dyadic <- createEffects(twoOrigins) |>
    addEffect(dyadic_covariate, origin = "originRegion", node.attribute = "sameRegionRegion")
  expect_equal(unname(getNetworkStatistics(dv, twoOrigins, twoOrigins$cache, comparisonEffects)[1]),
               unname(getNetworkStatistics(dv, twoOrigins, twoOrigins$cache, dyadic)))
})

test_that("with the same covariate for both, same and sim equal the stored matrices", {
  dv <- cppTestState$dep.var
  for (useCpp in c(TRUE, FALSE)) {
    old <- options(realm.useCpp = useCpp)
    for (effect in c("same_covariate", "sim_covariate")) {
      attribute <- if (effect == "same_covariate") "region" else "size"
      stored <- createEffectsObject(list(list(effect, origin = "origin", attribute.index = attribute)),
                                    checkProcessState = cppTestState)
      compared <- createEffectsObject(list(list(effect, origin = "origin", attribute.index = attribute,
                                                origin.attribute.index = attribute)),
                                      checkProcessState = cppTestState)
      expect_identical(unname(getNetworkStatistics(dv, cppTestState, cppTestState$cache, stored)),
                       unname(getNetworkStatistics(dv, cppTestState, cppTestState$cache, compared)))
      set.seed(4)
      a <- simulateNSteps(dv, cppTestState, cppTestState$cache, stored, 0.5, n = 2000)
      set.seed(4)
      b <- simulateNSteps(dv, cppTestState, cppTestState$cache, compared, 0.5, n = 2000)
      expect_identical(toListState(a$state, dv)[[dv]]$data, toListState(b$state, dv)[[dv]]$data)
    }
    options(old)
  }
})

test_that("effects that depend on an origin need a valid one", {
  expect_error(createEffects(myState) |> addEffect(loops), "specify the origin")
  expect_error(createEffects(myState) |> addEffect(loops, origin = "nowhere"),
               "not an origin of the state")
  expect_error(createEffects(myState) |> addEffect(in_weights_AC, origin = "origin"),
               "does not depend on an origin")
  # loops need an origin with the nodeset of the destinations
  expect_error(createEffects(twoOrigins) |> addEffect(loops, origin = "originRegion"),
               "needs an origin with the nodeset of the destinations")
  # same and sim covariates need a covariate of the origin's locations
  expect_error(createEffects(twoOrigins) |>
                 addEffect(same_covariate, origin = "originRegion", node.attribute = "size"),
               "specify the covariate of its locations with origin.attribute")
  expect_error(createEffects(twoOrigins) |>
                 addEffect(sim_covariate, origin = "originRegion", node.attribute = "size",
                           origin.attribute = "size"),
               "must be defined for regions")
  expect_error(createEffects(twoOrigins) |>
                 addEffect(sim_covariate, origin = "originRegion", node.attribute = "regionSize",
                           origin.attribute = "regionSize"),
               "must be defined for organisations")
})

test_that("loops can only be excluded for an origin with the nodeset of the destinations", {
  expect_error(createEffects(twoOrigins) |>
                 addEffect(concentration_AC, origin = "originRegion", exclude.loops = TRUE),
               "exclude.loops needs an origin with the nodeset of the destinations")
  expect_no_error(createEffects(twoOrigins) |>
                    addEffect(concentration_AC, origin = "origin", exclude.loops = TRUE))
})

test_that("covariates must have the nodesets of their role in the effect", {
  # an origin covariate of the organisations for an origin in regions
  expect_error(createEffects(twoOrigins) |>
                 addEffect(concentration_prop_orig_cov, origin = "originRegion",
                           node.attribute = "size"),
               "must be defined for regions")
  # a regions x organisations network for an origin in organisations
  expect_error(createEffects(twoOrigins) |>
                 addEffect(dyadic_covariate, origin = "origin",
                           node.attribute = "sameRegionRegion"),
               "must be defined for organisations x organisations")
  # a covariate of the origins where one of the destinations is needed
  expect_error(createEffects(twoOrigins) |> addEffect(alter_covariate, node.attribute = "regionSize"),
               "must be defined for organisations")
  expect_error(createEffects(twoOrigins) |>
                 addEffect(loops_resource_covar, origin = "origin",
                           edge.attribute = "size"),
               "must be defined for people")
})

test_that("every effect with a node attribute declares its role", {
  ns <- asNamespace("realm")
  withAttribute <- Filter(function(name) {
    f <- get(name, envir = ns)
    is.function(f) && "attribute.index" %in% names(formals(f)) &&
      "getTargetContribution" %in% names(formals(f))
  }, getNamespaceExports(ns))
  expect_setequal(withAttribute, names(attributeRoles))
})

test_that("the dependent variable and the origins are created from vectors", {
  expect_error(createDependent(mobilityEdgelist, c("organisations", "people")),
               "contains only the destinations")
  expect_error(createOrigin(mobilityEdgelist, c("organisations", "people")),
               "one location per person")
  transfers <- createDependent(mobilityEdgelist[, 2], nodes = "organisations", edges = "people")
  expect_s3_class(transfers, "dependent.realm")
  expect_equal(transfers$nodeSet, c("organisations", "people"))
  expect_error(createDependent(c(1, 2, NA), c("organisations", "people")), "missing values")
  expect_error(createDependent(c(0, 2), c("organisations", "people")), "numbered from one")
})

test_that("a state needs an origin of the same resources", {
  people <- createNodeSet(nrow(mobilityEdgelist))
  organisations <- createNodeSet(17)
  transfers <- createDependent(mobilityEdgelist[, 2], c("organisations", "people"))
  expect_error(realmDataCreate(transfers, people, organisations), "at least one origin")
  origin <- createOrigin(mobilityEdgelist[1:10, 1], c("organisations", "people"))
  expect_error(realmDataCreate(transfers, origin, people, organisations),
               "one location for each of the people")
  destinations <- createOrigin(mobilityEdgelist[, 1], c("organisations", "people"))
  expect_error(realmDataCreate(transfers, destinations, people, organisations),
               "cannot be used as the name of an origin")
})

test_that("loops_additional_origin is deprecated and equals loops with a second origin", {
  expect_warning(deprecated <- createEffects(myState) |>
                   addEffect(loops_additional_origin, edge.attribute = "second_or"),
                 "deprecated")
  # the additional origin of the example data as an origin
  local({
    secondOrigin <- myState
    secondOrigin$second_origin <- createOrigin(myState$second_or$data,
                                               nodeSet = c("organisations", "people"))
    secondOrigin$origins <- c("origin", "second_origin")
    secondOrigin$cache <- createInternalCache(
      secondOrigin, resourceCovariates = names(myState$cache[[myState$dep.var]]$destinations$resourceInflows)
    )
    assign("secondOrigin", secondOrigin, envir = globalenv())
  })
  on.exit(rm("secondOrigin", envir = globalenv()), add = TRUE)
  loopsSecond <- createEffects(secondOrigin) |> addEffect(loops, origin = "second_origin")
  dv <- myState$dep.var
  for (useCpp in c(TRUE, FALSE)) {
    old <- options(realm.useCpp = useCpp)
    expect_equal(unname(getNetworkStatistics(dv, secondOrigin, secondOrigin$cache, loopsSecond)),
                 unname(getNetworkStatistics(dv, myState, myState$cache, deprecated)))
    options(old)
  }
  set.seed(9)
  for (edge in sample(length(myState[[dv]]$data), 50)) {
    for (j in c(myState$second_or$data[edge], sample(17, 1))) {
      for (update in c(-1, 1)) {
        expect_equal(
          loops(dep.var = dv, origin = "second_origin", state = secondOrigin,
                cache = secondOrigin$cache, i = secondOrigin$second_origin$data[edge],
                j = j, edge = edge, update = update),
          loops_additional_origin(dep.var = dv, resource.attribute.index = "second_or",
                                  state = myState, cache = myState$cache, i = NA,
                                  j = j, edge = edge, update = update)
        )
      }
    }
  }
})
