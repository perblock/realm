# C++ versions of the exogenous, in-ties, present-relations and
# covariate-based effects (src/effects.cpp)

otherEffects <- createEffects(cppTestState) |>
  addEffect(same_covariate, origin = "origin", node.attribute = "region") |>
  addEffect(sim_covariate, origin = "origin", node.attribute = "size") |>
  addEffect(dyadic_covariate_resource_attribute, origin = "origin", node.attribute = "sameRegion",
            edge.attribute = "sex") |>
  addEffect(dyadic_covariate_resource_attribute, origin = "origin", node.attribute = "size_similarity",
            edge.attribute = "second_or") |>
  addEffect(in_ties_loops, origin = "origin") |>
  addEffect(present_relations, origin = "origin")

covariateEffects <- createEffects(cppTestState) |>
  addEffect(staying_by_prop_bin_inflow, origin = "origin", edge.attribute = "sex") |>
  addEffect(crowding_out_prop_covar_bin, origin = "origin", edge.attribute = "sex") |>
  addEffect(joining_similar_avoiding_dissimilar_covar_bin, origin = "origin", edge.attribute = "sex") |>
  addEffect(avoiding_dissimilar_covar_bin, origin = "origin", edge.attribute = "sex") |>
  addEffect(joining_similar_avoiding_dissimilar_covar_cont, origin = "origin", edge.attribute = "second_or") |>
  addEffect(avoiding_dissimilar_covar_cont, origin = "origin", edge.attribute = "second_or")

test_that("C++ exogenous, in-ties and present-relations effects give the same values as R", {
  expect_cpp_matches_r(otherEffects, state = cppTestState)
})

test_that("C++ covariate-based effects give the same values as R", {
  expect_cpp_matches_r(covariateEffects, state = cppTestState)
})

test_that("C++ and R chains are identical with these effects", {
  compare_chains(otherEffects, c(0.3, 0.2, 0.1, 0.01, 0.001, 0.2),
                 n = 5000, state = cppTestState)
  compare_chains(otherEffects, c(0.3, 0.2, 0.1, 0.01, 0.001, 0.2),
                 n = 500, state = cppTestState, multinomial = TRUE)
  compare_chains(covariateEffects, c(0.3, -0.2, 0.5, 0.3, 0.2, -0.05),
                 n = 5000, state = cppTestState)
  compare_chains(covariateEffects, c(0.3, -0.2, 0.5, 0.3, 0.2, -0.05),
                 n = 500, state = cppTestState, multinomial = TRUE)
})

# the effects that can exclude loops, with exclude.loops = TRUE
noLoopsEffects <- createEffects(cppTestState) |>
  addEffect(concentration_AC, origin = "origin", alpha = 2, exclude.loops = TRUE) |>
  addEffect(concentration_GW_dyad_covar, origin = "origin", node.attribute = "size_similarity",
            alpha = 2, exclude.loops = TRUE) |>
  addEffect(concentration_AC_resource_covar_bin, origin = "origin", edge.attribute = "sex",
            alpha = 2, exclude.loops = TRUE) |>
  addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion",
            exclude.loops = TRUE) |>
  addEffect(dyadic_covariate_resource_attribute, origin = "origin",
            node.attribute = "sameRegion", edge.attribute = "sex", exclude.loops = TRUE) |>
  addEffect(same_covariate, origin = "origin", node.attribute = "region", exclude.loops = TRUE) |>
  addEffect(present_relations, origin = "origin", exclude.loops = TRUE)

test_that("C++ effects that exclude loops give the same values as R", {
  expect_cpp_matches_r(noLoopsEffects, state = cppTestState)
  parameters <- c(0.3, 0.1, 0.2, 0.3, 0.1, 0.2, -0.2)
  compare_chains(noLoopsEffects, parameters, n = 5000, state = cppTestState)
  compare_chains(noLoopsEffects, parameters, n = 500, state = cppTestState, multinomial = TRUE)
})

test_that("loops count unless exclude.loops = TRUE", {
  dv <- myState$dep.var
  with <- createEffects(myState) |> addEffect(dyadic_covariate, origin = "origin",
                                              node.attribute = "sameRegion")
  without <- createEffects(myState) |> addEffect(dyadic_covariate, origin = "origin",
                                                 node.attribute = "sameRegion",
                                                 exclude.loops = TRUE)
  loops <- sum(diag(myState$cache[[dv]]$origin$valuedNetwork))
  expect_equal(unname(getNetworkStatistics(dv, myState, myState$cache, with)) -
                 unname(getNetworkStatistics(dv, myState, myState$cache, without)),
               loops)
})

test_that("all exported effects have a C++ version", {
  ns <- asNamespace("realm")
  isEffect <- vapply(getNamespaceExports(ns), function(name) {
    f <- get(name, envir = ns)
    is.function(f) && all(c("state", "cache", "getTargetContribution") %in% names(formals(f)))
  }, logical(1))
  expect_setequal(names(isEffect)[isEffect], cppEffects)
})
