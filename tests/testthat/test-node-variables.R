# createNodeVariables(): covariates from a data.frame, each column becoming an
# ordinary covariate of the state, or from a vector

baseObjects <- local({
  transfers <- createDependent(mobilityEdgelist[, 2], nodes = "organisations", edges = "people")
  origin <- createOrigin(mobilityEdgelist[, 1], nodes = "organisations", edges = "people")
  list(transfers = transfers, origin = origin, people = createNodeSet(nrow(mobilityEdgelist)),
       organisations = createNodeSet(length(orgRegion)))
})

test_that("a data.frame gives the same state as covariates created one by one", {
  oneByOne <- createProcessState(c(baseObjects, list(
    region = createNodeVariable(orgRegion, nodeSet = "organisations"),
    size = createNodeVariable(orgSize, nodeSet = "organisations"),
    sex = createNodeVariable(indSex, nodeSet = "people")
  )), "transfers")
  fromData <- createProcessState(c(baseObjects, list(
    orgCovars = monadicCovars(data.frame(region = orgRegion, size = orgSize), nodes = "organisations"),
    sex = monadicCovars(indSex, edges = "people")
  )), "transfers")
  expect_identical(fromData, oneByOne)
  # the same in chains, through the effects
  specs <- list(list("loops", origin = "origin"),
                list("alter_covariate", attribute.index = "size"),
                list("sim_covariate", origin = "origin", attribute.index = "size"),
                list("resource_covar_to_node_covar", attribute.index = "region",
                     resource.attribute.index = "sex"))
  effects <- createEffectsObject(specs, checkProcessState = fromData)
  compare_chains(effects, c(1, 0.1, 0.5, -0.5), n = 3000, state = fromData)
  compare_chains(effects, c(1, 0.1, 0.5, -0.5), n = 300, state = fromData, multinomial = TRUE)
})

test_that("realmDataCreate() names covariates by columns and vector objects", {
  transfers <- baseObjects$transfers
  origin <- baseObjects$origin
  people <- baseObjects$people
  organisations <- baseObjects$organisations
  orgCovars <- monadicCovars(data.frame(region = orgRegion, wage = orgSize), nodes = "organisations")
  oneColumn <- monadicCovars(data.frame(age = rnorm(742)), edges = "people")
  size <- monadicCovars(orgSize * 2, nodes = "organisations")
  state <- realmDataCreate(transfers, origin, people, organisations, orgCovars, oneColumn, size)
  expect_true(all(c("region", "wage", "age", "size") %in% names(state)))
  expect_false(any(c("orgCovars", "oneColumn") %in% names(state)))
  expect_s3_class(size, "nodeVar.realm")
  expect_equal(state$size$data, orgSize * 2)
  expect_equal(state$wage$data, as.vector(orgSize))
  expect_true("age" %in% names(state$cache[[state$dep.var]]$destinations$resourceInflows))
  expect_output(print(state), "wage")
  # a location column is no resource covariate
  expect_error(createEffectsObject(list(list("loops_resource_covar", origin = "origin",
                                             resource.attribute.index = "wage")),
                                   checkProcessState = state),
               "must be a covariate that applies to the nodeset 'people'")
})

test_that("covariates of the locations of another origin come from a data.frame", {
  objects <- c(baseObjects, list(
    regions = createNodeSet(2),
    originRegion = createOrigin(orgRegion[mobilityEdgelist[, 1]] + 1, nodes = "regions", edges = "people"),
    regionCovars = monadicCovars(data.frame(regionSize = c(5, 9)), nodes = "regions", addSame = FALSE)
  ))
  state <- createProcessState(objects, "transfers")
  effects <- createEffectsObject(list(list("concentration_prop_orig_cov", origin = "originRegion",
                                           attribute.index = "regionSize")),
                                 checkProcessState = state)
  expect_length(effects$name, 1)
})

test_that("matrices with column names work, and attributes of columns are dropped", {
  m <- monadicCovars(cbind(a = orgSize, b = orgRegion), nodes = "organisations", addSame = FALSE)
  expect_equal(names(m$variables), c("a", "b"))
  labelled <- structure(c(1, 2, 3), class = "someLabels", labels = c(low = 1))
  v <- monadicCovars(data.frame(x = I(labelled)), nodes = "n")
  expect_equal(v$variables$x$data, c(1, 2, 3))
})

test_that("wrong input to createNodeVariables stops", {
  expect_error(monadicCovars(data.frame(a = c("x", "y")), nodes = "n"),
               "must be numeric; recode them first. Input that caused the error: 'a'")
  expect_error(monadicCovars(data.frame(a = factor(c("x", "y"))), nodes = "n"), "must be numeric")
  expect_error(monadicCovars(factor(c("x", "y")), nodes = "n"), "must be numeric")
  expect_error(monadicCovars(matrix(1:4, 2), nodes = "n"), "must have a name")
  duplicated <- data.frame(a = 1:2, b = 1:2)
  names(duplicated) <- c("a", "a")
  expect_error(monadicCovars(duplicated, nodes = "n"), "the names must be unique")
  expect_error(monadicCovars(data.frame(a = c(1, NA)), nodes = "n"), "missing values")

  size <- monadicCovars(data.frame(people = orgSize), nodes = "organisations")
  expect_error(createProcessState(c(baseObjects, list(covars = size)), "transfers"),
               "must be unique; used more than once: people")
  wrongRows <- monadicCovars(data.frame(x = 1:3), nodes = "organisations")
  expect_error(createProcessState(c(baseObjects, list(covars = wrongRows)), "transfers"),
               "not of the same size")
  reserved <- monadicCovars(data.frame(cache = orgSize), nodes = "organisations")
  expect_error(createProcessState(c(baseObjects, list(covars = reserved)), "transfers"),
               "Input that caused the error: 'cache'")
})

test_that("addDummyEffects() only adds the dummies of realmDataCreate()", {
  transfers <- baseObjects$transfers
  origin <- baseObjects$origin
  people <- baseObjects$people
  organisations <- baseObjects$organisations
  orgCovars <- monadicCovars(data.frame(dummy_public = orgRegion), nodes = "organisations")
  dummyState <- realmDataCreate(transfers, origin, people, organisations, orgCovars,
                                fixedEffectDummies = TRUE)
  assign("dummyState", dummyState, envir = globalenv())
  on.exit(rm("dummyState", envir = globalenv()), add = TRUE)
  effects <- addDummyEffects(createEffects(dummyState), dummyState)
  expect_length(effects$name, 16)
  expect_false(any(grepl("dummy_public", effects$name)))
})
