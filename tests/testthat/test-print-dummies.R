# By default the printed effects and state leave out the location dummies and
# say how many there are.

# addEffect() looks up the state by name in the global environment
assign("printState", local({
  transfers <- myState$transfers
  origin <- myState$origin
  people <- myState$people
  organisations <- myState$organisations
  size <- myState$size
  realmDataCreate(transfers, origin, people, organisations, size,
                  fixedEffectDummies = TRUE)
}), envir = globalenv())
nDummies <- length(myState$organisations$ids) - 1

dummyEffects <- createEffects(printState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(alter_covariate, node.attribute = "size") |>
  addEffect(origin_to_node_covar, origin = "origin", node.attribute = "dummy2", dummy.value = 2) |>
  addFixedEffects(printState)

notShown <- "included but not shown. Use print(x, showDummies = TRUE) to show them."

test_that("printed effects leave out the fixed effects and say how many", {
  printed <- capture.output(print(dummyEffects))
  expect_length(printed, 2 + 3 + 1)
  expect_false(any(grepl("^ *alter_covariate +- +dummy", printed)))
  expect_match(printed[4], "^ *alter_covariate +- +size")
  expect_match(printed[5], "^ *origin_to_node_covar +origin +dummy2 +- +d 2")
  expect_equal(printed[6], paste(nDummies, "fixed effects are", notShown))
  # as when the object is printed at the console
  expect_equal(capture.output(print(dummyEffects)), printed)
})

test_that("printed effects show the fixed effects with showDummies = TRUE", {
  printed <- capture.output(print(dummyEffects, showDummies = TRUE))
  expect_length(printed, 2 + 3 + nDummies)
  expect_equal(sum(grepl("^ *alter_covariate +- +dummy[0-9]+", printed)), nDummies)
  expect_false(any(grepl("not shown", printed)))
})

test_that("printed effects handle models with only fixed effects, one, or none", {
  onlyFixed <- createEffects(printState) |> addFixedEffects(printState)
  expect_equal(capture.output(print(onlyFixed)),
               c("Effects", paste(nDummies, "fixed effects are", notShown)))

  oneFixed <- createEffects(printState) |>
    addEffect(alter_covariate, node.attribute = "dummy3")
  expect_equal(capture.output(print(oneFixed)),
               c("Effects", paste("1 fixed effect is", notShown)))

  empty <- createEffects(printState)
  expect_equal(capture.output(print(empty)),
               c("Effects", "The effects object does not contain any effects."))

  # nothing changes for effects without dummies
  noDummies <- createEffects(myState) |> addEffect(alter_covariate, node.attribute = "size")
  expect_false(any(grepl("not shown", capture.output(print(noDummies)))))
})

test_that("the printed state leaves out the dummies and says how many", {
  printed <- capture.output(print(printState))
  expect_false(any(grepl("^ *dummy[0-9]+", printed)))
  expect_true(any(grepl("^ *size ", printed)))
  line <- which(printed == paste(nDummies, "fixed-effect dummy variables are", notShown))
  expect_length(line, 1)
  # under the table of the nodeset of the dummies
  expect_lt(which(grepl("covariates of organisations", printed)), line)
  expect_lt(line, which(grepl("covariates of people", printed)))
})

test_that("the printed state shows the dummies with showDummies = TRUE", {
  printed <- capture.output(print(printState, showDummies = TRUE))
  expect_equal(sum(grepl("^ *dummy[0-9]+ ", printed)), nDummies)
  expect_false(any(grepl("not shown", printed)))
})

test_that("a nodeset with only dummies gives the line and no table", {
  onlyDummies <- printState
  onlyDummies$size <- NULL
  printed <- capture.output(print(onlyDummies))
  expect_false(any(grepl("no covariates available for nodeset organisations", printed)))
  expect_false(any(grepl("cov. name", printed)))
  expect_true(any(grepl("no covariates available for nodeset people", printed)))
  expect_equal(sum(grepl("not shown", printed)), 1)

  oneDummy <- onlyDummies
  oneDummy[paste0("dummy", 3:(nDummies + 1))] <- NULL
  expect_true(paste("1 fixed-effect dummy variable is", notShown) %in%
                capture.output(print(oneDummy)))
})

test_that("states without dummies print as before", {
  printed <- capture.output(print(myState))
  expect_false(any(grepl("not shown", printed)))
  expect_true(any(grepl("^ *size ", printed)))
})

rm("printState", envir = globalenv())
