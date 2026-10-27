# This script runs a simple example with the data from the realm package

library(realm)


##### create data objects from internal data files, 
# which are later combined to the process state #####

# extract number of individuals and organisations from the mobility data

N_ind <- nrow(mobilityEdgelist)
N_org <- length(unique(as.numeric(mobilityEdgelist)))

# Create a process state out of the mobility data objects:
# create objects (which are later combined to the process state).
# The dependent variable contains the destinations (second column of the
# mobility data), the origins are a separate object (first column)
transfers <- realmDependent(mobilityEdgelist[, 2],
                            nodes = "organisations",
                            edges = "people")
origin <- realmOrigin(mobilityEdgelist[, 1],
                      nodes = "organisations",
                      edges = "people")

people <- realmEdges(N_ind)
organisations <- realmNodes(N_org)

sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- dyadicCovar(sameRegion, nodes = "organisations")

region <- monadicCovar(orgRegion, nodes = "organisations")
size <- monadicCovar(orgSize, nodes = "organisations")
sex <- monadicCovar(indSex, edges = "people")


# the following lines create an artificial second origin used for illustration
# in the examples in the manual
other_origin <- sample(1:17, 742, replace = T)
resample <- as.logical(sample(0:1, 742, replace = T, prob = c(0.88, 0.12)))
other_origin[resample] <- transfers$data[resample]
second_or <- monadicCovar(other_origin, edges = "people")


# combine created objects to the process state
myState <- realmDataCreate(transfers,
                           origin,
                           people,
                           organisations,
                           sameRegion,
                           region,
                           size,
                           sex,
                           second_or,
                           fixedEffectDummies = TRUE)

# inspect the created object
myState

##### create effects object #####

# effects object
myEffects <- createEffects(myState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(concentration_AC, origin = "origin", alpha = 4) |>
  addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion") |>
  addEffect(alter_covariate, node.attribute = "size") |>
  addEffect(resource_covar_to_node_covar,
            node.attribute = "region",
            edge.attribute = "sex") |>
  addEffect(loops_resource_covar, origin = "origin", edge.attribute = "sex")

# inspect the created object
myEffects

# further effects object with fixed effects by location

myEffects_fe <- createEffects(myState) |>
  addEffect(loops, origin = "origin") |>
  addEffect(concentration_AC, origin = "origin", alpha = 4) |>
  addEffect(dyadic_covariate, origin = "origin", node.attribute = "sameRegion") |>
  addEffect(resource_covar_to_node_covar,
            node.attribute = "region",
            edge.attribute = "sex") |>
  addEffect(loops_resource_covar, origin = "origin", edge.attribute = "sex") |>
  addFixedEffects(myState)

myEffects_fe


##### create algorithm object #####

# define algorithm based on state and effects characteristics
myAlg <- realmAlgorithmCreate(myState, myEffects, nsubN2 = 3,
                              multinomialProposal = FALSE)


##### estimate mobility network model #####

# mobility network model 
myResDN <- realmEstimate(
  myState, myEffects, myAlg,
  parallel = TRUE, cpus = 4,
  verbose = TRUE,
  returnDeps = TRUE,
  fish = FALSE
)

# myResDN_fe <- realmEstimate(
#   myState, myEffects_fe, myAlg,
#   initialParameters = NULL,
#   # in case a pseudo-likelihood estimation was run, replace with
#   # initialParameters = initEst,
#   parallel = TRUE, cpus = 4,
#   verbose = TRUE,
#   returnDeps = TRUE,
#   fish = FALSE
# )

# check convergence
max(abs(myResDN$convergenceStatistics))

myResDN_old <- myResDN

# estimate mobility network model again based on previous results to improve convergence
# with an adjusted algorithm
myAlg <- realmAlgorithmCreate(myState, myEffects, multinomialProposal = TRUE, 
                              initialIterationsN2 = 100, nsubN2 = 1, initGain = 0.05, iterationsN3 = 1000)

myResDN <- realmEstimate(
  myState, myEffects, myAlg,
  prevAns = myResDN,
  parallel = TRUE, cpus = 4,
  verbose = TRUE,
  returnDeps = TRUE,
  fish = FALSE
)

# check convergence
max(abs(myResDN$convergenceStatistics))

# view results
myResDN


##### regression diagnostics #####

autoCorrelationTest(myResDN)

traces <- extractTraces(myResDN, myEffects)
plot(traces)


##### test whether other effects should be included #####

myEffects2 <- createEffects(myState) |>
  addEffect(in_weights_AC, alpha = 4)

test_ME.2 <- scoreTest(myResDN, myEffects2)
test_ME.2


##### goodness of fit #####

myGofIndegree <- realmGOF(ans = myResDN, 
                          gofFunction = getIndegree, 
                          lvls = 1:100)
plot(myGofIndegree,  lvls = 20:70)

myGofTieWeight <- realmGOF(ans = myResDN, 
                           gofFunction = getTieWeights, 
                           lvls = 1:30,
                           origin = "origin")
plot(myGofTieWeight, lvls = 1:15)


##### simulate mobility network #####

mySimDN <- realmSimulate(myState,
                         myEffects,
                         parameters = c(2, 1, 0.5, 0.1, -1, -0.5),
                         burnin = 45000,
                         thinning = 15000,
                         nSimulations = 10
)

mySimDN[[1]]
