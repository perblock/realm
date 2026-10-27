########## hiddenFunctions

# binarizeNetwork
binarizeNetwork <- function(network) {
  net <- 1 * (network >= 1)
  diag(net) <- 0
  net
}


# destinationNodeSet, resourceNodeSet
# the nodesets of the dependent variable: the locations (destinations) and the
# resources
destinationNodeSet <- function(state, dep.var = state$dep.var) {
  state[[dep.var]]$nodeSet[1]
}
resourceNodeSet <- function(state, dep.var = state$dep.var) {
  state[[dep.var]]$nodeSet[2]
}

# nonCovariates
# the names of the elements of a state that are not covariates
nonCovariates <- function(state) {
  dep.var <- state$dep.var
  nodesets <- unique(c(state[[dep.var]]$nodeSet,
                       unlist(lapply(state$origins, function(o) state[[o]]$nodeSet))))
  c(dep.var, state$origins, nodesets, "dep.var", "origins", "cache")
}


# checkProcessState
checkProcessState <- function(state) {
  
  # extract dependent variable
  dep.var <- state$dep.var
  
  # make sure that dependent variable has a corresponding object
  if (!(dep.var %in% names(state))) {
    stop("Dependent variable has no corresponding object in the state object.")
  }
  
  # is the dependent variable of class 'dependent.realm'?
  if (!(inherits(state[[dep.var]], "dependent.realm"))) {
    stop("Dependent variable is not of class 'dependent.realm'.")
  }
  
  # are there origins, and are they of class 'origin.realm'?
  origins <- state$origins
  if (length(origins) == 0) {
    stop("The state needs at least one origin (see createOrigin()).")
  }
  for (o in origins) {
    if (o %in% c("destinations", "dep.var", "origins", "cache")) {
      stop(paste0("'", o, "' cannot be used as the name of an origin."))
    }
    if (!(o %in% names(state)) || !inherits(state[[o]], "origin.realm")) {
      stop(paste(o, "is not an object of class 'origin.realm' in the state."))
    }
  }
  
  # do the nodeset names of the dependent variable and the origins have a
  # corresponding nodeset and are they of class 'nodeSet.realm'?
  nodesets <- unique(c(state[[dep.var]]$nodeSet,
                       unlist(lapply(origins, function(o) state[[o]]$nodeSet))))
  for (i in 1:length(nodesets)) {
    if (!(nodesets[i] %in% names(state))) {
      stop(paste(nodesets[i], "has no corresponding object.",
                 "Nodesets specified in the dependent variable and the origins",
                 "must have a corresponding object with the same name."))
    }
    if (!(inherits(state[[nodesets[i]]], "nodeSet.realm"))) {
      stop(paste(nodesets[i], "is not of class 'nodeSet.realm'."))
    }
  }
  
  # do dependent variable and individual nodeset have the same number of observations?
  resources <- resourceNodeSet(state)
  nResources <- length(state[[resources]]$ids)
  if (!(length(state[[dep.var]]$data) == nResources)) {
    stop("Dependent variable and individual nodeset do not have the same number of observations.")
  }
  if (max(state[[dep.var]]$data) > length(state[[destinationNodeSet(state)]]$ids)) {
    stop("Dependent variable contains locations that are not in its nodeset.")
  }
  
  # are the origins defined for the same resources, with locations of their nodeset?
  for (o in origins) {
    if (state[[o]]$nodeSet[2] != resources || length(state[[o]]$data) != nResources) {
      stop(paste("Origin", o, "must have one location for each of the", resources,
                 "of the dependent variable."))
    }
    if (max(state[[o]]$data) > length(state[[state[[o]]$nodeSet[1]]]$ids)) {
      stop(paste("Origin", o, "contains locations that are not in its nodeset."))
    }
  }
  
  # extract covar names
  covars <- names(state)[!(names(state) %in% nonCovariates(state))]
  
  # if there are covars, check them
  if (length(covars) != 0) {
    
    # are covars of class 'nodeVar.realm' or 'network.realm'? do they have correctly
    # specified nodeset names and are these nodesets of the correct size?
    for (i in 1:length(covars)) {
      #class
      if (!(inherits(state[[covars[i]]], "nodeVar.realm") || 
            inherits(state[[covars[i]]], "network.realm"))) {
        stop(paste(covars[i], "must either be of class 'nodeVar.realm' or 'network.realm'."))
      }
      # nodeset names
      if (!(all(state[[covars[i]]]$nodeSet %in% names(state)))) {
        stop(paste("The nodeset of covar '", covars[i], "' has no corresponding object.",
                   "Nodesets assigned when creating covars must have a corresponding object with the same name."           ))
      }
      # nodeset sizes
      if (!(state[[covars[i]]]$size[1] == length(state[[state[[covars[i]]]$nodeSet[1]]]$ids))) {
        stop(paste("The covar '", covars[i], "' is not of the same size as its assigned nodeset."))
      }
      if (inherits(state[[covars[i]]], "network.realm") &&
          !(state[[covars[i]]]$size[2] == length(state[[state[[covars[i]]]$nodeSet[2]]]$ids))) {
        stop(paste("The covar '", covars[i], "' is not of the same size as its assigned nodesets."))
      }
    }
  }
}


# createInternalCache
# The cache of the dependent variable: for each origin, the weighted network
# between the origin's locations (rows) and the destinations (columns), and
# the weighted resource networks (sums of each resource covariate); and for the
# destinations, the inflow and the sums of each resource covariate, used by
# effects that do not depend on an origin:
#   cache[[dep.var]]$destinations: inflow, resourceInflows
#   cache[[dep.var]][[origin]]: valuedNetwork, resourceNetworks, sameNodeSet
# sameNodeSet is TRUE when the origin has the same nodeset as the destinations,
# i.e. when a resource can stay at its origin (a loop).
createInternalCache <-
  function(processState,
           resourceCovariates = NULL) {
    cache <- list()
    name <- processState$dep.var
    if (!inherits(processState[[name]], "dependent.realm")) {
      stop(paste(name, "is not a dependent variable (see createDependent())."))
    }
    destinationSet <- destinationNodeSet(processState, name)
    nDestinations <- length(processState[[destinationSet]]$ids)
    nResources <- processState[[name]]$size[1]
    destinations <- processState[[name]]$data[seq_len(nResources)]

    # the values of the resource covariates; the values on each tie are summed
    # in the order of the resources, starting from 0, i.e. exactly as when
    # adding the resources one by one
    values <- lapply(resourceCovariates, function(ressCovar) {
      as.numeric(processState[[ressCovar]]$data[seq_len(nResources)])
    })
    names(values) <- resourceCovariates

    cache[[name]] <- list()
    cache[[name]]$destinations <- list(
      inflow = as.numeric(tabulate(destinations, nbins = nDestinations)),
      resourceInflows = lapply(values, function(v) {
        tieSumsCpp(as.integer(destinations), v, nDestinations)
      })
    )

    for (origin in processState$origins) {
      originSet <- processState[[origin]]$nodeSet[1]
      nOrigins <- length(processState[[originSet]]$ids)
      tie <- processState[[origin]]$data[seq_len(nResources)] +
        (destinations - 1) * nOrigins
      cache[[name]][[origin]] <- list(
        valuedNetwork = matrix(as.numeric(tabulate(tie, nbins = nOrigins * nDestinations)),
                               nOrigins, nDestinations),
        resourceNetworks = lapply(values, function(v) {
          matrix(tieSumsCpp(as.integer(tie), v, nOrigins * nDestinations),
                 nrow = nOrigins, ncol = nDestinations)
        }),
        sameNodeSet = originSet == destinationSet
      )
    }

    cache
  }


# getCovarianceMatrix
getCovarianceMatrix <- function(statistics) {
  meanStatistics <- colMeans(statistics)
  meanStatsMatrix <- meanStatistics %*% t(meanStatistics)

  nObs <- dim(statistics)[1]
  observationMatrices <-
    lapply(1:nObs, function(i) {
      statistics[i, ] %*% t(statistics[i, ])
    })
  meanObservationMatrix <- Reduce("+", observationMatrices) / nObs

  covMatrix <- meanObservationMatrix - meanStatsMatrix
  return(covMatrix)
}

# getNetStatsFromDeps
getNetStatsFromDeps <- function(dep.var, oneDep, ans, effects){
  depState <- ans$state
  depState[[dep.var]] <- oneDep
  
  resCovsInCache <- names(ans$cache[[dep.var]]$destinations$resourceInflows)
  if(is.null(resCovsInCache)){
    depCache <- createInternalCache(depState)
  } else {
    depCache <- createInternalCache(depState, resourceCovariates = resCovsInCache)
  }
  getNetworkStatistics(dep.var, depState, depCache, effects)
}

# effectOrigins
# the origin of each effect (NA for effects that do not depend on an origin)
effectOrigins <- function(effects) {
  vapply(effects$effectFormulas, function(f) {
    origin <- formals(f)$origin
    if (is.character(origin)) origin else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# effectSenders
# the origin of resource `edge` for each effect (NA for effects that do not
# depend on an origin), i.e. the i with which the effects are called
effectSenders <- function(state, origins, edge) {
  vapply(origins, function(origin) {
    if (is.na(origin)) NA_real_ else as.numeric(state[[origin]]$data[edge])
  }, numeric(1), USE.NAMES = FALSE)
}

# effectSenderTable
# effectSenders for all resources: a matrix with one row per resource and one
# column per effect
effectSenderTable <- function(state, origins) {
  nResources <- length(state[[state$dep.var]]$data)
  table <- matrix(NA_real_, nResources, length(origins))
  for (k in seq_along(origins)) {
    if (!is.na(origins[k])) {
      table[, k] <- as.numeric(state[[origins[k]]]$data)
    }
  }
  table
}

# getNetworkStatistics
# the target statistics: for effects with an origin, the sum of the target
# contributions over all cells (origin location i, destination j); for effects
# without, the sum over the destinations j
getNetworkStatistics <- function(dep.var, state, cache, effects) {
  if (useCppSampler(effects, state)) {
    return(cppNetworkStatistics(dep.var, state, cache, effects))
  }

  nDestinations <- length(state[[destinationNodeSet(state, dep.var)]]$ids)
  origins <- effectOrigins(effects)

  targetStatistics <- unlist(lapply(
    seq_along(effects$effectFormulas),
    function(k) {
      f <- effects$effectFormulas[[k]]
      if (is.na(origins[k])) {
        return(sum(vapply(1:nDestinations, function(j) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
            i = NA,
            j = j,
            update = NULL,
            edge = NULL,
            getTargetContribution = TRUE
          )
        }, numeric(1))))
      }
      nOrigins <- length(state[[state[[origins[k]]]$nodeSet[1]]]$ids)
      sum(apply(
        expand.grid(1:nOrigins, 1:nDestinations),
        1,
        function(v) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
            i = v[[1]],
            j = v[[2]],
            update = NULL,
            edge = NULL,
            getTargetContribution = TRUE
          )
        }
      ))
    }
  ))

  names(targetStatistics) <- effects$name

  return(targetStatistics)
}


# parallel computing
# One socket cluster per estimation (startCluster); the workers load the
# installed realm. Each call of applyChains runs one task per worker (chain)
# and gives every task its own L'Ecuyer-CMRG random number stream. The streams
# are derived from the random numbers of the master process, so results are
# reproducible with set.seed() and the chains are independent.
startCluster <- function(cpus) {
  cl <- makePSOCKcluster(cpus)
  tryCatch(clusterCall(cl, loadNamespace, "realm"), error = function(e) {
    stopCluster(cl)
    stop("realm could not be loaded on the parallel workers: ", conditionMessage(e))
  })
  cl
}

# n independent L'Ecuyer-CMRG streams, starting from a seed drawn with the
# random numbers of the master process (whose RNG kind is not changed)
chainSeeds <- function(n) {
  start <- sample.int(.Machine$integer.max, 1)
  oldSeed <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", oldSeed, envir = globalenv()))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(start)
  seeds <- list(get(".Random.seed", envir = globalenv()))
  for (k in seq_len(n - 1)) {
    seeds[[k + 1]] <- nextRNGStream(seeds[[k]])
  }
  seeds
}

runWithSeed <- function(.seed, .f, ...) {
  assign(".Random.seed", .seed, envir = globalenv())
  .f(...)
}

# run .f once per chain on the cluster cl; arguments in ... have one element
# per chain, arguments in MoreArgs are the same for all chains. The results
# are in the order of the chains.
applyChains <- function(cl, .f, ..., MoreArgs = list()) {
  clusterMap(cl, runWithSeed, .seed = chainSeeds(length(cl)), ...,
             MoreArgs = c(list(.f = .f), MoreArgs),
             SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

# iterations split as evenly as possible over n chains
splitIterations <- function(iterations, n) {
  perChain <- rep(iterations %/% n, n)
  rest <- iterations %% n
  if (rest > 0) {
    perChain[1:rest] <- perChain[1:rest] + 1
  }
  perChain
}


# runPhase1
runPhase1 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      initialParameters,
                      burnInN1,
                      iterationsN1,
                      thinningN1,
                      gainN1,
                      multinomialProposal = FALSE,
                      verbose = FALSE,
                      cl = NULL) {
  # simulate statistic matrix, with one chain per worker if cl is a cluster
  if (!is.null(cl)) {
    if (verbose) {
      cat(paste("Phase 1: \n burn-in", burnInN1, "steps\n", 
                iterationsN1, " iterations\n thinning", 
                thinningN1, "\n", length(cl), "cpus\n"))
    }
    statisticsMatrix <- applyChains(
      cl, simulateStatisticVectors,
      iterations = splitIterations(iterationsN1, length(cl)),
      MoreArgs = list(dep.var = dep.var, state = state, cache = cache,
                      effects = effects, initialParameters = initialParameters,
                      burnIn = burnInN1, thinning = thinningN1,
                      multinomialProposal = multinomialProposal)
    )
    statisticsMatrix <- Reduce("rbind", statisticsMatrix)
  } else {
  statisticsMatrix <-
    simulateStatisticVectors(
      dep.var,
      state,
      cache,
      effects,
      initialParameters,
      burnIn = burnInN1,
      iterations = iterationsN1,
      thinning = thinningN1,
      multinomialProposal = multinomialProposal,
      verbose = verbose
    )
  }

  # calculate covariance matric
  covarianceMatrix <- getCovarianceMatrix(statisticsMatrix)

  # calculate sensitivity
  sensitivityVector <- 1 / diag(covarianceMatrix)

  # initial update step if initial parameters are zero
  if (all(initialParameters == 0)) {
    observedStatistics <-
      getNetworkStatistics(dep.var, state, cache, effects)
    averageStatistics <- colMeans(statisticsMatrix)
    updatedParameters <-
      initialParameters - gainN1 * sensitivityVector * (averageStatistics - observedStatistics)
  } else {
    updatedParameters <- initialParameters
  }

  # returns initial estimates and sensitivity vector
  return(list(
    estimates = updatedParameters,
    sensitivityVector = sensitivityVector
  ))
}


# runPhase2
runPhase2 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      initialParameters,
                      sensitivityVector,
                      burnInN2,
                      nsubN2,
                      initGain,
                      thinningN2,
                      initialIterationsN2,
                      cl = NULL,
                      multinomialProposal = FALSE,
                      verbose = FALSE) {
  # calculate observed statistics
  observedStatistics <-
    getNetworkStatistics(dep.var, state, cache, effects)

  # with a cluster cl, one chain per worker; the parameters are averaged over
  # the chains after each sub phase
  parallel <- !is.null(cl)
  cpus <- if (parallel) length(cl) else 1

  # burn in
  if (parallel) {
    res <- applyChains(
      cl, simulateNSteps,
      n = rep(burnInN2, cpus),
      MoreArgs = list(dep.var = dep.var, state = state, cache = cache,
                      effects = effects, parameters = initialParameters,
                      multinomialProposal = multinomialProposal)
    )
  } else {
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        multinomialProposal = multinomialProposal,
        n = burnInN2
      )
    state <- res$state
    cache <- res$cache
  }

  # run n sub phases
  gain <- initGain
  iterations <- initialIterationsN2
  parameters <- initialParameters
  for (i in 1:nsubN2) {
    if (verbose) {
      cat(paste("\n"))
      cat(paste0("Sub phase", i, ":\n burn-in ", burnInN2, " steps\n ", 
                iterations, " iterations\n thinning ", 
                thinningN2, "\n ", cpus, " cpus\n"))
    }

    # with several chains, each chain continues from its own state and the
    # averaged parameters are passed to the next sub phase
    if (parallel) {
      res <- applyChains(
        cl, runSubphase2,
        state = lapply(res, "[[", "state"),
        cache = lapply(res, "[[", "cache"),
        MoreArgs = list(dep.var = dep.var, effects = effects,
                        initialParameters = parameters,
                        sensitivityVector = sensitivityVector,
                        observedStatistics = observedStatistics, gain = gain,
                        thinningN2 = thinningN2, iterations = iterations,
                        multinomialProposal = multinomialProposal,
                        verbose = FALSE)
      )
      parameterList <- lapply(res, "[[", "parameters")
      parameters <- colMeans(Reduce(rbind, parameterList))
      names(parameters) <- effects$name
    } else {
      resSub2 <- runSubphase2(
        dep.var,
        state,
        cache,
        effects,
        parameters,
        sensitivityVector,
        observedStatistics,
        gain,
        thinningN2,
        iterations,
        multinomialProposal = multinomialProposal,
        verbose = verbose
      )
      state <- resSub2$state
      cache <- resSub2$cache
      parameters <- resSub2$parameters
      names(parameters) <- effects$name
    }

    if (verbose) {
      cat(paste("\n"))
      cat(paste("New parameters:\n"))
    }
    if (verbose) {
      cat(paste(names(parameters), "\n", parameters, "\n"))
    }

    # determine number of iterations
    iterations <- round(iterations * 1.75)
    # determine gain
    gain <- gain / 2
  }

  return(parameters)
}


# runPhase2GT
# Alternative to runPhase2: Geyer-Thompson MCMC maximum likelihood estimation
# (Geyer & Thompson 1992; as in MCMC MLE for ERGMs). Networks are simulated at
# the current parameters theta0, and the parameters are updated with
# geyerThompsonStep. Iterations with iterationsGT draws continue until a step
# is an unrestricted Geyer-Thompson step, i.e. the observed statistics are
# within the range of the simulated ones; the estimate is then taken from a
# final step based on finalIterationsGT draws, which reduces Monte Carlo error.
runPhase2GT <- function(dep.var,
                        state,
                        cache,
                        effects,
                        initialParameters,
                        observedStatistics,
                        burnIn,
                        thinning,
                        iterationsGT,
                        maxIterationsGT,
                        finalIterationsGT = 3 * iterationsGT,
                        minESS = 0.5,
                        cl = NULL,
                        multinomialProposal = FALSE,
                        verbose = FALSE) {
  # with a cluster cl, one chain per worker
  cpus <- if (is.null(cl)) 1 else length(cl)

  runChain <- function(chain, dep.var, effects, parameters, burnIn, thinning,
                       multinomialProposal) {
    res <- simulateStatisticVectors(
      dep.var, chain$state, chain$cache, effects, parameters,
      burnIn = burnIn, iterations = chain$draws, thinning = thinning,
      multinomialProposal = multinomialProposal,
      returnState = TRUE
    )
    list(statistics = res$statistics, state = res$state, cache = res$cache)
  }

  # simulate a sample of statistics at the given parameters; chains continue
  # from where the previous sample ended
  chains <- lapply(seq_len(cpus), function(k) list(state = state, cache = cache))
  simulateSample <- function(parameters, n) {
    draws <- splitIterations(n, cpus)
    for (k in seq_len(cpus)) {
      chains[[k]]$draws <<- draws[k]
    }
    chainArgs <- list(dep.var = dep.var, effects = effects,
                      parameters = parameters, burnIn = burnIn,
                      thinning = thinning,
                      multinomialProposal = multinomialProposal)
    if (is.null(cl)) {
      chains <<- lapply(chains, function(chain) do.call(runChain, c(list(chain), chainArgs)))
    } else {
      chains <<- applyChains(cl, runChain, chain = chains, MoreArgs = chainArgs)
    }
    Reduce(rbind, lapply(chains, "[[", "statistics"))
  }

  report <- function(label, statistics, step, parameters) {
    if (!verbose) {
      return(invisible(NULL))
    }
    tRatios <- (colMeans(statistics) - observedStatistics) / apply(statistics, 2, sd)
    cat(paste0("\n", label, ": max |t-ratio| ", round(max(abs(tRatios)), 3),
               ", effective sample size ",
               if (is.na(step$ess)) "-" else round(step$ess), " of ", nrow(statistics),
               ", ", step$method, " step",
               if (step$restricted) " (restricted)" else "", "\n"))
    cat(paste(names(parameters), "\n", parameters, "\n"))
  }

  parameters <- initialParameters
  names(parameters) <- effects$name

  for (iter in 1:maxIterationsGT) {
    statistics <- simulateSample(parameters, iterationsGT)
    step <- geyerThompsonStep(statistics, observedStatistics, minESS)
    parameters <- parameters + step$delta
    report(paste("Geyer-Thompson iteration", iter), statistics, step, parameters)
    if (!step$restricted) {
      break
    }
  }

  # final step from a larger sample
  statistics <- simulateSample(parameters, finalIterationsGT)
  step <- geyerThompsonStep(statistics, observedStatistics, minESS)
  parameters <- parameters + step$delta
  report("Geyer-Thompson final step", statistics, step, parameters)

  parameters
}


# geyerThompsonStep
# One parameter update from a sample of statistics simulated at theta0.
# If the observed statistics lie within the range of the sample (measured by
# Mahalanobis distance), the Geyer-Thompson approximation of
# l(theta0 + delta) - l(theta0), -log mean_k exp(delta . (s_k - s_obs)), is
# maximised with Newton steps, keeping the effective sample size of the
# importance weights >= minESS * M. Otherwise, importance sampling is
# unreliable and a step based on the normal approximation of the statistics
# is taken towards a target between the sample mean and the observed
# statistics that lies within the sample range (similar to the step length
# of Hummel, Hunter & Handcock 2012).
geyerThompsonStep <- function(statistics, observedStatistics, minESS = 0.5,
                              maxNewton = 100, tol = 1e-8) {
  M <- nrow(statistics)
  p <- ncol(statistics)
  mu <- colMeans(statistics)
  C <- cov(statistics)
  solveC <- function(A, b) {
    tryCatch(solve(A, b), error = function(e) {
      solve(A + diag(1e-8 * max(1, max(diag(A))), p), b)
    })
  }

  # range of the sample: Mahalanobis distances of the draws from their mean
  centred <- sweep(statistics, 2, mu)
  sampleDistances <- sqrt(pmax(0, rowSums(centred * t(solveC(C, t(centred))))))
  range <- as.numeric(quantile(sampleDistances, 0.9))
  obsDistance <- sqrt(max(0, sum((observedStatistics - mu) *
                                   solveC(C, observedStatistics - mu))))

  if (obsDistance > range) {
    gamma <- range / obsDistance
    delta <- as.vector(solveC(C, gamma * (observedStatistics - mu)))
    return(list(delta = delta, ess = NA, restricted = TRUE, method = "normal"))
  }

  d <- sweep(statistics, 2, observedStatistics)
  weightsOf <- function(delta) {
    eta <- as.vector(d %*% delta)
    w <- exp(eta - max(eta))
    w / sum(w)
  }
  # objective, computed stably: -log mean exp(eta)
  objective <- function(delta) {
    eta <- as.vector(d %*% delta)
    -(max(eta) + log(mean(exp(eta - max(eta)))))
  }
  essOf <- function(w) 1 / sum(w^2)

  delta <- rep(0, p)
  restricted <- FALSE
  for (k in 1:maxNewton) {
    w <- weightsOf(delta)
    m <- colSums(w * d)
    Cw <- crossprod(d * sqrt(w)) - tcrossprod(m)
    newtonStep <- -solveC(Cw, m)

    # backtrack until the objective does not decrease and the effective
    # sample size is large enough
    a <- 1
    repeat {
      candidate <- delta + a * newtonStep
      if (objective(candidate) >= objective(delta) - 1e-12 &&
          essOf(weightsOf(candidate)) >= minESS * M) {
        break
      }
      a <- a / 2
      if (a < 1e-6) {
        candidate <- delta
        break
      }
    }
    if (a < 1) {
      restricted <- restricted ||
        essOf(weightsOf(candidate)) < minESS * M * 1.01
    }
    change <- max(abs(candidate - delta))
    delta <- candidate
    if (change < tol) {
      break
    }
  }

  list(delta = delta, ess = essOf(weightsOf(delta)), restricted = restricted,
       method = "Geyer-Thompson")
}


# getMPLE
# Maximum pseudo-likelihood estimate: a conditional logit of each resource's
# observed destination against all possible destinations, given the
# destinations of all other resources. Used as starting values for the
# estimation. Returns NULL if the estimate does not exist (e.g. separation) or
# the optimisation fails. The pseudo-likelihood is calculated in C++ if all
# effects have a C++ version (mpleObjectiveCpp), otherwise in R.
getMPLE <- function(dep.var, state, cache, effects,
                    maxNewton = 50, tol = 1e-8) {
  p <- length(effects$effectFormulas)
  evaluate <- if (useCppSampler(effects, state)) {
    mpleObjectiveCpp(dep.var, state, cache, effects)
  } else {
    mpleObjectiveR(dep.var, state, cache, effects)
  }

  theta <- rep(0, p)
  for (k in 1:maxNewton) {
    current <- evaluate(theta, derivatives = TRUE)
    if (!is.finite(current$logPL)) {
      return(NULL)
    }
    step <- tryCatch(-solve(current$hessian, current$gradient), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      return(NULL)
    }
    # step halving to keep the pseudo-likelihood increasing
    a <- 1
    repeat {
      candidate <- evaluate(theta + a * step, derivatives = FALSE)$logPL
      if ((is.finite(candidate) && candidate >= current$logPL - 1e-10) || a <= 1e-6) {
        break
      }
      a <- a / 2
    }
    theta <- theta + a * step
    if (max(abs(a * step)) < tol) {
      break
    }
  }
  if (any(!is.finite(theta)) || max(abs(theta)) > 1e3) {
    return(NULL)
  }
  names(theta) <- effects$name
  theta
}

# mpleObjectiveR
# returns a function of the parameters that gives the log pseudo-likelihood
# and, with derivatives = TRUE, its gradient and Hessian, using the R effect
# functions (see getMPLE); the change statistics are calculated once and stored
mpleObjectiveR <- function(dep.var, state, cache, effects) {
  state <- toMutableState(state, dep.var)
  cache <- toMutableCache(cache)
  dependent <- state[[dep.var]]
  nodeName <- destinationNodeSet(state, dep.var)
  resourceName <- resourceNodeSet(state, dep.var)
  destinations <- state[[nodeName]]$ids[state[[nodeName]]$considerWhenSampling]
  resources <- which(state[[resourceName]]$considerWhenSampling)
  effectFunctions <- effects$effectFormulas
  origins <- effectOrigins(effects)
  p <- length(effectFunctions)

  # for each resource: change statistics of adding it to every possible
  # destination after removing it, and which destination is observed
  choiceSets <- lapply(resources, function(e) {
    senders <- effectSenders(state, origins, e)
    j <- dependent$data[e]

    cells <- getCacheCells(cache[[dep.var]], e, j, state)
    updateCacheInPlace(cache[[dep.var]], e, j, state = state, update = -1)
    setEdgeInPlace(dependent, e, NA)
    X <- t(vapply(destinations, function(k) {
      vapply(seq_along(effectFunctions), function(f) {
        effectFunctions[[f]](dep.var = dep.var, state = state, cache = cache,
                             i = senders[f], j = k, update = 1, edge = e)
      }, numeric(1))
    }, numeric(p)))
    restoreCacheCells(cache[[dep.var]], e, j, state, cells)
    setEdgeInPlace(dependent, e, j)

    list(X = matrix(X, ncol = p), observed = match(j, destinations))
  })

  function(theta, derivatives = TRUE) {
    logPL <- sum(vapply(choiceSets, function(cs) {
      eta <- as.vector(cs$X %*% theta)
      eta[cs$observed] - max(eta) - log(sum(exp(eta - max(eta))))
    }, numeric(1)))
    if (!derivatives || !is.finite(logPL)) {
      return(list(logPL = logPL))
    }
    gradient <- rep(0, p)
    hessian <- matrix(0, p, p)
    for (cs in choiceSets) {
      eta <- as.vector(cs$X %*% theta)
      prob <- exp(eta - max(eta))
      prob <- prob / sum(prob)
      expected <- colSums(prob * cs$X)
      gradient <- gradient + cs$X[cs$observed, ] - expected
      hessian <- hessian - (crossprod(cs$X * sqrt(prob)) - tcrossprod(expected))
    }
    list(logPL = logPL, gradient = gradient, hessian = hessian)
  }
}


# runPhase3
runPhase3 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      parameters,
                      observedStatistics,
                      iterationsN3,
                      burnInN3,
                      thinningN3,
                      cl = NULL,
                      verbose = FALSE,
                      returnDeps = FALSE,
                      multinomialProposal = FALSE,
                      fish = fish) {
  # simulate statistic matrix
  # with a cluster cl, one chain per worker; the results are combined in the
  # order of the chains
  if (!is.null(cl)) {
    statsA <- applyChains(
      cl, simulateStatisticVectors,
      iterations = splitIterations(iterationsN3, length(cl)),
      MoreArgs = list(dep.var = dep.var, state = state, cache = cache,
                      effects = effects, initialParameters = parameters,
                      burnIn = burnInN3, thinning = thinningN3,
                      verbose = verbose,
                      returnDeps = returnDeps,
                      multinomialProposal = multinomialProposal, fish = FALSE)
    )


    if (returnDeps) {
      stats <- list()
      stats[[1]] <- lapply(statsA, function(x) {
        x[[1]]
      })
      stats[[2]] <-
        unlist(lapply(statsA, function(x) {
          x[[2]]
        }), recursive = FALSE)
      statisticsMatrix <- Reduce("rbind", stats[[1]])
    } else {
      stats <- statsA
      statisticsMatrix <- Reduce("rbind", stats)
    }
  } else {
    stats <- simulateStatisticVectors(
      dep.var,
      state,
      cache,
      effects,
      parameters,
      burnIn = burnInN3,
      iterations = iterationsN3,
      thinning = thinningN3,
      verbose = verbose,
      returnDeps = returnDeps,
      multinomialProposal = multinomialProposal,
      fish = fish
    )
    if (returnDeps) {
      statisticsMatrix <- stats[[1]]
    } else {
      statisticsMatrix <- stats
    }
  } # no parallel computing

  # calculate convergence statistics
  simulatedMeans <- colMeans(statisticsMatrix)
  simulatedSDs <- apply(statisticsMatrix, 2, sd)
  convergenceStatistic <-
    (simulatedMeans - observedStatistics) / simulatedSDs

  # calculate covariance matrix
  covarianceMatrix <- getCovarianceMatrix(statisticsMatrix)

  # calculate covariance matrix of parameters
  covarianceMatrixParameters <- solve(covarianceMatrix)

  # create empty object for the deps
  deps <- NULL
  
  # extract only the dependent var to return from deps
  if (returnDeps) {
    deps <- lapply(stats[[2]], function(x){
      x$state[[dep.var]]
    })
  }

  # returns covariance matrix and convergence statistics, and the statistics of
  # the simulated networks (in the order of deps)
  return(
    list(
      covarianceMatrix = covarianceMatrixParameters,
      convergenceStatistics = convergenceStatistic,
      deps = deps,
      statistics = statisticsMatrix
    )
  )
}


# runSubphase2
runSubphase2 <- function(dep.var,
                         state,
                         cache,
                         effects,
                         initialParameters,
                         sensitivityVector,
                         observedStatistics,
                         gain,
                         thinningN2,
                         iterations,
                         multinomialProposal = FALSE,
                         verbose = FALSE) {
  parameters <- c()

  for (i in 1:iterations) {
    if (verbose) {
      cat(".")
    }
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        multinomialProposal = multinomialProposal,
        n = (thinningN2 + 1)
      )
    state <- res$state
    cache <- res$cache

    # calculate statistics
    statistics <-
      getNetworkStatistics(dep.var, state, cache, effects)

    # calculate and store new parameters
    updatedParameters <-
      initialParameters - gain * sensitivityVector * (statistics - observedStatistics)
    initialParameters <- updatedParameters
    parameters <- rbind(parameters, updatedParameters)
  }
  if (verbose) {
    cat("\n")
  }

  ### return mean of the parameters discarding the 
  ### initial dk that are clearly before convergence
  # find the smaller of two numbers: 
  # 25% of the number of iterations OR
  # when the expected residual of the step-size is less than 0.05
  dk1 <- round(iterations/4)
  
  dk2 <- 1
  res <- (1-gain)
  
  repeat{
    if(dk2 > 100) break
    if(any(res < 0.05)) break
    dk2 <- dk2 + 1
    res <- res*(1-gain)
  }
  
  dk <- min(dk1, dk2)
  
  if(ncol(parameters) == 1){
    return(list(
      parameters = mean(parameters[dk:iterations,]),
      state = state,
      cache = cache
    ))
  }
  
  return(list(
    parameters = colMeans(parameters[dk:iterations,]),
    state = state,
    cache = cache
  ))
}


# simulateNSteps
# returns state and cache in their in-place form (see toMutableState);
# convert with toListState / toListCache before storing them
simulateNSteps <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           multinomialProposal = FALSE,
           debug = FALSE,
           n = 1) {
    if (!debug && useCppSampler(effects, state)) {
      res <- simulateChain(dep.var, state, cache, effects, parameters,
                           burnIn = n, iterations = 0, thinning = 0,
                           multinomialProposal = multinomialProposal)
      return(list(state = res$state, cache = res$cache))
    }

    state <- toMutableState(state, dep.var)
    cache <- toMutableCache(cache)
    # the origins of the resources do not change: their senders for each
    # effect are looked up in a table
    senderTable <- effectSenderTable(state, effectOrigins(effects))
    for (i in 1:n) {
      res <-
        simulateOneStep(
          dep.var,
          state,
          cache,
          effects,
          parameters,
          multinomialProposal,
          debug,
          senderTable = senderTable
        )
      state <- res$state
      cache <- res$cache
    }
    return(res)
  }


# simulateOneStep
# changes state and cache in place; both must be in their in-place form.
# A step moves one resource to another destination; its origins do not change.
simulateOneStep <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           multinomialProposal = FALSE,
           debug = FALSE,
           senderTable = effectSenderTable(state, effectOrigins(effects))) {
    resourceName <- resourceNodeSet(state, dep.var)
    nodeName <- destinationNodeSet(state, dep.var)
    resourceSet <-
      which(state[[resourceName]]$considerWhenSampling)

    # sample edges to swap
    randomEdge <- sample(resourceSet, size = 1)
    j <- state[[dep.var]]$data[randomEdge]
    # the origin of the resource for each effect (NA for effects that do not
    # depend on an origin)
    senders <- senderTable[randomEdge, ]

    nodeSet <-
      state[[nodeName]]$ids[state[[nodeName]]$considerWhenSampling]

    effectFunctions <- effects$effectFormulas
    nEffects <- length(effectFunctions)
    # the change statistics of all effects for moving the resource to or from
    # destination k
    changeStatistics <- function(k, update) {
      statistics <- numeric(nEffects)
      for (e in seq_len(nEffects)) {
        statistics[e] <- effectFunctions[[e]](
          dep.var = dep.var,
          state = state,
          cache = cache,
          i = senders[e],
          j = k,
          update = update,
          edge = randomEdge
        )
      }
      statistics
    }

    ## MULTINOMIAL PROPOSAL ## CONSIDER REFACTOR TO FUNCTION

    if (multinomialProposal) {
      if (debug) {
        print(paste("Proposed multinomial change of", randomEdge, "from", j))
      }

      # remove resource from process state
      statisticsDrop <- changeStatistics(j, -1)

      # update process state and cache
      updateCacheInPlace(cache[[dep.var]], randomEdge, j, state = state, update = -1)
      setEdgeInPlace(state[[dep.var]], randomEdge, NA)

      # for each possible receiver; one row per effect
      statisticsCreate <- matrix(
        vapply(nodeSet, function(k) changeStatistics(k, 1),
               numeric(length(effectFunctions))),
        nrow = length(effectFunctions)
      )
      objValues <-
        colSums((statisticsCreate + statisticsDrop) * parameters)
      pChange <- exp(objValues) / sum(exp(objValues))

      # draw receiver node
      draw <- which(cumsum(pChange) >= runif(1))[1]
      l <- nodeSet[draw]

      # update process state and cache
      updateCacheInPlace(cache[[dep.var]], randomEdge, l, state = state, update = 1)
      setEdgeInPlace(state[[dep.var]], randomEdge, l)
    }


    ## BINOMIAL PROPOSAL ## CONSIDER REFACTOR TO FUNCTION

    if (!multinomialProposal) {
      l <- sample(nodeSet, size = 1)

      if (debug) {
        print(paste("Proposed change of", randomEdge, "from", j, "to", l))
      }

      if (j == l) {
        return(list(state = state, cache = cache))
      }

      statisticsDrop <- changeStatistics(j, -1)

      # remove the resource from state and cache, keeping what is needed to
      # undo this if the change is rejected
      oldCells <- getCacheCells(cache[[dep.var]], randomEdge, j, state)
      updateCacheInPlace(cache[[dep.var]], randomEdge, j, state = state, update = -1)
      setEdgeInPlace(state[[dep.var]], randomEdge, NA)

      statisticsCreate <- changeStatistics(l, 1)

      changeContribution <-
        sum((statisticsDrop + statisticsCreate) * parameters)

      acceptanceProbability <- min(1, exp(changeContribution))

      if (acceptanceProbability >= runif(1)) {
        updateCacheInPlace(cache[[dep.var]], randomEdge, l, state = state, update = 1)
        setEdgeInPlace(state[[dep.var]], randomEdge, l)
      } else {
        restoreCacheCells(cache[[dep.var]], randomEdge, j, state, oldCells)
        setEdgeInPlace(state[[dep.var]], randomEdge, j)
      }
    }

    return(list(state = state, cache = cache))
  }


# simulateStatisticVectors
simulateStatisticVectors <- function(dep.var,
                                     state,
                                     cache,
                                     effects,
                                     initialParameters,
                                     burnIn,
                                     iterations,
                                     thinning,
                                     verbose = FALSE,
                                     multinomialProposal = FALSE,
                                     fish = FALSE,
                                     returnDeps = FALSE,
                                     returnState = FALSE) {
  # the whole chain in C++ if all effects have a C++ version
  if (useCppSampler(effects, state)) {
    if (verbose) {
      cat(paste("Starting burn-in with", burnIn, "steps\n"))
    }
    res <- simulateChain(dep.var, state, cache, effects, initialParameters,
                         burnIn = burnIn, iterations = iterations, thinning = thinning,
                         multinomialProposal = multinomialProposal,
                         returnDraws = returnDeps)
    statisticsMatrix <- res$statistics
    rownames(statisticsMatrix) <- rep("stats", nrow(statisticsMatrix))
    if (returnDeps) {
      return(list(statisticsMatrix, res$draws))
    }
    if (returnState) {
      return(list(statistics = statisticsMatrix, state = res$state, cache = res$cache))
    }
    return(statisticsMatrix)
  }

  statisticsMatrix <- c()

  if (returnDeps) {
    deps <- list()
  }

  # burn in
  if (verbose) {
    cat(paste("Starting burn-in with", burnIn, "steps\n"))
  }
  res <-
    simulateNSteps(
      dep.var,
      state,
      cache,
      effects,
      initialParameters,
      multinomialProposal = multinomialProposal,
      n = burnIn
    )
  cache <- res$cache
  state <- res$state

  # thinning and calculate statistics
  for (i in 1:iterations) {
    if (fish) {
      int <- i %% 11
      s <- "_,.-'``'-.,"
      cat(substr(s, int + 1, int + 1))
      if (runif(1) < 0.02) {
        cat("><(((A*>")
      }
    } else if (verbose) {
      cat(".")
    }
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        multinomialProposal = multinomialProposal,
        n = thinning + 1
      )
    state <- res$state
    cache <- res$cache

    # calculate and save statistics
    stats <-
      getNetworkStatistics(dep.var,
        state = state,
        cache = cache,
        effects = effects
      )
    statisticsMatrix <- rbind(statisticsMatrix, stats)

    if (returnDeps) {
      deps[[i]] <- list(state = toListState(state, dep.var),
                        cache = toListCache(cache))
    }
  }

  if (fish) {
    cat("\n")
  } else if (verbose) {
    cat("\n")
  }


  if (returnDeps) {
    return(list(statisticsMatrix, deps))
  } else if (returnState) {
    # the end of the chain (in-place form), to continue it later
    return(list(statistics = statisticsMatrix, state = state, cache = cache))
  } else {
    return(statisticsMatrix)
  }
}


# updateWeightedCache
# adds (update = 1) or removes (update = -1) resource `edge` at destination
# `receiver` in the cache of a dependent variable (as a list), i.e. in the
# networks of all origins and in the destinations
updateWeightedCache <- function(cacheDepVar,
                                edge,
                                receiver,
                                state,
                                update) {
  for (origin in names(cacheDepVar)) {
    x <- cacheDepVar[[origin]]
    if (origin == "destinations") {
      for (ressCovar in names(x$resourceInflows)) {
        x$resourceInflows[[ressCovar]][receiver] <-
          x$resourceInflows[[ressCovar]][receiver] + update * state[[ressCovar]]$data[edge]
      }
      x$inflow[receiver] <- x$inflow[receiver] + update
    } else {
      sender <- state[[origin]]$data[edge]
      # Updates of weighted resource caches
      for (ressCovar in names(x$resourceNetworks)) {
        v <- x$resourceNetworks[[ressCovar]][sender, receiver]
        x$resourceNetworks[[ressCovar]][sender, receiver] <-
          v + update * state[[ressCovar]]$data[edge]
      }
      # weighted updates
      x$valuedNetwork[sender, receiver] <- x$valuedNetwork[sender, receiver] + update
    }
    cacheDepVar[[origin]] <- x
  }
  cacheDepVar
}


########## in-place state and cache for MCMC chains
# Within one chain, the dependent variable of the state and the cache are
# stored in environments and updated in place. This avoids copying the
# destinations and all cache matrices in every simulation step. Effects read
# them as before (state[[dep.var]]$data, cache[[dep.var]][[origin]]$valuedNetwork[i, j],
# ...): the cache of a dependent variable is a list with one environment for
# the destinations and one for each origin.
# Updates are done by functions whose environment is that environment, using
# <<-; e$x[i, j] <- v would copy x, and evalq() is not byte-compiled.
# Anything that leaves a chain is converted back with toListState/toListCache.

# toMutableState
toMutableState <- function(state, dep.var) {
  if (is.environment(state[[dep.var]])) {
    return(state)
  }
  e <- list2env(unclass(state[[dep.var]]), envir = new.env(parent = baseenv()))
  e$.names <- names(state[[dep.var]])
  e$.class <- class(state[[dep.var]])

  # placeholder for the variable the function below updates in e (its
  # environment is set to e below)
  data <- NULL
  e$.setEdge <- function(edge, value) {
    data[edge] <<- value
  }
  environment(e$.setEdge) <- e

  state[[dep.var]] <- e
  state
}

# toListState
toListState <- function(state, dep.var) {
  e <- state[[dep.var]]
  if (!is.environment(e)) {
    return(state)
  }
  l <- mget(e$.names, envir = e)
  class(l) <- e$.class
  state[[dep.var]] <- l
  state
}

# toMutableCache
toMutableCache <- function(cache) {
  lapply(cache, function(cacheDepVar) {
    if (is.environment(cacheDepVar$destinations)) {
      return(cacheDepVar)
    }
    out <- lapply(names(cacheDepVar), function(name) {
      if (name == "destinations") {
        mutableDestinations(cacheDepVar[[name]])
      } else {
        mutableOrigin(cacheDepVar[[name]])
      }
    })
    names(out) <- names(cacheDepVar)
    out
  })
}

# mutableOrigin
# the cache of one origin as an environment, with functions that update it in
# place (same logic as updateWeightedCache)
mutableOrigin <- function(x) {
  e <- list2env(x, envir = new.env(parent = baseenv()))
  e$.names <- names(x)
  if (is.null(x$resourceNetworks)) {
    e$resourceNetworks <- NULL
  }

  # placeholders for the variables the functions below update in e (their
  # environment is set to e below)
  valuedNetwork <- resourceNetworks <- NULL

  e$.update <- function(sender, receiver, resourceID, state, update) {
    # Updates of weighted resource caches
    for (ressCovar in names(resourceNetworks)) {
      resourceNetworks[[ressCovar]][sender, receiver] <<-
        resourceNetworks[[ressCovar]][sender, receiver] +
        update * state[[ressCovar]]$data[resourceID]
    }

    # weighted updates
    valuedNetwork[sender, receiver] <<- valuedNetwork[sender, receiver] + update
  }

  # the cells that an update of sender-receiver changes, to restore them
  # exactly if a proposed change is rejected (adding back a numeric resource
  # covariate value is not always exact in floating point)
  e$.getCells <- function(sender, receiver) {
    list(
      valued = valuedNetwork[sender, receiver],
      resource = lapply(resourceNetworks, function(m) m[sender, receiver])
    )
  }

  e$.restoreCells <- function(sender, receiver, cells) {
    for (ressCovar in names(resourceNetworks)) {
      resourceNetworks[[ressCovar]][sender, receiver] <<- cells$resource[[ressCovar]]
    }
    valuedNetwork[sender, receiver] <<- cells$valued
  }

  environment(e$.update) <- e
  environment(e$.getCells) <- e
  environment(e$.restoreCells) <- e
  e
}

# mutableDestinations
# the cache of the destinations as an environment (as mutableOrigin)
mutableDestinations <- function(x) {
  e <- list2env(x, envir = new.env(parent = baseenv()))
  e$.names <- names(x)
  if (is.null(x$resourceInflows)) {
    e$resourceInflows <- NULL
  }

  inflow <- resourceInflows <- NULL

  e$.update <- function(receiver, resourceID, state, update) {
    for (ressCovar in names(resourceInflows)) {
      resourceInflows[[ressCovar]][receiver] <<-
        resourceInflows[[ressCovar]][receiver] +
        update * state[[ressCovar]]$data[resourceID]
    }
    inflow[receiver] <<- inflow[receiver] + update
  }

  e$.getCells <- function(receiver) {
    list(
      inflow = inflow[receiver],
      resource = lapply(resourceInflows, function(v) v[receiver])
    )
  }

  e$.restoreCells <- function(receiver, cells) {
    for (ressCovar in names(resourceInflows)) {
      resourceInflows[[ressCovar]][receiver] <<- cells$resource[[ressCovar]]
    }
    inflow[receiver] <<- cells$inflow
  }

  environment(e$.update) <- e
  environment(e$.getCells) <- e
  environment(e$.restoreCells) <- e
  e
}

# toListCache
toListCache <- function(cache) {
  lapply(cache, function(cacheDepVar) {
    lapply(cacheDepVar, function(e) {
      if (!is.environment(e)) {
        return(e)
      }
      mget(e$.names, envir = e)
    })
  })
}

# updateCacheInPlace
# same as updateWeightedCache, but changes the cache environments in place
updateCacheInPlace <- function(cacheDepVar, edge, receiver, state, update) {
  for (name in names(cacheDepVar)) {
    if (name == "destinations") {
      cacheDepVar[[name]]$.update(receiver, edge, state, update)
    } else {
      cacheDepVar[[name]]$.update(state[[name]]$data[edge], receiver, edge, state, update)
    }
  }
  invisible(NULL)
}

# getCacheCells
# the cells that adding or removing resource `edge` at destination `receiver`
# changes, for restoreCacheCells
getCacheCells <- function(cacheDepVar, edge, receiver, state) {
  cells <- list()
  for (name in names(cacheDepVar)) {
    cells[[name]] <- if (name == "destinations") {
      cacheDepVar[[name]]$.getCells(receiver)
    } else {
      cacheDepVar[[name]]$.getCells(state[[name]]$data[edge], receiver)
    }
  }
  cells
}

# restoreCacheCells
restoreCacheCells <- function(cacheDepVar, edge, receiver, state, cells) {
  for (name in names(cacheDepVar)) {
    if (name == "destinations") {
      cacheDepVar[[name]]$.restoreCells(receiver, cells[[name]])
    } else {
      cacheDepVar[[name]]$.restoreCells(state[[name]]$data[edge], receiver, cells[[name]])
    }
  }
  invisible(NULL)
}

# setEdgeInPlace
# sets the destination of one resource (NA while it is moved) in the
# environment of the dependent variable
setEdgeInPlace <- function(dependentEnv, edge, value) {
  dependentEnv$.setEdge(edge, value)
  invisible(NULL)
}
