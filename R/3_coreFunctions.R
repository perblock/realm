####### coreFunctions


#' createAlgorithm
#'
#' Specifies the algorithm used in the estimation based on characteristics
#' of the state and the effects.
#'
#' @param state An object of class "processState.realm" that contains all relevant information about
#' the outcome (the destinations), the origins, the nodesets, and covariates.
#' @param effects An object of class "effectsList.realm" that specifies the model.
#' @param multinomialProposal How should the next possible outcome in the simulation chains
#' be sampled? If TRUE, fewer simulation steps are needed, but each simulation
#' step takes considerably longer. Defaults to FALSE.
#' @param burnInN1 The number of simulation steps before the first draw in Phase 1.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' @param thinningN1 The number of simulation steps between two draws in Phase 1.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' @param iterationsN1 The number of draws taken in Phase 1.
#' A recommended value is at least 4 * n_effects which is set as default.
#' If the value is too low, there will be an error in Phase 1.
#' @param gainN1 The size of the updating step after Phase 1. A conservative
#' value is 0, values higher than 0.25 are courageous. Defaults to 0.1.
#' Only used if the estimation starts from zero, i.e., if no starting values
#' are given and maximum pseudo-likelihood estimates are not available.
#' @param burnInN2 The number of simulation steps before the first draw in Phase 1.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' @param thinningN2 The number of simulation steps between two draws in Phase 2.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' @param initialIterationsN2 The number of draws taken in subphase 1 of Phase 2.
#' For first estimations, a recommended value is around 50 (default to 50).
#' Note that in later subphases, the number of iterations increases.
#' If this is a further estimation to improve convergence, higher values (100+)
#' are recommended.
#' @param nsubN2 Number of subphases in Phase 2. In case this is the first
#' estimation, 4 subphases are recommended and set as default. If convergence
#' in a previous estimation was close, then 1-2 subphases should be enough.
#' @param initGain The magnitude of parameter updates in the first subphase of
#' Phase 2. Values of around 0.2 (default) are recommended.
#' @param burnInN3 The number of simulation steps before the first draw in Phase 3.
#' A recommended value is at least 3 * n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least 3 * n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' @param thinningN3 The number of simulation steps between two draws in Phase 3.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = FALSE, and at least 2 * n_Individuals if multinomialProposal = TRUE
#' which is set as default.
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param iterationsN3 Number of draws in Phase 3. Recommended are at the very
#' least 500 (default).
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param phase2Method How the parameters are estimated in Phase 2:
#' "RM" (default) uses the Robbins-Monro stochastic approximation
#' algorithm; "GT" uses Geyer-Thompson MCMC maximum likelihood estimation,
#' as used for exponential random graph models. With "GT", the arguments
#' nsubN2, initGain and initialIterationsN2 are not used.
#' @param iterationsGT With phase2Method = "GT": the number of draws in each
#' Geyer-Thompson iteration. Defaults to 10 * n_effects, but at least 100.
#' @param maxIterationsGT With phase2Method = "GT": the maximum number of
#' Geyer-Thompson iterations before the final step. Iterations stop earlier
#' once the observed statistics lie within the range of the simulated ones.
#' Defaults to 10.
#' @param finalIterationsGT With phase2Method = "GT": the number of draws for
#' the final Geyer-Thompson step, which determines the estimates.
#' Defaults to 3 * iterationsGT.
#'
#' @return An object of class "algorithm.realm".
#' @export
#'
#' @seealso [createProcessState()], [createEffectsObject()], [estimateMobilityNetwork()]
#'
#' @examples
#' # define algorithm based on state and effects characteristics
#' myAlg <- createAlgorithm(myState, myEffects, multinomialProposal = FALSE)
createAlgorithm <-
  function(state,
           effects,
           multinomialProposal = FALSE,
           burnInN1 = NULL,
           thinningN1 = NULL,
           iterationsN1 = NULL,
           gainN1 = 0.1,
           burnInN2 = NULL,
           thinningN2 = NULL,
           initialIterationsN2 = 50,
           nsubN2 = 4,
           initGain = 0.6,
           burnInN3 = NULL,
           thinningN3 = NULL,
           iterationsN3 = 500,
           phase2Method = c("RM", "GT"),
           iterationsGT = NULL,
           maxIterationsGT = 10,
           finalIterationsGT = NULL) {
    dep.var <- state$dep.var
    
    algorithm <- list()
    
    nameNodeSet1 <- destinationNodeSet(state, dep.var)
    nameNodeSet2 <- resourceNodeSet(state, dep.var)
    
    algorithm[["multinomialProposal"]] <- as.logical(multinomialProposal)

    if (is.null(burnInN1) & multinomialProposal == FALSE) {
      algorithm[["burnInN1"]] <- 
        length(state[[nameNodeSet1]][["ids"]]) * length(state[[nameNodeSet2]][["ids"]])
    }
    if (is.null(burnInN1) & multinomialProposal == TRUE) {
      algorithm[["burnInN1"]] <- length(state[[nameNodeSet2]][["ids"]])
    }
    if (!(is.null(burnInN1))) {
      algorithm[["burnInN1"]] <- burnInN1
    }

    if (is.null(iterationsN1)) {
      algorithm[["iterationsN1"]] <- length(effects[["effectFormulas"]]) * 8
    } else {
      algorithm[["iterationsN1"]] <- iterationsN1
    }

    if (is.null(thinningN1) & multinomialProposal == FALSE) {
      algorithm[["thinningN1"]] <- 
        length(state[[nameNodeSet2]][["ids"]]) * length(state[[nameNodeSet1]][["ids"]]) * 0.5
    }
    if (is.null(thinningN1) & multinomialProposal == TRUE) {
      algorithm[["thinningN1"]] <- length(state[[nameNodeSet2]][["ids"]])
    }
    if (!(is.null(thinningN1))) {
      algorithm[["thinningN1"]] <- thinningN1
    }

    algorithm[["gainN1"]] <- gainN1

    if (is.null(burnInN2) & multinomialProposal == FALSE) {
      algorithm[["burnInN2"]] <- 
        length(state[[nameNodeSet2]][["ids"]]) * length(state[[nameNodeSet1]][["ids"]])
    }
    if (is.null(burnInN2) & multinomialProposal == TRUE) {
      algorithm[["burnInN2"]] <- length(state[[nameNodeSet2]][["ids"]])
    }
    if (!(is.null(burnInN2))) {
      algorithm[["burnInN2"]] <- burnInN2
    }

    algorithm[["nsubN2"]] <- nsubN2

    algorithm[["initGain"]] <- initGain

    if (is.null(thinningN2) & multinomialProposal == FALSE) {
      algorithm[["thinningN2"]] <- 
        length(state[[nameNodeSet2]][["ids"]]) * length(state[[nameNodeSet1]][["ids"]]) * 0.5
    }
    if (is.null(thinningN2) & multinomialProposal == TRUE) {
      algorithm[["thinningN2"]] <- length(state[[nameNodeSet2]][["ids"]])
    }
    if (!(is.null(thinningN2))) {
      algorithm[["thinningN2"]] <- thinningN2
    }

    algorithm[["initialIterationsN2"]] <- initialIterationsN2

    algorithm[["iterationsN3"]] <- iterationsN3

    if (is.null(burnInN3) & multinomialProposal == FALSE) {
      algorithm[["burnInN3"]] <- 
        length(state[[nameNodeSet2]][["ids"]]) * length(state[[nameNodeSet1]][["ids"]]) * 3
    }
    if (is.null(burnInN3) & multinomialProposal == TRUE) {
      algorithm[["burnInN3"]] <- length(state[[nameNodeSet2]][["ids"]]) * 3
    }
    if (!(is.null(burnInN3))) {
      algorithm[["burnInN3"]] <- burnInN3
    }

    if (is.null(thinningN3) & multinomialProposal == FALSE) {
      algorithm[["thinningN3"]] <- 
        length(state[[nameNodeSet2]][["ids"]]) * length(state[[nameNodeSet1]][["ids"]])
    }
    if (is.null(thinningN3) & multinomialProposal == TRUE) {
      algorithm[["thinningN3"]] <- length(state[[nameNodeSet2]][["ids"]]) * 2
    }
    if (!(is.null(thinningN3))) {
      algorithm[["thinningN3"]] <- thinningN3
    }

    algorithm[["phase2Method"]] <- match.arg(phase2Method)
    if (is.null(iterationsGT)) {
      algorithm[["iterationsGT"]] <- max(100, 10 * length(effects[["effectFormulas"]]))
    } else {
      algorithm[["iterationsGT"]] <- iterationsGT
    }
    algorithm[["maxIterationsGT"]] <- maxIterationsGT
    if (is.null(finalIterationsGT)) {
      algorithm[["finalIterationsGT"]] <- 3 * algorithm[["iterationsGT"]]
    } else {
      algorithm[["finalIterationsGT"]] <- finalIterationsGT
    }

    class(algorithm) <- "algorithm.realm"
    algorithm
  }

#' realmAlgorithmCreate
#'
#' @rdname createAlgorithm
realmAlgorithmCreate <- createAlgorithm

#' createDependent
#'
#' Creates the dependent variable, the outcome to be modelled by realm: the
#' destination of each person/resource. The origins of the resources are
#' separate objects (see [createOrigin()]).
#'
#' @param dest A vector with the destination of each person/resource, numbered
#' from one to the number of locations.
#' @param nodeSet The nodesets of the dependent variable. This is a vector with
#' two entries referencing the names of the nodesets of locations and
#' individuals of the form c(locations, individuals).
#' @param nodes Alternative way to specify the nodeSet by naming nodes and edges: 
#' nodes denote the locations (destinations)
#' @param edges Alternative way to specify the nodeSet by naming nodes and edges: 
#' edges denote the individuals
#' 
#' @return An object of class "dependent.realm".
#' @export
#'
#' @seealso [createOrigin()], [createProcessState()]
#'
#' @examples
#' # the destinations of the example data
#' transfers <- createDependent(mobilityEdgelist[, 2], c("organisations", "people"))
createDependent <-
  function(dest, nodeSet = NULL, nodes = NULL, edges = NULL) {
    if (!is.null(dim(dest)) && ncol(dest) == 2) {
      stop(paste("The dependent variable contains only the destinations.",
                 "Create it from the second column of the edge list and the",
                 "origins from the first column with createOrigin()."))
    }
    locationVector(dest, nodeSet, nodes, edges, "dependent.realm")
  }

#' realmDependent
#'
#' @rdname createDependent
realmDependent <- createDependent


#' createOrigin
#'
#' Creates an origin: the location each person/resource comes from. A process
#' state can contain several origins, also of different nodesets of locations;
#' effects that depend on the origin name the origin they use.
#'
#' @param orig A vector with the origin of each person/resource, numbered
#' from one to the number of locations.
#' @param nodeSet The nodesets of the origin. This is a vector with two entries
#' referencing the names of the nodesets of locations and individuals of the
#' form c(locations, individuals).
#' @param nodes Alternative way to specify the nodeSet by naming nodes and edges: 
#' nodes denote the locations (origins)
#' @param edges Alternative way to specify the nodeSet by naming nodes and edges: 
#' edges denote the individuals
#' 
#' @return An object of class "origin.realm".
#' @export
#'
#' @seealso [createDependent()], [createProcessState()]
#'
#' @examples
#' # the origins of the example data
#' origin <- createOrigin(mobilityEdgelist[, 1], c("organisations", "people"))
createOrigin <-
  function(orig, nodeSet = NULL, nodes = NULL, edges = NULL) {
    if (!is.null(dim(orig)) && ncol(orig) != 1) {
      stop("An origin is a vector with one location per person/resource.")
    }
    locationVector(orig, nodeSet, nodes, edges, "origin.realm")
  }

#' realmOrigin
#'
#' @rdname createOrigin
realmOrigin <- createOrigin


# locationVector
# the dependent variable (destinations) and origins: one location per
# resource, with the nodesets c(locations, resources)
locationVector <- function(x, nodeSet, nodes, edges, class) {
  x <- as.vector(x)
  if (!is.numeric(x)) {
    x <- as.numeric(x)
  }
  if (any(is.na(x))) {
    stop(paste("Input data includes missing values or cannot be classified as numeric."))
  }
  if (any(x < 1) || any(x != round(x))) {
    stop("Locations should be numbered from one to the number of locations.")
  }
  # update nodeSet name according to which the user specified
  if (is.null(nodeSet)) {
    if (is.null(nodes) || is.null(edges)) {
      stop("Either nodeSet or nodes and edges need to be specified")
    }
    nodeSet <- c(nodes, edges)
  }
  if (length(nodeSet) != 2) {
    stop("Two nodesets need to be specified: nodes (locations) / edges (individuals)")
  }
  if (!inherits(nodeSet, "character")) {
    stop("nodeSet or nodes and edges need to be class character")
  }
  l <- list(
    data = x,
    nodeSet = nodeSet,
    size = length(x)
  )
  class(l) <- class
  l
}


#' createEffectsObject
#'
#' Specifies the model with endogenous and exogenous predictors.
#' The predictors in the model are called “effects”.
#'
#' @param effectInit A list of "effects", where each effect to be included is specified as a
#' further list that contains the effect name and the additional parameters it needs.
#' Effects without further parameters only contain the effect name (e.g., loops).
#' @param checkProcessState For internal use only.
#'
#' @return An object of class "effectsList.realm".
#' @export
#'
#' @examples
#' # create an effects object
#' myEffects <- createEffectsObject(
#'   list(
#'     list("loops", origin = "origin"),
#'     list("present_relations", origin = "origin"),
#'     list("dyadic_covariate", origin = "origin", attribute.index = "sameRegion"),
#'     list("alter_covariate", attribute.index = "size"),
#'     list("resource_covar_to_node_covar",
#'       attribute.index = "region",
#'       resource.attribute.index = "sex"
#'     ),
#'     list("loops_resource_covar", origin = "origin", resource.attribute.index = "sex")
#'   )
#' )
createEffectsObject <-
  function(effectInit, checkProcessState = NULL) {
    # TODO add default parameters to effect names
    effectNames <-
      lapply(effectInit, function(x) {
        ifelse(is.list(x), do.call(paste, x), x)
      })

    effects <-
      lapply(effectInit, function(x) {
        eval(parse(text = x[[1]]))
      })

    # Update signatures of the effects based on default parameters and above specified parameters
    signatures <- lapply(effects, formals)
    setParams <- lapply(effectInit, function(x) {
      x[-1]
    })

    # effects that depend on the origin of the resources need one, the others
    # cannot have one
    for (i in seq_along(effects)) {
      effectName <- if (is.list(effectInit[[i]])) effectInit[[i]][[1]] else effectInit[[i]]
      if (effectName %in% names(deprecatedEffects)) {
        warning(deprecatedEffects[[effectName]], call. = FALSE)
      }
      needsOrigin <- "origin" %in% names(signatures[[i]])
      origin <- if (is.list(setParams[[i]])) setParams[[i]]$origin else NULL
      if (needsOrigin && (!is.character(origin) || length(origin) != 1)) {
        stop(paste0("Effect ", effectNames[[i]], " depends on the origin of the resources; ",
                    "specify the origin, e.g. origin = \"origin\"."))
      }
      if (!needsOrigin && !is.null(origin)) {
        stop(paste0("Effect ", effectNames[[i]], " does not depend on an origin; ",
                    "remove the origin argument."))
      }
    }
    signatures <- mapply(
      function(s, p) {
        s[names(p)] <- p
        s["cache"] <- alist(cache = NULL)
        s
      },
      signatures, setParams,
      SIMPLIFY = FALSE
    )
    if ("matrix" %in% class(signatures)) {
      signatures <- apply(signatures, 2, invisible)
    }

    # Assign signatures with default values to generic functions
    for (i in 1:length(effects)) {
      formals(effects[[i]]) <- unlist(signatures[i], recursive = FALSE)
    }

    effectslist <- list(
      effectFormulas = effects,
      name = unlist(effectNames)
    )
    class(effectslist) <- "effectsList.realm"

    # checks of the data that effects refer to (see checkEffectInputs)
    if (!is.null(checkProcessState)) {
      checkEffectInputs(effectslist, checkProcessState)
    }

    effectslist
  }


#' createNetwork
#'
#' Defines a network between locations, generally to be used as a predictor in the model.
#' NOTE: The outcome variable of the model is not defined as a network, but as
#' the destinations of the resources (see [createDependent()])!
#'
#' @param m A matrix containing the network data: square for a network between
#' the locations of one nodeset, otherwise with the locations of the first
#' nodeset in rows and of the second in columns.
#' @param isSymmetric Currently not in use.
#' @param isBipartite Currently not in use.
#' @param nodeSet Which nodeset are the nodes of the network. Usually this will
#' be the locations in the data. Two nodesets, e.g. c(origin locations,
#' destination locations), for a network between the locations of two nodesets.
#' @param nodes Alternative way to specify the nodeSet by naming nodes: 
#' nodes denote the locations
#'
#' @return An object of class "network.realm".
#' @export
#'
#' @seealso [createProcessState()], [createDependent()]
#'
#' @examples
#' # create an object of class network.realm
#' sameRegion <- outer(orgRegion, orgRegion, "==") * 1
#' sameRegion <- createNetwork(sameRegion, nodeSet = c("organisations", "organisations"))
createNetwork <-
  function(m,
           isSymmetric = FALSE,
           isBipartite = FALSE,
           nodeSet = NULL,
           nodes = NULL) {
    if (!is.matrix(m)) {
      stop("Not a matrix.")
    }
    if (!is.numeric(m)) {
      m <- matrix(as.numeric(m), ncol = ncol(m))
    }
    if (any(is.na(m))) {
      stop(paste("Input matrix includes missing values or cannot be classified as numeric."))
    }
    # update nodeSet name according to which the user specified
    if(is.null(nodeSet)){
      if(is.null(nodes)){
        stop("Either nodeSet or nodes need to be specified")
      }
      nodeSet <- nodes
    }
    if (length(nodeSet) == 1) {
      nodeSet <- c(nodeSet, nodeSet)
    }
    # a network between the locations of one nodeset is square; between two
    # nodesets (e.g. origins in regions and destinations in organisations) it
    # has the first nodeset in rows and the second in columns
    if (nodeSet[1] == nodeSet[2] && nrow(m) != ncol(m)) {
      stop("Input matrix should have the same number of rows as columns.")
    }
    if(!inherits(nodeSet, "character")){
      stop("nodeSet or nodes need to be class character")
    }
    l <- list(
      data = m,
      isSymmetric = isSymmetric,
      isBipartite = isBipartite,
      nodeSet = nodeSet,
      size = dim(m)
    )
    class(l) <- "network.realm"
    l
  }

#' dyadicCovar
#' 
#' @rdname createNetwork
dyadicCovar <- createNetwork


#' createNodeSet
#'
#' Determines and names the nodesets of individuals and locations that make up the mobility network.
#'
#' @param x Either a single number indicating how many items are in this nodeset
#' or a vector from 1:n_items.
#' @param isPresent Currently not in use.
#' @param considerWhenSampling A boolean/logical vector of the length of the nodeset.
#' Only in use in special cases.
#' If the nodeset indicates a location, considerWhenSampling indicates whether the
#' location is a possible destination, or is only an origin (e.g. a training facility).
#' Entries in the vector of locations that cannot be a destination are FALSE.
#' If the nodeset indicates mobile individuals, considerWhenSampling indicates whether
#' their mobility should be modelled or whether it is structurally determined, that
#' is, their mobility is exogenously defined and does not follow the same logic as
#' the mobility of everybody else.
#'
#' @return An object of class "nodeSet.realm".
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # create an object of class nodeSet.realm
#' people <- createNodeSet(1:nrow(mobilityEdgelist))
#' organisations <- createNodeSet(length(orgRegion))
createNodeSet <-
  function(x = NULL,
           isPresent = NULL,
           considerWhenSampling = NULL) {
    if (is.null(x) && is.null(isPresent)) {
      stop("Only null parameters.")
    }
    if (is.null(x)) {
      x <- length(isPresent)
    }
    if (length(x) == 1 &&
      !is.numeric(x)) {
      stop("non-numeric value of single valued x.")
    }
    ids <- if (length(x) == 1) {
      1:x
    } else {
      x
    }
    if (is.null(isPresent)) {
      isPresent <- rep(T, length(ids))
    }
    if (is.null(considerWhenSampling)) {
      considerWhenSampling <- rep(T, length(ids))
    }
    if (length(isPresent) != length(ids)) {
      stop("Wrong length of presence vector.")
    }
    l <- list(
      isPresent = isPresent,
      ids = ids,
      considerWhenSampling = considerWhenSampling
    )
    class(l) <- "nodeSet.realm"
    l
  }

#' realmEdges
#' 
#' @rdname createNodeSet
realmEdges <- createNodeSet


#' realmNodes
#' 
#' @rdname createNodeSet
realmNodes <- createNodeSet


#' createNodeVariable
#'
#' Assigns a covariate to one nodeset, i.e., an exogenous characteristic of mobile
#' individuals/resources or locations.
#'
#' @param values A vector assigning the covariate value to each element of the nodeset.
#' @param range The range of possible values, currently not in use.
#' @param nodeSet The nodeset to which the covariate applies.
#' @param addSame Will the variable be used to model categorical homophily (e.g.,
#' with the same_covariate effect)? In this case, addSame needs to be set to TRUE.
#' @param addSim Will the variable be used to model continuous homophily (e.g.,
#' with the sim_covariate effect)? In this case, addSim needs to be set to TRUE.
#' @param nodes Alternative way to specify the nodeSet by naming nodes or edges: 
#' nodes denote the locations
#' @param edges Alternative way to specify the nodeSet by naming nodes or edges: 
#' edges denote the individuals
#' 
#' @return An object of class "nodeVar.realm".
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # create an object of class nodeVar.realm
#' region <- createNodeVariable(orgRegion, nodeSet = "organisations")
#' size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
#' sex <- createNodeVariable(indSex, nodeSet = "people")
createNodeVariable <-
  function(values,
           range = NULL,
           nodeSet = NULL,
           nodes = NULL,
           edges = NULL,
           addSame = NULL,
           addSim = NULL) {
    if (!(is.vector(values))) {
      stop("Input data should be of class 'vector'.")
    }
    if (!is.numeric(values)) {
      values <- as.numeric(values)
    }
    if (any(is.na(values))) {
      stop(paste("Input vector includes missing values or cannot be classified as numeric."))
    }
    
    # update nodeSet name according to which the user specified
    if(is.null(nodeSet)){
      if(is.null(nodes)){
        if(is.null(edges)){
          stop("Either nodeSet or nodes and edges need to be specified")
        } else {nodeSet <- edges}
      } else {nodeSet <- nodes}
    }
    
    if(!inherits(nodeSet, "character")){
      stop("nodeSet, nodes, or edges (whichever is used) need to be class character")
    }
    l <- list(
      data = values,
      range = range,
      nodeSet = nodeSet,
      size = length(values)
    )
    
    # add same and sim for nodesets that refer to loactions, usually 
    # smaller nodesets (i.e. < 500)
    
    if(is.null(addSame)){
      if(length(values) < 501){
        addSame <- TRUE
      } else {
        addSame <- FALSE
      }
    }
    
    if(is.null(addSim)){
      if(length(values) < 501){
        addSim <- TRUE
      } else {
        addSim <- FALSE
      }
    }

    if (addSame) {
      l[["same"]] <- outer(values, values, "==") * 1
    }
    if (addSim) {
      l[["sim"]] <-
        1 - abs(outer(values, values, "-")) / (range(values)[2] - range(values)[1])
    }

    class(l) <- "nodeVar.realm"
    l
  }

#' monadicCovar
#' 
#' @rdname createNodeVariable
monadicCovar <- createNodeVariable


#' createProcessState
#'
#' Creates the "Process state", i.e., a realm object that stores all information
#' about the data that will be used in the estimation. This includes the
#' outcome variable (the destinations), the origins, the nodesets, and all
#' covariates.
#'
#' @param elements A named list of the outcome variable (destinations), the
#' origins, the nodesets, and all covariates that contain the information about
#' the data that will be used in the estimation.
#' @param dependentVariable The name of the outcome variable (destinations) as
#' specified under "elements". This indicates what outcome
#' the researcher is interested in.
#'
#' @return An object of class "processState.realm". All elements of class
#' "origin.realm" are origins; their names are stored in the element `origins`.
#' @export
#'
#' @seealso [createDependent()], [createOrigin()], [createNodeSet()],
#' [createNodeVariable()], [createNetwork()]
#'
#' @examples
#' # Create a process state out of the mobility data objects:
#' # create objects (which are later combined to the process state)
#' transfers <- createDependent(mobilityEdgelist[, 2],
#'   nodeSet = c("organisations", "people")
#' )
#' origin <- createOrigin(mobilityEdgelist[, 1],
#'   nodeSet = c("organisations", "people")
#' )
#' people <- createNodeSet(1:nrow(mobilityEdgelist))
#' organisations <- createNodeSet(1:length(orgRegion))
#' sameRegion <- outer(orgRegion, orgRegion, "==") * 1
#' sameRegion <- createNetwork(sameRegion,
#'   nodeSet = c("organisations", "organisations")
#' )
#' region <- createNodeVariable(orgRegion, nodeSet = "organisations")
#' size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
#' sex <- createNodeVariable(indSex, nodeSet = "people")
#'
#' # combine created objects to the process state
#' myState <- createProcessState(list(
#'     transfers = transfers,
#'     origin = origin,
#'     people = people,
#'     organisations = organisations,
#'     sameRegion = sameRegion,
#'     region = region,
#'     size = size,
#'     sex = sex),
#'   dependentVariable = "transfers")
createProcessState <- function(elements, dependentVariable) {
  if (!is.list(elements)) {
    stop("Expecting a list.")
  }
  if (length(elements) == 0) {
    stop("List must have at least one element.")
  }
  if (is.null(names(elements))) {
    stop("List must be named.")
  }
  if (anyDuplicated(names(elements))) {
    stop(paste("The names of the elements must be unique; used more than once:",
               paste(unique(names(elements)[duplicated(names(elements))]), collapse = ", ")))
  }

  nodeSetIDs <- c()
  linkedElementIDs <- c()
  sizes <- c()
  for (i in 1:length(elements)) {
    e <- elements[[i]]
    if (!(
      class(e) %in% c(
        "dependent.realm",
        "origin.realm",
        "nodeSet.realm",
        "nodeVar.realm",
        "network.realm"
      )
    )) {
      stop(paste0("Unknown element of class '", class(e), 
           "'. Input objects should either be of classes 'dependent.realm', 'origin.realm', ",
           "'nodeSet.realm', 'nodeVar.realm or 'network.realm'."))
    }

    # TODO CLEANUP from here. What is necessary, hat should be extended?

    if (is(e, "nodeSet.realm")) {
      nodeSetIDs <- c(nodeSetIDs, i)
    }
    if (class(e) %in% c("network.realm", "nodeVar.realm")) {
      linkedElementIDs <- c(linkedElementIDs, i)
      sizes <- c(sizes, e$size)
    }
  }

  # if no node sets were found, create a default
  # TODO: the node set check should be based on the nodeSet names, not their size => done by checkProcessState?
  if (length(nodeSetIDs) == 0) {
    if (length(unique(sizes)) != 1) {
      stop("Differing element sizes without defining node sets.")
    }
    elements <-
      append(elements, list(actors = createNodeSet(sizes[1])))
  }

  # check if all elements in 'linkedElementIDs' have a corresponding node set
  # TODO implement => done by checkProcessState?
  
  elements$dep.var <- dependentVariable
  # the origins, in the order of the elements
  elements$origins <- names(elements)[vapply(elements, inherits, logical(1), "origin.realm")]

  class(elements) <- "processState.realm"
  
  checkProcessState(elements) # checks whether all input objects are correctly specified/valid
  
  
  # define resource covariates: the covariates of the resources (origins,
  # which are also defined for the resources, are not covariates)
  covars <- names(elements)[!(names(elements) %in% nonCovariates(elements))]
  resourceCovariates <- c()
  if(length(covars) > 0){
    for (i in 1:length(covars)) {
      if (resourceNodeSet(elements) %in% elements[[covars[i]]]$nodeSet) {
        resourceCovariates[length(resourceCovariates) + 1] <- covars[i]
      }
    }
  }
    
  # attach cache object
  elements$cache <- createInternalCache(elements, resourceCovariates = resourceCovariates)
  
  
  elements
}


#' createWeightedCache
#'
#' Since MoNAn version 1.0.0, this function no longer exists.
#'
#' @param processState Outdated.
#' @param resourceCovariates Outdated.
#'
#' @return Outdated.
#' @export
createWeightedCache <-
  function(processState,
           resourceCovariates = NULL) {
    cat(paste0("This function is not needed anymore since MoNAn version 1.0.0.", "\n", 
               "The cache is automatically created and included in the process state."))
  }


#' estimateMobilityNetwork
#'
#' The core function of the package in which the model for the analysis of
#' mobility tables is estimated.
#'
#' @param state An object of class "processState.realm" which contains all relevant information about
#' the outcome (the destinations), the origins, the nodesets, and covariates.
#' @param effects An object of class "effectsList.realm" that specifies the model.
#' @param algorithm An object of class "algorithm.realm" which specifies the algorithm used in the estimation.
#' @param initialParameters Starting values for the parameters. If neither
#' initialParameters nor prevAns are specified, maximum pseudo-likelihood
#' estimates are used, i.e., the estimates of the multinomial logit model for
#' which [getMultinomialStatistics()] provides the data. If these do not exist
#' (e.g., because of separation), the estimation starts from zero. For their
#' calculation, the change statistics of all possible destinations are stored,
#' once for each group of resources with the same origin, destination and
#' values of the resource covariates used by the effects. If this needs more
#' than `getOption("realm.mpleMemoryLimit")` bytes (default 1e9, i.e., 1 GB),
#' they are recalculated in each step instead, which is slower but needs little
#' memory.
#' @param prevAns If a previous estimation did not yield satisfactory convergence,
#' the outcome object of that estimation should be specified here to provide new
#' starting values for the estimation.
#' @param parallel Logical: computation on multiple cores? If TRUE, the
#' simulations run in `cpus` independent chains on a socket cluster (package
#' 'parallel'). Each chain gets its own random number stream derived from the
#' random numbers of the R session, so results are reproducible with
#' [set.seed()]; they differ from those of a serial estimation with the same
#' seed. The workers use the installed version of realm.
#' @param cpus Number of cores for computation in case parallel = TRUE.
#' @param verbose Logical: display information about estimation progress in the console?
#' @param returnDeps Logical: should the simulated values of Phase 3 be stored and returned?
#' This is necessary to run GoF tests.
#' Note that this might result in very large objects.
#' @param fish Logical: display a fish?
#' @param saveAlg Specify whether the algorithm object should be saved in the
#' results object. Defaults to FALSE.
#' @param cache Outdated parameter, no need to specify.
#'
#' @return The function `estimateMobilityNetwork` returns an object of class "result.realm" that contains the estimates, standard errors,
#' and convergence statistics. Furthermore, the covariance matrix used to calculate
#' the standard errors is included, which also shows collinearity between effects.
#' In case returnDeps = TRUE, the simulations of Phase 3 are included, too,
#' together with the statistics of the model's effects in the observed and the
#' simulated networks (used by [scoreTest()]).
#' @export
#'
#' @seealso [createProcessState()], [createEffectsObject()], [createAlgorithm()]
#'
#' @examples
#' \donttest{
#' # estimate mobility network model
#' 
#' myAlg_short <- createAlgorithm(myState, myEffects, multinomialProposal = FALSE,
#'                                nsubN2 = 1, iterationsN3 = 100)
#' 
#' myResDN <- estimateMobilityNetwork(myState, myEffects, myAlg_short,
#'                                    initialParameters = NULL,
#'                                    # in case a pseudo-likelihood estimation was run, replace with
#'                                    # initialParameters = initEst,
#'                                    parallel = TRUE, cpus = 2,
#'                                    verbose = TRUE,
#'                                    returnDeps = TRUE,
#'                                    fish = FALSE
#' )
#' 
#' # check convergence
#' max(abs(myResDN$convergenceStatistics))
#' 
#' # view results
#' myResDN
#' }
estimateMobilityNetwork <-
  function(state,
           effects,
           algorithm,
           initialParameters = NULL,
           prevAns = NULL,
           parallel = FALSE,
           cpus = 1,
           verbose = FALSE,
           returnDeps = FALSE,
           fish = FALSE,
           saveAlg = TRUE,
           cache = NULL) {
    
    if(!is.null(cache)){
      warning(paste("The cache object is automatically included in the process state", 
              "since MoNAn version 1.0.0 and does not need to be specified anymore."))
    }
    if(is.null(state$cache)){
      stop(paste("The cache object is automatically included in the process state", 
                    "since MoNAn version 1.0.0. Please recreate your process state."))
    }
    
    dep.var <- state$dep.var
    cache <- state$cache
    checkEffectInputs(effects, state)
    
    # extract parameters from algorithm object
    
    multinomialProposal <- algorithm$multinomialProposal
    burnInN1 <- algorithm$burnInN1
    iterationsN1 <- algorithm$iterationsN1
    thinningN1 <- algorithm$thinningN1
    gainN1 <- algorithm$gainN1
    burnInN2 <- algorithm$burnInN2
    nsubN2 <- algorithm$nsubN2
    initGain <- algorithm$initGain
    thinningN2 <- algorithm$thinningN2
    initialIterationsN2 <- algorithm$initialIterationsN2
    iterationsN3 <- algorithm$iterationsN3
    burnInN3 <- algorithm$burnInN3
    thinningN3 <- algorithm$thinningN3

    # set parameters to default values if not defined explicitly
    if (is.null(initialParameters)) {
      initialParameters <- rep(0, length(effects$name))
    }

    observedStatistics <-
      getNetworkStatistics(dep.var, state, cache, effects)

    # check wether the prevAns object has the same effects, if not cause error
    if (!is.null(prevAns)) {
      if (!(length(effects$name) == length(names(prevAns$estimates)))) {
        stop("prevAns object has different number of effects as the effects object")
      }
      if (!all(effects$name == names(prevAns$estimates))) {
        stop("effect names of prevAns object and effects object do not match")
      }
    }

    # starting values: if neither prevAns nor initialParameters are given,
    # use the maximum pseudo-likelihood estimate
    useGT <- identical(algorithm$phase2Method, "GT")
    if (is.null(prevAns) && all(initialParameters == 0)) {
      if (verbose) {
        cat("Calculating maximum pseudo-likelihood estimates as starting values\n")
      }
      mple <- getMPLE(dep.var, state, cache, effects)
      if (is.null(mple)) {
        cat("Maximum pseudo-likelihood estimates not available; starting from zero\n")
      } else {
        initialParameters <- mple
      }
    }

    # one cluster for all phases of a parallel estimation
    cl <- NULL
    if (parallel && cpus > 1) {
      cl <- startCluster(cpus)
      on.exit(stopCluster(cl), add = TRUE)
    }

    # decide whether to run phase 1: it is skipped with prevAns, and with
    # Geyer-Thompson estimation if there are starting values. With starting
    # values, phase 1 only determines the sensitivity vector for Robbins-Monro.
    if (!is.null(prevAns)) {
      # skip phase 1 as all information is already there
      cat("skipping phase 1 and taking values from prevAns\n")
      initialParameters <- prevAns$estimates
      sensitivityVector <- 1 / diag(solve(prevAns$covarianceMatrix))
    } else if (useGT && any(initialParameters != 0)) {
      # Geyer-Thompson does not need a sensitivity vector
    } else {
      # run phase 1 to get initial estimates and a sensitivity vector
      if (verbose) {
        cat("Starting phase 1\n")
      }
      resPhase1 <-
        runPhase1(
          dep.var,
          state,
          cache,
          effects,
          initialParameters,
          burnInN1,
          iterationsN1,
          thinningN1,
          gainN1,
          multinomialProposal,
          verbose,
          cl = cl
        )
      initialParameters <- resPhase1$estimates
      sensitivityVector <- resPhase1$sensitivityVector
    }

    # run phase 2 to get parameter estimates
    if (verbose) {
      cat("Starting phase 2\n")
    }
    # algorithm objects created before phase2Method existed use Robbins-Monro
    if (identical(algorithm$phase2Method, "GT")) {
      resPhase2 <- runPhase2GT(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        observedStatistics,
        burnIn = burnInN2,
        thinning = thinningN2,
        iterationsGT = algorithm$iterationsGT,
        maxIterationsGT = algorithm$maxIterationsGT,
        finalIterationsGT = algorithm$finalIterationsGT,
        cl = cl,
        multinomialProposal = multinomialProposal,
        verbose = verbose
      )
    } else {
      resPhase2 <- runPhase2(
        dep.var,
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
        cl = cl,
        multinomialProposal,
        verbose = verbose
      )
    }

    # run phase 3 to get final covariance matrix and convergence check
    if (verbose) {
      cat("Starting phase 3:\n")
      cat(paste(" burn-in", burnInN3, "steps\n", 
                iterationsN3, " iterations\n thinning", 
                thinningN3, "\n", cpus, "cpus\n"))
    }
    resPhase3 <- runPhase3(
      dep.var,
      state,
      cache,
      effects,
      resPhase2,
      observedStatistics,
      algorithm$iterationsN3,
      algorithm$burnInN3,
      algorithm$thinningN3,
      cl = cl,
      verbose = verbose,
      returnDeps = returnDeps,
      algorithm$multinomialProposal,
      fish = fish
    )

    if (!returnDeps) {
      result <- list(
        estimates = resPhase2,
        standardErrors = sqrt(diag(resPhase3$covarianceMatrix)),
        covarianceMatrix = resPhase3$covarianceMatrix,
        convergenceStatistics = resPhase3$convergenceStatistics,
        state = state,
        cache = cache
      )
    }

    if (returnDeps) {
      result <- list(
        estimates = resPhase2,
        standardErrors = sqrt(diag(resPhase3$covarianceMatrix)),
        covarianceMatrix = resPhase3$covarianceMatrix,
        convergenceStatistics = resPhase3$convergenceStatistics,
        state = state,
        cache = cache,
        deps = resPhase3$deps,
        # statistics of the model's effects in the observed and simulated
        # networks, used by scoreTest
        observedStatistics = observedStatistics,
        simulatedStatistics = `rownames<-`(resPhase3$statistics, NULL)
      )
    }
    
    if (saveAlg) {
      result[["algorithm"]] <- algorithm
    }

    class(result) <- "result.realm"

    return(result)
  }


#' estimateDistributionNetwork
#'
#' @rdname estimateMobilityNetwork
estimateDistributionNetwork <- estimateMobilityNetwork


#' realmEstimate
#'
#' @rdname estimateMobilityNetwork
realmEstimate <- estimateMobilityNetwork


#' simulateMobilityNetworks
#'
#' Simulates mobility networks for given data, effects, and parameters. This
#' function is mainly interesting to explore the behavior of the model or to
#' do counter-factual simulations.
#'
#' @param state An object of class "processState.realm" that contains all relevant information about
#' nodesets, and covariates. Further, the destinations of the dependent variable
#' need to be specified with the initial mobility network as starting value for the
#' simulation. For a large enough burn-in, any initial mobility network is allowed.
#' @param effects An object of class "effectsList.realm" that specifies the model.
#' @param parameters The parameters associated with the effects that shall be used
#' in the simulations.
#' @param burnin The number of simulation steps that are taken before the first draw of a
#' network is taken. A number too small will mean the first draw is influenced
#' by the initially specified network. A recommended value for the lower bound is 3 * n_Individuals *
#' n_locations.
#' @param thinning The number of simulation steps that are taken between two draws of a
#' network. A recommended value for the lower bound is n_Individuals * n_locations.
#' @param nSimulations The number of mobility networks to be simulated.
#' @param cache Outdated parameter, no need to specify.
#'
#' @return An object of class "sims.realm" with nSimulations entries, where each entry contains a further list with the
#' state and the cache of the current simulation stored.
#' @export
#'
#' @examples
#' \donttest{
#' # simulate a mobility network
#' # note that thinning and burn-in values are for this example only
#' # in real cases, choose values aprrox. times 10
#' mySimDN <- simulateMobilityNetworks(
#'   myState,
#'   myEffects,
#'   parameters = c(2, 1, 0.5, 0.1, -1, -0.5),
#'   burnin = 450,
#'   thinning = 150,
#'   nSimulations = 10
#' )
#'
#' mySimDN[[1]]
#' }
simulateMobilityNetworks <-
  function(state,
           effects,
           parameters,
           burnin,
           thinning,
           nSimulations,
           cache = NULL) {
    
    if(!is.null(cache)){
      warning(paste("The cache object is automatically included in the process state", 
                    "since MoNAn version 1.0.0 and does not need to be specified anymore."))
    }
    if(is.null(state$cache)){
      stop(paste("The cache object is automatically included in the process state", 
                    "since MoNAn version 1.0.0. Please recreate your process state."))
    }
    
    dep.var <- state$dep.var
    cache <- state$cache
    checkEffectInputs(effects, state)
    if (length(parameters) != length(effects$name)) {
      stop(paste("parameters must have one value per effect:", length(effects$name),
                 "effects, but", length(parameters), "parameters"))
    }
    
    # generate a state and cache after burnin with parameters
    r <-
      simulateNSteps(
        dep.var = dep.var,
        state = state,
        cache = cache,
        effects = effects,
        parameters = parameters,
        n = burnin
      )

    # extract state and cache
    simState <- r$state
    simCache <- r$cache

    # generate list in which simulated states and caches will be stored
    simulatedList <- list()

    for (i in 1:nSimulations) {
      r <-
        simulateNSteps(
          dep.var,
          simState,
          simCache,
          effects,
          n = thinning,
          parameters = parameters
        )
      simState <- r$state
      simCache <- r$cache
      simulatedList[[i]] <- list(state = toListState(simState, dep.var),
                                 cache = toListCache(simCache))
      int <- i %% 11
      s <- "_,.-'``'-.,"
      cat(substr(s, int + 1, int + 1))
      if (runif(1) < 0.02) {
        cat("<o_><")
      }
    }
    class(simulatedList) <- "sims.realm"
    return(simulatedList)
  }

#' simulateDistributionNetworks
#'
#' @rdname simulateMobilityNetworks
simulateDistributionNetworks <- simulateMobilityNetworks

#' realmSimulate
#'
#' @rdname simulateMobilityNetworks
realmSimulate <- simulateMobilityNetworks

