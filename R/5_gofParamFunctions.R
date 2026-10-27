########## gofParamFunctions


#' getIndegree
#' 
#' Calculates the weighted indegree distribution of all locations in the 
#' network. 
#' The weighted indegree is simply the column sum of the mobility table.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values. Level k is the
#' number of locations with a weighted indegree of k - 1, i.e., level 1 counts
#' locations with indegree 0. If NULL, all levels from 1 to the maximum
#' indegree + 1 are used.
#' @param ... Additional parameters.
#'
#' @return The degree distribution of simulated networks, named by the indegree
#' each count refers to.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'  
#' @keywords internal
getIndegree <- function(cache, dep.var, lvls, ...) {
  indegrees <- cache[[dep.var]]$destinations$inflow
  if (is.null(lvls)) {
    lvls <- 1:(max(indegrees) + 1)
  }
  counts <- tabulate(indegrees + 1, nbins = max(c(lvls, indegrees + 1)))
  setNames(counts[lvls], lvls - 1)
}


#' getTieWeights
#' 
#' Extracts the distribution of tie weights in the mobility network between
#' the locations of an origin and the destinations.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values. Level k is the
#' number of ties with weight k - 1, i.e., level 1 counts the empty dyads.
#' If NULL, all levels from 1 to the maximum tie weight + 1 are used.
#' @param origin The name of the origin (see [createOrigin()]) whose ties are
#' counted; passed from [gofMobilityNetwork()].
#' @param ... Additional parameters.
#'
#' @return The tie weight distribution of simulated networks, named by the
#' tie weight each count refers to.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'
#' @keywords internal
getTieWeights <- function(cache, dep.var, lvls, origin, ...) {
  if (missing(origin)) {
    stop("getTieWeights needs the origin of the ties, e.g. origin = \"origin\"")
  }
  weights <- as.vector(cache[[dep.var]][[origin]]$valuedNetwork)
  if (is.null(lvls)) {
    lvls <- 1:(max(weights) + 1)
  }
  counts <- tabulate(weights + 1, nbins = max(c(lvls, weights + 1)))
  setNames(counts[lvls], lvls - 1)
}
