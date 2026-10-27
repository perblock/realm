########## effectFunctions: exogenous effects

#' alter_covariate
#' 
#' Are locations higher on some attribute v more popular targets of mobility? 
#' E.g., do workers have a tendency to move to larger organisations?
#'
#' @param dep.var 
#' @param attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
alter_covariate <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    # the effect depends only on the destination: i is not used, and the target
    # contribution is the one of destination j
    if (getTargetContribution) {
      return(state[[attribute.index]]$data[j] * cache[[dep.var]]$destinations$inflow[j])
    }
    
    return(state[[attribute.index]]$data[j] * update)
  }


#' dyadic_covariate
#' 
#' Is mobility between locations predicted by the dyadic covariate u? E.g., is mobility 
#' likely between organisations that are in the same region? Note that in many cases 
#' dyadic covariates can convey the same information as
#' the ‘same covariate’ or the ‘covariate similarity’ effects.
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
#' @param exclude.loops Should loops (resources staying at their origin) be
#' excluded from the statistic? Defaults to FALSE, i.e. loops contribute.
#' Loops exist only for an origin with the nodeset of the destinations.
#' @param attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
dyadic_covariate <-
  function(dep.var = 1,
           origin,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           exclude.loops = FALSE,
           getTargetContribution = FALSE) {
    if (exclude.loops && i == j && cache[[dep.var]][[origin]]$sameNodeSet) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]][[origin]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i, j])
    }
    
    return(update * state[[attribute.index]]$data[i, j])
  }


#' dyadic_covariate_resource_attribute
#' 
#' Is mobility between locations predicted by the dyadic covariate u weighted by 
#' the individual covariate w? E.g., is mobility of women more likely between 
#' organisations that are in the same region? Note that this effect can be used 
#' to also model the interaction between the ‘same covariate’/‘covariate similarity’ 
#' effect and individual attributes, since sameness and similarity between locations 
#' can be translated into dyadic covariates.
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
#' @param exclude.loops Should loops (resources staying at their origin) be
#' excluded from the statistic? Defaults to FALSE, i.e. loops contribute.
#' Loops exist only for an origin with the nodeset of the destinations.
#' @param attribute.index 
#' @param resource.attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
dyadic_covariate_resource_attribute <-
  function(dep.var = 1,
           origin,
           attribute.index,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           exclude.loops = FALSE,
           getTargetContribution = FALSE) {
    if (exclude.loops && i == j && cache[[dep.var]][[origin]]$sameNodeSet) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]][[origin]]$resourceNetworks[[resource.attribute.index]][i, j] * state[[attribute.index]]$data[i, j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge] * state[[attribute.index]]$data[i, j])
  }

#' resource_covar_to_node_covar
#' 
#' Do individuals with some individual attribute w tend to move to locations with 
#' some location characteristic v? E.g., do women move to larger organisations than men?
#'
#' @param dep.var 
#' @param resource.attribute.index 
#' @param attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
resource_covar_to_node_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    # the effect depends only on the destination: i is not used, and the target
    # contribution is the one of destination j
    if (getTargetContribution) {
      return(state[[attribute.index]]$data[j] *
               cache[[dep.var]]$destinations$resourceInflows[[resource.attribute.index]][j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge] * state[[attribute.index]]$data[j])
  }


#' same_covariate
#' 
#' Is mobility more likely between locations that are identical on some attribute v? 
#' E.g., is mobility more likely between organisations that are located in the same region?
#' 
#' The attribute of the destination (node.attribute) is compared with the
#' attribute of the origin location. Without origin.attribute, this is the same
#' covariate (which needs the 'same' matrix, see [createNodeVariable()]); with
#' origin.attribute, a covariate of the origin's locations, which is needed for
#' an origin of another nodeset than the destinations (e.g. the region of
#' origin compared with the region of the destination organisation).
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
#' @param exclude.loops Should loops (resources staying at their origin) be
#' excluded from the statistic? Defaults to FALSE, i.e. loops contribute.
#' Loops exist only for an origin with the nodeset of the destinations.
#' @param attribute.index The covariate of the destinations (node.attribute in
#' [addEffect()]).
#' @param origin.attribute.index The covariate of the origin's locations
#' (origin.attribute in [addEffect()]); if NULL, attribute.index is used.
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
same_covariate <-
  function(dep.var = 1,
           origin,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           origin.attribute.index = NULL,
           exclude.loops = FALSE,
           getTargetContribution = FALSE) {
    if (exclude.loops && i == j) {
      return(0)
    }
    
    if (is.null(origin.attribute.index)) {
      if (is.null(state[[attribute.index]]$same)) {
        stop("Effect same_covariate expects covariates with a 'same' attribute.")
      }
      same <- state[[attribute.index]]$same[i, j]
    } else {
      same <- (state[[origin.attribute.index]]$data[i] == state[[attribute.index]]$data[j]) * 1
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]][[origin]]$valuedNetwork[i, j] * same)
    }
    
    return(update * same)
  }


#' sim_covariate
#' 
#' Is mobility more likely between locations that are similar on some attribute v? 
#' E.g., is mobility more likely between organisations that have a similar size?
#' 
#' The similarity of the attribute a of the origin location and the attribute b
#' of the destination is 1 - |a - b| / range, with the range of the values of
#' both. Without origin.attribute, both are the same covariate (which needs the
#' 'sim' matrix, see [createNodeVariable()]); with origin.attribute, a
#' covariate of the origin's locations, which is needed for an origin of
#' another nodeset than the destinations.
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
#' @param exclude.loops Should loops (resources staying at their origin) be
#' excluded from the statistic? Defaults to FALSE, i.e. loops contribute.
#' Loops exist only for an origin with the nodeset of the destinations.
#' @param attribute.index The covariate of the destinations (node.attribute in
#' [addEffect()]).
#' @param origin.attribute.index The covariate of the origin's locations
#' (origin.attribute in [addEffect()]); if NULL, attribute.index is used.
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
sim_covariate <-
  function(dep.var = 1,
           origin,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           origin.attribute.index = NULL,
           exclude.loops = FALSE,
           getTargetContribution = FALSE) {
    if (exclude.loops && i == j) {
      return(0)
    }
    
    if (is.null(origin.attribute.index)) {
      if (is.null(state[[attribute.index]]$sim)) {
        stop("Effect sim_covariate expects covariates with a 'sim' attribute.")
      }
      sim <- state[[attribute.index]]$sim[i, j]
    } else {
      # as the 'sim' matrix of createNodeVariable(), with the range of both
      # covariates
      a <- state[[origin.attribute.index]]$data
      b <- state[[attribute.index]]$data
      valueRange <- range(c(a, b))
      sim <- 1 - abs(a[i] - b[j]) / (valueRange[2] - valueRange[1])
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]][[origin]]$valuedNetwork[i, j] * sim)
    }
    
    return(update * sim)
  }
