########## effectFunctions


########## input checks of effects
# Checks of the data an effect refers to, run once when an effect is added
# (createEffectsObject) and when estimation or simulation starts, instead of
# in every call of the effect function: effects are called millions of times,
# and these checks take time proportional to the size of the data.

# checks per effect: functions of the state and the effect's arguments
effectInputChecks <- local({
  binaryResource <- function(effect) {
    function(state, args) {
      if (!all(state[[args$resource.attribute.index]]$data %in% c(0, 1))) {
        stop(paste("effect", effect, "only defined for binary covariates coded 0 1"))
      }
    }
  }
  list(
    crowding_out_prop_covar_bin = list(binaryResource("crowding_out_prop_covar_bin")),
    staying_by_prop_bin_inflow = list(binaryResource("staying_by_prop_bin_inflow"))
  )
})

# deprecated effects, with what to use instead; the warning is given when the
# effect is added to a model (createEffectsObject)
deprecatedEffects <- c(
  loops_additional_origin = paste(
    "loops_additional_origin is deprecated: create the additional origin as an",
    "origin (createOrigin()) and use addEffect(loops, origin = <its name>)."
  )
)

# the nodeset of the covariate given as attribute.index: a covariate of the
# destinations (used for j), of the origin's locations (used for i), or a
# dyadic covariate between both (used for [i, j])
attributeRoles <- c(
  alter_covariate = "destination",
  resource_covar_to_node_covar = "destination",
  same_covariate = "destination",
  sim_covariate = "destination",
  dyadic_covariate = "dyad",
  dyadic_covariate_resource_attribute = "dyad",
  concentration_AC_dyad_covar = "dyad",
  concentration_GW_dyad_covar = "dyad",
  loops_node_covar = "origin",
  loops_resource_covar_node_covar = "origin",
  concentration_prop_orig_cov = "origin"
)

# effects that need an origin with the same nodeset as the destinations:
# effects defined by loops (a resource staying at its origin)
sameNodeSetEffects <- c(
  "loops", "loops_AC", "loops_GW", "loops_node_covar", "loops_resource_covar",
  "loops_resource_covar_node_covar", "loops_x_loops_additional_origin",
  "in_ties_loops", "crowding_out_prop_covar_bin", "staying_by_prop_bin_inflow"
)

# effects that compare a covariate of the origin's locations
# (origin.attribute.index) with one of the destinations (attribute.index); for
# an origin of another nodeset than the destinations, the origin's covariate
# must be given
originAttributeEffects <- c("same_covariate", "sim_covariate")

# checkEffectNodeSets
# checks that the origin of an effect is in the state, and that the effect's
# origin and covariates have the nodesets the effect needs
checkEffectNodeSets <- function(effectName, args, state) {
  destinations <- destinationNodeSet(state)
  origin <- args$origin
  originSet <- NULL
  if (is.character(origin)) {
    if (!(origin %in% state$origins)) {
      stop(paste0("Effect ", effectName, ": '", origin, "' is not an origin of the state."))
    }
    originSet <- state[[origin]]$nodeSet[1]
  }
  if (effectName %in% sameNodeSetEffects && !identical(originSet, destinations)) {
    stop(paste0("Effect ", effectName, " needs an origin with the nodeset of the destinations (",
                destinations, ")."))
  }
  if (effectName %in% originAttributeEffects && !identical(originSet, destinations) &&
      !is.character(args$origin.attribute.index)) {
    stop(paste0("Effect ", effectName, ": the origin ", origin, " has another nodeset (",
                originSet, ") than the destinations; specify the covariate of its ",
                "locations with origin.attribute."))
  }
  if (is.character(args$origin.attribute.index)) {
    actual <- state[[args$origin.attribute.index]]$nodeSet
    if (!identical(as.character(actual), originSet)) {
      stop(paste0("Effect ", effectName, ": the covariate '", args$origin.attribute.index,
                  "' must be defined for ", originSet, ", not for ",
                  paste(actual, collapse = " x "), "."))
    }
  }
  if (isTRUE(args$exclude.loops) && !identical(originSet, destinations)) {
    stop(paste0("Effect ", effectName, ": exclude.loops needs an origin with the nodeset of ",
                "the destinations (", destinations, "); other origins have no loops."))
  }
  role <- attributeRoles[effectName]
  if (!is.na(role) && is.character(args$attribute.index)) {
    expected <- switch(role,
                       destination = destinations,
                       origin = originSet,
                       dyad = c(originSet, destinations))
    actual <- state[[args$attribute.index]]$nodeSet
    if (role == "dyad" && length(actual) == 1) {
      actual <- c(actual, actual)
    }
    if (!identical(as.character(actual), expected)) {
      stop(paste0("Effect ", effectName, ": the covariate '", args$attribute.index,
                  "' must be defined for ", paste(expected, collapse = " x "),
                  ", not for ", paste(actual, collapse = " x "), "."))
    }
  }
  if (is.character(args$resource.attribute.index)) {
    actual <- state[[args$resource.attribute.index]]$nodeSet
    if (!identical(as.character(actual), resourceNodeSet(state))) {
      stop(paste0("Effect ", effectName, ": the covariate '", args$resource.attribute.index,
                  "' must be defined for ", resourceNodeSet(state), "."))
    }
  }
  invisible(NULL)
}

# checkEffectInputs
checkEffectInputs <- function(effects, state) {
  for (k in seq_along(effects$effectFormulas)) {
    effectName <- strsplit(effects$name[k], " ")[[1]][1]
    args <- as.list(formals(effects$effectFormulas[[k]]))
    checkEffectNodeSets(effectName, args, state)
    for (check in effectInputChecks[[effectName]]) {
      check(state, args)
    }
  }
  invisible(NULL)
}



#' in_ties_loops
#' 
#' Are individuals that are in locations with a large inflow more likely to stay 
#' in their current location?
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
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
in_ties_loops <-
  function(dep.var = 1,
           origin,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (getTargetContribution) {
      return((cache[[dep.var]][[origin]]$valuedNetwork[i, j] * (i == j)) *
        (sum(cache[[dep.var]][[origin]]$valuedNetwork[, j]) - cache[[dep.var]][[origin]]$valuedNetwork[j, j]))
    }

    if (i == j) {
      return(update * (sum(cache[[dep.var]][[origin]]$valuedNetwork[, j]) - cache[[dep.var]][[origin]]$valuedNetwork[j, j]))
    }

    if (i != j) {
      return(update * cache[[dep.var]][[origin]]$valuedNetwork[j, j])
    }
  }




#' in_weights_AC
#' 
#' Is there a preferential attachment in the mobility network, i.e., do individuals 
#' move particularly to popular destinations?
#' The geometrically weighted version avoids degeneracy.
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
in_weights_AC <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 2) {
    if (alpha <= 0) {
      stop("Alpha parameter in in_weights_AC weights function must be positive")
    }
    
    # the effect depends only on the destination: i is not used, and the target
    # contribution is the one of destination j
    in_weight <- cache[[dep.var]]$destinations$inflow[j]
    
    if (getTargetContribution) {
      g_cum <- function(y, a){
        contr <- 0
        if(y>1){
          for(k in 1:(y-1)){
            contr <- contr + (1 - (1-1/a)^(k))
          }
        }
        return(contr)
      }
      
      return(g_cum(y = in_weight, a = alpha))
    }

    g_mar <- function(y, a){
      contr <- 0
      if(y>0) {
        contr <-  (1 - (1-1/a)^(y)) 
      } else {
        contr <- 0
      }
      return(contr)
    }
    
    if(update < 0){
      return(update * g_mar(y = (in_weight + update), a = alpha))
    }
    if(update > 0){
      return(update * g_mar(y = in_weight, a = alpha))
    }
  }


#' in_weights_GW
#' 
#' Is there a preferential attachment in the mobility network, i.e., do individuals 
#' move particularly to popular destinations?
#' The geometrically weighted version avoids degeneracy.
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
in_weights_GW <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 2) {
    if (alpha <= 0) {
      stop("Alpha parameter in in_weights_GW weights function must be positive")
    }
    
    # the effect depends only on the destination: i is not used, and the target
    # contribution is the one of destination j
    in_weight <- cache[[dep.var]]$destinations$inflow[j]
    
    if (getTargetContribution) {
      g_cum <- function(y, a){
        contr <- 0
        for(k in 0:y){
          contr <- contr + (y-k) * exp(-log(a)*k)
        }
        contr - y
      }
      
      return(g_cum(y = in_weight, a = alpha))
    }
    
    g_mar <- function(y, a){
      contr <- 0
      for(k in 0:y){
        contr <- contr + exp(-log(a)*k)
      }
      contr - 1
    }
    
    if(update < 0){
      return(update * g_mar(y = (in_weight + update), a = alpha))
    }
    if(update > 0){
      return(update * g_mar(y = in_weight, a = alpha))
    }
  }



#' present_relations
#' 
#' Do individuals move along many or few paths out of their origin? This models 
#' whether individuals have a tendency against being the only one
#' moving to a particular destination from their origin.
#'
#' @param dep.var 
#' @param origin The name of the origin (see [createOrigin()]) of the
#' resources.
#' @param exclude.loops Should loops (resources staying at their origin) be
#' excluded from the statistic? Defaults to FALSE, i.e. loops contribute.
#' Loops exist only for an origin with the nodeset of the destinations.
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
present_relations <-
  function(dep.var = 1,
           origin,
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
      return((cache[[dep.var]][[origin]]$valuedNetwork[i, j] > 0) * 1)
    }

    if (cache[[dep.var]][[origin]]$valuedNetwork[i, j] == 0 &&
      update > 0) {
      return(1)
    }
    if (cache[[dep.var]][[origin]]$valuedNetwork[i, j] == -update &&
      update < 0) {
      return(-1)
    }

    return(0)
  }



## Deprecated ##
## in_weights_exponent: the target statistic excludes loops and divides by the number of
## nodes, the change statistic does neither, so change and target statistics do not match.
## If de-commenting, fix this and put back the name in the namespace ##

#' #' in_weights_exponent
#' #' 
#' #' Is there a preferential attachment in the mobility network, i.e., do individuals 
#' #' move particularly to popular destinations?
#' #'
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #' @param exponent 
#' #'
#' #' 
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm.
#' #' @keywords internal
#' in_weights_exponent <-
#'   function(dep.var = 1,
#'            state,
#'            cache,
#'            i,
#'            j,
#'            edge,
#'            update,
#'            getTargetContribution = FALSE,
#'            exponent = .5) {
#'     if (i == j) {
#'       return(0)
#'     }
#'
#'     # a target contribution is calculated even for unconnected i-j pairs
#'     if (getTargetContribution) {
#'       return(sum(cache[[dep.var]]$valuedNetwork[, j])^exponent / length(cache[[dep.var]]$valuedNetwork[, j]))
#'     }
#'
#'     v <-
#'       (sum(cache[[dep.var]]$valuedNetwork[, j]) + update)^exponent -
#'       (sum(cache[[dep.var]]$valuedNetwork[, j]))^exponent
#'
#'     return(v)
#'   }
