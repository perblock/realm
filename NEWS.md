# realm 0.1.0

- First version, based on MoNAn 2.1.1.9000
  (<https://github.com/stocnet/MoNAn>). The package, its S3 classes
  (`*.realm`), the alias functions (`realmEstimate()`, `realmDataCreate()`,
  ...) and the options (`realm.useCpp`, `realm.mpleMemoryLimit`) are renamed;
  the alias `monan07()` is removed.
- Effects based on reciprocity (`reciprocity_*`) and transitivity
  (`transitivity_*`, `triad120*`) are removed, together with the net flow,
  min and two-path networks in the cache that only they used. The example
  data are re-estimated without `reciprocity_AC`.
- The dependent variable contains only the destinations of the resources
  (`createDependent()`/`realmDependent()`, class `dependent.realm`, which
  replace `createEdgelist()` and `edgelist.realm`). The origins are separate
  objects (`createOrigin()`/`realmOrigin()`, class `origin.realm`); a state
  can have several origins, also of another nodeset than the destinations
  (e.g. regions), whose networks are then rectangular.
- Effects that depend on the origin of the resources have the argument
  `origin`, which must be given in `addEffect()` (e.g.
  `addEffect(loops, origin = "origin")`) and becomes part of the effect name.
  `in_weights_AC`, `in_weights_GW`, `alter_covariate`,
  `resource_covar_to_node_covar` and `loops_additional_origin` depend only on
  the destinations and take no origin; their target statistics are sums over
  the destinations.
- The cache holds one network for each origin (`cache[[dep.var]][[origin]]`)
  and the inflow of the destinations (`cache[[dep.var]]$destinations`).
  `getTieWeights()` needs the origin, e.g.
  `gofMobilityNetwork(ans, getTieWeights, origin = "origin")`.
- The covariates of effects are checked against the nodesets they need
  (destinations, the origin's locations, or both for dyadic covariates), and
  effects about loops need an origin with the nodeset of the destinations.
- Loops are always allowed: the argument `allowLoops` of `createAlgorithm()`
  and `simulateMobilityNetworks()` and the argument `loop.contribution` of
  `alter_covariate` and `resource_covar_to_node_covar` are removed (loops
  always count, also in `resource_covar_to_node_covar`, which excluded them
  by default). The example data are re-estimated.
- `loops_additional_origin` is deprecated: create the additional origin as
  an origin and use `addEffect(loops, origin = <its name>)`, which gives the
  same statistics.
- Effects that excluded loops (`concentration_AC`, `concentration_GW` and
  their `_dyad_covar` and `_resource_covar_bin` versions, `dyadic_covariate`,
  `dyadic_covariate_resource_attribute`, `same_covariate`, `sim_covariate`,
  `present_relations`) now include them by default; the new argument
  `exclude.loops = TRUE` excludes them as before. The example data are
  re-estimated.
- `same_covariate` and `sim_covariate` have the argument `origin.attribute`
  (in `addEffect()`): a covariate of the origin's locations that is compared
  with the covariate of the destinations (`node.attribute`). It is needed for
  an origin of another nodeset than the destinations, e.g.
  `addEffect(same_covariate, origin = "originRegion", node.attribute =
  "orgRegion", origin.attribute = "regionId")`. The similarity uses the range
  of both covariates. Without it, the effects work as before.
- Effects check that their covariates are covariates of the state of the
  kind and size they read: a covariate (`createNodeVariable()`) of the
  destinations or the origin's locations, or a dyadic covariate
  (`createNetwork()`) between them, as the effect requires. Before, a node
  variable given as a dyadic covariate passed the checks and the C++ sampler
  read outside its data. The C++ effects also check the size of their
  covariates, and the elements of a process state must have unique names.
- `createOrigin()`/`realmOrigin()` have the argument `addDummies` (default
  `TRUE`), which adds the element `dummy` with a dummy for each origin value
  larger than 1 in the data.
- New effect `origin_to_node_covar`: do people from a category of an origin
  move to locations with high values of a covariate of the destinations?
  `addEffect(origin_to_node_covar, origin = "or_educat", node.attribute =
  "occ_covar")` adds one effect per dummy of the origin (value 1 is the
  reference), named e.g. `origin_to_node_covar or_educat 2 to occ_covar`, and
  reports the number of effects.
- New function `createNodeVariables()`, with the alias `monadicCovars()`, as
  the one function to create covariates: from a data.frame, each column
  becomes a covariate named by the column; from a vector, it gives one
  covariate as `createNodeVariable()`. Effects refer to the columns by name,
  e.g. `node.attribute = "size"`.
- `addDummyEffects()` only adds the dummies created by `realmDataCreate()`
  (`dummy2`, `dummy3`, ...), not every covariate with "dummy" in its name.
- Printing an effects object shows the origin of each effect in a column of
  its own, and the dummy of `origin_to_node_covar` under parameter (`d 2`).
- Printing an effects object or a state leaves out the fixed effects and the
  location dummies and gives their number in one line;
  `print(x, showDummies = TRUE)` shows them.
- `estimateMobilityNetwork()`/`realmEstimate()` have the argument
  `startingValues`: `"partial MPLE"` (the new default) sets endogenous effects
  to zero and takes the others from the MPLE of the model without them,
  `"full MPLE"` is the previous default and `"zero"` starts from zero.
  `initialParameters` are used as given: an all-zero vector no longer means
  MPLE. With the default, a model of only exogenous effects is a multinomial
  logit, so its MPLE is returned directly with standard errors from the
  Hessian, without phases 1 to 3; with `returnDeps = TRUE`, phase 3 simulates
  networks at these estimates.
