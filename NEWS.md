# realm 0.0.0.9000

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
