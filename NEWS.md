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
