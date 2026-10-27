# Adding an effect to realm

Note: This guide was written by Claude Code, which produced the code for the majority
of the C++ implementation of the MCMC. The C++ implementation is checked in all
steps and effects that it exactly reproduces the human-coded R-based implementation.

This guide is for developers of realm. It explains how to add a new effect so
that it works with both samplers:

- the **R sampler** (`R/1_hiddenFunctions.R`), which calls the effect's R function;
- the **C++ sampler** (`src/sampler.cpp`), which calls the effect's C++ class and
  is used whenever *all* effects of a model have a C++ version.

Every effect therefore exists twice: as an R function (the reference
implementation) and as a C++ class. The tests require both versions to give
**exactly** the same values, so that both samplers give identical results for
the same seed. If an effect only exists in R, models with it still work, but use
the (20–50 times slower) R sampler.

Adding an effect has five steps:

1. Write the R function.
2. Export and document it.
3. Add it to the change/target tests.
4. Write the C++ class.
5. Add it to the C++ tests.

The example used throughout is `alter_covariate_squared`: do individuals move to
locations with high (squared) values of a location covariate? Its target
statistic is the sum over all destinations of the number of individuals in the
destination times its squared covariate. It depends only on the destinations,
not on where the individuals come from; most effects depend on an origin (see
below).

---

## 1. The R function

Effects live in `R/2*_effectFunctions*.R`, grouped by type (exogenous, loops,
concentration, covariate-based). Every effect has this
signature:

```r
alter_covariate_squared <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (getTargetContribution) {
      return(cache[[dep.var]]$destinations$inflow[j] * state[[attribute.index]]$data[j]^2)
    }
    return(update * state[[attribute.index]]$data[j]^2)
  }
```

The dependent variable contains only the destinations of the resources. Their
origins are separate objects (`createOrigin()`), and a state can have several
origins, also of another nodeset than the destinations (e.g. the region an
individual comes from). **An effect that depends on the origin of the resources
has the argument `origin`** (without default, after `dep.var`), the name of the
origin the user chooses in `addEffect(..., origin = "origin")`; effects without
this argument depend only on the destinations.

The function is called in two ways:

- **Target contribution** (`getTargetContribution = TRUE`): for an effect with
  an origin, the contribution of the tie from location `i` of the origin to
  destination `j`; the target statistic is the sum over all `i`, `j`
  (`getNetworkStatistics()`). For an effect without an origin, the contribution
  of destination `j` (`i` is `NA`); the target statistic is the sum over all `j`.
- **Change statistic** (`getTargetContribution = FALSE`): how the target statistic
  changes when resource `edge` is added to (`update = 1`) or removed from
  (`update = -1`) destination `j`, coming from `i` (its location in the
  effect's origin). A move of a resource from `j` to `l` is calculated as
  removing it from `i -> j` and then adding it to `i -> l`.

**The change statistic must equal the difference in target statistics.** If it
does not, estimates are biased without any error; this is what the tests in
step 3 check.

What the function can use:

| object | content |
|---|---|
| `cache[[dep.var]][[origin]]$valuedNetwork` | number of resources on each tie (origin locations x destinations) |
| `cache[[dep.var]][[origin]]$resourceNetworks[[name]]` | sum of a resource covariate on each tie |
| `cache[[dep.var]][[origin]]$sameNodeSet` | whether the origin has the nodeset of the destinations, i.e. whether there are loops |
| `cache[[dep.var]]$destinations$inflow` | number of resources in each destination |
| `cache[[dep.var]]$destinations$resourceInflows[[name]]` | sum of a resource covariate in each destination |
| `state[[attribute.index]]$data` | a location covariate (vector) or dyadic covariate (matrix, origin locations x destinations) |
| `state[[resource.attribute.index]]$data` | a resource covariate (one value per resource) |
| `state[[dep.var]]$data` | the destination of each resource |
| `state[[origin]]$data` | the location of each resource in the origin |

Rules:

- **The cache already reflects the removal when a resource is added.** In a move,
  the "add" step is calculated after the resource has been removed from the
  cache, and its destination is `NA` at that point. Use `i` for its origin.
- **Loops exist only for an origin with the nodeset of the destinations.**
  Loops count by default. An effect for which excluding loops makes sense has
  the argument `exclude.loops = FALSE` and checks
  `exclude.loops && i == j && cache[[dep.var]][[origin]]$sameNodeSet`
  (`exclude.loops = TRUE` is only allowed for such origins). An effect that is
  about loops (a resource staying at its origin) goes into `sameNodeSetEffects`
  in `R/2_effectFunctions.R`, which allows it only for such origins.
- **Name arguments `attribute.index` (location or dyadic covariate) and
  `resource.attribute.index` (resource covariate).** Users specify them as
  `node.attribute` and `edge.attribute` in `addEffect()`, which renames them.
  Declare in `attributeRoles` (`R/2_effectFunctions.R`) whether
  `attribute.index` is a covariate of the destinations (used with `j`), of the
  origin's locations (used with `i`), or a dyadic covariate (`[i, j]`); the
  covariate's nodesets are checked against this. Further parameters (such as
  `alpha`) get defaults.
- **Don't check the input data in the function.** Effects are called millions of
  times. Checks that a covariate is symmetric, binary, etc. go into
  `effectInputChecks` in `R/2_effectFunctions.R`, which runs once when the effect
  is added and when estimation or simulation starts. Checks of parameters that
  cost O(1), such as `if (alpha < 1) stop(...)`, can stay in the function.
- **Use `cache[[dep.var]]` and `state[[dep.var]]`**, never a fixed name such as
  `state$transfers`.

## 2. Export and document

Add a roxygen block above the function, as for the other effects (description,
`@param` for each argument, `@return`, `@keywords internal`):

```r
#' alter_covariate_squared
#'
#' Do individuals move to locations with high values of a covariate, with the
#' covariate's effect increasing quadratically?
#'
#' @param dep.var
#' @param attribute.index
#' ...
#' @return Returns the change statistic or target statistic of the effect for
#' internal use by the estimation algorithm.
#' @keywords internal
```

`NAMESPACE` is written by hand (not by roxygen): add the effect's name to the
`export(...)` list, in the section for its group (e.g. `# effect functions:
exogenous`). Then run `devtools::document()` to create the help page.

## 3. Change/target tests

`tests/testthat/test-effects-change-target.R` checks for every exported effect
that the change statistic equals the difference in target statistics, for
random moves and moves into and out of loops. Add an entry to `effect_specs`
with the effect's arguments, using the covariates of the test data listed at the
top of the file:

```r
  alter_covariate_squared = list(attribute.index = "size"),
```

For further parameter values, add entries with `.effect`, e.g.
`concentration_AC_alpha4 = list(.effect = "concentration_AC", alpha = 4)`.
A test fails if an exported effect has no entry.

Run the tests: `devtools::test(filter = "effects-change-target")`. Set the
environment variable `REALM_FULL_EFFECT_TESTS=true` to move every resource
instead of a sample.

## 4. The C++ class

### The class

C++ effects live in `src/effects.cpp`. An effect implements the interface
`Effect` from `src/realm.h`:

```cpp
// alter_covariate_squared (R/2a_effectFunctions_exogenous.R)
class AlterCovariateSquared : public Effect {
  std::vector<double> attribute;
 public:
  explicit AlterCovariateSquared(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    return update * rpow(attribute[j], 2);
  }
  double target(const Chain& ch, int i, int j) const override {
    return ch.inflow[j] * rpow(attribute[j], 2);
  }
};
```

Put the class inside the unnamed `namespace { ... }` with the other effects.
Indices are **0-based** (`i`, `j`, `edge`; location ids in covariates stay as in
R, e.g. `ch.resourceValues[r][e] == j + 1`). The `Chain` and, for effects with
an origin, `net(ch)` (the cache of the effect's origin) give access to:

| C++ | R |
|---|---|
| `net(ch).v(i, j)` | `cache[[dep.var]][[origin]]$valuedNetwork[i, j]` |
| `net(ch).resourceNetworks[r][net(ch).idx(i, j)]` | `cache[[dep.var]][[origin]]$resourceNetworks[[name]][i, j]` |
| `net(ch).sameNodeSet` | `cache[[dep.var]][[origin]]$sameNodeSet` |
| `net(ch).colSum(m, j)`, `net(ch).rowSum(m, i)` | `sum(m[, j])`, `sum(m[i, ])` |
| `net(ch).origin[e]` | `state[[origin]]$data[e]` |
| `ch.inflow[j]`, `ch.resourceInflows[r][j]` | `cache[[dep.var]]$destinations$inflow[j]`, `...$resourceInflows[[name]][j]` |
| `ch.resourceValues[r][edge]` | `state[[resource.attribute.index]]$data[edge]` |
| `ch.receiver[e]` | `state[[dep.var]]$data[e]` (`-1` while `e` is moved) |
| `ch.n`, `net(ch).nRow` | number of destinations, of the origin's locations |

Matrices are column-major vectors of length `nRow * ch.n`, as in R. The effect's
origin is set from its specification (`o` in `Effect`, -1 for effects without
an origin), and the samplers pass its `i`.

### Creating the effect from its specification

Create the effect in `makeEffect()` at the end of `src/effects.cpp`:

```cpp
  if (name == "alter_covariate_squared") {
    return E(new AlterCovariateSquared(attributeOf(spec)));
  }
```

The specification `spec` comes from `cppEffectSpec()` in `R/1b_cppSampler.R`,
which translates the effect's arguments:

| R argument | in `spec` | read in C++ with |
|---|---|---|
| `alpha` | `spec["alpha"]` | `alphaOf(spec, "name")` (checks `alpha >= 1`) or `Rcpp::as<double>(spec["alpha"])` |
| `exclude.loops` | `spec["excludeLoops"]` | `excludeLoops` of `Effect`, set by `makeEffect()` |
| `attribute.index` | the covariate's data as a numeric vector | `attributeOf(spec)` |
| `resource.attribute.index` | the covariate's name | `chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]))` |

An effect with a new kind of argument needs an addition to `cppEffectSpec()`.
Finally, add the effect's name to `cppEffects` in `R/1b_cppSampler.R`; only
effects listed there use the C++ sampler.

Optionally, an effect can override `targetTotal()`, which calculates the whole
target statistic at once. The default sums `target()` over all cells (or over
the destinations for effects without an origin); override
it only if this is much slower than a direct calculation (see
`LoopsAdditionalOrigin` and `SimilarityCont`), and only if the result is exactly
the same.

### Getting exactly the same values as R

The tests compare C++ and R values with `identical()`, and chains only stay
identical if every acceptance decision is the same. Floating-point results must
therefore match R to the last bit:

- **Keep R's order of operations.** `a * b * c` is `(a * b) * c`;
  `update * g(y) * a` in R must be `update * g(y) * a` in C++, not
  `update * (g(y) * a)`.
- **Use `rpow(x, y)` for R's `x^y`.** R calculates `x^2` as `x * x`, which can
  differ from `std::pow(x, 2)`.
- **Accumulate sums in `long double`**, as R's `sum()`, `colSums()` and
  `cumsum()` do. `colSum()` and `rowSum()` of `OriginNet` in `src/realm.h` do
  this.
- **Reproduce R's sequences.** `for (k in 0:y)` counts *down* for negative `y`
  (use `rSeqFromZero()`), and `1:(y - 1)` goes up to `floor(y - 1)`.
- **Write loops in R's order.** A sum that R calculates as `contr <- contr + x`
  in a loop is a `double`, not a `long double`.
- **No fused multiply-add.** The pragma at the top of `src/realm.h` switches it
  off (the compiler would otherwise compute `a + b * c` with one rounding instead
  of two). Every C++ file must include `realm.h`.

Values that are integers in practice (counts of resources) are exact in any
order, so these rules matter most for effects with weights such as AC/GW
functions or continuous covariates.

## 5. C++ tests

Add the effect to a test file `tests/testthat/test-cpp-effects-*.R` for its
group, with the helpers from `tests/testthat/helper-cpp.R`:

```r
newEffects <- createEffects(myState) |>
  addEffect(alter_covariate_squared, node.attribute = "size")

test_that("C++ alter_covariate_squared gives the same values as R", {
  expect_cpp_matches_r(newEffects)
})

test_that("C++ and R chains are identical with alter_covariate_squared", {
  compare_chains(newEffects, 0.05, n = 5000)
  compare_chains(newEffects, 0.05, n = 500, multinomial = TRUE)
})
```

- `expect_cpp_matches_r()` compares the change statistic and target
  contribution of every effect with its R version for random moves, including
  moves into and out of loops, and the whole target statistics. On R
  built without `long double` (CRAN's "noLD" check), R sums in `double` while
  the C++ code still uses `long double`, so the comparison allows a relative
  difference of 1e-10 there (`expect_same()`); everywhere else it is exact.
- `compare_chains()` runs the R and the C++ sampler with the same seed and
  checks that the destinations, the cache and the random number state are identical
  afterwards.
- Use several parameter values (e.g. two values of `alpha`) and, for covariates,
  both a binary and a continuous one. `cppTestState` is `myState` with a
  continuous dyadic covariate (`size_similarity`) and a `same` matrix for
  `region`.

The test "all exported effects have a C++ version" in
`tests/testthat/test-cpp-effects-other.R` fails until the effect is in
`cppEffects`.

## Checklist

- [ ] R function in `R/2*_effectFunctions*.R`; data checks in `effectInputChecks`
- [ ] roxygen documentation; name in `NAMESPACE`; `devtools::document()`
- [ ] entry in `effect_specs` (`tests/testthat/test-effects-change-target.R`)
- [ ] C++ class and `makeEffect()` entry in `src/effects.cpp`; name in `cppEffects`
- [ ] C++ tests with `expect_cpp_matches_r()` and `compare_chains()`
- [ ] `devtools::test()` and `R CMD check` pass
- [ ] entry in `NEWS.md`

When an existing effect is changed, change **both** versions; the C++ tests fail
if they differ.

## Effects written only in R

Users can also write their own effects in R, with the signature from step 1
(with `origin` if they depend on the origin), and add them with `addEffect()`. The function must be defined in the global
environment, because `addEffect()` finds effects (and the state) by name there.
Models with such effects use the R sampler automatically. `options(realm.useCpp = FALSE)` makes all models use the R
sampler, e.g. to compare both.
