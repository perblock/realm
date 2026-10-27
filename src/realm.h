// Data structures of the C++ MCMC sampler: one chain (destinations, origins
// and cache) and the interface of effects. Indices are 0-based; matrices are
// column-major vectors, as in R.
#ifndef REALM_H
#define REALM_H

// No fused multiply-add: R rounds a + b * c in two steps, and the C++ effects
// must give exactly the same values as the R effects. (A pragma rather than
// -ffp-contract=off in Makevars, which R CMD check reports as non-portable.)
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include <Rcpp.h>
#include <Rmath.h>
#include <memory>
#include <string>
#include <vector>

// R's x^y (arithmetic.c): x * x for y == 2, R_pow() otherwise. R_pow() is
// R's own power function, which does not use pow() on all platforms (on
// 64-bit Windows, it uses the long double powl()).
inline double rpow(double x, double y) {
  if (y == 2.0) return x * x;
  return R_pow(x, y);
}

// the cache of one origin (as cache[[dep.var]][[origin]] in R): the network
// between the origin's locations (rows) and the destinations (columns)
struct OriginNet {
  int nRow, nCol;                         // origin locations, destinations
  bool sameNodeSet;                       // origins are destinations: loops exist
  std::vector<int> origin;                // per resource
  std::vector<double> valued;
  std::vector<std::vector<double>> resourceNetworks;

  inline int idx(int i, int j) const { return i + j * nRow; }
  inline double v(int i, int j) const { return valued[idx(i, j)]; }
  // sum(m[, j]) and sum(m[i, ]) as in R (long double accumulation)
  double colSum(const std::vector<double>& m, int j) const {
    long double total = 0;
    for (int i = 0; i < nRow; ++i) total += m[idx(i, j)];
    return (double)total;
  }
  double rowSum(const std::vector<double>& m, int i) const {
    long double total = 0;
    for (int j = 0; j < nCol; ++j) total += m[idx(i, j)];
    return (double)total;
  }
};

struct Chain {
  int n;                                  // number of destinations
  std::vector<int> receiver;              // per resource; -1 while it is moved
  std::vector<std::string> resourceNames; // resource covariates in the cache
  std::vector<std::vector<double>> resourceValues;  // per covariate, per resource
  std::vector<std::string> originNames;
  std::vector<OriginNet> origins;

  // cache of the destinations (as cache[[dep.var]]$destinations in R)
  std::vector<double> inflow;
  std::vector<std::vector<double>> resourceInflows;

  // same logic as updateWeightedCache / the in-place .update in R: resource
  // `edge` is added to (u = 1) or removed from (u = -1) destination r
  void update(int edge, int r, int u) {
    for (OriginNet& o : origins) {
      int s = o.origin[edge];
      for (size_t k = 0; k < o.resourceNetworks.size(); ++k) {
        o.resourceNetworks[k][o.idx(s, r)] += u * resourceValues[k][edge];
      }
      o.valued[o.idx(s, r)] += u;
    }
    for (size_t k = 0; k < resourceInflows.size(); ++k) {
      resourceInflows[k][r] += u * resourceValues[k][edge];
    }
    inflow[r] += u;
  }

  // cells changed by an update of resource `edge` at destination r, to
  // restore them exactly; getCells fills `c`, reusing its memory, as it is
  // called in every step
  struct Cells {
    std::vector<double> valued;    // per origin
    std::vector<double> resource;  // per origin and covariate
    double inflow;
    std::vector<double> resourceInflow;
  };
  void getCells(int edge, int r, Cells& c) const {
    c.valued.clear();
    c.resource.clear();
    c.resourceInflow.clear();
    for (const OriginNet& o : origins) {
      int s = o.origin[edge];
      c.valued.push_back(o.v(s, r));
      for (const auto& m : o.resourceNetworks) c.resource.push_back(m[o.idx(s, r)]);
    }
    c.inflow = inflow[r];
    for (const auto& m : resourceInflows) c.resourceInflow.push_back(m[r]);
  }
  void restoreCells(int edge, int r, const Cells& c) {
    size_t cell = 0;
    for (size_t o = 0; o < origins.size(); ++o) {
      OriginNet& net = origins[o];
      int s = net.origin[edge];
      for (size_t k = 0; k < net.resourceNetworks.size(); ++k) {
        net.resourceNetworks[k][net.idx(s, r)] = c.resource[cell++];
      }
      net.valued[net.idx(s, r)] = c.valued[o];
    }
    for (size_t k = 0; k < resourceInflows.size(); ++k) {
      resourceInflows[k][r] = c.resourceInflow[k];
    }
    inflow[r] = c.inflow;
  }

  int resourceIndex(const std::string& name) const {
    for (size_t k = 0; k < resourceNames.size(); ++k) {
      if (resourceNames[k] == name) return (int)k;
    }
    Rcpp::stop("resource covariate '%s' is not in the cache", name);
  }
};

// An effect: change statistic of adding (update = 1) or removing (update = -1)
// resource `edge` on i -> j, and target contribution of the i -> j cell.
// For effects that depend on an origin, i is the location of the resource in
// that origin; effects that depend only on the destinations get i = -1.
class Effect {
 public:
  int o = -1;  // the effect's origin (index in Chain::origins), -1 for none
  bool excludeLoops = false;  // exclude.loops in R (effects that can exclude loops)
  virtual ~Effect() {}
  virtual double change(const Chain& ch, int i, int j, int edge, int update) const = 0;
  virtual double target(const Chain& ch, int i, int j) const = 0;
  // the target statistic: the sum of the target contributions over all cells,
  // in the order of getNetworkStatistics() (i fastest), or over the
  // destinations for effects without an origin; effects can override this
  // with a faster calculation that gives exactly the same value
  virtual double targetTotal(const Chain& ch) const {
    long double total = 0;
    if (o < 0) {
      for (int j = 0; j < ch.n; ++j) total += target(ch, -1, j);
      return (double)total;
    }
    int nRow = ch.origins[o].nRow;
    for (int j = 0; j < ch.n; ++j) {
      for (int i = 0; i < nRow; ++i) total += target(ch, i, j);
    }
    return (double)total;
  }
  // the location of resource `edge` in the effect's origin (-1 for none)
  inline int senderOf(const Chain& ch, int edge) const {
    return o < 0 ? -1 : ch.origins[o].origin[edge];
  }
 protected:
  inline const OriginNet& net(const Chain& ch) const { return ch.origins[o]; }
};

// builds an effect from its specification (see cppEffectSpec in R)
std::unique_ptr<Effect> makeEffect(const Rcpp::List& spec, const Chain& chain);

#endif
