// C++ versions of the effects, translated line by line from the R effect
// functions (R/2*_effectFunctions*.R). Indices are 0-based.
#include "realm.h"
#include <algorithm>
#include <functional>

namespace {

// loops (R/2b_effectFunctions_loops.R)
class Loops : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return update;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    return net(ch).v(i, j) * (i == j);
  }
};

// marginal and cumulative contributions of the AC effects
inline double gMarAC(double y, double a) {
  if (y > 0) return 1 - rpow(1 - 1 / a, y);
  return 0;
}
inline double gCumAC(double y, double a) {
  double contr = 0;
  if (y > 1) {
    for (int k = 1; k <= (int)(y - 1); ++k) contr = contr + (1 - rpow(1 - 1 / a, k));
  }
  return contr;
}

double alphaOf(const Rcpp::List& spec, const char* effect) {
  double alpha = Rcpp::as<double>(spec["alpha"]);
  if (alpha < 1) Rcpp::stop("alpha parameter in %s function must be 1 or larger", effect);
  return alpha;
}

std::vector<double> attributeOf(const Rcpp::List& spec) {
  Rcpp::NumericVector a = spec["attribute"];
  return std::vector<double>(a.begin(), a.end());
}

// concentration_AC (R/2c_effectFunctions_concentration.R)
class ConcentrationAC : public Effect {
  double alpha;
 public:
  explicit ConcentrationAC(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    double tie = net(ch).v(i, j);
    if (update < 0) return update * gMarAC(tie + update, alpha);
    return update * gMarAC(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return gCumAC(net(ch).v(i, j), alpha);
  }
};

// dyadic_covariate (R/2a_effectFunctions_exogenous.R)
class DyadicCovariate : public Effect {
  std::vector<double> attribute;  // n x n, column-major
 public:
  explicit DyadicCovariate(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return update * attribute[net(ch).idx(i, j)];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return net(ch).v(i, j) * attribute[net(ch).idx(i, j)];
  }
};

// alter_covariate (R/2a_effectFunctions_exogenous.R)
class AlterCovariate : public Effect {
  std::vector<double> attribute;
 public:
  explicit AlterCovariate(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    return attribute[j] * update;
  }
  double target(const Chain& ch, int i, int j) const override {
    return attribute[j] * ch.inflow[j];
  }
};

// resource_covar_to_node_covar (R/2a_effectFunctions_exogenous.R)
class ResourceCovarToNodeCovar : public Effect {
  int resource;
  std::vector<double> attribute;
 public:
  ResourceCovarToNodeCovar(int r, std::vector<double> a)
    : resource(r), attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    return update * ch.resourceValues[resource][edge] * attribute[j];
  }
  double target(const Chain& ch, int i, int j) const override {
    return attribute[j] * ch.resourceInflows[resource][j];
  }
};

// loops_resource_covar (R/2b_effectFunctions_loops.R)
class LoopsResourceCovar : public Effect {
  int resource;
 public:
  explicit LoopsResourceCovar(int r) : resource(r) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    return update * ch.resourceValues[resource][edge];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i != j) return 0;
    return net(ch).resourceNetworks[resource][net(ch).idx(i, j)];
  }
};

// loops_AC (R/2b_effectFunctions_loops.R)
class LoopsAC : public Effect {
  double alpha;
 public:
  explicit LoopsAC(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    double tie = net(ch).v(i, i);
    if (update < 0) return update * gMarAC(tie + update, alpha);
    return update * gMarAC(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return gCumAC(net(ch).v(i, i), alpha);
    return 0;
  }
};

// R's 0:y: counts up to floor(y) for y >= 0, and down to ceiling(y) otherwise
inline std::vector<int> rSeqFromZero(double y) {
  std::vector<int> k;
  if (y >= 0) {
    for (int m = 0; m <= (int)std::floor(y); ++m) k.push_back(m);
  } else {
    for (int m = 0; m >= (int)std::ceil(y); --m) k.push_back(m);
  }
  return k;
}

// marginal and cumulative contributions of the GW effects, summed term by term
// over 0:y as in R
inline double gMarGW(double y, double a) {
  double contr = 0;
  for (int k : rSeqFromZero(y)) contr = contr + std::exp(-std::log(a) * k);
  return contr - 1;
}
inline double gCumGW(double y, double a) {
  double contr = 0;
  for (int k : rSeqFromZero(y)) contr = contr + (y - k) * std::exp(-std::log(a) * k);
  return contr - y;
}

// loops_GW (R/2b_effectFunctions_loops.R)
class LoopsGW : public Effect {
  double alpha;
 public:
  explicit LoopsGW(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    double tie = net(ch).v(i, i);
    if (update < 0) return update * gMarGW(tie + update, alpha);
    return update * gMarGW(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return gCumGW(net(ch).v(i, i), alpha);
    return 0;
  }
};

// loops_additional_origin and loops_x_loops_additional_origin
// (R/2b_effectFunctions_loops.R): resources moving to their additional origin
// (a resource covariate holding a location id); with onlyLoops, only if this
// is also their origin
class LoopsAdditionalOrigin : public Effect {
  int resource;
  bool onlyLoops;
  // whether resource e counts for destination j (and, with onlyLoops, the
  // origin location i)
  bool counts(const Chain& ch, int e, int i, int j) const {
    if (onlyLoops && net(ch).origin[e] != i) return false;
    return ch.receiver[e] == j && ch.resourceValues[resource][e] == j + 1;
  }
 public:
  LoopsAdditionalOrigin(int r, bool loops) : resource(r), onlyLoops(loops) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (j + 1 == ch.resourceValues[resource][edge] && (!onlyLoops || j == i)) return update;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (onlyLoops && i != j) return 0;
    double cont = 0;
    for (size_t e = 0; e < ch.receiver.size(); ++e) cont += counts(ch, (int)e, i, j);
    return cont;
  }
  // one pass over the resources instead of one per cell; the same count
  double targetTotal(const Chain& ch) const override {
    double total = 0;
    for (size_t e = 0; e < ch.receiver.size(); ++e) {
      int j = ch.receiver[e];
      if (j < 0 || (onlyLoops && net(ch).origin[e] != j)) continue;
      total += ch.resourceValues[resource][e] == j + 1;
    }
    return total;
  }
};

// loops_node_covar (R/2b_effectFunctions_loops.R)
class LoopsNodeCovar : public Effect {
  std::vector<double> attribute;
 public:
  explicit LoopsNodeCovar(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    return update * attribute[i];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i != j) return 0;
    return net(ch).v(i, j) * attribute[i];
  }
};

// loops_resource_covar_node_covar (R/2b_effectFunctions_loops.R)
class LoopsResourceCovarNodeCovar : public Effect {
  int resource;
  std::vector<double> attribute;
 public:
  LoopsResourceCovarNodeCovar(int r, std::vector<double> a) : resource(r), attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return update * ch.resourceValues[resource][edge] * attribute[j];
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    return (net(ch).resourceNetworks[resource][net(ch).idx(i, j)] * (i == j)) * attribute[i];
  }
};

// concentration_GW, concentration_AC_dyad_covar, concentration_GW_dyad_covar,
// concentration_AC_resource_covar_bin, concentration_GW_resource_covar_bin
// (R/2c_effectFunctions_concentration.R): the AC or GW function of the tie
// value, which is the number of resources (weighted by a dyadic covariate) or
// the sum of a resource covariate on the tie
class Concentration : public Effect {
 public:
  enum Weight { NONE, DYAD, RESOURCE };
 private:
  bool gw;
  Weight weight;
  double alpha;
  std::vector<double> dyad;  // DYAD: n x n covariate
  int resource;              // RESOURCE: resource covariate
  double tieValue(const Chain& ch, int i, int j) const {
    if (weight == RESOURCE) return net(ch).resourceNetworks[resource][net(ch).idx(i, j)];
    return net(ch).v(i, j);
  }
  double gMar(double y) const { return gw ? gMarGW(y, alpha) : gMarAC(y, alpha); }
  double gCum(double y) const { return gw ? gCumGW(y, alpha) : gCumAC(y, alpha); }
 public:
  Concentration(bool isGW, Weight w, double a, std::vector<double> d, int r)
    : gw(isGW), weight(w), alpha(a), dyad(std::move(d)), resource(r) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    if (weight == DYAD && dyad[net(ch).idx(i, j)] == 0) return 0;
    if (weight == RESOURCE && ch.resourceValues[resource][edge] == 0) return 0;
    double tie = tieValue(ch, i, j);
    double g = update < 0 ? update * gMar(tie + update) : update * gMar(tie);
    if (weight == DYAD) return g * dyad[net(ch).idx(i, j)];
    return g;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    double g = gCum(tieValue(ch, i, j));
    if (weight == DYAD) return g * dyad[net(ch).idx(i, j)];
    return g;
  }
};

// in_weights_AC, in_weights_GW (R/2_effectFunctions.R): the AC or GW function
// of the inflow of the destination
class InWeights : public Effect {
  bool gw;
  double alpha;
 public:
  InWeights(bool isGW, double a) : gw(isGW), alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double w = ch.inflow[j];
    double y = update < 0 ? w + update : w;
    return update * (gw ? gMarGW(y, alpha) : gMarAC(y, alpha));
  }
  double target(const Chain& ch, int i, int j) const override {
    double w = ch.inflow[j];
    return gw ? gCumGW(w, alpha) : gCumAC(w, alpha);
  }
};

// ---- remaining concentration effects (R/2c_effectFunctions_concentration.R)

// concentration_basic (power 2), concentration_basic_cube (power 3)
class ConcentrationBasic : public Effect {
  double power;
 public:
  explicit ConcentrationBasic(double p) : power(p) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    return rpow(net(ch).v(i, j) + update, power) - rpow(net(ch).v(i, j), power);
  }
  double target(const Chain& ch, int i, int j) const override {
    return rpow(net(ch).v(i, j), power);
  }
};

// concentration_prop, concentration_prop_AC: functions of the proportions of
// the resources of an origin going to each destination
class ConcentrationProp : public Effect {
  bool ac;
  double alpha;
  double g(double k) const { return 1 - rpow(1 - 1 / alpha, k); }
  static double sum(const std::vector<double>& row) {
    long double total = 0;
    for (double x : row) total += x;
    return (double)total;
  }
  // contribution of all ties from one origin; empty origins contribute 0
  double rowContribution(const std::vector<double>& row) const {
    double total = sum(row);
    if (total == 0) return 0;
    long double s = 0;
    if (ac) {
      for (double x : row) s += x * g(x / total);
      return (double)s;
    }
    for (double x : row) s += rpow(x, 2);
    return (double)s / total;
  }
  std::vector<double> rowOf(const Chain& ch, int i) const {
    std::vector<double> row(ch.n);
    for (int k = 0; k < ch.n; ++k) row[k] = net(ch).v(i, k);
    return row;
  }
 public:
  ConcentrationProp(bool isAC, double a) : ac(isAC), alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    std::vector<double> row = rowOf(ch, i);
    std::vector<double> newRow = row;
    newRow[j] = newRow[j] + update;
    return rowContribution(newRow) - rowContribution(row);
  }
  double target(const Chain& ch, int i, int j) const override {
    std::vector<double> row = rowOf(ch, i);
    double total = sum(row);
    if (total == 0) return 0;
    if (ac) return row[j] * g(row[j] / total);
    return rpow(row[j], 2) / total;
  }
};

// concentration_max: the largest number of resources on a tie from each origin
class ConcentrationMax : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double v = net(ch).v(i, j);
    // largest and second-largest value of row i (R: sort(row, decreasing = T)[2])
    double first = R_NegInf, second = R_NegInf;
    for (int k = 0; k < ch.n; ++k) {
      double x = net(ch).v(i, k);
      if (x > first) {
        second = first;
        first = x;
      } else if (x > second) {
        second = x;
      }
    }
    if (update == 1 && first == v) return update;
    if (update == -1 && v > second) return update;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    double rowMax = R_NegInf;
    for (int k = 0; k < ch.n; ++k) rowMax = std::max(rowMax, net(ch).v(i, k));
    return rowMax / ch.n;
  }
};

// ---- concentration effects restored from MoNAn 1.1.0

// row i of the valued network of an origin
inline std::vector<double> rowOf(const OriginNet& net, int i) {
  std::vector<double> row(net.nCol);
  for (int k = 0; k < net.nCol; ++k) row[k] = net.v(i, k);
  return row;
}

// R: sum(x^2) (long double accumulation)
inline double sumSquares(const std::vector<double>& x, double divisor = 1) {
  long double total = 0;
  for (double v : x) total += divisor == 1 ? rpow(v, 2) : rpow(v, 2) / divisor;
  return (double)total;
}

// concentration_basic_squared
class ConcentrationBasicSquared : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    std::vector<double> previous = rowOf(net(ch), i);
    std::vector<double> next = previous;
    next[j] = next[j] + update;
    return rpow(sumSquares(next), 2) - rpow(sumSquares(previous), 2);
  }
  double target(const Chain& ch, int i, int j) const override {
    double vij2 = rpow(net(ch).v(i, j), 2);
    long double total = 0;
    for (int k = 0; k < ch.n; ++k) total += vij2 * rpow(net(ch).v(i, k), 2);
    return (double)total;
  }
};

// concentration_norm, concentration_norm_squared
class ConcentrationNorm : public Effect {
  bool squared;
 public:
  explicit ConcentrationNorm(bool sq) : squared(sq) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    std::vector<double> previous = rowOf(net(ch), i);
    long double s = 0;
    for (double v : previous) s += v;
    double total = (double)s;
    double vij = net(ch).v(i, j);
    if (!squared) {
      if (update == -1) return (rpow(vij - 1, 2) - rpow(vij, 2)) / (total - 1);
      if (update == 1) return (rpow(vij + 1, 2) - rpow(vij, 2)) / total;
      Rcpp::stop("updates other than -1 and 1 are not supported");
    }
    std::vector<double> next = previous;
    if (update == -1) {
      next[j] = next[j] - 1;
      return rpow(sumSquares(next, total - 1), 2) - rpow(sumSquares(previous, total - 1), 2);
    }
    if (update == 1) {
      next[j] = next[j] + 1;
      return rpow(sumSquares(next, total), 2) - rpow(sumSquares(previous, total), 2);
    }
    Rcpp::stop("updates other than -1 and 1 are not supported");
  }
  double target(const Chain& ch, int i, int j) const override {
    long double s = 0;
    for (int k = 0; k < ch.n; ++k) s += net(ch).v(i, k);
    double total = (double)s;
    double vij2 = rpow(net(ch).v(i, j), 2);
    if (!squared) return vij2 / (total - 1);
    double denominator = rpow(total - 1, 2);
    long double sum = 0;
    for (int k = 0; k < ch.n; ++k) sum += vij2 * rpow(net(ch).v(i, k), 2) / denominator;
    return (double)sum;
  }
};

// concentration_prop_orig_cov (the proportion statistic of MoNAn 1.1.0,
// weighted by an attribute of the origin)
class ConcentrationPropOrigCov : public Effect {
  std::vector<double> attribute;
  double rowSumOf(const Chain& ch, int i) const {
    long double s = 0;
    for (int k = 0; k < ch.n; ++k) s += net(ch).v(i, k);
    return (double)s;
  }
 public:
  explicit ConcentrationPropOrigCov(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double v = net(ch).v(i, j), total = rowSumOf(ch, i), cont;
    if (update == -1) {
      cont = 2 * (-2 * v + 1 + total) / (total - 1);
    } else if (update == 1) {
      cont = 2 * (2 * v - total) / total;
    } else {
      Rcpp::stop("updates other than -1 and 1 are not supported");
    }
    return attribute[i] * cont;
  }
  double target(const Chain& ch, int i, int j) const override {
    double v = net(ch).v(i, j), total = rowSumOf(ch, i);
    double numerator = 2 * rpow(v, 2) - v * total - v;
    double denominator = total - 1;
    if (denominator <= 0) {
      numerator = 0;
      denominator = 1;
    }
    return attribute[i] * numerator / denominator;
  }
};

// concentration_rankGW: weights of the ties from an origin by their rank
// (R: rank(-row), ties get their average rank)
class ConcentrationRankGW : public Effect {
  double logWeight;  // log(1 / lambda)
  // average ranks of all values in decreasing order (rank 1 = largest)
  static std::vector<double> ranks(const std::vector<double>& row) {
    std::vector<double> sorted = row;
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());
    std::vector<double> r(row.size());
    for (size_t k = 0; k < row.size(); ++k) {
      // number of values larger than row[k], and equal to it
      size_t larger = std::lower_bound(sorted.begin(), sorted.end(), row[k],
                                       std::greater<double>()) - sorted.begin();
      size_t equal = std::upper_bound(sorted.begin(), sorted.end(), row[k],
                                      std::greater<double>()) - sorted.begin() - larger;
      r[k] = larger + (equal + 1) / 2.0;
    }
    return r;
  }
  double w(double rank) const { return std::exp(rank * logWeight); }
 public:
  explicit ConcentrationRankGW(double lambda) : logWeight(std::log(1 / lambda)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    std::vector<double> previous = rowOf(net(ch), i);
    std::vector<double> next = previous;
    next[j] = next[j] + update;
    std::vector<double> r1 = ranks(previous), r2 = ranks(next);
    if (r1[j] == r2[j]) return next[j] * w(r2[j]) - previous[j] * w(r1[j]);
    double v = 0;
    for (int k = 0; k < ch.n; ++k) {
      if ((r1[k] - r2[k]) != 0) v = v + next[k] * w(r2[k]) - previous[k] * w(r1[k]);
    }
    return v;
  }
  double target(const Chain& ch, int i, int j) const override {
    double n = net(ch).v(i, j);
    if (n == 0) return 0;
    return n * w(ranks(rowOf(net(ch), i))[j]);
  }
};

// ---- remaining exogenous and other effects (R/2a_effectFunctions_exogenous.R,
// R/2_effectFunctions.R)

// same_covariate, sim_covariate: a dyadic covariate (the 'same' or 'sim'
// matrix of a node covariate) times the number of resources
// (also used for dyadic_covariate, which has the same form)

// same_covariate, sim_covariate with origin.attribute: the covariate a of the
// origin's locations compared with the covariate b of the destinations
class CovariateComparison : public Effect {
  bool similarity;  // false: same_covariate
  std::vector<double> a, b;
  double valueRange;
  double value(int i, int j) const {
    if (!similarity) return (a[i] == b[j]) * 1;
    return 1 - std::abs(a[i] - b[j]) / valueRange;
  }
 public:
  CovariateComparison(bool sim, std::vector<double> origins, std::vector<double> destinations)
    : similarity(sim), a(std::move(origins)), b(std::move(destinations)) {
    // R: range(c(a, b))
    double lo = R_PosInf, hi = R_NegInf;
    for (const std::vector<double>* v : {&a, &b}) {
      for (double x : *v) {
        lo = std::min(lo, x);
        hi = std::max(hi, x);
      }
    }
    valueRange = hi - lo;
  }
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j) return 0;
    return update * value(i, j);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j) return 0;
    return net(ch).v(i, j) * value(i, j);
  }
};

// dyadic_covariate_resource_attribute
class DyadicCovariateResourceAttribute : public Effect {
  std::vector<double> dyad;
  int resource;
 public:
  DyadicCovariateResourceAttribute(std::vector<double> d, int r) : dyad(std::move(d)), resource(r) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return update * ch.resourceValues[resource][edge] * dyad[net(ch).idx(i, j)];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return net(ch).resourceNetworks[resource][net(ch).idx(i, j)] * dyad[net(ch).idx(i, j)];
  }
};

// in_ties_loops
class InTiesLoops : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return update * (net(ch).colSum(net(ch).valued, j) - net(ch).v(j, j));
    return update * net(ch).v(j, j);
  }
  double target(const Chain& ch, int i, int j) const override {
    return (net(ch).v(i, j) * (i == j)) * (net(ch).colSum(net(ch).valued, j) - net(ch).v(j, j));
  }
};

// present_relations
class PresentRelations : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    if (net(ch).v(i, j) == 0 && update > 0) return 1;
    if (net(ch).v(i, j) == -update && update < 0) return -1;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (excludeLoops && i == j && net(ch).sameNodeSet) return 0;
    return (net(ch).v(i, j) > 0) * 1;
  }
};

// ---- covariate-based effects (R/2f_effectFunctions_endogenousCovariateBased.R)

// staying_by_prop_bin_inflow, crowding_out_prop_covar_bin: loops weighted by
// the proportion of resources with the (binary) attribute among the inflow
class InflowProportion : public Effect {
  int resource;
  bool crowdingOut;
  double r(const Chain& ch, int i, int j) const { return net(ch).resourceNetworks[resource][net(ch).idx(i, j)]; }
  double propBefore(const Chain& ch, int j) const {
    double sumV = net(ch).colSum(net(ch).valued, j);
    if ((sumV - net(ch).v(j, j)) == 0) return 0.5;
    return (net(ch).colSum(net(ch).resourceNetworks[resource], j) - r(ch, j, j)) / (sumV - net(ch).v(j, j));
  }
 public:
  InflowProportion(int res, bool crowding) : resource(res), crowdingOut(crowding) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double before = propBefore(ch, j);
    if (i == j) {
      if (crowdingOut && ch.resourceValues[resource][edge] == 1) return 0;
      return before * update;
    }
    double loops = crowdingOut ? net(ch).v(j, j) - r(ch, j, j) : net(ch).v(j, j);
    double sumV = net(ch).colSum(net(ch).valued, j);
    double after;
    if (((sumV + update) - net(ch).v(j, j)) == 0) {
      after = 0.5;
    } else {
      after = ((net(ch).colSum(net(ch).resourceNetworks[resource], j) +
                update * ch.resourceValues[resource][edge]) - r(ch, j, j)) /
              ((sumV + update) - net(ch).v(j, j));
    }
    double change = after - before;
    return change * loops;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i != j) return 0;
    double loops = crowdingOut ? net(ch).v(i, j) - r(ch, i, j) : net(ch).v(i, j);
    double sumV = net(ch).colSum(net(ch).valued, j);
    double prop;
    if ((sumV - net(ch).v(i, j)) == 0) {
      prop = 0.5;
    } else {
      prop = (net(ch).colSum(net(ch).resourceNetworks[resource], j) - r(ch, i, j)) / (sumV - net(ch).v(i, j));
    }
    return loops * prop;
  }
};

// joining_similar_avoiding_dissimilar_covar_bin, avoiding_dissimilar_covar_bin
class SimilarityBin : public Effect {
  int resource;
  bool joining;  // false: avoiding_dissimilar
 public:
  SimilarityBin(int r, bool j) : resource(r), joining(j) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double n = net(ch).v(i, j);
    double n1 = net(ch).resourceNetworks[resource][net(ch).idx(i, j)];
    double attr = ch.resourceValues[resource][edge];
    if (update == -1) {
      if (n < 2) return 0;
      double originSize = net(ch).rowSum(net(ch).valued, i);
      if (joining) {
        if (attr == 1) return (-2 * (n1 - 1) + 2 * (n - n1)) / (originSize - 1);
        return (-2 * (n - n1 - 1) + 2 * n1) / (originSize - 1);
      }
      if (attr == 1) return (2 * (n - n1)) / (originSize - 1);
      return (2 * n1) / (originSize - 1);
    }
    if (update == 1) {
      if (n < 1) return 0;
      double originSize = net(ch).rowSum(net(ch).valued, i) + 1;
      if (joining) {
        if (attr == 1) return (2 * n1 - 2 * (n - n1)) / (originSize - 1);
        return (2 * (n - n1) - 2 * (n1)) / (originSize - 1);
      }
      if (attr == 1) return (-2 * (n - n1)) / (originSize - 1);
      return (-2 * (n1)) / (originSize - 1);
    }
    Rcpp::stop("updates other than -1 and 1 are not supported");
  }
  double target(const Chain& ch, int i, int j) const override {
    double n = net(ch).v(i, j);
    if (n < 2) return 0;
    double n1 = net(ch).resourceNetworks[resource][net(ch).idx(i, j)];
    double originSize = net(ch).rowSum(net(ch).valued, i);
    double numerator;
    if (joining) {
      numerator = n1 * (n1 - 1) + (n - n1) * (n - n1 - 1) - 2 * n1 * (n - n1);
    } else {
      numerator = -2 * n1 * (n - n1);
    }
    double denominator = originSize - 1;
    return numerator / denominator;
  }
};

// joining_similar_avoiding_dissimilar_covar_cont, avoiding_dissimilar_covar_cont:
// pairwise distances of a continuous attribute among the resources on a tie
class SimilarityCont : public Effect {
  int resource;
  bool joining;  // false: avoiding_dissimilar
  double attributeRange;
  // R: the pairwise term of the target statistic, for the distance d = |a - b|
  double term(double d) const {
    return joining ? (attributeRange - 2 * d) / attributeRange : d;
  }
  // the target contribution of one tie from the attributes of its resources
  // (in edgelist order), summing over the pairs in the order of upper.tri()
  double tieValue(const Chain& ch, int i, const std::vector<double>& x) const {
    long double numerator = 0;
    for (size_t c = 1; c < x.size(); ++c) {
      for (size_t r = 0; r < c; ++r) numerator += term(std::abs(x[r] - x[c]));
    }
    double denominator = net(ch).rowSum(net(ch).valued, i) - 1;
    return (double)numerator / denominator;
  }
 public:
  SimilarityCont(int r, bool j, double range) : resource(r), joining(j), attributeRange(range) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double n = net(ch).v(i, j);
    if (update == -1 && n < 2) return 0;
    if (update == 1 && n < 1) return 0;
    if (update != -1 && update != 1) Rcpp::stop("updates other than -1 and 1 are not supported");
    double attr = ch.resourceValues[resource][edge];
    long double total = 0;
    for (size_t e = 0; e < ch.receiver.size(); ++e) {
      if (net(ch).origin[e] == i && ch.receiver[e] == j && (int)e != edge) {
        total += term(std::abs(attr - ch.resourceValues[resource][e]));
      }
    }
    double originSize = net(ch).rowSum(net(ch).valued, i) + (update == 1 ? 1 : 0);
    if (update == -1) return -(double)total / (originSize - 1);
    return (double)total / (originSize - 1);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (net(ch).v(i, j) < 2) return 0;
    std::vector<double> x;
    for (size_t e = 0; e < ch.receiver.size(); ++e) {
      if (net(ch).origin[e] == i && ch.receiver[e] == j) x.push_back(ch.resourceValues[resource][e]);
    }
    return tieValue(ch, i, x);
  }
  // one pass over the resources instead of one per cell; the same values
  double targetTotal(const Chain& ch) const override {
    const OriginNet& N = net(ch);
    std::vector<std::vector<double>> byTie((size_t)N.nRow * ch.n);
    for (size_t e = 0; e < ch.receiver.size(); ++e) {
      if (ch.receiver[e] < 0) continue;
      byTie[N.idx(N.origin[e], ch.receiver[e])].push_back(ch.resourceValues[resource][e]);
    }
    long double total = 0;
    for (int j = 0; j < ch.n; ++j) {
      for (int i = 0; i < N.nRow; ++i) {
        total += net(ch).v(i, j) < 2 ? 0 : tieValue(ch, i, byTie[net(ch).idx(i, j)]);
      }
    }
    return (double)total;
  }
};

}  // namespace

std::unique_ptr<Effect> makeEffectFor(const Rcpp::List& spec, const Chain& chain) {
  std::string name = Rcpp::as<std::string>(spec["name"]);
  typedef std::unique_ptr<Effect> E;
  if (name == "loops") return E(new Loops());
  if (name == "concentration_AC") {
    return E(new ConcentrationAC(alphaOf(spec, "concentration_AC")));
  }
  if (name == "dyadic_covariate") return E(new DyadicCovariate(attributeOf(spec)));
  if (name == "alter_covariate") {
    return E(new AlterCovariate(attributeOf(spec)));
  }
  if (name == "resource_covar_to_node_covar") {
    return E(new ResourceCovarToNodeCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      attributeOf(spec)));
  }
  if (name == "loops_resource_covar") {
    return E(new LoopsResourceCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]))));
  }
  if (name == "loops_AC" || name == "loops_GW") {
    double alpha = Rcpp::as<double>(spec["alpha"]);
    if (alpha <= 0) Rcpp::stop("Alpha parameter in GW loops weights function must be positive");
    if (name == "loops_AC") return E(new LoopsAC(alpha));
    return E(new LoopsGW(alpha));
  }
  if (name == "loops_additional_origin" || name == "loops_x_loops_additional_origin") {
    return E(new LoopsAdditionalOrigin(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      name == "loops_x_loops_additional_origin"));
  }
  if (name == "loops_node_covar") return E(new LoopsNodeCovar(attributeOf(spec)));
  if (name == "loops_resource_covar_node_covar") {
    return E(new LoopsResourceCovarNodeCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      attributeOf(spec)));
  }
  if (name == "concentration_GW" || name == "concentration_AC_dyad_covar" ||
      name == "concentration_GW_dyad_covar" ||
      name == "concentration_AC_resource_covar_bin" ||
      name == "concentration_GW_resource_covar_bin") {
    bool gw = name.find("_GW") != std::string::npos;
    Concentration::Weight w = Concentration::NONE;
    std::vector<double> dyad;
    int resource = -1;
    // the R versions name the GW dyadic and unweighted effects
    // "concentration_GW" in their alpha error message
    std::string alphaName = name == "concentration_GW_dyad_covar" ? "concentration_GW" : name;
    double alpha = alphaOf(spec, alphaName.c_str());
    if (name.find("dyad_covar") != std::string::npos) {
      w = Concentration::DYAD;
      dyad = attributeOf(spec);
      int origin = Rcpp::as<int>(spec["origin"]);
      if (dyad.size() != (size_t)chain.origins[origin].nRow * chain.n) {
        Rcpp::stop("attribute.index in %s function must be a matrix of origins x destinations", name);
      }
    }
    if (name.find("resource_covar") != std::string::npos) {
      w = Concentration::RESOURCE;
      resource = chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]));
    }
    return E(new Concentration(gw, w, alpha, dyad, resource));
  }
  if (name == "in_weights_AC" || name == "in_weights_GW") {
    double alpha = Rcpp::as<double>(spec["alpha"]);
    if (alpha <= 0) {
      Rcpp::stop("Alpha parameter in %s weights function must be positive", name);
    }
    return E(new InWeights(name == "in_weights_GW", alpha));
  }
  if (name == "concentration_basic") return E(new ConcentrationBasic(2));
  if (name == "concentration_basic_cube") return E(new ConcentrationBasic(3));
  if (name == "concentration_max") return E(new ConcentrationMax());
  if (name == "concentration_basic_squared") return E(new ConcentrationBasicSquared());
  if (name == "concentration_norm") return E(new ConcentrationNorm(false));
  if (name == "concentration_norm_squared") return E(new ConcentrationNorm(true));
  if (name == "concentration_prop_orig_cov") {
    if (Rcpp::as<bool>(spec["attributeIsMatrix"])) {
      Rcpp::stop("attribute in concentration_prop_orig_cov must be monadic");
    }
    return E(new ConcentrationPropOrigCov(attributeOf(spec)));
  }
  if (name == "concentration_rankGW") {
    double lambda = Rcpp::as<double>(spec["lambda"]);
    if (lambda <= 0) Rcpp::stop("lambda parameter in concentration_rankGW function must be positive");
    return E(new ConcentrationRankGW(lambda));
  }
  if (name == "concentration_prop") return E(new ConcentrationProp(false, 0));
  if (name == "concentration_prop_AC") {
    return E(new ConcentrationProp(true, alphaOf(spec, "concentration_prop_AC")));
  }
  if (name == "same_covariate" || name == "sim_covariate") {
    if (spec.containsElementNamed("originAttribute")) {
      Rcpp::NumericVector a = spec["originAttribute"];
      return E(new CovariateComparison(name == "sim_covariate",
                                       std::vector<double>(a.begin(), a.end()),
                                       attributeOf(spec)));
    }
    return E(new DyadicCovariate(attributeOf(spec)));
  }
  if (name == "dyadic_covariate_resource_attribute") {
    return E(new DyadicCovariateResourceAttribute(
      attributeOf(spec), chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]))));
  }
  if (name == "in_ties_loops") return E(new InTiesLoops());
  if (name == "present_relations") return E(new PresentRelations());
  if (name == "staying_by_prop_bin_inflow" || name == "crowding_out_prop_covar_bin") {
    return E(new InflowProportion(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      name == "crowding_out_prop_covar_bin"));
  }
  if (name == "joining_similar_avoiding_dissimilar_covar_bin" ||
      name == "avoiding_dissimilar_covar_bin") {
    return E(new SimilarityBin(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      name == "joining_similar_avoiding_dissimilar_covar_bin"));
  }
  if (name == "joining_similar_avoiding_dissimilar_covar_cont" ||
      name == "avoiding_dissimilar_covar_cont") {
    int r = chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]));
    double lo = R_PosInf, hi = R_NegInf;
    for (double x : chain.resourceValues[r]) {
      lo = std::min(lo, x);
      hi = std::max(hi, x);
    }
    return E(new SimilarityCont(r, name == "joining_similar_avoiding_dissimilar_covar_cont", hi - lo));
  }
  Rcpp::stop("effect '%s' has no C++ version", name);
}

std::unique_ptr<Effect> makeEffect(const Rcpp::List& spec, const Chain& chain) {
  std::unique_ptr<Effect> effect = makeEffectFor(spec, chain);
  effect->o = Rcpp::as<int>(spec["origin"]);
  effect->excludeLoops = spec.containsElementNamed("excludeLoops") &&
                         Rcpp::as<bool>(spec["excludeLoops"]);
  return effect;
}
