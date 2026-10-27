// C++ implementation of the MCMC sampler (see R/1b_cppSampler.R).
// One step follows simulateOneStep() in R/1_hiddenFunctions.R exactly,
// including the order of random numbers, so that for the same seed the C++
// chain equals the R chain.
#include "realm.h"

#include <algorithm>
#include <map>

namespace {

// R's sample(x, 1): a single number x >= 1 means sample(1:x, 1)
int sampleLikeR(const std::vector<int>& x) {
  if (x.size() == 1 && x[0] >= 1) return (int)R_unif_index((double)x[0]) + 1;
  return x[(size_t)R_unif_index((double)x.size())];
}

std::vector<double> asVector(SEXP m) {
  Rcpp::NumericVector v(m);
  return std::vector<double>(v.begin(), v.end());
}

Rcpp::NumericMatrix asMatrix(const std::vector<double>& v, int nRow, int nCol) {
  Rcpp::NumericMatrix m(nRow, nCol);
  std::copy(v.begin(), v.end(), m.begin());
  return m;
}

class Sampler {
 public:
  Chain ch;
  std::vector<std::unique_ptr<Effect>> effects;
  std::vector<double> parameters;
  std::vector<int> resourceSet, nodeSet;
  bool multinomial;
  std::vector<int> senders;  // the origin of the moved resource for each effect
  std::vector<double> drop, create;
  Chain::Cells old;          // cells to restore if a proposal is rejected

  // the origin of resource `edge` for each effect, as effectSenders() in R
  void setSenders(int edge) {
    senders.resize(effects.size());
    for (size_t k = 0; k < effects.size(); ++k) senders[k] = effects[k]->senderOf(ch, edge);
  }

  void step() {
    int edge = sampleLikeR(resourceSet) - 1;
    int j = ch.receiver[edge];
    setSenders(edge);

    size_t K = effects.size();

    if (multinomial) {
      multinomialStep(edge, j);
      return;
    }

    int l = sampleLikeR(nodeSet) - 1;
    if (j == l) return;

    drop.resize(K);
    create.resize(K);
    for (size_t k = 0; k < K; ++k) drop[k] = effects[k]->change(ch, senders[k], j, edge, -1);

    // remove the resource, keeping what is needed to undo this if rejected
    ch.getCells(edge, j, old);
    ch.update(edge, j, -1);
    ch.receiver[edge] = -1;

    for (size_t k = 0; k < K; ++k) create[k] = effects[k]->change(ch, senders[k], l, edge, 1);

    long double contribution = 0;
    for (size_t k = 0; k < K; ++k) contribution += (drop[k] + create[k]) * parameters[k];
    double acceptance = std::min(1.0, std::exp((double)contribution));

    if (acceptance >= unif_rand()) {
      ch.update(edge, l, 1);
      ch.receiver[edge] = l;
    } else {
      ch.restoreCells(edge, j, old);
      ch.receiver[edge] = j;
    }
  }

  // multinomial proposal: the destination is drawn from all possible ones with
  // probabilities proportional to exp(parameters . change statistics), with the
  // same arithmetic as R (colSums, sum and cumsum accumulate in long double)
  void multinomialStep(int edge, int j) {
    size_t K = effects.size();
    drop.resize(K);
    for (size_t k = 0; k < K; ++k) drop[k] = effects[k]->change(ch, senders[k], j, edge, -1);

    ch.update(edge, j, -1);
    ch.receiver[edge] = -1;

    size_t D = nodeSet.size();
    std::vector<double> expValues(D);
    long double total = 0;
    for (size_t d = 0; d < D; ++d) {
      long double objective = 0;
      for (size_t k = 0; k < K; ++k) {
        double create = effects[k]->change(ch, senders[k], nodeSet[d] - 1, edge, 1);
        objective += (create + drop[k]) * parameters[k];
      }
      expValues[d] = std::exp((double)objective);
      total += expValues[d];
    }

    double u = unif_rand();
    long double cumulative = 0;
    size_t draw = D;
    for (size_t d = 0; d < D; ++d) {
      cumulative += expValues[d] / (double)total;
      if ((double)cumulative >= u) {
        draw = d;
        break;
      }
    }
    if (draw == D) Rcpp::stop("no destination drawn in multinomial proposal");
    int l = nodeSet[draw] - 1;

    ch.update(edge, l, 1);
    ch.receiver[edge] = l;
  }

  // as getNetworkStatistics(): sum of target contributions over all cells,
  // in the order of expand.grid (i fastest), or over the destinations
  std::vector<double> statistics() const {
    std::vector<double> s(effects.size());
    for (size_t k = 0; k < effects.size(); ++k) s[k] = effects[k]->targetTotal(ch);
    return s;
  }

  Rcpp::IntegerVector destinations() const {
    Rcpp::IntegerVector d(ch.receiver.size());
    for (size_t e = 0; e < ch.receiver.size(); ++e) d[e] = ch.receiver[e] + 1;
    return d;
  }

  // the cache as createInternalCache() in R
  Rcpp::List cache() const {
    Rcpp::List inflows(ch.resourceNames.size());
    for (size_t k = 0; k < ch.resourceNames.size(); ++k) {
      inflows[k] = Rcpp::NumericVector(ch.resourceInflows[k].begin(), ch.resourceInflows[k].end());
    }
    inflows.attr("names") = ch.resourceNames;
    Rcpp::List out(ch.origins.size() + 1);
    std::vector<std::string> names = {"destinations"};
    out[0] = Rcpp::List::create(
      Rcpp::Named("inflow") = Rcpp::NumericVector(ch.inflow.begin(), ch.inflow.end()),
      Rcpp::Named("resourceInflows") = inflows);
    for (size_t o = 0; o < ch.origins.size(); ++o) {
      const OriginNet& net = ch.origins[o];
      Rcpp::List resource(ch.resourceNames.size());
      for (size_t k = 0; k < ch.resourceNames.size(); ++k) {
        resource[k] = asMatrix(net.resourceNetworks[k], net.nRow, net.nCol);
      }
      resource.attr("names") = ch.resourceNames;
      out[o + 1] = Rcpp::List::create(
        Rcpp::Named("valuedNetwork") = asMatrix(net.valued, net.nRow, net.nCol),
        Rcpp::Named("resourceNetworks") = resource,
        Rcpp::Named("sameNodeSet") = net.sameNodeSet);
      names.push_back(ch.originNames[o]);
    }
    out.attr("names") = names;
    return out;
  }
};

// sets up chain and effects from the R objects (see cppChainInputs() in R):
// the destinations, the origins (a named list) and the cache of the
// dependent variable
void setupSampler(Sampler& s,
                  const Rcpp::IntegerVector& destinations,
                  const Rcpp::List& origins,
                  const Rcpp::List& cache,
                  const Rcpp::List& resourceCovariates,
                  const Rcpp::List& effectSpecs) {
  Rcpp::List dest = cache["destinations"];
  s.ch.inflow = asVector(dest["inflow"]);
  s.ch.n = (int)s.ch.inflow.size();
  for (R_xlen_t e = 0; e < destinations.size(); ++e) s.ch.receiver.push_back(destinations[e] - 1);
  Rcpp::List inflows = dest["resourceInflows"];
  if (inflows.size() > 0) {
    Rcpp::CharacterVector names = inflows.names();
    for (int k = 0; k < inflows.size(); ++k) {
      std::string name = Rcpp::as<std::string>(names[k]);
      s.ch.resourceNames.push_back(name);
      s.ch.resourceInflows.push_back(asVector(inflows[k]));
      s.ch.resourceValues.push_back(asVector(resourceCovariates[name]));
    }
  }
  if (origins.size() > 0) {
    Rcpp::CharacterVector names = origins.names();
    for (int o = 0; o < origins.size(); ++o) {
      std::string name = Rcpp::as<std::string>(names[o]);
      Rcpp::List oc = cache[name];
      Rcpp::NumericMatrix valued = oc["valuedNetwork"];
      OriginNet net;
      net.nRow = valued.nrow();
      net.nCol = s.ch.n;
      net.sameNodeSet = Rcpp::as<bool>(oc["sameNodeSet"]);
      Rcpp::IntegerVector origin = origins[o];
      for (R_xlen_t e = 0; e < origin.size(); ++e) net.origin.push_back(origin[e] - 1);
      net.valued = asVector(valued);
      Rcpp::List rn = oc["resourceNetworks"];
      for (const std::string& r : s.ch.resourceNames) net.resourceNetworks.push_back(asVector(rn[r]));
      s.ch.originNames.push_back(name);
      s.ch.origins.push_back(net);
    }
  }
  for (int k = 0; k < effectSpecs.size(); ++k) {
    s.effects.push_back(makeEffect(Rcpp::as<Rcpp::List>(effectSpecs[k]), s.ch));
  }
}

}  // namespace

// Runs one chain: burnIn steps, then `iterations` times (thinning + 1) steps
// followed by the calculation of the statistics, as simulateStatisticVectors().
// With iterations = 0, this is simulateNSteps(n = burnIn).
// [[Rcpp::export]]
Rcpp::List simulateChainCpp(Rcpp::IntegerVector destinations,
                            Rcpp::List origins,
                            Rcpp::List cache,
                            Rcpp::List resourceCovariates,
                            Rcpp::IntegerVector resourceSet,
                            Rcpp::IntegerVector nodeSet,
                            Rcpp::List effectSpecs,
                            Rcpp::NumericVector parameters,
                            int burnIn,
                            int iterations,
                            int thinning,
                            bool multinomialProposal,
                            bool returnDraws) {
  Sampler s;
  setupSampler(s, destinations, origins, cache, resourceCovariates, effectSpecs);
  s.parameters.assign(parameters.begin(), parameters.end());
  if (s.parameters.size() != s.effects.size()) {
    Rcpp::stop("number of parameters and effects differ");
  }
  s.resourceSet.assign(resourceSet.begin(), resourceSet.end());
  s.nodeSet.assign(nodeSet.begin(), nodeSet.end());
  s.multinomial = multinomialProposal;

  for (int t = 0; t < burnIn; ++t) s.step();

  Rcpp::NumericMatrix stats(iterations, (int)s.effects.size());
  Rcpp::List draws(returnDraws ? iterations : 0);
  for (int it = 0; it < iterations; ++it) {
    for (int t = 0; t < thinning + 1; ++t) s.step();
    std::vector<double> st = s.statistics();
    for (size_t k = 0; k < st.size(); ++k) stats(it, k) = st[k];
    if (returnDraws) {
      draws[it] = Rcpp::List::create(Rcpp::Named("destinations") = s.destinations(),
                                     Rcpp::Named("cache") = s.cache());
    }
  }

  return Rcpp::List::create(Rcpp::Named("statistics") = stats,
                            Rcpp::Named("destinations") = s.destinations(),
                            Rcpp::Named("cache") = s.cache(),
                            Rcpp::Named("draws") = draws);
}

// target statistics of a state, as getNetworkStatistics()
// [[Rcpp::export]]
Rcpp::NumericVector networkStatisticsCpp(Rcpp::IntegerVector destinations,
                                         Rcpp::List origins,
                                         Rcpp::List cache,
                                         Rcpp::List resourceCovariates,
                                         Rcpp::List effectSpecs) {
  Sampler s;
  setupSampler(s, destinations, origins, cache, resourceCovariates, effectSpecs);
  std::vector<double> st = s.statistics();
  return Rcpp::NumericVector(st.begin(), st.end());
}

// change statistic and target contribution of one effect, for tests; i is NA
// for effects without an origin
// [[Rcpp::export]]
Rcpp::NumericVector effectValuesCpp(Rcpp::IntegerVector destinations,
                                    Rcpp::List origins,
                                    Rcpp::List cache,
                                    Rcpp::List resourceCovariates,
                                    Rcpp::List effectSpec,
                                    int i, int j, int edge, int update) {
  Sampler s;
  setupSampler(s, destinations, origins, cache, resourceCovariates, Rcpp::List::create(effectSpec));
  int sender = i == NA_INTEGER ? -1 : i - 1;
  return Rcpp::NumericVector::create(
    Rcpp::Named("change") = s.effects[0]->change(s.ch, sender, j - 1, edge - 1, update),
    Rcpp::Named("target") = s.effects[0]->target(s.ch, sender, j - 1));
}

// sums of a resource covariate per tie (or per destination), adding the values
// in the order of the resources starting from 0 (for createInternalCache in R)
// [[Rcpp::export]]
Rcpp::NumericVector tieSumsCpp(Rcpp::IntegerVector tie, Rcpp::NumericVector values, int nCells) {
  Rcpp::NumericVector sums(nCells);
  for (R_xlen_t e = 0; e < tie.size(); ++e) sums[tie[e] - 1] += values[e];
  return sums;
}

// Log pseudo-likelihood of the conditional logit model of each resource's
// destination given all other resources (as getMPLE() in R), and optionally
// its gradient and Hessian, for parameters theta. For each resource, the
// change statistics of adding it to every possible destination are computed
// after removing it (as in the multinomial proposal); nothing is stored.
// [[Rcpp::export]]
Rcpp::List mpleNewtonCpp(Rcpp::IntegerVector destinations,
                         Rcpp::List origins,
                         Rcpp::List cache,
                         Rcpp::List resourceCovariates,
                         Rcpp::IntegerVector resourceSet,
                         Rcpp::IntegerVector nodeSet,
                         Rcpp::List effectSpecs,
                         Rcpp::NumericVector theta,
                         bool derivatives) {
  Sampler s;
  setupSampler(s, destinations, origins, cache, resourceCovariates, effectSpecs);
  size_t K = s.effects.size();
  if ((size_t)theta.size() != K) Rcpp::stop("number of parameters and effects differ");

  long double logPL = 0;
  std::vector<double> gradient(K, 0.0), hessian(K * K, 0.0);
  std::vector<double> X, eta, prob, expected(K);
  size_t D = nodeSet.size();

  for (int r : resourceSet) {
    int edge = r - 1;
    int j = s.ch.receiver[edge];
    s.setSenders(edge);
    int observed = -1;
    for (size_t d = 0; d < D; ++d) if (nodeSet[d] - 1 == j) observed = (int)d;
    if (observed < 0) {
      // the observed destination is not possible (as match() giving NA in R)
      return Rcpp::List::create(Rcpp::Named("logPL") = NA_REAL);
    }

    // change statistics of adding the resource to each destination
    s.ch.getCells(edge, j, s.old);
    s.ch.update(edge, j, -1);
    s.ch.receiver[edge] = -1;
    X.assign(D * K, 0.0);
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) {
        X[d * K + k] = s.effects[k]->change(s.ch, s.senders[k], nodeSet[d] - 1, edge, 1);
      }
    }
    s.ch.restoreCells(edge, j, s.old);
    s.ch.receiver[edge] = j;

    // conditional logit contributions (softmax with the maximum subtracted)
    eta.assign(D, 0.0);
    double maxEta = R_NegInf;
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) eta[d] += X[d * K + k] * theta[k];
      maxEta = std::max(maxEta, eta[d]);
    }
    double total = 0;
    prob.assign(D, 0.0);
    for (size_t d = 0; d < D; ++d) {
      prob[d] = std::exp(eta[d] - maxEta);
      total += prob[d];
    }
    logPL += eta[observed] - maxEta - std::log(total);
    if (!derivatives) continue;

    for (size_t d = 0; d < D; ++d) prob[d] /= total;
    std::fill(expected.begin(), expected.end(), 0.0);
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) expected[k] += prob[d] * X[d * K + k];
    }
    for (size_t k = 0; k < K; ++k) gradient[k] += X[observed * K + k] - expected[k];
    for (size_t a = 0; a < K; ++a) {
      for (size_t b = 0; b < K; ++b) {
        double second = 0;
        for (size_t d = 0; d < D; ++d) second += prob[d] * X[d * K + a] * X[d * K + b];
        hessian[a * K + b] -= second - expected[a] * expected[b];
      }
    }
  }

  Rcpp::NumericMatrix H(K, K);
  for (size_t a = 0; a < K; ++a) {
    for (size_t b = 0; b < K; ++b) H(a, b) = hessian[a * K + b];
  }
  return Rcpp::List::create(
    Rcpp::Named("logPL") = (double)logPL,
    Rcpp::Named("gradient") = Rcpp::NumericVector(gradient.begin(), gradient.end()),
    Rcpp::Named("hessian") = H);
}

// Choice sets of the conditional logit model of the MPLE (as mpleNewtonCpp),
// computed once and stored, each distinct choice set only once. Resources with
// the same origins, the same observed destination and the same values of the
// resource covariates in `keyCovariates` (those used by the model's effects)
// have the same change statistics, so one of them represents the group, with
// the group's size as weight. Returns possible = FALSE if an observed
// destination is not a possible destination (as match() giving NA in R), and
// tooLarge = TRUE, without computing the change statistics, if they would need
// more than maxBytes.
// [[Rcpp::export]]
Rcpp::List mpleChoiceSetsCpp(Rcpp::IntegerVector destinations,
                             Rcpp::List origins,
                             Rcpp::List cache,
                             Rcpp::List resourceCovariates,
                             Rcpp::IntegerVector resourceSet,
                             Rcpp::IntegerVector nodeSet,
                             Rcpp::List effectSpecs,
                             Rcpp::CharacterVector keyCovariates,
                             double maxBytes) {
  Sampler s;
  setupSampler(s, destinations, origins, cache, resourceCovariates, effectSpecs);
  size_t K = s.effects.size();
  std::vector<int> keyIndex;
  for (R_xlen_t k = 0; k < keyCovariates.size(); ++k) {
    keyIndex.push_back(s.ch.resourceIndex(Rcpp::as<std::string>(keyCovariates[k])));
  }

  // group the resources, in the order of resourceSet
  std::map<std::vector<double>, size_t> groupOf;
  std::vector<int> representative;
  std::vector<double> weight;
  for (int r : resourceSet) {
    int edge = r - 1;
    int j = s.ch.receiver[edge];
    bool possible = std::find(nodeSet.begin(), nodeSet.end(), j + 1) != nodeSet.end();
    if (!possible) {
      return Rcpp::List::create(Rcpp::Named("possible") = false);
    }
    std::vector<double> key;
    for (const OriginNet& net : s.ch.origins) key.push_back((double)net.origin[edge]);
    key.push_back((double)j);
    bool missing = false;
    for (int k : keyIndex) {
      double value = s.ch.resourceValues[k][edge];
      if (ISNAN(value)) missing = true;
      key.push_back(value);
    }
    if (!missing) {
      auto found = groupOf.find(key);
      if (found != groupOf.end()) {
        weight[found->second] += 1;
        continue;
      }
      groupOf[key] = representative.size();
    }
    // a new group (resources with missing values are not grouped)
    representative.push_back(edge);
    weight.push_back(1);
  }

  size_t G = representative.size();
  double rows = (double)G * nodeSet.size();
  if (rows * K * sizeof(double) > maxBytes) {
    return Rcpp::List::create(Rcpp::Named("possible") = true,
                              Rcpp::Named("tooLarge") = true);
  }

  // change statistics of adding each representative to each destination,
  // after removing it; one column per destination (effects in rows)
  Rcpp::NumericMatrix X((int)K, (int)rows);
  Rcpp::IntegerVector start((int)G + 1), observed((int)G);
  int row = 0;
  for (size_t g = 0; g < G; ++g) {
    int edge = representative[g];
    int j = s.ch.receiver[edge];
    s.setSenders(edge);
    start[g] = row;
    s.ch.getCells(edge, j, s.old);
    s.ch.update(edge, j, -1);
    s.ch.receiver[edge] = -1;
    for (int k : nodeSet) {
      int destination = k - 1;
      if (destination == j) observed[g] = row;
      for (size_t e = 0; e < K; ++e) {
        X(e, row) = s.effects[e]->change(s.ch, s.senders[e], destination, edge, 1);
      }
      ++row;
    }
    s.ch.restoreCells(edge, j, s.old);
    s.ch.receiver[edge] = j;
  }
  start[G] = row;

  return Rcpp::List::create(
    Rcpp::Named("possible") = true,
    Rcpp::Named("tooLarge") = false,
    Rcpp::Named("X") = X,
    Rcpp::Named("start") = start,
    Rcpp::Named("observed") = observed,
    Rcpp::Named("weight") = Rcpp::NumericVector(weight.begin(), weight.end()));
}

// Log pseudo-likelihood, and optionally gradient and Hessian, from the stored
// choice sets of mpleChoiceSetsCpp: as mpleNewtonCpp, with each choice set's
// contribution multiplied by its weight.
// [[Rcpp::export]]
Rcpp::List mpleEvaluateCpp(Rcpp::NumericMatrix X,
                           Rcpp::IntegerVector start,
                           Rcpp::IntegerVector observed,
                           Rcpp::NumericVector weight,
                           Rcpp::NumericVector theta,
                           bool derivatives) {
  size_t K = X.nrow();
  if ((size_t)theta.size() != K) Rcpp::stop("number of parameters and effects differ");
  const double* x = X.begin();
  long double logPL = 0;
  std::vector<double> gradient(K, 0.0), hessian(K * K, 0.0);
  std::vector<double> eta, prob, expected(K);

  for (R_xlen_t g = 0; g < weight.size(); ++g) {
    int first = start[g], D = start[g + 1] - start[g];
    int obs = observed[g] - first;
    double w = weight[g];
    const double* xg = x + (size_t)first * K;

    // conditional logit contributions (softmax with the maximum subtracted)
    eta.assign(D, 0.0);
    double maxEta = R_NegInf;
    for (int d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) eta[d] += xg[d * K + k] * theta[k];
      maxEta = std::max(maxEta, eta[d]);
    }
    double total = 0;
    prob.assign(D, 0.0);
    for (int d = 0; d < D; ++d) {
      prob[d] = std::exp(eta[d] - maxEta);
      total += prob[d];
    }
    logPL += w * (eta[obs] - maxEta - std::log(total));
    if (!derivatives) continue;

    for (int d = 0; d < D; ++d) prob[d] /= total;
    std::fill(expected.begin(), expected.end(), 0.0);
    for (int d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) expected[k] += prob[d] * xg[d * K + k];
    }
    for (size_t k = 0; k < K; ++k) gradient[k] += w * (xg[obs * K + k] - expected[k]);
    for (size_t a = 0; a < K; ++a) {
      for (size_t b = 0; b < K; ++b) {
        double second = 0;
        for (int d = 0; d < D; ++d) second += prob[d] * xg[d * K + a] * xg[d * K + b];
        hessian[a * K + b] -= w * (second - expected[a] * expected[b]);
      }
    }
  }

  if (!derivatives) {
    return Rcpp::List::create(Rcpp::Named("logPL") = (double)logPL);
  }
  Rcpp::NumericMatrix H(K, K);
  for (size_t a = 0; a < K; ++a) {
    for (size_t b = 0; b < K; ++b) H(a, b) = hessian[a * K + b];
  }
  return Rcpp::List::create(
    Rcpp::Named("logPL") = (double)logPL,
    Rcpp::Named("gradient") = Rcpp::NumericVector(gradient.begin(), gradient.end()),
    Rcpp::Named("hessian") = H);
}
