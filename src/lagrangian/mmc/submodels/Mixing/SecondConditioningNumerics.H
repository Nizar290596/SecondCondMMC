/* SPDX-License-Identifier: GPL-3.0-or-later
 * Numerical kernels shared by the solver and standalone regression tests.
 * tauMix is the pair-contrast relaxation time: variance decays at 2/tauMix.
 */
#ifndef SecondConditioningNumerics_H
#define SecondConditioningNumerics_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>

namespace Foam
{
namespace secondConditioningNumerics
{
inline double progressReaction
(
    double phi, double dt, double A, double Z, double maxSourceStep
)
{
    if (!(phi >= 0 && phi <= 1 && dt >= 0 && A >= 0 && Z >= 0
          && maxSourceStep > 0 && maxSourceStep <= 1)
        || !std::isfinite(phi + dt + A + Z + maxSourceStep))
        throw std::domain_error("Invalid second-conditioning reaction input");
    if (A == 0 || dt == 0 || phi == 1) return phi;
    // Exponential midpoint applied to y=1-phi. Each substep is bounded
    // without clipping and the method is second order as dt -> 0.
    const double count = std::max(1.0, std::ceil(A*dt/maxSourceStep));
    if (count > 100000)
        throw std::domain_error("Progress reaction requires >100000 substeps; reduce deltaT");
    const unsigned n = static_cast<unsigned>(count);
    const double h = dt/n;
    for (unsigned i=0; i<n; ++i)
    {
        const double k = A*std::exp(Z*(phi-1));
        const double mid = phi + (1-phi)*(-std::expm1(-0.5*h*k));
        phi += (1-phi)*(-std::expm1(-h*A*std::exp(Z*(mid-1))));
    }
    return phi;
}

inline double ouUpdate(double old, double dt, double tau, double normal)
{
    if (!(dt >= 0 && tau > 0))
        throw std::domain_error("OU requires dt >= 0 and tau > 0");
    return old*std::exp(-dt/tau) + std::sqrt(-std::expm1(-2*dt/tau))*normal;
}

// A pair-local stream, seeded from random values communicated with the
// Eulerian particle data. Every rank holding the pair obtains the same draw.
class PairRandom
{
    std::uint64_t state_;
public:
    PairRandom(double a, double b)
    : state_((std::uint64_t(a*4294967296.0) << 32)
             ^ std::uint64_t(b*4294967296.0)) {}
    double uniform()
    {
        std::uint64_t z = (state_ += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27))*UINT64_C(0x94d049bb133111eb);
        z ^= z >> 31;
        return (double(z >> 12) + 0.5)/4503599627370496.0;
    }
};

inline double mixingExtent
(
    double dt, double tau, bool modifiedCurl, double randomP, double randomQ
)
{
    if (!(dt >= 0 && tau > 0))
        throw std::domain_error("Mixing requires dt >= 0 and tau > 0");
    if (!modifiedCurl) return -std::expm1(-dt/tau);
    // Uniform-extent Curl events: E[(1-alpha)^2]=1/3. A Poisson rate
    // 3/tau therefore gives E[contrast^2]=exp(-2*dt/tau), at any dt.
    PairRandom random(randomP, randomQ);
    const double exposure = 3*dt/tau;
    double eventTime = 0, contrast = 1;
    while ((eventTime -= std::log(random.uniform())) < exposure)
    {
        contrast *= random.uniform();
        if (contrast < std::numeric_limits<double>::epsilon()) return 1;
    }
    return 1-contrast;
}

// A rotating greedy matching on the processor adjacency graph. Groups are
// disjoint, so no particle can be mixed independently by two rank groups.
inline std::vector<int> neighbourPartners
(
    int ranks, const std::vector<std::pair<int,int>>& edges, unsigned phase
)
{
    std::vector<int> partner(ranks, -1);
    if (edges.empty()) return partner;
    for (std::size_t i=0; i<edges.size(); ++i)
    {
        const auto edge = edges[(i+phase%edges.size())%edges.size()];
        if (edge.first < 0 || edge.second < 0 || edge.first >= ranks
            || edge.second >= ranks || edge.first == edge.second)
            throw std::domain_error("Invalid processor adjacency edge");
        if (partner[edge.first] < 0 && partner[edge.second] < 0)
        {
            partner[edge.first] = edge.second;
            partner[edge.second] = edge.first;
        }
    }
    return partner;
}

inline double modifiedProgress(double phi, double beta, double omega, double ageTerm=0)
{
    const double value = (phi + ageTerm)*std::exp(beta*omega);
    if (!(value >= 0) || !std::isfinite(value))
        throw std::domain_error("Non-finite modified progress variable");
    return value;
}
}
}
#endif
