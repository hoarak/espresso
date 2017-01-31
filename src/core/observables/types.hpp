#ifndef OBSERVABLES_TYPES_HPP
#define OBSERVABLES_TYPES_HPP

#include <array>

namespace Observables {
namespace ParticleProperties {
using std::size_t;
template <size_t dim> using extent_t = std::array<size_t, dim>;

template <size_t N, size_t M> struct Matrix {
  constexpr static size_t size() { return N * M; }
  constexpr static size_t rank() { return 2; }
  constexpr static extent_t<2> extent() { return {{N, M}}; }
};

template <size_t N> struct Vector {
  constexpr static size_t size() { return N; }
  constexpr static size_t rank() { return 1; }
  constexpr static extent_t<1> extent() { return {{N}}; }
};

struct Scalar {
  constexpr static size_t size() { return 1; }
  constexpr static size_t rank() { return 0; }
  constexpr static extent_t<0> extent() { return {}; };
};
}
}
#endif
